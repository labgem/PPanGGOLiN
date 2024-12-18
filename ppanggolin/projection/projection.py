#!/usr/bin/env python3

# default libraries
import argparse
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context, Value
import logging
import os
import time
from pathlib import Path
import tempfile
from typing import Tuple, Set, Dict, Optional, List, Iterable, Any
from collections import defaultdict
import csv
from itertools import chain

# installed libraries
from tqdm import tqdm
import networkx as nx
import yaml

# # local libraries
from ppanggolin.annotate.synta import get_contigs_from_fasta_file, get_dna_sequence
from ppanggolin.annotate.annotate import (
    init_contig_counter,
    read_anno_file,
    annotate_organism,
    local_identifiers_are_unique,
)
from ppanggolin.annotate import subparser as annotate_subparser
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import (
    detect_filetype,
    create_tmpdir,
    read_compressed_or_not,
    write_compressed_or_not,
    restricted_float,
    mk_outdir,
    get_config_args,
    parse_config_file,
    get_default_args,
    check_input_files,
    parse_input_paths_file,
)
from ppanggolin.align.alignOnPang import (
    write_gene_to_gene_family,
    get_input_seq_to_family_with_rep,
    get_input_seq_to_family_with_all,
    project_and_write_partition,
)
from ppanggolin.formats.writeSequences import write_gene_sequences_from_annotations
from ppanggolin.formats.readBinaries import check_pangenome_info
from ppanggolin.RGP.genomicIsland import naming_scheme, compute_org_rgp
from ppanggolin.RGP.spot import (
    make_spot_graph,
    check_sim,
    add_new_node_in_spot_graph,
    write_spot_graph,
)
from ppanggolin.genome import Organism
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.region import Region, Spot, Module
from ppanggolin.formats.writeFlatGenomes import (
    write_proksee_organism,
    manage_module_colors,
    write_gff_file,
    write_tsv_genome_file,
)
from ppanggolin.formats.writeFlatPangenome import (
    summarize_spots,
    summarize_genome,
    write_summaries_in_tsv,
    write_rgp_table,
)
from ppanggolin.formats.writeSequences import read_genome_file


class NewSpot(Spot):
    """
    This class represent a hotspot specifically
    created for the projected genome.
    """

    def __str__(self):
        return f"new_spot_{str(self.ID)}"


def check_pangenome_for_projection(pangenome: Pangenome, fast_aln: bool):
    """
    Check the status of a pangenome and determine whether projection is possible.

    :param pangenome: The pangenome to be checked.
    :param fast_aln: Whether to use the fast alignment option for gene projection.

    This function checks various attributes of a pangenome to determine whether it is suitable for projecting
    features into a provided genome.

    Returns:
        A tuple indicating whether RGP prediction, spot projection, and module projection
        are possible (True) or not (False) based on the pangenome's status.

    Raises:
        NameError: If the pangenome has not been partitioned.
        Exception: If the pangenome lacks gene sequences or gene family sequences, and fast alignment is not enabled.
    """

    project_modules = True
    predict_rgp = True
    project_spots = True

    if pangenome.status["partitioned"] not in ["Computed", "Loaded", "inFile"]:
        raise NameError(
            "The provided pangenome has not been partitioned. "
            "Annotation of an external genome is therefore not possible. "
            "See the 'partition' subcommands."
        )

    if pangenome.status["predictedRGP"] not in ["Computed", "Loaded", "inFile"]:
        logging.getLogger("PPanGGOLiN").info(
            "RGPs have not been predicted in the provided pangenome. "
            "Projection of RGPs and spots into the provided "
            "genome will not be performed."
        )
        predict_rgp = False
        project_spots = False

    elif pangenome.status["spots"] not in ["Computed", "Loaded", "inFile"]:
        logging.getLogger("PPanGGOLiN").info(
            "Spots have not been predicted in the provided pangenome. "
            "Projection of spots into the provided genome will not be performed."
        )
        project_spots = False

    if pangenome.status["modules"] not in ["Computed", "Loaded", "inFile"]:
        logging.getLogger("PPanGGOLiN").info(
            "Modules have not been predicted in the provided pangenome. "
            "Projection of modules into the provided genome will not be performed."
        )

        project_modules = False

    if (
        pangenome.status["geneSequences"] not in ["Loaded", "Computed", "inFile"]
        and not fast_aln
    ):
        raise Exception(
            "The provided pangenome has no gene sequences. "
            "Projection is still possible with the --fast option to use representative "
            "sequences rather than all genes to annotate input genes."
        )

    if pangenome.status["geneFamilySequences"] not in ["Loaded", "Computed", "inFile"]:
        raise Exception(
            "The provided pangenome has no gene families sequences. "
            "This is not possible to annotate an input genome to this pangenome."
        )

    return predict_rgp, project_spots, project_modules


def manage_input_genomes_annotation(
    pangenome,
    input_mode: str,
    anno: str,
    fasta: str,
    organism_name: str,
    circular_contigs: list,
    pangenome_params,
    cpu: int,
    use_pseudo: bool,
    disable_bar: bool,
    tmpdir: str,
    config: dict,
):
    """
    Manage the input genomes annotation based on the provided mode and parameters.

    :param pangenome: The pangenome object.
    :param input_mode: The input mode, either 'multiple' or 'single'.
    :param anno: The annotation file path or None.
    :param fasta: The FASTA file path or None.
    :param organism_name: The name of the organism.
    :param circular_contigs: List of circular contigs.
    :param pangenome_params: Parameters for pangenome processing.
    :param cpu: Number of CPUs to use.
    :param use_pseudo: Flag to use pseudo annotation.
    :param disable_bar: Flag to disable progress bar.
    :param tmpdir: Temporary directory path.
    :param config: Configuration dictionary.
    :return: A tuple of organisms, genome_name_to_path, and input_type.
    """

    genome_name_to_path = None
    input_type = None

    # Determine input type and parse paths based on the input mode
    if input_mode == "multiple":
        if anno:
            input_type = "annotation"
            genome_name_to_path = parse_input_paths_file(anno)
        elif fasta:
            input_type = "fasta"
            genome_name_to_path = parse_input_paths_file(fasta)
    elif input_mode == "single":
        circular_contigs = circular_contigs if circular_contigs else []
        if anno:
            input_type = "annotation"
            genome_name_to_path = {
                organism_name: {"path": anno, "circular_contigs": circular_contigs}
            }
        elif fasta:
            input_type = "fasta"
            genome_name_to_path = {
                organism_name: {"path": fasta, "circular_contigs": circular_contigs}
            }
    else:
        raise ValueError(
            f"Input mode '{input_mode}' is not valid. Expected 'multiple' or 'single'."
        )

    # Process annotation input type
    if input_type == "annotation":
        check_input_names(pangenome, genome_name_to_path)

        organisms, org_2_has_fasta = read_annotation_files(
            genome_name_to_path,
            cpu=cpu,
            pseudo=use_pseudo,
            translation_table=int(pangenome_params.cluster.translation_table),
            disable_bar=disable_bar,
        )

        # Check for genomes without associated sequence data
        if not all(has_fasta for has_fasta in org_2_has_fasta.values()):
            organisms_with_no_fasta = {
                org for org, has_fasta in org_2_has_fasta.items() if not has_fasta
            }
            if fasta:
                if input_mode == "multiple":
                    genome_name_to_fasta_path = parse_input_paths_file(fasta)
                else:
                    genome_name_to_fasta_path = {
                        organism_name: {
                            "path": fasta,
                            "circular_contigs": circular_contigs,
                        }
                    }
                get_gene_sequences_from_fasta_files(
                    organisms_with_no_fasta, genome_name_to_fasta_path
                )
            else:
                raise ValueError(
                    f"GFF files provided for {len(organisms_with_no_fasta)} (out of {len(organisms)}) genomes without "
                    "associated sequence data, and no FASTA sequences provided using the --fasta option. Cannot project "
                    "the pangenome onto these genomes. Genomes without sequence data: "
                    f"{', '.join(o.name for o in organisms_with_no_fasta)}"
                )

    # Process fasta input type
    elif input_type == "fasta":
        annotate_param_names = [
            "norna",
            "kingdom",
            "allow_overlap",
            "prodigal_procedure",
        ]
        annotate_params = manage_annotate_param(
            annotate_param_names, pangenome_params.annotate, config
        )

        check_input_names(pangenome, genome_name_to_path)
        organisms = annotate_fasta_files(
            genome_name_to_fasta_path=genome_name_to_path,
            tmpdir=tmpdir,
            cpu=cpu,
            translation_table=int(pangenome_params.cluster.translation_table),
            norna=annotate_params.norna,
            kingdom=annotate_params.kingdom,
            allow_overlap=annotate_params.allow_overlap,
            procedure=annotate_params.prodigal_procedure,
            disable_bar=disable_bar,
        )
    else:
        raise ValueError(
            f"Input type '{input_type}' is not valid. Expected 'fasta' or 'annotation'."
        )

    return organisms, genome_name_to_path, input_type


def write_projection_results(
    pangenome: Pangenome,
    organisms: Set[Organism],
    input_org_2_rgps: Dict[Organism, Set[Region]],
    input_org_to_spots: Dict[Organism, Set[Spot]],
    input_orgs_to_modules: Dict[Organism, Set[Module]],
    input_org_to_lonely_genes_count: Dict[Organism, int],
    write_proksee: bool,
    write_gff: bool,
    write_table: bool,
    add_sequences: bool,
    genome_name_to_path: Dict[str, dict],
    input_type: str,
    output_dir: Path,
    dup_margin: float,
    soft_core: float,
    metadata_sep: str,
    compress: bool,
    need_regions: bool,
    need_spots: bool,
    need_modules: bool,
):
    """
    Write the results of the projection of pangneome onto input genomes.

    :param pangenome: The pangenome onto which the projection is performed.
    :param organisms: A set of input organisms for projection.
    :param input_org_2_rgps: A dictionary mapping input organisms to sets of regions of genomic plasticity (RGPs).
    :param input_org_to_spots: A dictionary mapping input organisms to sets of spots.
    :param input_orgs_to_modules: A dictionary mapping input organisms to sets of modules.
    :param input_org_to_lonely_genes_count: A dictionary mapping input organisms to the count of lonely genes.
    :param write_proksee: Whether to write ProkSee JSON files.
    :param write_gff: Whether to write GFF files.
    :parama write_table: Whether to write table files.
    :param add_sequences: Whether to add sequences to the output files.
    :param genome_name_to_path: A dictionary mapping genome names to file paths.
    :param input_type: The type of input data (e.g., "annotation").
    :param output_dir: The directory where the output files will be written.
    :param dup_margin: The duplication margin used to compute completeness.
    :param soft_core: Soft core threshold


    Note:

    - If `write_proksee` is True and input organisms have modules, module colors for ProkSee are obtained.
    - The function calls other functions such as `summarize_projection`, `read_genome_file`, `write_proksee_organism`,
      `write_gff_file`, and `write_summaries` to generate various output files and summaries.
    """

    if write_proksee and input_orgs_to_modules:
        # get module color for proksee
        module_to_colors = manage_module_colors(set(pangenome.modules))

    # single_copy_families = get_single_copy_families(pangenome, dup_margin)
    # multigenics = pangenome.get_multigenics(pangenome_params.rgp.dup_margin)

    # dup margin value here is specified in argument and is used to compute completeness.
    # That means it can be different than dup margin used in spot and RGPS.

    pangenome_persistent_single_copy_families = (
        pangenome.get_single_copy_persistent_families(
            dup_margin=dup_margin, exclude_fragments=True
        )
    )
    pangenome_persistent_count = len(
        [fam for fam in pangenome.gene_families if fam.named_partition == "persistent"]
    )

    soft_core_families = pangenome.soft_core_families(soft_core)
    exact_core_families = pangenome.exact_core_families()

    summaries = []

    for organism in organisms:

        # summarize projection for all input organisms
        singleton_gene_count = input_org_to_lonely_genes_count[organism]

        org_summary = summarize_projected_genome(
            organism,
            pangenome_persistent_count,
            pangenome_persistent_single_copy_families,
            soft_core_families=soft_core_families,
            exact_core_families=exact_core_families,
            input_org_rgps=input_org_2_rgps.get(organism, None),
            input_org_spots=input_org_to_spots.get(organism, None),
            input_org_modules=input_orgs_to_modules.get(organism, None),
            pangenome_file=pangenome.file,
            singleton_gene_count=singleton_gene_count,
        )
        summaries.append(org_summary)

        yaml_outputfile = output_dir / organism.name / "projection_summary.yaml"
        write_summary_in_yaml(org_summary, yaml_outputfile)

        if (write_proksee or write_gff) and add_sequences:
            genome_sequences = read_genome_file(
                genome_name_to_path[organism.name]["path"], organism
            )
            genome_name_to_path[organism.name]["path"]
        else:
            genome_sequences = None

        if write_proksee:
            org_module_to_color = {
                org_mod: module_to_colors[org_mod]
                for org_mod in input_orgs_to_modules.get(organism, [])
            }

            output_file = output_dir / organism.name / f"{organism.name}_proksee.json"

            write_proksee_organism(
                organism,
                output_file,
                features="all",
                module_to_colors=org_module_to_color,
                genome_sequences=genome_sequences,
                compress=compress,
            )

        if write_gff:
            if (
                input_type == "annotation"
            ):  # if the genome has not been annotated by PPanGGOLiN
                annotation_sources = {
                    "rRNA": "external",
                    "tRNA": "external",
                    "CDS": "external",
                }
            else:
                annotation_sources = {}

            write_gff_file(
                organism,
                output_dir / organism.name,
                annotation_sources=annotation_sources,
                genome_sequences=genome_sequences,
                metadata_sep=metadata_sep,
                compress=compress,
            )

        if write_table:
            write_tsv_genome_file(
                organism,
                output_dir / organism.name,
                compress=compress,
                metadata_sep=metadata_sep,
                need_regions=need_regions,
                need_spots=need_spots,
                need_modules=need_modules,
            )

    output_file = output_dir / "summary_projection.tsv"
    write_summaries_in_tsv(
        summaries,
        output_file=output_file,
        dup_margin=dup_margin,
        soft_core=soft_core,
        compress=compress,
    )


def summarize_projected_genome(
    organism: Organism,
    pangenome_persistent_count: int,
    pangenome_persistent_single_copy_families: Set[GeneFamily],
    soft_core_families: Set[GeneFamily],
    exact_core_families: Set[GeneFamily],
    input_org_rgps: List[Region],
    input_org_spots: List[Spot],
    input_org_modules: List[Module],
    pangenome_file: str,
    singleton_gene_count: int,
) -> Dict[str, any]:
    """
    Summarizes the projected genome and generates an organism summary.

    :param organism: The Organism object for which the summary is generated.
    :param pangenome_persistent_count: The count of persistent genes in the pangenome.
    :param pangenome_persistent_single_copy_families: Set of single-copy persistent gene families.
    :param soft_core_families: Set of soft core families in the pangenome.
    :param exact_core_families: Set of exact core families in the pangenome.
    :param input_org_rgps: List of Region objects for the input organism.
    :param input_org_spots: List of Spot objects for the input organism.
    :param input_org_modules: List of Module objects for the input organism.
    :param pangenome_file: Filepath to the pangenome file.
    :param singleton_gene_count: Count of singleton genes in the organism.

    :return: A dictionary containing summarized information about the organism.
    """

    rgp_count = None if input_org_rgps is None else len(input_org_rgps)
    spot_count = None if input_org_spots is None else len(input_org_spots)
    module_count = None if input_org_modules is None else len(input_org_modules)

    organism_summary = summarize_genome(
        organism=organism,
        pangenome_persistent_count=pangenome_persistent_count,
        pangenome_persistent_single_copy_families=pangenome_persistent_single_copy_families,
        soft_core_families=soft_core_families,
        exact_core_families=exact_core_families,
        rgp_count=rgp_count,
        spot_count=spot_count,
        module_count=module_count,
    )

    # Add specific values for the projected genome
    organism_summary["Pangenome_file"] = pangenome_file
    cloud_without_specific_fams = (
        organism_summary["Cloud"]["families"] - singleton_gene_count
    )
    organism_summary["Cloud"]["families"] = cloud_without_specific_fams
    organism_summary["Cloud"]["specific families"] = singleton_gene_count

    input_org_spots = input_org_spots
    new_spot_count = (
        "Not computed"
        if input_org_spots is None
        else sum(1 for spot in input_org_spots if isinstance(spot, NewSpot))
    )
    organism_summary["New_spots"] = new_spot_count

    return organism_summary


def annotate_fasta_files(
    genome_name_to_fasta_path: Dict[str, dict],
    tmpdir: str,
    cpu: int = 1,
    translation_table: int = 11,
    kingdom: str = "bacteria",
    norna: bool = False,
    allow_overlap: bool = False,
    procedure: str = None,
    disable_bar: bool = False,
):
    """
    Main function to annotate a pangenome

    :param genome_name_to_fasta_path:
    :param fasta_list: List of fasta file containing sequences that will be base of pangenome
    :param tmpdir: Path to temporary directory
    :param cpu: number of CPU cores to use
    :param translation_table: Translation table (genetic code) to use.
    :param kingdom: Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation.
    :param norna: Use to avoid annotating RNA features.
    :param allow_overlap: Use to not remove genes overlapping with RNA features
    :param procedure: prodigal procedure used
    :param disable_bar: Disable the progresse bar
    """

    organisms = []
    arguments = []  # Argument given to annotate organism in same order than prototype
    for org_name, org_info in genome_name_to_fasta_path.items():
        arguments.append(
            (
                org_name,
                org_info["path"],
                org_info["circular_contigs"],
                tmpdir,
                translation_table,
                norna,
                kingdom,
                allow_overlap,
                procedure,
            )
        )

    logging.getLogger("PPanGGOLiN").info(
        f"Annotating {len(arguments)} genomes using {cpu} cpus..."
    )
    contig_counter = Value("i", 0)
    with ProcessPoolExecutor(
        mp_context=get_context("fork"),
        max_workers=cpu,
        initializer=init_contig_counter,
        initargs=(contig_counter,),
    ) as executor:
        with tqdm(total=len(arguments), unit="file", disable=disable_bar) as progress:
            futures = []

            for fn_args in arguments:
                future = executor.submit(annotate_organism, *fn_args)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            organisms.extend(future.result() for future in futures)

    return organisms


def read_annotation_files(
    genome_name_to_annot_path: Dict[str, dict],
    cpu: int = 1,
    pseudo: bool = False,
    translation_table: int = 11,
    disable_bar: bool = False,
) -> Tuple[List[Organism], Dict[Organism, bool]]:
    """
    Read the annotation from GBFF file

    :param pangenome: pangenome object
    :param organisms_file: List of GBFF files for each organism
    :param cpu: number of CPU cores to use
    :param pseudo: allow to read pseudogène
    :param translation_table: Translation table (genetic code) to use when /transl_table is missing from CDS tags.
    :param disable_bar: Disable the progresse bar
    """

    args = []
    organisms = []

    # we assume there are gene sequences in the annotation files,
    # unless a gff file without fasta is met (which is the only case where sequences can be absent)
    org_to_has_fasta_flag = {}

    args = [
        (
            org_name,
            org_info["path"],
            org_info["circular_contigs"],
            pseudo,
            translation_table,
        )
        for org_name, org_info in genome_name_to_annot_path.items()
    ]

    contig_counter = Value("i", 0)
    with ProcessPoolExecutor(
        mp_context=get_context("fork"),
        max_workers=cpu,
        initializer=init_contig_counter,
        initargs=(contig_counter,),
    ) as executor:
        with tqdm(total=len(args), unit="file", disable=disable_bar) as progress:
            futures = []

            for fn_args in args:
                future = executor.submit(read_anno_file, *fn_args)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                org, has_fasta = future.result()
                organisms.append(org)
                org_to_has_fasta_flag[org] = has_fasta

    genes = [gene for org in organisms for gene in org.genes]

    if local_identifiers_are_unique(genes):
        for gene in genes:
            gene.ID = (
                gene.local_identifier
            )  # Erase ppanggolin generated gene ids and replace with local identifiers
            gene.local_identifier = (
                ""  # this is now useless, setting it to default value
            )

        logging.getLogger("PPanGGOLiN").info(
            "Gene identifiers used in the provided annotation files were unique, "
            "PPanGGOLiN will use them."
        )
    else:
        logging.getLogger("PPanGGOLiN").info(
            "Gene identifiers used in the provided annotation files were not unique, "
            "PPanGGOLiN will use self-generated identifiers."
        )
    return organisms, org_to_has_fasta_flag


def get_gene_sequences_from_fasta_files(organisms, genome_name_to_annot_path):
    """
    Get gene sequences from fasta path file

    :param organisms: input pangenome
    :param fasta_file: list of fasta file
    """

    org_names = {org.name for org in organisms}

    if org_names & set(genome_name_to_annot_path) != org_names:
        missing = len(org_names - set(genome_name_to_annot_path))
        raise ValueError(
            f"You did not provided fasta for all the genomes found in annotation file. "
            f"{missing} are missing (out of {len(organisms)}). Missing genomes: {','.join(missing)}"
        )

    for org in organisms:

        org_fasta_file = genome_name_to_annot_path[org.name]["path"]

        with read_compressed_or_not(org_fasta_file) as currFastaFile:
            org_contig_to_seq = get_contigs_from_fasta_file(org, currFastaFile)

        for contig in org.contigs:
            try:
                contig_seq = org_contig_to_seq[contig.name]
            except KeyError:
                msg = (
                    f"Fasta file for genome {org.name} did not have the contig {contig.name} "
                    f"that was read from the annotation file. "
                )
                msg += (
                    f"The provided contigs in the fasta were : "
                    f"{', '.join(org_contig_to_seq)}."
                )
                raise KeyError(msg)

            for gene in contig.genes:
                gene.add_sequence(get_dna_sequence(contig_seq, gene))

            for rna in contig.RNAs:
                rna.add_sequence(get_dna_sequence(contig_seq, rna))


def check_input_names(pangenome, input_names):
    """
    Check if input organism names already exist in the pangenome.

    :param pangenome: The pangenome object.
    :param input_names: List of input organism names to check.
    :raises NameError: If duplicate organism names are found in the pangenome.
    """
    duplicated_names = set(input_names) & {org.name for org in pangenome.organisms}
    if len(duplicated_names) != 0:
        raise NameError(
            f"{len(duplicated_names)} provided genome name(s) already exist in the given pangenome: {' '.join(duplicated_names)}"
        )


def write_summary_in_yaml(summary_info: Dict[str, Any], output_file: Path):
    """
    Write summary information to a YAML file.

    This function takes a dictionary containing summary
    information. It writes this information to a YAML file.

    :param summary_info: A dictionary containing summary information.
    :param output_file: The file where the summary will be written.
    """

    yaml_string = yaml.dump(
        summary_info, default_flow_style=False, sort_keys=False, indent=4
    )

    with open(output_file, "w") as flout:
        flout.write("Projection_summary:")
        flout.write(yaml_string)


def annotate_input_genes_with_pangenome_families(
    pangenome: Pangenome,
    input_organisms: Iterable[Organism],
    output: Path,
    cpu: int,
    use_representatives: bool,
    no_defrag: bool,
    identity: float,
    coverage: float,
    tmpdir: Path,
    translation_table: int,
    keep_tmp: bool = False,
    disable_bar: bool = False,
):
    """
    Annotate input genes with pangenome gene families by associating them to a cluster.

    :param pangenome: Pangenome object.
    :param input_organisms: Iterable of input organism objects.
    :param output: Output directory for generated files.
    :param cpu: Number of CPU cores to use.
    :param no_defrag: Whether to use defragmentation.
    :param use_representatives: Use representative sequences of gene families rather than all sequence to align input genes
    :param identity: Minimum identity threshold for gene clustering.
    :param coverage: Minimum coverage threshold for gene clustering.
    :param tmpdir: Temporary directory for intermediate files.
    :param translation_table: Translation table ID for nucleotide sequences.
    :param keep_tmp: If True, keep temporary files.
    :param disable_bar: Whether to disable progress bar.

    :return: Number of genes that do not cluster with any of the gene families of the pangenome.
    """
    logging.getLogger("PPanGGOLiN").info("Writing gene sequences of input genomes.")

    input_genes = [gene for org in input_organisms for gene in org.genes]

    seq_fasta_file = output / "input_genes.fasta"

    write_gene_sequences_from_annotations(
        input_genes, seq_fasta_file, disable_bar=True, add="ppanggolin_"
    )

    with create_tmpdir(
        main_dir=tmpdir, basename="projection_tmp", keep_tmp=keep_tmp
    ) as new_tmpdir:
        if use_representatives:
            _, seqid_to_gene_family = get_input_seq_to_family_with_rep(
                pangenome=pangenome,
                sequence_files=seq_fasta_file,
                output=new_tmpdir,
                tmpdir=new_tmpdir,
                input_type="nucleotide",
                is_input_slf=True,
                cpu=cpu,
                no_defrag=no_defrag,
                identity=identity,
                coverage=coverage,
                translation_table=translation_table,
                disable_bar=disable_bar,
            )
        else:
            _, seqid_to_gene_family = get_input_seq_to_family_with_all(
                pangenome=pangenome,
                sequence_files=seq_fasta_file,
                output=new_tmpdir,
                tmpdir=new_tmpdir,
                input_type="nucleotide",
                is_input_slf=True,
                cpu=cpu,
                no_defrag=no_defrag,
                identity=identity,
                coverage=coverage,
                translation_table=translation_table,
                disable_bar=disable_bar,
            )

    input_org_to_lonely_genes_count = {}
    for input_organism in input_organisms:
        org_outdir = output / input_organism.name
        mk_outdir(org_outdir, force=True)

        seq_set = {
            gene.ID if gene.local_identifier == "" else gene.local_identifier
            for gene in input_organism.genes
        }

        project_and_write_partition(seqid_to_gene_family, seq_set, org_outdir)

        write_gene_to_gene_family(seqid_to_gene_family, seq_set, org_outdir)

        lonely_genes = set()
        for gene in input_organism.genes:
            gene_id = gene.ID
            try:
                gene_family = seqid_to_gene_family[gene_id]
            except KeyError:
                # the seqid is not in the dict so it does not align with any pangenome families
                # We consider it as cloud gene
                try:
                    # in some case a family exists already and has the same name of the gene id
                    # So gene id cannot be used
                    _ = pangenome.get_gene_family(gene_id)
                except KeyError:
                    gene_family = GeneFamily(pangenome.max_fam_id, gene_id)

                else:
                    # gene id already exists.
                    new_name = f"{input_organism.name}_{gene_id}"
                    logging.getLogger("PPanGGOLiN").warning(
                        "The input genome as a specific gene that does not align to any "
                        f"pangenome families with the same id ({gene_id}) than an existing gene family in the pangenome. "
                        f"The genome name is added to the family name: {new_name}"
                    )
                    gene_family = GeneFamily(pangenome.max_fam_id, new_name)

                pangenome.add_gene_family(gene_family)

                gene_family.partition = "Cloud"
                lonely_genes.add(gene)

            if gene_family.contains_gene_id(gene_id):
                new_name = f"{input_organism.name}_{gene_id}"
                logging.getLogger("PPanGGOLiN").warning(
                    "The input genome contains a gene that aligns to a pangenome family "
                    f"which already contains a gene with the same ID ({gene_id}). "
                    f"The genome name has been appended to the family name: {new_name}"
                )

                gene.ID = new_name

            # Add the gene to the gene family
            gene_family.add(gene)

        pangenome._mk_gene_getter()  # re-build the gene getter

        logging.getLogger("PPanGGOLiN").info(
            f"{input_organism.name} has {len(lonely_genes)}/{input_organism.number_of_genes()} "
            "specific genes that do not align to any gene of the pangenome."
        )
        # Write specific gene ids in a file
        with open(org_outdir / "specific_genes.tsv", "w") as fl:
            fl.write(
                "\n".join(
                    gene.ID if gene.local_identifier == "" else gene.local_identifier
                    for gene in lonely_genes
                )
                + "\n"
            )

        input_org_to_lonely_genes_count[input_organism] = len(lonely_genes)

    return input_org_to_lonely_genes_count


def predict_RGP(
    pangenome: Pangenome,
    input_organisms: List[Organism],
    persistent_penalty: int,
    variable_gain: int,
    min_length: int,
    min_score: int,
    multigenics: Set[GeneFamily],
    output_dir: Path,
    disable_bar: bool,
    compress: bool,
) -> Dict[Organism, Set[Region]]:
    """
    Compute Regions of Genomic Plasticity (RGP) for the given input organisms.

    :param pangenome: The pangenome object.
    :param input_organisms: List of the input organisms for which to compute RGPs.
    :param persistent_penalty: Penalty score to apply to persistent genes.
    :param variable_gain: Gain score to apply to variable genes.
    :param min_length: Minimum length (bp) of a region to be considered as RGP.
    :param min_score: Minimal score required for considering a region as RGP.
    :param multigenics: multigenic families.
    :param output_dir: Output directory where predicted rgps are going to be written.
    :param disable_bar: Flag to disable the progress bar.
    :param compress: Flag to compress the rgp table in gz.

    :return: Dictionary mapping organism with the set of predicted regions
    """

    logging.getLogger("PPanGGOLiN").info("Computing Regions of Genomic Plasticity...")

    name_scheme = naming_scheme(chain(pangenome.organisms, input_organisms))
    organism_to_rgps = {}

    for input_organism in input_organisms:
        rgps = compute_org_rgp(
            input_organism,
            multigenics,
            persistent_penalty,
            variable_gain,
            min_length,
            min_score,
            naming=name_scheme,
            disable_bar=disable_bar,
        )
        # turn on projected attribute in rgp objects
        # useful when associating spot to prevent failure when multiple spot are associated to a projected RGP
        for rgp in rgps:
            rgp.projected = True

        logging.getLogger("PPanGGOLiN").info(
            f"{len(rgps)} RGPs have been predicted in the input genomes."
        )

        org_outdir = output_dir / input_organism.name

        write_rgp_table(rgps, output=org_outdir, compress=compress)
        organism_to_rgps[input_organism] = rgps

    return organism_to_rgps


def write_rgp_to_spot_table(
    rgp_to_spots: Dict[Region, Set[str]],
    output: Path,
    filename: str,
    compress: bool = False,
):
    """
    Write a table mapping RGPs to corresponding spot IDs.

    :param rgp_to_spots: A dictionary mapping RGPs to spot IDs.
    :param output: Path to the output directory.
    :param filename: Name of the file to write.
    :param compress: Whether to compress the file.
    """
    fname = output / filename
    logging.getLogger("PPanGGOLiN").debug(f"Writing RGPs to spot table in {fname}")

    with write_compressed_or_not(fname, compress) as tab:
        fieldnames = ["region", "spot_id"]

        writer = csv.DictWriter(tab, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        regions = sorted(
            rgp_to_spots.keys(), key=lambda x: (x.organism.name, x.contig.name, x.ID)
        )
        for region in regions:
            row = {
                "region": region.name,
                "spot_id": ";".join(map(str, rgp_to_spots[region])),
            }

            writer.writerow(row)


def retrieve_gene_sequences_from_fasta_file(input_organism, fasta_file):
    """
    Get gene sequences from fastas

    :param pangenome: input pangenome
    :param fasta_file: list of fasta file
    """

    with read_compressed_or_not(fasta_file) as currFastaFile:
        contig_id2seq = get_contigs_from_fasta_file(input_organism, currFastaFile)

    for contig in input_organism.contigs:
        try:
            for gene in contig.genes:
                gene.add_dna(get_dna_sequence(contig_id2seq[contig.name], gene))

            for rna in contig.RNAs:
                rna.add_dna(get_dna_sequence(contig_id2seq[contig.name], rna))
        except KeyError:
            msg = (
                f"Fasta file for input genome {input_organism.name} did not have the contig {contig.name} "
                f"that was read from the annotation file. "
            )
            msg += (
                f"The provided contigs in the fasta were : "
                f"{', '.join(contig_id2seq.keys())}."
            )
            raise KeyError(msg)


def manage_annotate_param(
    annotate_param_names: List[str],
    pangenome_args: argparse.Namespace,
    config_file: Optional[str],
) -> argparse.Namespace:
    """
    Manage annotate parameters by collecting them from different sources and merging them.

    :param annotate_param_names: List of annotate parameter names to be managed.
    :param pangenome_args: Annotate arguments parsed from pangenomes parameters.
    :param config_file: Path to the config file, can be None if not provided.

    :return: An argparse.Namespace containing the merged annotate parameters with their values.
    """

    default_annotate_args = get_default_args("annotate", annotate_subparser)

    if config_file is None:
        config_annotate_args = argparse.Namespace()
    else:
        config = defaultdict(dict, parse_config_file(config_file))
        config_annotate_args = get_config_args(
            "annotate",
            annotate_subparser,
            config,
            "annotate",
            annotate_param_names,
            strict_config_check=False,
        )

    annotate_param_from_pangenome = {}
    annotate_param_from_config = {}
    annotate_param_from_default = {}

    annotate_params = argparse.Namespace()

    # Collecting annotate parameters from different sources
    # if they are found in pangenome param they are used
    # elif they are found in config they are used
    # else use the default value.
    for annotate_arg in annotate_param_names:
        if hasattr(pangenome_args, annotate_arg):
            param_val = getattr(pangenome_args, annotate_arg)
            annotate_param_from_pangenome[annotate_arg] = param_val
            setattr(annotate_params, annotate_arg, param_val)

        elif hasattr(config_annotate_args, annotate_arg):
            param_val = getattr(config_annotate_args, annotate_arg)
            annotate_param_from_config[annotate_arg] = param_val
            setattr(annotate_params, annotate_arg, param_val)

        else:
            param_val = getattr(default_annotate_args, annotate_arg)
            annotate_param_from_default[annotate_arg] = param_val
            setattr(annotate_params, annotate_arg, param_val)

    # Log the sources of the annotate parameters
    if len(annotate_param_from_pangenome) > 0:
        param_val_string = " ".join(
            [f"--{k} {v}" for k, v in annotate_param_from_pangenome.items()]
        )
        logging.getLogger("PPanGGOLiN").debug(
            f"{len(annotate_param_from_pangenome)}/{len(annotate_param_names)} annotate parameters extracted from pangenome parameters "
            f"(the parameters used to build the input pangenome): {param_val_string}"
        )

    if len(annotate_param_from_config) > 0:
        param_val_string = ";".join(
            [f" {k} : {v}" for k, v in annotate_param_from_config.items()]
        )
        logging.getLogger("PPanGGOLiN").debug(
            f"{len(annotate_param_from_config)}/{len(annotate_param_names)} annotate parameters were not found in pangenome internal parameters."
            f" They have been parsed from the annotate section in the config file: {param_val_string}"
        )

    if len(annotate_param_from_default) > 0:
        param_val_string = ";".join(
            [f" {k} : {v}" for k, v in annotate_param_from_default.items()]
        )
        logging.getLogger("PPanGGOLiN").debug(
            f"{len(annotate_param_from_default)}/{len(annotate_param_names)} annotate parameters were not found in the pangenome parameters "
            f"nor in the config file. Default values have been used: {param_val_string}"
        )

    return annotate_params


def check_spots_congruency(graph_spot: nx.Graph, spots: List[Spot]) -> None:
    """
    Check congruency of spots in the spot graph with the original spots.

    :param graph_spot: The spot graph containing the connected components representing the spots.
    :param spots: List of original spots in the pangenome.
    :return: None.
    """
    rgp_to_spot = {region: spot for spot in spots for region in spot.regions}

    spots = []
    for cc in nx.algorithms.components.connected_components(graph_spot):
        # one connected component is a spot
        regions_in_cc = set()
        for node in cc:
            regions_in_cc |= graph_spot.nodes[node]["rgp"]

        # check that region in cc are the regions of a spot
        spot_in_cc = {rgp_to_spot[rgp] for rgp in regions_in_cc}
        assert (
            len(spot_in_cc) == 1
        ), "More than one spot in a connected_components. Something went wrong when recomputing spots."
        current_spot = spot_in_cc.pop()
        # Add spot id to the graph
        for node in cc:
            graph_spot.nodes[node]["spot_id"] = str(current_spot)
            graph_spot.nodes[node]["spots"] = {current_spot}


def predict_spots_in_input_organisms(
    initial_spots: List[Spot],
    initial_regions: List[Region],
    input_org_2_rgps: Dict[Organism, Set[Region]],
    multigenics: Set[GeneFamily],
    output: Path,
    write_graph_flag: bool = False,
    graph_formats: List[str] = ["gexf"],
    overlapping_match: int = 2,
    set_size: int = 3,
    exact_match: int = 1,
    compress: bool = False,
) -> Dict[Organism, Set[Spot]]:
    """
    Create a spot graph from pangenome RGP and predict spots for input organism RGPs.

    :param initial_spots: List of original spots in the pangenome.
    :param initial_regions: List of original regions in the pangenome.
    :param input_org_2_rgps: Dictionary mapping input organisms to their RGPs.
    :param multigenics: Set of pangenome graph multigenic persistent families.
    :param output: Output directory to save the spot graph.
    :param write_graph_flag: If True, writes the spot graph in the specified formats. Default is False.
    :param graph_formats: List of graph formats to write (default is ['gexf']).
    :param overlapping_match: Number of missing persistent genes allowed when comparing flanking genes. Default is 2.
    :param set_size: Number of single copy markers to use as flanking genes for RGP during hotspot computation. Default is 3.
    :param exact_match: Number of perfectly matching flanking single copy markers required to associate RGPs. Default is 1.
    :param compress: Flag to compress output files

    :return: A dictionary mapping input organism RGPs to their predicted spots.
    """

    logging.getLogger("PPanGGOLiN").debug("Rebuilding original spot graph.")
    graph_spot = make_spot_graph(
        rgps=initial_regions,
        multigenics=multigenics,
        overlapping_match=overlapping_match,
        set_size=set_size,
        exact_match=exact_match,
    )

    original_nodes = set(graph_spot.nodes)

    # Check congruency with already computed spot and add spot id in node attributes
    check_spots_congruency(graph_spot, initial_spots)

    new_spot_id_counter = (
        max(s.ID for s in initial_spots) + 1 if len(initial_spots) != 0 else 1
    )

    input_org_to_spots = {}
    for input_organism, rgps in input_org_2_rgps.items():

        if len(rgps) == 0:
            logging.getLogger("PPanGGOLiN").debug(
                f"{input_organism.name}: No RGPs have been found. "
                "As a result, spot prediction and RGP output will be skipped."
            )

            input_org_to_spots[input_organism] = set()
            continue

        outdir_org = output / input_organism.name
        # Copy the graph spot, as each input organism are processed independently
        graph_spot_cp = graph_spot.copy()

        input_org_spots = predict_spot_in_one_organism(
            graph_spot_cp,
            input_org_rgps=rgps,
            original_nodes=original_nodes,
            new_spot_id_counter=new_spot_id_counter,
            multigenics=multigenics,
            organism_name=input_organism.name,
            output=outdir_org,
            write_graph_flag=write_graph_flag,
            graph_formats=graph_formats,
            overlapping_match=overlapping_match,
            set_size=set_size,
            exact_match=exact_match,
            compress=compress,
        )

        if len(input_org_spots) > 0:
            new_spot_id_counter = max(s.ID for s in input_org_spots) + 1

        input_org_to_spots[input_organism] = input_org_spots

    return input_org_to_spots


def predict_spot_in_one_organism(
    graph_spot: nx.Graph,
    input_org_rgps: List[Region],
    original_nodes: Set[int],
    new_spot_id_counter: int,
    multigenics: Set[GeneFamily],
    organism_name: str,
    output: Path,
    write_graph_flag: bool = False,
    graph_formats: List[str] = ["gexf"],
    overlapping_match: int = 2,
    set_size: int = 3,
    exact_match: int = 1,
    compress: bool = False,
) -> Set[Spot]:
    """
    Predict spots for input organism RGPs.

    :param graph_spot: The spot graph from the pangenome.
    :param input_org_rgps: List of RGPs from the input organism to be associated with spots.
    :param original_nodes: Set of original nodes in the spot graph.
    :param new_spot_id_counter: Counter for new spot IDs.
    :param multigenics: Set of pangenome graph multigenic persistent families.
    :param organism_name: Name of the input organism.
    :param output: Output directory to save the spot graph.
    :param write_graph_flag: If True, writes the spot graph in the specified formats. Default is False.
    :param graph_formats: List of graph formats to write (default is ['gexf']).
    :param overlapping_match: Number of missing persistent genes allowed when comparing flanking genes. Default is 2.
    :param set_size: Number of single copy markers to use as flanking genes for RGP during hotspot computation. Default is 3.
    :param exact_match: Number of perfectly matching flanking single copy markers required to associate RGPs. Default is 1.
    :param compress: Flag to compress output files

    Returns:
        Set[Spot]: The predicted spots for the input organism RGPs.
    """
    # Check which input RGP has a spot
    lost = 0
    used = 0

    input_org_node_to_rgps = defaultdict(set)

    for rgp in input_org_rgps:
        border = rgp.get_bordering_genes(set_size, multigenics)
        if len(border[0]) < set_size or len(border[1]) < set_size:
            lost += 1
        else:
            used += 1
            border_node = add_new_node_in_spot_graph(graph_spot, rgp, border)
            input_org_node_to_rgps[border_node].add(rgp)

    if len(input_org_node_to_rgps) == 0:
        logging.getLogger("PPanGGOLiN").debug(
            f"{organism_name}: no RGPs of the input genome will be associated with any spot of insertion "
            "as they are on a contig border (or have "
            f"less than {set_size} persistent gene families until the contig border). "
            "Projection of spots stops here"
        )
        return set()

    # remove node that were already in the graph
    new_nodes = set(input_org_node_to_rgps) - original_nodes

    logging.getLogger("PPanGGOLiN").debug(
        f"{organism_name}: {lost} RGPs were not used as they are on a contig border (or have"
        f"less than {set_size} persistent gene families until the contig border)"
    )

    logging.getLogger("PPanGGOLiN").debug(
        f"{organism_name}: {used} RGPs of the input genome will be associated to a spot of insertion"
    )

    # add potential edges from new nodes to the rest of the nodes
    all_nodes = list(graph_spot.nodes)
    for nodei in new_nodes:
        for nodej in all_nodes:
            if nodei == nodej:
                continue
            node_obj_i = graph_spot.nodes[nodei]
            node_obj_j = graph_spot.nodes[nodej]
            if check_sim(
                [node_obj_i["border0"], node_obj_i["border1"]],
                [node_obj_j["border0"], node_obj_j["border1"]],
                overlapping_match,
                set_size,
                exact_match,
            ):
                graph_spot.add_edge(nodei, nodej)

    input_rgp_to_spots = {}
    new_spots = []

    # determine spot ids of the new nodes and by extension to their rgps
    for comp in nx.algorithms.components.connected_components(graph_spot):
        # in very rare case one cc can have several original spots
        # that would mean a new nodes from the input organism have connected two old cc
        # in this case we report the two spots in the output
        spots_of_the_cc = set()
        for node in comp:
            if "spots" in graph_spot.nodes[node]:
                spots_of_the_cc |= set(graph_spot.nodes[node]["spots"])

        if len(spots_of_the_cc) == 0:
            # no spot associated with any node of the cc
            # that means this cc is only composed of new nodes
            # let's add a new spot id
            new_spot = NewSpot(new_spot_id_counter)
            new_spots.append(new_spot)
            spots_of_the_cc = {new_spot}  # {f"new_spot_{new_spot_id_counter}"}
            new_spot_id_counter += 1

        elif len(spots_of_the_cc) > 1:
            # more than one spot in the cc
            logging.getLogger("PPanGGOLiN").debug(
                f"{organism_name}: Some RGPs of the input genome "
                f"are connected to {len(spots_of_the_cc)} original spots of the pangenome."
            )

        input_rgps_of_the_cc = set()
        for node in comp:
            if node in input_org_node_to_rgps:
                input_rgps_of_the_cc |= input_org_node_to_rgps[node]

                if write_graph_flag:
                    graph_spot.nodes[node]["spots"] = spots_of_the_cc

                    graph_spot.nodes[node]["spot_id"] = ";".join(
                        str(spot) for spot in spots_of_the_cc
                    )
                    graph_spot.nodes[node]["includes_RGPs_from_the_input_genome"] = True

        for spot in spots_of_the_cc:
            for region in input_rgps_of_the_cc:
                spot.add(region)

        input_rgp_to_spots.update(
            {rgp: spots_of_the_cc for rgp in input_rgps_of_the_cc}
        )

    if write_graph_flag:
        # remove node that would not be writable in graph file
        for node in graph_spot.nodes:
            del graph_spot.nodes[node]["spots"]

        write_spot_graph(
            graph_spot, output, graph_formats, file_basename="projected_spotGraph"
        )

    write_rgp_to_spot_table(
        input_rgp_to_spots,
        output=output,
        filename="input_genome_rgp_to_spot.tsv",
        compress=compress,
    )

    input_org_spots = {spot for spots in input_rgp_to_spots.values() for spot in spots}
    new_spots = {spot for spot in input_org_spots if isinstance(spot, NewSpot)}

    logging.getLogger("PPanGGOLiN").debug(
        f"{organism_name}: {len(new_spots)} new spots have been created for the input genome."
    )

    if new_spots:
        summarize_spots(
            new_spots, output, compress=compress, file_name="new_spots_summary.tsv"
        )

    return input_org_spots


def project_and_write_modules(
    pangenome: Pangenome,
    input_organisms: Iterable[Organism],
    output: Path,
    compress: bool = False,
):
    """
    Write a tsv file providing association between modules and the input organism

    :param pangenome: Pangenome object
    :param input_organisms: iterable of the organisms that is being annotated
    :param output: Path to output directory
    :param compress: Compress the file in .gz
    """
    input_orgs_to_modules = {}
    for input_organism in input_organisms:
        output_file = output / input_organism.name / "modules_in_input_genome.tsv"

        input_organism_families = list(input_organism.families)
        counter = 0
        modules_in_input_org = []
        with write_compressed_or_not(output_file, compress) as fout:
            fout.write("module_id\tgenome\tcompletion\n")

            for mod in pangenome.modules:
                module_in_input_organism = any(
                    fam in input_organism_families for fam in mod.families
                )

                if module_in_input_organism:
                    counter += 1
                    modules_in_input_org.append(mod)

                    completion = round(
                        len(set(input_organism.families) & set(mod.families))
                        / len(set(mod.families)),
                        2,
                    )
                    fout.write(
                        f"module_{mod.ID}\t{input_organism.name}\t{completion}\n"
                    )

        logging.getLogger("PPanGGOLiN").debug(
            f"{input_organism.name}: {counter} modules have been projected to the input genomes."
        )

        logging.getLogger("PPanGGOLiN").debug(
            f"{input_organism.name}: Projected modules have been written in: '{output_file}'"
        )

        input_orgs_to_modules[input_organism] = modules_in_input_org

    return input_orgs_to_modules


def infer_input_mode(
    input_file: Path, expected_types: List[str], parser: argparse.ArgumentParser
) -> str:
    """
    Determine the input mode based on the provided input file and expected file types.

    :param input_file: A Path object representing the input file.
    :param expected_types: A list of expected file types (e.g., ['fasta', 'gff', 'gbff', 'tsv']).

    :return: A string indicating the input mode ('single' or 'multiple').
    """
    if not input_file.exists():
        parser.error(f"The provided file {input_file} does not exist.")

    try:
        filetype = detect_filetype(input_file)
    except Exception:
        parser.error(
            "Based on its content, the provided file is not recognized as a valid input file. Please ensure it is in one of the supported formats (FASTA, GFF/GBFF, or TSV)."
        )

    if filetype == "tsv":
        logging.getLogger("PPanGGOLiN").debug(
            f"The provided file ({input_file}) is detected as a TSV file."
        )
        mode = "multiple"
    elif filetype in expected_types:
        logging.getLogger("PPanGGOLiN").debug(
            f"The provided file ({input_file}) is detected as a single {'/'.join(expected_types)} file."
        )
        mode = "single"
    else:
        logging.getLogger("PPanGGOLiN").error(
            f"The provided file {input_file} is not recognized as a valid {'/'.join(expected_types)} file or a TSV file listing names and {'/'.join(expected_types)} files of genomes to annotate."
        )
        parser.error(
            f"The provided file {input_file} is not recognized as a valid {'/'.join(expected_types)} file or a TSV file listing names and files of genomes to annotate."
        )

    return mode


def check_projection_arguments(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> str:
    """
    Check the arguments provided for genome projection and raise errors if they are incompatible or missing.

    :param args: An argparse.Namespace object containing parsed command-line arguments.
    :param parser: parser of the command
    :return: A string indicating the input mode ('single' or 'multiple').
    """

    # Check if we annotate genomes from path files or only a single genome...
    if not args.anno and not args.fasta:
        parser.error(
            "Please provide either a FASTA file or a tab-separated file listing sequence files using the '--fasta' option, "
            "or an annotation file or a tab-separated file listing annotation files using the '--anno' option. "
            "You can specify these either through the command line or the configuration file."
        )

    mode_from_fasta, mode_from_anno = None, None
    if args.fasta:
        mode_from_fasta = infer_input_mode(args.fasta, ["fasta"], parser)
        input_mode = mode_from_fasta

    if args.anno:
        mode_from_anno = infer_input_mode(args.anno, ["gff", "gbff"], parser)
        input_mode = mode_from_anno

        logging.getLogger("PPanGGOLiN").debug("")

    if mode_from_fasta and mode_from_anno and mode_from_fasta != mode_from_anno:
        single_input, multiple_input = (
            ("fasta", "anno") if mode_from_fasta == "single" else ("anno", "fasta")
        )

        parser.error(
            f"You've provided both a single annotation/fasta file using the '--{single_input}' option and a list of files using "
            f"the '--{multiple_input}' option. Please choose either a single file or a tab-separated file listing genome files, but not both."
        )

    if input_mode == "multiple":
        # We are in paths file mode

        if args.circular_contigs:
            parser.error(
                "You provided a TSV file listing the files of genomes you wish to annotate. "
                "Therefore, the  argument '--circular_contigs' is incompatible with this multiple genomes file."
            )

        if args.fasta:
            check_input_files(args.fasta, True)

        if args.anno:
            check_input_files(args.anno, True)

    return input_mode


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """

    output_dir = Path(args.output)
    mk_outdir(output_dir, args.force)

    # For the moment these elements of the pangenome are predicted by default

    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)

    predict_rgp, project_spots, project_modules = check_pangenome_for_projection(
        pangenome, args.fast
    )

    need_graph = True if args.table else False

    check_pangenome_info(
        pangenome,
        need_annotations=True,
        need_families=True,
        disable_bar=args.disable_prog_bar,
        need_rgp=predict_rgp,
        need_modules=project_modules,
        need_gene_sequences=False,
        need_spots=project_spots,
        need_graph=need_graph,
        need_metadata=True,
    )

    logging.getLogger("PPanGGOLiN").info(
        "Retrieving parameters from the provided pangenome file."
    )
    pangenome_params = argparse.Namespace(
        **{
            step: argparse.Namespace(**k_v)
            for step, k_v in pangenome.parameters.items()
        }
    )

    if predict_rgp:
        # computing multigenics for rgp prediction first to have original family.number_of_genomes
        # and the same multigenics list as when rgp and spot were predicted
        multigenics = pangenome.get_multigenics(pangenome_params.rgp.dup_margin)

    organisms, genome_name_to_path, input_type = manage_input_genomes_annotation(
        pangenome=pangenome,
        input_mode=args.input_mode,
        anno=args.anno,
        fasta=args.fasta,
        organism_name=args.genome_name,
        circular_contigs=args.circular_contigs,
        pangenome_params=pangenome_params,
        cpu=args.cpu,
        use_pseudo=args.use_pseudo,
        disable_bar=args.disable_prog_bar,
        tmpdir=args.tmpdir,
        config=args.config,
    )

    input_org_to_lonely_genes_count = annotate_input_genes_with_pangenome_families(
        pangenome,
        input_organisms=organisms,
        output=output_dir,
        cpu=args.cpu,
        use_representatives=args.fast,
        no_defrag=args.no_defrag,
        identity=args.identity,
        coverage=args.coverage,
        tmpdir=args.tmpdir,
        translation_table=int(pangenome_params.cluster.translation_table),
        keep_tmp=args.keep_tmp,
        disable_bar=args.disable_prog_bar,
    )

    input_org_2_rgps, input_org_to_spots, input_orgs_to_modules = {}, {}, {}

    if predict_rgp:

        logging.getLogger("PPanGGOLiN").info("Detecting RGPs in input genomes.")

        input_org_2_rgps = predict_RGP(
            pangenome,
            organisms,
            persistent_penalty=pangenome_params.rgp.persistent_penalty,
            variable_gain=pangenome_params.rgp.variable_gain,
            min_length=pangenome_params.rgp.min_length,
            min_score=pangenome_params.rgp.min_score,
            multigenics=multigenics,
            output_dir=output_dir,
            disable_bar=args.disable_prog_bar,
            compress=args.compress,
        )

        if project_spots:
            logging.getLogger("PPanGGOLiN").info(
                "Predicting spot of insertion in input genomes."
            )
            input_org_to_spots = predict_spots_in_input_organisms(
                initial_spots=list(pangenome.spots),
                initial_regions=pangenome.regions,
                input_org_2_rgps=input_org_2_rgps,
                multigenics=multigenics,
                output=output_dir,
                write_graph_flag=args.spot_graph,
                graph_formats=args.graph_formats,
                overlapping_match=pangenome_params.spot.overlapping_match,
                set_size=pangenome_params.spot.set_size,
                exact_match=pangenome_params.spot.exact_match_size,
                compress=args.compress,
            )

    if project_modules:
        input_orgs_to_modules = project_and_write_modules(
            pangenome, organisms, output_dir, compress=args.compress
        )

    write_projection_results(
        pangenome,
        organisms,
        input_org_2_rgps,
        input_org_to_spots,
        input_orgs_to_modules,
        input_org_to_lonely_genes_count,
        write_proksee=args.proksee,
        write_gff=args.gff,
        write_table=args.table,
        add_sequences=args.add_sequences,
        genome_name_to_path=genome_name_to_path,
        input_type=input_type,
        output_dir=output_dir,
        dup_margin=args.dup_margin,
        soft_core=args.soft_core,
        metadata_sep=args.metadata_sep,
        compress=args.compress,
        need_modules=project_modules,
        need_spots=project_spots,
        need_regions=predict_rgp,
    )


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for projection command

    :return : parser arguments for projection command
    """
    parser = sub_parser.add_parser(
        "projection", formatter_class=argparse.RawTextHelpFormatter
    )
    parser_projection(parser)
    return parser


def parser_projection(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of projection command

    :param parser: parser for projection argument
    """
    required = parser.add_argument_group(title="Required arguments")

    required.add_argument(
        "-p", "--pangenome", required=False, type=Path, help="The pangenome.h5 file"
    )

    required.add_argument(
        "--fasta",
        required=False,
        type=Path,
        help="Specify a FASTA file containing the genomic sequences of the genome(s) you wish to annotate, "
        "or provide a tab-separated file listing genome names alongside their respective FASTA filepaths, with one line per genome.",
    )

    required.add_argument(
        "--anno",
        required=False,
        type=Path,
        help="Specify an annotation file in GFF/GBFF format for the genome you wish to annotate. "
        "Alternatively, you can provide a tab-separated file listing genome names alongside their respective annotation filepaths, "
        "with one line per genome. If both an annotation file and a FASTA file are provided, the annotation file will take precedence.",
    )

    required_single = parser.add_argument_group(
        title="Single Genome Arguments",
        description="Use these options when providing a single FASTA or annotation file:",
    )

    required_single.add_argument(
        "-n",
        "--genome_name",
        required=False,
        type=str,
        default="input_genome",
        help="Specify the name of the genome whose genome you want to annotate when providing a single FASTA or annotation file.",
    )

    required_single.add_argument(
        "--circular_contigs",
        nargs="+",
        required=False,
        type=tuple,
        help="Specify the contigs of the input genome that should be treated as circular when providing a single FASTA or annotation file.",
    )

    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument(
        "-o",
        "--output",
        required=False,
        type=Path,
        default="ppanggolin_projection"
        + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())
        + "_PID"
        + str(os.getpid()),
        help="Output directory",
    )

    optional.add_argument(
        "--no_defrag",
        required=False,
        action="store_true",
        help="DO NOT Realign gene families to link fragments with "
        "their non-fragmented gene family. (default: False)",
    )

    optional.add_argument(
        "--fast",
        required=False,
        action="store_true",
        help="Use representative sequences of gene families for input gene alignment. "
        "This option is faster but may be less sensitive. By default, all pangenome genes are used.",
    )

    optional.add_argument(
        "--identity",
        required=False,
        type=restricted_float,
        default=0.8,
        help="min identity percentage threshold",
    )

    optional.add_argument(
        "--coverage",
        required=False,
        type=restricted_float,
        default=0.8,
        help="min coverage percentage threshold",
    )

    optional.add_argument(
        "--use_pseudo",
        required=False,
        action="store_true",
        help="In the context of provided annotation, use this option to read pseudogenes. "
        "(Default behavior is to ignore them)",
    )

    optional.add_argument(
        "--dup_margin",
        required=False,
        type=restricted_float,
        default=0.05,
        help="minimum ratio of genomes in which the family must have multiple genes "
        "for it to be considered 'duplicated'. "
        "This metric is used to compute completeness and duplication of the input genomes",
    )

    optional.add_argument(
        "--soft_core",
        required=False,
        type=restricted_float,
        default=0.95,
        help="Soft core threshold used when generating general statistics on the projected genome. "
        "This threshold does not influence PPanGGOLiN's partitioning. "
        "The value determines the minimum fraction of genomes that must possess a gene family "
        "for it to be considered part of the soft core.",
    )

    optional.add_argument(
        "--spot_graph",
        required=False,
        action="store_true",
        help="Write the spot graph to a file, with pairs of blocks of single copy markers flanking RGPs "
        "as nodes. This graph can be used to visualize nodes that have RGPs from the input genome.",
    )

    optional.add_argument(
        "--graph_formats",
        required=False,
        type=str,
        choices=["gexf", "graphml"],
        nargs="+",
        default=["gexf"],
        help="Format of the output graph.",
    )

    optional.add_argument(
        "--gff",
        required=False,
        action="store_true",
        help="Generate GFF files with projected pangenome annotations for each input genome.",
    )

    optional.add_argument(
        "--proksee",
        required=False,
        action="store_true",
        help="Generate JSON map files for PROKSEE with projected pangenome annotations for each input genome.",
    )

    optional.add_argument(
        "--table",
        required=False,
        action="store_true",
        help="Generate a tsv file for each input genome with pangenome annotations.",
    )

    optional.add_argument(
        "--compress",
        required=False,
        action="store_true",
        help="Compress the files in .gz",
    )

    optional.add_argument(
        "--add_sequences",
        required=False,
        action="store_true",
        help="Include input genome DNA sequences in GFF and Proksee output.",
    )

    optional.add_argument(
        "-c",
        "--cpu",
        required=False,
        default=1,
        type=int,
        help="Number of available cpus",
    )

    optional.add_argument(
        "--tmpdir",
        required=False,
        type=Path,
        default=Path(tempfile.gettempdir()),
        help="directory for storing temporary files",
    )

    optional.add_argument(
        "--keep_tmp",
        required=False,
        default=False,
        action="store_true",
        help="Keeping temporary files (useful for debugging).",
    )

    optional.add_argument(
        "--add_metadata",
        required=False,
        action="store_true",
        help="Include metadata information in the output files "
        "if any have been added to pangenome elements (see ppanggolin metadata command).",
    )

    optional.add_argument(
        "--metadata_sources",
        default=None,
        nargs="+",
        help="Which source of metadata should be written. "
        "By default all metadata sources are included.",
    )

    optional.add_argument(
        "--metadata_sep",
        required=False,
        default="|",
        help="The separator used to join multiple metadata values for elements with multiple metadata"
        " values from the same source. This character should not appear in metadata values.",
    )
