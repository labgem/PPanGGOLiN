#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
from multiprocessing import get_context
import logging
import os
import time
from pathlib import Path
import tempfile
from typing import Tuple, Set, Dict, Iterator, Optional, List, Iterable, Any
from collections import defaultdict
import csv
from itertools import chain


# installed libraries
from tqdm import tqdm
import networkx as nx
import yaml
import pandas as pd


# # local libraries
from ppanggolin.annotate.synta import annotate_organism, read_fasta, get_dna_sequence
from ppanggolin.annotate.annotate import read_anno_file, launch_read_anno, launch_annotate_organism, local_identifiers_are_unique
from ppanggolin.annotate import subparser as annotate_subparser
from ppanggolin.pangenome import Pangenome
# from ppanggolin.genome import input_organism, Gene, RNA, Contig
from ppanggolin.utils import detect_filetype, create_tmpdir, read_compressed_or_not, write_compressed_or_not, restricted_float, mk_outdir, get_config_args, parse_config_file, get_default_args, check_input_files
from ppanggolin.align.alignOnPang import get_input_seq_to_family_with_rep,get_input_seq_to_family_with_all, project_and_write_partition
from ppanggolin.formats.writeSequences import write_gene_sequences_from_annotations
from ppanggolin.formats.readBinaries import check_pangenome_info
# from ppanggolin.formats import write_pangenome
from ppanggolin.RGP.genomicIsland import naming_scheme, compute_org_rgp
from ppanggolin.RGP.spot import make_spot_graph, check_sim, add_new_node_in_spot_graph, write_spot_graph
from ppanggolin.genome import Organism, Gene, RNA, Contig
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.region import Region, Spot, Module
from ppanggolin.formats.writeFlat import summarize_spots


class NewSpot(Spot):
    """
    This class represent a hotspot specifically 
    created for the projected genome.
    """

    def __str__(self):
        return f'new_spot_{str(self.ID)}'

def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """

    output_dir = Path(args.output)
    mk_outdir(output_dir, args.force)

    # For the moment these elements of the pangenome are predicted by default
    project_modules = True
    predict_rgp = True
    project_spots = True

    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)

    if pangenome.status["partitioned"] not in ["Computed", "Loaded", "inFile"]:
        raise NameError(f"The provided pangenome has not been partitioned. "
                        "Annotation of an external genome is therefore not possible. "
                        "See the 'partition' subcommands.")

    if pangenome.status["predictedRGP"] not in ["Computed", "Loaded", "inFile"]:
        logging.getLogger('PPanGGOLiN').info("RGPs have not been predicted in the provided pangenome. "
                                 "Projection of RGPs and spots into the provided genome will not be performed.")
        predict_rgp = False
        project_spots = False

    elif pangenome.status["spots"] not in ["Computed", "Loaded", "inFile"]:
        logging.getLogger('PPanGGOLiN').info("Spots have not been predicted in the provided pangenome. "
                                 "Projection of spots into the provided genome will not be performed.")
        project_spots = False

    if pangenome.status["modules"] not in ["Computed", "Loaded", "inFile"]:
        logging.getLogger('PPanGGOLiN').info("Modules have not been predicted in the provided pangenome. "
                                 "Projection of modules into the provided genome will not be performed.")

        project_modules = False
    
    if pangenome.status["geneSequences"] not in ["Loaded", "Computed", "inFile"] and not args.fast:
        raise Exception("The provided pangenome has no gene sequences. "
                        "Projection is still possible with the --fast option to use representative "
                        "sequences rather than all genes to annotate input genes.")

    if pangenome.status["geneFamilySequences"] not in ["Loaded", "Computed", "inFile"]:
        raise Exception("The provided pangenome has no gene families sequences. "
                        "This is not possible to annotate an input organism to this pangenome.")
    

    check_pangenome_info(pangenome, need_annotations=True, need_families=True, disable_bar=args.disable_prog_bar,
                         need_rgp=predict_rgp, need_modules=project_modules, need_gene_sequences=False,
                         need_spots=project_spots)


    logging.getLogger('PPanGGOLiN').info('Retrieving parameters from the provided pangenome file.')
    pangenome_params = argparse.Namespace(
        **{step: argparse.Namespace(**k_v) for step, k_v in pangenome.parameters.items()})
    

    genome_name_to_fasta_path, genome_name_to_annot_path = None, None

    if args.input_mode == "multiple":
        if args.anno:
            genome_name_to_annot_path = parse_input_paths_file(args.anno)
        
        if args.fasta:
            genome_name_to_fasta_path = parse_input_paths_file(args.fasta)

    else: #  args.input_mode == "single:

        circular_contigs = args.circular_contigs if args.circular_contigs else []
        if args.anno:
            genome_name_to_annot_path = {args.organism_name: {"path": args.annot,
                                                            "circular_contigs": circular_contigs}}
        
        if args.fasta:
            genome_name_to_fasta_path = {args.organism_name: {"path": args.fasta,
                                                            "circular_contigs": circular_contigs}}
    
    if genome_name_to_annot_path:
        check_input_names(pangenome, genome_name_to_annot_path)

        organisms, org_2_has_fasta = read_annotation_files(genome_name_to_annot_path, cpu=args.cpu, pseudo=args.use_pseudo,
                     disable_bar=args.disable_prog_bar)
        
        if not all((has_fasta for has_fasta in org_2_has_fasta.values())):
            organisms_with_no_fasta = {org for org, has_fasta in org_2_has_fasta.items() if not has_fasta}
            if args.fasta:
                get_gene_sequences_from_fasta_files(organisms_with_no_fasta, genome_name_to_fasta_path)
            else:
                raise ValueError(f"You provided GFF files for {len(organisms_with_no_fasta)} (out of {len(organisms)}) "
                                 "organisms without associated sequence data, and you did not provide "
                                "FASTA sequences using the --fasta or --single_fasta_file options. Therefore, it is impossible to project the pangenome onto the input genomes. "
                                f"The following organisms have no associated sequence data: {', '.join(o.name for o in organisms_with_no_fasta)}")

    elif genome_name_to_fasta_path:
        annotate_param_names = ["norna", "kingdom",
                                "allow_overlap", "prodigal_procedure"]

        annotate_params = manage_annotate_param(annotate_param_names, pangenome_params.annotate, args.config)

        
        check_input_names(pangenome, genome_name_to_fasta_path)
        organisms = annotate_fasta_files(genome_name_to_fasta_path=genome_name_to_fasta_path,  tmpdir=args.tmpdir, cpu=args.cpu,
                             translation_table=int(pangenome_params.cluster.translation_table), norna=annotate_params.norna, kingdom=annotate_params.kingdom,
                             allow_overlap=annotate_params.allow_overlap, procedure=annotate_params.prodigal_procedure, disable_bar=args.disable_prog_bar )


    input_org_to_lonely_genes_count = annotate_input_genes_with_pangenome_families(pangenome, input_organisms=organisms,
                                                                                output=output_dir, cpu=args.cpu, use_representatives=args.fast,
                                                                                no_defrag=args.no_defrag, identity=args.identity,
                                                                                coverage=args.coverage, tmpdir=args.tmpdir,
                                                                                translation_table=int(pangenome_params.cluster.translation_table),
                                                                                keep_tmp=args.keep_tmp, 
                                                                                disable_bar=args.disable_prog_bar)
            

    input_org_2_rgps, input_org_to_spots, input_orgs_to_modules = {}, {}, {}

    if predict_rgp:
        logging.getLogger('PPanGGOLiN').info('Detecting RGPs in input genomes.')

        logging.getLogger('PPanGGOLiN').debug("Detecting multigenic families...")
        multigenics = pangenome.get_multigenics(pangenome_params.rgp.dup_margin)

        input_org_2_rgps = predict_RGP(pangenome, organisms,  persistent_penalty=pangenome_params.rgp.persistent_penalty, variable_gain=pangenome_params.rgp.variable_gain,
                                     min_length=pangenome_params.rgp.min_length, min_score=pangenome_params.rgp.min_score, multigenics=multigenics, output_dir=output_dir,
                                     disable_bar=args.disable_prog_bar)


        if project_spots:
            logging.getLogger('PPanGGOLiN').info('Predicting spot of insertion in input genomes.')
            input_org_to_spots = predict_spots_in_input_organisms(initial_spots=list(pangenome.spots),
                                                            initial_regions=pangenome.regions,
                                                            input_org_2_rgps=input_org_2_rgps,
                                                            multigenics=multigenics,
                                                            output=output_dir,
                                                            write_graph_flag=args.spot_graph,
                                                            graph_formats=args.graph_formats,
                                                            overlapping_match=pangenome_params.spot.overlapping_match,
                                                            set_size=pangenome_params.spot.set_size,
                                                            exact_match=pangenome_params.spot.exact_match_size)

    if project_modules:
        input_orgs_to_modules = project_and_write_modules(pangenome, organisms, output_dir)

    organism_2_summary = {}
    for organism in organisms:
        # summarize projection for all input organisms
        organism_2_summary[organism] = summarize_projection(organism, pangenome,
                             input_org_2_rgps.get(organism, None),
                             input_org_to_spots.get(organism, None),
                             input_orgs_to_modules.get(organism, None),
                             input_org_to_lonely_genes_count[organism], output_dir)
        
        write_summaries(organism_2_summary, output_dir)


def annotate_fasta_files(genome_name_to_fasta_path: Dict[str,dict], tmpdir: str, cpu: int = 1, translation_table: int = 11,
                       kingdom: str = "bacteria", norna: bool = False, allow_overlap: bool = False, procedure: str = None,
                       disable_bar: bool = False):
    """
    Main function to annotate a pangenome

    :param genome_name_to_annot_path: 
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

        arguments.append((org_name, org_info['path'], org_info['circular_contigs'], tmpdir, translation_table,
                          norna, kingdom, allow_overlap, procedure))

    logging.getLogger("PPanGGOLiN").info(f"Annotating {len(arguments)} genomes using {cpu} cpus...")
    with get_context('fork').Pool(processes=cpu) as p:
        for organism in tqdm(p.imap_unordered(launch_annotate_organism, arguments), unit="genome",
                             total=len(arguments), disable=disable_bar):
            
            organisms.append(organism)
        p.close()
        p.join()

    return organisms


    
def read_annotation_files(genome_name_to_annot_path: Dict[str,dict], cpu: int = 1, pseudo: bool = False,
                     disable_bar: bool = False) -> Tuple[List[Organism], Dict[Organism,bool]]:
    """
    Read the annotation from GBFF file

    :param pangenome: pangenome object
    :param organisms_file: List of GBFF files for each organism
    :param cpu: number of CPU cores to use
    :param pseudo: allow to read pseudogène
    :param disable_bar: Disable the progresse bar
    """

    args = []
    organisms = []

    # we assume there are gene sequences in the annotation files,
    # unless a gff file without fasta is met (which is the only case where sequences can be absent)
    org_to_has_fasta_flag = {}


    for org_name, org_info in genome_name_to_annot_path.items():
        
        args.append((org_name, org_info['path'], org_info['circular_contigs'], pseudo))

    with get_context('fork').Pool(cpu) as p:
        for org, has_fasta in tqdm(p.imap_unordered(launch_read_anno, args), unit="file", total=len(args),
                              disable=disable_bar):
            
            organisms.append(org)
            org_to_has_fasta_flag[org] = has_fasta

    genes = (gene for org in organisms for gene in org.genes)

    if local_identifiers_are_unique(genes):
        for gene in genes:
            gene.ID = gene.local_identifier  # Erase ppanggolin generated gene ids and replace with local identifiers
            gene.local_identifier = ""  # this is now useless, setting it to default value

        logging.getLogger("PPanGGOLiN").info("Gene identifiers used in the provided annotation files were unique, "
                                             "PPanGGOLiN will use them.")
    else:
        logging.getLogger("PPanGGOLiN").info("Gene identifiers used in the provided annotation files were not unique, "
                                             "PPanGGOLiN will use self-generated identifiers.")
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
        raise ValueError(f"You did not provided fasta for all the organisms found in annotation file. "
                        f"{missing} are missing (out of {len(organisms)}). Missing organisms: {','.join(missing)}")
        
    for org in organisms:
        
        org_fasta_file = genome_name_to_annot_path[org.name]['path']

        with read_compressed_or_not(org_fasta_file) as currFastaFile:
            org_contig_to_seq, _ = read_fasta(org, currFastaFile)

        for contig in org.contigs:
            try:
                contig_seq = org_contig_to_seq[contig.name]
            except KeyError:
                msg = f"Fasta file for organism {org.name} did not have the contig {contig.name} " \
                      f"that was read from the annotation file. "
                msg += f"The provided contigs in the fasta were : " \
                       f"{', '.join([contig for contig in org_contig_to_seq])}."
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
        raise NameError(f"{len(duplicated_names)} provided organism names already exist in the given pangenome: {' '.join(duplicated_names)}")


def parse_input_paths_file(path_list_file: Path) -> Dict[str, Dict[str, List[str]]]:
    """
    Parse an input paths file to extract genome information.

    This function reads an input paths file, which is in TSV format, and extracts genome information
    including file paths and putative circular contigs.

    :param path_list_file: The path to the input paths file.
    :return: A dictionary where keys are genome names and values are dictionaries containing path information and
             putative circular contigs.
    :raises FileNotFoundError: If a specified genome file path does not exist.
    :raises Exception: If there are no genomes in the provided file.
    """
    logging.getLogger("PPanGGOLiN").info(f"Reading {path_list_file} to process organism files")
    genome_name_to_genome_path = {}

    for line in read_compressed_or_not(path_list_file):
        elements = [el.strip() for el in line.split("\t")]
        genome_file_path = Path(elements[1])
        genome_name = elements[0]
        putative_circular_contigs = elements[2:]

        if not genome_file_path.exists():  
            # Check if the file path doesn't exist and try an alternative path.
            genome_file_path_alt = path_list_file.parent.joinpath(genome_file_path)

            if not genome_file_path_alt.exists():
                raise FileNotFoundError(f"The file path '{genome_file_path}' for genome '{genome_name}' specified in '{path_list_file}' does not exist.")
            else:
                genome_file_path = genome_file_path_alt

        genome_name_to_genome_path[genome_name] = {
            "path": genome_file_path,
            "circular_contigs": putative_circular_contigs
        }

    if len(genome_name_to_genome_path) == 0:
        raise Exception(f"There are no genomes in the provided file: {path_list_file} ")
    
    return genome_name_to_genome_path



def write_summaries(organism_2_summary: Dict[Organism, Dict[str, Any]], output_dir: Path):
    """
    Write summary information to YAML files and create a summary projection in TSV format.

    This function takes a dictionary where keys are input organisms and values are dictionaries containing summary
    information. It writes this information to YAML files for each organism and creates a summary projection in TSV format.

    :param organism_2_summary: A dictionary where keys are input organisms and values are dictionaries containing
                               summary information.
    :param output_dir: The directory where the summary files will be written.
    """
    flat_summaries = []

    for input_organism, summary_info in organism_2_summary.items():
        yaml_string = yaml.dump(summary_info, default_flow_style=False, sort_keys=False, indent=4)

        with open(output_dir / input_organism.name / "projection_summary.yaml", 'w') as flout:
            flout.write('Projection_summary:')
            flout.write(yaml_string)

        flat_summary = {}
        for key, val in summary_info.items():
            if type(val) == dict:
                for nest_k, nest_v in val.items():
                    flat_summary[f"{key} {nest_k}"] =  nest_v
            else:
                flat_summary[key] = val

        flat_summaries.append(flat_summary)

    df_summary = pd.DataFrame(flat_summaries)

    df_summary.to_csv(output_dir / "summary_projection.tsv", sep='\t', index=False)


def summarize_projection(input_organism:Organism, pangenome:Pangenome, input_org_rgps:Region,
                         input_org_spots:Spot, input_org_modules:Module, singleton_gene_count:int, output_dir:Path):
    """

    :param singleton_gene_count: Number of genes that do not cluster with any of the gene families of the pangenome.

    """

    partition_to_gene = defaultdict(set)
    contigs_count = 0
    for contig in input_organism.contigs:
        contigs_count += 1
        for gene in contig.genes:
            partition_to_gene[gene.family.named_partition].add(gene)

    persistent_gene_count = len(partition_to_gene['persistent'])
    shell_gene_count = len(partition_to_gene['shell'])
    cloud_gene_count = len(partition_to_gene['cloud'])

    gene_count = persistent_gene_count + shell_gene_count + cloud_gene_count

    persistent_family_count = len({g.family for g in partition_to_gene['persistent']})
    shell_family_count = len({g.family for g in partition_to_gene['shell']})
    cloud_family_count = len({g.family for g in partition_to_gene['cloud']})

    families_count = persistent_family_count + shell_family_count + cloud_family_count

    rgp_count = "Not computed" if input_org_rgps is None else len(input_org_rgps)
    spot_count = "Not computed" if input_org_spots is None else len(input_org_spots)
    new_spot_count = "Not computed" if input_org_spots is None else sum(1 for spot in input_org_spots if isinstance(spot, NewSpot))
    module_count = "Not computed" if input_org_modules is None else len(input_org_modules)

    summary_info = {
        "Organism name": input_organism.name,
        "Pangenome file": pangenome.file,
        "Contigs": contigs_count,
        "Genes": gene_count,
        "Families": families_count,
        "Persistent": {"genes":persistent_gene_count, "families":persistent_family_count},
        "Shell": {"genes":shell_gene_count, "families":shell_family_count},
        "Cloud": {"genes":cloud_gene_count, "families":cloud_family_count - singleton_gene_count, "specific families":singleton_gene_count},
        "RGPs": rgp_count,
        "Spots": spot_count,
        "New spots": new_spot_count,
        "Modules": module_count
    }
    return summary_info

        
def annotate_input_genes_with_pangenome_families(pangenome: Pangenome, input_organisms: Iterable[Organism], output: Path,
                                                 cpu: int,use_representatives:bool, no_defrag: bool, 
                                                 identity: float, coverage: float, tmpdir: Path,
                                                 translation_table: int, keep_tmp:bool = False, disable_bar: bool =False):
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
    seq_fasta_files = []
    
    logging.getLogger('PPanGGOLiN').info(f'Writting gene sequences of input genomes.')

    for input_organism in input_organisms:

        seq_outdir = output / input_organism.name
        mk_outdir(seq_outdir, force=True)

        seq_fasta_file = seq_outdir / f"cds_sequences.fasta"

        with open(seq_fasta_file, "w") as fh_out_faa:
            write_gene_sequences_from_annotations(input_organism.genes, fh_out_faa, disable_bar=True, add=f"ppanggolin_")

        seq_fasta_files.append(seq_fasta_file)

    with create_tmpdir(main_dir=tmpdir, basename="align_input_seq_tmpdir", keep_tmp=keep_tmp) as new_tmpdir:
            
        if use_representatives:
            _, seqid_to_gene_family = get_input_seq_to_family_with_rep(pangenome, seq_fasta_files, output=new_tmpdir, tmpdir=new_tmpdir, is_input_seq_nt=True,
                                                        cpu=cpu, no_defrag=no_defrag, identity=identity, coverage=coverage, translation_table=translation_table)
        else:
            _, seqid_to_gene_family = get_input_seq_to_family_with_all(pangenome=pangenome, sequence_files=seq_fasta_files, 
                                                                                output=new_tmpdir, tmpdir=new_tmpdir, is_input_seq_nt=True,
                                                                                cpu=cpu, no_defrag=no_defrag, identity=identity, coverage=coverage,
                                                                                translation_table=translation_table, disable_bar=disable_bar)
    input_org_to_lonely_genes_count = {}
    for input_organism in input_organisms:
        
        org_outdir = output / input_organism.name

        seq_set = {gene.ID if gene.local_identifier == "" else gene.local_identifier for gene in input_organism.genes}

        project_and_write_partition(seqid_to_gene_family, seq_set, org_outdir)

        lonely_genes = set()
        for gene in input_organism.genes:
            gene_id = gene.ID if gene.local_identifier == "" else gene.local_identifier

            try:
                gene_family = seqid_to_gene_family[gene_id]
                gene_family.add(gene)
            except KeyError:
                # the seqid is not in the dict so it does not align with any pangenome families
                # We consider it as cloud gene
                try:
                    # in some case a family exists already and has the same name of the gene id
                    # So gene id cannot be used 
                    _ = pangenome.get_gene_family(gene_id)
                except KeyError:
                    new_gene_family = GeneFamily(pangenome.max_fam_id, gene_id)

                else:
                    # gene id already exists.
                    new_name=f"{input_organism.name}_{gene_id}"
                    logging.getLogger('PPanGGOLiN').warning('The input organism as a specific gene that does not align to any '
                                                            f'pangenome families with the same id ({gene_id}) than an existing gene family in the pangenome. '
                                                            f'The organism name is added to the family name: {new_name}')
                    new_gene_family = GeneFamily(pangenome.max_fam_id, new_name)

                pangenome.add_gene_family(new_gene_family)
                new_gene_family.add(gene)
                new_gene_family.partition = "Cloud"
                lonely_genes.add(gene)

        logging.getLogger('PPanGGOLiN').info(f"{input_organism.name} has {len(lonely_genes)}/{input_organism.number_of_genes()} "
                                "specific genes that do not align to any gene of the pangenome.")
        # Write specific gene ids in a file
        with open(org_outdir / "specific_genes.tsv", "w") as fl:
            fl.write('\n'.join((gene.ID if gene.local_identifier == "" else gene.local_identifier for gene in lonely_genes)) + '\n')

        input_org_to_lonely_genes_count[input_organism] = len(lonely_genes)

    return input_org_to_lonely_genes_count


def predict_RGP(pangenome: Pangenome, input_organisms: Organism, persistent_penalty: int, variable_gain: int,
                min_length: int, min_score: int, multigenics: float,
                output_dir:Path, disable_bar: bool) -> Dict[Organism, Set[Region]]:
    """
    Compute Regions of Genomic Plasticity (RGP) for the given input organisms.

    :param pangenome: The pangenome object.
    :param input_organisms: The input organism for which to compute RGPs.
    :param persistent_penalty: Penalty score to apply to persistent genes.
    :param variable_gain: Gain score to apply to variable genes.
    :param min_length: Minimum length (bp) of a region to be considered as RGP.
    :param min_score: Minimal score required for considering a region as RGP.
    :param multigenics: multigenic families.
    :param output_dir: Output directory where predicted rgps are going to be written.
    :param disable_bar: Flag to disable the progress bar.

    :return: Dictionary mapping organism with the set of predicted regions
    """

    logging.getLogger('PPanGGOLiN').info("Computing Regions of Genomic Plasticity...")

    name_scheme = naming_scheme(chain(pangenome.organisms, input_organisms))
    organism_to_rgps = {}

    for input_organism in input_organisms:
        rgps = compute_org_rgp(input_organism, multigenics, persistent_penalty, variable_gain, min_length,
                           min_score, naming=name_scheme, disable_bar=disable_bar)

        logging.getLogger('PPanGGOLiN').info(f"{len(rgps)} RGPs have been predicted in the input genomes.")
        
        
        org_outdir = output_dir / input_organism.name 
        
        write_predicted_regions(rgps, output=org_outdir, compress=False)
        organism_to_rgps[input_organism] = rgps

    return organism_to_rgps


def write_predicted_regions(regions: Set[Region],
                            output: Path, compress: bool = False):
    """
    Write the file providing information about predicted regions.

    :param regions: Set of Region objects representing predicted regions.
    :param output: Path to the output directory.
    :param compress: Whether to compress the file in .gz format.
    """
    fname = output / "plastic_regions.tsv"
    with write_compressed_or_not(fname, compress) as tab:
        fieldnames = ["region", "organism", "contig", "start",
                      "stop", "genes", "contigBorder", "wholeContig"]

        writer = csv.DictWriter(tab, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        regions = sorted(regions, key=lambda x: (
            x.organism.name, x.contig.name, x.ID))
        for region in regions:
            row = {
                "region": region.name,
                "organism": region.organism,
                "contig": region.contig,
                "start": region.starter,
                "stop": region.stopper,
                "genes": len(region),
                "contigBorder": region.is_contig_border,
                "wholeContig": region.is_whole_contig
            }

            writer.writerow(row)


def write_rgp_to_spot_table(rgp_to_spots: Dict[Region, Set[str]], output: Path, filename: str, compress: bool = False):
    """
    Write a table mapping RGPs to corresponding spot IDs.

    :param rgp_to_spots: A dictionary mapping RGPs to spot IDs.
    :param output: Path to the output directory.
    :param filename: Name of the file to write.
    :param compress: Whether to compress the file.
    """
    fname = output / filename
    logging.getLogger('PPanGGOLiN').debug(
        f'Writing RGPs to spot table in {fname}')

    with write_compressed_or_not(fname, compress) as tab:
        fieldnames = ["region", "spot_id"]

        writer = csv.DictWriter(tab, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        regions = sorted(rgp_to_spots.keys(), key=lambda x: (
            x.organism.name, x.contig.name, x.ID))
        for region in regions:
            row = {
                "region": region.name,
                "spot_id": ';'.join(map(str, rgp_to_spots[region]))
            }

            writer.writerow(row)


def retrieve_gene_sequences_from_fasta_file(input_organism, fasta_file):
    """
    Get gene sequences from fastas

    :param pangenome: input pangenome
    :param fasta_file: list of fasta file
    """

    with read_compressed_or_not(fasta_file) as currFastaFile:
        contig_id2deq, _ = read_fasta(input_organism, currFastaFile)

    for contig in input_organism.contigs:
        try:
            for gene in contig.genes:
                gene.add_dna(get_dna_sequence(
                    contig_id2deq[contig.name], gene))

            for rna in contig.RNAs:
                rna.add_dna(get_dna_sequence(contig_id2deq[contig.name], rna))
        except KeyError:
            msg = f"Fasta file for input_organism {input_organism.name} did not have the contig {contig.name} " \
                f"that was read from the annotation file. "
            msg += f"The provided contigs in the fasta were : " \
                f"{', '.join([contig for contig in contig_id2deq.keys()])}."
            raise KeyError(msg)


def manage_annotate_param(annotate_param_names: List[str], pangenome_args: argparse.Namespace,
                          config_file: Optional[str]) -> argparse.Namespace:
    """
    Manage annotate parameters by collecting them from different sources and merging them.

    :param annotate_param_names: List of annotate parameter names to be managed.
    :param pangenome_args: Annotate arguments parsed from pangenomes parameters.
    :param config_file: Path to the config file, can be None if not provided.

    :return: An argparse.Namespace containing the merged annotate parameters with their values.
    """

    default_annotate_args = get_default_args('annotate', annotate_subparser)

    if config_file is None:
        config_annotate_args = argparse.Namespace()
    else:
        config = defaultdict(dict, parse_config_file(config_file))
        config_annotate_args = get_config_args(
            'annotate', annotate_subparser, config, "annotate", annotate_param_names, strict_config_check=False)

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
        param_val_string = ' '.join(
            [f'--{k} {v}' for k, v in annotate_param_from_pangenome.items()])
        logging.getLogger("PPanGGOLiN").debug(f"{len(annotate_param_from_pangenome)}/{len(annotate_param_names)} annotate parameters extracted from pangenome parameters "
                                              f"(the parameters used to build the input pangenome): {param_val_string}")

    if len(annotate_param_from_config) > 0:
        param_val_string = ';'.join(
            [f' {k} : {v}' for k, v in annotate_param_from_config.items()])
        logging.getLogger("PPanGGOLiN").debug(f"{len(annotate_param_from_config)}/{len(annotate_param_names)} annotate parameters were not found in pangenome internal parameters."
                                              f" They have been parsed from the annotate section in the config file: {param_val_string}")

    if len(annotate_param_from_default) > 0:
        param_val_string = ';'.join(
            [f' {k} : {v}' for k, v in annotate_param_from_default.items()])
        logging.getLogger("PPanGGOLiN").debug(f"{len(annotate_param_from_default)}/{len(annotate_param_names)} annotate parameters were not found in the pangenome parameters "
                                              f"nor in the config file. Default values have been used: {param_val_string}")

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
        assert len(
            spot_in_cc) == 1, "More than one spot in a connected_components. Something went wrong when recomputing spots."
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
    output: str,
    write_graph_flag: bool = False, 
    graph_formats: List[str] = ['gexf'],
    overlapping_match: int = 2, 
    set_size: int = 3, 
    exact_match: int = 1 ) -> Dict[Organism, Set[Spot]]:
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

    :return: A dictionary mapping input organism RGPs to their predicted spots.
    """

    logging.getLogger("PPanGGOLiN").debug(f"Rebuilding original spot graph.")
    graph_spot = make_spot_graph(rgps=initial_regions, multigenics=multigenics,
                                 overlapping_match=overlapping_match, set_size=set_size, exact_match=exact_match)

    original_nodes = set(graph_spot.nodes)

    # Check congruency with already computed spot and add spot id in node attributes
    check_spots_congruency(graph_spot, initial_spots)
    
    new_spot_id_counter = max((s.ID for s in initial_spots)) + 1

    input_org_to_spots = {}
    for input_organism, rgps in input_org_2_rgps.items():
        
        if len(rgps) == 0:
            logging.getLogger('PPanGGOLiN').debug(f"{input_organism.name}: No RGPs have been found. "
                                             "As a result, spot prediction and RGP output will be skipped.")
        
            input_org_to_spots[input_organism] = set()
            continue
    
        outdir_org = output / input_organism.name
        # Copy the graph spot, as each input organism are processed independently
        graph_spot_cp = graph_spot.copy()
        
        input_org_spots = predict_spot_in_one_organism(graph_spot_cp, input_org_rgps=rgps, original_nodes=original_nodes, 
                                          new_spot_id_counter=new_spot_id_counter, multigenics=multigenics, organism_name=input_organism.name,
                                          output=outdir_org, write_graph_flag=write_graph_flag, graph_formats=graph_formats,
                                        overlapping_match=overlapping_match, set_size=set_size, exact_match=exact_match)
        
        new_spot_id_counter = max((s.ID for s in input_org_spots)) + 1
        
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
    graph_formats: List[str] = ['gexf'],
    overlapping_match: int = 2, 
    set_size: int = 3, 
    exact_match: int = 1 ) -> Set[Spot]:
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
        logging.getLogger("PPanGGOLiN").debug(f"{organism_name}: no RGPs of the input organism will be associated with any spot of insertion "
                                             "as they are on a contig border (or have "
                                             f"less than {set_size} persistent gene families until the contig border). "
                                             "Projection of spots stops here")
        return {}

    # remove node that were already in the graph
    new_nodes = set(input_org_node_to_rgps) - original_nodes

    logging.getLogger("PPanGGOLiN").debug(f"{organism_name}: {lost} RGPs were not used as they are on a contig border (or have"
                                         f"less than {set_size} persistent gene families until the contig border)")

    logging.getLogger("PPanGGOLiN").debug(
        f"{organism_name}: {used} RGPs of the input organism will be associated to a spot of insertion")

    # add potential edges from new nodes to the rest of the nodes
    all_nodes = list(graph_spot.nodes)
    for nodei in new_nodes:
        for nodej in all_nodes:
            if nodei == nodej:
                continue
            node_obj_i = graph_spot.nodes[nodei]
            node_obj_j = graph_spot.nodes[nodej]
            if check_sim([node_obj_i["border0"], node_obj_i["border1"]],
                         [node_obj_j["border0"], node_obj_j["border1"]],
                         overlapping_match, set_size, exact_match):
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
                spots_of_the_cc |= {
                    spot for spot in graph_spot.nodes[node]["spots"]}

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
            logging.getLogger("PPanGGOLiN").debug(f'{organism_name}: Some RGPs of the input organism '
                                                 f"are connected to {len(spots_of_the_cc)} original spots of the pangenome.")

        input_rgps_of_the_cc = set()
        for node in comp:
            if node in input_org_node_to_rgps:
                input_rgps_of_the_cc |= input_org_node_to_rgps[node]

                if write_graph_flag:
                    graph_spot.nodes[node]["spots"] = spots_of_the_cc

                    graph_spot.nodes[node]["spot_id"] = ';'.join(
                        (str(spot) for spot in spots_of_the_cc))
                    graph_spot.nodes[node]["includes_RGPs_from_the_input_organism"] = True

        for spot in spots_of_the_cc:
            for region in input_rgps_of_the_cc:
                spot.add(region)

        input_rgp_to_spots.update(
            {rgp: spots_of_the_cc for rgp in input_rgps_of_the_cc})

    if write_graph_flag:
        # remove node that would not be writable in graph file
        for node in graph_spot.nodes:
            del graph_spot.nodes[node]["spots"]

        write_spot_graph(graph_spot, output, graph_formats,
                         file_basename='projected_spotGraph')

    write_rgp_to_spot_table(input_rgp_to_spots, output=output,
                            filename='input_organism_rgp_to_spot.tsv')

    input_org_spots = {spot for spots in input_rgp_to_spots.values()
                 for spot in spots }
    new_spots = {spot for spot in input_org_spots if type(spot) == NewSpot}


    logging.getLogger('PPanGGOLiN').debug(
        f'{organism_name}: {len(new_spots)} new spots have been created for the input genome.')

    if new_spots:
        summarize_spots(new_spots, output, compress=False,
                        file_name="new_spots_summary.tsv")
        
    return input_org_spots

def project_and_write_modules(pangenome: Pangenome, input_organisms: Iterable[Organism],
                              output: Path, compress: bool = False):
    """
    Write a tsv file providing association between modules and the input organism

    :param pangenome: Pangenome object
    :param input_organisms: iterable of the organisms that is being annotated
    :param output: Path to output directory
    :param compress: Compress the file in .gz
    """
    input_orgs_to_modules = {}
    for input_organism in input_organisms:
        output_file = output / input_organism.name / "modules_in_input_organism.tsv"

        input_organism_families = list(input_organism.families)
        counter = 0
        modules_in_input_org = []
        with write_compressed_or_not(output_file, compress) as fout:
            fout.write("module_id\torganism\tcompletion\n")

            for mod in pangenome.modules:
                module_in_input_organism = any(
                    (fam in input_organism_families for fam in mod.families))

                if module_in_input_organism:
                    counter += 1
                    modules_in_input_org.append(mod)

                    completion = round(
                        len(set(input_organism.families) & set(mod.families)) / len(set(mod.families)), 2)
                    fout.write(
                        f"module_{mod.ID}\t{input_organism.name}\t{completion}\n")

        logging.getLogger('PPanGGOLiN').debug(
            f"{input_organism.name}: {counter} modules have been projected to the input genomes.")

        logging.getLogger('PPanGGOLiN').debug(
            f"{input_organism.name}: Projected modules have been written in: '{output_file}'")
        
        input_orgs_to_modules[input_organism] = modules_in_input_org
    
    return input_orgs_to_modules


def determine_input_mode(input_file: Path, expected_types: List[str], parser: argparse.ArgumentParser) -> str:
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
        parser.error("Based on its content, the provided file is not recognized as a valid input file. Please ensure it is in one of the supported formats (FASTA, GFF/GBFF, or TSV).")

    if filetype == "tsv":
        logging.getLogger('PPanGGOLiN').debug(f"The provided file ({input_file}) is detected as a TSV file.")
        mode = "multiple"
    elif filetype in expected_types:
        logging.getLogger('PPanGGOLiN').debug(f"The provided file ({input_file}) is detected as a single {'/'.join(expected_types)} file.")
        mode = "single"
    else:
        logging.getLogger('PPanGGOLiN').error(f"The provided file {input_file} is not recognized as a valid {'/'.join(expected_types)} file or a TSV file listing names and {'/'.join(expected_types)} files of genomes to annotate.")
        parser.error(f"The provided file {input_file} is not recognized as a valid {'/'.join(expected_types)} file or a TSV file listing names and files of genomes to annotate.")

    return mode


def check_projection_arguments(args: argparse.Namespace, parser: argparse.ArgumentParser ) -> str:
    """
    Check the arguments provided for genome projection and raise errors if they are incompatible or missing.

    :param args: An argparse.Namespace object containing parsed command-line arguments.
    :param parser : parser of the command
    :return: A string indicating the input mode ('single' or 'multiple').
    """
    
    # Check if we annotate genomes from path files or only a single genome...  
    if not args.anno and not args.fasta:
        parser.error("Please provide either a FASTA file or a tab-separated file listing sequence files using the '--fasta' option, "
                    "or an annotation file or a tab-separated file listing annotation files using the '--anno' option. "
                    "You can specify these either through the command line or the configuration file.")

    mode_from_fasta, mode_from_anno = None, None
    if args.fasta:
        mode_from_fasta = determine_input_mode(args.fasta, ['fasta'], parser)
        input_mode = mode_from_fasta

    if args.anno:
        mode_from_anno = determine_input_mode(args.anno, ['gff', "gbff"], parser)
        input_mode = mode_from_anno

        logging.getLogger('PPanGGOLiN').debug("")

    if mode_from_fasta and mode_from_anno and mode_from_fasta != mode_from_anno:
        single_input, multiple_input = ("fasta", "anno") if mode_from_fasta == "single" else ("anno", "fasta")

        parser.error(f"You've provided both a single annotation/fasta file using the '--{single_input}' option and a list of files using "
                    f"the '--{multiple_input}' option. Please choose either a single file or a tab-separated file listing genome files, but not both.")


    if input_mode == "multiple":
        # We are in paths file mode
        
        incompatible_args = ["organism_name", "circular_contigs"]
        for single_arg in incompatible_args:
            if getattr(args, single_arg) is not None:
                parser.error("You provided a TSV file listing the files of genomes you wish to annotate. "
                             f"Therefore, the single genome argument '--{single_arg}' is incompatible with this multiple genomes file.")

        if args.fasta:
            check_input_files(args.fasta, True)

        if args.anno:
            check_input_files(args.anno, True)

    elif input_mode == "single":
        # We are in single file mode
            
        if args.organism_name is None:
            parser.error("You directly provided a single FASTA/GBFF/GFF file. Please specify the name of the input organism you want to annotate. "
                        "You can use the --organism_name argument either through the command line or the config file.")
    
    return input_mode


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for projection command

    :return : parser arguments for projection command
    """
    parser = sub_parser.add_parser(
        "projection", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_projection(parser)
    return parser


def parser_projection(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of projection command

    :param parser: parser for projection argument
    """
    required = parser.add_argument_group(title="Required arguments")

    required.add_argument('-p', '--pangenome', required=False,
                          type=Path, help="The pangenome.h5 file")

    required.add_argument('--fasta', required=False, type=Path,
                                help="Specify a FASTA file containing the genomic sequences of the organism(s) you wish to annotate, "
                                "or provide a tab-separated file listing organism names alongside their respective FASTA filepaths, with one line per organism.")
    
    required.add_argument('--anno', required=False, type=Path,
                                    help="Specify an annotation file in GFF/GBFF format for the genome you wish to annotate. "
                                    "Alternatively, you can provide a tab-separated file listing organism names alongside their respective annotation filepaths, "
                                    "with one line per organism. If both an annotation file and a FASTA file are provided, the annotation file will take precedence.")

    required_single = parser.add_argument_group(title="Single Genome Arguments",
                                                description="Use these options when providing a single FASTA or annotation file:")

    required_single.add_argument("-n", '--organism_name', required=False, type=str,
                        help="Specify the name of the organism whose genome you want to annotate when providing a single FASTA or annotation file.")

    required_single.add_argument('--circular_contigs', nargs="+", required=False, type=tuple,
                                help="Specify the contigs of the input genome that should be treated as circular when providing a single FASTA or annotation file.")


    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument('-o', '--output', required=False, type=Path,
                          default="ppanggolin_projection" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                      time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory")

    optional.add_argument('--no_defrag', required=False, action="store_true",
                          help="DO NOT Realign gene families to link fragments with "
                               "their non-fragmented gene family. (default: False)")
    
    optional.add_argument("--fast", required=False, action="store_true",
                        help="Use representative sequences of gene families for input gene alignment. "
                            "This option is faster but may be less sensitive. By default, all pangenome genes are used.")

    optional.add_argument('--identity', required=False, type=restricted_float, default=0.8,
                          help="min identity percentage threshold")

    optional.add_argument('--coverage', required=False, type=restricted_float, default=0.8,
                          help="min coverage percentage threshold")

    optional.add_argument("--use_pseudo", required=False, action="store_true",
                          help="In the context of provided annotation, use this option to read pseudogenes. "
                               "(Default behavior is to ignore them)")

    optional.add_argument("--spot_graph", required=False, action="store_true",
                          help="Write the spot graph to a file, with pairs of blocks of single copy markers flanking RGPs "
                          "as nodes. This graph can be used to visualize nodes that have RGPs from the input organism.")

    optional.add_argument('--graph_formats', required=False, type=str, choices=['gexf', "graphml"], nargs="+",
                          default=['gexf'], help="Format of the output graph.")

    optional.add_argument("-c", "--cpu", required=False,
                          default=1, type=int, help="Number of available cpus")

    optional.add_argument("--tmpdir", required=False, type=Path, default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
        
    optional.add_argument("--keep_tmp", required=False, default=False, action="store_true",
                        help="Keeping temporary files (useful for debugging).")
