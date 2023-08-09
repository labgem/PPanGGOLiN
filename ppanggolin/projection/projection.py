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
from typing import Tuple, Set, Dict, Iterator, Optional, List
from collections import defaultdict
import csv

# installed libraries
from tqdm import tqdm
import networkx as nx

# # local libraries
from ppanggolin.annotate.synta import annotate_organism, read_fasta, get_dna_sequence
from ppanggolin.annotate.annotate import read_anno_file
from ppanggolin.annotate import subparser as annotate_subparser
from ppanggolin.pangenome import Pangenome
# from ppanggolin.genome import input_organism, Gene, RNA, Contig
from ppanggolin.utils import read_compressed_or_not, write_compressed_or_not, restricted_float, mk_outdir, get_config_args, parse_config_file, get_default_args
from ppanggolin.align.alignOnPang import get_seq2pang, project_and_write_partition
from ppanggolin.formats.writeSequences import write_gene_sequences_from_annotations
from ppanggolin.formats.readBinaries import check_pangenome_info
# from ppanggolin.formats import write_pangenome
from ppanggolin.RGP.genomicIsland import naming_scheme, compute_org_rgp
from ppanggolin.RGP.spot import make_spot_graph, check_sim, add_new_node_in_spot_graph, write_spot_graph
from ppanggolin.genome import Organism, Gene, RNA, Contig
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.region import Region, Spot
from ppanggolin.formats.writeFlat import summarize_spots



class NewSpot(Spot):
    """
    This class represent a hotspot specifically 
    created for the projected genome.
    """

    def __str__(self):
        return f'new_spot_{str(self.ID)}'


def annotate_input_genes_with_pangenome_families(pangenome: Pangenome, input_organism: Organism, output: Path, cpu: int, no_defrag: bool,
                                                  identity: float, coverage: float, tmpdir: Path,
                                                  translation_table: int):
    """
    Annotate input genes with pangenome gene families and perform clustering.

    :param pangenome: Pangenome object.
    :param input_organism: Input organism object.
    :param output: Output directory for generated files.
    :param cpu: Number of CPU cores to use.
    :param no_defrag: Whether to use defragmentation.
    :param identity: Minimum identity threshold for gene clustering.
    :param coverage: Minimum coverage threshold for gene clustering.
    :param tmpdir: Temporary directory for intermediate files.
    :param disable_bar: Whether to disable progress bar.
    :param translation_table: Translation table ID for nucleotide sequences.

    :return: None
    """

    seq_fasta_file = output / f"{input_organism.name}.fasta"

    with open(seq_fasta_file, "w") as fh_out_faa:
        write_gene_sequences_from_annotations(input_organism.genes, fh_out_faa, disable_bar=True)
    
    with tempfile.TemporaryDirectory(dir=tmpdir, prefix="seq_to_pang_tmpdir_") as new_tmpdir:
        seq_set, _, seqid_to_gene_family = get_seq2pang(pangenome, seq_fasta_file, output, Path(new_tmpdir),
                                                        cpu, no_defrag, identity=identity, coverage=coverage,
                                                        is_nucleotide=True, translation_table=translation_table)
    
    project_and_write_partition(seqid_to_gene_family, seq_set, output)

    lonely_gene = 0
    for gene in input_organism.genes:
        try:
            gene_family = seqid_to_gene_family[gene.ID]
            gene_family.add_gene(gene)
        except KeyError:
            new_gene_family = pangenome.add_gene_family(gene.ID)
            new_gene_family.add_gene(gene)
            new_gene_family.add_partition("Cloud")
            lonely_gene += 1

    logging.getLogger().info(f"The input organism has {lonely_gene}/{input_organism.number_of_genes()} " 
                             "genes that do not cluster with any of the gene families in the pangenome.")

def predict_RGP(pangenome: Pangenome, input_organism: Organism, persistent_penalty: int, variable_gain: int,
                min_length: int, min_score: int, multigenics: float,
                disable_bar: bool) -> None:
    """
    Compute Regions of Genomic Plasticity (RGP) for the given pangenome and input organism.

    :param pangenome: The pangenome object.
    :param input_organism: The input organism for which to compute RGPs.
    :param persistent_penalty: Penalty score to apply to persistent genes.
    :param variable_gain: Gain score to apply to variable genes.
    :param min_length: Minimum length (bp) of a region to be considered as RGP.
    :param min_score: Minimal score required for considering a region as RGP.
    :param multigenics: multigenic families.
    :param disable_bar: Flag to disable the progress bar.

    :return: None
    """

    logging.getLogger().info("Computing Regions of Genomic Plasticity...")
    name_scheme = naming_scheme(pangenome)

    rgps = compute_org_rgp(input_organism, multigenics, persistent_penalty, variable_gain, min_length,
                    min_score, naming=name_scheme, disable_bar=disable_bar)

    logging.getLogger().info(f"{len(rgps)} RGPs have been predicted in the input genomes.")
    return rgps


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
        fieldnames = ["region", "organism", "contig", "start", "stop", "genes", "contigBorder", "wholeContig"]

        writer = csv.DictWriter(tab, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        regions = sorted(regions, key=lambda x: (x.organism.name, x.contig.name, x.start))
        for region in regions:
            row = {
                "region": region.name,
                "organism": region.organism,
                "contig": region.contig,
                "start": region.start,
                "stop": region.stop,
                "genes": len(region.genes),
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
    logging.getLogger('PPanGGOLiN').info(f'Writing RGPs to spot table in {fname}')
    
    with write_compressed_or_not(fname, compress) as tab:
        fieldnames = ["region", "spot_id"]

        writer = csv.DictWriter(tab, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        regions = sorted(rgp_to_spots.keys(), key=lambda x: (x.organism.name, x.contig.name, x.start))
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
                gene.add_dna(get_dna_sequence(contig_id2deq[contig.name], gene))

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
        config_annotate_args = get_config_args('annotate', annotate_subparser, config, "annotate", annotate_param_names, strict_config_check=False)

    annotate_param_from_pangenome = {}
    annotate_param_from_config = {}
    annotate_param_from_default = {}

    annotate_params = argparse.Namespace()

    # Collecting annotate parameters from different sources
    for annotate_arg in annotate_param_names:
        if hasattr(pangenome_args, annotate_arg):
            param_val = getattr(pangenome_args, annotate_arg) 
            annotate_param_from_pangenome[annotate_arg] = param_val
            setattr(annotate_params, annotate_arg, param_val)

        elif hasattr(config_annotate_args, annotate_arg):
            param_val =  getattr(config_annotate_args, annotate_arg) 
            annotate_param_from_config[annotate_arg] = param_val
            setattr(annotate_params, annotate_arg, param_val)

        else:
            param_val =  getattr(default_annotate_args, annotate_arg) 
            annotate_param_from_default[annotate_arg] = param_val
            setattr(annotate_params, annotate_arg, param_val)

    # Log the sources of the annotate parameters
    if len(annotate_param_from_pangenome) > 0:
        param_val_string = ' '.join([f'--{k} {v}' for k, v in annotate_param_from_pangenome.items()])
        logging.getLogger("PPanGGOLiN").debug(f"{len(annotate_param_from_pangenome)}/{len(annotate_param_names)} annotate parameters extracted from pangenome parameters "
                                              f"(the parameters used to build the input pangenome): {param_val_string}")

    if len(annotate_param_from_config) > 0:
        param_val_string = ';'.join([f' {k} : {v}' for k, v in annotate_param_from_config.items()])
        logging.getLogger("PPanGGOLiN").debug(f"{len(annotate_param_from_config)}/{len(annotate_param_names)} annotate parameters were not found in pangenome internal parameters."
                                            f" They have been parsed from the annotate section in the config file: {param_val_string}")

    if len(annotate_param_from_default) > 0:
        param_val_string = ';'.join([f' {k} : {v}' for k, v in annotate_param_from_default.items()])
        logging.getLogger("PPanGGOLiN").debug(f"{len(annotate_param_from_default)}/{len(annotate_param_names)} annotate parameters were not found in the pangenome parameters "
                                            f"nor in the config file. Default values have been used: {param_val_string}")

    return annotate_params




def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """


    output_dir = Path(args.output)
    mk_outdir(output_dir, args.force)


    # For the moment this element of the pangenome are predicted by default
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
        logging.getLogger().info("RGPs have not been predicted in the provided pangenome. "
                                    "Projection of RGPs into the provided genome will not be performed.")
        predict_rgp = False

    if pangenome.status["spots"] not in ["Computed", "Loaded", "inFile"]:
        logging.getLogger().info("Spots have not been predicted in the provided pangenome. "
                                    "Projection of spots into the provided genome will not be performed.")
        project_spots = False


    if pangenome.status["modules"] not in ["Computed", "Loaded", "inFile"]:
        logging.getLogger().info("Modules have not been predicted in the provided pangenome. "
                                    "Projection of modules into the provided genome will not be performed.")
        
        project_modules = False


    check_pangenome_info(pangenome, need_annotations=True, need_families=True, disable_bar=args.disable_prog_bar, 
                         need_rgp=predict_rgp, need_modules=project_modules,
                         need_spots=project_spots)

    logging.getLogger().info('Retrieving parameters from the provided pangenome file.')
    pangenome_params = argparse.Namespace(**{step:argparse.Namespace(**k_v)  for step, k_v in pangenome.parameters.items()})

    
    if args.organism_name in [org.name for org in pangenome.organisms]:
        raise NameError(f"The provided organism name '{args.organism_name}' already exists in the given pangenome.")
    
    if args.annot_file is not None:
        # read_annotations(pangenome, args.anno, cpu=args.cpu, pseudo=args.use_pseudo, disable_bar=args.disable_prog_bar)
        input_organism, has_sequence = read_anno_file(organism_name = args.organism_name, 
                       filename=args.annot_file,
                        circular_contigs=[],
                        pseudo=args.use_pseudo)
        
        if not has_sequence:
            if args.fasta_file:
                retrieve_gene_sequences_from_fasta_file(input_organism, args.fasta_file)
            else:
                raise Exception("The gff/gbff provided did not have any sequence information, "
                            "Thus, we do not have the information we need to continue the projection.")

    elif args.fasta_file is not None:
        annotate_param_names = ["norna", "kingdom", "allow_overlap", "prodigal_procedure"]
                    
        annotate_params =  manage_annotate_param(annotate_param_names, pangenome_params.annotate, args.config)

        input_organism = annotate_organism(org_name=args.organism_name, file_name = args.fasta_file, circular_contigs=[], tmpdir=args.tmpdir,
                      code = args.translation_table, norna=annotate_params.norna, kingdom = annotate_params.kingdom,
                      overlap=annotate_params.allow_overlap, procedure=annotate_params.prodigal_procedure)

    else:
        raise Exception("At least one of --fasta_file or --anno_file must be given")


    # Add input organism in pangenome. This is temporary as the pangenome object is not going to be written.
    pangenome.add_organism(input_organism)

    annotate_input_genes_with_pangenome_families(pangenome, input_organism=input_organism, output=output_dir, cpu=args.cpu, 
                                                 no_defrag = args.no_defrag, identity = args.identity, coverage = args.coverage, tmpdir=args.tmpdir,
                                                 translation_table = args.translation_table)
    
    if predict_rgp:
        logging.getLogger().info('Detecting rgp in input genome.')


        logging.getLogger().info("Detecting multigenic families...")
        multigenics = pangenome.get_multigenics(pangenome_params.rgp.dup_margin)
        
        input_org_rgps = predict_RGP(pangenome, input_organism,  persistent_penalty=pangenome_params.rgp.persistent_penalty, variable_gain=pangenome_params.rgp.variable_gain,
                                    min_length=pangenome_params.rgp.min_length, min_score=pangenome_params.rgp.min_score, multigenics=multigenics, 
                                    disable_bar=args.disable_prog_bar)
        if len(input_org_rgps) == 0:

            logging.getLogger('PPanGGOLiN').info("No RGPs have been found in the input organisms. "
                                    "As a result, spot prediction and RGP output will be skipped.")
            
            
        else:

            write_predicted_regions(input_org_rgps, output=output_dir, compress=False)

            if project_spots and len(input_org_rgps) > 0:
                predict_spots_in_input_organism(initial_spots=pangenome.spots,
                                                                        initial_regions=pangenome.regions, 
                                                                        input_org_rgps=input_org_rgps,
                                                                        multigenics=multigenics, output=output_dir,
                                                                        write_graph_flag=args.spot_graph, graph_formats=args.graph_formats,
                                                                        overlapping_match=pangenome_params.spot.overlapping_match,
                                                                        set_size=pangenome_params.spot.set_size,
                                                                        exact_match=pangenome_params.spot.exact_match_size)

                    

            
        project_and_write_modules(pangenome, input_organism, output_dir)


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
        assert len(spot_in_cc) == 1, "More than one spot in a connected_components. Something went wrong when recomputing spots."
        current_spot = spot_in_cc.pop()
        # Add spot id to the graph
        for node in cc:
            graph_spot.nodes[node]["spot_id"] = str(current_spot)
            graph_spot.nodes[node]["spots"] = {current_spot}



def predict_spots_in_input_organism(initial_spots: List[Spot], initial_regions: List[Region], 
                                    input_org_rgps: List[Region], multigenics: Set[GeneFamily], output: str, 
                                    write_graph_flag: bool = False, graph_formats: List[str] = ['gexf'], 
                                    overlapping_match: int = 2, set_size: int = 3, exact_match: int = 1) -> Dict:
    """
    Create a spot graph from pangenome RGP and predict spots for input organism RGPs.

    :param initial_spots: List of original spots in the pangenome.
    :param initial_regions: List of original regions in the pangenome.
    :param input_org_rgps: List of RGPs from the input organism to be associated with spots.
    :param multigenics: Set of pangenome graph multigenic persistent families.
    :param output: Output directory to save the spot graph.
    :param write_graph_flag: If True, writes the spot graph in the specified formats.
    :param graph_formats: List of graph formats to write (default is ['gexf']).
    :param overlapping_match: Number of missing persistent genes allowed when comparing flanking genes.
    :param set_size: Number of single copy markers to use as flanking genes for RGP during hotspot computation.
    :param exact_match: Number of perfectly matching flanking single copy markers required to associate RGPs.

    :return: A dictionary mapping input organism RGPs to their predicted spots.
    """

    logging.getLogger("PPanGGOLiN").info(f"Rebuilding spot graph.")
    graph_spot = make_spot_graph(rgps=initial_regions, multigenics=multigenics,
                    overlapping_match=overlapping_match, set_size=set_size, exact_match=exact_match)
    
    original_nodes = set(graph_spot.nodes)

    # Check congruency with already computed spot and add spot id in node attributes
    check_spots_congruency(graph_spot, initial_spots)

    
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
        logging.getLogger("PPanGGOLiN").info(f"No RGPs of the input organism will be associated with any spot of insertion "
                                            "as they are on a contig border (or have "
                                            f"less than {set_size} persistent gene families until the contig border). "
                                            "Projection of spots stops here")
        return {}

    # remove node that were already in the graph 
    new_nodes = set(input_org_node_to_rgps) - original_nodes 
    
    logging.getLogger("PPanGGOLiN").info(f"{lost} RGPs of the input organism won't be associated with any spot of insertion "
                                         "as they are on a contig border (or have "
                                         f"less than {set_size} persistent gene families until the contig border)")
    
    logging.getLogger("PPanGGOLiN").info(f"{used} RGPs of the input organism will be associated to a spot of insertion")

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
    new_spot_id_counter = max((s.ID for s in initial_spots)) + 1 
    # determine spot ids of the new nodes and by extension to their rgps
    for comp in nx.algorithms.components.connected_components(graph_spot):
        # in very rare case one cc can have several original spots
        # that would mean a new nodes from the input organism have connected two old cc
        # in this case we report the two spots in the output
        spots_of_the_cc = set()
        for node in comp:
            if "spots" in graph_spot.nodes[node]:
                spots_of_the_cc |= {spot for spot in graph_spot.nodes[node]["spots"]}

        if len(spots_of_the_cc) == 0:
            # no spot associated with any node of the cc
            # that means this cc is only composed of new nodes
            # let's add a new spot id
            new_spot = NewSpot(new_spot_id_counter)
            new_spots.append(new_spot)
            spots_of_the_cc = {new_spot} # {f"new_spot_{new_spot_id_counter}"} 
            new_spot_id_counter += 1

        elif len(spots_of_the_cc) > 1:
            # more than one spot in the cc 
            logging.getLogger("PPanGGOLiN").info('Some RGPs of the input organism '
                                                 f"are connected to {len(spots_of_the_cc)} original spots of the pangenome.")
        
        input_rgps_of_the_cc = set()
        for node in comp:
            if node in input_org_node_to_rgps:
                input_rgps_of_the_cc |= input_org_node_to_rgps[node]
                
                if write_graph_flag:
                    graph_spot.nodes[node]["spots"] = spots_of_the_cc

                    graph_spot.nodes[node]["spot_id"] = ';'.join((str(spot) for spot in spots_of_the_cc))
                    graph_spot.nodes[node]["includes_RGPs_from_the_input_organism"] = True

        for spot in spots_of_the_cc:
            spot.add_regions(input_rgps_of_the_cc)

        input_rgp_to_spots.update({rgp:spots_of_the_cc for rgp in input_rgps_of_the_cc})
    
    if write_graph_flag:
        for node in graph_spot.nodes:
            del graph_spot.nodes[node]["spots"]

        write_spot_graph(graph_spot, output, graph_formats, file_basename='projected_spotGraph')

    
    write_rgp_to_spot_table(input_rgp_to_spots, output=output, filename='input_organism_rgp_to_spot.tsv')


    new_spots = {spot for spots in input_rgp_to_spots.values() for spot in spots if type(spot) == NewSpot}

    if new_spots:
        logging.getLogger('PPanGGOLiN').info(f'{len(new_spots)} new spots have been created for the input genome.')
        summarize_spots(new_spots, output, compress = False, file_name="new_spots_summary.tsv")


    return input_rgp_to_spots


def project_and_write_modules(pangenome:Pangenome, input_organism: Organism, output:Path, compress:bool=False):
    """
    Write a tsv file providing association between modules and the input organism

    :param output: Path to output directory
    :param compress: Compress the file in .gz
    """

    output_file = output / "modules_in_input_organism.tsv"

    input_organism_families = input_organism.families
    counter = 0
    with write_compressed_or_not(output_file, compress) as fout:
        fout.write("module_id\torganism\tcompletion\n")

        for mod in pangenome.modules:
            module_in_input_organism = any((fam in input_organism_families for fam in mod.families))

            if module_in_input_organism:
                counter += 1

                completion = round(len(input_organism.families & mod.families) / len(mod.families), 2)
                fout.write(f"module_{mod.ID}\t{input_organism.name}\t{completion}\n")

    logging.getLogger().info(f"{counter} modules have been projected to the input genomes.")

    logging.getLogger().info(
        f"Projected modules have been written in: '{output_file}'")


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("projection", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_projection(parser)
    return parser


def parser_projection(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of annotate command

    :param parser: parser for annotate argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=False, type=Path, help="The pangenome.h5 file")
    
    required.add_argument("-n", '--organism_name', required=False, type=str,
                        help="Name of the input_organism whose genome is being annotated with the provided pangenome.")
    
    required.add_argument('--fasta_file', required=False, type=Path,
                        help="The filepath of the genomic sequence(s) in FASTA format if the genome to annotate. "
                        "(Fasta file can be compressed with gzip)")

    required.add_argument('--annot_file', required=False, type=Path,
                        help="The filepath of the annotations in GFF/GBFF format for the genome to annotate with the provided pangenome. "
                        "(Annotation file can be compressed with gzip)")

    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument('-o', '--output', required=False, type=Path,
                          default="ppanggolin_output" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                      time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory")
    
    optional.add_argument("--tmpdir", required=False, type=Path, default=Path(tempfile.gettempdir()),
                        help="directory for storing temporary files")
    
    optional.add_argument('--no_defrag', required=False, action="store_true",
                          help="DO NOT Realign gene families to link fragments with"
                               "their non-fragmented gene family. (default: False)")
    
    optional.add_argument('--identity', required=False, type=restricted_float, default=0.5,
                          help="min identity percentage threshold")
    
    optional.add_argument('--coverage', required=False, type=restricted_float, default=0.8,
                          help="min coverage percentage threshold")
    
    optional.add_argument("--translation_table", required=False, default="11",
                          help="Translation table (genetic code) to use.")
    
    optional.add_argument("--use_pseudo", required=False, action="store_true",
                          help="In the context of provided annotation, use this option to read pseudogenes. "
                               "(Default behavior is to ignore them)")
    
    optional.add_argument("--spot_graph", required=False, action="store_true",
                        help="Write the spot graph to a file, with pairs of blocks of single copy markers flanking RGPs "
                            "as nodes. This graph can be used to visualize nodes that have RGPs from the input organism.")
    
    optional.add_argument('--graph_formats', required=False, type=str, choices=['gexf', "graphml"], nargs="+",
                          default=['gexf'], help="Format of the output graph.")  
    
    optional.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
