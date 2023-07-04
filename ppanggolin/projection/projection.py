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
from typing import Tuple, Set, Dict, Iterator

# installed libraries
from tqdm import tqdm

# # local libraries
from ppanggolin.annotate.synta import annotate_organism, read_fasta, get_dna_sequence
from ppanggolin.annotate.annotate import read_anno_file
from ppanggolin.pangenome import Pangenome
from ppanggolin.cluster.cluster import infer_singletons
# from ppanggolin.genome import input_organism, Gene, RNA, Contig
from ppanggolin.utils import read_compressed_or_not, write_compressed_or_not, mk_file_name, min_one, restricted_float, mk_outdir
from ppanggolin.align.alignOnPang import get_seq2pang, project_and_write_partition
from ppanggolin.formats.writeSequences import write_gene_sequences_from_annotations
from ppanggolin.formats import check_pangenome_info
# from ppanggolin.formats import write_pangenome
from ppanggolin.RGP.genomicIsland import naming_scheme, compute_org_rgp
# from ppanggolin.formats.readBinaries import retrieve_pangenome_parameters
from ppanggolin.genome import Organism, Gene, RNA, Contig
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.region import Region
from ppanggolin.formats.writeFlat import write_flat_files

def annotate_input_genes_with_pangenome_families(pangenome, input_organism,  output, cpu,  no_defrag, identity, coverage, tmpdir,
                                        disable_bar, translation_table, ):
    
    """

    """

    seq_fasta_file = output / f"{input_organism.name}.fasta"

    with open(seq_fasta_file, "w") as fh_out_faa:
        write_gene_sequences_from_annotations(input_organism, fh_out_faa, seq_attr_to_write="dna",
                                            disable_bar=disable_bar)
    
    # get corresponding gene families
    new_tmpdir = tempfile.TemporaryDirectory(dir=tmpdir, prefix="seq_to_pang_tmpdir_")
    seq_set, _, seqid_to_gene_family = get_seq2pang(pangenome, str(seq_fasta_file), str(output), new_tmpdir, cpu, no_defrag, identity=identity,
                                coverage=coverage, is_protein=False, translation_table=translation_table)
    
    # this function only write the seqid and partition associated in a file
    project_and_write_partition(seqid_to_gene_family, seq_set, str(output))

    # Add gene of the input organism in the associated gene family
    # when a gene is not associated with any gene family, a new family is created.
    lonely_gene = 0
    for gene in input_organism.genes:
        try:
            gene_family = seqid_to_gene_family[gene.ID]
            gene_family.add_gene(gene)

        except KeyError:
            # add a new gene family
            new_gene_family = pangenome.add_gene_family(gene.ID)
            new_gene_family.add_gene(gene)
            new_gene_family.add_partition("Cloud")
            lonely_gene += 1

    logging.getLogger().info(f"The input organisms have {lonely_gene}/{input_organism.number_of_genes()} " 
                             "genes that do not cluster with any of the gene families of the pangenome.")


def compute_RGP(pangenome: Pangenome, input_organism: Organism, persistent_penalty: int, variable_gain: int,
                min_length: int, min_score: int, dup_margin: float,
                disable_bar: bool) -> None:
    """
    Compute Regions of Genomic Plasticity (RGP) for the given pangenome and input organism.

    :param pangenome: The pangenome object.
    :param input_organism: The input organism for which to compute RGPs.
    :param persistent_penalty: Penalty score to apply to persistent genes.
    :param variable_gain: Gain score to apply to variable genes.
    :param min_length: Minimum length (bp) of a region to be considered as RGP.
    :param min_score: Minimal score required for considering a region as RGP.
    :param dup_margin: Minimum ratio of organisms in which a family must have multiple genes to be considered duplicated.
    :param disable_bar: Flag to disable the progress bar.

    :return: None
    """
    logging.getLogger().info("Detecting multigenic families...")
    multigenics = pangenome.get_multigenics(dup_margin)

    logging.getLogger().info("Computing Regions of Genomic Plasticity...")
    name_scheme = naming_scheme(pangenome)

    rgps = compute_org_rgp(input_organism, multigenics, persistent_penalty, variable_gain, min_length,
                    min_score, naming=name_scheme, disable_bar=disable_bar)

    print("Found RGPS ", len(rgps))
    return rgps


def write_predicted_regions(regions : Set[Region], output:Path, compress=False):
    """
    Write the file providing information about RGP content

    :param output: Path to output directory
    :param compress: Compress the file in .gz
    """
    fname = output / "plastic_regions.tsv"
    with write_compressed_or_not(fname, compress) as tab:
        tab.write("region\torganism\tcontig\tstart\tstop\tgenes\tcontigBorder\twholeContig\n")
        regions = sorted(regions, key=lambda x: (x.organism.name, x.contig.name, x.start))
        for region in regions:
            tab.write('\t'.join(map(str, [region.name, region.organism, region.contig, region.start, region.stop,
                                          len(region.genes), region.is_contig_border, region.is_whole_contig])) + "\n")
            

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
    



def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    

    output_dir = Path(args.output)
    mk_outdir(output_dir, args.force)
    
    # TODO check that the provided input_organism name is not found in pangenome 
    # if so add a warning or error

    # TODO some params are no keep in pangenome... like use_pseudo. what to do?
    logging.getLogger().info('Retrieving pangenome parameters from the provided pangenome file.')

    
    if args.annot_file is not None:
        # read_annotations(pangenome, args.anno, cpu=args.cpu, pseudo=args.use_pseudo, disable_bar=args.disable_prog_bar)
        input_organism, has_sequence = read_anno_file(organism_name = args.organism_name, 
                       filename=args.annot_file,
                        circular_contigs=[],
                        pseudo=False) 
        
        if not has_sequence:
            if args.fasta_file:
                retrieve_gene_sequences_from_fasta_file(input_organism, args.fasta_file)
            else:
                raise Exception("The gff/gbff provided did not have any sequence information, "
                            "Thus, we do not have the information we need to continue the projection.")

    elif args.fasta_file is not None:
        input_organism = annotate_organism(org_name=args.organism_name, file_name = args.fasta_file, circular_contigs=[], tmpdir=args.tmpdir,
                      code = args.annotate.translation_table, norna=args.annotate.norna, kingdom = args.annotate.kingdom,
                      overlap=args.annotate.allow_overlap, procedure=args.annotate.prodigal_procedure)

    else:
        raise Exception("At least one of --fasta_file or --anno_file must be given")


    # load pangenome

    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)

    check_pangenome_info(pangenome, need_annotations=True, need_families=True, disable_bar=args.disable_prog_bar, need_rgp=args.predict_rgp)

    # Add input organism in pangenome. This temporary as pangenome is not going to be written.
    pangenome.add_organism(input_organism)

    annotate_input_genes_with_pangenome_families(pangenome, input_organism=input_organism, output=output_dir, cpu=args.cluster.cpu, 
                                        no_defrag = args.cluster.no_defrag, identity = args.cluster.identity, coverage = args.cluster.coverage, tmpdir=args.tmpdir,
                                        disable_bar=args.disable_prog_bar, translation_table = args.annotate.translation_table )
    

    if args.predict_rgp:
        
        rgps = compute_RGP(pangenome, input_organism,  persistent_penalty=args.rgp.persistent_penalty, variable_gain=args.rgp.variable_gain,
                min_length=args.rgp.min_length, min_score=args.rgp.min_score, dup_margin=args.rgp.dup_margin,
                disable_bar=args.disable_prog_bar)
        
        write_predicted_regions(rgps, output=output_dir)

        
    # write_flat_files_for_input_genome(input_organism)


# def write_flat_files_for_input_genome(input_organism):
    
#     # create a pangenome object with only the input organism
#     pangenome = Pangenome()

#     pangenome.add_organism(input_organism)
    

#     write_flat_files(pangenome, regions=True)


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
    required.add_argument('-p', '--pangenome', required=False, type=str, help="The pangenome.h5 file")
    
    required.add_argument('--organism_name', required=False, type=str,
                        help="Name of the input_organism whose genome is being annotated with the provided pangenome.")
    
    required.add_argument('--fasta_file', required=False, type=str,
                        help="The filepath of the genomic sequence(s) in FASTA format for the projected genome. "
                        "(Fasta file can be compressed with gzip)")

    required.add_argument('--annot_file', required=False, type=str,
                        help="The filepath of the annotations in GFF/GBFF format for the projected genome. "
                        "(Annotation file can be compressed with gzip)")



    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('-o', '--output', required=False, type=str,
                          default="ppanggolin_output" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                      time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory")
    
    optional.add_argument('--predict_rgp', required=False, action='store_true', default=False,
                          help="Predict rgp on the input genome.")
    optional.add_argument('--project_modules', required=False, action='store_true', default=False,
                          help="Predict rgp on the input genome.")
    # optional.add_argument("--basename", required=False, default="pangenome", help="basename for the output file")
    
    # annotate = parser.add_argument_group(title="Annotation arguments")


    # annotate.add_argument('--allow_overlap', required=False, action='store_true', default=False,
    #                       help="Use to not remove genes overlapping with RNA features.")
    # annotate.add_argument("--norna", required=False, action="store_true", default=False,
    #                       help="Use to avoid annotating RNA features.")
    # annotate.add_argument("--kingdom", required=False, type=str.lower, default="bacteria",
    #                       choices=["bacteria", "archaea"],
    #                       help="Kingdom to which the prokaryota belongs to, "
    #                            "to know which models to use for rRNA annotation.")
    # annotate.add_argument("--translation_table", required=False, type=int, default=11,
    #                       help="Translation table (genetic code) to use.")

    # annotate.add_argument("--prodigal_procedure", required=False, type=str.lower, choices=["single", "meta"],
    #                       default=None, help="Allow to force the prodigal procedure. "
    #                                          "If nothing given, PPanGGOLiN will decide in function of contig length")
    # annotate.add_argument("--use_pseudo", required=False, action="store_true",
    #                     help="In the context of provided annotation, use this option to read pseudogenes. "
    #                         "(Default behavior is to ignore them)")

    # cluster = parser.add_argument_group(title="Clustering arguments")
    # cluster.add_argument('--no_defrag', required=False, action="store_true",
    #                       help="DO NOT Realign gene families to link fragments with"
    #                            "their non-fragmented gene family.")
    # cluster.add_argument('--identity', required=False, type=float, default=0.5,
    #                       help="min identity percentage threshold")
    # cluster.add_argument('--coverage', required=False, type=float, default=0.8,
    #                       help="min coverage percentage threshold")
    
