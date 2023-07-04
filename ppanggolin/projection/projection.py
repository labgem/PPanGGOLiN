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


# installed libraries
from tqdm import tqdm

# # local libraries
from ppanggolin.annotate.synta import annotate_organism, read_fasta, get_dna_sequence
from ppanggolin.annotate.annotate import read_anno_file
from ppanggolin.pangenome import Pangenome
# from ppanggolin.genome import input_organism, Gene, RNA, Contig
from ppanggolin.utils import read_compressed_or_not, mk_file_name, min_one, restricted_float, mk_outdir
from ppanggolin.align.alignOnPang import get_seq2pang, project_partition
from ppanggolin.formats.writeSequences import write_gene_sequences_from_annotations
from ppanggolin.formats import check_pangenome_info
# from ppanggolin.formats import write_pangenome
from ppanggolin.RGP.genomicIsland import naming_scheme, compute_org_rgp
from ppanggolin.formats.readBinaries import retrieve_pangenome_parameters

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
    seq_set, _, seq2pan = get_seq2pang(pangenome, str(seq_fasta_file), str(output), new_tmpdir, cpu, no_defrag, identity=identity,
                                coverage=coverage, is_protein=False, translation_table=translation_table)
    
    project_partition(seq2pan, seq_set, str(output))


def compute_RGP(pangenome, input_organism,  dup_margin, persistent_penalty, variable_gain, min_length, min_score):
    
    ## Computing RGPs ##   
    logging.getLogger().info("Detecting multigenic families...")
    multigenics = pangenome.get_multigenics(dup_margin)

    logging.getLogger().info("Compute Regions of Genomic Plasticity ...")
    name_scheme = naming_scheme(pangenome)

    compute_org_rgp(input_organism, multigenics, persistent_penalty, variable_gain, min_length,
                                                min_score, naming=name_scheme)
    



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

    step_to_params = retrieve_pangenome_parameters(args.pangenome)
    annotation_params_str = "  ".join([f"{param}={value}"  for param, value in step_to_params["annotation"].items()]) 
    logging.getLogger().debug(f'annotation params {annotation_params_str}' )
    
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
    check_pangenome_info(pangenome, need_annotations=True, need_families=True, disable_bar=args.disable_prog_bar)

    annotate_input_genes_with_pangenome_families(pangenome, input_organism=input_organism, output=output_dir, cpu=args.cluster.cpu, 
                                        no_defrag = args.cluster.no_defrag, identity = args.cluster.identity, coverage = args.cluster.coverage, tmpdir=args.tmpdir,
                                        disable_bar=args.disable_prog_bar, translation_table = args.annotate.translation_table )
    

    # compute_RGP(pangenome, input_organism,  dup_margin=0.05, persistent_penalty=3, variable_gain=1, min_length=3000, min_score=4)


    
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

    # required.add_argument('--fasta', required=False, type=str,
    #                       help="A tab-separated file listing the input_organism names, and the fasta filepath of its genomic "
    #                            "sequence(s) (the fastas can be compressed with gzip). One line per input_organism.")
    
    # required.add_argument('--anno', required=False, type=str,
    #                       help="A tab-separated file listing the input_organism names, and the gff/gbff filepath of its "
    #                            "annotations (the files can be compressed with gzip). One line per input_organism. "
    #                            "If this is provided, those annotations will be used.")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('-o', '--output', required=False, type=str,
                          default="ppanggolin_output" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                      time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory")
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
    
