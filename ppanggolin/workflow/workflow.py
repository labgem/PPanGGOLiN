#!/usr/bin/env python3
# coding:utf-8

# default libraries
import os
import time
import argparse
from collections import defaultdict
import logging

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mk_file_name, check_option_workflow, get_cmd_args_from_config, parse_config_file
from ppanggolin.annotate.annotate import annotate_pangenome, read_annotations, get_gene_sequences_from_fastas, parser_annot
from ppanggolin.cluster.cluster import clustering, read_clustering, parser_clust
from ppanggolin.graph.makeGraph import compute_neighbors_graph, parser_graph
from ppanggolin.nem.rarefaction import make_rarefaction_curve, parser_rarefaction
from ppanggolin.nem.partition import partition, parser_partition
from ppanggolin.formats.writeBinaries import write_pangenome
from ppanggolin.formats.writeFlat import write_flat_files, parser_flat
from ppanggolin.figures.ucurve import draw_ucurve
from ppanggolin.figures.tile_plot import draw_tile_plot
from ppanggolin.info.info import print_info


""" a global workflow that does everything in one go. """

def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """

    check_option_workflow(args)

    if args.config:
        config = parse_config_file(args.config)

        # convert config dict to defaultdict 
        config = defaultdict(dict, config)
    else:
        config = defaultdict(dict)

    general_params = ['help', 'fasta', 'clusters', 'anno', 'cpu', "output"]
    
    annotate_args = get_cmd_args_from_config("annotate", parser_annot, config['annotate'], general_params)
    cluster_args = get_cmd_args_from_config("cluster", parser_clust, config['cluster'], general_params)
    graph_args = get_cmd_args_from_config("graph", parser_graph, config['graph'], general_params)
    partition_args = get_cmd_args_from_config("partition", parser_partition, config['partition'], general_params)
    rarefaction_args = get_cmd_args_from_config("rarefaction", parser_rarefaction, config['rarefaction'], general_params)
    write_args = get_cmd_args_from_config("write", parser_flat, config['write'], general_params)

    pangenome = Pangenome()

    filename = mk_file_name(args.basename, args.output, args.force)

    if args.anno:  # if the annotations are provided, we read from it

        read_annotations(pangenome, args.anno, pseudo=annotate_args.use_pseudo,
                        cpu=args.cpu, disable_bar=args.disable_prog_bar)

        write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)

        if args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is None:
            raise Exception("The gff/gbff provided did not have any sequence informations, "
                            "you did not provide clusters and you did not provide fasta file. "
                            "Thus, we do not have the information we need to continue the analysis.")

        elif args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is not None:
            get_gene_sequences_from_fastas(pangenome, args.fasta)

        if args.clusters is not None:
            read_clustering(pangenome, args.clusters, disable_bar=args.disable_prog_bar, 
                            infer_singleton=cluster_args.infer_singleton)

        elif args.clusters is None:  # we should have the sequences here.

            clustering(pangenome, tmpdir=args.tmpdir, cpu=args.cpu, force=args.force, disable_bar=args.disable_prog_bar, 
                        defrag=not cluster_args.no_defrag, code=cluster_args.translation_table,
                        coverage=cluster_args.coverage, identity=cluster_args.identity, mode=cluster_args.mode) 
                        

    elif args.fasta is not None:
        pangenome = Pangenome()
        # annotate_pangenome(pangenome, args.fasta, args.tmpdir, args.cpu, contig_filter=args.contig_filter,
        #                    disable_bar=args.disable_prog_bar)
                           
        annotate_pangenome(pangenome, args.fasta, tmpdir=args.tmpdir, cpu=args.cpu, disable_bar=args.disable_prog_bar,
                        procedure=annotate_args.prodigal_procedure,
                        translation_table=annotate_args.translation_table, kingdom=annotate_args.kingdom, norna=annotate_args.norna,
                        overlap=annotate_args.overlap, contig_filter=annotate_args.contig_filter)

        write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)

        clustering(pangenome, tmpdir=args.tmpdir, cpu=args.cpu, force=args.force, disable_bar=args.disable_prog_bar, 
                    defrag=not cluster_args.no_defrag, code=cluster_args.translation_table,
                    coverage=cluster_args.coverage, identity=cluster_args.identity, mode=cluster_args.mode) 


    compute_neighbors_graph(pangenome, graph_args.remove_high_copy_number, args.force, disable_bar=args.disable_prog_bar)

    partition(pangenome,  
                tmpdir=args.tmpdir,
                outputdir = args.output,
                beta = partition_args.beta,
                sm_degree = partition_args.max_degree_smoothing,
                free_dispersion = partition_args.free_dispersion,
                chunk_size = partition_args.chunk_size,
                kval = partition_args.nb_of_partitions,
                krange = partition_args.krange,
                icl_margin = partition_args.ICL_margin,
                draw_icl = partition_args.draw_ICL,
                seed = partition_args.seed,
                keep_tmp_files = partition_args.keep_tmp_files,
                cpu = args.cpu,
                force = args.force,
                disable_bar=args.disable_prog_bar)

    write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)

    if not args.only_pangenome:
        if args.rarefaction:

            make_rarefaction_curve(pangenome=pangenome, output=args.output, tmpdir=args.tmpdir, 
                                    beta=rarefaction_args.beta,
                                    depth=rarefaction_args.depth, 
                                    min_sampling=rarefaction_args.min,
                                    max_sampling=rarefaction_args.max,
                                    sm_degree=rarefaction_args.max_degree_smoothing, 
                                    free_dispersion=rarefaction_args.free_dispersion,
                                    chunk_size=rarefaction_args.chunk_size,
                                    kval=rarefaction_args.nb_of_partitions,
                                    krange=rarefaction_args.krange,
                                    seed=rarefaction_args.seed,
                                    kestimate=rarefaction_args.reestimate_K,
                                    soft_core=rarefaction_args.soft_core, 
                                    cpu=args.cpu, disable_bar=args.disable_prog_bar)

        if 1 < len(pangenome.organisms) < 5000:
            draw_tile_plot(pangenome, args.output, nocloud=False if len(pangenome.organisms) < 500 else True)

        draw_ucurve(pangenome, args.output)

        # check that at least one output file is requested. if not write is not call.
        write_out_arguments = ["csv", "Rtab", "gexf", "light_gexf", "projection", "stats", 'json', "partitions"]

        if any(( getattr(write_args,arg) for arg in  write_out_arguments)):
            # some parameters are set to false because they have not been computed in this workflow
            write_flat_files(pangenome, args.output, cpu=args.cpu,  disable_bar=args.disable_prog_bar, 
                            soft_core=write_args.soft_core, dup_margin=write_args.dup_margin,
                            csv=write_args.csv, gene_pa=write_args.Rtab, gexf=write_args.gexf, light_gexf=write_args.light_gexf, projection=write_args.projection,
                            stats=write_args.stats, json=write_args.json, partitions=write_args.partitions, 
                            families_tsv=write_args.families_tsv, 
                            compress=write_args.compress, 
                            spot_modules=False, regions=False, modules=False, spots=False, borders=False)
        else:
            logging.getLogger().info(f'No flat file has been requested in config file. Writing output flat file is skipped.')


    print_info(filename, content=True)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("workflow", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title="Input arguments", description="The possible input arguments :")
    required.add_argument('--fasta', required=False, type=str,
                          help="A tab-separated file listing the organism names, "
                               "and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). "
                               "One line per organism. This option can be used alone.")
    required.add_argument('--anno', required=False, type=str,
                          help="A tab-separated file listing the organism names, and the gff filepath "
                               "of its annotations (the gffs can be compressed). One line per organism. "
                               "This option can be used alone IF the fasta sequences are in the gff files, "
                               "otherwise --fasta needs to be used.")
    required.add_argument("--clusters", required=False, type=str,
                          help="a tab-separated file listing the cluster names, the gene IDs, "
                               "and optionally whether they are a fragment or not.")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('-o', '--output', required=False, type=str,
                          default="ppanggolin_output" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                      time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory")
    optional.add_argument("--basename", required=False, default="pangenome", help="basename for the output file")
    
    optional.add_argument("--rarefaction", required=False, action="store_true",
                          help="Use to compute the rarefaction curves (WARNING: can be time consuming)")
    
    optional.add_argument("--only_pangenome", required=False, action="store_true",
                          help="Only generate the HDF5 pangenome file")
    
    optional.add_argument("--config", required=False, type=open, help="Config file.")
   
    
    return parser
