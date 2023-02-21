#!/usr/bin/env python3
# coding:utf-8

# default libraries
import os
import time
import argparse
import logging
from collections import defaultdict

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mk_file_name, mk_outdir, check_option_workflow, get_cmd_args_from_config, parse_config_file
from ppanggolin.annotate.annotate import annotate_pangenome, read_annotations, get_gene_sequences_from_fastas, parser_annot
from ppanggolin.cluster.cluster import clustering, read_clustering, parser_clust
from ppanggolin.graph.makeGraph import compute_neighbors_graph, parser_graph
from ppanggolin.nem.rarefaction import make_rarefaction_curve, parser_rarefaction
from ppanggolin.nem.partition import partition, parser_partition
from ppanggolin.formats.writeBinaries import write_pangenome
from ppanggolin.formats.writeFlat import write_flat_files, parser_flat
from ppanggolin.figures.ucurve import draw_ucurve
from ppanggolin.figures.tile_plot import draw_tile_plot
from ppanggolin.figures.draw_spot import draw_spots
from ppanggolin.info.info import print_info
from ppanggolin.RGP.genomicIsland import predict_rgp, parser_rgp
from ppanggolin.RGP.spot import predict_hotspots, parser_spot
from ppanggolin.mod.module import predict_modules, parser_module

"""a global workflow that does everything in one go."""


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

    rgp_args = get_cmd_args_from_config("rgp", parser_rgp, config['rgp'], general_params)
    spots_args = get_cmd_args_from_config("spot", parser_spot, config['spot'], general_params)
    module_args = get_cmd_args_from_config("module", parser_module, config['module'], general_params)
    

    pangenome = Pangenome()

    filename = mk_file_name(args.basename, args.output, args.force)

    writing_time, anno_time, clust_time, mod_time, desc_time = (None, None, None, None, None)

    if args.anno:  # if the annotations are provided, we read from it

        start_anno = time.time()
        read_annotations(pangenome, args.anno, pseudo=annotate_args.use_pseudo,
                        cpu=args.cpu, disable_bar=args.disable_prog_bar)
        anno_time = time.time() - start_anno

        start_writing = time.time()
        write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
        writing_time = time.time() - start_writing

        if args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is None:
            raise Exception("The gff/gbff provided did not have any sequence informations, "
                            "you did not provide clusters and you did not provide fasta file. "
                            "Thus, we do not have the information we need to continue the analysis.")

        elif args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is not None:
            get_gene_sequences_from_fastas(pangenome, args.fasta)

        start_clust = time.time()
        if args.clusters is not None:
            read_clustering(pangenome, args.clusters, disable_bar=args.disable_prog_bar, 
                            infer_singleton=cluster_args.infer_singleton)

        elif args.clusters is None:  # we should have the sequences here.

            clustering(pangenome, tmpdir=args.tmpdir, cpu=args.cpu, force=args.force, disable_bar=args.disable_prog_bar, 
                        defrag=not cluster_args.no_defrag, code=cluster_args.translation_table,
                        coverage=cluster_args.coverage, identity=cluster_args.identity, mode=cluster_args.mode) 
        clust_time = time.time() - start_clust

    elif args.fasta is not None:

        start_anno = time.time()
        annotate_pangenome(pangenome, args.fasta, tmpdir=args.tmpdir, cpu=args.cpu, disable_bar=args.disable_prog_bar,
                        procedure=annotate_args.prodigal_procedure,
                        translation_table=annotate_args.translation_table, kingdom=annotate_args.kingdom, norna=annotate_args.norna,
                        overlap=annotate_args.overlap, contig_filter=annotate_args.contig_filter)
        anno_time = time.time() - start_anno

        start_writing = time.time()
        write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
        writing_time = time.time() - start_writing

        start_clust = time.time()
        clustering(pangenome, tmpdir=args.tmpdir, cpu=args.cpu, force=args.force, disable_bar=args.disable_prog_bar, 
                    defrag=not cluster_args.no_defrag, code=cluster_args.translation_table,
                    coverage=cluster_args.coverage, identity=cluster_args.identity, mode=cluster_args.mode) 
        clust_time = time.time() - start_clust

    write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)

    start_graph = time.time()
    compute_neighbors_graph(pangenome, graph_args.remove_high_copy_number, args.force, disable_bar=args.disable_prog_bar)

    graph_time = time.time() - start_graph

    start_part = time.time()
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
    part_time = time.time() - start_part

    start_writing = time.time()
    write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
    writing_time = writing_time + time.time() - start_writing

    start_regions = time.time()
    # predict_rgp(pangenome, disable_bar=args.disable_prog_bar)
    predict_rgp(pangenome, persistent_penalty=rgp_args.persistent_penalty, variable_gain=rgp_args.variable_gain,
                min_length=rgp_args.min_length, min_score=rgp_args.min_score, dup_margin=rgp_args.dup_margin,
                force=args.force, disable_bar=args.disable_prog_bar)

    regions_time = time.time() - start_regions

    
    # predict_hotspots(pangenome, args.output, disable_bar=args.disable_prog_bar)
    start_spots = time.time()
    predict_hotspots(pangenome, args.output, force=args.force, spot_graph=spots_args.spot_graph,
                    overlapping_match=spots_args.overlapping_match, set_size=spots_args.set_size,
                    exact_match=spots_args.exact_match_size, disable_bar=args.disable_prog_bar)
    spot_time = time.time() - start_spots

    start_mods = time.time()
    # predict_modules(pangenome=pangenome, tmpdir=args.tmpdir, cpu=args.cpu, disable_bar=args.disable_prog_bar)
    predict_modules(pangenome=pangenome, tmpdir=args.tmpdir, cpu=args.cpu, dup_margin=module_args.dup_margin, size=module_args.size,
                min_presence=module_args.min_presence, transitive=module_args.transitive, jaccard=module_args.jaccard,
                force=args.force,
                disable_bar=args.disable_prog_bar)

    mod_time = time.time() - start_mods

    start_writing = time.time()
    write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
    writing_time = writing_time + time.time() - start_writing

    if not args.only_pangenome:

        start_spot_drawing = time.time()
        mk_outdir(args.output + '/spot_figures', force=True)
        draw_spots(pangenome=pangenome, output=args.output + '/spot_figures', spot_list='all',
                   disable_bar=args.disable_prog_bar)
        spot_time = spot_time + time.time() - start_spot_drawing

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

        start_desc = time.time()
        # write_flat_files(pangenome, args.output, args.cpu, csv=True, gene_pa=True, gexf=True, light_gexf=True,
        #                  projection=True, stats=True, json=True, partitions=True, regions=True, spots=True,
        #                  borders=True, modules=True, spot_modules=True)

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

        desc_time = time.time() - start_desc

    logging.getLogger().info(f"Annotation took : {round(anno_time, 2)} seconds")
    logging.getLogger().info(f"Clustering took : {round(clust_time, 2)} seconds")
    logging.getLogger().info(f"Building the graph took : {round(graph_time, 2)} seconds")
    logging.getLogger().info(f"Partitioning the pangenome took : {round(part_time, 2)} seconds")
    logging.getLogger().info(f"Predicting RGP took : {round(regions_time, 2)} seconds")
    logging.getLogger().info(f"Gathering RGP into spots took : {round(spot_time, 2)} seconds")
    logging.getLogger().info(f"Predicting modules took : {round(mod_time, 2)} seconds")
    logging.getLogger().info(f"Writing the pangenome data in HDF5 took : {round(writing_time, 2)} seconds")

    if not args.only_pangenome:
        logging.getLogger().info(f"Writing descriptive files for the pangenome took : {round(desc_time, 2)} seconds")
    
    print_info(filename, content=True)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("all", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title="Input arguments", description="The possible input arguments :")
    required.add_argument('--fasta', required=False, type=str,
                          help="A tab-separated file listing the organism names, "
                               "and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). "
                               "One line per organism. This option can be used alone.")
    required.add_argument('--anno', required=False, type=str,
                          help="A tab-separated file listing the organism names, and the gff filepath of "
                               "its annotations (the gffs can be compressed). One line per organism. "
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
    optional.add_argument("--basename", required=False, default="pangenome",
                            help="basename for the output file")
    optional.add_argument("--rarefaction", required=False, action="store_true",
                          help="Use to compute the rarefaction curves (WARNING: can be time consuming)")
    optional.add_argument("--only_pangenome", required=False, action="store_true",
                          help="Only generate the HDF5 pangenome file")
    optional.add_argument("--config", required=False, type=open, 
                        help="Config file in yaml format to launch the different step of "
                             "the workflow with specific arguments.")
   

    return parser
