#!/usr/bin/env python3
# coding:utf-8

# default libraries
import os
import time
import argparse
import logging
from collections import defaultdict
from collections.abc import Callable

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mk_file_name, mk_outdir, check_option_workflow, add_step_specific_args
from ppanggolin.annotate.annotate import annotate_pangenome, read_annotations, get_gene_sequences_from_fastas
from ppanggolin.cluster.cluster import clustering, read_clustering
from ppanggolin.graph.makeGraph import compute_neighbors_graph
from ppanggolin.nem.rarefaction import make_rarefaction_curve
from ppanggolin.nem.partition import partition
from ppanggolin.formats.writeBinaries import write_pangenome
from ppanggolin.formats.writeFlat import write_flat_files
from ppanggolin.figures.ucurve import draw_ucurve
from ppanggolin.figures.tile_plot import draw_tile_plot
from ppanggolin.figures.draw_spot import draw_spots
from ppanggolin.info.info import print_info
from ppanggolin.RGP.genomicIsland import predict_rgp
from ppanggolin.RGP.spot import predict_hotspots
from ppanggolin.mod.module import predict_modules

"""a global workflow that does everything in one go."""

def launch_workflow(args: argparse.Namespace, subcomamand_parser: Callable, panrgp: bool = True, panmodule: bool = True):
    """
    Unified function to launch ppanggolin workflow.

    :param args: All arguments provide by command line
    :param subcomamand_parser: Subparser function of the subcommand
    :param panrgp: Flag to launch RGP prediction
    :param panmodule: Flag to launch module prediction
    """

    check_option_workflow(args)
    
    pangenome = Pangenome()

    filename = mk_file_name(args.basename, args.output, args.force)

    writing_time, anno_time, clust_time, mod_time, desc_time = (None, None, None, None, None)

    if args.anno:  # if the annotations are provided, we read from it

        start_anno = time.time()
        read_annotations(pangenome, args.anno, pseudo=args.annotate.use_pseudo,
                        cpu=args.annotate.cpu, disable_bar=args.disable_prog_bar)
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
                            infer_singleton=args.cluster.infer_singletons)

        elif args.clusters is None:  # we should have the sequences here.

            clustering(pangenome, tmpdir=args.tmpdir, cpu=args.cluster.cpu, force=args.force, disable_bar=args.disable_prog_bar, 
                        defrag=not args.cluster.no_defrag, code=args.cluster.translation_table,
                        coverage=args.cluster.coverage, identity=args.cluster.identity, mode=args.cluster.mode) 
        clust_time = time.time() - start_clust

    elif args.fasta is not None:

        start_anno = time.time()
        annotate_pangenome(pangenome, args.fasta, tmpdir=args.tmpdir, cpu=args.annotate.cpu, disable_bar=args.disable_prog_bar,
                        procedure=args.annotate.prodigal_procedure,
                        translation_table=args.annotate.translation_table, kingdom=args.annotate.kingdom, norna=args.annotate.norna,
                        overlap=args.annotate.allow_overlap, contig_filter=args.annotate.contig_filter)
        anno_time = time.time() - start_anno

        start_writing = time.time()
        write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
        writing_time = time.time() - start_writing

        start_clust = time.time()
        clustering(pangenome, tmpdir=args.tmpdir, cpu=args.cluster.cpu, force=args.force, disable_bar=args.disable_prog_bar, 
                    defrag=not args.cluster.no_defrag, code=args.cluster.translation_table,
                    coverage=args.cluster.coverage, identity=args.cluster.identity, mode=args.cluster.mode) 
        clust_time = time.time() - start_clust

    write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)

    start_graph = time.time()
    compute_neighbors_graph(pangenome, args.graph.remove_high_copy_number, args.force, disable_bar=args.disable_prog_bar)

    graph_time = time.time() - start_graph

    start_part = time.time()
    partition(pangenome,  
                tmpdir=args.tmpdir,
                outputdir = args.output,
                beta = args.partition.beta,
                sm_degree = args.partition.max_degree_smoothing,
                free_dispersion = args.partition.free_dispersion,
                chunk_size = args.partition.chunk_size,
                kval = args.partition.nb_of_partitions,
                krange = args.partition.krange,
                icl_margin = args.partition.ICL_margin,
                draw_icl = args.partition.draw_ICL,
                seed = args.partition.seed,
                keep_tmp_files = args.partition.keep_tmp_files,
                cpu = args.partition.cpu,
                force = args.force,
                disable_bar=args.disable_prog_bar)
    part_time = time.time() - start_part

    start_writing = time.time()
    write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
    writing_time = writing_time + time.time() - start_writing

    if panrgp:
        start_regions = time.time()
        predict_rgp(pangenome, persistent_penalty=args.rgp.persistent_penalty, variable_gain=args.rgp.variable_gain,
                    min_length=args.rgp.min_length, min_score=args.rgp.min_score, dup_margin=args.rgp.dup_margin,
                    force=args.force, disable_bar=args.disable_prog_bar)

        regions_time = time.time() - start_regions

    
        start_spots = time.time()
        predict_hotspots(pangenome, args.output, force=args.force, spot_graph=args.spot.spot_graph,
                        overlapping_match=args.spot.overlapping_match, set_size=args.spot.set_size,
                        exact_match=args.spot.exact_match_size, disable_bar=args.disable_prog_bar)
        spot_time = time.time() - start_spots

    if panmodule:
        start_mods = time.time()
        predict_modules(pangenome=pangenome, tmpdir=args.tmpdir, cpu=args.module.cpu, dup_margin=args.module.dup_margin, size=args.module.size,
                    min_presence=args.module.min_presence, transitive=args.module.transitive, jaccard=args.module.jaccard,
                    force=args.force,
                    disable_bar=args.disable_prog_bar)

        mod_time = time.time() - start_mods

    start_writing = time.time()
    write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
    writing_time = writing_time + time.time() - start_writing

    if not args.only_pangenome:

        if panrgp:
            start_spot_drawing = time.time()
            mk_outdir(args.output + '/spot_figures', force=True)
            draw_spots(pangenome=pangenome, output=args.output + '/spot_figures', spot_list='all',
                    disable_bar=args.disable_prog_bar)
            spot_time = spot_time + time.time() - start_spot_drawing

        if args.rarefaction:
            make_rarefaction_curve(pangenome=pangenome, output=args.output, tmpdir=args.tmpdir, 
                                    beta=args.rarefaction.beta,
                                    depth=args.rarefaction.depth,
                                    min_sampling=args.rarefaction.min,
                                    max_sampling=args.rarefaction.max,
                                    sm_degree=args.rarefaction.max_degree_smoothing, 
                                    free_dispersion=args.rarefaction.free_dispersion,
                                    chunk_size=args.rarefaction.chunk_size,
                                    kval=args.rarefaction.nb_of_partitions,
                                    krange=args.rarefaction.krange,
                                    seed=args.rarefaction.seed,
                                    kestimate=args.rarefaction.reestimate_K,
                                    soft_core=args.rarefaction.soft_core, 
                                    cpu=args.rarefaction.cpu, disable_bar=args.disable_prog_bar)

        if 1 < len(pangenome.organisms) < 5000:
            draw_tile_plot(pangenome, args.output, nocloud=False if len(pangenome.organisms) < 500 else True)
        draw_ucurve(pangenome, args.output)

        start_desc = time.time()

        write_out_arguments = ["csv", "Rtab", "gexf", "light_gexf", "projection", "stats", 'json']

        # Check that we don't ask write to output something not computed.
        borders, spots, regions, spot_modules, modules = (False, False, False, False, False)

        if panmodule:
            modules = args.write.modules
            write_out_arguments.append('modules')

        if panrgp:
            borders, spots, regions = (args.write.borders, args.write.spots, args.write.regions)
            write_out_arguments +=  ["borders", "spots", 'regions']
        
        if panmodule and panrgp:
            spot_modules = args.write.spot_modules
            write_out_arguments.append('spot_modules')

        # check that at least one output file is requested. if not write is not call.
        if any(( getattr(args.write,arg) is True for arg in  write_out_arguments)):
            # some parameters are set to false because they have not been computed in this workflow
            write_flat_files(pangenome, args.output, cpu=args.write.cpu,  disable_bar=args.disable_prog_bar, 
                            soft_core=args.write.soft_core, dup_margin=args.write.dup_margin,
                            csv=args.write.csv, gene_pa=args.write.Rtab, gexf=args.write.gexf, light_gexf=args.write.light_gexf, projection=args.write.projection,
                            stats=args.write.stats, json=args.write.json, partitions=args.write.partitions, 
                            families_tsv=args.write.families_tsv, 
                            compress=args.write.compress, 
                            spot_modules=spot_modules, regions=regions, modules=modules, spots=spots, borders=borders)
        else:
            logging.getLogger().info(f'No flat file output has been requested in config file. Writing output flat file is skipped.')

        desc_time = time.time() - start_desc

    logging.getLogger().info(f"Annotation took : {round(anno_time, 2)} seconds")
    logging.getLogger().info(f"Clustering took : {round(clust_time, 2)} seconds")
    logging.getLogger().info(f"Building the graph took : {round(graph_time, 2)} seconds")
    logging.getLogger().info(f"Partitioning the pangenome took : {round(part_time, 2)} seconds")

    if panrgp:
        logging.getLogger().info(f"Predicting RGP took : {round(regions_time, 2)} seconds")
        logging.getLogger().info(f"Gathering RGP into spots took : {round(spot_time, 2)} seconds")
    if panmodule:
        logging.getLogger().info(f"Predicting modules took : {round(mod_time, 2)} seconds")
    
    logging.getLogger().info(f"Writing the pangenome data in HDF5 took : {round(writing_time, 2)} seconds")

    if not args.only_pangenome:
        logging.getLogger().info(f"Writing descriptive files for the pangenome took : {round(desc_time, 2)} seconds")
    
    print_info(filename, content=True)


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """

    launch_workflow(args, subparser,  panrgp=True, panmodule=True)


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

    optional.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")


    add_step_specific_args(optional)


    return parser


