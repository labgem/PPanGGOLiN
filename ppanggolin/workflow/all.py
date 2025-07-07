#!/usr/bin/env python3

# default libraries
import logging
import os
import time
import argparse
from pathlib import Path
import tempfile

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import (
    mk_file_name,
    mk_outdir,
    check_option_workflow,
    restricted_float,
)
from ppanggolin.annotate.annotate import (
    annotate_pangenome,
    read_annotations,
    get_gene_sequences_from_fastas,
    check_annotate_args,
)
from ppanggolin.cluster.cluster import clustering, read_clustering
from ppanggolin.graph.makeGraph import compute_neighbors_graph
from ppanggolin.nem.rarefaction import make_rarefaction_curve
from ppanggolin.nem.partition import partition
from ppanggolin.formats.writeBinaries import write_pangenome
from ppanggolin.formats.writeFlatPangenome import write_pangenome_flat_files
from ppanggolin.formats.writeFlatGenomes import write_flat_genome_files
from ppanggolin.figures.ucurve import draw_ucurve
from ppanggolin.figures.tile_plot import draw_tile_plot
from ppanggolin.figures.draw_spot import draw_spots
from ppanggolin.info.info import print_info
from ppanggolin.RGP.genomicIsland import predict_rgp
from ppanggolin.RGP.spot import predict_hotspots
from ppanggolin.mod.module import predict_modules

"""a global workflow that does everything in one go."""


def launch_workflow(
    args: argparse.Namespace, panrgp: bool = True, panmodule: bool = True
):
    """
    Unified function to launch ppanggolin workflow.

    :param args: All arguments provide by command line
    :param panrgp: Flag to launch RGP prediction
    :param panmodule: Flag to launch module prediction
    """

    check_option_workflow(args)
    check_annotate_args(args)
    pangenome = Pangenome()

    filename = mk_file_name(args.basename, args.output, args.force)

    writing_time, anno_time, clust_time, mod_time, desc_time = (
        None,
        None,
        None,
        None,
        None,
    )

    if args.anno:  # if the annotations are provided, we read from it

        start_anno = time.time()
        read_annotations(
            pangenome,
            args.anno,
            pseudo=args.annotate.use_pseudo,
            cpu=args.annotate.cpu,
            translation_table=args.annotate.translation_table,
            disable_bar=args.disable_prog_bar,
        )
        anno_time = time.time() - start_anno

        if args.clusters is not None:
            start_clust = time.time()
            read_clustering(
                pangenome,
                args.clusters,
                infer_singleton=args.cluster.infer_singletons,
                code=args.cluster.translation_table,
                cpu=args.cluster.cpu,
                tmpdir=args.tmpdir,
                keep_tmp=args.cluster.keep_tmp,
                force=args.force,
                disable_bar=args.disable_prog_bar,
            )
        else:  # args.cluster is None
            if pangenome.status["geneSequences"] == "No":
                if args.fasta is None:
                    raise Exception(
                        "The gff/gbff provided did not have any sequence information, "
                        "you did not provide clusters and you did not provide fasta file. "
                        "Thus, we do not have the information we need to continue the analysis."
                    )
                else:
                    get_gene_sequences_from_fastas(pangenome, args.fasta)
            start_clust = time.time()
            clustering(
                pangenome,
                tmpdir=args.tmpdir,
                cpu=args.cluster.cpu,
                force=args.force,
                disable_bar=args.disable_prog_bar,
                defrag=not args.cluster.no_defrag,
                code=args.cluster.translation_table,
                coverage=args.cluster.coverage,
                identity=args.cluster.identity,
                mode=args.cluster.mode,
                keep_tmp_files=args.cluster.keep_tmp,
            )
        clust_time = time.time() - start_clust

        start_writing = time.time()
        write_pangenome(
            pangenome, filename, args.force, disable_bar=args.disable_prog_bar
        )
        writing_time = time.time() - start_writing

    elif args.fasta is not None:
        if args.clusters is not None:
            message = """
            You provide a list of fasta file and a clustering, which is incompatible. 
            Please provide a list of annotation files and if needed the list of fasta file, to use your clustering.
            Or you can let PPanGGOLiN manage the clustering step from your fasta file.
            """
            raise argparse.ArgumentError(argument=None, message=message)

        start_anno = time.time()
        annotate_pangenome(
            pangenome,
            args.fasta,
            tmpdir=args.tmpdir,
            cpu=args.annotate.cpu,
            disable_bar=args.disable_prog_bar,
            procedure=args.annotate.prodigal_procedure,
            translation_table=args.annotate.translation_table,
            kingdom=args.annotate.kingdom,
            norna=args.annotate.norna,
            allow_overlap=args.annotate.allow_overlap,
        )
        anno_time = time.time() - start_anno

        start_writing = time.time()
        write_pangenome(
            pangenome, filename, args.force, disable_bar=args.disable_prog_bar
        )
        writing_time = time.time() - start_writing

        start_clust = time.time()
        clustering(
            pangenome,
            tmpdir=args.tmpdir,
            cpu=args.cluster.cpu,
            force=args.force,
            disable_bar=args.disable_prog_bar,
            defrag=not args.cluster.no_defrag,
            code=args.cluster.translation_table,
            coverage=args.cluster.coverage,
            identity=args.cluster.identity,
            mode=args.cluster.mode,
            keep_tmp_files=args.cluster.keep_tmp,
        )
        clust_time = time.time() - start_clust

    write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)

    start_graph = time.time()
    compute_neighbors_graph(
        pangenome,
        args.graph.remove_high_copy_number,
        args.force,
        disable_bar=args.disable_prog_bar,
    )

    graph_time = time.time() - start_graph

    start_part = time.time()
    partition(
        pangenome,
        tmpdir=args.tmpdir,
        output=args.output,
        beta=args.partition.beta,
        sm_degree=args.partition.max_degree_smoothing,
        free_dispersion=args.partition.free_dispersion,
        chunk_size=args.partition.chunk_size,
        kval=args.partition.nb_of_partitions,
        krange=args.partition.krange,
        icl_margin=args.partition.ICL_margin,
        draw_icl=args.partition.draw_ICL,
        seed=args.partition.seed,
        keep_tmp_files=args.partition.keep_tmp_files,
        cpu=args.partition.cpu,
        force=args.force,
        disable_bar=args.disable_prog_bar,
    )
    part_time = time.time() - start_part

    start_writing = time.time()
    write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
    writing_time = writing_time + time.time() - start_writing

    regions_time, spot_time = (0, 0)
    if panrgp:
        start_regions = time.time()
        predict_rgp(
            pangenome,
            persistent_penalty=args.rgp.persistent_penalty,
            variable_gain=args.rgp.variable_gain,
            min_length=args.rgp.min_length,
            min_score=args.rgp.min_score,
            dup_margin=args.rgp.dup_margin,
            force=args.force,
            disable_bar=args.disable_prog_bar,
        )

        regions_time = time.time() - start_regions

        start_spots = time.time()
        predict_hotspots(
            pangenome,
            args.output,
            force=args.force,
            spot_graph=args.spot.spot_graph,
            overlapping_match=args.spot.overlapping_match,
            set_size=args.spot.set_size,
            exact_match=args.spot.exact_match_size,
            disable_bar=args.disable_prog_bar,
        )
        spot_time = time.time() - start_spots

    if panmodule:
        start_mods = time.time()
        predict_modules(
            pangenome=pangenome,
            dup_margin=args.module.dup_margin,
            size=args.module.size,
            min_presence=args.module.min_presence,
            transitive=args.module.transitive,
            jaccard=args.module.jaccard,
            force=args.force,
            disable_bar=args.disable_prog_bar,
        )

        mod_time = time.time() - start_mods

    start_writing = time.time()
    write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
    writing_time = writing_time + time.time() - start_writing

    if args.rarefaction_flag:
        make_rarefaction_curve(
            pangenome=pangenome,
            output=args.output,
            tmpdir=args.tmpdir,
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
            cpu=args.rarefaction.cpu,
            disable_bar=args.disable_prog_bar,
        )

    if not args.no_flat_files:

        if panrgp and args.draw.draw_spots:
            start_spot_drawing = time.time()
            mk_outdir(args.output / "spot_figures", force=True)
            draw_spots(
                pangenome=pangenome,
                output=args.output / "spot_figures",
                spot_list=args.draw.spots,
                disable_bar=args.disable_prog_bar,
            )
            spot_time += time.time() - start_spot_drawing

        if args.draw.tile_plot:
            if (
                pangenome.number_of_organisms < 65000
                or pangenome.number_of_gene_families < 65000
            ):
                nocloud = (
                    args.draw.nocloud
                    if pangenome.number_of_organisms < 32767
                    or pangenome.number_of_gene_families < 32767
                    else True
                )
                draw_tile_plot(
                    pangenome,
                    args.output,
                    nocloud=nocloud,
                    disable_bar=args.disable_prog_bar,
                    draw_dendrogram=args.draw.add_dendrogram,
                    add_metadata=True,
                )
            else:
                logging.getLogger("PPanGGOLiN").warning(
                    "Tile plot output have been requested but there are too many genomes or families to produce a viewable tile plot."
                )

        if args.draw.ucurve:
            draw_ucurve(
                pangenome,
                args.output,
                disable_bar=args.disable_prog_bar,
                soft_core=args.draw.soft_core,
            )

        start_desc = time.time()

        write_pangenome_arguments = [
            "gexf",
            "light_gexf",
            "json",
            "csv",
            "Rtab",
            "stats",
            "partitions",
            "families_tsv",
        ]

        # Check that we don't ask write to output something not computed.
        borders, spots, spot_modules, modules, regions = (
            False,
            False,
            False,
            False,
            False,
        )

        if panmodule:
            modules = args.write_pangenome.modules
            write_pangenome_arguments.append("modules")

        if panrgp:
            borders, spots, regions = (
                args.write_pangenome.borders,
                args.write_pangenome.spots,
                args.write_pangenome.regions,
            )
            write_pangenome_arguments += ["borders", "spots", "regions"]

        if panmodule and panrgp:
            spot_modules = args.write_pangenome.spot_modules
            write_pangenome_arguments.append("spot_modules")

        # check that at least one output file is requested. if not write is not call.
        if any(
            getattr(args.write_pangenome, arg) is True
            for arg in write_pangenome_arguments
        ):
            # some parameters are set to false because they have not been computed in this workflow
            write_pangenome_flat_files(
                pangenome,
                args.output,
                cpu=args.write_pangenome.cpu,
                disable_bar=args.disable_prog_bar,
                soft_core=args.write_pangenome.soft_core,
                dup_margin=args.write_pangenome.dup_margin,
                csv=args.write_pangenome.csv,
                gene_pa=args.write_pangenome.Rtab,
                gexf=args.write_pangenome.gexf,
                light_gexf=args.write_pangenome.light_gexf,
                stats=args.write_pangenome.stats,
                json=args.write_pangenome.json,
                partitions=args.write_pangenome.partitions,
                families_tsv=args.write_pangenome.families_tsv,
                regions=regions,
                compress=args.write_pangenome.compress,
                spot_modules=spot_modules,
                modules=modules,
                spots=spots,
                borders=borders,
            )

        else:
            logging.getLogger("PPanGGOLiN").info(
                "No flat file describing the pangenome has been requested in config file. "
                "Writing output pangenome flat file is skipped."
            )

        write_genomes_arguments = ["proksee", "table", "gff"]
        if any(
            getattr(args.write_genomes, arg) is True for arg in write_genomes_arguments
        ):
            write_flat_genome_files(
                pangenome,
                args.output,
                proksee=args.write_genomes.proksee,
                table=args.write_genomes.table,
                gff=args.write_genomes.gff,
                add_metadata=True,
                compress=args.write_genomes.compress,
                disable_bar=args.disable_prog_bar,
                cpu=args.write_genomes.cpu,
            )
        else:
            logging.getLogger("PPanGGOLiN").info(
                "No flat file of genomes with pangenome annotation has been requested in config file. "
                "Writing output genomes flat file is skipped."
            )

        desc_time = time.time() - start_desc

    logging.getLogger("PPanGGOLiN").info(
        f"Annotation took : {round(anno_time, 2)} seconds"
    )
    logging.getLogger("PPanGGOLiN").info(
        f"Clustering took : {round(clust_time, 2)} seconds"
    )
    logging.getLogger("PPanGGOLiN").info(
        f"Building the graph took : {round(graph_time, 2)} seconds"
    )
    logging.getLogger("PPanGGOLiN").info(
        f"Partitioning the pangenome took : {round(part_time, 2)} seconds"
    )

    if panrgp:
        logging.getLogger("PPanGGOLiN").info(
            f"Predicting RGP took : {round(regions_time, 2)} seconds"
        )
        logging.getLogger("PPanGGOLiN").info(
            f"Gathering RGP into spots took : {round(spot_time, 2)} seconds"
        )

    if panmodule:
        logging.getLogger("PPanGGOLiN").info(
            f"Predicting modules took : {round(mod_time, 2)} seconds"
        )

    logging.getLogger("PPanGGOLiN").info(
        f"Writing the pangenome data in HDF5 took : {round(writing_time, 2)} seconds"
    )

    if not args.no_flat_files:
        logging.getLogger("PPanGGOLiN").info(
            f"Writing descriptive files for the pangenome took : {round(desc_time, 2)} seconds"
        )

    print_info(filename, content=True)


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """

    launch_workflow(args, panrgp=True, panmodule=True)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line
    :param sub_parser : sub_parser for all command
    :return : parser arguments for all command
    """
    parser = sub_parser.add_parser("all", formatter_class=argparse.RawTextHelpFormatter)

    add_workflow_args(parser)

    return parser


def add_workflow_args(parser: argparse.ArgumentParser):
    """
    Parser for important arguments that can be changed in CLI.
    Other (less important) arguments that are step specific can be changed in the config file.
    :param parser: parser for workflow argument
    """
    date = time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())
    required = parser.add_argument_group(
        title="Input arguments", description="The possible input arguments :"
    )

    required.add_argument(
        "--fasta",
        required=False,
        type=Path,
        help="A tab-separated file listing the genome names, "
        "and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). "
        "One line per genome. This option can be used alone.",
    )

    required.add_argument(
        "--anno",
        required=False,
        type=Path,
        help="A tab-separated file listing the genome names, and the gff filepath of "
        "its annotations (the gffs can be compressed). One line per genome. "
        "This option can be used alone IF the fasta sequences are in the gff files, "
        "otherwise --fasta needs to be used.",
    )

    required.add_argument(
        "--clusters",
        required=False,
        type=Path,
        help="a tab-separated file listing the cluster names, the gene IDs, "
        "and optionally whether they are a fragment or not.",
    )

    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument(
        "-o",
        "--output",
        required=False,
        type=Path,
        default=Path(f"ppanggolin_output{date}_PID{str(os.getpid())}"),
        help="Output directory",
    )

    optional.add_argument(
        "--basename",
        required=False,
        default="pangenome",
        help="basename for the output file",
    )

    optional.add_argument(
        "--rarefaction",
        required=False,
        action="store_true",
        dest="rarefaction_flag",
        help="Use to compute the rarefaction curves (WARNING: can be time consuming)",
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
        "--translation_table",
        required=False,
        type=int,
        default=11,
        help="Translation table (genetic code) to use.",
    )

    optional.add_argument(
        "--kingdom",
        required=False,
        type=str.lower,
        default="bacteria",
        choices=["bacteria", "archaea"],
        help="Kingdom to which the prokaryota belongs to, "
        "to know which models to use for rRNA annotation.",
    )

    optional.add_argument(
        "--mode",
        required=False,
        default="1",
        choices=["0", "1", "2", "3"],
        help="the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component),"
        " 2: CD-HIT-like, 3: CD-HIT-like (lowmem)",
    )

    optional.add_argument(
        "--coverage",
        required=False,
        type=restricted_float,
        default=0.8,
        help="Minimal coverage of the alignment for two proteins to be in the same cluster",
    )

    optional.add_argument(
        "--identity",
        required=False,
        type=restricted_float,
        default=0.8,
        help="Minimal identity percent for two proteins to be in the same cluster",
    )

    optional.add_argument(
        "--infer_singletons",
        required=False,
        action="store_true",
        help="Use this option together with --clusters. "
        "If a gene is not present in the provided clustering result file, "
        "it will be assigned to its own unique cluster as a singleton.",
    )

    optional.add_argument(
        "--use_pseudo",
        required=False,
        action="store_true",
        help="In the context of provided annotation, use this option to read pseudogenes. "
        "(Default behavior is to ignore them)",
    )

    optional.add_argument(
        "-K",
        "--nb_of_partitions",
        required=False,
        default=-1,
        type=int,
        help="Number of partitions to use. Must be at least 2. If under 2, "
        "it will be detected automatically.",
    )

    # This ensures compatibility with workflows built with the old option "defrag" when it was not the default
    optional.add_argument(
        "--no_defrag",
        required=False,
        action="store_true",
        help="DO NOT Realign gene families to link fragments with their non-fragmented gene family.",
    )

    optional.add_argument(
        "--no_flat_files",
        required=False,
        action="store_true",
        help="Generate only the HDF5 pangenome file.",
    )
    optional.add_argument(
        "--tmpdir",
        required=False,
        type=str,
        default=Path(tempfile.gettempdir()),
        help="directory for storing temporary files",
    )
