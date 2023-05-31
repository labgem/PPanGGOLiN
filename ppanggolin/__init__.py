import ppanggolin.nem.rarefaction
import ppanggolin.graph
import ppanggolin.annotate
import ppanggolin.cluster
import ppanggolin.figures
import ppanggolin.formats
# import ppanggolin.info
import ppanggolin.metrics
import ppanggolin.align
import ppanggolin.RGP
import ppanggolin.mod
import ppanggolin.context
import ppanggolin.workflow
# import ppanggolin.utility


SUBCOMMAND_TO_SUBPARSER = {
        "annotate":ppanggolin.annotate.subparser,
        "cluster":ppanggolin.cluster.subparser,
        "graph":ppanggolin.graph.subparser,
        "partition":ppanggolin.nem.partition.subparser,
        "rarefaction":ppanggolin.nem.rarefaction.subparser,
        "workflow":ppanggolin.workflow.workflow.subparser,
        "panrgp":ppanggolin.workflow.panRGP.subparser,
        "panModule":ppanggolin.workflow.panModule.subparser,
        "all":ppanggolin.workflow.all.subparser,
        "draw":ppanggolin.figures.subparser,
        "write":ppanggolin.formats.writeFlat.subparser,
        "fasta":ppanggolin.formats.writeSequences.subparser,
        "msa":ppanggolin.formats.writeMSA.subparser,
        "metrics":ppanggolin.metrics.metrics.subparser,
        "align":ppanggolin.align.subparser,
        "rgp":ppanggolin.RGP.genomicIsland.subparser,
        "spot":ppanggolin.RGP.spot.subparser,
        "rgp_cluster":ppanggolin.RGP.rgp_cluster.subparser,
        "module":ppanggolin.mod.subparser,
        "context":ppanggolin.context.subparser,# "info":ppanggolin.info.subparser, "default_config":ppanggolin.utility.default_config.subparser
        }
