import ppanggolin.nem.rarefaction
import ppanggolin.nem.partition
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
import ppanggolin.projection
# import ppanggolin.utility
import ppanggolin.meta

SUBCOMMAND_TO_SUBPARSER = {
    "annotate": ppanggolin.annotate.subparser,
    "cluster": ppanggolin.cluster.subparser,
    "graph": ppanggolin.graph.subparser,
    "partition": ppanggolin.nem.partition.subparser,
    "rarefaction": ppanggolin.nem.rarefaction.subparser,
    "workflow": ppanggolin.workflow.workflow.subparser,
    "panrgp": ppanggolin.workflow.panRGP.subparser,
    "panModule": ppanggolin.workflow.panModule.subparser,
    "all": ppanggolin.workflow.all.subparser,
    "draw": ppanggolin.figures.subparser,
    "write_pangenome": ppanggolin.formats.writeFlatPangenome.subparser,
    "write_genomes": ppanggolin.formats.writeFlatGenomes.subparser,
    "fasta": ppanggolin.formats.writeSequences.subparser,
    "msa": ppanggolin.formats.writeMSA.subparser,
    "metrics": ppanggolin.metrics.metrics.subparser,
    "align": ppanggolin.align.subparser,
    "rgp": ppanggolin.RGP.genomicIsland.subparser,
    "spot": ppanggolin.RGP.spot.subparser,
    "module": ppanggolin.mod.subparser,
    "context": ppanggolin.context.subparser,
    "projection":ppanggolin.projection.subparser,
    "rgp_cluster":ppanggolin.RGP.rgp_cluster.subparser,
    "metadata": ppanggolin.meta.subparser
}
