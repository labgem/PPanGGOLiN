from importlib.metadata import distribution
import ppanggolin.nem.rarefaction
import ppanggolin.nem.partition
import ppanggolin.graph
import ppanggolin.annotate
import ppanggolin.cluster
import ppanggolin.figures
import ppanggolin.formats
import ppanggolin.metrics
import ppanggolin.align
import ppanggolin.RGP
import ppanggolin.mod
import ppanggolin.context
import ppanggolin.workflow
import ppanggolin.projection
import ppanggolin.meta

SUBCOMMAND_TO_SUBPARSER = {
    "annotate": ppanggolin.annotate.subparser,
    "cluster": ppanggolin.cluster.subparser,
    "graph": ppanggolin.graph.subparser,
    "partition": ppanggolin.nem.partition.subparser,
    "rarefaction": ppanggolin.nem.rarefaction.subparser,
    "workflow": ppanggolin.workflow.workflow.subparser,
    "panrgp": ppanggolin.workflow.panRGP.subparser,
    "panmodule": ppanggolin.workflow.panModule.subparser,
    "all": ppanggolin.workflow.all.subparser,
    "draw": ppanggolin.figures.subparser,
    "write_pangenome": ppanggolin.formats.writeFlatPangenome.subparser,
    "write_genomes": ppanggolin.formats.writeFlatGenomes.subparser,
    "write_metadata": ppanggolin.formats.writeFlatMetadata.subparser,
    "fasta": ppanggolin.formats.writeSequences.subparser,
    "msa": ppanggolin.formats.writeMSA.subparser,
    "metrics": ppanggolin.metrics.metrics.subparser,
    "align": ppanggolin.align.subparser,
    "rgp": ppanggolin.RGP.genomicIsland.subparser,
    "spot": ppanggolin.RGP.spot.subparser,
    "module": ppanggolin.mod.subparser,
    "context": ppanggolin.context.subparser,
    "projection": ppanggolin.projection.subparser,
    "rgp_cluster": ppanggolin.RGP.rgp_cluster.subparser,
    "metadata": ppanggolin.meta.subparser,
}


version = distribution("ppanggolin").version

epilog = f"""
PPanGGOLiN ({version}) is an open-source bioinformatics tool developed by the LABGeM team, and distributed under the CeCILL Free Software License Agreement.
"""

pan_epilog = """
For pangenome analyses, please cite:
Gautreau G et al. (2020) PPanGGOLiN: Depicting microbial diversity via a partitioned pangenome graph.
PLOS Computational Biology 16(3): e1007732. https://doi.org/10.1371/journal.pcbi.1007732
"""
rgp_epilog = """
For genomic islands and spots of insertion detection, please cite:
Bazin et al., panRGP: a pangenome-based method to predict genomic islands and explore their diversity,
Bioinformatics, Volume 36, Issue Supplement_2, December 2020, Pages i651–i658, https://doi.org/10.1093/bioinformatics/btaa792
"""

mod_epilog = """
For module prediction, please cite:
Bazin et al., panModule: detecting conserved modules in the variable regions of a pangenome graph.
biorxiv. https://doi.org/10.1101/2021.12.06.471380
"""
