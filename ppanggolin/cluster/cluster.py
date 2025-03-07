#!/usr/bin/env python3

# default libraries
import logging
import tempfile
from collections import defaultdict
import os
import argparse
from typing import Tuple, Dict, Set
from pathlib import Path
import time
import gzip

# installed libraries
from networkx import Graph
from tqdm import tqdm
import pandas as pd

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Gene
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.utils import (
    is_compressed,
    restricted_float,
    run_subprocess,
    create_tmpdir,
    check_tools_availability,
)
from ppanggolin.formats.writeBinaries import write_pangenome, erase_pangenome
from ppanggolin.formats.readBinaries import (
    check_pangenome_info,
    write_gene_sequences_from_pangenome_file,
)
from ppanggolin.formats.writeSequences import (
    write_gene_sequences_from_annotations,
    translate_genes,
    create_mmseqs_db,
)


# Global functions
def check_pangenome_former_clustering(pangenome: Pangenome, force: bool = False):
    """
    Checks pangenome status and .h5 files for former clusterings, delete them if allowed or raise an error

    :param pangenome: Annotated Pangenome
    :param force: Force to write on existing pangenome information
    """
    if pangenome.status["genesClustered"] == "inFile" and not force:
        raise Exception(
            "You are trying to cluster genes that are already clustered together. If you REALLY want to "
            "do that, use --force (it will erase everything except annotation data in your HDF5 file!)"
        )
    elif pangenome.status["genesClustered"] == "inFile" and force:
        erase_pangenome(pangenome, gene_families=True)


# Clustering functions
def check_pangenome_for_clustering(
    pangenome: Pangenome,
    sequences: Path,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Check the pangenome statuses and write the gene sequences in the provided tmpFile.
    (whether they are written in the .h5 file or currently in memory)

    :param pangenome: Annotated Pangenome
    :param sequences: Path to write the sequences
    :param force: Force to write on existing pangenome information
    :param disable_bar: Allow to disable progress bar
    """
    check_pangenome_former_clustering(pangenome, force)
    if pangenome.status["geneSequences"] in ["Computed", "Loaded"]:
        logging.getLogger("PPanGGOLiN").debug(
            "Write sequences from annotation loaded in pangenome"
        )
        # we append the gene ids by 'ppanggolin' to avoid crashes from mmseqs when sequence IDs are only numeric.
        write_gene_sequences_from_annotations(
            pangenome.genes,
            sequences,
            add="ppanggolin_",
            compress=False,
            disable_bar=disable_bar,
        )
    elif pangenome.status["geneSequences"] == "inFile":
        logging.getLogger("PPanGGOLiN").debug("Write sequences from pangenome file")
        write_gene_sequences_from_pangenome_file(
            pangenome.file,
            sequences,
            add="ppanggolin_",
            compress=False,
            disable_bar=disable_bar,
        )  # write CDS sequences to the tmpFile
    else:
        raise Exception(
            "The pangenome does not include gene sequences, thus it is impossible to cluster "
            "the genes in gene families. Either provide clustering results (see --clusters), "
            "or provide a way to access the gene sequence during the annotation step "
            "(having the fasta in the gff files, or providing the fasta files through the --fasta option)"
        )


def first_clustering(
    sequences: Path,
    tmpdir: Path,
    cpu: int = 1,
    code: int = 11,
    coverage: float = 0.8,
    identity: float = 0.8,
    mode: int = 1,
) -> Tuple[Path, Path]:
    """
    Make a first clustering of all sequences in pangenome

    :param sequences: Sequence from pangenome
    :param tmpdir: Temporary directory
    :param cpu: number of CPU cores to use
    :param code: Genetic code used
    :param coverage: minimal coverage threshold for the alignment
    :param identity: minimal identity threshold for the alignment
    :param mode: MMseqs2 clustering mode

    :return: path to representative sequence file and path to tsv clustering result
    """

    seqdb = translate_genes(
        sequences=sequences,
        tmpdir=tmpdir,
        cpu=cpu,
        is_single_line_fasta=True,
        code=code,
    )
    logging.getLogger("PPanGGOLiN").info("Clustering sequences...")
    cludb = tmpdir / "cluster_db"
    cmd = list(
        map(
            str,
            [
                "mmseqs",
                "cluster",
                seqdb,
                cludb,
                tmpdir,
                "--cluster-mode",
                mode,
                "--min-seq-id",
                identity,
                "-c",
                coverage,
                "--threads",
                cpu,
                "--kmer-per-seq",
                80,
                "--max-seqs",
                300,
            ],
        )
    )
    run_subprocess(cmd, msg="MMSeqs2 cluster failed with the following error:\n")
    logging.getLogger("PPanGGOLiN").info("Extracting cluster representatives...")
    repdb = tmpdir / "representative_db"
    cmd = list(
        map(str, ["mmseqs", "result2repseq", seqdb, cludb, repdb, "--threads", cpu])
    )
    run_subprocess(cmd, msg="MMSeqs2 result2repseq failed with the following error:\n")
    reprfa = tmpdir / "representative_sequences.fasta"
    cmd = list(
        map(
            str,
            [
                "mmseqs",
                "result2flat",
                seqdb,
                seqdb,
                repdb,
                reprfa,
                "--use-fasta-header",
            ],
        )
    )
    run_subprocess(cmd, msg="MMSeqs2 result2flat failed with the following error:\n")
    logging.getLogger("PPanGGOLiN").info("Writing gene to family information")
    outtsv = tmpdir / "families_tsv"
    cmd = list(
        map(
            str,
            [
                "mmseqs",
                "createtsv",
                seqdb,
                seqdb,
                cludb,
                outtsv,
                "--threads",
                cpu,
                "--full-header",
            ],
        )
    )
    run_subprocess(cmd, msg="MMSeqs2 createtsv failed with the following error:\n")
    return reprfa, outtsv


def read_faa(faa_file_name: Path) -> Dict[str, str]:
    """
    Read a faa file to link pangenome families to sequences.

    :param faa_file_name: path to the faa file

    :return: dictionary with families ID as key and sequence as value
    """
    fam2seq = {}
    head = ""
    with open(faa_file_name) as faaFile:
        for line in faaFile:
            if line.startswith(">"):
                head = (
                    line[1:].strip().replace("ppanggolin_", "")
                )  # remove the eventual addition
            else:
                fam2seq[head] = line.strip()
    return fam2seq


def align_rep(
    faa_file: Path,
    tmpdir: Path,
    cpu: int = 1,
    coverage: float = 0.8,
    identity: float = 0.8,
) -> Path:
    """
    Align representative sequence

    :param faa_file: sequence of representative family
    :param tmpdir: Temporary directory
    :param cpu: number of CPU cores to use
    :param coverage: minimal coverage threshold for the alignment
    :param identity: minimal identity threshold for the alignment

    :return: Result of alignment
    """
    seqdb = create_mmseqs_db(
        [faa_file], "rep_sequence_db", tmpdir, db_mode=1, db_type=1
    )
    logging.getLogger("PPanGGOLiN").info("Aligning cluster representatives...")
    alndb = tmpdir / "rep_alignment_db"
    cmd = list(
        map(
            str,
            [
                "mmseqs",
                "search",
                seqdb,
                seqdb,
                alndb,
                tmpdir,
                "-a",
                "--min-seq-id",
                identity,
                "-c",
                coverage,
                "--cov-mode",
                1,
                "--threads",
                cpu,
            ],
        )
    )
    run_subprocess(cmd, msg="MMSeqs2 search failed with the following error:\n")
    logging.getLogger("PPanGGOLiN").info("Extracting alignments...")
    outfile = tmpdir / "rep_families.tsv"
    cmd = list(
        map(
            str,
            [
                "mmseqs",
                "convertalis",
                seqdb,
                seqdb,
                alndb,
                outfile,
                "--format-output",
                "query,target,qlen,tlen,bits",
            ],
        )
    )
    run_subprocess(cmd, msg="MMSeqs2 convertalis failed with the following error:\n")
    return outfile


def read_tsv(
    tsv_file_name: Path,
) -> Tuple[Dict[str, Tuple[str, bool]], Dict[str, Set[str]]]:
    """Reading tsv file

    :param tsv_file_name: path to the tsv

    :return: two dictionaries which link genes and families
    """
    genes2fam = {}
    fam2genes = defaultdict(set)
    with open(tsv_file_name) as tsvfile:
        for line in tsvfile:
            line = line.replace('"', "").replace("ppanggolin_", "").split()
            # remove the '"' char which protects the fields, and the eventual addition
            genes2fam[line[1]] = (
                line[0],
                False,
            )  # fam id, and it's a gene (and not a fragment)
            fam2genes[line[0]].add(line[1])
    return genes2fam, fam2genes


def refine_clustering(
    tsv: Path, aln_file: Path, fam_to_seq: dict
) -> Tuple[Dict[str, Tuple[str, bool]], Dict[str, str]]:
    """
    Refine clustering by removing fragment

    :param tsv: First clusterin result
    :param aln_file: Reprensentative alignment result
    :param fam_to_seq: Dictionary which link families to sequence

    :return: Two dictionary which link genes and families
    """
    simgraph = Graph()
    genes2fam, fam2genes = read_tsv(tsv)
    logging.getLogger("PPanGGOLiN").info(f"Starting with {len(fam_to_seq)} families")
    # create the nodes
    for fam, genes in fam2genes.items():
        simgraph.add_node(fam, nbgenes=len(genes))

    # add the edges
    with open(aln_file) as alnfile:
        for line in alnfile:
            line = (
                line.replace('"', "").replace("ppanggolin_", "").split()
            )  # remove the eventual addition

            if line[0] != line[1]:
                simgraph.add_edge(line[0], line[1], score=float(line[4]))
                simgraph.nodes[line[0]]["length"] = int(line[2])
                simgraph.nodes[line[1]]["length"] = int(line[3])

    for node, nodedata in sorted(simgraph.nodes(data=True)):
        choice = (None, 0, 0, 0)
        for neighbor in sorted(simgraph.neighbors(node)):
            nei = simgraph.nodes[neighbor]
            score = simgraph[neighbor][node]["score"]
            if (
                nei["length"] > nodedata["length"]
                and nei["nbgenes"] >= nodedata["nbgenes"]
                and choice[3] < score
            ):
                choice = (genes2fam[neighbor][0], nei["length"], nei["nbgenes"], score)
                # `genes2fam[neighbor]` instead of just neighbor in case that family has been assigned already
                # (this is for smaller fragments that are closer to other fragments than the actual gene family)

        if choice[0] is not None:
            genestochange = fam2genes[node]
            for gene in genestochange:
                genes2fam[gene] = (str(choice[0]), True)
                fam2genes[choice[0]].add(gene)
            del fam2genes[node]
    new_fam_to_seq = {}
    for fam in fam2genes:
        new_fam_to_seq[fam] = fam_to_seq[fam]
    logging.getLogger("PPanGGOLiN").info(
        f"Ending with {len(new_fam_to_seq)} gene families"
    )
    return genes2fam, new_fam_to_seq


def read_fam2seq(pangenome: Pangenome, fam_to_seq: Dict[str, str]):
    """
    Add gene family to pangenome and sequences to gene families

    :param pangenome: Annotated pangenome
    :param fam_to_seq: Dictionary which link families and sequences
    """
    logging.getLogger("PPanGGOLiN").info(
        "Adding protein sequences to the gene families"
    )
    for family, protein in fam_to_seq.items():
        fam = GeneFamily(pangenome.max_fam_id, family)
        fam.add_sequence(protein)
        pangenome.add_gene_family(fam)


def read_gene2fam(pangenome: Pangenome, gene_to_fam: dict, disable_bar: bool = False):
    """
    Add gene to pangenome families

    :param pangenome: Annotated Pangenome
    :param gene_to_fam: Dictionary which link gene to families
    :param disable_bar: Allow to disable progress bar
    """
    logging.getLogger("PPanGGOLiN").info(
        f"Adding {len(gene_to_fam)} genes to the gene families"
    )

    link = (
        True
        if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"]
        else False
    )
    if (
        link and len(gene_to_fam) != pangenome.number_of_genes
    ):  # then maybe there are genes with identical IDs
        logging.getLogger("PPanGGOLiN").debug(
            f"gene_to_fam size: {len(gene_to_fam)}, "
            f"Pangenome nb genes: {pangenome.number_of_genes}"
        )
        raise Exception(
            "Something unexpected happened during clustering (have less genes clustered than genes "
            "in the pangenome). A probable reason is that two genes in two different genomes have "
            "the same IDs; If you are sure that all of your genes have non identical IDs,  please post an "
            "issue at https://github.com/labgem/PPanGGOLiN/"
        )
    for gene, (family, is_frag) in tqdm(
        gene_to_fam.items(), unit="gene", total=len(gene_to_fam), disable=disable_bar
    ):
        try:
            fam = pangenome.get_gene_family(family)
        except KeyError:  # Family not found so create and add
            fam = GeneFamily(pangenome.max_fam_id, family)
            pangenome.add_gene_family(fam)
        if link:  # doing the linking if the annotations are loaded.
            gene_obj = pangenome.get_gene(gene)
        else:
            gene_obj = Gene(gene)
        gene_obj.is_fragment = is_frag
        fam.add(gene_obj)


def clustering(
    pangenome: Pangenome,
    tmpdir: Path,
    cpu: int = 1,
    defrag: bool = True,
    code: int = 11,
    coverage: float = 0.8,
    identity: float = 0.8,
    mode: int = 1,
    force: bool = False,
    disable_bar: bool = False,
    keep_tmp_files: bool = True,
):
    """
    Cluster gene sequences from an annotated pangenome into families.

    :param pangenome: Annotated Pangenome object.
    :param tmpdir: Path to a temporary directory for intermediate files.
    :param cpu: Number of CPU cores to use for clustering.
    :param defrag: Allow removal of fragmented sequences during clustering.
    :param code: Genetic code used for sequence translation.
    :param coverage: Minimum coverage threshold for sequence alignment during clustering.
    :param identity: Minimum identity threshold for sequence alignment during clustering.
    :param mode: Clustering mode (MMseqs2 mode).
    :param force: Force writing clustering results back to the pangenome.
    :param disable_bar: Disable the progress bar during clustering.
    :param keep_tmp_files: Keep temporary files (useful for debugging).
    """

    check_tools_availability(["mmseqs"])

    date = time.strftime("_%Y-%m-%d_%H-%M-%S", time.localtime())
    dir_name = f"clustering_tmpdir_{date}_PID{os.getpid()}"
    with create_tmpdir(tmpdir, basename=dir_name, keep_tmp=keep_tmp_files) as tmp_path:
        sequence_path = tmp_path / "nucleotide_sequences.fna"
        check_pangenome_for_clustering(
            pangenome, sequence_path, force, disable_bar=disable_bar
        )
        logging.getLogger("PPanGGOLiN").info("Clustering all of the genes sequences...")
        rep, tsv = first_clustering(
            sequence_path, tmp_path, cpu, code, coverage, identity, mode
        )

        fam2seq = read_faa(rep)
        if not defrag:
            logging.getLogger("PPanGGOLiN").debug("No defragmentation")
            genes2fam, _ = read_tsv(tsv)
        else:
            logging.getLogger("PPanGGOLiN").info(
                "Associating fragments to their original gene family..."
            )
            aln = align_rep(rep, tmp_path, cpu, coverage, identity)
            genes2fam, fam2seq = refine_clustering(tsv, aln, fam2seq)
            pangenome.status["defragmented"] = "Computed"
    read_fam2seq(pangenome, fam2seq)
    read_gene2fam(pangenome, genes2fam, disable_bar=disable_bar)

    pangenome.status["genesClustered"] = "Computed"
    pangenome.status["geneFamilySequences"] = "Computed"

    pangenome.parameters["cluster"] = {}
    pangenome.parameters["cluster"]["coverage"] = coverage
    pangenome.parameters["cluster"]["identity"] = identity
    pangenome.parameters["cluster"]["mode"] = mode
    pangenome.parameters["cluster"]["# defragmentation"] = defrag
    pangenome.parameters["cluster"]["no_defrag"] = not defrag

    pangenome.parameters["cluster"]["translation_table"] = code
    pangenome.parameters["cluster"]["# read_clustering_from_file"] = False


# Read clustering
def mk_local_to_gene(pangenome: Pangenome) -> dict:
    """Creates a dictionary that stores local identifiers, if all local identifiers are unique (and if they exist)

    :param pangenome: Input Pangenome

    :return: Dictionary with local identifiers
    """
    local_dict = {}
    for gene in pangenome.genes:
        old_len = len(local_dict)
        local_dict[gene.local_identifier] = gene
        if len(local_dict) == old_len:
            if (
                pangenome.parameters["annotate"]["# read_annotations_from_file"]
                and not pangenome.parameters["annotate"]["# used_local_identifiers"]
            ):
                raise Exception(
                    f"'{gene.local_identifier}' was found multiple times used as an identifier. "
                    f"The identifier of the genes (locus_tag, protein_id in gbff, ID in gff) were not "
                    f"unique throughout all of the files. It is thus impossible to differentiate the genes."
                    f" To use this function while importing annotate, all identifiers MUST be unique "
                    f"throughout all of your genomes"
                )
            return {}  # local identifiers are not unique.
    return local_dict


def infer_singletons(pangenome: Pangenome):
    """
    Creates a new family for each gene with no associated family.

    :param pangenome: Input pangenome object
    """
    singleton_counter = 0
    for gene in pangenome.genes:
        if gene.family is None:
            # Create a new family for the singleton gene
            fam = GeneFamily(family_id=pangenome.max_fam_id, name=gene.ID)
            fam.representative = gene
            fam.add(gene)

            # Try to add the new family
            try:
                pangenome.add_gene_family(fam)
            except KeyError:
                raise KeyError(
                    f"Cannot create singleton family with name='{fam.name}' for gene '{gene.ID}': "
                    f"A family with the same name already exists. Check the gene '{gene.ID}' in input cluster file."
                )

            singleton_counter += 1

    logging.getLogger("PPanGGOLiN").info(
        f"Inferred {singleton_counter} singleton families"
    )


def get_family_representative_sequences(
    pangenome: Pangenome,
    code: int = 11,
    cpu: int = 1,
    tmpdir: Path = None,
    keep_tmp: bool = False,
):

    logging.getLogger("PPanGGOLiN").info(
        "Retrieving protein sequences of family representatives."
    )

    tmpdir = Path(tempfile.gettempdir()) if tmpdir is None else tmpdir
    with create_tmpdir(tmpdir, "get_proteins_sequences", keep_tmp) as tmp:

        repres_path = tmp / "representative.fna.gz"

        with gzip.open(repres_path, mode="wt") as repres_seq:

            for family in pangenome.gene_families:

                if family.representative.dna is None:
                    raise ValueError(
                        f"DNA sequence of representative gene {family.representative} is None. "
                        "Sequence may not have been loaded correctly from the pangenome file or the pangenome has no gene sequences."
                    )

                repres_seq.write(f">{family.name}\n")
                repres_seq.write(f"{family.representative.dna}\n")

        translate_db = translate_genes(
            sequences=repres_path,
            tmpdir=tmp,
            cpu=cpu,
            is_single_line_fasta=True,
            code=code,
        )

        outpath = tmp / "representative_protein_genes.fna"
        cmd = list(map(str, ["mmseqs", "convert2fasta", translate_db, outpath]))
        run_subprocess(
            cmd, msg="MMSeqs convert2fasta failed with the following error:\n"
        )

        with open(outpath) as repres_prot:
            lines = repres_prot.readlines()
            while len(lines) > 0:
                family_name = lines.pop(0).strip()[1:]
                family_seq = lines.pop(0).strip()
                family = pangenome.get_gene_family(family_name)

                family.add_sequence(family_seq)


def read_clustering_file(families_tsv_path: Path) -> Tuple[pd.DataFrame, bool]:
    """
    Read and process a gene families clustering file.

    This function reads a tab-separated gene families file and processes it into a DataFrame with
    appropriate columns. It handles different formats of the input file based on the number of columns.

    The function expects the file to have 2, 3, or 4 columns:
    - 2 columns: ["family", "gene"]
    - 3 columns: ["family", "gene", "is_frag"] or ["family", "gene", "representative"]
    - 4 columns: ["family", "representative", "gene", "is_frag"]

    :param families_tsv_path: The path to the gene families file, which can be compressed or uncompressed.

    :raises ValueError: If the file has only one column or an unexpected number of columns.
    :raises Exception: If there are duplicated gene IDs in the clustering.

    :return: The processed DataFrame and a boolean indicating if any gene is marked as fragmented.
    """
    logging.getLogger("PPanGGOLiN").info(
        f"Reading clustering file to group genes into families: {families_tsv_path.as_posix()}"
    )

    # Detect compression type if any
    _, compress_type = is_compressed(families_tsv_path)

    # Read the file with inferred compression if necessary
    families_df = pd.read_csv(
        families_tsv_path,
        sep="\t",
        header=None,
        compression=compress_type if compress_type is not None else "infer",
        dtype=str,
    )

    # Process DataFrame based on the number of columns
    if families_df.shape[1] == 2:
        families_df.columns = ["family", "gene"]
        families_df["representative"] = families_df.groupby("family")["gene"].transform(
            "first"
        )
        families_df["is_frag"] = False

    elif families_df.shape[1] == 3:
        # Check if the third column is 'is_frag'
        if families_df[2].dropna().eq("F").all():
            families_df.columns = ["family", "gene", "is_frag"]
            families_df["is_frag"] = (
                families_df["is_frag"].replace("F", True).fillna(False)
            )
            families_df["representative"] = families_df.groupby("family")[
                "gene"
            ].transform("first")
        else:
            families_df.columns = ["family", "gene", "representative"]
            families_df["is_frag"] = False

    elif families_df.shape[1] == 4:
        families_df.columns = ["family", "representative", "gene", "is_frag"]

    else:
        raise ValueError(
            f"Unexpected number of columns ({families_df.shape[1]}). The file must have 2, 3, or 4 columns."
        )

    # Ensure columns are strings
    families_df["family"] = families_df["family"].astype(str)
    families_df["gene"] = families_df["gene"].astype(str)
    families_df["representative"] = families_df["representative"].astype(str)

    # Check for duplicate gene IDs
    duplicates = families_df[families_df["gene"].duplicated()]["gene"].unique()

    if len(duplicates) > 0:
        raise ValueError(
            f"Duplicate gene IDs found in your clustering: {', '.join(duplicates)}"
        )

    return (
        families_df[["family", "representative", "gene", "is_frag"]],
        families_df["is_frag"].any(),
    )


def read_clustering(
    pangenome: Pangenome,
    families_tsv_path: Path,
    infer_singleton: bool = False,
    code: int = 11,
    cpu: int = 1,
    tmpdir: Path = None,
    keep_tmp: bool = False,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Get the pangenome information, the gene families and the genes with an associated gene family.
    Reads a families tsv file from mmseqs2 output and adds the gene families and the genes to the pangenome.

    :param pangenome: Input Pangenome
    :param families_tsv_path: Clustering results path
    :param infer_singleton: creates a new family for each gene with no associated family
    :param code: Genetic code used for sequence translation.
    :param cpu: Number of CPU cores to use for clustering.
    :param tmpdir: Path to a temporary directory for intermediate files.
    :param keep_tmp: Keep temporary files (useful for debugging).
    :param force: force to write in the pangenome
    :param disable_bar: Allow to disable progress bar
    """
    check_pangenome_former_clustering(pangenome, force)

    if pangenome.status["geneSequences"] == "No":
        need_gene_sequences = False
    else:
        need_gene_sequences = True

    check_pangenome_info(
        pangenome,
        need_annotations=True,
        need_gene_sequences=need_gene_sequences,
        disable_bar=disable_bar,
    )

    families_df, frag = read_clustering_file(families_tsv_path)

    nb_gene_with_fam = 0
    local_dict = mk_local_to_gene(pangenome)

    def get_gene_obj(identifier):
        try:
            gene_obj = pangenome.get_gene(identifier)
        except KeyError:
            gene_obj = local_dict.get(identifier)
        return gene_obj

    for _, row in tqdm(
        families_df.iterrows(),
        total=families_df.shape[0],
        unit="line",
        disable=disable_bar,
    ):

        fam_id, reprez_id, gene_id, is_frag = (
            str(row["family"]),
            str(row["representative"]),
            str(row["gene"]),
            bool(row["is_frag"]),
        )

        gene = get_gene_obj(gene_id)

        if gene is not None:
            nb_gene_with_fam += 1

            try:
                fam = pangenome.get_gene_family(fam_id)

            except KeyError:  # Family not found so create and add
                fam = GeneFamily(pangenome.max_fam_id, fam_id)
                representative_gene = get_gene_obj(reprez_id)
                if representative_gene is None:
                    raise KeyError(
                        f"The gene {reprez_id} associated to family {fam_id} from the clustering file is not found in pangenome."
                    )

                fam.representative = representative_gene

                pangenome.add_gene_family(fam)
            gene.is_fragment = is_frag
            fam.add(gene)
        else:
            raise KeyError(
                f"The gene {gene_id} associated to family {fam_id} from the clustering file is not found in pangenome."
            )

    if (
        nb_gene_with_fam < pangenome.number_of_genes
    ):  # not all genes have an associated cluster
        if nb_gene_with_fam == 0:
            raise Exception(
                "No gene ID in the cluster file matched any gene ID from the annotation step."
                " Please ensure that the annotations that you loaded previously and the clustering results "
                "that you have used the same gene IDs. If you use .gff files it is the identifier stored in"
                " the field 'ID'. If you use .gbff files it is the identifier stored in 'locus_tag'."
            )
        else:
            if infer_singleton:
                infer_singletons(pangenome)
            else:
                raise Exception(
                    f"Some genes ({pangenome.number_of_genes - nb_gene_with_fam}) were not associated with a cluster. "
                    f"You can either update your cluster file to ensure each gene has a cluster assignment, "
                    f"or use the '--infer_singletons' option to automatically infer a cluster for each non-clustered gene."
                )
    if pangenome.status["geneSequences"] == "No":
        logging.getLogger("PPanGGOLiN").info(
            "The pangenome has no gene sequences so it is not possible to extract sequence of family representatives."
        )
    else:
        get_family_representative_sequences(pangenome, code, cpu, tmpdir, keep_tmp)

    pangenome.status["genesClustered"] = "Computed"
    if frag:  # if there was fragment information in the file.
        pangenome.status["defragmented"] = "Computed"
    pangenome.status["geneFamilySequences"] = "Computed"
    pangenome.parameters["cluster"] = {}
    pangenome.parameters["cluster"]["# read_clustering_from_file"] = True
    pangenome.parameters["cluster"]["infer_singletons"] = infer_singleton


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    if args.clusters is None:
        if args.infer_singletons is True:
            logging.getLogger("PPanGGOLiN").warning(
                "--infer_singletons option is not compatible with clustering "
                "creation. To infer singleton you should give a clustering"
            )

        clustering(
            pangenome,
            args.tmpdir,
            args.cpu,
            defrag=not args.no_defrag,
            code=args.translation_table,
            coverage=args.coverage,
            identity=args.identity,
            mode=args.mode,
            force=args.force,
            disable_bar=args.disable_prog_bar,
            keep_tmp_files=args.keep_tmp,
        )
        logging.getLogger("PPanGGOLiN").info("Done with the clustering")
    else:
        if None in [
            args.tmpdir,
            args.cpu,
            args.no_defrag,
            args.translation_table,
            args.coverage,
            args.identity,
            args.mode,
        ]:
            logging.getLogger("PPanGGOLiN").warning(
                "You are using an option compatible only with clustering creation."
            )
        read_clustering(
            pangenome,
            args.clusters,
            args.infer_singletons,
            args.translation_table,
            args.cpu,
            args.tmpdir,
            args.keep_tmp,
            args.force,
            disable_bar=args.disable_prog_bar,
        )
        logging.getLogger("PPanGGOLiN").info("Done reading the cluster file")
    write_pangenome(
        pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar
    )


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser(
        "cluster", formatter_class=argparse.RawTextHelpFormatter
    )
    parser_clust(parser)
    return parser


def parser_clust(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of cluster command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(
        title="Required arguments",
        description="One of the following arguments is required :",
    )
    required.add_argument(
        "-p", "--pangenome", required=False, type=Path, help="The pangenome .h5 file"
    )
    clust = parser.add_argument_group(title="Clustering arguments")
    clust.add_argument(
        "--identity",
        required=False,
        type=restricted_float,
        default=0.8,
        help="Minimal identity percent for two proteins to be in the same cluster",
    )
    clust.add_argument(
        "--coverage",
        required=False,
        type=restricted_float,
        default=0.8,
        help="Minimal coverage of the alignment for two proteins to be in the same cluster",
    )
    clust.add_argument(
        "--mode",
        required=False,
        default="1",
        choices=["0", "1", "2", "3"],
        help="the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component),"
        " 2: CD-HIT-like, 3: CD-HIT-like (lowmem)",
    )
    clust.add_argument(
        "--no_defrag",
        required=False,
        default=False,
        action="store_true",
        help="DO NOT Use the defragmentation strategy to link potential fragments "
        "with their original gene family.",
    )

    read = parser.add_argument_group(title="Read clustering arguments")
    read.add_argument(
        "--clusters",
        required=False,
        type=Path,
        help="A tab-separated list containing the result of a clustering. One line per gene. "
        "First column is cluster ID, and second is gene ID",
    )
    read.add_argument(
        "--infer_singletons",
        required=False,
        action="store_true",
        help="When reading a clustering result with --clusters, if a gene is not in the provided file"
        " it will be placed in a cluster where the gene is the only member.",
    )
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument(
        "--translation_table",
        required=False,
        default="11",
        help="Translation table (genetic code) to use.",
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
        "--tmpdir",
        required=False,
        type=Path,
        default=Path(tempfile.gettempdir()),
        help="directory for storing temporary files",
    )
    optional.add_argument(
        "--keep_tmp",
        required=False,
        default=False,
        action="store_true",
        help="Keeping temporary files (useful for debugging).",
    )


if __name__ == "__main__":
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser_clust(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
