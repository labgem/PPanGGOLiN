#!/usr/bin/env python3

# default libraries
import logging
from pathlib import Path
from typing import Dict, Any, Iterator, Set, List, Tuple, Optional
from collections import defaultdict

# installed libraries
from tqdm import tqdm
import tables

# local libraries
from ppanggolin.genome import Organism, Gene, RNA, Contig
from ppanggolin.pangenome import Pangenome
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.region import Region, Spot, Module
from ppanggolin.metadata import Metadata
from ppanggolin.utils import write_compressed_or_not


class Genedata:
    """
    This is a general class storing unique gene-related data to be written in a specific
    genedata table
    """

    def __init__(
        self,
        start: int,
        stop: int,
        strand: str,
        gene_type: str,
        position: int,
        name: str,
        product: str,
        genetic_code: int,
        coordinates: List[Tuple[int]] = None,
    ):
        """Constructor method

        :param start: Gene start position
        :param stop: Gene stop position
        :param strand: Associated strand
        :param gene_type: Gene type
        :param position: Position of the gene on its contig
        :param name: Name of the feature
        :param product: Associated product
        :param genetic_code: associated genetic code, if any
        """
        self.start = start
        self.stop = stop
        self.strand = strand
        self.gene_type = gene_type
        self.position = position
        self.name = name
        self.product = product
        self.genetic_code = genetic_code
        self.has_joined_coordinates = len(coordinates) > 1
        self.coordinates = coordinates

    def __eq__(self, other):
        return (
            self.start == other.start
            and self.stop == other.stop
            and self.strand == other.strand
            and self.gene_type == other.gene_type
            and self.position == other.position
            and self.name == other.name
            and self.product == other.product
            and self.genetic_code == other.genetic_code
            and self.coordinates == other.coordinates
        )

    def __hash__(self):
        return hash(
            (
                self.start,
                self.stop,
                self.strand,
                self.gene_type,
                self.position,
                self.name,
                self.product,
                self.genetic_code,
                tuple(self.coordinates),
            )
        )


def get_number_of_organisms(pangenome: Pangenome) -> int:
    """Standalone function to get the number of organisms in a pangenome

    :param pangenome: Annotated pangenome

    :return: Number of organisms in the pangenome
    """
    if hasattr(pangenome, "file"):
        filename = pangenome.file
    else:
        raise FileNotFoundError(
            "The provided pangenome does not have an associated .h5 file"
        )
    h5f = tables.open_file(filename, "r")
    annotations = h5f.root.annotations

    table = annotations.genes
    org_set = set()
    for org in read_chunks(table, column="genome"):
        org_set.add(org)
    h5f.close()
    return len(org_set)


def get_status(pangenome: Pangenome, pangenome_file: Path):
    """
    Checks which elements are already present in the file.

    :param pangenome: Blank pangenome
    :param pangenome_file: path to the pangenome file
    """
    h5f = tables.open_file(pangenome_file.absolute().as_posix(), "r")
    logging.getLogger("PPanGGOLiN").info("Getting the current pangenome status")
    status_group = h5f.root.status
    if status_group._v_attrs.genomesAnnotated:
        pangenome.status["genomesAnnotated"] = "inFile"
    if status_group._v_attrs.genesClustered:
        pangenome.status["genesClustered"] = "inFile"
    if status_group._v_attrs.geneSequences:
        pangenome.status["geneSequences"] = "inFile"
    if status_group._v_attrs.geneFamilySequences:
        pangenome.status["geneFamilySequences"] = "inFile"
    if status_group._v_attrs.NeighborsGraph:
        pangenome.status["neighborsGraph"] = "inFile"

    if hasattr(status_group._v_attrs, "version"):
        pangenome.status["ppanggolin_version"] = str(status_group._v_attrs.version)
    else:
        logging.getLogger("PPanGGOLiN").error(
            f"The provided pangenome file {pangenome_file} does not have a version stored in its status."
            " This issue may indicate that the file is corrupted."
        )
        pangenome.status["ppanggolin_version"] = None

    if status_group._v_attrs.Partitioned:
        pangenome.status["partitioned"] = "inFile"

    if (
        hasattr(status_group._v_attrs, "predictedRGP")
        and status_group._v_attrs.predictedRGP
    ):
        pangenome.status["predictedRGP"] = "inFile"

    if hasattr(status_group._v_attrs, "spots") and status_group._v_attrs.spots:
        pangenome.status["spots"] = "inFile"

    if hasattr(status_group._v_attrs, "modules") and status_group._v_attrs.modules:
        pangenome.status["modules"] = "inFile"

    if hasattr(status_group._v_attrs, "metadata") and status_group._v_attrs.metadata:
        metastatus = status_group.metastatus
        metasources = status_group.metasources
        for attr in metastatus._v_attrs._f_list():
            pangenome.status["metadata"][attr] = "inFile"
            pangenome.status["metasources"][attr] = metasources._v_attrs[attr]

    if "/info" in h5f:
        info_group = h5f.root.info
        pangenome.parameters = info_group._v_attrs.parameters
    h5f.close()


def read_chunks(table: tables.Table, column: str = None, chunk: int = 10000):
    """
    Reading entirely the provided table (or column if specified) chunk per chunk to limit RAM usage.

    :param table:
    :param column:
    :param chunk:
    """
    for i in range(0, table.nrows, chunk):
        yield from table.read(start=i, stop=i + chunk, field=column)


def read_genedata(h5f: tables.File) -> Dict[int, Genedata]:
    """
    Reads the genedata table and returns a genedata_id2genedata dictionary

    :param h5f: the hdf5 file handler

    :return: dictionary linking genedata to the genedata identifier

    :raises KeyError: If a Genedata entry with joined coordinates is not found in the annotations.joinCoordinates table.
    """

    genedata_id_to_coordinates = read_join_coordinates(h5f)

    table = h5f.root.annotations.genedata
    genedata_id2genedata = {}
    for row in read_chunks(table, chunk=20000):
        start = int(row["start"])
        stop = int(row["stop"])

        if (
            "has_joined_coordinates" in row.dtype.names
            and row["has_joined_coordinates"]
        ):
            # manage gene with joined coordinates if the info exists

            try:
                coordinates = genedata_id_to_coordinates[row["genedata_id"]]
            except KeyError:
                raise KeyError(
                    f'Genedata {row["genedata_id"]} is supposed to have joined '
                    "coordinates but is not found in annotations.joinCoordinates table"
                )
        else:
            coordinates = [(start, stop)]

        genedata = Genedata(
            start=start,
            stop=stop,
            strand=row["strand"].decode(),
            gene_type=row["gene_type"].decode(),
            position=int(row["position"]),
            name=row["name"].decode(),
            product=row["product"].decode(),
            genetic_code=int(row["genetic_code"]),
            coordinates=coordinates,
        )

        genedata_id = row["genedata_id"]
        genedata_id2genedata[genedata_id] = genedata

    return genedata_id2genedata


def read_join_coordinates(h5f: tables.File) -> Dict[str, List[Tuple[int, int]]]:
    """
    Read join coordinates from a HDF5 file and return a dictionary mapping genedata_id to coordinates.

    :param h5f: An HDF5 file object.
    :return: A dictionary mapping genedata_id to a list of tuples representing start and stop coordinates.
    """
    genedata_id_to_coordinates = defaultdict(list)

    if not hasattr(h5f.root.annotations, "joinedCoordinates"):
        # then the pangenome file has no joined annotations
        # or has been made before the joined annotations coordinates
        return {}

    table = h5f.root.annotations.joinedCoordinates

    for row in read_chunks(table, chunk=20000):
        genedata_id = row["genedata_id"]

        genedata_id_to_coordinates[genedata_id].append(
            (int(row["coordinate_rank"]), int(row["start"]), int(row["stop"]))
        )

    # sort coordinate by their rank
    genedata_id_to_sorted_coordinates = {}
    for genedata_id, coordinates in genedata_id_to_coordinates.items():
        sorted_coordinates = [
            (start, stop) for rank, start, stop in sorted(coordinates)
        ]
        genedata_id_to_sorted_coordinates[genedata_id] = sorted_coordinates

    return genedata_id_to_sorted_coordinates


def read_sequences(h5f: tables.File) -> dict:
    """
    Reads the sequences table and returns a sequence id to sequence dictionary
    :param h5f: the hdf5 file handler
    :return: dictionary linking sequences to the seq identifier
    """
    table = h5f.root.annotations.sequences
    seqid2seq = {}
    for row in read_chunks(table, chunk=20000):
        seqid2seq[row["seqid"]] = row["dna"].decode()
    return seqid2seq


def get_non_redundant_gene_sequences_from_file(
    pangenome_filename: str, output: Path, add: str = "", disable_bar: bool = False
):
    """
    Writes the non-redundant CDS sequences of the Pangenome object to a File object that can be filtered or not by a list of CDS,
    and adds the eventual str 'add' in front of the identifiers. Loads the sequences from a .h5 pangenome file.

    :param pangenome_filename: Name of the pangenome file
    :param output: Path to the output file
    :param add: Add a prefix to sequence header
    :param disable_bar: disable progress bar

    """

    logging.getLogger("PPanGGOLiN").info(
        f"Extracting and writing non redundant CDS sequences from {pangenome_filename}"
        f" to {output.absolute()}"
    )

    with tables.open_file(pangenome_filename, "r", driver_core_backing_store=0) as h5f:

        # get a dictionary mapping seqid to cds_name
        # seqid are uniq and can have multiple cds name.
        # We just want one of the cds name to have non-redundant fasta sequences
        seqid2cds_name = {}
        for row in read_chunks(h5f.root.annotations.geneSequences, chunk=20000):
            # Read the table chunk per chunk otherwise RAM dies on big pangenomes
            seqid2cds_name[row["seqid"]] = row["gene"].decode()

        table = h5f.root.annotations.sequences
        with open(output, "w") as file_obj:
            for row in tqdm(
                read_chunks(table, chunk=20000),
                total=table.nrows,
                unit="gene",
                disable=disable_bar,
            ):
                cds_name = seqid2cds_name[row["seqid"]]
                file_obj.write(f">{add}{cds_name}\n")
                file_obj.write(f'{row["dna"].decode()}\n')


def write_gene_sequences_from_pangenome_file(
    pangenome_filename: str,
    output: Path,
    list_cds: Optional[Iterator] = None,
    add: str = "",
    compress: bool = False,
    disable_bar: bool = False,
):
    """
    Writes the CDS sequences of the Pangenome object to a File object that can be filtered or not by a list of CDS,
    and adds the eventual str 'add' in front of the identifiers. Loads the sequences from a .h5 pangenome file.

    :param pangenome_filename: Name of the pangenome file
    :param output: Path to the sequences file
    :param list_cds: An iterable object of CDS
    :param add: Add a prefix to sequence header
    :param compress: Compress the output file
    :param disable_bar: Prevent to print disable progress bar
    """
    logging.getLogger("PPanGGOLiN").info(
        f"Extracting and writing CDS sequences from a {pangenome_filename} "
        "file to a fasta file..."
    )
    with tables.open_file(pangenome_filename, "r", driver_core_backing_store=0) as h5f:
        table = h5f.root.annotations.geneSequences
        list_cds = set(list_cds) if list_cds is not None else None
        seqid2seq = read_sequences(h5f)
        with write_compressed_or_not(output, compress) as file_obj:
            for row in tqdm(
                read_chunks(table, chunk=20000),
                total=table.nrows,
                unit="gene",
                disable=disable_bar,
            ):
                # Read the table chunk per chunk otherwise RAM dies on big pangenomes
                name_cds = row["gene"].decode()
                if row["type"] == b"CDS" and (list_cds is None or name_cds in list_cds):
                    file_obj.write(">" + add + name_cds + "\n")
                    file_obj.write(seqid2seq[row["seqid"]] + "\n")
    logging.getLogger("PPanGGOLiN").debug(
        "Gene sequences from pangenome file was written to "
        f"{output.absolute()}{'.gz' if compress else ''}"
    )


def read_rgp_genes_from_pangenome_file(h5f: tables.File) -> Set[bytes]:
    """
    Retrieves a list of RGP genes from the pangenome file.

    :param h5f: The open HDF5 pangenome file containing RGP gene data.
    :return: A list of gene names (as bytes) from the RGP.
    """
    rgp_genes = {row["gene"] for row in read_chunks(h5f.root.RGP, chunk=20000)}

    return rgp_genes


def get_families_from_genes(h5f: tables.File, genes: Set[bytes]) -> Set[bytes]:
    """
    Retrieves gene families associated with a specified set of genes from the pangenome file.

    :param h5f: The open HDF5 pangenome file containing gene family data.
    :param genes: A set of gene names (as bytes) for which to retrieve the associated families.
    :return: A set of gene family names (as bytes) associated with the specified genes.
    """

    families = set()
    for row in read_chunks(h5f.root.geneFamilies, chunk=20000):
        if row["gene"] in genes:
            families.add(row["geneFam"])

    return families


def read_module_families_from_pangenome_file(
    h5f: tables.File, module_name: str
) -> Set[bytes]:
    """
    Retrieves gene families associated with a specified module from the pangenome file.


    :param h5f: The open HDF5 pangenome file containing module data.
    :param module_name: The name of the module (as a string). The module ID is extracted from
                        the name by removing the "module_" prefix.
    :return: A set of gene family names (as bytes) associated with the specified module.
    """

    family_to_write = set()
    module_id = int(module_name[len("module_") :])
    module_table = h5f.root.modules

    for row in read_chunks(module_table, chunk=20000):
        if row["module"] == module_id:
            family_to_write.add(row["geneFam"])

    return family_to_write


def get_families_matching_partition(h5f: tables.File, partition: str) -> Set[bytes]:
    """
    Retrieves gene families that match the specified partition.

    :param h5f: The open HDF5 pangenome file containing gene family information.
    :param partition: The partition name (as a string). If "all", all gene families are included.
                      Otherwise, it filters by the first letter of the partition.
    :return: A set of gene family names (as bytes) that match the partition criteria.
    """

    family_to_write = set()

    gene_fam_info_table = h5f.root.geneFamiliesInfo
    parition_first_letter = partition[0].upper()

    for row in read_chunks(gene_fam_info_table, chunk=20000):

        if partition == "all" or row["partition"].decode().startswith(
            parition_first_letter
        ):
            family_to_write.add(row["name"])

    return family_to_write


def get_genes_from_families(h5f: tables.File, families: List[bytes]) -> Set[bytes]:
    """
    Retrieves a set of genes that belong to the specified families.

    This function reads the gene family data from an HDF5 pangenome file and returns
    a set of genes that are part of the given list of gene families.

    :param h5f: The open HDF5 pangenome file containing gene family data.
    :param families: A list of gene families (as bytes) to filter genes by.
    :return: A set of genes (as bytes) that belong to the specified families.
    """

    matching_genes = set()

    gene_fam_table = h5f.root.geneFamilies

    for row in read_chunks(gene_fam_table, chunk=20000):
        if row["geneFam"] in families:
            matching_genes.add(row["gene"])

    return matching_genes


def get_seqid_to_genes(
    h5f: tables.File,
    genes: Set[bytes],
    get_all_genes: bool = False,
    disable_bar: bool = False,
) -> Dict[int, List[str]]:
    """
    Creates a mapping of sequence IDs to gene names.

    :param h5f: The open HDF5 pangenome file containing gene sequence data.
    :param genes: A list of gene names to include in the mapping (if `get_all_genes` is False).
    :param get_all_genes: Boolean flag to indicate if all genes should be included in the mapping.
                          If set to True, all genes will be added regardless of the `genes` parameter.
    :param disable_bar: Boolean flag to disable the progress bar if set to True.
    :return: A dictionary mapping sequence IDs (integers) to lists of gene names (strings).
    """

    seq_id_to_genes = defaultdict(list)
    gene_seq_table = h5f.root.annotations.geneSequences
    match_count = 0
    for row in tqdm(
        read_chunks(gene_seq_table, chunk=20000),
        total=gene_seq_table.nrows,
        unit="gene",
        disable=disable_bar,
    ):

        if get_all_genes or row["gene"] in genes:
            seq_id_to_genes[row["seqid"]].append(row["gene"].decode())
            match_count += 1

    assert get_all_genes or match_count == len(
        genes
    ), f"Number of sequences found ({match_count}) does not match the number of expected genes {len(genes)}."

    return seq_id_to_genes


def write_genes_seq_from_pangenome_file(
    h5f: tables.File,
    outpath: Path,
    compress: bool,
    seq_id_to_genes: Dict[int, List[str]],
    disable_bar: bool,
):
    """
    Writes gene sequences from the pangenome file to an output file.

    Only sequences whose IDs match the ones in `seq_id_to_genes` will be written.

    :param h5f: The open HDF5 pangenome file containing sequence data.
    :param outpath: The path to the output file where sequences will be written.
    :param compress: Boolean flag to indicate whether output should be compressed.
    :param seq_id_to_genes: A dictionary mapping sequence IDs to lists of gene names.
    :param disable_bar: Boolean flag to disable the progress bar if set to True.
    """

    with write_compressed_or_not(file_path=outpath, compress=compress) as file_obj:

        seq_table = h5f.root.annotations.sequences

        with tqdm(
            total=len(seq_id_to_genes), unit="sequence", disable=disable_bar
        ) as pbar:

            for row in read_chunks(table=seq_table, chunk=20000):
                if row["seqid"] in seq_id_to_genes:

                    for seq_name in seq_id_to_genes[row["seqid"]]:
                        file_obj.write(f">{seq_name}\n")
                        file_obj.write(row["dna"].decode() + "\n")

                    pbar.update(1)


def get_gene_to_genome(h5f: tables.File) -> Dict[bytes, bytes]:
    """
    Generates a mapping between gene IDs and their corresponding genome.

    :param h5f: The open HDF5 pangenome file containing contig and gene annotations.
    :return: A dictionary mapping gene IDs to genome names.
    """

    contig_id_to_genome = {
        row["ID"]: row["genome"]
        for row in read_chunks(h5f.root.annotations.contigs, chunk=20000)
    }

    gene_to_genome = {
        row["ID"]: contig_id_to_genome[row["contig"]]
        for row in read_chunks(h5f.root.annotations.genes, chunk=20000)
    }

    return gene_to_genome


def get_family_to_genome_count(h5f: tables.File) -> Dict[bytes, int]:
    """
    Computes the number of unique genomes associated with each gene family.

    :param h5f: The open HDF5 pangenome file containing contig, gene, and gene family data.
    :return: A dictionary mapping gene family names (as bytes) to the count of unique genomes.
    """

    contig_id_to_genome = {
        row["ID"]: row["genome"]
        for row in read_chunks(h5f.root.annotations.contigs, chunk=20000)
    }

    gene_to_genome = {
        row["ID"]: contig_id_to_genome[row["contig"]]
        for row in read_chunks(h5f.root.annotations.genes, chunk=20000)
    }

    family_to_genomes = defaultdict(set)
    for row in read_chunks(h5f.root.geneFamilies, chunk=20000):
        family_to_genomes[row["geneFam"]].add(gene_to_genome[row["gene"]])

    family_to_genome_count = {
        fam: len(genomes) for fam, genomes in family_to_genomes.items()
    }

    return family_to_genome_count


def get_soft_core_families(h5f: tables.File, soft_core: float) -> Set[bytes]:
    """
    Identifies gene families that are present in at least a specified proportion of genomes.

    :param h5f: The open HDF5 pangenome file containing gene family and genome data.
    :param soft_core: The proportion of genomes (between 0 and 1) that a gene family must be present in
                      to be considered a soft core family.
    :return: A set of gene family names (as bytes) that are classified as soft core.
    """
    family_to_genome_count = get_family_to_genome_count(h5f)
    pangenome_info = read_info(h5f)
    genome_count = pangenome_info["Content"]["Genomes"]

    genome_count_threshold = genome_count * soft_core
    soft_core_families = {
        family
        for family, fam_genome_count in family_to_genome_count.items()
        if fam_genome_count >= genome_count_threshold
    }
    return soft_core_families


def write_fasta_gene_fam_from_pangenome_file(
    pangenome_filename: str,
    output: Path,
    family_filter: str,
    soft_core: float = 0.95,
    compress: bool = False,
    disable_bar=False,
):
    """
    Write representative nucleotide sequences of gene families

    :param pangenome: Pangenome object with gene families sequences
    :param output: Path to output directory
    :param gene_families: Selected partition of gene families
    :param soft_core: Soft core threshold to use
    :param compress: Compress the file in .gz
    :param disable_bar: Disable progress bar
    """

    outpath = output / f"{family_filter}_nucleotide_families.fasta"

    with tables.open_file(pangenome_filename, "r", driver_core_backing_store=0) as h5f:

        if family_filter in ["all", "persistent", "shell", "cloud"]:
            family_to_write = get_families_matching_partition(h5f, family_filter)

        elif family_filter.startswith("module_"):
            family_to_write = read_module_families_from_pangenome_file(
                h5f, module_name=family_filter
            )

        elif family_filter == "rgp":
            rgp_genes = read_rgp_genes_from_pangenome_file(h5f)
            family_to_write = get_families_from_genes(h5f, rgp_genes)

        elif family_filter in ["softcore", "core"]:
            if family_filter == "core":
                soft_core = 1.0

            family_to_write = get_soft_core_families(h5f, soft_core)

        if len(family_to_write) == 0:
            logging.getLogger("PPanGGOLiN").warning(
                f"No families matching filter {family_filter}."
            )
            return

        seq_id_to_genes = get_seqid_to_genes(h5f, set(family_to_write))

        write_genes_seq_from_pangenome_file(
            h5f, outpath, compress, seq_id_to_genes, disable_bar=disable_bar
        )

    logging.getLogger("PPanGGOLiN").info(
        "Done writing the representative nucleotide sequences "
        f"of the gene families : '{outpath}{'.gz' if compress else ''}"
    )


def write_fasta_prot_fam_from_pangenome_file(
    pangenome_filename: str,
    output: Path,
    family_filter: str,
    soft_core: float = 0.95,
    compress: bool = False,
    disable_bar=False,
):
    """
    Write representative amino acid sequences of gene families.

    :param pangenome: Pangenome object with gene families sequences
    :param output: Path to output directory
    :param prot_families: Selected partition of protein families
    :param soft_core: Soft core threshold to use
    :param compress: Compress the file in .gz
    :param disable_bar: Disable progress bar
    """

    outpath = output / f"{family_filter}_protein_families.faa"

    partition_filter = False
    family_to_write = []

    with (
        tables.open_file(pangenome_filename, "r", driver_core_backing_store=0) as h5f,
        write_compressed_or_not(outpath, compress) as fasta,
    ):

        if family_filter in ["all", "persistent", "shell", "cloud"]:
            partition_filter = True
            parition_first_letter = family_filter[0].upper()

        elif family_filter == "rgp":
            rgp_genes = read_rgp_genes_from_pangenome_file(h5f)
            family_to_write = get_families_from_genes(h5f, rgp_genes)

        elif family_filter.startswith("module_"):
            family_to_write = read_module_families_from_pangenome_file(
                h5f, module_name=family_filter
            )

        elif family_filter in ["softcore", "core"]:

            soft_core_to_apply = 1.0 if family_filter == "core" else soft_core

            family_to_write = get_soft_core_families(h5f, soft_core=soft_core_to_apply)

        gene_fam_info_table = h5f.root.geneFamiliesInfo

        for row in tqdm(
            read_chunks(gene_fam_info_table, chunk=20000),
            total=gene_fam_info_table.nrows,
            unit="family",
            disable=disable_bar,
        ):

            partition_match = partition_filter and (
                family_filter == "all"
                or row["partition"].decode().startswith(parition_first_letter)
            )
            family_match = row["name"] in family_to_write

            if partition_match or family_match:

                fasta.write(f">{row['name'].decode()}\n")
                fasta.write(row["protein"].decode() + "\n")

    logging.getLogger("PPanGGOLiN").info(
        f"Done writing the representative amino acid sequences of the gene families:"
        f"'{outpath}{'.gz' if compress else ''}'"
    )


def write_genes_from_pangenome_file(
    pangenome_filename: str,
    output: Path,
    gene_filter: str,
    soft_core: float = 0.95,
    compress: bool = False,
    disable_bar=False,
):
    """
    Write representative nucleotide sequences of gene families

    :param pangenome: Pangenome object with gene families sequences
    :param output: Path to output directory
    :param gene_families: Selected partition of gene families
    :param soft_core: Soft core threshold to use
    :param compress: Compress the file in .gz
    :param disable_bar: Disable progress bar
    """

    outpath = output / f"{gene_filter}_genes.fna"
    get_all_genes = False

    with tables.open_file(pangenome_filename, "r", driver_core_backing_store=0) as h5f:

        if gene_filter in [
            "persistent",
            "shell",
            "cloud",
            "softcore",
            "core",
        ] or gene_filter.startswith("module_"):

            if gene_filter.startswith("module_"):
                families = read_module_families_from_pangenome_file(
                    h5f, module_name=gene_filter
                )

            elif gene_filter in ["softcore", "core"]:

                soft_core_to_apply = 1.0 if gene_filter == "core" else soft_core
                families = get_soft_core_families(h5f, soft_core=soft_core_to_apply)

            else:
                families = get_families_matching_partition(h5f, gene_filter)

            genes_to_write = get_genes_from_families(h5f, families)

        elif gene_filter == "rgp":
            genes_to_write = read_rgp_genes_from_pangenome_file(h5f)

        elif gene_filter == "all":
            genes_to_write = set()
            get_all_genes = True

        seq_id_to_genes = get_seqid_to_genes(
            h5f, genes_to_write, get_all_genes=get_all_genes, disable_bar=disable_bar
        )

        write_genes_seq_from_pangenome_file(
            h5f, outpath, compress, seq_id_to_genes, disable_bar=disable_bar
        )

    logging.getLogger("PPanGGOLiN").info(
        "Done writing the representative nucleotide sequences "
        f"of the gene families : '{outpath}{'.gz' if compress else ''}"
    )


def read_graph(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read information about graph in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without graph information
    :param h5f: Pangenome HDF5 file with graph information
    :param disable_bar: Disable the progress bar
    """
    table = h5f.root.edges

    if pangenome.status["genomesAnnotated"] not in [
        "Computed",
        "Loaded",
    ] or pangenome.status["genesClustered"] not in ["Computed", "Loaded"]:
        raise Exception(
            "It's not possible to read the graph "
            "if the annotations and the gene families have not been loaded."
        )
    for row in tqdm(
        read_chunks(table, chunk=20000),
        total=table.nrows,
        unit="contig adjacency",
        disable=disable_bar,
    ):
        source = pangenome.get_gene(row["geneSource"].decode())
        target = pangenome.get_gene(row["geneTarget"].decode())
        pangenome.add_edge(source, target)
    pangenome.status["neighborsGraph"] = "Loaded"


def read_gene_families(
    pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False
):
    """
    Read gene families in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without gene families
    :param h5f: Pangenome HDF5 file with gene families information
    :param disable_bar: Disable the progress bar
    """
    table = h5f.root.geneFamilies

    link = (
        True
        if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"]
        else False
    )

    for row in tqdm(
        read_chunks(table, chunk=20000),
        total=table.nrows,
        unit="gene family",
        disable=disable_bar,
    ):
        try:
            fam = pangenome.get_gene_family(name=row["geneFam"].decode())
        except KeyError:
            fam = GeneFamily(
                family_id=pangenome.max_fam_id, name=row["geneFam"].decode()
            )
            pangenome.add_gene_family(fam)
        if link:  # linking if we have loaded the annotations
            gene_obj = pangenome.get_gene(row["gene"].decode())
        else:  # else, no
            gene_obj = Gene(row["gene"].decode())
        fam.add(gene_obj)
    pangenome.status["genesClustered"] = "Loaded"


def read_gene_families_info(
    pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False
):
    """
    Read information about gene families in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without gene families information
    :param h5f: Pangenome HDF5 file with gene families information
    :param disable_bar: Disable the progress bar
    """
    table = h5f.root.geneFamiliesInfo

    for row in tqdm(
        read_chunks(table, chunk=20000),
        total=table.nrows,
        unit="gene family",
        disable=disable_bar,
    ):
        fam = pangenome.get_gene_family(row["name"].decode())
        fam.partition = row["partition"].decode()
        fam.add_sequence(row["protein"].decode())

    if h5f.root.status._v_attrs.Partitioned:
        pangenome.status["partitioned"] = "Loaded"
    if h5f.root.status._v_attrs.geneFamilySequences:
        pangenome.status["geneFamilySequences"] = "Loaded"


def read_gene_sequences(
    pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False
):
    """
    Read gene sequences in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without gene sequence associate to gene
    :param h5f: Pangenome HDF5 file with gene sequence associate to gene
    :param disable_bar: Disable the progress bar
    """
    if pangenome.status["genomesAnnotated"] not in ["Computed", "Loaded"]:
        raise Exception(
            "It's not possible to read the pangenome gene dna sequences "
            "if the annotations have not been loaded."
        )
    table = h5f.root.annotations.geneSequences

    seqid2seq = read_sequences(h5f)
    for row in tqdm(
        read_chunks(table, chunk=20000),
        total=table.nrows,
        unit="gene",
        disable=disable_bar,
    ):
        gene = pangenome.get_gene(row["gene"].decode())
        gene.add_sequence(seqid2seq[row["seqid"]])
    pangenome.status["geneSequences"] = "Loaded"


def read_rgp(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read region of genomic plasticity in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without RGP
    :param h5f: Pangenome HDF5 file with RGP computed
    :param disable_bar: Disable the progress bar
    """
    if pangenome.status["genomesAnnotated"] not in [
        "Computed",
        "Loaded",
    ] or pangenome.status["genesClustered"] not in ["Computed", "Loaded"]:
        raise Exception(
            "It's not possible to read the RGP "
            "if the annotations and the gene families have not been loaded."
        )
    table = h5f.root.RGP

    for row in tqdm(
        read_chunks(table, chunk=20000),
        total=table.nrows,
        unit="region",
        disable=disable_bar,
    ):
        try:
            region = pangenome.get_region(row["RGP"].decode())
        except KeyError:
            region = Region(row["RGP"].decode())

            # starting from v2.2.1 score is part of RGP table in h5.
            if "score" in row.dtype.names:
                region.score = row["score"]

            pangenome.add_region(region)

        gene = pangenome.get_gene(row["gene"].decode())
        region.add(gene)
    pangenome.status["predictedRGP"] = "Loaded"


def read_spots(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read hotspots in the pangenome HDF5 file and add them to the pangenome object.

    :param pangenome: Pangenome object without spot
    :param h5f: Pangenome HDF5 file with spot computed
    :param disable_bar: Disable the progress bar
    """
    table = h5f.root.spots
    spots = {}
    curr_spot_id = None
    for row in tqdm(
        read_chunks(table, chunk=20000),
        total=table.nrows,
        unit="spot",
        disable=disable_bar,
    ):
        if curr_spot_id != int(row["spot"]):
            curr_spot_id = int(row["spot"])
            curr_spot = spots.get(curr_spot_id)
            if curr_spot is None:
                curr_spot = Spot(int(row["spot"]))
                spots[row["spot"]] = curr_spot
        region = pangenome.get_region(row["RGP"].decode())
        curr_spot.add(region)
    for spot in spots.values():
        spot.spot_2_families()
        pangenome.add_spot(spot)
    pangenome.status["spots"] = "Loaded"


def read_modules(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read modules in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without modules
    :param h5f: Pangenome HDF5 file with modules computed
    :param disable_bar: Disable the progress bar
    """
    if pangenome.status["genesClustered"] not in ["Computed", "Loaded"]:
        raise Exception(
            "It's not possible to read the modules if the gene families have not been loaded."
        )
    table = h5f.root.modules
    modules = {}  # id2mod
    for row in tqdm(
        read_chunks(table, chunk=20000),
        total=table.nrows,
        unit="module",
        disable=disable_bar,
    ):
        curr_module = modules.get(int(row["module"]))
        if curr_module is None:
            curr_module = Module(int(row["module"]))
            modules[row["module"]] = curr_module
        family = pangenome.get_gene_family(row["geneFam"].decode())
        curr_module.add(family)
    for module in modules.values():
        pangenome.add_module(module)
    pangenome.status["modules"] = "Loaded"


def read_organisms(
    pangenome: Pangenome,
    table: tables.Table,
    chunk_size: int = 20000,
    disable_bar: bool = False,
):
    """Read organism table in pangenome file to add them to the pangenome object

    :param pangenome: Pangenome object
    :param table: Organism table
    :param chunk_size: Size of the chunk reading
    :param disable_bar: Disable progress bar
    """
    for row in tqdm(
        read_chunks(table, chunk=chunk_size),
        total=table.nrows,
        unit="genome",
        disable=disable_bar,
    ):
        organism = Organism(row["name"].decode())
        pangenome.add_organism(organism)


def read_contigs(
    pangenome: Pangenome,
    table: tables.Table,
    chunk_size: int = 20000,
    disable_bar: bool = False,
):
    """Read contig table in pangenome file to add them to the pangenome object

    :param pangenome: Pangenome object
    :param table: Contig table
    :param chunk_size: Size of the chunk reading
    :param disable_bar: Disable progress bar
    """
    for row in tqdm(
        read_chunks(table, chunk=chunk_size),
        total=table.nrows,
        unit="contig",
        disable=disable_bar,
    ):
        contig = Contig(
            identifier=int(row["ID"]),
            name=row["name"].decode(),
            is_circular=row["is_circular"],
        )
        contig.length = int(row["length"])
        try:
            organism = pangenome.get_organism(row["genome"].decode())
        except KeyError:
            pass
        else:
            organism.add(contig)


def read_genes(
    pangenome: Pangenome,
    table: tables.Table,
    genedata_dict: Dict[int, Genedata],
    link: bool = True,
    chunk_size: int = 20000,
    disable_bar: bool = False,
):
    """Read genes in pangenome file to add them to the pangenome object

    :param pangenome: Pangenome object
    :param table: Genes table
    :param genedata_dict: Dictionary to link genedata with gene
    :param link: Allow to link gene to organism and contig
    :param chunk_size: Size of the chunk reading
    :param disable_bar: Disable progress bar
    """
    for row in tqdm(
        read_chunks(table, chunk=chunk_size),
        total=table.nrows,
        unit="gene",
        disable=disable_bar,
    ):
        gene = Gene(row["ID"].decode())
        genedata = genedata_dict[row["genedata_id"]]
        try:
            local = row["local"].decode()
        except ValueError:
            local = ""
        gene.fill_annotations(
            start=genedata.start,
            stop=genedata.stop,
            strand=genedata.strand,
            gene_type=genedata.gene_type,
            name=genedata.name,
            position=genedata.position,
            genetic_code=genedata.genetic_code,
            product=genedata.product,
            local_identifier=local,
            coordinates=genedata.coordinates,
        )
        gene.is_fragment = row["is_fragment"]
        if link:
            contig = pangenome.get_contig(identifier=int(row["contig"]))
            gene.fill_parents(contig.organism, contig)
            contig.add(gene)


def read_rnas(
    pangenome: Pangenome,
    table: tables.Table,
    genedata_dict: Dict[int, Genedata],
    link: bool = True,
    chunk_size: int = 20000,
    disable_bar: bool = False,
):
    """Read RNAs in pangenome file to add them to the pangenome object

    :param pangenome: Pangenome object
    :param table: RNAs table
    :param genedata_dict: Dictionary to link genedata with gene
    :param link: Allow to link gene to organism and contig
    :param chunk_size: Size of the chunk reading
    :param disable_bar: Disable progress bar
    """
    for row in tqdm(
        read_chunks(table, chunk=chunk_size),
        total=table.nrows,
        unit="gene",
        disable=disable_bar,
    ):
        rna = RNA(row["ID"].decode())
        genedata = genedata_dict[row["genedata_id"]]
        if genedata.start > genedata.stop:
            logging.warning(
                f"Wrong coordinates in RNA gene {genedata.name}: Start ({genedata.start}) should not be greater than stop ({genedata.stop}). This gene is ignored."
            )
            continue
        if genedata.start < 1 or genedata.stop < 1:
            logging.warning(
                f"Wrong coordinates in RNA gene {genedata.name}: Start ({genedata.start}) and stop ({genedata.stop}) should be greater than 0.  This gene is ignored."
            )
            continue

        rna.fill_annotations(
            start=genedata.start,
            stop=genedata.stop,
            strand=genedata.strand,
            gene_type=genedata.gene_type,
            name=genedata.name,
            product=genedata.product,
        )
        if link:
            contig = pangenome.get_contig(int(row["contig"]))
            rna.fill_parents(contig.organism, contig)
            contig.add_rna(rna)


def read_annotation(
    pangenome: Pangenome,
    h5f: tables.File,
    load_organisms: bool = True,
    load_contigs: bool = True,
    load_genes: bool = True,
    load_rnas: bool = True,
    chunk_size: int = 20000,
    disable_bar: bool = False,
):
    """
    Read annotation in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without annotation
    :param h5f: Pangenome HDF5 file with annotation
    :param load_organisms: Flag to load organisms
    :param load_contigs: Flag to load contigs
    :param load_genes: Flag to load genes
    :param load_rnas: Flag to load RNAs
    :param chunk_size: Size of chunks reading
    :param disable_bar: Disable the progress bar
    """
    annotations = h5f.root.annotations
    genedata_dict = None
    if load_organisms:
        read_organisms(
            pangenome,
            annotations.genomes,
            chunk_size=chunk_size,
            disable_bar=disable_bar,
        )

    if load_contigs:
        read_contigs(
            pangenome,
            annotations.contigs,
            chunk_size=chunk_size,
            disable_bar=disable_bar,
        )

    if load_genes:
        genedata_dict = read_genedata(h5f)
        read_genes(
            pangenome,
            annotations.genes,
            genedata_dict,
            all([load_organisms, load_contigs]),
            chunk_size=chunk_size,
            disable_bar=disable_bar,
        )
    if load_rnas:
        read_rnas(
            pangenome,
            annotations.RNAs,
            read_genedata(h5f) if genedata_dict is None else genedata_dict,
            all([load_organisms, load_contigs]),
            chunk_size=chunk_size,
            disable_bar=disable_bar,
        )
    pangenome.status["genomesAnnotated"] = "Loaded"


def create_info_dict(info_group: tables.group.Group):
    """
    Read the pangenome content

    :param info_group: group in pangenome HDF5 file containing information about pangenome
    """
    attributes = info_group._v_attrs._f_list()

    info_dict = {"Genes": int(info_group._v_attrs["numberOfGenes"])}

    if "numberOfGenomes" in attributes:
        info_dict["Genomes"] = int(info_group._v_attrs["numberOfGenomes"])

    if "numberOfClusters" in attributes:
        info_dict["Families"] = int(info_group._v_attrs["numberOfClusters"])

    if "numberOfEdges" in attributes:
        info_dict["Edges"] = int(info_group._v_attrs["numberOfEdges"])

    if "numberOfCloud" in attributes:  # then all the others are there

        persistent_stat = {
            "Family_count": int(info_group._v_attrs["numberOfPersistent"])
        }
        persistent_stat.update(info_group._v_attrs["persistentStats"])
        info_dict["Persistent"] = persistent_stat

        shell_stat = {"Family_count": int(info_group._v_attrs["numberOfShell"])}
        shell_stat.update(info_group._v_attrs["shellStats"])
        info_dict["Shell"] = shell_stat

        cloud_stat = {"Family_count": int(info_group._v_attrs["numberOfCloud"])}
        cloud_stat.update(info_group._v_attrs["cloudStats"])
        info_dict["Cloud"] = cloud_stat

        info_dict["Number_of_partitions"] = int(
            info_group._v_attrs["numberOfPartitions"]
        )

        if info_group._v_attrs["numberOfPartitions"] != 3:
            subpartition_stat = {
                f"Shell_{key}": int(val)
                for key, val in info_group._v_attrs["numberOfSubpartitions"].items()
            }
            info_dict.update(subpartition_stat)

    if "genomes_fluidity" in attributes:
        info_dict["Genomes_fluidity"] = {
            key: round(val, 3)
            for key, val in info_group._v_attrs["genomes_fluidity"].items()
        }

    if "family_fluidity" in attributes:
        info_dict["Family_fluidity"] = info_group._v_attrs["family_fluidity"]

    if "numberOfRGP" in attributes:
        info_dict["RGP"] = int(info_group._v_attrs["numberOfRGP"])

    if "numberOfSpots" in attributes:
        info_dict["Spots"] = int(info_group._v_attrs["numberOfSpots"])

    if "numberOfModules" in attributes:
        info_dict["Modules"] = {
            "Number_of_modules": int(info_group._v_attrs["numberOfModules"]),
            "Families_in_Modules": int(
                info_group._v_attrs["numberOfFamiliesInModules"]
            ),
            "Partition_composition": {
                "Persistent": info_group._v_attrs["PersistentSpecInModules"]["percent"],
                "Shell": info_group._v_attrs["ShellSpecInModules"]["percent"],
                "Cloud": info_group._v_attrs["CloudSpecInModules"]["percent"],
            },
        }
    return info_dict


def read_info(h5f):
    """
    Read the pangenome content

    :param h5f: Pangenome HDF5 file
    """
    if "/info" in h5f:
        info_group = h5f.root.info
        content = create_info_dict(info_group)
        return {"Content": content}


def read_metadata(
    pangenome: Pangenome,
    h5f: tables.File,
    metatype: str,
    sources: Set[str] = None,
    disable_bar: bool = False,
):
    """Read metadata to add them to the pangenome object

    :param pangenome: Pangenome object
    :param h5f: Pangenome file
    :param metatype: Object type to associate metadata
    :param sources: Source name of metadata
    :param disable_bar: Disable progress bar
    """
    metadata_group = h5f.root.metadata._f_get_child(metatype)
    for source in sources:
        source_table = metadata_group._f_get_child(source)
        for row in tqdm(
            read_chunks(source_table),
            total=source_table.nrows,
            unit="metadata",
            disable=disable_bar,
        ):
            meta_dict = {"source": source}
            try:
                meta_id = int(row["metadata_id"])
            except KeyError:
                meta_id = None
            identifier = (
                row["ID"].decode("utf-8")
                if isinstance(row["ID"], bytes)
                else int(row["ID"])
            )
            # else:
            #     identifier = row["name"].decode()
            if metatype == "families":
                element = pangenome.get_gene_family(identifier)
            elif metatype == "genomes":
                element = pangenome.get_organism(identifier)
            elif metatype == "genes":
                element = pangenome.get_gene(identifier)
            elif metatype == "RGPs":
                element = pangenome.get_region(identifier)
            elif metatype == "spots":
                element = pangenome.get_spot(identifier)
            elif metatype == "modules":
                element = pangenome.get_module(identifier)
            elif metatype == "contigs":
                element = pangenome.get_contig(identifier)
            else:
                expected_types = [
                    "families",
                    "genomes",
                    "contigs",
                    "genes",
                    "RGPs",
                    "spots",
                    "modules",
                ]
                raise KeyError(
                    f"The metatype {metatype} is unexpected. Object associated with metadata are {expected_types}"
                )
            for field in row.dtype.names:
                if field not in ["ID", "name"]:
                    meta_dict[field] = (
                        row[field].decode()
                        if isinstance(row[field], bytes)
                        else row[field]
                    )
            element.add_metadata(metadata=Metadata(**meta_dict), metadata_id=meta_id)
    pangenome.status["metadata"][metatype] = "Loaded"


def read_parameters(h5f: tables.File):
    """
    Read pangenome parameters

    :param h5f: Pangenome HDF5 file
    """
    step_to_parameters = get_pangenome_parameters(h5f)
    print("Parameters:")

    for step, param_name_to_value in step_to_parameters.items():
        print(f"    {step}:")
        for param_name, val in param_name_to_value.items():
            print(f"        {param_name}: {val}")
    print()
    # Cannot use yaml package because some of the parameters are yaml comment
    # yaml_output = yaml.dump({"Parameters":step_to_parameters}, default_flow_style=False, sort_keys=False, indent=4)
    # print(yaml_output)


def get_pangenome_parameters(h5f: tables.File) -> Dict[str, Dict[str, Any]]:
    """
    Read and return the pangenome parameters.

    :param h5f: Pangenome HDF5 file
    :return: A dictionary containing the name of the ppanggolin step as the key, and a dictionary of parameter names
             and their corresponding values used for that step.
    """
    if "/info" in h5f:
        info_group = h5f.root.info
        if "parameters" in info_group._v_attrs._f_list():
            return info_group._v_attrs["parameters"]


def read_pangenome(
    pangenome,
    annotation: bool = False,
    gene_families: bool = False,
    graph: bool = False,
    rgp: bool = False,
    spots: bool = False,
    gene_sequences: bool = False,
    modules: bool = False,
    metadata: bool = False,
    metatypes: Set[str] = None,
    sources: Set[str] = None,
    disable_bar: bool = False,
):
    """
    Reads a previously written pangenome, with all of its parts, depending on what is asked,
    with regard to what is filled in the 'status' field of the hdf5 file.

    :param pangenome: Pangenome object without some information
    :param annotation: get annotation
    :param gene_families: get gene families
    :param graph: get graph
    :param rgp: get RGP
    :param spots: get hotspot
    :param gene_sequences: get gene sequences
    :param modules: get modules
    :param metadata: get metadata
    :param metatypes: metatypes of the metadata to get
    :param sources: sources of the metadata to get (None means all sources)
    :param disable_bar: Allow to disable the progress bar
    """
    if pangenome.file is None:
        raise FileNotFoundError(
            "Your pangenome object has not been associated to any file."
        )
    filename = pangenome.file

    h5f = tables.open_file(filename, "r")

    if (
        annotation
    ):  # I place annotation here, to link gene to gene families if organism are not loaded
        if h5f.root.status._v_attrs.genomesAnnotated:
            logging.getLogger("PPanGGOLiN").info("Reading pangenome annotations...")
            read_annotation(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' has not been annotated, or has been improperly filled"
            )

    if gene_sequences:
        if h5f.root.status._v_attrs.geneSequences:
            logging.getLogger("PPanGGOLiN").info(
                "Reading pangenome gene dna sequences..."
            )
            read_gene_sequences(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have gene sequences, "
                f"or has been improperly filled"
            )

    if gene_families:
        if h5f.root.status._v_attrs.genesClustered:
            logging.getLogger("PPanGGOLiN").info("Reading pangenome gene families...")
            read_gene_families(pangenome, h5f, disable_bar=disable_bar)
            read_gene_families_info(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have gene families, or has been improperly filled"
            )

    if graph:
        if h5f.root.status._v_attrs.NeighborsGraph:
            logging.getLogger("PPanGGOLiN").info("Reading the neighbors graph edges...")
            read_graph(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have graph information, "
                f"or has been improperly filled"
            )

    if rgp:
        if h5f.root.status._v_attrs.predictedRGP:
            logging.getLogger("PPanGGOLiN").info("Reading the RGP...")
            read_rgp(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have RGP information, "
                f"or has been improperly filled"
            )

    if spots:
        if h5f.root.status._v_attrs.spots:
            logging.getLogger("PPanGGOLiN").info("Reading the spots...")
            read_spots(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have spots information, "
                f"or has been improperly filled"
            )
    if modules:
        if h5f.root.status._v_attrs.modules:
            logging.getLogger("PPanGGOLiN").info("Reading the modules...")
            read_modules(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise Exception(
                f"The pangenome in file '{filename}' does not have modules information, "
                f"or has been improperly filled"
            )

    if metadata:
        for metatype in metatypes:

            if h5f.root.status._v_attrs.metadata:
                metastatus = h5f.root.status._f_get_child("metastatus")
                metasources = h5f.root.status._f_get_child("metasources")

                metatype_sources = set(metasources._v_attrs[metatype]) & sources
                if metastatus._v_attrs[metatype] and len(metatype_sources) > 0:
                    logging.getLogger("PPanGGOLiN").info(
                        f"Reading the {metatype} metadata from sources {metatype_sources}..."
                    )
                    read_metadata(
                        pangenome,
                        h5f,
                        metatype,
                        metatype_sources,
                        disable_bar=disable_bar,
                    )
            else:
                raise KeyError(
                    f"The pangenome in file '{filename}' does not have metadata associated to {metatype}, "
                )

    h5f.close()


def get_need_info(
    pangenome,
    need_annotations: bool = False,
    need_families: bool = False,
    need_graph: bool = False,
    need_partitions: bool = False,
    need_rgp: bool = False,
    need_spots: bool = False,
    need_gene_sequences: bool = False,
    need_modules: bool = False,
    need_metadata: bool = False,
    metatypes: Set[str] = None,
    sources: Set[str] = None,
):
    need_info = {
        "annotation": False,
        "gene_families": False,
        "graph": False,
        "rgp": False,
        "spots": False,
        "gene_sequences": False,
        "modules": False,
        "metadata": False,
        "metatypes": metatypes,
        "sources": sources,
    }

    # TODO Automate call if one need another
    if need_annotations:
        if pangenome.status["genomesAnnotated"] == "inFile":
            need_info["annotation"] = True
        elif pangenome.status["genomesAnnotated"] not in ["Computed", "Loaded"]:
            raise Exception(
                "Your pangenome has no genes. See the 'annotate' subcommand."
            )
    if need_families:
        if pangenome.status["genesClustered"] == "inFile":
            need_info["gene_families"] = True
        elif pangenome.status["genesClustered"] not in ["Computed", "Loaded"]:
            raise Exception(
                "Your pangenome has no gene families. See the 'cluster' subcommand."
            )
    if need_graph:
        if pangenome.status["neighborsGraph"] == "inFile":
            need_info["graph"] = True
        elif pangenome.status["neighborsGraph"] not in ["Computed", "Loaded"]:
            raise Exception(
                "Your pangenome does not have a graph (no edges). See the 'graph' subcommand."
            )
    if need_partitions and pangenome.status["partitioned"] not in [
        "Computed",
        "Loaded",
        "inFile",
    ]:
        raise Exception(
            "Your pangenome has not been partitioned. See the 'partition' subcommand"
        )
    if need_rgp:
        if pangenome.status["predictedRGP"] == "inFile":
            need_info["rgp"] = True
        elif pangenome.status["predictedRGP"] not in ["Computed", "Loaded"]:
            raise Exception(
                "Your pangenome  regions of genomic plasticity have not been predicted. See the 'rgp' subcommand"
            )
    if need_spots:
        if pangenome.status["spots"] == "inFile":
            need_info["spots"] = True
        elif pangenome.status["spots"] not in ["Computed", "Loaded"]:
            raise Exception(
                "Your pangenome spots of insertion have not been predicted. See the 'spot' subcommand"
            )
    if need_gene_sequences:
        if pangenome.status["geneSequences"] == "inFile":
            need_info["gene_sequences"] = True
        elif pangenome.status["geneSequences"] not in ["Computed", "Loaded"]:
            raise Exception(
                "Your pangenome does not include gene sequences. "
                "This is possible only if you provided your own cluster file with the 'cluster' subcommand"
            )

    if need_modules:
        if pangenome.status["modules"] == "inFile":
            need_info["modules"] = True
        elif pangenome.status["modules"] not in ["Computed", "Loaded"]:
            raise Exception(
                "Your pangenome modules have not been predicted. See the 'module' subcommand"
            )

    metatypes_to_load = set()
    sources_to_load = set()
    if need_metadata:
        if metatypes is None:
            # load all metadata contained in the pangenome
            metatypes = [
                metatype
                for metatype, status in pangenome.status["metadata"].items()
                if status == "inFile"
            ]
        else:
            # check that specified types have metadata associated
            for metatype in metatypes:
                if pangenome.status["metadata"][metatype] not in [
                    "Computed",
                    "Loaded",
                    "inFile",
                ]:
                    logging.getLogger("PPanGGOLiN").warning(
                        "The pangenome does not have any metadata associated "
                        f"with {metatype}. See the 'metadata' subcommand"
                    )

        if sources is None:
            # load all metadata sources for each metatype
            for metatype in metatypes:
                sources_to_load |= set(pangenome.status["metasources"][metatype])
        else:
            # check that specified source exist for at least one metatype
            for source in set(sources):
                if any(
                    source in pangenome.status["metasources"][metatype]
                    for metatype in metatypes
                ):
                    sources_to_load.add(source)
                else:
                    logging.getLogger("PPanGGOLiN").warning(
                        f"There is no metadata assigned to any element of the pangenome with "
                        f"source={source}. This source is ignored"
                    )

        # select only metatypes that have a requested source .
        for metatype in metatypes:
            if set(pangenome.status["metasources"][metatype]) & sources_to_load:
                metatypes_to_load.add(metatype)
            else:
                logging.getLogger("PPanGGOLiN").debug(
                    f"There is no metadata assigned to {metatype} with specified sources:"
                    f" {', '.join(sources_to_load)} in the pangenome. This metatype is ignored."
                )
        if metatypes_to_load and sources_to_load:
            logging.getLogger("PPanGGOLiN").debug(
                f"metadata types to load: {', '.join(metatypes_to_load)}"
            )
            logging.getLogger("PPanGGOLiN").debug(
                f"metadata sources to load: {', '.join(sources_to_load)}"
            )
            need_info["metadata"] = True
            need_info["metatypes"] = metatypes_to_load
            need_info["sources"] = sources_to_load

    return need_info


def check_pangenome_info(
    pangenome,
    need_annotations: bool = False,
    need_families: bool = False,
    need_graph: bool = False,
    need_partitions: bool = False,
    need_rgp: bool = False,
    need_spots: bool = False,
    need_gene_sequences: bool = False,
    need_modules: bool = False,
    need_metadata: bool = False,
    metatypes: Optional[Set[str]] = None,
    sources: Optional[Set[str]] = None,
    disable_bar: bool = False,
):
    """
    Defines what needs to be read depending on what is needed, and automatically checks if the required elements
    have been computed with regard to the `pangenome.status`

    :param pangenome: Pangenome object without some information
    :param need_annotations: get annotation
    :param need_families: get gene families
    :param need_graph: get graph
    :param need_partitions: get partition
    :param need_rgp: get RGP
    :param need_spots: get hotspot
    :param need_gene_sequences: get gene sequences
    :param need_modules: get modules
    :param need_metadata: get metadata
    :param metatypes: metatypes of the metadata to get (None means all types with metadata)
    :param sources: sources of the metadata to get (None means all possible sources)
    :param disable_bar: Allow to disable the progress bar
    """
    need_info = get_need_info(
        pangenome,
        need_annotations,
        need_families,
        need_graph,
        need_partitions,
        need_rgp,
        need_spots,
        need_gene_sequences,
        need_modules,
        need_metadata,
        metatypes,
        sources,
    )
    if any([v for k, v in need_info.items() if k not in ["metatypes", "sources"]]):
        # if no flag is true, then nothing is needed.
        read_pangenome(pangenome, disable_bar=disable_bar, **need_info)
