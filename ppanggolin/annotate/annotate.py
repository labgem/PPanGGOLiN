#!/usr/bin/env python3

# default libraries
import argparse
import logging
from concurrent.futures import ProcessPoolExecutor

from multiprocessing import get_context
import os
from pathlib import Path
import tempfile
import time
from typing import List, Set, Tuple, Iterable, Dict, Generator, Union
import re
from collections import defaultdict
import warnings

# installed libraries
from tqdm import tqdm
from tables.path import check_name_validity, NaturalNameWarning

# local libraries
from ppanggolin.annotate.synta import (
    annotate_organism,
    get_contigs_from_fasta_file,
    get_dna_sequence,
    init_contig_counter,
    contig_counter,
)
from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Organism, Gene, RNA, Contig
from ppanggolin.utils import (
    read_compressed_or_not,
    mk_file_name,
    detect_filetype,
    check_input_files,
    has_non_ascii,
    replace_non_ascii,
    check_tools_availability,
)
from ppanggolin.formats import write_pangenome
from ppanggolin.metadata import Metadata

# ignore NaturalNameWarning
warnings.filterwarnings("ignore", category=NaturalNameWarning)

ctg_counter = contig_counter


def check_annotate_args(args: argparse.Namespace):
    """Check That the given arguments are usable

    :param args: All arguments provide by user

    :raise Exception:
    """
    if args.fasta is None and args.anno is None:
        raise argparse.ArgumentError(
            argument=None,
            message="You must provide at least a file with the --fasta option to annotate "
            "from sequences, or a file with the --anno option to load annotations from.",
        )

    if hasattr(args, "fasta") and args.fasta is not None:
        check_input_files(args.fasta, True)

    if hasattr(args, "anno") and args.anno is not None:
        check_input_files(args.anno, True)


def create_gene(
    org: Organism,
    contig: Contig,
    gene_counter: int,
    rna_counter: int,
    gene_id: str,
    dbxrefs: Set[str],
    coordinates: List[Tuple[int, int]],
    strand: str,
    gene_type: str,
    position: int = None,
    gene_name: str = "",
    product: str = "",
    genetic_code: int = 11,
    protein_id: str = "",
) -> Gene:
    """
    Create a Gene object and associate to contig and Organism

    :param org: Organism to add gene
    :param contig: Contig to add gene
    :param gene_counter: Gene counter to name gene
    :param rna_counter: RNA counter to name RNA
    :param gene_id: local identifier
    :param dbxrefs: cross-reference to external DB
    :param coordinates: Gene start and stop positions
    :param strand: gene strand association
    :param gene_type: gene type
    :param position: position in contig
    :param gene_name: Gene name
    :param product: Function of gene
    :param genetic_code: Genetic code used
    :param protein_id: Protein identifier
    """

    start, stop = coordinates[0][0], coordinates[-1][1]

    if any(
        dbxref.startswith("MaGe:") or dbxref.startswith("SEED:") for dbxref in dbxrefs
    ):
        if gene_name == "":
            gene_name = gene_id

        for dbxref in dbxrefs:
            if dbxref.startswith("MaGe:"):
                gene_id = dbxref.split(":")[1]
                break
            if dbxref.startswith("SEED:"):
                gene_id = dbxref.split(":")[1]
                break

    if gene_type == "CDS":
        if gene_id == "":
            gene_id = protein_id
            # on rare occasions, there are no 'locus_tag' from downloaded .gbk file.
            # So we use the protein_id field instead. (which is not supposed to be unique,
            # but was when cases like this were encountered)

        new_gene = Gene(org.name + "_CDS_" + str(gene_counter).zfill(4))
        new_gene.fill_annotations(
            start=start,
            stop=stop,
            strand=strand,
            coordinates=coordinates,
            gene_type=gene_type,
            name=gene_name,
            position=position,
            product=product,
            local_identifier=gene_id,
            genetic_code=genetic_code,
        )
        contig.add(new_gene)
    else:  # if not CDS, it is RNA
        new_gene = RNA(org.name + f"_{gene_type}_" + str(rna_counter).zfill(4))
        new_gene.fill_annotations(
            start=start,
            stop=stop,
            strand=strand,
            coordinates=coordinates,
            gene_type=gene_type,
            name=gene_name,
            product=product,
        )
        contig.add_rna(new_gene)
    new_gene.fill_parents(org, contig)
    return new_gene


def extract_positions(string: str) -> Tuple[List[Tuple[int, int]], bool, bool, bool]:
    """
    Extracts start and stop positions from a string and determines whether it is complement and pseudogene.

    Example of strings that the function is able to process:

    "join(190..7695,7695..12071)",
    "complement(join(4359800..4360707,4360707..4360962))",
    "join(6835405..6835731,1..1218)",
    "join(1375484..1375555,1375557..1376579)",
    "complement(6815492..6816265)",
    "6811501..6812109",
    "complement(6792573..>6795461)",
    "join(1038313,1..1016)"


    :param string: The input string containing position information.

    :return: A tuple containing a list of tuples representing start and stop positions,
             a boolean indicating whether it is complement,
             a boolean indicating whether it is a partial gene at start position and
             a boolean indicating whether it is a partial gene at end position.

    :raises ValueError: If the string is not formatted as expected or if positions cannot be parsed as integers.
    """
    complement = False
    coordinates = []
    has_partial_start = False
    has_partial_end = False

    # Check if 'complement' exists in the string
    if "complement" in string:
        complement = True

    if "(" in string:
        # Extract positions found inside the parenthesis
        inner_parentheses_regex = r"\(([^()]+)\)"
        inner_matches = re.findall(inner_parentheses_regex, string)

        try:
            positions = inner_matches[-1]
        except IndexError:
            raise ValueError(f"Gene position {string} is not formatted as expected.")
    else:
        positions = string.rstrip()

    # Check if '>' or '<' exists in the positions to identify partial genes

    if ">" in positions or "<" in positions:
        if "<" in positions.split(",")[0]:
            has_partial_start = True

        if ">" in positions.split(",")[-1]:
            has_partial_end = True

        inner_positions = ",".join(positions.split(",")[1:-1])

        if (
            ">" in inner_positions
            or "<" in inner_positions
            or (not has_partial_end and not has_partial_start)
        ):
            raise ValueError(
                f"Error parsing positions '{positions}' extracted from GBFF string '{string}'. "
                f"Chevrons are found in the inner position. This case is unexpected and not handle."
            )

    for position in positions.split(","):

        try:
            start, stop = position.replace(">", "").replace("<", "").split("..")
        except ValueError:
            # in some case there is only one position meaning that the gene is long of only one nt in this piece.
            # for instance : join(1038313,1..1016)
            start = position.replace(">", "").replace("<", "")
            stop = start
        try:
            start, stop = int(start), int(stop)
        except ValueError:
            raise ValueError(
                f"Error parsing position '{position}' extracted from GBFF string '{string}'. "
                f"Start position ({start}) and/or stop position ({stop}) are not valid integers."
            )

        coordinates.append((start, stop))

    return coordinates, complement, has_partial_start, has_partial_end


def parse_gbff_by_contig(
    gbff_file_path: Path,
) -> Generator[
    Tuple[Dict[str, str], Generator[Dict[str, Union[str, Set[str]]], None, None], str],
    None,
    None,
]:
    """
    Parse a GBFF file by contig and yield tuples containing header, feature, and sequence info for each contig.

    :param gbff_file_path: Path to the GBFF file.
    :return: A generator that yields tuples containing header lines, feature lines, and sequence info for each contig.
    """
    header_lines = []
    feature_lines = []
    sequence_lines = []

    current_section = None

    with read_compressed_or_not(gbff_file_path) as fl:
        for i, line in enumerate(fl):
            # Skip blank lines
            if not line.strip():
                continue

            if line.startswith("LOCUS") or line.startswith("CONTIG"):
                # CONTIG line are found between FEATURES and ORIGIN and are put in header section here for simplicity
                current_section = "header"

            elif line.startswith("FEATURES"):
                current_section = "feature"
                continue

            elif line.startswith("ORIGIN"):
                current_section = "sequence"
                continue

            if line.strip() == "//":
                # Check that each section has some lines
                assert header_lines and feature_lines and sequence_lines, (
                    "Missing section in GBFF file. "
                    f"Contig ending at line {i + 1} has an empty section. It has "
                    f"{len(header_lines)} header lines, "
                    f"{len(header_lines)} feature lines, "
                    f"and {len(sequence_lines)} sequence lines."
                )
                yield (
                    parse_contig_header_lines(header_lines),
                    parse_feature_lines(feature_lines),
                    parse_dna_seq_lines(sequence_lines),
                )

                header_lines = []
                feature_lines = []
                sequence_lines = []
                current_section = None
                continue

            if current_section == "header":
                header_lines.append(line)

            elif current_section == "feature":
                feature_lines.append(line)

            elif current_section == "sequence":
                sequence_lines.append(line)

            else:
                raise ValueError(
                    f"Unexpected structure in GBFF file: {gbff_file_path}. {line}"
                )

    # In case the last // is missing, return the last contig
    if header_lines or feature_lines or sequence_lines:
        yield (
            parse_contig_header_lines(header_lines),
            parse_feature_lines(feature_lines),
            parse_dna_seq_lines(sequence_lines),
        )


def parse_contig_header_lines(header_lines: List[str]) -> Dict[str, str]:
    """
    Parse required information from header lines of a contig from a GBFF file.

    :param header_lines: List of strings representing header lines of a contig from a GBFF file.
    :return: A dict with keys representing different fields and values representing their corresponding values joined by new line.
    """
    field = ""  # Initialize field
    field_to_value = defaultdict(
        list
    )  # Initialize defaultdict to store field-value pairs

    for line in header_lines:
        field_of_line = line[:12].strip()  # Extract field from the first 12 characters

        if len(field_of_line) > 1 and field_of_line.isupper():
            field = field_of_line  # Update current field

        # Append value to the current field in the defaultdict
        field_to_value[field].append(line[12:].strip())

    return {field: "\n".join(value) for field, value in field_to_value.items()}


def parse_feature_lines(
    feature_lines: List[str],
) -> Generator[Dict[str, Union[str, Set[str]]], None, None]:
    """
    Parse feature lines from a GBFF file and yield dictionaries representing each feature.

    :param feature_lines: List of strings representing feature lines from a GBFF file.
    :return: A generator that yields dictionaries, each representing a feature with its type, location, and qualifiers.
    """

    def stringify_feature_values(
        feature: Dict[str, List[str]]
    ) -> Dict[str, Union[str, Set[str]]]:
        """
        All value of the returned dict are str except for db_xref that is a list.
        When multiple values exist for the same tag only the first one is kept.
        """
        stringify_feature = {}
        for tag, val in feature.items():
            if tag == "db_xref":
                stringify_feature[tag] = set(val)
            elif isinstance(val, list):
                stringify_feature[tag] = val[0]
            else:
                stringify_feature[tag] = val
        return defaultdict(str, stringify_feature)

    current_feature = {}
    current_qualifier = None

    for line in feature_lines:
        # Check if the line starts a new feature
        if len(line[:21].strip()) > 0:
            if current_feature:
                # yield last feature
                yield stringify_feature_values(current_feature)

            current_feature = {
                "feature_type": line[:21].strip(),
                "location": [line[21:].strip()],
            }
            current_qualifier = "location"

        elif line.strip().startswith("/"):
            qualifier_line = line.strip()[1:]  # [1:] used to remove /

            if "=" in qualifier_line:
                current_qualifier, value = qualifier_line.split("=", 1)
            else:
                current_qualifier, value = qualifier_line, qualifier_line
            # clean value from quote
            value = value[1:] if value.startswith('"') else value
            value = value[:-1] if value.endswith('"') else value

            if current_qualifier in current_feature:
                current_feature[current_qualifier].append(value)

            else:
                current_feature[current_qualifier] = [value]

        else:
            # the line does not start a qualifier so it's the continuation of the last qualifier value.
            value = line.strip()
            value = value[:-1] if value.endswith('"') else value
            current_feature[current_qualifier][-1] += f" {value}"

    # Append the last feature
    if current_feature:
        yield stringify_feature_values(current_feature)


def parse_dna_seq_lines(sequence_lines: List[str]) -> str:
    """
    Parse sequence_lines from a GBFF file and return dna sequence

    :param sequence_lines: List of strings representing sequence lines from a GBFF file.
    :return: a string in upper case of the DNA sequences that have been cleaned
    """
    sequence = ""
    for line in sequence_lines:
        sequence += line[10:].replace(" ", "").strip().upper()
    return sequence


def combine_contigs_metadata(
    contig_to_metadata: Dict[Contig, Dict[str, str]]
) -> Tuple[Dict[str, str], Dict[Contig, Dict[str, str]]]:
    """
    Combine contig metadata to identify shared and unique metadata tags and values.

    :param contig_to_metadata: A dictionary mapping each contig to its associated metadata.
    :return: A tuple containing:
        - A dictionary of shared metadata tags and values present in all contigs.
        - A dictionary mapping each contig to its unique metadata tags and values.
    """
    # Flatten all metadata items and count their occurrences
    all_tag_to_value = [
        (tag, value)
        for source_info in contig_to_metadata.values()
        for tag, value in source_info.items()
        if isinstance(value, str)
    ]

    # Filter tags that would have a / as it is forbidden when writing the table in HDF5. Such tag can appear with db_xref formatting
    invalid_tag_names = []
    for tag, value in set(all_tag_to_value):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                check_name_validity(tag)
        except ValueError as err:
            logging.getLogger("PPanGGOLiN").debug(
                f"{err}. The tag {tag} is ignored for metadata."
            )
            invalid_tag_names.append(tag)

        if value == "":
            logging.getLogger("PPanGGOLiN").debug(
                f"Ignoring tag '{tag}' for metadata due to an empty value."
            )
            invalid_tag_names.append(tag)

    all_tag_to_value = [
        (tag, value) for tag, value in all_tag_to_value if tag not in invalid_tag_names
    ]

    contig_count = len(contig_to_metadata)

    # Identify tags and values shared by all contigs
    shared_tag_and_values = {
        tag_and_value
        for tag_and_value in all_tag_to_value
        if all_tag_to_value.count(tag_and_value) == contig_count
    }

    # Create a dictionary for shared metadata
    genome_metadata = dict(shared_tag_and_values)

    contig_to_uniq_metadata = {}
    for contig, tag_to_value in contig_to_metadata.items():
        # Identify unique metadata for each contig
        uniq_tag_to_value = {
            tag: value
            for tag, value in tag_to_value.items()
            if tag not in list(genome_metadata) + invalid_tag_names
            and isinstance(value, str)
        }
        if uniq_tag_to_value:
            contig_to_uniq_metadata[contig] = uniq_tag_to_value

    return genome_metadata, contig_to_uniq_metadata


def reverse_complement_coordinates(
    coordinates: List[Tuple[int, int]]
) -> List[Tuple[int, int]]:
    """
    Reverses and inverts the given list of coordinates. Each coordinate pair (start, end) is transformed into
    (-end, -start) and the order of the coordinates is reversed.

    :param coordinates: A list of tuples representing the coordinates to be reversed and inverted.
    :return: A list of reversed and inverted coordinates.
    """
    # Reverse the order of coordinates and invert each start and end
    return [(-end, -start) for start, end in coordinates[::-1]]


def shift_end_coordinates(
    coordinates: List[Tuple[int, int]], shift: int
) -> List[Tuple[int, int]]:
    """
    Shifts the end of a set of coordinates by a specified amount and then returns the final shifted coordinates.
    This involves reversing the coordinates twice, shifting the start, and then returning the original orientation.

    :param coordinates: A list of tuples representing the original coordinates.
    :param shift: The amount by which the end coordinate should be shifted.
    :return: The coordinates after the end shift and reverse complement transformations.
    """
    # Perform reverse complement, shift the start (corresponding to the original end), and reverse complement back
    rc_coordinates = reverse_complement_coordinates(coordinates)
    rc_coordinates_shifted = shift_start_coordinates(rc_coordinates, shift)
    final_coordinates = reverse_complement_coordinates(rc_coordinates_shifted)

    return final_coordinates


def shift_start_coordinates(
    coordinates: List[Tuple[int, int]], shift: int
) -> List[Tuple[int, int]]:
    """
    Shifts the start of the first coordinate in the list by the specified amount. If the shift results in
    a negative or zero-length interval for the first coordinate, this interval is removed, and the shift is
    propagated to the next coordinate if necessary.

    :param coordinates: A list of tuples representing the coordinates.
    :param shift: The amount by which the start coordinate should be shifted.
    :return: A new list of coordinates with the shifted start.
    """
    new_coordinates = coordinates.copy()

    # Shift the start of the first coordinate
    adjusted_start = new_coordinates[0][0] + shift

    adjusted_coordinate_part = (adjusted_start, new_coordinates[0][1])
    new_coordinates[0] = adjusted_coordinate_part

    # Calculate length of the adjusted part
    adjusted_part_length = adjusted_coordinate_part[1] - adjusted_coordinate_part[0] + 1

    if adjusted_part_length <= 0:
        # If the shift results in a zero or negative length, handle accordingly
        if len(new_coordinates) <= 1:
            raise ValueError(
                f"Shifting the start resulted in a gene with null or negative size. "
                f"Coordinates: {coordinates}, Shift: {shift}"
            )
        else:
            logging.getLogger("PPanGGOLiN").warning(
                f"Coordinate part {new_coordinates[0]} resulted in a non-positive length after shift. Removing it."
            )
            new_coordinates = new_coordinates[1:]

            # If length is negative, propagate the shift to the next coordinate
            if adjusted_part_length < 0:
                new_shift = -adjusted_part_length
                new_coordinates = shift_start_coordinates(new_coordinates, new_shift)

    return new_coordinates


def fix_partial_gene_coordinates(
    coordinates: List[Tuple[int, int]],
    is_complement: bool,
    start_shift: int,
    ensure_codon_multiple: bool = True,
) -> List[Tuple[int, int]]:
    """
    Adjusts gene coordinates if they have partial starts or ends, ensuring the gene length is a multiple of 3.

    If the gene is on the complement strand, the adjustments will be reversed (i.e., applied to the opposite ends).

    :param coordinates: List of coordinate tuples (start, stop) for the gene.
    :param is_complement: Flag indicating if the gene is on the complement strand.
    :param start_shift: The value by which the start coordinate should be shifted.
    :param ensure_codon_multiple: Flag to check that gene length is a multiple of 3.

    :return: A new list of adjusted coordinate tuples.
    """

    if not coordinates:
        raise ValueError(
            "No coordinates provided. Cannot fix partial gene coordinates."
        )

    # Non-complement strand adjustments
    if not is_complement:
        if start_shift != 0:
            coordinates = shift_start_coordinates(coordinates, start_shift)

        if ensure_codon_multiple:
            # Ensure the gene length is a multiple of 3 by adjusting the last end
            gene_length = sum([(stop - start + 1) for start, stop in coordinates])
            end_shift = gene_length % 3
            if end_shift != 0:
                coordinates = shift_end_coordinates(coordinates, end_shift)

    # Complement strand adjustments
    else:
        if start_shift != 0:
            coordinates = shift_end_coordinates(coordinates, start_shift)

        if ensure_codon_multiple:
            # Adjust first start for complement strand
            gene_length = sum([(stop - start + 1) for start, stop in coordinates])
            start_shift = gene_length % 3
            if start_shift != 0:
                coordinates = shift_start_coordinates(coordinates, start_shift)

    # Final length validation
    gene_length = sum([(stop - start + 1) for start, stop in coordinates])

    if gene_length % 3 != 0:
        logging.getLogger("PPanGGOLiN").warning(
            f"Gene with coordinates: {coordinates} has a length that is not a multiple of 3 after adjusting for partiality with new cordinates ({coordinates}): {gene_length}"
        )

    return coordinates


def read_org_gbff(
    organism_name: str,
    gbff_file_path: Path,
    circular_contigs: List[str],
    use_pseudogenes: bool = False,
    translation_table: int = 11,
) -> Tuple[Organism, bool]:
    """
    Read a GBFF file and fills Organism, Contig and Genes objects based on information contained in this file

    :param organism_name: Organism name
    :param gbff_file_path: Path to corresponding GBFF file
    :param circular_contigs: list of contigs
    :param use_pseudogenes: Allow to read pseudogenes
    :param translation_table: Translation table (genetic code) to use when /transl_table is missing from CDS tags.


    :return: Organism complete and true for sequence in file
    """
    used_transl_table_arg = 0
    global ctg_counter

    organism = Organism(organism_name)

    logging.getLogger("PPanGGOLiN").debug(
        f"Extracting genes information from the given gbff {gbff_file_path.name}"
    )
    gene_counter = 0
    rna_counter = 0
    contig_to_metadata = {}

    for header, features, sequence in parse_gbff_by_contig(gbff_file_path):
        if "LOCUS" not in header:
            raise ValueError("Missing LOCUS line in GBFF header.")

        if "VERSION" in header and header["VERSION"] != "":
            contig_id = header["VERSION"]
        else:
            # If contig_id is not specified in VERSION field like with Prokka, in that case we use the one in LOCUS
            contig_id = header["LOCUS"].split()[0]

        contig_len = int(header["LOCUS"].split()[1])

        if contig_len != len(sequence):
            logging.getLogger("PPanGGOLiN").warning(
                "Unable to determine if the contig is circular or linear in file "
                f"'{gbff_file_path}' from the LOCUS header information: {header['LOCUS']}. "
                "By default, the contig will be considered linear."
            )

        if "CIRCULAR" in header["LOCUS"].upper():
            # this line contains linear/circular word telling if the dna sequence is circularized or not
            is_circ = True

        elif "LINEAR" in header["LOCUS"].upper():
            is_circ = False

        else:
            is_circ = False
            logging.getLogger("PPanGGOLiN").warning(
                f"It's impossible to identify if contig {contig_id} is circular or linear."
                f"in file {gbff_file_path}."
            )

        try:
            contig = organism.get(contig_id)
        except KeyError:
            with contig_counter.get_lock():
                contig = Contig(
                    contig_counter.value,
                    contig_id,
                    True if contig_id in circular_contigs or is_circ else False,
                )
                contig_counter.value += 1
            organism.add(contig)
            contig.length = contig_len

        for feature in features:
            if feature["feature_type"] == "source":
                contig_to_metadata[contig] = {
                    tag: value
                    for tag, value in feature.items()
                    if tag not in ["feature_type", "location"]
                    and isinstance(value, str)
                }
                if "db_xref" in feature:
                    try:
                        db_xref_for_metadata = {
                            f"db_xref_{database}": identifier
                            for database_identifier in feature["db_xref"]
                            for database, identifier in [database_identifier.split(":")]
                        }
                        contig_to_metadata[contig].update(db_xref_for_metadata)
                    except ValueError:
                        logging.getLogger("PPanGGOLiN").warning(
                            f"db_xref values does not have the expected format. Expect 'db_xref=<database>:<identifier>' "
                            f"but got {feature['db_xref']} in file {gbff_file_path}. "
                            "db_xref tags is therefore not retrieved in contig/genomes metadata."
                        )
                    else:
                        contig_to_metadata[contig].update(db_xref_for_metadata)
            genetic_code = ""
            if feature["feature_type"] not in ["CDS", "rRNA", "tRNA"]:
                continue
            coordinates, is_complement, has_partial_start, has_partial_end = (
                extract_positions("".join(feature["location"]))
            )

            if "pseudo" in feature and not use_pseudogenes:
                continue

            elif "transl_except" in feature and not use_pseudogenes:
                # that's probably a 'stop' codon into selenocystein.
                logging.getLogger("PPanGGOLiN").info(
                    f"CDS '{feature['locus_tag']}' contains a 'transl_except' annotation ({feature['transl_except']}) "
                    f"in contig '{contig}' in file '{gbff_file_path}'. "
                    f"PPanGGOLiN does not handle 'transl_except' annotations. This gene's protein sequence "
                    "will likely contain an internal stop codon when translated with PPanGGOLiN."
                )

            for field in ["product", "gene", "db_xref"]:

                if field in feature and has_non_ascii(feature[field]):

                    logging.getLogger("PPanGGOLiN").warning(
                        f"In genome '{organism}', the '{field}' field of gene '{feature['locus_tag']}' contains non-ASCII characters: '{feature[field]}'. "
                        "These characters cannot be stored in the HDF5 file and will be replaced by underscores."
                    )
                    feature[field] = replace_non_ascii(feature[field])

            if feature["feature_type"] == "CDS":
                if feature["transl_table"] == "":
                    used_transl_table_arg += 1
                    genetic_code = translation_table
                else:
                    genetic_code = int(feature["transl_table"])

                if has_partial_start or has_partial_end:
                    start_shift = (
                        0
                        if "codon_start" not in feature
                        else int(feature["codon_start"]) - 1
                    )  # -1 is to be in zero based index.

                    coordinates = fix_partial_gene_coordinates(
                        coordinates,
                        is_complement=is_complement,
                        start_shift=start_shift,
                    )

            strand = "-" if is_complement else "+"

            gene = create_gene(
                org=organism,
                contig=contig,
                gene_counter=gene_counter,
                rna_counter=rna_counter,
                gene_id=feature["locus_tag"],
                dbxrefs=feature["db_xref"],
                coordinates=coordinates,
                strand=strand,
                gene_type=feature["feature_type"],
                position=contig.number_of_genes,
                gene_name=feature["gene"],
                product=feature["product"],
                genetic_code=genetic_code,
                protein_id=feature["protein_id"],
            )

            gene.add_sequence(get_dna_sequence(sequence, gene))

            if feature["feature_type"] == "CDS":
                gene_counter += 1
            else:
                rna_counter += 1

    genome_metadata, contig_to_uniq_metadata = combine_contigs_metadata(
        contig_to_metadata
    )
    organism.add_metadata(
        metadata=Metadata(source="annotation_file", **genome_metadata)
    )

    for contig, metadata_dict in contig_to_uniq_metadata.items():
        contig.add_metadata(Metadata(source="annotation_file", **metadata_dict))

    if used_transl_table_arg:
        logging.getLogger("PPanGGOLiN").debug(
            f"transl_table tag was not found for {used_transl_table_arg} CDS "
            f"in {gbff_file_path}. Provided translation_table argument value was used instead: {translation_table}."
        )

    return organism, True


def parse_db_xref_metadata(
    db_xref_values: List[str], annot_file_path: Path = ""
) -> Dict[str, str]:
    """
    Parses a list of db_xref values and returns a dictionary with formatted keys and identifiers.

    :param db_xref_values: List of db_xref strings in the format <database>:<identifier>.
    :param annot_file_path: Path to the annotation file being processed.
    :return: Dictionary with keys formatted as 'db_xref_<database>' and their corresponding identifiers.
    """
    db_xref_for_metadata = {}
    try:
        # Create a dictionary with keys formatted as 'db_xref_<database>' and their corresponding identifiers
        db_xref_for_metadata = {
            f"db_xref_{database}": identifier
            for database_identifier in db_xref_values
            for database, identifier in [database_identifier.split(":")]
        }
    except ValueError:
        logging.getLogger("PPanGGOLiN").warning(
            f"db_xref values do not have the expected format. Expected '<database>:<identifier>', "
            f"but got {db_xref_values} in file {annot_file_path}. "
            "db_xref tags are therefore not retrieved in contig/genome metadata."
        )
    return db_xref_for_metadata


def read_org_gff(
    organism: str,
    gff_file_path: Path,
    circular_contigs: List[str],
    pseudo: bool = False,
    translation_table: int = 11,
) -> Tuple[Organism, bool]:
    """
    Read annotation from GFF file

    :param organism: Organism name
    :param gff_file_path: Path corresponding to GFF file
    :param circular_contigs: List of circular contigs
    :param pseudo: Allow to read pseudogene
    :param translation_table: Translation table (genetic code) to use when transl_table is missing from CDS tags.


    :return: Organism object and if there are sequences associated or not
    """
    # TODO: This function would need some refactoring.
    used_transl_table_arg = 0
    global ctg_counter

    (
        gff_seqname,
        _,
        gff_type,
        gff_start,
        gff_end,
        _,
        gff_strand,
        frame,
        gff_attribute,
    ) = range(0, 9)

    # Missing values: source, score. They are unused.
    def get_gff_attributes(gff_fields: list) -> dict:
        """Parses the gff attribute's line and outputs the attributes_get in a dict structure.

        :param gff_fields: A gff line stored as a list. Each element of the list is a column of the gff.

        :return: Attributes get
        """
        attributes_field = [
            f for f in gff_fields[gff_attribute].strip().split(";") if len(f) > 0
        ]
        attributes_get = {}
        for att in attributes_field:
            try:
                (key, value) = att.strip().split("=")
                attributes_get[key.upper()] = value
            except ValueError:
                pass  # we assume that it is a strange, but useless field for our analysis
        return attributes_get

    def get_id_attribute(attributes_dict: dict) -> str:
        """
        Gets the ID of the element from which the provided attributes_get were extracted.
        Raises an error if no ID is found.

        :param attributes_dict: Attributes from one gff line

        :return: CDS identifier
        """
        element_id = attributes_dict.get("ID")
        if not element_id:
            raise Exception(
                f"Each CDS type of the gff files must own a unique ID attribute. "
                f"Not the case for file: {gff_file_path}"
            )
        return element_id

    def check_chevrons_in_start_and_stop(
        start: str, stop: str
    ) -> Tuple[int, int, bool]:
        """
        Checks for the presence of chevrons ('<' or '>') in the start and stop strings, removes them if present,
        and converts the remaining parts to integers.

        :param start: The start string which may contain chevrons.
        :param stop: The stop string which may contain chevrons.

        :return: A tuple containing the integer values of start and stop, and a boolean indicating if chevrons were present in either string.
        """
        chevrons_present = ">" in start or "<" in start or ">" in stop or "<" in stop

        if chevrons_present:
            start = int(start.replace("<", "").replace(">", ""))
            stop = int(stop.replace("<", "").replace(">", ""))
        else:
            start = int(start)
            stop = int(stop)

        return start, stop, chevrons_present

    contig = None  # initialize contig
    has_fasta = False
    fasta_string = ""
    org = Organism(organism)
    gene_counter = 0
    rna_counter = 0
    attr_prodigal = None
    contig_name_to_region_info = {}

    id_attr_to_gene_id = {}

    with read_compressed_or_not(gff_file_path) as gff_file:
        for line in gff_file:
            if has_fasta:
                fasta_string += line

            elif line.startswith("##", 0, 2):
                if line.startswith("FASTA", 2, 7):
                    has_fasta = True
                elif line.startswith("sequence-region", 2, 17):
                    fields = [el.strip() for el in line.split()]
                    if len(fields) != 4:
                        raise Exception(
                            "Pragma '##sequence-region' has an unexpected format. "
                            f"Expecting the format '##sequence-region seqid start stop', and got the following: '{line.strip()}'"
                        )
                    with contig_counter.get_lock():
                        contig = Contig(
                            contig_counter.value,
                            fields[1],
                            True if fields[1] in circular_contigs else False,
                        )
                        contig_counter.value += 1
                    org.add(contig)
                    contig.length = int(fields[-1]) - int(fields[2]) + 1
                else:
                    continue

            elif line.startswith("#"):
                if line.startswith("Sequence Data", 2, 15):  # GFF from prodigal
                    fields_prodigal = [
                        el.strip() for el in line.split(": ")[1].split(";")
                    ]
                    attr_prodigal = {
                        field.split("=")[0]: field.split("=")[1]
                        for field in fields_prodigal
                    }
                else:  # comment lines to be ignores by parsers
                    continue

            elif (
                line.rstrip() == ""
            ):  # empty lines are not expected, but they do not carry information, so we'll ignore them
                continue

            else:
                fields_gff = [el.strip() for el in line.split("\t")]
                attributes = get_gff_attributes(fields_gff)

                pseudogene = False

                gene_start, gene_stop, has_chevron = check_chevrons_in_start_and_stop(
                    start=fields_gff[gff_start], stop=fields_gff[gff_end]
                )

                for field in ["PRODUCT", "NAME", "DB_XREF", "DBXREF"]:
                    if field in attributes and has_non_ascii(attributes[field]):
                        logging.getLogger("PPanGGOLiN").warning(
                            f"In genome '{organism}', the '{field}' field of feature '{attributes['locus_tag']}' contains non-ASCII characters: '{attributes[field]}'. "
                            "These characters cannot be stored in the HDF5 file and will be replaced by underscores."
                        )
                        attributes[field] = replace_non_ascii(attributes[field])

                if fields_gff[gff_type] == "region":
                    # keep region attributes to add them as metadata of genome and contigs
                    # excluding some info as they are already contained in contig object.

                    contig_name_to_region_info[fields_gff[gff_seqname]] = {
                        tag.lower(): value
                        for tag, value in attributes.items()
                        if tag not in ["ID", "NAME", "IS_CIRCULAR", "DB_XREF", "DBXREF"]
                    }

                    if (
                        "DB_XREF" in attributes or "DBXREF" in attributes
                    ):  # db_xref can be written Dbxref and db_ref
                        dbxref_tag = "DB_XREF" if "DB_XREF" in attributes else "DBXREF"
                        dbxref_metadata = parse_db_xref_metadata(
                            attributes[dbxref_tag].split(","), gff_file_path
                        )
                        contig_name_to_region_info[fields_gff[gff_seqname]].update(
                            dbxref_metadata
                        )

                    if (
                        "IS_CIRCULAR" in attributes
                        and attributes["IS_CIRCULAR"] == "true"
                    ):
                        contig_name = fields_gff[gff_seqname]

                        if contig is not None:
                            logging.getLogger("PPanGGOLiN").debug(
                                f"Contig {contig.name} is circular."
                            )
                            contig.is_circular = True
                            assert contig.name == contig_name
                        else:
                            # contig object has not been initialized yet.
                            # let's keep the circularity info in the circular_contigs list
                            circular_contigs.append(contig_name)

                elif fields_gff[gff_type] == "CDS" or "RNA" in fields_gff[gff_type]:

                    id_attribute = get_id_attribute(attributes)
                    locus_tag = attributes.get("LOCUS_TAG")
                    protein_id = attributes.get("PROTEIN_ID")

                    if locus_tag is not None:
                        gene_id = locus_tag

                    elif protein_id is not None:
                        gene_id = protein_id

                    else:
                        gene_id = id_attribute

                    name = attributes.pop("NAME", attributes.pop("GENE", ""))

                    if "PSEUDO" in attributes or "PSEUDOGENE" in attributes:
                        pseudogene = True

                    if (
                        "PARTIAL" in attributes
                        and attributes["PARTIAL"].upper() == "TRUE"
                    ) or has_chevron:
                        is_partial = True
                    else:
                        is_partial = False

                    product = attributes.pop("PRODUCT", "")

                    if contig is None or contig.name != fields_gff[gff_seqname]:
                        # get the current contig
                        try:
                            contig = org.get(fields_gff[gff_seqname])
                        except KeyError:
                            with contig_counter.get_lock():
                                contig = Contig(
                                    contig_counter.value,
                                    fields_gff[gff_seqname],
                                    (
                                        True
                                        if fields_gff[gff_seqname] in circular_contigs
                                        else False
                                    ),
                                )
                                contig_counter.value += 1
                            org.add(contig)
                            if attr_prodigal is not None:
                                contig.length = int(attr_prodigal["seqlen"])

                    if fields_gff[gff_type] == "CDS" and (
                        not pseudogene or (pseudogene and pseudo)
                    ):

                        if "TRANSL_TABLE" in attributes:
                            genetic_code = int(attributes["TRANSL_TABLE"])
                        else:
                            used_transl_table_arg += 1
                            genetic_code = translation_table

                        gene_frame = 0
                        #  Get value of frame if valid
                        if fields_gff[frame] in ["1", "2", "0"]:
                            gene_frame = int(fields_gff[frame])

                        if (
                            id_attribute in id_attr_to_gene_id
                        ):  # the ID has already been seen at least once in this genome
                            existing_gene = id_attr_to_gene_id[id_attribute]
                            new_gene_info = {
                                "strand": fields_gff[gff_strand],
                                "type": fields_gff[gff_type],
                                "name": name,
                                "position": contig.number_of_genes,
                                "product": product,
                                "local_identifier": gene_id,
                                "start": gene_start,
                                "stop": gene_stop,
                                "frame": gene_frame,
                            }

                            check_and_add_extra_gene_part(existing_gene, new_gene_info)

                            continue

                        gene = Gene(org.name + "_CDS_" + str(gene_counter).zfill(4))

                        id_attr_to_gene_id[id_attribute] = gene

                        # here contig is filled in order, so position is the number of genes already stored in the contig.
                        gene.fill_annotations(
                            start=gene_start,
                            stop=gene_stop,
                            strand=fields_gff[gff_strand],
                            gene_type=fields_gff[gff_type],
                            name=name,
                            product=product,
                            position=contig.number_of_genes,
                            local_identifier=gene_id,
                            genetic_code=genetic_code,
                            is_partial=is_partial,
                            frame=gene_frame,
                        )

                        gene.fill_parents(org, contig)
                        gene_counter += 1
                        contig.add(gene)

                    elif "RNA" in fields_gff[gff_type]:

                        rna_type = fields_gff[gff_type]
                        rna = RNA(
                            org.name + f"_{rna_type}_" + str(rna_counter).zfill(4)
                        )

                        rna.fill_annotations(
                            start=gene_start,
                            stop=gene_stop,
                            strand=fields_gff[gff_strand],
                            gene_type=fields_gff[gff_type],
                            name=name,
                            product=product,
                            local_identifier=gene_id,
                        )
                        rna.fill_parents(org, contig)
                        rna_counter += 1
                        contig.add_rna(rna)

    # Fix partial genes coordinates
    for contig in org.contigs:
        for gene in contig.genes:
            if gene.is_partial:
                is_complement = gene.strand == "-"
                gene.coordinates = fix_partial_gene_coordinates(
                    gene.coordinates,
                    is_complement=is_complement,
                    start_shift=gene.frame,
                )

    # GET THE FASTA SEQUENCES OF THE GENES
    if fasta_string == "":
        has_fasta = False

    if has_fasta:
        contig_sequences = get_contigs_from_fasta_file(org, fasta_string.split("\n"))

        correct_putative_overlaps(org.contigs)

        for contig in org.contigs:

            for gene in contig.genes:
                gene.add_sequence(get_dna_sequence(contig_sequences[contig.name], gene))
            for rna in contig.RNAs:
                rna.add_sequence(get_dna_sequence(contig_sequences[contig.name], rna))

    # add metadata to genome and contigs
    if contig_name_to_region_info:
        add_metadata_from_gff_file(contig_name_to_region_info, org, gff_file_path)

    if used_transl_table_arg:
        logging.getLogger("PPanGGOLiN").debug(
            f"transl_table tag was not found for {used_transl_table_arg} CDS "
            f"in {gff_file_path}. Provided translation_table argument value was used instead: {translation_table}."
        )
    return org, has_fasta


def add_metadata_from_gff_file(
    contig_name_to_region_info: Dict[str, Dict[str, str]],
    org: Organism,
    gff_file_path: Path,
):
    """
    Add metadata to the organism object from a GFF file.

    :param contig_name_to_region_info: A dictionary mapping contig names to their corresponding region information.
    :param org: The organism object to which metadata will be added.
    :param gff_file_path: The path to the GFF file.
    """

    # Check if the number of contigs matches the expected number in the organism
    if len(contig_name_to_region_info) == org.number_of_contigs:
        contig_to_region_info = {
            org.get(name): region_info
            for name, region_info in contig_name_to_region_info.items()
        }
        genome_metadata, contig_to_uniq_metadata = combine_contigs_metadata(
            contig_to_region_info
        )
        if genome_metadata:
            org.add_metadata(Metadata(source="annotation_file", **genome_metadata))
    else:
        logging.getLogger("PPanGGOLiN").warning(
            f"Inconsistent data in GFF file {gff_file_path}: "
            f"expected {org.number_of_contigs} contigs but found {len(contig_name_to_region_info)} regions."
        )
        contig_to_uniq_metadata = {}

    for contig, metadata_dict in contig_to_uniq_metadata.items():
        contig.add_metadata(Metadata(source="annotation_file", **metadata_dict))


def check_and_add_extra_gene_part(
    gene: Gene, new_gene_info: Dict, max_separation: int = 10
):
    """
    Checks and potentially adds extra gene parts based on new gene information.
    This is done before checking for potential overlapping edge genes. Gene coordinates are expected to be in ascending order, and no circularity is taken into account here.

    :param gene: Gene object to be compared and potentially merged with new_gene_info.
    :param new_gene_info: Dictionary containing information about the new gene.
    :param max_separation: Maximum allowed separation between gene coordinates for merging. Default is 10.


    :raises AssertionError: If the start position is greater than the stop position in new_gene_info.
    :raises ValueError: If the coordinates of genes are too far apart to merge, or if the gene attributes do not match.
    """

    # Compare attributes of the existing gene with new_gene_info
    comparison = [
        gene.strand == new_gene_info["strand"],
        gene.type == new_gene_info["type"],
        gene.product == new_gene_info["product"],
        gene.name == new_gene_info["name"],
        gene.local_identifier == new_gene_info["local_identifier"],
    ]

    if all(comparison):
        # The new gene info seems concordant with the gene object. We can try to merge them
        assert (
            new_gene_info["start"] <= new_gene_info["stop"]
        ), "Start is greater than stop. Incorrect coordinates."

        new_gene_is_before = (
            gene.strand == "+" and new_gene_info["start"] < gene.start
        ) or (gene.strand == "-" and new_gene_info["start"] > gene.start)
        if new_gene_is_before:
            # new gene start if before the current gene
            # so its frame if used
            gene.frame = new_gene_info["frame"]

        # Add new coordinates to gene's coordinates
        gene.coordinates = sorted(
            gene.coordinates + [(new_gene_info["start"], new_gene_info["stop"])]
        )

        # Check if the coordinates are within the allowed maximum separation
        first_stop = gene.coordinates[0][1]
        for start, _ in gene.coordinates[1:]:
            if abs(start - first_stop) > max_separation:
                # This is maybe to restrictive but let's go with that first.
                raise ValueError(
                    f"The coordinates of genes are too far apart ({abs(start - first_stop)}nt). This is unexpected. "
                    f"Gene coordinates : {gene.coordinates}"
                )

        # Update start and stop positions based on new coordinates
        gene.start, gene.stop = gene.coordinates[0][0], gene.coordinates[-1][1]

        logging.getLogger("PPanGGOLiN").debug(
            f"Gene {new_gene_info['local_identifier']} is found in multiple parts. "
            "These parts are merged into one gene. "
            f"New gene coordinates: {gene.coordinates}"
        )

    else:
        detailed_comparison = {
            "same strand": f"{gene.strand} == {new_gene_info['strand']}: {gene.strand == new_gene_info['strand']}",
            "same type": f"{gene.type} == {new_gene_info['type']}: {gene.type == new_gene_info['type']}",
            "same product": f"{gene.product} == {new_gene_info['product']}: {gene.product == new_gene_info['product']}",
            "same name": f"{gene.name} == {new_gene_info['name']}: {gene.name == new_gene_info['name']}",
            "same local identifier": f"{gene.local_identifier} == {new_gene_info['local_identifier']}: {gene.local_identifier == new_gene_info['local_identifier']}",
        }

        raise ValueError(
            f"Two genes have the same id attributes but different info in some key attribute. {detailed_comparison}"
        )


def correct_putative_overlaps(contigs: Iterable[Contig]):
    """
    Corrects putative overlaps in gene coordinates for circular contigs.

    :param contigs: Iterable of Contig objects representing circular contigs.

    :raises ValueError: If a gene start position is higher than the length of the contig.
    """

    for contig in contigs:
        for gene in contig.genes:
            if gene.stop > len(contig):
                # Adjust gene coordinates to handle circular contig
                gene.start = 1  # Start gene at the beginning of the contig

                new_coordinates = []
                for start, stop in gene.coordinates:

                    if start > len(contig):
                        if len(new_coordinates) == 0:
                            raise ValueError(
                                f"First gene start position ({start}) is higher than contig "
                                f"length ({len(contig)}). This case is not handled."
                            )

                        new_start = start - len(contig)
                        new_stop = stop - len(contig)

                        new_coordinates.append((new_start, new_stop))

                        warn_msg = (
                            f"Start position ({start}) for gene {gene.name} is higher than contig {contig.name}"
                            f" length ({len(contig)}). New coordinate are {new_coordinates}"
                        )
                        logging.getLogger("PPanGGOLiN").warning(warn_msg)
                    elif stop > len(contig):
                        # Handle overlapping gene
                        new_stop = len(contig)
                        next_stop = stop - len(contig)
                        next_start = 1

                        new_coordinates.append((start, new_stop))
                        new_coordinates.append((next_start, next_stop))

                    else:
                        new_coordinates.append((start, stop))

                    logging.getLogger("PPanGGOLiN").debug(
                        f"Gene ({gene.ID} {gene.local_identifier}) coordinates ({gene.coordinates}) exceeded contig length ({len(contig)}). "
                        f"This is likely because the gene overlaps the edge of the contig. "
                        f"Adjusted gene coordinates: {new_coordinates}"
                    )

                gene.coordinates = new_coordinates


def read_anno_file(
    organism_name: str,
    filename: Path,
    circular_contigs: list,
    pseudo: bool = False,
    translation_table: int = 11,
) -> Tuple[Organism, bool]:
    """
    Read a GBFF file for one organism

    :param organism_name: Name of the organism
    :param filename: Path to the corresponding file
    :param circular_contigs: list of sequence in contig
    :param pseudo: allow to read pseudogene
    :param translation_table: Translation table (genetic code) to use when /transl_table is missing from CDS tags.

    :return: Annotated organism for pangenome and true for sequence in file
    """
    global ctg_counter
    filetype = detect_filetype(filename)
    if filetype == "gff":
        try:
            org, has_fasta = read_org_gff(
                organism_name, filename, circular_contigs, pseudo, translation_table
            )
        except Exception as err:
            raise Exception(
                f"Reading the gff3 file '{filename}' raised an error. {err}"
            )
        else:
            return org, has_fasta
    elif filetype == "gbff":
        try:
            org, has_fasta = read_org_gbff(
                organism_name, filename, circular_contigs, pseudo, translation_table
            )
        except Exception as err:
            raise Exception(
                f"Reading the gbff file '{filename}' raised an error. {err}"
            )
        else:
            return org, has_fasta
    elif filetype == "fasta":
        raise ValueError(
            f"Invalid file type provided for parameter '--anno'. "
            f"The file '{filename}' looks like a fasta file. "
            "Please use a .gff or .gbff file. You may be able to use --fasta instead of --anno."
        )

    else:
        raise ValueError(
            f"Invalid file type provided for parameter '--anno'. "
            f"The file '{filename}' appears to be of type '{filetype}'. Please use .gff or .gbff files."
        )


def chose_gene_identifiers(pangenome: Pangenome) -> bool:
    """
    Parses the pangenome genes to decide whether to use local_identifiers or ppanggolin generated gene identifiers.
    If the local identifiers are unique within the pangenome they are picked, otherwise ppanggolin ones are used.

    :param pangenome: input pangenome

    :return: Boolean stating True if local identifiers are used, and False otherwise
    """

    if local_identifiers_are_unique(pangenome.genes):

        for gene in pangenome.genes:
            gene.ID = (
                gene.local_identifier
            )  # Erase ppanggolin generated gene ids and replace with local identifiers
            gene.local_identifier = (
                ""  # this is now useless, setting it to default value
            )
        pangenome._mk_gene_getter()  # re-build the gene getter
        return True

    else:
        return False


def local_identifiers_are_unique(genes: Iterable[Gene]) -> bool:
    """
    Check if local_identifiers of genes are uniq in order to decide if they should be used as gene id.

    :param genes: Iterable of gene objects

    :return: Boolean stating True if local identifiers are uniq, and False otherwise
    """
    gene_id_2_local = {}
    local_to_gene_id = {}
    for gene in genes:
        gene_id_2_local[gene.ID] = gene.local_identifier
        local_to_gene_id[gene.local_identifier] = gene.ID
        if len(local_to_gene_id) != len(gene_id_2_local):
            # then, there are non unique local identifiers
            return False
    # if we reach this line, local identifiers are unique within the pangenome
    return True


def read_annotations(
    pangenome: Pangenome,
    organisms_file: Path,
    cpu: int = 1,
    pseudo: bool = False,
    translation_table: int = 11,
    disable_bar: bool = False,
):
    """
    Read the annotation from GBFF file

    :param pangenome: pangenome object
    :param organisms_file: List of GBFF files for each organism
    :param cpu: number of CPU cores to use
    :param pseudo: allow to read pseudogene
    :param translation_table: Translation table (genetic code) to use when /transl_table is missing from CDS tags.
    :param disable_bar: Disable the progress bar
    """

    logging.getLogger("PPanGGOLiN").info(
        f"Reading {organisms_file.name} the list of genome files ..."
    )

    pangenome.status["geneSequences"] = "Computed"
    # we assume there are gene sequences in the annotation files,
    # unless a gff file without fasta is met (which is the only case where sequences can be absent)
    args = []
    for line in read_compressed_or_not(organisms_file):
        if not line.strip() or line.strip().startswith("#"):
            continue
        elements = [el.strip() for el in line.split("\t")]
        org_path = Path(elements[1])
        name = elements[0]
        circular_contigs = elements[2:]
        if (
            not org_path.exists()
        ):  # Check tsv sanity test if it's not one it's the other
            org_path = organisms_file.parent.joinpath(org_path)

        args.append((name, org_path, circular_contigs, pseudo, translation_table))

    with ProcessPoolExecutor(
        mp_context=get_context("fork"),
        max_workers=cpu,
        initializer=init_contig_counter,
        initargs=(contig_counter,),
    ) as executor:
        with tqdm(total=len(args), unit="file", disable=disable_bar) as progress:
            futures = []

            for fn_args in args:
                future = executor.submit(read_anno_file, *fn_args)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                org, has_dna_sequence = future.result()
                pangenome.add_organism(org)

                if not has_dna_sequence:
                    pangenome.status["geneSequences"] = "No"

    # decide whether we use local ids or ppanggolin ids.
    used_local_identifiers = chose_gene_identifiers(pangenome)
    if used_local_identifiers:
        logging.getLogger("PPanGGOLiN").info(
            "gene identifiers used in the provided annotation files were unique, "
            "PPanGGOLiN will use them."
        )
    else:
        logging.getLogger("PPanGGOLiN").info(
            "gene identifiers used in the provided annotation files were not unique, "
            "PPanGGOLiN will use self-generated identifiers."
        )

    pangenome.status["genomesAnnotated"] = "Computed"
    pangenome.parameters["annotate"] = {}
    pangenome.parameters["annotate"][
        "# used_local_identifiers"
    ] = used_local_identifiers
    pangenome.parameters["annotate"]["use_pseudo"] = pseudo
    pangenome.parameters["annotate"]["# read_annotations_from_file"] = True

    if any(genome.has_metadata() for genome in pangenome.organisms):
        pangenome.status["metadata"]["genomes"] = "Computed"
        pangenome.status["metasources"]["genomes"].append("annotation_file")

    if any(contig.has_metadata() for contig in pangenome.contigs):
        pangenome.status["metadata"]["contigs"] = "Computed"
        pangenome.status["metasources"]["contigs"].append("annotation_file")


def get_gene_sequences_from_fastas(
    pangenome: Pangenome, fasta_files: Path, disable_bar: bool = False
):
    """
    Get gene sequences from fastas

    :param pangenome: Input pangenome
    :param fasta_files: list of fasta file
    :param disable_bar: Flag to disable progress bar
    """
    fasta_dict = {}
    for line in read_compressed_or_not(fasta_files):
        elements = [el.strip() for el in line.split("\t")]
        if len(elements) <= 1:
            logging.getLogger("PPanGGOLiN").error(
                "No tabulation separator found in genome file"
            )
            exit(1)
        try:
            org = pangenome.get_organism(elements[0])
        except KeyError:
            raise KeyError(
                f"One of the genome in your '{fasta_files}' was not found in the pan."
                f" This might mean that the genome names between your annotation file and "
                f"your fasta file are different."
            )
        with read_compressed_or_not(Path(elements[1])) as currFastaFile:
            fasta_dict[org] = get_contigs_from_fasta_file(org, currFastaFile)

            # When dealing with GFF files, some genes may have coordinates extending beyond the actual
            # length of contigs, especially when they overlap the edges. This usually needs to be split
            # into two parts to handle the circular genome wrapping.
            # If the GFF file lacks associated FASTA sequences and it was not possible to determine the
            # contig length from the GFF file, we must apply this correction while parsing the external FASTA file.

            correct_putative_overlaps(org.contigs)

    if set(pangenome.organisms) > set(fasta_dict.keys()):
        missing = pangenome.number_of_organisms - len(
            set(pangenome.organisms) & set(fasta_dict.keys())
        )
        raise KeyError(
            f"Not all of your pangenome genomes are present within the provided fasta file. "
            f"{missing} are missing (out of {pangenome.number_of_organisms})."
        )

    elif pangenome.number_of_organisms < len(fasta_dict):
        # Indicates that all organisms in the pangenome are present in the provided FASTA file,
        # but additional genomes are also detected in the file.
        diff_genomes = len(fasta_dict) - pangenome.number_of_organisms
        logging.getLogger("PPanGGOLiN").warning(
            f"The provided fasta file contains {diff_genomes} "
            "additional genomes compared to the pangenome."
        )

    with tqdm(
        total=pangenome.number_of_genes,
        unit="gene",
        disable=disable_bar,
        desc="Add sequences to genes",
    ) as bar:
        for org in pangenome.organisms:
            for contig in org.contigs:
                try:
                    for gene in contig.genes:
                        gene.add_sequence(
                            get_dna_sequence(fasta_dict[org][contig.name], gene)
                        )
                        bar.update()
                    # for rna in contig.RNAs:
                    #     rna.add_sequence(get_dna_sequence(fasta_dict[org][contig.name], rna))
                except KeyError:
                    msg = (
                        f"Fasta file for genome {org.name} did not have the contig {contig.name} "
                        f"that was read from the annotation file. "
                    )
                    msg += (
                        f"The provided contigs in the fasta were : "
                        f"{', '.join(fasta_dict[org].keys())}."
                    )
                    raise KeyError(msg)
    pangenome.status["geneSequences"] = "Computed"


def annotate_pangenome(
    pangenome: Pangenome,
    fasta_list: Path,
    tmpdir: str,
    cpu: int = 1,
    translation_table: int = 11,
    kingdom: str = "bacteria",
    norna: bool = False,
    allow_overlap: bool = False,
    procedure: str = None,
    disable_bar: bool = False,
):
    """
    Main function to annotate a pangenome

    :param pangenome: Pangenome with gene families to align with the given input sequences
    :param fasta_list: List of fasta file containing sequences that will be base of pangenome
    :param tmpdir: Path to temporary directory
    :param cpu: number of CPU cores to use
    :param translation_table: Translation table (genetic code) to use.
    :param kingdom: Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation.
    :param norna: Use to avoid annotating RNA features.
    :param allow_overlap: Use to not remove genes overlapping with RNA features
    :param procedure: prodigal procedure used
    :param disable_bar: Disable the progress bar
    """

    logging.getLogger("PPanGGOLiN").info(
        f"Reading {fasta_list} the list of genome files"
    )

    arguments = []  # Argument given to annotate organism in same order than prototype
    for line in read_compressed_or_not(fasta_list):

        elements = [el.strip() for el in line.split("\t")]
        org_path = Path(elements[1])

        if (
            not org_path.exists()
        ):  # Check tsv sanity test if it's not one it's the other
            org_path = fasta_list.parent.joinpath(org_path)

        arguments.append(
            (
                elements[0],
                org_path,
                elements[2:],
                tmpdir,
                translation_table,
                norna,
                kingdom,
                allow_overlap,
                procedure,
            )
        )

    if len(arguments) == 0:
        raise Exception("There are no genomes in the provided file")

    logging.getLogger("PPanGGOLiN").info(
        f"Annotating {len(arguments)} genomes using {cpu} cpus..."
    )
    with ProcessPoolExecutor(
        mp_context=get_context("fork"),
        max_workers=cpu,
        initializer=init_contig_counter,
        initargs=(contig_counter,),
    ) as executor:
        with tqdm(total=len(arguments), unit="file", disable=disable_bar) as progress:
            futures = []

            for fn_args in arguments:
                future = executor.submit(annotate_organism, *fn_args)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                pangenome.add_organism(future.result())

    logging.getLogger("PPanGGOLiN").info("Done annotating genomes")
    pangenome.status["genomesAnnotated"] = "Computed"  # the pangenome is now annotated.
    pangenome.status["geneSequences"] = (
        "Computed"  # the gene objects have their respective gene sequences.
    )
    pangenome.parameters["annotate"] = {}
    pangenome.parameters["annotate"]["norna"] = norna
    pangenome.parameters["annotate"]["kingdom"] = kingdom
    pangenome.parameters["annotate"]["translation_table"] = translation_table
    pangenome.parameters["annotate"]["prodigal_procedure"] = (
        None if procedure is None else procedure
    )
    pangenome.parameters["annotate"]["allow_overlap"] = allow_overlap
    pangenome.parameters["annotate"]["# read_annotations_from_file"] = False


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    check_annotate_args(args)
    filename = mk_file_name(args.basename, args.output, args.force)
    pangenome = Pangenome()
    if args.fasta is not None and args.anno is None:
        annotate_pangenome(
            pangenome,
            args.fasta,
            tmpdir=args.tmpdir,
            cpu=args.cpu,
            procedure=args.prodigal_procedure,
            translation_table=args.translation_table,
            kingdom=args.kingdom,
            norna=args.norna,
            allow_overlap=args.allow_overlap,
            disable_bar=args.disable_prog_bar,
        )
    elif args.anno is not None:
        # TODO add warning for option not compatible with read_annotations
        read_annotations(
            pangenome,
            args.anno,
            cpu=args.cpu,
            pseudo=args.use_pseudo,
            translation_table=args.translation_table,
            disable_bar=args.disable_prog_bar,
        )
        if pangenome.status["geneSequences"] == "No":
            if args.fasta:
                logging.getLogger("PPanGGOLiN").info(
                    f"Get sequences from FASTA file: {args.fasta}"
                )
                get_gene_sequences_from_fastas(
                    pangenome, args.fasta, disable_bar=args.disable_prog_bar
                )
            else:
                logging.getLogger("PPanGGOLiN").warning(
                    "You provided gff files without sequences, "
                    "and you did not provide fasta sequences. "
                    "Thus it was not possible to get the gene sequences."
                )
                logging.getLogger("PPanGGOLiN").warning(
                    "You will be able to proceed with your analysis "
                    "ONLY if you provide the clustering results in the next step."
                )

            if pangenome.contig_lengths_unavailable():
                raise ValueError(
                    "Unable to determine contig lengths from the provided GFF files. "
                    "Contig length must be specified using the ##sequence-region pragma. "
                    "Additionally, no FASTA sequences were provided. "
                    "As a result, contig lengths cannot be inferred.\n"
                    "To resolve this, please provide a FASTA file using the '--fasta' option, "
                    "or modify your GFF files to include sufficient information to deduce contig lengths."
                )

        else:
            if args.fasta:
                logging.getLogger("PPanGGOLiN").warning(
                    "You provided fasta sequences "
                    "but your gff files were already with sequences."
                    "PPanGGOLiN will use sequences in GFF and not from your fasta."
                )
    write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser(
        "annotate", formatter_class=argparse.RawTextHelpFormatter
    )
    parser_annot(parser)
    return parser


def parser_annot(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of annotate command

    :param parser: parser for annotate argument
    """
    date = time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())
    required = parser.add_argument_group(
        title="Required arguments",
        description="One of the following arguments is required :",
    )
    required.add_argument(
        "--fasta",
        required=False,
        type=Path,
        help="A tab-separated file listing the genome names, and the fasta filepath of its genomic "
        "sequence(s) (the fastas can be compressed with gzip). One line per genome.",
    )
    required.add_argument(
        "--anno",
        required=False,
        type=Path,
        help="A tab-separated file listing the genome names, and the gff/gbff filepath of its "
        "annotations (the files can be compressed with gzip). One line per genome. "
        "If this is provided, those annotations will be used.",
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
        "--allow_overlap",
        required=False,
        action="store_true",
        default=False,
        help="Use to not remove genes overlapping with RNA features.",
    )
    optional.add_argument(
        "--norna",
        required=False,
        action="store_true",
        default=False,
        help="Use to avoid annotating RNA features.",
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
        "--translation_table",
        required=False,
        type=int,
        default=11,
        help="Translation table (genetic code) to use.",
    )
    optional.add_argument(
        "--basename",
        required=False,
        default="pangenome",
        help="basename for the output file",
    )
    optional.add_argument(
        "--use_pseudo",
        required=False,
        action="store_true",
        help="In the context of provided annotation, use this option to read pseudogenes. "
        "(Default behavior is to ignore them)",
    )
    optional.add_argument(
        "-p",
        "--prodigal_procedure",
        required=False,
        type=str.lower,
        choices=["single", "meta"],
        default=None,
        help="Allow to force the prodigal procedure. "
        "If nothing given, PPanGGOLiN will decide in function of contig length",
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
        type=str,
        default=Path(tempfile.gettempdir()),
        help="directory for storing temporary files",
    )


if __name__ == "__main__":
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser_annot(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
