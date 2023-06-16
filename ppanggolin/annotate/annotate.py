#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from multiprocessing import get_context
import os
from pathlib import Path
import tempfile
import time

# installed libraries
from tqdm import tqdm

# local libraries
from ppanggolin.annotate.synta import annotate_organism, read_fasta, get_dna_sequence
from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Organism, Gene, RNA, Contig
from ppanggolin.utils import read_compressed_or_not, mk_file_name, detect_filetype
from ppanggolin.formats import write_pangenome


def check_annotate_args(args):
    """Check That the given arguments are usable

    :param args: All arguments provide by user

    :raise Exception:
    """
    if args.fasta is None and args.anno is None:
        raise Exception("You must provide at least a file with the --fasta option to annotate from sequences, "
                        "or a file with the --gff option to load annotations from.")


def create_gene(org: Organism, contig: Contig, gene_counter: int, rna_counter: int, gene_id: str, dbxref: set,
                start: int, stop: int, strand: str, gene_type: str, position: int = None, gene_name: str = "",
                product: str = "", genetic_code: int = 11, protein_id: str = ""):
    """
    Create a Gene object and associate to contig and Organism

    :param org: Organism to add gene
    :param contig: Contig to add gene
    :param gene_counter: Gene counter to name gene
    :param rna_counter: RNA counter to name RNA
    :param gene_id: local identifier
    :param dbxref: cross-reference to external DB
    :param start: Gene start position
    :param stop: Gene stop position
    :param strand: gene strand association
    :param gene_type: gene type
    :param position: position in contig
    :param gene_name: Gene name
    :param product: Function of gene
    :param genetic_code: Genetic code used
    :param protein_id: Protein identifier
    """
    if any('MaGe' or 'SEED' in dbref for dbref in dbxref):
        if gene_name == "":
            gene_name = gene_id
        for val in dbxref:
            if 'MaGe' in val:
                gene_id = val.split(':')[1]
                break
            if 'SEED' in val:
                gene_id = val.split(':')[1]
                break
    if gene_type == "CDS":
        if gene_id == "":
            gene_id = protein_id
            # on rare occasions, there are no 'locus_tag' from downloaded .gbk file.
            # So we use the protein_id field instead. (which is not supposed to be unique,
            # but was when cases like this were encountered)

        new_gene = Gene(org.name + "_CDS_" + str(gene_counter).zfill(4))
        new_gene.fill_annotations(start=start, stop=stop, strand=strand, gene_type=gene_type, name=gene_name,
                                  position=position, product=product, local_identifier=gene_id,
                                  genetic_code=genetic_code)
        contig.add_gene(new_gene)
    else:  # if not CDS, it is RNA
        new_gene = RNA(org.name + "_RNA_" + str(rna_counter).zfill(4))
        new_gene.fill_annotations(start=start, stop=stop, strand=strand, gene_type=gene_type, name=gene_name,
                                  product=product)
        contig.add_rna(new_gene)
    new_gene.fill_parents(org, contig)


def read_org_gbff(organism: str, gbff_file_path: Path, circular_contigs: list, pseudo: bool = False) -> (
        Organism, bool):
    """
    Read a GBFF file and fills Organism, Contig and Genes objects based on information contained in this file

    :param organism: Organism name
    :param gbff_file_path: Path to corresponding GFF file
    :param circular_contigs: list of contigs
    :param pseudo: Allow to read pseudogène

    :return: Organism complete and true for sequence in file
    """
    org = Organism(organism)

    logging.debug(f"Extracting genes informations from the given gbff {gbff_file_path.name}")
    # revert the order of the file, to read the first line first.
    lines = read_compressed_or_not(gbff_file_path).readlines()[::-1]
    gene_counter = 0
    rna_counter = 0
    while len(lines) != 0:
        line = lines.pop()
        # beginning of contig
        contig = None
        set_contig = False
        is_circ = False
        contig_locus_id = None
        if line.startswith('LOCUS'):
            if "CIRCULAR" in line.upper():
                # this line contains linear/circular word telling if the dna sequence is circularized or not
                is_circ = True
            contig_locus_id = line.split()[1]
            # If contig_id is not specified in VERSION afterward like with Prokka, in that case we use the one in LOCUS
            while not line.startswith('FEATURES'):
                if line.startswith('VERSION'):
                    contig_id = line[12:].strip()
                    if contig_id != "":
                        if contig_id in circular_contigs:
                            is_circ = True
                        contig = org.get_contig(contig_id, is_circ)
                        set_contig = True
                line = lines.pop()
        if not set_contig:
            # if no contig ids were filled after VERSION, we use what was found in LOCUS for the contig ID.
            # Should be unique in a dataset, but if there's an update the contig ID
            # might still be the same even though it should not(?)
            if contig_locus_id in circular_contigs:
                is_circ = True
            contig = org.get_contig(contig_locus_id, is_circ)
        # start of the feature object.
        dbxref = set()
        gene_name = ""
        product = ""
        locus_tag = ""
        obj_type = ""
        protein_id = ""
        genetic_code = ""
        useful_info = False
        start = None
        end = None
        strand = None
        line = lines.pop()
        while not line.startswith("ORIGIN"):
            curr_type = line[5:21].strip()
            if curr_type != "":
                if useful_info:
                    create_gene(org, contig, gene_counter, rna_counter, locus_tag, dbxref, start, end, strand, obj_type,
                                len(contig.genes), gene_name, product, genetic_code, protein_id)
                    if obj_type == "CDS":
                        gene_counter += 1
                    else:
                        rna_counter += 1
                useful_info = False
                obj_type = curr_type
                if obj_type in ['CDS', 'rRNA', 'tRNA']:
                    dbxref = set()
                    gene_name = ""
                    try:
                        if 'join' not in line[21:]:
                            useful_info = True
                            if line[21:].startswith('complement('):
                                strand = "-"
                                start, end = line[32:].strip().replace(
                                    ')', '').split("..")
                            else:
                                strand = "+"
                                start, end = line[21:].strip().split('..')
                            if '>' in start or '<' in start or '>' in end or '<' in end:
                                if not pseudo:
                                    # pseudogene likely
                                    useful_info = False
                                else:
                                    start = start.replace('>', '').replace('<', '')
                                    end = end.replace('>', '').replace('<', '')
                    except ValueError:
                        pass
                        # don't know what to do with that, ignoring for now.
                        # there is a protein with a frameshift mecanism.
            elif useful_info:  # current info goes to current objtype, if it's useful.
                if line[21:].startswith("/db_xref"):
                    dbxref.add(line.split("=")[1].replace('"', '').strip())
                elif line[21:].startswith("/locus_tag"):
                    locus_tag = line.split("=")[1].replace('"', '').strip()
                elif line[21:].startswith("/protein_id"):
                    protein_id = line.split("=")[1].replace('"', '').strip()
                elif line[21:].startswith('/gene'):  # gene name
                    gene_name = line.split("=")[1].replace('"', '').strip()
                elif line[21:].startswith('/transl_table'):
                    genetic_code = int(line.split("=")[1].replace('"', '').strip())
                elif line[21:].startswith('/product'):  # need to loop as it can be more than one line long
                    product = line.split('=')[1].replace('"', '').strip()
                    if line.count('"') == 1:  # then the product line is on multiple lines
                        line = lines.pop()
                        product += line.strip().replace('"', '')
                        while line.count('"') != 1:
                            line = lines.pop()
                            product += line.strip().replace('"', '')
                # if it's a pseudogene, we're not keeping it, unless pseudo
                elif line[21:].startswith("/pseudo") and not pseudo:
                    useful_info = False
                # that's probably a 'stop' codon into selenocystein.
                elif line[21:].startswith("/transl_except") and not pseudo:
                    useful_info = False
            line = lines.pop()
            # end of contig
        if useful_info:  # saving the last element...
            create_gene(org, contig, gene_counter, rna_counter, locus_tag, dbxref, start, end, strand, obj_type,
                        len(contig.genes), gene_name, product, genetic_code, protein_id)
            if obj_type == "CDS":
                gene_counter += 1
            else:
                rna_counter += 1

        # now extract the gene sequences
        line = lines.pop()  # first sequence line.
        # if the seq was to be gotten, it would be here.
        sequence = ""
        while not line.startswith('//'):
            sequence += line[10:].replace(" ", "").strip().upper()
            line = lines.pop()
        # get each gene's sequence.
        for gene in contig.genes:
            gene.add_dna(get_dna_sequence(sequence, gene))
    return org, True


def read_org_gff(organism: str, gff_file_path: Path, circular_contigs, pseudo: bool = False) -> (Organism, bool):
    """
    Read annotation from GFF file

    :param organism: Organism name
    :param gff_file_path: Path to corresponding GFF file
    :param circular_contigs: list of contigs
    :param pseudo: Allow to read pseudogène

    :return: Organism object and if there are sequences associate or not
    """
    (gff_seqname, _, gff_type, gff_start, gff_end, _, gff_strand, _, gff_attribute) = range(0, 9)

    # missing values : source, score, frame. They are unused.

    def get_gff_attributes(gff_fields: list) -> dict:
        """
        Parses the gff attribute's line and outputs the attributes_get in a dict structure.
        :param gff_fields: a gff line stored as a list. Each element of the list is a column of the gff.

        :return: attributes get
        """
        attributes_field = [f for f in gff_fields[gff_attribute].strip().split(';') if len(f) > 0]
        attributes_get = {}
        for att in attributes_field:
            try:
                (key, value) = att.strip().split('=')
                attributes_get[key.upper()] = value
            except ValueError:
                pass  # we assume that it is a strange, but useless field for our analysis
        return attributes_get

    def get_id_attribute(attributes_dict: dict) -> str:
        """
        Gets the ID of the element from which the provided attributes_get were extracted.
        Raises an error if no ID is found.
        :param attributes_dict: attributes from one gff line

        :return: CDS identifier
        """
        element_id = attributes_dict.get("ID")
        if not element_id:
            raise Exception(f"Each CDS type of the gff files must own a unique ID attribute. "
                            f"Not the case for file: {gff_file_path}")
        return element_id

    contig = None  # initialize contig
    has_fasta = False
    fasta_string = ""
    org = Organism(organism)
    gene_counter = 0
    rna_counter = 0
    with read_compressed_or_not(gff_file_path) as gff_file:
        for line in gff_file:
            if has_fasta:
                fasta_string += line
                continue
            elif line.startswith('##', 0, 2):
                if line.startswith('FASTA', 2, 7):
                    has_fasta = True
                elif line.startswith('sequence-region', 2, 17):
                    fields = [el.strip() for el in line.split()]
                    contig = org.get_contig(fields[1], True if fields[1] in circular_contigs else False)
                continue
            elif line.startswith('#'):  # comment lines to be ignores by parsers
                continue
            elif line == "":  # empty lines are not expected, but they do not carry information, so we'll ignore them
                continue
            fields_gff = [el.strip() for el in line.split('\t')]
            attributes = get_gff_attributes(fields_gff)
            pseudogene = False
            if fields_gff[gff_type] == 'region':
                if fields_gff[gff_seqname] in circular_contigs:
                    contig.is_circular = True
            elif fields_gff[gff_type] == 'CDS' or "RNA" in fields_gff[gff_type]:
                gene_id = attributes.get("PROTEIN_ID")
                # if there is a 'PROTEIN_ID' attribute, it's where the ncbi stores the actual gene ids, so we use that.
                if gene_id is None:
                    # if it's not found, we get the one under the 'ID' field which must exist
                    # (otherwise not a gff3 compliant file)
                    gene_id = get_id_attribute(attributes)
                try:
                    name = attributes.pop('NAME')
                except KeyError:
                    try:
                        name = attributes.pop('GENE')
                    except KeyError:
                        name = ""
                if "pseudo" in attributes or "pseudogene" in attributes:
                    pseudogene = True
                try:
                    product = attributes.pop('PRODUCT')
                except KeyError:
                    product = ""

                try:
                    genetic_code = int(attributes.pop("TRANSL_TABLE"))
                except KeyError:
                    genetic_code = 11
                if contig is None or contig.name != fields_gff[gff_seqname]:
                    # get the current contig
                    contig = org.get_contig(fields_gff[gff_seqname],
                                            True if fields_gff[gff_seqname] in circular_contigs else False)

                if fields_gff[gff_type] == "CDS" and (not pseudogene or (pseudogene and pseudo)):
                    gene = Gene(org.name + "_CDS_" + str(gene_counter).zfill(4))
                    # here contig is filled in order, so position is the number of genes already stored in the contig.
                    gene.fill_annotations(start=int(fields_gff[gff_start]), stop=int(fields_gff[gff_end]),
                                          strand=fields_gff[gff_strand], gene_type=fields_gff[gff_type], name=name,
                                          position=len(contig.genes), product=product, local_identifier=gene_id,
                                          genetic_code=genetic_code)
                    gene.fill_parents(org, contig)
                    contig.add_gene(gene)
                    gene_counter += 1
                elif "RNA" in fields_gff[gff_type]:
                    rna = RNA(org.name + "_CDS_" + str(rna_counter).zfill(4))
                    rna.fill_annotations(start=int(fields_gff[gff_start]), stop=int(fields_gff[gff_end]),
                                         strand=fields_gff[gff_strand], gene_type=fields_gff[gff_type], name=name,
                                         product=product, local_identifier=gene_id)
                    rna.fill_parents(org, contig)
                    contig.add_rna(rna)
                    rna_counter += 1

    # GET THE FASTA SEQUENCES OF THE GENES
    if has_fasta and fasta_string != "":
        contig_sequences, _ = read_fasta(org, fasta_string.split('\n'))  # _ is total contig length
        for contig in org.contigs:
            for gene in contig.genes:
                gene.add_dna(get_dna_sequence(contig_sequences[contig.name], gene))
            for rna in contig.RNAs:
                rna.add_dna(get_dna_sequence(contig_sequences[contig.name], rna))
    return org, has_fasta


def launch_read_anno(args: tuple) -> (Organism, bool):
    """ Allow to launch in multiprocessing the read of genome annotation

    :param args: Pack of argument for annotate_organism function

    :return: Organism object for pangenome
    """
    return read_anno_file(*args)


def read_anno_file(organism_name: str, filename: Path, circular_contigs: list, pseudo: bool = False) -> (
        Organism, bool):
    """
    Read a GBFF file for one organism

    :param organism_name: Name of the organism
    :param filename: Path to the corresponding file
    :param circular_contigs: list of sequence in contig
    :param pseudo: allow to read pseudogène

    :return: Annotated organism for pangenome
    """
    filetype = detect_filetype(filename)
    if filetype == "gff":
        try:
            return read_org_gff(organism_name, filename, circular_contigs, pseudo)
        except Exception:
            raise Exception(f"Reading the gff3 file '{filename}' raised an error.")
    elif filetype == "gbff":
        try:
            return read_org_gbff(organism_name, filename, circular_contigs, pseudo)
        except Exception:
            raise Exception(f"Reading the gbff file '{filename}' raised an error.")
    else:  # Fasta type obligatory because unknow raise an error in detect_filetype function
        raise Exception("Wrong file type provided. This looks like a fasta file. "
                        "You may be able to use --fasta instead.")


def chose_gene_identifiers(pangenome) -> bool:
    """
    Parses the pangenome genes to decide whether to use local_identifiers or ppanggolin generated gene identifiers.
    If the local identifiers are unique within the pangenome they are picked, otherwise ppanggolin ones are used.

    :param pangenome: input pangenome

    :return: Boolean stating True if local identifiers are used, and False otherwise
    """
    gene_id_2_local = {}
    local_to_gene_id = {}
    for gene in pangenome.genes:
        gene_id_2_local[gene.ID] = gene.local_identifier
        local_to_gene_id[gene.local_identifier] = gene.ID
        if len(local_to_gene_id) != len(gene_id_2_local):
            # then, there are non unique local identifiers
            return False
    # if we reach this line, local identifiers are unique within the pangenome
    for gene in pangenome.genes:
        gene.ID = gene.local_identifier  # Erase ppanggolin generated gene ids and replace with local identifiers
        gene.local_identifier = ""  # this is now useless, setting it to default value
    pangenome._mk_gene_getter()  # re-build the gene getter
    return True


def read_annotations(pangenome: Pangenome, organisms_file: Path, cpu: int = 1, pseudo: bool = False,
                     disable_bar: bool = False):
    """
    Read the annotation from GBFF file

    :param pangenome: pangenome object
    :param organisms_file: List of GBFF files for each organism
    :param cpu: number of CPU cores to use
    :param pseudo: allow to read pseudogène
    :param disable_bar: Disable the progresse bar
    """
    logging.info(f"Reading {organisms_file.name} the list of organism files ...")

    pangenome.status["geneSequences"] = "Computed"
    # we assume there are gene sequences in the annotation files,
    # unless a gff file without fasta is met (which is the only case where sequences can be absent)
    args = []
    for line in read_compressed_or_not(organisms_file):
        elements = [el.strip() for el in line.split("\t")]
        if len(elements) <= 1:
            raise Exception(f"No tabulation separator found in given --fasta file: '{organisms_file}'")
        org_path = Path(elements[1])
        if not org_path.exists():  # Check tsv sanity test if it's not one it's the other
            org_path = organisms_file.parent.joinpath(org_path)
        args.append((elements[0], org_path, elements[2:], pseudo))
    with get_context('fork').Pool(cpu) as p:
        for org, flag in tqdm(p.imap_unordered(launch_read_anno, args), unit="file", total=len(args),
                              disable=disable_bar):
            pangenome.add_organism(org)
            if not flag:
                pangenome.status["geneSequences"] = "No"

    # decide whether we use local ids or ppanggolin ids.
    used_local_identifiers = chose_gene_identifiers(pangenome)
    if used_local_identifiers:
        logging.info("gene identifiers used in the provided annotation files were unique, "
                     "PPanGGOLiN will use them.")
    else:
        logging.info("gene identifiers used in the provided annotation files were not unique, "
                     "PPanGGOLiN will use self-generated identifiers.")

    pangenome.status["genomesAnnotated"] = "Computed"
    pangenome.parameters["annotation"] = {}
    pangenome.parameters["annotation"]["used_local_identifiers"] = used_local_identifiers
    pangenome.parameters["annotation"]["read_pseudogenes"] = pseudo
    pangenome.parameters["annotation"]["read_annotations_from_file"] = True


def get_gene_sequences_from_fastas(pangenome, fasta_file):
    """
    Get gene sequences from fastas

    :param pangenome: input pangenome
    :param fasta_file: list of fasta file
    """
    fasta_dict = {}
    for line in read_compressed_or_not(fasta_file):
        elements = [el.strip() for el in line.split("\t")]
        if len(elements) <= 1:
            logging.error("No tabulation separator found in organisms file")
            exit(1)
        try:
            org = pangenome.get_organism(elements[0])
        except KeyError:
            raise KeyError(f"One of the genome in your '{fasta_file}' was not found in the pan."
                           f" This might mean that the genome names between your annotation file and "
                           f"your fasta file are different.")
        with read_compressed_or_not(elements[1]) as currFastaFile:
            fasta_dict[org], _ = read_fasta(org, currFastaFile)
    if set(pangenome.organisms) > set(fasta_dict.keys()):
        missing = len(pangenome.organisms) - len(set(pangenome.organisms) & set(fasta_dict.keys()))
        raise Exception(f"Not all of your pangenome organisms are present within the provided fasta file. "
                        f"{missing} are missing (out of {len(pangenome.organisms)}).")

    for org in pangenome.organisms:
        for contig in org.contigs:
            try:
                for gene in contig.genes:
                    gene.add_dna(get_dna_sequence(fasta_dict[org][contig.name], gene))
                for rna in contig.RNAs:
                    rna.add_dna(get_dna_sequence(fasta_dict[org][contig.name], rna))
            except KeyError:
                msg = f"Fasta file for organism {org.name} did not have the contig {contig.name} " \
                      f"that was read from the annotation file. "
                msg += f"The provided contigs in the fasta were : " \
                       f"{', '.join([contig for contig in fasta_dict[org].keys()])}."
                raise KeyError(msg)
    pangenome.status["geneSequences"] = "Computed"


def launch_annotate_organism(pack: tuple) -> Organism:
    """ Allow to launch in multiprocessing the genome annotation

    :param pack: Pack of argument for annotate_organism function

    :return: Organism object for pangenome
    """
    return annotate_organism(*pack)


def annotate_pangenome(pangenome: Pangenome, fasta_list: Path, tmpdir: str, cpu: int = 1, translation_table: int = 11,
                       kingdom: str = "bacteria", norna: bool = False, overlap: bool = True, procedure: str = None,
                       disable_bar: bool = False):
    """
    Main function to annotate a pangenome

    :param pangenome: Pangenome with gene families to align with the given input sequences
    :param fasta_list: List of fasta file containing sequences that will be base of pangenome
    :param tmpdir: Path to temporary directory
    :param cpu: number of CPU cores to use
    :param translation_table: Translation table (genetic code) to use.
    :param kingdom: Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation.
    :param norna: Use to avoid annotating RNA features.
    :param overlap: Use to not remove genes overlapping with RNA features
    :param procedure: prodigal procedure used
    :param disable_bar: Disable the progresse bar
    """

    logging.info(f"Reading {fasta_list} the list of organism files")

    arguments = []  # Argument given to annotate organism in same order than prototype
    for line in read_compressed_or_not(fasta_list):
        elements = [el.strip() for el in line.split("\t")]
        if len(elements) <= 1:  # TODO remove ? Already tested by check TSV sanity
            raise Exception("No tabulation separator found in organisms file")
        org_path = Path(elements[1])
        if not org_path.exists():  # Check tsv sanity test if it's not one it's the other
            org_path = fasta_list.parent.joinpath(org_path)
        arguments.append((elements[0], org_path, elements[2:], tmpdir, translation_table,
                          norna, kingdom, overlap, procedure))
    if len(arguments) == 0:
        raise Exception("There are no genomes in the provided file")
    logging.info(f"Annotating {len(arguments)} genomes using {cpu} cpus...")
    with get_context('fork').Pool(processes=cpu) as p:
        for organism in tqdm(p.imap_unordered(launch_annotate_organism, arguments), unit="genome",
                             total=len(arguments), disable=disable_bar):
            pangenome.add_organism(organism)
        p.close()
        p.join()

    logging.info("Done annotating genomes")
    pangenome.status["genomesAnnotated"] = "Computed"  # the pangenome is now annotated.
    pangenome.status["geneSequences"] = "Computed"  # the gene objects have their respective gene sequences.
    pangenome.parameters["annotation"] = {}
    pangenome.parameters["annotation"]["remove_Overlapping_CDS"] = overlap
    pangenome.parameters["annotation"]["annotate_RNA"] = True if not norna else False
    pangenome.parameters["annotation"]["kingdom"] = kingdom
    pangenome.parameters["annotation"]["translation_table"] = translation_table
    pangenome.parameters["annotation"]["read_annotations_from_file"] = False


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    check_annotate_args(args)
    filename = mk_file_name(args.basename, args.output, args.force)
    pangenome = Pangenome()
    if args.fasta is not None and args.anno is None:
        annotate_pangenome(pangenome, args.fasta, tmpdir=args.tmpdir, cpu=args.cpu, procedure=args.prodigal_procedure,
                           translation_table=args.translation_table, kingdom=args.kingdom, norna=args.norna,
                           overlap=args.allow_overlap, disable_bar=args.disable_prog_bar)
    elif args.anno is not None:
        read_annotations(pangenome, args.anno, cpu=args.cpu, pseudo=args.use_pseudo, disable_bar=args.disable_prog_bar)
        if pangenome.status["geneSequences"] == "No":
            if args.fasta:
                get_gene_sequences_from_fastas(pangenome, args.fasta)
            else:
                logging.warning("You provided gff files without sequences, and you did not provide "
                                "fasta sequences. Thus it was not possible to get the gene sequences.")
                logging.warning("You will be able to proceed with your analysis ONLY if you provide "
                                "the clustering results in the next step.")

    write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("annotate", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_annot(parser)
    return parser


def parser_annot(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of annotate command

    :param parser: parser for annotate argument
    """
    date = time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())
    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('--fasta', required=False, type=Path,
                          help="A tab-separated file listing the organism names, and the fasta filepath of its genomic "
                               "sequence(s) (the fastas can be compressed with gzip). One line per organism.")
    required.add_argument('--anno', required=False, type=Path,
                          help="A tab-separated file listing the organism names, and the gff/gbff filepath of its "
                               "annotations (the files can be compressed with gzip). One line per organism. "
                               "If this is provided, those annotations will be used.")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('-o', '--output', required=False, type=Path,
                          default=Path(f'ppanggolin_output{date}_PID{str(os.getpid())}'),
                          help="Output directory")
    optional.add_argument('--allow_overlap', required=False, action='store_true', default=False,
                          help="Use to not remove genes overlapping with RNA features.")
    optional.add_argument("--norna", required=False, action="store_true", default=False,
                          help="Use to avoid annotating RNA features.")
    optional.add_argument("--kingdom", required=False, type=str.lower, default="bacteria",
                          choices=["bacteria", "archaea"],
                          help="Kingdom to which the prokaryota belongs to, "
                               "to know which models to use for rRNA annotation.")
    optional.add_argument("--translation_table", required=False, type=int, default=11,
                          help="Translation table (genetic code) to use.")
    optional.add_argument("--basename", required=False, default="pangenome", help="basename for the output file")
    optional.add_argument("--use_pseudo", required=False, action="store_true",
                          help="In the context of provided annotation, use this option to read pseudogenes. "
                               "(Default behavior is to ignore them)")
    optional.add_argument("-p", "--prodigal_procedure", required=False, type=str.lower, choices=["single", "meta"],
                          default=None, help="Allow to force the prodigal procedure. "
                                             "If nothing given, PPanGGOLiN will decide in function of contig length")
    optional.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
    optional.add_argument("--tmpdir", required=False, type=str, default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    import tempfile

    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_annot(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
