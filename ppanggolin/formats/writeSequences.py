#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging

# installed libraries
from typing import TextIO

from tqdm import tqdm

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import write_compressed_or_not, mk_outdir, read_compressed_or_not, restricted_float
from ppanggolin.formats.readBinaries import check_pangenome_info, get_gene_sequences_from_file
from ppanggolin.annotate.annotate import detect_filetype

poss_values_log = "Possible values are 'all', 'persistent', 'shell', 'cloud', 'rgp', 'softcore', " \
                  "'core', 'module_X' with X being a module id."


def write_gene_sequences_from_annotations(pangenome: Pangenome, file_obj: TextIO, list_cds: list = None, add: str = '',
                                          disable_bar: bool = False):
    """
    Writes the CDS sequences given through list_CDS of the Pangenome object to a tmpFile object,
    and adds the str provided through add in front of it.
    Loads the sequences from previously computed or loaded annotations

    :param pangenome: Pangenome object with gene families sequences
    :param file_obj: Output file to write sequences
    :param list_cds: Selected genes
    :param add: Add prefix to gene ID
    :param disable_bar: Disable progress bar
    """
    counter = 0
    if list_cds is None:
        list_cds = pangenome.genes
    logging.getLogger().info("Writing all of the CDS sequences...")
    for gene in tqdm(sorted(list_cds, key=lambda x: x.ID), unit="gene", disable=disable_bar):
        if gene.type == "CDS":
            counter += 1
            file_obj.write('>' + add + gene.ID + "\n")
            file_obj.write(gene.dna + "\n")
    file_obj.flush()


def write_gene_sequences(pangenome: Pangenome, output: str, genes: str, soft_core: float = 0.95,
                         compress: bool = False, disable_bar: bool = False):
    """
    Write all nucleotide CDS sequences

    :param pangenome: Pangenome object with gene families sequences
    :param output: Path to output directory
    :param genes: Selected partition of gene
    :param soft_core: Soft core threshold to use
    :param compress: Compress the file in .gz
    :param disable_bar: Disable progress bar
    """
    logging.getLogger().info("Writing all the gene nucleotide sequences...")
    outname = output + f"/{genes}_genes.fna"

    genefams = select_families(pangenome, genes, "gene nucleotide sequences", soft_core)
    genes_to_write = []

    for fam in genefams:
        genes_to_write.extend(fam.genes)

    logging.getLogger().info(f"There are {len(genes_to_write)} genes to write")
    with write_compressed_or_not(outname, compress) as fasta:
        if pangenome.status["geneSequences"] in ["inFile"]:
            get_gene_sequences_from_file(pangenome.file, fasta, set([gene.ID for gene in genes_to_write]),
                                         disable_bar=disable_bar)
        elif pangenome.status["geneSequences"] in ["Computed", "Loaded"]:
            write_gene_sequences_from_annotations(pangenome, fasta, genes_to_write, disable_bar=disable_bar)
        else:
            # this should never happen if the pangenome has been properly checked before launching this function.
            raise Exception("The pangenome does not include gene sequences")
    logging.getLogger().info(f"Done writing the gene sequences : '{outname}'")


def select_families(pangenome: Pangenome, partition: str, type_name: str, soft_core: float = 0.95) -> set:
    """
    function used to filter down families to the given partition

    :param pangenome: Pangenome object
    :param partition: Selected partition
    :param type_name: Which type of sequence we want. Gene families, protein, gene
    :param soft_core: Soft core threshold to use

    :return: Selected gene families
    """
    genefams = set()
    if partition == 'all':
        logging.getLogger().info(f"Writing all of the {type_name}...")
        genefams = pangenome.gene_families
    elif partition in ['persistent', 'shell', 'cloud']:
        logging.getLogger().info(f"Writing the {type_name} of the {partition}...")
        for fam in pangenome.gene_families:
            if fam.named_partition == partition:
                genefams.add(fam)
    elif partition == "rgp":
        logging.getLogger().info(f"Writing the {type_name} in RGPs...")
        for region in pangenome.regions:
            genefams |= region.families
    elif partition == "softcore":
        logging.getLogger().info(
            f"Writing the {type_name} in {partition} genome, that are present in more than {soft_core} of genomes")
        threshold = pangenome.number_of_organisms() * soft_core
        for fam in pangenome.gene_families:
            if len(fam.organisms) >= threshold:
                genefams.add(fam)
    elif partition == "core":
        logging.getLogger().info(f"Writing the representative {type_name} of the {partition} gene families...")
        for fam in pangenome.gene_families:
            if len(fam.organisms) == pangenome.number_of_organisms():
                genefams.add(fam)
    elif "module_" in partition:
        logging.getLogger().info(f"Writing the representation {type_name} of {partition} gene families...")
        mod_id = int(partition.replace("module_", ""))
        for mod in pangenome.modules:
            # could be way more efficient with a dict structure instead of a set
            if mod.ID == mod_id:
                genefams |= mod.families
                break
    return genefams


def write_fasta_gene_fam(pangenome: Pangenome, output: str, gene_families: str, soft_core: float = 0.95,
                         compress: bool = False, disable_bar=False):
    """
    Write representative nucleotide sequences of gene families

    :param pangenome: Pangenome object with gene families sequences
    :param output: Path to output directory
    :param gene_families: Selected partition of gene families
    :param soft_core: Soft core threshold to use
    :param compress: Compress the file in .gz
    :param disable_bar: Disable progress bar
    """
    outname = output + f"/{gene_families}_nucleotide_families.fasta"

    genefams = select_families(pangenome, gene_families, "representative nucleotide sequences of the gene families",
                               soft_core)

    with write_compressed_or_not(outname, compress) as fasta:
        get_gene_sequences_from_file(pangenome.file, fasta, [fam.name for fam in genefams], disable_bar=disable_bar)

    logging.getLogger().info(f"Done writing the representative nucleotide sequences of the gene families : '{outname}'")


def write_fasta_prot_fam(pangenome: Pangenome, output: str, prot_families: str, soft_core: float = 0.95,
                         compress: bool = False, disable_bar: bool = False):
    """
    Write representative amino acid sequences of gene families.

    :param pangenome: Pangenome object with gene families sequences
    :param output: Path to output directory
    :param prot_families: Selected partition of protein families
    :param soft_core: Soft core threshold to use
    :param compress: Compress the file in .gz
    :param disable_bar: Disable progress bar
    """
    outname = output + f"/{prot_families}_protein_families.faa"

    genefams = select_families(pangenome, prot_families, "representative amino acid sequences of the gene families",
                               soft_core)

    with write_compressed_or_not(outname, compress) as fasta:
        for fam in tqdm(genefams, unit="prot families", disable=disable_bar):
            fasta.write('>' + fam.name + "\n")
            fasta.write(fam.sequence + "\n")
    logging.getLogger().info(f"Done writing the representative amino acid sequences of the gene families : '{outname}'")


def read_fasta_or_gff(filename: str) -> dict:
    """
    Read the genome file in fasta or gbff format

    :param filename: Path to genome file

    :return: Dictionary with all sequences associated to contig
    """
    sequence_dict = {}
    seqname = ""
    seq = ""
    z = False
    with read_compressed_or_not(filename) as f:
        for line in f:
            if line.startswith(">"):
                z = True
            if z:
                if line.startswith('>'):
                    if seq != "":
                        sequence_dict[seqname] = seq
                    seqname = line[1:].strip().split()[0]
                else:
                    seq += line.strip()
        if seq != "":
            sequence_dict[seqname] = seq
    return sequence_dict


def read_fasta_gbk(filename):
    """
    Read the genome file in gbk format

    :param filename: Path to genome file

    :return: Dictionary with all sequences associated to contig
    """
    # line.startswith("ORIGIN"):
    sequence_dict = {}
    lines = read_compressed_or_not(filename).readlines()[::-1]
    contig_id, contig_locus_id = ("", "")
    while len(lines) != 0:
        line = lines.pop()
        # beginning of contig
        if line.startswith('LOCUS'):
            contig_locus_id = line.split()[1]
            # If contig_id is not specified in VERSION afterwards like with Prokka,
            # in that case we use the one in LOCUS.
            while not line.startswith('FEATURES'):
                if line.startswith('VERSION'):
                    contig_id = line[12:].strip()
                line = lines.pop()
        if contig_id == "":
            contig_id = contig_locus_id
        while not line.startswith("ORIGIN"):
            line = lines.pop()  # stuff
        line = lines.pop()  # first sequence line.
        sequence = ""
        while not line.startswith('//'):
            sequence += line[10:].replace(" ", "").strip().upper()
            line = lines.pop()
        # get each gene's sequence.
        sequence_dict[contig_id] = sequence
        # end of contig
    return sequence_dict


def read_genome_file(file_dict: dict, genome_name: str) -> dict:
    """
    Read the genome file associated to organism

    :param file_dict: Dictionary given association between organism and fasta file
    :param genome_name: organism name

    :return: Dictionary with all sequences associated to contig
    """
    filetype = detect_filetype(file_dict[genome_name])
    if filetype in ["fasta", "gff"]:
        return read_fasta_or_gff(file_dict[genome_name])
    elif filetype == "gbff":
        return read_fasta_gbk(file_dict[genome_name])
    else:
        raise Exception(f"Unknown filetype detected: '{file_dict[genome_name]}'")


def write_spaced_fasta(sequence: str, space: int = 60):
    """Write a maximum of element per line

    :param sequence: sequence to write
    :param space: maximum of size for one line

    :return: a sequence of maximum space caracter
    """
    seq = ""
    j = 0
    while j < len(sequence):
        seq += sequence[j:j + space] + "\n"
        j += space
    return seq


def write_regions_sequences(pangenome: Pangenome, output: str, regions: str, fasta: str, anno: str,
                            compress: bool = False, disable_bar: bool = False):
    """
    Write representative amino acid sequences of gene families.

    :param pangenome: Pangenome object with gene families sequences
    :param output: Path to output directory
    :param regions: Write the RGP nucleotide sequences
    :param fasta: A tab-separated file listing the organism names, fasta filepath of its genomic sequences
    :param anno: A tab-separated file listing the organism names, and the gff/gbff filepath of its annotations
    :param compress: Compress the file in .gz
    :param disable_bar: Disable progress bar
    """
    organisms_file = fasta if fasta is not None else anno
    org_dict = {}
    for line in read_compressed_or_not(organisms_file):
        elements = [el.strip() for el in line.split("\t")]
        if len(elements) <= 1:
            raise Exception(f"No tabulation separator found in given --fasta or --anno file: '{organisms_file}'")
        org_dict[elements[0]] = elements[1]

    logging.getLogger().info(f"Writing {regions} rgp genomic sequences...")
    regions_to_write = []
    if regions == "complete":
        for region in pangenome.regions:
            if not region.is_contig_border:
                regions_to_write.append(region)
    else:
        regions_to_write = pangenome.regions

    regions_to_write = sorted(regions_to_write, key=lambda x: x.organism.name)
    # order regions by organism, so that we only have to read one genome at the time

    outname = output + f"/{regions}_rgp_genomic_sequences.fasta"
    with write_compressed_or_not(outname, compress) as fasta:
        loaded_genome = ""
        for region in tqdm(regions_to_write, unit="rgp", disable=disable_bar):
            if region.organism.name != loaded_genome:
                loaded_genome = region.organism.name
                genome_sequence = read_genome_file(org_dict, loaded_genome)
            fasta.write(f">{region.name}\n")
            fasta.write(write_spaced_fasta(genome_sequence[region.contig.name][region.start:region.stop], 60))
    logging.getLogger().info(f"Done writing the regions nucleotide sequences: '{outname}'")


def write_sequence_files(pangenome: Pangenome, output: str, fasta: str = None, anno: str = None,
                         soft_core: float = 0.95, regions: str = None, genes: str = None, gene_families: str = None,
                         prot_families: str = None, compress: bool = False, disable_bar: bool = False):
    """
    Main function to write sequence file from pangenome

    :param pangenome: Pangenome object containing sequences
    :param output: Path to output directory
    :param fasta: A tab-separated file listing the organism names, fasta filepath of its genomic sequences
    :param anno: A tab-separated file listing the organism names, and the gff/gbff filepath of its annotations
    :param soft_core: Soft core threshold to use
    :param regions: Write the RGP nucleotide sequences
    :param genes: Write all nucleotide CDS sequences
    :param gene_families: Write representative nucleotide sequences of gene families.
    :param prot_families: Write representative amino acid sequences of gene families.
    :param compress: Compress the file in .gz
    :param disable_bar: Disable progress bar
    """
    if not any(x for x in [regions, genes, prot_families, gene_families]):
        raise Exception("You did not indicate what file you wanted to write.")

    need_annotations = False
    need_families = False
    need_graph = False
    need_partitions = False
    need_spots = False
    need_regions = False
    need_modules = False

    if any(x is not None for x in [regions, genes, gene_families, prot_families]):
        need_annotations = True
        need_families = True
    if regions is not None or any(x == "rgp" for x in (genes, gene_families, prot_families)):
        need_regions = True
    if any(x in ["persistent", "shell", "cloud"] for x in (genes, gene_families, prot_families)):
        need_partitions = True
    for x in (genes, gene_families, prot_families):
        if x is not None and 'module_' in x:
            need_modules = True

    if not (need_annotations or need_families or need_graph or need_partitions or
            need_spots or need_regions or need_modules):
        # then nothing is needed, then something is wrong.
        # find which filter was provided
        provided_filter = ''
        if genes is not None:
            provided_filter = genes
        if gene_families is not None:
            provided_filter = gene_families
        if prot_families is not None:
            provided_filter = prot_families
        if regions is not None:
            provided_filter = regions
        raise Exception(f"The filter that you indicated '{provided_filter}' was not understood by PPanGGOLiN. "
                        f"{poss_values_log}")
    ex_gene_sequences = Exception("The provided pangenome has no gene sequences. "
                                  "This is not compatible with any of the following options : --genes, --gene_families")
    ex_gene_family_sequences = Exception("The provided pangenome has no gene families. "
                                         "This is not compatible with any of the following options : "
                                         "--prot_families, --gene_families")
    if not pangenome.status["geneSequences"] in ["inFile"] and (genes or gene_families):
        raise ex_gene_sequences
    if not pangenome.status["geneFamilySequences"] in ["Loaded", "Computed", "inFile"] and prot_families:
        raise ex_gene_family_sequences

    check_pangenome_info(pangenome, need_annotations=need_annotations, need_families=need_families,
                         need_graph=need_graph, need_partitions=need_partitions, need_rgp=need_regions,
                         need_spots=need_spots, need_modules=need_modules, disable_bar=disable_bar)

    if prot_families is not None:
        write_fasta_prot_fam(pangenome, output, prot_families, soft_core, compress, disable_bar)
    if gene_families is not None:
        write_fasta_gene_fam(pangenome, output, gene_families, soft_core, compress, disable_bar)
    if genes is not None:
        write_gene_sequences(pangenome, output, genes, soft_core, compress, disable_bar)
    if regions is not None:
        write_regions_sequences(pangenome, output, regions, fasta, anno, compress, disable_bar)


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    if args.regions is not None and args.fasta is None and args.anno is None:
        raise Exception("The --regions options requires the use of --anno or --fasta "
                        "(You need to provide the same file used to compute the pan)")
    mk_outdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    write_sequence_files(pangenome, args.output, fasta=args.fasta, anno=args.anno, soft_core=args.soft_core,
                         regions=args.regions, genes=args.genes, gene_families=args.gene_families,
                         prot_families=args.prot_families, compress=args.compress, disable_bar=args.disable_prog_bar)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("fasta", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_seq(parser)
    return parser


def parser_seq(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of fasta command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o', '--output', required=True, type=str,
                          help="Output directory where the file(s) will be written")

    context = parser.add_argument_group(title="Contextually required arguments",
                                        description="With --regions, the following arguments are required:")
    context.add_argument('--fasta', required=False, type=str,
                         help="A tab-separated file listing the organism names, and the fasta filepath of its genomic "
                              "sequence(s) (the fastas can be compressed with gzip). One line per organism.")
    context.add_argument('--anno', required=False, type=str,
                         help="A tab-separated file listing the organism names, and the gff/gbff filepath of its "
                              "annotations (the files can be compressed with gzip). One line per organism. "
                              "If this is provided, those annotations will be used.")
    onereq = parser.add_argument_group(title="Output file",
                                       description="At least one of the following argument is required. "
                                                   "Indicating 'all' writes all elements. Writing a partition "
                                                   "('persistent', 'shell' or 'cloud') write the elements associated "
                                                   "to said partition. Writing 'rgp' writes elements associated to RGPs"
                                       )
    onereq.add_argument("--genes", required=False, type=str,
                        help=f"Write all nucleotide CDS sequences. {poss_values_log}")
    onereq.add_argument("--prot_families", required=False, type=str,
                        help=f"Write representative amino acid sequences of gene families. {poss_values_log}")
    onereq.add_argument("--gene_families", required=False, type=str,
                        help=f"Write representative nucleotide sequences of gene families. {poss_values_log}")
    optional = parser.add_argument_group(title="Optional arguments")
    # could make choice to allow customization
    optional.add_argument("--regions", required=False, type=str, choices=["all", "complete"],
                          help="Write the RGP nucleotide sequences (requires --anno or --fasta used to compute "
                               "the pangenome to be given)")
    optional.add_argument("--soft_core", required=False, type=restricted_float, default=0.95,
                          help="Soft core threshold to use if 'softcore' partition is chosen")
    optional.add_argument("--compress", required=False, action="store_true", help="Compress the files in .gz")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_seq(main_parser)
    common = main_parser.add_argument_group(title="Common argument")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
    common.add_argument('-f', '--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
