#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging

# installed libraries
from tqdm import tqdm

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import write_compressed_or_not, mkOutdir, read_compressed_or_not, restricted_float
from ppanggolin.formats import checkPangenomeInfo, getGeneSequencesFromFile
from ppanggolin.annotate import detect_filetype

poss_values_log = "Possible values are 'all', 'persistent', 'shell', 'cloud', 'rgp', 'softcore', " \
                  "'core', 'module_X' with X being a module id."


def writeGeneSequencesFromAnnotations(pangenome, fileObj, list_CDS=None, disable_bar=False):
    """
    Writes the CDS sequences of the Pangenome object to a tmpFile object
    Loads the sequences from previously computed or loaded annotations
    """
    if list_CDS is None:
        list_CDS = pangenome.genes
    logging.getLogger().info("Writing all of the CDS sequences...")
    for gene in tqdm(list_CDS, unit="gene", disable=disable_bar):
        if gene.type == "CDS":
            fileObj.write('>' + gene.ID + "\n")
            fileObj.write(gene.dna + "\n")
    fileObj.flush()


def writeGeneSequences(pangenome, output, compress, genes, soft_core=0.95, disable_bar=False):
    logging.getLogger().info("Writing all the gene nucleotide sequences...")
    outname = output + f"/{genes}_genes.fna"

    genefams = selectFamilies(pangenome, genes, "gene nucleotide sequences", soft_core)
    genes_to_write = []

    for fam in genefams:
        genes_to_write.extend(fam.genes)

    logging.getLogger().info(f"There are {len(genes_to_write)} genes to write")
    with write_compressed_or_not(outname, compress) as fasta:
        if pangenome.status["geneSequences"] in ["inFile"]:
            getGeneSequencesFromFile(pangenome.file, fasta, set([gene.ID for gene in genes_to_write]),
                                     disable_bar=disable_bar)
        elif pangenome.status["geneSequences"] in ["Computed", "Loaded"]:
            writeGeneSequencesFromAnnotations(pangenome, fasta, genes_to_write, disable_bar=disable_bar)
        else:
            # this should never happen if the pangenome has been properly checked before launching this function.
            raise Exception("The pangenome does not include gene sequences")
    logging.getLogger().info(f"Done writing the gene sequences : '{outname}'")


def selectFamilies(pangenome, partition, type_name, soft_core):
    """ function used to filter down families to the given partition

    """
    genefams = set()
    if partition == 'all':
        logging.getLogger().info(f"Writing all of the {type_name}...")
        genefams = pangenome.geneFamilies
    elif partition in ['persistent', 'shell', 'cloud']:
        logging.getLogger().info(f"Writing the {type_name} of the {partition}...")
        for fam in pangenome.geneFamilies:
            if fam.namedPartition == partition:
                genefams.add(fam)
    elif partition == "rgp":
        logging.getLogger().info(f"Writing the {type_name} in RGPs...")
        for region in pangenome.regions:
            genefams |= region.families
    elif partition == "softcore":
        logging.getLogger().info(
            f"Writing the {type_name} in {partition} genome, that are present in more than {soft_core} of genomes")
        threshold = pangenome.number_of_organisms() * soft_core
        for fam in pangenome.geneFamilies:
            if len(fam.organisms) >= threshold:
                genefams.add(fam)
    elif partition == "core":
        logging.getLogger().info(f"Writing the representative {type_name} of the {partition} gene families...")
        for fam in pangenome.geneFamilies:
            if len(fam.organisms) == pangenome.number_of_organisms():
                genefams.add(fam)
    elif "module_" in partition:
        logging.getLogger().info(f"Writing the representation {type_name} of {partition} gene families...")
        modID = int(partition.replace("module_", ""))
        for mod in pangenome.modules:
            # could be way more efficient with a dict structure instead of a set
            if mod.ID == modID:
                genefams |= mod.families
                break
    return genefams


def writeFastaGeneFam(pangenome, output, compress, gene_families, soft_core=0.95, disable_bar=False):
    outname = output + f"/{gene_families}_nucleotide_families.fasta"

    genefams = selectFamilies(pangenome, gene_families, "representative nucleotide sequences of the gene families",
                              soft_core)

    with write_compressed_or_not(outname, compress) as fasta:
        getGeneSequencesFromFile(pangenome.file, fasta, [fam.name for fam in genefams], disable_bar=disable_bar)

    logging.getLogger().info(f"Done writing the representative nucleotide sequences of the gene families : '{outname}'")


def writeFastaProtFam(pangenome, output, compress, prot_families, soft_core=0.95, disable_bar=False):
    outname = output + f"/{prot_families}_protein_families.faa"

    genefams = selectFamilies(pangenome, prot_families, "representative amino acid sequences of the gene families",
                              soft_core)

    with write_compressed_or_not(outname, compress) as fasta:
        bar = tqdm(genefams, unit="prot families", disable=disable_bar)
        for fam in bar:
            fasta.write('>' + fam.name + "\n")
            fasta.write(fam.sequence + "\n")
        bar.close()
    logging.getLogger().info(f"Done writing the representative amino acid sequences of the gene families : '{outname}'")


def read_fasta_or_gff(filename):
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
    # line.startswith("ORIGIN"):
    sequence_dict = {}
    lines = read_compressed_or_not(filename).readlines()[::-1]
    while len(lines) != 0:
        line = lines.pop()
        # beginning of contig
        if line.startswith('LOCUS'):
            contigLocusID = line.split()[1]
            # If contigID is not specified in VERSION afterwards like with Prokka, in that case we use the one in LOCUS.
            while not line.startswith('FEATURES'):
                if line.startswith('VERSION'):
                    contigID = line[12:].strip()
                line = lines.pop()
        if contigID == "":
            contigID = contigLocusID
        while not line.startswith("ORIGIN"):
            line = lines.pop()  # stuff
        line = lines.pop()  # first sequence line.
        sequence = ""
        while not line.startswith('//'):
            sequence += line[10:].replace(" ", "").strip().upper()
            line = lines.pop()
        # get each gene's sequence.
        sequence_dict[contigID] = sequence
        # end of contig
    return sequence_dict


def read_genome_file(file_dict, genome_name):
    filetype = detect_filetype(file_dict[genome_name])
    if filetype in ["fasta", "gff"]:
        return read_fasta_or_gff(file_dict[genome_name])
    elif filetype == "gbff":
        return read_fasta_gbk(file_dict[genome_name])
    else:
        raise Exception(f"Unknown filetype detected: '{file_dict[genome_name]}'")


def write_spaced_fasta(sequence, space):
    seq = ""
    j = 0
    while j < len(sequence):
        seq += sequence[j:j + space] + "\n"
        j += space
    return seq


def writeRegionsSequences(pangenome, output, compress, regions, fasta, anno, disable_bar=False):
    organisms_file = fasta if fasta is not None else anno
    org_dict = {}
    for line in read_compressed_or_not(organisms_file):
        elements = [el.strip() for el in line.split("\t")]
        if len(elements) <= 1:
            logging.getLogger().error(
                f"No tabulation separator found in given --fasta or --anno file: '{organisms_file}'")
            exit(1)
        org_dict[elements[0]] = elements[1]

    logging.getLogger().info(f"Writing {regions} rgp genomic sequences...")
    regions_to_write = []
    if regions == "complete":
        for region in pangenome.regions:
            if not region.isContigBorder:
                regions_to_write.append(region)
    else:
        regions_to_write = pangenome.regions

    regions_to_write = sorted(regions_to_write, key=lambda x: x.organism.name)
    # order regions by organism, so that we only have to read one genome at the time

    outname = output + f"/{regions}_rgp_genomic_sequences.fasta"
    with write_compressed_or_not(outname, compress) as fasta:
        loaded_genome = ""
        bar = tqdm(regions_to_write, unit="rgp", disable=disable_bar)
        for region in bar:
            if region.organism.name != loaded_genome:
                loaded_genome = region.organism.name
                genome_sequence = read_genome_file(org_dict, loaded_genome)
            fasta.write(f">{region.name}\n")
            fasta.write(write_spaced_fasta(genome_sequence[region.contig.name][region.start:region.stop], 60))
        bar.close()
    logging.getLogger().info(f"Done writing the regions nucleotide sequences: '{outname}'")


def writeSequenceFiles(pangenome, output, fasta=None, anno=None, soft_core=0.95, regions=None, genes=None,
                       gene_families=None, prot_families=None, compress=False, disable_bar=False):
    if not any(x for x in [regions, genes, prot_families, gene_families]):
        raise Exception("You did not indicate what file you wanted to write.")

    needAnnotations = False
    needFamilies = False
    needGraph = False
    needPartitions = False
    needSpots = False
    needRegions = False
    needModules = False

    if any(x is not None for x in [regions, genes, gene_families, prot_families]):
        needAnnotations = True
        needFamilies = True
    if regions is not None or any(x == "rgp" for x in (genes, gene_families, prot_families)):
        needRegions = True
    if any(x in ["persistent", "shell", "cloud"] for x in (genes, gene_families, prot_families)):
        needPartitions = True
    for x in (genes, gene_families, prot_families):
        if x is not None and 'module_' in x:
            needModules = True

    if not (needAnnotations or needFamilies or needGraph or needPartitions or needSpots or needRegions or needModules):
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
        raise Exception(
            f"The filter that you indicated '{provided_filter}' was not understood by PPanGGOLiN. {poss_values_log}")
    ex_geneSequences = Exception("The provided pangenome has no gene sequences. "
                                 "This is not compatible with any of the following options : --genes, --gene_families")
    ex_geneFamilySequences = Exception("The provided pangenome has no gene families. "
                                       "This is not compatible with any of the following options : "
                                       "--prot_families, --gene_families")
    if not pangenome.status["geneSequences"] in ["inFile"] and (genes or gene_families):
        raise ex_geneSequences
    if not pangenome.status["geneFamilySequences"] in ["Loaded", "Computed", "inFile"] and prot_families:
        raise ex_geneFamilySequences

    checkPangenomeInfo(pangenome, needAnnotations=needAnnotations, needFamilies=needFamilies, needGraph=needGraph,
                       needPartitions=needPartitions, needRGP=needRegions, needSpots=needSpots, needModules=needModules,
                       disable_bar=disable_bar)

    if prot_families is not None:
        writeFastaProtFam(pangenome, output, compress, prot_families, soft_core=soft_core, disable_bar=disable_bar)
    if gene_families is not None:
        writeFastaGeneFam(pangenome, output, compress, gene_families, soft_core=soft_core, disable_bar=disable_bar)
    if genes is not None:
        writeGeneSequences(pangenome, output, compress, genes, soft_core=soft_core, disable_bar=disable_bar)
    if regions is not None:
        writeRegionsSequences(pangenome, output, compress, regions, fasta, anno, disable_bar=disable_bar)


def checkOptions(args):
    if hasattr(args, "regions") and args.regions is not None and args.fasta is None and args.anno is None:
        raise Exception("The --regions options requires the use of --anno or --fasta "
                        "(You need to provide the same file used to compute the pangenome)")


def launchSequences(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    checkOptions(args)
    writeSequenceFiles(pangenome, args.output, fasta=args.fasta, anno=args.anno, soft_core=args.soft_core,
                       regions=args.regions, genes=args.genes, gene_families=args.gene_families,
                       prot_families=args.prot_families, compress=args.compress, disable_bar=args.disable_prog_bar)


def writeSequenceSubparser(subparser):
    parser = subparser.add_parser("fasta", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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

    optional = parser.add_argument_group(
        title="Optional arguments. Indicating 'all' writes all elements. Writing a partition "
              "('persistent', 'shell' or 'cloud') write the elements associated to said partition. "
              "Writing 'rgp' writes elements associated to RGPs.")
    # could make choice to allow customization
    optional.add_argument("--regions", required=False, choices=["all", "complete"],
                          help="Write the RGP nucleotide sequences (requires --anno or --fasta used to compute "
                               "the pangenome to be given)")
    optional.add_argument("--genes", required=False, help=f"Write all nucleotide CDS sequences. {poss_values_log}")
    optional.add_argument("--prot_families", required=False,
                          help=f"Write representative amino acid sequences of gene families. {poss_values_log}")
    optional.add_argument("--gene_families", required=False,
                          help=f"Write representative nucleotide sequences of gene families. {poss_values_log}")
    optional.add_argument("--soft_core", required=False, type=restricted_float, default=0.95,
                          help="Soft core threshold to use if 'softcore' partition is chosen")
    optional.add_argument("--compress", required=False, action="store_true", help="Compress the files in .gz")
    return parser
