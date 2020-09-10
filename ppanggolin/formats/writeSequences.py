#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse
from multiprocessing import Pool
from collections import Counter, defaultdict
import logging
import pkg_resources
from statistics import median, mean, stdev
import os
from random import randint, shuffle

#installed libraries
from tqdm import tqdm

#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import write_compressed_or_not, mkOutdir, read_compressed_or_not
from ppanggolin.formats import checkPangenomeInfo, getGeneSequencesFromFile
from ppanggolin.annotate import detect_filetype

def writeGeneSequencesFromAnnotations(pangenome, fileObj, list_CDS=None, verbose = False):
    """
        Writes the CDS sequences of the Pangenome object to a tmpFile object
        Loads the sequences from previously computed or loaded annotations
    """
    if list_CDS is None:
        list_CDS = pangenome.genes
    logging.getLogger().info("Writing all of the CDS sequences...")
    bar =  tqdm(list_CDS, unit="gene", disable= not verbose)
    for gene in bar:
        if gene.type == "CDS":
            fileObj.write('>' + gene.ID + "\n")
            fileObj.write(gene.dna + "\n")
    fileObj.flush()
    bar.close()

def writeGeneSequences(pangenome, output, compress, genes):
    logging.getLogger().info("Writing all the gene nucleic sequences...")
    outname = output + f"/{genes}_genes.fna"

    genes_to_write = []
    if genes == 'all':
        logging.getLogger().info("Writing all of the gene sequences...")
        genes_to_write = pangenome.genes
    if genes in ['persistent','shell','cloud']:
        logging.getLogger().info(f"Writing all of the {genes} gene sequences...")
        for gene in pangenome.genes:
            if gene.family.namedPartition == genes:
                genes_to_write.append(gene)
    if genes == "rgp":
        logging.getLogger().info(f"Writing all of the gene sequences in RGP...")
        for region in pangenome.regions:
            genes_to_write.extend(region.genes)
    logging.getLogger().info(f"There are {len(genes_to_write)} genes to write")
    with write_compressed_or_not(outname,compress) as fasta:
        if pangenome.status["geneSequences"] in ["inFile"]:
            getGeneSequencesFromFile(pangenome,fasta, set([gene.ID for gene in genes_to_write]))
        elif pangenome.status["geneSequences"] in ["Computed","Loaded"]:
            writeGeneSequencesFromAnnotations(pangenome, fasta, genes_to_write)
        else:
            #this should never happen if the pangenome has been properly checked before launching this function.
            raise Exception("The pangenome does not include gene sequences")
    logging.getLogger().info(f"Done writing the gene sequences : '{outname}'")

def writeFastaGeneFam(pangenome, output, compress, gene_families):
    outname = output + f"/{gene_families}_nucleotide_families.fasta"

    genefams = set()
    if gene_families == 'all':
        logging.getLogger().info("Writing all of the representative nucleotide sequences of the gene families...")
        genefams = pangenome.geneFamilies
    if gene_families in ['persistent','shell','cloud']:
        logging.getLogger().info(f"Writing the representative nucleotide sequences of the {gene_families} gene families...")
        for fam in pangenome.geneFamilies:
            if fam.namedPartition == gene_families:
                genefams.add(fam)
    if gene_families == "rgp":
        logging.getLogger().info(f"Writing the representative nucleotide sequences of the gene families in RGPs...")
        for region in pangenome.regions:
            genefams |= region.families

    with write_compressed_or_not(outname,compress) as fasta:
        getGeneSequencesFromFile(pangenome,fasta,[fam.name for fam in genefams])

    logging.getLogger().info(f"Done writing the representative nucleotide sequences of the gene families : '{outname}'")

def writeFastaProtFam(pangenome, output, compress, prot_families):
    outname = output + f"/{prot_families}_protein_families.faa"

    genefams = set()
    if prot_families == 'all':
        logging.getLogger().info("Writing all of the representative amino acid sequences of the gene families...")
        genefams = pangenome.geneFamilies
    if prot_families in ['persistent','shell','cloud']:
        logging.getLogger().info(f"Writing the representative amino acid sequences of the {prot_families} gene families...")
        for fam in pangenome.geneFamilies:
            if fam.namedPartition == prot_families:
                genefams.add(fam)
    if prot_families == "rgp":
        logging.getLogger().info(f"Writing the representative amino acid sequences of the gene families in RGPs...")
        for region in pangenome.regions:
            genefams |= region.families

    with write_compressed_or_not(outname,compress) as fasta:
        bar = tqdm(genefams,unit="prot families")
        for fam in bar:
            fasta.write('>' +fam.name + "\n")
            fasta.write(fam.sequence + "\n")
        bar.close()
    logging.getLogger().info(f"Done writing the representative amino acid sequences of the gene families : '{outname}'")

def read_fasta_or_gff(filename):
    sequence_dict = {}
    line = ""
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
    #line.startswith("ORIGIN"):
    sequence_dict = {}
    line = ""
    lines = read_compressed_or_not(filename).readlines()[::-1]
    while len(lines) != 0:
        line = lines.pop()
        # beginning of contig
        if line.startswith('LOCUS'):
            contigLocusID = line.split()[1]#If contigID is not specified in VERSION afterwards like with Prokka, in that case we use the one in LOCUS.
            while not line.startswith('FEATURES'):
                if line.startswith('VERSION'):
                    contigID = line[12:].strip()
                line = lines.pop()
        if contigID == "":
            contigID = contigLocusID
        while not line.startswith("ORIGIN"):
            line = lines.pop()#stuff
        line = lines.pop()#first sequence line.
        sequence = ""
        while not line.startswith('//'):
            sequence += line[10:].replace(" ", "").strip().upper()
            line = lines.pop()
        #get each gene's sequence.
        sequence_dict[contigID] = sequence
        #end of contig
    return sequence_dict

def read_genome_file(file_dict, genome_name):
    filetype = detect_filetype(file_dict[genome_name])
    if filetype in ["fasta","gff"]:
        return read_fasta_or_gff(file_dict[genome_name])
    elif filetype == "gbff":
        return read_fasta_gbk(file_dict[genome_name])
    else:
        raise Exception(f"Unknown filetype detected: '{file_dict[genome_name]}'")

def write_spaced_fasta(sequence, space):
    seq = ""
    j = 0
    while j < len(sequence):
        seq += sequence[j:j+space] + "\n"
        j+=space
    return seq

def writeRegionsSequences(pangenome, output, compress, regions, fasta, anno):
    organisms_file = fasta if fasta is not None else anno
    org_dict = {}
    for line in read_compressed_or_not(organisms_file):
        elements = [el.strip() for el in line.split("\t")]
        if len(elements)<=1:
            logging.getLogger().error(f"No tabulation separator found in given --fasta or --anno file: '{organisms_file}'")
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
    
    regions_to_write = sorted(regions_to_write, key=lambda  x: x.organism.name)#order regions by organism, so that we only have to read one genome at the time

    outname = output + f"/{regions}_rgp_genomic_sequences.fasta"
    with write_compressed_or_not(outname,compress) as fasta:
        loaded_genome = ""
        bar = tqdm(regions_to_write, unit = "rgp")
        for region in bar:
            if region.organism.name != loaded_genome:
                loaded_genome = region.organism.name
                genome_sequence = read_genome_file(org_dict, loaded_genome)
            fasta.write(f">{region.name}\n")
            fasta.write(write_spaced_fasta(genome_sequence[region.contig.name][region.start:region.stop], 60))
        bar.close()
    logging.getLogger().info(f"Done writing the regions nucleotide sequences: '{outname}'")

def writeSequenceFiles(pangenome, output, fasta=None, anno=None, cpu=1, regions=None, genes=None, gene_families=None, prot_families=None, compress=False):
    if not any(x for x in [ regions, genes, prot_families, gene_families]):
        raise Exception("You did not indicate what file you wanted to write.")

    needAnnotations = False
    needFamilies = False
    needGraph = False
    needPartitions = False
    needSpots = False
    needRegions = False

    if any(x is not None for x in [regions, genes, gene_families, prot_families]):
        needAnnotations=True
    if regions is not None or any(x == "rgp" for x in (genes, gene_families, prot_families)):
        needRegions= True
    if regions is not None or gene_families is not None or prot_families is not None or genes in ["persistent","shell","cloud"]:
        needFamilies=True
    if any(x in ["persistent","shell","cloud"] for x in (genes, gene_families, prot_families)):
        needPartitions=True

    #need to deal with sequence-related flags outside of checkPangenomeInfo since 
    ex_geneSequences = Exception("The provided pangenome has no gene sequences. This is not compatible with any of the following options : --genes, --gene_families")
    ex_geneFamilySequences = Exception("The provided pangenome has no gene families. This is not compatible with any of the following options : --prot_families, --gene_families")
    if not pangenome.status["geneSequences"] in ["inFile"] and (genes or gene_families):
        raise ex_geneSequences
    if not pangenome.status["geneFamilySequences"] in ["Loaded","Computed","inFile"] and prot_families:
        raise ex_geneFamilySequences
    
    checkPangenomeInfo(pangenome, needAnnotations=needAnnotations, needFamilies=needFamilies, needGraph=needGraph, needPartitions= needPartitions, needRGP = needRegions, needSpots = needSpots)

    if prot_families is not None:
        writeFastaProtFam(pangenome, output, compress, prot_families)
    if gene_families is not None:
        writeFastaGeneFam(pangenome, output, compress, gene_families)
    if genes is not None:
        writeGeneSequences(pangenome, output, compress, genes)
    if regions is not None:
        writeRegionsSequences(pangenome, output, compress, regions, fasta, anno)

def checkOptions(args):
    if hasattr(args,"regions"):
        if args.regions is not None:
            if args.fasta is None and args.anno is None:
                raise Exception("The --regions options requires the use of --anno or --fasta (You need to provide the same file used to compute the pangenome)")

def launchSequences(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    checkOptions(args)
    writeSequenceFiles(pangenome, args.output, fasta=args.fasta, anno=args.anno, cpu=args.cpu, regions=args.regions, genes=args.genes, prot_families=args.prot_families, gene_families= args.gene_families, compress=args.compress)

def writeSequenceSubparser(subparser):
    parser = subparser.add_parser("fasta", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o','--output', required=True, type=str, help="Output directory where the file(s) will be written")

    context = parser.add_argument_group(title = "Contextually required arguments", description = "With --regions, the following arguments are required:")
    context.add_argument('--fasta',  required=False, type=str, help="A tab-separated file listing the organism names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed with gzip). One line per organism.")
    context.add_argument('--anno', required=False, type=str, help="A tab-separated file listing the organism names, and the gff/gbff filepath of its annotations (the files can be compressed with gzip). One line per organism. If this is provided, those annotations will be used.")

    optional = parser.add_argument_group(title = "Optional arguments. Indicating 'all' writes all elements. Writing a partition ('persistent', 'shell' or 'cloud') write the elements associated to said partition. Writing 'rgp' writes elements associated to RGPs.")
    ##could make choice to allow customization
    optional.add_argument("--regions", required=False, choices=["all","complete"], help = "Write the RGP nucleotide sequences (requires --anno or --fasta used to compute the pangenome to be given)")
    optional.add_argument("--genes", required=False,choices=["all","persistent","shell","cloud","rgp"], help = "Write all nucleotide CDS sequences")
    optional.add_argument("--prot_families", required=False, choices=["all","persistent","shell","cloud","rgp"], help = "Write representative amino acid sequences of gene families")
    optional.add_argument("--gene_families", required=False, choices=["all","persistent","shell","cloud","rgp"], help = "Write representative nucleotide sequences of gene families")
    optional.add_argument("--compress",required=False, action="store_true",help="Compress the files in .gz")
    return parser
