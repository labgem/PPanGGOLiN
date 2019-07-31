#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse
from multiprocessing import Pool
import logging
from pathlib import Path
import os
import time

#installed libraries
from tqdm import tqdm

#local libraries
from ppanggolin.annotate import  annotate_organism, read_fasta, get_dna_sequence
from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Organism, Gene, RNA
from ppanggolin.utils import read_compressed_or_not, getCurrentRAM, mkFilename, get_num_lines
from ppanggolin.formats import writePangenome, readPangenome

def read_org_line(pangenome, organism, gff_file_path, circular_contigs):
    (GFF_seqname, _, GFF_type, GFF_start, GFF_end, _, GFF_strand, _, GFF_attribute) = range(0,9)#missing values : source, score, frame. They are unused.
    def getGffAttributes(gff_fields):
        """
            Parses the gff attribute's line and outputs the attributes in a dict structure.
            :param gff_fields: a gff line stored as a list. Each element of the list is a column of the gff.
            :type list:
            :return: attributes:
            :rtype: dict
        """
        attributes_field = [f for f in gff_fields[GFF_attribute].strip().split(';') if len(f)>0]
        attributes = {}
        for att in attributes_field:
            (key, value) = att.strip().split('=')
            attributes[key.upper()]=value
        return attributes

    def getIDAttribute(attributes):
        """
            Gets the ID of the element from which the provided attributes were extracted. Raises an error if no ID is found.
            :param attribute:
            :type dict:
            :return: ElementID:
            :rtype: string
        """
        ElementID = attributes.get("ID")
        if not ElementID:
            logging.getLogger().error("Each CDS type of the gff files must own a unique ID attribute. Not the case for file: "+gff_file_path)
            exit(1)
        return ElementID

    hasFasta = False
    fastaString = ""
    org = Organism(organism)
    with read_compressed_or_not(gff_file_path) as gff_file:
        for line in gff_file:
            if hasFasta:
                fastaString += line
                continue
            elif line.startswith('##',0,2):
                if line.startswith('FASTA',2,7):
                    hasFasta = True
                elif line.startswith('sequence-region',2,17):
                    fields = [el.strip() for el in line.split()]
                    contig = org.addContig(fields[1], True if fields[1] in circular_contigs else False)
                continue
            elif line.startswith('#!',0,2):## special refseq comment lines for versionning softs, assemblies and annotations.
                continue
            gff_fields = [el.strip() for el in line.split('\t')]
            attributes = getGffAttributes(gff_fields)

            if gff_fields[GFF_type] == 'region':
                if gff_fields[GFF_seqname] in circular_contigs:
                    contig.is_circular = True
            elif gff_fields[GFF_type] == 'CDS' or "RNA" in gff_fields[GFF_type]:
                geneID = attributes.get("PROTEIN_ID")#if there is a 'PROTEIN_ID' attribute, it's where the ncbi stores the actual gene ids, so we use that.
                if geneID is None:#if its not found, we get the one under the 'ID' field which must exist (otherwise not a gff3 compliant file)
                    geneID = getIDAttribute(attributes)
                try:
                    name = attributes.pop('NAME')
                except KeyError:
                    try:
                        name = attributes.pop('GENE')
                    except KeyError:
                        name = ""

                try:
                    product = attributes.pop('PRODUCT')
                except KeyError:
                    product = ""

                try:
                    genetic_code = attributes.pop("TRANSL_TABLE")
                except KeyError:
                    genetic_code = "11"
                if contig.name != gff_fields[GFF_seqname]:
                    contig = org.addContig(gff_fields[GFF_seqname])#get the current contig
                if gff_fields[GFF_type] == "CDS":
                    gene = Gene(geneID)
                    #here contig is filled in order, so position is the number of genes already stored in the contig.
                    gene.fill_annotations(start = int(gff_fields[GFF_start]),
                                        stop = int(gff_fields[GFF_end]),
                                        strand =gff_fields[GFF_strand],
                                        geneType = gff_fields[GFF_type],
                                        position = len(contig.genes),
                                        name = name,
                                        product=product,
                                        genetic_code = genetic_code )
                    gene.fill_parents(org, contig)
                    contig.addGene(gene)
                elif "RNA" in gff_fields[GFF_type]:
                    rna = RNA(geneID)
                    rna.fill_annotations(start = int(gff_fields[GFF_start]),
                                        stop = int(gff_fields[GFF_end]),
                                        strand =gff_fields[GFF_strand],
                                        geneType = gff_fields[GFF_type],
                                        name = name,
                                        product=product)
                    rna.fill_parents(organism, contig)
                    contig.addRNA(rna)

    ### GET THE FASTA SEQUENCES OF THE GENES
    if hasFasta and fastaString != "":
        contigSequences = read_fasta( org, fastaString.split('\n'))
        for contig in org.contigs:
            for gene in contig.genes:
                gene.add_dna(get_dna_sequence(contigSequences[contig.name], gene))
            for rna in contig.RNAs:
                rna.add_dna(get_dna_sequence(contigSequences[contig.name], rna))
        pangenome.status["geneSequences"] = "Computed"
    elif hasFasta and fastaString == "":
        raise Exception(f"You have a combination of gff file with and without fasta sequences at the end, which is not accepted (as it makes everything difficult to process). Please remove the fasta sequences from all gff files and provide the sequences in separate files through the --fasta option, or add the fasta sequences at the end of ALL your gff files. Error was raised with gff file : '{gff_file_path}'")

    pangenome.addOrganism(org)

def readAnnotations(pangenome, organisms_file):
    logging.getLogger().info("Reading "+organisms_file+" the list of organism files ...")
    bar = tqdm(read_compressed_or_not(organisms_file),total=get_num_lines(organisms_file), unit = "gff file")
    for line in bar:
        elements = [el.strip() for el in line.split("\t")]
        if len(elements)<=1:
            logging.getLogger().error("No tabulation separator found in organisms file")
            exit(1)
        bar.set_description("Processing "+elements[1].split("/")[-1])
        bar.refresh()
        
        read_org_line(pangenome, elements[0], elements[1], elements[2:])
    bar.close()
    pangenome.status["genomesAnnotated"] = "Computed"

def getGeneSequencesFromFastas(pangenome, fasta_file):
    fastaDict = {}
    for line in read_compressed_or_not(fasta_file):
        elements = [el.strip() for el in line.split("\t")]
        if len(elements)<=1:
            logging.getLogger().error("No tabulation separator found in organisms file")
            exit(1)
        org = pangenome.addOrganism(elements[0])
        with read_compressed_or_not(elements[1]) as currFastaFile:
            fastaDict[org] = read_fasta(currFastaFile)
    
    for org in pangenome.organisms:
        for contig in org.contigs:
            for gene in contig.genes:
                gene.add_dna(get_dna_sequence(fastaDict[org][contig.name], gene))
            for rna in contig.RNAs:
                rna.add_dna(get_dna_sequence(fastaDict[org][contig.name], gene))
    pangenome.status["geneSequences"] = "Computed"

def launchAnnotateOrganism(pack):
    return annotate_organism(*pack)

def annotatePangenome(pangenome, fastaList, tmpdir,cpu, translation_table="11", kingdom = "bacteria", norna=False,  overlap=True):
    logging.getLogger().info(f"Reading {fastaList} the list of organism files")
    
    arguments = []
    for line in read_compressed_or_not(fastaList):
        elements = [el.strip() for el in line.split("\t")]
        if len(elements)<=1:
            logging.getLogger().error("No tabulation separator found in organisms file")
            exit(1)
        arguments.append((elements[0], elements[1], elements[2:], translation_table, kingdom, norna, tmpdir, overlap))
    
    logging.getLogger().info(f"Annotating {len(arguments)} genomes using {cpu} cpus...")
    with Pool(processes = cpu) as p:
        bar = tqdm(range(len(arguments)), unit = "genome")
        for organism in p.imap_unordered(launchAnnotateOrganism, arguments):
            bar.update()
            pangenome.addOrganism(organism)
    bar.close()
    logging.getLogger().info("Done annotating genomes")
    pangenome.status["genomesAnnotated"] = "Computed"#the pangenome is now annotated.
    pangenome.status["geneSequences"] = "Computed"#the gene objects have their respective gene sequences.
   
    return pangenome

def launch(args):
    filename = mkFilename(args.basename, args.output, args.force)
    pangenome = Pangenome()
    if args.fasta is not None:
        annotatePangenome(pangenome, args.fasta, args.tmpdir, args.cpu,  args.translation_table, args.kingdom, args.norna, args.overlap)
        
    elif args.gff is not None:
        readAnnotations(pangenome, args.gff)
        if pangenome.status["geneSequences"] == "No":
            if args.fasta:
                getGeneSequencesFromFastas(pangenome, args.fasta)
            else:
                logging.getLogger().warning("You provided gff files without sequences, and you did not provide fasta sequences. Thus it was not possible to get the gene sequences.")
                logging.getLogger().warning("You will be able to proceed with your analysis ONLY if you provide the clustering results in the next step.")
    writePangenome(pangenome, filename, args.force)

def syntaSubparser(subparser):
    parser = subparser.add_parser("annotate",help = "Annotate genomes")
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('-o','--output', required=False, type=str, default="ppanggolin_output"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory")
    optional.add_argument('--overlap', required=False, action='store_false',default=True, help="Use to not remove genes overlapping with RNA features.")
    optional.add_argument("--norna", required=False, action="store_true", default=False, help="Use to avoid annotating RNA features.")
    optional.add_argument("--kingdom",required = False, type = str.lower, default = "bacteria", choices = ["bacteria","archaea"], help = "Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation.")
    optional.add_argument("--translation_table",required=False, default="11", help = "Translation table (genetic code) to use.")
    optional.add_argument("--basename",required = False, default = "pangenome", help = "basename for the output file")
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('--fasta',  required=False, type=str, help="A tab-separated file listing the organism names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per organism.")
    required.add_argument('--gff', required=False, type=str, help="A tab-separated file listing the organism names, and the gff filepath of its annotations (the gffs can be compressed). One line per organism. If this is provided, those annotations will be used.")
    return parser