#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse
import os
import tempfile
from subprocess import Popen, PIPE
import ast
from collections import defaultdict

#local libraries
from ppanggolin.genome import Organism, Contig, Gene
from ppanggolin.utils import is_compressed, read_compressed_or_not
from ppanggolin.annotate import genetic_codes

def translate(seq, code):
    """ translates the given dna sequence with table code 11 of the ncbi (bacteria)"""
    # code:  https: //www.bioinformatics.org/sms/iupac.html
    start_table = code["start_table"]
    table = code["trans_table"]

    protein = ""
    if len(seq) % 3 == 0:
        protein = start_table[seq[0: 3]]
        for i in range(3, len(seq), 3):
            codon = seq[i: i + 3]
            try:
                protein += table[codon]
            except KeyError:  # codon was not planned for. Probably can't determine it.
                # print(codon)
                protein += 'X'  # X is for unknown
    else:
        print(len(seq))
        raise IndexError("Given sequence length modulo 3 was different than 0, which is unexpected.")
    return protein

def reverse_complement(seq):
    """ reverse complement the given dna sequence """
    complement = {'A':  'T', 'C':  'G', 'G':  'C', 'T':  'A', 'N': 'N', 'R': 'Y', 'Y': 'R',
                  'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}
    # see https://www.bioinformatics.org/sms/iupac.html for the code.
    # complement = {'A':  'T', 'C':  'G', 'G':  'C', 'T':  'A', 'N': 'N' } ## basic
    rcseq = ""
    for i in reversed(seq):
        rcseq += complement[i]
    return rcseq

def launch_aragorn(fnaFile, org):
    """ 
        launches Aragorn to annotate tRNAs. Takes a fna file name and a locustag to give an ID to the found genes.
        returns the annotated genes in a list of gene objects.
    """
    locustag = org.name
    cmd = ["aragorn", "-t", "-gcbact", "-l", "-w", fnaFile]
    p = Popen(cmd, stdout=PIPE)
    # loading the whole thing, reverting it to 'pop' in order.
    fileData = p.communicate()[0].decode().split("\n")[:: -1]
    geneObjs = defaultdict(set)
    c = 0
    while len(fileData) != 0:
        line = fileData.pop()
        if line.startswith(">"):
            header = line.replace(">", "").split()[0]
            fileData.pop()  # then next line must be removed too.
        elif len(line) > 0:  # if the line isn't empty, there's data to get.
            lineData = line.split()
            start, stop = ast.literal_eval(lineData[2].replace("c", ""))
            c += 1
            gene = Gene(ID = locustag+'_tRNA_'+str(c).zfill(3))
            gene.fill_annotations(start=start,
                                 stop=stop,
                                 strand="-" if lineData[2].startswith(
                                     "c") else "+",
                                 geneType="tRNA",
                                 product=lineData[1] + lineData[4])
            geneObjs[header].add(gene)
    return geneObjs

def launch_prodigal(fnaFile, org, code):
    """ 
        launches Prodigal to annotate CDS. Takes a fna file name and a locustag to give an ID to the found genes.
        returns the annotated genes in a list of gene objects.
    """
    locustag = org.name
    cmd = ["prodigal", "-f", "sco","-g",code, "-m", "-c", "-i", fnaFile, "-p", "single", "-q"]
    p = Popen(cmd, stdout=PIPE)

    geneObjs = defaultdict(set)
    c = 0
    for line in p.communicate()[0].decode().split("\n"):
        if line.startswith("# Sequence Data: "):
            for data in line.split(";"):
                if data.startswith("seqhdr"):
                    header = data.split("=")[1].replace('"', "").split()[0]
                    # print(header)

        elif line.startswith(">"):
            c += 1
            lineData = line[1:].split("_")  # not considering the '>'
            gene = Gene(ID = locustag + "_CDS_" + str(c).zfill(4))
            gene.fill_annotations(start=lineData[1],
                                 stop=lineData[2],
                                 strand=lineData[3],
                                 geneType="CDS",
                                 genetic_code=code)
            geneObjs[header].add(gene)

    return geneObjs

def launch_infernal(fnaFile, org, kingdom):
    """ 
        launches Infernal in hmmer-only mode to annotate rRNAs. Takes a fna file name and a locustag to give an ID to the found genes.
        returns the annotated genes in a list of gene objects.
    """
    locustag = org.name
    if kingdom == "bacteria":
        modelfile = os.path.dirname(os.path.realpath(__file__)) + "/rRNA_DB/rRNA_bact.cm"
    elif kingdom == "archaea":
        modelfile = os.path.dirname(os.path.realpath(__file__)) + "/rRNA_DB/rRNA_arch.cm"

    tmpFile = tempfile.NamedTemporaryFile(mode="r", dir = "/dev/shm/")
    cmd = ["cmscan", "--tblout", tmpFile.name, "--hmmonly", "--cpu",str(1), "--noali", modelfile, fnaFile]
    p = Popen(cmd, stdout=open(os.devnull, "w"), stderr=PIPE)
    err = p.communicate()[1].decode().split()
    if err != []:
        if err[0] == 'Error: ':
            raise Exception(f"Infernal (cmscan) failed with error:  '{ ' '.join(err) }'. If you never used this script, you should press the .cm file using cmpress executable from Infernal. You should find the file in '{os.path.dirname(os.path.realpath(__file__))}/rRNA_DB/'.")
        raise Exception(f"An error occurred with Infernal. Error is:  '{ ' '.join(err) }'.")
    # never managed to test what happens if the .cm files are compressed with a 'bad' version of infernal, so if that happens you are on your own.

    geneObjs = defaultdict(set)
    c = 0
    for line in tmpFile:
        if not line.startswith("#"):
            c += 1
            lineData = line.split()
            strand = lineData[9]
            if strand == "-":
                start = lineData[8]
                stop = lineData[7]
            else:
                start = lineData[7]
                stop = lineData[8]
            gene = Gene(ID = locustag + "_rRNA_" + str(c).zfill(3))
            gene.fill_annotations(start=start,
                                 stop=stop,
                                 strand=strand,
                                 geneType="rRNA",
                                 product=" ".join(lineData[17:]))
            geneObjs[lineData[2]].add(gene)

    return geneObjs

def read_fasta(org, fnaFile):
    """
        Reads a fna file and stores it in a dictionnary with contigs as key and sequence as value.
    """
    contigs = {}
    contig_seq = ""
    for line in fnaFile:
        if line.startswith('>'):
            if contig_seq != "":
                contigs[contig.name] = contig_seq
            contig_seq = ""
            contig = Contig(line.split()[0][1:])
            org.addContig(contig.name)
        else:
            contig_seq += line.strip()
    # processing the last contig
    if contig_seq != "":
        contigs[contig.name] = contig_seq
    return contigs

def write_tmp_fasta(contigs, tmpdir = "/dev/shm"):
    """
        Writes a temporary fna formated file, and returns the file-like object.

        This is for the cases where the given file is compressed, then we write a temporary file for the annotation tools to read from. The file will be deleted when close() is called.
    """
    tmpFile = tempfile.NamedTemporaryFile(mode="w", dir = tmpdir)
    for header in contigs.keys():
        tmpFile.write(f">{header}\n")
        j = 0
        while j < len(contigs[header]):
            tmpFile.write(contigs[header][j: j+60]+"\n")
            j += 60
    tmpFile.flush()  # force write what remains in the buffer.
    return tmpFile

def syntaxic_annotation(org, fastaFile, norna, kingdom, code):
    """
        Runs the different softwares for the syntaxic annotation.

        Takes in the file-like object containing the uncompressed fasta sequences to annotate
        the number of cpus that we can use.
        whether to annotate rna or not
        the locustag to give gene IDs.
    """
    # launching tools for syntaxic annotation
    genes = defaultdict(list)
    for key, items in launch_prodigal(fastaFile.name, org, code).items():
        genes[key].extend(items)
    if not norna:
        for key, items in launch_aragorn(fastaFile.name, org).items():
            genes[key].extend(items)
        for key, items in launch_infernal(fastaFile.name, org, kingdom).items():
            genes[key].extend(items)
    fastaFile.close()#closing either tmp file or original fasta file.
    return genes

def overlap_filter(allGenes, contigs, overlap):
    """
        Removes the CDS that overlap with RNA genes.
    """
    sortedGenes = defaultdict(set)
    for key, genes in allGenes.items():
        tmpGenes = sorted(genes, key=lambda x: x.start)
        rmGenes = set()
        if overlap:
            for i, gene_i in enumerate(tmpGenes):
                if i+1 < len(tmpGenes):
                    gene_j = tmpGenes[i+1]
                    if gene_i.type != "CDS" and gene_j.type == "CDS" and gene_i.stop > gene_j.start:
                        rmGenes.add(gene_j)
                    elif gene_i.type == "CDS" and gene_j.type != "CDS" and gene_i.stop > gene_j.start:
                        rmGenes.add(gene_i)

        for gene in rmGenes:
            tmpGenes.remove(gene)
        CDScounter = 0
        for gene in tmpGenes:
            if gene.type == "CDS":
                gene.position = CDScounter
                CDScounter+=1
        sortedGenes[key] = tmpGenes
    return sortedGenes

def get_protein_sequence(contigSeq, gene):
    if gene.strand == "+":
        return translate(contigSeq[gene.start-1:gene.stop], genetic_codes(gene.genetic_code))
    elif gene.strand == "-":
        return translate(reverse_complement(contigSeq[gene.start-1:gene.stop]), genetic_codes(gene.genetic_code))

def get_dna_sequence(contigSeq, gene):
    if gene.strand == "+":
        return contigSeq[gene.start-1:gene.stop]
    elif gene.strand == "-":
        return reverse_complement(contigSeq[gene.start-1:gene.stop])

def annotate_organism(orgName, fileName, circular_contigs, code, kingdom, norna, tmpdir, overlap):
    """
        Function to annotate a single organism
    """
    org = Organism(orgName)

    fastaFile = read_compressed_or_not(fileName)
    contigSequences = read_fasta(org, fastaFile)
    if is_compressed(fileName):
        fastaFile = write_tmp_fasta(contigSequences)

    genes = syntaxic_annotation(org, fastaFile, norna, kingdom, code)
    genes = overlap_filter(genes, contigSequences, overlap)

    for contigName, genes in genes.items():
        contig = org.addContig(contigName)
        if contig.name in circular_contigs:
            contig.is_circular = True
        for gene in genes:
            gene.add_dna(get_dna_sequence(contigSequences[contig.name], gene))
            if gene.type == "CDS":
                contig.addGene(gene)
            else:
                contig.addRNA(gene)
                
    return org

