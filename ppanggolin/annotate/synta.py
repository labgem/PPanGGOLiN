#!/usr/bin/env python3
#coding:utf-8

#default libraries
import os
import tempfile
from subprocess import Popen, PIPE
import ast
from collections import defaultdict

#local libraries
from ppanggolin.genome import Organism, Gene, RNA
from ppanggolin.utils import is_compressed, read_compressed_or_not

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
            gene = RNA(ID = locustag+'_tRNA_'+str(c).zfill(3))
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

def launch_infernal(fnaFile, org, kingdom, tmpdir):
    """
        launches Infernal in hmmer-only mode to annotate rRNAs. Takes a fna file name and a locustag to give an ID to the found genes.
        returns the annotated genes in a list of gene objects.
    """
    locustag = org.name
    if kingdom == "bacteria":
        modelfile = os.path.dirname(os.path.realpath(__file__)) + "/rRNA_DB/rRNA_bact.cm"
    elif kingdom == "archaea":
        modelfile = os.path.dirname(os.path.realpath(__file__)) + "/rRNA_DB/rRNA_arch.cm"

    tmpFile = tempfile.NamedTemporaryFile(mode="r", dir = tmpdir)
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
            gene = RNA(ID = locustag + "_rRNA_" + str(c).zfill(3))
            gene.fill_annotations(start=start,
                                 stop=stop,
                                 strand=strand,
                                 geneType="rRNA",
                                 product=" ".join(lineData[17:]))
            geneObjs[lineData[2]].add(gene)

    return geneObjs

def read_fasta(org, fnaFile, contig_filter):
    """
        Reads a fna file  (or stream, or string) and stores it in a dictionnary with contigs as key and sequence as value.
    """
    try:
        contigs = {}
        contig_seq = ""
        contig = None
        for line in fnaFile:
            if line.startswith('>'):
                if len(contig_seq) >= contig_filter:
                    contigs[contig.name] = contig_seq.upper()
                contig_seq = ""
                contig = org.getOrAddContig(line.split()[0][1:])
            else:
                contig_seq += line.strip()
        # processing the last contig
        if len(contig_seq) >= contig_filter:
            contigs[contig.name] = contig_seq.upper()
    except AttributeError as e:
        raise AttributeError(f"{e}\nAn error was raised when reading file: '{fnaFile.name}'. One possibility for this error is that the file did not start with a '>' as it would be expected from a fna file.")
    return contigs

def write_tmp_fasta(contigs, tmpdir ):
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

def syntaxic_annotation(org, fastaFile, norna, kingdom, code, tmpdir):
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
        for key, items in launch_infernal(fastaFile.name, org, kingdom, tmpdir).items():
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

def get_dna_sequence(contigSeq, gene):
    if gene.strand == "+":
        return contigSeq[gene.start-1:gene.stop]
    elif gene.strand == "-":
        return reverse_complement(contigSeq[gene.start-1:gene.stop])

def annotate_organism(orgName, fileName, circular_contigs, code, kingdom, norna, tmpdir, overlap, contig_filter):
    """
        Function to annotate a single organism
    """
    org = Organism(orgName)

    fastaFile = read_compressed_or_not(fileName)
    contigSequences = read_fasta(org, fastaFile, contig_filter)
    if is_compressed(fileName):
        fastaFile = write_tmp_fasta(contigSequences, tmpdir)

    genes = syntaxic_annotation(org, fastaFile, norna, kingdom, code, tmpdir)
    genes = overlap_filter(genes, contigSequences, overlap)

    for contigName, genes in genes.items():
        contig = org.getOrAddContig(contigName)
        if contig.name in circular_contigs:
            contig.is_circular = True
        for gene in genes:
            gene.add_dna(get_dna_sequence(contigSequences[contig.name], gene))
            gene.fill_parents(org, contig)
            if isinstance(gene, Gene):
                contig.addGene(gene)
            elif isinstance(gene, RNA):
                contig.addRNA(gene)
    return org
