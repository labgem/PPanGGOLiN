#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import os
import tempfile
from subprocess import Popen, PIPE
import ast
from collections import defaultdict

# local libraries
from ppanggolin.genome import Organism, Gene, RNA
from ppanggolin.utils import is_compressed, read_compressed_or_not


def reverse_complement(seq):
    """ reverse complement the given dna sequence """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'R': 'Y', 'Y': 'R',
                  'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}
    # see https://www.bioinformatics.org/sms/iupac.html for the code.
    # complement = {'A':  'T', 'C':  'G', 'G':  'C', 'T':  'A', 'N': 'N' } ## basic
    rcseq = ""
    for i in reversed(seq):
        rcseq += complement[i]
    return rcseq


def launch_aragorn(fna_file, org):
    """
        launches Aragorn to annotate tRNAs. Takes a fna file name and a locustag to give an ID to the found genes.
        returns the annotated genes in a list of gene objects.
    """
    locustag = org.name
    cmd = ["aragorn", "-t", "-gcbact", "-l", "-w", fna_file]
    p = Popen(cmd, stdout=PIPE)
    # loading the whole thing, reverting it to 'pop' in order.
    file_data = p.communicate()[0].decode().split("\n")[:: -1]
    gene_objs = defaultdict(set)
    c = 0
    header = ""
    while len(file_data) != 0:
        line = file_data.pop()
        if line.startswith(">"):
            header = line.replace(">", "").split()[0]
            file_data.pop()  # then next line must be removed too.
        elif len(line) > 0:  # if the line isn't empty, there's data to get.
            line_data = line.split()
            start, stop = ast.literal_eval(line_data[2].replace("c", ""))
            c += 1
            gene = RNA(identifier=locustag + '_tRNA_' + str(c).zfill(3))
            gene.fill_annotations(start=start, stop=stop, strand="-" if line_data[2].startswith(
                "c") else "+", gene_type="tRNA", product=line_data[1] + line_data[4])
            gene_objs[header].add(gene)
    return gene_objs


def launch_prodigal(fna_file, org, code, procedure):
    """
        launches Prodigal to annotate CDS. Takes a fna file name and a locustag to give an ID to the found genes.
        returns the annotated genes in a list of gene objects.
    """
    locustag = org.name
    cmd = list(map(str, ["prodigal", "-f", "sco", "-g", code, "-m", "-c", "-i", fna_file, "-p", procedure, "-q"]))
    logging.getLogger().debug(f"prodigal command : {' '.join(cmd)}")
    p = Popen(cmd, stdout=PIPE)

    gene_objs = defaultdict(set)
    c = 0
    header = ""
    for line in p.communicate()[0].decode().split("\n"):
        if line.startswith("# Sequence Data: "):
            for data in line.split(";"):
                if data.startswith("seqhdr"):
                    header = data.split("=")[1].replace('"', "").split()[0]
                    # print(header)

        elif line.startswith(">"):
            c += 1
            line_data = line[1:].split("_")  # not considering the '>'
            gene = Gene(gene_id=locustag + "_CDS_" + str(c).zfill(4))
            gene.fill_annotations(start=line_data[1], stop=line_data[2], strand=line_data[3],
                                  gene_type="CDS", genetic_code=code)
            gene_objs[header].add(gene)

    return gene_objs


def launch_infernal(fna_file, org, kingdom, tmpdir):
    """
        launches Infernal in hmmer-only mode to annotate rRNAs.
        Takes a fna file name and a locustag to give an ID to the found genes.
        returns the annotated genes in a list of gene objects.
    """
    locustag = org.name
    modelfile = ""
    if kingdom == "bacteria":
        modelfile = os.path.dirname(os.path.realpath(__file__)) + "/rRNA_DB/rRNA_bact.cm"
    elif kingdom == "archaea":
        modelfile = os.path.dirname(os.path.realpath(__file__)) + "/rRNA_DB/rRNA_arch.cm"

    tmp_file = tempfile.NamedTemporaryFile(mode="r", dir=tmpdir)
    cmd = ["cmscan", "--tblout", tmp_file.name, "--hmmonly", "--cpu", str(1), "--noali", modelfile, fna_file]
    p = Popen(cmd, stdout=open(os.devnull, "w"), stderr=PIPE)
    err = p.communicate()[1].decode().split()
    if err:
        if err[0] == 'Error: ':
            raise Exception(f"Infernal (cmscan) failed with error:  '{' '.join(err)}'. If you never used this script,"
                            f" you should press the .cm file using cmpress executable from Infernal. "
                            f"You should find the file in '{os.path.dirname(os.path.realpath(__file__))}/rRNA_DB/'.")
        raise Exception(f"An error occurred with Infernal. Error is:  '{' '.join(err)}'.")
    # never managed to test what happens if the .cm files are compressed with a 'bad' version of infernal,
    # so if that happens you are on your own.

    gene_objs = defaultdict(set)
    c = 0
    for line in tmp_file:
        if not line.startswith("#"):
            c += 1
            line_data = line.split()
            strand = line_data[9]
            if strand == "-":
                start = line_data[8]
                stop = line_data[7]
            else:
                start = line_data[7]
                stop = line_data[8]
            gene = RNA(identifier=locustag + "_rRNA_" + str(c).zfill(3))
            gene.fill_annotations(start=start, stop=stop, strand=strand, gene_type="rRNA",
                                  product=" ".join(line_data[17:]))
            gene_objs[line_data[2]].add(gene)

    return gene_objs


def read_fasta(org, fna_file, contig_filter=1):
    """
        Reads a fna file (or stream, or string) and stores it in a dictionary with contigs as key and sequence as value.
    """
    try:
        contigs = {}
        contig_seq = ""
        all_contig_len = 0
        contig = None
        for line in fna_file:
            if line.startswith('>'):
                if len(contig_seq) >= contig_filter:
                    contigs[contig.name] = contig_seq.upper()
                    all_contig_len += len(contig_seq)
                contig_seq = ""
                contig = org.get_or_add_contig(line.split()[0][1:])
            else:
                contig_seq += line.strip()
        # processing the last contig
        if len(contig_seq) >= contig_filter:
            contigs[contig.name] = contig_seq.upper()
            all_contig_len += len(contig_seq)
    except AttributeError as e:
        raise AttributeError(f"{e}\nAn error was raised when reading file: '{fna_file.name}'. "
                             f"One possibility for this error is that the file did not start with a '>' "
                             f"as it would be expected from a fna file.")
    return contigs, all_contig_len


def write_tmp_fasta(contigs, tmpdir):
    """
        Writes a temporary fna formated file, and returns the file-like object.
        This is for the cases where the given file is compressed,
        then we write a temporary file for the annotation tools to read from.
        The file will be deleted when close() is called.
    """
    tmp_file = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir)
    for header in contigs.keys():
        tmp_file.write(f">{header}\n")
        j = 0
        while j < len(contigs[header]):
            tmp_file.write(contigs[header][j: j + 60] + "\n")
            j += 60
    tmp_file.flush()  # force write what remains in the buffer.
    return tmp_file


def syntaxic_annotation(org, fasta_file, norna, kingdom, code, procedure, tmpdir):
    """
        Runs the different software for the syntaxic annotation.

        Takes in the file-like object containing the uncompressed fasta sequences to annotate
        the number of cpus that we can use.
        whether to annotate rna or not
        the locustag to give gene IDs.
    """
    # launching tools for syntaxic annotation
    genes = defaultdict(list)
    for key, items in launch_prodigal(fasta_file.name, org, code, procedure).items():
        genes[key].extend(items)
    if not norna:
        for key, items in launch_aragorn(fasta_file.name, org).items():
            genes[key].extend(items)
        for key, items in launch_infernal(fasta_file.name, org, kingdom, tmpdir).items():
            genes[key].extend(items)
    fasta_file.close()  # closing either tmp file or original fasta file.
    return genes


def overlap_filter(all_genes, overlap):
    """
        Removes the CDS that overlap with RNA genes.
    """
    sorted_genes = defaultdict(list)
    for key, genes in all_genes.items():
        tmp_genes = sorted(genes, key=lambda x: x.start)
        rm_genes = set()
        if overlap:
            for i, gene_i in enumerate(tmp_genes):
                if i + 1 < len(tmp_genes):
                    gene_j = tmp_genes[i + 1]
                    if gene_i.type != "CDS" and gene_j.type == "CDS" and gene_i.stop > gene_j.start:
                        rm_genes.add(gene_j)
                    elif gene_i.type == "CDS" and gene_j.type != "CDS" and gene_i.stop > gene_j.start:
                        rm_genes.add(gene_i)

        for gene in rm_genes:
            tmp_genes.remove(gene)
        cds_counter = 0
        for gene in tmp_genes:
            if gene.type == "CDS":
                gene.position = cds_counter
                cds_counter += 1
        sorted_genes[key] = tmp_genes
    return sorted_genes


def get_dna_sequence(contig_seq, gene):
    if gene.strand == "+":
        return contig_seq[gene.start - 1:gene.stop]
    elif gene.strand == "-":
        return reverse_complement(contig_seq[gene.start - 1:gene.stop])


def annotate_organism(org_name, file_name, circular_contigs, code, kingdom, norna, tmpdir, overlap, contig_filter,
                      procedure):
    """
        Function to annotate a single organism
    """
    org = Organism(org_name)

    fasta_file = read_compressed_or_not(file_name)
    contig_sequences, all_contig_len = read_fasta(org, fasta_file, contig_filter)
    if is_compressed(file_name):
        fasta_file = write_tmp_fasta(contig_sequences, tmpdir)
    if procedure is None:  # prodigal procedure is not force by user
        if all_contig_len < 20000:  # case of short sequence
            procedure = "meta"
        else:
            procedure = "single"
    genes = syntaxic_annotation(org, fasta_file, norna, kingdom, code, procedure, tmpdir)
    genes = overlap_filter(genes, overlap)

    for contigName, genes in genes.items():
        contig = org.get_or_add_contig(contigName)
        if contig.name in circular_contigs:
            contig.is_circular = True
        for gene in genes:
            gene.add_dna(get_dna_sequence(contig_sequences[contig.name], gene))
            gene.fill_parents(org, contig)
            if isinstance(gene, Gene):
                contig.add_gene(gene)
            elif isinstance(gene, RNA):
                contig.add_rna(gene)
    return org
