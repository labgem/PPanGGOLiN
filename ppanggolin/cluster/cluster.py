#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse
import logging
import tempfile
import subprocess
from collections import defaultdict
import os


#installed libraries
import networkx
from tqdm import tqdm

#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Gene
from ppanggolin.utils import read_compressed_or_not, getCurrentRAM
from ppanggolin.formats import writePangenome, readPangenome, getGeneSequencesFromFile
from ppanggolin.annotate import genetic_codes, translate

def alignRep(faaFile, tmpdir, cpu):
    newtmpdir = tempfile.TemporaryDirectory(dir = tmpdir.name)#create a tmpdir in the tmpdir provided.
    seqdb =  tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)
    cmd = ["mmseqs","createdb",faaFile.name, seqdb.name]
    logging.getLogger().debug(" ".join(cmd))
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    alndb =  tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)
    cmd = ["mmseqs","search",seqdb.name , seqdb.name, alndb.name, newtmpdir.name, "-a","--min-seq-id", "0.8", "-c", "0.8", "--cov-mode", "1", "--threads", str(cpu)]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Aligning cluster representatives...")
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    outfile =  tempfile.NamedTemporaryFile(mode="w", dir = tmpdir.name)
    cmd = ["mmseqs","convertalis", seqdb.name ,seqdb.name, alndb.name, outfile.name,"--format-output","query,target,qlen,tlen,bits"]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Extracting alignments...")
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    seqdb.close()
    alndb.close()
    newtmpdir.cleanup()
    return outfile

def firstClustering(sequences, tmpdir, cpu, code ):
    newtmpdir = tempfile.TemporaryDirectory(dir = tmpdir.name)#create a tmpdir in the tmpdir provided.
    seqNucdb = tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)
    cmd = ["mmseqs","createdb"] 
    cmd.append(sequences.name)
    cmd.extend([seqNucdb.name,"--dont-shuffle","false"])
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Creating sequence database...")
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    seqdb = tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)
    cmd = ["mmseqs","translatenucs", seqNucdb.name, seqdb.name, "--threads", str(cpu), "--translation-table",code]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Translating sequences...")
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    cludb = tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)
    cmd = ["mmseqs","cluster",seqdb.name, cludb.name, newtmpdir.name , "--min-seq-id", "0.8", "-c", "0.8", "--threads", str(cpu), "--kmer-per-seq","80","--max-seqs","300"]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Clustering sequences...")
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    logging.getLogger().info("Extracting cluster representatives...")
    repdb = tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)
    cmd = ["mmseqs","result2repseq", seqdb.name, cludb.name, repdb.name]
    logging.getLogger().debug(" ".join(cmd))
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    reprfa = tempfile.NamedTemporaryFile(mode="w", dir = tmpdir.name)
    cmd = ["mmseqs","result2flat",seqdb.name, seqdb.name, repdb.name, reprfa.name]
    logging.getLogger().debug(" ".join(cmd))
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    outtsv = tempfile.NamedTemporaryFile(mode="w", dir = tmpdir.name)
    cmd = ["mmseqs","createtsv",seqdb.name, seqdb.name, cludb.name,outtsv.name]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Writing gene to family informations")
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    repdb.close()
    seqdb.close()
    cludb.close()
    seqNucdb.close()
    newtmpdir.cleanup()#deleting temporary directory.
    return reprfa, outtsv

def read_faa(faFileName):
    fam2seq = {}
    head = ""
    with open(faFileName.name,"r") as faFile:
        for line in faFile:
            if line.startswith('>'):
                head = line[1:].strip()
            else:
                fam2seq[head] = line.strip()
    return fam2seq

def read_tsv(tsvfileName):
    # reading tsv file
    genes2fam = {}
    fam2genes = defaultdict(set)
    with open(tsvfileName.name, "r") as tsvfile:
        for line in tsvfile:
            line = line.split()
            genes2fam[line[1]] = (line[0],False)#fam id, and its a gene (and not a fragment)
            fam2genes[line[0]].add(line[1])
    return genes2fam, fam2genes

def refineClustering(tsv, alnFile, fam2seq):
    simgraph = networkx.Graph()
    genes2fam, fam2genes = read_tsv(tsv)
    logging.getLogger().info(f"Starting with {len(fam2seq)} families")
    #create the nodes
    for fam, genes in fam2genes.items():
        simgraph.add_node(fam, nbgenes = len(genes) )
    #add the edges
    with open(alnFile.name,"r") as alnfile:
        for line in alnfile:
            line = line.split()
            
            if line[0] != line[1]:
                simgraph.add_edge(line[0],line[1], score = float(line[4]))
                simgraph.nodes[line[0]]["length"] = int(line[2])
                simgraph.nodes[line[1]]["length"] = int(line[3])
    for node, nodedata in simgraph.nodes(data = True):
        choice = (None, 0, 0, 0)
        for neighbor in simgraph.neighbors(node):
            nei = simgraph.nodes[neighbor]
            score = simgraph[neighbor][node]["score"]
            if nei["length"] > nodedata["length"] and nei["nbgenes"] >= nodedata["nbgenes"] and  choice[3] < score:
                choice = (genes2fam[neighbor][0], nei["length"] , nei["nbgenes"], score)#genes2fam[neighbor] instead of just neighbor in case that family has been assigned already (this is for smaller fragments that are closer to other fragments than the actual gene family)
        if choice[0] is not None:
            genestochange = fam2genes[node]
            for gene in genestochange:
                genes2fam[gene] = (choice[0], "F")
                fam2genes[choice[0]].add(gene)
            del fam2genes[node]
    newFam2seq = {}
    for fam in fam2genes:
        newFam2seq[fam] = fam2seq[fam]
    logging.getLogger().info(f"Ending with {len(newFam2seq)} gene families")
    return genes2fam, newFam2seq

def read_gene2fam(pangenome, gene2fam):
    logging.getLogger().info(f"Adding {len(gene2fam)} genes to the gene families")

    link = True if pangenome.status["genomesAnnotated"] in ["Computed","Loaded"] else False
    
    for gene, (family, is_frag) in gene2fam.items():
        fam = pangenome.addGeneFamily(family)
        if link:#doing the linking if the annotations are loaded.
            geneObj = pangenome.getGene(gene)
        else:
            geneObj = Gene(gene)
        geneObj.is_fragment = is_frag
        fam.addGene(geneObj)

def read_fam2seq(pangenome, fam2seq):
    logging.getLogger().info("Adding protein sequences to the gene families")
    for family, protein in fam2seq.items():
        fam = pangenome.addGeneFamily(family)
        fam.addSequence(protein)

def writeGeneSequencesFromAnnotations(pangenome, fileObj):
    """
        Writes the CDS sequences of the Pangenome object to a tmpFile object
        Loads the sequences from previously computed or loaded annotations
    """
    logging.getLogger().info("Writing all of the CDS sequences from a Pangenome file to a fasta file")
    bar =  tqdm(pangenome.genes, unit="gene")
    for gene in bar:#reading the table chunk per chunk otherwise RAM dies on big pangenomes
        if gene.type == "CDS":
            fileObj.write('>' + gene.ID + "\n")
            fileObj.write(gene.dna + "\n")
    fileObj.flush()
    bar.close()

def checkPangenomeForClustering(pangenome, tmpFile):
    """
        Check the pangenome statuses and write the gene sequences in the provided tmpFile. (whether they are written in the .h5 file or currently in memory)
    """
    if pangenome.status["geneSequences"] in ["Computed","Loaded"]:
        writeGeneSequencesFromAnnotations(pangenome, tmpFile)
    elif pangenome.status["geneSequences"] == "inFile":
        getGeneSequencesFromFile(pangenome, tmpFile)#write CDS sequences to the tmpFile

def clustering(pangenome, tmpdir, cpu , defrag = False, code = "11"):
    newtmpdir = tempfile.TemporaryDirectory(dir = tmpdir)
    tmpFile = tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)

    checkPangenomeForClustering(pangenome, tmpFile)

    logging.getLogger().info("Clustering all of the genes sequences...")
    rep, tsv = firstClustering(tmpFile, newtmpdir, cpu, code)
    fam2seq = read_faa(rep)
    if not defrag:
        genes2fam = read_tsv(tsv)[0]
        
    else:
        logging.getLogger().info("Associating fragments to their original gene family...")
        aln = alignRep(rep, newtmpdir, cpu)
        genes2fam, fam2seq = refineClustering(tsv, aln, fam2seq)
        aln.close()
        pangenome.status["defragmented"] = "Computed"
    tmpFile.close()
    tsv.close()
    rep.close()
    newtmpdir.cleanup()
    read_fam2seq(pangenome, fam2seq)
    read_gene2fam(pangenome, genes2fam)
    pangenome.status["genesClustered"] = "Computed"
    pangenome.status["geneFamilySequences"] = "Computed"

def checkPangenomeForReadClustering(pangenome):
    if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"]:
        pass#nothing to do
    elif pangenome.status["genomesAnnotated"] == "inFile":
        readPangenome(pangenome, annotation = True)#reading the annotations so we are sure that the IDs in the cluster file matches the previously provided annotations.
    else:
        raise NotImplementedError("You tried to do something unplanned for")

def readClustering(pangenome, families_tsv_file):
    """
        Creates the pangenome, the gene families and the genes with an associated gene family.
        Reads a families tsv file from mmseqs2 output and adds the gene families and the genes to the pangenome.
    """
    checkPangenomeForReadClustering(pangenome)

    logging.getLogger().info("Reading "+families_tsv_file+" the gene families file ...")
    families_tsv_file = read_compressed_or_not(families_tsv_file)
    frag = False
    gene2families = {}
    for line in families_tsv_file:
        elements = [el.strip() for el in line.split()] # 2 or 3 fields expected
        if len(elements)<=1:
            logging.getLogger().error("No tabulation separator found in gene families file")
            exit(1)
        (fam_id, gene_id, is_frag) = elements if len(elements) == 3 else elements+[None]
        gene2families[gene_id] = (fam_id, True if is_frag == "F" else False)
        if is_frag == "F":
            frag=True
    families_tsv_file.close()
    read_gene2fam(pangenome, gene2families)
    pangenome.status["genesClustered"] = "Computed"
    if frag:#if there was fragment informations in the file.
        pangenome.status["defragmented"] = "Computed"

def launch(args):
    """ launch the clustering step"""
    # pangenome.stats()
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    if args.clusters is None:
        clustering(pangenome, args.tmpdir, args.cpu, args.defrag, args.translation_table)
        logging.getLogger().info("Done with the clustering")
    else:
        readClustering(pangenome, args.clusters)
        logging.getLogger().info("Done reading the cluster file")
    writePangenome(pangenome, pangenome.file, args.force)

def clusterSubparser(subparser):
    parser = subparser.add_parser("cluster",help = "Cluster proteins in protein families")
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('--defrag', required=False,default=False, action="store_true", help = "Use the defragmentation strategy to associated potential fragments with their original gene family.")
    optional.add_argument("--translation_table",required=False, default="11", help = "Translation table (genetic code) to use.")
    optional.add_argument('--clusters', required = False, type = str, help = "A tab-separated list containing the result of a clustering. One line per gene. First column is cluster ID, and second is gene ID")
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="A tab-separated file listing the organism names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per organism.")
    return parser