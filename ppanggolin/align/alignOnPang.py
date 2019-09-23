#!/usr/bin/env python3
#coding:utf-8

#default libraries
import logging
import tempfile
import subprocess
import argparse

#local libraries
from ppanggolin.formats import checkPangenomeInfo
from ppanggolin.utils import mkOutdir, read_compressed_or_not
from ppanggolin.pangenome import Pangenome


def createdb(fileObj, tmpdir):
    seqdb =  tempfile.NamedTemporaryFile(mode="w", dir = tmpdir.name)
    cmd = ["mmseqs","createdb",fileObj.name, seqdb.name]
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    return seqdb

def alignProtToPang(pangFile, protFile,  output, tmpdir, cpu = 1, defrag=False, identity = 0.8, coverage = 0.8):
    pangdb = createdb(pangFile, tmpdir)
    protdb = createdb(protFile, tmpdir)
    covmode = "0"
    if defrag:
        covmode = "1"
    alndb =  tempfile.NamedTemporaryFile(mode="w", dir = tmpdir.name)
    cmd = ["mmseqs","search",protdb.name , pangdb.name, alndb.name, tmpdir.name, "-a","--min-seq-id", str(identity), "-c", str(coverage), "--cov-mode", covmode, "--threads", str(cpu)]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Aligning proteins to cluster representatives...")
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    outfile =  output + "/protein_to_pangenome_associations.blast-tab"
    cmd = ["mmseqs","convertalis", protdb.name ,pangdb.name, alndb.name, outfile,"--format-mode","2"]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Extracting alignments...")
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    pangdb.close()
    protdb.close()
    alndb.close()

    return outfile

def readAlignments(outfile, pangenome):
    prot2pang = {}
    with open(outfile,"r") as alnFile:
        for line in alnFile:
            line = line.split()
            if prot2pang.get(line[0]) is None:#if no results were found yet
                prot2pang[line[0]] = pangenome.getGeneFamily(line[1])#then the best hit is the first one we see.
    return prot2pang

def getProt(protFile):
    protset = set()
    for line in protFile:
        if line.startswith(">"):
            protset.add(line[1:])
    return protset

def writeGeneFamSequences(pangenome, fileObj):
    for fam in pangenome.geneFamilies:
        fileObj.write(">" + fam.name + "\n")
        fileObj.write(fam.sequence + "\n")
    fileObj.flush()

def projectPartition(prot2pang, protSet, output):
    partitionProj = output + "/proteins_partition_projection.tsv"
    with open(partitionProj, "w") as partProjFile:
        for key, pangFam in prot2pang.items():
            partProjFile.write(key + "\t" + pangFam.namedPartition + "\n")
        for remainingProt in (prot2pang.keys() & protSet):
            partProjFile.write(remainingProt + "\tcloud\n")#if there is no hit, it's going to be cloud genes.
    return partitionProj

def align(pangenome, proteinFile, output, tmpdir, identity = 0.8, coverage=0.8, defrag = False, cpu = 1):
    if pangenome.status["geneFamilySequences"] not in ["inFile","Loaded","Computed"]:
        raise Exception("Cannot use this function as your pangenome does not have gene families representatives associated to it. For now this works only if the clustering is realised by PPanGGOLiN.")
    checkPangenomeInfo(pangenome, needFamilies=True)

    newtmpdir = tempfile.TemporaryDirectory(dir = tmpdir)
    tmpPangFile = tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)

    writeGeneFamSequences(pangenome, tmpPangFile)

    with read_compressed_or_not(proteinFile) as protFileObj:
        protSet = getProt(protFileObj)
        alignFile = alignProtToPang(tmpPangFile, protFileObj, output, newtmpdir, cpu, defrag, identity, coverage)

    prot2pang = readAlignments(alignFile, pangenome)
    partProj = projectPartition(prot2pang, protSet, output)

    logging.getLogger().info(f"{len(prot2pang)} proteins over {len(protSet)} have at least one hit in the pangenome.")
    logging.getLogger().info(f"Blast-tab file of the alignment : '{alignFile}'")
    logging.getLogger().info(f"proteins partition projection : '{partProj}'")
    tmpPangFile.close()
    newtmpdir.cleanup()

def launch(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    align(pangenome, args.proteins, args.output, args.tmpdir, args.identity, args.coverage, args.defrag, args.cpu)

def alignSubparser(subparser):
    parser = subparser.add_parser("align", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title = "Required arguments", description = "All of the following arguments are required :")
    required.add_argument('--proteins', required = True, type = str, help = "proteins sequences to align on the pangenome gene families")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o','--output', required=True, type=str, help="Output directory where the file(s) will be written")

    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('--defrag', required=False,default=False, action="store_true", help = "Use the defragmentation strategy to associate potential fragments with their original gene family.")
    optional.add_argument('--identity', required = False, type = float, default=0.5, help = "min identity percentage threshold")
    optional.add_argument('--coverage', required = False, type = float, default=0.8, help = "min coverage percentage threshold")

    return parser