#!/usr/bin/env python3
#coding:utf-8

#default libraries
import logging
import tempfile
import subprocess
import argparse
from collections import defaultdict

#local libraries
from ppanggolin.formats import checkPangenomeInfo
from ppanggolin.utils import mkOutdir, read_compressed_or_not
from ppanggolin.pangenome import Pangenome
from ppanggolin.annotate import detect_filetype, read_org_gff, read_org_gbff
from ppanggolin.cluster import writeGeneSequencesFromAnnotations
from ppanggolin.RGP import get_multigenics, compute_org_rgp


def createdb(fileObj, tmpdir):
    seqdb =  tempfile.NamedTemporaryFile(mode="w", dir = tmpdir.name)
    cmd = ["mmseqs","createdb",fileObj.name, seqdb.name]
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    return seqdb

def alignSeqToPang(pangFile, seqFile, output, tmpdir, cpu = 1, defrag=False, identity = 0.8, coverage = 0.8, is_nucl = False, code = 11):
    pangdb = createdb(pangFile, tmpdir)
    seqdb = createdb(seqFile, tmpdir)
    if is_nucl:
        seqNucdb = tempfile.NamedTemporaryFile(mode="w", dir = tmpdir.name)
        cmd = ["mmseqs","translatenucs", seqdb.name, seqNucdb.name, "--threads", str(cpu), "--translation-table",str(code)]
        logging.getLogger().debug(" ".join(cmd))
        subprocess.run(cmd, stdout=subprocess.DEVNULL)
        seqdb = seqNucdb
    covmode = "0"
    if defrag:
        covmode = "1"
    alndb =  tempfile.NamedTemporaryFile(mode="w", dir = tmpdir.name)
    cmd = ["mmseqs","search",seqdb.name , pangdb.name, alndb.name, tmpdir.name, "-a","--min-seq-id", str(identity), "-c", str(coverage), "--cov-mode", covmode, "--threads", str(cpu)]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Aligning proteins to cluster representatives...")
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    outfile =  output + "/input_to_pangenome_associations.blast-tab"
    cmd = ["mmseqs","convertalis", seqdb.name ,pangdb.name, alndb.name, outfile,"--format-mode","2"]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Extracting alignments...")
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    pangdb.close()
    seqdb.close()
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


def linkNewGenomeFamilies(orgPangenome, formerPangenome, blastTab):
    with open(blastTab,"r") as alnFile:
        for line in alnFile:
            line = line.split()
            gene = orgPangenome.getGene(line[0])
            if gene.family is None:#if no results were found yet
                oldFam = formerPangenome.getGeneFamily(line[1])
                fam = orgPangenome.addGeneFamily(line[1])#then the best hit is the first one we see.
                fam.addGene(gene)
                fam.partition = oldFam.partition
    #the remaining genes with no hits are cloud genes.
    for gene in orgPangenome.genes:
        if gene.family is None:
            fam = orgPangenome.addGeneFamily(gene.ID)#create a new family
            fam.addGene(gene)
            fam.addPartition("C")

def linkMultigenicFamilies(pangenome, multigenics):
    panMulti = set()
    for fam in multigenics:
        panFam = pangenome.getGeneFamily(fam.name)
        if panFam is not None:
            panMulti.add(panFam)
    return panMulti


def writeGbffRegions(filename, regions, output):
    ContigRegions = defaultdict(set)
    for region in regions:
        ContigRegions[region.contig.name].add(region)

    for contigName in ContigRegions.keys():
        ContigRegions[contigName] = sorted(ContigRegions[region.contig.name], key= lambda x : x.start, reverse=True)

    foutfile = open(output + "/genome_annotation.gbff", "w")

    curr_contig = None
    with read_compressed_or_not(filename) as fannot:
        for line in fannot:
            if line.startswith("VERSION"):
                curr_contig = line.split()[1]
            if curr_contig in ContigRegions and len(ContigRegions[curr_contig]) > 0:
                if line[0:5].strip() == "" and line[0:20].strip() != "" and len(line[20:].split("..")) == 2:#should be a FEATURE with its position
                    start = line[20:].replace("complement(","").replace(")","").split("..")[0]
                    if int(start) == ContigRegions[curr_contig][-1].start:
                        reg = ContigRegions[curr_contig].pop()
                        foutfile.write("     misc_feature    " + str(reg.start) + ".." + str(reg.stop)+"\n")
                        foutfile.write('                     /note="Region of genomic plasticity"\n')

            foutfile.write(line)
    
    for val in ContigRegions.values():
        if len(val) != 0:
            logging.getLogger().warning("Somes regions were not written in the new gbff file for unknown reasons")
    logging.getLogger().info(f"RGP have been written in the following file : '{output + '/genome_annotation.gbff'}' ")

def writeGffRegions(filename, regions, output):
    ContigRegions = defaultdict(set)
    for region in  regions:
        ContigRegions[region.contig.name].add(region)

    for contigName in ContigRegions.keys():
        ContigRegions[contigName] = sorted(ContigRegions[region.contig.name], key= lambda x : x.start, reverse=True)
    
    foutfile = open(output + "/genome_annotation.gff", "w")

    with read_compressed_or_not(filename) as fannot:
        for line in fannot:
            if line[0] == "#":
                pass
            else:
                features = line.split("\t")
                if len(features) == 9:#gff annotation lines are supposed to be 8 columns long
                    start = int(features[3])
                    if features[0] in ContigRegions:
                        if len(ContigRegions[features[0]]) > 0 and ContigRegions[features[0]][-1].start == start:
                            reg = ContigRegions[features[0]].pop()
                            foutfile.write('\t'.join(map(str,[features[0], "panRGP","sequence_feature",reg.start, reg.stop, reg.score, '+', '.', f'ID={region.name};note=Region of genomic plasticity;gbkey=misc_feature'])) + "\n")
            foutfile.write(line)

    for val in ContigRegions.values():
        if len(val) != 0:
            logging.getLogger().warning("Somes regions were not written in the new gff file for unknown reasons")
    logging.getLogger().info(f"RGP have been written in the following file : '{output + '/genome_annotation.gff'}' ")

def projectRGP(pangenome, annotation, output, tmpdir, identity = 0.8, coverage=0.8, defrag = False, cpu = 1, translation_table = 11):
    if pangenome.status["geneFamilySequences"] not in ["inFile","Loaded","Computed"]:
        raise Exception("Cannot use this function as your pangenome does not have gene families representatives associated to it. For now this works only if the clustering is realised by PPanGGOLiN.")

    #read given file
    logging.getLogger().info("Retrieving the annotations from the given file")
    singleOrgPang = Pangenome()#need to create a new 'pangenome' as the annotation reading functions take a pangenome as input.
    filetype = detect_filetype(annotation)
    if filetype == "gff":
        singleOrgPang.status["geneSequences"] = "Computed"#if there are no sequences in the gff, this value will change to 'No'
        read_org_gff(singleOrgPang, 'myGenome', annotation, [], True)
        if singleOrgPang.status["geneSequences"] == "No":
            raise Exception(f"The given annotation file did not have a FASTA sequence included (expected '##FASTA' pragma followed by a fasta-like file format). This is required for computing the Regions of Genomic Plasticity of your organism")
    elif filetype == "gbff":
        read_org_gbff(singleOrgPang, 'myGenome', annotation, [], True)

    #check and read given pangenome
    checkPangenomeInfo(pangenome, needFamilies=True, needPartitions = True, needAnnotations = True)

    newtmpdir = tempfile.TemporaryDirectory(dir = tmpdir)
    tmpPangFile = tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)
    tmpGeneFile = tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)

    writeGeneSequencesFromAnnotations(singleOrgPang, tmpGeneFile)
    writeGeneFamSequences(pangenome, tmpPangFile)

    blastout = alignSeqToPang(tmpPangFile, tmpGeneFile, output, newtmpdir, cpu, defrag, identity, coverage, True, translation_table)

    tmpPangFile.close()
    tmpGeneFile.close()
    newtmpdir.cleanup()
    #artificially reconstruct the gene families and their partitions
    linkNewGenomeFamilies(singleOrgPang, pangenome, blastout)

    multigenics = get_multigenics(pangenome, pangenome.parameters["RGP"]["dup_margin"])
    genomeMultigenics = linkMultigenicFamilies(singleOrgPang, multigenics)

    logging.getLogger().info("Predicting RGP in your genome")
    for org in singleOrgPang.organisms:
        genomeRGP = compute_org_rgp(org, pangenome.parameters["RGP"]["persistent_penalty"], pangenome.parameters["RGP"]["variable_gain"], pangenome.parameters["RGP"]["min_length"], pangenome.parameters["RGP"]["min_score"], genomeMultigenics)

    if filetype == "gff":
        #reread the file and insert sequence_feature objects corresponding to the predicted regions
        logging.getLogger().info("Writing the RGP in a gff file...")
        writeGffRegions(annotation, genomeRGP, output)
    elif filetype == "gbff":
        logging.getLogger().info("Writing the RGP in a gbff file...")
        writeGbffRegions(annotation, genomeRGP, output)


def align(pangenome, proteinFile, output, tmpdir, identity = 0.8, coverage=0.8, defrag = False, cpu = 1):
    if pangenome.status["geneFamilySequences"] not in ["inFile","Loaded","Computed"]:
        raise Exception("Cannot use this function as your pangenome does not have gene families representatives associated to it. For now this works only if the clustering is realised by PPanGGOLiN.")
    checkPangenomeInfo(pangenome, needFamilies=True)

    newtmpdir = tempfile.TemporaryDirectory(dir = tmpdir)
    tmpPangFile = tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)

    writeGeneFamSequences(pangenome, tmpPangFile)

    with read_compressed_or_not(proteinFile) as protFileObj:
        protSet = getProt(protFileObj)
        alignFile = alignSeqToPang(tmpPangFile, protFileObj, output, newtmpdir, cpu, defrag, identity, coverage)

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
    if args.proteins is not None:
        align(pangenome, args.proteins, args.output, args.tmpdir, args.identity, args.coverage, args.defrag, args.cpu)

    if args.annotation is not None:
        projectRGP(pangenome, args.annotation, args.output, args.tmpdir, args.identity, args.coverage, args.defrag, args.cpu, args.translation_table)

def alignSubparser(subparser):
    parser = subparser.add_parser("align", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title = "Required arguments", description = "All of the following arguments are required :")
    onereq = parser.add_argument_group(title = "Input file", description = "One of the following argument is required :")
    onereq.add_argument('--proteins', required = False, type = str, help = "proteins sequences to align on the pangenome gene families")
    onereq.add_argument('--annotation', required = False, type = str, help = "annotation input file (gff or gbff) from which to predict RGPs and partitions")

    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o','--output', required=True, type=str, help="Output directory where the file(s) will be written")

    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('--defrag', required=False,default=False, action="store_true", help = "Use the defragmentation strategy to associate potential fragments with their original gene family.")
    optional.add_argument('--identity', required = False, type = float, default=0.5, help = "min identity percentage threshold")
    optional.add_argument('--coverage', required = False, type = float, default=0.8, help = "min coverage percentage threshold")
    optional.add_argument("--translation_table",required=False, default="11", help = "Translation table (genetic code) to use.")
    return parser