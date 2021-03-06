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
from ppanggolin.RGP.genomicIsland import compute_org_rgp
from ppanggolin.RGP.spot import makeSpotGraph, draw_spots, subgraph


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
        try:
            panFam = pangenome.getGeneFamily(fam.name)
            panMulti.add(panFam)
        except KeyError:#the family is not in the genome
            pass
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

def projectRGP(pangenome, annotation, output, tmpdir, identity = 0.8, coverage=0.8, defrag = False, cpu = 1, translation_table = 11, pseudo=False):
    if pangenome.status["geneFamilySequences"] not in ["inFile","Loaded","Computed"]:
        raise Exception("Cannot use this function as your pangenome does not have gene families representatives associated to it. For now this works only if the clustering is realised by PPanGGOLiN.")
    ## could be possible either by picking a representative somehow, or by aligning on genes rather than on families, if they are in the pangenome.


    #read given file
    logging.getLogger().info("Retrieving the annotations from the given file")
    singleOrgPang = Pangenome()#need to create a new 'pangenome' as the annotation reading functions take a pangenome as input.
    filetype = detect_filetype(annotation)
    if filetype == "gff":
        org, hasFasta = read_org_gff('myGenome', annotation, [], True, pseudo=pseudo)
        singleOrgPang.addOrganism(org)
    elif filetype == "gbff":
        org, hasFasta = read_org_gbff('myGenome', annotation, [], pseudo=pseudo)
        singleOrgPang.addOrganism(org)
    if hasFasta:
        singleOrgPang.status["geneSequences"] = "Computed"
    else:
        raise Exception(f"The given annotation file did not have a FASTA sequence included (expected '##FASTA' pragma followed by a fasta-like file format). This is required for computing the Regions of Genomic Plasticity of your organism")

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

    multigenics = pangenome.get_multigenics(pangenome.parameters["RGP"]["dup_margin"])
    genomeMultigenics = linkMultigenicFamilies(singleOrgPang, multigenics)

    logging.getLogger().info("Predicting RGP in your genome")
    for org in singleOrgPang.organisms:
        genomeRGP = compute_org_rgp(org, pangenome.parameters["RGP"]["persistent_penalty"], pangenome.parameters["RGP"]["variable_gain"], pangenome.parameters["RGP"]["min_length"], pangenome.parameters["RGP"]["min_score"], genomeMultigenics)

    #write the rgps in a single file
    write_RGPs_cgview(genomeRGP, output)
    write_partitions_cgview(org, output)

    if filetype == "gff":
        #reread the file and insert sequence_feature objects corresponding to the predicted regions
        logging.getLogger().info("Writing the RGP in a gff file...")
        writeGffRegions(annotation, genomeRGP, output)
    elif filetype == "gbff":
        logging.getLogger().info("Writing the RGP in a gbff file...")
        writeGbffRegions(annotation, genomeRGP, output)

def write_partitions_cgview(organism, output):
    """writes the partition of each gene in a table that is compatible with cgview server"""
    with open(output + "/partitions.cgview.tab","w") as fout:
        fout.write("name\ttype\tstart\tstop\tstrand\n")
        for gene in organism.genes:
            fout.write("\t".join(map(str,[gene.ID, gene.family.namedPartition,gene.start, gene.stop, gene.strand]))+ "\n")
        fout.close()
    logging.getLogger().info(f"Done writing all of the partitions of each gene in a tsv file: {output + '/partitions.cgview.tab'}")

def write_RGPs_cgview(regions, output):
    """write RGPs in a table that is compatible with cgview server"""
    with open(output + "/regions.cgview.tab","w") as fout:
        fout.write("name\ttype\tstart\tstop\tstrand\n")
        for region in regions:
            fout.write("\t".join(map(str,[region.name, "RGP",region.start, region.stop, "+"]))+ "\n")
        fout.close()
    logging.getLogger().info(f"Done writing all of the RGP in a tsv file: {output + '/regions.cgview.tab'}")

def getFam2RGP(pangenome, multigenics):
    """associates families to the RGP they belong to, and those they are bordering"""
    fam2rgp = defaultdict(list)
    for rgp in pangenome.regions:
        for fam in rgp.families:
            fam2rgp[fam].append(rgp.name)
        for fam in [ gene.family for border in rgp.getBorderingGenes(pangenome.parameters["spots"]["set_size"], multigenics) for gene in border ]:
            fam2rgp[fam].append(rgp.name)
    return fam2rgp

def getFam2spot(pangenome, output, multigenics):
    """ reads a pangenome object and returns a dictionnary of family to RGP and family to spot, that indicates where each family is"""
    ##those are to be replaced as spots should be stored in the pangenome, and in the h5.
    fam2spot = defaultdict(list)
    fam2border = defaultdict(list)
    for spot in pangenome.spots:
        fams = set()
        famsBorder = set()
        for rgp in spot.regions:
            fams |= rgp.families
            famsBorder |= set([ gene.family for border in rgp.getBorderingGenes(pangenome.parameters["spots"]["set_size"], multigenics) for gene in border ])
        for fam in fams:
            fam2spot[fam].append(spot)
        for fam in famsBorder:
            fam2border[fam].append(spot)
    return fam2spot, fam2border,  multigenics

def add_spot_str(a):
    return "spot_" + str(a.ID)

def draw_spot_gexf(spots, output, multigenics, set_size = 3):
    for spot in spots:
        fname = "spot_" + str(spot.ID) + ".gexf"
        subgraph(spot,output, fname, set_size=set_size, multigenics=multigenics)

def checkLabelPriorityLogic(priority):
    for p in priority.split(','):
        if p.lower() not in ["name","id","family"]:
            raise Exception(f"You have indicated a label which is not supported with --label_priority. You indicated '{p}'. Supported labels are 'name', 'id' and 'family'")

def getProtInfo(prot2pang, pangenome, output, cpu, draw_related, priority):
    logging.getLogger().info("Writing RGP and spot information related to hits in the pangenome")
    multigenics = pangenome.get_multigenics(pangenome.parameters["RGP"]["dup_margin"])

    finfo = open(output+"/info_input_prot.tsv","w")
    finfo.write("input\tfamily\tpartition\tspot_list_as_member\tspot_list_as_border\trgp_list\n")
    fam2rgp = getFam2RGP(pangenome, multigenics)
    fam2spot, fam2border, multigenics= getFam2spot(pangenome, output, multigenics)
    spot_list = set()
    for prot, panfam in prot2pang.items():
        finfo.write(prot +'\t' + panfam.name + "\t" + panfam.namedPartition + "\t" + ",".join(map(add_spot_str, fam2spot[panfam])) + "\t" + ",".join(map(add_spot_str, fam2border[panfam])) + "\t" + ','.join(fam2rgp[panfam]) + "\n")
        spot_list |= set(fam2spot[panfam])
        spot_list |= set(fam2border[panfam])
    finfo.close()
    if draw_related:
        drawn_spots = set()
        for spot in spot_list:
            if len(spot.getUniqOrderedSet()) > 1:
                drawn_spots.add(spot)
        logging.getLogger().info(f"Drawing the {len(drawn_spots)} spots with more than 1 organization related to hits of the input proteins...")
        draw_spots(drawn_spots, output, cpu, 2, 1, 3, multigenics, [], priority)
        draw_spot_gexf(drawn_spots, output, multigenics = multigenics)

    logging.getLogger().info(f"File listing RGP and spots where proteins of interest are located : '{output+'/info_input_prot.tsv'}'")

def align(pangenome, proteinFile, output, tmpdir, identity = 0.8, coverage=0.8, defrag = False, cpu = 1, getinfo = False, draw_related = False, priority='name,ID'):
    if pangenome.status["geneFamilySequences"] not in ["inFile","Loaded","Computed"]:
        raise Exception("Cannot use this function as your pangenome does not have gene families representatives associated to it. For now this works only if the clustering is realised by PPanGGOLiN.")
    ## could be possible either by picking a representative somehow, or by aligning on genes rather than on families, if they are in the pangenome.

    if getinfo:
        checkLabelPriorityLogic(priority)
        checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, needRGP=True, needPartitions=True, needSpots=True)
    else:
        checkPangenomeInfo(pangenome, needFamilies=True)

    newtmpdir = tempfile.TemporaryDirectory(dir = tmpdir)
    tmpPangFile = tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)

    writeGeneFamSequences(pangenome, tmpPangFile)

    with read_compressed_or_not(proteinFile) as protFileObj:
        protSet = getProt(protFileObj)
        alignFile = alignSeqToPang(tmpPangFile, protFileObj, output, newtmpdir, cpu, defrag, identity, coverage)

    prot2pang = readAlignments(alignFile, pangenome)
    

    if getinfo:
        getProtInfo(prot2pang, pangenome, output, cpu, draw_related, priority)
    else:
        partProj = projectPartition(prot2pang, protSet, output)#write the partition assignation only
        logging.getLogger().info(f"proteins partition projection : '{partProj}'")
    logging.getLogger().info(f"{len(prot2pang)} proteins over {len(protSet)} have at least one hit in the pangenome.")
    logging.getLogger().info(f"Blast-tab file of the alignment : '{alignFile}'")

    tmpPangFile.close()
    newtmpdir.cleanup()

def launch(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    if args.proteins is not None:
        align(pangenome = pangenome, proteinFile = args.proteins, output = args.output, tmpdir = args.tmpdir, identity = args.identity, coverage =args.coverage, defrag =args.defrag,cpu= args.cpu, getinfo =args.getinfo, draw_related = args.draw_related, priority=args.label_priority )

    if args.annotation is not None:
        projectRGP(pangenome, args.annotation, args.output, args.tmpdir, args.identity, args.coverage, args.defrag, args.cpu, args.translation_table, pseudo=args.use_pseudo)

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
    optional.add_argument("--getinfo", required=False, action="store_true",help = "Use this option to extract info related to the best hit of each query, such as the RGP it is in, or the spots.")
    optional.add_argument("--draw_related", required=False, action = "store_true",help = "Draw figures and provide graphs in a gexf format of the eventual spots associated to the input proteins")
    optional.add_argument("--use_pseudo",required=False, action="store_true",help = "In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them)")
    optional.add_argument("--label_priority", required=False, type=str, default='name,ID', help = "Option to use with --draw_hotspots. Will indicate what to write in the figure labels as a comma-separated list. Order gives priority. Possible values are: name (for gene names), ID (for the gene IDs), family (for the family IDs)")
    return parser