#!/usr/bin/env python3
#coding:utf-8

#default libraries
import logging
import argparse
import time
import os
import tempfile
import subprocess
from itertools import combinations
from statistics import mean
#installed libraries
from tqdm import tqdm
import networkx as nx
from gmpy2 import xmpz, popcount
#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.region import Region
from ppanggolin.formats import checkPangenomeInfo, writePangenome
from ppanggolin.utils import mkOutdir


class Module:
    def __init__(self, geneFamilies):
        if not all(isinstance(fam, GeneFamily) for fam in geneFamilies):
            raise Exception(f"You provided elements that were not GeneFamily objetcs. Modules are only made of GeneFamily")
        self.families = set(geneFamilies)

    @property
    def distance(self):
        if hasattr(self, "_distance"):
            return self._distance
        else:
            self._get_dist()
            return self._distance

    def _get_dist(self):
        union = xmpz(0)
        inter = xmpz(-1)
        for fam in self.families:
            union = union | fam.bitarray
            inter = inter & fam.bitarray
        self._distance = float(popcount(inter) / popcount(union))

def writeRGPList(pangenome, tmpdir, complete = True):
    tmpfile = open(tmpdir.name + "/rgp.list","w")
    multi = pangenome.get_multigenics(pangenome.parameters["RGP"]["dup_margin"])#get the multigenics, in order to keep them.
    for rgp in pangenome.regions:
        if (complete and not rgp.isContigBorder) or not complete:
            tmpfile.write(' '.join([str(fam.name) for fam in rgp.families if (fam.namedPartition != "persistent" or fam in multi) ]) + "\n")
    tmpfile.flush()
    return tmpfile

def processJIMOut(pangenome, out_filename):
    outfile = open(out_filename, "r")
    mods = set()
    for line in outfile:
        geneFams = line.strip().split()[:-1]
        mods.add(Module(geneFamilies = [ pangenome.getGeneFamily(FamName) for FamName in geneFams ]))
    outfile.close()
    return mods

def getFam2Mod(modules):
    fam2mod = {}
    for mod in modules:
        for fam in mod.families:
            try:
                fam2mod[fam].add(mod)
            except KeyError:
                fam2mod[fam] = set([mod])
    return fam2mod

def groupModules(modules):
    """ groups modules that share gene families"""
    mymods = getFam2Mod(modules)
    modGraph = nx.Graph()
    bar = tqdm(mymods.values(),unit="module")
    for mods in bar:
        for mod in mods:
            modGraph.add_node(mod)
        if len(mods) > 1:
            for m1, m2 in combinations(mods, 2):
                modGraph.add_edge(m1, m2)
    # nbConnComp = len([comp for comp in nx.connected_components(modGraph)])
    # logging.getLogger().info(f"The module graph has {nx.number_of_edges(modGraph)} edges, {nx.number_of_nodes(modGraph)} nodes that are in {nbConnComp} connected components")

    newModules = set()
    for comp in nx.connected_components(modGraph):
        fams = set()
        for mod in comp:
            fams |= mod.families
        newModules.add(Module(fams))
    return newModules

def runJIM(pangenome, tmpdir, distance=95, minsize=3, minsup=6, group = True, measure = "jaccard", complete=True):
    """
        Runs a Jaccard Itemset Mining algorithm on the pangenome's RGP to extract groups of genes conserved together
        :Pangenome pangenome: a pangenome object on which to run JIM
        :str tmpdir: the temp directory path to write temporary files to
        :float distance: the maximal distance for the gene families of the same module, defaults to 0.95
        :int minsize: the minimum number of gene families to have a module, defaults to 3
        :int minsup: the minimal number of occurrence genes families need to have to be in a module, defaults to 6
        :str measure: the type of distance used. Possibilities include 'jaccard', 'hamming' and 'kulczynski'. defaults to jaccard
        :bool complete: if True, will use only RGP that are complete, i. e. bordered by persistent genes on both ends (or circular). defaults to True.
    """
    newtmpdir = tempfile.TemporaryDirectory(dir = tmpdir)

    rgp_file = writeRGPList(pangenome, newtmpdir, complete)
    out_filename = newtmpdir.name+"/module.out"

    start = time.time()
    if measure == "jaccard":
        dst = "j"
    elif measure == "hamming":
        dst = "m"
    elif measure == "kulczynski":
        dst = "k"
        # I assume the user provides 1-percentage of "mistakes" he allows for the set to exist, so that the option keeps some logic between the distance measure used.
        try:
            distance = (1 / ((100 - distance)*0.01)) * 100
        except ZeroDivisionError:  # in case we want perfect matches, distance == 100, so ...
            # should not allow any mismatch ??
            distance = len(pangenome.organisms) * 100 + 1

    # logging.getLogger().info(f"Running Jaccard Itemset Mining")
    cmd = ["jim", f"-i{dst}","-tm", f"-m{minsize}",
               f"-s-{minsup}", f"-c{distance}", rgp_file.name, out_filename]
    p = subprocess.Popen(cmd,  stderr=subprocess.DEVNULL)
    p.communicate()
    logging.getLogger().info(f"Got the results from JIM in {round(time.time() - start,4)} seconds.")

    modules = processJIMOut(pangenome, out_filename)
    # mindist = min([ mod.distance for mod in modules])
    # logging.getLogger().info(f"Minimal module distance : {round(mindist,2)}")
    if group:
        modules = groupModules(modules)
    fam2mod = getFam2Mod(modules)
    # mindist = min([ mod.distance for mod in modules])

    # logging.getLogger().info(f"Minimal module distance : {round(mindist,2)}")
    logging.getLogger().info(f"There are {len(fam2mod)} gene families in {len(modules)} modules.")

    return modules, fam2mod

def read_kofam(fname):
    kos = set()
    fam2ko = {}
    with open(fname,"r") as f:
        for line in f:
            if line.startswith('*'):#then the hit has significance
                ko = line.split()[2].upper()
                fam = line.split()[1]
                kos.add(ko)
                try:
                    fam2ko[fam].add(ko)
                except KeyError:
                    fam2ko[fam] = set([ko])

    return fam2ko, kos

def read_x2kodb(fname):
    ko2x = {}
    x2ko = {}
    with open(fname,"r") as f:
        for line in f:
            if line.startswith("path:ko"):#for the KEGG PATHWAY file
                continue
            x, ko = line.split()
            x = x.split(':')[1].upper()
            ko = ko.split(':')[1].upper()
            try:
                x2ko[x].add(ko)
            except KeyError:
                x2ko[x] = set([ko])
            try:
                ko2x[ko].add(x)
            except KeyError:
                ko2x[ko] = set([x])
    return ko2x, x2ko

def check_stats_global(ko, ko2x, x2ko, dbname):
    resulting_x = set()
    ko_not_in_x = set()
    for k in ko:
        try:
            resulting_x |= ko2x[k]
        except KeyError:
            ko_not_in_x.add(k)
    print(f"{len(ko) - len(ko_not_in_x)} KO are in KEGG {dbname}.")
    print(f"There are elements from {len(resulting_x)} {dbname.lower()} present in the pangenome")

    complete = set()
    for x in resulting_x:
        kos = x2ko[x]
        if len(kos & ko) == len(kos):
            complete.add(x)
    print(f"Among those, {len(complete)} are complete (i.e. all of the KO of the {dbname.lower()} are present in the pangenome)")

def check_partitions(fam2ko, ko2mod, ko2path):
    part2mods = {}
    part2path = {}
    for fname in ["persistent.txt","shell.txt","cloud.txt"]:
        partition = fname.split('.')[0]
        part2mods[partition] = set()
        part2path[partition] = set()
        fpart = open(fname,"r")
        nbPartMod = 0
        nbPartPath = 0
        totPart = 0
        totKo = 0
        for line in fpart:
            totPart+=1
            mykosfam = fam2ko.get(line.strip())
            if mykosfam is not None:
                totKo+=1
                for ko in mykosfam:
                    if ko in ko2mod:
                        part2mods[partition] |= ko2mod[ko]
                        nbPartMod += 1
                    if ko in ko2path:
                        part2path[partition] |= ko2path[ko]
                        nbPartPath += 1
        print(f"{totPart} families in {partition} : {totKo} have KOs. {nbPartMod} of those families are in {len(part2mods[partition])} modules and {nbPartPath} are in {len(part2path[partition])} pathways.")

def mkFam2X(pangenome, fam2ko, ko2x):
    fam2x = {}
    for fam, kos in fam2ko.items():
        for ko in kos:
            x = ko2x.get(ko)
            if x is not None:
                family = pangenome.getGeneFamily(fam)
                try:
                    fam2x[family] |= set(x)
                except KeyError:
                    fam2x[family] = set(x)
    return fam2x

def getReferences(pangenome):
    FAM2KO, KO = read_kofam("/home/abazin/Work/Datasets/CompleteColi/pang/saved_hits_KOpred.out")
    # print(f"There are {len(KO)} KO in the pangenome")
    KO2MOD, MOD2KO = read_x2kodb("/home/abazin/Work/scratch/ko_modules/ko2modules.txt")
    # check_stats_global(KO, KO2MOD, MOD2KO, "MODULES")
    KO2PATH, PATH2KO = read_x2kodb("/home/abazin/Work/scratch/ko_modules/ko2pathways.txt")
    # check_stats_global(KO, KO2PATH,PATH2KO,"PATHWAYS")
    # check_partitions(FAM2KO, KO2MOD, KO2PATH)
    fam2paths = mkFam2X(pangenome, FAM2KO, KO2PATH)
    return fam2paths

def test_pairs(fam2ref, fam2mod, modules, minsup, dist, fout, group):
    #test the pairs of relations that you get in the end.
    commonFams = set(fam2ref.keys()) & set(fam2mod.keys())
    print(f"There are {len(commonFams)} families that can be used to evaluate")

    famPairs = set()
    for mod in modules:
        commonMod = mod.families & commonFams
        if len(commonMod) >=2:
            for grp in combinations(commonMod,2):
                famPairs.add(frozenset(grp))
    print(f"There are {len(famPairs)} pairs to compare.")
    pos = 0
    neg = 0
    for pair0, pair1 in famPairs:
        refs0 = fam2ref[pair0]
        refs1 = fam2ref[pair1]
        if len(refs1 & refs0) != 0:
            pos+=1
        else:
            neg+=1
    fout.write('\t'.join(map(str, [len(modules), len(fam2mod), len(fam2ref),len(commonFams),pos, pos+neg, minsup, dist])) + "\n")
    fout.flush()
    print(f"There are {pos} / {pos+neg} pairs of families that are in the same Kegg Pathway.")


def test(fam2ref, fam2mod, modules, minsup, dist, fout, group):
    ### Testing the groups as a whole
    commonFams = set(fam2ref.keys()) & set(fam2mod.keys())
    print(f"There are {len(commonFams)} families that can be used to evaluate")

    fullPredRef = set()
    halfPredRef = set()
    for mod in modules:
        if len(mod.families & commonFams) == len(mod.families):
            #all the families of the module have a group of reference !!!
            fullPredRef.add(mod)
        if len(mod.families & commonFams) >= 2:#(len(mod.families)/2):
            halfPredRef.add(mod)

    nbCommGroup = 0
    for mod in fullPredRef:
        refInter = None
        refList = []
        for fam in mod.families:
            refs = fam2ref.get(fam)
            if refs is not None:
                if refInter is None:
                    refInter = refs
                refInter = refInter & refs
                refs = []
            refList.append(refs)
        if len(refInter) > 0:
            nbCommGroup+=1
        # print(refList)
    print(f"{nbCommGroup} / {len(fullPredRef)} modules with fully predicted gene families have a common group")

    nbCommGroup2 = 0
    for mod in halfPredRef:
        refInter = None
        refList = []
        for fam in mod.families:
            refs = fam2ref.get(fam)
            if refs is not None:
                if refInter is None:
                    refInter = refs
                refInter = refInter & refs
                refs = []
            refList.append(refs)
        if len(refInter) > 0:
            nbCommGroup2+=1
        # print(refList)
    print(f"{nbCommGroup2} / {len(halfPredRef)} modules with at least 2 predicted gene families have a common group")

    fout.write('\t'.join(map(str, [len(modules), len(fam2mod), len(fam2ref),len(commonFams),nbCommGroup, len(fullPredRef), nbCommGroup2, len(halfPredRef), minsup, dist, group])) + "\n")
    fout.flush()

def predictModules(pangenome, output, cpu, tmpdir):
    #check statuses and load info
    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, needGraph=False, needPartitions = True, needRGP = True)
    ##do the module thing
    pangenome.computeFamilyBitarrays()#might need that
    distance = 95
    fout = open("test.Rtab","w")
    # fout.write("nbModules\tnbFamsMod\tnbFamsRef\tnbFamCommun\tnbWellPredFull\tnbTotalPredFull\tnbWellPredicted2\tnbTotalPredicted2\tminsup\tdistance\tgroup\n")
    fout.write("nbModubles\tnbFamsMod\tnbFamsRef\tnbFamCommun\tnbGoodPairs\tnbTotalPairs\tminsup\tdistance\n")
    # for minsup in range(0,20):
    fam2Group = getReferences(pangenome)
    minsup = 50
    for distance in range(100,60,-1):
        print(distance, minsup)
        #not grouped
        jimModules, fam2Jim = runJIM(pangenome, tmpdir, distance = distance, minsize = 3, minsup = minsup, group=False)
        # test(fam2Group, fam2Jim, jimModules, minsup, distance, fout, group = False)
        test_pairs(fam2Group, fam2Jim, jimModules, minsup, distance, fout, group = False)
        #grouped
        # jimModules = groupModules(jimModules)
        # fam2Jim = getFam2Mod(jimModules)
        # # jimModules, fam2Jim = runJIM(pangenome, tmpdir, distance = distance, minsize = 3, minsup = minsup, group=True)
        # # test(fam2Group, fam2Jim, jimModules, minsup, distance, fout, group = True)
        # test_pairs(fam2Group, fam2Jim, jimModules, minsup, distance, fout, group = False)

    fout.close()
    #save parameters

def launch(args):
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    # mkOutdir(args.output, args.force)
    predictModules(pangenome = pangenome, output=args.output, cpu = args.cpu, tmpdir = args.tmpdir)
    #write modules to the hdf5 file

def moduleSubparser(subparser):
    parser = subparser.add_parser("module", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument("--size", required=False, type=int, default=3, help = "Minimal number of gene family in a module")
    optional.add_argument("--support", required=False, type=int, default=6, help = "Minimum number of times the module needs to be present in the pangenome to be reported")

    optional.add_argument('-o','--output', required=False, type=str, default="ppanggolin_output"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory")
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    return parser