#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import tempfile
import subprocess
import argparse
from collections import defaultdict

# local libraries
from ppanggolin.formats import check_pangenome_info
from ppanggolin.utils import mkOutdir, read_compressed_or_not
from ppanggolin.pangenome import Pangenome
from ppanggolin.figures.draw_spot import drawSelectedSpots, subgraph


def createdb(fileObj, tmpdir):
    """
    Create a MMseqs2 sequence database with the given fasta file

    :param fileObj: Fasta file
    :type fileObj: _io.TextIOWrapper
    :param tmpdir: temporary directory
    :type tmpdir: tempfile.TemporaryDirectory

    :return: DB file
    :rtype: _io.TextIOWrapper
    """
    seqdb = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name)
    cmd = ["mmseqs", "createdb", fileObj.name, seqdb.name, '--dbtype', '0']
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    return seqdb


def alignSeqToPang(pangFile, seqFile, output, tmpdir, cpu=1, no_defrag=False, identity=0.8, coverage=0.8):
    pang_db = createdb(pangFile, tmpdir)
    seq_db = createdb(seqFile, tmpdir)
    cov_mode = "0"  # coverage of query and target
    if not no_defrag:
        cov_mode = "1"  # coverage of target
    aln_db = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name)
    cmd = ["mmseqs", "search", seq_db.name, pang_db.name, aln_db.name, tmpdir.name, "-a", "--min-seq-id", str(identity),
           "-c", str(coverage), "--cov-mode", cov_mode, "--threads", str(cpu)]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Aligning sequences to cluster representatives...")
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    outfile = output + "/input_to_pangenome_associations.blast-tab_tmp"#write a tmp file of the results
    cmd = ["mmseqs", "convertalis", seq_db.name, pang_db.name, aln_db.name, outfile, "--format-mode", "2"]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Extracting alignments...")
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    pang_db.close()
    seq_db.close()
    aln_db.close()

    return outfile


def readAlignments(outfile, pangenome):
    seq2pang = {}
    outname = open(outfile.replace("_tmp",""),"w")#write the actual result file
    with open(outfile, "r") as alnFile:
        for line in alnFile:
            line = line.replace("ppanggolin_","")#remove the 'ppanggolin_' bit of the id
            outname.write(line)
            line = line.split()
            if seq2pang.get(line[0]) is None:  # if no results were found yet
                seq2pang[line[0]] = pangenome.getGeneFamily(line[1])  # then the best hit is the first one we see.
    outname.close()
    return seq2pang, outname.name


def getSeq(seqFile):
    seqset = set()
    for line in seqFile:
        if line.startswith(">"):
            seqset.add(line[1:])
    return seqset


def writeGeneFamSequences(pangenome, fileObj, add=""):
    for fam in pangenome.geneFamilies:
        fileObj.write(">" + add + fam.name + "\n")
        fileObj.write(fam.sequence + "\n")
    fileObj.flush()


def projectPartition(seq2pang, seqSet, output):
    partitionProj = output + "/sequences_partition_projection.tsv"
    with open(partitionProj, "w") as partProjFile:
        for key, pangFam in seq2pang.items():
            partProjFile.write(key + "\t" + pangFam.namedPartition + "\n")
        for remainingSeq in (seq2pang.keys() & seqSet):
            partProjFile.write(remainingSeq + "\tcloud\n")  # if there is no hit, it's going to be cloud genes.
    return partitionProj


def getFam2RGP(pangenome, multigenics):
    """associates families to the RGP they belong to, and those they are bordering"""
    fam2rgp = defaultdict(list)
    for rgp in pangenome.regions:
        for fam in rgp.families:
            fam2rgp[fam].append(rgp.name)
        for fam in [gene.family for border in
                    rgp.getBorderingGenes(pangenome.parameters["spots"]["set_size"], multigenics) for gene in border]:
            fam2rgp[fam].append(rgp.name)
    return fam2rgp


def getFam2spot(pangenome, multigenics):
    """
    reads a pangenome object and returns a dictionary of family to RGP and family to spot,
    that indicates where each family is
    """
    # those are to be replaced as spots should be stored in the pangenome, and in the h5.
    fam2spot = defaultdict(list)
    fam2border = defaultdict(list)
    for spot in pangenome.spots:
        fams = set()
        famsBorder = set()
        for rgp in spot.regions:
            fams |= rgp.families
            famsBorder |= set(
                [gene.family for border in rgp.getBorderingGenes(pangenome.parameters["spots"]["set_size"], multigenics)
                 for gene in border])
        for fam in fams:
            fam2spot[fam].append(spot)
        for fam in famsBorder:
            fam2border[fam].append(spot)
    return fam2spot, fam2border, multigenics


def add_spot_str(a):
    return "spot_" + str(a.ID)


def draw_spot_gexf(spots, output, multigenics, fam2mod, set_size=3):
    for spot in spots:
        fname = output + "/spot_" + str(spot.ID) + ".gexf"
        subgraph(spot, fname, set_size=set_size, multigenics=multigenics, fam2mod=fam2mod)


def getSeqInfo(seq2pang, pangenome, output, draw_related, disable_bar=False):
    logging.getLogger().info("Writing RGP and spot information related to hits in the pangenome")
    multigenics = pangenome.get_multigenics(pangenome.parameters["RGP"]["dup_margin"])

    finfo = open(output + "/info_input_seq.tsv", "w")
    finfo.write("input\tfamily\tpartition\tspot_list_as_member\tspot_list_as_border\trgp_list\n")
    fam2rgp = getFam2RGP(pangenome, multigenics)
    fam2spot, fam2border, multigenics = getFam2spot(pangenome, multigenics)
    spot_list = set()
    for seq, panfam in seq2pang.items():
        finfo.write(seq + '\t' + panfam.name + "\t" + panfam.namedPartition + "\t" + ",".join(
            map(add_spot_str, fam2spot[panfam])) + "\t" + ",".join(
            map(add_spot_str, fam2border[panfam])) + "\t" + ','.join(fam2rgp[panfam]) + "\n")
        spot_list |= set(fam2spot[panfam])
        spot_list |= set(fam2border[panfam])
    finfo.close()
    if draw_related:
        drawn_spots = set()
        for spot in spot_list:
            if len(spot.getUniqOrderedSet()) > 1:
                drawn_spots.add(spot)
        logging.getLogger().info(
            f"Drawing the {len(drawn_spots)} spots with more than 1 organization "
            f"related to hits of the input sequences...")
        drawSelectedSpots(drawn_spots, pangenome, output, pangenome.parameters["spots"]["overlapping_match"],
                          pangenome.parameters["spots"]["exact_match"], pangenome.parameters["spots"]["set_size"],
                          disable_bar=disable_bar)
        # fam2module
        fam2mod = {}
        if pangenome.status["modules"] != "No":
            for mod in pangenome.modules:
                for fam in mod.families:
                    fam2mod[fam] = f"module_{mod.ID}"

        draw_spot_gexf(drawn_spots, output, multigenics=multigenics, fam2mod=fam2mod)

    logging.getLogger().info(
        f"File listing RGP and spots where sequences of interest are located : '{output + '/info_input_seq.tsv'}'")


def get_seq2pang(pangenome, sequenceFile, output, tmpdir, cpu=1, no_defrag=False, identity=0.8, coverage=0.8):
    """
    Assign a pangenome gene family to the input sequences.

    :param pangenome: Pangenome with gene families to align with the given input sequences
    :type pangenome: Pangenome
    :param sequenceFile: Sequences in a .fasta file to align with the given Pangenome
    :type sequenceFile: str
    :param output: Output directory
    :type output: str
    :param tmpdir: Temporary directory
    :type tmpdir: tempfile.TemporaryDirectory
    :param cpu: number of CPU cores to use 
    :type cpu: int
    :param no_defrag: do not use the defrag workflow if true
    :type no_defrag: Boolean
    :param identity: minimal identity threshold for the alignment
    :type identity: float
    :param coverage: minimal identity threshold for the alignment
    :type coverage: float

    :return: sequence set, blast-tab result file string, and sequences aligned with families
    :rtype: set, str, dic
    """
    tmpPangFile = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name)

    writeGeneFamSequences(pangenome, tmpPangFile, add="ppanggolin_")

    with read_compressed_or_not(sequenceFile) as seqFileObj:
        seqSet = getSeq(seqFileObj)
        alignFile = alignSeqToPang(tmpPangFile, seqFileObj, output, tmpdir, cpu, no_defrag, identity, coverage)

    seq2pang, alignFile = readAlignments(alignFile, pangenome)

    tmpPangFile.close()

    return seqSet, alignFile, seq2pang


def align(pangenome, sequenceFile, output, tmpdir, identity=0.8, coverage=0.8, no_defrag=False, cpu=1, getinfo=False,
          draw_related=False, disable_bar=False):
    if pangenome.status["geneFamilySequences"] not in ["inFile", "Loaded", "Computed"]:
        raise Exception("Cannot use this function as your pangenome does not have gene families representatives "
                        "associated to it. For now this works only if the clustering is realised by PPanGGOLiN.")
    # could be possible either by picking a representative somehow, or by aligning on genes rather than on
    # families, if they are in the pangenome.

    if getinfo or draw_related:
        need_mod = False
        if pangenome.status["modules"] != "No":
            # modules are not required to be loaded, but if they have been computed we load them.
            need_mod = True
        check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_partitions=True, need_rgp=True,
                             need_spots=True, need_modules=need_mod, disable_bar=disable_bar)
    else:
        check_pangenome_info(pangenome, need_families=True, disable_bar=disable_bar)

    new_tmpdir = tempfile.TemporaryDirectory(dir=tmpdir)

    seqSet, alignFile, seq2pang = get_seq2pang(pangenome, sequenceFile, output, new_tmpdir,
                                               cpu, no_defrag, identity, coverage)

    if getinfo or draw_related:
        getSeqInfo(seq2pang, pangenome, output, draw_related, disable_bar=disable_bar)
    partProj = projectPartition(seq2pang, seqSet, output)  # write the partition assignation only
    logging.getLogger().info(f"sequences partition projection : '{partProj}'")
    logging.getLogger().info(f"{len(seq2pang)} sequences over {len(seqSet)} have at least one hit in the pangenome.")
    logging.getLogger().info(f"Blast-tab file of the alignment : '{alignFile}'")

    new_tmpdir.cleanup()


def launch(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    if args.interest or args.fig_margin or args.label_priority:
        logging.getLogger().warning("Options --interest, --fig_margin and --label_priority are deprecated, "
                                    "and the actions they defined are now doable directly in the interactive figures "
                                    "that are drawn")
    align(pangenome=pangenome, sequenceFile=args.sequences, output=args.output, tmpdir=args.tmpdir, cpu=args.cpu,
            identity=args.identity, coverage=args.coverage, no_defrag=args.no_defrag, getinfo=args.getinfo,
            draw_related=args.draw_related, disable_bar=args.disable_prog_bar)


def alignSubparser(subparser):
    parser = subparser.add_parser("align", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-S', '--sequences', required=True, type=str,
                          help="sequences (nucleotides or amino acids) to align on the pangenome gene families")

    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o', '--output', required=True, type=str,
                          help="Output directory where the file(s) will be written")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--defrag", required=False, action="store_true",
                          help=argparse.SUPPRESS)  # This ensures compatibility with the old option "defrag"
    optional.add_argument('--no_defrag', required=False, action="store_true",
                          help="DO NOT Realign gene families to link fragments with"
                               "their non-fragmented gene family. (default: False)")
    optional.add_argument('--identity', required=False, type=float, default=0.5,
                          help="min identity percentage threshold")
    optional.add_argument('--coverage', required=False, type=float, default=0.8,
                          help="min coverage percentage threshold")
    optional.add_argument("--translation_table", required=False, default="11",
                          help="Translation table (genetic code) to use.")
    optional.add_argument("--getinfo", required=False, action="store_true",
                          help="Use this option to extract info related to the best hit of each query, "
                               "such as the RGP it is in, or the spots.")
    optional.add_argument("--draw_related", required=False, action="store_true",
                          help="Draw figures and provide graphs in a gexf format of the eventual spots"
                               " associated to the input sequences")
    optional.add_argument("--interest", required=False, action="store_true",
                          help=argparse.SUPPRESS)  # This ensures compatibility with the old API
    # but does not use the option
    optional.add_argument("--fig_margin", required=False, action="store_true",
                          help=argparse.SUPPRESS)  # This ensures compatibility with the old API
    # but does not use the option
    optional.add_argument("--label_priority", required=False, action="store_true",
                          help=argparse.SUPPRESS)  # This ensures compatibility with the old API
    # but does not use the option
    optional.add_argument("--use_pseudo", required=False, action="store_true",
                          help="In the context of provided annotation, use this option to read pseudogenes. "
                               "(Default behavior is to ignore them)")
    return parser
