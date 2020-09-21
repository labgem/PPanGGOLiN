#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse
import logging
import tempfile
import subprocess
import time
from multiprocessing import Pool

#installed libraries
from tqdm import tqdm

#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import write_compressed_or_not, mkOutdir
from ppanggolin.formats import checkPangenomeInfo
from ppanggolin.genetic_codes import genetic_codes


def getFamiliesToWrite(pangenome, partitionFilter):
    fams = set()
    if partitionFilter == "all":
        return set(pangenome.geneFamilies)
    if partitionFilter in ["persistent","shell","cloud"]:
        for fam in pangenome.geneFamilies:
            if fam.namedPartition == partitionFilter:
                fams.add(fam)
    elif partitionFilter in ["core","accessory"]:
        nb_org = pangenome.number_of_organisms()
        if partitionFilter == "core":
            for fam in pangenome.geneFamilies:
                if len(fam.organisms) == nb_org:
                    fams.add(fam)
        elif partitionFilter == "accessory":
            for fam in pangenome.geneFamilies:
                if len(fam.organisms) < nb_org:
                    fams.add(fam)
    return fams

def translate(seq, code):
    """ translates the given dna sequence with the given translation table"""
    # code:  https://www.bioinformatics.org/sms/iupac.html
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
                protein += 'X'  # X is for unknown
    else:
        raise IndexError(
            "Given sequence length modulo 3 was different than 0, which is unexpected.")
    return protein

def writeFastaFamilies(family, tmpdir, source, code_table):

    #have a directory for each gene family, to make deletion of tmp files simpler
    
    fname = tmpdir.name + "/" + family.name + ".fasta"
    fObj = open(fname,"w")
    
    for gene in family.genes:
        fObj.write('>' + gene.ID + "\n")
        if source == "dna":
            fObj.write(gene.dna + '\n')
        elif source == "protein":
            fObj.write(translate(gene.dna, code_table)+ "\n")
        else:
            raise Exception("Unknown sequence source given (expected 'dna' or 'protein')")
    fObj.flush()

    return fname

def launchMafft(fname, output, fam_name):
    outname = output + "/" +fam_name + ".aln"
    cmd = ["mafft", "--thread", "1", fname]
    logging.getLogger().debug("command: " + " ".join(cmd))
    subprocess.run(cmd, stdout=open(outname, "w"), stderr = subprocess.DEVNULL, check=True)#

def launchMultiMafft(args):
    launchMafft(*args)

def computeMSA(families, output, cpu, tmpdir, compress, source, code):

    newtmpdir = tempfile.TemporaryDirectory(dir = tmpdir)

    write_total = 0
    msa_total = 0
    args = []
    bar = tqdm(families, unit="family")
    code_table = genetic_codes(code)

    for family in bar:
        start_write = time.time()
        fname = writeFastaFamilies(family, newtmpdir, source, code_table)
        write_total = write_total + (time.time() - start_write)
        args.append((fname, output, family.name))
    bar.close()
    start_msa = time.time()


    bar = tqdm(range(len(families)), unit = "family")
    with Pool(cpu) as p:
        for _ in p.imap_unordered(launchMultiMafft, args):
            bar.update()
    bar.close()

    msa_total = msa_total + (time.time() - start_msa)
    

def writeMSAFiles(pangenome, output, cpu = 1, partition = "core",compress = False, tmpdir = "/tmp", source="protein", force=False):
    
    needPartitions = False
    if partition in ["persistent","shell","cloud"]:
        needPartitions = True

    outname = output + f"/msa_{partition}_{source}/"
    mkOutdir(outname, force=force)

    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, needPartitions= needPartitions, needGeneSequences=True)
    logging.getLogger().info(f"computing msa for {partition} families...")
    families = getFamiliesToWrite(pangenome, partitionFilter=partition)

    #this must exist since we loaded the pangenome and families are required
    code = pangenome.parameters["cluster"]["translation_table"]

    computeMSA(families, outname, cpu=cpu, compress=compress, tmpdir=tmpdir, source=source, code = code)
    logging.getLogger().info(f"Done writing all {partition} MSA in: {outname}")

def launchMSA(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    writeMSAFiles(pangenome, args.output, cpu=args.cpu, partition = args.partition, compress=args.compress, tmpdir=args.tmpdir, source=args.source, force=args.force)

def writeMSASubparser(subparser):
    parser = subparser.add_parser("msa", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group(title = "Required arguments", description = "The following arguments are required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o','--output', required=True, type=str, help="Output directory where the file(s) will be written")

    optional = parser.add_argument_group(title = "Optional arguments. Indicating 'all' writes all elements. Writing a partition ('persistent', 'shell', 'cloud', 'core' or 'accessory') write the elements associated to said partition.")
    ##could make choice to allow customization
    optional.add_argument("--partition", required=False, default="core", choices=["all","persistent","shell","cloud","core","accessory"], help = "compute Multiple Sequence Alignement of the gene families in the given partition")
    optional.add_argument("--compress",required=False, action="store_true",help="Compress the files in .gz")
    optional.add_argument("--source",required=False, default = "protein", choices = ["dna","protein"], help = "indicates whether to use protein or dna sequences to compute the msa")
    return parser
