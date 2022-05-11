#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
import tempfile
import subprocess
import time
from multiprocessing import get_context

# installed libraries
from tqdm import tqdm

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mkOutdir, restricted_float
from ppanggolin.formats import checkPangenomeInfo
from ppanggolin.genetic_codes import genetic_codes



def is_single_copy(fam, dup_margin):
    nb_multi = 0
    for gene_list in fam.getOrgDict().values():
        if len(gene_list) > 1:
            nb_multi += 1
    dup_ratio = nb_multi / len(fam.organisms)
    if dup_ratio < dup_margin:
        return True
    return False


def getFamiliesToWrite(pangenome, partitionFilter, soft_core=0.95, dup_margin=0.95, single_copy=True):
    fams = set()
    nb_org = pangenome.number_of_organisms()

    if partitionFilter == "all":
        return set(pangenome.geneFamilies)
    if partitionFilter in ["persistent", "shell", "cloud"]:
        for fam in pangenome.geneFamilies:
            if fam.namedPartition == partitionFilter:
                if single_copy and is_single_copy(fam, dup_margin):
                    fams.add(fam)
                elif not single_copy:
                    fams.add(fam)
    elif partitionFilter in ["core", "accessory", "softcore"]:
        if partitionFilter == "core":
            for fam in pangenome.geneFamilies:
                if len(fam.organisms) == nb_org:
                    if single_copy and is_single_copy(fam, dup_margin):
                        fams.add(fam)
                    elif not single_copy:
                        fams.add(fam)
        elif partitionFilter == "accessory":
            for fam in pangenome.geneFamilies:
                if len(fam.organisms) < nb_org:
                    if single_copy and is_single_copy(fam, dup_margin):
                        fams.add(fam)
                    elif not single_copy:
                        fams.add(fam)
        elif partitionFilter == "softcore":
            for fam in pangenome.geneFamilies:
                if len(fam.organisms) >= nb_org * soft_core:
                    if single_copy and is_single_copy(fam, dup_margin):
                        fams.add(fam)
                    elif not single_copy:
                        fams.add(fam)
    return fams


def translate(seq, code):
    """ translates the given dna sequence with the given translation table"""
    # code:  https://www.bioinformatics.org/sms/iupac.html
    start_table = code["start_table"]
    table = code["trans_table"]

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


def writeFastaFamilies(family, tmpdir, source, use_gene_id, code_table):
    fname = tmpdir.name + "/" + family.name + ".fasta"
    fObj = open(fname, "w")
    # get genes that are present in only one copy for our family in each organism.
    single_copy_genes = []
    for _, genes in family.getOrgDict().items():
        if len(genes) == 1:
            single_copy_genes.extend(genes)

    for gene in single_copy_genes:
        if use_gene_id:
            fObj.write('>' + gene.ID + "\n")
        else:
            fObj.write('>' + gene.organism.name + "\n")
        if source == "dna":
            fObj.write(gene.dna + '\n')
        elif source == "protein":
            fObj.write(translate(gene.dna, code_table) + "\n")
        else:
            raise Exception("Unknown sequence source given (expected 'dna' or 'protein')")
    fObj.flush()

    return fname


def launchMafft(fname, output, fam_name):
    outname = output + "/" + fam_name + ".aln"
    cmd = ["mafft", "--thread", "1", fname]
    logging.getLogger().debug("command: " + " ".join(cmd))
    subprocess.run(cmd, stdout=open(outname, "w"), stderr=subprocess.DEVNULL, check=True)


def launchMultiMafft(args):
    launchMafft(*args)


def computeMSA(families, output, cpu, tmpdir, source, use_gene_id, code, disable_bar=False):
    newtmpdir = tempfile.TemporaryDirectory(dir=tmpdir)

    write_total = 0
    args = []
    logging.getLogger().info("Preparing input files for MSA...")
    bar = tqdm(families, unit="family", disable=disable_bar)
    code_table = genetic_codes(code)

    for family in bar:
        start_write = time.time()
        fname = writeFastaFamilies(family, newtmpdir, source, use_gene_id, code_table)
        write_total = write_total + (time.time() - start_write)
        args.append((fname, output, family.name))
    bar.close()

    logging.getLogger().info("Computing the MSA ...")
    bar = tqdm(range(len(families)), unit="family", disable=disable_bar)
    with get_context('fork').Pool(cpu) as p:
        for _ in p.imap_unordered(launchMultiMafft, args):
            bar.update()
    bar.close()


def writeWholeGenomeMSA(pangenome, families, phylo_name, outname, use_gene_id=False):
    phyloDict = {}
    for org in pangenome.organisms:
        phyloDict[org.name] = ""
    for fam in families:
        missing_genomes = set(phyloDict.keys())
        fin = open(outname + "/" + fam.name + ".aln", "r")
        genome_id = ""
        seq = ""
        curr_len = 0
        dup_gene = 0
        curr_phyloDict = {}

        for line in fin:
            if line.startswith('>'):
                if genome_id != "":
                    if genome_id not in missing_genomes:
                        dup_gene += 1
                        # duplicated genes. Replacing them with gaps.
                        curr_phyloDict[genome_id] = "-" * curr_len
                    else:
                        curr_phyloDict[genome_id] = seq
                        missing_genomes -= {genome_id}
                        curr_len = len(seq)
                if use_gene_id:
                    genome_id = pangenome.getGene(line[1:].strip()).organism.name
                else:
                    genome_id = line[1:].strip()
                seq = ""
            else:
                seq += line.strip()
        if genome_id != "":
            if genome_id not in missing_genomes:
                # duplicated genes. Replacing them with gaps.
                curr_phyloDict[genome_id] = "-" * curr_len
            else:
                curr_phyloDict[genome_id] = seq
                curr_len = len(seq)
        fin.close()
        for genome in missing_genomes:
            curr_phyloDict[genome] = "-" * curr_len

        for key, val in curr_phyloDict.items():
            phyloDict[key] += val

    fout = open(phylo_name, "w")
    for key, val in phyloDict.items():
        fout.write(">" + key + "\n")
        fout.write(val + "\n")
    fout.close()


def writeMSAFiles(pangenome, output, cpu=1, partition="core", tmpdir="/tmp", source="protein", soft_core=0.95,
                  phylo=False, use_gene_id=False, translation_table="11", dup_margin = 0.95, single_copy=True, force=False, disable_bar=False):
    needPartitions = False
    if partition in ["persistent", "shell", "cloud"]:
        needPartitions = True

    outname = output + f"/msa_{partition}_{source}/"
    mkOutdir(outname, force=force)

    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, needPartitions=needPartitions,
                       needGeneSequences=True, disable_bar=disable_bar)
    logging.getLogger().info(f"Doing MSA for {partition} families...")
    families = getFamiliesToWrite(pangenome, partitionFilter=partition, soft_core=soft_core, dup_margin=dup_margin, single_copy=single_copy)

    #check that the code is similar than the one used previously, if there is one
    if 'translation_table' in pangenome.parameters["cluster"]:
        if pangenome.parameters["cluster"]["translation_table"] != translation_table:
            logging.getLogger().warning(f"The translation table used during clustering ('{pangenome.parameters['cluster']['translation_table']}') is different than the one provided now ('{translation_table}')")
    code = translation_table

    computeMSA(families, outname, cpu=cpu, tmpdir=tmpdir, source=source, use_gene_id=use_gene_id, code=code,
               disable_bar=disable_bar)
    logging.getLogger().info(f"Done writing all {partition} MSA in: {outname}")

    if phylo:
        logging.getLogger().info("Writing the whole genome msa file")
        if partition == "softcore":
            phylo_name = output + f"/{partition}_{soft_core}_genome_alignment.aln"
        else:
            phylo_name = output + f"/{partition}_genome_alignment.aln"
        writeWholeGenomeMSA(pangenome, families, phylo_name, outname, use_gene_id=use_gene_id)
        logging.getLogger().info(f"Done writing the {partition} genome alignment in: '{phylo_name}'")


def launchMSA(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    writeMSAFiles(pangenome, args.output, cpu=args.cpu, partition=args.partition, tmpdir=args.tmpdir,
                  source=args.source, soft_core=args.soft_core, phylo=args.phylo, use_gene_id=args.use_gene_id,
                  translation_table=args.translation_table, dup_margin=args.dup_margin,
                  single_copy=args.single_copy, force=args.force, disable_bar=args.disable_prog_bar)


def writeMSASubparser(subparser):
    parser = subparser.add_parser("msa", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group(title="Required arguments",
                                         description="The following arguments are required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o', '--output', required=True, type=str,
                          help="Output directory where the file(s) will be written")

    optional = parser.add_argument_group(title="Optional arguments. Indicating 'all' writes all elements. "
                                               "Writing a partition ('persistent', 'shell', 'cloud', 'core' or "
                                               "'accessory') write the elements associated to said partition.")
    # could make choice to allow customization
    optional.add_argument("--soft_core", required=False, type=restricted_float, default=0.95,
                          help="Soft core threshold to use if 'softcore' partition is chosen")
    optional.add_argument("--dup_margin", required=False, type=restricted_float, default=0.05,
                          help="minimum ratio of organisms in which the family must have multiple genes "
                               "for it to be considered 'duplicated'")
    optional.add_argument("--single_copy", required=False, action="store_true", default=False,
                          help="Use report gene families that are considered 'single copy', for details see "
                               "option --dup_margin")
    optional.add_argument("--partition", required=False, default="core",
                          choices=["all", "persistent", "shell", "cloud", "core", "accessory", 'softcore'],
                          help="compute Multiple Sequence Alignement of the gene families in the given partition")
    optional.add_argument("--source", required=False, default="protein", choices=["dna", "protein"],
                          help="indicates whether to use protein or dna sequences to compute the msa")
    optional.add_argument("--phylo", required=False, action='store_true',
                          help="Writes a whole genome msa file for additional phylogenetic analysis")
    optional.add_argument("--use_gene_id", required=False, action='store_true',
                          help="Use gene identifiers rather than organism names for sequences in the family MSA"
                               " (organism names are used by default)")
    optional.add_argument("--translation_table", required=False, default="11",
                          help="Translation table (genetic code) to use.")
    return parser
