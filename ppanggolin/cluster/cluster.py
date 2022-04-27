#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import tempfile
import subprocess
from collections import defaultdict
import os
import argparse

# installed libraries
from networkx import Graph
from tqdm import tqdm

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Gene
from ppanggolin.utils import read_compressed_or_not, restricted_float
from ppanggolin.formats.writeBinaries import write_pangenome, erase_pangenome
from ppanggolin.formats.readBinaries import check_pangenome_info, get_gene_sequences_from_file
from ppanggolin.formats.writeSequences import write_gene_sequences_from_annotations


def align_rep(faa_file, tmpdir, cpu, coverage, identity):
    seqdb = tmpdir.name + '/rep_sequence_db'
    cmd = ["mmseqs", "createdb", faa_file, seqdb]
    logging.getLogger().debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    alndb = tmpdir.name + '/rep_alignment_db'
    cmd = ["mmseqs", "search", seqdb, seqdb, alndb, tmpdir.name, "-a", "--min-seq-id", str(identity), "-c",
           str(coverage), "--cov-mode", "1", "--threads", str(cpu)]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Aligning cluster representatives...")
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    outfile = tmpdir.name + '/rep_families.tsv'
    cmd = ["mmseqs", "convertalis", seqdb, seqdb, alndb, outfile, "--format-output", "query,target,qlen,tlen,bits"]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Extracting alignments...")
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    return outfile


def first_clustering(sequences, tmpdir, cpu, code, coverage, identity, mode):
    seq_nucdb = tmpdir.name + '/nucleotid_sequences_db'
    cmd = ["mmseqs", "createdb", sequences.name, seq_nucdb]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Creating sequence database...")
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    seqdb = tmpdir.name + '/aa_db'
    cmd = ["mmseqs", "translatenucs", seq_nucdb, seqdb, "--threads", str(cpu), "--translation-table", code]
    logging.getLogger().debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    cludb = tmpdir.name + '/cluster_db'
    cmd = ["mmseqs", "cluster", seqdb, cludb, tmpdir.name, "--cluster-mode", mode, "--min-seq-id", str(identity), "-c",
           str(coverage), "--threads", str(cpu), "--kmer-per-seq", "80", "--max-seqs", "300"]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Clustering sequences...")
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    logging.getLogger().info("Extracting cluster representatives...")
    repdb = tmpdir.name + '/representative_db'
    cmd = ["mmseqs", "result2repseq", seqdb, cludb, repdb]
    logging.getLogger().debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    reprfa = tmpdir.name + '/representative_sequences.fasta'
    cmd = ["mmseqs", "result2flat", seqdb, seqdb, repdb, reprfa, "--use-fasta-header"]
    logging.getLogger().debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    outtsv = tmpdir.name + '/families_tsv'
    cmd = ["mmseqs", "createtsv", seqdb, seqdb, cludb, outtsv, "--threads", str(cpu), "--full-header"]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Writing gene to family informations")
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    return reprfa, outtsv


def read_faa(faa_file_name):
    fam2seq = {}
    head = ""
    with open(faa_file_name, "r") as faFile:
        for line in faFile:
            if line.startswith('>'):
                head = line[1:].strip().replace("ppanggolin_", "")  # remove the eventual addition
            else:
                fam2seq[head] = line.strip()
    return fam2seq


def read_tsv(tsv_file_name):
    """
    reading tsv file
    :param tsv_file_name: path to the tsv
    :type tsv_file_name: str

    :return: two dictionnary which link genes and families
    :rtype: tuple(dict, dict)
    """
    genes2fam = {}
    fam2genes = defaultdict(set)
    with open(tsv_file_name, "r") as tsvfile:
        for line in tsvfile:
            line = line.replace('"', '').replace("ppanggolin_", "").split()
            # remove the '"' char which protects the fields, and the eventual addition
            genes2fam[line[1]] = (line[0], False)  # fam id, and it's a gene (and not a fragment)
            fam2genes[line[0]].add(line[1])
    return genes2fam, fam2genes


def refine_clustering(tsv, aln_file, fam_to_seq):
    simgraph = Graph()
    genes2fam, fam2genes = read_tsv(tsv)
    logging.getLogger().info(f"Starting with {len(fam_to_seq)} families")
    # create the nodes
    for fam, genes in fam2genes.items():
        simgraph.add_node(fam, nbgenes=len(genes))
    # add the edges
    with open(aln_file, "r") as alnfile:
        for line in alnfile:
            line = line.replace('"', '').replace("ppanggolin_", "").split()  # remove the eventual addition

            if line[0] != line[1]:
                simgraph.add_edge(line[0], line[1], score=float(line[4]))
                simgraph.nodes[line[0]]["length"] = int(line[2])
                simgraph.nodes[line[1]]["length"] = int(line[3])
    for node, nodedata in simgraph.nodes(data=True):
        choice = (None, 0, 0, 0)
        for neighbor in simgraph.neighbors(node):
            nei = simgraph.nodes[neighbor]
            score = simgraph[neighbor][node]["score"]
            if nei["length"] > nodedata["length"] and nei["nbgenes"] >= nodedata["nbgenes"] and choice[3] < score:
                choice = (genes2fam[neighbor][0], nei["length"], nei["nbgenes"], score)
                # `genes2fam[neighbor]` instead of just neighbor in case that family has been assigned already
                # (this is for smaller fragments that are closer to other fragments than the actual gene family)
        if choice[0] is not None:
            genestochange = fam2genes[node]
            for gene in genestochange:
                genes2fam[gene] = (choice[0], True)
                fam2genes[choice[0]].add(gene)
            del fam2genes[node]
    new_fam_to_seq = {}
    for fam in fam2genes:
        new_fam_to_seq[fam] = fam_to_seq[fam]
    logging.getLogger().info(f"Ending with {len(new_fam_to_seq)} gene families")
    return genes2fam, new_fam_to_seq


def read_gene2fam(pangenome, gene2fam, disable_bar=False):
    logging.getLogger().info(f"Adding {len(gene2fam)} genes to the gene families")

    link = True if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"] else False
    if link:
        if len(gene2fam) != len(pangenome.genes):  # then maybe there are genes with identical IDs
            raise Exception("Something unexpected happened during clustering "
                            "(have less genes clustered than genes in the pangenome). "
                            "A probable reason is that two genes in two different organisms have the same IDs; "
                            "If you are sure that all of your genes have non identical IDs, "
                            "please post an issue at https://github.com/labgem/PPanGGOLiN/")
    bar = tqdm(gene2fam.items(), unit="gene", disable=disable_bar)
    for gene, (family, is_frag) in bar:
        fam = pangenome.add_gene_family(family)
        if link:  # doing the linking if the annotations are loaded.
            gene_obj = pangenome.get_gene(gene)
        else:
            gene_obj = Gene(gene)
        gene_obj.is_fragment = is_frag
        fam.add_gene(gene_obj)
    bar.close()


def read_fam2seq(pangenome, fam2seq):
    logging.getLogger().info("Adding protein sequences to the gene families")
    for family, protein in fam2seq.items():
        fam = pangenome.add_gene_family(family)
        fam.add_sequence(protein)


def check_pangenome_former_clustering(pangenome, force):
    """ checks pangenome status and .h5 files for former clusterings, delete them if allowed or raise an error """
    if pangenome.status["genesClustered"] == "inFile" and not force:
        raise Exception("You are trying to cluster genes that are already clustered together. If you REALLY want to "
                        "do that, use --force (it will erase everything except annotation data in your HDF5 file!)")
    elif pangenome.status["genesClustered"] == "inFile" and force:
        erase_pangenome(pangenome, gene_families=True)


def check_pangenome_for_clustering(pangenome, tmp_file, force, disable_bar=False):
    """
        Check the pangenome statuses and write the gene sequences in the provided tmpFile.
        (whether they are written in the .h5 file or currently in memory)
    """
    check_pangenome_former_clustering(pangenome, force)
    if pangenome.status["geneSequences"] in ["Computed", "Loaded"]:
        # we append the gene ids by 'ppanggolin' to avoid crashes from mmseqs when sequence IDs are only numeric.
        write_gene_sequences_from_annotations(pangenome, tmp_file, add="ppanggolin_", disable_bar=disable_bar)
    elif pangenome.status["geneSequences"] == "inFile":
        get_gene_sequences_from_file(pangenome.file, tmp_file, add="ppanggolin_",
                                     disable_bar=disable_bar)  # write CDS sequences to the tmpFile
    else:
        tmp_file.close()  # closing the tmp file since an exception will be raised.
        raise Exception("The pangenome does not include gene sequences, thus it is impossible to cluster "
                        "the genes in gene families. Either provide clustering results (see --clusters), "
                        "or provide a way to access the gene sequence during the annotation step "
                        "(having the fasta in the gff files, or providing the fasta files through the --fasta option)")


def infer_singletons(pangenome):
    """creates a new family for each gene with no associated family"""
    singleton_counter = 0
    for gene in pangenome.genes:
        if gene.family is None:
            pangenome.add_gene_family(gene.ID).add_gene(gene)
            singleton_counter += 1
    logging.getLogger().info(f"Inferred {singleton_counter} singleton families")


def clustering(pangenome, tmpdir, cpu, defrag=True, code="11", coverage=0.8, identity=0.8, mode="1", force=False,
               disable_bar=False):
    newtmpdir = tempfile.TemporaryDirectory(dir=tmpdir)
    sequence_file = open(newtmpdir.name + '/nucleotid_sequences', "w")

    check_pangenome_for_clustering(pangenome, sequence_file, force, disable_bar=disable_bar)
    logging.getLogger().info("Clustering all of the genes sequences...")
    rep, tsv = first_clustering(sequence_file, newtmpdir, cpu, code, coverage, identity, mode)

    sequence_file.close()
    fam2seq = read_faa(rep)
    if not defrag:
        logging.getLogger().debug("No defragmentation")
        genes2fam = read_tsv(tsv)[0]
    else:
        logging.getLogger().info("Associating fragments to their original gene family...")
        aln = align_rep(rep, newtmpdir, cpu, coverage, identity)
        genes2fam, fam2seq = refine_clustering(tsv, aln, fam2seq)
        pangenome.status["defragmented"] = "Computed"
    newtmpdir.cleanup()
    read_fam2seq(pangenome, fam2seq)
    read_gene2fam(pangenome, genes2fam, disable_bar=disable_bar)

    pangenome.status["genesClustered"] = "Computed"
    pangenome.status["geneFamilySequences"] = "Computed"

    pangenome.parameters["cluster"] = {}
    pangenome.parameters["cluster"]["coverage"] = coverage
    pangenome.parameters["cluster"]["identity"] = identity
    pangenome.parameters["cluster"]["defragmentation"] = defrag
    pangenome.parameters["cluster"]["translation_table"] = code
    pangenome.parameters["cluster"]["read_clustering_from_file"] = False


def mk_local_to_gene(pangenome):
    """
        Creates a dictionary that stores local identifiers, if all local identifiers are unique (and if they exist)
    """
    local_dict = {}
    for gene in pangenome.genes:
        old_len = len(local_dict)
        local_dict[gene.local_identifier] = gene
        if len(local_dict) == old_len:
            if pangenome.parameters["annotation"]["read_annotations_from_file"] and not \
                    pangenome.parameters["annotation"]["used_local_identifiers"]:
                raise Exception(f"'{gene.local_identifier}' was found multiple times used as an identifier. "
                                f"The identifier of the genes (locus_tag, protein_id in gbff, ID in gff) were not "
                                f"unique throughout all of the files. It is thus impossible to differentiate the genes."
                                f" To use this function while importing annotation, all identifiers MUST be unique "
                                f"throughout all of your genomes")
            return {}  # local identifiers are not unique.
    return local_dict


def read_clustering(pangenome, families_tsv_file, infer_singleton=False, force=False, disable_bar=False):
    """
        Creates the pan, the gene families and the genes with an associated gene family.
        Reads a families tsv file from mmseqs2 output and adds the gene families and the genes to the pan.
    """
    check_pangenome_former_clustering(pangenome, force)
    check_pangenome_info(pangenome, need_annotations=True, disable_bar=disable_bar)

    logging.getLogger().info("Reading " + families_tsv_file + " the gene families file ...")
    filesize = os.stat(families_tsv_file).st_size
    families_tsv_file = read_compressed_or_not(families_tsv_file)
    frag = False
    # the genome annotations are necessarily loaded.
    nb_gene_with_fam = 0
    local_dict = mk_local_to_gene(pangenome)
    bar = tqdm(total=filesize, unit="bytes", disable=disable_bar)
    line_counter = 0
    for line in families_tsv_file:
        line_counter += 1
        bar.update(len(line))
        try:
            elements = [el.strip() for el in line.split()]  # 2 or 3 fields expected
            if len(elements) <= 1:
                raise ValueError("No tabulation separator found in gene families file")
            (fam_id, gene_id, is_frag) = elements if len(elements) == 3 else elements + [None]
            try:
                gene_obj = pangenome.get_gene(gene_id)
            except KeyError:
                gene_obj = local_dict.get(gene_id)
            if gene_obj is not None:
                nb_gene_with_fam += 1
                fam = pangenome.add_gene_family(fam_id)
                gene_obj.is_fragment = True if is_frag == "F" else False
                fam.add_gene(gene_obj)
            if is_frag == "F":
                frag = True
        except Exception:
            raise Exception(f"line {line_counter} of the file '{families_tsv_file.name}' raised an error.")
    bar.close()
    families_tsv_file.close()
    if nb_gene_with_fam < len(pangenome.genes):  # not all genes have an associated cluster
        if nb_gene_with_fam == 0:
            raise Exception("No gene ID in the cluster file matched any gene ID from the annotation step."
                            " Please ensure that the annotations that you loaded previously and the clustering results "
                            "that you have used the same gene IDs. If you use .gff files it is the identifier stored in"
                            " the field 'ID'. If you use .gbff files it is the identifier stored in 'locus_tag'.")
        else:
            if infer_singleton:
                infer_singletons(pangenome)
            else:
                raise Exception(f"Some genes ({len(pangenome.genes) - nb_gene_with_fam}) did not have an associated "
                                f"cluster. Either change your cluster file so that each gene has a cluster, "
                                f"or use the --infer_singletons option to infer a cluster for each non-clustered gene.")
    pangenome.status["genesClustered"] = "Computed"
    if frag:  # if there was fragment information in the file.
        pangenome.status["defragmented"] = "Computed"
    pangenome.parameters["cluster"] = {}
    pangenome.parameters["cluster"]["read_clustering_from_file"] = True
    pangenome.parameters["cluster"]["infer_singletons"] = infer_singleton


def launch(args):
    """ launch the clustering step"""
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    if args.clusters is None:
        clustering(pangenome, args.tmpdir, args.cpu, defrag=not args.no_defrag, code=args.translation_table,
                   coverage=args.coverage, identity=args.identity, mode=args.mode, force=args.force,
                   disable_bar=args.disable_prog_bar)
        logging.getLogger().info("Done with the clustering")
    else:
        read_clustering(pangenome, args.clusters, args.infer_singletons, args.force, disable_bar=args.disable_prog_bar)
        logging.getLogger().info("Done reading the cluster file")
    write_pangenome(pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar)


def subparser(sub_parser):
    parser = sub_parser.add_parser("cluster", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_clust(parser)
    return parser


def parser_clust(parser):
    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--defrag", required=False, action="store_true",
                          help=argparse.SUPPRESS)  # This ensures compatibility with the old option "defrag"
    optional.add_argument('--no_defrag', required=False, default=False, action="store_true",
                          help="DO NOT Use the defragmentation strategy to link potential fragments "
                               "with their original gene family.")
    optional.add_argument("--translation_table", required=False, default="11",
                          help="Translation table (genetic code) to use.")
    optional.add_argument('--clusters', required=False, type=str,
                          help="A tab-separated list containing the result of a clustering. One line per gene. "
                               "First column is cluster ID, and second is gene ID")
    optional.add_argument("--infer_singletons", required=False, action="store_true",
                          help="When reading a clustering result with --clusters, if a gene is not in the provided file"
                               " it will be placed in a cluster where the gene is the only member.")
    optional.add_argument("--mode", required=False, default="1", choices=["0", "1", "2", "3"],
                          help="the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component),"
                               " 2: CD-HIT-like, 3: CD-HIT-like (lowmem)")
    optional.add_argument("--coverage", required=False, type=restricted_float, default=0.8,
                          help="Minimal coverage of the alignment for two proteins to be in the same cluster")
    optional.add_argument("--identity", required=False, type=restricted_float, default=0.8,
                          help="Minimal identity percent for two proteins to be in the same cluster")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_clust(main_parser)
    common = main_parser.add_argument_group(title="Common argument")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--tmpdir", required=False, type=str, default=tempfile.gettempdir(),
                        help="directory for storing temporary files")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
    common.add_argument('-f', '--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
