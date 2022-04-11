#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import tempfile
import time
import logging

# installed libraries
from tqdm import tqdm
import networkx as nx
import pandas as pd

# local libraries
from ppanggolin.formats import checkPangenomeInfo
from ppanggolin.utils import mkOutdir, restricted_float, add_gene, connected_components
from ppanggolin.pangenome import Pangenome
from ppanggolin.align.alignOnPang import get_seq2pang, projectPartition
from ppanggolin.geneFamily import GeneFamily


class GeneContext:
    """
        A class used to represent a gene context

        Attributes
        ----------
        gc_id : int
            ID of the Gene context
        families : set
            Gene families related to the GeneContext

        Methods
        -------
        """

    def __init__(self, gc_id, families=None):
        """ Initial methods

        :param gc_id: ID of the GeneContext
        :type gc_id: int
        :param families: Gene families related to the GeneContext
        :type families: set
        """
        self.ID = gc_id
        self.families = set()
        if families is not None:
            if not all(isinstance(fam, GeneFamily) for fam in families):
                raise Exception(f"You provided elements that were not GeneFamily object."
                                f" GeneContext are only made of GeneFamily")
            self.families |= set(families)

    def add_family(self, family):
        """
        Allow to add one family in the GeneContext
        :param family: family to add
        :type family: GeneFamily
        """
        if not isinstance(family, GeneFamily):
            raise Exception("You did not provide a GenFamily object. Modules are only made of GeneFamily")
        self.families.add(family)


def search_gene_context_in_pangenome(pangenome, output, tmpdir, sequences=None, families=None, transitive=4,
                                     identity=0.5, coverage=0.8, jaccard=0.85, no_defrag=False, cpu=1,
                                     disable_bar=True):
    """
    Main function to search common gene contexts between sequence set and pangenome families

    :param pangenome: Pangenome containing GeneFamilies to align with sequence set
    :type pangenome: Pangenome
    :param sequences: Path to file containing the sequences
    :type sequences: str
    :param families: Path to file containing families name
    :type families: str
    :param output: Path to output directory
    :type output: str
    :param tmpdir: Path to temporary directory
    :type tmpdir: str
    :param transitive: number of genes to check on both sides of a family aligned with an input sequence
    :type transitive: int
    :param identity: minimum identity threshold between sequences and gene families for the alignment
    :type identity: float
    :param coverage: minimum coverage threshold between sequences and gene families for the alignment
    :type coverage: float
    :param jaccard: Jaccard index to filter edges in graph
    :type jaccard: float
    :param no_defrag: do not use the defrag workflow if true
    :type no_defrag: Boolean
    :param cpu: Number of core used to process
    :type cpu: int
    :param disable_bar: Allow preventing bar progress print
    :type disable_bar: Boolean
    """

    # check statuses and load info
    if sequences is not None and pangenome.status["geneFamilySequences"] not in ["inFile", "Loaded", "Computed"]:
        raise Exception("Cannot use this function as your pangenome does not have gene families representatives "
                        "associated to it. For now this works only if the clustering is realised by PPanGGOLiN.")

    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, disable_bar=disable_bar)
    gene_families = {}
    fam_2_seq = None
    if sequences is not None:
        # Alignment of sequences on pangenome families
        new_tmpdir = tempfile.TemporaryDirectory(dir=tmpdir)
        seq_set, _, seq2pan = get_seq2pang(pangenome, sequences, output, new_tmpdir, cpu, no_defrag,
                                           identity, coverage)
        projectPartition(seq2pan, seq_set, output)
        new_tmpdir.cleanup()
        for k, v in seq2pan.items():
            gene_families[v.name] = v
        fam_2_seq = fam2seq(seq2pan)

    if families is not None:
        with open(families, 'r') as f:
            for fam_name in f.read().splitlines():
                gene_families[fam_name] = pangenome.getGeneFamily(fam_name)

    # Compute the graph with transitive closure size provided as parameter
    start_time = time.time()
    logging.getLogger().info("Building the graph...")
    g = compute_gene_context_graph(families=gene_families, t=transitive, disable_bar=disable_bar)
    logging.getLogger().info(
        f"Took {round(time.time() - start_time, 2)} seconds to build the graph to find common gene contexts")
    logging.getLogger().debug(f"There are {nx.number_of_nodes(g)} nodes and {nx.number_of_edges(g)} edges")

    # extract the modules from the graph
    common_components = compute_gene_context(g, jaccard)

    families = set()
    for gene_context in common_components:
        families |= gene_context.families

    if len(families) != 0:
        export_to_dataframe(families, common_components, fam_2_seq, output)
    else:
        logging.getLogger().info(f"No gene contexts were found")

    logging.getLogger().info(f"Computing gene contexts took {round(time.time() - start_time, 2)} seconds")


def compute_gene_context_graph(families, t, disable_bar=False):
    """
    Construct the graph of gene contexts between families of the pangenome

    :param families: Gene families of interest
    :type families: dict
    :param t: transitive value
    :type t: int
    :param disable_bar: Prevents progress bar printing
    :type disable_bar: Boolean

    :return: Graph of gene contexts between interesting gene families of the pangenome
    :rtype: nx.Graph
    """
    g = nx.Graph()
    for family in tqdm(families.values(), unit="families", disable=disable_bar):
        for gene in family.genes:
            contig = gene.contig.genes
            pos_left, in_context_left, pos_right, in_context_right = extract_gene_context(gene, contig, families, t)
            if in_context_left or in_context_right:
                for env_gene in contig[pos_left:pos_right + 1]:
                    _compute_gene_context_graph(g, env_gene, contig, pos_right)
    return g


def _compute_gene_context_graph(g, env_gene, contig, pos_r):
    """
    Compute graph of gene contexts between one gene and the other part of the contig

    :param: Graph of gene contexts between interesting gene families of the pangenome
    :type: nx.Graph
    :param env_gene: Gene of the current position
    :type env_gene: Gene
    :param contig: Current contig to search a gene context
    :type contig: Contig
    :param pos_r: Gene to search a gene context
    :type pos_r: Gene
    """
    g.add_node(env_gene.family)
    add_gene(g.nodes[env_gene.family], env_gene, fam_split=False)
    pos = env_gene.position + 1
    while pos <= pos_r:
        if env_gene.family != contig[pos].family:
            g.add_edge(env_gene.family, contig[pos].family)
            edge = g[env_gene.family][contig[pos].family]
            add_gene(edge, env_gene)
            add_gene(edge, contig[pos])
        pos += 1


def extract_gene_context(gene, contig, families, t=4):
    """
    Extract gene context and whether said gene context exists

    :param gene: Gene of interest
    :type gene: Gene
    :param contig: Gene's contig
    :type contig: Contig
    :param families: Alignment results
    :param families: dict
    :param t: transitive value
    :type t: int

    :return: Position of the context and if it exists for each side ('left' and 'right')
    :rtype: (int, Bool, int, Bool)
    """
    pos_left, pos_right = (max(0, gene.position - t),
                           min(gene.position + t, len(contig) - 1))  # Gene positions to compare family
    in_context_left, in_context_right = (False, False)
    while pos_left < gene.position and not in_context_left:
        if contig[pos_left].family in families.values():
            in_context_left = True
        else:
            pos_left += 1

    while pos_right > gene.position and not in_context_right:
        if contig[pos_right].family in families.values():
            in_context_right = True
        else:
            pos_right -= 1

    return pos_left, in_context_left, pos_right, in_context_right


def compute_gene_context(g, jaccard=0.85):
    """
    Compute the gene contexts in the graph

    :param g: Graph of gene contexts between interesting gene families of the pangenome
    :type g: nx.Graph
    :param jaccard: Jaccard index
    :type jaccard: float

    :return: Set of gene contexts find in graph
    :rtype: Set
    """
    gene_contexts = set()
    c = 1
    for comp in connected_components(g, removed=set(), weight=jaccard):
        gene_contexts.add(GeneContext(gc_id=c, families=comp))
        c += 1
    return gene_contexts


def fam2seq(seq2pan):
    """
    Create a dictionary with gene families as keys and list of sequences id as values

    :param seq2pan: Dictionary storing the sequence ids as keys and the gene families
                    to which they are assigned as values
    :param seq2pan: dict

    :return: Dictionary reversed
    :rtype: dict
    """
    fam_2_seq = {}
    for sequence, family in seq2pan.items():
        if family.ID in fam_2_seq.keys():
            fam_2_seq[family.ID].append(sequence)
        else:
            fam_2_seq[family.ID] = [sequence]
    return fam_2_seq


def export_to_dataframe(families, gene_contexts, fam_2_seq, output):
    """ Export the results into dataFrame

    :param families: Families related to the connected components
    :type families: set
    :param gene_contexts: connected components found in the pangenome
    :type gene_contexts: set
    :param fam_2_seq: Dictionary with gene families as keys and list of sequence ids as values
    :type fam_2_seq: dict
    :param output: output path
    :type output: str
    """
    logging.getLogger().debug(f"There are {len(families)} families among {len(gene_contexts)} gene contexts")

    lines = []
    for gene_context in gene_contexts:
        for family in gene_context.families:
            line = [gene_context.ID]
            if fam_2_seq is None or fam_2_seq.get(family.ID) is None:
                line += [family.name, None, len(family.organisms), family.namedPartition]
            else:
                line += [family.name, ','.join(fam_2_seq.get(family.ID)),
                         len(family.organisms), family.namedPartition]
            lines.append(line)
    df = pd.DataFrame(lines,
                      columns=["GeneContext ID", "Gene family name", "Sequence ID", "Nb Genomes", "Partition"]
                      ).set_index("GeneContext ID")
    df.sort_values(["GeneContext ID", "Sequence ID"], na_position='last').to_csv(
        path_or_buf=f"{output}/gene_contexts.tsv", sep="\t", na_rep='NA')
    logging.getLogger(f"detected gene context(s) are listed in: '{output}/gene_contexts.tsv'")


def launch(args):
    if not any([args.sequences, args.family]):
        raise Exception("At least one of --sequences or --family must be given")
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    search_gene_context_in_pangenome(pangenome=pangenome, output=args.output, tmpdir=args.tmpdir,
                                     sequences=args.sequences, families=args.family, transitive=args.transitive,
                                     identity=args.identity, coverage=args.coverage, jaccard=args.jaccard,
                                     no_defrag=args.no_defrag, cpu=args.cpu, disable_bar=args.disable_prog_bar)


def subparser(sub_parser):
    """
    Parser arguments specific to context command

    :param sub_parser : sub_parser for align command
    :type sub_parser : argparse._SubParsersAction

    :return : parser arguments for align command
    :rtype : argparse.ArgumentParser
    """
    parser = sub_parser.add_parser("context", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    return parser_context(parser)


def parser_context(parser):
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o', '--output', required=True, type=str,
                          help="Output directory where the file(s) will be written")
    onereq = parser.add_argument_group(title="Input file", description="One of the following argument is required :")
    onereq.add_argument('-S', '--sequences', required=False, type=str,
                        help="Fasta file with the sequences of interest")
    onereq.add_argument('-F', '--family', required=False, type=str,
                        help="List of family IDs of interest from the pangenome")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('--no_defrag', required=False, action="store_true",
                          help="DO NOT Realign gene families to link fragments with"
                               "their non-fragmented gene family.")
    optional.add_argument('--identity', required=False, type=float, default=0.5,
                          help="min identity percentage threshold")
    optional.add_argument('--coverage', required=False, type=float, default=0.8,
                          help="min coverage percentage threshold")
    optional.add_argument("-t", "--transitive", required=False, type=int, default=4,
                          help="Size of the transitive closure used to build the graph. This indicates the number of "
                               "non related genes allowed in-between two related genes. Increasing it will improve "
                               "precision but lower sensitivity a little.")
    optional.add_argument("-s", "--jaccard", required=False, type=restricted_float, default=0.85,
                          help="minimum jaccard similarity used to filter edges between gene families. Increasing it "
                               "will improve precision but lower sensitivity a lot.")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import check_log

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_context(main_parser)
    common = main_parser.add_argument_group(title="Common argument")
    common.add_argument("--tmpdir", required=False, type=str, default=tempfile.gettempdir(),
                        help="directory for storing temporary files")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
    common.add_argument('-f', '--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    launch(main_parser.parse_args())
