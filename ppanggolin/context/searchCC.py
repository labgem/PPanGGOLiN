#!/usr/bin/env python3
# coding:utf-8


# default libraries
import argparse
import tempfile
import time
import logging
from importlib import import_module  # Allow to load module when it's needed

# installed libraries
from tqdm import tqdm
import networkx as nx
import pandas as pd

# local libraries
from ppanggolin.formats import checkPangenomeInfo
from ppanggolin.utils import mkOutdir, restricted_float, add_gene, connected_components
from ppanggolin.pangenome import Pangenome
from ppanggolin.align.alignOnPang import get_prot2pang, projectPartition
from ppanggolin.geneFamily import GeneFamily


def _write_graph(g, out_name='graph'):
    """
    Allow to export the graph in png
    :param g: graph to export
    :type g: Graph
    :param out_name: name of the png file
    :type out_name: str
    """
    plt = import_module("matplotlib.pyplot")  # import only if the function is used

    label_nodes = {node: len(node.genes) for node in g.nodes}
    label_edges = {edge: f'({str(len(edge[0].genes))}, {str(len(edge[1].genes))})' for edge in g.edges}
    plt.figure(figsize=(24, 24))
    p = nx.spring_layout(g)
    nx.draw(g, p, with_labels=True, labels=label_nodes)
    nx.draw_networkx_edge_labels(g, p, edge_labels=label_edges, font_color='red')
    plt.savefig(out_name)


class CC:
    """
        A class used to represent a common component

        Attributes
        ----------
        ID : int
            ID of the common component
        families : set
            Gene's families related to the common component

        Methods
        -------
        add_family(family)
            Prints the animals name and what sound it makes
        """
    def __init__(self, ID, families=None):
        """ Initial methods

        :param ID: ID of the common component
        :type ID: int
        :param families: Gene's families related to the common component
        :type families: set
        """
        self.ID = ID
        self.families = set()
        if families is not None:
            if not all(isinstance(fam, GeneFamily) for fam in families):
                raise Exception(f"You provided elements that were not GeneFamily object."
                                f" CC are only made of GeneFamily")
            self.families |= set(families)

    def add_family(self, family):
        """
        Allow to add one family in the common component
        :param family: family ti add
        :type family: GeneFamily
        """
        if not isinstance(family, GeneFamily):
            raise Exception("You did not provide a GenFamily object. Modules are only made of GeneFamily")
        self.families.add(family)


def search_cc_in_pangenome(pangenome, proteins, output, tmpdir, transitive=4, identity=0.5, coverage=0.8, jaccard=0.85,
                           no_defrag=False, cpu=1, disable_bar=True, write_graph='graph'):
    """
    Main function to search common component between protein set and pangenome's families

    :param pangenome: Pangenome containing GeneFamilies to align with protein set
    :type pangenome: Pangenome
    :param proteins: Path to file containing the proteins
    :type proteins: str
    :param output: Path to output directory
    :type output: str
    :param tmpdir: Path to temporary directory
    :type tmpdir: str
    :param transitive: number of considering gene's families at each side of family align with protein
    :type transitive: int
    :param identity: minimum rate identity between proteins and gene's families
    :type identity: float
    :param coverage: minimum rate covery between proteins and gene's families
    :type coverage: float
    :param jaccard: Jaccard index to filter edges in graph
    :type jaccard: float
    :param no_defrag: Allow preventing defragmentation in alignment
    :type no_defrag: Boolean
    :param cpu: Number of core used to process
    :type cpu: int
    :param disable_bar: Allow preventing bar progress print
    :type disable_bar: Boolean
    :param write_graph: Allow writing graph of common component
    :type write_graph: str
    """
    # check statuses and load info
    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, disable_bar=disable_bar)

    # Alignment of proteins on pangenome's families
    new_tmpdir = tempfile.TemporaryDirectory(dir=tmpdir)
    prot_set, _, prot2pan = get_prot2pang(pangenome, proteins, output, new_tmpdir, cpu, no_defrag, identity, coverage)
    projectPartition(prot2pan, prot_set, output)
    new_tmpdir.cleanup()

    # Compute the graph with transitive closure size provided as parameter
    start_time = time.time()
    logging.getLogger().info("Building the graph...")
    g = compute_cc_graph(alignment=prot2pan, t=transitive, disable_bar=disable_bar)
    logging.getLogger().info(
        f"Took {round(time.time() - start_time, 2)} seconds to build the graph to find common component in")
    logging.getLogger().info(f"There are {nx.number_of_nodes(g)} nodes and {nx.number_of_edges(g)} edges")

    # extract the modules from the graph
    common_components = compute_cc(g, jaccard)

    if write_graph is not None:
        _write_graph(g, out_name=f"{output}/{write_graph}.png")
        logging.getLogger().info(f"Write graph took {round(time.time() - start_time, 2)} seconds")

    families = set()
    for cc in common_components:
        families |= cc.families

    if len(families) != 0:
        logging.getLogger().info(f"There are {len(families)} families among {len(common_components)} common components")

        fam_2_prot = fam2prot(prot2pan)

        lines = []
        for cc in common_components:
            for family in cc.families:
                line = [cc.ID]
                if fam_2_prot.get(family.ID) is None:
                    line += [family.name, None, len(family._genePerOrg.keys())]
                else:
                    line += [family.name, ','.join(fam_2_prot.get(family.ID)), len(family._genePerOrg.keys())]
                lines.append(line)
        df = pd.DataFrame(lines, columns=["CC ID", "Gene family", "Protein ID", "Nb Genomes"]).set_index(
            "CC ID").sort_index()
        df.to_csv(path_or_buf=f"{output}/Common_component.tsv", sep="\t", na_rep='NA')
    else:
        logging.getLogger().info(f"No common component were find")

    logging.getLogger().info(f"Computing common components took {round(time.time() - start_time, 2)} seconds")


def compute_cc_graph(alignment, t, disable_bar=False):
    """
    Construct the graph of common component between proteins and pangenome gene's families

    :param alignment: Result of alignment between proteins and pangenome gene's families
    :type alignment: dict
    :param t: transitive value
    :type t: int
    :param disable_bar: Allow preventing bar progress print
    :type disable_bar: Boolean

    :return: Graph of common component between proteins and pangenome gene's families
    :rtype: nx.Graph
    """
    g = nx.Graph()
    for protein, gene_family in tqdm(alignment.items(), unit="proteins", disable=disable_bar):
        for gene in gene_family.genes:
            contig = gene.contig._genes_position  # TODO create method to extract
            pos_left, in_context_left, pos_right, in_context_right = extract_gene_context(gene, contig, alignment, t)
            if in_context_left or in_context_right:
                for env_gene in contig[pos_left:pos_right + 1]:
                    _compute_cc_graph(g, env_gene, contig, pos_right)
    return g


def _compute_cc_graph(g, env_gene, contig, pos_r):
    """
    Compute graph of common component between one gene and the other part of the contig

    :param: Graph of common component between proteins and pangenome gene's families
    :type: nx.Graph
    :param env_gene: Gene of the current position
    :type env_gene: Gene
    :param contig: Current contig to search CC
    :type contig: Contig
    :param pos_r: Gene to search a common component
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


def extract_gene_context(gene, contig, alignment, t=4):
    """
    Extract gene context if exist

    :param gene: Gene of interest
    :type gene: Gene
    :param contig: Gene's contig
    :type contig: Contig
    :param alignment: Alignment results
    :param alignment: dict
    :param t: transitive value
    :type t: int

    :return: Position of the context and if exist
    :rtype: (int, Bool, int, Bool)
    """
    pos_left, pos_right = (max(0, gene.position - t),
                           min(gene.position + t, len(contig)-1))  # Gene positions to compare family
    in_context_left, in_context_right = (False, False)
    while pos_left < gene.position and not in_context_left:
        if contig[pos_left].family in alignment.values():
            in_context_left = True
        else:
            pos_left += 1

    while pos_right > gene.position and not in_context_right:
        if contig[pos_right].family in alignment.values():
            in_context_right = True
        else:
            pos_right -= 1

    return pos_left, in_context_left, pos_right, in_context_right


def compute_cc(g, jaccard=0.85):
    """
    Compute the common component in the graph

    :param: Graph of common component between proteins and pangenome gene's families
    :type: nx.Graph
    :param jaccard: Jaccard index
    :type jaccard: float

    :return: Set of common component find in graph
    :rtype: Set
    """
    cc = set()
    c = 1
    for comp in connected_components(g, removed=set(), weight=jaccard):
        cc.add(CC(ID=c, families=comp))
        c += 1
    return cc


def fam2prot(prot2pan):
    """
    Create a dictionary with gene's families in keys and list of proteins in values

    :param prot2pan: Dictionary to reverse
    :param prot2pan: dict

    :return: Dictionary reversed
    :rtype: dict
    """
    fam_2_prot = {}
    for protein, family in prot2pan.items():
        if family.ID in fam_2_prot.keys():
            fam_2_prot[family.ID].append(protein)
        else:
            fam_2_prot[family.ID] = [protein]
    return fam_2_prot


def launch(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    search_cc_in_pangenome(pangenome=pangenome, proteins=args.proteins, output=args.output, identity=args.identity,
                           coverage=args.coverage, jaccard=args.jaccard, transitive=args.transitive, tmpdir=args.tmpdir,
                           no_defrag=args.no_defrag, cpu=args.cpu, disable_bar=args.disable_prog_bar,
                           write_graph=args.write_graph)


def contextSubparser(sub_parser):
    """
    Parser arguments specific to context command

    :param sub_parser : sub_parser for align command
    :type sub_parser : argparse._SubParsersAction

    :return : parser arguments for align command
    :rtype : argparse.ArgumentParser
    """
    parser = sub_parser.add_parser("context", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o', '--output', required=True, type=str,
                          help="Output directory where the file(s) will be written")
    required.add_argument('-P', '--proteins', required=True, type=str, help="Fasta with all proteins of interest")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('--no_defrag', required=False, action="store_true",
                          help="DO NOT Realign gene families to link fragments with"
                               "their non-fragmented gene family. (default: False)")
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
    optional.add_argument("--write_graph", type=str, default='graph',
                          help=argparse.SUPPRESS)  # Allow to export the graph in png
    return parser
