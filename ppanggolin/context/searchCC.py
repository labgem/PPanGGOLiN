#!/usr/bin/env python3
# coding:utf-8


# default libraries
import argparse
import tempfile
import networkx as nx

# local libraries
from ppanggolin.formats import checkPangenomeInfo
from ppanggolin.utils import mkOutdir
from ppanggolin.pangenome import Pangenome
from ppanggolin.align.alignOnPang import get_prot2pang
from ppanggolin.mod.module import add_gene, connected_components

import matplotlib.pyplot as plt
# def get_gene_position(pangenome):
#     pos2gene = {}
#     for gene in pangenome.genes:
#         pos2gene[gene.position]=gene
#     return pos2gene


# def add_gene(gene):

#def add_family(conti):

def write_graph(g):
    labelnodes = {node: node.ID for node in g.nodes}
    labeledges = {edge: f'({str(edge[0].ID)}, {str(edge[1].ID)})' for edge in g.edges}
    plt.figure(figsize=(24, 24))
    p = nx.spring_layout(g)
    nx.draw(g, p, with_labels=True, labels=labelnodes)
    nx.draw_networkx_edge_labels(g, p, edge_labels=labeledges, font_color='red')
    plt.savefig("test.png")


def search_cc_in_pangenome(pangenome, proteins, output, tmpdir,  transitive=4, identity=0.8, coverage=0.8, defrag=False, cpu=1):
    new_tmpdir = tempfile.TemporaryDirectory(dir=tmpdir)

    checkPangenomeInfo(pangenome, needFamilies=True, needAnnotations=True)

    prot2pan = get_prot2pang(pangenome, proteins, output, new_tmpdir, cpu, defrag, identity, coverage)[-1]

    cc_set = set()

    g = nx.Graph()
    for protein, gene_family in prot2pan.items():
        # print(protein, gene_family)
        for gene in gene_family.genes:
            # print(gene.contig._genes_position[gene.position-transitive:gene.position+transitive+1])
            pos_left, pos_right = (max(0, gene.position-transitive),
                                   gene.position+transitive)  # Gene position to compare family
            print(f'pos_left : {pos_left}; gene_pos : {gene.position}; pos_right : {pos_right}')
            in_context_left, in_context_right = (False, False)
            contig = gene.contig._genes_position  # TODO create method to extract
            while pos_left < gene.position and not in_context_left:
                if contig[pos_left].family in prot2pan.values():
                    in_context_left = True
                else:
                    pos_left += 1

            while pos_right < gene.position and not in_context_right:
                if contig[pos_right].family in prot2pan.values():
                    in_context_right = True
                else:
                    pos_right -= 1

            if in_context_left or in_context_right:
                print(contig[pos_left:pos_right+1], len(contig[pos_left:pos_right+1]), gene.position, pos_left,
                      pos_right)
                for env_gene in contig[pos_left:pos_right+1]:
                    g.add_node(env_gene.family)
                    add_gene(g.nodes[env_gene.family], gene, fam_split=False)
                    pos = env_gene.position + 1
                    while pos <= pos_right and pos < len(contig):
                        if env_gene.family != contig[pos].family:
                            g.add_edge(env_gene.family, contig[pos].family)
                            edge = g[env_gene.family][contig[pos].family]
                            add_gene(edge, env_gene)
                            add_gene(edge, contig[pos])
                        pos += 1
                cc_set = cc_set.union(set(contig[pos_left:pos_right+1]))
    l = set()
    for comp in connected_components(g, removed=set(), weight=0.85):
        if not any(fam.namedPartition == "persistent" for fam in comp):
            l.add(comp)

    # write_graph(g)
    new_tmpdir.cleanup()


def launch(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    search_cc_in_pangenome(pangenome=pangenome, proteins=args.proteins, output=args.output, tmpdir=args.tmpdir)


def contextSubparser(sub_parser):
    """
    Parser arguments specific to align command

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
    return parser
