#!/usr/bin/env python3

# default libraries
import argparse
import logging
from multiprocessing import get_context
from collections import Counter, defaultdict
import logging
from typing import TextIO, List, Dict, Set, Any
from pathlib import Path
from typing import TextIO
from importlib.metadata import distribution
from statistics import median, mean, stdev
import os
import csv

# installed libraries
import pandas as pd
from tqdm import tqdm

# local libraries
from ppanggolin.edge import Edge
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.genome import Organism
from ppanggolin.region import Region
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import (
    write_compressed_or_not,
    mk_outdir,
    restricted_float,
    flatten_nested_dict,
)
from ppanggolin.formats.readBinaries import check_pangenome_info

# global variable to store the pangenome
pan = Pangenome()  # TODO change to pangenome:Pangenome = Pangenome=() ?
needAnnotations = False
needFamilies = False
needGraph = False
needPartitions = False
needSpots = False
needRegions = False
needModules = False
needMetadata = False
metatype = False
ignore_err = False


def write_json_header(json: TextIO):
    """Write the header of json file to save graph

    :param json: file-like object, compressed or not
    """
    json.write('{"directed": false, "multigraph": false,')
    json.write(' "graph": {')
    json.write(' "genomes": {')
    orgstr = []
    for org in pan.organisms:
        orgstr.append('"' + org.name + '": {')
        contigstr = []
        for contig in org.contigs:
            contigstr.append(
                f'"{contig.name}": '
                + '{"is_circular": '
                + ("true" if contig.is_circular else "false")
                + "}"
            )
        orgstr[-1] += ", ".join(contigstr) + "}"

    json.write(", ".join(orgstr) + "}")
    # if other things are to be written such as the parameters, write them here
    json.write("},")


def write_json_gene_fam(gene_fam: GeneFamily, json: TextIO):
    """Write the gene families corresponding to node graph in json file

    :param gene_fam: file-like object, compressed or not
    :param json: file-like object, compressed or not
    """
    json.write(
        "{" + f'"id": "{gene_fam.name}", "nb_genes": {len(gene_fam)}, '
        f'"partition": "{gene_fam.named_partition}", "subpartition": "{gene_fam.partition}"'
    )
    org_dict = {}
    name_counts = Counter()
    product_counts = Counter()
    length_counts = Counter()
    for gene in gene_fam.genes:
        name_counts[gene.name] += 1
        product_counts[gene.product] += 1
        length_counts[gene.stop - gene.start] += 1
        try:
            org_dict[gene.organism][gene.contig].append(gene)
        except KeyError:
            try:
                org_dict[gene.organism][gene.contig] = [gene]
            except KeyError:
                org_dict[gene.organism] = {gene.contig: [gene]}

    json.write(
        f', "name": "{name_counts.most_common(1)[0][0]}", "product": "{product_counts.most_common(1)[0][0]}", '
        f'"length": {length_counts.most_common(1)[0][0]}'
    )

    json.write(', "genomes": {')
    orgstr = []
    for org in org_dict:
        orgstr.append('"' + org.name + '": {')
        contigstr = []
        for contig in org_dict[org]:
            contigstr.append('"' + contig.name + '": {')
            genestr = []
            for gene in org_dict[org][contig]:
                identifier = (
                    gene.ID if gene.local_identifier == "" else gene.local_identifier
                )
                genestr.append(
                    '"'
                    + identifier
                    + '": {'
                    + f'"name": "{gene.name}", "product": "{gene.product}", '
                    f'"is_fragment": {"true" if gene.is_fragment else "false"},'
                    f' "position": {gene.position}, "strand": "{gene.strand}",'
                    f' "end": {gene.stop}, "start": {gene.start}' + "}"
                )
            contigstr[-1] += ", ".join(genestr) + "}"
        orgstr[-1] += ", ".join(contigstr) + "}"
    json.write(", ".join(orgstr) + "}}")


def write_json_nodes(json: TextIO):
    """Write the node graph in json file

    :param json: file-like object, compressed or not
    """
    json.write('"nodes": [')
    fam_list = list(pan.gene_families)
    first_fam = fam_list[0]
    write_json_gene_fam(first_fam, json)
    for family in fam_list[1:]:
        json.write(", ")
        write_json_gene_fam(family, json)
    json.write("]")


def write_json_edge(edge: Edge, json: TextIO):
    """Write the edge graph in json file

    :param edge: file-like object, compressed or not
    :param json: file-like object, compressed or not
    """
    json.write("{")
    json.write(
        f'"weight": {len(edge.gene_pairs)}, "source": "{edge.source.name}", "target": "{edge.target.name}"'
    )
    json.write(', "genomes": {')
    orgstr = []
    for org in edge.organisms:
        orgstr.append('"' + org.name + '": [')
        genepairstr = []
        for gene_pair in edge.get_organism_genes_pairs(org):
            genepairstr.append(
                '{"source": "'
                + gene_pair[0].ID
                + '", "target": "'
                + gene_pair[1].ID
                + f'", "length": {gene_pair[0].start - gene_pair[1].stop}'
                + "}"
            )
        orgstr[-1] += ", ".join(genepairstr) + "]"
    json.write(", ".join(orgstr) + "}}")


def write_json_edges(json):
    """Write the edge graph in json file

    :param json: file-like object, compressed or not
    """
    json.write(', "links": [')
    edge_list = list(pan.edges)
    write_json_edge(edge_list[0], json)
    for edge in edge_list[1:]:
        json.write(", ")
        write_json_edge(edge, json)
    json.write("]")


def write_json(output: Path, compress: bool = False):
    """Writes the graph in a json file format

    :param output: Path to output directory
    :param compress: Compress the file in .gz
    """
    logging.getLogger("PPanGGOLiN").info(
        "Writing the json file for the pangenome graph..."
    )
    outname = output / "pangenomeGraph.json"
    with write_compressed_or_not(outname, compress) as json:
        write_json_header(json)
        write_json_nodes(json)
        write_json_edges(json)
        json.write("}")
    logging.getLogger("PPanGGOLiN").info(
        f"Done writing the json file : '{outname.as_posix()}'"
    )


def write_gexf_header(gexf: TextIO, light: bool = True):
    """Write the header of gexf file to save graph

    :param gexf: file-like object, compressed or not
    :param light: save the light version of the pangenome graph
    """
    index = None
    if not light:
        index = pan.get_org_index()  # has been computed already
    gexf.write(
        '<?xml version="1.1" encoding="UTF-8"?>\n<gexf xmlns:viz="https://www.gexf.net/1.2draft/viz"'
        ' xmlns="https://www.gexf.net/1.2draft" version="1.2">\n'
    )  # TODO update link
    gexf.write('  <graph mode="static" defaultedgetype="undirected">\n')
    gexf.write('    <attributes class="node" mode="static">\n')
    gexf.write('      <attribute id="0" title="nb_genes" type="long" />\n')
    gexf.write('      <attribute id="1" title="name" type="string" />\n')
    gexf.write('      <attribute id="2" title="product" type="string" />\n')
    gexf.write('      <attribute id="3" title="type" type="string" />\n')
    gexf.write('      <attribute id="4" title="partition" type="string" />\n')
    gexf.write('      <attribute id="5" title="subpartition" type="string" />\n')
    gexf.write('      <attribute id="6" title="partition_exact" type="string" />\n')
    gexf.write('      <attribute id="7" title="partition_soft" type="string" />\n')
    gexf.write('      <attribute id="8" title="length_avg" type="double" />\n')
    gexf.write('      <attribute id="9" title="length_med" type="long" />\n')
    gexf.write('      <attribute id="10" title="nb_genomes" type="long" />\n')

    if pan.number_of_spots > 0:
        gexf.write('      <attribute id="12" title="spot" type="string" />\n')
    if pan.number_of_modules > 0:
        gexf.write('      <attribute id="13" title="module" type="string" />\n')
    shift = 14

    source_fields = {
        m.source: m.fields
        for f in pan.gene_families
        if len(list(f.metadata)) > 0
        for m in f.metadata
    }
    for source_metadata_families in pan.metadata_sources("families"):
        for field in source_fields[source_metadata_families]:
            gexf.write(
                f'      <attribute id="{shift}" title="{source_metadata_families}_{field}" type="string" />\n'
            )
            shift += 1
    if not light:
        for org, org_idx in index.items():
            gexf.write(
                f'      <attribute id="{org_idx + shift}" title="{org.name}" type="string" />\n'
            )
    gexf.write("    </attributes>\n")
    gexf.write('    <attributes class="edge" mode="static">\n')
    gexf.write('      <attribute id="11" title="nb_genes" type="long" />\n')
    if not light:
        for org, org_idx in index.items():
            gexf.write(
                f'      <attribute id="{org_idx + len(index) + shift}" title="{org.name}" type="long" />\n'
            )
    gexf.write("    </attributes>\n")
    gexf.write("    <meta>\n")
    gexf.write(
        f'      <creator>PPanGGOLiN {distribution("ppanggolin").version}</creator>\n'
    )
    gexf.write("    </meta>\n")


def write_gexf_nodes(gexf: TextIO, light: bool = True, soft_core: False = 0.95):
    """Write the node of pangenome graph in gexf file

    :param gexf: file-like object, compressed or not
    :param light: save the light version of the pangenome graph
    :param soft_core: Soft core threshold to use
    """
    index = None
    gexf.write("    <nodes>\n")
    colors = {
        "persistent": 'a="0" b="7" g="165" r="247"',
        "shell": 'a="0" b="96" g="216" r="0"',
        "cloud": 'a="0" b="255" g="222" r="121"',
    }
    if not light:
        index = pan.get_org_index()

    pan_metadata_sources = pan.metadata_sources("families")

    for fam in pan.gene_families:
        name = Counter()
        product = Counter()
        gtype = Counter()
        lis = []
        for gene in fam.genes:
            name[gene.name] += 1
            product[gene.product.replace("&", "and")] += 1
            gtype[gene.type] += 1
            lis.append(gene.stop - gene.start)

        gexf.write(f'      <node id="{fam.ID}" label="{fam.name}">\n')
        gexf.write(f"        <viz:color {colors[fam.named_partition]} />\n")
        gexf.write(f'        <viz:size value="{fam.number_of_organisms}" />\n')
        gexf.write("        <attvalues>\n")
        gexf.write(f'          <attvalue for="0" value="{len(fam)}" />\n')
        gexf.write(
            f'          <attvalue for="1" value="{name.most_common(1)[0][0]}" />\n'
        )
        gexf.write(
            f'          <attvalue for="2" value="{product.most_common(1)[0][0]}" />\n'
        )
        gexf.write(
            f'          <attvalue for="3" value="{gtype.most_common(1)[0][0]}" />\n'
        )
        gexf.write(f'          <attvalue for="4" value="{fam.named_partition}" />\n')
        gexf.write(f'          <attvalue for="5" value="{fam.partition}" />\n')
        gexf.write(
            f'          <attvalue for="6" value="'
            f'{"exact_accessory" if fam.number_of_organisms != pan.number_of_organisms else "exact_core"}" />\n'
        )
        gexf.write(
            f'          <attvalue for="7" value="'
            f'{"soft_core" if fam.number_of_organisms >= (pan.number_of_organisms * soft_core) else "soft_accessory"}"'
            f" />\n"
        )
        gexf.write(
            f'          <attvalue for="8" value="{round(sum(lis) / len(lis), 2)}" />\n'
        )
        gexf.write(f'          <attvalue for="9" value="{int(median(lis))}" />\n')
        gexf.write(
            f'          <attvalue for="10" value="{fam.number_of_organisms}" />\n'
        )
        if pan.number_of_spots > 0:
            str_spot = "|".join([str(s) for s in list(fam.spots)])
            gexf.write(f'          <attvalue for="12" value="{str_spot}"/>\n')
        if pan.number_of_modules > 0:
            str_module = str(fam.module) if fam.has_module else ""
            gexf.write(f'          <attvalue for="13" value="{str_module}"/>\n')
        shift = 14
        source_fields = {
            m.source: m.fields
            for f in pan.gene_families
            if len(list(f.metadata)) > 0
            for m in f.metadata
        }
        for source_metadata_families in pan_metadata_sources:
            to_concat = defaultdict(list)
            for m in fam.metadata:
                if m.source == source_metadata_families:
                    for field in m.fields:
                        to_concat[field].append(str(getattr(m, field)))
            for field in source_fields[source_metadata_families]:
                concatenated_fields = "|".join(to_concat[field])
                gexf.write(
                    f'      <attvalue for="{shift}" value="{concatenated_fields}"/>\n'
                )
                shift += 1
        if not light:
            for org, genes in fam.get_org_dict().items():
                gexf.write(
                    f'          <attvalue for="'
                    f'{index[org] + shift}" '
                    f'value="{"|".join([gene.ID if gene.local_identifier == "" else gene.local_identifier for gene in genes])}" />\n'
                )
        gexf.write("        </attvalues>\n")
        gexf.write("      </node>\n")
    gexf.write("    </nodes>\n")


def write_gexf_edges(gexf: TextIO, light: bool = True):
    """Write the edge of pangenome graph in gexf file

    :param gexf: file-like object, compressed or not
    :param light: save the light version of the pangenome graph
    """
    gexf.write("    <edges>\n")
    edgeids = 0
    index = pan.get_org_index()
    shift = 14
    metadata_count = len(pan.metadata_sources("families"))
    for edge in pan.edges:
        gexf.write(
            f'      <edge id="{edgeids}" source="'
            f'{edge.source.ID}" target="{edge.target.ID}" weight="{edge.number_of_organisms}">\n'
        )
        gexf.write(f'        <viz:thickness value="{edge.number_of_organisms}" />\n')
        gexf.write("        <attvalues>\n")
        gexf.write(f'          <attvalue for="11" value="{len(edge.gene_pairs)}" />\n')
        if not light:
            for org, genes_pairs in edge.get_organisms_dict().items():
                gexf.write(
                    f'          <attvalue for="{index[org] + len(index) + metadata_count + shift}" value="{len(genes_pairs)}" />\n'
                )
        gexf.write("        </attvalues>\n")
        gexf.write("      </edge>\n")
        edgeids += 1
    gexf.write("    </edges>\n")


def write_gexf_end(gexf: TextIO):
    """Write the end of gexf file to save pangenome

    :param gexf: file-like object, compressed or not
    """
    gexf.write("  </graph>")
    gexf.write("</gexf>")


def write_gexf(output: Path, light: bool = True, compress: bool = False):
    """Write the node of pangenome in gexf file

    :param output: Path to output directory
    :param light: save the light version of the pangenome graph
    :param compress: Compress the file in .gz
    """
    txt = "Writing the "
    txt += (
        "light gexf file for the pangenome graph..."
        if light
        else "gexf file for the pangenome graph..."
    )

    logging.getLogger("PPanGGOLiN").info(txt)
    outname = output / f"pangenomeGraph{'_light' if light else ''}.gexf"
    with write_compressed_or_not(outname, compress) as gexf:
        graph_type = "ligth gexf" if light else "gexf"
        logging.getLogger("PPanGGOLiN").debug(f"Writing the {graph_type} header...")
        write_gexf_header(gexf, light)
        logging.getLogger("PPanGGOLiN").debug(f"Writing the {graph_type} nodes...")
        write_gexf_nodes(gexf, light)
        logging.getLogger("PPanGGOLiN").debug(f"Writing the {graph_type} edges...")
        write_gexf_edges(gexf, light)
        logging.getLogger("PPanGGOLiN").debug(f"Writing the {graph_type} ends...")
        write_gexf_end(gexf)
        logging.getLogger("PPanGGOLiN").info(
            f"Done writing the gexf file : '{gexf.name}'"
        )


def write_matrix(
    output: Path,
    sep: str = ",",
    ext: str = "csv",
    compress: bool = False,
    gene_names: bool = False,
):
    """
    Write a csv file format as used by Roary, among others.
    The alternative gene ID will be the partition, if there is one

    :param sep: Column field separator
    :param ext: file extension
    :param output: Path to output directory
    :param compress: Compress the file in .gz
    :param gene_names: write the genes name if there are saved in  pangenome
    """
    logging.getLogger("PPanGGOLiN").info(f"Writing the .{ext} file ...")
    outname = output / f"matrix.{ext}"
    with write_compressed_or_not(outname, compress) as matrix:

        index_org = {}
        default_dat = []
        for index, org in enumerate(pan.organisms):
            default_dat.append("0")
            index_org[org] = index

        matrix.write(
            sep.join(
                [
                    '"Gene"',  # 1
                    '"Non-unique Gene name"',  # 2
                    '"Annotation"',  # 3
                    '"No. isolates"',  # 4
                    '"No. sequences"',  # 5
                    '"Avg sequences per isolate"',  # 6
                    '"Accessory Fragment"',  # 7
                    '"Genome Fragment"',  # 8
                    '"Order within Fragment"',  # 9
                    '"Accessory Order with Fragment"',  # 10
                    '"QC"',  # 11
                    '"Min group size nuc"',  # 12
                    '"Max group size nuc"',  # 13
                    '"Avg group size nuc"',
                ]  # 14
                + ['"' + str(org) + '"' for org in pan.organisms]
            )
            + "\n"
        )  # 15
        default_genes = (
            ['""'] * pan.number_of_organisms
            if gene_names
            else ["0"] * pan.number_of_organisms
        )
        org_index = pan.get_org_index()  # should just return things
        for fam in pan.gene_families:
            genes = default_genes.copy()
            lis = []
            genenames = Counter()
            product = Counter()
            for org, gene_list in fam.get_org_dict().items():
                genes[org_index[org]] = (
                    " ".join(['"' + str(gene) + '"' for gene in gene_list])
                    if gene_names
                    else str(len(gene_list))
                )
                for gene in gene_list:
                    lis.append(gene.stop - gene.start)
                    product[gene.product] += 1
                    genenames[gene.name] += 1

            if fam.partition != "":
                alt = fam.named_partition
            else:
                alt = str(product.most_common(1)[0][0])

            lis = [gene.stop - gene.start for gene in fam.genes]
            matrix.write(
                sep.join(
                    [
                        '"' + fam.name + '"',  # 1
                        '"' + alt + '"',  # 2
                        '"' + str(product.most_common(1)[0][0]) + '"',  # 3
                        '"' + str(fam.number_of_organisms) + '"',  # 4
                        '"' + str(len(fam)) + '"',  # 5
                        '"'
                        + str(round(len(fam) / fam.number_of_organisms, 2))
                        + '"',  # 6
                        '"NA"',  # 7
                        '"NA"',  # 8
                        '""',  # 9
                        '""',  # 10
                        '""',  # 11
                        '"' + str(min(lis)) + '"',  # 12
                        '"' + str(max(lis)) + '"',  # 13
                        '"' + str(round(sum(lis) / len(lis), 2)) + '"',
                    ]  # 14
                    + genes
                )
                + "\n"
            )  # 15
    logging.getLogger("PPanGGOLiN").info(
        f"Done writing the matrix : '{outname.as_posix()}'"
    )


def write_gene_presence_absence(output: Path, compress: bool = False):
    """
    Write the gene presence absence matrix

    :param output: Path to output directory
    :param compress: Compress the file in .gz
    """
    logging.getLogger("PPanGGOLiN").info("Writing the gene presence absence file ...")
    outname = output / "gene_presence_absence.Rtab"
    with write_compressed_or_not(outname, compress) as matrix:
        index_org = {}
        default_dat = []
        for index, org in enumerate(pan.organisms):
            default_dat.append("0")
            index_org[org] = index

        matrix.write(
            "\t".join(["Gene"] + [str(org) for org in pan.organisms]) + "\n"  # 14
        )  # 15
        default_genes = ["0"] * pan.number_of_organisms
        org_index = pan.get_org_index()  # should just return things
        for fam in pan.gene_families:
            genes = default_genes.copy()
            for org in fam.organisms:
                genes[org_index[org]] = "1"

            matrix.write("\t".join([fam.name] + genes) + "\n")  # 14  # 15
    logging.getLogger("PPanGGOLiN").info(
        f"Done writing the gene presence absence file : '{outname.as_posix()}'"
    )


def summarize_genome(
    organism: Organism,
    pangenome_persistent_count: int,
    pangenome_persistent_single_copy_families: Set[GeneFamily],
    soft_core_families: Set[GeneFamily],
    exact_core_families: Set[GeneFamily],
    rgp_count: int,
    spot_count: int,
    module_count: int,
) -> Dict[str, any]:
    """
    Summarizes genomic information of an organism.

    :param input_organism: The organism for which the genome is being summarized.
    :param pangenome_persistent_count: Count of persistent genes in the pangenome.
    :param pangenome_persistent_single_copy_families: Set of gene families considered as persistent single-copy in the pangenome.
    :param soft_core_families: soft core families of the pangenome
    :param exact_core_families: exact core families of the pangenome
    :param input_org_rgps:  Number of regions of genomic plasticity in the input organism. None if not computed.
    :param input_org_spots:  Number of spots in the input organism. None if not computed.
    :param input_org_modules: Number of modules in the input organism. None if not computed.

    :return: A dictionary containing various summary information about the genome.
    """

    partition_to_genes = organism.group_genes_by_partition()

    persistent_gene_count = len(partition_to_genes["persistent"])
    shell_gene_count = len(partition_to_genes["shell"])
    cloud_gene_count = len(partition_to_genes["cloud"])

    gene_count = persistent_gene_count + shell_gene_count + cloud_gene_count

    persistent_family_count = len({g.family for g in partition_to_genes["persistent"]})
    shell_family_count = len({g.family for g in partition_to_genes["shell"]})
    cloud_family_count = len({g.family for g in partition_to_genes["cloud"]})

    persistent_fragmented_genes = {
        g for g in partition_to_genes["persistent"] if g.is_fragment
    }
    shell_fragmented_genes = {g for g in partition_to_genes["shell"] if g.is_fragment}
    cloud_fragmented_genes = {g for g in partition_to_genes["cloud"] if g.is_fragment}

    persistent_fragmented_genes_count = len(persistent_fragmented_genes)
    shell_fragmented_genes_count = len(shell_fragmented_genes)
    cloud_fragmented_genes_count = len(cloud_fragmented_genes)

    fragmented_genes_count = (
        persistent_fragmented_genes_count
        + shell_fragmented_genes_count
        + cloud_fragmented_genes_count
    )

    persistent_fragmented_family_count = len(
        {g.family for g in persistent_fragmented_genes}
    )
    shell_fragmented_family_count = len({g.family for g in shell_fragmented_genes})
    cloud_fragmented_family_count = len({g.family for g in cloud_fragmented_genes})

    families_with_fragment_count = (
        persistent_fragmented_family_count
        + shell_fragmented_family_count
        + cloud_fragmented_family_count
    )

    families_count = persistent_family_count + shell_family_count + cloud_family_count

    completeness = "NA"
    if pangenome_persistent_count > 0:
        completeness = round(
            (persistent_family_count / pangenome_persistent_count) * 100, 2
        )

    orgs_families_in_multicopy_by_part = defaultdict(set)
    for family in organism.families:
        if len(family.get_org_dict()[organism]) > 1:
            # the family has more than one gene in the genome
            orgs_families_in_multicopy_by_part[family.named_partition].add(family)

    orgs_persistent_families_in_multicopy_count = len(
        orgs_families_in_multicopy_by_part["persistent"]
    )
    orgs_shell_families_in_multicopy_count = len(
        orgs_families_in_multicopy_by_part["shell"]
    )
    orgs_cloud_families_in_multicopy_count = len(
        orgs_families_in_multicopy_by_part["cloud"]
    )

    orgs_families_in_multicopy_count = (
        orgs_persistent_families_in_multicopy_count
        + orgs_shell_families_in_multicopy_count
        + orgs_cloud_families_in_multicopy_count
    )

    single_copy_families_found_in_multicopy_count = len(
        pangenome_persistent_single_copy_families
        & orgs_families_in_multicopy_by_part["persistent"]
    )
    contamination = "NA"
    if len(pangenome_persistent_single_copy_families) > 0:
        contamination = round(
            100
            * single_copy_families_found_in_multicopy_count
            / len(pangenome_persistent_single_copy_families),
            2,
        )

    fragmentation = "NA"
    if families_count > 0:
        fragmentation = round(100.0 * families_with_fragment_count / families_count, 2)

    soft_core_genes = {
        gene for gene in organism.genes if gene.family in soft_core_families
    }
    exact_core_genes = {
        gene for gene in organism.genes if gene.family in exact_core_families
    }

    soft_core_families_count = len({gene.family for gene in soft_core_genes})
    exact_core_families_count = len({gene.family for gene in exact_core_genes})

    rgp_count = "Not computed" if rgp_count is None else rgp_count
    spot_count = "Not computed" if spot_count is None else spot_count
    module_count = "Not computed" if module_count is None else module_count

    summary_info = {
        "Genome_name": organism.name,
        "Contigs": organism.number_of_contigs,
        "Genes": gene_count,
        "Fragmented_genes": fragmented_genes_count,
        "Families": families_count,
        "Families_with_fragments": families_with_fragment_count,
        "Families_in_multicopy": orgs_families_in_multicopy_count,
        "Soft_core": {
            "families": soft_core_families_count,
            "genes": len(soft_core_genes),
        },
        "Exact_core": {
            "families": exact_core_families_count,
            "genes": len(exact_core_genes),
        },
        "Persistent": {
            "genes": persistent_gene_count,
            "fragmented_genes": persistent_fragmented_genes_count,
            "families": persistent_family_count,
            "families_with_fragments": persistent_fragmented_family_count,
            "families_in_multicopy": orgs_persistent_families_in_multicopy_count,
        },
        "Shell": {
            "genes": shell_gene_count,
            "fragmented_genes": shell_fragmented_genes_count,
            "families": shell_family_count,
            "families_with_fragments": shell_fragmented_family_count,
            "families_in_multicopy": orgs_shell_families_in_multicopy_count,
        },
        "Cloud": {
            "genes": cloud_gene_count,
            "fragmented_genes": cloud_fragmented_genes_count,
            "families": cloud_family_count,
            "families_with_fragments": cloud_fragmented_family_count,
            "families_in_multicopy": orgs_cloud_families_in_multicopy_count,
        },
        "Completeness": completeness,
        "Contamination": contamination,
        "Fragmentation": fragmentation,
        "RGPs": rgp_count,
        "Spots": spot_count,
        "Modules": module_count,
    }
    return summary_info


def write_persistent_duplication_statistics(
    pangenome: Pangenome, output: Path, dup_margin: float, compress: bool
) -> Set[GeneFamily]:
    """
    Writes statistics on persistent duplications in gene families to a specified output file.

    :param pangenome: The Pangenome object containing gene families.
    :param output: The Path specifying the output file location.
    :param dup_margin: The duplication margin used for determining single copy markers.
    :param compress: A boolean indicating whether to compress the output file.

    :return :
    """
    logging.getLogger("PPanGGOLiN").info(
        "Writing statistics on persistent duplication..."
    )

    single_copy_persistent = set()  # Could use bitarrays if speed is needed

    with write_compressed_or_not(
        output / "mean_persistent_duplication.tsv", compress
    ) as outfile:
        fieldnames = [
            "persistent_family",
            "duplication_ratio",
            "mean_presence",
            "is_single_copy_marker",
        ]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for fam in pangenome.gene_families:
            if fam.named_partition == "persistent":
                mean_pres = len(fam) / fam.number_of_organisms
                dup_ratio = fam.duplication_ratio(exclude_fragment=True)
                is_scm = dup_ratio < dup_margin

                if is_scm:
                    single_copy_persistent.add(fam)

                writer.writerow(
                    {
                        "persistent_family": fam.name,
                        "duplication_ratio": round(dup_ratio, 3),
                        "mean_presence": round(mean_pres, 3),
                        "is_single_copy_marker": is_scm,
                    }
                )

    logging.getLogger("PPanGGOLiN").info("Done writing stats on persistent duplication")
    return single_copy_persistent


def write_summaries_in_tsv(
    summaries: List[Dict[str, Any]],
    output_file: Path,
    dup_margin: float,
    soft_core: float,
    compress: bool = False,
):
    """
    Writes summaries of organisms stored in a dictionary into a Tab-Separated Values (TSV) file.

    :param summaries: A list containing organism summaries.
    :param output_file: The Path specifying the output TSV file location.
    :param soft_core: Soft core threshold used
    :param dup_margin: minimum ratio of organisms in which family must have multiple genes to be considered duplicated
    :param compress: Compress the file in .gz
    """
    # Flatten the nested dictionaries within the summaries dictionary
    flat_summaries = [flatten_nested_dict(summary_info) for summary_info in summaries]

    # Create a DataFrame from the flattened summaries
    df_summary = pd.DataFrame(flat_summaries)

    with write_compressed_or_not(output_file, compress) as flout:
        flout.write(f"#soft_core={round(soft_core, 3)}\n")
        flout.write(f"#duplication_margin={round(dup_margin, 3)}\n")

        # Write the DataFrame to a TSV file
        df_summary.to_csv(flout, sep="\t", index=False)


def write_stats(
    output: Path,
    soft_core: float = 0.95,
    dup_margin: float = 0.05,
    compress: bool = False,
):
    """
    Write pangenome statistics for each genomes

    :param output: Path to output directory
    :param soft_core: Soft core threshold to use
    :param dup_margin: minimum ratio of organisms in which family must have multiple genes to be considered duplicated
    :param compress: Compress the file in .gz

    :return: A set containing gene families identified as single copy persistent markers.
    """
    logging.getLogger("PPanGGOLiN").info("Writing pangenome statistics...")

    single_copy_persistent = write_persistent_duplication_statistics(
        pangenome=pan, output=output, dup_margin=dup_margin, compress=compress
    )

    logging.getLogger("PPanGGOLiN").info(
        "Writing genome per genome statistics (completeness and counts)..."
    )

    soft_core_families = pan.soft_core_families(soft_core)
    exact_core_families = pan.exact_core_families()

    pangenome_persistent_single_copy_families = pan.get_single_copy_persistent_families(
        dup_margin=dup_margin, exclude_fragments=True
    )
    assert pangenome_persistent_single_copy_families == single_copy_persistent
    pangenome_persistent_count = len(
        [fam for fam in pan.gene_families if fam.named_partition == "persistent"]
    )
    summaries = []

    for organism in pan.organisms:

        rgp_count = (
            organism.number_of_regions if pan.status["predictedRGP"] != "No" else None
        )
        spot_count = organism.number_of_spots if pan.status["spots"] != "No" else None
        module_count = (
            organism.number_of_modules if pan.status["modules"] != "No" else None
        )

        organism_summary = summarize_genome(
            organism=organism,
            pangenome_persistent_count=pangenome_persistent_count,
            pangenome_persistent_single_copy_families=pangenome_persistent_single_copy_families,
            soft_core_families=soft_core_families,
            exact_core_families=exact_core_families,
            rgp_count=rgp_count,
            spot_count=spot_count,
            module_count=module_count,
        )

        summaries.append(organism_summary)

    write_summaries_in_tsv(
        summaries,
        output_file=output / "genomes_statistics.tsv",
        dup_margin=dup_margin,
        soft_core=soft_core,
        compress=compress,
    )

    logging.getLogger("PPanGGOLiN").info("Done writing genome per genome statistics")


def write_partitions(output: Path, soft_core: float = 0.95):
    """
    Write the list of gene families for each partition

    :param output: Path to output directory
    :param soft_core: Soft core threshold to use
    """
    logging.getLogger("PPanGGOLiN").info(
        "Writing the list of gene families for each partition ..."
    )
    if not os.path.exists(output / "partitions"):
        os.makedirs(output / "partitions")

    part_sets = defaultdict(set)
    # initializing key, value pairs so that files exist even if they are empty
    for needed_key in [
        "soft_core",
        "exact_core",
        "exact_accessory",
        "soft_accessory",
        "persistent",
        "shell",
        "cloud",
    ]:
        part_sets[needed_key] = set()

    for fam in pan.gene_families:
        part_sets[fam.named_partition].add(fam.name)

        # write sub shell partitions
        if fam.partition.startswith("S"):
            part_sets[fam.partition].add(fam.name)

        if fam.number_of_organisms >= pan.number_of_organisms * soft_core:
            part_sets["soft_core"].add(fam.name)
            if fam.number_of_organisms == pan.number_of_organisms:
                part_sets["exact_core"].add(fam.name)
            else:
                part_sets["exact_accessory"].add(fam.name)
        else:
            part_sets["soft_accessory"].add(fam.name)
            part_sets["exact_accessory"].add(fam.name)

    for key, val in part_sets.items():
        with open(output / f"partitions/{key}.txt", "w") as curr_key_file:
            if len(val) > 0:
                curr_key_file.write("\n".join(val) + "\n")

    logging.getLogger("PPanGGOLiN").info(
        "Done writing the list of gene families for each partition"
    )


def write_gene_families_tsv(
    output: Path, compress: bool = False, disable_bar: bool = False
):
    """
    Write the file providing the association between genes and gene families

    :param output: Path to output directory
    :param compress: Compress the file in .gz
    :param disable_bar: Flag to disable progress bar
    """
    logging.getLogger("PPanGGOLiN").info(
        "Writing the file providing the association between genes and gene families..."
    )
    outname = output / f"gene_families.tsv{'.gz' if compress else ''}"
    out_list = []
    for fam in tqdm(
        pan.gene_families,
        total=pan.number_of_gene_families,
        unit="family",
        disable=disable_bar,
    ):
        for gene in fam.genes:
            out_list.append(
                [
                    fam.name,
                    gene.ID,
                    gene.local_identifier,
                    "F" if gene.is_fragment else "",
                ]
            )
    out_df = pd.DataFrame(out_list, columns=["GeneFam", "Gene", "local_id", "is_frag"])
    out_df["count"] = out_df.groupby("GeneFam")["GeneFam"].transform("count")
    out_df = out_df.sort_values(
        by=["count", "Gene", "local_id", "is_frag"], ascending=[False, True, True, True]
    )
    out_df = out_df.drop(columns=["count"])
    out_df.to_csv(
        outname,
        sep="\t",
        index=False,
        header=False,
        compression="infer" if compress else None,
    )
    logging.getLogger("PPanGGOLiN").info(
        "Done writing the file providing the association between genes and "
        f"gene families: '{outname}'"
    )


def summarize_spots(
    spots: set, output: Path, compress: bool = False, file_name="summarize_spots.tsv"
):
    """
    Write a file providing summarize information about hotspots

    :param spots: set of spots in pangenome
    :param output: Path to output directory
    :param compress: Compress the file in .gz
    :patam file_name: Name of the output file
    """

    def r_and_s(value: float):
        """rounds to dp figures and returns a str of the provided value"""
        return str(round(value, 3)) if isinstance(value, float) else str(value)

    file_path = output / file_name

    with write_compressed_or_not(file_path, compress) as fout:
        fout.write(
            "spot\tnb_rgp\tnb_families\tnb_unique_family_sets\tmean_nb_genes\t"
            "stdev_nb_genes\tmax_nb_genes\tmin_nb_genes\n"
        )
        for spot in sorted(spots, key=lambda x: len(x), reverse=True):
            tot_fams = set()
            len_uniq_content = len(spot.get_uniq_content())
            size_list = []
            for rgp in spot.regions:
                tot_fams |= set(rgp.families)
                size_list.append(len(rgp))
            mean_size = mean(size_list)
            stdev_size = stdev(size_list) if len(size_list) > 1 else 0
            max_size = max(size_list)
            min_size = min(size_list)
            fout.write(
                "\t".join(
                    map(
                        r_and_s,
                        [
                            f"{str(spot)}",
                            len(spot),
                            len(tot_fams),
                            len_uniq_content,
                            mean_size,
                            stdev_size,
                            max_size,
                            min_size,
                        ],
                    )
                )
                + "\n"
            )
    logging.getLogger("PPanGGOLiN").info(f"Done writing spots in '{file_path}'")


def write_regions(output: Path, compress: bool = False):
    """
    Write the file providing information about RGP content

    :param output: Path to output directory
    :param compress: Compress the file in .gz
    """

    write_rgp_table(pan.regions, output, compress)


def write_rgp_table(regions: Set[Region], output: Path, compress: bool = False):
    """
    Write the file providing information about regions of genomic plasticity.

    :param regions: Set of Region objects representing regions.
    :param output: Path to the output directory.
    :param compress: Whether to compress the file in .gz format.
    """
    fname = output / "regions_of_genomic_plasticity.tsv"
    with write_compressed_or_not(fname, compress) as tab:
        fieldnames = [
            "region",
            "genome",
            "contig",
            "genes",
            "first_gene",
            "last_gene",
            "start",
            "stop",
            "length",
            "coordinates",
            "score",
            "contigBorder",
            "wholeContig",
        ]

        writer = csv.DictWriter(tab, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        regions = sorted(regions, key=lambda x: (x.organism.name, x.contig.name, x.ID))

        for region in regions:
            row = {
                "region": region.name,
                "genome": region.organism,
                "contig": region.contig,
                "genes": len(region),
                "first_gene": region.starter,
                "last_gene": region.stopper,
                "start": region.start,
                "stop": region.stop,
                "length": region.length,
                "coordinates": region.string_coordinates(),
                "score": region.score,
                "contigBorder": region.is_contig_border,
                "wholeContig": region.is_whole_contig,
            }
            writer.writerow(row)


def spot2rgp(spots: set, output: Path, compress: bool = False):
    """Write a tsv file providing association between spot and rgp

    :param spots: set of spots in pangenome
    :param output: Path to output directory
    :param compress: Compress the file in .gz
    """
    with write_compressed_or_not(output / "spots.tsv", compress) as fout:
        fout.write("spot_id\trgp_id\n")
        for spot in spots:
            for rgp in spot.regions:
                fout.write(f"spot_{spot.ID}\t{rgp.name}\n")


def write_spots(output: Path, compress: bool = False):
    """Write tsv files providing spots information and association with RGP

    :param output: Path to output directory
    :param compress: Compress the file in .gz
    """
    if pan.number_of_spots > 0:
        spot2rgp(pan.spots, output, compress)
        summarize_spots(pan.spots, output, compress)


def write_borders(output: Path, dup_margin: float = 0.05, compress: bool = False):
    """Write all gene families bordering each spot

    :param output: Path to output directory
    :param compress: Compress the file in .gz
    :param dup_margin: minimum ratio of organisms in which family must have multiple genes to be considered duplicated
    """
    multigenics = pan.get_multigenics(dup_margin=dup_margin)
    all_fams = set()
    with write_compressed_or_not(output / "spot_borders.tsv", compress) as fout:
        fout.write("spot_id\tnumber\tborder1\tborder2\n")
        for spot in sorted(pan.spots, key=lambda x: len(x), reverse=True):
            curr_borders = spot.borders(pan.parameters["spot"]["set_size"], multigenics)
            for c, border in curr_borders:
                famstring1 = ",".join([fam.name for fam in border[0]])
                famstring2 = ",".join([fam.name for fam in border[1]])
                all_fams |= set(border[0])
                all_fams |= set(border[1])
                fout.write(f"{spot.ID}\t{c}\t{famstring1}\t{famstring2}\n")

    with write_compressed_or_not(
        output / "border_protein_genes.fasta", compress
    ) as fout:
        for fam in all_fams:
            fout.write(f">{fam.name}\n")
            fout.write(f"{fam.sequence}\n")


def write_module_summary(output: Path, compress: bool = False):
    """
    Write a file providing summarize information about modules

    :param output: Path to output directory
    :param compress: Compress the file in .gz
    """
    logging.getLogger("PPanGGOLiN").info("Writing functional modules summary...")
    with write_compressed_or_not(output / "modules_summary.tsv", compress) as fout:
        fout.write(
            "module_id\tnb_families\tnb_genomes\tpartition\tmean_number_of_occurrence\n"
        )
        for mod in pan.modules:
            org_dict = defaultdict(set)
            partition_counter = Counter()
            for family in mod.families:
                partition_counter[family.named_partition] += 1
                for gene in family.genes:
                    org_dict[gene.organism].add(gene)
            fout.write(
                f"module_{mod.ID}\t{len(mod)}\t{len(org_dict)}\t{partition_counter.most_common(1)[0][0]}\t"
                f"{round((sum([len(genes) for genes in org_dict.values()]) / len(org_dict)) / len(mod), 3)}\n"
            )
        fout.close()

    logging.getLogger("PPanGGOLiN").info(
        f"Done writing module summary: '{output.as_posix() + '/modules_summary.tsv'}'"
    )


def write_modules(output: Path, compress: bool = False):
    """Write a tsv file providing association between modules and gene families

    :param output: Path to output directory
    :param compress: Compress the file in .gz
    """
    logging.getLogger("PPanGGOLiN").info("Writing functional modules...")
    with write_compressed_or_not(output / "functional_modules.tsv", compress) as fout:
        fout.write("module_id\tfamily_id\n")
        for mod in pan.modules:
            for family in mod.families:
                fout.write(f"module_{mod.ID}\t{family.name}\n")
        fout.close()

    logging.getLogger("PPanGGOLiN").info(
        f"Done writing functional modules to: '{output.as_posix() + '/functional_modules.tsv'}'"
    )


def write_org_modules(output: Path, compress: bool = False):
    """Write a tsv file providing association between modules and organisms

    :param output: Path to output directory
    :param compress: Compress the file in .gz
    """
    logging.getLogger("PPanGGOLiN").info("Writing modules to genomes associations...")
    with write_compressed_or_not(output / "modules_in_genomes.tsv", compress) as fout:
        fout.write("module_id\tgenome\tcompletion\n")
        for mod in pan.modules:
            mod_orgs = set()
            for fam in mod.families:
                mod_orgs |= set(fam.organisms)
            for org in mod_orgs:
                completion = len(set(org.families) & set(mod.families)) / len(mod)
                fout.write(f"module_{mod.ID}\t{org.name}\t{completion:.2}\n")
        fout.close()
    logging.getLogger("PPanGGOLiN").info(
        f"Done writing modules to genomes associations to: '{output.as_posix() + '/modules_in_genomes.tsv'}'"
    )


def write_spot_modules(output: Path, compress: bool = False):
    """Write a tsv file providing association between modules and spots

    :param output: Path to output directory
    :param compress: Compress the file in .gz
    """
    logging.getLogger("PPanGGOLiN").info("Writing modules to spot associations...")

    with write_compressed_or_not(output / "modules_spots.tsv", compress) as fout:
        fout.write("module_id\tspot_id\n")

        for spot in pan.spots:
            curr_mods = defaultdict(set)
            for rgp in spot.get_uniq_content():
                for fam in rgp.families:
                    if fam.module is not None:
                        curr_mods[fam.module].add(fam)

            for module, mod_families_in_spot in curr_mods.items():

                if mod_families_in_spot == set(module.families):
                    # if all the families in the module are found in the spot, write the association
                    fout.write(f"module_{module.ID}\tspot_{spot.ID}\n")

    logging.getLogger("PPanGGOLiN").info(
        f"Done writing module to spot associations to: {output.as_posix() + '/modules_spots.tsv'}"
    )


def write_rgp_modules(output: Path, compress: bool = False):
    """Write a tsv file providing association between modules and RGP

    :param output: Path to output directory
    :param compress: Compress the file in .gz
    """
    logging.getLogger("PPanGGOLiN").info("Clustering RGPs based on module content...")

    lists = write_compressed_or_not(output / "modules_RGP_lists.tsv", compress)
    lists.write("representative_RGP\tnb_spots\tmod_list\tRGP_list\n")
    fam2mod = {}
    for mod in pan.modules:
        for fam in mod.families:
            fam2mod[fam] = mod

    region2spot = {}
    for spot in pan.spots:
        for region in spot.regions:
            region2spot[region] = spot

    mod_group2rgps = defaultdict(list)

    for region in pan.regions:
        curr_mod_list = set()
        for fam in region.families:
            mod = fam2mod.get(fam)
            if mod is not None:
                curr_mod_list.add(mod)
        if curr_mod_list != set():
            mod_group2rgps[frozenset(curr_mod_list)].append(region)

    for mod_list, regions in mod_group2rgps.items():
        spot_list = set()
        for region in regions:
            myspot = region2spot.get(region)
            if myspot is not None:
                spot_list.add(region2spot[region])
        lists.write(
            f"{regions[0].name}\t{len(spot_list)}\t{','.join(['module_' + str(mod.ID) for mod in mod_list])}\t"
            f"{','.join([reg.name for reg in regions])}\n"
        )
    lists.close()

    logging.getLogger("PPanGGOLiN").info(
        f"RGP and associated modules are listed in : {output.as_posix() + '/modules_RGP_lists.tsv'}"
    )


def write_pangenome_flat_files(
    pangenome: Pangenome,
    output: Path,
    cpu: int = 1,
    soft_core: float = 0.95,
    dup_margin: float = 0.05,
    csv: bool = False,
    gene_pa: bool = False,
    gexf: bool = False,
    light_gexf: bool = False,
    stats: bool = False,
    json: bool = False,
    partitions: bool = False,
    families_tsv: bool = False,
    regions: bool = False,
    spots: bool = False,
    borders: bool = False,
    modules: bool = False,
    spot_modules: bool = False,
    compress: bool = False,
    disable_bar: bool = False,
):
    """
    Main function to write flat files from pangenome

    :param pangenome: Pangenome object
    :param output: Path to output directory
    :param cpu: Number of available core
    :param soft_core: Soft core threshold to use
    :param dup_margin: minimum ratio of organisms in which family must have multiple genes to be considered duplicated
    :param csv: write csv file format as used by Roary
    :param gene_pa: write gene presence absence matrix
    :param gexf: write pangenome graph in gexf format
    :param light_gexf: write pangenome graph with only gene families
    :param stats: write statistics about pangenome
    :param json: write pangenome graph in json file
    :param partitions: write the gene families for each partition
    :param families_tsv: write gene families information
    :param regions: write RGP information
    :param spots: write information on spots
    :param borders: write gene families bordering spots
    :param modules: write information about modules
    :param spot_modules: write association between modules and RGP and modules and spots
    :param compress: Compress the file in .gz
    :param disable_bar: Disable progress bar
    """
    # TODO Add force parameter to check if output already exist
    if not any(
        x
        for x in [
            csv,
            gene_pa,
            gexf,
            light_gexf,
            stats,
            json,
            partitions,
            spots,
            borders,
            families_tsv,
            modules,
            spot_modules,
            regions,
        ]
    ):
        raise Exception("You did not indicate what file you wanted to write.")

    processes = []
    global pan
    global needAnnotations
    global needFamilies
    global needGraph
    global needPartitions
    global needSpots
    global needRegions
    global needModules
    global needMetadata
    global metatype
    global ignore_err

    pan = pangenome

    if (
        csv
        or gene_pa
        or gexf
        or light_gexf
        or stats
        or json
        or partitions
        or spots
        or families_tsv
        or borders
        or modules
        or spot_modules
        or regions
    ):
        needAnnotations = True
        needFamilies = True
    if stats or partitions or spots or borders:
        needPartitions = True
    if gexf or light_gexf or json or stats:
        needGraph = True
        needRegions = True if pan.status["predictedRGP"] == "inFile" else False
        needSpots = True if pan.status["spots"] == "inFile" else False
        needModules = True if pan.status["modules"] == "inFile" else False
        if pangenome.status["metadata"]["families"] == "inFile":
            needMetadata = True
            metatype = "families"
        else:
            needMetadata = False
    if spots or borders or spot_modules or regions:
        needRegions = True
    if spots or borders or spot_modules:  # or projection:
        needSpots = True
    if modules or spot_modules:  # or projection:
        needModules = True

    check_pangenome_info(
        pangenome,
        need_annotations=needAnnotations,
        need_families=needFamilies,
        need_graph=needGraph,
        need_partitions=needPartitions,
        need_rgp=needRegions,
        need_spots=needSpots,
        need_modules=needModules,
        need_metadata=needMetadata,
        metatypes=[metatype],
        sources=None,
        disable_bar=disable_bar,
    )
    pan.get_org_index()  # make the index because it will be used most likely
    with get_context("fork").Pool(processes=cpu) as p:
        if csv:
            processes.append(
                p.apply_async(
                    func=write_matrix, args=(output, ",", "csv", compress, True)
                )
            )
        if gene_pa:
            processes.append(
                p.apply_async(func=write_gene_presence_absence, args=(output, compress))
            )
        if gexf:
            processes.append(
                p.apply_async(func=write_gexf, args=(output, False, compress))
            )
        if light_gexf:
            processes.append(
                p.apply_async(func=write_gexf, args=(output, True, compress))
            )
        if stats:
            processes.append(
                p.apply_async(
                    func=write_stats, args=(output, soft_core, dup_margin, compress)
                )
            )
        if json:
            processes.append(p.apply_async(func=write_json, args=(output, compress)))
        if partitions:
            processes.append(
                p.apply_async(func=write_partitions, args=(output, soft_core))
            )
        if families_tsv:
            processes.append(
                p.apply_async(
                    func=write_gene_families_tsv, args=(output, compress, disable_bar)
                )
            )
        if spots:
            processes.append(p.apply_async(func=write_spots, args=(output, compress)))
        if regions:
            processes.append(p.apply_async(func=write_regions, args=(output, compress)))
        if borders:
            processes.append(
                p.apply_async(func=write_borders, args=(output, dup_margin, compress))
            )
        if modules:
            processes.append(p.apply_async(func=write_modules, args=(output, compress)))
            processes.append(
                p.apply_async(func=write_module_summary, args=(output, compress))
            )
            processes.append(
                p.apply_async(func=write_org_modules, args=(output, compress))
            )
        if spot_modules:
            processes.append(
                p.apply_async(func=write_spot_modules, args=(output, compress))
            )
            processes.append(
                p.apply_async(func=write_rgp_modules, args=(output, compress))
            )

        for process in processes:
            process.get()  # get all the results


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    mk_outdir(args.output, args.force)
    global pan
    pan.add_file(args.pangenome)
    write_pangenome_flat_files(
        pan,
        args.output,
        cpu=args.cpu,
        soft_core=args.soft_core,
        dup_margin=args.dup_margin,
        csv=args.csv,
        gene_pa=args.Rtab,
        gexf=args.gexf,
        light_gexf=args.light_gexf,
        stats=args.stats,
        json=args.json,
        partitions=args.partitions,
        families_tsv=args.families_tsv,
        regions=args.regions,
        spots=args.spots,
        borders=args.borders,
        modules=args.modules,
        spot_modules=args.spot_modules,
        compress=args.compress,
        disable_bar=args.disable_prog_bar,
    )


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser(
        "write_pangenome", formatter_class=argparse.RawTextHelpFormatter
    )
    parser_flat(parser)
    return parser


def parser_flat(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of write command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(
        title="Required arguments",
        description="One of the following arguments is required :",
    )
    required.add_argument(
        "-p", "--pangenome", required=False, type=Path, help="The pangenome .h5 file"
    )
    required.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output directory where the file(s) will be written",
    )
    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument(
        "--soft_core",
        required=False,
        type=restricted_float,
        default=0.95,
        help="Soft core threshold to use",
    )

    optional.add_argument(
        "--dup_margin",
        required=False,
        type=restricted_float,
        default=0.05,
        help="minimum ratio of genomes in which the family must have multiple genes "
        "for it to be considered 'duplicated'",
    )

    optional.add_argument(
        "--gexf",
        required=False,
        action="store_true",
        help="write a gexf file with all the annotations and all the genes of each gene family",
    )
    optional.add_argument(
        "--light_gexf",
        required=False,
        action="store_true",
        help="write a gexf file with the gene families and basic information about them",
    )

    optional.add_argument(
        "--json",
        required=False,
        action="store_true",
        help="Writes the graph in a json file format",
    )

    optional.add_argument(
        "--csv",
        required=False,
        action="store_true",
        help="csv file format as used by Roary, among others. "
        "The alternative gene ID will be the partition, if there is one",
    )
    optional.add_argument(
        "--Rtab",
        required=False,
        action="store_true",
        help="tabular file for the gene binary presence absence matrix",
    )

    optional.add_argument(
        "--stats",
        required=False,
        action="store_true",
        help="tsv files with some statistics for each each gene family",
    )

    optional.add_argument(
        "--partitions",
        required=False,
        action="store_true",
        help="list of families belonging to each partition, with one file per partitions and "
        "one family per line",
    )

    optional.add_argument(
        "--families_tsv",
        required=False,
        action="store_true",
        help="Write a tsv file providing the association between genes and gene families",
    )

    optional.add_argument(
        "--regions",
        required=False,
        action="store_true",
        help="Writes the predicted RGP and descriptive metrics in 'plastic_regions.tsv'",
    )
    optional.add_argument(
        "--spots",
        required=False,
        action="store_true",
        help="Write spot summary and a list of all RGP in each spot",
    )
    optional.add_argument(
        "--borders",
        required=False,
        action="store_true",
        help="List all borders of each spot",
    )
    optional.add_argument(
        "--modules",
        required=False,
        action="store_true",
        help="Write a tsv file listing functional modules and the families that belong to them",
    )
    optional.add_argument(
        "--spot_modules",
        required=False,
        action="store_true",
        help="writes 2 files comparing the presence of modules within spots",
    )

    optional.add_argument(
        "--compress",
        required=False,
        action="store_true",
        help="Compress the files in .gz",
    )
    optional.add_argument(
        "-c",
        "--cpu",
        required=False,
        default=1,
        type=int,
        help="Number of available cpus",
    )


if __name__ == "__main__":
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser_flat(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
