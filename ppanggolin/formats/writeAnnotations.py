#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from typing import Dict, Tuple, Union

# installed libraries
from tqdm import tqdm
import tables

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Gene, RNA
from ppanggolin.formats.readBinaries import Genedata


genedata_counter = 0


def get_max_len_annotations(pangenome: Pangenome) -> Tuple[int, int, int, int, int]:
    """
    Get the maximum size of each annotation information to optimize disk space

    :param pangenome: Annotated pangenome

    :return: Maximum size of each annotation
    """
    max_org_len, max_contig_len, max_gene_id_len, max_rna_id_len, max_gene_local_id = 1, 1, 1, 1, 1
    for org in pangenome.organisms:
        if len(org.name) > max_org_len:
            max_org_len = len(org.name)
        for contig in org.contigs:
            if len(contig.name) > max_contig_len:
                max_contig_len = len(contig.name)
            for gene in contig.genes:
                if len(gene.ID) > max_gene_id_len:
                    max_gene_id_len = len(gene.ID)
                if len(gene.local_identifier) > max_gene_local_id:
                    max_gene_local_id = len(gene.local_identifier)
            for rna in contig.RNAs:
                if len(rna.ID) > max_rna_id_len:
                    max_rna_id_len = len(rna.ID)

    return max_org_len, max_contig_len, max_gene_id_len, max_rna_id_len, max_gene_local_id


def organism_desc(org_len: int, contig_len: int) -> Dict[str, tables.StringCol]:
    """
    Table description to save organism-related information

    :param org_len: Maximum size of organism name.
    :param contig_len: Maximum size of contigs name

    :return: Formatted table
    """
    return {'name': tables.StringCol(itemsize=org_len),
            "contig": tables.StringCol(itemsize=contig_len)}


def write_organisms(pangenome: Pangenome, h5f: tables.File, annotation: tables.Group,
                    organism_desc: Dict[str, tables.StringCol], disable_bar=False):
    """Write organisms information in the pangenome file

    :param pangenome: Annotated pangenome object
    :param h5f: Pangenome file
    :param annotation: Annotation table group
    :param organism_desc: Organisms table description.
    :param disable_bar: Allow disabling progress bar
    """
    organism_table = h5f.create_table(annotation, "genomes", organism_desc,
                                      expectedrows=pangenome.number_of_organisms)
    logging.getLogger("PPanGGOLiN").debug(f"Writing {pangenome.number_of_organisms} genomes")
    organism_row = organism_table.row
    for org in tqdm(pangenome.organisms, total=pangenome.number_of_organisms, unit="genome", disable=disable_bar):
        for contig in org.contigs:
            organism_row["name"] = org.name
            organism_row["contig"] = contig.name
            organism_row.append()
    organism_table.flush()


def contig_desc(contig_len: int, max_gene_id_len: int,
                max_rna_id_len: int) -> Dict[str, Union[tables.StringCol, tables.BoolCol, tables.UInt32Col]]:
    """Table description to save contig-related information

    :param contig_len: Maximum size of contig name
    :param max_gene_id_len: Maximum size of gene name
    :param max_rna_id_len: Maximum size of rna name

    :return: Formatted table
    """
    return {'name': tables.StringCol(itemsize=contig_len),
            "is_circular": tables.BoolCol(dflt=False),
            'length': tables.UInt32Col(),
            "gene": tables.StringCol(itemsize=max_gene_id_len),
            "rna": tables.StringCol(itemsize=max_rna_id_len)}


def write_contigs(pangenome: Pangenome, h5f: tables.File, annotation: tables.Group,
                  contig_desc: Dict[str, Union[tables.StringCol, tables.BoolCol, tables.UInt32Col]],
                  disable_bar=False):
    """Write contigs information in the pangenome file
    :param pangenome: Annotated pangenome object
    :param h5f: Pangenome file
    :param annotation: Annotation table group
    :param contig_desc: Contigs table description
    :param disable_bar: Allow disabling progress bar
    """
    contig_table = h5f.create_table(annotation, "contigs", contig_desc, expectedrows=pangenome.number_of_contigs)
    logging.getLogger("PPanGGOLiN").debug(f"Writing {pangenome.number_of_contigs} contigs")
    contig_row = contig_table.row
    for contig in tqdm(pangenome.contigs, total=pangenome.number_of_contigs, unit="contigs", disable=disable_bar):
        if contig.number_of_genes >= contig.number_of_rnas:
            rna_list = list(contig.RNAs)
            for index, gene in enumerate(contig.genes):
                contig_row["name"] = contig.name
                contig_row["is_circular"] = contig.is_circular
                contig_row["length"] = len(contig)
                contig_row["gene"] = gene.ID
                if index < len(rna_list):
                    contig_row["rna"] = rna_list[index].ID
                contig_row.append()
        else:
            gene_list = list(contig.genes)
            for index, rna in enumerate(contig.RNAs):
                contig_row["name"] = contig.name
                contig_row["is_circular"] = contig.is_circular
                contig_row["rna"] = rna.ID
                if index < len(gene_list):
                    contig_row["gene"] = gene_list[index].ID
                contig_row.append()
    contig_table.flush()


def gene_desc(id_len, max_local_id) -> Dict[str, Union[tables.StringCol, tables.UInt32Col, tables.BoolCol]]:
    """Table description to save gene-related information

    :param id_len: Maximum size of gene name
    :param max_local_id: Maximum size of gene local identifier

    :return: Formatted table
    """
    return {'ID': tables.StringCol(itemsize=id_len),
            'genedata_id': tables.UInt32Col(),
            'local': tables.StringCol(itemsize=max_local_id),
            'is_fragment': tables.BoolCol(dflt=False)}


def write_genes(pangenome: Pangenome,  h5f: tables.File, annotation: tables.Group,
                gene_desc: Dict[str, Union[tables.StringCol, tables.UInt32Col, tables.BoolCol]],
                disable_bar=False) -> Dict[Genedata, int]:
    """Write genes information in the pangenome file

    :param pangenome: Annotated pangenome object
    :param h5f: Pangenome file
    :param annotation: Annotation table group
    :param gene_desc: Genes table description
    :param disable_bar: Allow to disable progress bar

    :returns: Dictionnary linking genedata to gene identifier
    """
    global genedata_counter
    genedata2gene = {}
    gene_table = h5f.create_table(annotation, "genes", gene_desc, expectedrows=pangenome.number_of_genes)
    logging.getLogger("PPanGGOLiN").debug(f"Writing {pangenome.number_of_genes} genes")
    gene_row = gene_table.row
    for gene in tqdm(pangenome.genes, total=pangenome.number_of_genes, unit="gene", disable=disable_bar):
        gene_row["ID"] = gene.ID
        gene_row["is_fragment"] = gene.is_fragment
        gene_row["local"] = gene.local_identifier
        genedata = get_genedata(gene)
        genedata_id = genedata2gene.get(genedata)
        if genedata_id is None:
            genedata_id = genedata_counter
            genedata2gene[genedata] = genedata_id
            genedata_counter += 1
        gene_row["genedata_id"] = genedata_id
        gene_row.append()
    gene_table.flush()
    return genedata2gene


def rna_desc(id_len) -> Dict[str, Union[tables.StringCol, tables.UInt32Col]]:
    """Table description to save rna-related information

    :param id_len: Maximum size of RNA identifier

    :return: Formatted table
    """
    return {'ID': tables.StringCol(itemsize=id_len),
            'genedata_id': tables.UInt32Col()}


def write_rnas(pangenome: Pangenome,  h5f: tables.File, annotation: tables.Group,
               rna_desc: Dict[str, Union[tables.StringCol, tables.UInt32Col]],
               disable_bar=False) -> Dict[Genedata, int]:
    """Write RNAs information in the pangenome file

    :param pangenome: Annotated pangenome object
    :param h5f: Pangenome file
    :param annotation: Annotation table group
    :param rna_desc: RNAs table description
    :param disable_bar: Allow to disable progress bar

    :returns: Dictionnary linking genedata to RNA identifier
    """
    global genedata_counter
    genedata2rna = {}
    rna_table = h5f.create_table(annotation, "RNAs", rna_desc, expectedrows=pangenome.number_of_genes)
    logging.getLogger("PPanGGOLiN").debug(f"Writing {pangenome.number_of_genes} genes")
    rna_row = rna_table.row
    for rna in tqdm(pangenome.RNAs, total=pangenome.number_of_rnas, unit="RNA", disable=disable_bar):
        rna_row["ID"] = rna.ID
        genedata = get_genedata(rna)
        genedata_id = genedata2rna.get(genedata)
        if genedata_id is None:
            genedata_id = genedata_counter
            genedata2rna[genedata] = genedata_id
            genedata_counter += 1
        rna_row["genedata_id"] = genedata_id
        rna_row.append()
    rna_table.flush()
    return genedata2rna


def genedata_desc(type_len, name_len, product_len):
    """
    Creates a table for gene-related data

    :param type_len: Maximum size of gene Type.
    :param name_len: Maximum size of gene name
    :param product_len: Maximum size of gene product
    :return: Formatted table for gene metadata
    """
    return {
        'genedata_id': tables.UInt32Col(),
        'start': tables.UInt32Col(),
        'stop': tables.UInt32Col(),
        'strand': tables.StringCol(itemsize=1),
        'gene_type': tables.StringCol(itemsize=type_len),
        'position': tables.UInt32Col(),
        'name': tables.StringCol(itemsize=name_len),
        'product': tables.StringCol(itemsize=product_len),
        'genetic_code': tables.UInt32Col(dflt=11),
    }


def get_max_len_genedata(pangenome: Pangenome) -> Tuple[int, int, int]:
    """
    Get the maximum size of each gene data information to optimize disk space

    :param pangenome: Annotated pangenome
    :return: maximum size of each annotation
    """
    max_name_len = 1
    max_product_len = 1
    max_type_len = 1
    for org in pangenome.organisms:
        for contig in org.contigs:
            for gene in contig.genes:
                if len(gene.name) > max_name_len:
                    max_name_len = len(gene.name)
                if len(gene.product) > max_product_len:
                    max_product_len = len(gene.product)
                if len(gene.type) > max_type_len:
                    max_type_len = len(gene.type)
            for gene in contig.RNAs:
                if len(gene.name) > max_name_len:
                    max_name_len = len(gene.name)
                if len(gene.product) > max_product_len:
                    max_product_len = len(gene.product)
                if len(gene.type) > max_type_len:
                    max_type_len = len(gene.type)

    return max_type_len, max_name_len, max_product_len


def get_genedata(feature: Union[Gene, RNA]) -> Genedata:
    """
    Gets the genedata type of Feature

    :param feature: Gene or RNA object

    :return: Tuple with a Feature associated data
    """
    position = None
    genetic_code = 11
    if isinstance(feature, Gene):
        position = feature.position
        genetic_code = feature.genetic_code
    return Genedata(feature.start, feature.stop, feature.strand, feature.type, position, feature.name,
                    feature.product, genetic_code)


def write_genedata(pangenome: Pangenome,  h5f: tables.File, annotation: tables.Group,
                   genedata2gene: Dict[Genedata, int], disable_bar=False):
    """Writting genedata information in pangenome file

    :param pangenome: Pangenome object filled with annotation.
    :param h5f: Pangenome file
    :param annotation: Annotation group in Table
    :param genedata2gene: Dictionnary linking genedata to gene identifier.
    :param disable_bar: Allow disabling progress bar
    """
    try:
        genedata_table = annotation.genedata
    except tables.exceptions.NoSuchNodeError:
        genedata_table = h5f.create_table(annotation, "genedata", genedata_desc(*get_max_len_genedata(pangenome)),
                                          expectedrows=len(genedata2gene))
    logging.getLogger("PPanGGOLiN").debug(f"Writing {len(genedata2gene)} gene-related data "
                                          "(can be lower than the number of genes)")
    genedata_row = genedata_table.row
    for genedata, genedata_id in tqdm(genedata2gene.items(), unit="genedata", disable=disable_bar):
        genedata_row["genedata_id"] = genedata_id
        genedata_row["start"] = genedata.start
        genedata_row["stop"] = genedata.stop
        genedata_row["strand"] = genedata.strand
        genedata_row["gene_type"] = genedata.gene_type
        if genedata.gene_type == "CDS":
            genedata_row["position"] = genedata.position
            genedata_row["genetic_code"] = genedata.genetic_code
        genedata_row["name"] = genedata.name
        genedata_row["product"] = genedata.product
        genedata_row.append()
    genedata_table.flush()


def write_annotations(pangenome: Pangenome, h5f: tables.File, rec_organisms: bool = True,
                      rec_contigs: bool = True, rec_genes: bool = True, rec_rnas: bool = True, disable_bar: bool = False):
    """Function writing all the pangenome annotations

    :param pangenome: Annotated pangenome
    :param h5f: Pangenome HDF5 file
    :param rec_organisms: Allow writing organisms in pangenomes
    :param rec_contigs: Allow writing contigs in pangenomes
    :param rec_genes: Allow writing genes in pangenomes
    :param rec_rnas: Allow writing RNAs in pangenomes
    :param disable_bar: Alow to disable progress bar
    """
    annotation = h5f.create_group("/", "annotations", "Annotations of the pangenome organisms")

    org_len, contig_len, gene_id_len, rna_id_len, gene_local_id = get_max_len_annotations(pangenome)

    # I add these boolean in case we would one day only load organism, contig or genes, without the other.

    if rec_organisms:
        desc = organism_desc(org_len, contig_len)
        write_organisms(pangenome, h5f, annotation, desc, disable_bar)
    if rec_contigs:
        desc = contig_desc(contig_len, gene_id_len, rna_id_len)
        write_contigs(pangenome, h5f, annotation, desc, disable_bar)
    if rec_genes:
        desc = gene_desc(gene_id_len, gene_local_id)
        genedata2gene = write_genes(pangenome, h5f, annotation, desc, disable_bar)
        write_genedata(pangenome, h5f, annotation, genedata2gene, disable_bar)
    if rec_rnas:
        desc = rna_desc(rna_id_len)
        genedata2rna = write_rnas(pangenome, h5f, annotation, desc, disable_bar)
        write_genedata(pangenome, h5f, annotation, genedata2rna, disable_bar)


def get_gene_sequences_len(pangenome: Pangenome) -> Tuple[int, int]:
    """
    Get the maximum size of gene sequences to optimize disk space
    :param pangenome: Annotated pangenome
    :return: maximum size of each annotation
    """
    max_gene_id_len = 1
    max_gene_type = 1
    for gene in pangenome.genes:
        if len(gene.ID) > max_gene_id_len:
            max_gene_id_len = len(gene.ID)
        if len(gene.type) > max_gene_type:
            max_gene_type = len(gene.type)
    return max_gene_id_len, max_gene_type


def gene_sequences_desc(gene_id_len, gene_type_len) -> dict:
    """
    Create table to save gene sequences
    :param gene_id_len: Maximum size of gene sequence identifier
    :param gene_type_len: Maximum size of gene type
    :return: Formated table
    """
    return {
        "gene": tables.StringCol(itemsize=gene_id_len),
        "seqid": tables.UInt32Col(),
        "type": tables.StringCol(itemsize=gene_type_len)
    }


def get_sequence_len(pangenome: Pangenome) -> int:
    """
    Get the maximum size of gene sequences to optimize disk space
    :param pangenome: Annotated pangenome
    :return: maximum size of each annotation
    """
    max_seq_len = 1
    for gene in pangenome.genes:
        if len(gene.dna) > max_seq_len:
            max_seq_len = len(gene.dna)
    return max_seq_len


def sequence_desc(max_seq_len: int) -> dict:
    """
    Table description to save sequences
    :param max_seq_len: Maximum size of gene type
    :return: Formated table
    """
    return {
        "seqid": tables.UInt32Col(),
        "dna": tables.StringCol(itemsize=max_seq_len)
    }


def write_gene_sequences(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Function writing all the pangenome gene sequences
    :param pangenome: Pangenome with gene sequences
    :param h5f: Pangenome HDF5 file without sequences
    :param disable_bar: Disable progress bar
    """
    gene_seq = h5f.create_table("/annotations", "geneSequences", gene_sequences_desc(*get_gene_sequences_len(pangenome)),
                                expectedrows=pangenome.number_of_genes)
    # process sequences to save them only once
    seq2seqid = {}
    id_counter = 0
    gene_row = gene_seq.row
    for gene in tqdm(sorted(pangenome.genes, key=lambda x: x.ID), total=pangenome.number_of_genes, unit="gene",
                     disable=disable_bar):
        curr_seq_id = seq2seqid.get(gene.dna)
        if curr_seq_id is None:
            curr_seq_id = id_counter
            seq2seqid[gene.dna] = id_counter
            id_counter += 1
        gene_row["gene"] = gene.ID
        gene_row["seqid"] = curr_seq_id
        gene_row["type"] = gene.type
        gene_row.append()
    gene_seq.flush()

    seq_table = h5f.create_table("/annotations", "sequences", sequence_desc(get_sequence_len(pangenome)),
                                 expectedrows=len(seq2seqid))

    seq_row = seq_table.row
    for seq, seqid in seq2seqid.items():
        seq_row["dna"] = seq
        seq_row["seqid"] = seqid
        seq_row.append()
    seq_table.flush()
