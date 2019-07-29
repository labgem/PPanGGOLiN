#!/usr/bin/env python3
#coding:utf-8

# default libraries
import logging
from collections import Counter

#external libraries
from tqdm import tqdm

#local libraries
from ppanggolin.utils import read_compressed_or_not, get_num_lines, getCurrentRAM
from ppanggolin.genome import Gene, Organism

def linkGeneFamsToGenes(gene2families, pangenome):
    """
        Reads a dictionnary with gene IDs as key and gene families objects as values, and a pangenome whose annotations have already been filled, and link the annotations with the gene families.
    """
    if len(pangenome.organisms) > 0:#check if annotations are actually there
        logging.getLogger().info("Linking genes to gene families...")
        bar = tqdm(pangenome.organisms, desc = "Linking genes to gene families", unit = "organism")
        for org in bar:
            for contig in org.contigs:
                for gene in contig.genes:
                    fam, is_frag = gene2families[gene.ID]
                    gene.fill_family(fam, is_frag)
        bar.close()
    else:
        raise Exception("The pangenome does not seem to have annotations.")

def read_org_line(pangenome, organism, gff_file_path, circular_contigs, lim_occurence, infer_singletons, add_rna_to_the_pangenome, gene2families ):
    (GFF_seqname, _, GFF_type, GFF_start, GFF_end, _, GFF_strand, _, GFF_attribute) = range(0,9)#missing values : source, score, frame. They are unused.
    def getGffAttributes(gff_fields):
        """
            Parses the gff attribute's line and outputs the attributes in a dict structure.
            :param gff_fields: a gff line stored as a list. Each element of the list is a column of the gff.
            :type list:
            :return: attributes:
            :rtype: dict
        """
        attributes_field = [f for f in gff_fields[GFF_attribute].strip().split(';') if len(f)>0]
        attributes = {}
        for att in attributes_field:
            (key, value) = att.strip().split('=')
            attributes[key.upper()]=value
        return attributes

    def getIDAttribute(attributes):
        """
            Gets the ID of the element from which the provided attributes were extracted. Raises an error if no ID is found.
            :param attribute:
            :type dict:
            :return: ElementID:
            :rtype: string
        """
        ElementID = attributes.get("ID")
        if not ElementID:
            logging.getLogger().error("Each CDS type of the gff files must own a unique ID attribute. Not the case for file: "+gff_file_path)
            exit(1)
        return ElementID
    
    organism = pangenome.addOrganism(organism)
    cpt_fam_occ = Counter()
    with read_compressed_or_not(gff_file_path) as gff_file:
        for line in gff_file:
            if line.startswith('##',0,2):
                if line.startswith('FASTA',2,7):
                    break
                elif line.startswith('sequence-region',2,17):
                    fields = [el.strip() for el in line.split()]
                    
                    contig = organism.addContig(fields[1], True if fields[1] in circular_contigs else False)
                continue
            if line.startswith('#!',0,2):## special refseq comment lines for versionning softs, assemblies and annotations.
                continue
            gff_fields = [el.strip() for el in line.split('\t')]
            attributes = getGffAttributes(gff_fields)

            if gff_fields[GFF_type] == 'region':
                if gff_fields[GFF_seqname] in circular_contigs:
                    contig.is_circular = True
            elif gff_fields[GFF_type] == 'CDS' or (add_rna_to_the_pangenome and gff_fields[GFF_type].find("RNA")):
                protein = getIDAttribute(attributes)
                family = gene2families.get(protein)
                if not family:## the protein id was not found under "ID", searching elsewhere
                    proteinID = attributes.get("PROTEIN_ID")
                    if proteinID:
                        protein = proteinID
                        family = gene2families.get(protein)

                if not family:
                    if infer_singletons:## if we did not find the associated family at some point above, and if infer_singletons is set
                        family = pangenome.addGeneFamily(protein)
                        gene2families[protein] = (family,False)
                        logging.getLogger().info("infered singleton: "+protein)
                    else:
                        raise KeyError("Unknown families:"+protein+ ", check your families file or run again the program using the option to infer singleton")

                try:
                    name = attributes.pop('NAME')
                except KeyError:
                    try:
                        name = attributes.pop('GENE')
                    except KeyError:
                        name = ""

                try:
                    product = attributes.pop('PRODUCT')
                except KeyError:
                    product = ""
                gene = Gene(protein)
                gene.fill_family(family[0], is_fragment=family[1])
                #here contig is filled in order, so position is the number of genes already stored in the contig.
                gene.fill_annotations(start = int(gff_fields[GFF_start]),
                                    stop = int(gff_fields[GFF_end]),
                                    strand =gff_fields[GFF_strand],
                                    geneType = gff_fields[GFF_type],
                                    position = len(contig.genes),
                                    name = name,
                                    product=product )
                contig.addGene(gene)
                cpt_fam_occ[family]+=1

    if (lim_occurence > 0):
        fam_to_remove = [key for key, val in cpt_fam_occ if val > lim_occurence ]
        logging.getLogger().debug("highly repeated families found (>"+str(lim_occurence)+" in "+organism+"): "+" ".join(fam_to_remove))
        pangenome.removeFamsFrom(fam_to_remove)

def readAnnotations(pangenome, organisms_file, lim_occurence, infer_singletons, add_rna_to_the_pangenome, gene2families):
    logging.getLogger().info("Reading "+organisms_file+" the list of organism files ...")
    bar = tqdm(read_compressed_or_not(organisms_file),total=get_num_lines(organisms_file), unit = "gff file")

    for line in bar:
        elements = [el.strip() for el in line.split("\t")]
        if len(elements)<=1:
            logging.getLogger().error("No tabulation separator found in organisms file")
            exit(1)
        bar.set_description("Processing "+elements[1].split("/")[-1])
        bar.refresh()
        
        #multi? => If pangenome is not given only(so no infer_singletons) (maybe)
        read_org_line(pangenome, elements[0], elements[1], elements[2:], lim_occurence, infer_singletons, add_rna_to_the_pangenome, gene2families)
    bar.close()
    logging.getLogger().debug(f"RAM at the end of readAnnotations: {getCurrentRAM()}")

def readClustering(pangenome, families_tsv_file):
    """
        Creates the pangenome, the gene families and the genes with an associated gene family.
        Reads a families tsv file from mmseqs2 output and adds the gene families and the genes to the pangenome.
    """
    logging.getLogger().info("Reading "+families_tsv_file+" the gene families file ...")

    families_tsv_file = read_compressed_or_not(families_tsv_file)

    gene2families = {}
    for line in families_tsv_file:
        elements = [el.strip() for el in line.split()] # 2 or 3 fields expected
        if len(elements)<=1:
            logging.getLogger().error("No tabulation separator found in gene families file")
            exit(1)
        (fam_id, gene_id, is_frag) = elements if len(elements) == 3 else elements+[None]
        geneFam = pangenome.addGeneFamily(fam_id)
        gene2families[gene_id] = (geneFam, True if is_frag == "F" else False)
    families_tsv_file.close()

    return gene2families

def computeNeighborsGraph(pangenome):
    """
        Creates the Pangenome Graph with the loaded information.
    """
    bar = tqdm(pangenome.organisms, total = len(pangenome.organisms), unit = "organism")
    for org in bar:
        bar.set_description(f"Processing {org.name}")
        bar.refresh()
        for contig in org.contigs:
            prev = None
            for gene in contig.genes:
                if prev is not None and not gene.family.removed:
                    if not (prev.family == gene.family and (prev.is_fragment or gene.is_fragment)):
                        pangenome.addEdge(gene, prev, org)
                prev = gene
            if contig.is_circular:
                pangenome.addEdge(contig.genes[0],prev, org)

