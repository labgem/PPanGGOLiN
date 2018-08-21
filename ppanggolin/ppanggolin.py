#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-

import warnings
warnings.filterwarnings("ignore")
from collections import defaultdict, OrderedDict, Counter, deque
from bidict import bidict
from ordered_set import OrderedSet
import networkx as nx
import logging
import sys
import math
from time import time, sleep
import os
import shutil
import gzip
import tempfile
from tqdm import tqdm
from random import sample
from multiprocessing import Pool, Semaphore
import contextlib
from nem import *
from .utils import *
import pdb
from fa2 import ForceAtlas2
import plotly.plotly as py
import plotly.offline as out_plotly
import plotly.graph_objs as go
from ascii_graph import Pyasciigraph

(TYPE, FAMILY, START, END, STRAND, NAME, PRODUCT) = range(0, 7)#data index in annotation
(ORGANISM_INDEX,CONTIG_INDEX,POSITION_INDEX) = range(0, 3)#index
(ORGANISM_ID, ORGANISM_GFF_FILE) = range(0, 2)#data index in the file listing organisms 
(GFF_seqname, GFF_source, GFF_feature, GFF_start, GFF_end, GFF_score, GFF_strand, GFF_frame, GFF_attribute) = range(0,9) 
(MU,EPSILON,PROPORTION) = range(0, 3)
(FAMILIES_PARTITION,PARTITION_PARAMETERS) = range(0, 2)
RESERVED_WORDS = set(["id", "label", "name", "weight", "partition", "partition_exact", "partition_soft", "length", "length_min", "length_max", "length_avg", "length_med", "product", 'nb_genes','subpartition_shell',"viz"])
SHORT_TO_LONG = {'EA':'exact_accessory','EC':'exact_core','SA':'soft_accessory','SC':'soft_core','P':'persistent','S':'shell','C':'cloud','U':'undefined'}
COLORS = {"pangenome":"black", "exact_accessory":"#EB37ED", "exact_core" :"#FF2828", "soft_core":"#e6e600", "soft_accessory":"#996633","shell": "#00D860", "persistent":"#F7A507", "cloud":"#79DEFF", "undefined":"#828282"}
COLORS_RGB = {"pangenome":{'r': 0, 'g': 0, 'b': 0, 'a': 0}, "exact_accessory":{'r': 235, 'g': 55, 'b': 237, 'a': 0}, "exact_core" :{'r': 255, 'g': 40, 'b': 40, 'a': 0},  "soft_core":{'r': 255, 'g': 255, 'b': 0, 'a': 0}, "soft_accessory": {'r': 153, 'g': 102, 'b': 51, 'a': 0},"shell": {'r': 0, 'g': 216, 'b': 96, 'a': 0}, "persistent":{'r': 247, 'g': 165, 'b': 7, 'a': 0}, "cloud":{'r': 121, 'g': 222, 'b': 255, 'a': 0}, "undefined":{'r': 130, 'g': 130, 'b': 130, 'a': 0}}

"""
    :mod:`ppanggolin` -- Depict microbial diversity
===================================

.. module:: ppanggolin
   :platform: Unix
   :synopsis: Depict microbial diversity via a partionned pangenome graph.
    .. moduleauthor:: Guillaume GAUTREAU (LABGeM, Genoscope, France) ggautrea@genoscope.cns.fr

    Description
    -------------------
    Pangenomes are generally stored in a binary matrix denoting the presence or absence of each gene family across organisms. However, this structure does not handle the genomic organization of gene families in each organism. We propose a graph model where nodes represent families and edges chromosomal neighborhood information. Indeed, it is known that core gene families share conserved organizations whereas variable regions are rather randomly distributed along genomes. Moreover, our method classifies gene families through an Expectation/Maximization algorithm based on Bernoulli mixture model. This approach splits pangenomes in three groups: (1) persistent genome, equivalent to a relaxed core genome (genes conserved in all but a few genomes); (2) shell genome, genes having intermediate frequencies corresponding to moderately conserved genes potentially associated to environmental adaptation capabilities; (3) cloud genome, genes found at very low frequency.
""" 

class PPanGGOLiN:

    """  
          The ``PPanGGOLiN`` class
        ======================
        .. class:: PPanGGOLiN

            Pangenomes are generally stored in a binary matrix denoting the presence or absence of each gene family across organisms. 
            However, this structure does not handle the genomic organization of gene families in each organism. 
            We propose a graph model where nodes represent families and edges chromosomal neighborhood information. 
            Indeed, it is known that core gene families share conserved organizations whereas variable regions are rather randomly distributed along genomes.
            The PPanGGOLiN class models the genemic diversity of a pangenome, this modelisation organize the genemic diveristy via a pangenome graph of chromosomal neigborhood.
            Moreover, our method classifies gene families through an Expectation/Maximization algorithm based on Bernoulli mixture model. 
            This approach splits pangenomes in three groups: 
                1. *persistent genome*, equivalent to a relaxed core genome (genes conserved in all but a few genomes); 
                2. *shell genome*, genes having intermediate frequencies corresponding to moderately conserved genes potentially associated to environmental adaptation capabilities; 
                3. *cloud genome*, genes found at very low frequency.

            .. attribute:: annotations

                multilevel dictionnaries containing a dictionary of contig for each organism, and a dictionary of lists containing annotations for each contig

            .. attribute:: neighbors_graph

                a networkx graph. Node correspond to gene families and edges to chromosomal colocalization beween families. Organisms supporting each edge are stored in edge attribute as weel as the edge weight (number of organism coverinf each edge).
                Nodes attributes contains the gene identifiers of each organism supporting this node.

            .. attribute:: organisms

                an ordored-set contains the imported organisms 

            .. attribute:: nb_organisms

                a int giving the number of imported organisms 

            .. attribute:: circular_contig_size

                a dict containing the size of the contigs (contigs identifiers as keys) which are both well assembled and circular contigs (used to circularize the graph). The contigs wich are not circular are not in this dictionnaries.

            .. attribute:: families_repeted_th
   
                a int containing the threshold of the maximum number copy of each families. Families exceeding this threshold are removed and are listed in the next attribute.

            .. attribute:: families_repeted

                a set containing the family identifiers of ones having a maximum number of copy in at least one organism above the families_repeted_th threshold attribute.

            .. attribute:: pan_size

                The number of nodes into the graph.
                .. warning:: this number is not necesserally equal to the number of imported gene families. Indeed, if imported families are not present into organism, superfluous will be discard.

            .. attribute:: is_partitionned
            
                a boolean specifying if the pangenome graph has been partitionned or not

            .. attribute:: nem_intermediate_files

                a str storing the path to the nem nem intermediate_files

            .. attribute:: partitions

                a dict providing the families present in each partion:
                    * partitions["exact_core"] contains the list of core families (exact core)
                    * partitions["exact_accessory"] contains the list of families not in the exact core 
                    * partitions["soft_core"] contains the list of soft core families (core exact)
                    * partitions["soft_accessory"] contains the list of families not in the soft core
                    * partitions["persistent"] contains the list of persistent families
                    * partitions["shell"] contains the list of shell families
                    * partitions["cloud"] contains the list of cloud families
                    * partitions["undefined"] contains the list of families unable to be classified (probably because the number of organisms is too small)

            .. attribute:: BIC

                a float providing the Bayesian Information Criterion. This Criterion give an estimation of the quality of the partionning (a low value means a good one)
                . seealso:: https://en.wikipedia.org/wiki/Bayesian_information_criterion

            .. attribute:: soft_core_th

                a float between 0 and 1 providing the threshold ratio of presence to attribute a gene families to the soft core genome
    """ 
    def __init__(self, init_from = "args", *args):
        """ 
            :param init_from: specified the excepted input (can be "file", "args", "database")
            :param *args: depending on the previous paramter, args can take multiple forms
            :type init_from: str
            :type *args: list

            :Example:

            >>>pan = PPanGGOLiN("file", organisms, gene_families, remove_high_copy_number_families)
            >>>pan = PPanGGOLiN("args", annotations, organisms, circular_contig_size, families_repeted, directed, distance_CDS_fragments)# load direclty the main attributes
        """ 
        self.directed                       = False
        self.annotations                    = dict()
        self.neighbors_graph                = None
        self.untangled_neighbors_graph      = None
        self.index                          = bidict()
        self.organisms                      = OrderedSet()
        self.nb_organisms                   = 0
        self.circular_contig_size           = dict()
        self.families_repeted_th            = 0
        self.families_repeted               = set()
        self.pan_size                       = 0
        self.is_partitionned                = False
        self.nem_intermediate_files         = None
        self.partitions                     = {}
        for p in SHORT_TO_LONG.values():
            self.partitions[p]              = list()
        self.BIC                            = None # Bayesian Index Criterion
        self.partitions_by_organism         = dict()
        self.subpartitions_shell_parameters = {}
        self.subpartition_shell             = {}
        self.distance_CDS_fragments         = 0
        self.CDS_fragments                  = {}
        self.soft_core_th                   = None
        if init_from == "file":
            self.__initialize_from_files(*args)
        elif init_from == "args":
            (self.annotations,
             self.organisms,
             self.circular_contig_size,
             self.families_repeted,
             self.directed,
             self.distance_CDS_fragments) = args 
        elif init_from == "database":
            logging.getLogger().error("database is not yet implemented")
            pass
        else:
            raise ValueError("init_from parameter is required")
        self.nb_organisms = len(self.organisms)

        logging.getLogger().info("Computing gene neighborhood ...")
        self.__neighborhood_computation()

    def __initialize_from_files(self, organisms_file, families_tsv_file, lim_occurence = 0, infer_singletons = False, distance_CDS_fragments = 0, directed = False):
        """ 
            :param organisms_file: a file listing organims by compute, first column is organism name, second is path to gff file and optionnally other other to provide the name of circular contig
            :param families_tsv_file: a file listing families. The first element is the family identifier (by convention, we advice to use the identifier of the average gene of the family) and then the next elements are the identifiers of the genes belonging to this family.
            :param lim_occurence: a int containing the threshold of the maximum number copy of each families. Families exceeding this threshold are removed and are listed in the families_repeted attribute.
            :param infer_singletons: a bool specifying if singleton must be explicitely present in the families_tsv_file (False) or if single gene in gff files must be automatically infered as a singleton family (True)
             :param distance_CDS_fragments: an integer specifying the distance between two consecutive genes belonging to the same gene families to consider them as CDS fragments (and so to not add the reflexive links)            param directed: a bool specifying if the pangenome graph is directed or undirected
            :type file: 
            :type file: 
            :type int: 
            :type bool: 
            :type int:
            :type bool: 
        """ 
        self.directed = directed
        self.distance_CDS_fragments = distance_CDS_fragments
        logging.getLogger().info("Reading "+families_tsv_file.name+" the gene families file ...")

        families_tsv_file = read_compressed_or_not(families_tsv_file)
        organisms_file    = read_compressed_or_not(organisms_file)
        families    = dict()
        first_iter  = True
        for line in families_tsv_file:
            elements = [el.strip() for el in line.split()]
            for gene in elements[1:]:
                families[gene]=elements[0]

        self.circular_contig_size = {}

        logging.getLogger().info("Reading "+organisms_file.name+" the list of organism files ...")

        bar = tqdm(organisms_file,total=get_num_lines(organisms_file), unit = "gff file")

        for line in bar:
            elements = [el.strip() for el in line.split("\t")]
            bar.set_description("Processing "+elements[ORGANISM_GFF_FILE])
            bar.refresh()
            if len(elements)>2:
                self.circular_contig_size.update({contig_id: None for contig_id in elements[2:len(elements)]})  # size of the circular contig is initialized to None (waiting to read the gff files to fill the dictionnaries with the correct values)
            self.annotations[elements[0]] = self.__load_gff(elements[ORGANISM_GFF_FILE], families, elements[ORGANISM_ID], lim_occurence, infer_singletons)
        check_circular_contigs = {contig: size for contig, size in self.circular_contig_size.items() if size == None }
        if len(check_circular_contigs) > 0:
            logging.getLogger().error("""
                The following identifiers of circular contigs in the file listing organisms have not been found in any region feature of the gff files: '"""+"'\t'".join(check_circular_contigs.keys())+"'")
            exit()

    def __load_gff(self, gff_file_path, families, organism, lim_occurence = 0, infer_singletons = False):
        """
            Load the content of a gff file
            :param gff_file_path: a valid gff file path where only feature of the type 'CDS' will be imported as genes. Each 'CDS' feature must have a uniq ID as attribute (afterall called gene id).
            :param families: a dictionary having the gene as key and the identifier of the associated family as value. Depending on the infer_singletons attribute, singleton must be explicetly present on the dictionnary or not
            :param organism: a str containing the organim name
            :param lim_occurence: a int containing the threshold of the maximum number copy of each families. Families exceeding this threshold are removed and are listed in the next attribute.
            :param infer_singletons: a bool specifying if singleton must be explicitely present in the families parameter (False) or if single gene automatically infered as a singleton family (True)
            :type str: 
            :type dict: 
            :type str: 
            :type int: 
            :type bool: 
            :return: annot: 
            :rtype: dict 
        """ 

        logging.getLogger().debug("Reading "+gff_file_path+" file ...")
        if organism not in self.organisms:
            self.organisms.add(organism)
            annot = defaultdict(OrderedDict)

            ctp_prev = 1
            cpt_fam_occ = defaultdict(int)

            gene_id_auto = False

            with read_compressed_or_not(gff_file_path) as gff_file:
                for line in gff_file:
                    if line.startswith('##',0,2):
                        if line.startswith('FASTA',2,7):
                            break
                        elif line.startswith('sequence-region',2,17):
                            fields = [el.strip() for el in line.split()]
                            if fields[1] in self.circular_contig_size:
                                self.circular_contig_size[fields[1]] = int(fields[3])
                            else:
                                logging.getLogger().debug(fields[1]+" is not circular")
                        continue
                    gff_fields = [el.strip() for el in line.split('\t')]
                    if GFF_feature == 'region':
                        if GFF_seqname in self.circular_contig_size:
                            self.circular_contig_size = int(GFF_end)
                            continue

                    elif gff_fields[GFF_feature] == 'CDS':
                        attributes_field = [f for f in gff_fields[GFF_attribute].strip().split(';') if len(f)>0]
                        attributes = {}
                        for att in attributes_field:
                            (key, value) = att.strip().split('=')
                            attributes[key.upper()]=value
                        try:
                            protein = attributes["ID"]
                        except:
                            logging.getLogger().error("Each CDS feature of the gff files must own a unique ID attribute. Not the case for file: "+gff_file_path)
                            exit(1)
                        try:
                            family = families[protein]
                        except KeyError:
                            if infer_singletons:
                                families[protein] = protein
                                family            = families[protein]
                                logging.getLogger().info("infered singleton: "+protein)
                            else:
                                raise KeyError("Unknown families:"+protein, ", check your families file or run again the program using the option to infer singleton")

                        cpt_fam_occ[family]+=1
                        prev = families[protein]

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

                        annot[gff_fields[GFF_seqname]][protein] = ["CDS",family,int(gff_fields[GFF_start]),int(gff_fields[GFF_end]),gff_fields[GFF_strand], name, product]

            for seq_id in list(annot):#sort genes by annotation start coordinate
                annot[seq_id] = OrderedDict(sorted(annot[seq_id].items(), key = lambda item: item[1][START]))
            if (lim_occurence > 0):
                fam_to_remove =[fam for fam, occ in cpt_fam_occ.items() if occ > lim_occurence]
                logging.getLogger().debug("highly repeted families found (>"+str(lim_occurence)+" in "+organism+"): "+" ".join(fam_to_remove))
                self.families_repeted = self.families_repeted.union(set(fam_to_remove))

            return(annot)
        else:
            raise KeyError("Redondant organism names was found ("+organism+")")

    def __str__(self):
        """ Return an overview of the statistics of the pangenome as a formated string """ 

        def str_histogram(title, values, force_max_value=25):
            ret = "\n".join([l for l in Pyasciigraph(force_max_value=force_max_value).graph(str(title), [("degree "+str(i),v) for i, v in enumerate(values[:force_max_value+1])])])
            if len(values) > force_max_value:
                ret+="\n"
                ret+="And "+str(sum(values[50:]))+" nodes having degree above "+str(force_max_value)+"..."
            return("↵"+ret)

        pan_str ="\n"
        pan_str += "----------- Statistics -----------\n"
        pan_str += "Number of organisms: "+str(self.nb_organisms)+"\n"
        pan_str += "Pangenome size:"+str(self.pan_size)+"\n"
        pan_str += "\n"
        if self.is_partitionned:
            pan_str += "Exact core genome size:"+str(len(self.partitions["exact_core"]))+"\n"
            pan_str += "Exact accessory genome size:"+str(len(self.partitions["exact_accessory"]))+"\n"
            pan_str += "\n"
            pan_str += "Soft core (>="+str(self.soft_core_th*100)+"%) genome size:"+str(len(self.partitions["soft_core"]))+"\n"
            pan_str += "Soft accessory (<"+str(self.soft_core_th*100)+"%) genome size:"+str(len(self.partitions["soft_accessory"]))+"\n"
            pan_str += "\n"
            pan_str += "Persistent genome size:"+str(len(self.partitions["persistent"]))+"\n"
            pan_str += "Shell genome size:"+str(len(self.partitions["shell"]))+"\n"
            pan_str += "Cloud genome cloud:"+str(len(self.partitions["cloud"]))+"\n"
            pan_str += "\n"
            pan_str += "Gene families with undefined partition:"+str(len(self.partitions["undefined"]))+"\n"
            # pan_str += "\n"
            # pan_str += str_histogram("Degree distribution of the core genome partition: ",nx.degree_histogram(self.neighbors_graph.subgraph(self.partitions["exact_core"])))+"\n"
            # pan_str += str_histogram("Degree distribution of the exact_accessory genome partition: ",nx.degree_histogram(self.neighbors_graph.subgraph(self.partitions["shell"])))+"\n"
            # pan_str += "\n"
            # pan_str += str_histogram("Degree distribution of the persistent partition: ",nx.degree_histogram(self.neighbors_graph.subgraph(self.partitions["persistent"])))+"\n"
            # pan_str += str_histogram("Degree distribution of the shell partition: ",nx.degree_histogram(self.neighbors_graph.subgraph(self.partitions["shell"])))+"\n"
            # pan_str += str_histogram("Degree distribution of the cloud partition: ",nx.degree_histogram(self.neighbors_graph.subgraph(self.partitions["cloud"])))+"\n"
            # pan_str += "\n"
        else:
            pan_str += "No partitioning have been performed on this PPanGGOLiN instance\n"
            pan_str += "Run the partition() method to obtain more detailled statistics...\n"
            pan_str += "\n"
        pan_str += "Number of edges in the pangenome graph: "+str(nx.number_of_edges(self.neighbors_graph))+"\n"
        pan_str += str_histogram("Degree distribution of the pangenome graph: ",nx.degree_histogram(self.neighbors_graph))+"\n"
        pan_str += "\n"
        pan_str += "----------------------------------"

        return(pan_str)

    def add_organism(self, new_orgs, new_annotations, new_circular_contig_size, new_families_repeted):
        self.annotations.update(new_annotations)
        self.index                     = bidict()
        self.organisms = self.organisms + new_orgs
        self.nb_organisms = len(self.organisms)
        self.circular_contig_size.update(new_circular_contig_size)
        self.families_repeted = self.families_repeted + new_families_repeted
        self.delete_nem_intermediate_files()
        self.partitions                    = {}
        for p in SHORT_TO_LONG.values():
            self.partitions[p] = list()
        self.BIC                           = None 
        self.__neighborhood_computation(update=new_orgs)

    def __repr__(self):
        return(self.__str__())

    def __len__(self):
        """ return the number of gene families into this pangenome """
        self.pan_size


    def __iadd__(self, another_pan):
        """ add a pangenome to this pangenome (reset the partionning) """

        self.annotations.update(another_pan.annotations)
        self.neighbors_graph               = None
        self.organisms                     = self.organisms.union(another_pan.organisms)
        self.nb_organisms                  = len(self.organisms)
        self.circular_contig_size          = self.circular_contig_size.update(another_pan.circular_contig_size)
        self.families_repeted              = self.families_repeted.union(another_pan.families_repeted)
        self.pan_size                      = 0
        self.is_partitionned               = False
        self.partitions                    = {}
        for p in SHORT_TO_LONG.values():
            self.partitions[p] = list()
        self.BIC                           = None
        
        self.delete_nem_intermediate_files()
        self.__neighborhood_computation()

        return(self)

    def __add_gene(self, fam_id, org, gene, name, length, product, graph_type = "neighbors_graph"):
        """
            Add gene to the pangenome graph
            :param fam_id: The family identifier
            :param org: The organism name
            :param gene : The gene identifier
            :param name: The biological name of the gene
            :param length: The number of nucleotide of the gene
            :param product: The name of the protein function
            :param length: The number of nucleotide of the gene
            :type str: 
            :type str:
            :type str: 
            :type str: 
            :type str: 
            :type str: 
            :type str: 
        """ 
        graph = self.neighbors_graph
        if graph_type == "untangled_neighbors_graph":
            graph = self.untangled_neighbors_graph

        graph.add_node(fam_id)

        try: 
            graph.node[fam_id]["nb_genes"]+=1
        except KeyError:
            graph.node[fam_id]["nb_genes"]=1
        try:
            graph.node[fam_id][org].add(gene)
        except KeyError:
            graph.node[fam_id][org] = set([gene])

        for attribute in ["name","length","product"]:
            try:
                graph.node[fam_id][attribute].add(locals()[attribute])
            except KeyError:
                graph.node[fam_id][attribute]=set([locals()[attribute]])

    def __add_link(self, fam_id, fam_id_nei, org, length, graph_type = "neighbors_graph"):
        """
            Add line between families of a the pangenome graph
            :param fam_id: The family identifier the first node (need to be have at least one gene belonging to this family already added to the graph via the method __add_gene())
            :param fam_id_nei: The family identifier the second node (need to be have at least one gene belonging to this family already added to the graph via the method __add_gene())
            :param org : The identifier of the organism supporting this link
            :param length : The distance in number of base between the genes adding this link
            :type str: 
            :type str:
            :type str: 
        """ 
        graph = self.neighbors_graph
        if graph_type == "untangled_neighbors_graph":
            graph = self.untangled_neighbors_graph

        if not self.neighbors_graph.has_edge(fam_id,fam_id_nei):
            graph.add_edge(fam_id, fam_id_nei)
            # logging.getLogger().debug([str(i) for i in [fam_id, fam_id_nei, org]])
        try:
            graph[fam_id][fam_id_nei][org]+=1
        except KeyError:
            graph[fam_id][fam_id_nei][org]=1
            try:
                graph[fam_id][fam_id_nei]["weight"]+=1.0
            except KeyError:
                graph[fam_id][fam_id_nei]["weight"]=1.0
        try:
            graph[fam_id][fam_id_nei]["length"].add(length)
        except KeyError:
            graph[fam_id][fam_id_nei]["length"]=set([length])

    def get_gene_info(self, gene):
        """ return annotation info about a gene
            :param gene: a gene identifier
            :type str: 
            :return: gene_info: a dict of info about the gene
            :rtype: dict 
        """ 
        gene_index = self.index[gene]
        info = self.annotations[gene_index[ORGANISM_INDEX]][gene_index[CONTIG_INDEX]][gene]
        return(dict(zip(("ORGANISM","CONTIG","POSITION","TYPE", "FAMILY", "START", "END", "STRAND", "NAME", "PRODUCT"),list(gene_index)+info)))

    def __neighborhood_computation(self, update=False):#,light = False, 
        """ Use the information already loaded (annotation) to build the pangenome graph
            :param update: an optional list of organism to update a previous graph
            :type int: 
            :type bool: 
            :type list:
        """ 
        #:param light: a bool specifying is the annotation attribute must be detroyed at each step to save memory
        if self.neighbors_graph is None:
            if self.directed:
                self.neighbors_graph = nx.DiGraph()
            else:
                self.neighbors_graph = nx.Graph()
        if update:
            orgs = update
        else:
            orgs = tqdm(list(self.annotations),total=len(self.annotations),unit = "organism")
        for organism in orgs:
            if not update:
                orgs.set_description("Processing "+organism)
                orgs.refresh()
            for contig, contig_annot in self.annotations[organism].items():
                try:
                    (gene_start, gene_info_start) = contig_annot.popitem(last=False)
                    while (gene_info_start[FAMILY] in self.families_repeted):
                            (gene_start, gene_info_start) = contig_annot.popitem(last=False)
                except KeyError:
                    continue

                self.__add_gene(gene_info_start[FAMILY],
                                organism,
                                gene_start,
                                gene_info_start[NAME],
                                gene_info_start[END]-gene_info_start[START],
                                gene_info_start[PRODUCT])
                self.index[gene_start]=(organism,contig,0)

                gene_nei, gene_info_nei = gene_start, gene_info_start
                logging.getLogger().debug(gene_info_start)
                for pos, (gene, gene_info) in enumerate(contig_annot.items()):
                    logging.getLogger().debug(gene_info)
                    logging.getLogger().debug(gene)
                    if gene_info[FAMILY] not in self.families_repeted:
                        self.__add_gene(gene_info[FAMILY],
                                        organism,
                                        gene,
                                        gene_info[NAME],
                                        gene_info[END]-gene_info[START],
                                        gene_info[PRODUCT])
                        self.index[gene]=(organism,contig,pos+1)
                        self.neighbors_graph.add_node(gene_info_nei[FAMILY])
                        if not (gene_info[FAMILY] == gene_info_nei[FAMILY] and
                                gene_info[STRAND] == gene_info_nei[STRAND] and
                                (gene_info[START] - gene_info_nei[END]) <= self.distance_CDS_fragments):# to avoid reflexive links with gene fragments
                            self.__add_link(gene_info[FAMILY],gene_info_nei[FAMILY],organism, gene_info[START] - gene_info_nei[END])
                        else:
                            self.CDS_fragments[gene]=gene_info[FAMILY]
                            self.CDS_fragments[gene_nei]=gene_info[FAMILY]#or gene_info_nei[FAMILY]
                        gene_nei, gene_info_nei = gene, gene_info
                
                if contig in self.circular_contig_size:#circularization
                    if not (gene_info_start[FAMILY] == gene_info_nei[FAMILY] and \
                            gene_info_start[STRAND] == gene_info_nei[STRAND] and \
                            (self.circular_contig_size[contig] - gene_info_nei[END] + gene_info_start[START]) <= self.distance_CDS_fragments):# to avoid reflexive links with gene fragments
                        self.__add_link(gene_info_start[FAMILY],gene_info_nei[FAMILY],organism, (self.circular_contig_size[contig] - gene_info_nei[END]) + gene_info_start[START])
                    else:
                        self.CDS_fragments[gene]=gene_info[FAMILY]
                        self.CDS_fragments[gene_nei]=gene_info[FAMILY]
                if sys.version_info < (3,):
                    ordered_dict_prepend(contig_annot,gene_start,gene_info_start)#insert at the top
                else:
                    contig_annot[gene_start]=gene_info_start
                    contig_annot.move_to_end(gene_start, last=False)#move to the beginning
            # if light:
            #     del self.annotations[organism]
        self.pan_size = nx.number_of_nodes(self.neighbors_graph)

    def untangle_neighbors_graph(self, K = 3):
        
        self.untangled_neighbors_graph = self.neighbors_graph.copy()
        
        separation_tree = defaultdict(set)

        def absolute_orientation(a_list):
                orientation1=tuple(a_list)
                orientation2=tuple(reversed(a_list))
                if (hash(orientation1)>hash(orientation2)):
                    return(orientation1)
                else:
                    return(orientation2)

        def update_seed_path(seed_path):
            new_seed_paths = set()
            try:
                for i, element in enumerate(seed_path):
                    
                        if element in separation_tree:
                            print(element)
                            for child in separation_tree[element]:
                                new_seed_path = list(seed_path)
                                new_seed_path[i]=child
                                new_seed_paths.add(tuple(new_seed_path))
            except TypeError as e:
                pdb.set_trace()
            return(new_seed_paths)

        def extends_seeds(all_path_k):
            all_path_k_p_1 = set()
            all_path_k = deque(all_path_k)
            while all_path_k:
                p = all_path_k.pop()
                if p is not None:
                    try:
                        for nei in self.untangled_neighbors_graph.neighbors(p[0]):
                            path = (nei,)+p
                            
                            #absolute_orientation(path)# give an absolute orientation
                            all_path_k_p_1.add(absolute_orientation(path))
                    except nx.exception.NetworkXError as e:
                        print(path)
                        all_path_k.extendleft(update_seed_path(p))
                        continue
            return(all_path_k_p_1)

        def merge_overlapping_extremities(data): #inspired from http://stackoverflow.com/a/9114443/7500030s
            sets = (set(e) for e in data if e)
            try:
                results = [next(sets)]
                for e_set in sets:
                    to_update = []
                    for i,res in enumerate(results):
                        if not e_set.isdisjoint(res):
                            to_update.insert(0,i)
                    if not to_update:
                        results.append(e_set)
                    else:
                        last = results[to_update.pop(-1)]
                        for i in to_update:
                            last |= results[i]
                            del results[i]
                        last |= e_set
                return results
            except StopIteration:
                return {}

        # def filter_seed_paths(validated_seed_paths):
        #     print("validated_seed_paths: "+str(len(validated_seed_paths)))
        #     filtered_validated_seed_paths = set()
        #     for validated_seed_path in validated_seed_paths:
        #         for fam in validated_seed_path:
        #             if fam in self.untangled_neighbors_graph:
        #                 filtered_validated_seed_paths.add(validated_seed_path)

        #     print("filtered_validated_seed_paths: "+str(len(filtered_validated_seed_paths)))
        #     return(filtered_validated_seed_paths)

        def align_on_graph(path):
            extremities_seed_path = defaultdict(lambda: defaultdict(set))

            path_complete=defaultdict(lambda : False)

            #print(path)
            #pdb.set_trace()
            for org in [o for o in self.untangled_neighbors_graph.nodes[path[0]] if o not in RESERVED_WORDS]:
                for gene in self.untangled_neighbors_graph.nodes[path[0]][org]:
                    logging.getLogger().debug("here"+gene)
                    (pos,contig)   = (self.index[gene][POSITION_INDEX],self.index[gene][CONTIG_INDEX])
                    orientation    = None
                    path_exist     = True
                    tmp_path       = list()                                
                    tmp_path.append(self.index.inv[(org,contig,pos)])
                    circular = sys.maxsize if contig not in self.circular_contig_size else len(self.annotations[org][contig])
                    for i, fam in enumerate(path[1:]):
                        try:
                            if orientation is None:
                                if fam == self.annotations[org][contig][self.index.inv[(org,contig,(pos+(i+1))%circular)]][FAMILY]:
                                    orientation = 1
                                elif fam == self.annotations[org][contig][self.index.inv[(org,contig,(pos-(i+1))%circular)]][FAMILY]:
                                    orientation = -1
                                else:
                                    path_exist = False
                                    break
                            gene_i = self.index.inv[(org,contig,(pos+((i+1)*orientation))%circular)]
                            
                            #logging.getLogger().debug("fam="+fam+"  gene_i="+gene_i+"    self.annotations[org][contig][gene_i][FAMILY]="+self.annotations[org][contig][gene_i][FAMILY])
                            if fam == self.annotations[org][contig][gene_i][FAMILY]:
                                tmp_path.append(gene_i)
                            else:
                                path_exist = False
                                break
                        except KeyError:
                            break
                    if path_exist:
                        if len(tmp_path) == len(path):
                            path_complete[frozenset([path[0],path[len(path)-1]])] = True
                            extremities_seed_path[frozenset([path[0],path[len(path)-1]])][org].add(absolute_orientation(tuple(tmp_path)))
                        else:
                            path_complete[frozenset([path[0],path[len(path)-1]])] = path_complete[frozenset([path[0],path[len(path)-1]])]
                            extremities_seed_path[frozenset([path[0],path[len(path)-1]])][org].add(tuple(tmp_path))
                        logging.getLogger().debug(tmp_path)
            for org in [o for o in self.untangled_neighbors_graph.nodes[path[len(path)-1]] if o not in RESERVED_WORDS]:
                # if org == "org1" and path == ('fam1', 'fam2_5', 'fam2_5_bis', 'fam3',):
                #     pdb.set_trace()
                for gene in self.untangled_neighbors_graph.nodes[path[len(path)-1]][org]:
                    logging.getLogger().debug("here"+gene)
                    (pos,contig)   = (self.index[gene][POSITION_INDEX],self.index[gene][CONTIG_INDEX])
                    orientation    = None
                    path_exist     = True
                    tmp_path       = list()
                    tmp_path.append(self.index.inv[(org,contig,pos)])
                    circular = sys.maxsize if contig not in self.circular_contig_size else len(self.annotations[org][contig])
                    for i, fam in enumerate(reversed(path[:len(path)-1])):
                        try:
                            if orientation is None:
                                if fam == self.annotations[org][contig][self.index.inv[(org,contig,(pos+(i+1))%circular)]][FAMILY]:
                                    orientation = 1
                                elif fam == self.annotations[org][contig][self.index.inv[(org,contig,(pos-(i+1))%circular)]][FAMILY]:
                                    orientation = -1
                                else:
                                    path_exist = False
                                    break
                            gene_i = self.index.inv[(org,contig,(pos+((i+1)*orientation))%circular)]
                            logging.getLogger().debug("fam="+fam+"  gene_i="+gene_i+"    self.annotations[org][contig][gene_i][FAMILY]="+self.annotations[org][contig][gene_i][FAMILY])
                            if fam == self.annotations[org][contig][gene_i][FAMILY]:
                                tmp_path.append(gene_i)
                            else:
                                path_exist = False
                                break
                        except KeyError:
                            break
                    if path_exist:
                        if len(tmp_path) == len(path):
                            path_complete[frozenset([path[0],path[len(path)-1]])] = True
                            extremities_seed_path[frozenset([path[0],path[len(path)-1]])][org].add(absolute_orientation(tuple(tmp_path)))
                        else:
                            path_complete[frozenset([path[0],path[len(path)-1]])] = path_complete[frozenset([path[0],path[len(path)-1]])]
                            extremities_seed_path[frozenset([path[0],path[len(path)-1]])][org].add(tuple(tmp_path))
                        
                        logging.getLogger().debug(tmp_path)
            #print(path_complete)
            for ext,complete in path_complete.items():
                if not complete:
                    del extremities_seed_path[ext]
            return(extremities_seed_path)

        validated_seed_paths = set()
        
        for k in range(1,K+1):
            logging.getLogger().debug("k="+str(k))
            if k == 1:
                LILO_seed_paths = deque([(p,) for p in self.untangled_neighbors_graph.nodes()])
            else:
                logging.getLogger().debug(validated_seed_paths)
                LILO_seed_paths = deque(extends_seeds(validated_seed_paths))
                validated_seed_paths = set()
            #pdb.set_trace()
            with tqdm(total=len(LILO_seed_paths)) as pbar:
                while LILO_seed_paths:
                    seed_path = LILO_seed_paths.popleft()
                    pbar.update(1)
                    
                    neighbors=set()
                    all_extremities_seed_path = defaultdict(lambda: defaultdict(set))
                    extremity_groups=None
                    #print(seed_path)
                    if seed_path is not None:
                        new_seed_paths = update_seed_path(seed_path)
                        if len(new_seed_paths)>0:
                            LILO_seed_paths.extendleft(new_seed_paths)
                            continue
                    else:
                            continue
                    #print(separation_tree)
                    tested_extremities=set()
                    neighbors = set()
                    set_seed_path = set(seed_path)
                    try :
                        #if nx.degree(self.untangled_neighbors_graph,seed_path[0])>2 and nx.degree(self.untangled_neighbors_graph,seed_path[k-1])>2:
                        for nei_left in self.untangled_neighbors_graph.neighbors(seed_path[0]):
                            if nei_left in set_seed_path and nei_left != seed_path[0] and nei_left != seed_path[k-1]:
                                continue
                            for nei_right in self.untangled_neighbors_graph.neighbors(seed_path[k-1]):
                                if nei_right in set_seed_path and nei_left != seed_path[0] and nei_left != seed_path[k-1]:
                                    continue
                                extremities = frozenset([nei_right,nei_left])
                                if len(extremities)>1 and extremities not in tested_extremities:
                                    flanked_seed_path = (nei_left,)+seed_path+(nei_right,)
                                    logging.getLogger().debug(flanked_seed_path)
                                    res = align_on_graph(flanked_seed_path)
                                    all_extremities_seed_path.update(res)
                                    # pdb.set_trace()
                                    if len(res)>0:
                                        # validated_neighbors.add(nei_left)
                                        # validated_neighbors.add(nei_rigth)
                                        validated_seed_paths.add(seed_path)# way to extends seeds in a optimized way
                                        logging.getLogger().debug("validated: "+str(seed_path))
                                    tested_extremities.add(extremities)
                    except:
                        logging.getLogger().debug(seed_path[0]+"not in graph")
                        pass


                    #not_validated_neighbors = neighbors - validated_seed_paths
                    #logging.getLogger().debug(all_extremities_seed_path)
                    extremity_groups = merge_overlapping_extremities(all_extremities_seed_path)
                    #print(extremity_groups)
                    #pdb.set_trace()
                    # if inc>10000:
                    #     pdb.set_trace()
                    if len(extremity_groups) > 1:
                        #separation
                        for nb, extremitiy_group in enumerate(extremity_groups):
                            suffix = "$k"+str(k)+"_"+str(nb)
                            new_seed_path = [None]*(k+2)
                            for extremity,org_genes in all_extremities_seed_path.items():
                                if extremity & extremitiy_group:
                                    gene_info_prec = None
                                    for org, gene_series in org_genes.items():
                                        # if org == "org1":
                                        #     pdb.set_trace()
                                        for genes in gene_series:
                                            gene_info_prec = self.annotations[org][self.index[genes[0]][CONTIG_INDEX]][genes[0]]
                                            new_seed_path[0]= gene_info_prec[FAMILY]
                                            for i, gene in enumerate(genes[1:]):
                                                gene_info = self.annotations[org][self.index[gene][CONTIG_INDEX]][gene]
                                                logging.getLogger().debug(gene_info[FAMILY])
                                                new_family_name = gene_info[FAMILY]+suffix
                                                new_seed_path[i+1]=new_family_name
                                                if i < k:                                                   
                                                    if gene_info[FAMILY] in self.untangled_neighbors_graph:
                                                        #some things to do                                                        
                                                        self.untangled_neighbors_graph.remove_node(gene_info[FAMILY])
                                                    separation_tree[gene_info[FAMILY]].add(new_family_name)
                                                    
                                                    self.__add_gene(new_family_name,org,gene,gene_info[NAME],gene_info[END]-gene_info[START],gene_info[PRODUCT],"untangled_neighbors_graph")
                                                    self.annotations[org][self.index[gene][CONTIG_INDEX]][gene][FAMILY]=new_family_name
                                                    if self.is_partitionned:
                                                        self.untangled_neighbors_graph.nodes[new_family_name]["partition"]="undefined"

                                                    #self.untangled_neighbors_graph.nodes[new_family_name]["untangled"]=k
                                                # circular cases
                                                
                                                length = (gene_info[START] - gene_info_prec[END]) if (gene_info_prec[START] < gene_info[START]) else (gene_info_prec[START] - gene_info[END])
                                                self.__add_link(self.annotations[org][self.index[gene][CONTIG_INDEX]][gene][FAMILY], gene_info_prec[FAMILY], org, length,"untangled_neighbors_graph")
                                                
                                                gene_info_prec = gene_info

                                            LILO_seed_paths.extend(tuple(new_seed_path[1:len(new_seed_path)-1]))
                                            #pdb.set_trace()
                                            validated_seed_paths.add(tuple(new_seed_path[1:len(new_seed_path)-1]))
                                            validated_seed_paths.add(tuple(new_seed_path[0:len(new_seed_path)-2]))
                                            validated_seed_paths.add(tuple(new_seed_path[2:len(new_seed_path)]))
                        
                    else:
                        pass#validated_seed_paths
                                                #add new seed to LILO queue 
                                                # think to non_valided_neibors
                                                #refine validated_seed_paths
                    all_extremities_seed_path = None

    def get_gene_families_related_to_metadata(self, metadata, output_path = None):
        exclusively_in = dict()
        never_in       = dict()
        variables_and_possible_values = defaultdict(set)
        for org, metaD in metadata.items():
            for variable, value in metaD.items():
                exclusively_in[variable] = defaultdict(set)
                never_in[variable]       = defaultdict(set)
                variables_and_possible_values[variable].add(value)
        for node_name, node_organisms in self.neighbors_graph.nodes(data=True):
            for variable, possible_values in variables_and_possible_values.items():
                if output_path is not None:
                    try:
                        os.makedirs(output_path+"/exclusively_in/"+variable)
                        os.makedirs(output_path+"/never_in/"+variable)
                    except:
                        pass
                res = set([metadata[org][variable] for org in metadata if org in node_organisms])
                opposite_res = (possible_values - res)

                if len(res)==1:
                    value = res.pop()
                    if output_path is not None:
                        with open(output_path+"/exclusively_in/"+variable+"/"+str(value),"a") as out_file:
                            out_file.write(node_name+"\n")
                    exclusively_in[variable][value].add(node_name)
                elif len(opposite_res) == 1:
                    value = opposite_res.pop()
                    if output_path is not None:
                        with open(output_path+"/never_in/"+variable+"/"+str(value),"a") as out_file:
                            out_file.write(node_name+"\n")
                    never_in[variable][value].add(node_name)
        return ({"exclusively_in":exclusively_in,"never_in":never_in})

    def __write_nem_input_files(self, nem_dir_path, select_organisms, init = "default", low_disp=0.1, filter_by_partition = None):
        if len(select_organisms)<=10:# below 10 organisms a statistical computation do not make any sence
            logging.getLogger().warning("The number of selected organisms is too low ("+str(len(select_organisms))+" organisms used) to partition the pangenome graph in persistent, shell and cloud genome. Add new organisms to obtain more robust metrics.")

        if not os.path.exists(nem_dir_path):
            #NEM requires 5 files: nem_file.index, nem_file.str, nem_file.dat, nem_file.m and nem_file.nei
            os.makedirs(nem_dir_path)

        logging.getLogger().debug("Writing nem_file.str nem_file.index nem_file.nei nem_file.dat and nem_file.m files")
        with open(nem_dir_path+"/nem_file.str", "w") as str_file,\
             open(nem_dir_path+"/nem_file.index", "w") as index_file,\
             open(nem_dir_path+"/column_org_file", "w") as org_file,\
             open(nem_dir_path+"/nem_file.nei", "w") as nei_file,\
             open(nem_dir_path+"/nem_file.dat", "w") as dat_file,\
             open(nem_dir_path+"/nem_file.m", "w") as m_file:

            nei_file.write("1\n")
            
            org_file.write(" ".join(["\""+org+"\"" for org in select_organisms])+"\n")
            org_file.close()

            index_fam = OrderedDict()
            for node_name, node_organisms in self.neighbors_graph.nodes(data=True):
                if filter_by_partition is not None and "partition" in node_organisms and node_organisms["partition"] != filter_by_partition:
                    continue
                logging.getLogger().debug(node_organisms)
                logging.getLogger().debug(select_organisms)
                
                if not select_organisms.isdisjoint(node_organisms): # if at least one commun organism
                    dat_file.write("\t".join(["1" if org in node_organisms else "0" for org in select_organisms])+"\n")
                    index_fam[node_name] = len(index_fam)+1
                    index_file.write(str(len(index_fam))+"\t"+str(node_name)+"\n")
            for node_name, index in index_fam.items():
                row_fam         = []
                row_dist_score  = []
                neighbor_number = 0
                try:
                    for neighbor in set(nx.all_neighbors(self.neighbors_graph, node_name)):
                        if filter_by_partition is not None and "partition" in self.neighbors_graph.node[neighbor] and self.neighbors_graph.node[neighbor]["partition"] == filter_by_partition:
                            continue
                        coverage = 0
                        if self.neighbors_graph.is_directed():
                            cov_sens, cov_antisens = (0,0)
                            try:
                                cov_sens = sum([pre_abs for org, pre_abs in self.neighbors_graph[node_name][neighbor].items() if ((org in select_organisms) and (org not in RESERVED_WORDS))])
                            except KeyError:
                                pass
                            try:
                                cov_antisens = sum([pre_abs for org, pre_abs in self.neighbors_graph[neighbor][node_name].items() if ((org in select_organisms) and (org not in RESERVED_WORDS))])
                            except KeyError:
                                pass
                            coverage = cov_sens + cov_antisens
                        else:
                            coverage = sum([pre_abs for org, pre_abs in self.neighbors_graph[node_name][neighbor].items() if ((org in select_organisms) and (org not in RESERVED_WORDS))])

                        if coverage==0:
                            continue
                        distance_score = coverage#/len(select_organisms)
                        row_fam.append(str(index_fam[neighbor]))
                        row_dist_score.append(str(round(distance_score,4)))
                        neighbor_number += 1
                    if neighbor_number>0:
                        nei_file.write("\t".join([str(item) for sublist in [[index_fam[node_name]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
                    else:
                        nei_file.write(str(index_fam[node_name])+"\t0\n")
                        logging.getLogger().debug("The family: "+node_name+" is an isolated family in the selected organisms")
                except nx.exception.NetworkXError as nxe:
                    print(nxe)
                    logging.getLogger().debug("The family: "+node_name+" is an isolated family")
                    nei_file.write(str(index_fam[node_name])+"\t0\n")

            if init is not None:
                m_file.write("1 ")# 1 to initialize parameter,
                if init == "default":
                    m_file.write("0.33333 0.33333 ")# 0.333 and 0.333 for to give one third of initial proportition to each class (last 0.33 is automaticaly determined by substraction)
                    m_file.write(" ".join(["1"]*len(select_organisms))+" ") # persistent binary vector
                    m_file.write(" ".join(["0.5"]*len(select_organisms))+" ") # shell binary vector (1 ou 0, whatever because dispersion will be of 0.5)
                    m_file.write(" ".join(["0"]*len(select_organisms))+" ") # cloud binary vector
                    m_file.write(" ".join([str(low_disp)]*len(select_organisms))+" ") # persistent dispersition vector (low)
                    m_file.write(" ".join(["0.5"]*len(select_organisms))+" ") # shell dispersition vector (high)
                    m_file.write(" ".join([str(low_disp)]*len(select_organisms))) # cloud dispersition vector (low)
                elif isinstance(init,dict):
                    all_orgs_in_groups = set([org for orgs in init.values() for org in orgs])
                    (a,b,c)=(0,0,0)
                    for p in range(0,len(init)):
                        a+=1
                        m_file.write(str(round(float(1)/len(init),4))+" ")

                    for org_groups in init.values():
                        b+=1
                        m_file.write(" ".join(["1" if org in org_groups else "0" if org in all_orgs_in_groups else "0.5" for org in select_organisms])+" ")
                    m_file.write(" ".join(["0.5"]*len(select_organisms))+" ")
                    for org_groups in init.values():
                        c+=1
                        m_file.write(" ".join([str(low_disp) if org in org_groups else str(low_disp) if org in all_orgs_in_groups else "0.5" for org in select_organisms])+" ")
                    m_file.write(" ".join(["0.5"]*len(select_organisms)))

                    print("a="+str(a)+"    b="+str(b)+"      c"+str(c))
                elif isinstance(init,list):
                    m_file.write("0.33333 0.33333 ")
                    (positive,negative) = init
                    m_file.write(" ".join(["1" if org in positive else "0" if org in negative else "0.5" for org in select_organisms])+" ")
                    m_file.write(" ".join(["0" if org in positive else "1" if org in negative else "0.5" for org in select_organisms])+" ")
                    m_file.write(" ".join(["0.5"]*len(select_organisms))+" ")
                    m_file.write(" ".join([str(low_disp) if org in positive else str(low_disp) if org in negative else "0.5" for org in select_organisms])+" ")
                    m_file.write(" ".join([str(low_disp) if org in positive else str(low_disp) if org in negative else "0.5" for org in select_organisms])+" ")
                    m_file.write(" ".join(["0.5"]*len(select_organisms)))

            str_file.write("S\t"+str(len(index_fam))+"\t"+
                                 str(len(select_organisms))+"\n")

    def partition(self, nem_dir_path     = tempfile.mkdtemp(),
                        select_organisms = None,
                        beta             = 0.5,
                        free_dispersion  = False,
                        chunck_size      = 500,
                        soft_core_th     = 0.95,
                        inplace          = True,
                        just_stats       = False,
                        nb_threads       = 1):
        """
            Use the graph topology and the presence or absence of genes from each organism into families to partition the pangenome in three groups ('persistent', 'shell' and 'cloud')
            . seealso:: Read the Mo Dang's thesis to understand NEM, a summary is available here : http://www.kybernetika.cz/content/1998/4/393/paper.pdf
            :param nem_dir_path: a str containing a path to store temporary file of the NEM program
            :param select_organisms: a list of organism to used to obtain the partition (must be included in the organism attributes of the object) or None to used all organisms in the object
            :param beta: a float containing the spatial coefficient of smoothing of the clustering results using the weighted graph topology (0.00 turn off the spatial clustering)
            :param free_dispersion: a bool specyfing if the dispersion around the centroid vector of each paritition is the same for all the organisms or if the dispersion is free
            :param chunck_size: an int specifying the size of the chunks
            :param soft_core_th a float between 0 and 1 providing the threshold ratio of presence to attribute a gene families to the soft core genome
            :param inplace: a boolean specifying if the partition must be stored in the object of returned (throw an error if inplace is true and organisms parameter i not None)
            :param just_stats: a boolean specifying if the partitions must be returned or just stats about them (number of families in each partition)
            :param nb_threads: an integer specifying the number of threads to use (works only if the number of organisms is higher than the chunck_size)
            :type str: 
            :type list: 
            :type float: 
            :type bool:
            :type float
            :type int: 
            :type bool: 
            :type bool: 
            :type int: 
        """ 
        
        if select_organisms is None:
            select_organisms = self.organisms
        else:
            select_organisms = OrderedSet(select_organisms)
            if len(select_organisms - self.organisms)>0:
                raise Exception("select_organisms parameter must be included in the organisms attribute of the objet")
            if inplace:
                raise Exception("inplace can't be true if the select_organisms parameter has not the same size than organisms attribute")

        if soft_core_th > 1 or soft_core_th < 0:
            raise Exception("soft_core_th parameter must be a float between 0 and 1")
        if inplace:
            self.soft_core_th = soft_core_th
        # if self.neighbors_graph is None:
        #     raise Exception("The neighbors_graph is not built, please use the function neighborhood_computation before")
        if self.is_partitionned and inplace:
            logging.getLogger().warning("The pangenome was already partionned, inplace=true parameter will erase previous nem file intermediate files, partitions and subpartitions")
            self.delete_nem_intermediate_files()

        if inplace:
            self.nem_intermediate_files = nem_dir_path

        stats = defaultdict(int)
        partitionning_results = {}
        
        #core first
        families = []
        for node_name, data_organisms in self.neighbors_graph.nodes(data=True):
            families.append(node_name)
            pres_abs_vector = [True if org in data_organisms else False for org in select_organisms]
            
            nb_true  = pres_abs_vector.count(True)
            if nb_true == 0:
                continue # family absent of the selected_organisms

            if nb_true==len(select_organisms):
                stats["exact_core"]+=1
            else:
                stats["exact_accessory"]+=1

            if nb_true>=len(select_organisms)*soft_core_th:
                stats["soft_core"]+=1
            else:
                stats["soft_accessory"]+=1
        BIC = 0
        
        if len(select_organisms) > chunck_size:

            cpt_partition = OrderedDict()
            for fam in families:
                cpt_partition[fam]= {"P":0,"S":0,"C":0,"U":0}
            
            total_BIC = 0

            @contextlib.contextmanager
            def empty_cm():
                yield None

            validated = set()
            cpt=0

            if inplace:
                bar = tqdm(total = stats["exact_accessory"]+stats["exact_core"], unit = "families partitionned")

            sem = Semaphore(nb_threads)

            def validate_family(result):
                #nonlocal total_BIC
                try:
                    partitions = result
                    #total_BIC += BIC
                    for node,nem_class in partitions[FAMILIES_PARTITION].items():
                        cpt_partition[node][nem_class]+=1
                        sum_partionning = sum(cpt_partition[node].values())
                        if (sum_partionning > len(select_organisms)/chunck_size and max(cpt_partition[node].values()) >= sum_partionning*0.5) or (sum_partionning > len(select_organisms)):
                            if node not in validated:
                                if inplace:
                                    bar.update()
                                if max(cpt_partition[node].values()) < sum_partionning*0.5:
                                    cpt_partition[node]["U"] = sys.maxsize #if despite len(select_organisms) partionning, the abosolute majority is found, then the families is set to undefined 
                                validated.add(node)
                                # if max(cpt_partition[node], key=cpt_partition[node].get) == "P" and cpt_partition[node]["S"]==0 and cpt_partition[node]["C"]==0:
                                        #     validated[node]="P"
                                        # elif cpt_partition[node]["S"]==0:
                                        #     validated[node]="C"
                                        # else:
                                        #     validated[node]="S" 
                finally:
                    sem.release()

            with contextlib.closing(Pool(processes = nb_threads)) if nb_threads>1 else empty_cm() as pool:
            
                #proba_sample = OrderedDict(zip(select_organisms,[len(select_organisms)]*len(select_organisms)))

                pan_size = stats["exact_accessory"]+stats["exact_core"]
                
                while len(validated)<pan_size:
                    if sem.acquire() if nb_threads>1 else True:#
                        # print(select_organisms)
                        # print(chunck_size)
                        # print(proba_sample.values())
                        # print(len(proba_sample.values()))
                        # min_o = min(proba_sample.values()) 
                        # max_o = max(proba_sample.values()) 
                        # range_o = max_o-min_o
                        # if min_o != max_o:
                        #     p = [(p-min_o/range_o)/len(select_organisms) for p in proba_sample.values()]
                        # else:
                        #     p = list(proba_sample.values())
                        # print(p)
                        #s = sum(proba_sample.values())
                        
                        #orgs = np.random.choice(select_organisms, size = chunck_size, replace = False, p = [p/s for p in proba_sample.values()])#
                        orgs = sample(select_organisms,chunck_size)
                        orgs = OrderedSet(orgs)

                        # for org, p in proba_sample.items():
                        #     if org in orgs:
                        #         proba_sample[org] = p - len(select_organisms)/chunck_size if p >1 else 1
                        #     else:
                        #         proba_sample[org] = p + len(select_organisms)/chunck_size

                        index = self.__write_nem_input_files(nem_dir_path+"/"+str(cpt)+"/",
                                                             orgs)
                        if nb_threads>1:
                            res = pool.apply_async(run_partitioning,
                                                   args = (nem_dir_path+"/"+str(cpt)+"/",#nem_dir_path
                                                           len(orgs),
                                                           beta,
                                                           free_dispersion),                                                        
                                                   callback = validate_family)
                        else:
                            res = run_partitioning(nem_dir_path+"/"+str(cpt)+"/",#nem_dir_path
                                                   len(orgs),
                                                   beta,
                                                   free_dispersion)
                            validate_family(res)
                        cpt +=1
                    else:
                        sleep(0.01)

                    # if inplace:
                    #     bar.update()    
                if nb_threads>1:
                    sleep(1)
                    pool.close()
                    pool.join() 
                #BIC = total_BIC/cpt
                BIC = 0
            # if just_stats:
            #     print('len(validated)= '+str(len(validated)))
            #     print('len(cpt_partition)= '+str(len(cpt_partition)))

            for fam, data in cpt_partition.items():
                partitionning_results[fam]=max(data, key=data.get)

            # if just_stats:
            #     print("stat")
            #     #print(partitions)
            #     c = Counter(partitions)      
            #     print('stats["P"] '+str(c["P"]))
            #     print('stats["S"] '+str(c["S"]))
            #     print('stats["C"] '+str(c["C"])) 
            #     print('stats["U"] '+str(c["U"]))           
            #     print('total '+str(c["P"]+c["S"]+c["C"]+c["U"]))


            #     print('stats["exact_accessory"] '+str(stats["exact_accessory"]))
            #     print('stats["exact_core"] '+str(stats["exact_core"]))
            #     print('total '+str(stats["exact_accessory"]+stats["exact_core"]))
            #     print(' ')
        else:
            self.__write_nem_input_files(nem_dir_path+"/",
                                         select_organisms)
            partitionning_results = run_partitioning(nem_dir_path, len(select_organisms), beta, free_dispersion)[FAMILIES_PARTITION]
            
        if inplace:
            self.BIC = BIC
            if self.is_partitionned:
                for p in SHORT_TO_LONG.values():
                    self.partitions[p] = list()# erase older values
            for node, nem_class in partitionning_results.items():
                nb_orgs=0
                for key in list(self.neighbors_graph.node[node].keys()):
                    if key not in RESERVED_WORDS:
                        #self.partitions_by_organisms[key][partition[int(nem_class)]].add(self.neighbors_graph.node[node][key])
                        nb_orgs+=1

                self.neighbors_graph.node[node]["partition"]=SHORT_TO_LONG[nem_class]
                
                self.partitions[SHORT_TO_LONG[nem_class]].append(node)

                if nb_orgs == self.nb_organisms:
                    self.partitions["exact_core"].append(node)#EXACT CORE
                    self.neighbors_graph.node[node]["partition_exact"]="exact_core"
                elif nb_orgs < self.nb_organisms:
                    self.partitions["exact_accessory"].append(node)#EXACT ACCESSORY
                    self.neighbors_graph.node[node]["partition_exact"]="exact_accessory"
                else:
                    logging.getLogger().error("nb_orgs can't be > to self.nb_organisms")
                    exit(1)

                if nb_orgs >= (self.nb_organisms)*self.soft_core_th:
                    self.partitions["soft_core"].append(node)#SOFT CORE
                    self.neighbors_graph.node[node]["partition_soft"]="soft_core"
                elif nb_orgs < self.nb_organisms:
                    self.partitions["soft_accessory"].append(node)#SOFT ACCESSORY
                    self.neighbors_graph.node[node]["partition_soft"]="soft_accessory"
                else:
                    logging.getLogger().error("nb_orgs can't be > to self.nb_organisms")
                    exit(1)

                self.neighbors_graph.nodes[node]["viz"]={}
                if nem_class != "U":
                    self.neighbors_graph.nodes[node]["viz"]['color']=COLORS_RGB[self.neighbors_graph.node[node]["partition"]]
                else:
                    self.neighbors_graph.nodes[node]["viz"]['color']=COLORS_RGB[self.neighbors_graph.node[node]["partition_exact"]]
                self.neighbors_graph.nodes[node]["viz"]['size']=nb_orgs
            if self.families_repeted_th > 0:
                if len(self.families_repeted)>0:
                    logging.getLogger().info("Gene families that have been discarded because there are repeated:\t"+" ".join(self.families_repeted))
                else:
                    logging.getLogger().info("No gene families have been discarded because there are repeated")

            logging.getLogger().debug(nx.number_of_edges(self.neighbors_graph))

            self.is_partitionned=True
        else:
            if just_stats:
                
                for node_name, nem_class in partitionning_results.items():
                    stats[SHORT_TO_LONG[nem_class]]+=1
                return stats
            else:
                return partitionning_results

    def partition_shell(self, nem_dir_path = tempfile.mkdtemp(),
                        subpart_name    = "subpartition_shell",
                        beta            = 0.5,
                        free_dispersion = False,
                        Q               = "auto",
                        exclusity_th    = 0.1,
                        init_using_qual = None):
        """

        """ 
        if not self.is_partitionned:
            logging.getLogger().warning("The pangenome must be already partionned to subpartition the shell genome")
        else:
            if Q == "auto":
                if init_using_qual is None:
                    stats = self.projection(nem_dir_path,self.organisms)
                    Q = float(len(self.partitions["shell"]))/stats[1]
                    Q = int(round(Q, 0))+1 # +1 to store unexclusive families
                else:
                    if isinstance(init_using_qual,dict):
                        Q = len(init_using_qual)+1
                    elif isinstance(init_using_qual,list):
                        Q = 3
                    else:
                        print("else")
            elif Q <= 1:
                logging.getLogger().error('Q must be above 1 or equals to "auto"')
                return()
            self.__write_nem_input_files(nem_dir_path+"/",
                                        self.organisms,
                                        init=init_using_qual,
                                        filter_by_partition="shell")
            subpartitions = run_partitioning(nem_dir_path, self.nb_organisms, beta, free_dispersion, Q = Q, init="random" if init_using_qual is None else "param_file")
            self.subpartitions_shell_parameters = {} 
            self.organisms_subpartitions_shell = defaultdict(set)
            proportion_exclusive = 0
            labels = {}
            for k, parameters in subpartitions[PARTITION_PARAMETERS].items():
                #pdb.set_trace()
                label = str(k)

                if mean(parameters[EPSILON])<exclusity_th:
                    label = label+"_exclusive:"+str(round(mean(parameters[EPSILON]),2))
                    proportion_exclusive+=parameters[PROPORTION]
                else:
                    label = label+"_shared:"+str(round(mean(parameters[EPSILON]),2))
                if isinstance(init_using_qual,dict):
                    coverages = {}
                    for group_name, set_org in init_using_qual.items():
                        k_set = set([org for (org, boolean) in zip(self.organisms, parameters[MU]) if boolean])
                        coverage=0
                        if len(k_set)>0:
                            coverage = float(len(k_set & set_org))/float(len(set_org))
                            if coverage >= 0.5:
                                coverages[group_name]=coverage
                        print(group_name+"size:"+str(len(set_org))+"  "+label+"  size: "+str(len(k_set))+"      inter:"+str(len(k_set & set_org))+"    coverage:"+str(coverage)) 
                    label = label+"_"+"|".join(sorted(coverages, key=coverages.get)) if len(coverages)>0 else label

                self.subpartitions_shell_parameters[label]=([org for (org, boolean) in zip(self.organisms, parameters[MU]) if boolean],
                                                       mean(parameters[EPSILON]),
                                                       parameters[PROPORTION],)
                for org in self.subpartitions_shell_parameters[label][0]:
                    self.organisms_subpartitions_shell[org].add(label)
                labels[k]=label
                print(parameters[EPSILON])
                print(parameters[PROPORTION])
                
            self.subpartition_shell = defaultdict(list)
            nx.set_node_attributes(self.neighbors_graph,nx.get_node_attributes(self.neighbors_graph, "partition"),subpart_name)
            for node, nem_class in subpartitions[FAMILIES_PARTITION].items():
                nb_orgs = 0
                self.neighbors_graph.node[node][subpart_name]=labels[nem_class]
                self.subpartition_shell[labels[nem_class]].append(node)
            return(Q)
            
    def compute_layout(self,
                       iterations = 500,
                       graph_type = "neighbors_graph",
                       outboundAttractionDistribution=True,  
                       linLogMode=False,  
                       adjustSizes=False,  
                       edgeWeightInfluence=1.0,
                       jitterTolerance=1.0, 
                       barnesHutOptimize=True,
                       barnesHutTheta=1.2,
                       multiThreaded=False,
                       scalingRatio=50000,
                       strongGravityMode=True,
                       gravity=1.0,
                       verbose=False):
        G=None
        if graph_type == "untangled_neighbors_graph":
            G = self.untangled_neighbors_graph
        else:
            G = self.neighbors_graph

        forceatlas2 = ForceAtlas2(
                          outboundAttractionDistribution=outboundAttractionDistribution,
                          linLogMode=linLogMode,
                          adjustSizes=adjustSizes, 
                          edgeWeightInfluence=edgeWeightInfluence,
                          jitterTolerance=jitterTolerance,
                          barnesHutOptimize=barnesHutOptimize,
                          barnesHutTheta=barnesHutTheta,
                          multiThreaded=multiThreaded,
                          scalingRatio=scalingRatio,
                          strongGravityMode=strongGravityMode,
                          gravity=gravity,
                          verbose=verbose)

        for node, pos_x_y in forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=iterations).items():
            z=(0,)
            if self.is_partitionned:
                if self.neighbors_graph.nodes[node]["partition"]=="persistent":
                    z=(2,)
                elif self.neighbors_graph.nodes[node]["partition"]=="shell":
                    z=(1,)
            self.neighbors_graph.nodes[node]["viz"]['position']=dict(zip(["x","y","z"],pos_x_y+z))

    def export_to_GEXF(self, graph_output_path, compressed=False, metadata = None, all_node_attributes = True, all_edge_attributes = True, graph_type = "neighbors_graph"):
        """
            Export the Partionned Pangenome Graph Of Linked Neighbors to a GEXF file  
            :param graph_output_path: a str containing the path of the GEXF output file
            :param compressed: a bool specifying if the file must be compressed in gzip or not
            :param metadata: a dict having the organisms as key emcompassing a dict having metadata as keys of its value as value 
            :param all_node_attributes: a bool specifying if organisms and genes attributes of each family node must be in the file or not.
            :param all_edge_attributes: a bool specifying if organisms count of each edge must be in the file or not.
            :type str: 
            :type bool: 
            :type bool: 
            :type dict: 
        """
        graph_to_save=None
        if "neighbors_graph":
            graph_to_save = self.neighbors_graph.copy()
        elif "untangled_neighbors_graph":
            graph_to_save = self.untangled_neighbors_graph.copy()
            
        for node in self.neighbors_graph.nodes():            
            for key in list(self.neighbors_graph.node[node].keys()):
                if key == "viz":
                    continue
                if not all_node_attributes and key in self.organisms:
                    del graph_to_save.node[node][key]
                else:
                    try:
                        if not isinstance(self.neighbors_graph.node[node][key], str):
                            graph_to_save.node[node][key]="|".join(self.neighbors_graph.node[node][key])#because networkx and gephi do not support list type in gexf despite it is possible according to the specification using liststring (https://gephi.org/gexf/1.2draft/data.xsd)
                    except TypeError:
                        if key == "length":
                            l = list(self.neighbors_graph.node[node][key])
                            graph_to_save.node[node]["length_avg"] = float(mean(l))
                            graph_to_save.node[node]["length_med"] = float(median(l))
                            graph_to_save.node[node]["length_min"] = min(l)
                            graph_to_save.node[node]["length_max"] = max(l)
                            del graph_to_save.node[node]["length"]
        for node_i, node_j, data in self.neighbors_graph.edges(data = True):
            l = list(data["length"])
            graph_to_save[node_i][node_j]["length_avg"] = float(mean(l))
            graph_to_save[node_i][node_j]["length_med"] = float(median(l))
            graph_to_save[node_i][node_j]["length_min"] = min(l)
            graph_to_save[node_i][node_j]["length_max"] = max(l)
            del graph_to_save[node_i][node_j]["length"]
            
            atts = set()
            for key in data.keys():
                
                if key in self.organisms:
                    if metadata:
                        for att, value in metadata[key].items():
                            atts.add(att)
                            try:
                                graph_to_save[node_i][node_j][att].add(value)
                            except KeyError:
                                graph_to_save[node_i][node_j][att]=set([value])
                    if not all_edge_attributes:
                        del graph_to_save[node_i][node_j][key] 

            for att in atts:
                graph_to_save[node_i][node_j][att]="|".join(sorted(graph_to_save[node_i][node_j][att]))

            graph_to_save[node_i][node_j]["viz"]={"thickness":graph_to_save[node_i][node_j]["weight"]}

        graph_output_path = graph_output_path+".gexf"
        if compressed:
            graph_output_path = gzip.open(graph_output_path+".gz","w")
        #pdb.set_trace()
        nx.write_gexf(graph_to_save, graph_output_path)

    # def import_from_GEXF(self, path_graph_to_update):
    #     """
    #         Import an already built Partionned Pangenome Graph Of Linked Neighbors from a GEXF file  
    #         :param nem_dir_path: a str containing the path of the GEXF input file (compressed or not)
    #         :type str: 

    #         .. warning:: please import a full GEXF input file having all the attributes
    #     """ 
    #     file = gzip.open(path_graph_to_update.name,"r")
    #     try:
    #         file.readline()
    #         file.seek(0)
    #         self.neighbors_graph = nx.read_gexf(file)
    #     except IOError as e:
    #         if e.message == "Not a gzipped file":
    #             self.neighbors_graph = nx.read_gexf(path_graph_to_update)
    #         else:
    #             logging.getLogger().error("Unable to open file "+path_graph_to_update)

    #     for node, data in self.neighbors_graph.nodes(data=True):
    #         for key, value in data.items():
    #             new_value = set(value.split('|'))
    #             self.neighbors_graph.node[node][key]=new_value
    #             if key not in RESERVED_WORDS:
    #                 self.organisms.add(key)
    #         if len(self.organisms) == 0:
    #             logging.getLogger().error("No attributes containing organisms names found on this node: "+str(node))
    #     logging.getLogger().debug(self.organisms)
    #     self.nb_organisms = len(self.organisms)

    #     for source, target, data in self.neighbors_graph.edges(data=True):
    #         try:
    #             del self.neighbors_graph[source][target]['id']
    #         except KeyError:
    #             logging.getLogger().warnings("No previous edge id found in gexf input file for edge: "+source+" <-> "+target)

    def write_matrix(self, path, header=True, csv = True, Rtab = True):
        """
            Export the pangenome as a csv_matrix similar to the csv et Rtab matrix exported by Roary (https://sanger-pathogens.github.io/Roary/)
            :param nem_dir_path: a str containing the path of the out files (csv+Rtab)
            :param header: a bool specifying if the header must be added to the file or not
            :type str: 
            :type bool: 
        """ 
        if self.is_partitionned:
            def write_file(ext, gene_or_not, sep):
                with open(path+"."+ext,"w") as matrix:
                    if header:
                        matrix.write(sep.join(['"Gene"',#1
                                               '"Non-unique Gene name"',#2
                                               '"Annotation"',#3
                                               '"No. isolates"',#4
                                               '"No. sequences"',#5
                                               '"Avg sequences per isolate"',#6
                                               '"Accessory Fragment"',#7
                                               '"Genome Fragment"',#8
                                               '"Order within Fragment"',#9
                                               '"Accessory Order with Fragment"',#10
                                               '"QC"',#11
                                               '"Min group size nuc"',#12
                                               '"Max group size nuc"',#13
                                               '"Avg group size nuc"']#14
                                               +['"'+org+'"' for org in list(self.organisms)])+"\n")#15

                    for node, data in self.neighbors_graph.nodes(data=True):
                        genes  = [('"'+"|".join(data[org])+'"' if gene_or_not else str(len(data[org]))) if org in data else ('""' if gene_or_not else "0") for org in self.organisms]
                        nb_org = len([gene for gene in genes if gene != ('""' if gene_or_not else "0")])
                        l = list(data["length"])
                        matrix.write(sep.join(['"'+node+'"',#1
                                               '"'+data["partition"]+'"',#2
                                               '"'+"|".join(data["product"])+'"',#3
                                               str(nb_org),#4
                                               str(data["nb_genes"]),#5
                                               str(round(data["nb_genes"]/nb_org,2)),#6
                                               '""',#data["subpartition_shell"],#7
                                               '""',#8
                                               '""',#9
                                               '""',#10
                                               '""',#11
                                               str(min(l)),#12
                                               str(max(l)),#13
                                               str(round(mean(l),2))]#14
                                               +genes)+"\n")#15
            if csv:
                logging.getLogger().info("Writing csv matrix")
                write_file("csv",True,",")
            if Rtab:
                logging.getLogger().info("Writing Rtab matrix")
                write_file("Rtab",False,"\t")
        else:
            logging.getLogger().error("The pangenome need to be partionned before being exported to a file matrix")
    # def delete_pangenome_graph(self, delete_NEM_files = False):
    #     """
    #         Delete all the pangenome graph of eventuelly the statistic of the partionning process (including the temporary file)
    #     """ 
    #     if delete_NEM_files:
    #         self.delete_nem_intermediate_files()

    #     self.nem_output               = None
    #     self.neighbors_graph          = None
    #     self.pan_size                 = 0
    #     self.nem_output               = None
    #     self.is_partitionned          = False
    #     self.partitions               = {}
    #     self.partitions["undefined"]  = list()
    #     self.partitions["persistent"] = list()
    #     self.partitions["shell"]      = list()
    #     self.partitions["cloud"]      = list()
    #     self.partitions["exact_core"] = list()
    #     self.partitions["exact_accessory"]  = list()
    #     self.BIC                      = None
    #     #self.partitions_by_organisms  = defaultdict(lambda: defaultdict(set))

    def delete_nem_intermediate_files(self):
        """
            Delete all the tempory files used to partion the pangenome
        """ 
        if self.nem_intermediate_files is not None:
            logging.getLogger().info("delete "+self.nem_intermediate_files)
            shutil.rmtree(self.nem_intermediate_files)
            self.nem_intermediate_files = None

    def ushaped_plot(self, out_file):
        """
            generate ushaped representation
            :param out_file: a str containing the path of the output file
            :type str: 
        """ 
        max_bar = 0
        count = defaultdict(lambda : defaultdict(int))
        for node, data in self.neighbors_graph.nodes(data=True):
            nb_org  = len([True for org in self.organisms if org in data])
            if self.is_partitionned:
                count[nb_org][data["partition"]]+=1
            count[nb_org]["pangenome"]+=1
            max_bar = count[nb_org]["pangenome"] if count[nb_org]["pangenome"] > max_bar else max_bar
        data_plot = []

        if self.is_partitionned and len(self.partitions["undefined"]) == 0:
            persistent_values = []
            shell_values      = []
            cloud_values      = []
            for nb_org in range(1,self.nb_organisms+1):
                persistent_values.append(count[nb_org]["persistent"])
                shell_values.append(count[nb_org]["shell"])
                cloud_values.append(count[nb_org]["cloud"])
            data_plot.append(go.Bar(x=list(range(1,self.nb_organisms+1)),y=persistent_values,name='persistent', marker=dict(color = COLORS["persistent"])))
            data_plot.append(go.Bar(x=list(range(1,self.nb_organisms+1)),y=shell_values,name='shell', marker=dict(color = COLORS["shell"])))
            data_plot.append(go.Bar(x=list(range(1,self.nb_organisms+1)),y=cloud_values,name='cloud', marker=dict(color = COLORS["cloud"])))
        else:
            text = 'undefined' if len(self.partitions["undefined"]) else "pangenome"
            undefined_values = []
            for nb_org in range(1,self.nb_organisms+1):
                undefined_values.append(count[nb_org][text])
            data_plot.append(go.Bar(x=list(range(1,self.nb_organisms+1)),y=undefined_values,name=text, marker=dict(color = COLORS[text])))
        layout = None
        if self.soft_core_th:
            x = self.nb_organisms*self.soft_core_th
            print(str(x))
            layout =  go.Layout(barmode='stack', shapes=[dict(type='line', x0=x, x1=x, y0=0, y1=max_bar, line = dict(dict(width=5, dash='dashdot', color="grey")))])
        else:
            layout = go.Layout(barmode='stack')

        fig = go.Figure(data=data_plot, layout=layout)
        out_plotly.plot(fig, filename = out_file+".html", auto_open=False)

    def tile_plot(self, outdir):
        """
            generate tile plot representation (not work for the moment)
            :param outdir: a str containing the path of the output file
            :type str: 
        """ 

        binary_data = []
        fam_order = []
        cpt = 1
        for node, data in self.neighbors_graph.nodes(data=True):
            fam_order.append(node)
            binary_data.append([1 if org in data else 0 for org in self.organisms])
            cpt+=1
            if cpt>1500:
                break
        binary_data = []
        fam_order=[]
        for org in self.organisms:
            l = []
            fam_order = []
            for node, data in self.neighbors_graph.nodes(data=True):
                fam_order.append(node)
                if org in data:
                    l.append(1)
                else:
                    l.append(0)
            binary_data.append(l)
        data_plot = go.Heatmap(z=binary_data,
                               x=list(self.organisms),
                               y=fam_order,
                               colorscale=[[0, 'rgb(0, 0, 0)'],[0, 'rgb(1, 1, 1)']],
                               colorbar={"tick0":0,"dtick":1})
        out_plotly.plot([data_plot], filename = outdir+"/tile_plot.html", auto_open=False)

        ##########

    # def identify_communities_in_each_partition(self):
    #     """
    #         Use the Louvain's algorithm to label the nodes by their community in each partition
    #     """ 
    #     size_communities=defaultdict(lambda : defaultdict(int))
    #     for partition in ["persistent","shell", "cloud"]:
    #         subgraph = self.neighbors_graph.subgraph([nodes for nodes,data in self.neighbors_graph.nodes(data=True) if data['partition']==partition])
    #         comm = community.best_partition(subgraph)# = nx.algorithms.community.asyn_fluidc(subgraph, 100)
    #         for node, id_com in comm.items():
    #             self.neighbors_graph.node[node]['community'] = partition+"_"+str(id_com)
    #             size_communities[partition][id_com]+=1

    # def identify_shell_subpaths(self, k_range = range(2,10),nem_dir_path = tempfile.mkdtemp()):
        
    #     subgraph = self.neighbors_graph.subgraph([nodes for nodes,data in self.neighbors_graph.nodes(data=True) if data['partition']=='Shell'])

    #     if not os.path.exists(nem_dir_path):
    #         #NEM requires 5 files: nem_file.index, nem_file.str, nem_file.dat, nem_file.m and nem_file.nei
    #         os.makedirs(nem_dir_path)
    #     self.nem_intermediate_files = nem_dir_path

    #     logging.getLogger().info("Writing nem_file.str nem_file.index nem_file.nei nem_file.dat and nem_file.m files")
    #     str_file   = open(nem_dir_path+"/nem_file.str", "w")
    #     index_file = open(nem_dir_path+"/nem_file.index", "w")
    #     org_file   = open(nem_dir_path+"/column_org_file", "w")
    #     nei_file   = open(nem_dir_path+"/nem_file.nei", "w")
    #     dat_file   = open(nem_dir_path+"/nem_file.dat", "w")

    #     str_file.write("S\t"+str(nx.number_of_nodes(subgraph))+"\t"+str(self.nb_organisms)+"\n")
    #     str_file.close()

    #     nei_file.write("1\n")#to enable weigthed partionning
        
    #     index = {node: index+1 for index, node in enumerate(subgraph.nodes(data=False))}
    #     index_inv = {i: node for node, i in index.items()}
    #     org_file.write(" ".join([org for org in self.organisms])+"\n")
    #     org_file.close()

    #     for node_name, node_organisms in subgraph.nodes(data=True):

    #         index_file.write(str(index[node_name])+"\t"+str(node_name)+"\n")
    #         logging.getLogger().debug(node_organisms)
    #         logging.getLogger().debug(self.organisms)
    #         dat_file.write("\t".join(["1" if org in node_organisms else "0" for org in self.organisms])+"\n")

    #         row_fam         = []
    #         row_dist_score  = []
    #         neighbor_number = 0
    #         try:
    #             for neighbor in nx.all_neighbors(subgraph, node_name):
    #                 #nb_presences = sum([pre_abs for org, pre_abs in self.neighbors_graph[node_name][neighbor].items() if org not in RESERVED_WORDS])
    #                 #self.neighbors_graph[node_name][neighbor]["weight"]= nb_presences
    #                 distance_score = subgraph[node_name][neighbor]["weight"]/self.nb_organisms
    #                 row_fam.append(str(index[neighbor]))
    #                 row_dist_score.append(str(round(distance_score,4)))
    #                 neighbor_number += 1
    #             if neighbor_number>0:
    #                 nei_file.write("\t".join([str(item) for sublist in [[index[node_name]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
    #                 #logging.getLogger().debug("\t".join([str(item) for sublist in [[[index[node_name]]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
    #             else:
    #                 raise nx.exception.NetworkXError("no neighbors in selected organismss")
    #         except nx.exception.NetworkXError as nxe:
    #             logging.getLogger().debug("The family: "+node_name+" is an isolated family")
    #             nei_file.write(str(index[node_name])+"\t0\n")

    #     index_file.close()
    #     nei_file.close()
    #     dat_file.close()

    #     for k in k_range:

    #         logging.getLogger().info("Running NEM uing "+str(k)+" class")
    #         # weighted_degree = sum(list(self.neighbors_graph.degree(weight="weight")).values())/nx.number_of_edges(self.neighbors_graph)
    #         # logging.getLogger().debug("weighted_degree: "+str(weighted_degree))
    #         # logging.getLogger().debug("org/weighted_degree: "+str(self.nb_organisms/weighted_degree))    
    #         #weighted_degree = sum(self.neighbors_graph.degree(weight="weight").values())/nx.number_of_edges(self.neighbors_graph)


    #         ALGO           = "ncem" #fuzzy classification by mean field approximation
    #         ITERMAX        = 100 # number of iteration max 
    #         MODEL          = "bern" # multivariate Bernoulli mixture model
    #         PROPORTION     = "pk" #equal proportion :  "p_"     varying proportion : "pk"
    #         VARIANCE_MODEL = "sk_" #one variance per partition and organism : "sdk"      one variance per partition, same in all organisms : "sd_"   one variance per organism, same in all partion : "s_d"    same variance in organisms and partitions : "s__" 
    #         NEIGHBOUR_SPEC = "f"# "f" specify to use all neighbors, orther argument is "4" to specify to use only the 4 neighbors with the higher weight (4 because for historic reason due to the 4 pixel neighbors of each pixel)
    #         CONVERGENCE_TH = "clas "+str(0.000001)

    #         HEURISTIC      = "heu_d"# "psgrad" = pseudo-likelihood gradient ascent, "heu_d" = heuristic using drop of fuzzy within cluster inertia, "heu_l" = heuristic using drop of mixture likelihood
    #         STEP_HEURISTIC = 1 # step of beta increase
    #         BETA_MAX       = float(self.nb_organisms) #maximal value of beta to test,
    #         DDROP          = 0.8 #threshold of allowed D drop (higher = less detection)
    #         DLOSS          = 0.5 #threshold of allowed D loss (higher = less detection)
    #         LLOSS          = 0.02 #threshold of allowed L loss (higher = less detection)
            
    #         BETA = ["-B",HEURISTIC,"-H",str(STEP_HEURISTIC),str(BETA_MAX),str(DDROP),str(DLOSS),str(LLOSS)]

    #         command = " ".join([NEM_LOCATION, 
    #                             nem_dir_path+"/nem_file",
    #                             str(k),
    #                             "-a", ALGO,
    #                             "-i", str(ITERMAX),
    #                             "-m", MODEL, PROPORTION, VARIANCE_MODEL,
    #                             "-s r 5",
    #                             *BETA,
    #                             "-n", NEIGHBOUR_SPEC,
    #                             "-c", CONVERGENCE_TH,
    #                             "-f fuzzy",
    #                             "-l y"])
         
    #         logging.getLogger().info(command)
    #         proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #         (out,err) = proc.communicate()  
    #         logging.getLogger().debug(out)
    #         logging.getLogger().debug(err)

    #         if os.path.isfile(nem_dir_path+"/nem_file.uf"):
    #             logging.getLogger().info("Reading NEM results...")
    #         else:
    #             logging.getLogger().error("No NEM output file found")
            
    #         with open(nem_dir_path+"/nem_file.uf","r") as classification_nem_file, open(nem_dir_path+"/nem_file.mf","r") as parameter_nem_file:
    #             classification = []
                
    #             parameter = parameter_nem_file.readlines()
    #             M = float(parameter[6].split()[3]) # M is markov ps-like
    #             BIC = -2 * M - (k * self.nb_organisms * 2 + k - 1) * math.log(len(self.partitions["Shell"]))

    #             logging.getLogger().info("The Bayesian Criterion Index of the partionning for "+str(k)+" is "+str(BIC))

    #             for i, line in enumerate(classification_nem_file):
    #                 elements = [float(el) for el in line.split()]
    #                 max_prob = max([float(el) for el in elements])
    #                 classes = [pos for pos, prob in enumerate(elements) if prob == max_prob]

    #                 self.neighbors_graph.node[index_inv[i+1]]["subshell"]=str(classes[0])

    def projection(self, out_dir, organisms_to_project):
        """
            generate files about the projection of the partition of the graph on the organisms
            return statistics about the number of genes in each organism to project for each partition in the file out_dir/nb_genes.csv
            :param outdir: a str containing the path of the output directory (name of files will be the name of organisms)
            :param organisms_to_project: a list of str containing the name of the organism
            :type str:
            :type list
            :return: stats: 
            :rtype: dict 
        """ 
        sep=","
        if self.is_partitionned:
            with open(out_dir+"/nb_genes.csv","w") as nb_genes_file:
                nb_genes_file.write("org\tpersistent\tshell\tcloud\texact_core\texact_accessory\tpangenome\n")
                for organism in organisms_to_project:
                    nb_genes_by_partition = defaultdict(int)
                    with open(out_dir+"/"+organism+".csv","w") as out_file:
                        out_file.write(sep.join(["gene","contig","coord_start","coord_end","strand","ori","family","nb_copy_in_org","partition","persistent","shell","cloud"])+"\n")
                        for contig, contig_annot in self.annotations[organism].items():
                            for gene, gene_info in contig_annot.items():
                                if gene_info[FAMILY] not in self.families_repeted:
                                    nb_genes_by_partition[self.neighbors_graph.node[gene_info[FAMILY]]["partition"]]+=1
                                    nb_genes_by_partition[self.neighbors_graph.node[gene_info[FAMILY]]["partition_exact"]]+=1
                                    nb_genes_by_partition[self.neighbors_graph.node[gene_info[FAMILY]]["partition_soft"]]+=1
                                    nb_genes_by_partition["pangenome"]+=1
                                    nei_partitions = [self.neighbors_graph.node[nei]["partition"] for nei in nx.all_neighbors(self.neighbors_graph,gene_info[FAMILY])]
                                    out_file.write(sep.join([gene,
                                                              contig,
                                                              str(gene_info[START]),
                                                              str(gene_info[END]),
                                                              gene_info[STRAND],
                                                              "T" if (gene_info[NAME].upper() == "DNAA" or gene_info[PRODUCT].upper() == "DNAA") else "F",
                                                              gene_info[FAMILY],
                                                              str(len(self.neighbors_graph.node[gene_info[FAMILY]][organism])),
                                                              self.neighbors_graph.node[gene_info[FAMILY]]["partition"],
                                                              str(nei_partitions.count("persistent")),
                                                              str(nei_partitions.count("shell")),
                                                              str(nei_partitions.count("cloud"))])+"\n")
                    self.partitions_by_organism[organism]=nb_genes_by_partition
                    nb_genes_file.write("\t".join([organism,
                                                  str(nb_genes_by_partition["persistent"]),
                                                  str(nb_genes_by_partition["shell"]),
                                                  str(nb_genes_by_partition["cloud"]),
                                                  str(nb_genes_by_partition["exact_core"]),
                                                  str(nb_genes_by_partition["exact_accessory"]),
                                                  str(nb_genes_by_partition["pangenome"])])+"\n")
        else:
            logging.getLogger().warning("The pangenome must be partionned before using this method (projection)")
        persistent_stats = []
        shell_stats = []
        cloud_stats = []
            
        for org, part in self.partitions_by_organism.items():
            persistent_stats.append(part["persistent"])
            shell_stats.append(part["shell"])
            cloud_stats.append(part["cloud"])
            
        return((mean(persistent_stats),mean(shell_stats),mean(cloud_stats),))
    
################ END OF CLASS PPanGGOLiN ################

################ FUNCTION run_partitioning ################
""" """
def run_partitioning(nem_dir_path, nb_org, beta, free_dispersion, Q = 3, init="param_file_default"):
    
    logging.getLogger().debug("Running NEM...")
    # weighted_degree = sum(list(self.neighbors_graph.degree(weight="weight")).values())/nx.number_of_edges(self.neighbors_graph)
    # logging.getLogger().debug("weighted_degree: "+str(weighted_degree))
    # logging.getLogger().debug("org/weighted_degree: "+str(self.nb_organisms/weighted_degree))    
    #weighted_degree = sum(self.neighbors_graph.degree(weight="weight").values())/nx.number_of_edges(self.neighbors_graph)

    ALGO           = b"ncem" #fuzzy classification by mean field approximation
    ITERMAX        = 100 # number of iteration max 
    MODEL          = b"bern" # multivariate Bernoulli mixture model
    PROPORTION     = b"pk" #equal proportion :  "p_"     varying proportion : "pk"
    VARIANCE_MODEL = b"skd" if free_dispersion else b"sk_"#one variance per partition and organism : "sdk"      one variance per partition, same in all organisms : "sd_"   one variance per organism, same in all partion : "s_d"    same variance in organisms and partitions : "s__" 
    #NEIGHBOUR_SPEC = "f"# "f" specify to use all neighbors, orther argument is "4" to specify to use only the 4 neighbors with the higher weight (4 because for historic reason due to the 4 pixel neighbors of each pixel)
    CONVERGENCE    = b"clas"
    CONVERGENCE_TH = 0.00000001
    (INIT_SORT, INIT_RANDOM, INIT_PARAM_FILE, INIT_FILE, INIT_LABEL, INIT_NB) = range(0,6)

    # HEURISTIC      = "heu_d"# "psgrad" = pseudo-likelihood gradient ascent, "heu_d" = heuristic using drop of fuzzy within cluster inertia, "heu_l" = heuristic using drop of mixture likelihood
    # STEP_HEURISTIC = 0.5 # step of beta increase
    # BETA_MAX       = float(len(organisms)) #maximal value of beta to test,
    # DDROP          = 0.8 #threshold of allowed D drop (higher = less detection)
    # DLOSS          = 0.5 #threshold of allowed D loss (higher = less detection)
    # LLOSS          = 0.02 #threshold of allowed L loss (higher = less detection)
    
    # BETA           = ["-B",HEURISTIC,"-H",str(STEP_HEURISTIC),str(BETA_MAX),str(DDROP),str(DLOSS),str(LLOSS)] if beta == float('Inf') else ["-b "+str(beta)]

    #WEIGHTED_BETA = beta#*nb_org
    # command = " ".join([NEM_LOCATION, 
    #                     nem_dir_path+"/nem_file",
    #                     str(Q),
    #                     "-a", ALGO,
    #                     "-i", str(ITERMAX),
    #                     "-m", MODEL, PROPORTION, VARIANCE_MODEL,
    #                     "-s m "+ nem_dir_path+"/nem_file.m",
    #                     "-b "+str(WEIGHTED_BETA),
    #                     "-n", NEIGHBOUR_SPEC,
    #                     "-c", CONVERGENCE_TH,
    #                     "-f fuzzy",
    #                     "-l y"])
 
    # logging.getLogger().debug(command)

    # proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE if logging.getLogger().getEffectiveLevel() == logging.INFO else None,
    #                                             stderr=subprocess.PIPE if logging.getLogger().getEffectiveLevel() == logging.INFO else None)
    # (out,err) = proc.communicate()
    # if logging.getLogger().getEffectiveLevel() == logging.INFO:
    #     with open(nem_dir_path+"/out.txt", "wb") as file_out, open(nem_dir_path+"/err.txt", "wb") as file_err:
    #         file_out.write(out)
    #         file_err.write(err)

    # logging.getLogger().debug(out)
    #logging.getLogger().debug(err)
    nem(Fname          = nem_dir_path.encode('ascii')+b"/nem_file",
        nk             = Q,
        algo           = ALGO,
        beta           = beta,
        convergence    = CONVERGENCE,
        convergence_th = CONVERGENCE_TH,
        format         = b"fuzzy",
        it_max         = ITERMAX,
        dolog          = True,
        model_family   = MODEL,
        proportion     = PROPORTION,
        dispersion     = VARIANCE_MODEL,
        init_mode      = INIT_PARAM_FILE if init.startswith("param_file") else INIT_RANDOM)
    # arguments_nem = [str.encode(s) for s in ["nem", 
    #                  nem_dir_path+"/nem_file",
    #                  str(Q),
    #                  "-a", ALGO,
    #                  "-i", str(ITERMAX),
    #                  "-m", MODEL, PROPORTION, VARIANCE_MODEL,
    #                  "-s m "+ nem_dir_path+"/nem_file.m",
    #                  "-b "+str(WEIGHTED_BETA),
    #                  "-n", NEIGHBOUR_SPEC,
    #                  "-c", CONVERGENCE_TH,
    #                  "-f fuzzy",
    #                  "-l y"]]

    # command = " ".join([NEM_LOCATION, 
    #                     nem_dir_path+"/nem_file",
    #                     str(Q),
    #                     "-a", ALGO,
    #                     "-i", str(ITERMAX),
    #                     "-m", MODEL, PROPORTION, VARIANCE_MODEL,
    #                     "-s m "+ nem_dir_path+"/nem_file.m",
    #                     "-b "+str(WEIGHTED_BETA),
    #                     "-n", NEIGHBOUR_SPEC,
    #                     "-c", CONVERGENCE_TH,
    #                     "-f fuzzy",
    #                     "-l y"])

    # logging.getLogger().debug(command)

    # proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE if logging.getLogger().getEffectiveLevel() == logging.INFO else None,
    #                                             stderr=subprocess.PIPE if logging.getLogger().getEffectiveLevel() == logging.INFO else None)
    # (out,err) = proc.communicate()
    # if logging.getLogger().getEffectiveLevel() == logging.INFO:
    #     with open(nem_dir_path+"/out.txt", "wb") as file_out, open(nem_dir_path+"/err.txt", "wb") as file_err:
    #         file_out.write(out)
    #         file_err.write(err)

    # logging.getLogger().debug(out)
    # logging.getLogger().debug(err)

    # array = (c_char_p * len(arguments_nem))()
    # array[:] = arguments_nem
    # nemfunctions.mainfunc(c_int(len(arguments_nem)), array)

    # if beta == float('Inf'):
    #     starting_heuristic = False
    #     with open(nem_dir_path+"/beta_evol.txt", "w") as beta_evol_file:
    #         for line in str(err).split("\\n"):
    #             if not starting_heuristic and line.startswith("* * Starting heuristic * *"):
    #                 beta_evol_file.write("beta\tLmix\n")
    #                 starting_heuristic = True
    #                 continue
    #             if starting_heuristic:
    #                 elements = line.split("=")
    #                 if elements[0] == " * * Testing beta ":
    #                    beta = float(elements[1].split("*")[0].strip())
    #                 elif elements[0] == "  criterion NEM ":
    #                     Lmix = float(elements[len(elements)-1].strip())
    #                     beta_evol_file.write(str(beta)+"\t"+str(Lmix)+"\n")
    
    if os.path.isfile(nem_dir_path+"/nem_file.uf"):
        logging.getLogger().debug("Reading NEM results...")
    else:
        logging.getLogger().warning("No NEM output file found: "+ nem_dir_path+"/nem_file.uf")
    index_fam = []
    with open(nem_dir_path+"/nem_file.index","r") as index_nem_file:
        for line in index_nem_file:
            index_fam.append(line.split("\t")[1].strip())
    
    partitions_list = ["U"] * len(index_fam)
    all_parameters = {}
    try:
        with open(nem_dir_path+"/nem_file.uf","r") as partitions_nem_file, open(nem_dir_path+"/nem_file.mf","r") as parameter_nem_file:
            parameter = parameter_nem_file.readlines()
            M = float(parameter[2].split()[3]) # M is markov ps-like
            BIC = -2 * M - (Q * nb_org * 2 + Q - 1) * math.log(len(index_fam))
            
            sum_mu_k = []
            sum_epsilon_k = []

            logging.getLogger().debug("The Bayesian Criterion Index of the partionning is "+str(BIC))
            for k, line in enumerate(parameter[-Q:]):
                logging.getLogger().debug(line)
                vector = line.split() 
                mu_k = [bool(float(mu_kj)) for mu_kj in vector[0:nb_org]]
                logging.getLogger().debug(mu_k)
                logging.getLogger().debug(len(mu_k))

                epsilon_k = [float(epsilon_kj) for epsilon_kj in vector[nb_org+1:]]
                logging.getLogger().debug(epsilon_k)
                logging.getLogger().debug(len(epsilon_k))
                proportion = float(vector[nb_org])
                logging.getLogger().debug(proportion)
                sum_mu_k.append(sum(mu_k))
                logging.getLogger().debug(sum(mu_k))
                sum_epsilon_k.append(sum(epsilon_k))
                logging.getLogger().debug(sum(epsilon_k))
                all_parameters[k]=(mu_k,epsilon_k,proportion)

            if init=="param_file_default":
                
                #persistent is defined by a sum of mu near of nb_organism and a low sum of epsilon
                max_mu_k     = max(sum_mu_k)
                persistent_k = sum_mu_k.index(max_mu_k)

                #shell is defined by an higher sum of epsilon_k
                max_epsilon_k = max(sum_epsilon_k)
                shell_k       = sum_epsilon_k.index(max_epsilon_k)

                # the other one should be cloud (basicaly with low sum_mu_k and low epsilon_k)
                cloud_k = set([0,1,2]) - set([persistent_k, shell_k])
                cloud_k = list(cloud_k)[0]

                # but if the difference between epsilon_k of shell and cloud is tiny, we check using the sum of mu_k which basicaly must be lower in cloud
                if ((sum_epsilon_k[shell_k]-sum_epsilon_k[cloud_k])/nb_org)<0.1 and sum_mu_k[shell_k]<sum_mu_k[cloud_k]:
                     # otherwise we permutate
                     (shell_k, cloud_k) = (cloud_k, shell_k)

                logging.getLogger().debug(sum_mu_k)
                logging.getLogger().debug(sum_epsilon_k)

                logging.getLogger().debug(persistent_k)
                logging.getLogger().debug(shell_k)
                logging.getLogger().debug(cloud_k)

                partition               = {}
                partition[persistent_k] = "P"#PERSISTENT
                partition[shell_k]      = "S"#SHELL
                partition[cloud_k]      = "C"#CLOUD

                if partition[0] != "P" or partition[1] != "S" or partition[2] != "C":
                    raise ValueError("vector mu_k and epsilon_k value in the mf file are not consistent with the initialisation value in the .m file")

            for i, line in enumerate(partitions_nem_file):
                elements = [float(el) for el in line.split()]
                max_prob = max([float(el) for el in elements])
                positions_max_prob = [pos for pos, prob in enumerate(elements) if prob == max_prob]
                logging.getLogger().debug(positions_max_prob)
                logging.getLogger().debug(i)
                
                if init=="param_file_default":
                    if (len(positions_max_prob)>1):
                        partitions_list[i]="S"#SHELL in case of doubt (equiprobable partition), gene families is attributed to shell
                    else:
                        partitions_list[i]=partition[positions_max_prob.pop()]
                else:
                    partitions_list[i]=positions_max_prob.pop()

            #logging.getLogger().debug(index.keys())
    except IOError:
        logging.getLogger().warning("Statistical partitioning do not works (the number of organisms used is probably too low), see logs here to obtain more details "+nem_dir_path+"/nem_file.log")
    except ValueError:
        ## return the default partitions_list which correspond to undefined
        pass
    return((dict(zip(index_fam, partitions_list)),all_parameters))

################ END OF FILE ################
