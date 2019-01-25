#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-

import warnings
warnings.filterwarnings("ignore")
import pandas
import numpy
import scipy
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
import json
import tempfile
from tqdm import tqdm
import random
from multiprocessing import Pool, Semaphore
from nem import *
from .utils import *
import pdb
from fa2 import ForceAtlas2
import plotly.plotly as py
import plotly.offline as out_plotly
import plotly.graph_objs as go
import plotly.figure_factory as ff
from ascii_graph import Pyasciigraph
from scipy.spatial.distance import squareform, pdist
from scipy.sparse import csr_matrix, bsr_matrix, csc_matrix
from scipy.stats import iqr, linregress, pearsonr
from scipy.cluster.hierarchy import linkage, dendrogram
import markov_clustering as mc
import io
from contextlib import redirect_stdout

import glob

(TYPE, FAMILY, START, END, STRAND, NAME, PRODUCT) = range(0, 7)#data index in annotation
(ORGANISM_INDEX,CONTIG_INDEX,POSITION_INDEX) = range(0, 3)#index
(ORGANISM_ID, ORGANISM_GFF_FILE) = range(0, 2)#data index in the file listing organisms
(GFF_seqname, GFF_source, GFF_type, GFF_start, GFF_end, GFF_score, GFF_strand, GFF_frame, GFF_attribute) = range(0,9)
(MU,EPSILON,PROPORTION) = range(0, 3)
(FAMILIES_PARTITION,PARTITION_PARAMETERS, LOG_LIKELIHOOD) = range(0, 3)
RESERVED_WORDS = set(["id", "label", "name", "weight", "partition", "partition_exact", "partition_soft", "length", "length_min", "length_max", "length_avg", "length_med", "product", 'nb_genes','subpartition',"viz","type","path","correlated_paths"])
SHORT_TO_LONG = {'EA':'exact_accessory','EC':'exact_core','SA':'soft_accessory','SC':'soft_core','P':'persistent','S':'shell','C':'cloud','U':'undefined'}
COLORS = {"pangenome":"black", "exact_accessory":"#EB37ED", "exact_core" :"#FF2828", "soft_core":"#e6e600", "soft_accessory":"#996633","shell": "#00D860", "persistent":"#F7A507", "cloud":"#79DEFF", "undefined":"#828282"}
COLORS_RGB = {"pangenome":{'r': 0, 'g': 0, 'b': 0, 'a': 0}, "exact_accessory":{'r': 235, 'g': 55, 'b': 237, 'a': 0}, "exact_core" :{'r': 255, 'g': 40, 'b': 40, 'a': 0},  "soft_core":{'r': 255, 'g': 255, 'b': 0, 'a': 0}, "soft_accessory": {'r': 153, 'g': 102, 'b': 51, 'a': 0},"shell": {'r': 0, 'g': 216, 'b': 96, 'a': 0}, "persistent":{'r': 247, 'g': 165, 'b': 7, 'a': 0}, "cloud":{'r': 121, 'g': 222, 'b': 255, 'a': 0}, "undefined":{'r': 130, 'g': 130, 'b': 130, 'a': 0}}
MAX_Q = 20

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
            >>>pan = PPanGGOLiN("args", annotations, organisms, circular_contig_size, families_repeted, directed)# load direclty the main attributes
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
        self.Q                              = None
        self.beta                           = None
        self.free_dispersion                = None
        self.th_degree                      = None
        self.chunk_size                     = None
        self.subpartitions_shell             = defaultdict(list)
        self.partition_parameters           = {}
        self.CDS_fragments                  = {}
        self.soft_core_th                   = None
        self.path_groups_vectors            = {}## key : id, value : vecteur numpy de moyenne de présence / absence des organismes pour chaque famille.
        self.path_vectors                   = {}## key : id, value : vecteur numpy de moyenne de présence / absence des organismes pour chaque famille.

        if init_from == "file":
            self.__initialize_from_files(*args)
        elif init_from == "args":
            (self.annotations,
             self.organisms,
             self.circular_contig_size,
             self.families_repeted,
             self.directed)
        elif init_from == "database":
            logging.getLogger().error("database is not yet implemented")
            pass
        elif init_from == "json":
            self.import_from_json(*args)
        else:
            raise ValueError("init_from parameter is required")
        self.nb_organisms = len(self.organisms)

        if init_from != "json":
            logging.getLogger().info("Computing gene neighborhood ...")
            self.__neighborhood_computation()

    def __initialize_from_files(self, organisms_file, families_tsv_file, lim_occurence = 0, infer_singletons = False, add_rna_to_the_pangenome = False, directed = False):
        """ 
            :param organisms_file: a file listing organims by compute, first column is organism name, second is path to gff file and optionnally other other to provide the name of circular contig
            :param families_tsv_file: a file listing families. The first element is the family identifier (by convention, we advice to use the identifier of the average gene of the family) and then the next elements are the identifiers of the genes belonging to this family.
            :param lim_occurence: a int containing the threshold of the maximum number copy of each families. Families exceeding this threshold are removed and are listed in the families_repeted attribute.
            :param add_rna_to_the_pangenome: a bool specifying if the rna genes must be added to the pangenome or not.
            :param infer_singletons: a bool specifying if singleton must be explicitely present in the families_tsv_file (False) or if single gene in gff files must be automatically infered as a singleton family (True)
            :type file: 
            :type file: 
            :type int: 
            :type bool: 
            :type int:
            :type bool: 
            :type bool: 
        """ 
        self.directed = directed
        logging.getLogger().info("Reading "+families_tsv_file.name+" the gene families file ...")

        families_tsv_file = read_compressed_or_not(families_tsv_file)
        organisms_file    = read_compressed_or_not(organisms_file)
        families    = dict()
        for line in families_tsv_file:
            elements = [el.strip() for el in line.split()] # 2 or 3 fields expected
            if len(elements)<=1:
                logging.getLogger().error("No tabulation separator found in gene families file")
                exit(1)
            (fam_id, gene_id, is_frag) = elements if len(elements) == 3 else elements+[None]
            families[gene_id]          = fam_id
            if is_frag == "F":
                self.CDS_fragments[gene_id] = fam_id

        self.circular_contig_size = {}
        logging.getLogger().info("Reading "+organisms_file.name+" the list of organism files ...")
        bar = tqdm(organisms_file,total=get_num_lines(organisms_file), unit = "gff file")

        for line in bar:
            elements = [el.strip() for el in line.split("\t")]
            if len(elements)<=1:
                logging.getLogger().error("No tabulation separator found in organisms file")
                exit(1)
            try:
                bar.set_description("Processing "+elements[ORGANISM_GFF_FILE])
                bar.refresh()
            except:
                if len(elements)>2:
                    self.circular_contig_size.update({contig_id: None for contig_id in elements[2:len(elements)]})  # size of the circular contig is initialized to None (waiting to read the gff files to fill the dictionnaries with the correct values)
            self.annotations[elements[0]] = self.__load_gff(elements[ORGANISM_GFF_FILE], families, elements[ORGANISM_ID], lim_occurence, infer_singletons, add_rna_to_the_pangenome)
        check_circular_contigs = {contig: size for contig, size in self.circular_contig_size.items() if size == None }
        if len(check_circular_contigs) > 0:
            logging.getLogger().error("""
                The following identifiers of circular contigs in the file listing organisms have not been found in any region type of the gff files: '"""+"'\t'".join(check_circular_contigs.keys())+"'")
            exit(1)

    def __load_gff(self, gff_file_path, families, organism, lim_occurence = 0, infer_singletons = False, add_rna_to_the_pangenome = False):
        """
            Load the content of a gff file
            :param gff_file_path: a valid gff file path where only type 'CDS' will be imported as genes. Each 'CDS' type must have a uniq ID as attribute (afterall called gene id).
            :param families: a dictionary having the gene as key and the identifier of the associated family as value. Depending on the infer_singletons attribute, singleton must be explicetly present on the dictionnary or not
            :param organism: a str containing the organim name
            :param lim_occurence: a int containing the threshold of the maximum number copy of each families. Families exceeding this threshold are removed and are listed in the next attribute.
            :param infer_singletons: a bool specifying if singleton must be explicitely present in the families parameter (False) or if single gene automatically infered as a singleton family (True)
            :param add_rna_to_the_pangenome: a bool specifying if the rna genes must be added to the pangenome or not.
            :type str:
            :type dict:
            :type str:
            :type int:
            :type bool:
            :type bool:
            :return: annot:
            :rtype: dict
        """

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

        def getIDAttribute(attribute):
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

        logging.getLogger().debug("Reading "+gff_file_path+" file ...")
        if organism not in self.organisms:
            self.organisms.add(organism)
            annot = defaultdict(OrderedDict)

            ctp_prev = 1
            cpt_fam_occ = defaultdict(int)
            gene_id_auto = False
            with read_compressed_or_not(gff_file_path) as gff_file:
                prevPseudoID = ""## default value. Can be anything but None, and should not be found in the ID field, for the default value.
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
                    if line.startswith('#!',0,2):## special refseq comment lines for versionning softs, assemblies and annotations.
                        continue
                    gff_fields = [el.strip() for el in line.split('\t')]

                    attributes = getGffAttributes(gff_fields)


                    if GFF_type == 'region':
                        if GFF_seqname in self.circular_contig_size:
                            self.circular_contig_size = int(GFF_end)
                            continue
                    elif gff_fields[GFF_type] == 'CDS' or (add_rna_to_the_pangenome and gff_fields[GFF_type].find("RNA")):
                        if "PSEUDO" in attributes:## if it is a pseudogene, this CDS is not actually translated.
                            continue
                        parent = attributes.get("PARENT")
                        if parent == prevPseudoID:## if the CDS has a Parent attribute and this parent corresponds to the last pseudogene element that was seen, this CDS is not actually translated.
                            continue
                        protein = getIDAttribute(attributes)
                        family = families.get(protein)
                        if not family:## the protein id was not found under "ID", searching elsewhere
                            proteinID = attributes.get("PROTEIN_ID")
                            if proteinID:
                                protein = proteinID
                                family = families.get(protein)

                        if not family:
                            if infer_singletons:## if we did not find the associated family at some point above, and if infer_singletons is set
                                families[protein] = protein
                                family            = families[protein]
                                logging.getLogger().info("infered singleton: "+protein)
                            else:
                                raise KeyError("Unknown families:"+protein+ ", check your families file or run again the program using the option to infer singleton")

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

                        annot[gff_fields[GFF_seqname]][protein] = [gff_fields[GFF_type],family,int(gff_fields[GFF_start]),int(gff_fields[GFF_end]),gff_fields[GFF_strand], name, product]

                    if attributes.get("PSEUDO") or attributes.get("PSEUDOGENE"):## if the element has the attribute pseudo, or pseudogene.
                        prevPseudoID = getIDAttribute(attributes)


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
            ret = "\n".join([l for l in Pyasciigraph(force_max_value=force_max_value, graphsymbol='*').graph(str(title), [("node(s) having degree "+str(i),v) for i, v in enumerate(values[:force_max_value+1])])])
            if len(values) > force_max_value:
                ret+="\n"
                ret+="And "+str(sum(values[50:]))+" nodes having degree above "+str(force_max_value)+"..."
            return(ret)

        pan_str ="\n"
        pan_str += "----------- Statistics -----------\n"
        pan_str += "Number of organisms:"+str(self.nb_organisms)+"\n"
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
            pan_str += "Q:"+str(self.Q)+"\n"
            pan_str += "beta:"+str(self.beta)+"\n"
            pan_str += "free dispersion:"+str(self.free_dispersion)+"\n"
            pan_str += "max node degree for smoothing:"+str(self.th_degree)+"\n"
            pan_str += "chunk size:"+"not partitionned by chunks\n " if self.chunk_size is None else str(self.chunk_size)+"\n"

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
        weight = sum(nx.get_edge_attributes(self.neighbors_graph, "weight").values())
        pan_str += "Sum of edges weight in the pangenome graph: "+str(weight)+"\n"
        pan_str += "Sum of normalized edges weight (divided by # of genomes) in the pangenome graph: "+str(round(weight/self.nb_organisms, 2))+"\n"
        try:
            pan_str += str_histogram("Degree distribution of the pangenome graph: ",nx.degree_histogram(self.neighbors_graph))+"\n"
        except:
            pan_str += "Degree distribution of the pangenome graph: "+str(nx.degree_histogram(self.neighbors_graph))+"\n"
        pan_str += "\n"
        pan_str += "----------------------------------"

        return(pan_str)

    def add_organisms(self, new_orgs, new_annotations, new_circular_contig_size, new_families_repeted):
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

    def __add_gene(self, fam_id, org, gene, name, length, product, type="CDS", graph_type = "neighbors_graph"):
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

        for attribute in ["name","length","product","type"]:
            try:
                graph.node[fam_id][attribute].add(locals()[attribute])
            except KeyError:
                graph.node[fam_id][attribute]=set([locals()[attribute]])

    def __add_link(self, fam_id, id, fam_id_nei, id_nei, organism, length, graph_type = "neighbors_graph"):
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
            graph[fam_id][fam_id_nei][organism].append({"source":id_nei,"target":id,"length":length})
        except:
            graph[fam_id][fam_id_nei][organism]= [ { "source" : id_nei,
                                                "target" : id,
                                                "length" : length}]
                                                             
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
                                gene_info_start[PRODUCT],
                                gene_info_start[TYPE])
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
                                        gene_info[PRODUCT],
                                        gene_info[TYPE])
                        self.index[gene]=(organism,contig,pos+1)
                        self.neighbors_graph.add_node(gene_info_nei[FAMILY])
                        if not (gene_info[FAMILY] == gene_info_nei[FAMILY] and
                                gene_info[STRAND] == gene_info_nei[STRAND] and
                                (gene in self.CDS_fragments or gene_nei in self.CDS_fragments)):# to avoid reflexive links with gene fragments
                            self.__add_link(fam_id      = gene_info[FAMILY],
                                            id          = gene,
                                            fam_id_nei  = gene_info_nei[FAMILY],
                                            id_nei      = gene_nei,
                                            organism   = organism,
                                            length      = gene_info[START] - gene_info_nei[END])
                        gene_nei, gene_info_nei = gene, gene_info
                
                if contig in self.circular_contig_size:#circularization
                    if not (gene_info_start[FAMILY] == gene_info_nei[FAMILY] and 
                            gene_info_start[STRAND] == gene_info_nei[STRAND] and 
                            (gene in self.CDS_fragments or gene_start in self.CDS_fragments)):# to avoid reflexive links with gene fragments
                        self.__add_link(fam_id      = gene_info_start[FAMILY],
                                        id          = gene_start,
                                        fam_id_nei  = gene_info_nei[FAMILY],
                                        id_nei      = gene_nei,
                                        organism   = organism,
                                        length      = (self.circular_contig_size[contig] - gene_info_nei[END]) + gene_info_start[START])
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
                                                self.__add_link(fam_id = self.annotations[org][self.index[gene][CONTIG_INDEX]][gene][FAMILY],
                                                                id = gene,
                                                                fam_id_nei = gene_info_prec[FAMILY],
                                                                id_nei = gene_prec,
                                                                organism = org,
                                                                length = length,
                                                                graph_type = "untangled_neighbors_graph")
                                                
                                                gene_prec, gene_info_prec = gene, gene_info

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
                if metaD is not None:
                    exclusively_in[variable] = defaultdict(set)
                    never_in[variable]       = defaultdict(set)
                    if value is not None:
                        variables_and_possible_values[variable].add(value)
        for node_name, node_organisms in self.neighbors_graph.nodes(data=True):
            for variable, possible_values in variables_and_possible_values.items():
                if output_path is not None:
                    try:
                        os.makedirs(output_path+"/exclusively_in/"+variable)
                        os.makedirs(output_path+"/never_in/"+variable)
                    except:
                        pass
                res = set([metadata[org][variable] for org in metadata if org in node_organisms]) - set([None])
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

    def __write_nem_input_files(self, nem_dir_path, select_organisms, old_nem_dir = None, th_degree = 20 ):
        
        if len(select_organisms)<=10:# below 10 organisms a statistical computation do not make any sence
            logging.getLogger().warning("The number of selected organisms is too low ("+str(len(select_organisms))+" organisms used) to partition the pangenome graph in persistent, shell and cloud genome. Add new organisms to obtain more robust metrics.")

        if not os.path.exists(nem_dir_path):
            #NEM requires 5 files: nem_file.index, nem_file.str, nem_file.dat, nem_file.m and nem_file.nei
            os.makedirs(nem_dir_path)

        # nx.set_edge_attributes(self.neighbors_graph,'tmp',nx.get_edge_attributes(self.neighbors_graph,"weight"))
        
        total_edges_weight = 0

        logging.getLogger().debug("Writing nem_file.str nem_file.index nem_file.nei and nem_file.dat files")
        with open(nem_dir_path+"/nem_file.str", "w") as str_file,\
            open(nem_dir_path+"/nem_file.index", "w") as index_file,\
            open(nem_dir_path+"/column_org_file", "w") as org_file,\
            open(nem_dir_path+"/nem_file.nei", "w") as nei_file,\
            open(nem_dir_path+"/nem_file.dat", "w") as dat_file:
            
            nei_file.write("1\n")
            if old_nem_dir:
                select_organisms = self.__write_nem_param_from_file(nem_dir_path, select_organisms, old_nem_dir) ## Writing NEM file to run NEM using former nem output.
            
            org_file.write(" ".join(["\""+org+"\"" for org in select_organisms])+"\n")
            org_file.close()

            index_fam = OrderedDict()
            for node_name, node_organisms in self.neighbors_graph.nodes(data=True):
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
                # if not self.is_partitionned or self.neighbors_graph.node["partition"]!="shell":
                try:
                    if nx.degree(self.neighbors_graph,node_name) < th_degree:
                        for neighbor in set(nx.all_neighbors(self.neighbors_graph, node_name)):
                            coverage = 0
                            if self.neighbors_graph.is_directed():
                                cov_sens, cov_antisens = (0,0)
                                try:
                                    cov_sens = sum([len(pre_abs) for org, pre_abs in self.neighbors_graph[node_name][neighbor].items() if ((org in select_organisms) and (org not in RESERVED_WORDS))])
                                except KeyError:
                                    pass
                                try:
                                    cov_antisens = sum([len(pre_abs) for org, pre_abs in self.neighbors_graph[neighbor][node_name].items() if ((org in select_organisms) and (org not in RESERVED_WORDS))])
                                except KeyError:
                                    pass
                                coverage = cov_sens + cov_antisens
                            else:
                                coverage = sum([len(pre_abs) for org, pre_abs in self.neighbors_graph[node_name][neighbor].items() if ((org in select_organisms) and (org not in RESERVED_WORDS))])

                            if coverage==0:
                                continue
                            
                            #distance_score = coverage/len(select_organisms)# len((set(self.neighbors_graph.node[node_name]) & set(self.neighbors_graph.node[neighbor])) - RESERVED_WORDS) & select_organisms)                          
                            distance_score = coverage/len(((set(self.neighbors_graph.node[node_name]) | set(self.neighbors_graph.node[neighbor])) - RESERVED_WORDS) & select_organisms)
                            total_edges_weight+=distance_score
                            row_fam.append(str(index_fam[neighbor]))
                            row_dist_score.append(str(round(distance_score,4)))
                            neighbor_number += 1
                        else:
                            logging.getLogger().debug("The family: "+node_name+" has a too high degree > "+str(th_degree)+" and will not be included in the spatial smoothing")
                # else:
                #     try:
                #         all_neighbors   = set(nx.all_neighbors(self.neighbors_graph, node_name))
                #         already_browsed = set()
                #         for neighbor_a in all_neighbors:
                #             already_browsed.add(neighbor_a)
                #             if neighbor_a self.neighbors_graph.node[neighbor_a]["partition"]!="persistent":
                #                 pass
                #             if len(all_neighbors - already_browsed):
                #                 for neighbor_b in all_neighbors - already_browsed:
                #                     path_cost = self.neighbors_graph[node_name][neighbor_a]['weight'] + self.neighbors_graph[node_name][neighbor_b]['weight']
                #                     alternative_path_cost = 0
                #                     total_alternative_path_cost = 0
                #                     removed_nodes_for_view = set()
                #                     removed_nodes_for_view.add(node_name)
                #                     def rule_filter_nodes(n):
                #                         if n in removed_nodes_for_view:
                #                             return(False)
                #                         else:
                #                             return(True)
                #                     while alternative_path_cost < path_cost:
                #                         total_alternative_path_cost+=alternative_path_cost
                #                         v = nx.classes.graphviews.subgraph_view(self.neighbors_graph, rule_filter_nodes)
                #                         try:
                #                             alternative_path = nx.dijkstra_path(v,neighbor_a,neighbor_b, weight = "weight")
                #                             alternative_path_cost = sum([d['weight'] for u,v,d in nx.subgraph(v,alternative_path).edges(data=True)])
                #                             removed_nodes_for_view.add(alternative_path[1])
                #                             removed_nodes_for_view.add(alternative_path[len(alternative_path)-1])
                #                         except nx.NetworkXNoPath:
                #                             break
                #                     self.neighbors_graph[node_name][neighbor_b]['tmp'] = self.neighbors_graph[node_name][neighbor_a]['tmp']-total_alternative_path_cost/2
                #                     self.neighbors_graph[node_name][neighbor_b]['tmp'] = self.neighbors_graph[node_name][neighbor_b]['tmp']-total_alternative_path_cost/2
                            
                                
                        # if self.neighbors_graph.is_directed():
                        #     cov_sens, cov_antisens = (0,0)
                        #     try:
                        #         cov_sens = sum([pre_abs for org, pre_abs in self.neighbors_graph[node_name][neighbor].items() if ((org in select_organisms) and (org not in RESERVED_WORDS))])
                        #     except KeyError:
                        #         pass
                        #     try:
                        #         cov_antisens = sum([pre_abs for org, pre_abs in self.neighbors_graph[neighbor][node_name].items() if ((org in select_organisms) and (org not in RESERVED_WORDS))])
                        #     except KeyError:
                        #         pass
                        #     coverage = cov_sens + cov_antisens
                        # else:
                        #     coverage = sum([pre_abs for org, pre_abs in self.neighbors_graph[node_name][neighbor].items() if ((org in select_organisms) and (org not in RESERVED_WORDS))])

                        # if coverage==0:
                        #     continue
                            
                    if neighbor_number>0:
                        nei_file.write("\t".join([str(item) for sublist in [[index_fam[node_name]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
                    else:
                        nei_file.write(str(index_fam[node_name])+"\t0\n")
                        logging.getLogger().debug("The family: "+node_name+" has a too high degree or is an isolated family in the selected organisms")
                except nx.exception.NetworkXError as nxe:
                    print(nxe)
                    logging.getLogger().debug("The family: "+node_name+" is an isolated family")
                    nei_file.write(str(index_fam[node_name])+"\t0\n")

            str_file.write("S\t"+str(len(index_fam))+"\t"+
                            str(len(select_organisms))+"\n")
        
        return(total_edges_weight)

    def __evaluate_nb_partitions(self, nem_dir_path = tempfile.mkdtemp(),
                                 old_nem_dir        = None,
                                 select_organisms   = None,
                                 free_dispersion    = False,
                                 Qmin_Qmax          = [3,MAX_Q],
                                 ICL_th             = 0.01,
                                 nb_iter            = 10,
                                 nb_threads         = 1,
                                 seed               = 42):
        if select_organisms is None:
            select_organisms = self.organisms
        else:
            select_organisms = OrderedSet(select_organisms)
            if len(select_organisms - self.organisms)>0:
                raise Exception("select_organisms parameter must be included in the organisms attribute of the objet")
        if old_nem_dir:
            init = "init_from_old"
        else:
            init = "param_file"
        
        edges_weight = self.__write_nem_input_files(nem_dir_path,select_organisms, old_nem_dir = old_nem_dir)

        def run_several_quick_partitioning (all_Q_to_partition):
            all_log_likelihood = []
            with contextlib.closing(Pool(processes = nb_threads)) if nb_threads>1 else empty_cm() as pool:
                    # logging.disable(logging.INFO)# disable INFO message
                    # logging.disable(logging.WARNING)# disable WARNING message
                if nb_threads>1: 
                    all_log_likelihood.extend(pool.starmap(run_partitioning, [(nem_dir_path,
                                                                               len(select_organisms),
                                                                               0,#quick, beta=0
                                                                               free_dispersion,
                                                                               q,
                                                                               seed,
                                                                               init,
                                                                               nb_iter,#quick
                                                                               True) for q in all_Q_to_partition]))
                else:
                    for q in all_Q_to_partition:
                        all_log_likelihood.append(run_partitioning(nem_dir_path,
                                                                len(select_organisms),
                                                                0, #quick, beta=0
                                                                free_dispersion,
                                                                q,
                                                                seed,
                                                                init,
                                                                nb_iter,#quick, only 10 iterations
                                                                True))
            logging.disable(logging.NOTSET)# restaure message
            all_BICs = defaultdict(float)
            all_ICLs = defaultdict(float)
            all_LLs  = defaultdict(float)
            for Q_candidate, log_likelihood, entropy in all_log_likelihood:
                if log_likelihood is not None:
                    all_BICs[Q_candidate] = calculate_BIC(log_likelihood,Q_candidate * (len(select_organisms) + 1 + (len(select_organisms) if free_dispersion else 1)),self.pan_size)
                    all_ICLs[Q_candidate] = all_BICs[Q_candidate] - entropy
                    all_LLs[Q_candidate]  = log_likelihood
            return(tuple([all_BICs,all_ICLs,all_LLs]))
        #kneedle    = None
        best_icl_Q     = None
        max_icl_Q      = 0
        slope      = 0
        intercept  = 0
        BICs,ICLs,LLs = run_several_quick_partitioning(list(range(Qmin_Qmax[0]-1,Qmin_Qmax[1]+1)))
        #mean_BICs = {q: mean(q_BICs) for q, q_BICs in BICs.items()}
        slope, intercept, r_value, p_value, std_err = (0,0,None,None,None)
        if len(BICs)>3:
            slope, intercept, r_value, p_value, std_err = linregress(list(ICLs.keys()), list(ICLs.values()))
            if slope > 0:
                max_icl_Q  = max(ICLs, key=ICLs.get)
                delta_ICL  = (ICLs[max_icl_Q]-min(ICLs.values()))*ICL_th
                best_Q = min({q for q, icl in ICLs.items() if icl>=ICLs[max_icl_Q]-delta_ICL})                
                # try:
                #     with redirect_stdout(io.StringIO()):# to capture error message from KneeLocator
                #         kneedle = KneeLocator(Qs, BICs, direction='decreasing',curve='convex')
                #     if kneedle.knee is None or kneedle.knee < 3:
                #         best_Q = max_Q if self.Q is None else self.Q# if no knee was found, its means that no discrete gain is obtain by increasing Q
                #     else:
                #         best_Q = kneedle.knee
                #         logging.getLogger().debug(" ".join([str(Qs[Q_idx]) for Q_idx in kneedle.xmx_idx]))
                # except ValueError:
                #     best_Q = 3# if KneeLocator raise a ValueError it means that overpartionning do not works, so that a unique central shell partition is enough
            else:
                best_Q = 3# if slope is increasing it means that increasing the number of partitions is useless, so 3 partition is enough
        else:
            best_Q = 3# not enough data, so 3 partition is enough
        best_Q = 3 if best_Q < 3 else best_Q
        # otherwise add new point ?
        # putative_best_Q = kneedle.knee
        #all_log_likelihood = run_several_quick_partitioning([putative_best_Q-1,putative_best_Q+1])
        # valid_Q , all_BICs = ([],[])
        # for Q, log_likelihood in zip(all_Q,all_log_likelihood):
        #     if log_likelihood is not None:
        #         valid_Q.append(Q)
        #         all_BICs.append(calculate_BIC(log_likelihood,Q * (len(select_organisms) + 1 + (len(select_organisms) if free_dispersion else 1)),self.pan_size))
        # print(valid_Q)
        # print(all_BICs)
        # kneedle = KneeLocator(valid_Q, all_BICs,direction='decreasing',curve='convex')
        # # otherwise add new point ?
        if len(BICs)>0:
            traces = []
            # traces.append(go.Scatter(x=[q for q, bics in BICs.items() for bic in bics],
            #                          y=[bic for q, bics in BICs.items() for bic in bics],
            #                          name = "all BICs",
            #                          mode = "markers"))
            traces.append(go.Scatter(x=list(BICs.keys()),
                                     y=list(BICs.values()),
                                     name = "BIC",
                                     mode = "lines+markers"))
            traces.append(go.Scatter(x=list(ICLs.keys()),
                                     y=list(ICLs.values()),
                                     name = "ICL",
                                     mode = "lines+markers"))
            traces.append(go.Scatter(x=list(LLs.keys()),
                                     y=list(LLs.values()),
                                     name = "log likelihood",
                                     mode = "lines+markers"))
            layout = go.Layout(title = 'ICL curve (best Q is '+str(best_Q)+', ICL_th= is '+str(ICL_th)+")",#, "+ ("y = "+str(round(slope,2))+"x + "+str(round(intercept,2))+", r = "+str(round(r_value,2)) if r_value else ""
                               titlefont = dict(size = 20),
                               xaxis = dict(title='number of overpartitions'),
                               yaxis = dict(title='ICL, BIC, log likelihood'),
                               shapes=[dict(type='line', x0=best_Q, x1=best_Q, y0=0, y1=ICLs[best_Q], line = dict(dict(width=1, dash='dashdot', color="black"))),
                                       dict(type='line', x0=max_icl_Q, x1=max_icl_Q, y0=0, y1=ICLs[max_icl_Q], line = dict(dict(width=1, dash='dashdot', color="black"))),
                                       dict(type='line', x0=best_Q, x1=max_icl_Q, y0=ICLs[max_icl_Q], y1=ICLs[max_icl_Q], line = dict(dict(width=1, dash='dashdot', color="black"))),
                                       dict(type='line', x0=2, x1=Qmin_Qmax[1], y0=ICLs[best_Q], y1=ICLs[best_Q], line = dict(dict(width=1, dash='dashdot', color="black")))])
                                       #dict(type='line', x0=2, x1=max_icl_Q, y0=(3*slope)+intercept, y1=(max_icl_Q*slope)+intercept, line = dict(dict(width=1, dash='dashdot', color="grey")))])
            fig = go.Figure(data=traces, layout=layout)
            out_plotly.plot(fig, filename=nem_dir_path+"/ICL_curve_Q"+str(best_Q)+".html", auto_open=False)
        return(best_Q, edges_weight)

    def __write_nem_param_from_file(self, nem_dir_path, select_organisms, old_nem_dir):
        
        old_orgs = open(old_nem_dir + "/column_org_summary","r").readline().replace('"','').split() ## former organisms used in precedent partitionning.
        new_orgs = [ orgs for orgs in select_organisms if orgs not in old_orgs ]## listing the new organisms
        ## listing removed organisms index (to not write the variables from the nem files corresponding to this organism)
        rmOrgsIndex = []
        kept_old = []
        for org in range(len(old_orgs)):
            if old_orgs[org] not in select_organisms:
                rmOrgsIndex.append(org)## saving the index of the organism not present in selected_organisms
            else:
                kept_old.append(old_orgs[org])## else saving the organism.
        select_organisms = OrderedSet(kept_old + new_orgs)## ordered of the kept organisms is known and the same than former NEM files.
        oldOrgsLen = len(old_orgs)
        newOrgsLen = len(new_orgs)
        
        fname = glob.glob(old_nem_dir + "*summary.mf", recursive = False)[0]
        with open(fname,"r") as parameter_nem_file:
            parameter = parameter_nem_file.readlines()
            
            MuLists = []
            EpsLists = []
            PkList = []
            for line in parameter:
                currParams = line.split()
                
                currMu = currParams[0:oldOrgsLen]
                PkList.append(currParams[oldOrgsLen])
                currEps = currParams[oldOrgsLen+1:]

                MuList = [currMu[index] for index in range(oldOrgsLen) if index not in rmOrgsIndex ]
                EpsList = [currEps[index] for index in range(oldOrgsLen) if index not in rmOrgsIndex ]
                
                nbs = []
                for val in MuList:
                    if val == '0.5':## cause int('0.5') raises an error.
                        nbs.append(0.5)## setting to 1 at random.
                    else:
                        nbs.append(int(val))
                
                MuList.extend([str(int(median(nbs)))] * newOrgsLen)
                EpsList.extend([str(median(map(float, EpsList)))] * newOrgsLen)
                
                ### NEM does not allow epsilon to be 0, so in case of 0s we fix epsilon to 1e10-6
                for n, i in enumerate(EpsList):
                    if float(i) == 0:
                        EpsList[n] = "0.000001"
                
                MuLists.append(MuList)
                EpsLists.append(EpsList)
                
        with open(nem_dir_path+"/nem_file_init_"+str(len(MuLists))+".m", "w") as m_file:
            m_file.write("1 ")# 1 to initialize parameter, 
            m_file.write(" ".join(PkList[:-1]) + " ")
            for MuList in MuLists:
                m_file.write( " ".join(MuList) + " ")
            for EpsList in EpsLists:
                m_file.write( " ".join(EpsList) + " ")
                
        return select_organisms
                                    
    def partition(self, nem_dir_path     = tempfile.mkdtemp(),
                        old_nem_dir      = None,
                        select_organisms = None,
                        Q                = -1,
                        Qmin_Qmax        = [3,MAX_Q],
                        beta             = 0.5,
                        th_degree        = float("inf"),
                        free_dispersion  = False,
                        chunck_size      = 500,
                        soft_core_th     = 0.95,
                        ICL_th           = 0.01,
                        inplace          = True,
                        just_stats       = False,
                        nb_threads       = 1,
                        seed             = 42):
        """
            Use the graph topology and the presence or absence of genes from each organism into families to partition the pangenome in three groups ('persistent', 'shell' and 'cloud')
            . seealso:: Read the Mo Dang's thesis to understand NEM, a summary is available here : http://www.kybernetika.cz/content/1998/4/393/paper.pdf
            :param nem_dir_path: a str containing a path to store temporary file of the NEM program.
            :param old_nem_dir: a str containing a path where nem files from another run are stored.
            :param select_organisms: a list of organism to used to obtain the partition (must be included in the organism attributes of the object) or None to used all organisms in the object
            :param beta: a float containing the spatial coefficient of smoothing of the clustering results using the weighted graph topology (0.00 turn off the spatial clustering)
            :param free_dispersion: a bool specyfing if the dispersion around the centroid vector of each paritition is the same for all the organisms or if the dispersion is free
            :param chunck_size: an int specifying the size of the chunks
            :param soft_core_th a float between 0 and 1 providing the threshold ratio of presence to attribute a gene families to the soft core genome
            :param inplace: a boolean specifying if the partition must be stored in the object of returned (throw an error if inplace is true and organisms parameter i not None)
            :param just_stats: a boolean specifying if the partitions must be returned or just stats about them (number of families in each partition)
            :param nb_threads: an integer specifying the number of threads to use (works only if the number of organisms is higher than the chunck_size)
            :type str: 
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
        
        if old_nem_dir:## then the NEM files will be written before run_partitionning.
            init = "init_from_old"
        else:## else, write generic nem files
            init = "param_file"
        
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
        
        def run_evaluate_nb_partitions(orgs, Q):
            if Q == -1:
                if self.Q == 3:
                    Q = 3
                    edges_weight = self.__write_nem_input_files(nem_dir_path,orgs, old_nem_dir = old_nem_dir, th_degree = th_degree)
                else:
                    if inplace:
                        logging.getLogger().info("Estimating the optimal number of partitions...")
                    (Q,edges_weight) = self.__evaluate_nb_partitions(nem_dir_path     = nem_dir_path,
                                                                     old_nem_dir      = old_nem_dir,
                                                                     select_organisms = orgs,
                                                                     free_dispersion  = free_dispersion,
                                                                     Qmin_Qmax        = Qmin_Qmax if self.Q is None else [3,self.Q],
                                                                     ICL_th           = ICL_th,
                                                                     nb_iter          = 10,
                                                                     nb_threads       = nb_threads,
                                                                     seed             = seed)
                    if inplace:
                        logging.getLogger().info("Best Q is "+str(Q))
            else:
                edges_weight = self.__write_nem_input_files(nem_dir_path,orgs, old_nem_dir = old_nem_dir, th_degree = th_degree)
            return(Q, edges_weight)
        
        if len(select_organisms) > chunck_size:
            Q,edges_weight = run_evaluate_nb_partitions(OrderedSet(random.sample(select_organisms,chunck_size)),Q)
            if inplace:
                logging.getLogger().info("Partitioning...")
                self.chunck_size = chunck_size
            cpt_partition = OrderedDict()
            for fam in families:
                cpt_partition[fam]= {"P":0,"S":0,"C":0,"U":0}
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
                        cpt_partition[node][nem_class[0]]+=1
                        sum_partionning = sum(cpt_partition[node].values())
                        if (sum_partionning > len(select_organisms)/chunck_size and max(cpt_partition[node].values()) >= sum_partionning*0.5) or (sum_partionning > len(select_organisms)):
                            if node not in validated:
                                if inplace:
                                    bar.update()
                                if max(cpt_partition[node].values()) < sum_partionning*0.5:
                                    cpt_partition[node]["U"] = sys.maxsize #if despite len(select_organisms) partionning, an abosolute majority is not found then the families is set to undefined 
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
                    if (sem.acquire() if nb_threads>1 else True):#
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
                        orgs = random.sample(select_organisms,chunck_size)
                        orgs = OrderedSet(orgs)

                        # for org, p in proba_sample.items():
                        #     if org in orgs:
                        #         proba_sample[org] = p - len(select_organisms)/chunck_size if p >1 else 1
                        #     else:
                        #         proba_sample[org] = p + len(select_organisms)/chunck_size

                        edges_weight = self.__write_nem_input_files(nem_dir_path+"/"+str(cpt)+"/",orgs, old_nem_dir = old_nem_dir)
                        if nb_threads>1:
                            res = pool.apply_async(run_partitioning,
                                                   args = (nem_dir_path+"/"+str(cpt)+"/",#nem_dir_path
                                                           len(orgs),
                                                           beta,#*((stats["exact_accessory"]+stats["exact_core"])/edges_weight),
                                                           free_dispersion,
                                                           Q,
                                                           seed,
                                                           init),                                                        
                                                   callback = validate_family)
                        else:
                            res = run_partitioning(nem_dir_path+"/"+str(cpt)+"/",#nem_dir_path
                                                   len(orgs),
                                                   beta, #*((stats["exact_accessory"]+stats["exact_core"])/edges_weight),
                                                   free_dispersion,
                                                   Q,
                                                   seed,
                                                   init)
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

            partitionning_results = [partitionning_results,[]]

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
            Q,edges_weight = run_evaluate_nb_partitions(select_organisms,Q)
            if inplace:
                logging.getLogger().info("Partitioning...")

            partitionning_results = run_partitioning(nem_dir_path, len(select_organisms), beta , free_dispersion, Q = Q, seed = seed, init = init)# * ((stats["exact_accessory"]+stats["exact_core"])/edges_weight)

            #partitionning_results = partitionning_results[FAMILIES_PARTITION]

            # all_Q = []
            # all_BIC = []
            # pool.map(f, range(10))
            # for i in range(3,20):
            #     print("here (__write_nem_input_files) "+str(i)+ "   "+nem_dir_path+"Q_"+str(i)+"/")
            #     self.__write_nem_input_files(nem_dir_path+"Q_"+str(i)+"/",select_organisms, Q=i)
            #     print("here (run_partitioning) "+str(i))
            #     partitionning_results = run_partitioning(nem_dir_path+"Q_"+str(i)+"/", len(select_organisms), beta, free_dispersion, Q=i)
            #     BIC = partitionning_results[2]

            #     all_Q.append(i)
            #     all_BIC.append(BIC)

            #     print("here (after run_partitioning) "+str(i)+"   BIC="+str(BIC))
            # kneedle = KneeLocator(all_Q, all_BIC, curve='concave', direction='increasing')
            # best_Q = kneedle.knee
            # print("here (__write_nem_input_files) "+str(i)+ "   "+nem_dir_path+"best_Q_"+str(best_Q)+"/")
            # print("here (run_partitioning) "+str(best_Q))

            
            # s=0
            # p=0
            # #
            # try:
            #     for k, (mu,epsilon,prop) in partitionning_results[PARTITION_PARAMETERS].items():
            #         if k!= 0 and k!= i:
            #             s += sum([e*prop for e in epsilon])
            #             p += prop
            #     #print(s/p/len(select_organisms))
            # except:
            #     pass

            #= s*i
            
            #pdb.set_trace()
        
            # #if partitionning_results[]

            
            # #partitionning_results = partitionning_results[PARTITION_PARAMETERS]                
                
        if inplace:
            #self.BIC = BIC
            self.Q                    = Q
            self.beta                 = beta
            self.free_dispersion      = free_dispersion
            self.th_degree            = th_degree
            self.partition_parameters = partitionning_results[PARTITION_PARAMETERS]
            if self.is_partitionned:
                for p in SHORT_TO_LONG.values():
                    self.partitions[p] = list()# erase older values
                self.subpartitions_shell = defaultdict(list)
                # for p in self.subpartitions_shell:
                #     self.subpartitions_shell[p] = list()# erase older values
            for node, nem_class in partitionning_results[FAMILIES_PARTITION].items():
                nb_orgs=0
                for key in list(self.neighbors_graph.node[node].keys()):
                    if key not in RESERVED_WORDS:
                        #self.partitions_by_organisms[key][partition[int(nem_class)]].add(self.neighbors_graph.node[node][key])
                        nb_orgs+=1

                self.neighbors_graph.node[node]["partition"]=SHORT_TO_LONG[nem_class[0]]
                self.partitions[SHORT_TO_LONG[nem_class[0]]].append(node)
                if nem_class[0] == "S":
                    self.subpartitions_shell[nem_class].append(node)
                self.neighbors_graph.node[node]["subpartition"]=nem_class
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
                for node_name, nem_class in partitionning_results[FAMILIES_PARTITION].items():
                    stats[SHORT_TO_LONG[nem_class[0]]]+=1
                return stats
            else:
                return partitionning_results

    def compute_layout(self,
                       iterations = 5,
                       graph_type = "neighbors_graph",
                       outboundAttractionDistribution=True,  
                       linLogMode=False,  
                       adjustSizes=False,  
                       edgeWeightInfluence=1.0,
                       jitterTolerance=1.0, 
                       barnesHutOptimize=True,
                       barnesHutTheta=1.2,
                       multiThreaded=False,
                       scalingRatio=500,
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

    def _setOrgsAsNodeAttr(self, graph):
        """
            Creates an organisms set and stores it as a graph attribute. Stores in nodes every protein ID from every organism belonging to the gene family.
            Assumes the ppanggolin organisation of the pangenome graph with networkx.
            Saves also the annotation data for each gene beloning to the processed node.
        """		
        # ~ start = time.time()
        def fillGene(annotation_data, position):
            geneDict = { "position" : position, "strand" : annotation_data[STRAND], "end" : annotation_data[END], "start" : annotation_data[START], "name" : annotation_data[NAME], "product" : annotation_data[PRODUCT] }
            return geneDict

        def fillContig(annotation_data, position, gene):
            contigDict = { gene : fillGene(annotation_data, position)}
            return contigDict

        def fillOrg(annotation_data, position, gene, contig_name):
            orgDict = { contig_name : fillContig(annotation_data, position, gene)}
            return orgDict
        ## prepare annotation object to parse the data faster than using self.annotations directly.
        annotation = self.annotations.copy()
        for org in self.annotations.keys():
            for contig, genes in self.annotations[org].items():
                annotation[org][contig] = list(genes.values())

        for node in graph.nodes.data():
            
            nodeObj = graph.nodes[node[0]]
            ## Counters to establish the most commun name, product and length of the gene family represented by this node.
            names_count = Counter()
            product_count = Counter()
            length_count = Counter()

            attrs = node[1].copy()
            nodeObj["organisms"] = {}
            for att in attrs:
                if att in self.organisms:## then att is an organism.
                    

                    if type(attrs[att]) == str:
                        nodeObj["organisms"][att].add(attrs[att])
                    elif type(attrs[att]) in [ list, set ]:
                        for gene in attrs[att]:
                            gene_data = self.index[gene]
                            annotation_data = annotation[gene_data[0]][gene_data[1]][gene_data[2]]
                            names_count[annotation_data[NAME]] += 1
                            product_count[annotation_data[PRODUCT]] += 1
                            length_count[annotation_data[END] - annotation_data[START]] +=1

                            ## Process the structure and add stuff depending on what's already present.
                            org = nodeObj["organisms"].get(gene_data[0])
                            if not org :## There's no gene from that organism for now, we add the organism, with the contig and the gene
                                nodeObj["organisms"][gene_data[0]] = fillOrg(annotation_data, gene_data[2], gene, gene_data[1])
                            else:## there's already a gene from this gene_data[0] organism in this family
                                contig = org.get(gene_data[1])
                                if not contig:## no gene from that contig in this gene family for now, we add the contig with the gene
                                    org[gene_data[1]] = fillContig(annotation_data, gene_data[2], gene)
                                else:## there's already a gene from this gene_data[1] contig in this organism in this family, we add the gene.
                                    contig[gene] = fillGene(annotation_data, gene_data[2])

                            # nodeObj["organisms"][att].add(gene)    
                    else:
                        raise TypeError("Unexpected type in a node attribute, expecting only str or list.")
                    
                    del nodeObj[att]## deleting the attribute as it is be in the node's organisms attribute.
            
            nodeObj["name"] = names_count.most_common(1)[0][0]
            nodeObj["product"] = product_count.most_common(1)[0][0]
            nodeObj["length"] = length_count.most_common(1)[0][0]
            ## assuming that a gene family cannot really have different type?
            if len(nodeObj["type"]) > 1:## this is assuming the data is stored in a set. Remove the warning from code if you change that in the data structure.
                logging.getLogger().warning("For the gene family " + str(node[0]) + " multiple sequence types were found ( " + ",".join(nodeObj["type"]) + "). This is quite unexpected from a gene family. Only the first one was kept.")
            nodeObj["type"] = list(nodeObj["type"])[0]
    
    def _setOrgsAsEdgeAttr(self, graph):
        """
            Stores in edge attribute 'organisms' the presence of an edge in an organism.
        """
        # ~ start = time.time()            
        for node_i, node_j, data in graph.edges(data = True):
            orgsPA = {}
            currData = data.copy()
            length_count = Counter()
            for att in currData:
                if att in self.organisms:
                    orgsPA[att] = currData[att]
                    for link in currData[att]:## usually there is only one, but there can be more.
                        length_count[link["length"]] +=1
                    
                    del graph[node_i][node_j][att]

            graph[node_i][node_j]["length"] = length_count.most_common(1)[0][0]
            graph[node_i][node_j]["organisms"] = orgsPA

    def _setGraphAttr(self, graph):
        """
            Fills the graph with attributes and parameters to save.
        """

        ## organisms names, their affiliated contigs, and whether they are circular or not
        graph.graph["organisms"] = {}
        for org in self.annotations.keys():
            graph.graph["organisms"][org] = {}
            for contig in self.annotations[org].keys():
                ## filling the graph's 'organisms' attribute
                if contig in self.circular_contig_size.keys():
                    graph.graph["organisms"][org][contig] = { "is_circular" : True }
                else:
                    graph.graph["organisms"][org][contig] = { "is_circular" : False }

        ## output parameters of each nem partition.
        graph.graph["params"] = {}
        # print(self.partition_parameters)
        ## self.partition_parameters is inconsistent and empty when multithreaded...
        ## For now saving only if it has the expected organisation, which is a dictionnary with the partition ID as key and a 3-values tuple as value.
        if isinstance(self.partition_parameters, dict):
            for key, val in self.partition_parameters.items():
                graph.graph["params"][key] = { "mu":[ int(mu) for mu in val[0] ] , "eps":val[1],"pk":val[2] }
        
        graph.graph["soft_core_threshold"] = self.soft_core_th 
        graph.graph["number_of_partitions"] = self.Q
        graph.graph["beta"] = self.beta
        graph.graph["free_dispersion"] = self.free_dispersion
        graph.graph["pan_size"] = self.pan_size
        graph.graph["is_partitionned"] = self.is_partitionned
        graph.graph["th_degree"] = self.th_degree
        graph.graph["families_repeted_th"] = self.families_repeted_th
    
    def export_to_json(self, graph_output_path, compressed = False, metadata = None, all_node_attributes = True, all_edge_attributes = True, graph_type = "neighbors_graph"):
        """
            Export the Partionned Pangenome Graph Of Linked Neighbors to a JSON file  
            :param graph_output_path: a str containing the path of the JSON output file
            :param compressed: a bool specifying if the file must be compressed in gzip or not
            :param metadata: a dict having the organisms as key emcompassing a dict having metadata as keys of its value as value 
            :param all_node_attributes: a bool specifying if organisms and genes attributes of each family node must be in the file or not.
            :param all_edge_attributes: a bool specifying if organisms count of each edge must be in the file or not.
            :type str: 
            :type bool: 
            :type bool: 
            :type dict: 
        """

        graph = None
        if graph_type =="neighbors_graph":
            graph = self.neighbors_graph.copy()
        elif graph_type == "untangled_neighbors_graph":
            graph = self.neighbors_graph.copy()

        logging.getLogger().debug("Setting orgs as node attributes")
        self._setOrgsAsNodeAttr(graph)
        logging.getLogger().debug("Setting orgs as edge attributes")
        self._setOrgsAsEdgeAttr(graph)
        logging.getLogger().debug("setting the graph attributes")
        self._setGraphAttr(graph)
        logging.getLogger().debug("Turning graph into JSON-friendly datastructure and writing JSON file.")
        
        graph_output_path = graph_output_path+".json"
        if compressed:
            graph_output_path = gzip.open(graph_output_path+".gz","wt")
        else:
            graph_output_path = open(graph_output_path, "w")

        json.dump(nx.node_link_data(graph), graph_output_path, cls = PanEncoder, indent = None)

    def import_from_json(self, graph_input_path, lim_occurence = 0, infer_singletons = False, add_rna_to_the_pangenome = False, directed = False):
        """
            Import the Partionned Pangenome Graph Of Linked Neighbors from a JSON file  
            :param graph_output_path: a str containing the path of the JSON input file
            :type str: 
        """
        if type(graph_input_path) == list:
            if len(graph_input_path) ==1:
                graph_input_path = graph_input_path[0]

        def readNodesAttr(graph):
            #self.annotations # Dict [key = organism], of dict [ key = contig ID] of OrderedDict [ keys = gene IDs, values = tuple of annotation_data ordered as (TYPE, FAMILY, START, END, STRAND, NAME, PRODUCT)  ]
            #self.index, Dict [key = gene ID], values is tuple of ( organism, contig, position on contig)
            for node in graph.nodes:
                gene_type = graph.nodes[node]["type"]
                for org in graph.nodes[node]["organisms"]:
                    graph.nodes[node][org] = [] ##
                    for contig in graph.nodes[node]["organisms"][org]:
                        for gene in graph.nodes[node]["organisms"][org][contig]:
                            gene_data = graph.nodes[node]["organisms"][org][contig][gene]
                            self.annotations[org][contig][gene] = (gene_type, node, gene_data["start"], gene_data["end"], gene_data["strand"], gene_data["name"], gene_data["product"])
                            self.index[gene] = (org, contig, gene_data["position"])
                            graph.nodes[node][org].append(gene)
                
                del graph.nodes[node]["organisms"]


            ## need to sort the OrderedDict of self.annotations as the order they are supposed to be in :
            logging.getLogger().debug("Done with storing data in self.annotations and self.index. Now need to order self.annotations[orgs][contig] OrderedDicts properly.")
            for org in self.annotations.keys():
                for contig in self.annotations[org].keys():
                    self.annotations[org][contig] = OrderedDict(sorted(self.annotations[org][contig].items(), key = lambda t : self.index[t[0]][2] ))## ordered by the gene position stored in index !
                    

            for node in graph.nodes:
                for attr in graph.nodes[node].keys():
                    if attr in ["name","length","product","type"]:
                        graph.nodes[node][attr] = set([graph.nodes[node][attr]])## storing as set since it is initialized as such when using gff files and protein clusters.
                    if attr in SHORT_TO_LONG.keys():
                        self.partitions[SHORT_TO_LONG[attr]].append(node)
                
                if graph.nodes[node]["partition"] == "shell":
                    self.subpartitions_shell[graph.nodes[node]["subpartition"]].append(node)

        def readEdgesAttr(graph):
            for node_i, node_j in graph.edges():

                length_list = set()## stored as a set currently.
                for org in graph[node_i][node_j]["organisms"]:
                    ## changing to the organisation used in ppanggolin currently.
                    graph[node_i][node_j][org] = graph[node_i][node_j]["organisms"][org].copy()
                    length_list.update([ contig["length"] for contig in graph[node_i][node_j]["organisms"][org] ])
                del graph[node_i][node_j]["organisms"]

                graph[node_i][node_j]["length"] = length_list

        def readGraphAttr(graph):
            for org in graph.graph["organisms"].keys():
                self.organisms.add(org)
                self.annotations[org] = {}
                for contig in graph.graph["organisms"][org].keys():
                    self.annotations[org][contig] = OrderedDict()
                    if graph.graph["organisms"][org][contig]["is_circular"] == True:
                        self.circular_contig_size[contig] = None### value of the intergenic region of the genes at the contig's "extremities" is unknown from the graph attributes alone. Saving those informations for later.

            for key in graph.graph["params"].keys():
                vals = graph.graph["params"][key]
                self.partition_parameters[key] = ([str(mu) for mu in vals["mu"] ], vals["eps"], vals["pk"])
            
            self.soft_core_th  = graph.graph["soft_core_threshold"]
            self.Q = graph.graph["number_of_partitions"]
            self.beta = graph.graph["beta"]
            self.free_dispersion = graph.graph["free_dispersion"]
            self.pan_size = graph.graph["pan_size"]
            self.is_partitionned = graph.graph["is_partitionned"]
            self.th_degree = graph.graph["th_degree"]
            self.families_repeted_th = graph.graph["families_repeted_th"]

        logging.getLogger().info("Loading the pangenome graph from the provided json file :" + graph_input_path)
        logging.getLogger().debug("Getting the graph object and metadata from json, uncompressing if need be")
        graph = nx.node_link_graph(json.load(read_compressed_or_not(graph_input_path)))
        logging.getLogger().debug("Reading graph attributes")
        readGraphAttr(graph)
        logging.getLogger().debug("Reading edges attributes")
        readEdgesAttr(graph)
        logging.getLogger().debug("Reading nodes attributes")
        readNodesAttr(graph)
        self.neighbors_graph = graph

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
                    # if metadata and key in metadata:
                    #     for att, value in metadata[key].items():
                    #         if value is not None:
                    #             atts.add(att)
                    #             try:
                    #                 graph_to_save[node_i][node_j][att].add(value)
                    #             except KeyError:
                    #                 graph_to_save[node_i][node_j][att]=set([value])
                    if not all_edge_attributes:
                        del graph_to_save[node_i][node_j][key]
                    else:
                        graph_to_save[node_i][node_j][key] = len(graph_to_save[node_i][node_j][key])
                # if key == "path_group":
                #     if metadata:
                #         path_group = self.neighbors_graph[node_i][node_j]["path_group"]
                #         orgs = set([org for org, v in self.path_groups_vectors[path_group].items() if round(v)>0]) & data.keys()
                #         for o in orgs:
                #             if o in metadata:
                #                 for att, value in metadata[o].items():
                #                     if value is not None:
                #                         try:
                #                             graph_to_save[node_i][node_j]["path_group_"+att].add(value)
                #                         except KeyError:
                #                             graph_to_save[node_i][node_j]["path_group_"+att]=set([value])
            for att in atts:
                graph_to_save[node_i][node_j][att]="|".join(sorted(graph_to_save[node_i][node_j][att]))
                # if "path_group_"+att in graph_to_save[node_i][node_j]:
                #     graph_to_save[node_i][node_j]["path_group_"+att]="|".join(sorted(graph_to_save[node_i][node_j]["path_group_"+att]))
            try:
                graph_to_save[node_i][node_j]["viz"]["thickness"] = graph_to_save[node_i][node_j]["weight"]
            except:
                graph_to_save[node_i][node_j]["viz"] = {"thickness" : graph_to_save[node_i][node_j]["weight"]}
            graph_to_save[node_i][node_j]["viz"]["color"] = average_color(graph_to_save.node[node_i]["viz"]["color"],graph_to_save.node[node_j]["viz"]["color"])
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

    def write_parameters(self, path):
        if self.is_partitionned:
            with open(path,"w") as parameters_file:
                for partition, (mu, ep, pi) in self.partition_parameters.items():
                    parameters_file.write(partition+":\n")
                    parameters_file.write("mu_qj:"+",".join([str(int(m)) for m in mu])+"\n")
                    parameters_file.write("epsilon_qj:"+",".join([str(e) for e in ep])+"\n")
                    parameters_file.write("pi_q:"+str(pi)+"\n")
                    parameters_file.write("===\n")

    def write_melted_matrix(self, path):
        if self.is_partitionned:
            with open(path+".csv","w") as melted_matrix:
                for node, data in self.neighbors_graph.nodes(data=True):
                    for key, d in data.items():
                        if key in self.organisms:
                            melted_matrix.write('"'+'"'.join([node,key]+list(d))+'"\n')

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
                                               '"'+(data["path_group"] if "path_group" in data else "NA")+'"',#7
                                               '"'+(data["path"] if "path" in data else "NA")+'"',#8
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

    def delete_nem_intermediate_files(self):
        """
            Delete all the tempory files used to partion the pangenome
        """ 
        if self.nem_intermediate_files is not None:
            logging.getLogger().info("delete "+self.nem_intermediate_files)
            shutil.rmtree(self.nem_intermediate_files)
            self.nem_intermediate_files = None
            
    def keep_nem_intermediate_files(self, outdir):
        """
            Saves the temporary files used to partition the pangenome, and creates summary files to reuse for futur partitionning.
            :param out_file: a str containing the path of the output directory for nem files
            :type str: 
        """
        if self.nem_intermediate_files:
            shutil.copytree(self.nem_intermediate_files, outdir)
            ## checking if the analysis has outputs in chunks
            subdirs = [ x[0] for x in os.walk(self.nem_intermediate_files) if x[0] != self.nem_intermediate_files ]
            
            if len(subdirs) > 0:
                ## then create a summary matrice... ?
                self.make_nem_matrix_summary(subdirs, outdir)
            # ~ for subdir in subdirs:
                # ~ basename = os.path.basename(subdir)
                # ~ shutil.copytree(subdir, outdir + "/" + basename)
            ## copying files of the corresponding partition
            else:
                self.make_nem_matrix_summary([self.nem_intermediate_files], outdir)
            
    def make_nem_matrix_summary(self, subdirs, outdir):
        """
            Saves a resulting global matrix from NEM submatrices.
            :param outdir: a str containing the path of the output directory for nem files
            :param subdirs: a list containing the subdirectories of the temporary directory where the nem files are stored
            :type str: 
            :type list:
        """
        matrice = dict()
        OrgSet = OrderedSet()
        PkLists = []
        
        for i in range(self.Q):
            PkLists.append([])
        ## gathering data
        for subdir in subdirs:
            fname = glob.glob(subdir + "/*_" + str(self.Q) + ".mf", recursive = False)[0]
            orgs = open(subdir + "/column_org_file","r").readline().replace('"','').split() ## orgs present in this subpartitionning.
            
            for org in orgs:
                if org not in OrgSet:
                    matrice[org] = {"mu":[], "eps":[]}
                    for x in range(self.Q):
                        matrice[org]["mu"].append([])
                        matrice[org]["eps"].append([])
                    OrgSet.add(org)
            
            
                            
            with open(fname,"r") as parameter_nem_file:
                parameter = parameter_nem_file.readlines()
                MuLists = []
                EpsLists = []
                currPkList = []
                lineNb=0
                for line in parameter[8:]:
                    currParams = line.split()
                    for paramIndex in range(len(orgs)):
                        matrice[orgs[paramIndex]]["mu"][lineNb].append(currParams[paramIndex])
                        matrice[orgs[paramIndex]]["eps"][lineNb].append(currParams[len(orgs)+ 1 + paramIndex])
                        
                    currPkList.append(currParams[len(orgs)])
                    lineNb+=1
                
                for pk in range(len(currPkList)):
                    PkLists[pk].append(currPkList[pk])
        PkList = []
        
        org_file = open(outdir + "/column_org_summary","w")
        org_file.write(" ".join(["\""+org+"\"" for org in OrgSet])+"\n")
        org_file.close()
        mf_file = open(outdir+"/nem_file_summary.mf", "w")
        
        for i in range(self.Q):
            currMu = []
            currEps = []
            PkList.append(str(median(map(float, PkLists[i]))))
            
            for org in OrgSet:
                nbs = []
                for val in matrice[org]["mu"][i]:
                    if val == '0.5':## cause int('0.5') raises an error.
                        nbs.append(0.5)## setting to 1 at random.
                    else:
                        nbs.append(int(val))
                
                currMu.append(str(int(median(nbs))))
                currEps.append(str(round(median(map(float, matrice[org]["eps"][i])), 6)))
               
            mf_file.write( " ".join(currMu) + " ")
            mf_file.write(str(median(map(float, PkLists[i]))) + " ")
            mf_file.write( " ".join(currEps) + "\n")
            
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
        chao = "NA"
        if count[1]["pangenome"] > 0:
            chao = round(self.pan_size + ((count[0]["pangenome"]^2)/(count[1]["pangenome"]*2)),2)

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
            layout =  go.Layout(title = "Gene families frequency distribution (U shape), chao="+str(chao),
                                xaxis = dict(title='Occurring in x genomes'),
                                yaxis = dict(title='# of gene families (F)'),
                                barmode='stack', shapes=[dict(type='line', x0=x, x1=x, y0=0, y1=max_bar, line = dict(dict(width=5, dash='dashdot', color="grey")))])
        else:
            layout = go.Layout(title = "Gene frequency distribution (U shape), chao="+str(chao),
                               xaxis = dict(title='Occurring in x genomes'),
                               yaxis = dict(title='# of gene families (F)'),
                               barmode='stack')

        fig = go.Figure(data=data_plot, layout=layout)
        out_plotly.plot(fig, filename = out_file+".html", auto_open=False)

    def tile_plot(self, outdir, shell_persistent_only = False):
        """
            generate tile plot representation
            :param outdir: a str containing the path of the output file
            :param outdir: a bool specifying if only the persistent and shell genome is plotted or all the pangenome
            :type str: 
        """ 
        data        = []
        all_indexes = []
        all_columns = []
        nodes_order = bidict()
        if shell_persistent_only:
            graph       = nx.Graph.subgraph(self.neighbors_graph,self.partitions["shell"]+self.partitions["persistent"])
        else:
            graph       = self.neighbors_graph
        for row, (node_name, node_organisms) in enumerate(graph.nodes(data=True)):
            new_col=[i for i, org in enumerate(self.organisms) if org in node_organisms]
            all_indexes.extend([row]*len(new_col))
            all_columns.extend(new_col)
            data.extend([1.0]*len(new_col))
            nodes_order[row]=node_name

        mat_p_a = csc_matrix((data, (all_indexes,all_columns)), shape = (len(graph),self.nb_organisms), dtype='float')
        dist    = pdist(1 - jaccard_similarities(mat_p_a,0).todense())
        hc      = linkage(dist, 'single')

        dendro = dendrogram(hc,no_plot=True)
        order_organisms = self.organisms[dendro['leaves']]

        binary_data = []
        text_data   = []
        fam_order   = []

        ordored_nodes_p = sorted(self.partitions["persistent"], key=lambda n:len(graph.nodes[n]), reverse=True)
        ordored_nodes_s = sorted(self.partitions["shell"], key=lambda n:len(graph.nodes[n]), reverse=True)
        ordored_nodes_c = sorted(self.partitions["cloud"], key=lambda n:len(graph.nodes[n]), reverse=True)

        for node in ordored_nodes_p+ordored_nodes_s+ordored_nodes_c:
            fam_order.append(node)
            data = graph.nodes[node]
            binary_data.append([len(data[org]) if org in data else numpy.nan for org in order_organisms])
            text_data.append([("\n".join(data[org])) if org in data else numpy.nan for org in order_organisms])
        # fam_order=[]
        # for org in self.organisms:
        #     l = []
        #     fam_order = []
        #     for node, data in self.neighbors_graph.nodes(data=True):
        #         fam_order.append(node)
        #         if org in data:
        #             l.append(1)
        #         else:
        #             l.append(0)
        #     binary_data.append(l)

        heatmap = go.Heatmap(z              = binary_data,
                             x              = list(self.organisms),
                             y              = fam_order,
                             text           = text_data,
                             zauto          = False,
                             zmin           = 1,
                             zmax           = 2,
                             autocolorscale = False,
                             colorscale     = [[0.50, 'rgb(100, 15, 78)'],[1, 'rgb(59, 157, 50)']],
                             colorbar       = dict(title     = 'Presence/Absence',
                                                   titleside = 'top',
                                                   tickmode  = 'array',
                                                   tickvals  = [1,2],
                                                   ticktext  = ['Presence','Multicopy'],
                                          #    tick0    = 0,
                                          #    dtick    = 0.333,
                                          #    nticks   = 3,
                                                   ticks     = 'outside'))

        sep1 = len(ordored_nodes_p)-0.5
        sep2 = sep1 + len(ordored_nodes_s)
        sep3 = sep2 + len(ordored_nodes_c)

        layout = go.Layout(title  = "presence/absence matrix",
                           xaxis  = dict(title='organisms'),
                           yaxis  = dict(title='gene families'),
                           shapes = [dict(type='line', x0=-1, x1=-1, y0=0, y1=sep1, line = dict(dict(width=10, color=COLORS["persistent"]))),
                                     dict(type='line', x0=self.nb_organisms, x1=self.nb_organisms, y0=0, y1=sep1, line = dict(dict(width=10, color=COLORS["persistent"]))),
                                     dict(type='line', x0=-1, x1=self.nb_organisms, y0=sep1, y1=sep1, line = dict(dict(width=1, color=COLORS["persistent"]))),
                                     dict(type='line', x0=-1, x1=-1, y0=sep1, y1=sep2, line = dict(dict(width=10, color=COLORS["shell"]))),
                                     dict(type='line', x0=self.nb_organisms, x1=self.nb_organisms, y0=sep1, y1=sep2, line = dict(dict(width=10, color=COLORS["shell"]))),
                                     dict(type='line', x0=-1, x1=self.nb_organisms, y0=sep2, y1=sep2, line = dict(dict(width=1, color=COLORS["shell"]))),
                                     dict(type='line', x0=-1, x1=-1, y0=sep2, y1=sep3, line = dict(dict(width=10, color=COLORS["cloud"]))),
                                     dict(type='line', x0=self.nb_organisms, x1=self.nb_organisms, y0=sep2, y1=sep3, line = dict(dict(width=10, color=COLORS["cloud"])))])

        # figure = ff.create_dendrogram(mat_p_a, orientation='bottom', labels=self.organisms,distfun = lambda x: dist, linkagefun = lambda x: hc)
        # for i in range(len(figure['data'])):
        #     figure['data'][i]['yaxis'] = 'y2'
        # pdb.set_trace()
        # heatmap[0]['x'] = figure['layout']['xaxis']['tickvals']
        # figure.add_traces(heatmap)
        
        out_plotly.plot(go.Figure(data=[heatmap], layout=layout), filename = outdir+"/tile_plot.html", auto_open=False)

        binary_data = []
        path_order  = []
        if self.path_vectors:
            layout = go.Layout(title  = "presence/absence matrix",
                               xaxis  = dict(title='organisms'),
                               yaxis  = dict(title='paths'))
            for path, vector in self.path_vectors.items():
                path_order.append("p"+path)
                binary_data.append(vector)
            heatmap = go.Heatmap(z   = binary_data,
                                x    = list(self.organisms),
                                y    = path_order,
                                zauto = False,
                                zmin = 0,
                                zmax = 1,
                                autocolorscale = False,
                                colorscale=[[0.50, 'rgb(100, 15, 78)'],[1, 'rgb(59, 157, 50)']],
                                colorbar = dict(title     = 'Presence ratio',
                                                titleside = 'top',
                                                # tickmode  = 'array',
                                                # tickvals  = [0,1],
                                                # ticktext  = ['Presence','Multicopy'],
                                            #    tick0    = 0,
                                            #    dtick    = 0.333,
                                            #    nticks   = 3,
                                                ticks     = 'outside'))
            
            out_plotly.plot(go.Figure(data=[heatmap], layout=layout), filename = outdir+"/tile_plot_paths.html", auto_open=False)

    def extract_shell_paths(self,jaccard_similarity_th = 0.7, inflation_mcl_path_groups = 1.5, breaker_th = 1.0/3.0):
        #graph   = nx.Graph.subgraph(self.neighbors_graph,self.partitions["shell"]).copy()
        #graph   = self.neighbors_graph.copy()
        #mat_p_a = pandas.DataFrame(False, columns=graph.nodes(data=False), index=self.organisms, dtype=numpy.bool)
        
        for subpartition in self.subpartitions_shell:
            data        = []
            all_indexes = []
            all_columns = []
            graph       = nx.Graph.subgraph(self.neighbors_graph,self.subpartitions_shell[subpartition])
            nodes_order = bidict()
            for col, (node_name, node_organisms) in enumerate(graph.nodes(data=True)):
                new_ind          = [i for i, org in enumerate(self.organisms) if org in node_organisms]
                nodes_order[col] = node_name
                all_indexes.extend(new_ind)
                all_columns.extend([col]*len(new_ind))
                data.extend([1.0]*len(new_ind))

            mat_p_a     = csc_matrix((data, (all_indexes,all_columns)), shape = (self.nb_organisms,len(graph)), dtype='float')#, dtype=numpy.bool_
            path_groups = mc.get_clusters(mc.run_mcl(jaccard_similarities(mat_p_a,jaccard_similarity_th), inflation = inflation_mcl_path_groups))

            # pdb.set_trace()
            # mat_hamming = pandas.DataFrame(squareform(pdist(mat_p_a.values, metric='jaccard')), index = mat_p_a.index, columns = mat_p_a.index)
            # mat_similarity_hamming  = mat_hamming.apply(lambda row: [-x+1 for x in row])
            # mat_similarity_hamming[mat_similarity_hamming < hamming_similarity] = 0
            # numpy.fill_diagonal(mat_similarity_hamming.values, 0)
            # nodes_indexes = mat_similarity_hamming.index
            self.path_groups_vectors = {}
            self.path_vectors = {}
            #path_groups = mc.get_clusters(mc.run_mcl(csr_matrix(mat_similarity_hamming.values), inflation= inflation_mcl_path_groups))
            #bar = tqdm(path_groups,total=len(path_groups), unit = "path groups")
            for i, path_group_index in enumerate(path_groups):
                #bar.set_description("Extracting path groups "+str(i))
                #bar.refresh()
                path_group = [nodes_order[ind] for ind in path_group_index]
                self.path_groups_vectors[subpartition+"#"+str(i)]=numpy.asarray(mat_p_a[:,path_group_index].sum(axis=1)/len(path_group)).flatten()
                if len(path_group)>1:
                    subg        = nx.Graph.subgraph(graph,path_group).copy()
                    all_weights = list(nx.get_edge_attributes(subg,"weight").values())
                    subg.remove_edges_from([(u,v) for u,v,d in subg.edges(data=True) if d["weight"] < numpy.median(all_weights)*breaker_th])
                    subg.remove_edges_from([(u,v) for u,v,d in subg.edges(data=True) if subg.degree(u)>2 or subg.degree(v)>2])
                    for j, path in enumerate(nx.algorithms.components.connected_components(subg)):
                        path_index = [nodes_order.inv[n] for n in path]
                        self.path_vectors[subpartition+"#"+str(i)+"#"+str(j)]=numpy.asarray(mat_p_a[:,path_index].sum(axis=1)/len(path_index)).flatten()
                        for node in path:
                            self.neighbors_graph.nodes[node]["path"]=str(i)+"#"+str(j)
                            self.neighbors_graph.nodes[node]["path_group"] = str(i)
                            subg.node[node]["path"]=subpartition+"#"+str(i)+"#"+str(j)
                        path_graph = nx.Graph.subgraph(subg,path)
                        for u, v in path_graph.edges():
                            self.neighbors_graph[u][v]["path"]=subpartition+"#"+str(i)+"#"+str(j)
                            self.neighbors_graph[u][v]["path_group"] = subpartition+"#"+str(i)
                else:
                    node = path_group.pop()
                    self.neighbors_graph.nodes[node]["path"]=subpartition+"#"+str(i)+"#1"
                    self.neighbors_graph.nodes[node]["path_group"] = subpartition+"#"+str(i)
                    self.path_vectors[subpartition+"#"+str(i)+"#1"]=numpy.asarray(mat_p_a[:,nodes_order.inv[node]].todense()).flatten()
            all_paths = list(set(nx.get_node_attributes(self.neighbors_graph,"path").values()))
            needed_col = colors(len(all_paths), except_list = [COLORS_RGB["persistent"],COLORS_RGB["shell"],COLORS_RGB["cloud"]])
            all_paths_colors = dict(zip(all_paths,needed_col))
            for n, d in self.neighbors_graph.nodes(data=True):
                if "path" in d:
                    self.neighbors_graph.node[n]["viz"]['color']=all_paths_colors[d["path"]]
            
            # for u, v, d in self.neighbors_graph.edges(data=True):
            #     if "path" in d:
            #         try:
            #             self.neighbors_graph[u][v]["viz"]['color']=all_paths_colors[d["path"]]
            #         except:
            #             self.neighbors_graph[u][v]["viz"]={'color':all_paths_colors[d["path"]]}
        if self.is_partitionned:
            return(self.path_groups_vectors,self.path_vectors)
        else:
            logging.getLogger().error("The pangenome is not partitionned")

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

def run_partitioning(nem_dir_path, nb_org, beta, free_dispersion, Q = 3, seed = 42, init="param_file", itermax=100, just_log_likelihood=False):
    logging.getLogger().debug("Running NEM...")
    if (Q<3 and not just_log_likelihood) or Q<2:
        logging.getLogger().error("Bad usage, Q must be higher or equal to 3 except for testing just_log_likelihood of a 2 paritions model")

    if init=="param_file":
        with open(nem_dir_path+"/nem_file_init_"+str(Q)+".m", "w") as m_file:
            m_file.write("1 ")# 1 to initialize parameter, 
            try:           
                m_file.write(" ".join([str(round(1/float(Q),2)) for q in range(Q-1)])+" ")# 1/Q give the initial proportition to each class (the last proportion is automaticaly determined by substraction in nem)
            except TypeError:
                pdb.set_trace()
            mu=[]
            epsilon=[]
            # mu = ["1"]*len(select_organisms)+ ["0"]*len(select_organisms)# persistent binary vector and cloud binary vector
            # epsilon = [str(low_disp)]*2*len(select_organisms)# persistent dispersition vector and cloud dispersition vector

            # for q in range(1,Q-1):
            #     mu += numpy.random.choice(["0","1"],len(select_organisms)).tolist()  # shell binary vector (1 ou 0, whatever because dispersion will be of 0.5)
            #     epsilon += ["0.5"]*len(select_organisms) # shell dispersition vector (high)
            step = 0.5/(math.ceil(Q/2))
            for q in range(1,Q+1):
                if q <= Q/2:
                    mu += ["1"]*nb_org
                    epsilon += [str(step*q)]*nb_org
                else:
                    mu += ["0"]*nb_org
                    epsilon += [str(step*(Q-q+1))]*nb_org
            
            m_file.write(" ".join(mu)+" "+" ".join(epsilon))

    # weighted_degree = sum(list(self.neighbors_graph.degree(weight="weight")).values())/nx.number_of_edges(self.neighbors_graph)
    # logging.getLogger().debug("weighted_degree: "+str(weighted_degree))
    # logging.getLogger().debug("org/weighted_degree: "+str(self.nb_organisms/weighted_degree))    
    # weighted_degree = sum(self.neighbors_graph.degree(weight="weight").values())/nx.number_of_edges(self.neighbors_graph)

    ALGO           = b"nem" #fuzzy classification by mean field approximation
    MODEL          = b"bern" # multivariate Bernoulli mixture model
    PROPORTION     = b"pk" #equal proportion :  "p_"     varying proportion : "pk"
    VARIANCE_MODEL = b"skd" if free_dispersion else b"sk_"#one variance per partition and organism : "sdk"      one variance per partition, same in all organisms : "sd_"   one variance per organism, same in all partion : "s_d"    same variance in organisms and partitions : "s__" 
    #NEIGHBOUR_SPEC = "f"# "f" specify to use all neighbors, orther argument is "4" to specify to use only the 4 neighbors with the higher weight (4 because for historic reason due to the 4 pixel neighbors of each pixel)
    CONVERGENCE    = b"clas"
    CONVERGENCE_TH = 0.01
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
    nem(Fname           = nem_dir_path.encode('ascii')+b"/nem_file",
        nk              = Q,
        algo            = ALGO,
        beta            = beta,
        convergence     = CONVERGENCE,
        convergence_th  = CONVERGENCE_TH,
        format          = b"fuzzy",
        it_max          = itermax,
        dolog           = True,
        model_family    = MODEL,
        proportion      = PROPORTION,
        dispersion      = VARIANCE_MODEL,
        init_mode       = INIT_PARAM_FILE if init in ["param_file","init_from_old"] else INIT_RANDOM,
        init_file       = nem_dir_path.encode('ascii')+b"/nem_file_init_"+str(Q).encode('ascii')+b".m",
        out_file_prefix = nem_dir_path.encode('ascii')+b"/nem_file_"+str(Q).encode('ascii'),
        seed            = seed)
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
    
    if os.path.isfile(nem_dir_path+"/nem_file_"+str(Q)+".uf"):
        logging.getLogger().debug("Reading NEM results...")
    elif not just_log_likelihood:
        logging.getLogger().warning("No NEM output file found: "+ nem_dir_path+"/nem_file_"+str(Q)+".uf")
    else:
        logging.getLogger().debug("No NEM output file found: "+ nem_dir_path+"/nem_file_"+str(Q)+".uf")
    index_fam = []
    all_parameters = {}
    # if just_log_likelihood:
    #     try:
    #         with open(nem_dir_path+"/nem_file_"+str(Q)+".mf","r") as parameter_nem_file:
    #             parameter = parameter_nem_file.readlines()
    #             log_likelihood = float(parameter[2].split()[3]) # M
    #             return(tuple([Q,log_likelihood]))
    #     except:
    #         return tuple([Q,None])
    
    with open(nem_dir_path+"/nem_file.index","r") as index_nem_file:
        for line in index_nem_file:
            index_fam.append(line.split("\t")[1].strip())

    partitions_list = ["U"] * len(index_fam)
    all_parameters  = {}
    U,D,L,log_likelihood,Z,error  = [None]*6
    entropy         = None
    try:
        with open(nem_dir_path+"/nem_file_"+str(Q)+".uf","r") as partitions_nem_file, open(nem_dir_path+"/nem_file_"+str(Q)+".mf","r") as parameters_nem_file:
            parameters      = parameters_nem_file.readlines()
            U,D,L,log_likelihood,Z,error = [float(p) for p in parameters[2].split()]
            print("U="+str(U))
            print("D="+str(D))
            print("log_likelihood="+str(log_likelihood))
            sum_mu_k       = []
            sum_epsilon_k  = []

            for k, line in enumerate(parameters[-Q:]):
                logging.getLogger().debug(line)
                vector = line.split() 
                mu_k = [bool(float(mu_kj)) for mu_kj in vector[0:nb_org]]
                logging.getLogger().debug(mu_k)
                logging.getLogger().debug(len(mu_k))
                epsilon_k = [round(float(epsilon_kj),2) for epsilon_kj in vector[nb_org+1:]]
                logging.getLogger().debug(epsilon_k)
                logging.getLogger().debug(len(epsilon_k))
                proportion = round(float(vector[nb_org]),2)
                logging.getLogger().debug(proportion)
                sum_mu_k.append(sum(mu_k))
                logging.getLogger().debug(sum(mu_k))
                sum_epsilon_k.append(sum(epsilon_k))
                logging.getLogger().debug(sum(epsilon_k))
                if k == 0:
                    all_parameters["persistent"]=(mu_k,epsilon_k,proportion)
                elif k == Q-1:
                    all_parameters["cloud"]=(mu_k,epsilon_k,proportion)
                else:
                    all_parameters["shell_"+str(k)]=(mu_k,epsilon_k,proportion)
                # for k2, (mu_k_2,epsilon_k_2,proportion_2) in all_parameters.items():
                #     for i in range(nb_org):
                #         inertia += (0.5-epsilon_k_2[i])(0.5-epsilon_k[i])(abs(mu_k_2[i]-mu_k[i]))
                
            # #persistent is defined by a sum of mu near of nb_organism and a low sum of epsilon
            # max_mu_k     = max(sum_mu_k)
            # persistent_k = sum_mu_k.index(max_mu_k)

            # #shell is defined by an higher sum of epsilon_k
            # max_epsilon_k = max(sum_epsilon_k)
            # shell_k       = sum_epsilon_k.index(max_epsilon_k)

            # # the other one should be cloud (basicaly with low sum_mu_k and low epsilon_k)
            # cloud_k = set([0,1,2]) - set([persistent_k, shell_k])
            # cloud_k = list(cloud_k)[0]

            # # but if the difference between epsilon_k of shell and cloud is tiny, we check using the sum of mu_k which basicaly must be lower in cloud
            # if ((sum_epsilon_k[shell_k]-sum_epsilon_k[cloud_k])/nb_org)<0.1 and sum_mu_k[shell_k]<sum_mu_k[cloud_k]:
            #     # otherwise we permutate
            #     (shell_k, cloud_k) = (cloud_k, shell_k)

            # logging.getLogger().debug(sum_mu_k)
            # logging.getLogger().debug(sum_epsilon_k)

            # logging.getLogger().debug(persistent_k)
            # logging.getLogger().debug(shell_k)
            # logging.getLogger().debug(cloud_k)

            partition               = {}
            # partition[persistent_k] = "P"#PERSISTENT
            # partition[shell_k]      = "S"#SHELL
            # partition[cloud_k]      = "C"#CLOUD
            # if partition[0] != "P" or partition[1] != "S" or partition[2] != "C":
            #     raise ValueError("vector mu_k and epsilon_k value in the mf file are not consistent with the initialisation value in the .m file")
            partition[0]   = "P"#PERSISTENT
            partition[Q-1] = "C"#CLOUD
            for i in range(1,Q-1):
                partition[i]="S"+str(i)
            entropy = 0

            for i, line in enumerate(partitions_nem_file):
                elements = [float(el) for el in line.split()]
                if just_log_likelihood:
                    entropy +=sum([math.log(float(el))*float(el) if float(el)>0 else 0 for el in elements])
                else:
                    max_prob = max([float(el) for el in elements])
                    positions_max_prob = [pos for pos, prob in enumerate(elements) if prob == max_prob]
                    logging.getLogger().debug(positions_max_prob)
                    logging.getLogger().debug(i)
                    if (len(positions_max_prob)>1 or max_prob<0.5):
                        partitions_list[i]="S_"#SHELL in case of doubt gene families is attributed to shell
                    else:
                        partitions_list[i]=partition[positions_max_prob.pop()]
            #logging.getLogger().debug(index.keys())
    except IOError:
        if not just_log_likelihood:
            logging.getLogger().warning("Statistical partitioning do not works (the number of organisms used is probably too low), see logs here to obtain more details "+nem_dir_path+"/nem_file.log")
        else:
            logging.getLogger().debug("Statistical partitioning do not works (the number of organisms used is probably too low), see logs here to obtain more details "+nem_dir_path+"/nem_file.log")
    except ValueError:
        ## return the default partitions_list which correspond to undefined
        pass
    if just_log_likelihood:
        return (tuple([Q,log_likelihood,entropy]))
    else:
        return((dict(zip(index_fam, partitions_list)),all_parameters,log_likelihood))

################ END OF FILE ################
