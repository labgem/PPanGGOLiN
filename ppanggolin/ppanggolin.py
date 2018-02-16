#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-

from collections import defaultdict, OrderedDict, Counter
from ordered_set import OrderedSet
import networkx as nx
import os
import sys
import math
import logging
import shutil
import gzip
import numpy as np
#import community
import tempfile
import subprocess
from tqdm import tqdm
import mmap
from random import sample
import time
from multiprocessing.pool import ThreadPool
from multiprocessing import Semaphore
from highcharts import Highchart
import contextlib

#import forceatlas2 

NEM_LOCATION  = os.path.dirname(os.path.abspath(__file__))+"/nem_exe"
(TYPE, FAMILY, START, END, STRAND, NAME, PRODUCT) = range(0, 7)#data index in annotation
(ORGANISM_ID, ORGANISM_GFF_FILE) = range(0, 2)#data index in the file listing organisms 

(GFF_seqname, GFF_source, GFF_feature, GFF_start, GFF_end, GFF_score, GFF_strand, GFF_frame, GFF_attribute) = range(0,9) 

RESERVED_WORDS = set(["id", "label", "name", "weight", "partition", "partition_exact", "length", "length_min", "length_max", "length_avg", "length_med", "product", 'nb_gene', 'community'])

SHORT_TO_LONG = {'A':'accessory','CE':'core_exact','P':'persistent','S':'shell','C':'cloud','U':'undefined'}

COLORS = {"pangenome":"black", "accessory":"#EB37ED", "core_exact" :"#FF2828", "shell": "#00D860", "persistent":"#F7A507", "cloud":"#79DEFF"}

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

                a networkx undirected graph. Node correspond to gene families and edges to chromosomal colocalization beween families. Organisms supporting each edge are stored in edge attribute as weel as the edge weight (number of organism coverinf each edge).
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
                    * partitions["core_exact"] contains the list of core families (core exact)
                    * partitions["accessory"] contains the list of families not in the core exact
                    * partitions["persistent"] contains the list of persistent families
                    * partitions["shell"] contains the list of shell families
                    * partitions["cloud"] contains the list of cloud families
                    * partitions["undefined"] contains the list of families unable to be classified (probably because the number of organisms is too small)

            .. attribute:: BIC

                a float providing the Bayesian Information Criterion. This Criterion give an estimation of the quality of the partionning (a low value means a good one)
                . seealso:: https://en.wikipedia.org/wiki/Bayesian_information_criterion
    """ 
    def __init__(self, init_from = "args", *args):
        """ 
            :param init_from: specified the excepted input (can be "file", "args", "database")
            :param *args: depending on the previous paramter, args can take multiple forms
            :type init_from: str
            :type *args: list

            :Example:

            >>>pan = PPanGGOLiN("file", organisms, gene_families, remove_high_copy_number_families)
            >>>pan = PPanGGOLiN("args", annotations, organisms, circular_contig_size, families_repeted)# load direclty the main attributes
        """ 
        self.undirected               = False
        self.annotations              = dict()
        self.neighbors_graph          = None
        self.organisms                = OrderedSet()
        self.nb_organisms             = 0
        self.circular_contig_size     = dict()
        self.families_repeted_th      = 0
        self.families_repeted         = set()
        self.pan_size                 = 0
        self.is_partitionned          = False
        self.nem_intermediate_files   = None
        self.partitions               = {}
        self.partitions["undefined"]  = list()
        self.partitions["persistent"] = list()
        self.partitions["shell"]      = list()
        self.partitions["cloud"]      = list()
        self.partitions["core_exact"] = list()
        self.partitions["accessory"]  = list()
        self.BIC                      = None # Bayesian Index Criterion
        #self.partitions_by_organisms  = defaultdict(lambda: defaultdict(set))

        if init_from == "file":
            self.__initialize_from_files(*args)
        elif init_from == "args":
            (self.annotations,
             self.organisms,
             self.circular_contig_size,
             self.families_repeted,
             self.undirected) = args 
        elif init_from == "database":
            logging.getLogger().error("database is not yet implemented")
            pass
        else:
            raise ValueError("init_from parameter is required")
        self.nb_organisms = len(self.organisms)

        logging.getLogger().info("Computing gene neighborhood ...")
        self.__neighborhood_computation(undirected = self.undirected)

    def __initialize_from_files(self, organisms_file, families_tsv_file, lim_occurence = 0, infer_singletons = False, undirected = False):
        """ 
            :param organisms_file: a file listing organims by compute, first column is organism name, second is path to gff file and optionnally other other to provide the name of circular contig
            :param families_tsv_file: a file listing families. The first element is the family identifier (by convention, we advice to use the identifier of the average gene of the family) and then the next elements are the identifiers of the genes belonging to this family.
            :param lim_occurence: a int containing the threshold of the maximum number copy of each families. Families exceeding this threshold are removed and are listed in the families_repeted attribute.
            :param infer_singletons: a bool specifying if singleton must be explicitely present in the families_tsv_file (False) or if single gene in gff files must be automatically infered as a singleton family (True)
            :param undirected: a bool specifying if the pangenome graph is undirected or directed
            :type file: 
            :type file: 
            :type int: 
            :type bool: 
            :type bool: 
        """ 
        self.undirected = undirected

        logging.getLogger().info("Reading "+families_tsv_file.name+" families file ...")
        families    = dict()
        first_iter  = True
        for line in families_tsv_file:
            elements = [el.strip() for el in line.split()]
            for gene in elements[1:]:
                families[gene]=elements[0]

        self.circular_contig_size = {}

        logging.getLogger().info("Reading "+organisms_file.name+" list of organism files ...")

        def get_num_lines(file_path):
            fp = open(file_path, "r+")
            buf = mmap.mmap(fp.fileno(), 0)
            lines = 0
            while buf.readline():
                lines += 1
            return lines

        bar = tqdm(organisms_file,total=get_num_lines(organisms_file.name), unit = "gff file")

        for line in bar:
            elements = [el.strip() for el in line.split("\t")]
            bar.set_description("Processing "+elements[ORGANISM_GFF_FILE])
            bar.refresh()
            if len(elements)>2:
                self.circular_contig_size.update({contig_id: None for contig_id in elements[2:len(elements)]})# size of the circular contig is initialized to None (waiting to read the gff files to fill the dictionnaries with the correct values)
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

            with open(gff_file_path,'r') as gff_file:
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
                annot[seq_id] = OrderedDict(sorted(annot[seq_id].items(), key= lambda item: item[1][START]))
                    
            if (lim_occurence > 0):
                fam_to_remove =[fam for fam, occ in cpt_fam_occ.items() if occ > lim_occurence]
                logging.getLogger().debug("highly repeted families found (>"+str(lim_occurence)+" in "+organism+"): "+" ".join(fam_to_remove))
                self.families_repeted = self.families_repeted.union(set(fam_to_remove))

            return(annot)
        else:
            raise KeyError("Redondant organism names was found ("+organism+")")

    def __str__(self):
        """ Return an overview of the statistics of the pangenome as a formated string """ 
        pan_str ="\n"
        pan_str += "----------- Statistics -----------\n"
        pan_str += "Number of organisms: "+str(self.nb_organisms)+"\n"

        if self.pan_size != 0:
            pan_str += "Pan-genome size:"+str(self.pan_size)+"\n"
            pan_str += "Exact core-genome size:"+str(len(self.partitions["core_exact"]))+"\n"
            pan_str += "Exact variable-genome size:"+str(self.pan_size-len(self.partitions["core_exact"]))+"\n"
            pan_str += "Persistent genome size:"+str(len(self.partitions["persistent"]))+"\n"
            pan_str += "Shell genome size:"+str(len(self.partitions["shell"]))+"\n"
            pan_str += "Cloud genome cloud:"+str(len(self.partitions["cloud"]))+"\n"
        else:
            pan_str += "No partitioning have been performed on this Pangenome instance\n"
            pan_str += "Run the partitioning function to obtain more detailled statistics...\n"
        pan_str += "----------------------------------"

        return(pan_str)    

    def __add_gene(self, fam_id, org, gene, name, length, product):#, multi_copy = None
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
        self.neighbors_graph.add_node(fam_id)

        try: 
            self.neighbors_graph.node[fam_id]["nb_gene"]+=1
        except KeyError:
            self.neighbors_graph.node[fam_id]["nb_gene"]=1
        try:
            self.neighbors_graph.node[fam_id][org].add(gene)
        except KeyError:
            self.neighbors_graph.node[fam_id][org] = set([gene])

        for attribute in ["name","length","product"]:
            try:
                self.neighbors_graph.node[fam_id][attribute].add(locals()[attribute])
            except KeyError:
                self.neighbors_graph.node[fam_id][attribute]=set([locals()[attribute]])

    def __add_link(self, fam_id, fam_id_nei, org, length):
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

        if not self.neighbors_graph.has_edge(fam_id,fam_id_nei):
            self.neighbors_graph.add_edge(fam_id, fam_id_nei)
            # logging.getLogger().debug([str(i) for i in [fam_id, fam_id_nei, org]])
        try:
            self.neighbors_graph[fam_id][fam_id_nei][org]+=1
        except KeyError:
            self.neighbors_graph[fam_id][fam_id_nei][org]=1
            try:
                self.neighbors_graph[fam_id][fam_id_nei]["weight"]+=1.0
            except KeyError:
                self.neighbors_graph[fam_id][fam_id_nei]["weight"]=1.0
        try:
            self.neighbors_graph[fam_id][fam_id_nei]["length"].add(length)
        except KeyError:
            self.neighbors_graph[fam_id][fam_id_nei]["length"]=set([length])

    def __neighborhood_computation(self, undirected = False):#,light = False, untangle_multi_copy_families = False
        """ Use the information already loaded (annotation) to build the pangenome graph
            :param undirected: a bool specifying if the graph is directed or undirected
            :type bool: 
        """ 
        #:param light: a bool specifying is the annotation attribute must be detroyed at each step to save memory
        if self.neighbors_graph is None:
            if undirected:
                self.neighbors_graph = nx.Graph()
            else:
                self.neighbors_graph = nx.DiGraph()

        bar = tqdm(list(self.annotations),total=len(self.annotations),unit = "organism")
        for organism in bar:

            bar.set_description("Processing "+organism)
            bar.refresh()

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

                family_id_nei, end_family_nei  = gene_info_start[FAMILY], gene_info_start[END]
                logging.getLogger().debug(gene_info_start)
                for gene, gene_info in contig_annot.items():
                    logging.getLogger().debug(gene_info)
                    logging.getLogger().debug(gene)
                    if gene_info[FAMILY] not in self.families_repeted:
                        self.__add_gene(gene_info[FAMILY],
                                        organism,
                                        gene,
                                        gene_info[NAME],
                                        gene_info[END]-gene_info[START],
                                        gene_info[PRODUCT])
                        self.neighbors_graph.add_node(family_id_nei)
                        self.__add_link(gene_info[FAMILY],family_id_nei,organism, gene_info[START] - end_family_nei)
                        family_id_nei  = gene_info[FAMILY]
                        end_family_nei = gene_info[END]
                
                if contig in self.circular_contig_size:#circularization
                    self.__add_link(gene_info_start[FAMILY],family_id_nei,organism, (self.circular_contig_size[contig] - end_family_nei) + gene_info_start[START])

                contig_annot[gene_start]=gene_info_start

            # if light:
            #     del self.annotations[organism]

        self.pan_size = nx.number_of_nodes(self.neighbors_graph)

    def partition(self, nem_dir_path    = tempfile.mkdtemp(),
                        organisms       = None,
                        beta            = 0.5,
                        free_dispersion = False,
                        chunck_size     = 200,
                        inplace         = True,
                        just_stats      = False,
                        nb_threads      = 1):
        print(nb_threads)
        """
            Use the graph topology and the presence or absence of genes from each organism into families to partition the pangenome in three groups ('persistent', 'shell' and 'cloud')
            . seealso:: Read the Mo Dang's thesis to understand NEM, a summary is available here : http://www.kybernetika.cz/content/1998/4/393/paper.pdf
            :param mode : a str specyfying if the pangenome will be parition in 2 partitions ("persistent" and "accessory") or 3 partitions ("persistent", "shell" and "accessory")
            :param nem_dir_path: a str containing a path to store temporary file of the NEM program
            :param beta: a float containing the spatial coefficient of smoothing of the clustering results using the weighted graph topology (0.00 turn off the spatial clustering)
            :param free_dispersion: a bool specyfing if the dispersion around the centroid vector of each paritition is the same for all the organisms or if the dispersion is free
            :param organisms: a list of organism to used to obtain the partition (must be included in the organism attributes of the object) or None to used all organisms in the object
            :param inplace: a boolean specifying if the partition must be stored in the object of returned (throw an error if inplace is true and organisms parameter i not None)
            :param chunck_size: an int specifying the size of the chunks
            :type str: 
            :type float: 
            :type bool: 
            :type list: 

            .. warning:: please use the function neighborhoodComputation before
        """ 
        
        if organisms is None:
            organisms = self.organisms
        else:
            organisms = OrderedSet(organisms)
            if len(organisms - self.organisms)>0:
                raise Exception("organisms parameter must be included in the organisms attribute of the objet")
            if inplace:
                raise Exception("inplace can't be true if the organisms parameter has not the same size than organisms attribute")

        # if self.neighbors_graph is None:
        #     raise Exception("The neighbors_graph is not built, please use the function neighborhood_computation before")
        if self.is_partitionned and inplace:
            logging.getLogger().warning("The pangenome was already partionned, inplace=true parameter will erase previous nem file intermediate files")
            self.delete_nem_intermediate_files()

        if not os.path.exists(nem_dir_path):
            #NEM requires 5 files: nem_file.index, nem_file.str, nem_file.dat, nem_file.m and nem_file.nei
            os.makedirs(nem_dir_path)
        if inplace:
            self.nem_intermediate_files = nem_dir_path

        stats = defaultdict(int)
        classification = []
        
        #core exact first
        families = []
        for node_name, data_organisms in self.neighbors_graph.nodes(data=True):
            compressed_vector = set([True if org in data_organisms else False for org in organisms])
            if len(compressed_vector)>1:
                families.append(node_name)
                stats["accessory"]+=1
            elif True in compressed_vector:# if size = 1 and contains just True, then core_exact
                families.append(node_name)
                stats["core_exact"]+=1

        BIC = 0

        if len(organisms) > chunck_size:

            cpt_partition = OrderedDict()
            for fam in families:
                cpt_partition[fam]= {"P":0,"S":0,"C":0,"U":0}
            
            total_BIC = 0
            #with contextlib.closing(ThreadPool(nb_threads)) as pool:
            sem = Semaphore(nb_threads)

            validated = set()
            cpt=0

            proba_sample = OrderedDict(zip(organisms,[len(organisms)]*len(organisms)))

            pan_size = stats["accessory"]+stats["core_exact"]
            if inplace:
                bar = tqdm(total = stats["accessory"]+stats["core_exact"], unit = "families partitionned")

            def validate_family(result):                    
                #nonlocal total_BIC
                try :
                    (BIC, partitions) = result

                    #total_BIC += BIC
                    for node,nem_class in partitions.items():
                        cpt_partition[node][nem_class]+=1
                        sum_partionning = sum(cpt_partition[node].values()) 

                        if sum_partionning > len(organisms)/chunck_size and max(cpt_partition[node].values()) > sum_partionning*0.5:
                            if node not in validated:
                                if inplace:
                                    bar.update()
                                validated.add(node)
                                # if max(cpt_partition[node], key=cpt_partition[node].get) == "P" and cpt_partition[node]["S"]==0 and cpt_partition[node]["C"]==0:
                                #     validated[node]="P"
                                # elif cpt_partition[node]["S"]==0:
                                #     validated[node]="C"
                                # else:
                                #     validated[node]="S" 
                finally:
                    sem.release()
            
            while len(validated)<pan_size:
                if sem.acquire():
                    # print(organisms)
                    # print(chunck_size)
                    # print(proba_sample.values())
                    # print(len(proba_sample.values()))
                    # min_o = min(proba_sample.values()) 
                    # max_o = max(proba_sample.values()) 
                    # range_o = max_o-min_o
                    # if min_o != max_o:
                    #     p = [(p-min_o/range_o)/len(organisms) for p in proba_sample.values()]
                    # else:
                    #     p = list(proba_sample.values())
                    # print(p)
                    #s = sum(proba_sample.values())
                    
                    #orgs = np.random.choice(organisms, size = chunck_size, replace = False, p = [p/s for p in proba_sample.values()])#
                    orgs = sample(organisms,chunck_size)
                    orgs = OrderedSet(orgs)
                    for org, p in proba_sample.items():
                        if org in orgs:
                            proba_sample[org] = p - len(organisms)/chunck_size if p >1 else 1
                        else:
                            proba_sample[org] = p + len(organisms)/chunck_size
                    # res = pool.apply_async(self.partition,
                    #                        args = (nem_dir_path+"/"+str(cpt)+"/",#nem_dir_path
                    #                                orgs,#organisms
                    #                                beta,#beta
                    #                                free_dispersion,#free dispersion
                    #                                chunck_size,#chunck_size
                    #                                False,#inplace
                    #                                False,#just_stats
                    #                                1),#nb_threads
                    #                        callback = validate_family)

                    res = self.partition(nem_dir_path+"/"+str(cpt)+"/",#nem_dir_path
                                                   orgs,#organisms
                                                   beta,#beta
                                                   free_dispersion,#free dispersion
                                                   chunck_size,#chunck_size
                                                   False,#inplace
                                                   False,#just_stats
                                                   1)
                    validate_family(res)

                    cpt +=1
                else:
                    time.sleep(0.01)

                if inplace:
                    bar.update()
                #pool.terminate()                

                #BIC = total_BIC/cpt
                BIC = 0
            classification = list()

            # if just_stats:
            #     print('len(validated)= '+str(len(validated)))
            #     print('len(cpt_partition)= '+str(len(cpt_partition)))

            for fam, data in cpt_partition.items():
                classification.append(max(data, key=data.get))

            # if just_stats:
            #     print("stat")
            #     #print(classification)
            #     c = Counter(classification)      
            #     print('stats["P"] '+str(c["P"]))
            #     print('stats["S"] '+str(c["S"]))
            #     print('stats["C"] '+str(c["C"])) 
            #     print('stats["U"] '+str(c["U"]))           
            #     print('total '+str(c["P"]+c["S"]+c["C"]+c["U"]))


            #     print('stats["accessory"] '+str(stats["accessory"]))
            #     print('stats["core_exact"] '+str(stats["core_exact"]))
            #     print('total '+str(stats["accessory"]+stats["core_exact"]))
            #     print(' ')
        else:
            classification = ["U"] * (stats["core_exact"]+stats["accessory"])

            if len(organisms)<=10:# below 10 organisms a statistical computation do not make any sence
                logging.getLogger().warning("The number of organisms is too low ("+str(len(organisms))+" organisms used) to partition the pangenome graph in persistent, shell, cloud partition, traditional partitions only (Core and Accessory genome) will be provided")
            else:
                logging.getLogger().debug("Writing nem_file.str nem_file.index nem_file.nei nem_file.dat and nem_file.m files")
                with open(nem_dir_path+"/nem_file.str", "w") as str_file,\
                     open(nem_dir_path+"/nem_file.index", "w") as index_file,\
                     open(nem_dir_path+"/column_org_file", "w") as org_file,\
                     open(nem_dir_path+"/nem_file.nei", "w") as nei_file,\
                     open(nem_dir_path+"/nem_file.dat", "w") as dat_file,\
                     open(nem_dir_path+"/nem_file.m", "w") as m_file:

                    nei_file.write("1\n")
                    
                    org_file.write(" ".join([org for org in organisms])+"\n")
                    org_file.close()

                    
                    index_fam = OrderedDict()
                    for node_name, node_organisms in self.neighbors_graph.nodes(data=True):
                        logging.getLogger().debug(node_organisms)
                        logging.getLogger().debug(organisms)
                        
                        if not organisms.isdisjoint(node_organisms): # if at least one commun organism
                            dat_file.write("\t".join(["1" if org in node_organisms else "0" for org in organisms])+"\n")
                            index_fam[node_name] = len(index_fam)+1
                            index_file.write(str(len(index_fam))+"\t"+str(node_name)+"\n")
                    for node_name, index in index_fam.items():
                        row_fam         = []
                        row_dist_score  = []
                        neighbor_number = 0
                        try:
                            for neighbor in set(nx.all_neighbors(self.neighbors_graph, node_name)):
                                coverage = 0
                                if self.neighbors_graph.is_directed():
                                    cov_sens, cov_antisens = (0,0)
                                    try:
                                        cov_sens = sum([pre_abs for org, pre_abs in self.neighbors_graph[node_name][neighbor].items() if ((org in organisms) and (org not in RESERVED_WORDS))])
                                    except KeyError:
                                        pass
                                    try:
                                        cov_antisens = sum([pre_abs for org, pre_abs in self.neighbors_graph[neighbor][node_name].items() if ((org in organisms) and (org not in RESERVED_WORDS))])
                                    except KeyError:
                                        pass
                                    coverage = cov_sens + cov_antisens
                                else:
                                    coverage = sum([pre_abs for org, pre_abs in self.neighbors_graph[node_name][neighbor].items() if ((org in organisms) and (org not in RESERVED_WORDS))])

                                if coverage==0:
                                    continue
                                distance_score = coverage/self.nb_organisms
                                row_fam.append(str(index_fam[neighbor]))
                                row_dist_score.append(str(round(distance_score,4)))
                                neighbor_number += 1
                            if neighbor_number>0:
                                nei_file.write("\t".join([str(item) for sublist in [[index_fam[node_name]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
                            else:
                                nei_file.write(str(len(index_fam))+"\t0\n")
                                logging.getLogger().debug("The family: "+node_name+" is an isolated family in the selected organisms")
                        except nx.exception.NetworkXError as nxe:
                            logging.getLogger().debug("The family: "+node_name+" is an isolated family")
                            nei_file.write(str(len(index_fam))+"\t0\n")

                    m_file.write("1 0.33333 0.33333 ") # 1 to initialize parameter, 0.333 and 0.333 for to give one third of initial proportition to each class (last 0.33 is automaticaly determined by substraction)
                    m_file.write(" ".join(["1"]*len(organisms))+" ") # persistent binary vector
                    m_file.write(" ".join(["1"]*len(organisms))+" ") # shell binary vector (1 ou 0, whatever because dispersion will be of 0.5)
                    m_file.write(" ".join(["0"]*len(organisms))+" ") # cloud binary vector
                    m_file.write(" ".join(["0.1"]*len(organisms))+" ") # persistent dispersition vector (low)
                    m_file.write(" ".join(["0.5"]*len(organisms))+" ") # shell dispersition vector (high)
                    m_file.write(" ".join(["0.1"]*len(organisms))) # cloud dispersition vector (low)

                    str_file.write("S\t"+str(len(index_fam))+"\t"+
                                         str(len(organisms))+"\n")

                logging.getLogger().debug("Running NEM...")
                # weighted_degree = sum(list(self.neighbors_graph.degree(weight="weight")).values())/nx.number_of_edges(self.neighbors_graph)
                # logging.getLogger().debug("weighted_degree: "+str(weighted_degree))
                # logging.getLogger().debug("org/weighted_degree: "+str(self.nb_organisms/weighted_degree))    
                #weighted_degree = sum(self.neighbors_graph.degree(weight="weight").values())/nx.number_of_edges(self.neighbors_graph)

                Q              = 3 # number of partitions
                ALGO           = "ncem" #fuzzy classification by mean field approximation
                ITERMAX        = 100 # number of iteration max 
                MODEL          = "bern" # multivariate Bernoulli mixture model
                PROPORTION     = "pk" #equal proportion :  "p_"     varying proportion : "pk"
                VARIANCE_MODEL = "skd" if free_dispersion else "sk_"#one variance per partition and organism : "sdk"      one variance per partition, same in all organisms : "sd_"   one variance per organism, same in all partion : "s_d"    same variance in organisms and partitions : "s__" 
                NEIGHBOUR_SPEC = "f"# "f" specify to use all neighbors, orther argument is "4" to specify to use only the 4 neighbors with the higher weight (4 because for historic reason due to the 4 pixel neighbors of each pixel)
                CONVERGENCE_TH = "clas "+str(0.000001)

                # HEURISTIC      = "heu_d"# "psgrad" = pseudo-likelihood gradient ascent, "heu_d" = heuristic using drop of fuzzy within cluster inertia, "heu_l" = heuristic using drop of mixture likelihood
                # STEP_HEURISTIC = 0.5 # step of beta increase
                # BETA_MAX       = float(len(organisms)) #maximal value of beta to test,
                # DDROP          = 0.8 #threshold of allowed D drop (higher = less detection)
                # DLOSS          = 0.5 #threshold of allowed D loss (higher = less detection)
                # LLOSS          = 0.02 #threshold of allowed L loss (higher = less detection)
                
                # BETA           = ["-B",HEURISTIC,"-H",str(STEP_HEURISTIC),str(BETA_MAX),str(DDROP),str(DLOSS),str(LLOSS)] if beta == float('Inf') else ["-b "+str(beta)]

                WEIGHTED_BETA = beta*len(organisms)
                command = " ".join([NEM_LOCATION, 
                                    nem_dir_path+"/nem_file",
                                    str(Q),
                                    "-a", ALGO,
                                    "-i", str(ITERMAX),
                                    "-m", MODEL, PROPORTION, VARIANCE_MODEL,
                                    "-s m "+ nem_dir_path+"/nem_file.m",
                                    "-b "+str(WEIGHTED_BETA),
                                    "-n", NEIGHBOUR_SPEC,
                                    "-c", CONVERGENCE_TH,
                                    "-f fuzzy",
                                    "-l y"])
             
                logging.getLogger().debug(command)

                proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE if logging.getLogger().getEffectiveLevel() == logging.INFO else None,
                                                            stderr=subprocess.PIPE if logging.getLogger().getEffectiveLevel() == logging.INFO else None)
                (out,err) = proc.communicate()
                if logging.getLogger().getEffectiveLevel() == logging.INFO:
                    with open(nem_dir_path+"/out.txt", "wb") as file_out, open(nem_dir_path+"/err.txt", "wb") as file_err:
                        file_out.write(out)
                        file_err.write(err)

                logging.getLogger().debug(out)
                logging.getLogger().debug(err)

                if beta == float('Inf'):
                    starting_heuristic = False
                    with open(nem_dir_path+"/beta_evol.txt", "w") as beta_evol_file:
                        for line in str(err).split("\\n"):
                            if not starting_heuristic and line.startswith("* * Starting heuristic * *"):
                                beta_evol_file.write("beta\tLmix\n")
                                starting_heuristic = True
                                continue
                            if starting_heuristic:
                                elements = line.split("=")
                                if elements[0] == " * * Testing beta ":
                                   beta = float(elements[1].split("*")[0].strip())
                                elif elements[0] == "  criterion NEM ":
                                    Lmix = float(elements[len(elements)-1].strip())
                                    beta_evol_file.write(str(beta)+"\t"+str(Lmix)+"\n")

                if os.path.isfile(nem_dir_path+"/nem_file.uf"):
                    logging.getLogger().debug("Reading NEM results...")
                else:
                    logging.getLogger().error("No NEM output file found: "+ nem_dir_path+"/nem_file.uf")
                
                
                try:
                    with open(nem_dir_path+"/nem_file.uf","r") as classification_nem_file, open(nem_dir_path+"/nem_file.mf","r") as parameter_nem_file:
                        sum_mu_k = []
                        sum_epsilon_k = []
                        proportion = []

                        parameter = parameter_nem_file.readlines()
                        M = float(parameter[6].split()[3]) # M is markov ps-like
                        BIC = -2 * M - (Q * len(organisms) * 2 + Q - 1) * math.log(len(index_fam))
                        logging.getLogger().debug("The Bayesian Criterion Index of the partionning is "+str(BIC))

                        for k, line in enumerate(parameter[-3:]):
                            logging.getLogger().debug(line)
                            vector = line.split() 
                            mu_k = [bool(float(mu_kj)) for mu_kj in vector[0:len(organisms)]]
                            logging.getLogger().debug(mu_k)
                            logging.getLogger().debug(len(mu_k))

                            epsilon_k = [float(epsilon_kj) for epsilon_kj in vector[len(organisms)+1:]]
                            logging.getLogger().debug(epsilon_k)
                            logging.getLogger().debug(len(epsilon_k))
                            proportion = float(vector[len(organisms)])
                            logging.getLogger().debug(proportion)
                            sum_mu_k.append(sum(mu_k))
                            logging.getLogger().debug(sum(mu_k))
                            sum_epsilon_k.append(sum(epsilon_k))
                            logging.getLogger().debug(sum(epsilon_k))


                        
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
                        # if (sum_epsilon_k[shell_k]-sum_epsilon_k[cloud_k])<0.05 and sum_mu_k[shell_k]<sum_mu_k[cloud_k]:
                        #     # otherwise we permutate
                        #     (shell_k, cloud_k) = (cloud_k, shell_k)

                        logging.getLogger().debug(sum_mu_k)
                        logging.getLogger().debug(sum_epsilon_k)

                        logging.getLogger().debug(persistent_k)
                        logging.getLogger().debug(shell_k)
                        logging.getLogger().debug(cloud_k)

                        partition                 = {}
                        partition[persistent_k] = "P"#PERSISTENT
                        partition[shell_k]      = "S"#SHELL
                        partition[cloud_k]      = "C"#CLOUD

                        for i, line in enumerate(classification_nem_file):
                            elements = [float(el) for el in line.split()]
                            max_prob = max([float(el) for el in elements])
                            positions_max_prob = [pos for pos, prob in enumerate(elements) if prob == max_prob]
                            logging.getLogger().debug(positions_max_prob)
                            logging.getLogger().debug(i)

                            if (len(positions_max_prob)>1):
                                classification[i]="S"#SHELL in case of doubt (equiprobable partition), gene families is attributed to shell
                            else:
                                classification[i] = partition[positions_max_prob.pop()]

                        logging.getLogger().debug(partition)
                        #logging.getLogger().debug(index.keys())
                except FileNotFoundError:
                    logging.getLogger().warning("Statistical partitioning do not works, see log here to obtain more details "+nem_dir_path+"/nem_file.log")
        if inplace:

            self.BIC = BIC

            if self.is_partitionned:
                for p in SHORT_TO_LONG.values():
                    self.partitions[p] = list()# erase older values

            for node, nem_class in zip(self.neighbors_graph.nodes(), classification):
                nb_orgs=0
                for key in list(self.neighbors_graph.node[node].keys()):
                    if key not in RESERVED_WORDS:
                        #self.partitions_by_organisms[key][partition[int(nem_class)]].add(self.neighbors_graph.node[node][key])
                        nb_orgs+=1

                self.neighbors_graph.node[node]["partition"]=SHORT_TO_LONG[nem_class]
                self.partitions[SHORT_TO_LONG[nem_class]].append(node)

                if nb_orgs == self.nb_organisms:
                    self.partitions["core_exact"].append(node)#CORE EXACT
                    self.neighbors_graph.node[node]["partition_exact"]="core_exact"
                elif nb_orgs < self.nb_organisms:
                    self.partitions["accessory"].append(node)#ACCESSORY
                    self.neighbors_graph.node[node]["partition_exact"]="accessory"
                else:
                    logging.getLogger().error("nb_orgs can't be > to self.nb_organisms")
                    exit(1)

            if len(self.families_repeted)>0:
                logging.getLogger().info("Gene families that have been discarded because there are repeated:\t"+" ".join(self.families_repeted))
            else:
                logging.getLogger().info("No gene families have been discarded because there are repeated")

            logging.getLogger().debug(nx.number_of_edges(self.neighbors_graph))

            self.is_partitionned=True
        else:
            if just_stats:
                for node_name, nem_class in zip(self.neighbors_graph.nodes(data=False), classification):
                    stats[SHORT_TO_LONG[nem_class]]+=1
                return stats
            else:
                return (BIC,dict(zip(index_fam.keys(), classification)))

    def export_to_GEXF(self, graph_output_path, compressed=False, metadata = None, all_node_attributes = True, all_edge_attributes = True):
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
        graph_to_save = self.neighbors_graph.copy()

        for node in self.neighbors_graph.nodes():            
            for key in list(self.neighbors_graph.node[node].keys()):
                if not all_node_attributes and key in self.organisms:
                    del graph_to_save.node[node][key]
                else:
                    try:
                        if not isinstance(self.neighbors_graph.node[node][key], str):
                            graph_to_save.node[node][key]="|".join(self.neighbors_graph.node[node][key])#because networkx and gephi do not support list type in gexf despite it is possible according to the specification using liststring (https://gephi.org/gexf/1.2draft/data.xsd)
                    except TypeError:
                        if key == "length":
                            l = list(self.neighbors_graph.node[node][key])
                            graph_to_save.node[node]["length_avg"] = float(np.mean(l))
                            graph_to_save.node[node]["length_med"] = float(np.median(l))
                            graph_to_save.node[node]["length_min"] = min(l)
                            graph_to_save.node[node]["length_max"] = max(l)
                            del graph_to_save.node[node]["length"]

        for node_i, node_j, data in self.neighbors_graph.edges(data = True):
            l = list(data["length"])
            graph_to_save[node_i][node_j]["length_avg"] = float(np.mean(l))
            graph_to_save[node_i][node_j]["length_med"] = float(np.median(l))
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

        graph_output_path = graph_output_path+".gexf"
        if compressed:
            graph_output_path = gzip.open(graph_output_path+".gz","w")

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
            Exported the pangenome as a csv_matrix similar to the csv et Rtab matrix exported by Roary (https://sanger-pathogens.github.io/Roary/)
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
                        genes  = [('"'+"|".join(data[org])+'"' if gene_or_not else "1") if org in data else ('""' if gene_or_not else "0") for org in self.organisms]
                        nb_org = len([gene for gene in genes if gene != ('""' if gene_or_not else "0")])
                        l = list(data["length"])
                        matrix.write(sep.join(['"'+node+'"',#1
                                               '"'+data["partition"]+'"',#2
                                               '"'+"|".join(data["product"])+'"',#3
                                               str(nb_org),#4
                                               str(data["nb_gene"]),#5
                                               str(round(data["nb_gene"]/nb_org,2)),#6
                                               '""',#7
                                               '""',#8
                                               '""',#9
                                               '""',#10
                                               '""',#11
                                               str(min(l)),#12
                                               str(max(l)),#13
                                               str(round(np.mean(l),2))]#14
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
    #     self.partitions["core_exact"] = list()
    #     self.partitions["accessory"]  = list()
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

    def ushaped_plot(self, outdir):
        """
            generate ushaped representation
            :param outdir: a str containing the path of the output file
            :type str: 
        """ 
        ushaped_plot = Highchart(width = 1800, height = 800)

        count = defaultdict(lambda : defaultdict(int))

        for node, data in self.neighbors_graph.nodes(data=True):
            nb_org  = len([True for org in self.organisms if org in data])
            count[nb_org][data["partition"]]+=1

        persistent_values = []
        shell_values      = []
        cloud_values      = []

        for nb_org in range(1,self.nb_organisms+1):
            persistent_values.append(count[nb_org]["persistent"])
            shell_values.append(count[nb_org]["shell"])
            cloud_values.append(count[nb_org]["cloud"])
        
        options_ushaped_plot={
        'title': {'text':'Distribution of gene families frequency in the pangenome'},
        'xAxis': {'tickInterval': 1, 'categories': list(range(1,self.nb_organisms+1)), 'title':{'text':'# of organisms in which each family is present'}},
        'yAxis': {'allowDecimals': False, 'title' : {'text':'# of families'}},
        'tooltip': {'headerFormat': '<span style="font-size:11px"># of orgs: <b>{point.x}</b></span><br>',
                    'pointFormat': '<span style="color:{point.color}">{series.name}</span>: {point.y}<br/>',
                    'shared': True},
        'plotOptions': {'column': {'stacking': 'normal'}}
        }
        ushaped_plot.set_dict_options(options_ushaped_plot)
        ushaped_plot.add_data_set(persistent_values,'column','Persistent', color = COLORS["persistent"])
        ushaped_plot.add_data_set(shell_values,'column','Shell', color = COLORS["shell"])
        ushaped_plot.add_data_set(cloud_values,'column','Cloud', color = COLORS["cloud"])

        ushaped_plot.save_file(filename = outdir+"/ushaped_plot")

    # def tile_plot(self, outdir):
    #     """
    #         generate tile plot representation (not work for the moment)
    #         :param outdir: a str containing the path of the output file
    #         :type str: 
    #     """ 

    #     tile_plot = Highchart(width = 1600, height = 1280)

    #     binary_data = []
    #     fam_order = []
    #     cpt = 1
    #     for node, data in self.neighbors_graph.nodes(data=True):
    #         fam_order.append(node)
    #         v = [[1,2,3] if org in data else [0,0,0] for org in self.organisms]
    #         binary_data.append(v)
    #         cpt+=1
    #         if cpt>50:
    #             break
     
    #     print(binary_data)
    #     options_tile_plot={
    #     'chart': {'type': 'heatmap', 'plotBorderWidth': 1},
    #     'title': {'text':'Presence/Absence matrix'},
    #     'xAxis': {'categories': fam_order},
    #     'yAxis': {'categories': list(self.organisms)},
    #     'colorAxis': {'min': 0, 'max': 1, 'minColor': '#FFFFFF', 'maxColor': '#7CB5EC'}
    #     # ,
    #     # 'legend', {'align': 'right', 'layout': 'vertical', 'margin': 0, 'verticalAlign': 'top', 'y': 25, 'symbolHeight': 280}   
    #     }
    #     print(options_tile_plot)
    #     tile_plot.set_dict_options(options_tile_plot)
    #     tile_plot.add_data_set(binary_data)

    #     tile_plot.save_file(filename = outdir+"/tile_plot")

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


    def projection_polar_histogram(self, out_dir, organisms_to_project):
        """
            generate tile projection_polar_histogram representation
            :param outdir: a str containing the path of the output directory (name of files will be the name of organisms)
            :organisms_to_project: a list of str containing the name of the organism
            :type str:
            :type list
        """ 
        for organism in organisms_to_project:
            with open(out_dir+"/"+organism+".csv","w") as out_file:
                out_file.write("gene\tcontig\tori\tfamily\tpartition\tpersistent\tshell\tcloud\n")
                for contig, contig_annot in self.annotations[organism].items():
                    for gene, gene_info in contig_annot.items():
                        if gene_info[FAMILY] not in self.families_repeted:
                            nei_partitions = [self.neighbors_graph.node[nei]["partition"] for nei in nx.all_neighbors(self.neighbors_graph,gene_info[FAMILY])]
                            out_file.write("\t".join([gene,
                                                      contig,
                                                      "T" if (gene_info[NAME].upper() == "DNAA" or gene_info[PRODUCT].upper() == "DNAA") else "F",
                                                      gene_info[FAMILY],
                                                      self.neighbors_graph.node[gene_info[FAMILY]]["partition"],
                                                      str(nei_partitions.count("persistent")),
                                                      str(nei_partitions.count("shell")),
                                                      str(nei_partitions.count("cloud"))])+"\n")

################ END OF CLASS PPanGGOLiN ################
