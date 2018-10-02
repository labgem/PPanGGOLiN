#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-

import warnings
warnings.filterwarnings("ignore")
import numpy
import scipy.optimize as optimization
import pandas
from collections import defaultdict, OrderedDict
from ordered_set import OrderedSet
import networkx as nx
import logging
import sys
import os
import argparse
from random import shuffle, sample
from tqdm import tqdm
tqdm.monitor_interval = 0
from concurrent.futures import ProcessPoolExecutor, as_completed
from time import strftime, time, localtime
import subprocess
import pkg_resources
import traceback
import shutil
from .ppanggolin import *
from .utils import *

### PATH AND FILE NAME
OUTPUTDIR                   = None 
TMP_DIR                     = None
NEM_DIR                     = "/NEM_results/"
FIGURE_DIR                  = "/figures/"
PROJECTION_DIR              = "/projections/"
METADATA_DIR                = "/metadata/"
EVOLUTION_DIR               = "/evolutions/"
PARTITION_DIR               = "/partitions/"
PATH_DIR                    = "/paths/"
GRAPH_FILE_PREFIX           = "/graph"
MATRIX_FILES_PREFIX         = "/matrix"
USHAPE_PLOT_PREFIX          = "/Ushaped_plot"
MATRIX_PLOT_PREFIX          = "/presence_absence_matrix_plot"
EVOLUTION_CURVE_PREFIX      = "/evolution_curve"
EVOLUTION_STATS_FILE_PREFIX = "/evol_stats"
EVOLUTION_PARAM_FILE_PREFIX = "/evol_param"
SUMMARY_STATS_FILE_PREFIX   = "/summary_stats"
CORRELATED_PATHS_PREFIX     = "/correlated_paths_prefix"
SCRIPT_R_FIGURE             = "/generate_plots.R"

def plot_Rscript(script_outfile, verbose=True):
    """
    """
    rscript = "#!/usr/bin/env R\n"+("options(warn=-1)\n" if not verbose else "")+"""
options(show.error.locations = TRUE)

library("ggplot2")
library("reshape2")

color_chart = c(pangenome="black", "exact_accessory"="#EB37ED", "exact_core" ="#FF2828",  "soft_core"="#e6e600", "soft_accessory"="#996633", shell = "#00D860", persistent="#F7A507", cloud = "#79DEFF","soft_core"="#e6e600", "soft_accessory"="#996633","undefined"="#828282")

########################### START U SHAPED PLOT #################################

binary_matrix         <- read.table('"""+OUTPUTDIR+MATRIX_FILES_PREFIX+""".Rtab', header=TRUE, sep='\\t', check.names = FALSE)
data_header           <- c("Gene","Non-unique Gene name","Annotation","No. isolates","No. sequences","Avg sequences per isolate","Accessory Fragment","Genome Fragment","Order within Fragment","Accessory Order with Fragment","QC","Min group size nuc","Max group size nuc","Avg group size nuc") 
family_data           <- binary_matrix[,colnames(binary_matrix) %in% data_header]
binary_matrix         <- binary_matrix[,!(colnames(binary_matrix) %in% data_header)]
binary_matrix[binary_matrix>1] <- 1
occurences            <- rowSums(binary_matrix)
classification_vector <- family_data[,2]

c = data.frame(nb_org = occurences, partition = classification_vector)

plot <- ggplot(data = c) + 
    geom_bar(aes_string(x = "nb_org", fill = "partition")) +
    scale_fill_manual(name = "partition", values = color_chart, breaks=c("persistent","shell","cloud")) +
    scale_x_discrete(limits = seq(1, ncol(binary_matrix))) +
    xlab("# of organisms in which each familly is present")+
    ylab("# of families")+
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave('"""+OUTPUTDIR+FIGURE_DIR+USHAPE_PLOT_PREFIX+""".pdf', device = "pdf", height= (par("din")[2]*1.5),plot)

########################### END U SHAPED PLOT #################################

########################### START PRESENCE/ABSENCE MATRIX #################################

nb_org                  <- ncol(binary_matrix)

binary_matrix_hclust    <- hclust(dist(t(binary_matrix), method="binary"))
binary_matrix           <- data.frame(binary_matrix,"NEM partitions" = classification_vector, occurences = occurences, check.names=FALSE)

binary_matrix[occurences == nb_org, "Former partitions"] <- "exact_core"
binary_matrix[occurences != nb_org, "Former partitions"] <- "exact_accessory"
binary_matrix = binary_matrix[order(match(binary_matrix$"NEM partitions",c("persistent", "shell", "cloud")),
                                    match(binary_matrix$"Former partitions",c("exact_core", "exact_accessory")),
                                    -binary_matrix$occurences),
                              c(binary_matrix_hclust$label[binary_matrix_hclust$order],"NEM partitions","Former partitions")]

binary_matrix$familles <- seq(1,nrow(binary_matrix))
data = melt(binary_matrix, id.vars=c("familles"))

colnames(data) = c("fam","org","value")

data$value <- factor(data$value, levels = c(1,0,"persistent", "shell", "cloud", "exact_core", "exact_accessory"), labels = c("presence","absence","persistent", "shell", "cloud", "exact_core", "exact_accessory"))

plot <- ggplot(data = data)+
        geom_raster(aes_string(x="org",y="fam", fill="value"))+
        scale_fill_manual(values = c("presence"="green","absence"="grey80",color_chart)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=1), panel.border = element_blank(), panel.background = element_blank())

ggsave('"""+OUTPUTDIR+FIGURE_DIR+MATRIX_PLOT_PREFIX+""".pdf', device = "pdf", plot)

########################### END PRESENCE/ABSENCE MATRIX #################################

########################### START EVOLUTION CURVE #################################

library("ggrepel")
library("data.table")
library("minpack.lm")

if (file.exists('"""+OUTPUTDIR+EVOLUTION_DIR+EVOLUTION_STATS_FILE_PREFIX+""".csv')){
    data <- read.csv('"""+OUTPUTDIR+EVOLUTION_DIR+EVOLUTION_STATS_FILE_PREFIX+""".csv', header = TRUE)
    data <- melt(data, id = "nb_org")
    colnames(data) <- c("nb_org","partition","value")

    final_state = data[data$nb_org == max(data$nb_org,na.rm=T),]
    final_state = final_state[!duplicated(final_state[,c("nb_org","partition")]), ]
    final <- structure(names = as.character(final_state$partition), as.integer(final_state$value))

    #gamma and kappa are calculated according to the Tettelin et al. 2008 approach
    median_by_nb_org <- setDT(data)[,list(med=as.numeric(median(value))), by=c("nb_org","partition")]
    colnames(median_by_nb_org) <- c("nb_org_comb","partition","med")

    for (part in as.character(final_state$partition)){
        regression  <- nlsLM(med~kapa*(nb_org_comb^gama),median_by_nb_org[which(median_by_nb_org$partition == part),],start=list(kapa=1000,gama=0))
        coefficient <- coef(regression)
        final_state[final_state$partition == part,"formula" ] <- paste0("F == ", format(coefficient["kapa"],decimal.mark = ",",digits =2),"~N^{",format(coefficient["gama"],digits =2),"}")
    }

    plot <- ggplot(data = data, aes_string(x="nb_org",y="value", colour = "partition"))+
            ggtitle(bquote(list("Rarefaction curve. Heap's law parameters based on Tettelin et al. 2008 approach", n == kappa~N^gamma)))+
            geom_smooth(data        = median_by_nb_org[median_by_nb_org$partition %in% c("pangenome","shell","cloud","exact_accessory", "persistent", "exact_core") ,],# 
                        mapping     = aes_string(x="nb_org_comb",y="med",colour = "partition"),
                        method      = "nlsLM",
                        formula     = y~kapa*(x^gama),method.args =list(start=c(kapa=1000,gama=0)),
                        linetype    ="twodash",
                        size        = 1.5,
                        se          = FALSE,
                        show.legend = FALSE)+
            stat_summary(fun.ymin = function(z) { quantile(z,0.25) },  fun.ymax = function(z) { quantile(z,0.75) }, geom="ribbon", alpha=0.1,size=0.1, linetype="dashed", show.legend = FALSE)+
            stat_summary(fun.y=median, geom="line",size=0.5)+
            stat_summary(fun.y=median, geom="point",shape=4,size=1, show.legend = FALSE)+
            stat_summary(fun.ymax=max,fun.ymin=min,geom="errorbar",linetype="dotted",size=0.1,width=0.2)+
            scale_x_continuous(breaks = as.numeric(unique(data$nb_org)))+
            scale_y_continuous(limits=c(0,max(data$value,na.rm=T)), breaks = seq(0,max(data$value,na.rm=T),1000))+
            scale_colour_manual(name = "NEM partitioning", values = color_chart, breaks=names(sort(final, decreasing = TRUE)))+
            geom_label_repel(data = final_state, aes_string(x="nb_org", y="value", colour = "partition", label = "value"), show.legend = FALSE,
                      fontface = 'bold', fill = 'white',
                      box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.5, "lines"),
                      segment.color = 'grey50',
                      nudge_x = 45) +
            geom_label_repel(data = final_state, aes(x = nb_org*0.9, y = value, label = formula), size = 2, parse = TRUE, show.legend = FALSE, segment.color = NA) + 
            xlab("# of organisms")+
            ylab("# of families")+
            ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.minor = element_blank())
    
    ggsave('"""+OUTPUTDIR+FIGURE_DIR+EVOLUTION_CURVE_PREFIX+""".pdf', device = "pdf", width = (par("din")[1]*2) ,plot)

}
########################### END EVOLUTION CURVE #################################

########################### START PROJECTION #################################

for (org_csv in list.files(path = '"""+OUTPUTDIR+PROJECTION_DIR+"""', pattern = "*.csv$", full.names = T)){
    org_name <- tools::file_path_sans_ext(basename(org_csv))
    data <- read.csv(org_csv, header = T)
    if(org_name=="nb_genes"){
        data.melted = melt(data,id.var="org")
        colnames(data.melted) <- c("org","partition","nb_geness")
        is_outlier <- function(x) {
            return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
        }
        data.melted_persistent = data.melted[data.melted$partition =="persistent",]
        data.melted_persistent <- data.melted_persistent[is_outlier(data.melted_persistent$nb_geness),]
        data.melted_shell = data.melted[data.melted$partition =="shell",]
        data.melted_shell <- data.melted_shell[is_outlier(data.melted_shell$nb_geness),]
        data.melted_cloud = data.melted[data.melted$partition =="cloud",]
        data.melted_cloud <- data.melted_cloud[is_outlier(data.melted_cloud$nb_geness),]
        data.melted_exact_core = data.melted[data.melted$partition =="exact_core",]
        data.melted_exact_core <- data.melted_exact_core[is_outlier(data.melted_exact_core$nb_geness),]
        data.melted_exact_accessory = data.melted[data.melted$partition =="exact_accessory",]
        data.melted_exact_accessory <- data.melted_exact_accessory[is_outlier(data.melted_exact_accessory$nb_geness),]
        data.melted_pangenome = data.melted[data.melted$partition =="cloud",]
        data.melted_pangenome <- data.melted_pangenome[is_outlier(data.melted_pangenome$nb_geness),]

        plot = ggplot(data.melted)+
               ggtitle(paste0("number of genes resulting of the projection of the partionning on each organism (nb organism=", nrow(data),")"))+
               geom_boxplot(aes(x = partition,y = nb_geness, fill = partition))+
               geom_text_repel(data = rbind(data.melted_persistent,data.melted_shell,data.melted_cloud, data.melted_exact_core, data.melted_exact_accessory, data.melted_pangenome),
                               aes(x=partition, y=nb_geness, label = org))+
               scale_fill_manual(name = "partitioning", values = color_chart)
    }
    else{
        data <- cbind(data, pos = seq(nrow(data)))

        max_degree_log2p1 <- max(apply(data,1,FUN = function(x){
                sum(log2(as.numeric(x[10:12])+1))
            }))

        ori <- which(data$ori == T, arr.ind=T)
        data$ori <- NULL

        duplicated_fam     <- unique(data[duplicated(data$family),"family"])
        data$family <- ifelse(data$family %in% duplicated_fam, data$family, NA)
        data$family = as.factor(data$family)
        colors_duplicated_fam <- rainbow(length(duplicated_fam))
        names(colors_duplicated_fam) <- duplicated_fam

        data_melted <- melt(data, id.var=c("contig", "gene","family","nb_copy_in_org","partition","pos","strand","coord_start","coord_end"))
        data_melted$variable <- factor(data_melted$variable, levels = rev(c("persistent","shell","cloud")), ordered=TRUE)

        contig <- unique(data_melted$contig)
        contig_color <-  rainbow(length(contig))
        names(contig_color) <- contig

        data_melted$value <- log2(data_melted$value+1)
        plot = ggplot(data = data_melted)+
        ggtitle(paste0("plot corresponding to the file", org_name))+
        geom_bar(aes_string(x = "gene", y = "value", fill = "variable"),stat="identity", show.legend = FALSE)+
        scale_y_continuous(limits = c(-30, max_degree_log2p1), breaks = seq(0,ceiling(max_degree_log2p1)))+
        geom_hline(yintercept = 0)+
        geom_rect(aes_string(xmin ="pos-1/2", xmax = "pos+1/2", fill = "partition"), ymin = -10, ymax=-1, color = NA, show.legend = FALSE)+
        geom_hline(yintercept = -10)+
        geom_rect(aes_string(xmin ="pos-1/2", xmax = "pos+1/2", fill = "family"), ymin = -20, ymax=-11,  color = NA, show.legend = FALSE)+
        geom_hline(yintercept = -20)+
        geom_rect(aes_string(xmin ="pos-1/2", xmax = "pos+1/2", fill = "contig"), ymin = -30, ymax=-21,  color = NA)+
        geom_vline(xintercept = ori)+
        scale_fill_manual(values = c(color_chart,colors_duplicated_fam, contig_color), na.value = "grey80")+
        coord_polar()+
        ylab("log2(degree+1) of the families in wich each gene is")+
        theme(axis.line        = ggplot2::element_blank(),
              axis.text.x      = ggplot2::element_blank(),
              axis.ticks.x     = ggplot2::element_blank(),
              axis.title.x     = ggplot2::element_blank(),
              panel.background = ggplot2::element_blank(),
              panel.border     = ggplot2::element_blank(),
              panel.grid.major.x = ggplot2::element_blank(),
              panel.grid.minor.x = ggplot2::element_blank(),
              plot.background  = ggplot2::element_blank(),
              plot.margin      = grid::unit(c(0,0,0,0), "cm"),
              panel.spacing    = grid::unit(c(0,0,0,0), "cm"))
    }
    ggsave(paste0('"""+OUTPUTDIR+FIGURE_DIR+"""',org_name,'.pdf'), device = "pdf", height= 40, width = 49, plot)

}
########################### END PROJECTION #################################

    """
    logging.getLogger().info("Writing R script generating plot")
    with open(script_outfile,"w") as script_file:
        script_file.write(rscript)

#### START - NEED TO BE AT THE HIGHEST LEVEL OF THE MODULE TO ALLOW MULTIPROCESSING

shuffled_comb = []
evol = None
pan = None
options = None
EVOLUTION = None

def resample(index):
    global shuffled_comb
    nem_dir_path    = TMP_DIR+EVOLUTION_DIR+"/nborg"+str(len(shuffled_comb[index]))+"_"+str(index)
    stats = pan.partition(nem_dir_path     = nem_dir_path,
                          select_organisms = shuffled_comb[index],
                          beta             = options.beta_smoothing[0],
                          free_dispersion  = options.free_dispersion,
                          chunck_size      = options.chunck_size[0],
                          soft_core_th     = options.soft_core_threshold[0],
                          inplace          = False,
                          just_stats       = True,
                          nb_threads       = 1)
    if not options.keep_nem_temporary_files:
        shutil.rmtree(nem_dir_path)
    evol.write(",".join([str(len(shuffled_comb[index])),
                          str(stats["persistent"]) if stats["undefined"] == 0 else "NA",
                          str(stats["shell"]) if stats["undefined"] == 0 else "NA",
                          str(stats["cloud"]) if stats["undefined"] == 0 else "NA",
                          str(stats["undefined"]),
                          str(stats["exact_core"]),
                          str(stats["exact_accessory"]),
                          str(stats["soft_core"]),
                          str(stats["soft_accessory"]),
                          str(stats["exact_core"]+stats["exact_accessory"])])+"\n")
    evol.flush()

# def replication(index):
#     subset = random.sample(pan.organisms, 2)
#     old_res = None
#     before_pangenome = None
#     before_shell = None
#     for i in range(pan.nb_organisms):
#         res = pan.partition(subset,just_stats=False,inplace=False)
#         if old_res is not None:
#             pangenome    = set(res["exact_accessory"]+res["exact_core"])
#             new_families = pangenome - before_pangenome
#             new_shell    = set(res["shell"]) - before_shell
#     subset.add(random.sample(pan.organisms-subset,STEP))
#     evol.write("\t".join([str(len(shuffled_comb[index])),
#                           str(stats["persistent"]) if stats["undefined"] == 0 else "NA",
#                           str(stats["shell"]) if stats["undefined"] == 0 else "NA",
#                           str(stats["cloud"]) if stats["undefined"] == 0 else "NA",
#                           str(stats["exact_core"]),
#                           str(stats["exact_accessory"]),
#                           str(stats["exact_core"]+stats["exact_accessory"])])+"\n")
#     evol.flush()


#### END - NEED TO BE AT THE HIGHEST LEVEL OF THE MODULE TO ALLOW MULTIPROCESSING

def __main__():
    parser = argparse.ArgumentParser(prog = "ppanggolin",
                                     description='Build a partitioned pangenome graph from annotated genomes (GFF files) and gene families (TSV files). Reserved words are: '+' '.join(RESERVED_WORDS), 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-?', '--version', action='version', version=pkg_resources.get_distribution("ppanggolin").version)
    parser.add_argument('-o', '--organisms', type=argparse.FileType('r'), nargs=1, metavar=('ORGANISMS_FILE'), help="""
    File: A tab-delimited file containing at least 2 mandatory fields by row and as many optional fields as the number of well assembled circular contigs. 
    Each row corresponds to an organism to be added to the pangenome.
    The first field is the organism ID.
    The organism ID can be any string but must be unique and can't contain any space, quote, double quote and reserved word.
    The second field is the gff file containing the annotations associated to the organism. 
    This path can be absolute or relative. 
    The gff file must contain an ID for each line.
    Accepted types are CDS or xRNA (rRNA,tRNA,tmRNA) in the type column.
    The contig ID and gene ID can be any string but must be unique and can't contain any space, quote, double quote, pipe and reserved words.
    (optional): The next fields contain the name of perfectly assembled circular contigs (to take in account the link between the first and the last gene in the graph). 
    In this case, it is mandatory to provide the contig size in the gff files either by adding a "region" type having the correct contig ID attribute or using a '##sequence-region' pragma.
    """, required=True)
    parser.add_argument('-gf', '--gene_families', type=argparse.FileType('r'), nargs=1, metavar=('FAMILIES_FILE'), help="""
    File: A tab-delimited file containing the gene families. Each row contains 2 or 3 fields.
    The first field is the family ID. 
    The second field is the gene IDs associated with this family. 
    The third field (optional) is a flag "F" to specify if the gene is a gene fragment (empty otherwise).
    If several consecutive genes belonging to the same gene families have the flag, then the no reflexive links are drawn between them in the graph.
    The family ID can be any string but must be unique and can't contain any space, quote, double quote and reserved word. 
    The family ID is often the name of the most representative sequence of the family.
    Gene IDs can be any string corresponding to the IDs of the gff files. They must be uniques and can't contain any spaces, quote, double quote, pipe and reserved words.
    """,  required=True)
    parser.add_argument('-od', '--output_directory', type=str, nargs=1, default=["PPanGGOLiN_outputdir"+strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", localtime())+"_PID"+str(os.getpid())], metavar=('OUTPUT_DIR'), help="""
    Dir: The output directory""")
    parser.add_argument('-td', '--temporary_directory', type=str, nargs=1, default=["/tmp/PPanGGOLiN_tempdir"+strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", localtime())+"_PID"+str(os.getpid())], metavar=('TMP_DIR'), help="""
    Dir: Temporary directory to store nem intermediate files""")
    parser.add_argument('-f', '--force', action="store_true", help="""
    Flag: Force overwriting existing output directory""")
    parser.add_argument('-r', '--remove_high_copy_number_families', type=int, nargs=1, default=[0], metavar=('REPETITION_THRESHOLD'), help="""
    Positive Number: Remove families having a number of copy of gene in a single family above or equal to this threshold in at least one organism (0 or negative values are ignored). 
    """)#When -update is set, only work on new organisms added
    parser.add_argument('-s', '--infer_singletons', default=False, action="store_true", help="""
    Flag: If a gene id found in a gff file is absent of the gene families file, the singleton will be automatically infered as a gene families having a single element. 
    if this argument is not set, the program will raise KeyError exception if a gene id found in a gff file is absent of the gene families file.""")
    #    parser.add_argument("-up", "--update", default = None, type=argparse.FileType('r'), nargs=1, help="""
    # Pangenome Graph to be updated (in gexf format)""")
    parser.add_argument("-u", "--untangle", type=int, default = 0, nargs=1, help="""
    Flag: (in test) max size of the untangled paths to be untangled""")
    parser.add_argument("-b", "--beta_smoothing", default = [float("0.5")], type=float, nargs=1, metavar=('BETA_VALUE'), help = """
    Positive decimal number: This option determines the strength of the smoothing (:math:beta) of the partitioning based on the graph topology (using a Markov Random Field). 
    b must be a positive float, b = 0.0 means to discard spatial smoothing and 1.00 means strong smoothing (can be above 1.00 but it is not advised).
    0.5 is generally advised as a good trade off.
    """)
    parser.add_argument("-th", "--soft_core_threshold", type = float, nargs = 1, default = [0.95], metavar=('SOFT_CORE_THRESHOLD'), help = """
    Postitive decimal number: A value between 0 and 1 providing the threshold ratio of presence to attribute a gene families to the soft core genome""")
    parser.add_argument("-fd", "--free_dispersion", default = False, action="store_true", help = """
    Flag: Specify if the dispersion around the centroid vector of each partition is the same for all the organisms or if the dispersion is free
    """)
    parser.add_argument("-kf", "--keep_nem_temporary_files", default=False, action="store_true", help="""
    Flag: Delete temporary files used by NEM""")
    parser.add_argument("-cg", "--compress_graph", default=False, action="store_true", help="""
    Flag: Compress (using gzip) the files containing the partionned pangenome graph""")
    parser.add_argument("-c", "--cpu", default=[1],  type=int, nargs=1, metavar=('NB_CPU'), help="""
    Positive Number: Number of cpu to use (several cpu will be used only if the option -e is set or/and if the -ck option is below the number of organisms provided)""")
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="""
    Flag: Show all messages including debugging ones""")
    # parser.add_argument("-as", "--already_sorted", default=False, action="store_true", help="""
    # Accelerate loading of gff files if there are sorted by the coordinate of gene annotations (starting point) for each contig""")
    #parser.add_argument("-l", "--freemem", default=False, action="store_true", help="""
    #Free the memory elements which are no longer used""")
    parser.add_argument("-p", "--plots", default=False, action="store_true", help="""
    Flag: Run the Rscript generating the plots (required: R in the path and the packages ggplot2, ggrepel, data.table, minpack.lm and reshape2 to be installed).""")
    # parser.add_argument("-di", "--directed", default=False, action="store_true", help="""
    # generate directed graph
    # """)
    parser.add_argument("-e", "--evolution", default=False, action="store_true", help="""
    Flag: Partition the pangenome using multiple subsamples of a croissant number of organisms in order to obtain a curve of the evolution of the pangenome metrics
    """)
    parser.add_argument("-ep", "--evolution_resampling_param", nargs=5, default=[0.1,10,10,1,float("Inf")], metavar=('RESAMPLING_RATIO','MINIMUM_RESAMPLING','MAXIMUM_RESAMPLING','STEP','LIMIT'), help="""
    5 Positive Numbers (or Inf for the last one):
    1st argument is the resampling ratio (FLOAT)
    2st argument is the minimum number of resampling for each number of organisms (INTEGER)
    3nd argument is the maximum number of resampling for each number of organisms (INTEGER or Inf)
    4rd argument is the step between each number of organisms (INTEGER)
    5rd argument is the limit of the size of the samples (INTEGER or Inf)
    """)
    parser.add_argument("-pr", "--projection", type = int, nargs = "+", metavar=('LINE_NUMBER_OR_ZERO'), help="""
    Positive Number: project the graph as a circos plot on each organism.
    Expected parameters are the line number (1 based) of each organism on which the graph will be projected.
    It provides a circular plot (well-assembled representative organisms must be prefered).
    0 means all organisms (it is discouraged to use -p and -pr 0 in the same time because the projection of the graph on all the organisms can take a long time).
    """)
    parser.add_argument("-ck", "--chunck_size", type = int, nargs = 1, default = [500], metavar=('SIZE'), help="""
    Positive Number: Size of the chunks to perform the partioning by chunks.
    If the number of organisms used is higher than SIZE, the partioning will be performed by chunks of size SIZE
    """)
    parser.add_argument("-mt", "--metadata", type=argparse.FileType('r'), default = [None], nargs=1, metavar=('METADATA_FILE'), help="""
    File: It is possible to add metainformation to the pangenome graph. These information must be associated to each organism via a METADATA_FILE. During the construction of the graph, metainformation about the organisms are used to label the covered edges.
    METADATA_FILE is a tab-delimitated file. The first line contains the names of the attributes and the following lines contain associated information for each organism.
    Ex:
    org metadata1   metadata2
    org1    x1  y1
    org2    x2  y1
    ...
    Metadata can't contain reserved word or exact organism name.
    """)
    parser.add_argument("-ra", "--add_rna_to_the_pangenome", default = False, action="store_true", help = """
    Flag: If the specified the xRNA (rRNA,tRNA,tmRNA...) are added to the pangenome graph, (as for proteins, RNA genes need to be clustered)""")
    #parser.add_argument("-dc", "--distance_CDS_fragments", type = int, nargs = 1, default = [-9], metavar=('DISTANCE_CDS_FRAGMENTS'), help = """
    #Number: When several consecutive genes belonging to the same gene families separated by less or equals than DISTANCE_CDS_FRAGMENTS (in nucleotides) are found, there are considered as CDS fragments and then no reflexive links are generated between the associated gene families. Negative values are considered as overlapping CDS fragments (short gene overlapping seems to be frequent in bacteria https://doi.org/10.1186/1471-2164-15-721 ).""")
    # parser.add_argument("-ss", "--subpartition_shell", default = 0, type=int, nargs=1, help = """
    # Number: (in test) Subpartition the shell genome in k subpartitions, k can be detected automatically using k = -1, if k = 0 the partioning will used the first column of metadata to subpartition the shell""")
    parser.add_argument("-l", "--compute_layout", default = False, action="store_true", help = """
    Flag: (in test) precalculated the ForceAtlas2 layout""")

    global options
    options = parser.parse_args()

    level = logging.INFO
    if options.verbose:
        level = logging.DEBUG

    logging.basicConfig(stream=sys.stdout, level = level, format = '\n%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    logging.getLogger().info("Command: "+" ".join([arg for arg in sys.argv]))
    logging.getLogger().info("PPanGGOLiN version: "+pkg_resources.get_distribution("ppanggolin").version)
    logging.getLogger().info("Python version: "+sys.version)
    logging.getLogger().info("Networkx version: "+nx.__version__)
    global OUTPUTDIR
    global TMP_DIR

    
    OUTPUTDIR = options.output_directory[0]
    TMP_DIR   = options.temporary_directory[0]

    logging.getLogger().info("Output directory is: "+OUTPUTDIR)
    logging.getLogger().info("Temporary directory is: "+TMP_DIR)

    list_dir = ["",FIGURE_DIR,PARTITION_DIR,PATH_DIR]
    if options.projection:
        list_dir.append(PROJECTION_DIR)
    if options.metadata[0]:
        list_dir.append(METADATA_DIR)
    if options.evolution:
        list_dir.append(EVOLUTION_DIR)
        (RESAMPLING_RATIO, RESAMPLING_MIN, RESAMPLING_MAX, STEP, LIMIT) = options.evolution_resampling_param
        (RESAMPLING_RATIO, RESAMPLING_MIN, RESAMPLING_MAX, STEP, LIMIT) = (float(RESAMPLING_RATIO), int(RESAMPLING_MIN), int(RESAMPLING_MAX) if str(RESAMPLING_MAX).upper() != "Inf" else sys.maxsize, int(STEP), int(LIMIT) if str(LIMIT).upper() != "INF" else sys.maxsize)
    for directory in list_dir:
        if not os.path.exists(OUTPUTDIR+directory):
            os.makedirs(OUTPUTDIR+directory)
        elif not options.force:
            logging.getLogger().error(OUTPUTDIR+directory+" already exists")
            exit(1)

    #-------------
    metadata = None

    if options.metadata[0]:
        metadata = list()
        attribute_names = list()
        metadata = OrderedDict()
        for num, line in enumerate(options.metadata[0]):
            elements = [el.strip() for el in line.split("\t")]
            if num == 0:
                attribute_names = elements
            else:
                o = elements[0]
                if o == "":
                    logging.getLogger().error("Organism name empty in metadata file is empty")
                    exit(1)
                elements = [e if e != "" else None for e in elements]
                metadata[o]=dict(zip(attribute_names,elements))
    start_loading = time()
    global pan

    pan = PPanGGOLiN("file",
                     options.organisms[0],
                     options.gene_families[0],
                     options.remove_high_copy_number_families[0],
                     options.infer_singletons,
                     options.add_rna_to_the_pangenome,
                     False)#options.directed)

    # with open(OUTPUTDIR+CDS_FRAGMENTS_FILE_PREFIX+".csv","w") as CDS_fragments_file:
    #     CDS_fragments_file.write(",".join(["CDS_fragment","CDS_fragment_length","corresponding_gene_family","gene_family_median_length"])+"\n")
    #     for gene, family in pan.CDS_fragments.items():
    #         info = pan.get_gene_info(gene)
    #         CDS_fragments_file.write(",".join([gene,
    #                                            str(info["END"]-info["START"]),
    #                                            family,
    #                                            str(median(pan.neighbors_graph.node[family]["length"]))])+"\n")

    # # if options.update is not None:
    # #     pan.import_from_GEXF(options.update[0])
    # end_loading_file = time.time()
    # #-------------

    # #-------------
    
    # start_neighborhood_computation = time.time()
    end_loading = time()
    #-------------

    #-------------
    logging.getLogger().info("Partitioning...")

    start_partitioning = time()
    pan.partition(nem_dir_path    = TMP_DIR+NEM_DIR,
                  beta            = options.beta_smoothing[0],
                  free_dispersion = options.free_dispersion,
                  chunck_size     = options.chunck_size[0],
                  soft_core_th    = options.soft_core_threshold[0],
                  inplace         = True,
                  just_stats      = False,
                  nb_threads      = options.cpu[0])
    end_partitioning = time()
    #-------------
    if options.metadata[0]:
        if options.force:
            shutil.rmtree(OUTPUTDIR+METADATA_DIR)
        pan.get_gene_families_related_to_metadata(metadata,OUTPUTDIR+METADATA_DIR)

    # if options.subpartition_shell:
    #     if options.subpartition_shell[0] <0:
    #         Q = pan.partition_shell(Q = "auto")
    #         logging.getLogger().info(str(Q)+" subpartitions has been used to subpartition the shell genome...")
    #     elif options.subpartition_shell[0]==0:
    #         init=defaultdict(set)
    #         for orgs, metad in metadata.items():
    #             init[metad].add(orgs)
    #         pan.partition_shell(init_using_qual=init)
    #     else:
    #         pan.partition_shell(options.subpartition_shell[0])


    #-------------
    # th = 100

    # cpt_partition = {}
    # for fam in pan.neighbors_graph.node:
    #     cpt_partition[fam]= {"persistent":0,"shell":0,"cloud":0}

    # cpt = 0
    # validated = set()
    # while(len(validated)<pan.pan_size):
    #     sample = pan.sample(n=100)
    #     sample.neighborhood_computation(options.undirected, light=True)
    #     sample.partition(EVOLUTION+"/"+str(cpt), float(50), options.free_dispersion)#options.beta_smoothing[0]
    #     cpt+=1
    #     for node,data in pan.neighbors_graph.nodes(data=True):
    #         cpt_partition[node][data["partition"]]+=1
    #         if sum(cpt_partition[node].values()) > th:
    #             validated.add(node)

    # for fam, data in cpt_partition.items():
    #     pan.neighbors_graph.nodes[fam]["partition_bis"]= max(data, key=data.get)


    # print(cpt_partition)
    #-------------

    #-------------
    start_writing_output_file = time()

    logging.getLogger().info("Extract and label paths")
    start_paths = time()
    correlated_paths = pan.extract_shell_paths()
    end_paths = time()

    with open(OUTPUTDIR+"/"+PATH_DIR+"/"+CORRELATED_PATHS_PREFIX+"_vectors.csv","w") as correlated_paths_vectors, open(OUTPUTDIR+"/"+PATH_DIR+"/"+CORRELATED_PATHS_PREFIX+"_confidences.csv","w") as correlated_paths_confidences:
        header = []
        for i, (path, vector) in enumerate(correlated_paths.items()):
            if i==0:
                header = list(vector.keys())
                correlated_paths_vectors.write(",".join(["Gene","Non-unique Gene name","Annotation","No. isolates","No. sequences","Avg sequences per isolate","Accessory Fragment","Genome Fragment","Order within Fragment","Accessory Order with Fragment","QC","Min group size nuc","Max group size nuc","Avg group size nuc"]+header)+"\n")
                correlated_paths_confidences.write(",".join(["correlated_paths"]+header)+"\n")
            correlated_paths_vectors.write(",".join([path]+["","",str(sum(vector.values())),str(sum(vector.values())),"","","","","","","","",""]+[str(int(round(v))) for v in vector.values()])+("\n" if i < len(correlated_paths)-1 else ""))
            correlated_paths_confidences.write(",".join([path]+["","",str(sum(vector.values())),str(sum(vector.values())),"","","","","","","","",""]+[str(v) for v in vector.values()])+("\n" if i < len(correlated_paths)-1 else ""))

    if options.compute_layout:
        logging.getLogger().info("Computing layout")
        start_layout = time()
        pan.compute_layout()#multiThreaded=options.cpu[0])
        end_layout = time()

    logging.getLogger().info("Writing GEXF file")
    pan.export_to_GEXF(OUTPUTDIR+GRAPH_FILE_PREFIX+(".gz" if options.compress_graph else ""), options.compress_graph, metadata)
    logging.getLogger().info("Writing GEXF light file")
    pan.export_to_GEXF(OUTPUTDIR+GRAPH_FILE_PREFIX+"_light"+(".gz" if options.compress_graph else ""), options.compress_graph, metadata, False,False)
    with open(OUTPUTDIR+"/pangenome.txt","w") as pan_text:
        for partition, families in pan.partitions.items():
            file = open(OUTPUTDIR+PARTITION_DIR+"/"+partition+".txt","w")
            if len(families):
                file.write("\n".join(families)+"\n")
                if partition in set(["exact_core","exact_accessory"]):
                    pan_text.write("\n".join(families)+"\n")
            file.close()
    pan.write_matrix(OUTPUTDIR+MATRIX_FILES_PREFIX)
    if options.projection:
        logging.getLogger().info("Projection...")
        start_projection = time()
        pan.projection(OUTPUTDIR+PROJECTION_DIR, [pan.organisms.__getitem__(index-1) for index in options.projection] if options.projection[0] > 0 else list(pan.organisms))
        end_projection = time()
    end_writing_output_file = time()
    pan.ushaped_plot(OUTPUTDIR+FIGURE_DIR+"/"+USHAPE_PLOT_PREFIX)
    #pan.tile_plot(OUTPUTDIR+FIGURE_DIR)
    del pan.annotations # no more required for the following process
    # print(pan.partitions_by_organisms)
    # partitions_by_organisms_file = open(OUTPUTDIR+"/partitions_by_organisms.txt","w")
    # exact_by_organisms_file = open(OUTPUTDIR+"/exacte_by_organisms.txt","w")
    # for org, partitions in pan.partitions_by_organisms.items(): 
    #     partitions_by_organisms_file.write(org+"\t"+str(len(partitions["persistent"]))+
    #                                            "\t"+str(len(partitions["shell"]))+
    #                                            "\t"+str(len(partitions["cloud"]))+"\n")
    #     exact_by_organisms_file.write(org+"\t"+str(len(partitions["exact_core"]))+
    #                                       "\t"+str(len(partitions["exact_accessory"]))+"\n")
    # partitions_by_organisms_file.close()
    # exact_by_organisms_file.close()
    #-------------

    logging.getLogger().info(pan)
    with open(OUTPUTDIR+"/"+SUMMARY_STATS_FILE_PREFIX+".txt","w") as file_stats:
        file_stats.write(str(pan))
    #-------------

    if options.untangle>0:
        pan.untangle_neighbors_graph(options.untangle[0])
        pan.export_to_GEXF(OUTPUTDIR+GRAPH_FILE_PREFIX+(".gz" if options.compress_graph else ""), options.compress_graph, metadata,"untangled_neighbors_graph" )

    

    plot_Rscript(script_outfile = OUTPUTDIR+"/"+SCRIPT_R_FIGURE, verbose=options.verbose)

    if options.evolution:

        logging.getLogger().info("Evolution...")

        start_evolution = time()
        if not options.verbose:
            logging.disable(logging.INFO)# disable INFO message to not disturb the progess bar
            logging.disable(logging.WARNING)# disable WARNING message to not disturb the progess bar
        combinations = samplingCombinations(list(pan.organisms), sample_ratio=RESAMPLING_RATIO, sample_min=RESAMPLING_MIN, sample_max=RESAMPLING_MAX)
        
        global shuffled_comb
        shuffled_comb = combinations
        shuffled_comb = [OrderedSet(comb) for nb_org, combs in combinations.items() for comb in combs if nb_org%STEP == 0 and nb_org<=LIMIT]
        shuffle(shuffled_comb)

        global evol
        evol =  open(OUTPUTDIR+EVOLUTION_DIR+EVOLUTION_STATS_FILE_PREFIX+".csv","w")

        evol.write(",".join(["nb_org","persistent","shell","cloud","undefined","exact_core","exact_accessory","soft_core","soft_accessory","pangenome"])+"\n")
        evol.write(",".join([str(pan.nb_organisms),    
                              str(len(pan.partitions["persistent"])),
                              str(len(pan.partitions["shell"])),
                              str(len(pan.partitions["cloud"])),
                              str(len(pan.partitions["undefined"])),
                              str(len(pan.partitions["exact_core"])),
                              str(len(pan.partitions["exact_accessory"])),
                              str(len(pan.partitions["soft_core"])),
                              str(len(pan.partitions["soft_accessory"])),
                              str(len(pan.partitions["exact_accessory"])+len(pan.partitions["exact_core"]))])+"\n")
        evol.flush()
        with ProcessPoolExecutor(options.cpu[0]) as executor:
            futures = [executor.submit(resample,i) for i in range(len(shuffled_comb))]
            for f in tqdm(as_completed(futures), total = len(shuffled_comb), unit = 'pangenome resampled'):
                ex = f.exception()
                if ex:
                    executor.shutdown(wait=False)
                    raise ex
        evol.close()

        def heap_law(N, kappa, gamma):
            return kappa*N**(gamma)

        annotations=[]
        traces = []
        data_evol = pandas.read_csv(OUTPUTDIR+EVOLUTION_DIR+EVOLUTION_STATS_FILE_PREFIX+".csv",index_col=False)
        params_file = open(OUTPUTDIR+EVOLUTION_DIR+EVOLUTION_PARAM_FILE_PREFIX+".csv","w")
        params_file.write("partition,kappa,gamma,kappa_std_error,gamma_std_error\n")
        for partition in list(pan.partitions.keys())+["pangenome"]:
            half_stds = pandas.Series({i:numpy.std(data_evol[data_evol["nb_org"]==i][partition])/2 for i in range(1,pan.nb_organisms+1)}).dropna()
            mins = pandas.Series({i:numpy.min(data_evol[data_evol["nb_org"]==i][partition]) for i in range(1,pan.nb_organisms+1)}).dropna()
            maxs = pandas.Series({i:numpy.max(data_evol[data_evol["nb_org"]==i][partition]) for i in range(1,pan.nb_organisms+1)}).dropna()
            medians = pandas.Series({i:numpy.median(data_evol[data_evol["nb_org"]==i][partition]) for i in range(1,pan.nb_organisms+1)}).dropna()
            initial_kappa_gamma = numpy.array([0.0, 0.0])
            res = optimization.curve_fit(heap_law, medians.index, medians, initial_kappa_gamma)
            kappa, gamma = res[0]
            error_k,error_o = numpy.sqrt(numpy.diag(res[1])) # to calculate the fitting error. The variance of parameters are the diagonal elements of the variance-co variance matrix, and the standard error is the square root of it. source https://stackoverflow.com/questions/25234996/getting-standard-error-associated-with-parameter-estimates-from-scipy-optimize-c
            if numpy.isinf(error_k) and numpy.isinf(error_o):
                params_file.write(",".join([partition,"NA","NA","NA","NA"])+"\n")
            else:
                params_file.write(",".join([partition,str(kappa),str(gamma),str(error_k),str(error_o)])+"\n")
                regression = numpy.apply_along_axis(heap_law, 0, range(1,pan.nb_organisms+1), kappa, gamma)
                traces.append(go.Scatter(x=list(range(1,pan.nb_organisms+1)), 
                                    y=regression, 
                                    name = partition,
                                    line = dict(color = COLORS[partition],
                                                width = 4,
                                                dash = 'dash'),
                                    visible = "legendonly" if partition == "undefined" else True))
                annotations.append(dict(x=pan.nb_organisms,
                                        y=heap_law(pan.nb_organisms,kappa, gamma),
                                        ay=0,
                                        ax=100,
                                        text="F="+str(round(kappa,0))+"N"+"<sup>"+str(round(gamma,5))+"</sup>",
                                        showarrow=True,
                                        arrowhead=7,
                                        font=dict(size=16,color='white'),
                                        align='center',
                                        arrowcolor=COLORS[partition],
                                        bordercolor='#c7c7c7',
                                        borderwidth=2,
                                        borderpad=4,
                                        bgcolor=COLORS[partition],
                                        opacity=0.8))
            
            traces.append(go.Scatter(x=medians.index, 
                                     y=medians, 
                                     name = "medians",
                                     mode="lines+markers",
                                     error_y=dict(type='data',
                                                     symmetric=False,
                                                     array=maxs.subtract(medians),
                                                     arrayminus=medians.subtract(mins),
                                                     visible=True,
                                                     color = COLORS[partition],
                                                     thickness =1),
                                     line = dict(color = COLORS[partition],
                                                 width = 1),
                                     marker=dict(color = COLORS[partition]),
                                     visible = "legendonly" if partition == "undefined" else True))
            up = medians.add(half_stds)
            down = medians.subtract(half_stds)
            down[down < 0] = 0
            sd_area = up.append(down[::-1])
            traces.append(go.Scatter(x=sd_area.index, 
                                     y=sd_area, 
                                     name = "sd",
                                     fill='toself',
                                     mode="lines",
                                     hoveron="points",
                                     #hovertext=[str(round(e)) for e in half_stds.multiply(2)],
                                     line=dict(color=COLORS[partition]),
                                     marker=dict(color = COLORS[partition]),
                                     visible = "legendonly" if partition == "undefined" else True))
        layout = go.Layout(
            title = "Evolution curve ",
            titlefont = dict(size = 20),
            xaxis = dict(title='size of genome subsets (N)'),
            yaxis = dict(title='# of gene families (F)'),
            annotations=annotations)
        fig = go.Figure(data=traces, layout=layout)
        out_plotly.plot(fig, filename=OUTPUTDIR+"/"+FIGURE_DIR+"/"+EVOLUTION_CURVE_PREFIX+".html", auto_open=False)
        params_file.close()

        # evolution_curve = Highchart(width = 1800, height = 800)
        # options_evolution_curve_plot={
        # 'title': {'text':'Evolution curve of the pangenome metrics with a growing number of organisms'},
        # 'xAxis': {'tickInterval': 1, 'categories': list(range(1,min(pan.nb_organisms,LIMIT)+1), 'title':{'text':'size of the subsets'}},
        # 'yAxis': {'allowDecimals': False, 'title' : {'text':'# of families'}},
        # 'tooltip': {'headerFormat': '<span style="font-size:11px"># of orgs: <b>{point.x}</b></span><br>',
        #             'pointFormat': '<span style="color:{point.color}">{series.name}</span>: {point.y}<br/>',
        #             'shared': True}
        # }
        # evolution_curve.set_dict_options(options_evolution_curve_plot)

        # for i, line in enumerate(open(OUTPUTDIR+EVOLUTION_DIR+EVOLUTION_STATS_FILE_PREFIX+".txt","r")):
        #     if i == 0:
        #         continue
        #     elements = [int(e.strip()) for e in line.split(",")]


        
        # pandas.read_csv(file_evol_path+"/projections/nb_genes.csv", sep = "\t").dropna()

        # ushaped_plot.add_data_set(persistent_values,'column','Persistent', color = COLORS["persistent"])
        end_evolution = time()
        logging.disable(logging.NOTSET)#restaure info and warning messages 

    # if options.new_genes_evolution:
    #     logging.getLogger().info("New genes evolution...")
    #     start_evolution = time()
    #     if not options.verbose:
    #         logging.disable(logging.INFO)# disable INFO message to not disturb the progess bar
    #         logging.disable(logging.WARNING)# disable WARNING message to not disturb the progess bar
        
        
    #         with ProcessPoolExecutor(options.cpu[0]) as executor:
    #             futures = [executor.submit(resample,i) for i in range(len(shuffled_comb))]

    #         for f in tqdm(as_completed(futures), total = len(shuffled_comb), unit = 'pangenome resampled'):
    #             ex = f.exception()
    #             if ex:
    #                 executor.shutdown(wait=False)
    #                 raise ex
    #     evol.close()

    #     end_evolution = time()
    #     logging.disable(logging.NOTSET)#restaure info and warning messages 
    #-------------

    logging.getLogger().info("\n"+
    "Execution time of loading and neighborhood computation: """ +str(round(end_loading-start_loading, 2))+" s\n"+
    "Execution time of partitioning: " +str(round(end_partitioning-start_partitioning, 2))+" s\n"+
    "Execution time of path detection: " +str(round(end_paths-start_paths, 2))+" s\n"+
    (("Execution time of layout computation: " +str(round(end_layout-start_layout, 2))+" s\n") if options.compute_layout else "")+
    "Execution time of writing output files: " +str(round(end_writing_output_file-start_writing_output_file, 2))+" s\n"+
    (("Execution time of evolution: " +str(round(end_evolution-start_evolution, 2))+" s\n") if options.evolution else "")+
    "Total execution time: " +str(round(time()-start_loading, 2))+" s\n")

    logging.getLogger().info("""The pangenome computation is complete.""")

    if options.plots:
        logging.getLogger().info("Running R script generating plot")
        cmd = "Rscript "+OUTPUTDIR+SCRIPT_R_FIGURE
        logging.getLogger().info("""Several plots will be generated using R (in the directory: """+OUTPUTDIR+FIGURE_DIR+""").
    If R and the required package (ggplot2, reshape2, ggrepel(>0.6.6), data.table, minpack.lm) are not installed don't worry, the R script is saved in the directory. To generate the figures later, just use the following command :
    """+cmd)
        
        logging.getLogger().info(cmd)
        proc = subprocess.Popen(cmd, shell=True)
        proc.communicate()

    if not options.keep_nem_temporary_files:
        pan.delete_nem_intermediate_files()  

    logging.getLogger().info("Finished !")
    exit(0)

if __name__ == "__main__":
    __main__()
