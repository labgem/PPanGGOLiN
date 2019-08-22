#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse
from multiprocessing import Pool
import time
from collections import Counter
import logging
import pkg_resources
from statistics import median
#installed libraries
from tqdm import tqdm

#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import write_compressed_or_not, mkOutdir
from ppanggolin.formats import checkPangenomeInfo

#global variable to store the pangenome
pan = None

def writeGEXFheader(gexf, light):
    if not light:
        index = pan.getIndex()#has been computed already
    gexf.write('<?xml version="1.1" encoding="UTF-8"?>\n<gexf xmlns:viz="http://www.gexf.net/1.3draft/viz" xmlns="http://www.gexf.net/1.3draft" version="1.3">\n')
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
    gexf.write('      <attribute id="10" title="nb_organisms" type="long" />\n')
    if not light:
        for org, orgIndex in index.items():
            gexf.write(f'      <attribute id="{orgIndex + 12}" title="{org.name}" type="string" />\n')
    
    gexf.write('    </attributes>\n')
    gexf.write('    <attributes class="edge" mode="static">\n')
    gexf.write('      <attribute id="11" title="nb_genes" type="long" />\n')
    if not light:
        for org, orgIndex in index.items():
            gexf.write(f'      <attribute id="{orgIndex + len(index) + 12}" title="{org.name}" type="long" />\n')
    # gexf.write('      <attribute id="12" title="nb_organisms" type="long" />\n')#useless because it's the weight of the edge
    gexf.write('    </attributes>\n')
    gexf.write('    <meta>\n')
    gexf.write(f'      <creator>PPanGGOLiN {pkg_resources.get_distribution("ppanggolin").version}</creator>\n')
    gexf.write('    </meta>\n')

def writeGEXFnodes(gexf, light, soft_core = 0.95):
    gexf.write('    <nodes>\n')
    colors = {"persistent":'a="0" b="7" g="165" r="247"','shell':'a="0" b="96" g="216" r="0"', 'cloud':'a="0" b="255" g="222" r="121"'}
    if not light:
        index = pan.getIndex()

    for fam in pan.geneFamilies:
        name = Counter()
        product = Counter()
        gtype = Counter()
        l = []
        for gene in fam.genes:
            name[gene.name] +=1
            product[gene.product] += 1
            gtype[gene.type] += 1
            l.append(gene.stop - gene.start)
        
        gexf.write(f'      <node id="{fam.ID}" label="{fam.name}">\n')
        gexf.write(f'        <viz:color {colors[fam.namedPartition]} />\n')
        gexf.write(f'        <viz:size value="{len(fam.organisms)}" />\n')
        gexf.write(f'        <attvalues>\n')
        gexf.write(f'          <attvalue for="0" value="{len(fam.genes)}" />\n')
        gexf.write(f'          <attvalue for="1" value="{name.most_common(1)[0][0]}" />\n')
        gexf.write(f'          <attvalue for="2" value="{product.most_common(1)[0][0]}" />\n')
        gexf.write(f'          <attvalue for="3" value="{gtype.most_common(1)[0][0]}" />\n')
        gexf.write(f'          <attvalue for="4" value="{fam.namedPartition}" />\n')
        gexf.write(f'          <attvalue for="5" value="{fam.partition}" />\n')
        gexf.write(f'          <attvalue for="6" value="{"exact_accessory" if len(fam.organisms) != len(pan.organisms) else "exact_core"}" />\n')
        gexf.write(f'          <attvalue for="7" value="{"soft_core" if len(fam.organisms) > (len(pan.organisms)*soft_core) else "soft_accessory"}" />\n')
        gexf.write(f'          <attvalue for="8" value="{round(sum(l) / len(l),2)}" />\n')
        gexf.write(f'          <attvalue for="9" value="{ int(median(l))}" />\n')
        gexf.write(f'          <attvalue for="10" value="{len(fam.organisms)}" />\n')
        if not light:
            for org, genes in fam.getOrgDict().items():
                gexf.write(f'          <attvalue for="{index[org]+12}" value="{"|".join([ gene.ID for gene in genes])}" />\n')
        gexf.write(f'        </attvalues>\n')
        gexf.write(f'      </node>\n')
    gexf.write('    </nodes>\n')

def writeGEXFedges(gexf, light):
    gexf.write('    <edges>\n')
    edgeids = 0
    index = pan.getIndex()

    for edge in pan.edges: 
        gexf.write(f'      <edge id="{edgeids}" source="{edge.source.ID}" target="{edge.target.ID}" weight="{len(edge.organisms)}">\n')
        gexf.write(f'        <viz:thickness value="{len(edge.organisms)}" />\n')
        gexf.write('        <attvalues>\n')
        gexf.write(f'          <attribute id="11" value="{len(edge.genePairs)}" />\n')
        if not light:
            for org, genes in edge.getOrgDict().items():
                gexf.write(f'          <attvalue for="{index[org]+len(index)+12}" value="{len(genes)}" />\n')
        gexf.write('        </attvalues>\n')
        gexf.write('      </edge>\n')
        edgeids+=1
    gexf.write('    </edges>\n')

def writeGEXFend(gexf):
    gexf.write("  </graph>")
    gexf.write("</gexf>")

def writeGEXF(output, light = True, soft_core = 0.95, compress=False):
    txt = "Writing the gexf file for the pangenome graph..."
    if light:
        txt = "Writing the light gexf file for the pangenome graph..."
    logging.getLogger().info(txt)
    outname = output + "/pangenomeGraph"
    outname += "_light" if light else ""
    outname += ".gexf"
    with write_compressed_or_not(outname,compress) as gexf:
        writeGEXFheader(gexf, light)
        writeGEXFnodes(gexf, light)
        writeGEXFedges(gexf, light)
        writeGEXFend(gexf)
    logging.getLogger().info(f"Done writing the gexf file : '{outname}'")

def writeMatrix(sep, ext, output, compress=False, geneNames = False):
    logging.getLogger().info(f"Writing the .{ext} file ...")
    outname = output + "/matrix." + ext
    with write_compressed_or_not(outname,compress) as matrix:

        index_org = {}
        default_dat = []
        for index, org in enumerate(pan.organisms):
            default_dat.append('0')
            index_org[org] = index

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
                                +['"'+str(org)+'"' for org in pan.organisms])+"\n")#15
        default_genes = ['""'] * len(pan.organisms) if geneNames else ["0"] * len(pan.organisms)
        org_index = pan.getIndex()#should just return things
        for fam in pan.geneFamilies:
            genes = default_genes.copy()
            l = []
            product = Counter()
            for org, gene_list in fam.getOrgDict().items():
                genes[org_index[org]] = " ".join([ '"' + str(gene) + '"' for gene in gene_list]) if geneNames else str(len(gene_list))
                for gene in gene_list:
                    l.append(gene.stop - gene.start)
                    product[gene.product] +=1

            l = [ gene.stop - gene.start for gene in fam.genes ]
            matrix.write(sep.join(['"'+fam.name+'"',#1
                                    '"'+fam.namedPartition+'"',#2
                                    '"'+ str(product.most_common(1)[0][0])  +'"',#3
                                    '"' + str(len(fam.organisms)) + '"',#4
                                    '"' + str(len(fam.genes)) + '"',#5
                                    '"' + str(round(len(fam.genes)/len(fam.organisms),2)) + '"',#6
                                    '"NA"',#7
                                    '"NA"',#8
                                    '""',#9
                                    '""',#10
                                    '""',#11
                                    '"' + str(min(l)) + '"',#12
                                    '"' + str(max(l)) + '"',#13
                                    '"' + str(round(sum(l)/len(l),2)) + '"']#14
                                    +genes)+"\n")#15
    logging.getLogger().info(f"Done writing the matrix : '{outname}'")

def writeGenePresenceAbsence(output, compress=False):
    logging.getLogger().info(f"Writing the gene presence absence file ...")
    outname = output + "/gene_presence_absence.Rtab"
    with write_compressed_or_not(outname,compress) as matrix:
        index_org = {}
        default_dat = []
        for index, org in enumerate(pan.organisms):
            default_dat.append('0')
            index_org[org] = index

        matrix.write('\t'.join(['Gene']#14
                                +[str(org) for org in pan.organisms])+"\n")#15
        default_genes =  ["0"] * len(pan.organisms)
        org_index = pan.getIndex()#should just return things
        for fam in pan.geneFamilies:
            genes = default_genes.copy()
            for org in fam.organisms:
                genes[org_index[org]] = "1"

            matrix.write('\t'.join([fam.name]#14
                                    +genes)+"\n")#15
    logging.getLogger().info(f"Done writing the gene presence absence file : '{outname}'")

def writeFlatFiles(pangenome, output, cpu = 1, soft_core = 0.95, csv=False, genePA = False, gexf = False, light_gexf = False, projection = False, compress = False):
    global pan
    pan = pangenome
    processes = []
    if any(x for x in [csv, genePA, gexf, light_gexf, projection]):
        #then it's useful to load the pangenome.
        checkPangenomeInfo(pan, needAnnotations=True, needFamilies=True, needGraph=True)
        if not pan.status["partitionned"] == "Loaded" and (light_gexf or gexf or csv) :
            raise Exception("The provided pangenome has not been partitionned. This is not compatible with any of the following options : --light_gexf, --gexf, --csv")
        pan.getIndex()#make the index because it will be used most likely
        with Pool(processes = cpu) as p:
            if csv:
                processes.append(p.apply_async(func = writeMatrix, args = (',', "csv", output, compress, True)))
            if genePA:
                processes.append(p.apply_async(func = writeGenePresenceAbsence, args = (output, compress)))
            if gexf:
                processes.append(p.apply_async(func = writeGEXF, args = (output, False, soft_core, compress)))
            if light_gexf:
                processes.append(p.apply_async(func = writeGEXF, args = (output, True, soft_core, compress)))
            for process in processes:
                process.get()#get all the results

def launch(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    writeFlatFiles(pangenome, args.output, args.cpu, args.soft_core, args.csv, args.Rtab, args.gexf, args.light_gexf, args.projection, args.compress)


def writeFlatSubparser(subparser):
    parser = subparser.add_parser("write",help = "Writes 'flat' files representing the pangenome that can be used with other softwares", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument("--soft_core",required=False, default = 0.95, help = "Soft core threshold to use")
    optional.add_argument("--gexf",required = False, action = "store_true", help = "write a gexf file with all the annotations and all the genes of each gene family")
    optional.add_argument("--light_gexf",required = False, action="store_true",help = "write a gexf file with the gene families and basic informations about them")
    optional.add_argument("--csv", required=False, action = "store_true",help = "csv file format as used by Roary, among others. The alternative gene ID will be the partition, if there is one")
    optional.add_argument("--Rtab", required=False, action = "store_true",help = "tabular file for the gene binary presence absence matrix")
    optional.add_argument("--projection", required=False, action = "store_true",help = "a csv file for each organism providing informations on the projection of the graph on the organism")
    optional.add_argument("--compress",required=False, action="store_true",help="Compress the files in .gz")
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o','--output', required=True, type=str, help="Output directory where the file(s) will be written")
    return parser