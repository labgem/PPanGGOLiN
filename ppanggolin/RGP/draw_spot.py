#!/usr/bin/env python3
#coding:utf-8


#default libraries
import argparse
import time
import os
import logging
from collections import defaultdict
import random
from math import pi

#local libraries
from ppanggolin.utils import mkOutdir, jaccard_similarities
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import checkPangenomeInfo
from ppanggolin.RGP.spot import compBorder

#installed libraries
from scipy.spatial.distance import pdist
from scipy.sparse import csc_matrix
from scipy.cluster.hierarchy import linkage, dendrogram
from numpy import percentile

from tqdm import tqdm
from bokeh.plotting import ColumnDataSource, figure, show, save
from bokeh.io import output_file
from bokeh.layouts import column, row
from bokeh.models import WheelZoomTool, LabelSet, Slider, CustomJS, HoverTool, RadioGroup, Div

def checkPredictedSpots(pangenome):
    """ checks pangenome status and .h5 files for predicted spots, raises an error if they were not predicted"""
    if pangenome.status["spots"] == "No":
        raise Exception("You are trying to draw spots for a pangenome that does not have spots predicted. Please see the 'spot' subcommand.")


def makeColorsForFams(fams):
    """randomly picks 256 colors for gene families"""
    famcol = {}
    for fam in fams:
        col =  list(random.choices(range(256), k=3))
        famcol[fam] = '#%02x%02x%02x' % (col[0], col[1], col[2])
    return famcol

def orderGeneLists(geneLists, ordered_counts, overlapping_match, exact_match, set_size):
    geneLists = lineOrderGeneLists(geneLists, overlapping_match, exact_match, set_size)
    return rowOrderGeneLists(geneLists, ordered_counts)

def rowOrderGeneLists(geneLists, ordered_counts):
    famDict = defaultdict(set)

    for index, genelist in enumerate([genelist[0] for genelist in geneLists]):
        for gene in genelist:
            if hasattr(gene,"family"):
                famDict[gene.family].add(index)
    all_indexes = []
    all_columns = []
    data = []
    for famIndex, RGPindexes in enumerate(famDict.values()):
        all_indexes.extend([famIndex] * len(RGPindexes))
        all_columns.extend(RGPindexes)
        data.extend([1.0]*len(RGPindexes))

    mat_p_a = csc_matrix((data, (all_indexes,all_columns)), shape = (len(famDict),len(geneLists)), dtype='float')
    dist    = pdist(1 - jaccard_similarities(mat_p_a,0).todense())
    hc      = linkage(dist, 'single')

    dendro = dendrogram(hc,no_plot=True)

    new_geneLists = [ geneLists[index] for index in dendro["leaves"]]
    new_ordered_counts = [ ordered_counts[index] for index in dendro["leaves"] ]
    return new_geneLists, new_ordered_counts

def lineOrderGeneLists(geneLists, overlapping_match, exact_match, set_size):
    classified = set([0])#first gene list has the right order
    new_classify = set()

    to_classify = set(range(1, len(geneLists)))#the others may (or may not) have it

    while len(to_classify) != 0:
        for classIndex in classified:
            base_border1 = [ gene.family for gene in geneLists[classIndex][1][0] ]
            base_border2 = [ gene.family for gene in geneLists[classIndex][1][1] ]
            for unclassIndex in list(to_classify):
                border1 = [ gene.family for gene in geneLists[unclassIndex][1][0] ]
                border2 = [ gene.family for gene in geneLists[unclassIndex][1][1] ]
                if compBorder(base_border1, border1, overlapping_match, exact_match, set_size) and compBorder(base_border2, border2, overlapping_match, exact_match, set_size):
                    to_classify.discard(unclassIndex)
                    new_classify.add(unclassIndex)
                elif compBorder(base_border2, border1, overlapping_match, exact_match, set_size) and compBorder(base_border1, border2, overlapping_match, exact_match, set_size):
                    geneLists[unclassIndex][0] = geneLists[unclassIndex][0][::-1]#reverse the order of the genes to match the 'reference'
                    #inverse the borders
                    former_border_1 = geneLists[unclassIndex][1][0]
                    former_border_2 = geneLists[unclassIndex][1][1]
                    geneLists[unclassIndex][1][0] = former_border_2
                    geneLists[unclassIndex][1][1] = former_border_1

                    #specify the new 'classified' and remove from unclassified
                    to_classify.discard(unclassIndex)
                    new_classify.add(unclassIndex)
        classified = new_classify#the newly classified will help to check the unclassified, the formerly classified are not useful for what remains (if something remains)
        new_classify = set()
    return geneLists


def mkSourceData(genelists, famCol):

    partitionColors = {"shell": "#00D860", "persistent":"#F7A507", "cloud":"#79DEFF"}

    df = {'name':[],'ordered':[],'strand':[],"start":[],"stop":[],'x':[],'y':[],'width':[], 'family_color':[], 'partition_color':[], 'partition':[], "family":[], "product":[],"x_label":[],"y_label":[], "label":[], "gene_type":[],'gene_ID':[], "gene_local_ID":[]}

    for index, GeneList in enumerate(genelists):
        genelist = GeneList[0]

        if genelist[0].start < genelist[1].start:
            ordered=True
            start = genelist[0].start
        else:
            ordered = False
            start = genelist[0].stop

        for gene in genelist:
            df["ordered"].append(str(ordered))
            df["strand"].append(gene.strand)
            df["start"].append(gene.start)
            df["stop"].append(gene.stop)
            df["gene_type"].append(gene.type)
            df["product"].append(gene.product)
            df["gene_local_ID"].append(gene.local_identifier)
            df['gene_ID'].append(gene.ID)
            if "RNA" in gene.type:
                df["name"].append(gene.product)
                df["family"].append(gene.type)
                df["partition"].append("none")
                df["family_color"].append("#A109A7")
                df["partition_color"].append("#A109A7")
                #mod and mod color
            else:
                df["name"].append(gene.name)
                df["family"].append(gene.family.name)
                df["partition"].append(gene.family.namedPartition)
                df["family_color"].append(famCol[gene.family])
                df["partition_color"].append(partitionColors[gene.family.namedPartition])
                #mod and mod color
            df["x"].append((abs(gene.start - start) + abs(gene.stop - start)) / 2)
            df["width"].append(gene.stop - gene.start)
            df["x_label"].append(df["x"][-1] - int(df["width"][-1]/2))
            if ordered:
                if gene.strand == "+":
                    df["y"].append((index * 10) + 1)
                else:
                    
                    df["y"].append((index * 10) - 1)
            else:
                if gene.strand == "+":
                    df["y"].append((index * 10) - 1)
                else:
                    df["y"].append((index * 10) + 1)
            df["y_label"].append(df["y"][-1] + 1.5)
    df["label"] = df["name"]
    df["line_color"] = df["partition_color"]
    df["fill_color"] = df["family_color"]

    TOOLTIPS = [
        ("start","@start"),
        ("stop","@stop"),
        ("name", "@name"),
        ("product","@product"),
        ("family","@family"),
        ("partition","@partition"),
        ("local identifier","@gene_local_ID"),
        ("gene ID","@gene_ID"),
        ("ordered","@ordered"),
        ("strand","@strand"),
    ]

    return ColumnDataSource(data=df), TOOLTIPS

def changeColors(recs, sourceData):
    def colorStr(color_element):
        """javascript code to switch between partition, family and module color for the given 'color_element'"""
        return f"""
            if(this.active == 0){{
                source.data['{color_element}'] = source.data['partition_color'];
            }}else if(this.active == 1){{
                source.data['{color_element}'] = source.data['family_color'];
            }}else if(this.active == 2){{
                source.data['{color_element}'] = source.data['module_color'];
            }}
            recs.{color_element} = source.data['{color_element}'];
            source.change.emit();
        """

    radio_line_color = RadioGroup(labels=["partition", "family", "module"], active=0)
    radio_fill_color = RadioGroup(labels=["partition", "family", "module"], active=1)

    radio_line_color.js_on_click(CustomJS(args=dict(recs=recs,source=sourceData),
        code=colorStr("line_color")))
    
    radio_fill_color.js_on_click(CustomJS(args=dict(recs=recs,source=sourceData),
        code=colorStr("fill_color")))

    color_header = Div(text="<b>Colors:</b>")
    line_title =  Div(text="""Color to use for gene outlines:""",
        width=200, height=100)
    fill_title = Div(text="""Color to fill genes with:""",
        width=200, height=100)
    color_block = column(color_header, row(column(line_title, radio_line_color), column(fill_title, radio_fill_color)))

    return color_block

def addLabels(fig, sourceData):

    labels = LabelSet(x='x_label', y='y_label', text='label', source=sourceData, render_mode='canvas', text_font_size="18px")
    slider_font = Slider(start=0, end=64, value=16, step=1, title="Label font size in px")
    slider_angle = Slider(start=0, end=pi/2, value=0, step=0.01, title="Label angle in radian")

    radio_label_type = RadioGroup(labels=["name", "product", "family","local identifier","gene ID", "none"], active=0)

    slider_angle.js_link('value',labels,'angle')

    slider_font.js_on_change('value',
        CustomJS(args=dict(other=labels),
                code="other.text_font_size = this.value+'px'"
        )
    )

    radio_label_type.js_on_click(CustomJS(args=dict(other=labels, source=sourceData),
                code="""
                if(this.active == 5){
                    source.data['label'] = [];
                    for(var i=0;i<source.data['name'].length;i++){
                        source.data['label'].push('');
                    }
                }else if(this.active == 3){
                    source.data['label'] = source.data['gene_local_ID'];
                }else if(this.active == 4){
                    source.data['label'] = source.data['gene_ID'];
                }
                else{
                    source.data['label'] = source.data[this.labels[this.active]];
                }
                other.source = source;
                source.change.emit();
                """
        ))

    label_header = Div(text="<b>Labels:</b>")
    radio_title =  Div(text="""Labels to use:""",
    width=200, height=100)
    labels_block = column(label_header, row(slider_font, slider_angle), column(radio_title, radio_label_type))
    
    fig.add_layout(labels)
    
    return labels_block

def mkGenomes(geneLists, ordered_counts):
    df = {"name":[], "width":[], "occurrences":[],'x':[],'y':[]}

    for index, GeneList in enumerate(geneLists):
        genelist = GeneList[0]
        df["occurrences"].append(ordered_counts[index])

        df["y"].append(index*10)
        df["width"].append(genelist[-1].stop - genelist[0].start)
        if genelist[0].start < genelist[1].start:
            df["x"].append((df["width"][-1])/2)
        else:
            df["x"].append((df["width"][-1])/2)
        df["name"].append(genelist[0].organism.name)
    TOOLTIP = [
        ("name","@name"),
        ("occurrences","@occurrences"),
    ]
    return ColumnDataSource(data=df), TOOLTIP

def drawCurrSpot(genelists, ordered_counts, famCol, filename):
    #prepare the source data

    output_file(filename + "_test.html")

    #generate the figure and add some tools to it
    wheel_zoom = WheelZoomTool()
    fig = figure(title="spot graphic", plot_width=1600, plot_height=600, tools=["pan","box_zoom","reset","save",wheel_zoom,"ywheel_zoom","xwheel_zoom"])
    fig.axis.visible = True
    fig.toolbar.active_scroll = wheel_zoom

    #genome rectangles
    genomeSource, genomeTooltip = mkGenomes(genelists, ordered_counts)
    genomeRecs = fig.rect(x='x',y='y',fill_color="dimgray",width="width", height=0.5, source=genomeSource)
    genomeRecs_hover = HoverTool(renderers=[genomeRecs], tooltips=genomeTooltip)
    fig.add_tools(genomeRecs_hover)
    #genome labels?

    #gene rectanges
    GeneSource, GeneTooltips = mkSourceData(genelists, famCol)
    recs = fig.rect(x='x', y='y',line_color='line_color', fill_color='fill_color', width='width',height = 2,line_width=5, source=GeneSource)
    recs_hover = HoverTool(renderers=[recs], tooltips=GeneTooltips)
    fig.add_tools(recs_hover)
    #color modification tools
    color_change_tools = changeColors(recs, GeneSource)

    #gene labels and label modification tools
    labels_change_tools = addLabels(fig, GeneSource)

    save(column(fig, row(labels_change_tools, color_change_tools)))

def drawSelectedSpots(selected_spots, pangenome, output, overlapping_match, exact_match, set_size, disable_bar):
    logging.getLogger().info("Selecting and ordering genes among regions...")

    multigenics = pangenome.get_multigenics(pangenome.parameters["RGP"]["dup_margin"])
    #bar = tqdm(range(len(selected_spots)), unit = "spot", disable = disable_bar)
    print(pangenome.modules)
    fam2mod = { (fam, f"module_{mod.ID}") for mod in pangenome.modules for fam in mod.families }
    print(fam2mod)
    for spot in selected_spots:

        fname = output + '/spot_' + str(spot.ID)

        ##write identical rgps and the rgps they are identical to
        uniqRGPS = set()
        out_struc = open(fname + '_identical_rgps.tsv','w')
        out_struc.write('representative_rgp\trepresentative_rgp_organism\tidentical_rgp\tidentical_rgp_organism\n')
        for keyRGP, otherRGPs in spot.getUniq2RGP().items():
            uniqRGPS.add(keyRGP)
            for rgp in otherRGPs:
                out_struc.write(f"{keyRGP.name}\t{keyRGP.organism.name}\t{rgp.name}\t{rgp.organism.name}\n")
        out_struc.close()

        Fams = set()
        GeneLists = []

        countUniq = spot.countUniqOrderedSet()

        #order unique rgps by occurrences
        sortedUniqRGPs = sorted(uniqRGPS, key = lambda x : countUniq[x], reverse=True)
        for rgp in sortedUniqRGPs:
            borders = rgp.getBorderingGenes(set_size, multigenics)
            minpos = min([ gene.position for border in borders for gene in border ])
            maxpos = max([ gene.position for border in borders for gene in border ])
            GeneList = rgp.contig.genes[minpos:maxpos+1]
            minstart = min([ gene.start for border in borders for gene in border ])
            maxstop = max([ gene.stop for border in borders for gene in border ])
            RNAstoadd = set()
            for rna in rgp.contig.RNAs:
                if rna.start > minstart and rna.start < maxstop:
                    RNAstoadd.add(rna)
            GeneList.extend(RNAstoadd)
            GeneList = sorted(GeneList, key = lambda x : x.start)

            Fams |= { gene.family for gene in GeneList if gene.type == "CDS"}

            GeneLists.append([GeneList, borders, rgp])
        famcolors = makeColorsForFams(Fams)
        ordered_counts = sorted(countUniq.values(), reverse = True)
        GeneLists, ordered_counts = orderGeneLists(GeneLists, ordered_counts, overlapping_match, exact_match, set_size)
        fname = output + '/spot_' + str(spot.ID)

        drawCurrSpot(GeneLists, ordered_counts, famcolors, fname)


def drawSpots(pangenome, output, spot_list, disable_bar):
    #check that the pangenome has spots
    checkPredictedSpots(pangenome)

    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, needGraph=False, needPartitions = True,
    needRGP=True, needSpots=True, disable_bar=disable_bar)

    selected_spots = set()
    curated_spot_list = ['spot_' + str(s) if 'spot' not in s else str(s) for s in spot_list.split(',')]
    
    if spot_list == 'all' or any(x == 'all' for x in curated_spot_list):
        selected_spots = [s for s in pangenome.spots if len(s.getUniqOrderedSet()) > 1]
    else:
        selected_spots = [ s for s in pangenome.spots if "spot_" + str(s.ID) in curated_spot_list]
    if len(selected_spots) < 10:
        logging.getLogger().info(f"Drawing the following spots: {','.join(['spot_' + str(s.ID) for s in selected_spots])}")
    else:
        logging.getLogger().info(f"Drawing {len(selected_spots)} spots")

    multi = pangenome.get_multigenics(pangenome.parameters["RGP"]["dup_margin"])

    drawSelectedSpots(selected_spots, pangenome, output,
                overlapping_match = pangenome.parameters["spots"]["overlapping_match"],
                exact_match = pangenome.parameters["spots"]["exact_match"],
                set_size = pangenome.parameters["spots"]["set_size"],
                disable_bar = disable_bar)


def launch(args):
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    if args.spots == "":
        raise Exception("You did not provide any spot to draw. see the '--spots' options.")
    mkOutdir(args.output, args.force)
    drawSpots(pangenome=pangenome, output = args.output, spot_list=args.spots, disable_bar=args.disable_prog_bar)

def drawSpotSubparser(subparser):
    parser = subparser.add_parser("drawspot", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('-o','--output', required=False, type=str, default="ppanggolin_output"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory")
    optional.add_argument("--spots", required=False, type=str, default='', help =  "a comma-separated list of spots to draw (or 'all' to draw all spots)")
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    return parser