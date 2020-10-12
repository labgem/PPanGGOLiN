#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse
from multiprocessing import Pool
from collections import Counter, defaultdict
import logging
import pkg_resources
from statistics import median, mean, stdev
import os

#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import write_compressed_or_not, mkOutdir, restricted_float
from ppanggolin.formats import checkPangenomeInfo

#global variable to store the pangenome
pan = None

def writeJSONheader(json):
    json.write('{"directed": false, "multigraph": false,')
    json.write(' "graph": {')
    json.write(' "organisms": {')
    orgstr = []
    for org in pan.organisms:
        orgstr.append('"'+org.name+'": {')
        contigstr = []
        for contig in org.contigs:
            contigstr.append('"' + contig.name + '": {"is_circular": ' + ('true' if contig.is_circular else 'false')+'}')
        orgstr[-1] += ', '.join(contigstr) + "}"

    json.write(', '.join(orgstr) + "}")
    ##if other things are to be written such as the parameters, write them here
    json.write('},')

def writeJSONGeneFam(geneFam, json):
        json.write('{'+ f'"id": "{geneFam.name}", "nb_genes": {len(geneFam.genes)}, "partition": "{geneFam.namedPartition}", "subpartition": "{geneFam.partition}"')
        orgDict = {}
        name_counts = Counter()
        product_counts = Counter()
        length_counts = Counter()
        for gene in geneFam.genes:
            name_counts[gene.name] += 1
            product_counts[gene.product] += 1
            length_counts[gene.stop - gene.start] += 1
            try:
                orgDict[gene.organism][gene.contig].append(gene)
            except KeyError:
                try:
                    orgDict[gene.organism][gene.contig] = [gene]
                except KeyError:
                    orgDict[gene.organism] = {gene.contig : [gene]}
        json.write(f', "name": "{name_counts.most_common(1)[0][0]}", "product": "{product_counts.most_common(1)[0][0]}", "length": {length_counts.most_common(1)[0][0]}')
        json.write(', "organisms": {')
        orgstr = []
        for org in orgDict:
            orgstr.append('"' + org.name + '": {')
            contigstr = []
            for contig in orgDict[org]:
                contigstr.append('"' + contig.name + '": {')
                genestr = []
                for gene in orgDict[org][contig]:
                    identifier = gene.ID if gene.local_identifier == "" else gene.local_identifier
                    genestr.append('"' + identifier + '": {' + f'"name": "{gene.name}", "product": "{gene.product}", "is_fragment": {"true" if gene.is_fragment else "false"}, "position": {gene.position}, "strand": "{gene.strand}", "end": {gene.stop}, "start": {gene.start}'+'}')
                contigstr[-1] += ", ".join(genestr) + "}"
            orgstr[-1] += ", ".join(contigstr) + "}"
        json.write(", ".join(orgstr) + "}}")

def writeJSONnodes(json):
    json.write('"nodes": [')
    famList = list(pan.geneFamilies)
    firstFam = famList[0]
    writeJSONGeneFam(firstFam, json)
    for geneFam in famList[1:]:
        json.write(', ')
        writeJSONGeneFam(geneFam, json)
    json.write(']')

def writeJSONedge(edge, json):
    json.write("{")
    json.write(f'"weight": {len(edge.genePairs)}, "source": "{edge.source.name}", "target": "{edge.target.name}"')
    json.write(', "organisms": {')
    orgstr = []
    for org in edge.getOrgDict():
        orgstr.append('"' + org.name + '": [')
        genepairstr = []
        for genepair in  edge.getOrgDict()[org]:
            genepairstr.append('{"source": "' + genepair[0].ID + '", "target": "' + genepair[1].ID + f'", "length": {genepair[0].start - genepair[1].stop}' + '}')
        orgstr[-1] += ', '.join(genepairstr) + ']'
    json.write(', '.join(orgstr) + "}}")

def writeJSONedges(json):
    json.write(', "links": [')
    edgelist = list(pan.edges)
    writeJSONedge(edgelist[0], json)
    for edge in edgelist[1:]:
        json.write(", ")
        writeJSONedge(edge, json)
    json.write(']')

def writeJSON(output, compress):
    logging.getLogger().info("Writing the json file for the pangenome graph...")
    outname = output + "/pangenomeGraph.json"
    with write_compressed_or_not(outname, compress) as json:
        writeJSONheader(json)
        writeJSONnodes(json)
        writeJSONedges(json)
        json.write("}")
    logging.getLogger().info(f"Done writing the json file : '{outname}'")

def writeGEXFheader(gexf, light):
    if not light:
        index = pan.getIndex()#has been computed already
    gexf.write('<?xml version="1.1" encoding="UTF-8"?>\n<gexf xmlns:viz="http://www.gexf.net/1.2draft/viz" xmlns="http://www.gexf.net/1.2draft" version="1.2">\n')
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
                gexf.write(f'          <attvalue for="{index[org]+12}" value="{"|".join([ gene.ID if gene.local_identifier == "" else gene.local_identifier for gene in genes])}" />\n')
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
            alt = fam.namedPartition if fam.partition != "" else False
            genenames = Counter()
            product = Counter()
            for org, gene_list in fam.getOrgDict().items():
                genes[org_index[org]] = " ".join([ '"' + str(gene) + '"' for gene in gene_list]) if geneNames else str(len(gene_list))
                for gene in gene_list:
                    l.append(gene.stop - gene.start)
                    product[gene.product] +=1
                    genenames[gene.name] += 1

            if fam.partition != "":
                alt = fam.namedPartition
            else:
                alt = str(product.most_common(1)[0][0])

            l = [ gene.stop - gene.start for gene in fam.genes ]
            matrix.write(sep.join(['"'+fam.name+'"',#1
                                    '"'+alt+'"',#2
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

def writeStats(output, soft_core, dup_margin, compress=False):
    logging.getLogger().info("Writing pangenome statistics...")
    logging.getLogger().info("Writing statistics on persistent duplication...")
    single_copy_markers = set()#could use bitarrays if speed is needed
    with write_compressed_or_not(output + "/mean_persistent_duplication.tsv", compress) as outfile:
        outfile.write(f"#duplication_margin={round(dup_margin,3)}\n")
        outfile.write("\t".join(["persistent_family","duplication_ratio","mean_presence","is_single_copy_marker"]) + "\n")
        for fam in pan.geneFamilies:
            if fam.namedPartition == "persistent":
                mean_pres = len(fam.genes) / len(fam.organisms)
                nb_multi = 0
                for gene_list in fam.getOrgDict().values():
                    if len(gene_list) > 1:
                        nb_multi +=1
                dup_ratio = nb_multi / len(fam.organisms)
                is_SCM = False
                if dup_ratio < dup_margin:
                    is_SCM = True
                    single_copy_markers.add(fam)
                outfile.write("\t".join([fam.name,
                                         str(round(dup_ratio,3)),
                                         str(round(mean_pres,3)),
                                         str(is_SCM)]) + "\n")
    logging.getLogger().info("Done writing stats on persistent duplication")
    logging.getLogger().info("Writing genome per genome statistics (completeness and counts)...")
    soft = set()#could use bitarrays if speed is needed
    core = set()
    for fam in pan.geneFamilies:
        if len(fam.organisms) >= pan.number_of_organisms() * soft_core:
            soft.add(fam)
        if len(fam.organisms) == pan.number_of_organisms():
            core.add(fam)

    with write_compressed_or_not(output + "/organisms_statistics.tsv", compress) as outfile:
        outfile.write(f"#soft_core={round(soft_core,3)}\n")
        outfile.write(f"#duplication_margin={round(dup_margin,3)}\n")
        outfile.write("\t".join(["organism","nb_families","nb_persistent_families","nb_shell_families","nb_cloud_families","nb_exact_core","nb_soft_core","nb_genes","nb_persistent_genes","nb_shell_genes","nb_cloud_genes","nb_exact_core_genes","nb_soft_core_genes","completeness","nb_single_copy_markers"]) + "\n")

        for org in pan.organisms:
            fams = org.families
            nb_pers = 0
            nb_shell = 0
            nb_cloud = 0
            for fam in fams:
                if fam.namedPartition == "persistent":
                    nb_pers+=1
                elif fam.namedPartition == "shell":
                    nb_shell+=1
                else:
                    nb_cloud+=1

            nb_gene_pers = 0
            nb_gene_shell = 0
            nb_gene_soft = 0
            nb_gene_cloud = 0
            nb_gene_core = 0
            for gene in org.genes:
                if gene.family.namedPartition == "persistent":
                    nb_gene_pers +=1
                elif gene.family.namedPartition == "shell":
                    nb_gene_shell +=1
                else:
                    nb_gene_cloud += 1
                if gene.family in soft:
                    nb_gene_soft+=1
                    if gene.family in core:
                        nb_gene_core+=1
            completeness = "NA"
            if len(single_copy_markers) > 0:
                completeness = round((len(fams & single_copy_markers) / len(single_copy_markers))*100,2)
            outfile.write("\t".join(map(str,[org.name,
                                    len(fams),
                                    nb_pers,
                                    nb_shell,
                                    nb_cloud,
                                    len(core & fams),
                                    len(soft & fams),
                                    org.number_of_genes(),
                                    nb_gene_pers,
                                    nb_gene_shell,
                                    nb_gene_cloud,
                                    nb_gene_core,
                                    nb_gene_soft,
                                    completeness,
                                    len(fams & single_copy_markers)])) + "\n")

    logging.getLogger().info("Done writing genome per genome statistics")

def writeOrgFile(org, output, compress=False):
    with write_compressed_or_not(output + "/" + org.name + ".tsv",compress) as outfile:
        outfile.write("\t".join(["gene","contig","start","stop","strand","ori","family","nb_copy_in_org","partition","persistent_neighbors","shell_neighbors","cloud_neighbors"]) + "\n")
        for contig in org.contigs:
            for gene in contig.genes:
                nb_pers = 0
                nb_shell = 0
                nb_cloud = 0
                for neighbor in gene.family.neighbors:
                    if neighbor.namedPartition == "persistent":
                        nb_pers+=1
                    elif neighbor.namedPartition == "shell":
                        nb_shell+=1
                    else:
                        nb_cloud+=1
                outfile.write("\t".join(map(str,[ gene.ID if gene.local_identifier == "" else gene.local_identifier,
                                        contig.name,
                                        gene.start,
                                        gene.stop,
                                        gene.strand,
                                        "T" if (gene.name.upper() == "DNAA" or gene.product.upper() == "DNAA") else "F",
                                        gene.family.name,
                                        len(gene.family.getGenesPerOrg(org)),
                                        gene.family.namedPartition,
                                        nb_pers,
                                        nb_shell,
                                        nb_cloud
                                        ])) + "\n")

def writeProjections(output, compress=False):
    logging.getLogger().info("Writing the projection files...")
    outdir = output+"/projection"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for org in pan.organisms:
        writeOrgFile(org, outdir, compress)
    logging.getLogger().info("Done writing the projection files")

def writeParts(output, soft_core, compress=False):
    logging.getLogger().info("Writing the list of gene families for each partitions...")
    if not os.path.exists(output + "/partitions"):
        os.makedirs(output + "/partitions")
    partSets = defaultdict(set)
    #initializing key, value pairs so that files exist even if they are empty
    for neededKey in ["undefined","soft_core","exact_core","exact_accessory","soft_accessory","persistent","shell","cloud"]:
        partSets[neededKey] = set()
    for fam in pan.geneFamilies:
        partSets[fam.namedPartition].add(fam.name)
        if fam.partition.startswith("S"):
            partSets[fam.partition].add(fam.name)
        if len(fam.organisms) >= len(pan.organisms) * soft_core:
            partSets["soft_core"].add(fam.name)
            if len(fam.organisms) == len(pan.organisms):
                partSets["exact_core"].add(fam.name)
            else:
                partSets["exact_accessory"].add(fam.name)
        else:
            partSets["soft_accessory"].add(fam.name)
            partSets["exact_accessory"].add(fam.name)

    for key, val in partSets.items():
        currKeyFile = open(output + "/partitions/" + key + ".txt","w")
        if len(val) > 0:
            currKeyFile.write('\n'.join(val) + "\n")
        currKeyFile.close()
    logging.getLogger().info("Done writing the list of gene families for each partition")

def writeGeneFamiliesTSV(output, compress=False):
    logging.getLogger().info("Writing the file providing the association between genes and gene families...")
    outname = output + "/gene_families.tsv"
    with write_compressed_or_not(outname,compress) as tsv:
        for fam in pan.geneFamilies:
            for gene in fam.genes:
                tsv.write("\t".join([fam.name, gene.ID if gene.local_identifier == "" else gene.local_identifier, "F" if gene.is_fragment else ""])+"\n")
    logging.getLogger().info(f"Done writing the file providing the association between genes and gene families : '{outname}'")

def writeRegions(output, compress = False):
    fname = output + "/plastic_regions.tsv"
    with write_compressed_or_not(fname, compress) as tab:
        tab.write("region\torganism\tcontig\tstart\tstop\tgenes\tcontigBorder\twholeContig\n")
        regions = sorted(pan.regions, key = lambda x : (x.organism.name, x.contig.name, x.start))
        for region in regions:
            tab.write('\t'.join(map(str,[region.name, region.organism, region.contig, region.start, region.stop, len(region.genes), region.isContigBorder, region.isWholeContig]))+"\n")

def summarize_spots(spots, output, compress):

    def r_and_s(value):
        """ rounds to dp figures and returns a str of the provided value"""
        if isinstance(value, float):
            return str(round(value,3))
        else:
            return str(value)

    with write_compressed_or_not(output + "/summarize_spots.tsv", compress) as fout:
        fout.write("spot\tnb_rgp\tnb_families\tnb_unique_family_sets\tmean_nb_genes\tstdev_nb_genes\tmax_nb_genes\tmin_nb_genes\n")
        for spot in sorted(spots, key=lambda x : len(x.regions), reverse=True):
            tot_fams = set()
            rgp_list = list(spot.regions)
            len_uniq_content = len(spot.getUniqContent())
            size_list = []
            for rgp in spot.regions:
                tot_fams |= rgp.families
                size_list.append(len(rgp.genes))
            mean_size = mean(size_list)
            stdev_size = stdev(size_list) if len(size_list) > 1 else 0
            max_size = max(size_list)
            min_size = min(size_list)
            fout.write("\t".join(map(r_and_s,[f"spot_{spot.ID}", len(rgp_list), len(tot_fams), len_uniq_content, mean_size,stdev_size,max_size, min_size])) + "\n")
    logging.getLogger().info(f"Done writing spots in : '{output + '/summarize_spots.tsv'}'")

def spot2rgp(spots, output, compress):
    with write_compressed_or_not(output + "/spots.tsv", compress) as fout:
        fout.write("spot_id\trgp_id\n")
        n_spot = 0
        for spot in spots:
            for rgp in spot.regions:
                fout.write(f"spot_{spot.ID}\t{rgp.name}\n")
            n_spot+=1

def writeSpots(output, compress):
    if len(pan.spots) > 0:
        spot2rgp(pan.spots, output, compress)
        summarize_spots(pan.spots, output, compress)

def writeBorders(output, dup_margin, compress):
    multigenics = pan.get_multigenics(dup_margin=dup_margin)
    all_fams = set()
    with write_compressed_or_not(output+"/spot_borders.tsv",compress) as fout:
        fout.write("spot_id\tnumber\tborder1\tborder2\n")
        for spot in sorted(pan.spots, key= lambda x: len(x.regions), reverse=True):
            curr_borders=spot.borders(pan.parameters["spots"]["set_size"], multigenics)
            for c, border in curr_borders:
                famstring1 = ",".join([ fam.name for fam in border[0] ])
                famstring2 = ",".join([ fam.name for fam in border[1]])
                all_fams |= set(border[0])
                all_fams |= set(border[1])
                fout.write(f"{spot.ID}\t{c}\t{famstring1}\t{famstring2}\n")

    with write_compressed_or_not(output + "/border_protein_genes.fasta",compress) as fout:
        for fam in all_fams:
            fout.write(f">{fam.name}\n")
            fout.write(f"{fam.sequence}\n")

def writeFlatFiles(pangenome, output, cpu = 1, soft_core = 0.95, dup_margin = 0.05, csv=False, genePA = False, gexf = False, light_gexf = False, projection = False, stats = False, json = False, partitions=False,regions = False, families_tsv = False, spots = False, borders=False, compress = False):

    if not any(x for x in [csv, genePA, gexf, light_gexf, projection, stats, json, partitions, regions, spots, borders, families_tsv]):
        raise Exception("You did not indicate what file you wanted to write.")

    global pan
    pan = pangenome
    processes = []
    needAnnotations = False
    needFamilies = False
    needGraph = False
    needPartitions = False
    needSpots = False
    needRegions = False

    if csv or genePA or gexf or light_gexf or projection or stats or json or partitions or regions or spots or families_tsv or borders:
        needAnnotations = True 
    if csv or genePA or gexf or light_gexf or projection or stats or json or partitions or regions or spots or families_tsv  or borders:
        needFamilies = True
    if projection or stats or partitions or regions or spots or borders:
        needPartitions = True
    if gexf or light_gexf or json:
        needGraph = True
    if regions or spots or borders:
        needRegions = True
    if spots or borders:
        needSpots = True

    checkPangenomeInfo(pan, needAnnotations=needAnnotations, needFamilies=needFamilies, needGraph=needGraph, needPartitions= needPartitions, needRGP = needRegions, needSpots = needSpots)
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
        if projection:
            processes.append(p.apply_async(func = writeProjections, args = (output, compress)))
        if stats:
            processes.append(p.apply_async(func = writeStats, args = (output, soft_core, dup_margin, compress)))
        if json:
            processes.append(p.apply_async(func = writeJSON, args = (output, compress)))
        if partitions:
            processes.append(p.apply_async(func = writeParts, args = (output, soft_core, compress)))
        if families_tsv:
            processes.append(p.apply_async(func = writeGeneFamiliesTSV, args = (output, compress)))
        if regions:
            processes.append(p.apply_async(func = writeRegions, args = (output, compress)))
        if spots:
            processes.append(p.apply_async(func = writeSpots, args=(output, compress)))
        if borders:
            processes.append(p.apply_async(func=writeBorders, args=(output, dup_margin, compress)))

        for process in processes:
            process.get()#get all the results

def launchFlat(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    writeFlatFiles(pangenome, args.output,cpu= args.cpu, soft_core=args.soft_core,dup_margin= args.dup_margin,csv= args.csv,genePA= args.Rtab, gexf=args.gexf, light_gexf=args.light_gexf, projection=args.projection, stats=args.stats, json=args.json, partitions=args.partitions, regions=args.regions, families_tsv=args.families_tsv, spots=args.spots, borders=args.borders, compress=args.compress)

def writeFlatSubparser(subparser):
    parser = subparser.add_parser("write", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o','--output', required=True, type=str, help="Output directory where the file(s) will be written")
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument("--soft_core",required=False, type=restricted_float, default = 0.95, help = "Soft core threshold to use")
    optional.add_argument("--dup_margin", required=False, type=restricted_float, default=0.05, help = "minimum ratio of organisms in which the family must have multiple genes for it to be considered 'duplicated'")
    optional.add_argument("--gexf",required = False, action = "store_true", help = "write a gexf file with all the annotations and all the genes of each gene family")
    optional.add_argument("--light_gexf",required = False, action="store_true",help = "write a gexf file with the gene families and basic informations about them")
    optional.add_argument("--csv", required=False, action = "store_true",help = "csv file format as used by Roary, among others. The alternative gene ID will be the partition, if there is one")
    optional.add_argument("--Rtab", required=False, action = "store_true",help = "tabular file for the gene binary presence absence matrix")
    optional.add_argument("--projection", required=False, action = "store_true",help = "a csv file for each organism providing informations on the projection of the graph on the organism")
    optional.add_argument("--stats",required=False, action = "store_true",help = "tsv files with some statistics for each organism and for each gene family")
    optional.add_argument("--partitions", required=False, action = "store_true", help = "list of families belonging to each partition, with one file per partitions and one family per line")
    optional.add_argument("--compress",required=False, action="store_true",help="Compress the files in .gz")
    optional.add_argument("--json", required=False, action = "store_true", help = "Writes the graph in a json file format")
    optional.add_argument("--regions", required=False, action = "store_true", help = "Write the RGP in a tab format, one file per genome")
    optional.add_argument("--spots", required=False, action = "store_true", help = "Write spot summary and a list of all rgp in each spot")
    optional.add_argument("--borders", required=False, action = "store_true", help = "List all borders of each spot")
    optional.add_argument("--families_tsv", required=False, action = "store_true", help = "Write a tsv file providing the association between genes and gene families")
    return parser
