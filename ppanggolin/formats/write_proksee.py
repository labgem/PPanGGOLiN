#!/usr/bin/env python3
# coding:utf-8

# default libraries
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import json
import logging
from pathlib import Path
from random import randint
from tqdm import tqdm
from typing import Dict, List, Tuple, Set
from uuid import uuid4
from itertools import cycle
from collections import defaultdict

# installed libraries
from bokeh.palettes import Category20
from plotly.express.colors import qualitative

from ppanggolin.genome import Organism, Contig, Gene
from ppanggolin.region import Spot, Module

# local libraries
from ppanggolin.pangenome import Pangenome


def palette() -> List[Tuple[int]]:
    palette = qualitative.Vivid + qualitative.Pastel2 + qualitative.Pastel1 + qualitative.Antique  + qualitative.Safe +  qualitative.Bold
    palette = cycle(palette)

    return palette


def read_settings(settings_data: dict):
    if "format" not in settings_data:
        settings_data["format"] = "circular"
    if "geneticCode" not in settings_data:
        # TODO Manage genetic code
        settings_data["geneticCode"] = "11"


def write_legend_items(features: List[str], module_to_color: Dict[Module, str]):#, modules: Set[Module]): #, sources: List[str]):
    """
    
    """
    # use https://medialab.github.io/iwanthue/ to find nice colors 
    # that associate well with established partition colors (orange, light green, light blue)    
    main_colors = {
                    "orange": "#e59c04",
                    "light green": "#00d860" ,
                    "light blue": "#79deff",
                    "purple": "#a567bb",
                    "dark green": "#7a9a4c",
                    "dark red": "#ca5c55",
                    }

    legend_data = {"items" : [
                            {"name": "persistent", "swatchColor": main_colors['orange'], "decoration": "arrow"},
                            {"name": "shell", "swatchColor": main_colors['light green'], "decoration": "arrow"},
                            {"name": "cloud", "swatchColor": main_colors['light blue'], "decoration": "arrow"},
                            {"name": "RNA", "swatchColor": main_colors['purple'], "decoration": "arrow"},
                        ]
                 }
    if "rgp" in features or "all" in features:
        legend_data["items"].append({"name": "RGP", "swatchColor": main_colors['dark green'], "decoration": "arc"}),

    if "modules" in features or "all" in features:
        # if modules is None:
        # legend_data["items"].append({"name": "Module", "swatchColor": main_colors['dark red'], "decoration": "arc"})
        # else:
        for mod, color in sorted(module_to_color.items(), key=lambda x: x[0].ID):
            legend_data["items"].append({"name": f"module_{mod.ID}", "decoration": "arc", "swatchColor": color, "visible":False})

    return legend_data

def write_tracks(features: List[str]):

    tracks = [{"name": "Gene", "separateFeaturesBy": "strand", "position": "outside", "thicknessRatio": 1,
               "dataType": "feature", "dataMethod": "source", "dataKeys": "Gene"}
               ]
    
    if "rgp" in features or "all" in features:
        tracks.append({"name": "RGP", "separateFeaturesBy": "None", "position": "inside", "thicknessRatio": 1,
                       "dataType": "feature", "dataMethod": "source", "dataKeys": "RGP"}),
    # if "spots" in features or "all" in features:
    #     tracks.append({"name": "Spots", "separateFeaturesBy": "None", "position": "inside", "thicknessRatio": 1,
    #                    "dataType": "feature", "dataMethod": "source", "dataKeys": "Spot"})
    if "modules" in features or "all" in features:
        tracks.append({"name": "Module", "separateFeaturesBy": "None", "position": "inside", "thicknessRatio": 1,
                       "dataType": "feature", "dataMethod": "source", "dataKeys": "Module"})

    return tracks


def read_data(template: Path, features: List[str], modules: List[str] = None) -> dict:
    """
    """
    
    with open(template, "r") as template_file:
        proksee_data = json.load(template_file)

    now = datetime.now()

    if "created" in proksee_data["cgview"]:
        proksee_data["cgview"]["updated"] = now.strftime("%Y-%m-%d %H:%M:%S")
        last_version = proksee_data["cgview"]["version"].split('.')
        proksee_data["cgview"]["version"] = ".".join(last_version[:-1] + last_version[-1] + 1)
    else:
        proksee_data["cgview"]["created"] = now.strftime("%Y-%m-%d %H:%M:%S")
        proksee_data["cgview"]["version"] = "1.0"

    if "name" not in proksee_data["cgview"]:
        proksee_data["cgview"]["name"] = "PPanGGOLiN annotations at genome levels"
        proksee_data["cgview"]["id"] = uuid4().hex

    read_settings(proksee_data["cgview"]["settings"])

    if "items" not in proksee_data["cgview"]["legend"]:
        write_legend_items(proksee_data["cgview"]["legend"], features, modules)

    if "tracks" not in proksee_data["cgview"]:
        proksee_data["cgview"]["tracks"] = write_tracks(features)
    return proksee_data

def initiate_proksee_data(features, org_name, module_to_color):
    """

    """

    proksee_legends =  write_legend_items(features, module_to_color)

    proksee_tracks = write_tracks(features)

    proksee_captions =  {
                        "name": f"{org_name} annotated with PPanGGOLiN",
                        "position": "bottom-center",
                        "font": "sans-serif,plain,18",
                        "backgroundColor": "rgba(255,255,255,0.4)"
                        }

    cgview_data = {"name": "PPanGGOLiN annotations at genome levels",
                   "version": "1.5.0",
                   'settings':{},
                   "legend":proksee_legends,
                   "tracks":proksee_tracks,
                   "sequence":{},
                   'captions':[proksee_captions],
                  }

    return {"cgview":cgview_data}


def write_contig(organism: Organism, genome_sequences):
    """
    """
    contigs_data_list = []
    for contig in tqdm(organism.contigs, unit="contig", disable=True):
        contig_info = {"name": contig.name,
                        "length": contig.length,
                        "orientation": "+",
                        }
        if genome_sequences:
            contig_info['seq'] = genome_sequences[contig.name]
        contigs_data_list.append(contig_info)
    return contigs_data_list


def write_genes(organism: Organism, disable_bar=True):
    genes_data_list = []
    gf2gene = defaultdict(list)

    for gene in tqdm(organism.genes, total=organism.number_of_genes(), unit="genes", disable=disable_bar):

        gf = gene.family
        gf2gene[gf.name].append(gene)

        genes_data_list.append({"name": gene.name,
                                "type": "Gene",
                                "contig": gene.contig.name,
                                "start": gene.start,
                                "stop": gene.stop,
                                "strand": 1 if gene.strand == "+" else -1,
                                "product": gene.product,
                                "tags": [gene.family.named_partition, gene.family.name],
                                "source": "Gene",
                                "legend": gene.family.named_partition,
                                "meta": ""#annotations
                                })
        
    for gene in tqdm(organism.rna_genes, total=organism.number_of_rnas(), unit="rnas", disable=disable_bar):

        genes_data_list.append({"name": gene.name,
                                "type": "Gene",
                                "contig": gene.contig.name,
                                "start": gene.start,
                                "stop": gene.stop,
                                "strand": 1 if gene.strand == "+" else -1,
                                "product": gene.product,
                                "tags": [],
                                "source": "Gene",
                                "legend": "RNA",
                                "meta": ""#annotations
                                })
        
    return genes_data_list, gf2gene


def write_partition(organism: Organism):
    partition_data_list = []
    for gene in tqdm(organism.genes, total=organism.number_of_genes(), unit="genes", disable=True):
        partition_data_list.append({"name": gene.family.name,
                                    "presence": gene.family.named_partition,
                                    "contig": gene.contig.name,
                                    "start": gene.start,
                                    "stop": gene.stop,
                                    "source": "partition",
                                    "legend": gene.family.named_partition,
                                    "tags": ["partition"]})
        
    return partition_data_list


def write_rgp(pangenome: Pangenome, organism: Organism):
    rgp_data_list = []
    for rgp in tqdm(pangenome.regions, unit="RGP", disable=True):
        if rgp.organism == organism:
            rgp_data_list.append({"name": rgp.name,
                                  "contig": rgp.contig.name,
                                  "start": rgp.start,
                                  "stop": rgp.stop,
                                  "legend": "RGP",
                                  'source':"RGP",
                                  "tags": []})
    return rgp_data_list


def write_spots(pangenome: Pangenome, organism: Organism, gf2genes: Dict[str, List[Gene]]):
    spots_data_list = []
    for spot in tqdm(pangenome.spots, unit="Spot", disable=True):
        spot: Spot
        spot_orgs = set()

        for gf in spot.families:
            spot_orgs |= set(gf.organisms)

        if organism in spot_orgs:
            gf_intersection = set(organism.families) & set(spot.families)
            completion = round(len(gf_intersection) / spot.number_of_families, 2)

            for gf in gf_intersection:
                for gene in gf2genes[gf.name]:
                    spots_data_list.append({"name": f"Spot_{spot.ID}",
                                            "start": gene.start,
                                            "stop": gene.stop,
                                            "contig": gene.contig.name,
                                            "legend": "Spot",
                                            "source":"Spot",
                                            "tags": [],
                                            "meta": {
                                                "completion": completion
                                            }})
    return spots_data_list


def write_modules(modules: List[Module], organism: Organism, gf2genes: Dict[str, List[Gene]]):
    modules_data_list = []
    for module in modules:
        gf_intersection = set(organism.families) & set(module.families)
        if gf_intersection:
            completion = round(len(gf_intersection) / len(set(module.families)), 2)
            for gf in gf_intersection:
                for gene in gf2genes[gf.name]:
                    modules_data_list.append({"name": f"Module_{module.ID}", 
                                                "presence": "Module",
                                                "start": gene.start,
                                                "stop": gene.stop,
                                                "contig": gene.contig.name,
                                                "legend":f"module_{module.ID}",
                                                "source": "Module",
                                                "tags": [],
                                                "meta": {
                                                    "completion": completion
                                                }})
    return modules_data_list


def write_proksee_organism(pangenome: Pangenome, organism: Organism, output: Path,
                           features: List[str] = None, module_to_colors: Dict[Module,str] = None, genome_sequences= None):
    
    proksee_data = initiate_proksee_data(features, organism.name, module_to_colors)


    proksee_data["cgview"]["sequence"]["contigs"] = write_contig(organism, genome_sequences)


    genes_features, gf2genes = write_genes(organism)

    proksee_data["cgview"]["features"] = genes_features
    # proksee_data["cgview"]["features"] += write_partition(organism)

    if "rgp" in features or "all" in features:
        proksee_data["cgview"]["features"] += write_rgp(pangenome=pangenome, organism=organism)

    # if "spots" in features or "all" in features:
    #     proksee_data["cgview"]["features"] += write_spots(pangenome=pangenome, organism=organism, gf2genes=gf2genes)

    if "modules" in features or "all" in features:
        proksee_data["cgview"]["features"] += write_modules(modules=module_to_colors, organism=organism, gf2genes=gf2genes)

    logging.debug(f"Write proksee for {organism.name}")
    with open(output.joinpath(organism.name).with_suffix(".json"), "w") as out_json:
        json.dump(proksee_data, out_json, indent=2)


# def write_proksee(pangenome: Pangenome, output: Path, features: List[str] = None, sources: List[str] = None,
#                   template: Path = None, organisms_list: List[str] = None, threads: int = 1, disable_bar: bool = False):
#     assert features is not None
#     if template is None:
#         template = Path(__file__).parent.joinpath("proksee_template").with_suffix(".json")
#     if organisms_list is not None:
#         organisms = [organism for organism in pangenome.organisms if organism.name in organisms_list]
#     else:
#         organisms = pangenome.organisms
#     with ThreadPoolExecutor(max_workers=threads) as executor:
#         with tqdm(total=len(organisms), unit='organism', disable=disable_bar) as progress:
#             futures = []
#             for organism in organisms:
#                 future = executor.submit(write_proksee_organism, pangenome, organism, output,
#                                          template, features, sources)
#                 future.add_done_callback(lambda p: progress.update())
#                 futures.append(future)

#             for future in futures:
#                 future.result()
