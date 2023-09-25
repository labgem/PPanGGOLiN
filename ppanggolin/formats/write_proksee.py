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
from typing import Dict, List, Tuple
from uuid import uuid4

# installed libraries
from bokeh.palettes import Category20
from ppanggolin.genome import Organism, Contig, Gene
from ppanggolin.region import Spot

# local libraries
from ppanggolin.pangenome import Pangenome


def palette() -> List[Tuple[int]]:
    palette = []
    for hex_color in list(Category20[20]):
        palette.append(tuple(int(hex_color.strip('#')[i:i + 2], 16) for i in (0, 2, 4)))
    return palette


def read_settings(settings_data: dict):
    if "format" not in settings_data:
        settings_data["format"] = "circular"
    if "geneticCode" not in settings_data:
        # TODO Manage genetic code
        settings_data["geneticCode"] = "11"


def write_legend_items(legend_data: dict, features: List[str]): #, sources: List[str]):
    colors = palette()
    legend_data["items"] = [#{"name": "CDS", "swatchColor": f"rgba({','.join(map(str, colors.pop(1)))},0.5)", "decoration": "arrow"},
                            {"name": "persistent", "swatchColor": "rgba(229,156,4,1)", "decoration": "arrow"},
                            {"name": "shell", "swatchColor": "rgba(60,254,91,1)", "decoration": "arrow"},
                            {"name": "cloud", "swatchColor": f"rgba({','.join(map(str, colors.pop(17)))},1)", "decoration": "arrow"},
                             {"name": "RNA", "swatchColor": "rgba(137,23,207,0.5)", "decoration": "arrow"},]
    if "rgp" in features or "all" in features:
        legend_data["items"].append({"name": "RGP", "swatchColor": f"rgba({','.join(map(str, colors.pop(6)))}, 1)", "decoration": "arc"}),
    if "spots" in features or "all" in features:
        legend_data["items"].append({"name": "Spot", "swatchColor": f"rgba({','.join(map(str, colors.pop(5)))}, 1)", "decoration": "arc"})
    if "modules" in features or "all" in features:
        legend_data["items"].append({"name": "Module", "swatchColor": f"rgba({','.join(map(str, colors.pop(3)))},1)", "decoration": "arc"})
    # if "systems" in features or "all" in features:
    #     for source in sources:
    #         color = ','.join(map(str, colors.pop(randint(0, len(colors) - 1))))
    #         legend_data["items"].append({"name": source, "decoration": "arc", "swatchColor": f"rgba({color},1)"})


def write_tracks(features: List[str]):

    tracks = [{"name": "Gene", "separateFeaturesBy": "None", "position": "outside", "thicknessRatio": 1,
               "dataType": "feature", "dataMethod": "source", "dataKeys": "Gene"},
              {"name": "Partition", "separateFeaturesBy": "strand", "position": "outside", "thicknessRatio": 1,
               "dataType": "feature", "dataMethod": "source", "dataKeys": "partition"}]
    
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


def read_data(template: Path, features: List[str], sources: List[str] = None) -> dict:
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
        write_legend_items(proksee_data["cgview"]["legend"], features)

    if "tracks" not in proksee_data["cgview"]:
        proksee_data["cgview"]["tracks"] = write_tracks(features)
    return proksee_data


def write_contig(organism: Organism):
    contigs_data_list = []
    for contig in tqdm(organism.contigs, unit="contig", disable=True):

        contigs_data_list.append({"name": contig.name,
                                  "length": contig.length,
                                  "orientation": "+", # "seq": "".join([gene.dna for gene in contig.genes])
                                  })
    return contigs_data_list


def write_genes(organism: Organism, sources: List[str]=None):
    genes_data_list = []
    gf2gene = {}

    for gene in tqdm(organism.genes, total=organism.number_of_genes(), unit="genes", disable=True):

        gf = gene.family
        if gf.name in gf2gene:
            gf2gene[gf.name].append(gene)
        else:
            gf2gene[gf.name] = [gene]
        # annotations = {source: "|".join(list(map(str, gf.get_source(source)))) for source in gf.sources if
        #                source in sources}
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
        
    for gene in tqdm(organism.rna_genes, total=organism.number_of_rnas(), unit="rnas", disable=True):

        if gf.name in gf2gene:
            gf2gene[gf.name].append(gene)
        else:
            gf2gene[gf.name] = [gene]
        # annotations = {source: "|".join(list(map(str, gf.get_source(source)))) for source in gf.sources if
        #                source in sources}
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
    c=0
    for gene in tqdm(organism.genes, total=organism.number_of_genes(), unit="genes", disable=True):
        c += 1
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


def write_modules(pangenome: Pangenome, organism: Organism, gf2genes: Dict[str, List[Gene]]):
    modules_data_list = []
    for module in tqdm(pangenome.modules, unit="Module", disable=True):
        mod_orgs = set()
        for gf in module.families:
            mod_orgs |= set(gf.organisms)
        if organism in mod_orgs:
            gf_intersection = set(organism.families) & set(module.families)
            completion = round(len(gf_intersection) / len(set(module.families)), 2)
            for gf in gf_intersection:
                for gene in gf2genes[gf.name]:
                    modules_data_list.append({"name": f"Module_{module.ID}", 
                                              "presence": "Module",
                                              "start": gene.start,
                                              "stop": gene.stop,
                                              "contig": gene.contig.name,
                                              "legend": "Module",
                                              "source": "Module",
                                              "tags": [],
                                              "meta": {
                                                  "completion": completion
                                              }})
    return modules_data_list


def write_proksee_organism(pangenome: Pangenome, organism: Organism, output: Path, template: Path,
                           features: List[str] = None, sources: List[str] = None):
    
    proksee_data = read_data(template=template, features=features, sources=sources)
    
    if "name" not in proksee_data["cgview"]["captions"]:
        proksee_data["cgview"]["captions"][0]["name"] = f"{organism.name} annotated with PPanGGOLiN"

    proksee_data["cgview"]["sequence"]["contigs"] = write_contig(organism)

    if "features" not in proksee_data["cgview"]:
        proksee_data["cgview"]["features"] = []

    genes_features, gf2genes = write_genes(organism, sources=sources)
    print(len(genes_features))
    proksee_data["cgview"]["features"] += genes_features
    proksee_data["cgview"]["features"] += write_partition(organism)

    if "rgp" in features or "all" in features:
        proksee_data["cgview"]["features"] += write_rgp(pangenome=pangenome, organism=organism)
    if "spots" in features or "all" in features:
        proksee_data["cgview"]["features"] += write_spots(pangenome=pangenome, organism=organism, gf2genes=gf2genes)
    if "modules" in features or "all" in features:
        proksee_data["cgview"]["features"] += write_modules(pangenome=pangenome, organism=organism, gf2genes=gf2genes)

    logging.debug(f"Write proksee for {organism.name}")
    with open(output.joinpath(organism.name).with_suffix(".json"), "w") as out_json:
        json.dump(proksee_data, out_json, indent=2)


def write_proksee(pangenome: Pangenome, output: Path, features: List[str] = None, sources: List[str] = None,
                  template: Path = None, organisms_list: List[str] = None, threads: int = 1, disable_bar: bool = False):
    assert features is not None
    if template is None:
        template = Path(__file__).parent.joinpath("proksee_template").with_suffix(".json")
    if organisms_list is not None:
        organisms = [organism for organism in pangenome.organisms if organism.name in organisms_list]
    else:
        organisms = pangenome.organisms
    with ThreadPoolExecutor(max_workers=threads) as executor:
        with tqdm(total=len(organisms), unit='organism', disable=disable_bar) as progress:
            futures = []
            for organism in organisms:
                future = executor.submit(write_proksee_organism, pangenome, organism, output,
                                         template, features, sources)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                future.result()
