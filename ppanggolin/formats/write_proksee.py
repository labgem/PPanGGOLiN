#!/usr/bin/env python3
# coding:utf-8

# default libraries
import json
import logging
from pathlib import Path
from tqdm import tqdm
from typing import Dict, List, Tuple
from collections import defaultdict

# installed libraries


# local libraries
from ppanggolin.genome import Organism, Gene
from ppanggolin.region import Module, Region
from ppanggolin.pangenome import Pangenome



def write_legend_items(features: List[str], module_to_color: Dict[Module, str]):
    """
    Generates legend items based on the selected features and module-to-color mapping.

    :param features: A list of features to include in the legend.
    :param module_to_color: A dictionary mapping modules to their assigned colors.

    :return: A data structure containing legend items based on the selected features and module colors.
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
        for mod, color in sorted(module_to_color.items(), key=lambda x: x[0].ID):
            legend_data["items"].append({"name": f"module_{mod.ID}", "decoration": "arc", "swatchColor": color, "visible":False})

    return legend_data

def write_tracks(features: List[str]):
    """
    Generates track information based on the selected features.

    :param features: A list of features to include in the ProkSee data.

    :return: A list of track configurations based on the selected features.
    """
    tracks = [
        {
            "name": "Gene",
            "separateFeaturesBy": "strand",
            "position": "outside",
            "thicknessRatio": 1,
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": "Gene"
        }
    ]

    if "rgp" in features or "all" in features:
        tracks.append({
            "name": "RGP",
            "separateFeaturesBy": "None",
            "position": "inside",
            "thicknessRatio": 1,
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": "RGP"
        })

    if "modules" in features or "all" in features:
        tracks.append({
            "name": "Module",
            "separateFeaturesBy": "None",
            "position": "inside",
            "thicknessRatio": 1,
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": "Module"
        })

    return tracks


def initiate_proksee_data(features: List[str], org_name: str, module_to_color: Dict[Module, str]):
    """
    Initializes ProkSee data structure with legends, tracks, and captions.

    :param features: A list of features to include in the ProkSee data.
    :param org_name: The name of the organism for which the ProkSee data is being generated.
    :param module_to_color: A dictionary mapping modules to their assigned colors.

    :return: ProkSee data structure containing legends, tracks, and captions.
    """
    proksee_legends = write_legend_items(features, module_to_color)
    proksee_tracks = write_tracks(features)

    proksee_captions = {
        "name": f"{org_name} annotated with PPanGGOLiN",
        "position": "bottom-center",
        "font": "sans-serif,plain,18",
        "backgroundColor": "rgba(255,255,255,0.4)"
    }

    cgview_data = {
        "name": "PPanGGOLiN annotation at genome level",
        "version": "1.5.0",
        'settings': {},
        "legend": proksee_legends,
        "tracks": proksee_tracks,
        "sequence": {},
        'captions': [proksee_captions],
    }

    return {"cgview": cgview_data}


def write_contig(organism: Organism, genome_sequences: Dict[str, str] = None) -> List[Dict]:
    """
    Writes contig data for a given organism in proksee format.

    :param organism: The organism for which contig data will be written.
    :param genome_sequences: A dictionary mapping contig names to their DNA sequences (default: None).

    :return: A list of contig data in a structured format.
    """
    contigs_data_list = []

    for contig in tqdm(organism.contigs, unit="contig", disable=True):
        contig_info = {
            "name": contig.name,
            "length": contig.length,
            "orientation": "+",
        }

        if genome_sequences:
            contig_info['seq'] = genome_sequences.get(contig.name, "")

        contigs_data_list.append(contig_info)

    return contigs_data_list



def write_genes(organism: Organism, disable_bar: bool = True) -> Tuple[List[Dict], Dict[str, List[Gene]]]:
    """
    Writes gene data for a given organism, including both protein-coding genes and RNA genes.

    :param organism: The organism for which gene data will be written.
    :param disable_bar: A flag to disable the progress bar when processing genes (default: True).
    
    :return: A tuple containing a list of gene data in a structured format and a dictionary mapping gene families to genes.
    """
    genes_data_list = []
    gf2gene = defaultdict(list)

    # Process protein-coding genes
    for gene in tqdm(organism.genes, total=organism.number_of_genes(), unit="genes", disable=disable_bar):
        gf = gene.family
        gf2gene[gf.name].append(gene)

        genes_data_list.append({
            "name": gene.name,
            "type": "Gene",
            "contig": gene.contig.name,
            "start": gene.start,
            "stop": gene.stop,
            "strand": 1 if gene.strand == "+" else -1,
            "product": gene.product,
            "tags": [gene.family.named_partition, gene.family.name],
            "source": "Gene",
            "legend": gene.family.named_partition,
            "meta": ""  # annotations
        })

    # Process RNA genes
    for gene in tqdm(organism.rna_genes, total=organism.number_of_rnas(), unit="rnas", disable=disable_bar):
        genes_data_list.append({
            "name": gene.name,
            "type": "Gene",
            "contig": gene.contig.name,
            "start": gene.start,
            "stop": gene.stop,
            "strand": 1 if gene.strand == "+" else -1,
            "product": gene.product,
            "tags": [],
            "source": "Gene",
            "legend": "RNA",
            "meta": ""  # annotations
        })

    return genes_data_list, gf2gene


def write_rgp(rgps: Pangenome, organism: Organism):
    """
    Writes RGP (Region of Genomic Plasticity) data for a given organism in proksee format.

    :param pangenome: The pangenome containing information about RGPs.
    :param organism: The specific organism for which RGP data will be written.

    :return: A list of RGP data in a structured format.
    """
    rgp_data_list = []

    # Iterate through each RGP in the pangenome
    for rgp in tqdm(rgps, unit="RGP", disable=True):
        if rgp.organism == organism:
            # Create an entry for the RGP in the data list
            rgp_data_list.append({
                "name": rgp.name,
                "contig": rgp.contig.name,
                "start": rgp.start,
                "stop": rgp.stop,
                "legend": "RGP",
                "source": "RGP",
                "tags": []
            })

    return rgp_data_list


def write_modules(modules: List[Module], organism: Organism, gf2genes: Dict[str, List[Gene]]):
    """
    Writes module data in proksee format for a list of modules associated with a given organism.

    :param modules: A list of modules for which data will be written.
    :param organism: The organism to which the modules are associated.
    :param gf2genes: A dictionary that maps gene families to the genes they contain.

    :return: A list of module data in a structured format.
    """
    modules_data_list = []

    # Iterate through each module and find intersecting gene families
    for module in modules:
        gf_intersection = set(organism.families) & set(module.families)

        if gf_intersection:
            # Calculate the completion percentage
            completion = round(len(gf_intersection) / len(set(module.families)), 2)

            # Create module data entries for genes within intersecting gene families
            for gf in gf_intersection:
                for gene in gf2genes[gf.name]:
                    modules_data_list.append({
                        "name": f"Module_{module.ID}",
                        "presence": "Module",
                        "start": gene.start,
                        "stop": gene.stop,
                        "contig": gene.contig.name,
                        "legend": f"module_{module.ID}",
                        "source": "Module",
                        "tags": [],
                        "meta": {
                            "completion": completion
                        }
                    })

    return modules_data_list


def write_proksee_organism(organism: Organism, output_file: Path,
                           features: List[str] = None,
                           module_to_colors: Dict[Module, str] = None,
                           rgps:List[Region] = None,
                           genome_sequences: Dict[str,str] = None):
    """
    Write ProkSee data for a given organism.

    :param organism: The organism for which ProkSee data will be written.
    :param output_file: The output file where ProkSee data will be written.
    :param features: A list of features to include in the ProkSee data, e.g., ["rgp", "modules", "all"].
    :param module_to_colors: A dictionary mapping modules to their assigned colors.
    :patram rgps: list of RGPs that belong to the organisms
    :param genome_sequences: The genome sequences for the organism.

    This function writes ProkSee data for a given organism, including contig information, genes colored by partition, RGPs,
    and modules. The resulting data is saved as a JSON file in the specified output file.
    """
    proksee_data = initiate_proksee_data(features, organism.name, module_to_colors)

    proksee_data["cgview"]["sequence"]["contigs"] = write_contig(organism, genome_sequences)

    genes_features, gf2genes = write_genes(organism)

    proksee_data["cgview"]["features"] = genes_features

    if "rgp" in features or "all" in features:
        proksee_data["cgview"]["features"] += write_rgp(rgps, organism=organism)

    if "modules" in features or "all" in features:
        proksee_data["cgview"]["features"] += write_modules(modules=module_to_colors, organism=organism, gf2genes=gf2genes)

    logging.debug(f"Write ProkSee for {organism.name}")
    with open(output_file, "w") as out_json:
        json.dump(proksee_data, out_json, indent=2)
