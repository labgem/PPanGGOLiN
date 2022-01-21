#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging

# installed libraries
from gmpy2 import popcount
from itertools import combinations
from tqdm import tqdm

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import checkPangenomeInfo


def genomes_fluidity(pangenome, disable_bar=False):
    """ Compute the genomes fluidity from the pangenome

    :param pangenome: pangenome which will be used to compute the genomes fluidity
    :type pangenome: Pangenome
    :param disable_bar: Disable the progress bar
    :type disable_bar: bool

    :return: Genomes fluidity value from the pangenome
    :rtype:float
    """

    # check statuses and load info
    logging.getLogger().info("Check information in pangenome")
    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, disable_bar=disable_bar)
    logging.getLogger().debug("Compute binaries sequences corresponding to presence / absence of families in organisms")
    pangenome.compute_org_bitarrays()  # Compute binaries corresponding to presence / absence of families in organisms
    g_sum = 0
    logging.getLogger().debug("Get number of families in each organisms")
    org2_nb_fam = nb_fam_per_org(pangenome, disable_bar)
    logging.getLogger().info("Compute rate of unique family for each genome combination")
    for c_organisms in tqdm(list(combinations(pangenome.organisms, 2)), unit="combination", disable=disable_bar):
        tot_fam = org2_nb_fam.get(c_organisms[0].name) + org2_nb_fam.get(c_organisms[1].name)
        common_fam = popcount(c_organisms[0].bitarray & c_organisms[1].bitarray) - 1
        g_sum += (tot_fam - 2 * common_fam) / tot_fam
    return (2 / (pangenome.number_of_organisms() * (pangenome.number_of_organisms() - 1))) * g_sum


def nb_fam_per_org(pangenome, disable_bar=False):
    """
    Create a dictionary with for each organism the number of gene families

    :param pangenome: Pangenome which contain the organisms and gene families
    :type pangenome: Pangenome
    :param disable_bar: Disable the progress bar
    :type disable_bar: bool

    :return: Dictionary with organisms as key and number of families as value
    :rtype: dict
    """
    org2_nb_fam = dict()
    for org in tqdm(pangenome.organisms, unit='organism', disable=disable_bar):
        org2_nb_fam[org.name] = popcount(org.bitarray)
    return org2_nb_fam


# TODO Function to normalize genome fluidity
# def genomes_fluidity_norm(p_pangenome: Pangenome, disable_bar: bool = False) -> float:

# TODO Function to compute mash distance between genome for normalization
