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
    fluidity_dict = {'all': None, 'shell': None, 'cloud': None, 'accessory': None}
    for subset in fluidity_dict.keys():
        logging.getLogger().debug(f"Compute binaries for {subset} partition")
        pangenome.compute_org_bitarrays(part=subset)
        # Compute binaries corresponding to presence / absence of families in organisms
        g_sum = 0
        logging.getLogger().debug("Get number of families in each organisms")
        org2_nb_fam = nb_fam_per_org(pangenome, disable_bar)
        logging.getLogger().info("Compute rate of unique family for each genome combination")
        for c_organisms in tqdm(list(combinations(pangenome.organisms, 2)), unit="combination", disable=disable_bar):
            tot_fam = org2_nb_fam.get(c_organisms[0].name) + org2_nb_fam.get(c_organisms[1].name)
            common_fam = popcount(c_organisms[0].bitarray & c_organisms[1].bitarray) - 1
            if tot_fam > 0 and common_fam > 0:
                g_sum += (tot_fam - 2 * common_fam) / tot_fam
        fluidity_dict[subset] = (2 / (pangenome.number_of_organisms() * (pangenome.number_of_organisms() - 1))) * g_sum
    return fluidity_dict


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

# TODO Create function to compute module fluidity

# TODO Function to compute mash distance between genome for normalization

def fam_fluidity(pangenome, disable_bar=False):
    """ Compute the family fluidity from the pangenome

    :param pangenome: pangenome which will be used to compute the genomes fluidity
    :type pangenome: Pangenome
    :param disable_bar: Disable the progress bar
    :type disable_bar: bool

    :return: family fluidity value from the pangenome
    :rtype:float
    """
    # check statuses and load info
    logging.getLogger().info("Check information in pangenome")
    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, disable_bar=disable_bar)
    fluidity_dict = {'all': None, 'shell': None, 'cloud': None, 'accessory': None}
    for subset in fluidity_dict.keys():
        logging.getLogger().debug(f"Compute binaries for {subset} partition")
        pangenome.computeFamilyBitarrays(part=subset)
        # Compute binaries corresponding to presence / absence of families in organisms
        f_sum = 0
        logging.getLogger().debug("Get number of families in each organisms")
        fam_2_nb_org = nb_org_per_fam(pangenome, disable_bar)
        logging.getLogger().info("Compute rate of unique organism for each family combination")
        for c_fam in tqdm(list(combinations(pangenome.geneFamilies, 2)), unit="combination", disable=disable_bar):
            tot_org = fam_2_nb_org.get(c_fam[0].name) + fam_2_nb_org.get(c_fam[1].name)
            common_fam = popcount(c_fam[0].bitarray & c_fam[1].bitarray) - 1
            if tot_org > 0 and common_fam > 0:
                f_sum += (tot_org - 2 * common_fam) / tot_org
        fluidity_dict[subset] = (2 / (pangenome.number_of_geneFamilies() *
                                      (pangenome.number_of_geneFamilies() - 1))) * f_sum
    return fluidity_dict


def nb_org_per_fam(pangenome, disable_bar=False):
    """
    Create a dictionary with for each gene families the number of organism

    :param pangenome: Pangenome which contain the organisms and gene families
    :type pangenome: Pangenome
    :param disable_bar: Disable the progress bar
    :type disable_bar: bool

    :return: Dictionary with organisms as key and number of families as value
    :rtype: dict
    """
    fam_2_nb_org = dict()
    for fam in tqdm(pangenome.geneFamilies, unit='gene families', disable=disable_bar):
        fam_2_nb_org[fam.name] = popcount(fam.bitarray)
    return fam_2_nb_org
