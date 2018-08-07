#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-

import sys
import gzip
from decimal import Decimal
from collections import defaultdict, OrderedDict
import math
from random import sample
from io import TextIOWrapper
import mmap

""" argument can be a file descriptor (compressed or not) or file path (compressed or not) and return a readable file descriptor"""
def read_compressed_or_not(file_or_file_path):
    file = file_or_file_path
    if type(file) == str:
        file = open(file,"rb")
    else:
        try:
            file = open(file.name,"rb")
        except:
            return(file)
    if file.read(2).startswith(b'\x1f\x8b'):
        file.seek(0)
        if sys.version_info < (3,):# if python2
            return(gzip.open(filename=file.name, mode = "r"))
        else:# if python3
            return(TextIOWrapper(gzip.open(filename=file, mode = "r")))
    else:
        file.close()
        file = open(file.name,"r")
        return(file)

def get_num_lines(file):
    fp = open(file.name, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

""" The number of combinations of n things taken k at a time."""
def comb_k_n(k,n):
    if (k == 0):
        return 1
    if (k > n):
        return 0
    result = 1
    for i in range(0, k):
        result *= float(n - i)/(i + 1);
    return sys.maxsize if result==float("Inf") else int(round(result)) 

# proportional sampling
def samplingCombinations(items, sample_ratio, sample_min, sample_max=100, step = 1):
    samplingCombinationList = defaultdict(list)
    item_size = len(items)
    combTotNb = pow(2,item_size)-1
    for k in range(1, item_size, step):
        tmp_comb = []
        combNb = Decimal(comb_k_n(item_size, k))
        combNb = sys.float_info.max if combNb>sys.float_info.max else combNb # to avoid to reach infinit values
        combNb_sample = math.ceil(Decimal(combNb)/Decimal(sample_ratio))
        # Plus petit echantillonage possible pour un k donn<C3><A9> = sample_min
        if ((combNb_sample < sample_min) and k != item_size):
            combNb_sample = sample_min
        # Plus grand echantillonage possible
        if (sample_max != None and (combNb_sample > sample_max)):
            combNb_sample = sample_max
        
        i = 0;
        while(i < combNb_sample):
            comb_sub = sample(items,k)
            # Echantillonnage sans remise
            if (comb_sub not in tmp_comb):
                tmp_comb.append(comb_sub)
                samplingCombinationList[len(comb_sub)].append(comb_sub)
            i+=1
    return samplingCombinationList

"""simple arithmetic mean"""
def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

"""simple median"""
def median(numbers):
    numbers = sorted(numbers)
    n = len(numbers)
    if n == 0:
        return(None)
    if n%2 == 1:
        return numbers[n//2]
    else:
        i = n//2
        return (numbers[i - 1] + numbers[i])/2

def standard_deviation(lst, population=True):
    """Calculates the standard deviation for a list of numbers.
    from https://codeselfstudy.com/blogs/how-to-calculate-standard-deviation-in-python"""
    num_items = len(lst)
    mean = sum(lst) / num_items
    differences = [x - mean for x in lst]
    sq_differences = [d ** 2 for d in differences]
    ssd = sum(sq_differences)
 
    # Note: it would be better to return a value and then print it outside
    # the function, but this is just a quick way to print out the values along
    # the way.
    if population is True:
        variance = ssd / num_items
    else:
        variance = ssd / (num_items - 1)
    sd = math.sqrt(variance)

    return(sd)
"""insertion of element at the top of an OrderedDict (for compatibility with python 2.7)
from : https://stackoverflow.com/questions/16664874/how-can-i-add-an-element-at-the-top-of-an-ordereddict-in-python/18326914
"""
def ordered_dict_prepend(dct, key, value, dict_setitem=dict.__setitem__):
    root = dct._OrderedDict__root
    first = root[1]

    if key in dct:
        link = dct._OrderedDict__map[key]
        link_prev, link_next, _ = link
        link_prev[1] = link_next
        link_next[0] = link_prev
        link[0] = root
        link[1] = first
        root[1] = first[0] = link
    else:
        root[1] = first[0] = dct._OrderedDict__map[key] = [root, first, key]
        dict_setitem(dct, key, value)