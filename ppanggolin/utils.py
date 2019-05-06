#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-

import sys
import gzip
import json
from decimal import Decimal
import contextlib
from collections import defaultdict, OrderedDict
import math
import random
from io import TextIOWrapper
import mmap
import numpy
from scipy.stats import chi2_contingency
import pdb

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
        result *= float(n - i)/(i + 1)
    return sys.maxsize if result==float("Inf") else int(round(result)) 

# proportional sampling
def samplingCombinations(items, sample_ratio, sample_min, sample_max=100, step = 1, seed=42):
    random.seed(seed)
    samplingCombinationList = defaultdict(list)
    item_size = len(items)
    for k in range(1, item_size, step):
        tmp_comb = []
        combNb = Decimal(comb_k_n(k,item_size))
        combNb = sys.float_info.max if combNb>sys.float_info.max else combNb # to avoid to reach infinit values
        combNb_sample = math.ceil(Decimal(combNb)/Decimal(sample_ratio))
        # Plus petit echantillonage possible pour un k donné = sample_min
        if ((combNb_sample < sample_min) and k != item_size):
            combNb_sample = sample_min
        # Plus grand echantillonage possible
        if (sample_max != None and (combNb_sample > sample_max)):
            combNb_sample = sample_max
        i = 0
        while(i < combNb_sample):
            comb_sub = random.sample(items,k)
            # Echantillonnage sans remise
            if (comb_sub not in tmp_comb):
                tmp_comb.append(comb_sub)
                samplingCombinationList[len(comb_sub)].append(comb_sub)
            i+=1
    return samplingCombinationList

"""simple arithmetic mean"""
def mean(numbers):
    try:
        return float(sum(numbers)) / len(numbers)
    except ZeroDivisionError:
        return(float("nan"))
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
    if num_items>0:
        mean = sum(lst) / num_items
        differences = [x - mean for x in lst]
        sq_differences = [d ** 2 for d in differences]
        ssd = sum(sq_differences)
        if population is True:
            variance = ssd / num_items
        else:
            variance = ssd / (num_items - 1)
        sd = math.sqrt(variance)
        return(sd)
    else:
        return(float("nan"))

def seq(start, stop, step=1):
    n = int(round((stop - start)/float(step)))
    if n > 1:
        return([start + step*i for i in range(n+1)])
    elif n == 1:
        return([start])
    else:
        return([])

def cramers_corrected_stat(confusion_matrix):
    """ calculate Cramers V statistic for categorical-categorical association.
        uses correction from Bergsma Wicher, 
        Journal of the Korean Statistical Society 42 (2013): 323-328
        source: https://stackoverflow.com/questions/51859894/how-to-plot-a-cramer-s-v-heatmap-for-categorical-features
    """
    chi2 = chi2_contingency(confusion_matrix)
    n = confusion_matrix.sum().sum()
    phi2 = chi2[0]/n
    r,k = confusion_matrix.shape
    phi2corr = max(0, phi2 - ((k-1)*(r-1))/(n-1))
    rcorr = r - ((r-1)**2)/(n-1)
    kcorr = k - ((k-1)**2)/(n-1)
    return (tuple([chi2[1], numpy.sqrt(phi2corr / min( (kcorr-1), (rcorr-1)))]))

def calculate_BIC(log_likelihood,nb_params,nb_points):
    return( log_likelihood - 0.5 *(math.log(nb_points) * nb_params))

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

def jaccard_similarities(mat,jaccard_similarity_th):
    cols_sum = mat.getnnz(axis=0)
    ab = mat.T * mat
    # for rows
    aa = numpy.repeat(cols_sum, ab.getnnz(axis=0))
    # for columns
    bb = cols_sum[ab.indices]
    similarities = ab.copy()
    similarities.data /= (aa + bb - ab.data)
    similarities.data[similarities.data<jaccard_similarity_th] = 0
    similarities.eliminate_zeros()
    return similarities

@contextlib.contextmanager
def empty_cm():
    yield None

def diff_col(c1,c2,min_diff=1):
    r_diff = abs(c1["r"]-c2["r"])
    g_diff = abs(c1["g"]-c2["g"])
    b_diff = abs(c1["b"]-c2["b"])
    if (r_diff+g_diff+b_diff)<min_diff:
        return(False)
    else:
        return(True)
def colors(n,min_diff=10, maxiter=float("inf"), except_list = []):
    ret = []
    iter = 0
    while len(ret) < n and iter<maxiter:
        r = int(random.random() * 256)
        g = int(random.random() * 256)
        b = int(random.random() * 256)
        new_element = {'r': r, 'g': g, 'b': b, 'a': 0}
        valid = True
        for previous_element in ret:
            if not diff_col(new_element,previous_element, min_diff):
                valid = False
                break
        for previous_element in except_list:
            if not diff_col(new_element,previous_element, min_diff):
                valid = False
                break
        if valid:
            ret.append({'r': r, 'g': g, 'b': b, 'a': 0}) 
        iter += 1
    return ret

def average_color(c1,c2):
    r_average = int(round(math.sqrt((c1["r"]^2 + c2["r"]^2)/2),0))
    g_average = int(round(math.sqrt((c1["g"]^2 + c2["g"]^2)/2),0))
    b_average = int(round(math.sqrt((c1["b"]^2 + c2["b"]^2)/2),0))
    return ({'r': r_average, 'g': g_average, 'b': b_average, 'a': 0})

class PanEncoder(json.JSONEncoder):
	"""
		json encoder class for specific list-like or dict-like python classes.
	"""
	def default(self, obj):
		if isinstance(obj, set):
			return list(obj)# turns sets into a list.
			# Let the base class default method raise the TypeError in case there is a problem
		return json.JSONEncoder.default(self, obj)
