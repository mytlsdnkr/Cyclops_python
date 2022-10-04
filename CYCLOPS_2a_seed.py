#!/usr/bin/env python
# coding: utf-8

# In[154]:


import pandas as pd
import math
import numpy as np


# ! means in-place function that induced data is changed at the result
def clean_data_ip(data, bluntpercent = 0.99):
    ngenes = len(data.index)
    nsamples = len(data.columns)
    nfloor = 1+math.floor((1-bluntpercent)*nsamples)
    nceiling = math.ceil(bluntpercent*nsamples)
    for row in range(0, ngenes):
        sorted_data = sorted(list(data.iloc[row,:]))
        vfloor = sorted_data[nfloor]
        vceil = sorted_data[nceiling]
        for sample in range(0, nsamples):
            data.iloc[row, sample] = max(vfloor, data.iloc[row, sample])
            data.iloc[row, sample] = min(vceil, data.iloc[row, sample])
    return data

def getseed(data, symbol_list, maxcv = 0.75, mincv = 0.07, minmean = 500, blunt = 0.99):
    data_symbols = data.iloc[0:, 1]
    data_data = data.iloc[0:, 3:]
    data_data = clean_data_ip(data_data, blunt)
    ngenes = len(data_data.index)
    nsamples = len(data_data.columns)
    
    # axis = 1 means by row
    gene_means = data_data.mean(axis=1)
    gene_sds = data_data.std(axis=1)
    gene_cvs = gene_sds / gene_means
    
    # store index
    data_symbols = data_symbols.values.tolist()
    symbol_list = symbol_list.values.tolist()
    criteria1 = sorted([data_symbols.index(common) for common in set(data_symbols).intersection(symbol_list)])
    criteria2 = [true for true in range(0, len(gene_means)) if gene_means.iloc[true] > minmean]
    criteria3 = [true for true in range(0, len(gene_cvs)) if gene_cvs.iloc[true] > mincv]
    criteria4 = [true for true in range(0, len(gene_cvs)) if gene_cvs.iloc[true] < maxcv]

    allcriteria = set(criteria1).intersection(criteria2)
    allcriteria = set(allcriteria).intersection(criteria3)
    allcriteria = set(allcriteria).intersection(criteria4)
    
    seed_data = pd.DataFrame([data_data.iloc[i,:] for i in allcriteria])
    seed_symbols = pd.DataFrame([data_symbols[i] for i in allcriteria])
    
    return seed_symbols, seed_data

def dispersion_ip (data):
    ngenes = len(data.index)
    nsamples = len(data.columns)
    for gene in range(0, ngenes):
        genemean = data.iloc[gene, :].mean()
        for sample in range(0, nsamples):
            data.iloc[gene, sample] = (data.iloc[gene, sample] - genemean) / genemean
    return data

