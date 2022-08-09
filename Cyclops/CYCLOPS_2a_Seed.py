import polars as pl
import pandas as pd
import pyarrow
import math

def clean_data(expression_data,bluntpercent=0.99):

    expression_data=expression_data.to_pandas()
    ngenes,nsamples=expression_data.shape
    nfloor=1+math.floor((1-bluntpercent)*nsamples)
    nceiling=math.ceil(bluntpercent*nsamples)

    for row in range(0,ngenes):

        temp=expression_data.iloc[row].to_list()
        temp.sort()
        vfloor=int(temp[nfloor])
        vceil=int(temp[nceiling])

        for sample in range(0,nsamples):
            expression_data.iloc[row,sample]=max(vfloor,expression_data.iloc[row,sample])
            expression_data.iloc[row,sample]=min(vceil ,expression_data.iloc[row,sample])

    return expression_data


def getseed(data,colnames,symbol_list,maxcv=.75,mincv=.07,minmean=500,blunt=.99):
    data_symbols=data.get_column(colnames[0]).to_list()
    expression_data=data.drop(colnames[0])
    expression_data=expression_data.with_columns(pl.all().cast(pl.Float64,strict=False))
    expression_data=clean_data(expression_data,blunt)
    ngenes,nsamples=expression_data.shape

    gene_means=expression_data.mean(axis=1).to_list()
    gene_sds=expression_data.std(axis=1).to_list()

    gene_cvs = [i / j for i, j in zip(gene_sds, gene_means)]

    both = set(data_symbols).intersection(symbol_list)

    criteria1=[data_symbols.index(x) for x in both]


    print(data_symbols)
    print(symbol_list)
    print(criteria1)

    





    




# function getseed(data::Array{Any,2},symbol_list,maxcv=.75,mincv=.07,minmean=500,blunt=.99)
# 	data_symbols=data[2:end ,2]
# 	data_data=data[2:end , 4:end];
# 	data_data=float64(data_data)
# 	data_data=clean_data!(data_data,blunt)
# 	ngenes,namples=size(data_data)

# 	gene_means=mean(data_data,2)
# 	gene_sds=std(data_data,2)
# 	gene_cvs= gene_sds ./ gene_means

# 	criteria1=findin(data_symbols,symbol_list)
# 	criteria2=findin((gene_means .> minmean),true)
# 	criteria3=findin((gene_cvs .> mincv),true)
# 	criteria4=findin((gene_cvs .< maxcv),true)

# 	allcriteria=intersect(criteria1,criteria2,criteria3,criteria4)
# 	seed_data=data_data[allcriteria,:]
# 	seed_symbols=data_symbols[allcriteria,:]
# 	seed_symbols, seed_data
# end