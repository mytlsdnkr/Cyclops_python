import pandas as pd
import math

def clean_data(expression_data,bluntpercent=0.99):

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


def getseed(data,symbol_list,maxcv=.75,mincv=.07,minmean=500,blunt=.99):
    print(symbol_list)
    expression_data=data.iloc[0:,1:]
    expression_data.astype('float64').dtypes
    expression_data=clean_data(expression_data,blunt)


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