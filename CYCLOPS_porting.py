#!/usr/bin/env python
# coding: utf-8

# In[13]:


# runCYCLOPS_Eigen.jl

# parse the input parameter

import argparse
import pandas as pd
from CYCLOPS_2a_seed import getseed, dispersion_ip
from CYCLOPS_2a_PreNPostprocessModule import GetEigenGenes

## -c, -i, -s, -o required -> need to edit for final version
parser = argparse.ArgumentParser(description = "CYCLOPS porting into python")
#parser.add_argument("-c", "--cycdir", dest = "cycdir", action = "store", required = True, help = "the complete path of CYCLOPS code")
parser.add_argument("-i", "--infile", dest = "infile", action = "store", default = "/disk1_ssd/bikyw/2.Circadian/1.Cyclops_tutorial/Oslops/example.ToRunCYCLOPS.csv", help = "the input expression dataset, csv file path")
parser.add_argument("-s", "--seedfile", dest = "seedfile", action = "store", default = "/disk1_ssd/bikyw/2.Circadian/1.Cyclops_tutorial/Oslops/clockList.csv",  help = "the seed gene list, csv file path")
parser.add_argument("-o", "--outdir", dest = "outdir", action = "store", default = "/disk1_ssd/bikyw/2.Circadian/3.Cyclops_python/result/", help = "the directory store output csv files")
parser.add_argument("--Frac_Var", dest = "Frac_Var", action = "store", type = float, default = 0.99, help = "set number of dimensions of SVD to maintain this fraction of variance")
parser.add_argument("--DFrac_Var", dest = "DFrac_Var", action = "store", type = float, default = 0.01, help = "set number of dimensions of SVD to so that incremetal fraction of variance of var is at least this much")
parser.add_argument("--Seed_MinCV", dest = "Seed_MinCV", action = "store", type = float, default = 0.14, help = "set the minimal CV for filtering the genes in the seed list, genes with CV below this value will be removed")
parser.add_argument("--Seed_MaxCV", dest = "Seed_MaxCV", action = "store", type = float, default = 0.85, help = "set the maximal CV for filtering the genes in the seed list, genes with CV above this value will be removed")
parser.add_argument("--Seed_Blunt", dest = "Seed_Blunt", action = "store", type = float, default = 0.975, help = "set the blunt number for removing outliers")
parser.add_argument("--Seed_MinMean", dest = "Seed_MinMean", action = "store", type = float, default = -100000000000000.00, help = "set the minimal mean expression cutoff for seed list, genes with mean expression level below this value will be removed")
parser.add_argument("--Out_Symbol", dest = "Out_Symbol", action = "store", help = "the symbol used in the output file")

args = parser.parse_args()

Seed_MinCV   = args.Seed_MinCV
Seed_MaxCV   = args.Seed_MaxCV
Seed_Blunt   = args.Seed_Blunt
Seed_MinMean = args.Seed_MinMean
Frac_Var     = args.Frac_Var
DFrac_Var    = args.DFrac_Var

# read the seed gene list
seedfile              = str(args.seedfile)
fullnonseed_data_BHTC = pd.read_csv(seedfile)
bhtc_seeds            = fullnonseed_data_BHTC.iloc[0:,0]

## if you want to slice df, should use .iloc function
# get the eigen gene expression profile
## in python version, colnames are not included in dataframe
## in julia version, colnames are included in matrix
fullnonseed_data_merge            = pd.read_csv(str(args.infile))
fullnonseed_data2                 = pd.concat([fullnonseed_data_merge.iloc[:,0],fullnonseed_data_merge.iloc[:,0],fullnonseed_data_merge], axis = 1) # hcat means merging matrix by column
#alldata_samples                   = fullnonseed_data2.iloc[0,3:]
alldata_samples                   = list(fullnonseed_data2.columns)[3:]
seed_symbols_bhtc, seed_data_bhtc = getseed(fullnonseed_data2, bhtc_seeds, Seed_MaxCV, Seed_MinCV, Seed_MinMean, Seed_Blunt)
seed_data_bhtc                    = dispersion_ip(seed_data_bhtc)

outs_bhtc, norm_seed_data_bhtc, eigen_sig_fraction, eigen_sum_fraction = GetEigenGenes(seed_data_bhtc, Frac_Var, DFrac_Var, 300)

# output the dispersed seed gene expression and the eigen expression matrix
seed_out = pd.concat([seed_symbols_bhtc, seed_data_bhtc], axis = 1)
seed_out = seed_out.rename(columns = {seed_out.columns[0] : "geneName"}) # already sample col names are assigned
seed_out.to_csv(str(args.outdir) + str(args.Out_Symbol) + "_SeedgeneExp.csv", index = False)
eigen_out = pd.DataFrame(["eigen_" + str(i+1) + "_" + str(round(eigen_sig_fraction[i], 4)) for i in range(int(outs_bhtc))])
eigen_out = pd.concat([eigen_out, pd.DataFrame(norm_seed_data_bhtc)], axis = 1)
eigen_out.columns = ["eigenName"] + alldata_samples
eigen_out.to_csv(str(args.outdir) + str(args.Out_Symbol) + "_EigengeneExp.csv", index = False)

# write the file parameters out
args_dict = vars(args)
outdir = str(args.outdir)
parafile = str(args.Out_Symbol) + "_para.txt"
f = open(outdir + parafile, "a")
for key, values in args_dict.items():
    f.write(str(key))
    f.write(" : ")
    f.write(str(values)) 
    f.write("\n")
f.write("\nThe important intermediate parameters in ordering with each eigen cluster is listed below:\n")
f.close()

### should change datatype into float64 due to difference between julia version and python version

