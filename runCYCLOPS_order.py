#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse
import datetime
import pandas as pd
import numpy as np
import random


# parsing
parser = argparse.ArgumentParser(description = "CYCLOPS porting into Python")
parser = argparse.ArgumentParser(description = "CYCLOPS porting into python")
#parser.add_argument("-c", "--cycdir", dest = "cycdir", action = "store", required = True, help = "the complete path of CYCLOPS code")
parser.add_argument("-i", "--infile", dest = "infile", action = "store", default = "/disk1_ssd/bikyw/2.Circadian/1.Cyclops_tutorial/Oslops/example.ToRunCYCLOPS.csv", help = "the input expression dataset, csv file path")
parser.add_argument("-n", "--normseed", dest = "normseed", action = "store", default = "/disk1_ssd/bikyw/2.Circadian/3.Cyclops_python/result/None_SeedgeneExp.csv", help = "the complete file name of seed expression data outputed from runCYCLOPS_Eigen.jl")
parser.add_argument("-e", "--eigenfile", dest = "eigenfile", action = "store", default = "/disk1_ssd/bikyw/2.Circadian/3.Cyclops_python/result/None_EigengeneExp.csv", help = "the complete file name of eigen expression data outputed from runCYCLOPS_Eigen.jl")
parser.add_argument("-p", "--oscopecluster", dest = "oscopecluster", action = "store", default = "/disk1_ssd/bikyw/2.Circadian/3.Cyclops_python/result/None_EigengeneExp.csv", help = "the complete file name of eigen expression data outputed from runCYCLOPS_Eigen.jl")
parser.add_argument("-o", "--outdir", dest = "outdir", action = "store", default = "/disk1_ssd/bikyw/2.Circadian/3.Cyclops_python/result/", help = "the directory store output csv files")
parser.add_argument("--Seed_Random", dest = "Seed_Random", action = "store", default = int(12345), help = "the random seed number used for CYCLOPS")
parser.add_argument("--N_best", dest = "N_best", action = "store", default = int(40), help = "Number of random initial conditions to try for each optimization")
parser.add_argument("--N_background", dest = "N_background", action = "store", default = int(40), help = "Number of background runs for global background distribution")
parser.add_argument("--Seed_MinMean", dest = "Seed_MinMean", action = "store", default = float(-100000000000000.00), help = "Set the mimimal mean expression cutoff for seed list, genes with mean expression level below this value will be removed")
parser.add_argument("--Out_Symbol", dest = "Out_Symbol", action = "store", help = "the symbol used in the output file")

args = parser.parse_args()

N_best               = args.N_best
total_background_num = args.N_background
Seed_MinMean         = args.Seed_MinMean
Seed_Random          = args.Seed_Random

# write the file parameters out
args_dict = vars(args)
outdir = str(args.outdir)
parafile = str(args.Out_Symbol) + "_para_from_runCYCLOPS_order.txt"
f = open(outdir + parafile, "a")
f.write("CYCLOPS StartTime : ")
f.write(str(datetime.datetime.now()))
f.write("\nThe input file parameters for running CYCLOPS is as below: \n")
for key, values in args_dict.items():
    f.write(str(key))
    f.write(" : ")
    f.write(str(values)) 
    f.write("\n")
f.write("\nThe important intermediate parameters in ordering with each eigen cluster is listed below:\n")
f.close()

# read in the expression data
fullnonseed_data_merge            = pd.read_csv(str(args.infile))
fullnonseed_data2                 = pd.concat([fullnonseed_data_merge.iloc[:,0],fullnonseed_data_merge.iloc[:,0],fullnonseed_data_merge], axis = 1) # hcat means merging matrix by column
alldata_samples                   = list(fullnonseed_data2.columns)[3:]

# read in the data output from 'runCYCLPOS_Eigen.jl'
seed_data_bhtc = pd.read_csv(str(args.normseed))
seed_symbols_bhtc = seed_data_bhtc.iloc[:, 0]
seed_data_bhtc = seed_data_bhtc.iloc[:,1:].astype(np.float64)
norm_seed_data_bhtc = pd.read_csv(str(args.eigenfile))
norm_seed_data_bhtc = norm_seed_data_bhtc.iloc[:,1:].astype(np.float64)

random.seed(Seed_Random)
oscope_eigen_group = pd.read_csv(str(args.oscopecluster))

print(oscope_eigen_group)
exit()

