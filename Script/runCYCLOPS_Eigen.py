import argparse
import os
import sys
import pandas as pd

parser = argparse.ArgumentParser(description='Cyclops argparser')


#### Add parser options ####
parser.add_argument("--cycdir","-c",
help="the complete path of CYCLOPS code. eg. /TOOLS/code/JuliaCYCLOPS2a",
required=True)

parser.add_argument("--infile","-i",
help="the input expression dataset, csv file (including path)",
required=True)

parser.add_argument("--seedfile","-s",
help="the seed gene list, csv file (including path)",
required=True)

parser.add_argument("--outdir","-o",
help="the directory store ouput csv files",
required=True)

parser.add_argument("--Frac_Var",
help="Set Number of Dimensions of SVD to maintain this fraction of variance",
type=float,
default=0.99)

parser.add_argument("--DFrac_Var",
help="Set Number of Dimensions of SVD to so that incremetal fraction of variance of var is at least this much",
type=float,
default=0.01)


parser.add_argument("--Seed_MinCV",
help="Set the minimal CV for filtering the genes in the seed list, genes with CV below this value will be removed",
type=float,
default=0.14)

parser.add_argument("--Seed_MaxCV",
help="Set the maximal CV for filtering the genes in the seed list, genes with CV above this value will be removed",
type=float,
default=0.85)

parser.add_argument("--Seed_Blunt",
help="Set the blunt number for removing outliers",
type=float,
default=0.975)

parser.add_argument("--Seed_MinMean",
help="Set the mimimal mean expression cutoff for seed list, genes with mean expression level below this value will be removed",
type=float,
default=-100000000000000.00)

parser.add_argument("--Out_Symbol",
help="The symbol used in the output file")


args=parser.parse_args()

# Variable initilize #

print("Below is your options\n\n")
cyc_dir=args.cycdir
print("Cyclops directory:",cyc_dir)
infile=args.infile
print("Gene expression file:",infile)
seedfile=args.seedfile
print("Seed gene list:",seedfile)
outdir=args.outdir
print("Output directory",outdir)


Seed_MinCV=args.Seed_MinCV
print("Seed Min CV:",Seed_MinCV)
Seed_MaxCV=args.Seed_MaxCV
print("Seed Max CV:",Seed_MaxCV)
Seed_Blunt=args.Seed_Blunt
print("Seed Blunt:",Seed_Blunt)
Seed_MinMean=args.Seed_MinMean
print("Seed MinMean:",Seed_MinMean)
Frac_Var=args.Frac_Var
print("Frac Var:",Frac_Var)
DFrac_Var=args.DFrac_Var
print("DFrac Var:",DFrac_Var)


sys.path.insert(0,cyc_dir)

from CYCLOPS_2a_Seed import *


fullnonseed_data_BHTC=pd.read_csv(seedfile)
bhtc_seeds=fullnonseed_data_BHTC['geneSymbol'].tolist()


fullnonseed_data_merge=pd.read_csv(infile)
# geneSymbol=fullnonseed_data_merge['Gene.Symbol'].tolist()

# fullnonseed_geneSymbol=pd.DataFrame({'GeneSymbol0':geneSymbol,
# 'GeneSymbol1':geneSymbol
# })
#fullnonseed_data2=pd.concat([fullnonseed_geneSymbol,fullnonseed_data_merge],axis=1)

sampleName=fullnonseed_data_merge.columns[1:].tolist()


getseed(fullnonseed_data_merge,bhtc_seeds,Seed_MaxCV,Seed_MinCV,Seed_MinMean,Seed_Blunt)

#seed_symbols_bhtc, seed_data_bhtc=getseed(fullnonseed_data_merge,bhtc_seeds,Seed_MaxCV,Seed_MinCV,Seed_MinMean,Seed_Blunt)




