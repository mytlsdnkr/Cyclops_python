import argparse
import os
import sys
import polars as pl
#import pandas as pd

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
from CYCLOPS_2a_PreNPostprocessModule import *

fullnonseed_data_BHTC=pl.read_csv(seedfile)

colnames=fullnonseed_data_BHTC.columns

bhtc_seeds=fullnonseed_data_BHTC.get_column(colnames[0]).to_list()

fullnonseed_data_merge=pl.read_csv(infile)
colnames=fullnonseed_data_merge.columns

geneSymbol=fullnonseed_data_merge.get_column(colnames[0]).to_list()

seed_symbols_bhtc,seed_data_bhtc=getseed(fullnonseed_data_merge,colnames,bhtc_seeds,Seed_MaxCV,Seed_MinCV,Seed_MinMean,Seed_Blunt)
seed_data_bhtc=dispersion(seed_data_bhtc)


print(seed_data_bhtc)



outs_bhtc, norm_seed_data_bhtc, eigen_sig_fraction, eigen_sum_fraction=GetEigenGenes(seed_data_bhtc,Frac_Var,DFrac_Var,300)


print("outs_bhtc:",outs_bhtc)

print("norm_seed_data_bhtc:",norm_seed_data_bhtc)

exit()


seed_data_bhtc=pl.DataFrame(seed_data_bhtc)

seed_out=seed_data_bhtc.select([pl.Series(name="GeneSymbol",values=seed_symbols_bhtc),pl.all()])


Seed_gene_path=outdir+"SeedgeneExp.csv"
seed_out.write_csv(Seed_gene_path,sep=",")

eigen_list=[]

eigen_sig_fraction=eigen_sig_fraction[0:outs_bhtc]

for i in range(0,outs_bhtc):
    eigen_list.append("eigen_"+str(i))

norm_seed_data_bhtc=pl.DataFrame(norm_seed_data_bhtc)

for i in range(0,len(colnames)-1):
    norm_seed_data_bhtc=norm_seed_data_bhtc.rename({"column_"+str(i):colnames[i+1]})

norm_seed_data_bhtc=norm_seed_data_bhtc.select([pl.Series(name="Sig_fraction",values=eigen_sig_fraction),pl.all()])
norm_seed_data_bhtc=norm_seed_data_bhtc.select([pl.Series(name="Eigen_name",values=eigen_list),pl.all()])

Eigen_gene_filename="EigengeneExp.csv"
norm_seed_data_bhtc.write_csv(outdir+Eigen_gene_filename,sep=",")
args_dict=vars(args)
f=open(outdir+"Parameter.txt","w")

f.write("The important intermediate parameters in ordering with each eigen cluster is listed below:\n")
for key in args_dict:
    f.write(key+":"+str(args_dict[key])+"\n")

f.close()


#### Run Oscope ####

cmd="Rscript runOscope.R --infile "+Eigen_gene_filename+" --outdir "+outdir

os.system(cmd)
