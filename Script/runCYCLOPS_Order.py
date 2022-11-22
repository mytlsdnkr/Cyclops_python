import argparse
import os
import sys
import polars as pl
import datetime
import random
#import pandas as pd

parser = argparse.ArgumentParser(description='Cyclops ordering argparser')

#### Add parser options ####
parser.add_argument("--cycdir","-c",
help="the complete path of CYCLOPS code. eg. /TOOLS/code/Cyclops",
required=True)

parser.add_argument("--infile","-i",
help="the input expression dataset, csv file (including path)",
required=True)

parser.add_argument("--normseed","-n",
help = "the complete file name (including path) of seed expression data outputed from 'runCYCLOPS_Eigen.py'",
required=True)

parser.add_argument("--eigenfile","-e",
help = "the complete file name (including path) of eigen expression data outputed from 'runCYCLOPS_Eigen.jl' ",
required=True)

parser.add_argument("--oscopecluster", "-p",
        help = "the complete file name (including path) of eigen expression data outputed from 'runCYCLOPS_Eigen.jl' ",
required=True)

parser.add_argument("--outdir", "-o",
        help = "the directory store CYCLOPS' ouput csv files",
required=True)

parser.add_argument(
    "--Seed_Random",
	help = "the random seed number used for CYCLOPS",
    type = int,
	default  = 12345,
)

parser.add_argument(
    "--N_best",
	help = "Number of random initial conditions to try for each optimization",
    type = int,
	default  = 40,
)

parser.add_argument(
    "--N_background",
	help = "Number of background runs for global background distribution",
    type = int,
	default  = 40,
)

parser.add_argument(
    "--Seed_MinMean",
	help = "Set the mimimal mean expression cutoff for seed list, genes with mean expression level below this value will be removed",
    type = float,
	default = -100000000000000.00,
)

parser.add_argument(
    "--Out_Symbol",
	help = "The symbol used in the output file to differentiate multiple runs for the same input file",
)

args=parser.parse_args()


# Variable initilize #

print("Below is your options\n\n")
cyc_dir=args.cycdir
print("Cyclops directory:",cyc_dir)

infile=args.infile
print("Gene expression file:",infile)

seedfile=args.normseed
print("Seed gene list:",seedfile)

normseed=args.normseed
print("Seed gene expression file: ",normseed)

eigenfile=args.eigenfile
print("Eigen gene expression file: ",eigenfile)

oscope_cluster=args.oscopecluster
print("Oscope cluster file: ",oscope_cluster)


seed_random=args.Seed_Random
print("Seed Random: ",seed_random)

N_best=args.N_best
print("Number of random initial conditions: ",N_best)

total_background_num=args.N_background
print("Number of background runs: ",total_background_num)

seed_minmean=args.Seed_MinMean
print("Minimal mean expression cutoff for seed list: ",seed_minmean)

outdir=args.outdir
print("Output directory",outdir)

parafile = "Parameter.txt"

f = open(outdir+parafile, "a")
f.write("CYCLOPS StartTime:"+str(datetime.datetime.now())+ "\nThe input file parameters for running CYCLOPS is as below: \n")
f.close()

sys.path.insert(0,cyc_dir)




def filter_ebox(eboxphases_bhtc):
    if gene_means!=0:
        return x>gene_means

    if gene_cvs_min!=0:
        return x>gene_cvs_min
    
    if gene_cvs_max!=0:
        return x<gene_cvs_max
    


from CYCLOPS_2a_Seed import *
from CYCLOPS_2a_PreNPostprocessModule import *
from CYCLOPS_2a_AutoEncoderModule import *
from CYCLOPS_2a_MultiCoreModule_Smooth import *
from CYCLOPS_2a_MCA import *
###########################################################
#read in the expression data
fullnonseed_data_merge=pl.read_csv(infile)


#fullnonseed_data2=hcat(fullnonseed_data_merge[:,1], fullnonseed_data_merge[:,1], fullnonseed_data_merge) ## the default get seed function assumes first column=probe, 2nd symbol, 3rd entrez (or just text)
#alldata_samples=fullnonseed_data2[1,4:end]

alldata_samples=fullnonseed_data_merge.columns[1:]


#read in the data output from 'runCYCLOPS_Eigen.jl'
seed_data_bhtc=pl.read_csv(seedfile)

seed_colnames=seed_data_bhtc.columns
seed_symbols_bhtc=seed_data_bhtc.get_column(seed_colnames[0]).to_list()


seed_data_bhtc=seed_data_bhtc[1:,1:]
seed_data_bhtc=seed_data_bhtc.with_columns(pl.all().cast(pl.Float64,strict=False))
norm_seed_data_bhtc = pl.read_csv(eigenfile)
norm_seed_data_bhtc = norm_seed_data_bhtc[0:,2:].with_columns(pl.all().cast(pl.Float64,strict=False))
random.seed(seed_random)

oscope_eigen_group=pl.read_csv(oscope_cluster)

nrow=oscope_eigen_group.shape[0]
ncol=oscope_eigen_group.shape[1]


def condition(x,pval=0,ptr=0,gene_mean=0):
    if pval!=0:
        return x<pval


    if ptr!=0:
        return x>ptr

    if gene_mean!=0:
        return x>gene_mean




for i in range(0,nrow):
    tind=oscope_eigen_group[i,ncol-1].split("|")
    ind=list(map(int,tind))

    norm_seed_data_oscope=norm_seed_data_bhtc[ind,:]
    outs_bhtc_oscope=len(ind)

    N_best=1

    estimated_phaselist_bhtc,bestnet_bhtc,global_var_metrics_bhtc = CYCLOPS_Order(outs_bhtc_oscope,norm_seed_data_oscope,N_best)
    #estimated_phaselist_bhtc=[value + (2*math.pi) for value in estimated_phaselist]
    #estimated_phaselist_bhtc=[value % (2*math.pi) for value in estimated_phaselist]
    global_smooth_metrics_bhtc=smoothness_measures(seed_data_bhtc,norm_seed_data_bhtc,estimated_phaselist_bhtc)

    global_metrics_bhtc=global_var_metrics_bhtc

    pvals=multicore_backgroundstatistics_global_eigen(seed_data_bhtc,outs_bhtc_oscope,N_best,total_background_num,global_metrics_bhtc)

    cosinor_skin_bhtc=Compile_MultiCore_Cosinor_Statistics(fullnonseed_data_merge,estimated_phaselist_bhtc,4,48)

    eboxgenes=["DBP","HLF","TEF","PER1","NR1D2","CIART","C1orf51","PER3"]

    index_list=[]

    for i in range(0,len(eboxgenes)):
        index_list.append(cosinor_skin_bhtc.index[cosinor_skin_bhtc['GeneSymbol']==eboxgenes[i]].to_list())



    index_list=list(filter(None,index_list))
    real_index=[]

    for i in range(0,len(index_list)):
        real_index.append(index_list[i][0])

    eboxphases_bhtc=cosinor_skin_bhtc.iloc[real_index,:]


    pvalue_list=eboxphases_bhtc["pval"].to_list()
    ptr_list=eboxphases_bhtc["ptr"].to_list()
    mean_list=eboxphases_bhtc["fitMean"].to_list()


    criteria1=[idx for idx, element in enumerate(pvalue_list) if condition(element,pval=0.1)]
    criteria2=[idx for idx, element in enumerate(ptr_list) if condition(element,ptr=1.25)]
    criteria3=[idx for idx, element in enumerate(mean_list) if condition(element,gene_mean=seed_minmean)]

    print("Criteria1:",criteria1)
    print("Criteria2:",criteria2)
    print("Criteria3:",criteria3)
    
    temp=set(criteria1).intersection(criteria2)
    all_criteria=set(temp).intersection(criteria3)

    print("all_criteria:",all_criteria)

    


    
    







    exit()


    



###########################################################



'''
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

outs_bhtc, norm_seed_data_bhtc, eigen_sig_fraction, eigen_sum_fraction=GetEigenGenes(seed_data_bhtc,Frac_Var,DFrac_Var,300)

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
'''
