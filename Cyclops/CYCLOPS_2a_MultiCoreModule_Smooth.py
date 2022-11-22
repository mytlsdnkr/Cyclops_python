import polars as pl
import pandas as pd
import pyarrow
import math
import numpy as np
import random
from multiprocessing import Pool
from multiprocessing import Process
from numpy import array,matlib
from scipy.linalg import svd
from CYCLOPS_2a_AutoEncoderModule import *
from CYCLOPS_2a_MultiCoreModule_Smooth import *
from CYCLOPS_2a_PreNPostprocessModule import *



def circ_diff(data):
    diff=data.diff(axis=1)
    diff.iloc[:,0]=data.iloc[:,0]-data.iloc[:,-1]

    circd=diff

    return circd




def smoothness_measures(l_seeddata,l_eigendata,estimated_phaselist):

    l_eigendata=l_eigendata.to_pandas()
    l_seeddata=l_seeddata.to_pandas()
    estimated_phaselist1s=[value % (2*math.pi) for value in estimated_phaselist]
    
    estimated_phaselist1s=np.array(estimated_phaselist1s)


    use_order_c=np.argsort(estimated_phaselist1s,kind='mergesort')

    l_eigen_temp=np.array(l_eigendata.loc[0].to_list())

    use_order_l=np.argsort(l_eigen_temp,kind='mergesort')

    new_phaselist=estimated_phaselist1s[use_order_c]

    print("new phase list:",new_phaselist)

    print("use_order_c:",use_order_c)

    print("seed_Data:",l_seeddata)

    circ_seeddata=l_seeddata.iloc[:,use_order_c]
    circ_eigendata=l_eigendata.iloc[:,use_order_c]

    lin_seeddata=l_seeddata.iloc[:,use_order_l]
    lin_eigendata=l_eigendata.iloc[:,use_order_l]

    nsamp=l_seeddata.shape[1]


######################################################

    gdiffs=circ_diff(circ_seeddata)
    gdiffs=gdiffs.to_numpy()
    gdiffs2 =gdiffs*gdiffs
    gdiffsm=np.sqrt(np.sum(gdiffs2,axis=0))


    num=np.sum(gdiffsm)/nsamp



    gdiffs=lin_seeddata.diff(axis=1)
    gdiffs.iloc[:,0]=lin_seeddata.iloc[:,0]-lin_seeddata.iloc[:,-1] # Remove NaN column
    gdiffs=gdiffs.to_numpy()
    gdiffs2 =gdiffs*gdiffs
    gdiffsm=np.sqrt(np.sum(gdiffs2,axis=0))
    denom=sum(gdiffsm)/(nsamp-1)

    measure_raw=num/denom

######################################################

    gdiffs=circ_diff(circ_seeddata)
    gdiffs=gdiffs.to_numpy()

    gdiffs2 =gdiffs*gdiffs

    gdiffsm=np.sum(gdiffs2,axis=0)

    num=np.sum(gdiffsm)/nsamp


    gdiffs=lin_seeddata.diff(axis=1)
    gdiffs.iloc[:,0]=lin_seeddata.iloc[:,0]-lin_seeddata.iloc[:,-1] # Remove NaN column
    gdiffs=gdiffs.to_numpy()

    gdiffs2 =gdiffs*gdiffs
    gdiffsm=np.sum(gdiffs2,axis=0)
    denom=np.sum(gdiffsm)/(nsamp-1)

    measure_raw_sq=num/denom

    
######################################################

    gdiffs=circ_diff(circ_eigendata)
    gdiffs=gdiffs.to_numpy()

    gdiffs2 =gdiffs*gdiffs
    gdiffsm=np.sqrt(np.sum(gdiffs2,axis=0))
    num=np.sum(gdiffsm)/nsamp


    gdiffs=lin_eigendata.diff(axis=1)
    gdiffs.iloc[:,0]=lin_eigendata.iloc[:,0]-lin_eigendata.iloc[:,-1] # Remove NaN column
    gdiffs=gdiffs.to_numpy()

    gdiffs2 =gdiffs*gdiffs
    gdiffsm=np.sqrt(np.sum(gdiffs2,axis=0))
    denom=np.sum(gdiffsm)/(nsamp-1)

    measure_eigen=num/denom

#####################################################

    gdiffs=circ_diff(circ_eigendata)
    gdiffs=gdiffs.to_numpy()

    gdiffs2 =gdiffs*gdiffs

    gdiffsm=np.sum(gdiffs2,axis=0)

    num=np.sum(gdiffsm)/nsamp


    gdiffs=circ_eigendata.diff(axis=1)
    gdiffs.iloc[:,0]=circ_eigendata.iloc[:,0]-circ_eigendata.iloc[:,-1] # Remove NaN column
    gdiffs=gdiffs.to_numpy()

    gdiffs2 =gdiffs*gdiffs

    gdiffsm=np.sum(gdiffs2,axis=0)
    denom=np.sum(gdiffsm)/(nsamp-1)

    measure_eigen_sq=num/denom
    temp= [measure_eigen,measure_eigen_sq,measure_raw,measure_raw_sq]


    return temp

def backgroundmetrics_global_eigen(seed_ldata,ESize,N_best,N_runs):
    #background_global_metrics=[[0]*3 for _ in range(N_best)]
    background_global_metrics=np.zeros((N_runs,3))

    netsize=ESize

    for count in range(0,N_runs):
        best_error=10E40
        shuffled_seed_data=Row_Shuffle(seed_ldata)
        EigenShuffledData=GetNEigenGenes(shuffled_seed_data, ESize)
        EigenShuffledData=pl.DataFrame(EigenShuffledData).transpose()
        estimated_phaselist,bestnet,variance_global_matrics=CYCLOPS_Order(ESize,EigenShuffledData,N_best)
        background_global_metrics[count,:]=variance_global_matrics


    return background_global_metrics


def multicore_backgroundmetrics_global_eigen(seed_ldata,ESize,N_best,N_trials):
    '''   
    num_cores=5

    pool=Pool(num_cores)
    a1=pool.starmap(backgroundmetrics_global_eigen,zip(seed_ldata,ESize,N_best,(N_trials/5)))


    a1=backgroundmetrics_global_eigen(seed_ldata,ESize,N_best,int(N_trials/5))
    '''

    a1=backgroundmetrics_global_eigen(seed_ldata,ESize,N_best,int(N_trials/5))
    a2=backgroundmetrics_global_eigen(seed_ldata,ESize,N_best,int(N_trials/5))
    a3=backgroundmetrics_global_eigen(seed_ldata,ESize,N_best,int(N_trials/5))



    return [a1,a2,a3]


def ret_num(gm,testmetrics,N_best):

    count=0
    for row in range(0,N_best):
        for col in range(0,3):
            if gm[row,col]<=testmetrics:
                count=count+1

    return count
    
def multicore_backgroundstatistics_global_eigen(seed_ldata,ESize,N_best,N_trials,testmetrics):
    backgroundmetrics=multicore_backgroundmetrics_global_eigen(seed_ldata, ESize, N_best, N_trials)

    print("backgroundmetrics:",backgroundmetrics)

    gm1=backgroundmetrics[0]
    gm2=backgroundmetrics[1]
    gm3=backgroundmetrics[2]

    nless1=ret_num(gm1, testmetrics[0], N_best)
    nless2=ret_num(gm2, testmetrics[1], N_best)
    nless3=ret_num(gm3, testmetrics[2], N_best)


    p1=nless1/N_trials
    p2=nless2/N_trials
    p3=nless3/N_trials


    return [p1,p2,p3]


    





    









