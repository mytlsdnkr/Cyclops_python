import polars as pl
import pandas as pd
import pyarrow
import math
import numpy as np
from numpy import array
from scipy.linalg import svd

def condition(x,frac_var=0,dfrac_var=0):
    if frac_var!=0:
        return x<frac_var

    if dfrac_var!=0:
        return x>dfrac_var


def diff(x):
    list=[]
    for i in range(0,len(x)):
        list.append(1-x[i])

    return list
        

def GetEigenGenes(numeric_data,fraction_var=.75, dfrac_var=0.05,maxeig=50):

    numeric_data=numeric_data.to_numpy()

    print(numeric_data)
    u,w,v=svd(numeric_data)

    expvar=np.cumsum(np.array(w)**2)/np.sum(np.array(w)**2)


    eigen_sum_fraction=expvar
    eigen_sig_fraction=(np.array(w)**2)/np.sum(np.array(w)**2)

    print("Expvar:",expvar)
    
    ReductionDim1=1+len([idx for idx, element in enumerate(expvar) if condition(element,frac_var=fraction_var)])

    print(ReductionDim1)

    vardif=diff(expvar)

    ReductionDim2=1+len([idx for idx, element in enumerate(vardif) if condition(element,dfrac_var=dfrac_var)])

    ReductionDim=min(ReductionDim1,ReductionDim2)
    

    Transform=v[0:ReductionDim,:]


    return ReductionDim,10*Transform,eigen_sig_fraction,eigen_sum_fraction
