import polars as pl
import pandas as pd
import pyarrow
import math
from numpy import array
from scipy.linalg import svd


def GetEigenGenes(numeric_data,fraction_var=.75, dfrac_var=0.05,maxeig=50):

    numeric_data=numeric_data.to_numpy()

    print(numeric_data)
    u,w,v=svd(numeric_data)
    print("u:",u)
    print("v:",v)
    print("w:",w)


    # expvar=cumsum(w.^2)/sum(w.^2)
    # println("expvar:",expvar)
    # eigen_sum_fraction = expvar
    # eigen_sig_fraction = (w.^2)/sum(w.^2)  #### 전체 Eigengene 중에서 특정 eigengene의 비율을 본다. 즉 특정 eigengene이 얼마나 중요한지에 대해 확인할 수 있다.                                                        ##modification here
    
    # println("eigen_sig_fraction:",eigen_sig_fraction)
    
    # ReductionDim1=1+length(expvar[expvar .<= fraction_var]);
    # println(ReductionDim1)

    # println("expvar:",expvar)
    # vardif=diff(expvar)
    # println("vardif:",vardif)
    # ReductionDim2=1+length(vardif[vardif .>= dfrac_var]);
    # ReductionDim=minimum([ReductionDim1,ReductionDim2])
    # println("ReductionDim:",ReductionDim)
    # Transform=v[:,1:ReductionDim]'
    # println("Transform:",Transform)

    # ReductionDim, 10*Transform, eigen_sig_fraction, eigen_sum_fraction    