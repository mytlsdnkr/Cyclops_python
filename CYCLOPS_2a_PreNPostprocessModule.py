#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np

def GetEigenGenes(numeric_data, fraction_var = 0.75, dfrac_var = 0.05, maxeig = 50):
    numeric_data = np.array(numeric_data).astype(np.float64)
    u, w, v = np.linalg.svd(numeric_data, full_matrices= False) # in julia, default of svd function is thin svd #(17,17)(17,)(533,17)
    v = v.T
    expvar = np.cumsum(w**2)/np.sum(w**2)
    eigen_sum_fraction = expvar
    eigen_sig_fraction = (w**2)/np.sum(w**2)
    ReductionDim1 = 1 + len(expvar[expvar <= fraction_var])
    vardif = np.diff(expvar) # infinite difference function
    ReductionDim2 = 1 + len(vardif[vardif >= dfrac_var])
    ReductionDim = min(ReductionDim1, ReductionDim2)
    Transform = v[:,0:ReductionDim].T
    
    return ReductionDim, 10*Transform, eigen_sig_fraction, eigen_sum_fraction

    


# a = np.array([[0.5600838385894501,0.671481430808863,0.7510562341172665,0.8012195266889552,0.8435022872072205,0.8769554355957212,0.9016200464679558,0.920145075916983,0.9356067445584285,0.9489062810629535,0.9609903219760676,0.9707664280728471,0.978807707384975,0.9858120476886983,0.9914702983368832,0.9961424175635639,1.0000000000000004]])
# b = np.diff(a)
# 
# ReductionDim1 = 1 + len(a[a <= 0.99])
# ReductionDim2 = 1 + len(b[b >= 0.01])
# 
# print(ReductionDim1)
# print(ReductionDim2)
# ReductionDim = min(ReductionDim1, ReductionDim2)
# 10*a
