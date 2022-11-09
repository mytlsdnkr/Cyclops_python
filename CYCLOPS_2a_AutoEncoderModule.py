#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np

def linr(z :np.float64, dummy = 1.0):
    return z

def circ(z :np.float64, zstar :np.float64):
    result = z/(math.sqrt(z^2+zstar^2))
    return result

def linr_deriv(z :np.float64, dummy = 1.0):
    1.0

def Find_Partner(x :int):
    grp = (x-1)//2
    elm = (x-1)%2
    pelm = (1-elm)
    partner = 1+(grp*2)+pelm
    return partner

##############################################################################
def Input_Layer(a):
    a = np.float64(a)
    return a

def Create_InputLayer(in_dim :int):
    l = Input_Layer(np.array([0]*in_dim))
    return l
##############################################################################
def BottleNeck_Layer(z, a, jstar):
    z = np.float64(z)
    a = np.float64(a)
    
    jstar = np.int8(jstar)
    
    return z, a, jstar

# in julia, repmat function contain original data shape
# in python, np.tile doesn't contain original data shape

def Create_BottleNeckLayer(layer_size :int, n_circ :int):
    z = np.array([0]*layer_size)
    a = np.array([0]*layer_size)
    
    a_func = np.append(np.tile(np.array(circ),n_circ), np.tile(np.array(linr),layer_size - n_circ))
    
    jstar = np.array([0]*layer_size)
    for i in range(0, n_circ):
        jstar[i] = Find_Partner(i+1)
    for i in range(n_circ,layer_size):
        jstar[i] = i+1
    z, a, jstar = BottleNeck_Layer(z, a, jstar)
    
    return z, a, a_func, jstar
##############################################################################          
def Output_Layer(z, a):
    z = np.float64(z)
    a = np.float64(a)
    
    return z, a
    
def Create_OutputLayer(out_dim :int, activ_fun, activ_deriv):
    z       = np.array([0]*out_dim)
    a       = np.array([0]*out_dim)
    a_func  = np.tile(np.array([activ_fun]), out_dim)
    a_deriv = np.tile(np.array([activ_deriv]), out_dim)
    
    z, a = Output_Layer(z,a)
    
    return z,a,a_func,a_deriv
##############################################################################     
def Create_Network(size :int, bottle_size :int, n_circ :int):
    ilayer = Create_InputLayer(size)
    blayer = Create_BottleNeckLayer(bottle_size, n_circ)
    olayer = Create_OutputLayer(size, linr, linr_deriv)
    
    return ilayer,blayer,olayer
    
tmpNN = Create_Network(4,2,2)

print(tmpNN[0])
print(tmpNN[1])
print(tmpNN[2])

