#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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

class Layer_Connections:
    w = np.array([])
    b = np.array([])

def Initialize_Layer_Connections(layer_dim, in_dim):
    w = np.random.normal(size = (layer_dim, in_dim))
    b = np.random.normal(size = layer_dim)/100
    
    Layer_Connections.w = w
    Layer_Connections.b = b
##############################################################################
class Input_Layer:
    a = np.array([])

def Create_InputLayer(in_dim :int):
    Input_Layer.a = np.array([0]*in_dim)
##############################################################################
class BottleNeck_Layer:
    z = np.array([])
    a = np.array([])
    
    a_func = np.array([])
    jstar  = np.int8()

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
    
    BottleNeck_Layer.z      = z
    BottleNeck_Layer.a      = a
    BottleNeck_Layer.a_func = a_func
    BottleNeck_Layer.jstar  = jstar
##############################################################################          
class Output_Layer:
    z = np.array([])
    a = np.array([])
    
    a_func  = np.array([])
    a_deriv = np.array([])
    
def Create_OutputLayer(out_dim :int, activ_fun, activ_deriv):
    z       = np.array([0]*out_dim)
    a       = np.array([0]*out_dim)
    a_func  = np.tile(np.array([activ_fun]), out_dim)
    a_deriv = np.tile(np.array([activ_deriv]), out_dim)
    
    Output_Layer.z        = z
    Output_Layer.a        = a
    Output_Layer.a_func   = a_func
    Output_Layer.a_deriv  = a_deriv
##############################################################################     
class NeuralNetwork:
    dim     = int()
    nbottle = int()
    ncirc   = int()
    ilayer  = Input_Layer()
    blayer  = BottleNeck_Layer()
    olayer  = Output_Layer()
    c2      = Layer_Connections()
    c3      = Layer_Connections()
        
def Create_Network(size :int, bottle_size :int, n_circ :int):
    Create_InputLayer(size)
    Create_BottleNeckLayer(bottle_size, n_circ)
    Create_OutputLayer(size, linr, linr_deriv)
    
    NeuralNetwork.dim     = size
    NeuralNetwork.nbottle = bottle_size
    NeuralNetwork.ncirc   = n_circ
    NeuralNetwork.ilayer  = Input_Layer
    NeuralNetwork.blayer  = BottleNeck_Layer
    NeuralNetwork.olayer  = Output_Layer
    
    Initialize_Layer_Connections(bottle_size, size)
    NeuralNetwork.c2      = Layer_Connections
    
    Initialize_Layer_Connections(size, bottle_size)
    NeuralNetwork.c3      = Layer_Connections    
############################################################################## 
def Feed_Forward_ip(data, NN :NeuralNetwork):
    NN.ilayer.a = data
    NN.blayer.z = (NN.c2.w)*(NN.ilayer.a)+(NN.c2.b)
    
    for j in range(0, NN.ncirc):
        jstar = Find_Partner(j+1)
        NN.blayer.a[j] = NN.blayer.a_func[j](NN.blayer.z[j], NN.blayer.z[jstar])
        
    for j in range(NN.ncirc, NN.nbottle):
        NN.blayer.a[j] = NN.blayer.a_func[j](NN.blayer.z[j])
        
    NN.olayer.z = (NN.c3.w)*(NN.blayer.a)+(NN.c3.b)
    
    for j in range(0, NN.dim):
        NN.olayer.a[j] = NN.olayer.a_func[j](NN.olayer.z[j])

