import polars as pl
import pandas as pd
import pyarrow
import math
import numpy as np
import random
from numpy import array,matlib
from scipy.linalg import svd
from sklearn.decomposition import PCA


class input_layer:
    def __init__(self,l):
        self.l=l

class BottleNeck_Layer:
    def __init__(self,z,a,a_func,jstar):
        self.z=z
        self.a=a
        self.a_func=a_func
        self.jstar=jstar
class Output_Layer:
    def __init__(self,z,a,a_func,a_deriv):
        self.z=z
        self.a=a
        self.a_func=a_func
        self.a_deriv=a_deriv

class Layer_Connections:
    def __init__(self,w,b):
        self.w=w
        self.b=b

class Create_Network:
    def __init__(self,w,b):
        self.w=w
        self.b=b

class NeuralNetwork:
    def __init__(self,dim,nbottle,ncirc,ilayer,blayer,olayer,c2,c3):
        self.dim=dim
        self.nbottle=nbottle
        self.ncirc=ncirc
        self.ilayer=ilayer
        self.blayer=blayer
        self.olayer=olayer
        self.c2=c2
        self.c3=c3



def linr(z):
    return z

def linr_deriv(z,dummy=1.0):
    return 1.

def circ(z,zstar):


    return z/(np.sqrt((z ** 2)+(zstar ** 2)))


def Find_Partner(x):

    grp=int((x)/2)
    elm=int((x)%2)

    pelm=(1-elm)
    partner=1+(grp)*2+pelm
    return partner-1

def Create_BottleNeckLayer(layer_size,n_circ):

    print("Bottleneck layer size:",layer_size)

    print("Number of circle nodes in bottleneck layer:",n_circ)

    z=[0.0 for i in range(0,layer_size)]
    a=[0.0 for i in range(0,layer_size)]

    jstar=[0 for i in range(0,layer_size)]


    a_func_circ=[circ for i in range(n_circ)]

    print("Number of activation function in circ nodes:",len(a_func_circ))

    a_func_linr=[linr for i in range(layer_size-n_circ)]

    print("Number of activation function in linear nodes:",len(a_func_linr))
    a_func=a_func_circ+a_func_linr

    print("Activation function in bottleneck layer:",a_func)

    for i in range(0,n_circ):
        jstar[i]=Find_Partner(i)

    for i in range(n_circ+1,layer_size):
        jstar[i]=i

    print("partner index in bottleneck layer:",jstar)

    return BottleNeck_Layer(z,a,a_func,jstar)

def Create_OutputLayer(out_dim,activ_func,activ_deriv):
    z=[0.0 for i in range(0,out_dim)]
    a=[0.0 for i in range(0,out_dim)]

    a_func=[activ_func for i in range(out_dim)]

    print("Activation function in output layer:",a_func)
    a_deriv=[activ_deriv for i in range(out_dim)]

    print("A_driv in output layer:",a_deriv)

    return Output_Layer(z, a, a_func, a_deriv)

def Initialize_Layer_Connections(layer_dim,in_dim):

    print("Layer_dimension:",layer_dim)
    print("Input dimension:",in_dim)

    w=np.random.randn(layer_dim, in_dim)
    b=np.random.randn(layer_dim)/100


    return Layer_Connections(w,b)


def Create_InputLayer(in_dim):
    l=[0.0 for i in range(0,in_dim)]

    return input_layer(l)

def Create_Network(size,bottle_size,n_circ):
    print("Create Network Size:",size)
    print("Bottlenectk Size:",bottle_size)
    print("Number of circular nodes:",n_circ)
    ilayer=Create_InputLayer(size)
    print(ilayer.__dict__)
    blayer=Create_BottleNeckLayer(bottle_size,n_circ)
    print(blayer.__dict__)
    olayer=Create_OutputLayer(size,linr,linr_deriv)
    print(olayer.__dict__)
    print("C2!!!")
    c2=Initialize_Layer_Connections(bottle_size,size)
    print(c2.__dict__)
    print("C3!!")
    c3=Initialize_Layer_Connections(size,bottle_size)
    print(c3.__dict__)

    return NeuralNetwork(size,bottle_size,n_circ,ilayer,blayer,olayer,c2,c3)


def Feed_Forward(data,NN):
    data=data.to_list()
    NN.ilayer.l=data # set input nodes=data

    NN.ilayer.l=np.array(NN.ilayer.l)
    NN.blayer.z=((NN.c2.w)@(NN.ilayer.l))+(NN.c2.b)

    print("Input:",NN.ilayer.l)

    print("NN.blayer.z:",NN.blayer.z)




    for i in range(0,NN.ncirc):
        jstar=Find_Partner(i)
        NN.blayer.a[i]=NN.blayer.a_func[i](NN.blayer.z[i],NN.blayer.z[jstar])

        print("Partner nodes::",jstar)

        print("Activation:",NN.blayer.a[i])



    

    for i in range(NN.ncirc+1,NN.nbottle):
        NN.blayer.a[i]=NN.blayer.a_func[i](NN.blayer.z[i])

    NN.olayer.z=(NN.c3.w)@(NN.blayer.a)+(NN.c3.b)


    for i in range(0,NN.dim):
        NN.olayer.a[i]=NN.olayer.a_func[i](NN.olayer.z[i])




def Find_Gradients(data,NN):

    print("Gradients data:",data)
    Feed_Forward(data,NN)
    c2grads=Layer_Connections(0* NN.c2.w, 0*NN.c2.b)
    c3grads=Layer_Connections(0* NN.c3.w, 0*NN.c3.b)


    print("C2grads:",c2grads.__dict__)

    print("C3grads:",c2grads.__dict__)



    data=data.to_list()
    data=[value * -1 for value in data]


    #del_o=(NN.olayer.a-data)  #simple difference between layer and goal

    del_o=[ai+bi for ai,bi in zip(data,NN.olayer.a)]


    err=(0.5*(sum(del_o)**2))

    
    for i in range(0,NN.dim):
        print(i)
        del_o[i]=NN.olayer.a_deriv[i](NN.olayer.z[i])*del_o[i]
        c3grads.b[i]=del_o[i]

        print("c3grads_b:",c3grads.b[i])

        for j in range(0,NN.nbottle):
            c3grads.w[i,j]=(del_o[i])*(NN.blayer.a[j])

            print("c3grads.w:",c3grads.w[i,j])

    r=[0.0 for i in range(0,NN.ncirc)]

    for i in range(0,NN.ncirc):
        jstar=Find_Partner(i)
        r[i]=math.sqrt((NN.blayer.z[i])**2+(NN.blayer.z[jstar])**2)

        print("r[i]:",r[i])


    del_b=[0.0 for i in range(0,NN.nbottle)]

    for i in range(0,NN.ncirc):
        jstar=Find_Partner(i)
        dsum=0
        for j in range(0,NN.dim):
            d=del_o[j]*1/(r[i]**3)*(NN.c3.w[j,i]*(NN.blayer.z[jstar])**2-NN.c3.w[j,jstar]*(NN.blayer.z[i])*(NN.blayer.z[jstar]))

            print("d:",d)
            dsum=dsum+d

            print("dsum:",dsum)

        del_b[i]=dsum
        c2grads.b[i]=del_b[i]
        for j in range(0,NN.dim):
            c2grads.w[i,j]=(del_b[i])*(NN.ilayer.l[j])
            print("c2grads.w:",c2grads.w[i,j])




    print(NN.ncirc)
    print(NN.nbottle)
    for i in range(NN.ncirc,NN.nbottle):
        print("shouldn't evaluate if only circular bottle")
        dsum=0
        for j in range(0,NN.dim):
            d=NN.c3.w[j,i]*del_o[j]
            dsum=dsum+d

        del_b[i]=dsum
        c2grads.b[i]=del_b[i]
        
        for j in range(0,NN.dim):
            c2grads.w[i,j]=(del_b[i])*(NN.ilayer.l[j])



    return c2grads,c3grads,err




def Find_Total_Gradients(data_matrix,NN):
    
    c2wt=np.zeros(np.shape(NN.c2.w))
    c3wt=np.zeros(np.shape(NN.c3.w))
    c2bt=np.zeros(np.shape(NN.c2.b))
    c3bt=np.zeros(np.shape(NN.c2.b))
    terr=0
    ntimepoints=data_matrix.shape[1]
    
    print("timepoints:",ntimepoints)


    for n in range(0,ntimepoints):
        (c2gradients,c3gradients, err)=Find_Gradients(data_matrix.iloc[:,n],NN)
        c2wt=c2wt+(c2gradients.w)/ntimepoints
        c3wt=c3wt+(c3gradients.w)/ntimepoints
        c2bt=c2bt+(c2gradients.b)/ntimepoints
        c3bt=c3bt+(c3gradients.b)/ntimepoints
        terr=terr+err/ntimepoints


    c2tgrads=Layer_Connections(c2wt, c2bt)
    c3tgrads=Layer_Connections(c3wt, c3bt)

    return c2tgrads,c3tgrads,terr



def Train_Momentum_Stochastic(data_matrix,NN,batchsize=10,rate=.3,momentum=.5,tol=0.0005,epochs=10):


    oc2wt=np.zeros(np.shape(NN.c2.w))
    oc3wt=np.zeros(np.shape(NN.c3.w))
    oc2bt=np.zeros(np.shape(NN.c2.b))
    oc3bt=np.zeros(np.shape(NN.c3.b))



    o2changes=Layer_Connections(oc2wt, oc2bt)
    o3changes=Layer_Connections(oc3wt, oc3bt)

    new2changes=Layer_Connections(oc2wt,oc2bt)
    new3changes=Layer_Connections(oc3wt,oc3bt)
    (c2tgrads,c3tgrads,e0)=Find_Total_Gradients(data_matrix,NN)

    if NN.nbottle != NN.ncirc:
        rate=rate/NN.dim
    
    max=1
    n=0
    e1=1.0


    rowsize,trainsize=data_matrix.shape

    while (max/e1)>tol and (n<epochs):
        n=n+1
        new2changes.b=-rate*(c2tgrads.b)+momentum*o2changes.b
        new2changes.w=-rate*(c2tgrads.w)+momentum*o2changes.b
        new3changes.b=-rate*(c3tgrads.b)+momentum*o3changes.b
        new3changes.w=-rate*(c3tgrads.w)+momentum*o3changes.b

        NN.c2.b=NN.c2.b+new2changes.b
        NN.c2.w=NN.c2.w+new2changes.w
        NN.c3.b=NN.c3.b+new3changes.b
        NN.c3.w=NN.c3.w+new3changes.w

        e0=e1

        (c2tgrads,c3tgrads,e1)=Find_Total_Gradients(data_matrix.sample(n=batchsize,axis='columns'), NN)

        o2changes.b=new2changes.b
        o2changes.w=new2changes.w
        o3changes.b=new3changes.b
        o3changes.w=new3changes.w


        max=np.max([np.max(abs(c2tgrads.b)),np.max(abs(c2tgrads.w)),np.max(abs(c3tgrads.b)),np.max(abs(c3tgrads.w))])

       # (c2tgrads,c3tgrads,e1)=Find_Total_Gradients(data_matrix[0:random.sample(range(trainsize),batchsize)], NN)


def Train_bold(data_matrix,NN,rate=0.3,tol=0.0005,epochs=500):
    
    c2tgrads,c3tgrads,e0=Find_Total_Gradients(data_matrix, NN)

    max=1
    n=0
    e1=1.0


    if NN.nbottle != NN.ncirc:
        rate=rate/NN.dim

    while (max/e1>tol) and (n<epochs) and (rate>0.00001):
        n=n+1
        ONN=NN
        oc2tgrads,oc3tgrads=c2tgrads,c3tgrads
        NN.c2.b=NN.c2.b-rate*(c2tgrads.b)
        NN.c3.b=NN.c3.b-rate*(c3tgrads.b)
        NN.c2.w=NN.c2.w-rate*(c2tgrads.w)
        NN.c3.w=NN.c3.w-rate*(c3tgrads.w)

        c2tgrads,c3tgrads,e1=Find_Total_Gradients(data_matrix, NN)
        
        if e1>e0:
            NN=ONN
            rate=rate*0.75
            c2tgrads,c3tgrads,e1=oc2tgrads,oc3tgrads,e0

        elif e1<e0:
            rate=rate*1.05
            oc2tgrads,oc3tgrads,e0=c2tgrads,c3tgrads,e1
        max=np.max([np.max(abs(c2tgrads.b)),np.max(abs(c2tgrads.w)),np.max(abs(c3tgrads.b)),np.max(abs(c3tgrads.w))])

        return e1


def Default_Metrics(data,outsize):

    print("data:",data)
    print("shape:",data.shape)

    data=data.T
    
    totalsse=0

    pca=PCA()
    print("pca:",pca)

    PC=pca.fit_transform(data)

    print("PC:",PC)


    pca_result=pd.DataFrame(data=PC,columns=['PC1','PC2'])


    usesize=pca_result.shape[1]

    totalsse=totalsse+pca_result['PC1'].var()
    totalsse=totalsse+pca_result['PC2'].var()

    linearsse=totalsse-pca_result['PC1'].var()

    print("totalsse:",totalsse)

    print("Linear sse:",linearsse)


    return totalsse,linearsse





def ExtractPhase(data_matrix,NN):

    print("data_matrix:",data_matrix)
    

    points=data_matrix.shape[1]

    print("points:",points)
    phases=[0.0 for i in range(0,points)]

    print("phases:",phases)
    for n in range(0,points):
        Feed_Forward(data_matrix.iloc[:,n], NN)
        phases[n]=math.atan2(NN.blayer.a[0],NN.blayer.a[1])
        print("phases[n]:",phases[n])
    

    return phases
    
    

def CYCLOPS_Order(l_outs,l_norm_seed_data,N_trials):

    l_norm_seed_data=l_norm_seed_data.to_pandas()
    scalefactor=1;
    besterror=10E20
    bestnet=0;
    NET=0;

    print("l_outs:",l_outs)
    print("l_norm_seed_data:",l_norm_seed_data)
    print("n_Trials:",N_trials)


    for trial in range(0,N_trials):
        NET=Create_Network(l_outs,2,2);


        print("NET:",NET.__dict__)

        Train_Momentum_Stochastic(l_norm_seed_data,NET)  
        global_circ_sse=Train_bold(l_norm_seed_data,NET)

        if global_circ_sse<besterror:
            besterror=global_circ_sse
            bestnet=NET


        print("global_circ_sse:",global_circ_sse)



    total_sse,onedlinear_sse=Default_Metrics(l_norm_seed_data,l_outs)

    besterror=besterror*2
    global_metrics=[besterror,besterror/(total_sse-besterror),besterror/(onedlinear_sse-besterror)]

    print("global_metrics:",global_metrics)

    estimated_phaselist=ExtractPhase(l_norm_seed_data,bestnet);
    estimated_phaselist=[value + (2*math.pi) for value in estimated_phaselist]
    estimated_phaselist=[value % (2*math.pi) for value in estimated_phaselist]

    return estimated_phaselist,bestnet,global_metrics

