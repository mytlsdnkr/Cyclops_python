import numpy as np
from scipy.stats import f
from scipy.stats import *
import pandas as pd


def Corrected_Cosinor_Statistics_Faster(expression,o_PREDICTED_PHASELIST,n_shifts):

    print("expression:",expression)
    length=len(o_PREDICTED_PHASELIST)

    shift_list=np.linspace(2*np.pi/n_shifts,2*np.pi,n_shifts)

    expression=expression.to_list()

    min_error=np.var(expression)*10E20


    best_shift="error"


    # Find Best linear model #

    for shift in shift_list:
        l_PREDICTED_PHASELIST= (o_PREDICTED_PHASELIST+shift)%(2*np.pi)

        XLINEAR=np.array(l_PREDICTED_PHASELIST)

        

        A=np.vstack([XLINEAR,np.ones(len(XLINEAR))]).T

        m,c=np.linalg.lstsq(A,expression,rcond=None)[0]

        predict_linear=XLINEAR*m+c
        
        sse_linear=np.sum(np.square(expression-predict_linear))
        

        if sse_linear < min_error:
            min_error=sse_linear
            best_shift=shift

    l_PREDICTED_PHASELIST= (o_PREDICTED_PHASELIST+best_shift)%(2*np.pi)

    SIN_l_PREDICTED_PHASELIST=np.sin(l_PREDICTED_PHASELIST)
    COS_l_PREDICTED_PHASELIST=np.cos(l_PREDICTED_PHASELIST)


    XLINEAR=np.array(l_PREDICTED_PHASELIST)
    A=np.vstack([XLINEAR,np.ones(len(XLINEAR))]).T
    m_linear,b_linear=np.linalg.lstsq(A,expression,rcond=None)[0]
    predict_linear=XLINEAR*m_linear+b_linear

    print("SIN_l_PREDICTED_PHASELIST:",SIN_l_PREDICTED_PHASELIST)
    print("SIN_l_PREDICTED_PHASELIST shape:",np.shape(SIN_l_PREDICTED_PHASELIST))
    print("COS_l_PREDICTED_PHASELIST:",COS_l_PREDICTED_PHASELIST)
    print("COS_l_PREDICTED_shape:",np.shape(COS_l_PREDICTED_PHASELIST))
    print("XLINEAR:",XLINEAR)
    print("XLINEAR.sahpe:",np.shape(XLINEAR))

    print("predicted_linear:",predict_linear)
    print("predicted_linear.shape:",np.shape(predict_linear))
    

    # Full model #
    XFULL=np.vstack((l_PREDICTED_PHASELIST,SIN_l_PREDICTED_PHASELIST,COS_l_PREDICTED_PHASELIST))

    print("XFULL:",XFULL)
    print("XFULL.shape:",np.shape(XFULL))
    A=XFULL.T

    print("A:",A)
    print("A.shape:",np.shape(A))
    mod_full=np.linalg.lstsq(A,expression,rcond=None)
    m_full=mod_full[0]
    b_full=mod_full[1]
    predict_full=A@m_full+b_full

    print("predict_full:",predict_full)

    # Find Residual of two model #
    sse_linear=np.sum(np.square(expression-predict_linear))
    sse_full=np.sum(np.square(expression-predict_full))

    print("sse_linear:",sse_linear)
    print("sse_full:",sse_full)



    # F-Test #
    f_metric=((sse_linear-sse_full)/2)/((sse_full)/(length-5))

    print("length:",length)

    print("f_metric:",f_metric)

    pval=1-f.cdf(f_metric,2,length-5)

    print("pval:",pval)



    SIN_o_PREDICTED_PHASELIST=np.sin(o_PREDICTED_PHASELIST)
    COS_o_PREDICTED_PHASELIST=np.cos(o_PREDICTED_PHASELIST)
    
    XCOSINOR=np.vstack((l_PREDICTED_PHASELIST,SIN_l_PREDICTED_PHASELIST,COS_l_PREDICTED_PHASELIST))
    
    A=XCOSINOR.T
    mod_cosinor=np.linalg.lstsq(A,expression,rcond=None)
    m_cosinor=mod_full[0]
    b_cosinor=mod_full[1]
    predict_cosinor=A@m_cosinor+b_cosinor

    sse_cosinor=np.sum(np.square(expression-predict_cosinor))
    expression=np.array(expression)
    mean_expression=np.mean(expression)
    sse_base=np.sum(np.square(expression-mean_expression))

    r2=1-(sse_cosinor/sse_base)

    print("r2:",r2)

    sinterm,costerm=m_cosinor[0],b_cosinor[0]
    
    print("sinterm,costerm",m_cosinor,b_cosinor)
    acrophase=np.arctan2(sinterm,costerm)
    print("acrophase:",acrophase)
    amp=np.sqrt(np.square(sinterm)+np.square(costerm))

    print("amp:",amp)
    fitavg=b_cosinor[0]
    print("fitavg:",fitavg)
    avg=np.mean(expression)
    
    return pval,acrophase,amp,fitavg,avg,r2





def MultiCore_Cosinor_Statistics(data,PREDICTED_PHASELIST,n_shifts=20):


    ngenes,nsamples=data.shape

    PrbPval=np.zeros(ngenes)
    PrbPhase=np.zeros(ngenes)
    PrbAmp=np.zeros(ngenes)
    PrbFitMean=np.zeros(ngenes)
    PrbMean=np.zeros(ngenes)
    PrbR2=np.zeros(ngenes)
    
    for row in range(0,ngenes):
        PrbPval[row],PrbPhase[row],PrbAmp[row],PrbFitMean[row],PrbMean[row],PrbR2[row]=Corrected_Cosinor_Statistics_Faster(data.iloc[row,:],PREDICTED_PHASELIST,n_shifts)


    return PrbPval,PrbPhase,PrbAmp,PrbFitMean,PrbMean,PrbR2





    

'''
def MultiCore_Cosinor_Statistics(data::Array{Float64,2},PREDICTED_PHASELIST::Array{Float64,1},n_shifts=20) 
    ngenes,nsamples=size(data)

    for row=1:ngenes
        PrbPval[row],PrbPhase[row],PrbAmp[row],PrbFitMean[row],PrbMean[row],PrbR2[row]=Corrected_Cosinor_Statistics_Faster(data[row,:],PREDICTED_PHASELIST,n_shifts)
    end
    PrbPval,PrbPhase,PrbAmp,PrbFitMean,PrbMean,PrbR2
end


'''


def Bonferroni_Adjust(pvalues):
    n=pvalues.shape[0]

    adjusted=np.zeros(n)



    for cnt in range(0,n):
        adjusted[cnt]=min([n*pvalues[cnt],1.])



    return adjusted



def Compile_MultiCore_Cosinor_Statistics(annotated_data,PREDICTED_PHASELIST,NumericStartCol=4,n_shifts=20):

    annotated_data=annotated_data.to_pandas()

    annotated=annotated_data.iloc[:,0].to_list()

    data=annotated_data.iloc[:,1:]



    ngenes=annotated_data.shape

    estimated_phaselist=PREDICTED_PHASELIST


    PrbPval,PrbPhase,PrbAmp,PrbFitMean,PrbMean,PrbRsq=MultiCore_Cosinor_Statistics(data,estimated_phaselist)

    PrbPtr=(PrbAmp+PrbFitMean)/(PrbFitMean-PrbAmp)

    PrbBon=Bonferroni_Adjust(PrbPval)


    cosinor_output=pd.DataFrame({
        "GeneSymbol":annotated,
        "pval":PrbPval,
        "bon_pval":PrbBon,
        "phase":PrbPhase,
        "amp":PrbAmp,
        "fitMean":PrbFitMean,
        "mean":PrbMean,
        "rsq":PrbRsq,
        "ptr":PrbPtr
        })




    print("cosinor_output:",cosinor_output)


    return cosinor_output
    
