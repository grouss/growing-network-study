import pickle
import numpy as np
import warnings
from scipy.special import zeta
import scipy.optimize as opt

from . import config

warnings.filterwarnings("ignore")

libname="mle"

# Copyright 2025 Université Paris Cité, France
# Author: [Guillaume Rousseau](https://www.linkedin.com/in/grouss/), Physics Department, Paris, France
# This is supplemental materials and replication package associated with the preprint available on 
# - arXiv (https://arxiv.org/abs/2501.10145)
# - SSRN  (http://ssrn.com/abstract=5191689
# Current version of python scripts and ressources are available on [github's author page](https://github.com/grouss/growing-network-study)
# This work is currently licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)



# Functions implementing method to estimate alpha
# based on Clauset et al. 2009



def GetXmin(x,AllValues=False):
    """
    This function returns the list of threshold values x_c corresponding to 
    the beginning of the distribution's tail, over which the power-law exponent 
    alpha[x, x_c] will be estimated.
        Parameters:
            x (array) : Observed values from the distribution.
            AllValues : bool, optional (default=False).
                        If True, return all values of x for which the empirical distribution y[x] != 0.
        Returns:
            x_min (array) : Threshold values x_c to be used for estimating the power-law exponent.
    """
    if AllValues:
        return x[x!=0][:-1] # excluding last one and x==0
    else:
        missing=0
        x_min=[]
        # Here we want to avoid the calculation of alpha for all values
        xlogspace=np.unique(np.logspace(0, np.log10(max(x)+0.1),1000,dtype=int))
        xlogspace=xlogspace[xlogspace>=1]
        for xm in xlogspace:
            if xm in x:
                x_min.append(xm)
            else:
                # as soon as we have reached the first missing value of x
                break
        # and we add all larger values of x that exist
        for xm in x[x>max(x_min)]:
            if xm!=x[-1]:
                # exclude last value (it leads to warning because it leads to 1/ln(1))
                x_min.append(xm)
        return np.array(x_min)

def f_exposant_hat(x,y,shift=0.5, AllValues=False,WKS="weightedKS",debug=False):
    """
    Partial implementation (step 1 and 2) of the method proposed 
    by Clauset et al. (2009) for estimating the exponent of a power-law 
    distribution from empirical data.
        Parameters :
            x : array of int or float
                Array of observed values for which the number of occurrences is known.
            y : array of int
                Array of counts such that y[i] corresponds to the number of occurrences of value x[i].
            Shift : float, optional (default=0.5)
                    Constant used to adapt the exponent estimation depending on the data type:
                        - shift = 0 for continuous data
                        - shift = 0.5 for discrete data
            AllValues : bool, optional (default=False)
                        Used to return all non zero values of x
            WKS : bool (optional, default=True)
                  Parameter to set if we use KS, weightedKS, or Kuiper

        Returns :
            x_min : array of int
                    List of threshold values x_c corresponding to the beginning of the distribution's tail.            
            e_x_min : array of float
                      Estimated power-law exponent alpha[x, x_c] for each candidate threshold value.          
            D_max : array of float
                    Kolmogorov–Smirnov distance between the empirical distribution and the best-fit power-law for each candidate x_c.
    """
    # List of threshold values x_c corresponding 
    # to the beginning of the distribution's tail.
    x_min=GetXmin(x,AllValues)
    
    # Estimated power-law exponent alpha[x, x_c] 
    # for each candidate threshold value. 
    e_x_min = np.zeros_like(x_min,dtype='float')
    for i,xm in enumerate(x_min):
        e_x_min[i]=exposant_hat(xm,x,y,shift)
    
    # Kolmogorov–Smirnov distance between the empirical distribution 
    # and the best-fit power-law for each candidate x_c.
    D_max=np.zeros_like(x_min,dtype='float')
    debugD=[]
    for i,xm in enumerate(x_min): 
        if shift==0:
            # perfect power-law distribution for a continuous random variable
            pc=pcdf_continuous(xm,e_x_min[i],x[x>=xm]) 
        elif shift==0.5:
            # Perfect power-law distribution for a discrete random variable
            pc=pcdf_discrete(xm,e_x_min[i],x[x>=xm]) # "perfect (power law)
        else:
            print("ERROR unknown shift value")
        s=y[x>=xm]
        sc=1-np.cumsum(s[::-1])[::-1]/np.sum(s) # observed complementary cumulative distribution function for  x>=xm
        if WKS=="weightedKS": 
            # weighted D_KS
            mask=np.logical_and(pc!=0,pc!=1)
            D=np.abs(sc-pc)[mask]/np.sqrt(pc*(1-pc))[mask] 
            D_max[i]=D.max()
            debugD.append({"x":x[x>=xm][mask],"D":D,"pc":pc[mask],"sc":sc[mask]})
        elif WKS=="KS":
            # KS
            D=np.abs(sc-pc)
            D_max[i]=D.max()
            debugD.append({"x":x[x>=xm],"D":D,"pc":pc,"sc":sc})
        elif WKS=="Kuiper":
            # Kuiper test
            D_max[i]=(sc-pc).max()+(pc-sc).max()
        else:
            print("ERROR unknown WKS")
            return 0
    if debug:
        return x_min,e_x_min,D_max,debugD
    else:
        return x_min,e_x_min,D_max

def exposant_hat(xmin,x,y,shift=0.5):
    """
    This function return the estimator of the scaling-factor 
    based on ne maximum of the likelyhood of the sample. 
    Two expressiosn are used corresponding to the case of continuous 
    random variable, and an approximation for discrete random variable,
    which can be used for large values of xmin according to clauset et al. 2009 
    In practice we use this approximation for xmin>=100, and keep the explicit form for xmin<100.
        Parameters:
            xmin : float or int 
                   threshold values x_c corresponding to the beginning of the distribution's tail
            x : array of int or float
                Array of observed values for which the number of occurrences is known.
            y : array of int
                Array of counts such that y[i] corresponds to the number of occurrences of value x[i].
            Shift : float, optional (default=0.5)
                    Constant used to adapt the exponent estimation depending on the data type:
                        - shift = 0 for continuous data
                        - shift = 0.5 for discrete data
        Returns:
            e_hat : float 
                    Value of the estimate
    """
    sy = y[x>=xmin]
    sx = x[x>=xmin]
    n=sy.sum()
    sumv=np.sum(sy*np.log(sx))
    #e_hat = 1 + sy.sum()/ ( sy*np.log( sx / (xmin-shift) ) ).sum()
    e_hat = 1 + n/ ( sumv-n*np.log(xmin-shift))
    if xmin<=100 and shift==0.5:
        alphastart=e_hat # we use the approximation as initial value
        result= opt.minimize(LLDiscretePL,x0=alphastart,args=(xmin,n,sumv),tol=alphastart/100,options={"disp":False}) 
        e_hat=result.x[0]
    return e_hat

def LLDiscretePL(alpha,xmin,n,sumv):
    """
    (- LogLikelihood) function to minimize 
        Parameters:
            alpha : float.
                    Value of the scaling factor for which the LL is calculated
            n : int
                sum_i(y_i)
            xmin : int or float.
                   cut off value 
            sumv : float
                sum_i(y_i*log(x_i))
        Returns:
            LLDiscretePL : float
                   n*np.log(zeta(alpha,xmin))+alpha*sumv
    """
    return n*np.log(zeta(alpha,xmin))+alpha*sumv

# continuous case
def pcdf_continuous(xmin,e,s):
    """
    Cumulative distribution fonction of a power-law distribution for a continuous random variable x>xmin and an exponent e 
        Parameters:
            xmin : float (xmin!=0)
                   Minimal value for which the pdf is not zero
            e : float
                Scaling factor (i.e. exponent of the power-law distribution)
            s : float,int, array of float, or aray of int
                Value for which the 
        Returns:
            pcdf : float or array of float 
                   Cumulative Distribution Function 
    """
    pcdf=1-(s/xmin)**(1-e)
    return pcdf

# dicrete case
def pcdf_discrete(xmin,e,s):
    """
    Cumulative distribution fonction of a power-law distribution for a continuous random variable x>=xmin and an exponent e 
         Parameters:
            xmin : float or int (xmin!=0)
                   Minimal value for which the pdf is not zero
            e : float
                Scaling factor (i.e. exponent of the power-law distribution)
            s : float,int, array of float, or aray of int
                Value for which the 
        Returns:
            pcdf : float or array of float 
                   Cumulative Distribution Function
    """
    pcdf=1-zeta(e,s)/zeta(e,xmin)
    return pcdf
    

def Get_yfit_yc_y_x(x,y,yc,Fall=False,Verbose=False,shift=0.5,WKS="weightedKS",XCmin=0):
    """
    This function call f_exposant_hat, determine the best value for the cut-off value 
    associated with the min of the max of the KS distance, and return a fitted function, 
    proportional to the pdf of the best power-law distribution according 
    to step 1 and 2 of Clauset et al. 2009. 
        Parameters:
            x : array of int or float
                Array of observed values for which the number of occurrences is known.
            y : array of int
                Array of counts such that y[i] corresponds to the number of occurrences of value x[i].
            yc : Array of int or float
                 complementary cumulative distribution function
            Fall : bool (optional, default=False)
                   Return only yfit, xfit (False) or  
                   yfit,xfit,imin,x_min_shift,e_x_min_shift,D_max_shift (True)
            Verbose : bool (optional, default=False)
                      Print extra information for debug/verbose mode
            Shift : float, optional (default=0.5)
                    Constant used to adapt the exponent estimation depending on the data type:
                        - shift = 0 for continuous data
                        - shift = 0.5 for discrete data (bias correction)
            AllValues : bool, optional (default=False)
                        Used to return all non zero values of x
            WKS : bool (optional, default="weightedKS")
                  Parameter to set if we use KS, weightedKS, or Kuiper
            XCmin : Integer or float (optional, default=0)
                    Minimal value of x for the estimate of D.max()
        Returns:
            tuple:
                yfit,xfit : if Fall==False  
                yfit,xfit,imin,x_min_shift,e_x_min_shift,D_max_shift : if Fall==True  
    """
    x_min_shift,e_x_min_shift,D_max_shift=f_exposant_hat(x,y,shift,WKS=WKS)
    mask=~np.isnan(D_max_shift) # nan in extreme case calling zeta function
    x_min_shift=x_min_shift[mask]
    e_x_min_shift=e_x_min_shift[mask]
    # excluding value cc<XCmin
    D_max_shift=D_max_shift[mask]
    imin=np.argmin(np.where(x_min_shift<XCmin,np.max(D_max_shift),D_max_shift)) 
    #print("XCmin",XCmin,x_min_shift[imin])
    minDmax=D_max_shift[imin]
    xfit=x[x!=0]
    cte=np.sum(y[x>=x_min_shift[imin]]) 
    if shift==0:
        yfit=(x_min_shift[imin]/xfit)**(e_x_min_shift[imin]-1)*cte
    else:
        cte=cte/zeta(e_x_min_shift[imin],x_min_shift[imin])
        yfit=zeta(e_x_min_shift[imin],xfit)*cte
        
    # cte is such that yfit[imin]=yc[imin] 
    # i.e. (xfit,yfit) cross the cdf of the dataset (x,yc) 
    # for the best xmin value
    # using continuous random variable
    if Fall:
        if Verbose: print("e hat",e_x_min_shift[imin])
        return yfit,xfit,imin,x_min_shift,e_x_min_shift,D_max_shift
    else:
        return yfit,xfit



def power_law(x, a, b):
    return np.log(a)+b*np.log(x)