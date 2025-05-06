import time
import numpy as np
import warnings
import datetime
from scipy.special import zeta


from . import config

warnings.filterwarnings("ignore")
libname="model"

# Copyright 2025 Université Paris Cité, France
# Author: [Guillaume Rousseau](https://www.linkedin.com/in/grouss/), Physics Department, Paris, France
# This is supplemental materials and replication package associated with the preprint available on 
# - arXiv (https://arxiv.org/abs/2501.10145)
# - SSRN  (http://ssrn.com/abstract=5191689
# Current version of python scripts and ressources are available on [github's author page](https://github.com/grouss/growing-network-study)
# This work is currently licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)

def BarabasiAlbertGraph(n = 10,m = 2,seed=None,Verbose=False,SelfLoop=False):
    """
    Implementation of the Barabasi-Albert Model, 
    using Compressed Sparse Row format to store nodes and edges. 
    Initial graph is a complete graph with m nodes.
        Parameters:
            n : int. 
                Number of nodes in the graph.
            m : int.
                Number of new edges per new nodes.
            seed : int. (Optional, Default=None)
                   Seed of the random generator
            Verbose : bool. (Optional, Default=False)
                      Display elapsed time information 
            SelfLoop : bool. (Optional, Default=False)
                       if True, nodes link to themself.
        Returns:
            nodes : array of int. 
            edges : array of int.
            Nnodes : int.
                Number of nodes
            Nedges : int.
                Number of edges
            SelfLoop : bool.
                       True if nodes have a selfloop (False otherwise) 
    """
    if n<m+1:
        print("Wrong n wrt n0",n,m+1)
        return -1
    n0=m+1
    nodes=np.arange(0,n+1,dtype='uint32')*m
    m0=m*n0
    edges=np.zeros(m*n,dtype='uint32')

    # initial condition
    k=0
    for i in range(n0):
        for j in range(n0):
            if i!=j:
                edges[k]=j
                k+=1
                
    if seed!=None:
        np.random.seed(seed)
        
    ti=time.time()
    # performance of this part could be improved
    # if needed.
    for i in range(n0,n):
        localset=set([])
        while len(localset)!=m:
            r=np.random.random()
            if r<0.5:
                j=int(i*r*2)
            else:
                j=edges[int(m*i*(r-0.5)*2)]
            localset.add(j)
        edges[m*i:m*i+m]=np.array(list(localset),dtype='uint32')
    tf=time.time()
    if Verbose:
        print("Elapse",tf-ti)
    nodes[-1]=m*n
    
    Nnodes=n
    Nedges=n*m
    
        
    if SelfLoop:
        nodes+=np.arange(0,n+1,dtype='uint32')
        edgestmp=np.zeros((n,m+1),dtype='uint32')
        edgestmp[:,-1]=np.arange(0,n,dtype='uint32')
        edgestmp[:,:-1]=edges.reshape((n,m))
        Nedges=n*(m+1)
        edges=edgestmp.flatten()
    return nodes,edges,Nnodes,Nedges,SelfLoop

def zipf_distribution(alpha,xmin,N,Debug=False,Fccdf=False):
    """
    Provide a discrete random sample of size N, scaling factor alpha, 
    and minimal value xmin following a disceret powerlaw (i.e. a zipf distribution)
    Use zipf function provided by numpy if xmin==0.
    """
    ti=time.time()
    rng = np.random.default_rng()
    if xmin==0:
        s = rng.zipf(alpha, size=N)
        y=np.bincount(s)
        x=np.arange(len(y))
        mask=y!=0
        x=x[mask]
        y=y[mask]
    else:
        s = rng.uniform(0,1,N)
        x=[]
        y=[]
        xi=xmin
        zetamin=zeta(alpha,xmin)
        while len(s)!=0:
            p=1-zeta(alpha,xi+1)/zetamin
            mask=s<=p
            ni=np.sum(mask)
            if ni!=0:
                x.append(xi)
                y.append(ni)
                s=s[np.logical_not(mask)]
            xi+=1
        x=np.array(x)
        y=np.array(y)
    tf=time.time()
    if Debug: print("Elpase ",tf-ti,"(s)")
    if Fccdf:
        return x,y,np.cumsum(y[::-1])[::-1]
    else:
        return x,y

def Sum_Two_Zipf_distribution(alpha_1,xmin_1,N_1,alpha_2,xmin_2,N_2,Fccdf=False):
    """
    Merge 2 Zipf distribution.
    """
    x_1,y_1=zipf_distribution(alpha_1,xmin_1,N_1)
    x_2,y_2=zipf_distribution(alpha_2,xmin_2,N_2)

    # merge the 2 distributions 

    x=np.arange(max(max(x_1),max(x_2))+1)
    y=np.zeros_like(x)

    y[x_1]+=y_1
    y[x_2]+=y_2
    
    mask=y!=0
    x=x[mask]
    y=y[mask]
    
    if Fccdf:
        return x,y,np.cumsum(y[::-1])[::-1]
    else:
        return x,y