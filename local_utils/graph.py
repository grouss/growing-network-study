import pickle
import time
import json
import numpy as np



import warnings
from . import config
from .dataset import *

warnings.filterwarnings("ignore")
libname="graph"

# Copyright 2025 Université Paris Cité, France
# Author: [Guillaume Rousseau](https://www.linkedin.com/in/grouss/), Physics Department, Paris, France
# This is supplemental materials and replication package associated with the preprint available on 
# - arXiv (https://arxiv.org/abs/2501.10145)
# - SSRN  (http://ssrn.com/abstract=5191689
# Current version of python scripts and ressources are available on [github's author page](https://github.com/grouss/growing-network-study)
# This work is currently licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)

def GetNodesSonsId(index,nodes,edges):
    """
    Returns an array of target node indices for the edges originating from the given source node index,
    assuming a Compressed Sparse Row (CSR) format.
    See: https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)

    Parameters:
        index: Index (int) of the source node.
        nodes: Array of nodes in CSR format.
        edges: Array of edges in CSR format.

    Returns:
        array of indices of the target nodes connected to the source node.
    """
    return edges[nodes[index]:nodes[index+1]]

def GetNodesDescedants(index,nodes,edges):
    """
    Returns a set of node indices reachable by traversing the graph downstream 
    from the given source node index, assuming a Compressed Sparse Row (CSR) format.

    Parameters:
        index: Index (int) of the source node.
        nodes: Array of nodes in CSR format.
        edges: Array of edges in CSR format.

    Returns:
        knownindexset: Set of indices of the nodes reachable by downstream traversal of the graph.
                       NB: The starting node's index is not included in the returned set.
    """
    # basic implementation without using recursive call 
    knownindexset=set([])
    alreadyKnown=set([])
    currentset=set([index])    
    while len(currentset)!=0: 
        sonlist=[]
        for i in currentset: 
            sonlist+=edges[nodes[i]:nodes[i+1]].tolist()  
        sonset=set(sonlist)
        currentset=sonset-knownindexset
        knownindexset.update(currentset)            
    return knownindexset   


def GetNodesType(index,dRVindexMin=156799786,dOindexMax=139524532):
        """
        This function return node's type providing the node index, for the specific graph studied here
            Parameters:
                index (int): node index
                dRVindexMin (int): RV node min index
                dOindexMax  (int): O node max index
            Returns:
                string: The node type.
        """
        if index>=dRVindexMin:
            return "RV"
        elif index <= dOindexMax:
            return "O"
        else:
            return "RL"

def GetNodesTypesArrayTSL(nodes,edges,d=GetD(),depth=None, Debug=False):
    """
    Returns the TSL node type and the corresponding encoding.

    Parameters:
        nodes: Array of nodes in CSR format.
        edges: Array of edges in CSR format.
        d: Unused parameter (kept for consistency with the parameter list in related functions).
        depth (int or None): Maximum depth used for the TSL encoding.
        Debug (bool): If True, performs additional integrity checks.

    Returns:
        nodestype: Array of integers indicating the TSL type of each node, 
                   based on the associated encoding.
        encoding: Array of TSL-encoded strings corresponding to each type 
                  of node.
    """
    # TSL
        # L=1 if selfloop exiss, =0 otherwise
        # T=np.min(depth,nb of incoming edges)
        # S=np.min(depth,nb of outgoing edges)
    if type(edges)==type(None):
        print("ERROR edges==None")
        return None,None
    sourceEdges=GetSourceEdge(nodes)
    mask=edges==sourceEdges
    L=np.bincount(edges[mask],minlength=len(nodes)-1)
    # ! selfloop add 1 to S and L
    S=np.bincount(sourceEdges[np.logical_not(mask)],minlength=len(nodes)-1)
    T=np.bincount(edges[np.logical_not(mask)],minlength=len(nodes)-1)
    if Debug:
        if np.sum(S)!=np.sum(T):
            print("ERROR np.sum(S)!=np.sum(T)")
            return None,None
        if np.sum(S)+np.sum(T)+2*np.sum(L)!=2*len(edges):
            print("ERROR np.sum(S)+np.sum(T)+np.sum(L)!=len(edges)")
            print(np.sum(S),np.sum(T),2*np.sum(L),2*len(edges))
            print(np.sum(S)+np.sum(T)+2*np.sum(L),2*len(edges))
            print(np.sum(S)+np.sum(T)+2*np.sum(L)-2*len(edges))
            return None,None
    S[S>depth]=depth
    T[T>depth]=depth
    length=len(str(depth))
    encoding=[(2*length+1)*'X']*(depth+1)**3
    for t in range(depth+1):
        for s in range(depth+1):
            for l in range(2):
                #for debug 
                i=t*(depth+1)**2+s*(depth+1)+l
                encoding[i]=str(t).zfill(length)+str(s).zfill(length)+str(l)
    nodestype=T*(depth+1)**2+S*(depth+1)+L
    return nodestype,encoding
    
def GetEdgesTypesArrayTSL(nodes,edges,d=GetD(),depth=None,Debug=True):
    """
    Returns the TSL edge type and the corresponding encoding.

    Parameters:
        nodes: Array of nodes in CSR format.
        edges: Array of edges in CSR format.
        d: Unused parameter (kept for consistency with the parameter list in related functions).
        depth (int or None): Maximum depth used for the TSL encoding.
        Debug (bool): If True, performs additional integrity checks.

    Returns:
        nodestype: Array of integers indicating the TSL type of each edge, 
                   based on the associated encoding.
        encoding: Array of TSL-encoded strings corresponding to each type 
                  of edge.
    """
    if type(edges)==type(None):
        print("ERROR edges==None")
        return None,None
    sourceEdges=GetSourceEdge(nodes)
    nodestype,nodesencoding=GetNodesTypesArrayTSL(nodes,edges,d=d,depth=depth,Debug=Debug)
    edgestype=nodestype[sourceEdges]*len(nodesencoding)+nodestype[edges]
    edgesencoding=[]
    for se in nodesencoding:
        for te in nodesencoding:
            edgesencoding.append(se+">"+te)
    return edgestype,edgesencoding
    
def GetNodesTypesArray(nodes,edges,d=GetD(),depth=None,Debug=True):
    """
    Returns the node type and the corresponding encoding.

    Parameters:
        nodes: Array of nodes in CSR format.
        edges: Array of edges in CSR format.
        d: Dictionnary containing min and max index of each type of nodes.
        depth (int or None): Maximum depth used for the TSL encoding.
        Debug (bool): If True, performs additional integrity checks.

    Returns:
        nodestype: Array of integers indicating the type of each node, 
                   based on the associated encoding.
        encoding: Array of encoded strings corresponding to each type 
                  of node.
    """
    if depth!=None:
        return GetNodesTypesArrayTSL(nodes,edges,d=d,depth=depth,Debug=Debug)
    # depth!=None is used to return to switch to the TSL type of the O-(RV/RL)-O derived graph
    
    # the following one, which the default one, works for the full graph with
    # O/RV/RL nodes, but also for the O-(RV/RL)-O derived graph which only contains 
    # O nodes
    
    encoding=["O","RL","RV"]

    # 0 => "O"
    # 1 => "RL"
    # 2 => "RV"
    
    nodestype=2*np.ones(len(nodes)-1,dtype='int8')
    
    try:
        nodestype[d["OindexMin"]:d["OindexMax"]+1]=0
    except:
        pass
    try:
        nodestype[d["RLindexMin"]:d["RLindexMax"]+1]=1
    except:
        pass
    
    if Debug:
        for key,value in d.items():
            if "index" in key:
                if value<len(nodestype) and encoding[nodestype[value]]!=key.split("index")[0]:
                    print("ERROR GetNodesTypesArray",key,value,encoding[nodestype[value]])
                    return -1,-1
    return nodestype,encoding

def GetEdgesTypesArray(nodes,edges,d=GetD(),depth=None,Debug=False):
    """
    Returns the edge type and the corresponding encoding.

    Parameters:
        nodes: Array of nodes in CSR format.
        edges: Array of edges in CSR format.
        d: Dictionnary containing min and max index of each type of nodes.
        depth (int or None): Maximum depth used for the TSL encoding
        Debug (bool): If True, performs additional integrity checks.

    Returns:
        nodestype: Array of integers indicating the type of each edge, 
                   based on the associated encoding.
        encoding: Array of encoded strings corresponding to each type 
                  of edge.
    """
    if depth!=None:
        return GetEdgesTypesArrayTSL(nodes,edges,d=d,depth=depth,Debug=Debug)
    # depth!=None is used to switch to the TSL type of the O-(RV/RL)-O derived graph
    # and call the right function

    # # edges 3* source node type encoding + target node type encoding
    # 0 => "O>O"   # does not exist
    # 1 => "O>RL"
    # 2 => "O>RV"
    # 3 => "RL>O"  # does not exist
    # 4 => "RL>RL"
    # 5 => "RL>RV"
    # 6 => "RV>O"  # does not exist
    # 7 => "RV>RL" # does not exist
    # 8 => "RV>RV"
    encoding=["O>O","O>RL","O>RV","RL>O","RL>RL","RL>RV","RV>O","RV>RL","RV>RV"]
    
    # set source type bits (could use more generic approach, see below)
    # we assume edge type "RV>RV" by default
    edgestype=8*np.ones_like(edges,dtype='int8') 
    try:
        edgestype[nodes[d["RLindexMin"]]:nodes[d["RLindexMax"]+1]]-=3
    except:
        pass
    try:
        edgestype[nodes[d["OindexMin"]]:nodes[d["OindexMax"]+1]]-=6
    except:
        pass
    edgestype[edges<=d["RLindexMax"]]-=1
    edgestype[edges<=d["OindexMax"]]-=1
    if Debug:
        if np.min(edgestype)<0 or np.max(edgestype)>8:
            print("ERROR GetEdgesTypesArray",np.min(edgestype),np.max(edgestype))
            return -1,-1
    return edgestype,encoding

def GetTSL(nodetype):
    """
    Returns the integer values associated with the T, S, and L components in the TSL encoding.
    """
    length=(len(nodetype)-1)//2 # default length is 3, but can be larger.
    if nodetype[0]=='X':
        return -1,-1,-1
    else:
        return int(nodetype[:length]),int(nodetype[length:2*length]),int(nodetype[-1])
    
def GetNodesTypeException(nodetype,depth):
    """
    Check if nodetype exists and return 
         False (no exception,nodetype exists)
         True (exception, nodetype does not exist)
    """    
    if depth==None:
        return False 
     # return False if this node type does not exist
    #length=(len(nodetype)-1)/2
    T,S,L=GetTSL(nodetype)
    if T==-1: #"XXX" return type exception
        return True
    if T>0 and L==0: # can not be a target if L==0
        return True
    return False
    
def GetEdgesTypeException(edgetype,depth):
    """
    Check if edgetype exists and return 
         False (no exception,nodetype exists)
         True (exception, nodetype does not exist)
    """    
    # return False if this edge type does not exist 
    # assuming selfloops are excluded
    if depth==None:
        return False 
    length=int((len(edgetype)-1)/2)
    sourcenodetype=edgetype[:length]
    targetnodetype=edgetype[length+1:]
    if GetNodesTypeException(sourcenodetype,depth) or GetNodesTypeException(targetnodetype,depth):
        return True
    sT,sS,sL=GetTSL(sourcenodetype)
    tT,tS,tL=GetTSL(targetnodetype)
    if sT==-1 or tT==-1:
        # XXX => exception (particular cases for  type encoding 
        return True
    if tL==0:
        return True
    if sS==0:
        return True
    if tT==0:
        return True
        
    return False

            
def GetSourceEdge(nodes,Nedges=None):
    """
    Build the array containing source node index of edges 
    """
    if Nedges==None:
        Nedges=nodes[-1]
    SourceEdges=np.cumsum(np.bincount(nodes[1:-1],minlength=Nedges),dtype='uint32')
    if SourceEdges.shape[0]>Nedges: # take into account cases where last nodes do not have edges
        SourceEdges=SourceEdges[:Nedges] 
    return SourceEdges
    
def GetEdgesTranspose(nodes,argsortedges,Nedges=None):
    """
    Build edge graph transpose using numpy lib
    """
    if Nedges==None:
        Nedges=nodes[-1]
    ti=time.time()
    #print("Edges transpose Start")
    #edges_transpose=np.zeros(Nedges,dtype='uint32')
    #print("Target",Nnodes//int(1e7))
    #nodesistart=nodes[0]
    ##i=np.uint32(0)
    #oneuint32=np.uint32(1)
    #for nodeinext in nodes[1:]:
    #    if i%int(1e7)==0: print(i//int(1e7),end=" ")
    #    edges_transpose[nodesistart:nodeinext]=i
    #    nodesistart=nodeinext
    #    i+=oneuint32
    print("Edges transpose Start")
    edges_transpose=GetSourceEdge(nodes,Nedges) 
    tf=time.time()
    print("Edges transpose Step 1/2 elapse(s)",np.round(tf-ti),"(s)")
    
    ti=time.time()
    edges_transpose=edges_transpose[argsortedges]
    tf=time.time()
    print("Edges transpose Step 2/2 elapse(s)",np.round(tf-ti),"(s)")
    return edges_transpose

def GetNodesTranspose(edges,Nnodes):
    """
    Build node graph transpose using numpy lib
    """
    print("Nodes transpose Start")
    ti=time.time()
    nodes_transpose=np.cumsum(np.bincount(edges+1,minlength=Nnodes+1),dtype='uint32') # +1 to be sure that cumsum start at 0
    tf=time.time()
    print("Nodes transpose elapse(s)",np.round(tf-ti),"(s)")    
    return nodes_transpose

def ArgSort(array,kind="quicksort"):
    """
    wrapping of the numpy argsort function
    """
    #kind{‘quicksort’, ‘mergesort’, ‘heapsort’, ‘stable’}, optional
    #Sorting algorithm. The default is ‘quicksort’. Note that both ‘stable’ and ‘mergesort’ use timsort 
    # under the covers and, in general, the actual implementation will vary with data type. 
    # The ‘mergesort’ option is retained for backwards compatibility.
    return np.argsort(array,kind=kind)

def timestampsarray2yearmonth(arrayTS):
    """
    Convert timestamp array (in seconds since EPOCH) in months since EPOCH
    """
    return (np.array(arrayTS, dtype='datetime64[s]') - np.datetime64('1970-01-01')).astype('timedelta64[M]').astype('int')


    
def GetSourceEdgeTimeStamp(nodes,edges,nodesTS,d):
    """
    Builds an array of timestamps corresponding to the source nodes of each edge.

    Parameters:
        nodes: Node array in CSR format.
        edges: Edge array in CSR format.
        nodesTS: Array of timestamps associated with each node.
        d : not used (deprecated)

    Returns:
        Array of timestamps for the source nodes of all edges.
    """
    return nodesTS[GetSourceEdge(nodes,edges.shape[0])]
    
def GetTargetEdgeTimeStamp(nodes,edges,nodesTS,d):
    """
    Builds an array of timestamps corresponding to the target nodes of each edge.

    Parameters:
        nodes: Node array in CSR format.
        edges: Edge array in CSR format.
        nodesTS: Array of timestamps associated with each node.
        d : not used (deprecated)

    Returns:
        Array of timestamps for the target nodes of all edges.
    """
    return nodesTS[edges]
    

def GetEdgeTs(nodes,edges,nodesad,d,Nnodes,Nedges):
    """
    Build timestamp array of all source nodes and target nodes
    """
    print(80*"-")
    ti=time.time()
    sourcearrayTS=GetSourceEdgeTimeStamp(nodes,edges,nodesad,d)
    targetarrayTS=GetTargetEdgeTimeStamp(nodes,edges,nodesad,d)
    tf=time.time()
    print("Edge Source/Target TS Building ",np.round(tf-ti,2),"(s)")

    # Exclude edges where timestamp==0 or timestamp==2**32-1
    ti=time.time()
    maskEdgeTS=np.logical_and(
        np.logical_and(sourcearrayTS!=0,sourcearrayTS!=2**32-1),
        np.logical_and(targetarrayTS!=0,targetarrayTS!=2**32-1))
    tf=time.time()
    print("Mask Building ",np.round(tf-ti,2),"(s)")
    valid_edges=np.sum(maskEdgeTS)
    print(f'{valid_edges:,} valid edges over a total of {edges.shape[0]:,} i.e. {np.round(valid_edges/edges.shape[0]*100,2)} % of valid edges')
    print(80*"-")

    # ! Timestamp are uint32 (must be converted to int before doing the diff)
    # deltaTS still in second (negative or positive)
    ti=time.time()
    sourcearrayTS=sourcearrayTS[maskEdgeTS]
    targetarrayTS=targetarrayTS[maskEdgeTS]
    deltaTS=sourcearrayTS.astype('int')-targetarrayTS.astype('int')
    tf=time.time()
    print("Building DeltaTS ",np.round(tf-ti,2),"(s)")

    # convert s since EPOCH to Months since EPOCH
    ti=time.time()
    sourcearrayTS=timestampsarray2yearmonth(sourcearrayTS)
    targetarrayTS=timestampsarray2yearmonth(targetarrayTS)
    tf=time.time()
    print("Applying Mask ",np.round(tf-ti,2),"(s)")

    ti=time.time()
    arraytype,encoding=GetEdgesTypesArray(nodes,edges,d)
    arraytype=arraytype[maskEdgeTS]
    tf=time.time()
    print("Building edge type array and applying mask",np.round(tf-ti,2),"(s)")
    print(80*"-")
    return arraytype,encoding,sourcearrayTS,targetarrayTS,deltaTS
