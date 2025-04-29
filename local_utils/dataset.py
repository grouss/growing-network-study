import pickle
import time
import json
import numpy as np



import warnings
from . import config

warnings.filterwarnings("ignore")
libname="dataset"

# Copyright 2025 Université Paris Cité, France
# Author: [Guillaume Rousseau](https://www.linkedin.com/in/grouss/), Physics Department, Paris, France
# This is supplemental materials and replication package associated with the preprint available on 
# - arXiv (https://arxiv.org/abs/2501.10145)
# - SSRN  (http://ssrn.com/abstract=5191689
# Current version of python scripts and ressources are available on [github's author page](https://github.com/grouss/growing-network-study)
# This work is currently licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)
    
def LoadAllArray(transpose=False):
    """
    Return nodes, edges and nodesad array using Compressed Sparse Row format. 
    edges[Nodes[index]:Nodes[index+1]] returns the list of source node indexes for a given target node index.
        Parameters:
            transpose : bool, optional (default=Falss)
                  if equals True will load the transpose.
        Returns:
            nodes : array (1D, len=Nnodes+1)
            edges : array (1D, len=Nedges)
            nodesad : array (1D, len=Nodes)
            Nnodes : int 
                     Number of nodes.
            Nedges : int
                     Number of edges.
    """
    if transpose:
        transpose="transpose_"
    else:
        transpose=""
    try:
        filename_nodes=config.graphpath+"nodes_"+transpose+"20240310.pkl"
        nodes=pickle.load(open(filename_nodes,"rb"))
        print("Loaded :", filename_nodes)

        filename_edges=config.graphpath+"edges_"+transpose+"20240310.pkl"
        edges=pickle.load(open(filename_edges,"rb"))
        print("Loaded :", filename_edges)
        
        filename_nodesad=config.graphpath+"nodesad_20240310.pkl"
        nodesad=pickle.load(open(filename_nodesad,"rb"))
        print("Loaded :", filename_nodesad)

    except:
        print("Did you download dataset files from https://zenodo.org/records/15260640")
        print("If not, download it and put all files in your import data directory")
        print("that can be defined in ./local_utils/config.py")
        if transpose=="transpose_":
            print("Transpose graph is not provided in the zenodo deposit")
            print("use SM03_BuildingTransposeGraph.ipynb to build it")
        return -1
        
    d=GetD()
    Nnodes=nodes.shape[0]-1
    Nedges=edges.shape[0]
    return nodes,edges,nodesad,d,Nnodes,Nedges    
 
def LoadAllArray_OO(keypath="BigO"):
    """
    Same as LoadAllArray, but for the derived O-(RV/RL)-O network.
    Return nodes, edges and nodesad array using Compressed Sparse Row format. 
    edges[Nodes[index]:Nodes[index+1]] returns the list of source node indexes for a given target node index.
        Parameters:
            keypath : string, optional (default="BigO")
                  "BigO" => with inheritance. "L0" => without inheritance.
        Returns:
            nodes : array (1D, len=Nnodes+1)
            edges : array (1D, len=Nedges)
            nodesad : array (1D, len=Nodes)
            Nnodes : int 
                     Number of nodes.
            Nedges : int
                     Number of edges.
    """
    try:
        
        filename_nodes="nodes_o_derived_O-RVRL-O_"+keypath+"_20240429.pkl"
        nodes=pickle.load(open(config.graphpath+filename_nodes,"rb"))
        print("Loaded :", filename_nodes)
        
        filename_edges="edges_o_derived_O-RVRL-O_"+keypath+"_20240429.pkl"
        edges=pickle.load(open(config.graphpath+filename_edges,"rb"))
        print("Loaded :", filename_edges)
        
        # nodesadderived only depends on temporal partitionning, not on the derived/inheritance policy
        filename_nodesad="nodesadderived_derived_O-RVRL-O_BigO_20240429.pkl"
        nodesad=pickle.load(open(config.graphpath+filename_nodesad,"rb"))
        print("Loaded :", filename_nodesad)
        
    except:
        print("Did you download dataset files from https://zenodo.org/records/15260640")
        print("If not, download it and put all files in your import dayta directory")
        print("that can be defined in ./local_utils/config.py")
        return -1
    
    d=GetD()
    Nnodes=nodes.shape[0]-1
    Nedges=edges.shape[0]
    return nodes,edges,nodesad,d,Nnodes,Nedges  

def GetD(Verbose=False):
    """
    This function return index min and index max for each node types, and total number of nodes of each types. 
    This is specific to the graph representation and dataset used because nodes'type encoding assumes node ordering by type.
        Parameters:
            Verbose: Default=False. 
                     Display info 
        Returns:
            dict: key:int. 
                  keys are of the forme TYPEindexMin, TYPEindexMax, TYPE
    """
    
    d={
       "RVindexMin": 156799786,
       "RVindexMax": 2224378839,
       "RV": 2067579054,
       "OindexMin": 0,
       "OindexMax": 139524532,
       "O": 139524533,
       "RLindexMin": 139524533,
       "RLindexMax": 156799785,
       "RL": 17275253
        }
    if Verbose:
        print(json.dump(d,depth=3))
    # Origins            139 524 533
    # Revisions        2 067 579 054
    # Releases            17 275 253
    return d


def SetParamsBA_d(NnodesBA):
    """
    Set Params similar to the network studyied here
    """
    dBA={}
    for ntype in ["O","RV","RL"]:
        dBA[ntype]=0
        dBA[ntype+"indexMin"]=0
        dBA[ntype+"indexMax"]=0

    dBA["O"]=NnodesBA
    dBA["OindexMin"]=0
    dBA["OindexMax"]=NnodesBA-1


    for key in ["RV","RL"]:
        dBA[key]=0
        dBA[key+"indexMax"]=NnodesBA
        dBA[key+"indexMin"]=NnodesBA
    return dBA

def SetParamsBA_nodesad(NnodesBA):
    """
    Set Params similar to the network studyied here
    """
    nodesadBA=1+(2021-1970)*365*24*3600*np.log(np.arange(1,NnodesBA+1,dtype='float'))/np.log(NnodesBA)
    nodesadBA=nodesadBA.astype('uint32')
    return nodesadBA