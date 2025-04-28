import pickle
import time
import json
import numpy as np


import warnings
#from . import config
from .dataset import *
from .graph import *

warnings.filterwarnings("ignore")
libname="stat"

# Copyright 2025 Université Paris Cité, France
# Author: [Guillaume Rousseau](https://www.linkedin.com/in/grouss/), Physics Department, Paris, France
# This is supplemental materials and replication package associated with the preprint available on 
# - arXiv (https://arxiv.org/abs/2501.10145)
# - SSRN  (http://ssrn.com/abstract=5191689
# Current version of python scripts and ressources are available on [github's author page](https://github.com/grouss/growing-network-study)
# This work is currently licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)



def GetXYhist(edgesindexarray,Nnodes):
    """
    Computes the histogram and complementary cumulative distribution function (CCDF)
    of node degrees based on an edge index array.
        Parameters:
            edgesindexarray: Array of node indices (e.g. edge sources or targets).
            Nnodes (int): Total number of nodes in the graph.

        Returns:
            dict: A dictionary with the following keys:
                - 'x': numpy array of unique degree values.
                - 'y': numpy array  of the frequency of each degree (i.e. histogram).
                - 'ccdf': numpy array  of the complementary cumulative distribution function,
                          computed as the reverse cumulative sum of 'y'.
    """
    bincountarray=np.bincount(np.bincount(edgesindexarray,minlength=Nnodes))
    x=np.arange(bincountarray.shape[0])
    mask=bincountarray!=0
    y=bincountarray[mask]
    return {"x":x[mask],"y":y,"ccdf":np.cumsum(y[::-1])[::-1]}


def convertdict(dic):
    """
    This function convert an histogram using dictionnary to numpy array
        Parameters:
            dic: dictionnary {x_i:y_i} (not necessary sorted)
        Returns:
            x  : numpy array of sorted x_i/key from dic
            y  : numpy array of associated y_i/value from dic
            yc : complementary cumulative distribution fonction
            
    """
    x=[]
    y=[]
    for key in sorted(dic.keys()):
        x.append(key)
        y.append(dic[key])
    x=np.array(x)
    y=np.array(y)
    yc=np.cumsum(y[::-1])[::-1]
    
    # similar to / to be updated with 
    # xy=np.array(list(dic.items()))
    # x=xy[:,0]
    # y=xy[:,1]
    # mask=y!=0
    # x=x[mask]
    # y=y[mask]
    # reordering=np.argsort(x)
    # x=x[reordering]
    # y=y[reordering]
    # yc=np.cumsum(y[::-1])[::-1]
    
    return x,y,yc

def BinArray2dict(array):
    """
    Converts an array into a dictionary by mapping each non-zero value to its index.
        Parameters:
            array: Input array of values.
        Returns:
            dict: Dictionary where keys are indices and values are the corresponding non-zero elements of the array.
    """
    dic={i:array[i] for i in range(len(array)) if array[i]!=0}
    return dic


def DisplayTypeStats(nodes,edges,d,depth=None,Debug=False,FilePath=None):
    """
    Computes and displays the frequency statistics of types for nodes and edges.
    Depending on the presence of the `depth` parameter, the function either computes simple type
    statistics or TSL-encoded types with depth-based context.
    For edge types, it separates standard edges from self-loops, and optionally exports the results.
        Parameters:
            nodes: Node array in CSR format.
            edges: Edge array in CSR format.
            d: dictionnary used for type encoding .
            depth (int, optional): If provided, enables TSL encoding with depth context. Default is None.
            Debug (bool, optional): If True, prints additional information such as ignored types. Default is False.
            FilePath (str, optional): If provided, saves the output dictionary as a pickle file to this path.
        Returns:
            dict: Dictionary mapping each TSL type (or equivalent) to its frequency in the graph.
        Side effects:
            - Prints a summary of type frequencies and their relative percentages.
            - Performs integrity checks on type counts.
            - Prints timing information.
            - If Debug is True, prints the list of excluded types.
            - If FilePath is set, saves the resulting stats dictionary as a pickle file.
    """
    
    def GetStat(typearray,encoding):
        stats=np.bincount(typearray,minlength=len(encoding))
        #return {encoding[i]:int(stats[i]) for i in range(len(encoding)) if stats[i]!=0}
        return {encoding[i]:int(stats[i]) for i in range(len(encoding)) if encoding[i]!='X'*len(encoding[0])}
    
    statsoutput={}
    # We use the same function 
    if depth==None:
        Liste=[
         (f'GetNodesTypesArray',GetNodesTypesArray,(nodes,edges,d),len(nodes)-1,GetNodesTypeException,False)
        ,(f'GetEdgesTypesArray',GetEdgesTypesArray,(nodes,edges,d),len(edges),GetEdgesTypeException,False)
                        ]
    else:
        Liste=[
         (f'GetNodesTypesArray',GetNodesTypesArray,(nodes,edges,d,depth,Debug),len(nodes)-1,GetNodesTypeException,False)
        ,(f'GetEdgesTypesArray',GetEdgesTypesArray,(nodes,edges,d,depth,Debug),len(edges),GetEdgesTypeException,True)
        ]
    exception=[]                   
    for field,f,para,total,TypeException,Fmask in Liste:
        totaltmp=0
        statsselfloop={}
        
        ti=time.time()
        
        typearray,encoding=f(*para)
        
        if Fmask and field==f'GetEdgesTypesArray':
            # exception for selfloop management
            mask=GetSourceEdge(nodes)==edges
            statsselfloop=GetStat(typearray[mask],encoding)
            stats=GetStat(typearray[np.logical_not(mask)],encoding)
        else:
            stats=GetStat(typearray,encoding)
    
        tf=time.time()
        print(field+f' [Elapse time : {np.round(tf-ti)} (s)]')
        for key,value in stats.items():
            if not TypeException(key,depth):
                totaltmp+=value
                print(f'___ {key:5} : {value:15,} ({np.round(100*value/np.sum(list(stats.values())),2)}%)')
                statsoutput[key]=value
            else:
                exception.append((key,value))
        print(f'_'*30) 
        print(f'___ Total : {totaltmp:15,} ({np.round(100*totaltmp/np.sum(list(stats.values())),2)}%)')        
        if Fmask and field==f'GetEdgesTypesArray':
            print()
            print(field+f' Self Loop [Elapse time : {np.round(tf-ti)} (s)]')
            totalselflooptmp=0
            for key,value in statsselfloop.items():
                if value!=0:
                    totalselflooptmp+=value
                    print(f'{key.replace(">","=")} {value:15,} ({np.round(100*value/np.sum(list(statsselfloop.values())),2)}%)')
                    statsoutput[key.replace(">","=")]=value
            print(f'_'*30) 
            print(f'___ Total : {totalselflooptmp:15,} ({np.round(100*totalselflooptmp/np.sum(list(statsselfloop.values())),2)}%)') 
            
            # integrity check
            if totaltmp+totalselflooptmp!=total:
                print("ERROR total",field,np.sum(list(stats.values())),
                      totaltmp,np.sum(list(statsselfloop.values())),totaltmp,totalselflooptmp,total)
        else:
            
            # integrity check
            if totaltmp!=total:
                print("ERROR total",field,np.sum(list(stats.values())),
                      totaltmp,np.sum(list(statsselfloop.values())),totaltmp,total)

        if Debug:
            print("Debug",exception)
            
        print()
        
        if FilePath!=None:
            pickle.dump(statsoutput,open(FilePath,"wb"))
            
    return statsoutput
            

def BuildNodesTimeStampHisto(nodes,edges,nodesTS,d,stat={},depth=None):
    """
    Computes an histogram of node types based on their timestamps (months since EPOCH).
        Parameters:
            nodes: Nodes in CSR format.
            edges: Edges in CSR format.
            nodesTS: Node timestamps.
            d (dict): Parameter for type classification.
            stat (dict, optional): Dictionary to fill with results (default: empty).
            depth (int, optional): Depth for TSL encoding (if applicable).
        Returns:
            dict: Dict with counts per months since EPOCH and type string.
    """
    nodestype,encoding=GetNodesTypesArray(nodes,edges,d,depth=depth)
    nodesTSM=timestampsarray2yearmonth(nodesTS)
    TSM2stat(nodesTSM,nodestype,encoding,stat=stat)
    return stat
        
def TSM2stat(TSM,arraytype,encoding,stat={}):
    """
    Aggregates counts of node types over time into a flat histogram.
        Parameters:
            TSM: Timestamps mapped to months since EPOCH integers.
            arraytype: Array of type indices (corresponding to positions in `encoding`).
            encoding: Array of type labels.
            stat (dict, optional): Dictionary to fill with results. Default is empty.
        Returns:
            dict: Dictionary mapping each type label to an array of counts over time.
    """
    TSmax=1634 # =1+max(nodesTS)=1+1633
    TSM=TSM+TSmax*arraytype
    TSM=np.bincount(TSM,minlength=TSmax*len(encoding))
    for i,Ntype in enumerate(encoding):
        stat[Ntype]=TSM[TSmax*i:TSmax*(i+1)]
    return stat
    
def BuildEdgesTimeStampHisto(nodes,edges,nodesTS,d,stat={},depth=None,Verbose=False):
    """
    Builds a histogram of edge types based on the timestamps of their source and target nodes.
        Parameters:
            nodes: Nodes in CSR format.
            edges: Edges in CSR format.
            nodesTS: Timestamps associated with each node.
            d (dic): Parameter used for edge type classification.
            stat (dict, optional): Dictionary to accumulate time-based type statistics. Default is empty.
            depth (int, optional): Enables TSL-based edge encoding. If set, self-loops are excluded.
            Verbose (bool, optional): If True, prints timing information for each step.
        Returns:
            dict: Dictionary with edge type counts per months since EPOCH.
    """

    # build edges timestamp
    # we use edgesTSM for all temporary array
    if Verbose: ti=time.time()
    edgesTSM=GetSourceEdgeTimeStamp(nodes,edges,nodesTS,d)
    if Verbose:
        tf=time.time()
        print("edgesTSM 1/2 ",np.round(tf-ti),"(s)")
        ti=tf
    #source time stamp (max TS from source and edge)
    edgesTSM=timestampsarray2yearmonth(np.maximum(
        GetSourceEdgeTimeStamp(nodes,edges,nodesTS,d),
        GetTargetEdgeTimeStamp(nodes,edges,nodesTS,d)
        ))
    if Verbose:
        tf=time.time()
        print("edgesTSM 2/2 ",np.round(tf-ti),"(s)")
        ti=tf
    # edges type
    edgestype,encoding=GetEdgesTypesArray(nodes,edges,d,depth=depth)
    if type(depth)!=type(None):
        # excluding self loop
        mask=edges!=GetSourceEdge(nodes)
        edgesTSM=edgesTSM[mask]
        edgestype=edgestype[mask]
    if Verbose:
        tf=time.time()
        print("edgestype 1/1 ",np.round(tf-ti),"(s)")
        ti=tf
    TSM2stat(edgesTSM,edgestype,encoding,stat=stat)
    if Verbose:
        tf=time.time()
        print("TSM2stat 1/1 ",np.round(tf-ti),"(s)")
        ti=tf
    return stat

def DisplayTimestampException(arraytype,encoding,arrayTS,TypeString,Verbose=False):
    """
    Displays timestamp anomalies (zero or max value) for each type in a typed array.
    Specifically checks for:
    - Timestamps equal to 0 (interpreted as 1970-01-01).
    - Timestamps equal to 2^32 - 1 (interpreted as 2106-02-07).
        Parameters:
            arraytype: Array of type indices (corresponding to positions in `encoding`).
            encoding: Array of type labels.
            arrayTS: Array of timestamps to check.
            TypeString (str): Label to prefix in the printed output (e.g., "Node" or "Edge").
            Verbose (bool, optional): If True, prints elapsed time. Default is False.
        Returns:
            None
    """

    ti=time.time()
    mask_zero=arrayTS==0
    mask_max=arrayTS==2**32-1
    for i,Ntype in enumerate(encoding):
        mask_type=arraytype==i
        Ntot=np.sum(mask_type)
        if Ntot!=0:
            Nzero=np.sum(np.logical_and(mask_type,mask_zero))
            Nmax=np.sum(np.logical_and(mask_type,mask_max))
            print(f'{TypeString:12} Timestamp Exceptions {Ntype:5} | Total={Ntot:15,} | Zero(1970/1/1)={Nzero:15,} | Ts(2106/2/7)={Nmax:15,}')
    tf=time.time()
    
    if Verbose: print("Elapse time ",np.round(tf-ti),"(s)")



def GetDegreeStats(sourceedges,targetedges,TSY,s,dout,din,Nnodes,
                   #YearBegin=1990,YearEnd=2030,YearSlice=40,FlagMonth=False,MonthSlice=1):
                   YearBegin=1980,YearEnd=2025,YearSlice=1,FlagMonth=False,TSM=None,MonthSlice=1,Verbose=False):
    """
    Computes and stores time-sliced in-degree and out-degree histograms for a given entity.

    Depending on the `FlagMonth` setting, the function aggregates degree statistics either by year 
    or by months since EPOCH, using timestamp-based filtering on the edge arrays.
        Parameters:
            sourceedges: Array of source node indices of edges.
            targetedges: Array of target node indices of edges.
            TSY: Timestamps in years relative to 1970 (used when FlagMonth=False).
            s: Key string used to store results in the output dictionaries.
            dout: Output dictionary for out-degree statistics (modified in place).
            din: Output dictionary for in-degree statistics (modified in place).
            Nnodes: Total number of nodes in the graph.
            YearBegin (int, optional): Starting year for aggregation (default: 1980).
            YearEnd (int, optional): Ending year for aggregation (default: 2025).
            YearSlice (int, optional): Step between years (default: 1).
            FlagMonth (bool, optional): If True, aggregates by months since EPOCH instead of by year.
            TSM (array-like, optional): Timestamps in months since Jan 1970 (required if FlagMonth=True).
            MonthSlice (int, optional): Step between months (default: 1).
            Verbose (bool, optional): If True, prints progress and timing info.

        Returns:
            None: Results are stored in-place in `dout[s]` and `din[s]`.
    """
    ti=time.time()
    dout[s]={}
    din[s]={}
    if FlagMonth:
        Nstep=max(1,((YearEnd-YearBegin)*12)//10)
        for count,Month in enumerate(range(12*YearBegin,12*(YearEnd+1),MonthSlice)):
            tf=time.time()
            if Verbose: 
                print(Month//12,Month%12,"(",int(tf-ti),"s)",end=" ")
            else:
                if count%Nstep==0: 
                    print(Month//12,Month%12,end=" ")
            mask=TSM<(Month-1970*12) # <= 1/Month/YEAR
            dout[s][Month]=GetXYhist(sourceedges[mask],Nnodes)
            din[s][Month]=GetXYhist(targetedges[mask],Nnodes)   
        # Year=3000 => no constrainte
        dout[s][3000*12]=GetXYhist(sourceedges,Nnodes)
        din[s][3000*12]=GetXYhist(targetedges,Nnodes)   
    else:
        Nstep=max(1,(YearEnd-YearBegin)//10)
        for count,Year in enumerate(range(YearBegin,YearEnd+1,YearSlice)):
            tf=time.time()
            if Verbose: 
                print(Year,"(",int(tf-ti),"s)",end=" ")
            else:
                if count%Nstep==0: 
                    print(Year,end=" ")
            mask=TSY<(Year-1970) # < 1/1/YEAR
            dout[s][Year]=GetXYhist(sourceedges[mask],Nnodes)
            din[s][Year]=GetXYhist(targetedges[mask],Nnodes)  
        # Year=3000 => no constrainte
        dout[s][3000]=GetXYhist(sourceedges,Nnodes)
        din[s][3000]=GetXYhist(targetedges,Nnodes)  
    print()
    
def GetAllTypesDegreeStats(sourceedges,targetedges,TSY,dout,din,Nnodes,FlagMonth=False,TSM=None):
    """
    Computes degree statistics over time for all node types combined.
    Calls `GetDegreeStats` with a fixed label ("All types") and stores the resulting
    in-degree and out-degree histograms into the provided `dout` and `din` dictionaries.
        Parameters:
            sourceedges: Array of source node indices of edges.
            targetedges: Array of target node indices of edges.
            TSY: Timestamps in years relative to 1970 (used when FlagMonth=False).
            s: Key string used to store results in the output dictionaries.
            dout: Output dictionary for out-degree statistics (modified in place).
            din: Output dictionary for in-degree statistics (modified in place).
            Nnodes: Total number of nodes in the graph.
            FlagMonth (bool, optional): If True, aggregates by months since EPOCH instead of by year.
            TSM (array-like, optional): Timestamps in months since Jan 1970 (required if FlagMonth=True).

        Returns:
            None: Results are stored in-place in `dout["All types"]` and `din["All types"]`.
    """
    s="All types"
    print("Start "+s)
    ti=time.time()
    GetDegreeStats(sourceedges,targetedges,TSY,s,dout,din,Nnodes,FlagMonth=FlagMonth,TSM=TSM)
    tf=time.time()
    print("All types elapse : ",np.round(tf-ti,2),"(s)")
    
def GetPerTypesDegreeStats(sourceedges,targetedges,TSY,Edgestype,Edgesencoding,dout,din,Nnodes,
                           PerFlag="EDGE",FlagMonth=False,TSM=None):
    """
    Computes time-sliced degree statistics for each edge type (or source/target projection).
    For each unique edge type (or aggregated type if PerFlag is "SOURCE" or "TARGET"), the function
    filters the corresponding edges and computes in-degree and out-degree histograms over time.
        Parameters:
            sourceedges: Array of source node indices of edges.
            targetedges: Array of target node indices of edges.
            TSY: Timestamps in years relative to 1970 (used when FlagMonth=False).
            arraytype: Array of type indices (corresponding to positions in `encoding`).
            encoding: Array of type labels.
            dout: Output dictionary for out-degree statistics (modified in place).
            din: Output dictionary for in-degree statistics (modified in place).
            Nnodes: Total number of nodes in the graph.
            PerFlag (str, optional): Aggregation mode. One of {"EDGE", "SOURCE", "TARGET"}.
            FlagMonth (bool, optional): If True, aggregates by months since EPOCH instead of by year.
            TSM (array-like, optional): Timestamps in months since Jan 1970 (required if FlagMonth=True).
        Returns:
            None: Updates `dout` and `din` in-place with time-sliced stats for each edge type or group.
    """
    s="Per "+PerFlag+" types "
    print("Start "+s)
    ti=time.time()
    alreadydone=[]
    for indextype,edgetype in enumerate(Edgesencoding):
        if edgetype not in alreadydone:
            if PerFlag=="EDGE":
                stmp=edgetype
                masktype=Edgestype==indextype
                alreadydone.append(edgetype)
            elif PerFlag=="SOURCE":
                stmp=edgetype.split('>')[0]+">"
                # on agrège tous les types correspondant
                masktype=np.zeros_like(Edgestype,dtype='bool')
                for indextype2,edgetype2 in enumerate(Edgesencoding):
                    stmp2=edgetype2.split('>')[0]+">"
                    if stmp2==stmp:
                        masktype=np.logical_or(masktype,Edgestype==indextype2)
                        alreadydone.append(edgetype2)
            elif PerFlag=="TARGET":
                stmp=">"+edgetype.split('>')[1]
                # on agrège tous les types correspondant
                masktype=np.zeros_like(Edgestype,dtype='bool')
                for indextype2,edgetype2 in enumerate(Edgesencoding):
                    stmp2=">"+edgetype2.split('>')[1]
                    if stmp2==stmp:
                        masktype=np.logical_or(masktype,Edgestype==indextype2)
                        alreadydone.append(edgetype2)            
            else:
                print("ERROR PerFlag not supported")    
                return
            ncount=np.sum(masktype)
            if ncount!=0:
                print(f'Start type {stmp} / {ncount:,} edges',end=" | ")
                sourceedgestmp=sourceedges[masktype]
                targetedgestmp=targetedges[masktype]
                if not FlagMonth:
                    TSYtmp=TSY[masktype]
                    GetDegreeStats(sourceedgestmp,targetedgestmp,TSYtmp,stmp,dout,din,Nnodes,FlagMonth=False,TSM=None)
                else:
                    TSMtmp=TSM[masktype]
                    GetDegreeStats(sourceedgestmp,targetedgestmp,None,stmp,dout,din,Nnodes,FlagMonth=True,TSM=TSMtmp)
    tf=time.time()
    print("Per edges types elapse : ",np.round(tf-ti,2),"(s)")

def GetNodeEdgeStat(nodesTM,nodesDegree,debug=False):
    """ 
    Computes the number of nodes and outgoing edges per month.
    Assumes that `nodesTM` and `nodesDegree` are sorted by timestamp (in months since the epoch).
        Parameters:
            nodesTM (array-like): Array of node timestamps expressed in months since epoch.
            nodesDegree (array-like): Out-degree of each node, in the same order as `nodesTM`.
            debug (bool, optional): If True, performs consistency checks and prints errors if detected.
        Returns:
            tuple:
                - nodesStat (np.ndarray): Array of node counts per month (length = 1634).
                - edgesStat (np.ndarray): Array of total outgoing edges per month.
    """
    # nodesTM,nodesDegree has to be sorted according to nodesTM
    # nodesTM is the timestamp converted to number of month since EPOCH 
    # max using timestamp in second since EPOCH and uiint32 => [0,1633]
    nodesStat=np.bincount(nodesTM,minlength=1634)
    degreeoutcumsum=np.zeros(1+len(nodesDegree),dtype='int')
    #argsortarray=np.argsort(nodesTS)
    #degreeoutcumsum[1:]=np.cumsum(nodesDegree[np.argsort(nodesTS)])
    degreeoutcumsum[1:]=np.cumsum(nodesDegree)
    indexcumsum=np.cumsum(nodesStat)
    degreeoutcumsum=degreeoutcumsum[indexcumsum]
    edgesStat=np.zeros_like(degreeoutcumsum)
    edgesStat[0]=degreeoutcumsum[0]
    edgesStat[1:]=degreeoutcumsum[1:]-degreeoutcumsum[:-1]
    if debug:
        if np.sum(nodesDegree)!=np.sum(edgesStat):
            print("ERROR 1",np.sum(nodesDegree),np.sum(edgesStat))
        if nodesTM.shape[0]!=np.sum(nodesStat):
            print("ERROR 2",nodesTM.shape[0],np.sum(nodesStat))
        # check boundary XXX
        # check min and max XXX
    return nodesStat,edgesStat

