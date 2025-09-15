import os,gc

from . import config

from .graph import *
from .stat import *
from .plot import *

def GetArrayFromSandBox(name,Verbose=False,FilterException=False):
    filename=config.exportpath+"sandbox/"+name+".pkl"
    if Verbose:
        print("Load ",name,"from sandbox")
    if not FilterException:
        return pickle.load(open(filename,"rb"))
    elif "node" in name:
        return pickle.load(open(filename,"rb"))[GetTS_Exception("node",Verbose=Verbose)]
    elif "edge" in name:
        return pickle.load(open(filename,"rb"))[GetTS_Exception("edge",Verbose=Verbose)]
        
    print("ERROR GetArrayFromSandBox",name,"Verbose",Verbose,"FilterException",FilterException)
    
def SaveArrayToSandBox(array,name,Verbose=False):
    filename=config.exportpath+"sandbox/"+name+".pkl"
    if Verbose:
        print("Save ",name,"from sandbox")
    return pickle.dump(array,open(filename,"wb"))
    
def GetTSMFromArrayTS(name,EPOCH,FilterException=False):
    return timestampsarray2yearmonth(GetArrayFromSandBox(name+"TS",FilterException=FilterException),EPOCH=EPOCH)
        
def GetTS_Exception(name,Verbose=False):
    return GetArrayFromSandBox(name+"TS_Exception",Verbose=Verbose)

def ReportMain_A(nodes,edges,nodesad,d,EPOCH='1970-01-01',FSaveTemporary=True,Verbose=False):

    DisplayTypeStats(nodes,edges,d)   
    lstring=140
    
    if FSaveTemporary:
        sandboxpath=config.exportpath+"sandbox/"
        if not os.path.exists(sandboxpath):
            os.makedirs(sandboxpath)
            print(f"Create '{sandboxpath}'")
        
    print(f'_'*lstring)

    # Exception Timestamp nodes per type
    arraytype,encoding=GetNodesTypesArray(nodes,edges,d)
    
    
    if FSaveTemporary:
        SaveArrayToSandBox(arraytype,"nodetype",Verbose=Verbose)
        SaveArrayToSandBox(encoding,"nodeencoding",Verbose=Verbose)
        SaveArrayToSandBox(nodesad,"nodeTS",Verbose=Verbose)
        
        #pickle.dump(arraytype,open(sandboxpath+"nodetype.pkl","wb"))
        #pickle.dump(encoding,open(sandboxpath+"nodeencoding.pkl","wb"))
        #pickle.dump(nodesad,open(sandboxpath+"nodeTS.pkl","wb"))
        #nodeTSM=timestampsarray2yearmonth(nodesad,EPOCH=EPOCH)
        #pickle.dump(nodeTSM,open(sandboxpath+"nodeTSM.pkl","wb"))


    DisplayTimestampException(arraytype,encoding,nodesad,f'Nodes',EPOCH=EPOCH)
    del arraytype,encoding;gc.collect()
    #del nodeTSM
    
        
    print(f'_'*lstring)

    # edge source timestamp exception
    arraytype,encoding=GetEdgesTypesArray(nodes,edges,d)
    if FSaveTemporary:
        SaveArrayToSandBox(arraytype,"edgetype",Verbose=Verbose)
        SaveArrayToSandBox(encoding,"edgeencoding",Verbose=Verbose)
        #pickle.dump(arraytype,open(sandboxpath+"edgetype.pkl","wb"))
        #pickle.dump(encoding,open(sandboxpath+"edgeencoding.pkl","wb"))    
    
    arrayTS=GetSourceEdgeTimeStamp(nodes,edges,nodesad,d)
    
    if FSaveTemporary:
        SaveArrayToSandBox(arrayTS,"sourceedgeTS",Verbose=Verbose)

        #pickle.dump(arrayTS,open(sandboxpath+"sourceedgeTS.pkl","wb"))
        #sourceedgeTSM=timestampsarray2yearmonth(arrayTS,EPOCH=EPOCH)
        #pickle.dump(sourceedgeTSM,open(sandboxpath+"sourceedgeTSM.pkl","wb"))
        #del sourceedgeTSM
        
    DisplayTimestampException(arraytype,encoding,arrayTS,f'Source Edges',EPOCH=EPOCH)
    del arrayTS;gc.collect()
    
    print(f'_'*lstring)

    # edge target timestamp exception
    # arraytype already known
    #arraytype,encoding=GetEdgesTypesArray(nodes,edges,d)    
    arrayTS=GetTargetEdgeTimeStamp(nodes,edges,nodesad,d)
    
    if FSaveTemporary:
        SaveArrayToSandBox(arrayTS,"targetedgeTS",Verbose=Verbose)

        #pickle.dump(arrayTS,open(sandboxpath+"targetedgeTS.pkl","wb"))
        #targetedgeTSM=timestampsarray2yearmonth(arrayTS,EPOCH=EPOCH)
        #pickle.dump(targetedgeTSM,open(sandboxpath+"targetedgeTSM.pkl","wb"))
        #del targetedgeTSM
    
    DisplayTimestampException(arraytype,encoding,arrayTS,f'Target Edges',EPOCH=EPOCH)
    del arrayTS,arraytype,encoding;gc.collect()
    
    if FSaveTemporary:
        # split the building to reduce memory usage
        sourceedgeTS=GetArrayFromSandBox("sourceedgeTS",Verbose=Verbose)
        #sourceedgeTS=pickle.load(open(config.exportpath+"sandbox/sourceedgeTS.pkl","rb"))
        edgeDeltaTS=sourceedgeTS.astype('int')
        del sourceedgeTS;gc.collect()
        targetedgeTS=GetArrayFromSandBox("targetedgeTS",Verbose=Verbose)
        #targetedgeTS=pickle.load(open(config.exportpath+"sandbox/targetedgeTS.pkl","rb"))
        edgeDeltaTS-=targetedgeTS.astype('int')
        del targetedgeTS;gc.collect()
        SaveArrayToSandBox(edgeDeltaTS,"edgeDeltaTS",Verbose=Verbose)
        #pickle.dump(edgeDeltaTS,open(sandboxpath+"edgeDeltaTS.pkl","wb"))
        del edgeDeltaTS;gc.collect()
    if FSaveTemporary:
        # export nodes and edges list to filter timestamp exception (by default 0 and 2**32-1)
        #nodeTS=GetArrayFromSandBox("nodeTS")
        NodeTS_Exception=np.ones(len(nodes)-1,dtype='bool')
        NodeTS_Exception[nodesad==0]=False
        NodeTS_Exception[nodesad==np.uint(2**32-1)]=False
        SaveArrayToSandBox(NodeTS_Exception,"nodeTS_Exception",Verbose=Verbose)
        del NodeTS_Exception; gc.collect()

        edgeTS_Exception=np.ones(len(edges),dtype='bool')

        sourceedgeTS=GetArrayFromSandBox("sourceedgeTS",Verbose=Verbose)
        edgeTS_Exception[sourceedgeTS==0]=False
        edgeTS_Exception[sourceedgeTS==np.uint(2**32-1)]=False
        del sourceedgeTS; gc.collect()

        targetedgeTS=GetArrayFromSandBox("targetedgeTS",Verbose=Verbose)
        edgeTS_Exception[targetedgeTS==0]=False
        edgeTS_Exception[targetedgeTS==np.uint(2**32-1)]=False
        del targetedgeTS; gc.collect()       
        
        SaveArrayToSandBox(edgeTS_Exception,"edgeTS_Exception",Verbose=Verbose)
        del edgeTS_Exception; gc.collect()

        
    print(f'_'*lstring)

def ReportMain_B(Nnodes,Nedges,
                 EPOCH='1970-01-01',
                 yearplain=[],yeardash=[],yearlongdash=[],
                 yscale="log",figsize=(18,4),FOnly=False,Verbose=False,FilterException=True):    
    
    #nodetype,nodeencoding=pickle.load(open(config.exportpath+"sandbox/nodetype.pkl","rb"))
    #nodetype=pickle.load(open(config.exportpath+"sandbox/nodetype.pkl","rb"))
    nodeencoding=GetArrayFromSandBox("nodeencoding")
    edgeencoding=GetArrayFromSandBox("edgeencoding")
        
        
    if len(nodeencoding)==1:
        print("Encoding of the single type of nodes",nodeencoding)
        todolistNodes=list(enumerate(nodeencoding))
        todolistEdges=list(enumerate(edgeencoding))
    elif FOnly:
        # aggregated stats only
        print("Encoding of the different type of nodes",nodeencoding)
        print("Aggregated stat ONLY")
        todolistNodes=[(-1,'"All Node types"')]
        todolistEdges=[(-1,'"All Edge types"')]
    
    else:
        # aggregated stats first
        print("Encoding of the different type of nodes",nodeencoding)
        print("Aggregated stat first")
        todolistNodes=[(-1,'"All Node types"')]+list(enumerate(nodeencoding))
        todolistEdges=[(-1,'"All Edge types"')]+list(enumerate(edgeencoding))
        
    print("-"*80)
    print("Report will display the different types of nodes in the following order")
    print(todolistNodes)
    print("-"*80)

    #nodetype=GetArrayFromSandBox("nodetype")
    
    if Verbose: print("Start todolistNodes")
    for i,NodeType in todolistNodes:
        if i==-1:
            nodeTSM=GetTSMFromArrayTS("node",EPOCH,FilterException=FilterException)
        else:
            nodeTSM=GetTSMFromArrayTS("node",EPOCH,FilterException=FilterException)[GetArrayFromSandBox("nodetype",FilterException=FilterException)==i]
            
        if len(nodeTSM)!=0:   
            print("Start for nodes of type :",NodeType)
            print(80*"-")
            PlotTSoverTimeNodes(nodeTSM,
                            NodeType,
                            EPOCH=EPOCH,Xmin=int(EPOCH[:4]),
                            yearplain=yearplain,yeardash=yeardash,yearlongdash=yearlongdash,
                            yscale=yscale,figsize=figsize)
        del nodeTSM 
        gc.collect()
        if Verbose: print("End NodeType in todolistNodes")
    
  
  
    print("-"*80)
    print("Report will display the different types of edges in the following order")
    print(todolistEdges)
    print("-"*80)
    
    if Verbose: print("Start todolistEdges")

    for i,EdgeType in todolistEdges:
        if i!=-1:
            
            NodeType=EdgeType.split(">")[0]
            sourceedgeTSM=GetTSMFromArrayTS("sourceedge",EPOCH,FilterException=FilterException)[GetArrayFromSandBox("edgetype",FilterException=FilterException)==i ]
            
            iNodeType=np.where(np.array(nodeencoding)==NodeType)[0][0]
            nodeTSM=GetTSMFromArrayTS("node",EPOCH,FilterException=FilterException)[GetArrayFromSandBox("nodetype",FilterException=FilterException)==iNodeType]
            
            
        else:
            nodeTSM=GetTSMFromArrayTS("node",EPOCH,FilterException=FilterException)
            sourceedgeTSM=GetTSMFromArrayTS("sourceedge",EPOCH,FilterException=FilterException)         
            
        if len(sourceedgeTSM)==0: # only for nodes types that exist
            del sourceedgeTSM,nodeTSM
            gc.collect()
            print(EdgeType,"doesn't exist, FilterException=",FilterException)
            print(80*"-")
            if Verbose: print("End",EdgeType," in todolistEdges")
        else:
            print("Start for edges of type :",EdgeType)
            print(80*"-")
                
            PlotTSoverTimeEdgesNodes(nodeTSM,NodeType,sourceedgeTSM,EdgeType,
                        EPOCH=EPOCH,Xmin=int(EPOCH[:4]),
                            yearplain=yearplain,yeardash=yeardash,yearlongdash=yearlongdash,
                            yscale=yscale,figsize=figsize)
            del nodeTSM
            gc.collect()
            if Verbose: print("End PlotTSoverTimeEdgesNodes")
            
            
            if Verbose: print("Start PlotTSoverTimeEdges")
            if i==-1:
                targetedgeTSM=GetTSMFromArrayTS("targetedge",EPOCH,FilterException=FilterException)
            else:
                targetedgeTSM=GetTSMFromArrayTS("targetedge",EPOCH,FilterException=FilterException)[GetArrayFromSandBox("edgetype",FilterException=FilterException)==i ]
                
            PlotTSoverTimeEdges(sourceedgeTSM,targetedgeTSM,EdgeType,
                            EPOCH=EPOCH,Xmin=int(EPOCH[:4]),
                            yearplain=yearplain,yeardash=yeardash,yearlongdash=yearlongdash,
                            yscale=yscale,figsize=figsize)
            del targetedgeTSM,sourceedgeTSM
            gc.collect()
            if Verbose: print("End PlotTSoverTimeEdges")      
                
            if Verbose: print("Start load deltaTS")
            if i==-1:
                deltaTS=GetArrayFromSandBox("edgeDeltaTS",FilterException=FilterException)
            else:
                deltaTS=GetArrayFromSandBox("edgeDeltaTS",FilterException=FilterException)[GetArrayFromSandBox("edgetype",FilterException=FilterException)==i ]
            
            if Verbose: print("Start PlotTSHisto")
            PlotTSHisto(deltaTS,EdgeType,yscale=yscale,figsize=figsize)
            
            if Verbose: print("End PlotTSHisto")      
            
            if Verbose: print("Start DisplayTSstat")
            DisplayTSstat(deltaTS,EdgeType)
            if Verbose: print("End DisplayTSstat")
            del deltaTS
            gc.collect()
            
    
def ReportMain_C(nodes,edges,nodesad,d,Nnodes,Nedges,
                 EPOCH='1970-01-01',
                 yscale="log",xscale="log",
                 figsize=(18,4),YearBegin=1970,YearEnd=2025,YearSlice=2,YearList=None,FAllTypes=False,FilterException=False):    
    

     
    # building/preprocessing
    
    Edgestype=GetArrayFromSandBox("edgetype")
    Edgesencoding=GetArrayFromSandBox("edgeencoding")
    Nodesencoding=GetArrayFromSandBox("nodeencoding")
    Nodestype=GetArrayFromSandBox("nodetype")
    sourceedges=GetSourceEdge(nodes,len(edges))
    sourceType=Nodestype[sourceedges]
    sourceTSY=GetTSMFromArrayTS("sourceedge",EPOCH,)//12
    #sourceTSY=timestampsarray2yearmonth(GetSourceEdgeTimeStamp(nodes,edges,nodesad,d),EPOCH=EPOCH)//12

    if FilterException:
        node_mask=GetTS_Exception("node")
        edge_mask=GetTS_Exception("edge")
        
        Edgestype=Edgestype[edge_mask]
        Nodestype=Nodestype[node_mask]
        sourceedges=sourceedges[edge_mask]
        sourceType=sourceType[edge_mask]
        sourceTSY=sourceTSY[edge_mask]
        edgestmp=edges[edge_mask]
    else:
        edgestmp=edges

    dout={}
    din={}
    
    if FAllTypes:
        GetAllTypesDegreeStats(sourceedges,edgestmp,sourceTSY,dout,din,Nnodes,
                               YearBegin=YearBegin,YearEnd=YearEnd,EPOCH=EPOCH,YearSlice=YearSlice,YearList=YearList) 

    else:
        GetPerTypesDegreeStats(sourceedges,edgestmp,sourceTSY,Edgestype,Edgesencoding, dout,din,Nnodes,
                           PerFlag="EDGE",
                           YearBegin=YearBegin,YearEnd=YearEnd,EPOCH=EPOCH,YearSlice=YearSlice,YearList=YearList) 

       
    Epoch=EPOCH2Epoch(EPOCH)
    nodesstat={}
    
    nodesTSY=GetTSMFromArrayTS("node",EPOCH,FilterException=FilterException)//12
    
    def GetSourceNodeType(EdgeType,Nodesencoding):
        if EdgeType=="All types":
            SourceNodeType=EdgeType
            SourceNodeTypeIndex=-1
        else:
            SourceNodeTypeIndex=np.where(np.array(Nodesencoding)==EdgeType.split(">")[0])[0][0]
            SourceNodeType=Nodesencoding[SourceNodeTypeIndex]
        return SourceNodeTypeIndex,SourceNodeType
    
    for EdgeType in din.keys():
        SourceNodeTypeIndex,SourceNodeType=GetSourceNodeType(EdgeType,Nodesencoding)
        print("Start SourceNodeType",SourceNodeType)
        if SourceNodeType not in nodesstat:
            nodesstat[SourceNodeType]={}
        #YearListTodo=list(din[EdgeType].keys())
        #if YearList!=None:
        #    YearListTodo+=YearList
        #    YearListTodo=list(set(YearListTodo))
        for year in din[EdgeType].keys():
            if year not in nodesstat[SourceNodeType]:
                if year==3000:
                    if SourceNodeType=="All types":
                        mask=np.ones_like(Nodestype,dtype='bool')
                    else:
                        mask=Nodestype==SourceNodeTypeIndex
                else:
                    if SourceNodeType=="All types":
                        mask=nodesTSY<(year-Epoch)
                    else:
                        mask=np.logical_and(Nodestype==SourceNodeTypeIndex,nodesTSY<(year-Epoch)) # < YYY/01/01 (ie up to (YYYY-1)/12/31) included)
                nodesstat[SourceNodeType][year]=np.sum(mask)
    # normalize din/dout wrt nodestat
    for key,value in din.items():
        SourceNodeTypeIndex,SourceNodeType=GetSourceNodeType(key,Nodesencoding)
        #for year in sorted(value.keys()):      
        for year in sorted(value.keys()):      
            if dout[key][year]["ccdf"][0]!=nodesstat[SourceNodeType][year]:
                #print(key,year,'ToBeFixed')
                delta=dout[key][year]["ccdf"][0]-nodesstat[SourceNodeType][year]
                if delta<0 or dout[key][year]["y"][0]<delta:
                    print("ERROR","key",key,"year",year,"delta",delta,dout[key][year]["ccdf"][0],nodesstat[SourceNodeType][year])
                dout[key][year]["y"][0]-=delta # juste le premier
                dout[key][year]["ccdf"][0]-=delta # tous
            if din[key][year]["ccdf"][0]!=nodesstat[SourceNodeType][year]:
                #print(key,year,'ToBeFixed')
                delta=din[key][year]["ccdf"][0]-nodesstat[SourceNodeType][year]
                din[key][year]["y"][0]-=delta # juste le premier
                din[key][year]["ccdf"][0]-=delta # tous


    del node_mask,edge_mask,Edgestype,Nodestype,sourceedges,sourceType,sourceTSY
    if FilterException:
        del edgestmp
    gc.collect()
        
    nrow=4
    ncol=4
    c="k"
    lw=2
    markersize=lw
    
    if YearList==None:
        YearList=np.arange(YearBegin,YearEnd,YearSlice*((YearEnd-YearBegin)//YearSlice//4))
        YearList=YearList[1:]
        print("Default YearList",YearList)
    
    
    def GetNormalize(x,pdf,ccdf):
        cx_ccdf=1
        cx_pdf=1

        if 0 in x:
            cy_ccdf=ccdf[np.where(x==0)[0][0]]
            cy_pdf=cy_ccdf
        else:
            cy_ccdf=ccdf[0]
            cy_pdf=ccdf[0]

        return cx_pdf,cy_pdf,cx_ccdf,cy_ccdf

    def GetRescaledMedian(x,pdf,ccdf):
        #return GetNormalize(x,pdf,ccdf)
        #cy from normalized
        _,cy_pdf,_,cy_ccdf=GetNormalize(x,pdf,ccdf)

        # look for the mediane (excluding 0)
        if ccdf[-1]/cy_ccdf<=0.5:
            # more precise value could be extrapolated
            tmp=np.abs(ccdf/cy_ccdf-0.5)
            cx_ccdf=x[np.argmin(tmp)]
            #cx_ccdf=x[np.where(ccdf/cy_ccdf<=0.5)[0][0]]
        else:
            cx_ccdf=1
        cx_pdf=cx_ccdf  
        return cx_pdf,cy_pdf,cx_ccdf,cy_ccdf

    def GetRescaledMean(x,pdf,ccdf):
        #return GetNormalize(x,pdf,ccdf)
        #cy from normalized
        _,cy_pdf,_,cy_ccdf=GetNormalize(x,pdf,ccdf)

        # look for the average (excluding 0)
        cx_ccdf=np.sum(x*pdf)/np.sum(pdf)
        cx_pdf=cx_ccdf  
        return cx_pdf,cy_pdf,cx_ccdf,cy_ccdf


    for s in din.keys():

        
        #YearList=[1940,1950,1960,1970,1980,1990,2000,2010,2020]
        ColorList=["tab:orange","tab:green","tab:red","tab:blue","k","tab:purple","tab:brown","tab:pink","tab:gray","tab:olive","tab:cyan"]
        ColorList=ColorList[:len(YearList)]

        
        fig, axes = plt.subplots(nrow, ncol, figsize=(18, nrow*4))
        fig.suptitle('Edge type '+s)
        for i in range(nrow):
            for j in range(ncol):
                axes[i,j].grid()
                axes[i,j].set_yscale(yscale)
                axes[i,j].set_xscale(xscale)
        axes[0,0].title.set_text(r'$pdf~(~\delta_{out}(<YYYY)~)$')
        axes[0,1].title.set_text(r'$ccdf~(~\delta_{out}(<YYYY)~)$')
        axes[0,2].title.set_text(r'$pdf~(~\delta_{in}(<YYYY)~)$')
        axes[0,3].title.set_text(r'$ccdf~(~\delta_{in}(<YYYY)~)$')

        axes[0,0].set_ylabel("Cumulated")
        axes[1,0].set_ylabel("Normalized")
        axes[2,0].set_ylabel("Rescaled/Median")
        if nrow==4:
            axes[3,0].set_ylabel("Rescaled/Mean")
        for year in sorted(dout[s].keys()):
            for i in range(nrow):
                if i==0:
                    # cumulative over time
                    cx_pdf_out,cy_pdf_out,cx_ccdf_out,cy_ccdf_out=1,1,1,1
                    cx_pdf_in,cy_pdf_in,cx_ccdf_in,cy_ccdf_in=1,1,1,1
                elif i==1:
                    # normalized
                    cx_pdf_out,cy_pdf_out,cx_ccdf_out,cy_ccdf_out=GetNormalize(dout[s][year]["x"],dout[s][year]["y"],dout[s][year]["ccdf"])
                    cx_pdf_in,cy_pdf_in,cx_ccdf_in,cy_ccdf_in=GetNormalize(din[s][year]["x"],din[s][year]["y"],din[s][year]["ccdf"])
                elif i==2:
                    # rescaled / median
                    cx_pdf_out,cy_pdf_out,cx_ccdf_out,cy_ccdf_out=GetRescaledMedian(dout[s][year]["x"],dout[s][year]["y"],dout[s][year]["ccdf"])
                    cx_pdf_in,cy_pdf_in,cx_ccdf_in,cy_ccdf_in=GetRescaledMedian(din[s][year]["x"],din[s][year]["y"],din[s][year]["ccdf"])
                    #cx_pdf_in,cy_pdf_in,cx_ccdf_in,cy_ccdf_in=cx_pdf_out,cy_pdf_out,cx_ccdf_out,cy_ccdf_out

                elif i==3:
                    # rescaled / mean
                    cx_pdf_out,cy_pdf_out,cx_ccdf_out,cy_ccdf_out=GetRescaledMean(dout[s][year]["x"],dout[s][year]["y"],dout[s][year]["ccdf"])
                    cx_pdf_in,cy_pdf_in,cx_ccdf_in,cy_ccdf_in=GetRescaledMean(din[s][year]["x"],din[s][year]["y"],din[s][year]["ccdf"])
                    #cx_pdf_in,cy_pdf_in,cx_ccdf_in,cy_ccdf_in=cx_pdf_out,cy_pdf_out,cx_ccdf_out,cy_ccdf_out
                mask_out=dout[s][year]["x"]!=0
                mask_in=din[s][year]["x"]!=0
                if year in YearList:
                    c=ColorList[np.where(np.array(YearList)==year)[0][0]]
                    axes[i,0].plot(dout[s][year]["x"][mask_out]/cx_pdf_out,dout[s][year]["y"][mask_out]/cy_pdf_out,'o-',markersize=markersize,c=c,lw=lw,label=str(year),alpha=0.6)
                    axes[i,1].plot(dout[s][year]["x"][mask_out]/cx_ccdf_out,dout[s][year]["ccdf"][mask_out]/cy_ccdf_out,'o-',markersize=markersize,c=c,lw=lw,label=str(year),alpha=0.6)
                    axes[i,2].plot(din[s][year]["x"][mask_in]/cx_pdf_in,din[s][year]["y"][mask_in]/cy_pdf_in,'o-',markersize=markersize,c=c,lw=lw,label=str(year),alpha=0.6)
                    axes[i,3].plot(din[s][year]["x"][mask_in]/cx_ccdf_in,din[s][year]["ccdf"][mask_in]/cy_ccdf_in,'o-',markersize=markersize,c=c,lw=lw,label=str(year),alpha=0.6)
                elif year!=3000: # exclude all time disitribution if not in YearList
                    c='k'
                    axes[i,0].plot(dout[s][year]["x"][mask_out]/cx_pdf_out,dout[s][year]["y"][mask_out]/cy_pdf_out,markersize=markersize,c=c,lw=1,alpha=0.25,zorder=0)
                    axes[i,1].plot(dout[s][year]["x"][mask_out]/cx_ccdf_out,dout[s][year]["ccdf"][mask_out]/cy_ccdf_out,markersize=markersize,c=c,lw=1,alpha=0.25,zorder=0)
                    axes[i,2].plot(din[s][year]["x"][mask_in]/cx_pdf_in,din[s][year]["y"][mask_in]/cy_pdf_in,markersize=markersize,c=c,lw=1,alpha=0.25,zorder=0)
                    axes[i,3].plot(din[s][year]["x"][mask_in]/cx_ccdf_in,din[s][year]["ccdf"][mask_in]/cy_ccdf_in,markersize=markersize,c=c,lw=1,alpha=0.25,zorder=0)
        axes[0,ncol-1].legend()
        plt.tight_layout()
        plt.show()

    del din,dout
    gc.collect()