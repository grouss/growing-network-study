import pickle,json
import numpy as np
import matplotlib.pyplot as plt
import warnings

import networkx as nx

from . import config
from .mle import *
from .graph import *

warnings.filterwarnings("ignore")
libname="plot"

# This module provides specific yet reusable plotting functions developed for this study. 
# These functions are currently undocumented, but will be documented in the near future.


# Copyright 2025 Université Paris Cité, France
# Author: [Guillaume Rousseau](https://www.linkedin.com/in/grouss/), Physics Department, Paris, France
# This is supplemental materials and replication package associated with the preprint available on 
# - arXiv (https://arxiv.org/abs/2501.10145)
# - SSRN  (http://ssrn.com/abstract=5191689
# Current version of python scripts and ressources are available on [github's author page](https://github.com/grouss/growing-network-study)
# This work is currently licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)



#pathexportfigure="../paper-growing-network-2024/figures/"

# fontsize according to number of subfigures in the figure
FontSize={}
FontSize[1]=28
FontSize[2]=FontSize[1]+4
FontSize[3]=FontSize[1]+8


def PlotTSoverTimeNodes(nodearrayTSM,Ntype,
                        EPOCH='1970-01-01',Xmin=1970,Xmax=2025,
                        yearplain=[2017,2018,2019,2020,2021],
                        yeardash=[2008,2011,2014],
                        yearlongdash=[2021.167],
                       yscale="log",figsize=(18,4)):
    x=np.arange(1634)/12+EPOCH2Epoch(EPOCH) # switch from month to year scale
    plt.figure(figsize=figsize)
    plt.title("New Nodes of Type "+Ntype+" per month")
    plt.fill_between(x,np.bincount(nodearrayTSM,minlength=1634),color="tab:blue",label=Ntype+r' $Nodes$',alpha=0.5)
    verticalline(yearplain=yearplain,yeardash=yeardash,yearlongdash=yearlongdash)
    plt.yscale(yscale) 
    plt.xlim(Xmin,Xmax)
    plt.legend(loc='upper left')
    plt.grid()
    plt.show()
    
def PlotTSoverTimeEdges(sourcearrayTSM,targetarrayTSM,Etype,
                        EPOCH='1970-01-01',Xmin=1970,Xmax=2025,
                        yearplain=[2017,2018,2019,2020,2021],
                        yeardash=[2008,2011,2014],
                        yearlongdash=[2021.167],
                       yscale="log",figsize=(18,4)):
    x=np.arange(1634)/12+EPOCH2Epoch(EPOCH)
    plt.figure(figsize=figsize)
    plt.title("New Edges of Type "+Etype+" per month, using Source Node or Target Node TimeStamps")
    plt.fill_between(x,np.bincount(sourcearrayTSM,minlength=1634),color="tab:blue",label=Etype+r' $Source$',alpha=0.5)
    plt.fill_between(x,np.bincount(targetarrayTSM,minlength=1634),color="tab:orange",label=Etype+r' $Target$',alpha=0.5)
    verticalline(yearplain=yearplain,yeardash=yeardash,yearlongdash=yearlongdash)
    plt.yscale(yscale) 
    plt.xlim(Xmin,Xmax)
    plt.legend(loc='upper left')
    plt.grid()
    plt.show()
    
def PlotTSoverTimeEdgesNodes(nodeTSM,Ntype,edgeTSM,Etype,
                        EPOCH='1970-01-01',Xmin=1970,Xmax=2025,
                        yearplain=[2017,2018,2019,2020,2021],
                        yeardash=[2008,2011,2014],
                        yearlongdash=[2021.167],
                       yscale="log",figsize=(18,4)):
    x=np.arange(1634)/12+EPOCH2Epoch(EPOCH)
    plt.figure(figsize=figsize)
    plt.title("Average Edges/Nodes  : New Edges of Type "+Etype+" / New Nodes of Type "+Ntype+" (per month)")
    plt.scatter(x,np.bincount(edgeTSM,minlength=1634)/np.bincount(nodeTSM,minlength=1634),
                color="tab:blue",label=Etype+'/'+Ntype,alpha=0.5)
    verticalline(yearplain=yearplain,yeardash=yeardash,yearlongdash=yearlongdash)
    plt.yscale(yscale) 
    plt.xlim(Xmin,Xmax)
    plt.legend(loc='upper left')
    plt.grid()
    plt.show() 
    
def verticalline(alpha=0.5,lw=1,yearplain=[2017,2018,2019,2020,2021],yeardash=[2008,2011,2014],yearlongdash=[2021.167]):
    #print("Plain vertical lines :",yearplain)
    #print("Dash  vertical lines :",yeardash)
    for year in yearplain:
        plt.axvline(x=year,color='k',alpha=alpha,lw=lw)
    for year in yearlongdash:
        plt.axvline(x=year,color='k',linestyle='--',alpha=1,lw=lw)
    for year in yeardash:
        plt.axvline(x=year,color='k',linestyle=':',alpha=1,lw=lw)
        
def DisplayTSstat(deltaTS,Etype):
    #typed_deltaTS=deltaTS[mask_type]
    mn=60;hour=60*60;day=24*hour;week=7*day;year=365*day;month=year//12
    TSmax=2**32-1
    for tsmin,tsmax,field in [
        (-TSmax,-10*year, f'             TS <= -10 Years'),
        (-10*year,-year,  f'-10 Years <  TS <= - 1 Year '),
        (-year,-month,    f'- 1 Year  <  TS <= - 1 Month'),
        (-month,-week,    f'- 1 Month <  TS <= - 1 Week '),
        (-week,-day,      f'- 1 Week  <  TS <= - 1 Day  '),
        (-day,-hour,      f'- 1 Day   <  TS <= - 1 Hour '),
        (-hour,-mn,       f'- 1 Hour  <  TS <= - 1 Mn   '),
        (-mn,-1,          f'- 1 Mn    <  TS <= - 1 S    '),
        (None,None,    None), # skip line
        (-TSmax,-1,     Etype+f'     Total TS <  0 S (negative)'),
        (None,None,    None), # skip line
        (0,0,          Etype+f'     ***** TS == 0 S (null)    '),
        (None,None,    None), # skip line
        (1,TSmax,      Etype+f'     Total TS >  0 S (positive)'),
        (None,None,    None), # skip line
        (1,mn,            f'  1 S     <= TS <    1 Mn   '),
        (mn,hour,         f'  1 Mn    <= TS <    1 Hour '),
        (hour,day,        f'  1 Hour  <= TS <    1 Day  '),
        (day,week,        f'  1 Day   <= TS <    1 Week '),
        (week,month,      f'  1 Week  <= TS <    1 Month'),
        (month,year,      f'  1 Month <= TS <    1 Year '),
        (year,10*year,    f'  1 Year  <= TS <   10 Years'),
        (10*year,TSmax,   f' 10 Years <= TS             '),
    ]:
        stat={}
        if field==None:
            print('.'*80)
        else:
            if tsmin==0 and tsmax==0:
                a=deltaTS==0
                count=np.count_nonzero(a)
                del a
            elif tsmin<0: 
                a=deltaTS>tsmin
                b=deltaTS<=tsmax
                a=np.logical_and(a,b)
                count=np.count_nonzero(a)
                del a,b
            else:
                a=deltaTS>=tsmin
                b=deltaTS<tsmax
                a=np.logical_and(a,b)
                count=np.count_nonzero(a)
                del a,b
            print(field+f' {count:15,} i.e {np.round(100*count/len(deltaTS),2):>6.2f} %')
        
def PlotTSHisto(deltaTS,Etype,EPOCH='1970-01-01',yscale="log",figsize=(18,4)):
    mn=60;hour=60*60;day=24*hour;week=7*day;year=365*day;month=year//12
    #x=np.arange(1634)/12+EPOCH2Epoch(EPOCH)
    #h=np.sign(deltaTS)*np.log10(np.abs(deltaTS))
    #mask=np.where(deltaTS==0)[0]
    #h[np.where(deltaTS==0)[0]]=0
    h=np.where(deltaTS==0,0,np.sign(deltaTS)*np.log10(np.abs(deltaTS)))
    plt.figure(figsize=figsize)
    #plt.title("Edge Type :"+Etype+"\n Vertical black lines : -1Year/-1month/-1Week/-1Days/-1Hour/-1mn/0s/1mn/1Hour/1Day/1Week/1Month/1Year")
    # vertical line 1mn, 1H, 1D,1W,1Y
    xticks=[-np.log10(10*year),-np.log10(year),-np.log10(month),-np.log10(week),-np.log10(day),-np.log10(hour),-np.log10(mn),0,
            np.log10(mn),np.log10(hour),np.log10(day),np.log10(week),np.log10(month),np.log10(year),np.log10(10*year)]
   
    ax = plt.gca() 
    xtick_labels = ['-10 Year','-1 Year','-1 Month','-1 Week','-1 Day','-1 Hour','-1 Mn','0 s','1 Mn','1 Hour','1 Day','1 Week','1 Month','1 Year','10 Year']
    plt.xticks(xticks)    
    ax.set_xticklabels(xtick_labels,rotation=30,ha='right')    
    
    for xv in [mn,day,hour,week,month,year,10*year]:
        plt.axvline(x=np.log10(xv),color='k',alpha=0.5,lw=2)
        plt.axvline(x=-np.log10(xv),color='k',alpha=0.5,lw=2)
   
    h=plt.hist(h,bins=1001,range=[-10,10],density=False,label=r'$\operatorname{sgn}(\Delta TS){\rm Log}_{10}(|\Delta TS|)$',alpha=0.8)
    plt.legend(loc="upper left")
    plt.yscale("log")
    plt.grid()
    plt.show()
    del h
    
    print("End for edge type :",Etype)   


def Plot_Figure_Degree(x,y,yc,Fall=True,fontsize=28,year=2100,filename='ehat_indegree_RVRV_',
                       height=6,WKS="weightedKS",XCmin=0,
                       xrangeCCDF=None,xrangeDKS=None,xrangeAlpha=None,
                       yrangeCCDF=None,yrangeDKS=None,yrangeAlpha=None,
                       yticksCCDF=None,yticksDKS=None,yticksAlpha=None,
                       xtics=[1e0,1e1,1e2,1e3,1e4,1e5,1e6],
                      ):
    yfit,xfit,imin,x_min_shift,e_x_min_shift,D_max_shift=Get_yfit_yc_y_x(x,y,yc,Fall=Fall,WKS=WKS,XCmin=XCmin)
    fig, axs = plt.subplots(3,figsize=(18,2*height),sharex=True)   
    plt.tight_layout()
    axs[0].loglog(x,yc,"-",color="tab:blue",lw=6,markersize=10,alpha=0.8)
    axs[0].loglog(x,y,"o",color="tab:blue",markersize=8,alpha=0.8)
    axs[0].loglog(xfit,yfit,'-',color="k",lw=3,zorder=-10)
    axs[2].semilogx(x_min_shift,e_x_min_shift,"o",markersize=10,color="k")
    axs[2].axhline(y=e_x_min_shift[imin],ls="--",color="k")
    #axs[2].axhline(y=3,ls="-",color="k",alpha=0.5)
    #axs[2].axhline(y=2,ls="-",color="k",alpha=0.5)
    axs[1].loglog(x_min_shift,D_max_shift,"o",markersize=10,color="k")

    

    xtext=1.5e5
    if type(yrangeCCDF)==type(None): 
        ytext=10**(np.max(np.log10(yc))/3)
    else:
        ytext=10**(np.min(yrangeCCDF)+(np.max(np.log10(yrangeCCDF))-np.min(np.log10(yrangeCCDF))/3))
                   
    if year!=2100:        
        axs[0].text(xtext,ytext,r'ccdf$(\delta_{in})$'+'\n'+str(year),fontdict={"fontsize":fontsize+5,'backgroundcolor':'w'})

    else:
        axs[0].text(1.5e5,10**(np.max(np.log10(yc))/3),
                    r'ccdf$(\delta_{in})$',fontdict={"fontsize":fontsize+5,'backgroundcolor':'w'})

    if type(yrangeAlpha)==type(None): 
        ytext=np.min(e_x_min_shift)+(np.max(e_x_min_shift)-np.min(e_x_min_shift))/3
    else:
        ytext=np.min(yrangeAlpha)+(np.max(yrangeAlpha)-np.min(yrangeAlpha))/3
    axs[2].text(xtext,ytext,r'$\hat{\alpha}(x_c)$'+'\n'+r'$\alpha\simeq$'+str(round(e_x_min_shift[imin],1)),fontdict={"fontsize":fontsize+5,'backgroundcolor':'w'})
    
    if type(yrangeDKS)==type(None): 
        ytext=10**(np.min(np.log10(D_max_shift))+(np.max(np.log10(D_max_shift))-np.min(np.log10(D_max_shift)))/3)
    else:
        ytext=10**(np.min(np.log10(yrangeDKS))+(np.max(np.log10(yrangeDKS))-np.min(np.log10(yrangeDKS)))/3)
    
    Dname={
        "KS":r'$D_{KS}(x_c)$',
        "weightedKS":r'$D_{KSw}(x_c)$',
        "Kuiper":r'$D_{KU}(x_c)$'
    }
    try:
        name=Dname[WKS]
    except:
        name="D_{unknown}"
        print("WARNING unknown distance in Plot_Figure_Degree")
    axs[1].text(xtext,ytext,name+'\n'+r'$x_c^{min}=$'+f'{x_min_shift[imin]}',fontdict={"fontsize":fontsize+5,'backgroundcolor':'w'})
    
    
    xrange=[xrangeCCDF,xrangeDKS,xrangeAlpha]
    yrange=[yrangeCCDF,yrangeDKS,yrangeAlpha]
    yticks=[yticksCCDF,yticksDKS,yticksAlpha]
    
    for i in range(3):
        axs[i].set_xticks(xtics)
        axs[i].axvline(x=x_min_shift[imin],ls="--",color="k",markersize=3,lw=3,zorder=-10,alpha=0.9)
        #axs[i].set_xlim(max(0.5,min(x)/2),max(x)*2)
        if type(xrange[i])!=type(None):
            axs[i].set_xlim(xrange[i][0],xrange[i][1])
        if type(yrange[i])!=type(None):
            axs[i].set_ylim(yrange[i][0],yrange[i][1])
        if type(yticks[i])!=type(None):
            axs[i].set_yticks(yticks[i])



    #axs[1].set_ylim(min(0.5*np.min(D_max_shift),0.005),2)
    
    for i in range(3):
        axs[i].grid()
        axs[i].tick_params(axis='both', which='major', labelsize=fontsize)
    #for xpath in ["./",pathexportfigure]:
    #    plt.savefig(xpath+filename+str(year)+'.png', format='png', dpi=300,bbox_inches='tight')

    plt.show()



def DisplayDegreeMetricsPerOrigin(name,xMinMax,yMinMax,FlagGlobalOnly=True):
    inputfile=pickle.load(open(config.exportpath+name+"_degree_dt_20240930.pkl","rb"))
    nfig=3
    fontsize=28
    fig, axs = plt.subplots(nfig,figsize=(18,4*nfig))   
    #fig.title("Origin's URL: "+inputfile["url"])
    #plt.tight_layout()
    for i in range(nfig):
        #axs[i].set_xlim(1990,2022)
        axs[i].grid()
        if xMinMax[i]!=None:
            axs[i].set_xlim(xMinMax[i][0],xMinMax[i][1])
        if yMinMax[i]!=None:
            axs[i].set_ylim(yMinMax[i][0],yMinMax[i][1])
    # degree in 
    maxy=0
    if FlagGlobalOnly:
        DisplayList=[("din",r'$\delta_{in}$',"k")]
    else:
        DisplayList=[("din",r'$\delta_{in}$',"k"),("din_inner",r'$\delta_{in}^{inner}$',"tab:orange")]
        
    for key,label,color in DisplayList:
        x=inputfile[key]["x"]
        y=inputfile[key]["y"]
        ccdf=inputfile[key]["ccdf"]
        axs[0].loglog(x,y,"o",markersize=3,c=color,alpha=0.7)
        axs[0].loglog(x,ccdf,lw=3,label=r'ccdf('+label+')',c=color,alpha=0.7)
        maxy=max(maxy,np.max(y))
    axs[0].legend(loc='best',fontsize=fontsize,ncol=1,markerscale=1)
    axs[0].tick_params(axis='both',which='major',labelsize=fontsize)

    axs[0].text(1,maxy*10,inputfile["url"],fontdict={"fontsize":fontsize+4,'backgroundcolor':'w'})

    if FlagGlobalOnly:
        ytmp=np.bincount(inputfile["adsource"])
        xtmp=1970.+np.arange(len(ytmp))/12
        axs[1].fill_between(xtmp,ytmp,label=r'$t_{s}$',alpha=0.5)
    else:
        ytmp=np.bincount(inputfile["adsource_outer"])
        xtmp=1970.+np.arange(len(ytmp))/12
        axs[1].fill_between(xtmp,ytmp,label=r'$t_{s}^{outer}$',alpha=0.5)

    
    ytmp=np.bincount(inputfile["adtarget"])
    xtmp=1970.+np.arange(len(ytmp))/12
    axs[1].fill_between(xtmp,ytmp,label=r'$t_{t}$',alpha=0.5)
    
    axs[1].legend(loc='best',fontsize=fontsize+4,ncol=1,markerscale=3)
    axs[1].set_yscale("log") 
    #axs[1].set_xlim(xmin,xmax)
    axs[1].tick_params(axis='both',which='major',labelsize=fontsize-2)
    #verticalline(fig,axs[1])

    rmax=max(abs(min(inputfile["dt"])),abs(max(inputfile["dt"])))
    brange=(-rmax,rmax)
    bins=int(2*len(inputfile["dt"])**(1/3)) #rice formulae
    
    if FlagGlobalOnly:
        axs[2].hist(inputfile["dt"],label=r'$\operatorname{sgn}(t_s-t_t){\rm Log}_{10}(|t_s-t_t|)$',range=brange,bins=bins,color="k",log=True,alpha=0.7)
    else:
        axs[2].hist(inputfile["dt_outer"],label=r'$\operatorname{sgn}(t_s^{outer}-t_t){\rm Log}_{10}(|t_s^{outer}-t_t|)$',range=brange,bins=bins,color="tab:blue",log=True,alpha=0.5)
        axs[2].hist(inputfile["dt_inner"],label=r'$\operatorname{sgn}(t_s^{inner}-t_t){\rm Log}_{10}(|t_s^{inner}-t_t|)$',range=brange,bins=bins,color="tab:orange",log=True,alpha=0.5)
        
    axs[2].tick_params(axis='both',which='major',labelsize=fontsize-2)
    axs[2].legend(loc='best',fontsize=fontsize,ncol=1,markerscale=2)
    #verticalline(fig,axs[2])
    #for xpath in ['./',pathexportfigure]:
    #    plt.savefig(xpath+name+'origin_RVRV_stats.png', format='png', dpi=300,bbox_inches='tight')

    plt.show()

def DisplayTSLGraph_Delta_1(statsoutput):
    depth=1
    try:
        statsoutputtmp=statsoutput[depth]
    except:
        return -1
        
    D = nx.DiGraph()
    
    def SizeValueNode(value):
        return 1+value
    
    def SizeValueEdge(value):
        return 200*int(6+np.log10(max(1,value)))
    
    SizeTotNodes=0
    for key, value in statsoutput[depth].items():
        if ">" not in key and '=' not in key:
            SizeTotNodes+=value
    SizeTotEdges=0
    
    for key, value in statsoutput[depth].items():
        if ">" in key: #  or "=" in key:
            SizeTotEdges+=value
    
    selfloops=[]
    regularedges=[]
    
    
                
    for key, value in statsoutput[depth].items():
        if ">" in key:
            s,t=key.split(">")
            if s!=t and value!=0:
                regularedges.append((s,t))
                D.add_edge(s,t,size=1+value,labels=f'{np.round(100*value/SizeTotEdges,depth):,}%')
            #else:
            #    selfloops.append((s,t))
            #    D.add_edge(s,t,size=10+value,labels=f'{np.round(100*value/SizeTotEdges,depth):,}%')
                
    for key, value in statsoutput[depth].items():
        if value!=0:
            if ">" not in key and '=' not in key:
                if key+">"+key in statsoutput[depth].keys() and statsoutput[depth][key+">"+key]!=0:
                    D.add_node(key,size=SizeValueEdge(value),
                               labels=key+f'\n{np.round(100*value/SizeTotNodes,depth):,}%\n[{np.round(100*statsoutput[depth][key+">"+key]/SizeTotEdges,depth):,}%]')
                else:
                    D.add_node(key,size=SizeValueEdge(value),
                               labels=key+f'\n{np.round(100*value/SizeTotNodes,depth):,}%')

   
    plt.figure(figsize=(6, 4))
    scale=10
    pos={
        '111':np.array([0,0])*scale,
        '101':np.array([3,0])*scale,
        '000':np.array([3,2])*scale,
        '001':np.array([3,-2])*scale,
        '010':np.array([-3,2])*scale,
        '011':np.array([-3,-2])*scale,
    }
    #pos=nx.kamada_kawai_layout(D,scale=10,weight="size")
    nx.draw(D
            ,pos=pos
            ,node_size=list(nx.get_node_attributes(D, 'size').values())
            ,with_labels=False
            ,alpha=0.5,
           )
    
    nx.draw_networkx_labels(D, pos, labels=nx.get_node_attributes(D, 'labels'), font_size=10)
    edgelabels=nx.get_edge_attributes(D, 'labels')
   
    #edge_labels=nx.get_edge_attributes(D,'labels')
    #print(edge_labels)
    #for u,v in list(edge_labels):
    #    if u==v:
    #        del edge_labels[(u,v)] # do not add label for self loop
            
    nx.draw_networkx_edge_labels(D, pos, edge_labels=edgelabels,rotate=True)
    #print(pos)
    Xscale=6*1.4
    Yscale=4*1.4
    plt.xlim(-Xscale*scale/2, Xscale*scale/2)
    plt.ylim(-Yscale*scale/2, Yscale*scale/2)
    #plt.axis('equal')
    plt.show()

def DisplayTSLGraph_Delta_1_2(statsoutput):
    if 1 not in statsoutput or 2 not in statsoutput:
        return -1
        
    D = nx.DiGraph()
    Dsub = nx.DiGraph()

    depth=1
    SizeTotNodes=0
    for key, value in statsoutput[depth].items():
        if ">" not in key and '=' not in key:
            SizeTotNodes+=value
    SizeTotEdges=0
    
    for key, value in statsoutput[depth].items():
        if ">" in key: #  or "=" in key:
            SizeTotEdges+=value
    
    selfloops=[]
    regularedges=[]

    def SizeValueNode(value):
        return 1+value
    def SizeValueEdge(value):
        return 200*int(6+np.log10(max(1,value)))
    
                
    for key, value in statsoutput[depth].items():
        if ">" in key and value!=0:
            s,t=key.split(">")
            if s!=t:
                regularedges.append((s+'(1)',t+'(1)'))
                D.add_edge(s+'(1)',t+'(1)',size=SizeValueNode(value),labels=f'{np.round(100*value/SizeTotEdges,1):,}%')
            #else:
            #    selfloops.append((s,t))
            #    D.add_edge(s,t,size=10+value,labels=f'{np.round(100*value/SizeTotEdges,depth):,}%')
                
    for key, value in statsoutput[depth].items():
        if value!=0:
            if ">" not in key and '=' not in key:
                if key+">"+key in statsoutput[depth].keys() and statsoutput[depth][key+">"+key]!=0:
                    D.add_node(key+'(1)',size=SizeValueEdge(value),
                               labels=key+f'({depth})'+f'\n{np.round(100*value/SizeTotNodes,1):,}%\n[{np.round(100*statsoutput[depth][key+">"+key]/SizeTotEdges,1):,}%]')
                else:
                    D.add_node(key+'(1)',size=SizeValueEdge(value),
                               labels=key+f'({depth})'+f'\n{np.round(100*value/SizeTotNodes,1):,}%')

    


    t="010";
    depth=2
    for s in ["010","020"]:
        if statsoutput[2][s]!=0:
            D.add_node(s+"(2)",size=SizeValueEdge(statsoutput[2][s])
                           ,labels=s+f'(2)\n{np.round(100*statsoutput[2][s]/SizeTotNodes,1):,}%'
                  )    
        value=0
        for stmp in ["111","211","121","221","101","201"]:
                value+=statsoutput[2][s+">"+stmp] 
        if value!=0:
            D.add_edge(s+"(2)",t+"(1)",size=SizeValueNode(value),labels=f'{np.round(100*value/SizeTotEdges,1):,}%')

    t="011";
    depth=2
    for s in ["011","021"]:
        if statsoutput[2][s]!=0:
            D.add_node(s+"(2)",size=SizeValueEdge(statsoutput[2][s])
                           ,labels=s+f'(2)\n{np.round(100*statsoutput[2][s]/SizeTotNodes,1):,}%'
                  )    
        value=0
        for stmp in ["111","211","121","221","101","201"]:
            value+=statsoutput[2][s+">"+stmp]    
        if value!=0:
            D.add_edge(s+"(2)",t+"(1)",size=SizeValueNode(value),labels=f'{np.round(100*value/SizeTotEdges,1):,}%')

    s="101";
    depth=2
    for t in ["101","201"]:
        if statsoutput[2][t]!=0:
            D.add_node(t+"(2)",size=SizeValueEdge(statsoutput[2][t])
                       ,labels=s+f'(2)\n{np.round(100*statsoutput[2][s]/SizeTotNodes,1):,}%'
              )    
        value=0
        for ttmp in ["111","211","121","221","011","021","010","020"]:
            value+=statsoutput[2][ttmp+">"+t]    
        if value!=0:
            D.add_edge(s+"(1)",t+"(2)",size=SizeValueNode(value),labels=f'{np.round(100*value/SizeTotEdges,1):,}%')
    
    for s in ["111","211","121","221"]:
        if s+">"+t not in statsoutput[2]:
            Dsub.add_node(s+"(2)",size=SizeValueEdge(statsoutput[2][s])
                       ,labels=s+f'(2)\n{np.round(100*statsoutput[2][s]/SizeTotNodes,2):,}%'
              )
        else:
            Dsub.add_node(s+"(2)",size=SizeValueEdge(statsoutput[2][s])
                       ,labels=s+f'(2)\n{np.round(100*statsoutput[2][s]/SizeTotNodes,2):,}%\n[{np.round(100*statsoutput[2][s+">"+s]/SizeTotEdges,2):,}%]'
              )
            
        for t in ["111","211","121","221"]:
            if s!=t:
                    value=statsoutput[2][s+">"+t]
                    Dsub.add_edge(s+'(2)',t+'(2)',size=SizeValueNode(value),labels=f'{np.round(100*value/SizeTotEdges,2):,}%')
                
    scale=10

    pos={
        '111(1)':np.array([0,0])*scale,
        '101(1)':np.array([3,0])*scale,
        '000(1)':np.array([3,2])*scale,
        '001(1)':np.array([3,-2])*scale,
        '010(1)':np.array([-3,2])*scale,
        '011(1)':np.array([-3,-2])*scale,


        '010(2)':np.array([-6,2+1])*scale,
        '020(2)':np.array([-6,2-1])*scale,
        
        '011(2)':np.array([-6,-2+1])*scale,
        '021(2)':np.array([-6,-2-1])*scale,

        '101(2)':np.array([6,0+1])*scale,
        '201(2)':np.array([6,0-1])*scale,
    }
    pos2={

        '111(2)':np.array([-1.5,-1.5])*scale,
        '211(2)':np.array([1.5,-1.5])*scale,
        '121(2)':np.array([-1.5,+1.5])*scale,
        '221(2)':np.array([1.5,+1.5])*scale,
    }
    #pos=nx.kamada_kawai_layout(D,scale=10,weight="size")
    plt.figure(figsize=(10, 5))
    nx.draw(D
            ,pos=pos
            ,node_size=list(nx.get_node_attributes(D, 'size').values())
            ,with_labels=False
            ,alpha=0.5,
           )
    
    nx.draw_networkx_labels(D, pos, labels=nx.get_node_attributes(D, 'labels'), font_size=10)
    edgelabels=nx.get_edge_attributes(D, 'labels')
    nx.draw_networkx_edge_labels(D, pos, edge_labels=edgelabels,rotate=True)
    Xscale=16
    Yscale=8
    plt.xlim(-Xscale/2*scale, Xscale/2*scale)
    plt.ylim(-Yscale*scale/2, Yscale*scale/2)
    plt.show()
    
    #print(pos)
    if True:
        plt.figure(figsize=(4, 4))

        nx.draw(Dsub
                ,pos=pos2
                ,node_size=list(nx.get_node_attributes(Dsub, 'size').values())
                ,with_labels=False
                ,alpha=0.5,
               )
        nx.draw_networkx_labels(Dsub, pos2, labels=nx.get_node_attributes(Dsub, 'labels'), font_size=8)
        
        edgelabels2=nx.get_edge_attributes(Dsub, 'labels')
        nx.draw_networkx_edge_labels(Dsub, pos2, edge_labels=edgelabels2,label_pos=0.3)
        Xscale=6
        Yscale=6
        plt.xlim(-Xscale/2*scale, Xscale/2*scale)
        plt.ylim(-Yscale*scale/2, Yscale*scale/2)
        plt.show()
    # subgraph 111(1) <=. 111(2), 211(2),121(2), 221(2)
    return 

def DisplayTSLGraph_Delta_1(statsoutput):
    depth=1
    try:
        statsoutputtmp=statsoutput[depth]
    except:
        return -1
        
    D = nx.DiGraph()
    def SizeValueNode(value):
        return 1+value
    def SizeValueEdge(value):
        return 200*int(6+np.log10(max(1,value)))
    SizeTotNodes=0
    for key, value in statsoutput[depth].items():
        if ">" not in key and '=' not in key:
            SizeTotNodes+=value
    SizeTotEdges=0
    
    for key, value in statsoutput[depth].items():
        if ">" in key: #  or "=" in key:
            SizeTotEdges+=value
    
    selfloops=[]
    regularedges=[]
    
    
                
    for key, value in statsoutput[depth].items():
        if ">" in key:
            s,t=key.split(">")
            if s!=t and value!=0:
                regularedges.append((s,t))
                D.add_edge(s,t,size=1+value,labels=f'{np.round(100*value/SizeTotEdges,depth):,}%')
            #else:
            #    selfloops.append((s,t))
            #    D.add_edge(s,t,size=10+value,labels=f'{np.round(100*value/SizeTotEdges,depth):,}%')
                
    for key, value in statsoutput[depth].items():
        if value!=0:
            if ">" not in key and '=' not in key:
                if key+">"+key in statsoutput[depth].keys() and statsoutput[depth][key+">"+key]!=0:
                    D.add_node(key,size=SizeValueEdge(value),
                               labels=key+f'\n{np.round(100*value/SizeTotNodes,depth):,}%\n[{np.round(100*statsoutput[depth][key+">"+key]/SizeTotEdges,depth):,}%]')
                else:
                    D.add_node(key,size=SizeValueEdge(value),
                               labels=key+f'\n{np.round(100*value/SizeTotNodes,depth):,}%')

   
    plt.figure(figsize=(6, 4))
    scale=10
    pos={
        '111':np.array([0,0])*scale,
        '101':np.array([3,0])*scale,
        '000':np.array([3,2])*scale,
        '001':np.array([3,-2])*scale,
        '010':np.array([-3,2])*scale,
        '011':np.array([-3,-2])*scale,
    }
    #pos=nx.kamada_kawai_layout(D,scale=10,weight="size")
    nx.draw(D
            ,pos=pos
            ,node_size=list(nx.get_node_attributes(D, 'size').values())
            ,with_labels=False
            ,alpha=0.5,
           )
    
    nx.draw_networkx_labels(D, pos, labels=nx.get_node_attributes(D, 'labels'), font_size=10)
    edgelabels=nx.get_edge_attributes(D, 'labels')
   
    #edge_labels=nx.get_edge_attributes(D,'labels')
    #print(edge_labels)
    #for u,v in list(edge_labels):
    #    if u==v:
    #        del edge_labels[(u,v)] # do not add label for self loop
            
    nx.draw_networkx_edge_labels(D, pos, edge_labels=edgelabels,rotate=True)
    #print(pos)
    Xscale=6*1.4
    Yscale=4*1.4
    plt.xlim(-Xscale*scale/2, Xscale*scale/2)
    plt.ylim(-Yscale*scale/2, Yscale*scale/2)
    #plt.axis('equal')
    plt.show()

def DisplayTSLGraph_Delta_1_2(statsoutput):
    if 1 not in statsoutput or 2 not in statsoutput:
        return -1
        
    D = nx.DiGraph()
    Dsub = nx.DiGraph()

    depth=1
    SizeTotNodes=0
    for key, value in statsoutput[depth].items():
        if ">" not in key and '=' not in key:
            SizeTotNodes+=value
    SizeTotEdges=0
    
    for key, value in statsoutput[depth].items():
        if ">" in key: #  or "=" in key:
            SizeTotEdges+=value
    
    selfloops=[]
    regularedges=[]

    def SizeValueNode(value):
        return 1+value
    def SizeValueEdge(value):
        return 200*int(6+np.log10(max(1,value)))
    
                
    for key, value in statsoutput[depth].items():
        if ">" in key and value!=0:
            s,t=key.split(">")
            if s!=t:
                regularedges.append((s+'(1)',t+'(1)'))
                D.add_edge(s+'(1)',t+'(1)',size=SizeValueNode(value),labels=f'{np.round(100*value/SizeTotEdges,1):,}%')
            #else:
            #    selfloops.append((s,t))
            #    D.add_edge(s,t,size=10+value,labels=f'{np.round(100*value/SizeTotEdges,depth):,}%')
                
    for key, value in statsoutput[depth].items():
        if value!=0:
            if ">" not in key and '=' not in key:
                if key+">"+key in statsoutput[depth].keys() and statsoutput[depth][key+">"+key]!=0:
                    D.add_node(key+'(1)',size=SizeValueEdge(value),
                               labels=key+f'({depth})'+f'\n{np.round(100*value/SizeTotNodes,1):,}%\n[{np.round(100*statsoutput[depth][key+">"+key]/SizeTotEdges,1):,}%]')
                else:
                    D.add_node(key+'(1)',size=SizeValueEdge(value),
                               labels=key+f'({depth})'+f'\n{np.round(100*value/SizeTotNodes,1):,}%')

    


    t="010";
    depth=2
    for s in ["010","020"]:
        if statsoutput[2][s]!=0:
            D.add_node(s+"(2)",size=SizeValueEdge(statsoutput[2][s])
                           ,labels=s+f'(2)\n{np.round(100*statsoutput[2][s]/SizeTotNodes,1):,}%'
                  )    
        value=0
        for stmp in ["111","211","121","221","101","201"]:
                value+=statsoutput[2][s+">"+stmp] 
        if value!=0:
            D.add_edge(s+"(2)",t+"(1)",size=SizeValueNode(value),labels=f'{np.round(100*value/SizeTotEdges,1):,}%')

    t="011";
    depth=2
    for s in ["011","021"]:
        if statsoutput[2][s]!=0:
            D.add_node(s+"(2)",size=SizeValueEdge(statsoutput[2][s])
                           ,labels=s+f'(2)\n{np.round(100*statsoutput[2][s]/SizeTotNodes,1):,}%'
                  )    
        value=0
        for stmp in ["111","211","121","221","101","201"]:
            value+=statsoutput[2][s+">"+stmp]    
        if value!=0:
            D.add_edge(s+"(2)",t+"(1)",size=SizeValueNode(value),labels=f'{np.round(100*value/SizeTotEdges,1):,}%')

    s="101";
    depth=2
    for t in ["101","201"]:
        if statsoutput[2][t]!=0:
            D.add_node(t+"(2)",size=SizeValueEdge(statsoutput[2][t])
                       ,labels=s+f'(2)\n{np.round(100*statsoutput[2][s]/SizeTotNodes,1):,}%'
              )    
        value=0
        for ttmp in ["111","211","121","221","011","021","010","020"]:
            value+=statsoutput[2][ttmp+">"+t]    
        if value!=0:
            D.add_edge(s+"(1)",t+"(2)",size=SizeValueNode(value),labels=f'{np.round(100*value/SizeTotEdges,1):,}%')
    
    for s in ["111","211","121","221"]:
        if s+">"+t not in statsoutput[2]:
            Dsub.add_node(s+"(2)",size=SizeValueEdge(statsoutput[2][s])
                       ,labels=s+f'(2)\n{np.round(100*statsoutput[2][s]/SizeTotNodes,2):,}%'
              )
        else:
            Dsub.add_node(s+"(2)",size=SizeValueEdge(statsoutput[2][s])
                       ,labels=s+f'(2)\n{np.round(100*statsoutput[2][s]/SizeTotNodes,2):,}%\n[{np.round(100*statsoutput[2][s+">"+s]/SizeTotEdges,2):,}%]'
              )
            
        for t in ["111","211","121","221"]:
            if s!=t:
                    value=statsoutput[2][s+">"+t]
                    Dsub.add_edge(s+'(2)',t+'(2)',size=SizeValueNode(value),labels=f'{np.round(100*value/SizeTotEdges,2):,}%')
                
    scale=10

    pos={
        '111(1)':np.array([0,0])*scale,
        '101(1)':np.array([3,0])*scale,
        '000(1)':np.array([3,2])*scale,
        '001(1)':np.array([3,-2])*scale,
        '010(1)':np.array([-3,2])*scale,
        '011(1)':np.array([-3,-2])*scale,


        '010(2)':np.array([-6,2+1])*scale,
        '020(2)':np.array([-6,2-1])*scale,
        
        '011(2)':np.array([-6,-2+1])*scale,
        '021(2)':np.array([-6,-2-1])*scale,

        '101(2)':np.array([6,0+1])*scale,
        '201(2)':np.array([6,0-1])*scale,
    }
    pos2={

        '111(2)':np.array([-1.5,-1.5])*scale,
        '211(2)':np.array([1.5,-1.5])*scale,
        '121(2)':np.array([-1.5,+1.5])*scale,
        '221(2)':np.array([1.5,+1.5])*scale,
    }
    #pos=nx.kamada_kawai_layout(D,scale=10,weight="size")
    plt.figure(figsize=(10, 5))
    nx.draw(D
            ,pos=pos
            ,node_size=list(nx.get_node_attributes(D, 'size').values())
            ,with_labels=False
            ,alpha=0.5,
           )
    
    nx.draw_networkx_labels(D, pos, labels=nx.get_node_attributes(D, 'labels'), font_size=10)
    edgelabels=nx.get_edge_attributes(D, 'labels')
    nx.draw_networkx_edge_labels(D, pos, edge_labels=edgelabels,rotate=True)
    Xscale=16
    Yscale=8
    plt.xlim(-Xscale/2*scale, Xscale/2*scale)
    plt.ylim(-Yscale*scale/2, Yscale*scale/2)
    plt.show()
    
    #print(pos)
    if True:
        plt.figure(figsize=(4, 4))

        nx.draw(Dsub
                ,pos=pos2
                ,node_size=list(nx.get_node_attributes(Dsub, 'size').values())
                ,with_labels=False
                ,alpha=0.5,
               )
        nx.draw_networkx_labels(Dsub, pos2, labels=nx.get_node_attributes(Dsub, 'labels'), font_size=8)
        
        edgelabels2=nx.get_edge_attributes(Dsub, 'labels')
        nx.draw_networkx_edge_labels(Dsub, pos2, edge_labels=edgelabels2,label_pos=0.3)
        Xscale=6
        Yscale=6
        plt.xlim(-Xscale/2*scale, Xscale/2*scale)
        plt.ylim(-Yscale*scale/2, Yscale*scale/2)
        plt.show()
    # subgraph 111(1) <=. 111(2), 211(2),121(2), 221(2)
    return 

