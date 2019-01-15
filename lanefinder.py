#class gel_lane:
#    def __init__(self):
import numpy as np
from copy import copy,deepcopy
from typing import List
import skimage
import networkx as nx
from sklearn.linear_model import LinearRegression
from operator import attrgetter
from scipy import ndimage as ndi

class gel:
    def __init__(self,gelimg:np.ndarray,peaklblimg,geltype:str="protein"):
        '''constructor requires processed gelimage and gelbands determined elsewhere'''
        self.lanes=[]
        self.intensityimage=gelimg
        self.lblpkimage=peaklblimg
        self.all_gelbands=[]
        self.all_regionprops=skimage.measure.regionprops(peaklblimg,intensity_image=gelimg,coordinates='rc')
        for rp in self.all_regionprops:
            self.all_gelbands.append(gel_band(rp,gelimg))
        self.sizegood_gelbands=[]
        self.size2small_gelbands=[]
        self.size2big_gelbands=[]
        self.sizegood_orientgood_gelbands=[]

        if geltype=='protein':
            for gb in self.all_gelbands:
                if gb.rp.inertia_tensor[1,1]>10 and gb.rp.inertia_tensor[1,1]<50:
                    self.sizegood_gelbands.append(gb)
                elif gb.rp.inertia_tensor[1,1]<=10:
                    self.size2small_gelbands.append(gb)
                elif gb.rp.inertia_tensor[1,1]>50:
                    self.size2big_gelbands.append(gb)            
            self.sizegood_orientgood_gelbands=[gb for gb in self.sizegood_gelbands \
                                            if abs(gb.rp.inertia_tensor[0,0]/gb.rp.inertia_tensor[1,1])<0.2 ]

class LaneCluster:
    def __init__(self,gelobj:gel,gbs:List):
        allcoords=[]
        centers=[]
        self.gel=gelobj
        self.gbs=gbs
        self.csize=len(self.gbs)
        self.pband=[]
        self.lcprobimage=None
        for gb in self.gbs:
            centers.append(gb.rp.weighted_centroid)
            for coord in gb.rp.coords:
                allcoords.append(coord)
        allcoords=np.array(allcoords)
        centers=np.array(centers)
    
        self.row_max=max(allcoords[:,0])
        self.row_min=min(allcoords[:,0])
        self.col_max=max(allcoords[:,1])
        self.col_min=min(allcoords[:,1])
        xarray=np.array([np.array([x]) for x in centers[:,0]])

        reg=LinearRegression().fit(X=xarray,y=centers[:,1].T)
        self.centerline_tuples=[(r,int(reg.intercept_+reg.coef_[0]*r)) for r in range(self.gel.intensityimage.shape[0])]

    def score_band_alignment(self,othergb):
#        print(self.centerline_tuples)
        centerline_intersect=set(self.centerline_tuples).intersection(othergb.coord_tuples)
        if len(centerline_intersect)>0:
#            print('hi',len(centerline_intersect))
            centerdist=np.inf
            for oicoord in centerline_intersect:
                curdist=(np.sqrt(
                             (othergb.rp.weighted_centroid[0]-oicoord[0])**2. \
                            +(othergb.rp.weighted_centroid[1]-oicoord[1])**2.))
                centerdist=min(curdist,centerdist)
            othergb_axlength=othergb.rp.major_axis_length
            overlapvalue=(othergb_axlength-centerdist)/othergb_axlength
            return overlapvalue
        else:
            return None

def get_bandspan(gb,prior_width=18,exspan=32,rdim=3,cdim=3):
    gb.bprobspan=[]

    ctr_row=int(gb.rp.weighted_centroid[0])
    ctr_col=int(gb.rp.weighted_centroid[1])
    ctrwindow=gb.intensityimage[ctr_row-int(rdim/2):ctr_row+int(rdim/2)+1,\
                        ctr_col-int(cdim/2):ctr_col+int(cdim/2)+1]
    ctrwindow_intensity=np.mean(ctrwindow)
    osteps=np.array(range(-int(exspan/2),int(exspan/2)+1) )
    prior_logistic=1-np.power(1.7,np.abs(osteps)-prior_width/2)/(1+np.power(1.7,np.abs(osteps)-prior_width/2))
    for spos,ostep in enumerate(osteps):
        deltarc=gb.norm_ict*ostep
        rcpair=np.array(gb.rp.weighted_centroid)+deltarc
        imgwindow=gb.intensityimage[int(rcpair[0])-int(rdim/2):int(rcpair[0])+int(rdim/2)+1,\
                        int(rcpair[1])-int(cdim/2):int(rcpair[1])+int(cdim/2)+1]
        wdw_mean=np.mean(imgwindow)
        pdata_inlane=min(1.0,wdw_mean/ctrwindow_intensity)
        numy=pdata_inlane*(0.05+0.95*prior_logistic[spos])
        gb.bprobspan.append([rcpair,numy])
        

def infer_bandwidths(gb1,gb2):
    region_image=np.zeros_like(gb1.intensityimage)
    prob_image=np.zeros_like(gb1.intensityimage)
    for bprob1,bprob2 in zip(gb1.bprobspan,gb2.bprobspan):
        rc1int=bprob1[0].astype(int)
        rc2int=bprob2[0].astype(int)
        connect_reg=LinearRegression().fit(np.reshape(np.array([rc1int[0],rc2int[0]]),(2,1)),
                                           np.reshape(np.array([rc1int[1],rc2int[1]]),(2,1)))
        for rval in range(rc1int[0]+1,rc2int[0]):
            cval=int(connect_reg.intercept_+connect_reg.coef_[0]*rval)
            region_image[rval,cval]=1
        region_image[rc1int[0],rc1int[1]]=1
        region_image[rc2int[0],rc2int[1]]=1
    #to fill in holes...
    region_image=ndi.binary_fill_holes(region_image)
    for bprob1,bprob2 in zip(gb1.bprobspan,gb2.bprobspan):
        rc1int=bprob1[0].astype(int)
        rc2int=bprob2[0].astype(int)
        region_image[rc1int[0],rc1int[1]]=0
        region_image[rc2int[0],rc2int[1]]=0
    labeled_peaks,_=ndi.label(region_image,structure=[[1,1,1],[1,1,1],[1,1,1]])
    rirp=skimage.measure.regionprops(labeled_peaks)
    #now fill in prop. to distance to each prob_band
    gb1pband_coords=[x[0] for x in gb1.bprobspan]
    gb2pband_coords=[x[0] for x in gb2.bprobspan]
    for rc in rirp[0].coords:
        ds2pband1=[np.sqrt( (rc[0]-x[0])**2. + (rc[1]-x[1])**2.) for x in gb1pband_coords]
        mindist_topband1=min(ds2pband1)
        pband1_contribution=gb1.bprobspan[ds2pband1.index(mindist_topband1)][1]

        ds2pband2=[np.sqrt( (rc[0]-x[0])**2. + (rc[1]-x[1])**2.) for x in gb2pband_coords]
        mindist_topband2=min(ds2pband2)
        pband2_contribution=gb2.bprobspan[ds2pband2.index(mindist_topband2)][1]

        probval=pband1_contribution*(mindist_topband1/(mindist_topband1+mindist_topband2))
        probval+=pband2_contribution*(mindist_topband2/(mindist_topband1+mindist_topband2))
        prob_image[rc[0],rc[1]]=probval
    return prob_image

def new_def_lw(lcluster):
    #order bands by row pos
    lcbands=sorted(list(lcluster.gbs),key=lambda x:x.rp.weighted_centroid[0])#key=attrgetter('rp.weighted_centroid[0]'))
    window_size=3 #3x3 windows
    #get probability axis for band 1
    get_bandspan(lcbands[0])
    #for lcbidx in range(1,len(lcbands)):
    allprobimage=np.zeros_like(lcbands[0].intensityimage)
    for lcbidx in range(1,len(lcbands)):
        get_bandspan(lcbands[lcbidx])
        curprobimage=infer_bandwidths(lcbands[lcbidx-1],lcbands[lcbidx])
        allprobimage+=curprobimage#.max(a,probimage)
    for lcband in lcbands:
        for rcp in lcband.bprobspan:
            rcint=rcp[0].astype(int)
            allprobimage[rcint[0],rcint[1]]=rcp[1]
    return allprobimage

    #for bands 2,...n, get probab axis, then infer band axis btwn


def new_def_lws(lclusters):
    probimage=np.zeros_like(list(lclusters[0].gbs)[0].intensityimage)
    for lcluster in lclusters:
        lcprobimage=new_def_lw(lcluster)
        lcluster.lcprobimage=lcprobimage
        probimage+=lcprobimage
    probimage[probimage>1]=1
    return probimage


import scipy

def lane_model(rvals,c0):#0):
    #return c0#+(rvals-232)*crslope
    return np.array([c0 for x in range(len(rvals))])#)#+(rvals-232)*crslope

def l2_loss(tup, rvals, pbandimg):
    x0 = tup[0]
    predicted_cvals_int = lane_model(rvals,x0).astype(int)
    
    observed_prob_band=np.array([pbandimg[r,c] for r,c in zip(rvals,predicted_cvals_int) ])
#    print(observed_prob_band)
#   print(predicted_cvals_int)
    delta= np.ones((len(observed_prob_band))) - observed_prob_band
    #print(np.dot(delta,delta))
    return np.dot(delta, delta)
#    print(predicted_cvals_int)
#    predicted_rcpairs=[np.array(r,c) for r,c in zip(rvals,predicted_cvals_int)]
#    predicted_prob_band=np.array([1.0 for r,c in zip(rvals,predicted_cvals_int)])
#    print(predicted_prob_band)
 #   print([[rcpair[0],rcpair[1]] for rcpair in zip(rvals,predicted_cvals_int)])

def lane_model_b(rvals,cstart,r0,crslope):#0):
    #return c0#+(rvals-232)*crslope
    predcols=cstart+(rvals-r0)*crslope
    return predcols
#    return np.array([c0 for x in range(len(rvals))])#)#+(rvals-232)*crslope

def l2_loss_b(tup, all_rvals, pbandimg):
    c0,r0,cstep,crslope= tup
    all_deltas=[]
    for lnum,lanerows in enumerate(all_rvals):
        cstart=c0+lnum*cstep
        predicted_cvals_int = lane_model_b(lanerows,cstart,r0,crslope).astype(int)
        observed_prob_band=np.array([pbandimg[r,c] if c<775 else 0 for r,c in zip(lanerows,predicted_cvals_int) ])
        delta= np.ones((len(observed_prob_band))) - observed_prob_band
        all_deltas.extend(delta)
    #print(np.dot(all_deltas,all_deltas))
    #print(np.dot(delta,delta))
    return np.dot(all_deltas, all_deltas)
#

#this functino can probably be deleted! 1/12/19. Not used, use new_def_lws instead
def define_lane_widths(lclusters):
    bandsize=3
    allrow_min=min(x.row_min for x in lclusters)
    allrow_max=max(x.row_max for x in lclusters)
    gimg=lclusters[0].gel.intensityimage
    stuffs=[]
    allbands=[]
    for lc in lclusters[:1]:
        centerseg_tuples=[x for x in lc.centerline_tuples if (x[0]>=allrow_min and x[0]<=allrow_max)]
        centerband=np.zeros((len(centerseg_tuples),bandsize))
        for rcpos,rc in enumerate(centerseg_tuples):
            centerband[rcpos,:]=[   gimg[centerseg_tuples[rcpos][0] ,centerseg_tuples[rcpos][1] + x]     
                                for x in range(-int((bandsize)/2),int(bandsize/2)+1)]
        inlane=True
        offset=0
        while(inlane):
            offset+=1
            shiftband=np.zeros((len(centerseg_tuples),3))
            for rcpos,rc in enumerate(centerseg_tuples):
                shiftband[rcpos,:]=[   gimg[centerseg_tuples[rcpos][0] ,centerseg_tuples[rcpos][1] + x + offset]     
                                     for x in range(-int((bandsize)/2),int(bandsize/2)+1)]
#            fracband=
            U,s,Vh=np.linalg.svd(shiftband)
            stuffs.append([U,s,Vh])
            #Udiff=U-Uc
            #Udiff_abs=np.abs(Udiff)
            diffband=np.array([centerband[x]-shiftband[x] for x in range(len(centerband))])
            avgdiffband=[np.mean(x) for x in diffband]
            allbands.append(avgdiffband)
#            Udiff_abs=np.abs(avgdiffband)
#            U,s,Vh=np.linalg.svd(diffband)
#            stuffs.append([U,s,Vh])
#            print(np.mean(Udiff_abs))
 
        
        allbands.append(np.array([np.mean(x) for x in centerband[:,]]))
        Uc,sc,Vhc=np.linalg.svd(centerband)
        stuffs.append([Uc,sc,Vhc])
        #move to the right!
        inlane=True
        offset=0
        while(inlane):
            offset+=1
            shiftband=np.zeros((len(centerseg_tuples),3))
            for rcpos,rc in enumerate(centerseg_tuples):
                shiftband[rcpos,:]=[   gimg[centerseg_tuples[rcpos][0] ,centerseg_tuples[rcpos][1] + x + offset]     
                                     for x in range(-int((bandsize)/2),int(bandsize/2)+1)]
            U,s,Vh=np.linalg.svd(shiftband)
            stuffs.append([U,s,Vh])
            #Udiff=U-Uc
            #Udiff_abs=np.abs(Udiff)
            diffband=np.array([centerband[x]-shiftband[x] for x in range(len(centerband))])
            avgdiffband=[np.mean(x) for x in diffband]
            allbands.append(avgdiffband)
#            Udiff_abs=np.abs(avgdiffband)
#            U,s,Vh=np.linalg.svd(diffband)
#            stuffs.append([U,s,Vh])
#            print(np.mean(Udiff_abs))
            if offset>50:
                inlane=False
    return stuffs,allbands

class GelDG:
    def __init__(self,gelobj:gel):
        self.gel=gelobj
        self.band_DG=None
        self.build_band_DG()
        self.lane_clusters=[]
        self.unclustered=[]
    def build_band_DG(self):
        self.band_DG=nx.DiGraph()
        for gb in self.gel.sizegood_orientgood_gelbands:
            self.band_DG.add_node(gb)
        print("graph built, adding edges")
        for n1num,node1 in enumerate(list(self.band_DG.nodes)):
            if n1num>0 and n1num%100==0:print(f'added {n1num} nodes...')
            for node2 in list(self.band_DG.nodes):
                if node1==node2:continue
                ol=node1.score_band_alignment(node2)
                ndistance=np.sqrt( (node1.rp.weighted_centroid[0]-node2.rp.weighted_centroid[0])**2.
                           +(node1.rp.weighted_centroid[1]-node2.rp.weighted_centroid[1])**2.)
                if ol is not None:
                    self.band_DG.add_edge(node1,node2,alignscore=ol,distance=ndistance)        
        print(f'final graph contains {self.band_DG.number_of_nodes()} nodes and {self.band_DG.number_of_edges()} edges')

    def make_lane_clusters(self):
        self.lane_clusters=[]
        prevclusterbands=[]
        for band in list(self.band_DG.nodes):
            curbands=[]
            if band not in prevclusterbands:
#            if True not in [band in x for x in self.lane_clusters]:
                curbands=self._findothers(band,set())
            if len(curbands)>1:
                self.lane_clusters.append(LaneCluster(self.gel,curbands))
                prevclusterbands.extend(curbands)
            else:
                self.unclustered.append(band)
        print([len(x.gbs) for x in self.lane_clusters])

    def _findothers(self,gb,lanebands):
        lanebands.add(gb)
        neighbors=set(self.band_DG.successors(gb)).intersection(self.band_DG.predecessors(gb))
        strong_neighbors=[n for n in neighbors if (self.band_DG[gb][n]['alignscore']>0.8 
                                                   and self.band_DG[n][gb]['alignscore']>0.8) 
                                                   and self.band_DG[n][gb]['distance']<50 ]    
        for sn in strong_neighbors:
            if sn not in lanebands:
                lanebands=self._findothers(sn,lanebands)
        return lanebands

    def extend_lane_clusters(self):
        lcs=sorted(self.lane_clusters,key=attrgetter('csize'))        
        for lcpos in range(len(lcs)-1,-1,-1):
            curlc=lcs[lcpos]
            for ucgbidx in range(len(self.unclustered)-1,-1,-1):
                curgb=self.unclustered[ucgbidx]
                ol=curlc.score_band_alignment(curgb)
                rowdistance=min(abs(curlc.row_max-curgb.rp.weighted_centroid[0]),
                                abs(curlc.row_min-curgb.rp.weighted_centroid[0]))
                if ol is not None:
                    if ol>0.2 and rowdistance<0.75*(curlc.row_max-curlc.row_min):
                        curlc.gbs.add(curgb)
                        self.unclustered.pop(ucgbidx)
            curlc=LaneCluster(self.gel,curlc.gbs)
            self.lane_clusters[lcpos]=curlc

    def _trymerge(self):
        lcs=sorted(self.lane_clusters,key=attrgetter('csize'))        
        for lcpos in range(len(lcs)-1,0,-1):
            curlc=lcs[lcpos]
            for olcpos in range(lcpos-1,-1,-1):
                domerge=False
                olc=lcs[olcpos]
                for curgb in olc.gbs:
                    ol=curlc.score_band_alignment(curgb)
                    rowdistance=min(abs(curlc.row_max-curgb.rp.weighted_centroid[0]),
                                abs(curlc.row_min-curgb.rp.weighted_centroid[0]))
                    if ol is not None:
                        if ol>0.2 and rowdistance<(curlc.row_max-curlc.row_min):
                            domerge=True
                if domerge:
                    mergegbs=curlc.gbs.union(olc.gbs)
                    curlc=LaneCluster(self.gel,mergegbs)
                    lcs.pop(olcpos)
                    self.lane_clusters=lcs
                    return True
        return False
 

    def merge_lane_clusters(self):
        didmerge=True
        while(didmerge):
            didmerge=self._trymerge()


class gel_lane:
    def __init__(self):
        self.bands=[]

class gel_band:
    def __init__(self,rp,gelimg):
        self.rp=rp
        self.posterior_laneprob=[]
        self.prior_laneprob=[]
        self.current_lane=None
        #these get added as necessary
        self.orthogline=None
        self.coord_tuples=None
        self.orthogline_tuples=None
        self.intensityimage=gelimg
        self.bprobspan=[]
#        itensor=self.rp.inertia_tensor
        ict_notnorm=self.rp.inertia_tensor[:,1]
        self.norm_ict=ict_notnorm/np.linalg.norm(ict_notnorm)
        self.norm_ict=np.array([-self.norm_ict[0],self.norm_ict[1]])
        dval=ict_notnorm[0]/(-ict_notnorm[1])
        self.ict_orthog_vector=np.array([1,dval])
        self.norm_ict_orthog=self.ict_orthog_vector/np.linalg.norm(self.ict_orthog_vector)

        center_point=self.rp.weighted_centroid
        numrows=self.intensityimage.shape[0]
        self.orthogline=[center_point[1]+(x-center_point[0])*(self.ict_orthog_vector[1]/self.ict_orthog_vector[0]) \
                for x in range(numrows)]    
        self.coord_tuples=[tuple(x) for x in self.rp.coords]
        self.orthogline_tuples=[ (x, int(self.orthogline[x])) for x in range(numrows)]

    def score_band_alignment(self,othergb):
        orthog_intersect=set(self.orthogline_tuples).intersection(othergb.coord_tuples)
        if len(orthog_intersect)>0:
            centerdist=np.inf
            for oicoord in orthog_intersect:
                curdist=(np.sqrt(
                             (othergb.rp.weighted_centroid[0]-oicoord[0])**2. \
                            +(othergb.rp.weighted_centroid[1]-oicoord[1])**2.))
                centerdist=min(curdist,centerdist)
            othergb_axlength=othergb.rp.major_axis_length
            overlapvalue=(othergb_axlength-centerdist)/othergb_axlength
            return overlapvalue
        else:
            return None

