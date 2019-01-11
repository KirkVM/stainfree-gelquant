#class gel_lane:
#    def __init__(self):
import numpy as np
from copy import copy,deepcopy
from typing import List
import skimage
import networkx as nx
from sklearn.linear_model import LinearRegression
from operator import attrgetter

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

def new_def_lw(lcluster):
    itensor=list(lc0.gbs)[0].rp.inertia_tensor
    ict=np.linalg.norm(it[:,1]
print(np.linalg.norm(ict/np.linalg.norm(ict)))
    
    logistic=1-np.power(1.7,x-9)/(1+np.power(1.7,x-9))



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
            fracband=
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
                    if ol>0.2 and rowdistance<(curlc.row_max-curlc.row_min):
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

        itensor=self.rp.inertia_tensor
        dval=itensor[0,1]/(-itensor[1,1])
        ic_orthog_vector=np.array([1,dval])
        center_point=self.rp.weighted_centroid
        numrows=gelimg.shape[0]
        self.orthogline=[center_point[1]+(x-center_point[0])*(ic_orthog_vector[1]/ic_orthog_vector[0]) \
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


def get_likelihood(lanes):
    rate_term=2.0
    sum_likelihood=0.0
    for lidx,lane in enumerate(lanes):
        for bidx,band in enumerate(lane.bands):
            for olidx,olane in enumerate(lanes):
                inlane_overlaps=0.1
                outlane_overlaps=0.1
                for obidx,oband in enumerate(olane.bands):
                    if olidx==lidx and obidx==bidx:continue
                    obandcoord_tuple=((x,y) for x,y in oband.regionprops.coords)
                    bandcoord_tuple=((x,y) for x,y in band.regionprops.coords)
                    frac_overlap1=len(band.band_extrap_set.intersection(obandcoord_tuple))/len(oband.regionprops.coords)
                    frac_overlap2=len(oband.band_extrap_set.intersection(bandcoord_tuple))/len(band.regionprops.coords)
                    frac_overlap=max(frac_overlap1,frac_overlap2)
                    if lidx==olidx:
                        inlane_overlaps+=frac_overlap
                    else:
                        outlane_overlaps+=frac_overlap
                likelihood_band=1.0-np.exp(-rate_term*inlane_overlaps/(inlane_overlaps+outlane_overlaps))
                sum_likelihood+=likelihood_band
    return sum_likelihood
#            
#        for band in lane.bands:
#