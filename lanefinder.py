#class gel_lane:
#    def __init__(self):
import numpy as np
from copy import copy,deepcopy

class gel:
    def __init__(self,gelimg):
        self.lanes=[]
    def assign_lane(self,gb):
        likelihoods=[]
        for lidx,lane in enumerate(self.lanes):
            lanes_copy=self.lanes[:]#deepcopy(self.lanes)
#            print(lanes_copy)
#            print(len(lanes_copy))
            lanes_copy[lidx].bands.append(gb)
            curlikelihood=get_likelihood(lanes_copy)
            likelihoods.append(curlikelihood)
        if len(likelihoods)==0:
            new_lane=gel_lane()
            new_lane.bands.append(gb)
            self.lanes.append(new_lane)
        else:
        #now for case where band is added to new lane:
            lanes_copy=self.lanes[:]
            new_lane=gel_lane()
            new_lane.bands.append(gb)
            lanes_copy.append(new_lane)
            best_current_likelihood=max(likelihoods)
            best_current_lane=likelihoods.index(best_current_likelihood)
            newlane_likelihood=get_likelihood(lanes_copy)
            print(likelihoods,best_current_likelihood,best_current_lane)
            print(newlane_likelihood)
            if newlane_likelihood>best_current_likelihood:
                self.lanes=lanes_copy
            else:
                self.lanes[best_current_lane].bands.append(gb)
        print(len(self.lanes))
class gel_lane:
    def __init__(self):
        self.bands=[]

class gel_band:
    def __init__(self,rp,gimg):
        self.regionprops=rp
        self.posterior_laneprob=[]
        self.prior_laneprob=[]
        self.current_lane=None
        itensor=self.regionprops.inertia_tensor
        dval=itensor[0,1]/(-itensor[1,1])
        ic_orthog_vector=np.array([1,dval])
        center_point=self.regionprops.weighted_centroid
        numrows=gimg.shape[0]
        numcols=gimg.shape[1]
        self.orthogline=[center_point[1]+(x-center_point[0])*(ic_orthog_vector[1]/ic_orthog_vector[0]) \
                for x in range(numrows)]    

        self.coord_tuples=[tuple(x) for x in self.regionprops.coords]
        self.orthogline_tuples=[ (x, int(self.orthogline[x])) for x in range(numrows)]

        #band_extrap_set=set()
        #for coord in self.regionprops.coords:
        #    rcvals=[(x,int(coord[1]+(x-coord[0])*(ic_orthog_vector[1]/ic_orthog_vector[0]))) \
        #        for x in range(numrows)]
        #    for rcval in rcvals:
        #        if rcval[1]>0 and rcval[1]<numcols:
        #            band_extrap_set.add(rcval)
        #self.band_extrap_set=band_extrap_set
    def score_band_connection(self,othergb):
        orthog_intersect=set(self.orthogline_tuples).intersection(othergb.coord_tuples)
        #print(orthog_intersect)
        if len(orthog_intersect)>0:
            centerdist=np.inf
            for oicoord in orthog_intersect:
                curdist=(np.sqrt(
                             (othergb.regionprops.weighted_centroid[0]-oicoord[0])**2. \
                            +(othergb.regionprops.weighted_centroid[1]-oicoord[1])**2.))
                centerdist=min(curdist,centerdist)
#            othergb_axlength=max(othergb.regionprops.major_axis_length,centerdist+0.1)
            othergb_axlength=othergb.regionprops.major_axis_length
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