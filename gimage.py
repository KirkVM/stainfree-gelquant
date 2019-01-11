import skimage
from operator import attrgetter
import numpy as np
from scipy import ndimage as ndi
from scipy.optimize import curve_fit
def sharpnegpeak(x,bpoint,slope,offset):
    return offset+ slope*np.abs(x-bpoint)

def get_optrotate(gelimg):
    """returns rotation angle that optimally reduces information content 
    according to SVD U-dot-s
    
    Arguments:
        gelimg
    """
    fltimg=skimage.img_as_float(gelimg)
    gross_angles=[-30,-20,-10,-1,1,10,20,30]
    gross_content=[]
    for ga in gross_angles:
        rimg=skimage.transform.rotate(fltimg,ga,mode='edge')
        U,s,Vh=np.linalg.svd(rimg,full_matrices=False)
        gross_content.append(np.sum(
                [np.linalg.norm(U@np.diag(s)[:,x]) for x in range(len(s))]
                ))
    popt,pcov=curve_fit(sharpnegpeak,gross_angles,gross_content,p0=[0,10,100])
    fg_angles=[popt[0]-4.0+0.75*x for x in range(13) if abs(popt[0]-4.0+0.75*x)>1]
    fg_content=[]
    for fga in fg_angles:
        rimg=skimage.transform.rotate(fltimg,fga,mode='edge')
        U,s,Vh=np.linalg.svd(rimg,full_matrices=False)
        fg_content.append(np.sum(
                [np.linalg.norm(U@np.diag(s)[:,x]) for x in range(len(s))]
                ))
    popt,pcov=curve_fit(sharpnegpeak,fg_angles,fg_content,p0=popt)
    return popt[0]

#what does this do??
def get_optrscutoff(gelimg):
    """returns lower cutoff in a rescaling that reduces U-dot-S less than 1%

    Arguments:
        gelimg
    """
    fltimg=skimage.img_as_float(gelimg)
    U,s,Vh=np.linalg.svd(fltimg,full_matrices=False)
    start_energy=np.sum([np.linalg.norm(U@np.diag(s)[:,x]) for x in range(len(s))])
    lower_co=0
    rel_energy=1
    while(rel_energy>0.99):
        lower_co+=0.01
        stimg = skimage.exposure.rescale_intensity(fltimg, in_range=(lower_co,1.0))
        U,s,Vh=np.linalg.svd(stimg,full_matrices=False)
        new_energy=np.sum([np.linalg.norm(U@np.diag(s)[:,x]) for x in range(len(s))])
        rel_energy=new_energy/start_energy
    return lower_co


def get_gelbands(gelimg):
    el_map=skimage.filters.sobel_h(gelimg)
    el_map=skimage.img_as_float(el_map)/el_map.max()
    el_map=skimage.exposure.rescale_intensity(el_map,in_range=(0,0.1))
    markers=np.zeros_like(gelimg)
    markers[el_map>0.3]=2
    markers[el_map<0.1]=1
    seg=skimage.morphology.watershed(el_map,markers)
    peakseg=np.zeros_like(gelimg)
    peakseg[seg==2]=1
    peakseg=ndi.binary_fill_holes(peakseg)
    labeled_peaks,_=ndi.label(peakseg)#ndi.binary_fill_holes(h))
    labeled_peaks=skimage.morphology.remove_small_objects(labeled_peaks,min_size=10)
    #rps=skimage.measure.regionprops(simple,coordinates='rc')
    #labeled_peaks=skimage.morphology.remove_small_holes(labeled_peaks,min_size=10)
    return labeled_peaks

def peakcheck(rprops,labeled_peaks):
    if rprops.bbox_area<100:
        return False
    return True

def add_wspeaks(plblmtrx,add_plblmtrx,gelimg):
    pk_props=skimage.measure.regionprops(plblmtrx)
    addpk_props=skimage.measure.regionprops(add_plblmtrx)
    sorted_addpkprops=sorted(addpk_props,key=attrgetter('bbox_area'))
    for addpkidx in range(len(sorted_addpkprops)-1,-1,-1):
        candidate_pkprop=sorted_addpkprops[addpkidx]
        candidate_pklblmtrx=np.zeros_like(plblmtrx)
        #merge_pklblmtrx=np.zeros_like(plblmtrx)
        merge_pklblmtrx=plblmtrx.copy()
        for coord in candidate_pkprop.coords:
            candidate_pklblmtrx[coord[0],coord[1]]=1
            merge_pklblmtrx[coord[0],coord[1]]=-1
        addboundaries=skimage.segmentation.find_boundaries(candidate_pklblmtrx,mode='outer')
        prevboundaries=skimage.segmentation.find_boundaries(plblmtrx,mode='outer')
        mergeboundaries=skimage.segmentation.find_boundaries(merge_pklblmtrx,mode='outer')
        print(addboundaries.max(),addboundaries.min())
        edges=(prevboundaries | addboundaries)

#        merge_edges=
        #^mergeboundaries
#        edges=edges^mergeboundaries
#        edges=skimage.img_as_uint(edges)
#        edgelblmtrx=np.zeros_like(plblmt)
        edgelblmtrx[addboundaries is False]=0
        edgelblmtrx[addboundaries is False]=0
        print(edges.max(),edges.min())
        print(edgelblmtrx.shape)
        print(edgelblmtrx.max(),edgelblmtrx.min())
#        merge_pklblmtrx[edges!=0]+=1
        #merge_pklblmtrx-=1
#        merge_pklblmtrx[merge_pklblmtrx%1==0.0]=0
        return edgelblmtrx,candidate_pklblmtrx,plblmtrx


def get_labeled_peaks(gelimg):
    el_map=skimage.filters.sobel(gelimg)
    markers=np.zeros_like(gelimg)
    markers[gelimg<0.1]=1
    markers[gelimg>.1]=2
    markers[gelimg>.2]=3
    markers[gelimg>.3]=4
    markers[gelimg>.4]=5
    markers[gelimg>.5]=6
    markers[gelimg>.6]=7
    markers[gelimg>.7]=8
    markers[gelimg>.8]=9
    markers[gelimg>.9]=10
    seg=skimage.morphology.watershed(el_map,markers)

    prevseg=np.zeros_like(seg)
    curseg=np.zeros_like(seg)
    curlabeled_peaks=np.zeros_like(seg)
    prevlabeled_peaks=np.zeros_like(seg)
    curpk_props=[]
    tempseg=np.zeros_like(seg)
    tempseg[seg==10]=1
    labeled_peaks,_=ndi.label(ndi.binary_fill_holes(tempseg))
#    prevpk_props=skimage.measure.regionprops(labeled_peaks)

    for mrkrval in range(9,8,-1):
        curseg=np.zeros_like(seg)
        curseg[seg==mrkrval]=1
        curseg=ndi.binary_fill_holes(curseg)
        curlabeled_peaks,_=ndi.label(curseg)
        #now sort these in ascending order
        plblmtrx,candidate_pklblmtrx,edges=add_wspeaks(labeled_peaks,curlabeled_peaks,gelimg)



    return plblmtrx,candidate_pklblmtrx,edges
#        for curpk_prop in curpk_props:
#            curbbox=curpk_prop.bbox
#            insides=[]
#            for prevpk_prop in prevpk_props:
#                prevbbox=prevpk_prop.bbox
#                insidepoints=0
#                if curbbox[0]>prevbbox[0]:insidepoints+=1
#                if curbbox[1]>prevbbox[1]:insidepoints+=1
#                if curbbox[2]<prevbbox[2]:insidepoints+=1
#                if curbbox[3]<prevbbox[3]:insidepoints+=1
#                if insidepoints>=4:
#                    insides.append(prevpk_prop.label)
#            if len(insides)>0:
#                if len(insides)>1:
#                    print(f'skipping at {newcoord[0]},{newcoord[1]}')
#                else:
#                    for newcoord in curpk_prop.coords:
#                        labeled_peaks[newcoord[0],newcoord[1]]=insides[0]
        #elif len(connecteds)>0:
        #    print('hi')
