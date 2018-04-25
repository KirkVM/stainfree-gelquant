import yaml,argparse,os
import pandas as pd
from operator import itemgetter


def clean_infodf(infodf):
    #needs 5 column names
    misnamedcol=False
    #reqcols_=['lane','sample','volume','mw','wcount']
    reqcols_=['sample','lane','volume','mw','wcount']
    reqcolorder_=[0,1,2,3,4]
    renameHT={}
    for x,cname in enumerate(list(infodf.columns.values)):
        if cname.lower().strip() in reqcols_:
            renameHT[cname]=reqcols_[reqcols_.index(cname.lower().strip())]
        else:
            print('error cant identify column labeled {}'.format(cname))
    newinfodf=infodf.rename(renameHT,axis='columns')
    return newinfodf

def read_gelbandfile(lbfpath,linfodf):
    lbdf_all=pd.read_excel(lbfpath,skiprows=[0])
    assert(lbdf_all.shape[1]%11 == 0)
    colnames_=['Lane','Band No.','Band Label','Mol. Wt. (KDa)',
               'Relative Front','Adj. Volume (Int)','Volume (Int)',
               'Abs. Quant.','Rel. Quant.','Band %','Lane %']
    lbdf=lbdf_all[colnames_]
    lbdf=lbdf.dropna(how='all')
    lbdfHT={}
    lbdfHT[lbdf['Lane'].iat[0]]=lbdf
    for rnumidx in range(int(lbdf_all.shape[1]/11))[1:]:
        rcolnames_=[x+'.{}'.format(rnumidx) for x in colnames_]
        renameHT={x+'.{}'.format(rnumidx):x for x in colnames_}
        lbdf=lbdf_all[rcolnames_]
        lbdf=lbdf.dropna(how='all')
        if lbdf.shape[0]==0:
            continue
        lbdf=lbdf.rename(renameHT,axis='columns')
        lbdfHT[lbdf['Lane'].iat[0]]=lbdf
    return lbdfHT

def readin():
    with open("fileconfigs.yml","r") as f:
        xlfpathinfo=yaml.load(f)
    lbfpath=os.path.join(xlfpathinfo['laneband_basepath'],
            xlfpathinfo['laneband_fldr'],xlfpathinfo['laneband_fname'])
    linfofpath=os.path.join(xlfpathinfo['lanedescription_basepath'],
            xlfpathinfo['lanedescription_fldr'],xlfpathinfo['lanedescription_fname'])
    linfodf_raw=pd.read_excel(linfofpath)
    linfodf=clean_infodf(linfodf_raw)
    lbdfHT=read_gelbandfile(lbfpath,linfodf)
    return lbdfHT,linfodf


#    lbdf=pd.read_excel(lbfpath,skiprows=[0],usecols=[0,1,2,3,4,5,6,7,8,9,10])
#    return linfodf
    #pd.read_excel()
