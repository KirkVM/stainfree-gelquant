import yaml,argparse,os
import pandas as pd



def readin():
    with open("fileconfigs.yaml","r") as f:
        xlfpathinfo=yaml.load(f)
    lbfpath=os.path.join(xlfpathinfo['laneband_basepath'],
            xlfpathinfo['laneband_fldr'],xlfpathinfo['laneband_fname'])
    linfofpath=os.path.join(xlfpathinfo['lanedescription_basepath'],
            xlfpathinfo['lanedescription_fldr'],xlfpathinfo['lanedescription_fname'])
    lbdf=pd.read_excel(lbfpath,skiprows=[0],usecols=[0,1,2,3,4,5,6,7,8,9,10])
    return lbdf
    #pd.read_excel()
