
# ## How to combine multiple STIPS outputs into one observation and create the final image 

import matplotlib.pyplot as plt 
from stips.scene_module import SceneModule
import numpy as np 
from stips.observation_module import ObservationModule
import glob
import copy 
import pickle

with open('obm_1001_F129.save', 'rb') as file: 
    obm1 = pickle.load(file)
    
obm = copy.deepcopy(obm1)
obm.__dict__['imgbase'] =  '/Users/tumlinson/Dropbox/FOGGIE/collab/Roman/roman_images/stream_2457_51/elvis'

obs_list = [] # empty list to contain all the unpickled observations 
for i in 1000+np.arange(18)+1: 
    with open('obm_'+str(i)+'_F129.save', 'rb') as file: 
        obs_list.append(pickle.load(file))    


for i in np.arange(18): 
    #for each SCA, add up all the observations by lopering over the obs_list
    print("Combining SCA ", i) 
    for j in np.arange(18): 
        obm.__dict__['instrument'].__dict__['detectors'][i].data = \
            obm.__dict__['instrument'].__dict__['detectors'][i].data + \
            obs_list[i].__dict__['instrument'].__dict__['detectors'][j].data

obm.finalize(mosaic=False)

