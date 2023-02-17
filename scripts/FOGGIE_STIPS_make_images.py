import yt 
from foggie.utils.foggie_load import *
import numpy as np 
import argparse
from astropy.table import Table, vstack
import glob

import cProfile 

from stips.scene_module import SceneModule
from stips.observation_module import ObservationModule
import matplotlib.pyplot as plt 

def gather_tables(sca_number):

    wildcard = 'SCA'+str(sca_number)+'/master_star_table_sca'+sca_number+'.fits'
    print(wildcard)
    filelist = glob.glob(wildcard)
    print("There are : ", len(filelist), " files in the pile.")
    print(filelist)

    t0 = Table.read(filelist[0])
    print(t0)

    i = 0
    for file in filelist[1:]:
        print(file)
        t_new = Table.read(file)
        print(t_new)
        t0 = vstack([t0, t_new])
        print("Merged in Table : ", i, " name of : ", file)
        i = i+1

    t0.write('SCA'+scanumber+'/sim_stars_merged_SCA'+str(sca_number)+'.fits', overwrite=True)
    print("Saved file : ", 'SCA'+str(scanumber)+'sim_stars_merged_SCA'+str(sca_number)+'.fits')

    return t0


def run_one_sca(sca_number, filter, suffix, nstars):

    wildcard = 'SCA'+str(sca_number) + '/sim_stars_????????_SCA'+str(sca_number)+'.fits'
    print(wildcard)

    filelist = glob.glob(wildcard)
    print("There are : ", len(filelist), " files in the pile.")
    print(filelist)

    nn = np.min([int(nstars), len(filelist)])

    offset_ra=800.
    offset_dec=1000.

    obs = { 
     'instrument': 'WFI',
     'filters': [filter],
     'detectors': 18,
     'suffix':'_'+filter,
     'scastart':1,
     'distortion': False,
     'oversample': 5,
     'pupil_mask': '',
     'background': 'avg',
     'observations_id': 1000+int(sca_number),
     'exptime': 10000,
     'offsets': [{'offset_id': 1, 'offset_centre': False, 'offset_ra': offset_ra + 180.*3600., 'offset_dec': offset_dec, 'offset_pa':0.0}]
    }  # making offset_dec more negative moves the stars UP in the FOV


    obm = ObservationModule(obs, parallel_enable=True, parallel_ncores=2)
    obm.nextObservation()

    i = 0
    for file in filelist[0:nn]:
        print(file)
        if ('conv' not in file):
            output_stellar_catalogues = obm.addCatalogue(file, parallel=True, parallel_enable=True, parallel_ncores=2)
            i = i + 1
            print("We are at file:",i+1," out of ", nn)

    psf_file = obm.addError()
    fits_file, mosaic_file, params = obm.finalize(mosaic=False)




def parse_args():
    parser = argparse.ArgumentParser(description="drives STIPS for FOGGIE") 

    parser.add_argument('--scanumber', dest='scanumber', 
                        help='which SCA do you want to run?')
    parser.set_defaults(scanumber=1)

    parser.add_argument('--filter', dest='filter', 
                        help='which band do you want to run?')
    parser.set_defaults(filter='F129')

    parser.add_argument('--suffix', dest='suffix', 
                        help='add this string to final fits')
    parser.set_defaults(filter='')

    parser.add_argument('--nstars', dest='nstars', 
                        help='how many stars to compute?')
    parser.set_defaults(nstars=1000) 
    
    parser.add_argument('--perftest', dest='perftest', 
                        help='run cProfile?')
    parser.set_defaults(perftest = False) 
     
    args = parser.parse_args()
    return args

if __name__ == "__main__":

    args = parse_args()

    print("We are going to run the STIPS driver for ", args.scanumber, " with ", args.nstars) 

    if (args.perftest): 
        cProfile.runctx('run_one_sca(args.scanumber, args.filter, args.suffix, args.nstars)', globals(), locals()) 
    else: 
        run_one_sca(args.scanumber, args.filter, args.suffix, args.nstars)





