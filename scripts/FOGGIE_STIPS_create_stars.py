import yt
from foggie.utils.foggie_load import *
import numpy as np
from astropy.table import Table, vstack 
import argparse
import glob
import copy 
import os 

from stips.scene_module import SceneModule
from stips.observation_module import ObservationModule
import matplotlib.pyplot as plt


def rot_star_table(stars, theta_in_degrees):

    rotated_stars = copy.deepcopy(stars)
    
    theta = (theta_in_degrees / 180 * np.pi)
    new_xx = stars['xx'] * np.cos(theta) + stars['yy'] * np.sin(theta)
    new_yy = -1.*stars['xx'] * np.sin(theta) + stars['yy'] * np.cos(theta)
    rotated_stars['xx'] = new_xx
    rotated_stars['yy'] = new_yy

    return rotated_stars


def prescreen_sca(stars, scanumber, obm): 

    screened_stars = copy.deepcopy(stars)
    
    ra0_sca0 = (obm.__dict__['instrument'].__dict__['detectors'][scanumber-1].__dict__['header']['RA_APER']) % 360.
    dec0_sca0 = obm.__dict__['instrument'].__dict__['detectors'][scanumber-1].__dict__['header']['DEC_APER']
    
    padding = 10. # how far outside each SCA to pick up particles - in arcsec 

    screened_stars = screened_stars[screened_stars['xx'] > (ra0_sca0 - (2044 * 0.11 + padding) / 3600.)]     # left/west edge 
    screened_stars = screened_stars[screened_stars['xx'] < (ra0_sca0 + (2044 * 0.11 + padding) / 3600.)]     # right/east edge 
    screened_stars = screened_stars[screened_stars['yy'] > (dec0_sca0 - (2044 * 0.11 + padding) / 3600.)]    # bottom/south edge
    screened_stars = screened_stars[screened_stars['yy'] < (dec0_sca0 + (2044 * 0.11 + padding) / 3600.)]    # top/north edge 
    
    return screened_stars 

def create_catalogs(streamfile, scanumber=1, rotangle=0., offset_ra=0., offset_dec=0., offset_pa=0., exptime=100000): 

    obs = {
     'instrument': 'WFI',
     'filters': ['F129'],
     'detectors': 18,
     'suffix':'',
     'scastart':1,
     'distortion': False,
     'oversample': 5,
     'pupil_mask': '',
     'background': 'avg',
     'observations_id': 1,
     'exptime': exptime,
     'offsets': [{'offset_id': 1, 'offset_centre': False, 'offset_ra': offset_ra + 180.*3600., 'offset_dec': offset_dec, 'offset_pa':0.0}]    
    }  # making offset_dec more negative moves the stars UP in the FOV

    obm = ObservationModule(obs)
    obm.nextObservation() 
    
    stars = Table.read(streamfile) 
    stars = rot_star_table(stars, rotangle)
    stars['xx'] = stars['xx'] + 180. # This 180 degree offset is here to place streams in the "middle of the sky"
    stars['yy'] = stars['yy'] 

    stars = prescreen_sca(stars, scanumber, obm)
    print(stars) 

    scm = SceneModule()

    i = 0
    for line_of_stars in stars:    
        print("line = ", i, " of ", len(stars)) 
        print(line_of_stars)
        stellar = {
                    'n_stars': line_of_stars['mass'],
                    'age_low': line_of_stars['age'], 'age_high': line_of_stars['age'],
                    'z_low': line_of_stars['metallicity'], 'z_high': line_of_stars['metallicity'],
                    'imf': 'salpeter', 'alpha': -2.35,
                    'binary_fraction': 0.1,
                    'distribution': 'uniform', 'clustered': False,
                    'radius': line_of_stars['r_split'], 'radius_units': 'arcsec',
                    'distance_low': 0.1, 'distance_high': 0.11, #<---- these distances are kpc
                    'offset_ra': line_of_stars['xx']*3600., 'offset_dec': line_of_stars['yy']*3600.
                   }

        stellar_cat_file, t_new = scm.CreatePopulation(stellar, id=int(line_of_stars['index']))
        oldfile = 'sim_stars_'+str(int(line_of_stars['index']))+'.fits'
        newfile = 'SCA'+str(scanumber)+'/sim_stars_'+str(int(line_of_stars['index']))+'_SCA'+str(scanumber)+'.fits'
        os.rename(oldfile, newfile) 

        if (i<1): # for the first iteration create a copy of the table 
            t0 = t_new 

        t0 = vstack([t0, t_new])
        del t_new 

        i = i+1 

    return t0, obm

def plot_view(streamfile, scanumber, angle, observation, table): 

    st = rot_star_table(Table.read(streamfile), angle)

    plt.scatter(st['xx']+180., st['yy'])    # the 180 is there because we have moved to the middle of the sky 
    plt.scatter(table['ra'], table['dec'])

    padding = 10. # how far outside each SCA to pick up particles - in arcsec 
    
    for sca in [1,2,3,4, 5, 6, 7,8,9,10,11,12,13,14,15,16,17,18 ]: 
            ra0_sca0 = (observation.__dict__['instrument'].__dict__['detectors'][sca-1].__dict__['header']['RA_APER']) % 360.
            dec0_sca0 = observation.__dict__['instrument'].__dict__['detectors'][sca-1].__dict__['header']['DEC_APER']
            dx = (2044 * 0.11 + padding) / 3600.
            plt.plot([ra0_sca0-dx, ra0_sca0 + dx], [dec0_sca0+dx, dec0_sca0+dx], color='brown')
            plt.plot([ra0_sca0-dx, ra0_sca0 + dx], [dec0_sca0-dx, dec0_sca0-dx], color='brown')
            plt.plot([ra0_sca0-dx, ra0_sca0 - dx], [dec0_sca0-dx, dec0_sca0+dx], color='brown')
            plt.plot([ra0_sca0+dx, ra0_sca0 + dx], [dec0_sca0-dx, dec0_sca0+dx], color='brown')

    plt.savefig('SCA'+str(scanumber)+'/SCA'+str(scanumber)) 



def parse_args():
    parser = argparse.ArgumentParser(description="creates star catalogs for all start particles on a given SCA") 

    parser.add_argument('--scanumber', dest='scanumber', type=int, 
                        help='which SCA do you want to run?')
    parser.set_defaults(scanumber=1)

    parser.add_argument('--streamfile', dest='streamfile', type=str, 
                        help='file containing stream particles?')
    parser.set_defaults(streamfile='table.fits')

    args = parser.parse_args()
    return args





if __name__ == "__main__":

    args = parse_args()

    print("We are going to run star creation for: ", args.scanumber) 

    master_table, obs = create_catalogs(args.streamfile, scanumber=args.scanumber, rotangle=60., 
                    offset_ra=800., offset_dec=1000., offset_pa=0., exptime=10000)

    print(master_table) 
    plot_view(args.streamfile, args.scanumber, 60., obs, master_table)  
    
    master_table.write('SCA'+str(args.scanumber)+'/master_star_table_sca'+str(args.scanumber)+'.fits', overwrite=True) 

