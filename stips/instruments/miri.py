from __future__ import absolute_import, division
__filetype__ = "detector"

#External Modules
import os

import numpy as np

from astropy.io import fits as pyfits

#Local Modules
from ..astro_image import AstroImage
from .instrument import Instrument
from .jwst_instrument import JwstInstrument
from ..utilities import OffsetPosition

class MIRI(JwstInstrument):
    __classtype__ = "detector"
    """
    The MIRI class contains the necessary constants and modifications to run MIRI
    observations.
        
        detectors : array of detectors, each an AstroImage, and each with its own RA/DEC
        instrument: string, which instrument this is
        filter    : string, what filter of the instrument is being observed
    """

    def __init__(self, **kwargs):
        """
        Still only a base class. Init is super init.
        """
        #Initialize superclass
        super(MIRI,self).__init__(**kwargs)

        #Set oversampling
        self.oversample = kwargs.get('oversample', 1)
        
        self.k = self.K[kwargs.get('miri_mods', 'fast')]
        self.a = self.A[kwargs.get('miri_mods', 'fast')]

    def resetPSF(self):
        import webbpsf
        if self.filter not in self.FILTERS:
            raise ValueError("Filter %s is not a valid MIRI filter" % (self.filter))
        have_psf = False
        if os.path.exists(os.path.join(self.out_path, "psf_cache")):
            if os.path.exists(os.path.join(self.out_path, "psf_cache", "psf_{}_{}_{}.fits".format("MIRI", self.filter, self.oversample))):
                with pyfits.open(os.path.join(self.out_path, "psf_cache", "psf_{}_{}_{}.fits".format("MIRI", self.filter, self.oversample))) as psf:
                    if psf[0].header['VERSION'] >= webbpsf.__version__ and (self.psf_commands is None or self.psf_commands == ''):
                        self.psf = AstroImage(data=psf[0].data, detname="MIRI {} PSF".format(self.filter), logger=self.logger)
                        have_psf = True
        if not have_psf:
            base_state = self.getState()
            self.updateState(base_state+"<br /><span class='indented'>Generating PSF</span>")
            self._log("info", "Creating PSF")
            ins = webbpsf.MIRI()
            self._log("info", "Setting PSF attributes")
            if self.psf_commands is not None and self.psf_commands != '':
                for attribute,value in self.psf_commands.iteritems():
                    self._log("info", "Setting PSF attribute {} to {}".format(attribute, value))
                    setattr(ins,attribute,value)
            self._log("info", "Setting PSF filter to '{}'".format(self.filter))
            ins.filter = self.filter
            max_safe_size = int(np.floor(30. * self.PHOTPLAM[self.filter] / (2. * self.SCALE[0])))
            max_ins_size = max(self.DETECTOR_SIZE) * self.oversample
            max_conv_size = int(np.floor(self.convolve_size / (2*self.oversample)))
            psf_size = min(max_safe_size, max_ins_size, max_conv_size)
            self._log("info", "PSF choosing between {}, {}, and {}, chose {}".format(max_safe_size, max_ins_size, max_conv_size, psf_size))
            psf = ins.calcPSF(oversample=self.oversample,fov_pixels=psf_size)
            psf[0].header['VERSION'] = webbpsf.__version__
            if os.path.exists(os.path.join(self.out_path, "psf_cache")):
                dest = os.path.join(self.out_path, "psf_cache", "psf_{}_{}_{}.fits".format("MIRI", self.filter, self.oversample))
                pyfits.writeto(dest, psf[0].data, header=psf[0].header, overwrite=True)
            self.psf = AstroImage(data=psf[0].data,detname="MIRI %s PSF" % (self.filter),logger=self.logger)
            self.updateState(base_state)
        
    def generateReadnoise(self):
        """
        Readnoise formula that is similar to JWST ETC.
        """
        k = self.k
        a = self.a

        if self.exptime > 1000:
            numGroups = np.ceil(self.exptime/1000.0)
            timePerGroup = self.exptime/numGroups
        else:
            numGroups = 1
            timePerGroup = self.exptime

        rdns = np.sqrt(numGroups) * k * timePerGroup**(a)
        return rdns
    
    @classmethod
    def handleDithers(cls,form):
        """
        Handles MIRI-specific dither patterns
        """
        dither_pattern = form['dither_type']
        if dither_pattern == "CYCLING":
            dither_points = int(form['dither_points'])
            dither_size = form['dither_size']
            initial_dithers = []
            while len(initial_dithers) <= dither_points:
                initial_dithers += cls.DITHER_OFFSETS[dither_pattern][dither_size]
            initial_dithers = initial_dithers[:dither_points]
        elif dither_pattern == "REULEAUX" or dither_pattern == "GAUSSIAN":
            dither_size = form['dither_size']
            initial_dithers = cls.DITHER_OFFSETS[dither_pattern][dither_size]
        else: #no dither pattern
            initial_dithers = [(0.,0.)]
        dither_subpixel = form['dither_subpixel']
        if dither_subpixel == 'WITH':
            subpixel_dithers = cls.DITHER_OFFSETS['SUBPIXEL']
        else:
            subpixel_dithers = [(0.,0.)]
        return cls.doSubpixel(initial_dithers,subpixel_dithers)

    
    INSTRUMENT = "MIRI"
    DETECTOR = "MIRI"
    K = {'fast':43.78, 'slow':27.84}
    A = {'fast':-0.257, 'slow':-0.257}

    # Offsets are in (arcseconds_ra,arcseconds_dec,degrees_angle)
    # Instrument Offset is from JWST V2,V3 centre
    INSTRUMENT_OFFSET = (-433.68,-375.6,4.55)
    # Detector Offsets are RELATIVE TO INSTRUMENT OFFSET
    DETECTOR_OFFSETS = ((0.,0.,0.),)
    # N_DETECTORS is a set of options on how many of the instrument's detectors you want to use
    N_DETECTORS = [1]
    OFFSET_NAMES = ("MIRI",)
    DETECTOR_SIZE = (1024,1024) #pixels
    PIXEL_SIZE = 25.0 #um
    SCALE = [0.11,0.11]
    DISTORTION =    {
                        'MIRI': {
                                    'DIST_A':  [[   0.,             0.,             1.163e-7,       -2.86545e-9,    1.99e-12],
                                                [   0.,             -1.74755e-8,    2.52e-11,       -1.56e-13,      0.],
                                                [   -8.13208e-7,    -1.599e-9,      1.99e-12,       0.,             0.],
                                                [   -2.38e-12,      1.07e-13,       0.,             0.,             0.],
                                                [   2.95e-13,       0.,             0.,             0.,             0.]],
                                    'DIST_B':  [[   0.,             0.,             2.346e-8,       5.28e-11,       -1.50e-13],
                                                [   0.,             -1.63366e-6,    -8.07e-10,      1.37e-12,       0.],
                                                [   -2.848e-8,      4.81e-11,       -1.14e-13,      0.,             0.],
                                                [   5.88e-11,       1.10e-12,       0.,             0.,             0.],
                                                [   9.91e-14,       0.,             0.,             0.,             0.]]
                                }
                    }
    FILTERS = ('F560W','F770W','F1000W','F1130W','F1280W','F1500W','F1800W','F2100W','F2550W')
    DEFAULT_FILTER = 'F1000W'
    FLATFILE = 'err_flat_test.fits'
    DARKFILE = 'err_rdrk_miri_im.fits'  # ETC short
    BACKGROUND = {  'none': {  'F560W': 0., 'F770W': 0., 'F1000W' :0., 'F1130W': 0., 'F1280W': 0., 'F1500W': 0., 'F1800W':0., 
                                'F2100W': 0., 'F2550W': 0.},
                    'avgt': {   'F560W': 2.466E+00,  'F770W': 2.372E+01,  'F1000W': 3.584E+01, 'F1130W': 1.978E+01, 
                                'F1280W': 6.272E+01, 'F1500W': 1.645E+02, 'F1800W': 1.796E+02, 'F2100W': 5.656E+02, 
                                'F2550W': 1.286E+03}
                 }
    BACKGROUNDS_V = ['none', 'avg', 'med', 'max', 'min']
    BACKGROUNDS = ['None', 'Average zodiacal and thermal background', 'Median zodiacal and thermal background', 'Maximum zodiacal and thermal background', 'Minimum zodiacal and thermal background']
    BGTEXT = {'none': 'None', 'avg': 'Average zodiacal and thermal background', 'med': 'Median zodiacal and thermal background', 'max': 'Maximum zodiacal and thermal background', 'min': 'Minimum zodiacal and thermal background'}
    #PHOTFNU has units of Jy
    PHOTFNU = { 'F560W':1.872E-07,  'F770W':8.608E-08, 'F1000W':1.406E-07, 'F1130W':3.707E-07, 
                'F1280W':1.648E-07, 'F1500W':1.083E-07, 'F1800W':1.996E-07, 'F2100W':1.761E-07, 
                'F2550W':3.091E-07
              }
    #PHOTPLAM has units of um
    PHOTPLAM = {    'F560W':5.6270, 'F770W':7.6329, 'F1000W':9.9417, 'F1130W':11.3493, 
                    'F1280W':12.8919, 'F1500W':15.1981, 'F1800W':17.9315, 'F2100W':20.8882, 
                    'F2550W':25.2396
               }
    DITHERS = ("SUBPIXEL ONLY", "CYCLING", "REULEAUX", "GAUSSIAN")
    DITHER_POINTS = {
                        "SUBPIXEL ONLY": ["0"],
                        "CYCLING": [str(x+1) for x in range(312)], #if people want multiple cycles, do multiple observations
                        "REULEAUX": ["12"],
                        "GAUSSIAN": ["5"]
                     }
    DITHER_SIZE = {
                    "SUBPIXEL ONLY": ["STANDARD"],
                    "CYCLING": ["SMALL","MEDIUM","LARGE"],
                    "REULEAUX": ["SMALL","MEDIUM","LARGE"],
                    "GAUSSIAN": ["SMALL","MEDIUM","LARGE"]
                  }
    DITHER_SUBPIXEL = {
                        "SUBPIXEL ONLY": ["WITH"],
                        "CYCLING": ["WITH","WITHOUT"],
                        "REULEAUX": ["WITH","WITHOUT"],
                        "GAUSSIAN": ["WITH","WITHOUT"]
                      }
    DITHER_OFFSETS = {
                        "CYCLING":  {
                                        "SMALL":    [(-0.024,0.036),(0.127,0.030),(0.880,-0.275),(0.165,-0.770),(0.550,-0.440),
                                                     (-0.935,-0.715),(-1.210,-0.275),(-0.275,-0.770),(-0.440,-1.210),
                                                     (0.055,1.045),(1.100,-0.275),(-0.605,0.110),(-0.770,-0.770),(0.715,-0.605),
                                                     (0.880,0.055),(-0.055,-0.550),(0.770,0.770),(-0.825,-0.825),
                                                     (-1.100,-0.275),(-0.495,-0.770),(-0.770,0.770),(0.825,0.275),
                                                     (-0.770,0.275),(0.495,-0.220),(-1.210,0.110),(-0.385,-0.825),
                                                     (0.000,-0.715),(-0.825,0.880),(-0.990,-0.770),(1.045,0.605),
                                                     (-0.880,0.055),(-1.155,0.330),(0.110,-0.550),(1.155,0.055),(0.770,-0.495),
                                                     (0.165,-0.880),(-0.880,0.000),(-1.155,0.385),(0.110,-0.495),(1.045,0.110),
                                                     (-0.990,-0.440),(0.165,1.045),(-0.660,1.045),(-0.715,0.440),(0.550,-0.110),
                                                     (0.605,-0.495),(-0.990,-0.495),(-0.385,0.660),(0.440,0.770),
                                                     (-0.165,-0.825),(-0.330,1.155),(0.495,-0.770),(0.660,0.220),
                                                     (1.045,-0.495),(0.660,0.825),(-1.045,0.440),(-0.550,0.880),
                                                     (-0.715,-0.275),(-0.550,0.275),(0.825,-0.660),(0.110,0.440),
                                                     (-0.385,0.605),(0.440,-0.605),(0.495,0.990),(-0.330,0.000),
                                                     (-0.605,0.275),(-0.220,0.935),(-1.045,-0.330),(0.000,-0.220),
                                                     (0.715,0.275),(-0.550,-0.165),(0.495,0.000),(1.100,0.440),(-1.045,0.605),
                                                     (-0.660,0.715),(0.605,0.440),(-0.440,0.000),(-0.495,-0.605),(0.990,0.605),
                                                     (-0.935,0.330),(0.000,0.330),(-0.055,-0.605),(-0.220,0.825),
                                                     (-0.385,-0.330),(0.220,0.990),(0.715,-0.935),(0.550,-0.165),
                                                     (-0.715,0.330),(-0.770,-0.660),(-0.385,-0.275),(0.330,-0.715),
                                                     (0.385,0.440),(0.660,-0.660),(0.935,-0.055),(-0.990,-0.495),
                                                     (-0.715,-0.110),(-0.220,-0.550),(0.165,-0.825),(0.000,0.935),
                                                     (-0.495,0.770),(0.000,-0.220),(-0.825,-0.385),(-0.880,-0.055),
                                                     (0.165,0.990),(0.220,-0.770),(0.935,-0.495),(-0.770,0.495),(0.055,-0.330),
                                                     (-0.770,-0.550),(1.155,0.495),(0.660,0.275),(0.715,0.880),(-0.110,0.440),
                                                     (0.275,0.055),(0.110,-0.825),(-0.825,-0.330),(-0.550,0.550),(0.605,0.715),
                                                     (-0.330,-1.155),(0.055,-0.110),(0.000,-0.660),(0.825,-0.825),(0.770,0.165),
                                                     (0.715,0.660),(-0.660,0.550),(-0.055,0.275),(-0.990,0.385),(-0.055,-1.210),
                                                     (0.770,0.000),(0.605,0.935),(0.220,0.495),(-0.495,-0.110),(-0.880,0.000),
                                                     (0.825,0.165),(0.110,-1.045),(0.055,-0.660),(0.440,0.880),(0.385,0.495),
                                                     (-0.330,0.275),(-0.935,-0.110),(1.210,0.110),(-0.715,0.825),
                                                     (-0.330,-0.275),(0.715,0.110),(0.660,0.880),(0.275,-0.935),(1.100,0.165),
                                                     (0.275,0.330),(-0.660,-0.330),(0.055,1.045),(-1.100,0.385),(-1.155,-0.220),
                                                     (0.330,-0.880),(-0.275,1.045),(0.880,-0.605),(0.825,0.440),(-0.990,0.440),
                                                     (-0.935,-0.385),(0.220,0.495),(0.605,0.000),(-0.550,-0.880),(0.495,-1.045),
                                                     (-0.440,-0.055),(0.935,0.330),(0.440,-0.550),(0.165,0.935),(-0.660,0.495),
                                                     (-0.605,-0.440),(-0.990,-0.440),(0.055,-0.605),(-0.880,0.275),
                                                     (0.165,0.990),(0.550,-0.440),(-0.605,-0.385),(-0.330,0.605),
                                                     (-0.605,-1.100),(0.990,0.110),(0.275,-0.275),(0.110,0.495),(-0.825,-0.660),
                                                     (-0.110,0.990),(0.715,-0.495),(0.770,0.055),(-1.045,-0.110),(0.330,-0.990),
                                                     (-0.825,-0.605),(-0.330,1.045),(0.275,-0.220),(0.770,-0.550),(0.385,0.275),
                                                     (0.220,-0.825),(-0.165,-0.660),(1.100,-0.220),(-0.385,0.495),
                                                     (-1.100,0.495),(-0.715,0.550),(0.220,0.660),(0.055,-0.605),(0.330,-1.045),
                                                     (0.605,-0.770),(-0.660,-0.660),(0.715,0.055),(-0.880,0.825),(-0.605,0.330),
                                                     (0.440,0.550),(0.715,-0.275),(0.220,-0.715),(-0.385,0.990),(-0.880,-0.880),
                                                     (-0.605,0.715),(-1.100,-0.495),(-0.165,0.550),(0.220,-0.990),
                                                     (-1.155,0.055),(0.990,0.275),(-0.935,0.440),(-0.770,-0.220),(-0.715,0.825),
                                                     (0.220,-1.045),(0.165,0.770),(0.660,-0.440),(-1.155,0.275),(-0.990,-0.495),
                                                     (1.045,-0.220),(0.550,-0.880),(0.055,0.495),(0.000,-0.935),(0.715,-0.110),
                                                     (-0.880,-0.440),(1.045,0.385),(0.770,-0.715),(0.275,-0.220),
                                                     (-0.440,-0.880),(0.385,0.605),(0.110,-0.715),(-0.825,-0.330),
                                                     (0.880,-0.660),(0.495,-0.825),(0.770,0.165),(0.275,0.990),(-0.770,-0.990),
                                                     (-0.825,0.275),(0.110,-0.605),(-0.825,-0.220),(0.660,0.000),(0.055,1.045),
                                                     (0.330,-1.045),(0.715,0.550),(-0.880,-0.770),(0.935,-0.275),
                                                     (-0.220,-0.495),(-0.935,-0.330),(-0.220,0.220),(-1.155,-0.055),
                                                     (-0.990,-0.715),(0.825,0.330),(-0.220,0.550),(-0.825,-0.165),(0.330,0.385),
                                                     (-0.605,-0.440),(0.660,0.770),(-0.825,0.275),(-0.220,0.715),(-0.275,0.330),
                                                     (-0.220,-0.990),(0.825,-0.825),(0.330,-0.495),(-0.715,-0.770),
                                                     (0.880,-0.440),(0.605,0.935),(0.660,0.385),(-0.495,0.440),(0.440,-0.550),
                                                     (1.155,-0.165),(0.770,-0.715),(-0.605,-0.550),(-0.440,0.220),(0.165,0.385),
                                                     (-0.330,-0.825),(-0.715,0.440),(0.990,-0.110),(0.825,-0.605),
                                                     (-0.660,-0.825),(-0.275,0.440),(0.220,-0.330),(0.165,1.155),(0.330,0.825),
                                                     (-0.825,0.000),(-0.990,0.440),(0.055,-1.045),(0.440,0.165),(0.935,-0.550),
                                                     (0.000,0.770),(0.715,-0.825),(-0.440,1.045),(0.055,-0.220),(0.660,0.000),
                                                     (0.715,-0.385),(-0.330,0.055),(1.045,-0.660),(1.100,-0.220),(-0.165,0.715),
                                                     (0.550,-0.055),(0.495,-0.660),(-0.220,0.220),(-0.935,0.055),(0.880,-0.165),
                                                     (0.935,-0.660),(-0.880,-0.550),(-0.385,-0.495),(-0.550,-0.825)],
                                        "MEDIUM":   [(0.387,-0.484),(0.054,0.024),(1.210,4.565),(-3.795,9.625),(3.850,-3.190),
                                                     (0.935,3.190),(-5.500,0.495),(1.155,2.145),(1.100,-3.520),(-5.775,1.100),
                                                     (3.850,-1.375),(-3.905,-6.985),(8.140,-1.870),(-5.115,6.820),(1.760,1.265),
                                                     (3.355,-6.545),(-3.850,2.310),(3.795,0.000),(1.980,2.035),(2.255,-1.595),
                                                     (0.330,2.750),(1.485,0.000),(-3.190,0.605),(-0.495,-1.705),(0.000,0.660),
                                                     (-6.875,0.000),(2.530,6.325),(1.485,-1.815),(-1.650,3.630),(-5.445,-0.440),
                                                     (0.000,-0.275),(-3.465,-1.705),(-0.330,-1.100),(4.345,-0.880),
                                                     (0.880,0.385),(0.825,2.475),(2.420,4.510),(-0.165,7.590),(0.440,-3.355),
                                                     (-1.925,2.585),(-5.170,-3.520),(-1.925,-2.970),(2.200,-1.265),
                                                     (3.575,-1.155),(0.880,-5.060),(-0.715,2.420),(-4.070,-1.815),(1.485,3.795),
                                                     (3.630,6.930),(-1.045,-9.460),(5.830,1.595),(-11.055,-1.375),(0.110,5.060),
                                                     (-2.585,-5.280),(0.880,-4.235),(-0.825,-1.265),(3.300,-0.110),
                                                     (-2.695,-3.740),(-2.420,2.805),(5.115,-1.045),(2.750,-1.980),
                                                     (0.055,-3.630),(0.110,-3.355),(8.745,0.935),(-1.980,4.620),(-0.055,1.980),
                                                     (4.290,-0.385),(4.015,0.935),(0.770,4.180),(-4.125,2.420),(-1.650,-0.055),
                                                     (-4.235,0.715),(2.310,1.760),(5.665,1.430),(2.640,0.715),(-2.365,1.485),
                                                     (-0.550,4.950),(5.445,-0.220),(-0.770,-1.595),(-4.785,0.055),(1.210,1.430),
                                                     (7.315,1.320),(-2.310,-3.575),(3.465,-1.265),(-1.540,3.410),(-0.275,6.820),
                                                     (0.110,3.135),(5.885,-11.715),(-2.640,0.550),(-0.385,0.440),(3.630,0.605),
                                                     (3.465,0.055),(-4.620,-5.280),(1.155,-0.330),(4.840,3.575),(-0.605,-4.895),
                                                     (-1.100,-0.550),(-2.915,-1.210),(-1.980,-1.265),(-3.245,0.935),
                                                     (2.200,3.520),(0.165,-0.550),(2.860,2.365),(-0.495,5.335),(-3.850,0.000),
                                                     (8.415,-0.110),(4.070,5.335),(2.145,-5.775),(-2.090,-3.630),(-0.165,3.300),
                                                     (0.000,2.805),(0.165,2.145),(3.080,4.620),(-3.795,-10.780),(0.880,5.885),
                                                     (-1.705,4.015),(1.210,-0.220),(-1.045,6.600),(0.770,4.015),(5.225,4.345),
                                                     (-3.300,-0.990),(3.245,2.970),(-2.640,-6.875),(2.255,4.015),(2.750,-0.550),
                                                     (0.605,-3.190),(-5.060,-2.475),(0.605,1.595),(2.640,1.430),(3.355,-2.530),
                                                     (-2.750,1.815),(4.785,-1.595),(-2.860,1.100),(2.365,2.530),(4.070,1.815),
                                                     (3.355,-4.235),(-1.540,-0.660),(0.715,0.440),(0.550,2.255),(1.485,-0.605),
                                                     (-1.760,4.290),(2.805,-2.860),(-3.300,-3.575),(-7.865,4.675),
                                                     (-1.980,7.040),(-0.605,-0.330),(-0.550,3.245),(2.915,-0.495),
                                                     (0.770,-1.430),(5.225,8.470),(-4.510,-2.585),(-2.695,-1.925),(1.210,2.090),
                                                     (-1.045,-3.630),(2.310,-4.675),(-2.475,1.705),(1.320,-2.420),
                                                     (-4.015,5.610),(-2.640,-0.165),(-1.375,-0.825),(-1.100,1.430),
                                                     (-5.335,-0.440),(-0.990,2.365),(1.595,3.245),(-0.220,7.810),(-0.385,2.200),
                                                     (3.190,5.555),(6.765,-2.915),(3.850,-4.840),(-1.045,3.080),(2.310,-0.605),
                                                     (2.915,2.805),(5.280,-0.660),(-0.935,-4.620),(3.630,2.695),(-0.605,-0.165),
                                                     (2.530,-6.710),(-6.435,-0.110),(-2.970,4.565),(2.475,-1.595),
                                                     (-3.850,-1.650),(-3.025,-0.110),(1.650,0.825),(0.825,5.665),
                                                     (-5.500,-0.770),(-1.045,2.090),(-2.640,5.555),(-0.275,-10.175),
                                                     (0.000,1.760),(-5.005,-2.860),(-1.100,-1.045),(-0.605,-8.525),
                                                     (4.180,-2.750),(-6.875,3.740),(4.620,6.105),(4.235,-2.035),(5.940,-1.540),
                                                     (2.145,1.650),(-0.660,-3.135),(8.195,-3.135),(1.320,-1.540),(-2.585,4.180),
                                                     (5.830,2.695),(-3.465,-1.925),(6.160,3.850),(-1.375,-1.980),
                                                     (-4.620,-1.265),(-0.385,-2.145),(3.850,3.410),(2.365,-3.630),
                                                     (3.630,-2.035),(-2.695,-4.565),(2.530,6.710),(-1.155,-1.210),
                                                     (3.190,-0.825),(1.265,4.235),(0.000,-0.660),(1.265,5.500),(3.850,1.595),
                                                     (6.325,-2.035),(2.090,-0.990),(6.875,0.880),(-2.310,-3.795),(2.255,1.595),
                                                     (1.870,0.440),(2.805,1.210),(-2.640,1.595),(-1.925,2.475),(-3.740,1.320),
                                                     (-0.055,1.430),(-0.880,-2.475),(4.125,3.795),(1.760,-3.520),(-3.905,4.730),
                                                     (0.660,1.485),(-2.805,3.465),(-2.530,2.860),(2.585,-5.500),(7.920,0.495),
                                                     (3.135,-0.605),(-0.330,3.410),(-4.565,-6.050),(6.050,-1.705),
                                                     (9.295,-4.455),(2.530,-1.540),(-2.365,-1.650),(0.440,-0.165),
                                                     (-1.265,0.385),(1.100,5.940),(3.575,5.060),(7.700,3.685),(-1.155,5.665),
                                                     (-0.990,-4.950),(2.585,-3.740),(-2.090,-0.605),(-0.715,-0.715),
                                                     (-8.910,-3.190),(7.865,-5.280),(5.060,-5.225),(-7.095,-0.495),
                                                     (-3.080,3.520),(-0.715,-6.600),(3.850,-1.925),(-2.475,1.925),
                                                     (-1.540,-0.330),(5.995,-2.090),(-8.470,2.255),(3.025,3.795),
                                                     (-0.440,-3.300),(-2.255,-6.820),(1.540,-2.695),(3.465,1.375),
                                                     (2.420,-2.970),(7.425,-0.330),(-0.990,-2.805),(6.325,0.715),(4.290,2.860),
                                                     (4.125,0.550),(-0.440,0.605),(-2.585,0.385),(0.000,3.960),(0.495,-2.750),
                                                     (-3.630,0.275),(-1.375,3.905),(-2.970,-1.980),(0.935,-4.180),
                                                     (0.330,-2.255),(-1.925,4.015),(3.960,-4.290),(-6.215,6.270),(2.530,5.225),
                                                     (-3.465,3.355),(-3.080,-3.630),(1.375,3.410),(-0.770,4.675),(5.555,-0.385),
                                                     (5.500,-2.200),(-1.815,-1.320),(-1.320,-0.385),(-0.825,-6.765),
                                                     (4.180,2.090),(-4.015,3.080),(9.130,1.045),(3.905,-2.255),(-0.330,-1.430),
                                                     (-4.235,-1.540),(1.650,2.255),(0.385,-4.125),(-0.880,-1.980),
                                                     (0.495,-2.420),(-2.310,4.235)],
                                        "LARGE":    [(-0.194,-0.448),(0.103,0.048),(2.530,9.075),(7.865,6.875),(7.700,-6.380),
                                                     (-8.525,-3.190),(-11.110,0.935),(2.365,4.345),(3.960,3.960),
                                                     (-11.495,2.310),(7.810,-2.805),(0.055,5.665),(-0.660,9.900),
                                                     (-1.265,-0.770),(8.580,-1.925),(7.755,-4.455),(-7.810,4.620),(7.535,0.110),
                                                     (2.200,-7.205),(-0.935,-3.465),(0.660,5.500),(2.915,0.000),(-6.490,1.155),
                                                     (4.565,-3.135),(0.000,1.320),(-13.695,0.110),(5.170,12.705),(2.915,-3.685),
                                                     (-3.300,7.260),(-10.835,-0.990),(0.000,-0.605),(-6.985,-3.355),
                                                     (-0.660,-2.310),(3.575,2.420),(1.760,0.825),(1.705,4.895),(4.840,9.020),
                                                     (-0.825,3.190),(0.880,-6.765),(-3.795,5.225),(4.510,-2.530),
                                                     (-3.795,-6.050),(-10.340,-7.205),(7.095,-2.255),(1.760,-10.120),
                                                     (0.165,-7.260),(-8.140,-3.685),(2.915,7.535),(7.370,13.970),
                                                     (-1.705,-2.530),(11.770,3.135),(3.465,-1.815),(0.220,10.120),
                                                     (-5.225,-10.670),(1.760,-8.415),(10.175,2.365),(6.600,-0.330),
                                                     (1.595,8.470),(-4.840,5.555),(10.285,-2.035),(5.500,-4.070),(-1.485,4.840),
                                                     (-5.280,-7.645),(8.525,6.655),(11.330,2.860),(-0.055,3.960),(8.580,-0.825),
                                                     (8.085,1.815),(0.220,-6.710),(-8.305,4.840),(-3.300,-0.055),(-8.415,1.375),
                                                     (4.730,3.520),(-4.125,9.240),(5.280,1.485),(-4.785,2.915),(-1.100,9.900),
                                                     (-1.705,-3.190),(10.780,-0.495),(-9.515,0.055),(2.530,2.860),(6.435,8.800),
                                                     (-4.620,-7.205),(6.875,-2.475),(-3.190,6.930),(-0.495,13.750),
                                                     (0.330,6.325),(7.315,1.265),(-5.390,1.210),(-0.715,0.880),(4.180,-11.605),
                                                     (6.985,0.165),(-0.330,6.600),(2.365,-0.660),(9.680,7.205),(-7.755,-14.025),
                                                     (-2.310,-1.210),(-5.885,-2.420),(-3.960,-2.585),(-6.435,1.815),
                                                     (4.400,7.040),(0.275,-1.100),(5.830,4.785),(-1.045,10.615),(-7.810,0.000),
                                                     (-1.815,11.660),(8.250,10.615),(-4.235,-7.425),(-3.190,6.710),
                                                     (-9.405,-10.560),(-1.100,-9.735),(0.385,4.235),(6.270,9.350),(6.325,2.200),
                                                     (4.400,7.975),(-6.765,-2.035),(2.530,-0.440),(-2.035,13.200),(1.540,8.085),
                                                     (3.025,-1.155),(-3.410,7.920),(6.435,6.050),(-5.280,-13.695),
                                                     (1.815,11.825),(5.610,-1.100),(1.155,-6.380),(-10.230,-4.895),
                                                     (1.155,3.245),(5.390,2.970),(6.765,-5.060),(-5.610,3.575),(9.625,-3.245),
                                                     (-5.830,2.310),(4.785,5.060),(8.250,3.685),(6.765,-8.525),(-3.190,-1.320),
                                                     (1.485,0.990),(1.100,4.565),(10.395,8.635),(-3.630,8.580),(5.665,-5.830),
                                                     (-6.600,-7.205),(-2.915,6.435),(5.280,-8.250),(-10.175,13.640),
                                                     (-1.210,6.435),(5.885,-1.045),(-2.640,-1.540),(-4.015,3.740),
                                                     (-9.020,-5.225),(-5.335,-3.795),(2.420,4.290),(-2.145,-7.260),
                                                     (4.730,-9.295),(-4.895,3.355),(2.750,-4.950),(-8.085,11.330),
                                                     (-5.280,-0.385),(1.595,-3.025),(-2.310,2.970),(-10.615,-0.880),
                                                     (-2.090,4.785),(-2.145,6.215),(-0.330,-4.070),(4.785,-1.210),
                                                     (6.490,11.165),(13.475,-5.885),(7.700,-9.680),(3.135,6.380),(-0.770,4.455),
                                                     (5.775,5.665),(10.560,-1.320),(-11.165,-1.650),(7.370,5.335),
                                                     (-1.265,-0.275),(5.060,-13.530),(-12.815,-0.330),(-5.940,9.075),
                                                     (4.895,-3.135),(-6.050,-0.330),(-7.755,-3.410),(0.000,3.685),
                                                     (1.595,11.275),(-1.870,-9.240),(-2.145,4.290),(-5.390,11.165),
                                                     (-0.055,-5.115),(-6.930,-3.740),(-10.065,-5.720),(-2.310,-2.145),
                                                     (-0.385,8.525),(8.470,-5.500),(-13.695,7.480),(9.350,12.155),
                                                     (8.525,-4.015),(11.990,-3.190),(4.235,3.410),(-1.430,-6.215),(6.105,4.235),
                                                     (-5.060,8.470),(2.695,-3.080),(11.660,5.445),(3.355,1.705),(12.320,7.700),
                                                     (-2.695,-4.070),(-9.350,-2.585),(-0.825,-4.235),(-4.840,3.300),
                                                     (4.675,-7.370),(7.260,-4.125),(-5.445,-9.185),(5.060,13.530),
                                                     (-2.365,-2.420),(6.490,-1.595),(2.585,8.415),(0.000,-1.320),
                                                     (12.595,-3.960),(7.700,3.245),(3.685,-7.095),(4.180,-1.980),(8.305,7.480),
                                                     (-4.730,-7.645),(4.455,3.245),(3.740,0.880),(5.665,2.420),(12.210,-3.465),
                                                     (-3.795,4.895),(-7.590,2.750),(-0.055,2.860),(-1.870,-5.005),
                                                     (13.805,1.815),(2.530,11.110),(-7.865,9.570),(1.320,3.025),(-5.555,6.875),
                                                     (-5.170,5.720),(5.115,-11.000),(9.240,-4.125),(6.325,-1.155),
                                                     (-0.770,6.820),(-9.185,-12.100),(-5.390,3.245),(-2.035,-9.955),
                                                     (11.990,-4.180),(-4.675,-3.410),(0.880,-0.275),(-2.585,0.715),
                                                     (2.200,11.990),(7.205,10.120),(6.050,-12.925),(-2.365,11.275),
                                                     (-7.920,2.860),(5.115,-7.480),(-4.290,-1.155),(-1.375,-1.485),
                                                     (-9.020,-2.200),(-0.275,13.530),(10.230,-10.505),(-3.245,-5.665),
                                                     (-6.270,7.040),(-1.485,-13.310),(7.700,-3.905),(-5.005,3.905),
                                                     (-3.080,-0.660),(5.225,-3.080),(-8.580,7.315),(6.105,7.535),
                                                     (-0.990,-6.710),(-4.565,-13.750),(3.190,-5.445),(6.875,2.805),
                                                     (4.840,-5.940),(-3.025,6.930),(-2.090,-5.665),(12.595,1.375),(8.690,5.830),
                                                     (8.195,1.100),(-0.990,1.265),(-5.225,0.825),(0.000,8.030),(1.045,-5.500),
                                                     (-7.370,0.605),(-2.695,7.755),(-5.940,-3.960),(1.925,-8.470),
                                                     (0.660,-4.455),(-3.905,7.975),(8.030,-8.690),(-12.485,12.540),
                                                     (5.060,10.395),(-6.875,6.765),(-4.620,8.360),(2.695,6.820),(-1.540,9.295),
                                                     (11.165,-0.715),(11.110,-4.400),(-3.685,-2.750),(-2.640,-0.825),
                                                     (-1.595,-13.475),(8.470,4.180),(-7.975,6.270),(-2.310,-2.365),
                                                     (6.765,-13.145),(-0.660,-2.860),(0.935,-4.950),(3.410,4.565),
                                                     (0.825,-8.305),(6.930,-8.910),(1.925,6.490),(-6.160,-7.315)]
                                    },
                        "REULEAUX": {
                                        "SMALL": [(0.660,-1.320),(0.990,-0.715),(1.155,0.000),(1.045,0.715),(0.660,1.320),
                                                  (0.000,1.265),(-0.605,0.880),(-1.155,0.495),(-1.430,0.000),(-1.100,-0.495),
                                                  (-0.605,-0.880),(0.055,-1.265)],
                                        "MEDIUM": [(1.430,-2.640),(2.035,-1.320),(2.200,0.055),(2.035,1.375),(1.430,2.640),
                                                   (0.165,2.420),(-1.100,1.925),(-2.255,1.045),(-2.970,0.000),(-2.255,-0.990),
                                                   (-1.100,-1.925),(0.165,-2.475)],
                                        "LARGE": [(2.970,-5.280),(4.125,-2.640),(4.400,0.055),(4.125,2.695),(2.970,5.280),
                                                  (0.275,4.840),(-2.200,3.905),(-4.455,2.145),(-6.050,0.000),(-4.455,-2.090),
                                                  (-2.200,-3.905),(0.275,-4.895)]
                                    },
                        "GAUSSIAN": {
                                        "SMALL": [(0.000,-0.110),(-0.275,-0.990),(-1.210,0.385),(0.165,0.165),(1.238,0.578)],
                                        "MEDIUM": [(0.000,-0.220),(-0.605,-1.980),(-2.310,0.825),(0.275,0.275),(2.547,1.238)],
                                        "LARGE": [(0.000,-0.440),(-1.155,-4.070),(-4.730,1.595),(0.495,0.605),(5.198,2.448)]
                                    },
                        "SUBPIXEL": [(0.000,0.000),(0.275,0.000),(0.440,0.275),(0.165,0.275)]
                     }
