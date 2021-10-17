# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 19:55:16 2021

@author: ar-h1
"""

import matplotlib.pyplot as plt 
import numpy as np
import astropy.io.fits as pyfits 

#hdu = pyfits.open('C:/Users/ar-h1/research/TESS/WASP-4b/WASP-4b_Sector2_Exptime120s/MAST_2020-12-02T1330/TESS/tess2018234235059-s0002-0000000402026209-0121-s/tess2018235142541-s0002-s0002-0000000402026209-00109_dvt.fits')

#hdu = pyfits.open('C:/Users/ar-h1/research/TESS/WASP-4b/WASP-4b_Sector28_Exptime120s/MAST_2020-12-02T1334/TESS/tess2020212050318-s0028-0000000402026209-0190-s/tess2020213081515-s0028-s0028-0000000402026209-00364_dvt.fits')

hdu = pyfits.open('C:/Users/ar-h1/research/TESS/WASP-4b/WASP-4b_Sector29_Exptime120s/MAST_2020-12-02T1335/TESS/tess2020238165205-s0029-0000000402026209-0193-s/tess2020239173514-s0029-s0029-0000000402026209-00382_dvt.fits')


DATA_REL = hdu[0].header['DATA_REL']

LenMedDetrender_hr = hdu[1].header['MEDDETR']

print('Data release version number: %d'%(DATA_REL))
print('Running median filter length: %.2f hour'%(LenMedDetrender_hr))