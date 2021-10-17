# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 14:18:30 2020

@author: ar-h1
"""

import numpy as np
import lightkurve as lk
import matplotlib.pyplot as plt
import pickle 
import astropy.io.fits as pyfits 
import os 
 
### For XO-6b
# FileToLoad = 'XO-6b/MAST_2020-04-14T1401/MAST_2020-04-14T1401/TESS/tess2019331140908-s0019-0000000138168780-0164-s/tess2019332134924-s0019-s0019-0000000138168780-00278_dvt.fits'
# InitialLCFileToLoad = 'XO-6b/MAST_2020-04-14T1401/MAST_2020-04-14T1401/TESS/tess2019331140908-s0019-0000000138168780-0164-s/tess2019331140908-s0019-0000000138168780-0164-s_lc.fits'

# FileToLoad = 'XO-6b/MAST_2020-04-14T1411/MAST_2020-04-14T1411/TESS/tess2019357164649-s0020-0000000138168780-0165-s/tess2019358235523-s0020-s0020-0000000138168780-00285_dvt.fits'
# InitialLCFileToLoad = 'XO-6b/MAST_2020-04-14T1411/MAST_2020-04-14T1411/TESS/tess2019357164649-s0020-0000000138168780-0165-s/tess2019357164649-s0020-0000000138168780-0165-s_lc.fits'

# FileToLoad = 'XO-6b/Sector26/MAST_2020-08-07T1535/TESS/tess2020160202036-s0026-0000000138168780-0188-s/tess2020161181517-s0026-s0026-0000000138168780-00350_dvt.fits'
# InitialLCFileToLoad = 'XO-6b/Sector26/MAST_2020-08-07T1535/TESS/tess2020160202036-s0026-0000000138168780-0188-s/tess2020160202036-s0026-0000000138168780-0188-s_lc.fits'

### For WASP-12b
#FileToLoad = 'WASP12bMAST_2020-09-03T1338/MAST_2020-09-03T1338/TESS/tess2019357164649-s0020-0000000086396382-0165-s/tess2019358235523-s0020-s0020-0000000086396382-00285_dvt.fits'
#InitialLCFileToLoad = 'WASP12bMAST_2020-09-03T1338/MAST_2020-09-03T1338/TESS/tess2019357164649-s0020-0000000086396382-0165-s/tess2019357164649-s0020-0000000086396382-0165-s_lc.fits'

# FileToLoad = 'WASP-12b_RedownloadedToCheckVersion_MAST_2020-11-10T1302/MAST_2020-11-10T1302/TESS/tess2019357164649-s0020-0000000086396382-0165-s/tess2019358235523-s0020-s0020-0000000086396382-00285_dvt.fits'
# InitialLCFileToLoad = 'WASP-12b_RedownloadedToCheckVersion_MAST_2020-11-10T1302/MAST_2020-11-10T1302/TESS/tess2019357164649-s0020-0000000086396382-0165-s/tess2019357164649-s0020-0000000086396382-0165-s_lc.fits'

# sector = 20

#### WASP-4b


# name = 'WASP-4b'

# ###Exp time 120 s
# exptime = 120
# # sector = 2
# # path = '%s/%s_Sector%d_Exptime120s/MAST_2020-12-02T1330/TESS/tess2018234235059-s0002-0000000402026209-0121-s'%(name,name,sector)
# # InitialLCFileToLoad = '%s/tess2018234235059-s0002-0000000402026209-0121-s_lc.fits'%(path)
# # FileToLoad = '%s/tess2018235142541-s0002-s0002-0000000402026209-00109_dvt.fits'%(path)

# # sector = 28
# # path = '%s/%s_Sector%d_Exptime120s/MAST_2020-12-02T1334/TESS/tess2020212050318-s0028-0000000402026209-0190-s'%(name,name,sector)
# # InitialLCFileToLoad = '%s/tess2020212050318-s0028-0000000402026209-0190-s_lc.fits'%(path)
# # FileToLoad = '%s/tess2020213081515-s0028-s0028-0000000402026209-00364_dvt.fits'%(path)

# sector = 29
# path = '%s/%s_Sector%d_Exptime120s/MAST_2020-12-02T1335/TESS/tess2020238165205-s0029-0000000402026209-0193-s'%(name,name,sector)
# InitialLCFileToLoad = '%s/tess2020238165205-s0029-0000000402026209-0193-s_lc.fits'%(path)
# FileToLoad = '%s/tess2020239173514-s0029-s0029-0000000402026209-00382_dvt.fits'%(path)

### Exp time 20 s 
# exptime = 20
# sector = 28
# path = '%s/%s_Sector%d_Exptime20s/MAST_2020-12-02T1333/TESS/tess2020212050318-s0028-0000000402026209-0190-a_fast'%(name,name,sector)
# InitialLCFileToLoad = '%s/tess2020212050318-s0028-0000000402026209-0190-a_fast-lc.fits'%(path)
#### FileToLoad = '%s/tess2020213081515-s0028-s0028-0000000402026209-00364_dvt.fits'%(path)

# exptime = 20
# sector = 29
# path = '%s/%s_Sector%d_Exptime20s/MAST_2020-12-02T1334/TESS/tess2020238165205-s0029-0000000402026209-0193-a_fast'%(name,name,sector)
# InitialLCFileToLoad = '%s/tess2020238165205-s0029-0000000402026209-0193-a_fast-lc.fits'%(path)


# name = 'KELT-22Ab'
# sector = 29

# #exptime = 120
# ##path = '%s/%s_Sector%d_Exptime20s/MAST_2020-12-02T1334/TESS/tess2020238165205-s0029-0000000402026209-0193-a_fast'%(name,name,sector)
# # path = 'KELT-22Ab/KELT-22Ab_sector29/MAST_2021-01-20T2036_fast/MAST_2021-01-20T2036/TESS/tess2020238165205-s0029-0000000077031414-0193-a_fast/'
# # InitialLCFileToLoad = '%s/tess2020238165205-s0029-0000000077031414-0193-a_fast-lc.fits'%(path)


# ### What seems to be KELT-22Ab (target 77031414)
# # path = 'KELT-22Ab/KELT-22Ab_sector29/MAST_2021-01-20T2037/MAST_2021-01-20T2037/TESS/tess2020238165205-s0029-0000000077031414-0193-s'
# # InitialLCFileToLoad = '%s/tess2020238165205-s0029-0000000077031414-0193-s_lc.fits'%(path)
# # FileToLoad = '%s/tess2020239173514-s0029-s0029-0000000077031414-00382_dvt.fits'%(path)

# ### Unknown target 77031413 (thought the B companion but it shows transits)
# # path = 'KELT-22Ab/KELT-22Ab_sector29/MAST_2021-01-20T2037_target77031413_LikelyB/MAST_2021-01-20T2037/TESS/tess2020238165205-s0029-0000000077031413-0193-s/'
# # InitialLCFileToLoad = '%s/tess2020238165205-s0029-0000000077031413-0193-s_lc.fits'%(path)
# # FileToLoad = '%s/tess2020239173514-s0029-s0029-0000000077031413-00382_dvt.fits'%(path)

# ## For the 20s exposure time 
# exptime = 20
# path = 'KELT-22Ab/KELT-22Ab_sector29/MAST_2021-01-20T2036_fast/MAST_2021-01-20T2036/TESS/tess2020238165205-s0029-0000000077031414-0193-a_fast/'
# InitialLCFileToLoad = '%s/tess2020238165205-s0029-0000000077031414-0193-a_fast-lc.fits'%(path)
# FileToLoad = '%s/tess2020239173514-s0029-s0029-0000000077031413-00382_dvt.fits'%(path)


# name = 'KELT-11b'
# sector = 9
# exptime = 120
# path = 'C:/Users/ar-h1/research/TESS/KELT-11b/KELT_11b_MAST_2021-02-19T1434/MAST_2021-02-19T1434/TESS/tess2019058134432-s0009-0000000055092869-0139-s'
# InitialLCFileToLoad = '%s/tess2019058134432-s0009-0000000055092869-0139-s_lc.fits'%(path)
# FileToLoad = '%s/tess2019059170935-s0009-s0009-0000000055092869-00198_dvt.fits'%(path)



### For WASP-49b
name = 'WASP-49b'
# sector = 33
# exptime = 120
# path = 'C:/Users/ar-h1/research/TESS/WASP-49b/Sector33/MAST_2021-09-23T1753/TESS/tess2020351194500-s0033-0000000306362738-0203-s'
# InitialLCFileToLoad = '%s/tess2020351194500-s0033-0000000306362738-0203-s_lc.fits'%(path)
# FileToLoad = '%s/tess2020353052510-s0033-s0033-0000000306362738-00430_dvt.fits'%(path)

sector = 33
exptime = 20
path = 'C:/Users/ar-h1/research/TESS/WASP-49b/Sector33_fast/MAST_2021-09-23T1753/TESS/tess2020351194500-s0033-0000000306362738-0203-a_fast'
InitialLCFileToLoad = '%s/tess2020351194500-s0033-0000000306362738-0203-a_fast-lc.fits'%(path)
FileToLoad = '%s/tess2020353052510-s0033-s0033-0000000306362738-00430_dvt.fits'%(path)




InitialLCopenedFile = pyfits.open(InitialLCFileToLoad)

tbdataLC = InitialLCopenedFile[1].data

tLC = tbdataLC['TIME']
bgLC = tbdataLC['SAP_BKG']
lcLC = tbdataLC['SAP_FLUX']
lcLC_Err = tbdataLC['SAP_FLUX_ERR']
pdclcLC = tbdataLC['PDCSAP_FLUX']
PDCSAP_FLUX_ERR = tbdataLC['PDCSAP_FLUX_ERR']
##pdclcLC = ((pdclcLC/18100)-1)*1000  ## Rough conversion to ppt 



LCNanIndices = np.isnan(pdclcLC)

tLC_NoNans =  tLC[~LCNanIndices]
lcLC_NoNans = lcLC[~LCNanIndices]
lcLC_Err_NoNans = lcLC_Err[~LCNanIndices]
pdclcLC_NoNans = pdclcLC[~LCNanIndices]
PDCSAP_FLUX_ERR_NoNans = PDCSAP_FLUX_ERR[~LCNanIndices]

SAPflux = np.zeros((len(tLC_NoNans),3))
SAPflux[:,0] = tLC_NoNans
SAPflux[:,1] = lcLC_NoNans
SAPflux[:,2] = lcLC_Err_NoNans

plt.figure()
plt.errorbar(SAPflux[:,0],SAPflux[:,1],SAPflux[:,2],capsize=3)
plt.title('SAP flux (not PDC correct)')

header_SAP_FLUX = 'Time - 2457000.0 (Barycentric Dynamical Time, TDB), SAP_FLUX (counts), SAP_FLUX_ERR (counts)'

### Line to save: 
np.savetxt('%s/%s_Sector%s_SAP_FLUX_exptime%s.txt'%(name,name,sector,exptime),SAPflux,header=header_SAP_FLUX)


resultsLC = np.zeros((len(tLC_NoNans),3)) ### For the PDC_SAP_FLUX

#resultsLC[:,0] = tLC
resultsLC[:,0] = tLC_NoNans
resultsLC[:,1] = pdclcLC_NoNans
resultsLC[:,2] = PDCSAP_FLUX_ERR_NoNans

header_text_lc = 'Time - 2457000.0 (Barycentric Dynamical Time, TDB), PDCSAP_FLUX (counts), PDCSAP_FLUX_ERR (counts)'

## 'KELT-22Ab_Target77031413'

PDCSAPFLUX_SaveFileName = '%s/%s_Sector%d_PDCSAP_FLUX_exptime%ds.txt'%(name,name,sector,exptime)
##PDCSAPFLUX_SaveFileName = '%s/%s_Sector%d_PDCSAP_FLUX_exptime%ds.txt'%(name,'KELT-22Ab_Target77031414',sector,exptime)

### Line to save: 
np.savetxt(PDCSAPFLUX_SaveFileName,resultsLC,delimiter=',',header=header_text_lc)

plt.figure()
plt.title('PDCSAP_FLUX')
plt.errorbar(resultsLC[:,0],resultsLC[:,1],resultsLC[:,2],capsize=3)

# plt.figure()

# #################

openedFile = pyfits.open(FileToLoad)

h1 = openedFile[1].header
h2 = openedFile[2].header

tb1 = openedFile[1].data
tb2 = openedFile[2].data

detrended_lc1 = tb1['LC_DETREND']*1000  ##ppt 
detrended_Whitenedlc1 = tb1['LC_WHITE']  ##ppt 

LC_INIT = tb1['LC_INIT']*1000 ## column title: Detrended initial light curve (ga

LC_INIT_ERR = tb1['LC_INIT_ERR']*1000 ## column title: Error in the detrended initial li
time1 = tb1['time']

# ###############
nanindices = np.isnan(detrended_lc1)

tNoNans = time1[~nanindices]
fluxNoNans = detrended_lc1[~nanindices]
uncNoNans = LC_INIT_ERR[~nanindices]

# dnoNans = np.zeros((len(fluxNoNans),3))

# dnoNans[:,0] = tNoNans
# dnoNans[:,1] = fluxNoNans
# dnoNans[:,2] = uncNoNans


# ##################

# results = np.zeros((len(time1),3))

# results[:,0] = time1
# results[:,1] = detrended_lc1
# results[:,2] = LC_INIT_ERR

### DVT light curve in results 
results = np.zeros((len(tNoNans),3))

results[:,0] = tNoNans
results[:,1] = fluxNoNans
results[:,2] = uncNoNans

# ##############

dvt_header_text = 'Time - 2457000.0 (Barycentric Dynamical Time, TDB), DVT flux (ppt), DVT flux error (ppt)'

SaveFileName = '%s/%s_Sector%d_dvt_detrended_exptime120s.txt'%(name,name,sector)
###SaveFileName = '%s/%s_Sector%d_dvt_detrended_exptime120s.txt'%(name,'KELT-22Ab_Target77031413',sector)

### Line to save the DVT light curve: 
np.savetxt(SaveFileName,results,delimiter=',',header=dvt_header_text)

plt.figure()
plt.errorbar(results[:,0],results[:,1],results[:,2],capsize=3)
plt.title('DVT light curve')

###################################


# relative_time_error = h1['TIERRELA']

# PDCSAP_FLUX = tb2['PDCSAP_FLUX']
# PDCSAP_FLUX_ERR = tb2['PDCSAP_FLUX_ERR']
# time2 = tb2['TIME']

# PDCSAP_FLUX_TransitsCut_indices = np.where(((PDCSAP_FLUX>6000)&(6075)))
# PDCSAP_FLUX_TransitsCut = PDCSAP_FLUX[PDCSAP_FLUX_TransitsCut_indices]

# Normalization_PDCSAP_FLUX = np.mean(PDCSAP_FLUX_TransitsCut)
# Normalized_PDCSAP_FLUX = PDCSAP_FLUX/Normalization_PDCSAP_FLUX

# Normalized_PDCSAP_FLUX = (Normalized_PDCSAP_FLUX - 1)*1000

# #### Do the subplots of detreneded and not detrended 
# plt.figure()
# plt.errorbar(time1,detrended_lc1,LC_INIT_ERR,capsize=3)
# plt.title('detrended. Errors are called \'Error in the detrended initial li\'\n so not sure if they are applicable')
# plt.xlabel('Time (BJD - 2457000.0)')
# plt.ylabel('Normalized flux (ppt)')

# plt.figure(figsize=(10, 10))

# plt.subplot(2,1,1)
# #plt.title('Sector %d TESS pipeline PDCSAP_FLUX'%(sector))
# #plt.plot(tLC,pdclcLC,'k')
# plt.plot(tLC,lcLC,'k.',markersize=2)
# plt.ylabel('Raw flux (counts)')

# plt.xlim(tLC[0],tLC[-1])
# #plt.ylim((4950,5220))

# plt.subplot(2,1,2)
# #plt.plot(time1,detrended_lc1,'k.',markersize=2)
# #plt.ylabel('Normalized flux (ppt)')

# plt.plot(time1,Normalized_PDCSAP_FLUX,'k.',markersize=2)
# plt.ylabel('Normalized flux (ppt)')



# # plt.plot(time2,PDCSAP_FLUX,'k.',markersize=2)
# # plt.ylabel('Detrended flux (counts)')
# #plt.title('Sector%d_DVT_detrended'%(sector))
# plt.xlabel('Time (BJD - 2457000.0)')

# #plt.xlim(time1[0],time1[-1])
# #plt.savefig('WASP-12b_Sector%s_LC_and_PDCSAP_FLUX.png'%(sector),dpi=400)

# #plt.savefig('Sector26FromDVT_detrended.png',dpi=400)




# #plt.plot(time1,-np.ones_like(time1)*0.014437159587233598*1000)

# #plt.savefig('KELT-22Ab_sector2_dvt_file_detrended_lc.png',dpi=400)

# plt.figure()
# plt.errorbar(time1,LC_INIT,LC_INIT_ERR,capsize=3)
# plt.title('LC_INIT \'Detrended initial light curve (ga\'')
# plt.xlabel('Time (BJD - 2457000.0)')
# plt.ylabel('Normalized flux (ppt)')
# plt.savefig('KELT-22Ab_sector2_dvt_file_LC_INIT.png',dpi=400)


# plt.figure()

# plt.errorbar(time2,PDCSAP_FLUX,PDCSAP_FLUX_ERR,capsize=3)
# plt.title('PDCSAP_FLUX')
# plt.xlabel('Time (BJD - 2457000.0)')
# plt.ylabel('Flux')

# plt.figure()
# plt.plot(time2,PDCSAP_FLUX)
# plt.figure()
# plt.plot(PDCSAP_FLUX_ERR)