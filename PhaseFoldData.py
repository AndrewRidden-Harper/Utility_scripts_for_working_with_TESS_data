# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 13:13:16 2020

@author: ar-h1
"""


import numpy as np
import matplotlib.pyplot as plt 
import lightkurve as lk
from scipy.stats import binned_statistic


def RoughBinnedUnc(vect):
    
    meanunc = np.mean(vect)
    
    NumPoints = len(vect)
    
    return meanunc/np.sqrt(NumPoints)




# LoadedData = np.loadtxt('WASP-12b_Sector20_dvt_file_detrended_lc.txt',delimiter=',')
# DetrendMethodStr = 'DVT'

# LoadedData = np.loadtxt('WASP-12b_Sector20_PDCSAP_FLUX.txt',delimiter=',')
# DetrendMethodStr = 'PDCSAP_FLUX'


# LoadedData = np.loadtxt('WASP-12b_TESS_sector20_WorkshopCodeNewErrors.txt',delimiter=',')
# DetrendMethodStr = 'WorkshopCode'


# planet = 'WASP-4b'

# P = 1.33823152840
# T0 = 2458355.1848698 + P/2


# LoadedData = np.loadtxt('WASP-12b_TESS_sector20_WorkshopCodeNewErrors.txt',delimiter=',')
# DetrendMethodStr = 'WorkshopCode'


# planet = 'WASP-4b'

# P = 1.33823152840
# T0 = 2458355.1848698 + P/2


#LoadedData = np.loadtxt('%s/%s_Sector2_dvt_detrended_exptime120s.txt'%(planet,planet),delimiter=',')

# LoadedDataP1 = np.loadtxt('%s/%s_Sector28_dvt_detrended_exptime120s.txt'%(planet,planet),delimiter=',')
# LoadedDataP2 = np.loadtxt('%s/%s_Sector29_dvt_detrended_exptime120s.txt'%(planet,planet),delimiter=',')

#LoadedData = np.loadtxt('%s/%s_Sector2_PDCSAP_FLUX_exptime120s.txt'%(planet,planet),delimiter=',')

planet = 'KELT-22Ab'

P = 1.3866529043360933
T0 = 2457288.8596280282363296 + P/2
#T0 = 2457288.8585 + P/2
	

LoadedData = np.loadtxt('KELT-22Ab_sector2_dvt_file_detrended_lc.txt',delimiter=',')

#LoadedData = np.loadtxt('KELT-22Ab/KELT-22Ab_Target77031414_Sector29_dvt_detrended_exptime120s.txt',delimiter=',')
#LoadedData = np.loadtxt('KELT-22Ab/KELT-22Ab_Target77031414_Sector29_PDCSAP_FLUX_exptime120s.txt',delimiter=',')




# LoadedDataP1 = np.loadtxt('%s/%s_Sector28_PDCSAP_FLUX_exptime120s.txt'%(planet,planet),delimiter=',')
# LoadedDataP2 = np.loadtxt('%s/%s_Sector29_PDCSAP_FLUX_exptime120s.txt'%(planet,planet),delimiter=',')

# LoadedData = np.vstack((LoadedDataP1,LoadedDataP2))





DetrendMethodStr = 'DVT'
#DetrendMethodStr = 'PDCSAP'

SaveDir = '%s/PhaseFolded'%(planet)

Sector = '2'
#Sector = '28_and_29'
#Sector = '29'

ExpTime = 120

############ BAD error 
# LoadedData = np.loadtxt('WASP-12b_TESS_sector20_WorkshopCode.txt',delimiter=',')
# DetrendMethodStr = 'WorkshopCode'

#LoadedData[:,2] = LoadedData[:,2]*((LoadedData[:,1]+1000)/1000)
#LoadedData[:,2] = LoadedData[:,2]*((LoadedData[:,1]+1000))


Descriptor = 'AllData'
#Descriptor = 'FirstHalfOfData'
#Descriptor = 'LastHalfOfData'


NumBins = 150
#NumBins = 50


#t = LoadedData[:,0] - np.min(LoadedData[:,0]) #+ 2457000.0
t = LoadedData[:,0] + 2457000.0
f = LoadedData[:,1]
e = LoadedData[:,2]

#MiddleIndex = int(len(t)/2)
MiddleIndex = 8795

if Descriptor == 'FirstHalfOfData':
    t = t[:MiddleIndex]
    f = f[:MiddleIndex]
    e = e[:MiddleIndex]
    
if Descriptor == 'LastHalfOfData':
    t = t[MiddleIndex:]
    f = f[MiddleIndex:]
    e = e[MiddleIndex:] 


plt.plot(t,f,color='grey',lw=1)

lc = lk.lightcurve.LightCurve(time=t, flux=f,flux_err=e)
lc_err = lk.lightcurve.LightCurve(time=t, flux=e)
#lc_err = lk.lightcurve.LightCurve(time=t, flux=(LoadedData[:,2]*LoadedData[:,1]))


folded_lc = lc.fold(period=P, t0=T0)

folded_lc.plot()

flcp = folded_lc.phase
flcf = folded_lc.flux 

ferror = lc_err.fold(period=P, t0=T0)
fep = ferror.phase
fef = ferror.flux

#flce = folded_lc.error 

binnedflux,binedges,bin_number = binned_statistic(flcp,flcf,'mean',bins=NumBins)

roughbinnedunc, errbinedges,errbin_number = binned_statistic(fep,fef,RoughBinnedUnc,bins=NumBins)


#PlotBinValues = binedges[:-1] + np.diff(binedges)
PlotBinValues = (binedges[0:-1] + binedges[1:])/2

plt.errorbar(PlotBinValues,binnedflux,yerr=roughbinnedunc)

output = np.empty((len(flcp),3))
binnedoutput = np.empty((len(binnedflux),3))

if DetrendMethodStr == 'DVT':
    outputheadertext1 = 'phase, flux (ppt), error (ppt)'
    outputheadertext2 = 'phase, mean binned flux (ppt), rough binned error (ppt)'
    
if DetrendMethodStr == 'PDCSAP':
    outputheadertext1 = 'phase, flux, error'
    outputheadertext2 = 'phase, mean binned flux, rough binned error'


output[:,0] = flcp
output[:,1] = flcf
output[:,2] = fef

np.savetxt('%s/%s_TESS_ExpTime%d_Sector%s_%s_PhaseFolded_OccultationAtZero.txt'%(SaveDir,planet,ExpTime,Sector,DetrendMethodStr),output,delimiter=',',header=outputheadertext1)

outputheadertext = 'phase, mean binned flux (ppt), rough binned error (ppt)'
binnedoutput[:,0] = PlotBinValues
binnedoutput[:,1] = binnedflux
binnedoutput[:,2] = roughbinnedunc

np.savetxt('%s/%s_TESS_ExpTime%ds_Sector%s_%s_%d_Bins_PhaseFolded_OccultationAtZero.txt'%(SaveDir,planet,ExpTime,Sector,DetrendMethodStr,NumBins),binnedoutput,delimiter=',',header=outputheadertext2)

plt.figure()
plt.errorbar(PlotBinValues,binnedflux,yerr=roughbinnedunc)
#plt.ylim((-1.5,1.5)) 
plt.xlabel('phase')
plt.ylabel('binnedflux (ppt)')
#plt.savefig('WASP-12b_binned_secondary.png',dpi=400)

plt.figure()
plt.plot(flcp,flcf)

plt.xlabel('phase')
plt.ylabel('flux (ppt)')




