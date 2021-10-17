# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 23:34:47 2020

@author: ar-h1
"""

import numpy as np
import matplotlib.pyplot as plt 
from brokenaxes import brokenaxes
from matplotlib.gridspec import GridSpec

s19 = np.loadtxt('XO6b_Sector19_dvt_file_detrended_lc.txt',delimiter=',')
s20 = np.loadtxt('XO6b_Sector20_dvt_file_detrended_lc.txt',delimiter=',')
s26 = np.loadtxt('XO6b_Sector26_dvt_file_detrended_lc.txt',delimiter=',')


s19LC = np.loadtxt('TESS_LC_XO6b_Sector19.txt',delimiter=',')
s20LC = np.loadtxt('TESS_LC_XO6b_Sector20.txt',delimiter=',')
s26LC = np.loadtxt('TESS_LC_XO6b_Sector26.txt',delimiter=',')

s19PDCSAPFLUX = np.loadtxt('TESS_PDCSAP_FLUX_XO-6b_Sector19.txt',delimiter=',')
s20PDCSAPFLUX = np.loadtxt('TESS_PDCSAP_FLUX_XO-6b_Sector20.txt',delimiter=',')
s26PDCSAPFLUX = np.loadtxt('TESS_PDCSAP_FLUX_XO-6b_Sector26.txt',delimiter=',')

s19LC_shifted = np.copy(s19LC)
s20LC_shifted = np.copy(s20LC)
s26LC_shifted = np.copy(s26LC)

s19LC_shifted[:,1] = s19LC[:,1] - np.nanmedian(s19LC[:,1]) + 1000 
s20LC_shifted[:,1] = s20LC[:,1] - np.nanmedian(s20LC[:,1]) + 1000 
s26LC_shifted[:,1] = s26LC[:,1] - np.nanmedian(s26LC[:,1]) + 1000 


d = np.vstack((s19,s20,s26))
dLC = np.vstack((s19LC_shifted,s20LC_shifted,s26LC_shifted))

d_PDCSAPFLUX = np.vstack((s19PDCSAPFLUX,s20PDCSAPFLUX,s26PDCSAPFLUX))

OnlyPDCSAPFLUXTime = d_PDCSAPFLUX[:,0]
OnlyPDCSAPFLUX = d_PDCSAPFLUX[:,1]
OnlyPDCSAPFLUXErr = d_PDCSAPFLUX[:,2]

pdcsapflux_nanind = np.isnan(OnlyPDCSAPFLUX)

OnlyPDCSAPFLUXTime_NoNans = OnlyPDCSAPFLUXTime[~pdcsapflux_nanind]
OnlyPDCSAPFLUX_NoNans = OnlyPDCSAPFLUX[~pdcsapflux_nanind]
OnlyPDCSAPFLUXErr_NoNans = OnlyPDCSAPFLUXErr[~pdcsapflux_nanind]

AllPDCSAPFLUX_NoNans = np.empty((len(OnlyPDCSAPFLUXTime_NoNans),3))
AllPDCSAPFLUX_NoNans[:,0] = OnlyPDCSAPFLUXTime_NoNans
AllPDCSAPFLUX_NoNans[:,1] = OnlyPDCSAPFLUX_NoNans
AllPDCSAPFLUX_NoNans[:,2] = OnlyPDCSAPFLUXErr_NoNans

FileHeader = 'time - 2457000.0 (Barycentric Dynamical Time, TDB), PDCSAP_FLUX (counts), PDCSAP_FLUX_ERR'
np.savetxt('XO-6b_PDCSAP_FLUX_file_S19_20_26.txt',AllPDCSAPFLUX_NoNans,delimiter=',',header=FileHeader)


plt.figure()
plt.title('All sectors of PDCSAP_FLUX')
plt.errorbar(AllPDCSAPFLUX_NoNans[:,0],AllPDCSAPFLUX_NoNans[:,1],AllPDCSAPFLUX_NoNans[:,2],capsize=10)

plt.figure()

dLCflux = dLC[:,1]
idxs = np.where((dLCflux<600) | (dLCflux>1200))

dLCflux[idxs] = np.nan
dLCflux[idxs] = np.nan

plt.plot(dLC[:,0],dLCflux,'k.',markersize=1)
#plt.ylim((600,1200))

nanindices = np.isnan(d[:,1])

t = d[:,0]
flux = d[:,1]
unc = d[:,2]

tNoNans = t[~nanindices]
fluxNoNans = flux[~nanindices]
uncNoNans = unc[~nanindices]

dnoNans = np.zeros((len(fluxNoNans),3))

dnoNans[:,0] = tNoNans
dnoNans[:,1] = fluxNoNans
dnoNans[:,2] = uncNoNans

#plt.plot(d[:,0],d[:,1],'k.',markersize=1)
plt.plot(dnoNans[:,0],dnoNans[:,1],'k.',markersize=1)


fig = plt.figure(figsize=(6,2))
#bax = brokenaxes(xlims=((1815, 1870), (2009, 2036)), ylims=((-20,10)))
bax = brokenaxes(xlims=((1814, 1870), (2009, 2036)))

bax.plot(d[:,0], d[:,1], 'k.', markersize=1)
bax.set_xlabel('Time (BJD - 2457000.0)')
bax.set_ylabel('Normalized flux (ppt)')

plt.savefig('XO-6b_dvt_S19_20_26.png',dpi=400)

FileHeader = 'time - 2457000.0 (Barycentric Dynamical Time, TDB), detrended flux (ppt), error (ppt)'
#np.savetxt('XO-6b_dvt_file_S19_20_26.txt',dnoNans,delimiter=',',header=FileHeader)
   

##################

fig = plt.figure(figsize=(10,10))

sps1, sps2 = GridSpec(2,1)

bax = brokenaxes(xlims=((1814, 1870), (2009, 2036)),subplot_spec=sps1)
bax.plot(dLC[:,0], dLCflux, 'k.', markersize=1)
bax.set_ylabel('Shifted raw flux (counts)')

bax = brokenaxes(xlims=((1814, 1870), (2009, 2036)),subplot_spec=sps2)
bax.plot(d[:,0], d[:,1], 'k.', markersize=1)
bax.set_xlabel('Time (BJD - 2457000.0)')
bax.set_ylabel('Normalized flux (ppt)')

#plt.savefig('XO-6b_LC_and_DVT_S19_20_26.png',dpi=400)










