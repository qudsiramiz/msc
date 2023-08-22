#!/usr/bin/env python3
import os
import glob
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import TurbAn.Analysis.Simulations.AnalysisFunctions as af

from pylab import rcParams
rcParams['figure.figsize'] = 4.83,6.67
rcParams['xtick.direction'] = 'in'
rcParams['font.size'] = 10

plotpoints=400
basedir='/scratch/PSP/DATA/Derived_Products/E1/'
RAW_DIR=basedir+   '/8hr_chunks/'
TRB_DIR=basedir+'/tq_8hr_chunks/'
tfiles=sorted(glob.glob(TRB_DIR+'e1_8hr_08.p'))
rfiles=sorted(glob.glob(TRB_DIR+'df_e1_8hr_08.p'))
noise_seconds=7

f,ax = plt.subplots(4,1,sharex='col')
for flt,flr in zip(tfiles,rfiles):
   print(tfiles.index(flt),flt)
   d=pickle.load(open(flt,'rb'))
   dr=pd.read_pickle(flr)
   k  =d['zp']['k']*1e3 # 1/m -> 1/km
   ekp=d['zpr']['ek']+d['zpt']['ek']+d['zpn']['ek']
   ekm=d['zmr']['ek']+d['zmt']['ek']+d['zmn']['ek']
   ekv=d['vp']['ek']
   ekb=d['bv']['ek']
   di = 228e3/np.sqrt(d['par']['np'])

   sc = (ekp - ekm)/(ekp + ekm)
   sr = (ekv - ekb)/(ekv + ekb)
   ct = sc/np.sqrt(1 - sr**2)
   ra = ekv/ekb
   xx = np.logspace(-6.,-4.,100)

   ax[0].cla()
   ax[0].plot(k,sc)
   ax[0].set_ylabel(r'$\sigma_c$',size='x-large')
   ax[0].axvline(1/di/100,linestyle='--',color='g')#,label=r'$(10^2 d_i)^{-1}$')
   ax[0].axvline(1/di/1e3,linestyle='--',color='m')#,label=r'$(10^3 d_i)^{-1}$')
   ax[0].axvline(2*np.pi*1000/(noise_seconds*d['par']['vsw']),linestyle=':',label='SPC Noise Floor?')
#  ax[0].semilogx(xx,0.21-0.1*np.log10(xx),'-.',label=r'$\sigma_c \sim \log_{10}(f^{-0.1})$')
   ax[0].semilogx(ax[0].get_xlim(),[0.75,0.75],'-.',label=r'$\sigma_c = 0.75$')
   ax[0].legend(loc='lower left')
   ax[0].text(7e-4,0.36,r'${1/100  d_i}$')
   ax[0].text(7e-5,0.36,r'${1/1000 d_i}$')
   ax[0].set_title('Bin t='+str(dr.index[len(dr.index)//2].round('T'))[:-6])

   ax[1].cla()
   ax[1].plot(k,sr)
   ax[1].set_ylabel(r'$\sigma_r$',size='x-large')
   ax[1].axvline(1/di/100,linestyle='--',color='g')
   ax[1].axvline(1/di/1e3,linestyle='--',color='m')
   ax[1].axvline(2*np.pi*1000/(noise_seconds*d['par']['vsw']),linestyle=':')

   ax[2].cla()
   ax[2].plot(k,ct)
   ax[2].set_ylabel(r'$\cos\ \theta_{vb}$',size='x-large')
   ax[2].axvline(1/di/100,linestyle='--',color='g')
   ax[2].axvline(1/di/1e3,linestyle='--',color='m')
   ax[2].axvline(2*np.pi*1000/(noise_seconds*d['par']['vsw']),linestyle=':')

   ax[3].cla()
   ax[3].semilogx(k,ra)
   ax[3].set_xlabel(r'lag(1/km)')
   ax[3].set_ylabel(r'$r_A$',size='x-large')
   ax[3].set_xlim(3e-7,4e-3)
   ax[3].axvline(1/di/100,linestyle='--',color='g')
   ax[3].axvline(1/di/1e3,linestyle='--',color='m')
   ax[3].axvline(2*np.pi*1000/(noise_seconds*d['par']['vsw']),linestyle=':')
   
f.align_ylabels([i for i in ax])
plt.subplots_adjust(hspace=0.03,wspace=0)
plt.savefig('sc_sr_ct_ra_scale08.pdf',bbox_inches='tight', pad_inches=0.01)
