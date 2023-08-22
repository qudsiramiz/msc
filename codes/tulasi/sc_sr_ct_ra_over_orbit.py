#!/usr/bin/env python3
import os
import glob
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#import TurbAn.Analysis.Simulations.AnalysisFunctions as af

from pylab import rcParams
rcParams['figure.figsize'] = 4.83,6.67
rcParams['xtick.direction'] = 'in'
rcParams['font.size'] = 10
#rcParams['font.family'] = 'FreeMono'

plotpoints=400
#basedir='/scratch/PSP/DATA/Derived_Products/E1_ENHANCED_DATAFRAME/without_zp_zm_sc/'
basedir='/scratch/PSP/DATA/Derived_Products/E1/'
RAW_DIR=basedir+   '8hr_chunks/'
TRB_DIR=basedir+'tq_8hr_chunks/'
tfiles=sorted(glob.glob(TRB_DIR+'e1_8hr_*.p'))
rfiles=sorted(glob.glob(TRB_DIR+'df_e1_8hr_*.p'))

lags=[10000,3000,1000]
linestyles=['-','--']
colors=['r','k']

outd={}
outd['t']=[]
for i in ['sc','sr','ct','ra']:
   outd[i]={}
   for j in range(len(lags)):
      outd[i][lags[j]]=[]

for flt,flr in zip(tfiles,rfiles):
   print(tfiles.index(flt),flt)
   d=pickle.load(open(flt,'rb'))
   dr=pd.read_pickle(flr)
   k  =d['zpr']['k']*1e3 # 1/m -> 1/km
   ekp=d['zpr']['ek']+d['zpt']['ek']+d['zpn']['ek']
   ekm=d['zmr']['ek']+d['zmt']['ek']+d['zmn']['ek']
   ekv=d['vp']['ek']
   ekb=d['bv']['ek']
   di = 228e3/np.sqrt(d['par']['np'])

   ind={}
   ind['sc'] = (ekp - ekm)/(ekp + ekm)
   ind['sr'] = (ekv - ekb)/(ekv + ekb)
   ind['ct'] = ind['sc']/np.sqrt(1 - ind['sr']**2)
   ind['ra'] = ekv/ekb

   outd['t'].append(dr.index[len(dr.index)//2])
   for i in ['sc','sr','ct','ra']:
      for j in lags:
         idx = np.argmin(np.abs(k-1/di/j))
         outd[i][j].append(ind[i][idx])


f,ax = plt.subplots(4,1,sharex='col')
for i in range(len(lags) -1 ):
   ax[0].plot(outd['t'],outd['sc'][lags[i]],linestyle=linestyles[i],color=colors[i])
   ax[1].plot(outd['t'],outd['sr'][lags[i]],linestyle=linestyles[i],color=colors[i],label=r'$\Delta={0:03d}d_i$'.format(lags[i]))
   ax[2].plot(outd['t'],outd['ct'][lags[i]],linestyle=linestyles[i],color=colors[i])
   ax[3].plot(outd['t'],outd['ra'][lags[i]],linestyle=linestyles[i],color=colors[i])

ax[1].legend(ncol=2)#,loc='upper center',bbox_to_anchor=[0.5,1.26])#,prop={'size':8})
ax[0].set_ylabel(r'$\sigma_c$')
ax[1].set_ylabel(r'$\sigma_r$')
ax[2].set_ylabel(r'$\cos\theta$')
ax[3].set_ylabel(r'$r_A$')


ax[3].set_xlabel('date')

f.align_ylabels([i for i in ax])
plt.xticks(rotation=23)
plt.tight_layout()
plt.subplots_adjust(hspace=0.03,wspace=0)
f.canvas.draw()
labels = [item.get_text() for item in ax[3].get_xticklabels()]
labels[0] = ''; labels[-1] = ''; labels[2]=''; labels[4]=''
ax[3].set_xticklabels(labels)
plt.savefig('sc_sr_ct_ra_orbit.pdf',bbox_inches='tight',pad_inches=0.01)

