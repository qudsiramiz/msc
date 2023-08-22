#!/usr/bin/env python3
import os
import glob
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import TurbAn.Analysis.Simulations.AnalysisFunctions as af
from smooth import SmoothySpec

from pylab import rcParams
rcParams['figure.figsize'] = 4.83,6.67
rcParams['xtick.direction'] = 'in'
rcParams['font.size'] = 10
rcParams['axes.grid']       = True
rcParams['grid.alpha']      = 0.3
rcParams['axes.grid.which'] = 'major'
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

plotpoints=400
basedir='/scratch/PSP/DATA/Derived_Products/E1/'
RAW_DIR=basedir+   '/8hr_chunks/'
TRB_DIR=basedir+'/tq_8hr_chunks/'
tfiles=sorted(glob.glob(TRB_DIR+'e1_8hr_08.p'))
rfiles=sorted(glob.glob(TRB_DIR+'df_e1_8hr_08.p'))

f,ax = plt.subplots(4,1,sharex='col')
for flt,flr in zip(tfiles,rfiles):
   print(tfiles.index(flt),flt)
   d=pickle.load(open(flt,'rb'))
   dr=pd.read_pickle(flr)
   fq =d['zp']['fq']
   ekp=d['zp']['ef']
   ekm=d['zm']['ef']
   ekv=d['vp']['ef']
   ekb=d['bv']['ef']

   di = 228e3/np.sqrt(d['par']['np'])

   sc = (ekp - ekm)/(ekp + ekm)
   sr = (ekv - ekb)/(ekv + ekb)
   ct = sc/np.sqrt(1 - sr**2)
   ra = ekv/ekb
   xx = np.logspace(-3,-1,100)

   ax[0].cla()
   ax[0].plot(fq[1:],SmoothySpec(sc[1:]))
   ax[0].set_ylabel(r'$\sigma_c$',size='x-large')
   ax[0].semilogx(xx,0.40-0.12*np.log10(xx),':',label=r'$\sigma_c \sim \log_{10}({f^{-0.12}})$')
   ax[0].semilogx(ax[0].get_xlim(),[0.75,0.75],':',label=r'$\sigma_c = 0.75$')
   ax[0].axvline(0.06,linestyle='--',label='SPC Noise Floor?')
   ax[0].set_ylim(0.15,0.82)
   ax[0].legend()
   ax[0].set_title('Bin t='+str(dr.index[len(dr.index)//2].round('T'))[:-6])

   ax[1].cla()
   ax[1].plot(fq[1:],SmoothySpec(sr[1:]))
   ax[1].set_ylabel(r'$\sigma_r$',size='x-large')
#  ax[1].semilogx(xx,0.5+0.30*np.log10(xx),':',label=r'$\sigma_r \sim \log_{10}({f^{0.30}})$')
   ax[1].semilogx(xx,0.3+0.20*np.log10(xx),':',label=r'$\sigma_r \sim \log_{10}({f^{0.20}})$')
   ax[1].axvline(0.06,linestyle='--')
   ax[1].set_ylim(-0.45,0.60)
   ax[1].legend()

   ax[2].cla()
   ax[2].plot(fq[1:],SmoothySpec(ct[1:]))
   ax[2].semilogx(xx,0.40-0.12*np.log10(xx),':',label=r'$\cos \theta_{vb} \sim \log_{10}({f^{-0.12}})$')
   ax[2].semilogx(ax[2].get_xlim(),[0.75,0.75],':',label=r'$\cos \theta_{vb} = 0.75$')
   ax[2].set_ylabel(r'$\cos \theta_{vb}$',size='x-large')
   ax[2].axvline(0.06,linestyle='--')
   ax[2].set_ylim(0.15,0.85)
   ax[2].legend()

   ax[3].cla()
   ax[3].semilogx(fq[1:],SmoothySpec(ra[1:]))
   ax[3].axvline(0.06,linestyle='--')
#  ax[3].semilogx(xx,2.30+0.55*np.log10(xx),':',label=r'$r_A \sim \log_{10}({f^{0.55}})$')
   ax[3].semilogx(xx,1.60+0.35*np.log10(xx),':',label=r'$r_A \sim \log_{10}({f^{0.35}})$')
   ax[3].set_ylim(0.4,3.3)
   ax[3].legend()
   ax[3].set_xlabel(r'f(Hz)')
   ax[3].set_ylabel(r'$r_A$',size='x-large')
   ax[3].set_xlim(7.5e-4,0.55)
f.align_ylabels([i for i in ax])
plt.subplots_adjust(hspace=0.03,wspace=0)
plt.savefig(os.environ['HOME']+'/WorkSpace/PSP/ApJ/plots/sc_sr_ct_ra_fft08.pdf',bbox_inches='tight', pad_inches=0.01)
