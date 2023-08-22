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
rcParams['figure.figsize'] = 10.15,5.5
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['font.size'] = 10

plotpoints=400
basedir='/scratch/PSP/DATA/Derived_Products/E1/'
RAW_DIR=basedir+'tq_8hr_chunks/'
rfiles=sorted(glob.glob(RAW_DIR+'df_e1_8hr_07.p'))
#!~ RAW_DIR=basedir+'8hr_chunks/'
#!~ rfiles=sorted(glob.glob(RAW_DIR+'e1_8hr_08.p'))

f,ax = plt.subplots(4,1,sharex='col')
for flr in rfiles:
   print(rfiles.index(flr),flr)
   df=pd.read_pickle(flr)
   npts=len(df)//plotpoints
   plt.clf() 
   ax3=plt.subplot2grid((4,7),(3,0),colspan=6            )
   ax0=plt.subplot2grid((4,7),(0,0),colspan=6, sharex=ax3)
   ax1=plt.subplot2grid((4,7),(1,0),colspan=6, sharex=ax3)
   ax2=plt.subplot2grid((4,7),(2,0),colspan=6, sharex=ax3)

   hx0=plt.subplot2grid((4,7),(0,6), sharey=ax0)
   hx1=plt.subplot2grid((4,7),(1,6), sharey=ax1)
   hx2=plt.subplot2grid((4,7),(2,6), sharey=ax2)
   hx3=plt.subplot2grid((4,7),(3,6), sharey=ax3)

   sc = df.sc
   print('Computing ekv,ekb,sr,ct,ra')
   vmg      = np.sqrt(df.vpr**2+df.vpt**2+df.vpn**2)
   bmg      = np.sqrt(df.bvr**2+df.bvt**2+df.bvn**2)
   df['ekv']= 0.5*vmg**2
   df['ekb']= 0.5*bmg**2
   df['ssc']= (df.vpr*df.bvr+df.vpt*df.bvt+df.vpn*df.bvn)/(df.ekv+df.ekb)
   df['sr' ]= (df.ekv - df.ekb)/(df.ekv + df.ekb)
   df['ct' ]= df.ssc/np.sqrt(1 - df.sr**2)
   df['cct']= (df.vpr*df.bvr+df.vpt*df.bvt+df.vpn*df.bvn)/(vmg*bmg)
   df['ra' ]= df.ekv/df.ekb

#  ax0.plot(k,sc)
   print('Plotting sc ')
   df.ssc[::npts].plot(ax=ax0,label='1 NYs cadence',legend=True)
   tmp=df.ssc.rolling(1429,min_periods=1).mean()
   tmp[::npts].plot(ax=ax0,legend=True,label='1250s rolling')
   df.ssc.hist(ax=hx0,bins=32,histtype=u'step',density=False,orientation='horizontal',\
      label=df[['ssc']].describe(exclude='dtype').to_string())
   ax0.text(df.index[6000],-0.4,r'$\langle \sigma_c \rangle$ = {0:0.3f}'.format(\
   df.ssc.mean()),bbox=dict(facecolor='w',edgecolor='r',boxstyle='round,pad=1',alpha=0.55))
   ax0.set_ylabel(r'$\sigma_c$',size='x-large')
   ax0.set_title('Bin t='+str(df.index[len(df.index)//2].round('T'))[:-6])
   ax0.legend(loc='lower right')

   print('Plotting sr ')
   df.sr[::npts].plot(ax=ax1,legend=False)
   tmp=df.sr.rolling(1429,min_periods=1).mean()
   tmp[::npts].plot(ax=ax1,legend=False,label='1250s rolling')
   df.sr.hist(ax=hx1,bins=32,histtype=u'step',density=False,orientation='horizontal',\
      label=df[['sr']].describe(exclude='dtype').to_string())
   ax1.text(df.index[6000],0.2,r'$\langle \sigma_r \rangle$ = {0:0.3f}'.format(\
   df.sr.mean()),bbox=dict(facecolor='w',edgecolor='r',boxstyle='round,pad=1',alpha=0.55))
   ax1.set_ylabel(r'$\sigma_r$',size='x-large')
#  ax1.legend(loc='upper right')

   try:
      print('Plotting ct ')
      df.cct[::npts].plot(ax=ax2,legend=False)
      tmp=df.cct.rolling(1429,min_periods=1).mean()
      tmp[::npts].plot(ax=ax2,legend=False,label='1250s rolling')
      df.cct.hist(ax=hx2,bins=32,histtype=u'step',density=False,orientation='horizontal',\
         label=df[['cct']].describe(exclude='dtype').to_string())
      ax2.text(df.index[6000],-0.4,r'$\langle \cos \theta \rangle$ = {0:0.3f}'.format(\
      df.cct.mean()),bbox=dict(facecolor='w',edgecolor='r',boxstyle='round,pad=1',alpha=0.55))
      ax2.set_ylabel(r'$\cos \theta_{vb}$',size='x-large')
#     ax2.legend(loc='lower right')
   except:
      pass

   print('Plotting ra ')
   ra=df.ra.clip(upper=1.5)#,inplace=True)
   ra[::npts].plot(ax=ax3,legend=False)
   tmp=ra.rolling(1429,min_periods=1).mean()
   tmp[::npts].plot(ax=ax3,legend=False,label='1250s rolling')
   ra.hist(ax=hx3,bins=32,histtype=u'step',density=False,orientation='horizontal',\
      label=df[['ra']].describe(exclude='dtype').to_string())
   ax3.text(df.index[6000],.9,r'$\langle r_A \rangle$ = {0:0.3f}'.format(\
   df.ra.mean()),bbox=dict(facecolor='w',edgecolor='r',boxstyle='round,pad=1',alpha=0.55))
   ax3.set_ylim(-0.1,1.5)
   #ax3.set_yscale('log')
   ax3.set_ylabel(r'$r_A$',size='x-large')
   ax3.set_xlabel(r'datetime',size='x-large')
#  ax3.legend(loc='upper right')

   for i in [hx0,hx1,hx2,hx3]:
      i.grid(False)

#  hx0.set_title('density')
f.align_ylabels([ax0,ax1,ax2,ax3])
plt.savefig(os.environ['HOME']+'/WorkSpace/PSP/ApJ/plots/sc_sr_ct_ra_hist07.pdf',bbox_inches='tight',pad_inches=0.01)
