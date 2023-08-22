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
rcParams['figure.figsize'] = 13,8.5
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['axes.grid']       = True
rcParams['grid.alpha']      = 0.3
rcParams['axes.grid.which'] = 'major'
#rcParams['font.size'] = 10

def plotz(tfiles,rfiles,name,plotpoints):
   for flt,flr in zip(tfiles,rfiles):
      print(tfiles.index(flt),flt)
      d=pickle.load(open(flt,'rb'))
      dr=pd.read_pickle(flr)
      jumps=len(dr)//plotpoints

      if name in ['b','vp','zp','zm','b_v']:
         names=[name+'r',name+'t',name+'n']

      labels={'zpr':'z^+_r','zpt':'z^+_t','zpn':'z^+_n',
              'zmr':'z^-_r','zmt':'z^-_t','zmn':'z^-_n'}

      plt.clf() 
      ax0=plt.subplot2grid((2,4),(0,0),colspan=4)
      ax1=plt.subplot2grid((2,4),(1,0))
      ax2=plt.subplot2grid((2,4),(1,1))
      ax3=plt.subplot2grid((2,4),(1,2))
      ax4=plt.subplot2grid((2,4),(1,3))
      
      if len(names) == 3:
          for l in range(3):
            ax0.plot(dr.index[::jumps],dr[names[l]][::jumps],'.-',label=r'$'+labels[names[l]]+'$')
      else:
        ax0.plot(dr.index[::jumps],dr[names[0]][::jumps],'.-',label=names[0])
      ax0.set_xlabel('date')
      ax0.set_ylabel(r'$z_{r,t,n}$ (m/s)')
      ax0.set_title('Bin t='+str(dr.index[len(dr.index)//2].round('T'))[:-6])
      ax0.legend(prop={'size':8})

      ax0.set_ylim(-3.1e5,6.1e5)
      ax0.ticklabel_format(axis='y',style='sci',scilimits=(3,5))

      if len(names) == 3:
         ax1.plot(d[name]['times'],d[name]['cr']/3.,label='Cr')
         tcor=d[name]['times'][np.argmin(np.abs(d[name]['cr']/3.-1./np.e))]
      else:
         ax1.plot(d[name]['times'],d[name]['cr'],label='Cr')
         tcor=d[name]['times'][np.argmin(np.abs(d[name]['cr']-1./np.e))]
      ax1.axvline(tcor,linestyle='--',color='r')
      ax1.text(tcor*1.3, 0.8,r'$t_{{corr}}$ = {0:.02f} s'.format(tcor))
      ax1.set_xlabel('$\Delta$t (s)')
      ax1.set_ylabel(r'Cr($\Delta$t)')
      ax1.legend()

      ax2.loglog(d[name]['times'],d[name]['sfn'][0,:],label='S$^{(2)}$')
      yyy=d[name]['sfn'][0,:][58]
#      af.pltpwrl(50,yyy*2.0,xi=1,xf=500,alpha=2./3,label='2/3',ax=ax2)
      yy2=d[name]['sfn'][0,:][1]*1.03
      ax2.text(3,yy2,r'S$^{{(2)}}(dt_{{max}})$={0:.01e} $m^2/s^2$'.format(d[name]['sfn'][0,-1]))
      ax2.set_xlabel('$\Delta$t (s)')
      ax2.set_ylabel(r'S$^{(2)}$($\Delta$t)')
      for tick in ax2.get_yticklabels():
         tick.set_rotation(90)
      ax2.legend()

      if len(names) == 3:
         for l in range(3):
            ax3.semilogx(d[names[l]]['tau'],d[names[l]]['sdk'],label='$\kappa_{'+labels[names[l]]+'}$')
      else:
        ax3.semilogx(d[name]['tau'],d[name]['sdk'],label='$\kappa_{'+names[0]+'}$')
      ax3.set_xlabel('$\Delta$t (s)')
      ax3.set_ylabel(r'$\kappa$($\Delta$t)')
      ax3.legend()

      if len(names) == 3:
         ax4.semilogy(d[names[0]]['bn1'   ],d[names[0]]['pdf1'   ],label=r'$\Delta = 1 dt$')
         ax4.semilogy(d[names[0]]['bn10'  ],d[names[0]]['pdf10'  ],label=r'$\Delta = 10 dt$')
         ax4.semilogy(d[names[0]]['bn100' ],d[names[0]]['pdf100' ],label=r'$\Delta = 100 dt$')
         ax4.semilogy(d[names[0]]['bn1000'],d[names[0]]['pdf1000'],label=r'$\Delta = 1000 dt$')
      else:
         ax4.semilogy(d[name]['bn1'   ],d[name]['pdf1'   ],label=r'$\Delta = 1 dt$')
         ax4.semilogy(d[name]['bn10'  ],d[name]['pdf10'  ],label=r'$\Delta = 10 dt$')
         ax4.semilogy(d[name]['bn100' ],d[name]['pdf100' ],label=r'$\Delta = 100 dt$')
         ax4.semilogy(d[name]['bn1000'],d[name]['pdf1000'],label=r'$\Delta = 1000 dt$')
      ax4.legend()
      ax4.set_xlabel(r'$\Delta '+labels[names[0]]+r'/\sigma_{\Delta '+labels[names[0]]+'}$')
      ax4.set_ylabel('PDF')
      for tick in ax4.get_yticklabels():
         tick.set_rotation(90)
      
   #  plt.tight_layout()
      plt.savefig(name+'-overview.pdf',bbox_inches='tight',pad_inches=0.01)

plotpoints=400
#basedir='/scratch/PSP/DATA/Derived_Products/E1_ENHANCED_DATAFRAME/without_zp_zm_sc/'
basedir='/scratch/PSP/DATA/Derived_Products/E1'
RAW_DIR=basedir+   '/8hr_chunks/'
TRB_DIR=basedir+'/tq_8hr_chunks/'
tfiles=sorted(glob.glob(TRB_DIR+'e1_8hr_07.p'))
rfiles=sorted(glob.glob(TRB_DIR+'df_e1_8hr_07.p'))

allnames=['zp','zm']
for name in allnames:
   plotz(tfiles,rfiles,name,plotpoints)
