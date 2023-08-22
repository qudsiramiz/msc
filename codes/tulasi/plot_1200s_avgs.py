#!/usr/bin/env python
import pandas as pd
from pylab import *
rcParams['axes.grid']       = True
rcParams['grid.alpha']      = 0.3
rcParams['axes.grid.which'] = 'major'
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
#chunks=[2,4,8,16,24]
#marks =['P','*','o','s','^']
#colors=['r','g','b','c','m']
chunks=[4,8,24]
marks =['o','s','^']
colors=['b','c','m']
basedir='/scratch/PSP/DATA/Derived_Products/mag_products/MAG_CHUNKS_SWEAP_CADENCE'
df=pd.read_pickle('/scratch/PSP/DATA/Derived_Products/psp_swe_fld_merged_swe_cadence_20181031_20181111_FirstPerihelion_Hampel.p')
dfr=df.resample('0.874S').median()
den=dfr.np_moment.rolling(1429,min_periods=1).mean()
di = 228./sqrt(den)
vpr=dfr.vpr.rolling(1429,min_periods=1).mean()
bmag=sqrt(dfr.br**2+dfr.bt**2+dfr.bn**2)
bm=bmag.rolling(1429,min_periods=1).mean()
va=bm*1e-9/sqrt(4*np.pi*1e-7*1.67*1e-27*den*1e6)/1e3
d={}
for i in chunks:
   d[i] = pd.read_pickle(basedir+'/tq_'+str(i)+'hr_chunks/tq_bin_avg_'+str(i)+'hr.p')

f,ax=subplots(5,1,sharex='col')
f.set_figheight(10)
f.set_figwidth(6)
#f.suptitle('1250s Averages')

for i in chunks:
   ax[0].plot(d[i].t,d[i].tcor,markersize=7,marker=marks[chunks.index(i)],linewidth=0.3,label=str(i)+'hr bins',alpha=0.1,color=colors[chunks.index(i)])
   ax[0].plot(d[i].t,d[i].tcor.rolling(24//i).mean(),label=str(i)+'hr @24hr avg',color=colors[chunks.index(i)])
ax[0].set_ylim(0,1200)
ax[0].set_ylabel('$t_{corr}$ (s)')
ax[0].legend(ncol=3,prop={'size':8})

ax[1].plot(dfr.index,dfr.np_moment,alpha=0.1,label='$n_p$')
ax[1].plot(dfr.index,den,label='$n_p$ 1250s avg')
ax[1].set_xlim(datetime.datetime(2018,11,3),datetime.datetime(2018,11,10))
ax[1].set_ylim(0,800)
plt.setp(ax[1].get_yticklabels()[-1], visible=False)
ax[1].legend(ncol=2)
ax[1].set_ylabel(r'$n_p$ (cm$^{-3}$)' )

ax[2].plot(dfr.index,di,label='$d_i$ from 1250s avg $n_p$')
ax[2].set_xlim(datetime.datetime(2018,11,3),datetime.datetime(2018,11,10))
ax[2].legend()
ax[2].set_ylabel('$d_i$ (Km)')

ax[3].plot(dfr.index,dfr.vpr,alpha=0.1,label='$V_{sw}$')
ax[3].plot(dfr.index,vpr,label='$V_{sw}$ 1250s avg')
ax[3].plot(dfr.index,va ,label='$V_a$ from 1250s avgs')
ax[3].set_xlim(datetime.datetime(2018,11,3),datetime.datetime(2018,11,10))
ax[3].set_ylabel('$V_{sw}$, $V_a$ (km/s)')
ax[3].legend(ncol=3)

ax[3].set_xlim(datetime.datetime(2018,11,3),datetime.datetime(2018,11,10))

ax[4].plot(dfr.index,vpr/va,label='ratio of 1250s avgs')
ax[4].legend()
ax[4].set_ylim(0,6)
ax[4].set_ylabel('$M_A$')
ax[4].set_xlabel('date')
plt.setp(ax[4].get_yticklabels()[-1], visible=False)
plt.xticks(rotation=45)
f.subplots_adjust(hspace=0.02)
f.align_ylabels([i for i in ax])
#savefig('1250_rolling.png',dpi=300,bbox_inches='tight')
import os
savefig('1250_rolling.pdf',bbox_inches='tight',pad_inches=0.01)
#plt.show()
