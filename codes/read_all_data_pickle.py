from spacepy.pycdf import CDF as cdf

from numpy import *
from glob import glob
import pandas as pd
import datetime
import pytz
import h5py as hf

#
read_files = input( 'Read the files? ==> ' )

if read_files :

    fnames = sort( glob('/home/ahmadr/SpacePlasFest/msc/data/*.p') )

    ldf = []
    msc_n = []

    for f in fnames :

        ldf.append( pd.read_pickle( f ) )
        msc_n.append( f[36:-46] )

plt.close( 'all' )

#plt.figure( )

marker = [ 'd', 'o', '^', 'v', 'D', '*' ]
c = [ 'r', 'b', 'g', 'm', 'k', 'c' ]

for xx in range( len( ldf ) ) :

    if( xx == 0 ) :
        fig = plt.figure( num=None, figsize=( 7, 15 ), dpi=200, facecolor='w', 
                           edgecolor='gray' )

        fig.subplots_adjust( left=0.15, right=0.95, top=0.98, bottom=0.10,
                             wspace=0.02, hspace=0. )

        # Define the axes in the figure

        axs1 = fig.add_subplot( 2, 1, 1, sharex=None )

    axs1.scatter( ldf[xx].sc_r_median, ldf[xx].np_median, c=c[xx], marker=marker[xx], s=5, label=msc_n[xx] )

#	axs1.set_xlim( 0, 100 )
    axs1.set_ylim( 1e-3, 1e2 )

    axs1.set_yscale( 'log' )

    axs1.legend( loc=1, fontsize=18 )

    axs1.set_ylabel( r'$n_p (/cm^{-3})$', fontsize=18 )

    if( xx == 0 ) :

        axs2 = fig.add_subplot( 2, 1, 2, sharex=axs1 )

    axs2.scatter( ldf[xx].sc_r_median, ldf[xx].vp_median, c=c[xx], marker=marker[xx], s=5, label=msc_n[xx] )

    axs2.set_xlim( 0, 100 )
    axs2.set_ylim( 7e1, 1e3 )

    axs2.set_yscale( 'log' )

    axs2.legend( loc=3, fontsize=18 )

    axs2.set_ylabel( r'$v_p (km/s)$', fontsize=18 )

    axs2.set_xlabel( r'$R_s (AU)$', fontsize=18 )

plt.savefig( 'np_vp_all_sc_median.png' )

# Plot magnetic Fields

for xx in range( len( ldf ) ) :

    if( xx == 0 ) :
        fig = plt.figure( num=None, figsize=( 7, 15 ), dpi=200, facecolor='w', 
                           edgecolor='gray' )

        fig.subplots_adjust( left=0.15, right=0.95, top=0.98, bottom=0.10,
                             wspace=0.02, hspace=0. )

        # Define the axes in the figure

        axs1 = fig.add_subplot( 4, 1, 1, sharex=None )

    axs1.scatter( ldf[xx].sc_r_median, ldf[xx].br_median, c=c[xx], marker=marker[xx], s=5, label=msc_n[xx] )

#	axs1.set_xlim( 0, 100 )
    axs1.set_ylim( -5, 5 )

#	axs1.set_yscale( 'log' )

    axs1.legend( loc=1, fontsize=18 )

    axs1.set_ylabel( r'$B_R (nT)$', fontsize=18 )

    if( xx == 0 ) :

        axs2 = fig.add_subplot( 4, 1, 2, sharex=axs1 )

    axs2.scatter( ldf[xx].sc_r_median, ldf[xx].br_median, c=c[xx], marker=marker[xx], s=5, label=msc_n[xx] )

#	axs2.set_xlim( 0, 100 )
    axs2.set_ylim( -5, 5 )

#	axs2.set_yscale( 'log' )

    axs2.legend( loc=3, fontsize=18 )

    axs2.set_ylabel( r'$B_T (nT)$', fontsize=18 )

    if( xx == 0 ) :

        axs3 = fig.add_subplot( 4, 1, 3, sharex=axs1 )

    axs3.scatter( ldf[xx].sc_r_median, ldf[xx].bn_median, c=c[xx], marker=marker[xx], s=5, label=msc_n[xx] )

#	axs3.set_xlim( 0, 100 )
    axs3.set_ylim( -5, 5 )

#	axs3.set_yscale( 'log' )

    axs3.legend( loc=3, fontsize=18 )

    axs3.set_ylabel( r'$B_N (nT)$', fontsize=18 )

    if( xx == 0 ) :

        axs4 = fig.add_subplot( 4, 1, 4, sharex=axs1 )

    axs4.scatter( ldf[xx].sc_r_median, ldf[xx].bm_median, c=c[xx], marker=marker[xx], s=5, label=msc_n[xx] )

    axs4.set_xlim( 0, 100 )
    axs4.set_ylim( 1e-2, 2e1 )

    axs4.set_yscale( 'log' )

    axs4.legend( loc=3, fontsize=18 )

    axs4.set_ylabel( r'$B_M (nT)$', fontsize=18 )

    axs4.set_xlabel( r'$R_s (AU)$', fontsize=18 )

plt.savefig( 'b_mag_all_sc_median.png' )
