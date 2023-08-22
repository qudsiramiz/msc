# Import and configure the necessary modules.

from numpy import *

import h5py

import matplotlib
import matplotlib.colors
import matplotlib.lines
import matplotlib.patches
import matplotlib.pyplot

matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{starfont}'

# Define the function for generating a plot.

def plot_indiv( key=None, fl=None, ylim=None, yscale=None, ylabel=None, keycorner=1 ) :

    # Define the location of the plot within the figure as well as the size of each.

    inch_plot_x      = 5.00
    inch_plot_y      = 3.00

    inch_hist_x      = inch_plot_x
    inch_hist_y      = 1.00

    inch_mrgn_left   = 0.65
    inch_mrgn_right  = 0.05
    inch_mrgn_bottom = 0.45
    inch_mrgn_top    = 0.20

    inch_fig_x = ( inch_mrgn_left   + inch_plot_x + inch_mrgn_right )
    inch_fig_y = ( inch_mrgn_bottom + inch_plot_y + inch_hist_y + inch_mrgn_top   )

    figsize = [ inch_fig_x, inch_fig_y ]

    pos_plot = [ inch_mrgn_left   / inch_fig_x,
                 inch_mrgn_bottom / inch_fig_y,
                 inch_plot_x      / inch_fig_x,
                 inch_plot_y      / inch_fig_y  ]

    pos_hist = [ inch_mrgn_left   / inch_fig_x,
                 ( inch_mrgn_bottom + inch_plot_y ) / inch_fig_y,
                 inch_hist_x      / inch_fig_x,
                 inch_hist_y      / inch_fig_y  ]

    # Define the settings for the plot.

    xlim   = [ 0.01 , 150. ]
    xscale = 'log'
    xlabel = r'\rm Distance from Sun, $r$ [au]'

    hlim   = [ 3.E1, 3.E5 ]
    hscale = 'log'
    hlabel = r'\rm Intervals'

    # Initialize the figure.

    fig  = matplotlib.pyplot.figure( figsize=figsize )

    matplotlib.rc( 'text', usetex=True )

    # Initialize the plot.

    fig_plot = fig.add_subplot( 211, position=pos_plot )

    fig_plot.set_xlim( xlim )
    fig_plot.set_ylim( ylim )

    fig_plot.set_xscale( xscale )
    fig_plot.set_yscale( yscale )

    fig_plot.set_xlabel( xlabel, size='small' )
    fig_plot.set_ylabel( ylabel, size='small' )

    # Initialize the histogram.

    fig_hist = fig.add_subplot( 212, position=pos_hist )

    fig_hist.set_xlim( xlim )
    fig_hist.set_ylim( hlim )

    fig_hist.set_xscale( xscale )
    fig_hist.set_yscale( hscale )

    fig_hist.set_ylabel( hlabel, size='small' )

    fig_hist.set_xticklabels( [] )

    # Initialize the key.

    if ( keycorner == 2 ) :
        fig_key = fig.add_axes( [ 0.13, 0.57, 0.50, 0.15 ] )
    else :
        fig_key = fig.add_axes( [ 0.13, 0.11, 0.50, 0.15 ] )

    fig_key.set_xlim( [ 0., 3. ] )
    fig_key.set_ylim( [ 0., 3. ] )

    fig_key.set_xticks( [] )
    fig_key.set_yticks( [] )

    # Generate the solar-system object markers.

    d_mercury =  0.3871
    d_venus   =  0.7233
    d_earth   =  1.000
    d_mars    =  1.524
    d_ceres   =  2.767
    d_jupiter =  5.204
    d_saturn  =  9.583
    d_uranus  = 19.22
    d_neptune = 30.11
    d_pluto   = 39.48
    d_eris    = 67.781

    lst_d = [ d_mercury, d_venus, d_earth, d_mars, d_ceres,
              d_jupiter, d_saturn, d_uranus, d_neptune, d_pluto, d_eris ]

    lst_lab = [ r'\starfontserif \Mercury', r'\starfontserif \Venus'  ,
                r'\starfontserif \Terra'  , r'\starfontserif \Mars'   ,
                r'\starfontserif \Ceres'  , r'\starfontserif \Jupiter',
                r'\starfontserif \Saturn' , r'\starfontserif \Uranus' ,
                r'\starfontserif \Neptune', r'\starfontserif \Pluto'  ,
                r'\starfontserif \Mars'                                 ]

    y = ( ( inch_mrgn_bottom + inch_plot_y + inch_hist_y ) / inch_fig_y ) + 0.005

    for i in range( len( lst_d ) ) :

        d = lst_d[i]

        if ( ( d == d_ceres ) or ( d == d_pluto ) or ( d == d_eris ) ) :
            c = 'gray'
        else :
            c = 'black'

        if ( d == d_eris ) :
            shift_x  =  -1.3
            shift_y  =  0.000
            rotation = 221
        else :
            shift_x  = 0
            shift_y  = 0
            rotation = 0

        fig_hist.annotate( lst_lab[i],
                           ( d + shift_x, y + shift_y ),
                           xycoords=( 'data', 'figure fraction' ),
                           color=c, size='medium',
                           horizontalalignment='center',
                           verticalalignment='bottom',
                           rotation=rotation                       )

        fig_hist.plot( ( d, d ), hlim,
                       color=c, linewidth=0.5, linestyle=':' )

        fig_plot.plot( ( d, d ), ylim,
                       color=c, linewidth=0.5, linestyle=':' )

    # Load and plot the data from each spacecraft

    n_sc = 10

    sc_name = [ 'PSP'         , 'Helios 1'  , 'Helios 2'  ,
                'Cassini'     , 'Pioneer 10', 'Pioneer 11',
                'New Horizons', 'Voyager 1' , 'Voyager 2', 'Ulysses' ]

    sc_color = [ 'tab:brown', 'tab:orange', 'tab:olive',
                 'tab:green', 'tab:purple', 'tab:pink',
                 'tab:red'  , 'tab:blue'  , 'tab:cyan','tab:cyan'    ]

    sc_marker = [ 'X', '^', 'v', 'o', '<', '>', '*', 's', 'd', 'D' ]

    p1 = '/mnt/cephadrius/udel_research/msc/data/binned/scaled/v16/'
    p2 = '_coho1hr_merged_mag_plasma_'
    p3 = '_v06_80_binned_scaled_80_binned_v16.hf'

    sc_fl = [ p1 + 'psp'       + p2 + '20180101_20210701' + p3,
              p1 + 'helios1'      + p2 + '19740101_19811201' + p3,
              p1 + 'helios2'      + p2 + '19760101_19801201' + p3,
              p1 + 'cassini'      + p2 + '20000101_20040101' + p3,
              p1 + 'pioneer10'    + p2 + '19720101_19950901' + p3,
              p1 + 'pioneer11'    + p2 + '19730101_19941201' + p3,
              p1 + 'new_horizons' + p2 + '20081010210700' + p3,
              p1 + 'voyager1'     + p2 + '19770101_20181201' + p3,
              p1 + 'voyager2'     + p2 + '19770101_20181201' + p3,
              p1 + 'uy'      + p2 + '19900101_19920201' + p3 ]

    for sc in range( n_sc ) :

        # Load the data.

        #print(sc_fl[sc])

        dat = h5py.File( sc_fl[sc], 'r' )
        #print(list(dat.keys()))
        dat_r = array( dat['sc_r_iqr_50'] )
        dat_k = array( dat[key+'_iqr_50'] )
        #dat_n = array( dat[key+'_num']    )

        tk = where( ( isfinite( dat_r ) ) & ( isfinite( dat_k ) ) )[0]
#		            ( isfinite( dat_n ) )                           )[0]

        # Populate the plot.

        fig_plot.plot( dat_r[tk], dat_k[tk],
                       markersize=3.0, markeredgewidth=0,
                       color=sc_color[sc], marker=sc_marker[sc],
                       linestyle=''                              )

        # Populate the histogram.

        #fig_hist.plot( dat_r[tk], dat_n[tk],
        #               markersize=3.0, markeredgewidth=0,
        #               color=sc_color[sc], marker=sc_marker[sc],
        #               linestyle=''                              )

    # Populate the key.

    for sc in range( n_sc ) :

        x_sc = sc / 3
        y_sc = sc % 3

        fig_key.plot( [ 0.9*x_sc + 0.2 ], [ y_sc + 0.5 ],
                      markersize=3.0, markeredgewidth=0,
                      color=sc_color[sc], marker=sc_marker[sc],
                      linestyle=''                              )

        fig_key.annotate( sc_name[sc],
                          ( 0.9*x_sc + 0.3, y_sc + 0.5 ),
                          xycoords='data',
                          color='k', size='small',
                          horizontalalignment='left',
                          verticalalignment='center'  )

    # Output the figure.

    fname = '/mnt/cephadrius/udel_research/msc/figures/indiv/indiv_' + fl + '_v16.png'

    matplotlib.pyplot.savefig( fname )


# Generate the plots.

plot_indiv( key='bm', fl='b',
            ylim=[ 0.01, 200. ], yscale='log',
            ylabel=r'\rm Magnetic Field Strength, $B$ [nT]' )

plot_indiv( key='np', fl='np',
            ylim=[ 0.0005, 500. ], yscale='log',
            ylabel=r'\rm Proton Density, $n_p$ [${\rm cm}^{-3}$]' )

plot_indiv( key='vp_m', fl='vp',
            ylim=[ 0., 750. ], yscale='linear',
            ylabel=r'\rm Proton Speed, $v_p$ [${\rm km}/{\rm s}$]' )

plot_indiv( key='vp_r', fl='vpr',
            ylim=[ 0., 750. ], yscale='linear',
            ylabel=r'\rm Proton Radial Speed, $v_{pr}$ [${\rm km}/{\rm s}$]' )

plot_indiv( key='Tp', fl='tp',
            ylim=[ 1.E3, 2.E6 ], yscale='log',
            ylabel=r'\rm Proton Temperature, $T_p$ [K]' )

plot_indiv( key='vA', fl='va', keycorner=2,
            ylim=[ 0., 3.E5 ], yscale='linear',
            ylabel=r"\rm Alf\'en Speed, $v_A$ [${\rm m}/{\rm s}$]" )

plot_indiv( key='alfven_ratio', fl='na',
            ylim=[ 0., 40. ], yscale='linear',
            ylabel=r"\rm Alf\'en Number, $N_A = ???$" )

plot_indiv( key='parker_angle', fl='ang',
            ylim=[ 0., 3.14159/2. ], yscale='linear',
            ylabel=r"\rm Parker Angle, ??? [???]" )

plot_indiv( key='particle_flux', fl='flux',
            ylim=[ 1.E34, 1.E36 ], yscale='log',
            ylabel=r"\rm Proton Flux, $\Phi_p = n_p V_p r^2$ [/(s $\Omega$)]" )

plot_indiv( key='sig_c', fl='sigma',
                ylim=[ -1., 1. ], yscale='linear',
                ylabel=r"$\sigma_c$ [???]" )

plot_indiv( key='sc_r', fl='heliographicLatitude', keycorner=2,
            ylim=[ -90, 90 ], yscale='linear',
            ylabel=r"Heliographic Latitude (degrees)" )

plot_indiv( key='sc_r', fl='heliographicLongitude', keycorner=2,
            ylim=[ -90, 90 ], yscale='linear',
            ylabel=r"Heliographic Longitude (degrees)" )

plt.close('all')

