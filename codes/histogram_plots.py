def brzl_hist( x=None, y=None, z=[], xlim=None, ylim=None, ztype='', spc='', 
              date='', bins=50, n_min=5, dens=False, vmin=None, vmax=None, mode=None) :

    # Check to make sure all required values are passed. If not, abort.
    # Note: 'z = None' produces standard bins of number of data

    if ( x is None ) :

        raise ValueError( 'x data not provided' )
        return

    if ( y is None ) :

        raise ValueError( 'y data not provided' )
        return

    # Change x, y, and z into numpy arrays

    x = np.array( x )
    y = np.array( y )
    z = np.array( z )

    # Determine if 1 or 2 bin values are given and separate values
    # if necessary.

    if ( hasattr( bins, 'len' ) ) :

        nx = bins[0]
        ny = bins[1]

    else :

        nx = bins
        ny = bins

    if( xlim is None ) :

        xlim = [ np.nanmin( x ), np.nanmax( x ) ]

    if( ylim is None ) :

        ylim = [ np.nanmin( y ), np.nanmax( y ) ]

    # If no z data have been provided, generate a standard histogram.
    # Otherwise, manually generate bins and bin the data.

    if len(z) == 0 :

        # Note: bins variables are lists of all leading edges plus
        # the lagging edge of the final bin.

        xbins = np.linspace( xlim[0], xlim[1], nx+1 )
        ybins = np.linspace( ylim[0], ylim[1], ny+1 )

        hst = plt.hist2d( x, y, bins=bins, cmin=n_min, range=[ xlim, ylim ] )

        plt.close('all')
        # number of data points
        N = np.nansum( np.nansum( hst[0], axis=0 ) )

        zplot1 = np.transpose(hst[0])
    else :

        # set number of data points to zero for running tally
        N = 0

        # Note: bins variables are lists of all leading edges plus
        # the lagging edge of the final bin.

        xbins = np.linspace( xlim[0], xlim[1], nx+1 )
        ybins = np.linspace( ylim[0], ylim[1], ny+1 )

        # For each bin, find all the z data whose corresponding
        # x and y values fall within that bin and save the medians
        # as an array to be plotted.

        zplot1 = np.array( [ [ float('nan') for n in range( nx ) ]
                                        for m in range( ny ) ] )

        tkx = []
        tky = []
        tkz = []

        # For each bin along the x-axis...

        for i in range( nx ) :

            # find the x data with beta values within that bin

            tkx = np.where( ( x >= xbins[i]   ) &
                     ( x <  xbins[i+1] )   )[0]

            # If no valid beta values were found,
            # move to the next bin.

            if len( tkx ) == 0 :

                continue

            # For each bin along the y-axis...

            for j in range( ny ) :

                # find the y data with anisotropy values
                # within that bin.

                tkxy = np.where( ( y[tkx] >= ybins[j]   ) &
                              ( y[tkx] <  ybins[j+1] )   )[0]

                # If any valid data were found, take the
                # median of those data (ignoring 'Nan's) and
                # assign that value to the corresponding
                # location in 'zplot' for plotting.

                if len( tkxy ) >= n_min :
                    if mode == 'median' :
                        zplot1[i,j] = np.nanmedian( z[tkx][tkxy])
                    elif mode == 'mean' :
                        zplot1[i,j] = np.nanmean( z[tkx][tkxy])
                    # running tally of number of data points
                    N += len( np.where( z[tkx][tkxy] != float('nan') )[0] )

        # Generate the histogram (not actually a matplotlib histogram)
        # Note: imshow() plots according to ( y, x )

        zplot1 = np.transpose( zplot1 )

    return xbins, ybins, zplot1

def brz_plotting_routine(data=None, ex=None, ey=None, xmin=None, xmax=None, ymin=None, ymax=None,
                         bins=None, nbins=None, cmin=None, cmax=None, norm=None, cmap=None,
                         vmin=None, vmax=None, spc=None, save_fig=False, plot_show=False, mode=None) :

    pad = 0.02
    clabelpad = 5.0
    labelsize = 30
    ticklabelsize = 30
    cticklabelsize = 20
    clabelsize = 25
    ticklength = 10
    mticklength = 4
    tickwidth = 1.5
    tick_width = 1.5
    tick_length = 10
    ctick_length = 5
    mcticklength = 4
    aspect = 'auto'
    labelrotation = 0
    alpha = 0.9

    if None in [xmin, xmax, ymin, ymax] :
        xmin = min(ex)
        xmax = max(ex)
        ymin = min(ey)
        ymax = max(ey)

    if bins is None :
        if (nbins is None):
            nbins = 50
        bins = (np.logspace(np.log10(xmin), np.log10(xmax),nbins),
                np.logspace(np.log10(ymin), np.log10(ymax),nbins))

    if cmap is None :
        cmap = plt.cm.Blues

    if cmin is None :
        cmin = 5
    if cmax is None :
        cmax = 1e2

    if norm is None :
        norm = mpl.colors.LogNorm(vmin, vmax)

    fig = plt.figure(num=None, figsize=( 24, 7.5), dpi=300, facecolor='w', edgecolor='gray')
    fig.subplots_adjust(left=0.01, right=0.95, top=0.99, bottom=0.01, wspace=0.0, hspace=0.0)

    axs_xlabel = r'$\beta_{\parallel \mathrm{p}} = 2 \mu_0 n_\mathrm{p} k_\mathrm{B}$' +\
                 r'$T_{\parallel \mathrm{p}}/B^2$'

    axs_ylabel = r'$R_{\mathrm{p}} = T_{\perp \mathrm{p}}/T_{\parallel \mathrm{p}}$'

    cbar_label = r'$\tilde{p}( \beta_{\parallel \mathrm{p}}, R_\mathrm{p}) = n/( N \Delta \beta_{\parallel \mathrm{p}} \Delta R_\mathrm{p} )$'

    # Figure 1
    axs1 = fig.add_subplot(1, 4, 1, sharex=None)

    im1 = axs1.pcolormesh(ex, ey, data[0], alpha=alpha, shading='auto', cmap=cmap, norm=norm)

    axs1.axhline( 0, lw=1.0, c='k' )

    axs1.set_xlabel( axs_xlabel, fontsize=labelsize )
    axs1.set_ylabel( axs_ylabel, fontsize=labelsize )

    axs1.set_xlim(xmin, xmax)
    axs1.set_ylim(ymin, ymax)
    axs1.set_xscale('log')
    axs1.set_yscale('log')

    im1.axes.tick_params(which='major', axis='both', direction='in', labelbottom=True, bottom=True,
                         labeltop=False, top=True, labelleft=True, left=True, labelright=False,
                         right=True, width=tickwidth, length=ticklength, labelsize=ticklabelsize,
                         labelrotation=labelrotation)
    im1.axes.tick_params(which='minor', axis='both', direction='in', labelbottom=False, bottom=True,
                         labeltop=False, top=True, labelleft=False, left=True, labelright=False,
                         right=True, width=tickwidth, length=mticklength,
                         labelsize=ticklabelsize, labelrotation=labelrotation)

    divider1 = make_axes_locatable(axs1)
    cax1 = divider1.append_axes("top", size="5%", pad=pad )
    cbar1 = plt.colorbar(im1, cax=cax1, orientation='horizontal', ticks=None, fraction=0.05,
                         pad=0.0)

    cbar1.ax.tick_params(axis='x', which='both', direction='in', labeltop=True, top=True, 
                         labelbottom=False, bottom=False, width=tick_width, length=ctick_length, 
                         labelsize=clabelsize, labelrotation=0, pad=0)

    cbar1.ax.xaxis.set_label_position('top')

    cbar1.set_label( r'$N$', fontsize=clabelsize, labelpad=clabelpad )


    # Figure 2
    axs2 = fig.add_subplot(1, 4, 2, sharex=None)

    im2 = axs2.pcolormesh(ex, ey, data[1], alpha=alpha, shading='auto',
                          norm=mpl.colors.LogNorm(8e-5, 6e-1), cmap=cmap)

    axs2.axhline( 0, lw=1.0, c='k' )

    axs2.set_xlabel( axs_xlabel, fontsize=labelsize )

    axs2.set_xlim(xmin, xmax)
    axs2.set_ylim(ymin, ymax)
    axs2.set_xscale('log')
    axs2.set_yscale('log')

    im2.axes.tick_params(which='major', axis='both', direction='in', labelbottom=True, bottom=True,
                         labeltop=False, top=True, labelleft=False, left=True, labelright=False,
                         right=True, width=tickwidth, length=ticklength, labelsize=ticklabelsize,
                         labelrotation=labelrotation)
    im2.axes.tick_params(which='minor', axis='both', direction='in', labelbottom=False, bottom=True,
                         labeltop=False, top=True, labelleft=False, left=True, labelright=False,
                         right=True, width=tickwidth, length=mticklength,
                         labelsize=ticklabelsize, labelrotation=labelrotation)

    #fmtx = FuncFormatter( lambda x, pos: tickformatx( x ) )
    #axs2.xaxis.set_major_formatter( fmtx )
    #axs2.yaxis.set_major_formatter( fmtx )

    divider2 = make_axes_locatable(axs2)
    cax2 = divider2.append_axes("top", size="5%", pad=pad )
    cbar2 = plt.colorbar(im2, cax=cax2, orientation='horizontal', ticks=None, fraction=0.05,
                         pad=0.0)

    cbar2.ax.tick_params(axis='x', which='both', direction='in', labeltop=True, top=True, 
                         labelbottom=False, bottom=False, width=tick_width, length=ctick_length,
                         labelsize=clabelsize, labelrotation=0, pad=0)

    cbar2.ax.xaxis.set_label_position('top')
    cbar2.set_label( r'$\Gamma_{\max}/\Omega_{\rm cp}$', fontsize=clabelsize, labelpad=clabelpad )

    # Figure 3
    axs3 = fig.add_subplot(1, 4, 3, sharex=axs1)

    im3 = axs3.pcolormesh(ex, ey, data[2], alpha=alpha, shading='auto',
                          norm=mpl.colors.LogNorm(8e-5, 6e-1), cmap=cmap)

    axs3.axhline( 0, lw=1.0, c='k' )

    axs3.set_xlabel( axs_xlabel, fontsize=labelsize )

    axs3.set_xlim(xmin, xmax)
    axs3.set_ylim(ymin, ymax)
    axs3.set_xscale('log')
    axs3.set_yscale('log')

    im3.axes.tick_params(which='major', axis='both', direction='in', labelbottom=True, bottom=True,
                         labeltop=False, top=True, labelleft=False, left=True, labelright=False,
                         right=True, width=tickwidth, length=ticklength, labelsize=ticklabelsize,
                         labelrotation=labelrotation)
    im3.axes.tick_params(which='minor', axis='both', direction='in', labelbottom=False, bottom=True,
                         labeltop=False, top=True, labelleft=False, left=True, labelright=False,
                         right=True, width=tickwidth, length=mticklength,
                         labelsize=ticklabelsize, labelrotation=labelrotation)

    #fmtx = FuncFormatter( lambda x, pos: tickformatx( x ) )
    #axs3.xaxis.set_major_formatter( fmtx )
    #axs3.yaxis.set_major_formatter( fmtx )

    divider3 = make_axes_locatable(axs3)
    cax3 = divider3.append_axes("top", size="5%", pad=pad )
    cbar3 = plt.colorbar(im3, cax=cax3, orientation='horizontal', ticks=None, fraction=0.05,
                         pad=0.0)

    cbar3.ax.tick_params(axis='x', which='both', direction='in', labeltop=True, top=True, 
                         labelbottom=False, bottom=False, width=tick_width, length=ctick_length, 
                         labelsize=clabelsize, labelrotation=0, pad=0)

    cbar3.ax.xaxis.set_label_position('top')

    cbar3.set_label( r'$\omega_{\rm nl}/\Omega_{\rm cp}$', fontsize=clabelsize, labelpad=clabelpad )


    # Figure 4
    axs4 = fig.add_subplot(1, 4, 4, sharex=None)

    im4 = axs4.pcolormesh(ex, ey, data[3], alpha=alpha, shading='auto',
                          norm=mpl.colors.LogNorm(1e-2, 1e2), cmap=plt.cm.bwr)

    axs4.axhline( 0, lw=1.0, c='k' )

    axs4.set_xlabel( axs_xlabel, fontsize=labelsize )

    axs4.set_xlim(xmin, xmax)
    axs4.set_ylim(ymin, ymax)
    axs4.set_xscale('log')
    axs4.set_yscale('log')

    im4.axes.tick_params(which='major', axis='both', direction='in', labelbottom=True, bottom=True,
                         labeltop=False, top=True, labelleft=False, left=True, labelright=True,
                         right=True, width=tickwidth, length=ticklength, labelsize=ticklabelsize,
                         labelrotation=labelrotation)
    im4.axes.tick_params(which='minor', axis='both', direction='in', labelbottom=False, bottom=True,
                         labeltop=False, top=True, labelleft=False, left=True, labelright=False,
                         right=True, width=tickwidth, length=mticklength,
                         labelsize=ticklabelsize, labelrotation=labelrotation)

    #fmtx = FuncFormatter( lambda x, pos: tickformatx( x ) )
    #axs4.xaxis.set_major_formatter( fmtx )
    #axs4.yaxis.set_major_formatter( fmtx )

    divider4 = make_axes_locatable(axs4)
    cax4 = divider4.append_axes("top", size="5%", pad=pad )
    cbar4 = plt.colorbar(im4, cax=cax4, orientation='horizontal', ticks=None, fraction=0.05,
                         pad=0.0)

    cbar4.ax.tick_params(axis='x', which='both', direction='in', labeltop=True, top=True, 
                         labelbottom=False, bottom=False, width=tick_width, length=ctick_length, 
                         labelsize=clabelsize, labelrotation=0, pad=0)

    cbar4.ax.xaxis.set_label_position('top')

    cbar4.set_label( r'$\Gamma_{\max}/\omega_{\rm nl}$', fontsize=clabelsize, labelpad=clabelpad )

    xticks2 = im2.axes.xaxis.get_majorticklabels()
    xticks2[1].set_visible(False)
    xticks2[-2].set_visible(False)

    xticks3 = im3.axes.xaxis.get_majorticklabels()
    xticks3[-2].set_visible(False)

    #yticks1 = axs1.xaxis.get_minor_ticks()
    #yticks1[-1].set_visible(False)

    #xticks2 = axs2.xaxis.get_major_ticks()
    #xticks2[0].set_visible(False)

    #xticks3 = axs3.xaxis.get_major_ticks()
    #xticks3[0].set_visible(False)

    #xticks4 = axs4.xaxis.get_major_ticks()
    #xticks4[0].set_visible(False)

    #yticks4 = axs1.xaxis.get_major_ticks()
    #yticks4[-1].set_visible(False)

    #cxticks2 = cax2.xaxis.get_major_ticks()
    #cxticks2[0].set_visible(False)

    cxticks3 = cax3.xaxis.get_major_ticks()
    cxticks3[-1].set_visible(False)

    #cxticks4 = cax4.xaxis.get_major_ticks()
    #cxticks4[0].set_visible(False)

    if (save_fig) :
        plt.savefig(f'../figures/{spc}_brz_omega_gamma_{mode}_v000.pdf', bbox_inches='tight',
                    pad_inches = 0.05, format='pdf', dpi=300)

        print(f'Figure saved for spacecraft {spc}')

    if (plot_show==True) :
        plt.show()
    #plt.close('all')