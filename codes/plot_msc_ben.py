# Import and configure the necessary modules.
import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Suppress warnings
import warnings
warnings.filterwarnings('ignore')

# matplotlib.rcParams['text.latex.preamble'] = '\usepackage{starfont} \usepackage{amsmath}'
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{starfont}'
# matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
# set latex use to false
matplotlib.rc('text', usetex=True)

# Define the model functions for fitting:
try:
    plt.close("all")
except Exception:
    pass


def mod_power(x, c, a):
    return c * (x**a)


def mod_broken(x, c, x0, a1, a2, force=0):

    y = []

    for xx in x:

        if (((force == 0) and (xx <= x0)) or (force == 1)):
            y.append(c * xx**a1)
        else:
            y.append(c * (x0**(a1 - a2)) * (xx**a2))

    return y


def mod_broken2(x, c, xa, xb, a1, a2, a3, force=0):

    y = []

    for xx in x:

        if (((force == 0) and (xx <= xa)) or (force == 1)):
            y.append(c * xx**a1)
        elif (((force == 0) and (xx <= xb)) or
                (force == 2)):
            y.append(c * (xa**(a1 - a2)) * (xx**a2))
        else:
            y.append(c * (xa**(a1 - a2)) * (xb**(a2 - a3)) * (xx**a3))

    return y


# Define the function for generating a plot.

def plot_msc(key, scaled=True):

    # Define the location of the plot within the figure as well as the size of each.

    inch_plot_x = 5.00
    inch_plot_y = 3.00

    inch_hist_x = inch_plot_x
    inch_hist_y = 1.00

    inch_mrgn_left = 0.65
    inch_mrgn_right = 0.05
    inch_mrgn_bottom = 0.45
    inch_mrgn_top = 0.20

    inch_fig_x = (inch_mrgn_left + inch_plot_x + inch_mrgn_right)
    inch_fig_y = (inch_mrgn_bottom + 2. * inch_plot_y + inch_hist_y + inch_mrgn_top)

    figsize = [inch_fig_x, inch_fig_y]

    pos_comp = [inch_mrgn_left / inch_fig_x,
                inch_mrgn_bottom / inch_fig_y,
                inch_plot_x / inch_fig_x,
                inch_plot_y / inch_fig_y]

    pos_indiv = [inch_mrgn_left / inch_fig_x,
                 (inch_mrgn_bottom + inch_plot_y) / inch_fig_y,
                 inch_plot_x / inch_fig_x,
                 inch_plot_y / inch_fig_y]

    pos_hist = [inch_mrgn_left / inch_fig_x,
                (inch_mrgn_bottom + 2. * inch_plot_y) / inch_fig_y,
                inch_hist_x / inch_fig_x,
                inch_hist_y / inch_fig_y]

    # Generate the plot settings specific for this key.

    if (key == 'b'):
        dkey = 'bm'
        lab_name = r'Magnetic Field Strength'
        lab_sym = r'B'
        lab_unit = r'{\rm nT}'
        corner = 3
        yscale = 'log'
        ylim_scaled = [0.003, 300.]
        ylim_unscaled = [0.008, 8000.]
        unit = 1.
    elif (key == 'np'):
        dkey = 'np'
        lab_name = r'Proton Density'
        lab_sym = r'n_{\rm p}'
        lab_unit = r'{\rm cm}^{-3}'
        corner = 3
        yscale = 'log'
        ylim_scaled = [0.00008, 5000.]
        ylim_unscaled = [0.0005, 20000.]
        unit = 1.
    elif (key == 'vp'):
        dkey = 'vp_m'
        lab_name = r'Proton Speed'
        lab_sym = r'v_{\rm p}'
        lab_unit = r'{\rm km}/{\rm s}'
        corner = 3
        yscale = 'linear'
        ylim_scaled = [0.0, 1.5]
        ylim_unscaled = [1., 850.]
        unit = 1.
    elif (key == 'vpr'):
        dkey = 'vp_r'
        lab_name = r'Proton Radial Speed'
        lab_sym = r'v_{{\rm p}r}'
        lab_unit = r'{\rm km}/{\rm s}'
        corner = 3
        yscale = 'linear'
        ylim_scaled = [0.0, 1.5]
        ylim_unscaled = [1., 850.]
        unit = 1.
    elif (key == 'tp'):
        dkey = 'Tp'
        lab_name = r'Proton Temperature'
        lab_sym = r'T_{\rm p}'
        lab_unit = r'{\rm K}'
        corner = 1
        yscale = 'log'
        ylim_scaled = [0.01, 100.]
        ylim_unscaled = [8.E2, 2.E6]
        unit = 1.
    elif (key == 'va'):
        dkey = 'vA'
        lab_name = r'Alfv\'en Speed'
        lab_sym = r'v_{\rm A}'
        lab_unit = r'{\rm km}/{\rm s}'
        corner = 3
        yscale = 'log'
        ylim_scaled = [0.07, 8.]
        ylim_unscaled = [8., 200.]
        unit = 1.E-3
    elif (key == 'na'):
        dkey = 'alfven_ratio'
        lab_name = r'Alfv\'en Number'
        lab_sym = r'(v_{\rm p}/v_{\rm A})'
        lab_unit = None
        corner = 4
        yscale = 'log'
        ylim_scaled = [0.1, 5.0]
        ylim_unscaled = [0.2, 50.]
        unit = 1.
    elif (key == 'ang'):
        dkey = 'parker_angle'
        lab_name = r'Parker-Angle'
        # lab_sym = r'\arccos\lvert\hat{\bf B}\cdot\hat{\bf r}\rvert'
        lab_sym = r'\arccos \lvert'
        lab_unit = r'{\rm deg}'
        corner = 4
        yscale = 'linear'
        ylim_scaled = [0., 2.]
        ylim_unscaled = [0, 90.]
        unit = 180. / np.pi
    elif (key == 'loss'):
        dkey = 'particle_flux'
        lab_name = r'Proton Loss Rate'
        lab_sym = r'4 \pi r^2 n_{\rm p} v_{\rm p}'
        lab_unit = r'10^{36}/{\rm s}'
        corner = 2
        yscale = 'linear'
        ylim_scaled = [0., 3.9]
        ylim_unscaled = [0., 2.]
        unit = 4. * np.pi / 1.E36
    elif (key == 'lat'):
        dkey = 'heliographicLatitude'
        lab_name = r'Heliographic-Latitude'
        lab_sym = r'\lambda'
        lab_unit = r'{\rm deg}'
        corner = 3
        yscale = 'linear'
        ylim_scaled = [-180., 540.]
        ylim_unscaled = [-180., 540.] 
        unit = 1.
    elif (key == 'long'):
        dkey = 'heliographicLongitude'
        lab_name = r'Heliographic-Longitude'
        lab_sym = r'\lambda'
        lab_unit = r'{\rm deg}'
        corner = 3
        yscale = 'linear'
        ylim_scaled = [-180., 540.]
        ylim_unscaled = [-180., 540.]
        unit = 1.
    elif (key == 'proton_beta'):
        dkey = 'proton_beta'
        lab_name = r'Proton-Beta'
        lab_sym = r'\beta_{\rm p}'
        lab_unit = None
        corner = 3
        yscale = 'log'
        ylim_scaled = [0.008, 6.0]
        ylim_unscaled = [0.008, 6.0]
        unit = 1.
    elif (key == 'alfven_ratio'):
        dkey = 'proton_beta'
        lab_name = r'Alfven-Ratio'
        lab_sym = r'v_{\rm mp}/v_{\rm A}'
        lab_unit = None
        corner = 3
        yscale = 'log'
        ylim_scaled = [0.008, 2.0]
        ylim_unscaled = [1, 50]
        unit = 1.

    # Generate the y-axis label for the individual plot.
    ylabel_indiv = r'\rm '

    if (lab_unit is None):
        ylabel_indiv += (lab_name + r', $' + lab_sym + '$')
    else:
        ylabel_indiv += (lab_name + r', $' + lab_sym + r'$ [$' + lab_unit + r'$]')

    # Generate the y-axis label for the composite plot.

    ylabel_comp = r'\rm '

    if (scaled):
        ylabel_comp += (r'Scaled ' + lab_name + r', $' + lab_sym + r'/\langle ' + lab_sym + r'\rangle_{\oplus}$')
    else:
        if (lab_unit is None):
            ylabel_comp += (lab_name + r', $' + lab_sym + '$')
        else:
            ylabel_comp += (lab_name + r', $' + lab_sym + r'$ [$' + lab_unit + r'$]')

    # Generate the other plot settings.

    xlim = [0.05, 150.]
    xscale = 'log'
    xlabel = r'\rm Distance from Sun, $r$ [au]'

    hlim = [10.**0, 10.**4.2]  # [10.**1.1, 10.**4.2]
    hscale = 'log'
    hlabel = r'\rm Intervals'

    ylim_indiv = ylim_scaled if scaled else ylim_unscaled
    ylim_comp = ylim_scaled if scaled else ylim_unscaled

    # Initialize the figure.

    fig = plt.figure(figsize=figsize)

    matplotlib.rc('text', usetex=True)

    # Initialize the histogram.

    fig_hist = fig.add_axes(pos_hist)

    fig_hist.set_xlim(xlim)
    fig_hist.set_ylim(hlim)

    fig_hist.set_xscale(xscale)
    fig_hist.set_yscale(hscale)

    fig_hist.set_ylabel(hlabel, size='small')

    fig_hist.set_xticklabels([])

    # Initialize the individual plot.

    fig_indiv = fig.add_axes(pos_indiv)

    fig_indiv.set_xlim(xlim)
    fig_indiv.set_ylim(ylim_indiv)

    fig_indiv.set_xscale(xscale)
    fig_indiv.set_yscale(yscale)

    fig_indiv.set_ylabel(ylabel_indiv, size='small')

    fig_indiv.set_xticklabels([])

    # Initialize the composite plot.

    fig_comp = fig.add_axes(pos_comp)

    fig_comp.set_xlim(xlim)
    fig_comp.set_ylim(ylim_comp)

    fig_comp.set_xscale(xscale)
    fig_comp.set_yscale(yscale)

    fig_comp.set_xlabel(xlabel, size='small')
    fig_comp.set_ylabel(ylabel_comp, size='small')

    # Initialize the key.
    if (corner == 1):
        pos_key = [0.43, 0.77, 0.55, 0.06]
    elif (corner == 2):
        pos_key = [0.13, 0.77, 0.55, 0.06]
    elif (corner == 3):
        pos_key = [0.13, 0.46, 0.55, 0.06]
    elif (corner == 4):
        pos_key = [0.43, 0.46, 0.55, 0.06]
    else:
        pos_key = [0.13, 0.46, 0.55, 0.06]

    fig_key = fig.add_axes(pos_key)

    fig_key.set_xlim([0., 4.05])  # Changed from 4. for New Horizons
    fig_key.set_ylim([0., 3.])

    fig_key.set_xticks([])
    fig_key.set_yticks([])

    # Generate the solar-system object markers.

    d_mercury = 0.3871
    d_venus = 0.7233
    d_earth = 1.000
    d_mars = 1.524
    d_ceres = 2.767
    d_jupiter = 5.204
    d_saturn = 9.583
    d_uranus = 19.22
    d_neptune = 30.11
    d_pluto = 39.48
    d_eris = 67.781

    lst_d = [d_mercury, d_venus, d_earth, d_mars, d_ceres,
             d_jupiter, d_saturn, d_uranus, d_neptune, d_pluto, d_eris]

    lst_lab = [r'\starfontserif \Mercury', r'\starfontserif \Venus',
               r'\starfontserif \Terra', r'\starfontserif \Mars',
               r'\starfontserif \Ceres', r'\starfontserif \Jupiter',
               r'\starfontserif \Saturn', r'\starfontserif \Uranus',
               r'\starfontserif \Neptune', r'\starfontserif \Pluto',
               r'\starfontserif \Mars']

    y = ((inch_mrgn_bottom + 2. * inch_plot_y + inch_hist_y) / inch_fig_y) + 0.005

    for i in range(len(lst_d)):

        d = lst_d[i]

        if ((d == d_ceres) or (d == d_pluto) or (d == d_eris)):
            c = 'gray'
        else:
            c = 'black'

        if (d == d_eris):
            shift_x = -1.3
            shift_y = 0.000
            rotation = 221
        else:
            shift_x = 0
            shift_y = 0
            rotation = 0

        fig_hist.annotate(lst_lab[i], (d + shift_x, y + shift_y),
                          xycoords=('data', 'figure fraction'), color=c, size='medium', horizontalalignment='center', verticalalignment='bottom', rotation=rotation)

        fig_hist.plot((d, d), hlim, color=c, linewidth=0.5, linestyle=':')
        fig_indiv.plot((d, d), ylim_indiv, color=c, linewidth=0.5, linestyle=':')
        fig_comp.plot((d, d), ylim_comp, color=c, linewidth=0.5, linestyle=':')

    # If relevant, add the unit-scale axis.

    if (scaled):

        fig_comp.plot(xlim, [1., 1.], color='k', linewidth=0.5, linestyle=':')

    # Load data for and populate the individual plot and histogram.

    n_sc = 12

    sc_name = [r'$\rm PSP$', r'$\rm Helios~1$', r'$\rm Helios~2$', r'$\rm Ulysses$',
               r'$\rm Mariner~2$', r'$\rm Mariner~10$', r'$\rm Cassini$', r'$\rm Pioneer~10$',
               r'$\rm Pioneer~11$', r'$\rm New~Horizons$', r'$\rm Voyager~1$', r'$\rm Voyager~2$']

    sc_color = ['tab:brown', 'tab:orange', 'tab:olive', 'w', 'k', 'tab:gray', 'tab:green',
                'tab:purple', 'tab:pink', 'tab:red', 'tab:blue', 'tab:cyan']

    sc_marker = ['d', '^', 'v', 'o', 'P', 'X', 'p', '<', '>', '*', 's', 'D']

    sc_markeredgewidth = [0, 0, 0,
                          0.5, 0, 0,
                          0, 0, 0,
                          0, 0, 0]

    sc_markersize = [3, 3, 3,
                     2, 3, 3,
                     3, 3, 3,
                     3, 3, 3]

    sc_markeredgecolor = [None, None, None,
                          'k', None, None,
                          None, None, None,
                          None, None, None]

    if scaled:
        p1 = '/mnt/cephadrius/udel_research/msc/data/binned/scaled/v2024.1/'
        p2 = '_coho1hr_merged_mag_plasma_'
        p3 = "_80_binned_scaled_v2024.1.hf"
    else:
        p1 = "/mnt/cephadrius/udel_research/msc/data/binned/unscaled/v2024.1/"
        p2 = '_coho1hr_merged_mag_plasma_'
        p3 = "_v2024.1.hf"

    sc_fl = [
        p1 + "psp" + p2 + "20180101_20231001" + p3,
        p1 + "helios1" + p2 + "19740101_19811201" + p3,
        p1 + "helios2" + p2 + "19760101_19801201" + p3,
        p1 + "uy" + p2 + "19900101_19920201" + p3,
        p1 + "mariner2" + p2 + "19620830_19621116" + p3,
        p1 + "mariner10" + p2 + "19731103_19740918" + p3,
        p1 + "cassini" + p2 + "20000101_20040101" + p3,
        p1 + "pioneer10" + p2 + "19720101_19950901" + p3,
        p1 + "pioneer11" + p2 + "19730101_19941201" + p3,
        p1 + "new_horizons" + p2 + "20081010_20230731" + p3,
        p1 + "voyager1" + p2 + "19770101_20181201" + p3,
        p1 + "voyager2" + p2 + "19770101_20181201" + p3,
    ]

    for sc in range(n_sc):

        # Load and select the data.

        dat = h5py.File(sc_fl[sc], 'r')

        dat_r = np.array(dat['sc_r_iqr_50'])
        dat_k = np.array(dat[dkey + '_iqr_50'])
        dat_n = np.array(dat[dkey + '_num'])

        tk = np.where((np.isfinite(dat_r)) &
                      (np.isfinite(dat_k)) &
                      (np.isfinite(dat_n)) &
                      (dat_n >= 1))[0]

        sel_r = dat_r[tk]
        sel_k = unit * dat_k[tk]
        sel_n = dat_n[tk]

        # Populate the individual plot.

        fig_indiv.plot(sel_r, sel_k,
                       linestyle='', color=sc_color[sc],
                       marker=sc_marker[sc],
                       markersize=sc_markersize[sc],
                       markeredgewidth=sc_markeredgewidth[sc],
                       markeredgecolor=sc_markeredgecolor[sc])

        # Populate the histogram.

        fig_hist.plot(sel_r, sel_n,
                      linestyle='', color=sc_color[sc],
                      marker=sc_marker[sc],
                      markersize=sc_markersize[sc],
                      markeredgewidth=sc_markeredgewidth[sc],
                      markeredgecolor=sc_markeredgecolor[sc])

        # Make a horizontal line for the median.
        fig_hist.axhline(25, color='k', linestyle='--')

        # Populate the key.

        x_sc = sc // 3
        y_sc = 2 - (sc % 3)

        opt_sc = 0. if (x_sc == 0) else -0.15

        fig_key.plot([x_sc + 0.1 + opt_sc], [y_sc + 0.5],
                     linestyle='', color=sc_color[sc],
                     marker=sc_marker[sc],
                     markersize=sc_markersize[sc],
                     markeredgewidth=sc_markeredgewidth[sc],
                     markeredgecolor=sc_markeredgecolor[sc])

        fig_key.annotate(sc_name[sc],
                         (x_sc + 0.175 + opt_sc, y_sc + 0.5),
                         xycoords='data',
                         color='k', size='x-small',
                         horizontalalignment='left',
                         verticalalignment='center')

    # Load data for and populate the composite.

    p = "/mnt/cephadrius/udel_research/msc/data/all_data/v2024.1/all_spacecraft_data_"

    if (scaled):
        dat = h5py.File(p + 'scaled_80_binned_v2024.1.hf', 'r')
    else:
        dat = h5py.File(p + "_80_binned_v2024.1.hf", "r")

    dat_r = np.array(dat['sc_r_iqr_50'])
    dat_10 = np.array(dat[dkey + '_iqr_10'])
    dat_25 = np.array(dat[dkey + '_iqr_25'])
    dat_50 = np.array(dat[dkey + '_iqr_50'])
    dat_75 = np.array(dat[dkey + '_iqr_75'])
    dat_90 = np.array(dat[dkey + '_iqr_90'])
    dat_num = np.array(dat[dkey + '_num'])

    # if (key == 'b'):
    #     print(transpose([dat_r, dat_num]))

    tk = np.where((np.isfinite(dat_r)) & (np.isfinite(dat_25)) &
                  (np.isfinite(dat_50)) & (np.isfinite(dat_75)) &
                  (np.isfinite(dat_90)) & (np.isfinite(dat_num)) &
                  (dat_num >= 1))[0]

    n_sel = len(tk)

    sel_r = dat_r[tk]
    sel_10 = dat_10[tk]
    sel_25 = dat_25[tk]
    sel_50 = dat_50[tk]
    sel_75 = dat_75[tk]
    sel_90 = dat_90[tk]
    sel_num = dat_num[tk]

    if (not scaled):
        sel_10 = unit * sel_10
        sel_25 = unit * sel_25
        sel_50 = unit * sel_50
        sel_75 = unit * sel_75
        sel_90 = unit * sel_90

    for i in range(n_sel):

        if (scaled):
            color = 'm'
        else:
            color = 'b'

        fig_comp.plot([sel_r[i], sel_r[i]],
                      [sel_10[i], sel_90[i]],
                      color=color, linestyle='-',
                      linewidth=0.5)

        fig_comp.plot([sel_r[i], sel_r[i]],
                      [sel_25[i], sel_75[i]],
                      color=color, linestyle='-',
                      linewidth=1.0)

        fig_comp.plot(sel_r[i], sel_50[i],
                      color=color, linestyle='', marker='o',
                      markersize=3.0, markeredgewidth=0)

    # Analyze the composite data and add the analysis to the composite plot.

    if (key in ['loss', 'betap']):

        fit = curve_fit(mod_power, sel_r, sel_50,
                        p0=(1., 0.),
                        sigma=sel_50 / np.sqrt(sel_num))

        crv_r = np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), 300)
        crv_k = mod_power(crv_r, fit[0][0], fit[0][1])

        fig_comp.plot(crv_r, crv_k, color='tab:blue', linestyle='-', linewidth=1.0)

        fig_fit = fig.add_axes([0.13, 0.07, 0.25, 0.03])

        fig_fit.set_xlim([0., 1.])
        fig_fit.set_ylim([0.5, 1.5])

        fig_fit.set_xticks([])
        fig_fit.set_yticks([])

        txt = r'$\alpha = {:.4f} \pm {:.4f}$'.format(fit[0][1], np.sqrt(fit[1][1][1]))

        fig_fit.annotate(txt,
                         (0.05, 1.0),
                         xycoords='data',
                         color='tab:blue',
                         size='medium',
                         horizontalalignment='left',
                         verticalalignment='center')

    elif (key in ['np', 'vp', 'vpr']):

        if (key in ['np']):	
            p0 = [1., 4., -2., -2.]
        elif (key in ['vp', 'vpr']):	
            p0 = [1., 4., 0., 0.]
        elif (key in ['va']):
            p0 = [1., 3.0, -0.6, 0.]
        elif (key in ['na']):	
            p0 = [1., 3.0, 0.6, 0.2]
        elif (key in ['betap']):
            p0 = [1., 8., -0.4, -0.2]
        else:
            p0 = [1., 4., 0., 0.]

        fit = curve_fit(mod_broken, sel_r, sel_50,
                        p0=p0,
                        sigma=sel_50 / np.sqrt(sel_num))

        crv_r = np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), 300)
        crv_k = mod_broken(crv_r, fit[0][0], fit[0][1], fit[0][2], fit[0][3])
        crv_k1 = mod_broken(crv_r, fit[0][0], fit[0][1], fit[0][2], fit[0][3], force=1)
        crv_k2 = mod_broken(crv_r, fit[0][0], fit[0][1], fit[0][2], fit[0][3], force=2)

        fig_comp.plot(crv_r, crv_k2,
                      color='tab:blue', linestyle=':', linewidth=1.0)

        fig_comp.plot(crv_r, crv_k1,
                      color='tab:blue', linestyle=':', linewidth=1.0)

        fig_comp.plot(crv_r, crv_k,
                      color='tab:blue', linestyle='-', linewidth=1.0)

        fig_comp.plot([fit[0][1], fit[0][1]], ylim_comp,
                      color='tab:blue', linestyle='--', linewidth=1.0)

        fig_fit = fig.add_axes([0.13, 0.07, 0.25, 0.09])

        fig_fit.set_xlim([0., 1.])
        fig_fit.set_ylim([0.5, 3.5])

        fig_fit.set_xticks([])
        fig_fit.set_yticks([])

        txt_a1 = r'$\alpha_1={:.3f} \pm {:.3f}$'.format(fit[0][2], np.sqrt(fit[1][2][2]))
        txt_a2 = r'$\alpha_2={:.3f} \pm {:.3f}$'.format(fit[0][3], np.sqrt(fit[1][3][3]))

        txt_rc_a = r'$r_c=\left({:.2f}\pm {:.2f}'.format(fit[0][1], np.sqrt(fit[1][1][1]))
        txt_rc_b = r'\right)\,{\rm au}$'

        arr_txt = [txt_rc_a + txt_rc_b, txt_a2, txt_a1]

        for t in range(3):

            txt = arr_txt[t]

            fig_fit.annotate(txt,
                             (0.05, 1.0 + t),
                             xycoords='data',
                             color='tab:blue',
                             size='medium',
                             horizontalalignment='left',
                             verticalalignment='center')

    elif (key in ['b', 'tp', 'va', 'na']):

        bounds = (-np.inf, np.inf)

        if (key in ['va']):
            p0 = [1., 0.6, 5.0, -0.8, -0.4, 0.]
        elif (key in ['na']):
            p0 = [1., 0.4, 5.0, 0.6, 0.4, 0.]
            bounds = ((0, 0.2, 1., -np.inf, -np.inf, -np.inf),
                       (np.inf, 1., 10., np.inf, np.inf, np.inf))
        elif (key in ['tp']):
            p0 = [1., 0.2, 30., 0.2, -0.7, 0.2]
        else:
            p0 = [1., 0.3, 4., -2., -1., -0.5]

        fit = curve_fit(mod_broken2, sel_r, sel_50,
                        p0=p0, bounds=bounds,
                        sigma=sel_50 / np.sqrt(sel_num))

        if (key == 'na'):
            print('    Median Alfven point:', (1. / fit[0][0])**(1. / fit[0][3]), 'au')

        crv_r = np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), 300)
        crv_k = mod_broken2(crv_r, fit[0][0], fit[0][1], fit[0][2],
                            fit[0][3], fit[0][4], fit[0][5])
        crv_k1 = mod_broken2(crv_r, fit[0][0], fit[0][1], fit[0][2],
                             fit[0][3], fit[0][4], fit[0][5], force=1)
        crv_k2 = mod_broken2(crv_r, fit[0][0], fit[0][1], fit[0][2],
                             fit[0][3], fit[0][4], fit[0][5], force=2)
        crv_k3 = mod_broken2(crv_r, fit[0][0], fit[0][1], fit[0][2],
                             fit[0][3], fit[0][4], fit[0][5], force=3)

        fig_comp.plot(crv_r, crv_k3, color='tab:blue',
                      linestyle=':', linewidth=1.0)
        fig_comp.plot(crv_r, crv_k2,
                      color='tab:blue',
                      linestyle=':', linewidth=1.0)
        fig_comp.plot(crv_r, crv_k1,
                      color='tab:blue',
                      linestyle=':', linewidth=1.0)
        fig_comp.plot(crv_r, crv_k,
                      color='tab:blue',
                      linestyle='-', linewidth=1.0)

        for i in [1, 2]:
            fig_comp.plot([fit[0][i], fit[0][i]], ylim_comp,
                          color='tab:blue', linestyle='--', linewidth=1.0)

        fig_fit = fig.add_axes([0.13, 0.07, 0.25, 0.15])

        fig_fit.set_xlim([0., 1.])
        fig_fit.set_ylim([0.5, 5.5])

        fig_fit.set_xticks([])
        fig_fit.set_yticks([])

        txt_a1 = r'$\alpha_1={:.3f} \pm {:.3f}$'.format(fit[0][3], np.sqrt(fit[1][3][3]))
        txt_a2 = r'$\alpha_2={:.3f} \pm {:.3f}$'.format(fit[0][4], np.sqrt(fit[1][4][4]))
        txt_a3 = r'$\alpha_3={:.3f} \pm {:.3f}$'.format(fit[0][5], np.sqrt(fit[1][5][5]))

        txt_ra = r'$r_a=\left({:.2f}\pm {:.2f}'.format(fit[0][1], np.sqrt(fit[1][1][1]))
        txt_rb = r'$r_b=\left({:.2f}\pm {:.2f}'.format(fit[0][2], np.sqrt(fit[1][2][2]))

        txt_ra += r'\right)\,{\rm au}$'
        txt_rb += r'\right)\,{\rm au}$'

        arr_txt = [txt_rb, txt_ra, txt_a3, txt_a2, txt_a1]

        for t in range(5):

            txt = arr_txt[t]

            fig_fit.annotate(txt,
                             (0.05, 1.0 + t),
                             xycoords='data',
                             color='tab:blue',
                             size='medium',
                             horizontalalignment='left',
                             verticalalignment='center')

    # Output the figure.

    if (scaled):
        # fname = 'svg/plot_scaled_' + key + '.svg'
        fname = '/home/cephadrius/Dropbox/v2024.1/plot_scaled_' + key + '.pdf'
    else:
        # fname = 'svg/plot_unscaled_' + key + '.svg'
        fname = '/home/cephadrius/Dropbox/v2024.1/plot_unscaled_' + key + '.pdf'

    plt.savefig(fname, dpi=300, bbox_inches='tight', pad_inches=0.02, transparent=False)


# Generate the figures.

n_arr = 10
arr_key = ['b', 'np', 'vp', 'vpr', 'tp', 'va', 'na', 'loss', 'ang', 'proton_beta', 'alfven_ratio']
r_min = -1.2
r_max = 2
n_bin = 80
r_bin = np.logspace(r_min, r_max, n_bin)

# arr_scaled = [True, True, True, True, True, True, False, True, False, False]

for i in range(8, 9):
    print(arr_key[i])
    for scaled in [False, True]:
        plot_msc(arr_key[i], scaled=scaled)
