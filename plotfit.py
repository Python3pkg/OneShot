#!/usr/bin/env python
import numpy as _np
import mytools as _mt


def plotfit(
        x                      ,
        y                      ,
        beta                   ,
        X                      ,
        top      = None        ,
        bottom   = 'Arbitrary' ,
        figpath  = None        ,
        error    = None        ,
        figlabel = None        ,
        axes     = None        ,
        fontsize = None        ,
        **kwargs
        ):
    # ======================================
    # Set up PDF saving
    # ======================================
    # print axes
    # if axes is None:
    #     if not (figlabel == None):
    #         fig = _mt.figure(figlabel)
    #     else:
    #         fig = _plt.figure()
    #     axes=fig.add_subplot(1, 1, 1)
    # else:
    #     pass
    #     # fig = axes.get_figure()

    # Convert from meters to um
    y_data_mm_sq = y * 1e6
    y_fit_mm_sq = _np.dot(X, beta) * 1e6
    # error_mm_sq = error * 1e6

    # =================================================
    # VARIANCE FITS
    # =================================================
    # Plot data with error bars
    # axes.errorbar(x, y_data_mm_sq, error_mm_sq, fmt='.-')
    # axes.plot(x, y_data_mm_sq, 'o-', **kwargs)

    # Plot fits
    # axes.plot(x, y_fit_mm_sq, '-', **kwargs)

    # =================================================
    # SIGMA FITS
    # =================================================
    # Plot data with error bars
    # axes.errorbar(x, y_data_mm_sq, error_mm_sq, fmt='.-')
    axes.plot(x, _np.sqrt(y_data_mm_sq)*1e3, 'o-', **kwargs)

    # Plot fits
    axes.plot(x, _np.sqrt(y_fit_mm_sq)*1e3, '-', **kwargs)

    # _mt.addlabel(top, bottom, '$\sigma_x^2$ [mm$^2$]')
    if fontsize is None:
        axes.legend(['Measured Slices', 'Fit to Measurement'])
    else:
        axes.legend(['Measured Slices', 'Fit to Measurement'], fontsize=fontsize)
    _mt.addlabel(axes=axes, xlabel='Slice Energy [GeV]', ylabel='Slice Spot Size $\\sigma_x$ [$\\mu$m]', toplabel=top)
    # if fig is not None:
    #         fig.tight_layout()

    # if not (figpath == None):
    #         _mt.graphics.savefig(figlabel, figpath)

    return axes
