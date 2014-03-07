#!/usr/bin/env python
import numpy as _np
import matplotlib.pyplot as _plt
import mytools as _mt
from matplotlib.backends.backend_pdf import PdfPages

def plotfit(x,y,beta,X,top,bottom='Arbitrary',figpath=None,error=None,figlabel=None,axes=None):
	# ======================================
	# Set up PDF saving
	# ======================================
	if axes is None:
		if not (figlabel == None):
			fig = _mt.figure(figlabel)
		else:
			fig = _plt.figure()
		axes=_plt.gca()

	# Convert from meters to um
	y_data_mm_sq = y * 1e6
	y_fit_mm_sq = _np.dot(X,beta) * 1e6
	error_mm_sq = error * 1e6

	# Plot data with error bars
	axes.errorbar(x,y_data_mm_sq,error_mm_sq,fmt='.-')
	# axes.plot(x,y_data_mm_sq,'-')

	# _mt.figure('splitme')

	# Plot fits
	axes.plot(x,y_fit_mm_sq,'.-')

	_mt.addlabel(top,bottom,'$\sigma_x^2$ [mm$^2$]')
	axes.legend(['Measured Spot Variance','Fit Spot Variance'])

	if not (figpath == None):
		_mt.graphics.savefig(top,figpath)

	# return fig
