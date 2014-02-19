#!/usr/bin/env python
import numpy as _np
import matplotlib.pyplot as _plt
import mytools as _mt

def plotfit(x,y,beta,X,top,bottom='Arbitrary',elepath=None,error=None,figlabel=None):
	# Convert from meters to um
	y_data_mm_sq = y * 1e6
	y_fit_mm_sq = _np.dot(X,beta) * 1e6
	error_mm_sq = error * 1e6

	if not (figlabel == None):
		fig = _mt.figure(figlabel)
	else:
		fig = _plt.figure()

	# Plot data with error bars
	_plt.errorbar(x,y_data_mm_sq,error_mm_sq,fmt='.-')
	# _plt.plot(x,y_data_mm_sq,'-')

	# _mt.figure('splitme')

	# Plot fits
	_plt.plot(x,y_fit_mm_sq,'.-')

	_mt.addlabel(top,bottom,'$\sigma_x^2$ [mm$^2$]')
	_plt.legend(['Simulated Data','Fit to Data'])

	if not (elepath == None):
		_mt.graphics.savefig(top,elepath)

	return fig
