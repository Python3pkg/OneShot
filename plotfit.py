#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import mytools as mt

def plotfit(x,y,beta,X,spotexpected,top,elepath='figs',error=None):
	# Recast as fractional energy
	e_frac=x+1
	# Plot data with error bars
	plt.errorbar(e_frac,y*(1e3**2),error*(1e3**2),fmt='.-')

	# Plot fits
	plt.plot(e_frac,np.dot(X,beta)*(1e3**2),'-')

	mt.addlabel(top,'$E/E_0$','$\sigma_x^2$ [mm$^2$]')
	plt.legend(['Simulated Data','Fit to Data'])
	mt.graphics.savefig(top,elepath)
