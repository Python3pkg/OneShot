#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import mytools as mt

def plotfit(filt,davg,stddev,beta,X,spotexpected,top,elepath):
	plt.plot(davg[filt]+1,stddev[filt]*1e6,'.',davg[filt]+1,np.sqrt(np.dot(X,beta))*1e6,'-',davg[filt]+1,spotexpected*1e6,'-')
	mt.addlabel(top,'$E/E_0$','$\sigma_x$ [$\mu$m]')
	# plt.legend(['Simulated Data','Fit to Data','Expected (via\nInitial Conditions)'])
	mt.graphics.savefig(top,elepath)
