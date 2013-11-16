#!/usr/bin/env python
import numpy as np

# Energy spot size strips {{{
def histenergy(x,d,res):
	h,xe,ye=np.histogram2d(x,d,res)
	xval = (xe[1]-xe[0])/2. + xe
	xval = xval[0:-1]
	
	davg=(ye[1]-ye[0])/2. + ye
	davg=davg[0:-1]
	
	# filt=(davg>-0.01) & (davg < 0.01)

	out = [h,xval,davg]
	return out
#}}}
