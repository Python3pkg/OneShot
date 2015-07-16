#!/usr/bin/env python
import numpy as np


# Cherenkov Spot size strips {{{
def histcher(x, y, res):
    h, xe, ye = np.histogram2d(x, y, res)
    xval      = (xe[1]-xe[0])/2. + xe
    xval      = xval[0:-1]
    
    dely = (ye[1]-ye[0])
    yavg = dely/2. + ye
    yavg = yavg[0:-1]
    
    etay = 0.070673884756039224
    davg = (etay/(etay-yavg)) - 1
    # davg = yavg/etay
    
    # filt=(davg>-0.01) & (davg < 0.01)

    out = [h, xval, davg]
    return out
# }}}
