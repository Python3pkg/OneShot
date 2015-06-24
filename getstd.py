#!/usr/bin/env python
import numpy as np


# Get std dev (spot size) {{{
def getstd(res, h, xval):
    stddevsq  = np.zeros(res)
    indivbool = False
    if indivbool:
        figscan = plt.figure()  # noqa
    
    def gauss(x, A, mu, sig):
        return A*np.exp(-np.power(x-mu, 2)/(2*np.power(sig, 2)))
    
    for i, row in enumerate(np.transpose(h)):
        # A = max(row)
        mean = np.sum(xval*row)/row.sum()
        var = np.sum(np.power(xval-mean, 2)*row)/row.sum()
        # root = np.sqrt(var)
        # pguess = [A, mean, root]
        # popt = pguess
        # popt, pcov = spopt.curve_fit(gauss, xval, row, pguess)
        # # print "A: {}, mean: {}, sig: {}".format(popt[0], popt[1], popt[2])
        # # print "Percent diff: {}%".format(100*(popt[2]-root)/root)
        # fit = gauss(xval, popt[0], popt[1], popt[2])
        # unchangedroot = gauss(xval, popt[0], popt[1], root)
        # if indivbool: plt.plot(xval, row, xval, fit, xval, unchangedroot)
    
        # # plt.plot(xval, row)
        # if indivbool: raw_input("Any key.")
        # if indivbool: figscan.clf()
        # # stddevsq[i] = np.power(popt[2], 2)
        stddevsq[i] = var
    
    # stddev=np.sqrt(stddevsq)
    return stddevsq
# }}}
