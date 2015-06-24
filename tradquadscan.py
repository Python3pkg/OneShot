#!/usr/bin/env python
import numpy as np
import copy
import collections as col
import mytools as mt


# Fit bowtie {{{
def tradquadscan(beamline, y, twiss, emitx, error=None, verbose=False):
    # beamline_manip = copy.deepcopy(beamline)
    numsteps = beamline.size
    betax          = twiss.beta
    alphax         = twiss.alpha
    gammax         = twiss.gamma
    y              = y[np.newaxis]
    # gamma          = (1+x)*40000
    X              = np.zeros([numsteps, 3])
    spotexpected   = np.zeros(numsteps)
    R              = np.zeros([6, 6, numsteps])
    betaf          = np.zeros(numsteps)
    
    for i, bl in enumerate(beamline):
        R11             = bl.R[0, 0]
        R12             = bl.R[0, 1]
        R[:, :, i]        = bl.R
        X[i, 0]          = R11*R11
        X[i, 1]          = 2*R11*R12
        X[i, 2]          = R12*R12
        T2              = np.dot(np.dot(R[0:2, 0:2, i], twiss.T), np.transpose(R[0:2, 0:2, i]))
        betaf[i]        = T2[0, 0]
        spotexpected[i] = np.sqrt((R11*R11*betax - 2*R11*R12*alphax + R12*R12*gammax)*emitx)
    
    if verbose:
        # mt.figure('Expected')
        # xax=np.linspace(1, numsteps, numsteps)
        # plt.plot(xax, np.sqrt(y.transpose()), xax, spotexpected)
        # plt.show()
        pass

    # beta is the best solution of beam parameters:
    # sig^2 = R11^2 <x^2> + 2 R11 R12 <xx'> + R12^2 <x'^2>
    # beta[0, 0] = <x^2>
    # beta[1, 0] = <xx'>
    # beta[2, 0] = <x'^2>

    X_unweighted = copy.deepcopy(X)
    # if ( False ):
    if error is not None:
        for i, el in enumerate(X):
            X[i, :] = el/error[i]
        y_err = y/error

    # This is the linear least squares matrix formalism
    y_err = y_err.transpose()
    beta  = np.dot(np.linalg.pinv(X) , y_err)
    covar = np.linalg.inv(np.dot(np.transpose(X), X))
    
    emit = np.sqrt( beta[0, 0] * beta[2, 0] - np.square(beta[1, 0]) )
    del_emit_sq = np.power(1/(2*emit), 2) * \
        (
            np.power(beta[2, 0] * covar[0, 0]   , 2) +
            np.power(2*beta[1, 0] * covar[1, 1] , 2) +
            np.power(beta[0, 0] * covar[2, 2]   , 2)  # +
            # 2 * (-beta[2, 0] * 2 * beta[1, 0]) * covar[0, 1] +
            # 2 * (beta[2, 0] * beta[0, 0]) * covar[0, 2] +
            # 2 * (-2*beta[1, 0] * beta[0, 0]) * covar[1, 2]
        )

    chisq_red = mt.chisquare(y.transpose(), np.dot(X_unweighted, beta), error, ddof=3, verbose=verbose)
# def chisquare(observe, expect, error, ddof, verbose=True):

    if verbose:
        print('Emittance error is:\t\t{}.'.format(np.sqrt(del_emit_sq)))
        beta0 = beta[0, 0]/emit
        gamma0 = beta[2, 0]/emit
        alpha0 = -np.sign(beta[1, 0])*np.sqrt(beta0*gamma0-1)
        print('Emittance fit:\t\t\t{}.'.format(emit))
        print('Normalized emittance fit:\t{}.'.format(emit*40000))
        print('Initial beta fit:\t\t{}.'.format(beta0))
        print('Initial alpha fit:\t\t{}.'.format(alpha0))
        print('Initial gamma fit:\t\t{}.'.format(gamma0))
        print('Initial spot from fit:\t\t{}.'.format(np.sqrt(beta[0, 0])))

    output = col.namedtuple('fitbowtie_output', ['spotexpected', 'X', 'X_unweighted', 'beta', 'covar', 'chisq_red'])
    out = output(spotexpected, X, X_unweighted, beta, covar, chisq_red)

    return out
# }}}
