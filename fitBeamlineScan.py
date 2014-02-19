#!/usr/bin/env python
import numpy as _np
import copy as _copy
import collections as _col
import mytools as _mt
import matplotlib.pyplot as _plt
import pdb as _pdb

# Fit bowtie {{{
def fitBeamlineScan(beamline,y,emitx,error=None,verbose=False):
	beamline_manip = _copy.deepcopy(beamline)
	numsteps = beamline.size
	y              = y[_np.newaxis]
	X              = _np.zeros([numsteps,3])
	spotexpected   = _np.zeros(numsteps)
	
	for i,bl in enumerate(beamline):
		R11             = bl.R[0,0]
		R12             = bl.R[0,1]
		X[i,0]          = R11*R11
		X[i,1]          = 2*R11*R12
		X[i,2]          = R12*R12
		# spotexpected[i] = _np.sqrt((R11*R11*twiss.beta - 2*R11*R12*twiss.alpha + R12*R12*twiss.gamma)*emitx)
		# twiss=beamline.twiss_x
		# spotexpected[i] = bl.twiss.transport(bl.R[0:2,0:2]).spotsize(emitx)
		spotexpected[i] = bl.spotsize_x_end(emitx)
	
	if verbose:
		# _mt.figure('Expected')
		# xax=_np.linspace(1,numsteps,numsteps)
		# plt.plot(xax,_np.sqrt(y.transpose()),xax,spotexpected)
		# plt.show()
		pass

	# beta is the best solution of beam parameters:
	# sig^2 = R11^2 <x^2> + 2 R11 R12 <xx'> + R12^2 <x'^2>
	# beta[0,0] = <x^2>
	# beta[1,0] = <xx'>
	# beta[2,0] = <x'^2>

	X_unweighted = _copy.deepcopy(X)
	y_unweighted = _copy.deepcopy(y)
	if ( error != None ):
	# if ( False ):
		for i,el in enumerate(X):
			X[i,:] = el/error[i]
		y_err = y/error

	# This is the linear least squares matrix formalism
	y_err = y_err.transpose()
	# _pdb.set_trace()
	beta  = _np.dot(_np.linalg.pinv(X) , y_err)
	covar = _np.linalg.inv(_np.dot(_np.transpose(X),X))
	
	emit = _np.sqrt( beta[0,0] * beta[2,0] - _np.square(beta[1,0]) )
	del_emit_sq = _np.power(1/(2*emit),2) * \
		(
			_np.power(beta[2,0] * covar[0,0]   ,2) +
			_np.power(2*beta[1,0] * covar[1,1] ,2) +
			_np.power(beta[0,0] * covar[2,2]   ,2) # +
			# 2 * (-beta[2,0] * 2 * beta[1,0]) * covar[0,1] +
			# 2 * (beta[2,0] * beta[0,0]) * covar[0,2] +
			# 2 * (-2*beta[1,0] * beta[0,0]) * covar[1,2]
		)

	chisq_red = _mt.chisquare(y.transpose(),_np.dot(X_unweighted,beta),error,ddof=3,verbose=verbose)
# def chisquare(observe,expect,error,ddof,verbose=True):

	if verbose:
		print 'Emittance error is:\t\t{}.'.format(_np.sqrt(del_emit_sq))
		beta0 = beta[0,0]/emit
		gamma0 = beta[2,0]/emit
		alpha0 = -_np.sign(beta[1,0])*_np.sqrt(beta0*gamma0-1)
		print 'Emittance fit:\t\t\t{}.'.format(emit)
		print 'Normalized emittance fit:\t{}.'.format(emit*40000)
		print 'Initial beta fit:\t\t{}.'.format(beta0)
		print 'Initial alpha fit:\t\t{}.'.format(alpha0)
		print 'Initial gamma fit:\t\t{}.'.format(gamma0)
		print 'Initial spot from fit:\t\t{}.'.format(_np.sqrt(beta[0,0]))

	output = _col.namedtuple('fitbowtie_output',['spotexpected','X','X_unweighted','beta','covar','chisq_red'])
	out = output(spotexpected,X,X_unweighted,beta,covar,chisq_red)

	return out
#}}}
