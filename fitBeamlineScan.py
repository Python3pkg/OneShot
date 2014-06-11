#!/usr/bin/env python
import numpy as _np
import copy as _copy
import collections as _col
import mytools as _mt
import matplotlib.pyplot as _plt
import pdb as _pdb
import mytools.slactrac as _sltr

class LinLsqFit(object):
	def __init__(self,y_unweighted,X_unweighted,y_error=None):
		self._force_recalc()

		# print y_unweighted.shape
		# print X_unweighted.shape
		# print y_error.shape
		self.y_unweighted=y_unweighted
		self.y_error=y_error
		self.X_unweighted = X_unweighted

	# ======================================
	# Resets stored values for calculated
	# quantities
	# ======================================
	def _force_recalc(self):
		self._X = None
		self._y = None
		self._beta = None
		self._covar = None
		self._chisq_red = None
		self._emit = None
		self._twiss=None

	# ======================================
	# y_unweighted
	# ======================================
	def _get_y_unweighted(self):
		return self._y_unweighted
	def _set_y_unweighted(self,val):
		self._force_recalc()
		self._y_unweighted=val
	y_unweighted = property(_get_y_unweighted,_set_y_unweighted)

	# ======================================
	# y_error
	# ======================================
	def _get_y_error(self):
		return self._y_error
	def _set_y_error(self,val):
		self._force_recalc()
		self._y_error=val
	y_error = property(_get_y_error,_set_y_error)
	
	# ======================================
	# X_unweighted
	# ======================================
	def _get_X_unweighted(self):
		return self._X_unweighted
	def _set_X_unweighted(self,val):
		self._force_recalc()
		self._X_unweighted=val
	X_unweighted = property(_get_X_unweighted,_set_X_unweighted)

	# ======================================
	# X (calculated)
	# ======================================
	def _get_X(self):
		if self._X == None:
			X = _copy.deepcopy(self.X_unweighted)
			# print 'X shape is {}'.format(X.shape)
			for i,el in enumerate(X):
				X[i,:]=el/self.y_error[i]
			# print 'New X shape is {}'.format(X.shape)
			self._X = X
		return self._X
	X = property(_get_X)

	# ======================================
	# y (calculated)
	# ======================================
	def _get_y(self):
		if self._y==None:
			self._y = self.y_unweighted/self.y_error
		return self._y
	y = property(_get_y)

	# ======================================
	# beta (calculated)
	# ======================================
	def _get_beta(self):
		if self._beta==None:
			# This is the linear least squares matrix formalism
			self._beta = _np.dot(_np.linalg.pinv(self.X) , self.y)
		return self._beta
	beta = property(_get_beta)

	# ======================================
	# covar (calculated)
	# ======================================
	def _get_covar(self):
		if self._covar==None:
			self._covar=_np.linalg.inv(_np.dot(_np.transpose(self.X),self.X))
		return self._covar
	covar = property(_get_covar)
	
	# ======================================
	# chisq_red (calculated)
	# ======================================
	def _get_chisq_red(self):
		if self._chisq_red==None:
			self._chisq_red = _mt.chisquare(self.y_unweighted.transpose(),_np.dot(self.X_unweighted,self.beta),self.y_error,ddof=3,verbose=False)
		return self._chisq_red
	chisq_red=property(_get_chisq_red)

	# ======================================
	# emit (calculated)
	# ======================================
	def _get_emit(self):
		if self._emit==None:
			self._emit = _np.sqrt( self.beta[0,0] * self.beta[2,0] - _np.square(self.beta[1,0]) )
		return self._emit
	emit=property(_get_emit)

	# ======================================
	# twiss (calculated)
	# ======================================
	def _get_twiss(self):
		if self._twiss==None:
			beta0 = self.beta[0,0]/self.emit
			gamma0 = self.beta[2,0]/self.emit
			alpha0 = -_np.sign(self.beta[1,0])*_np.sqrt(beta0*gamma0-1)
			self._twiss = _sltr.Twiss(beta=beta0,alpha=alpha0)
		return self._twiss
	twiss=property(_get_twiss)


class BeamlineScanFit(object):
	def __init__(self,fitresults,spotexpected):
		self.fitresults=fitresults
		self.spotexpected=spotexpected

# Fit bowtie {{{
def fitBeamlineScan(beamline,y,emitx,error=None,verbose=False,plot=False,eaxis=None):
	beamline_manip = _copy.deepcopy(beamline)
	numsteps = beamline.size
	y              = y[_np.newaxis]
	error          = error[_np.newaxis].transpose()
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
	
	if plot:
		_mt.figure('Expected')
		xax=_np.linspace(1,numsteps,numsteps)
		_plt.plot(xax,_np.sqrt(y.transpose()),xax,spotexpected)
		_plt.show()
		pass

	# beta is the best solution of beam parameters:
	# sig^2 = R11^2 <x^2> + 2 R11 R12 <xx'> + R12^2 <x'^2>
	# beta[0,0] = <x^2>
	# beta[1,0] = <xx'>
	# beta[2,0] = <x'^2>

	myfit = LinLsqFit(y_unweighted=y.transpose(),X_unweighted=X,y_error=error)
	beta = myfit.beta
	covar = myfit.covar
	chisq_red = myfit.chisq_red
	emit = myfit.emit

	# emit = _np.sqrt( beta[0,0] * beta[2,0] - _np.square(beta[1,0]) )
	del_emit_sq = _np.power(1/(2*emit),2) * \
		(
			_np.power(beta[2,0] * covar[0,0]   ,2) +
			_np.power(2*beta[1,0] * covar[1,1] ,2) +
			_np.power(beta[0,0] * covar[2,2]   ,2) # +
			# 2 * (-beta[2,0] * 2 * beta[1,0]) * covar[0,1] +
			# 2 * (beta[2,0] * beta[0,0]) * covar[0,2] +
			# 2 * (-2*beta[1,0] * beta[0,0]) * covar[1,2]
		)

	# chisq_red = _mt.chisquare(y.transpose(),_np.dot(X_unweighted,beta),error,ddof=3,verbose=verbose)

	# beta0 = beta[0,0]/emit
	# gamma0 = beta[2,0]/emit
	# alpha0 = -_np.sign(beta[1,0])*_np.sqrt(beta0*gamma0-1)
	argmin = _np.argmin(_np.dot(myfit.X,myfit.beta))
	e_gamma = eaxis[argmin]/(0.5109989e-3)
	print eaxis[argmin]
	print e_gamma

	if verbose:
		print 'Emittance error is:\t\t{}.'.format(_np.sqrt(del_emit_sq))
		print 'Emittance fit:\t\t\t{}.'.format(emit)
		print 'Normalized emittance fit:\t{}.'.format(myfit.emit*e_gamma)
		print 'Initial beta fit:\t\t{}.'.format(myfit.twiss.beta)
		print 'Initial alpha fit:\t\t{}.'.format(myfit.twiss.alpha)
		print 'Initial gamma fit:\t\t{}.'.format(myfit.twiss.gamma)
		print 'Initial spot from fit:\t\t{}.'.format(_np.sqrt(myfit.beta[0,0]))
		print 'Min spot size (gauss fit): \t{}.'.format(min(_np.sqrt(y[0])))
		print 'Min spot size (emit fit): \t{}.'.format(myfit.twiss.minspotsize(myfit.emit))

		print 'Beta* is: \t\t\t{}.'.format(myfit.twiss.betastar)
		print 's* is: \t\t\t\t{}.'.format(myfit.twiss.sstar)

	# output = _col.namedtuple('fitbowtie_output',['emit','twiss','spotexpected','X','X_unweighted','beta','covar','chisq_red'])
	# out = output(emit,,spotexpected,myfit.X,myfit.X_unweighted,myfit.beta,myfit.covar,myfit.chisq_red)

	out = BeamlineScanFit(fitresults=myfit,spotexpected=spotexpected)

	return out
#}}}
