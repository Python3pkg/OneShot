#!/usr/bin/env python
import numpy as np

# Fit bowtie {{{
def fitbowtie(beamline,davg,stddevsq,filt,T,twiss,emitx):
	betax=twiss[0]
	alphax=twiss[1]
	gammax=twiss[2]
	y=stddevsq[filt]
	y=y[np.newaxis]
	gamma = (1+davg[filt])*40000
	X            = np.zeros([len(gamma),3])
	# spot         = np.zeros(len(gamma))
	spotexpected = np.zeros(len(gamma))
	R = np.zeros([6,6,len(gamma)])
	betaf = np.zeros(len(gamma))
	
	for i,g in enumerate(gamma):
		beamline.change_energy(g)
		beamline.calc_mat()
		R11 = beamline.R[0,0]
		R12 = beamline.R[0,1]
		R[:,:,i] = beamline.R
		X[i,0] = R11*R11
		X[i,1] = 2*R11*R12
		X[i,2] = R12*R12
		T2 = np.dot(np.dot(R[0:2,0:2,i],T),np.transpose(R[0:2,0:2,i]))
		betaf[i] = T2[0,0]
		spotexpected[i] = np.sqrt((R11*R11*betax - 2*R11*R12*alphax + R12*R12*gammax)*emitx)
	
	beta = np.dot(np.linalg.pinv(X) , y.transpose())
	# print beta
	
	emit = np.sqrt( beta[0,0] * beta[2,0] - np.square(beta[1,0]) )
	beta0 = beta[0,0]/emit
	gamma0 = beta[2,0]/emit
	alpha0 = -np.sign(beta[1,0])*np.sqrt(beta0*gamma0-1)
	print 'Emittance fit:\t{}.'.format(emit)
	print 'Normalized emittance fit:\t{}.'.format(emit*40000)
	print 'Initial beta fit:\t{}.'.format(beta0)
	print 'Initial alpha fit:\t{}.'.format(alpha0)
	print 'Initial gamma fit:\t{}.'.format(gamma0)
	print 'Initial spot from fit:\t{}.'.format(np.sqrt(beta[0,0]))
	out = [spotexpected,X,beta]
	return out
#}}}
