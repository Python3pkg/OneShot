import mytools.slactrac as _sltr

gamma_default  = 39824
QS1_K1_default = 3.077225846087095e-01
QS2_K1_default = -2.337527121004531e-01

def IP_to_lanex(twiss_x,twiss_y,
		gamma  = gamma_default,
		QS1_K1 = QS1_K1_default,
		QS2_K1 = QS2_K1_default
		):
	# Beamline elements
	IP2QS1    = _sltr.Drift(length = 5.4217)
	QS1       = _sltr.Quad(length = 5.000000000E-01,K1 = QS1_K1)
	LQS12QS2  = _sltr.Drift(length = 4.00E+00)
	QS2       = _sltr.Quad(length = 5.000000000E-01,K1 = QS2_K1)
	LQS22BEND = _sltr.Drift(length = 0.7428E+00)
	B5D36     = _sltr.Bend(
			length= 2*4.889500000E-01,
			angle= 6.0E-03, 	     
			order=1,
			rotate=0
		       )
	LBEND2ELANEX = _sltr.Drift(length = 8.792573)
	beamline     = _sltr.Beamline(
			element_list=[
				IP2QS1        ,
				QS1           ,
				QS1           ,
				LQS12QS2      ,
				QS2           ,
				QS2           ,
				LQS22BEND     ,
				B5D36         ,
				LBEND2ELANEX	
				],
			gamma = gamma,
			twiss_x = twiss_x,
			twiss_y = twiss_y
			)
	return beamline

def IP_to_lanex_nobend(twiss_x,twiss_y,
		gamma=39824,
		QS1_K1 = 3.077225846087095e-01,
		QS2_K1 = -2.337527121004531e-01
		):
	# Beamline elements
	IP2QS1    = _sltr.Drift(length = 5.4217)
	QS1       = _sltr.Quad(length = 5.000000000E-01,K1 = QS1_K1)
	LQS12QS2  = _sltr.Drift(length = 4.00E+00)
	QS2       = _sltr.Quad(length = 5.000000000E-01,K1 = QS2_K1)
	LQS22BEND = _sltr.Drift(length = 0.7428E+00)
	B5D36     = _sltr.Drift(length= 2*4.889500000E-01)
	LBEND2ELANEX = _sltr.Drift(length = 8.792573)
	beamline     = _sltr.Beamline(
			element_list=[
				IP2QS1        ,
				QS1           ,
				QS1           ,
				LQS12QS2      ,
				QS2           ,
				QS2           ,
				LQS22BEND     ,
				B5D36         ,
				LBEND2ELANEX	
				],
			gamma = gamma,
			twiss_x = twiss_x,
			twiss_y = twiss_y
			)
	return beamline

def IP_to_cherfar(twiss_x,twiss_y,
		gamma  = gamma_default,
		QS1_K1 = QS1_K1_default,
		QS2_K1 = QS2_K1_default
		):
	beamline = IP_to_lanex(twiss_x,twiss_y,
		gamma  = gamma,
		QS1_K1 = QS1_K1,
		QS2_K1 = QS2_K1
		)

	beamline.elements[8].length = beamline.elements[8].length + 0.8198

	return beamline
