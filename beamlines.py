import mytools.slactrac as _sltr

gamma_default  = 39824
# QS1_K1_default = 3.077225846087095e-01
# QS2_K1_default = -2.337527121004531e-01
QS1_K1_default = 3.8743331090707228e-1
QS2_K1_default = -2.5439067538354171e-1
PEXT_Z = 1994.97
QS1_Z = 1998.71
# IP2QS1_length = 5.4217
IP2QS1_length = QS1_Z-PEXT_Z

def IP_to_lanex(beam_x,beam_y,
		gamma  = gamma_default,
		QS1_K1 = QS1_K1_default,
		QS2_K1 = QS2_K1_default
		):
	# Beamline elements
	IP2QS1    = _sltr.Drift(length = IP2QS1_length)
	QS1       = _sltr.Quad(length = 5.000000000E-01,K1 = QS1_K1)
	LQS12QS2  = _sltr.Drift(length = 4.00E+00)
	QS2       = _sltr.Quad(length = 5.000000000E-01,K1 = QS2_K1)
	LQS22BEND = _sltr.Drift(length = 0.7428E+00)
	B5D36     = _sltr.Bend(
			length= 2*4.889500000E-01,
			angle= 6.0E-03, 	     
			order=1,
			rotate=90
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
			beam_x = beam_x,
			beam_y = beam_y
			)
	return beamline

def IP_to_lanex_nobend(beam_x,beam_y,
		gamma  = gamma_default,
		QS1_K1 = QS1_K1_default,
		QS2_K1 = QS2_K1_default
		):

	beamline = IP_to_lanex(
			beam_x=beam_x,
			beam_y=beam_y,
			gamma  = gamma,
			QS1_K1 = QS1_K1,
			QS2_K1 = QS2_K1
			)

	# Replace bend with drift
	B5D36_drift = _sltr.Drift(length= 2*4.889500000E-01)
	beamline.elements[7] = B5D36_drift

	return beamline

def IP_to_cherfar(beam_x,beam_y,
		gamma  = gamma_default,
		QS1_K1 = QS1_K1_default,
		QS2_K1 = QS2_K1_default
		):
	beamline = IP_to_lanex(beam_x,beam_y,
		gamma  = gamma,
		QS1_K1 = QS1_K1,
		QS2_K1 = QS2_K1
		)

	beamline.elements[8].length = beamline.elements[8].length + 0.8198

	return beamline
