/*
 * lumpModelParam.cpp
 *
 *  Created on: Jan 29, 2019
 *      Author: brian
 */

#include "lumpActiveModelParam.hpp"



void ParamLumpActive::paramPreloopComputations(double Tc)
{

	hillMaxVal_ventricle = computeHillMax( m1_ventricle, m2_ventricle, tau1_ventricle, tau2_ventricle, Tc );

	hillMaxVal_atrium = computeHillMax( m1_atrium, m2_atrium, tau1_atrium, tau2_atrium, Tc );


}


void ParamLumpActive::setDefaultParams()
{
	Nw = 2;

	//----------- 1. Lump models -----------//
	// fixed values
	Ppv_fixed = 1.0;
	Psp_fixed = 10.66;

	// Compliances (converted from ml/mmHg to ml/kPa)
	Csa = 15;

	//  Hydraulic resistance (converted from  mmHg*s/ml to kPa*s/ml)
	Rsa = 0.05;
	Rpv = 0.001;
	Rpao = 0.001;

	// Timing and computation
	Tc = 1.0;
	Ncyc = 10;
	dtMin = 0.001;
	dtMax = 0.001;

	// newton iteration parameters
	xMinDiff = 1e-3;
	fNewtMin = 1e-3;
	dX = 1e-7;
	maxNewtIter = 3;


	// valves
	Rmvo = 0.0001;
	Rmvc = 10;
	Raovo = 0.001;
	Raovc = 10;
	beta = 250;


	// Left atrium
	V0_la = 20;
	Vw_la = 20;
	Ta0_la = 20;
	Tp0_la = 0.5;
	cp_la = 8;
	cv_la = 0;
	v0_la = 10;
	ls0_la = 1.8;
	lsmax_la = 2.3;
	lsw_la = 0.35;


	//----------- 5. Activation model -----------//
	// ventricle activation
	Ts_ventricle = 0.07;
	m1_ventricle = 1.32;
	m2_ventricle = 27.4;
	tau1_ventricle = 0.2;
	tau2_ventricle = 0.4;
	activeDurationAdjustment = 0;

	// atrial activation
	Ts_atrium = -0.00;
	m1_atrium = 1.32;
	m2_atrium = 13.1;
	tau1_atrium = 0.07;
	tau2_atrium = 0.12;


}



void ParamLumpActive::importParams( string prmPrefix, DataContainer & prmData )
{
	// sets parameters to defaults
	setDefaultParams();


	//----------- 1. Lump models -----------//
	Csa = readPrmValue( "Csa", prmPrefix, prmData);
	Rsa = readPrmValue( "Rsa", prmPrefix, prmData);
	Rpv = readPrmValue( "Rpv", prmPrefix, prmData);
	Rpao = readPrmValue( "Rpao", prmPrefix, prmData);

	// fixed compartments
	Ppv_fixed = readPrmValue( "Ppv_fixed", prmPrefix, prmData);
	Psp_fixed = readPrmValue( "Psp_fixed", prmPrefix, prmData);


	//----------- 2. Single fiber model chambers -----------//

	// Left atrium
	V0_la = readPrmValue( "V0_la", prmPrefix, prmData);
	Vw_la = readPrmValue( "Vw_la", prmPrefix, prmData);
	Ta0_la = readPrmValue( "Ta0_la", prmPrefix, prmData);
	Tp0_la = readPrmValue( "Tp0_la", prmPrefix, prmData);
	cp_la = readPrmValue( "cp_la", prmPrefix, prmData);
	cv_la = readPrmValue( "cv_la", prmPrefix, prmData);
	v0_la = readPrmValue( "v0_la", prmPrefix, prmData);
	ls0_la = readPrmValue( "ls0_la", prmPrefix, prmData);
	lsmax_la = readPrmValue( "lsmax_la", prmPrefix, prmData);
	lsw_la = readPrmValue( "lsw_la", prmPrefix, prmData);



	Rmvo = readPrmValue( "Rmvo", prmPrefix, prmData );
	Rmvc = readPrmValue( "Rmvc", prmPrefix, prmData );
	Raovo = readPrmValue( "Raovo", prmPrefix, prmData );
	Raovc = readPrmValue( "Raovc", prmPrefix, prmData );
	beta = readPrmValue( "beta", prmPrefix, prmData );


	//----------- 4. Others -----------//
	Tc = readPrmValue( "Tc", prmPrefix, prmData);
	Ncyc = (int) readPrmValue( "Ncyc", prmPrefix, prmData);
	dtMin = readPrmValue( "dtMin", prmPrefix, prmData);
	dtMax = readPrmValue( "dtMax", prmPrefix, prmData);


	// newton iteration parameters
	xMinDiff = readPrmValue( "xMinDiff", prmPrefix, prmData);
	fNewtMin = readPrmValue( "fNewtMin", prmPrefix, prmData);
	dX = readPrmValue( "dX", prmPrefix, prmData);
	maxNewtIter = readPrmValue( "maxNewtIter", prmPrefix, prmData);


	//----------- 5. Activation model -----------//

	// ventricle activation
	Ts_ventricle = readPrmValue( "Ts_ventricle", prmPrefix, prmData);
	m1_ventricle = readPrmValue( "m1_ventricle", prmPrefix, prmData);
	m2_ventricle = readPrmValue( "m2_ventricle", prmPrefix, prmData);
	tau1_ventricle = readPrmValue( "tau1_ventricle", prmPrefix, prmData);
	tau2_ventricle = readPrmValue( "tau2_ventricle", prmPrefix, prmData);

	activeDurationAdjustment = readPrmValue( "activeDurationAdjustment", prmPrefix, prmData);

	// atrial activation
	Ts_atrium = readPrmValue( "Ts_atrium", prmPrefix, prmData);
	m1_atrium = readPrmValue( "m1_atrium", prmPrefix, prmData);
	m2_atrium = readPrmValue( "m2_atrium", prmPrefix, prmData);
	tau1_atrium = readPrmValue( "tau1_atrium", prmPrefix, prmData);
	tau2_atrium = readPrmValue( "tau2_atrium", prmPrefix, prmData);


}
