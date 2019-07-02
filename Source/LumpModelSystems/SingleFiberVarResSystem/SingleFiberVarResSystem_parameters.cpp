/*
 * SingleFiberVarResSystem_parameters.cpp
 *
 *  Created on: Dec 11, 2018
 *      Author: brian
 */


#include "SingleFiberVarResSystem.hpp"


void ParamSingleFiberVarRes::paramPreloopComputations()
{

	dt = Tc/( (double) NtStepsCycle );

	hillMaxVal_ventricle = computeHillMax( m1_ventricle, m2_ventricle, tau1_ventricle, tau2_ventricle, Tc );

	hillMaxVal_atrium = computeHillMax( m1_atrium, m2_atrium, tau1_atrium, tau2_atrium, Tc );


}


void ParamSingleFiberVarRes::setDefaultParams()
{

	//----------- 1. Lump models -----------//
	// Compliances (converted from ml/mmHg to ml/kPa)
	Csa = 15;

	//  Hydraulic resistance (converted from  mmHg*s/ml to kPa*s/ml)
	Rsa = 0.05;
	Rpv = 0.001;
	Rpao = 0.001;

	// fixed values
	Ppv_fixed = 1.0;
	Psp_fixed = 10.66;


	//----------- 2. Single fiber model chambers -----------//
	// Left ventricle
	V0_lv = 50;
	Vw_lv = 100;
	Ta0_lv = 120;
	Tp0_lv = 1.0;
	cp_lv = 8;
	cv_lv = 0;
	v0_lv = 10;
	ls0_lv = 1.8;
	lsmax_lv = 2.3;
	lsw_lv = 0.35;

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


	//----------- 4. Others -----------//

	// Timing and computation
	Tc = 1.0;
	Ncyc = 10;
	NtStepsCycle = 2000;

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


	//----------- 5. Activation model -----------//
	// ventricle activation
	Ts_ventricle = 0.07;
	m1_ventricle = 1.32;
	m2_ventricle = 27.4;
	tau1_ventricle = 0.2;
	tau2_ventricle = 0.4;

	// atrial activation
	Ts_atrium = -0.1;
	m1_atrium = 1.32;
	m2_atrium = 13.1;
	tau1_atrium = 0.07;
	tau2_atrium = 0.12;

	activeDurationAdjustment = 0.0;

	outputFile = "SingleFiberVarRes_output.m";
	paramOutputFile = "SingleFiberVarRes_params.m";
	initialValuesFileName = "Input/initialValuesSingleFiberVarRes.txt";


}


void ParamSingleFiberVarRes::importInitialValues( vector <double> & w,  DataContainer & prmData )
{

	string prmPrefix = "initValuesSingleFiberVarRes";

	// pull values from file
	readPrmVector( "w", prmPrefix, prmData,  w.data(), Nvars );


	MPI_Bcast(w.data(), Nvars, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);


}


void setDefaultInitialValuesSingleFiberVarRes( double * w )
{
	// default initial values satisfy the equilibrium at the default simple lump
	w[0] = 64.7;
	w[1] = 113.5;
	w[2] = 13.6;


}

void ParamSingleFiberVarRes::importParams( DataContainer & prmData )
{
	// sets parameters to defaults
	setDefaultParams();

	string prmPrefix = "SingleFiberVarRes";

	//----------- 1. Lump models -----------//
	Csa = readPrmValue( "Csa", prmPrefix, prmData);
	Rsa = readPrmValue( "Rsa", prmPrefix, prmData);
	Rpv = readPrmValue( "Rpv", prmPrefix, prmData);
	Rpao = readPrmValue( "Rpao", prmPrefix, prmData);

	// fixed compartments
	Ppv_fixed = readPrmValue( "Ppv_fixed", prmPrefix, prmData);
	Psp_fixed = readPrmValue( "Psp_fixed", prmPrefix, prmData);


	//----------- 2. Single fiber model chambers -----------//
	// Left ventricle
	V0_lv = readPrmValue( "V0_lv", prmPrefix, prmData);
	Vw_lv = readPrmValue( "Vw_lv", prmPrefix, prmData);
	Ta0_lv = readPrmValue( "Ta0_lv", prmPrefix, prmData);
	Tp0_lv = readPrmValue( "Tp0_lv", prmPrefix, prmData);
	cp_lv = readPrmValue( "cp_lv", prmPrefix, prmData);
	cv_lv = readPrmValue( "cv_lv", prmPrefix, prmData);
	v0_lv = readPrmValue( "v0_lv", prmPrefix, prmData);
	ls0_lv = readPrmValue( "ls0_lv", prmPrefix, prmData);
	lsmax_lv = readPrmValue( "lsmax_lv", prmPrefix, prmData);
	lsw_lv = readPrmValue( "lsw_lv", prmPrefix, prmData);

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
	NtStepsCycle = readPrmValue( "NtStepsCycle", prmPrefix, prmData);;

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

	// atrial activation
	Ts_atrium = readPrmValue( "Ts_atrium", prmPrefix, prmData);
	m1_atrium = readPrmValue( "m1_atrium", prmPrefix, prmData);
	m2_atrium = readPrmValue( "m2_atrium", prmPrefix, prmData);
	tau1_atrium = readPrmValue( "tau1_atrium", prmPrefix, prmData);
	tau2_atrium = readPrmValue( "tau2_atrium", prmPrefix, prmData);
	activeDurationAdjustment = readPrmValue( "activeDurationAdjustment", prmPrefix, prmData);

	outputFile = readPrmString( "outputFile", prmPrefix, prmData);
	paramOutputFile = readPrmString( "paramOutputFile", prmPrefix, prmData);
	initialValuesFileName = readPrmString( "initialValuesFileName", prmPrefix, prmData);


}
