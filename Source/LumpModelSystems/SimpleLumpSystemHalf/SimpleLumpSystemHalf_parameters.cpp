/*
 * SimpleLumpSystemHalf_parameters.cpp
 *
 *  Created on: Sep 14, 2018
 *      Author: brian
 */


#include "SimpleLumpSystemHalf.hpp"


void ParamSimpleLumpHalf::paramPreloopComputations()
{

	dt = Tc/( (double) NtStepsCycle );

	hillMaxVal_ventricle = computeHillMax( m1_ventricle, m2_ventricle, tau1_ventricle, tau2_ventricle, Tc );

	hillMaxVal_atrium = computeHillMax( m1_atrium, m2_atrium, tau1_atrium, tau2_atrium, Tc );


}


void ParamSimpleLumpHalf::setDefaultParams()
{

	//----------- 1. Lump models -----------//
	// Compliances (converted from ml/mmHg to ml/kPa)
	Csa = 4;

	//  Hydraulic resistance (converted from  mmHg*s/ml to kPa*s/ml)
	Rsa = 0.01;
	Rpv = 0.001;
	Rpao = 0.001;

	// fixed values
	Ppv_fixed = 1.0;
	Psp_fixed = 10.66;



	//----------- 2. Varying elastance chambers -----------//
	// Left ventricle
	Emin_lv = 0.0107;
	Emax_lv = 0.4;
	Ks_lv = 4e-9;
	V0_lv = 20;

	// Left atrium
	Emin_la = 0.0107;
	Emax_la = 0.0227;
	Ks_la = 10e-9;
	V0_la = 3;

	//----------- 4. Others -----------//

	// Timing and computation
	Tc = 1.0;
	Ncyc = 5;
	NtStepsCycle = 1000;

	// newton iteration parameters
	xMinDiff = 1e-3;
	fNewtMin = 1e-3;
	dX = 1e-7;
	maxNewtIter = 3;


	// valves
	Rmvo = 0.001;
	Rmvc = 10;
	Raovo = 0.001;
	Raovc = 10;
	beta = 100;


	//----------- 5. Activation model -----------//
	// ventricle activation
	m1_ventricle = 1.32;
	m2_ventricle = 27.4;
	tau1_ventricle = 0.269;
	tau2_ventricle = 0.452;
	Ts_ventricle = 0.1;

	// atrial activation
	Ts_atrium = -0.15;
	m1_atrium = 1.32;
	m2_atrium = 13.1;
	tau1_atrium = 0.11;
	tau2_atrium = 0.18;


	outputFile = "Output/simpleLumpSystemHalf_output.m";
	paramOutputFile = "Output/SimpleLumpSystemHalf_params.m";
	initialValuesFileName = "Input/initialValuesSimpleLumpSystemHalf.txt";

}


void ParamSimpleLumpHalf::importInitialValues( vector <double> & w,  DataContainer & prmData )
{

	string prmPrefix = "initValuesSimpleLumpHalf";

	// pull values from file
	readPrmVector( "w", prmPrefix, prmData,  w.data(), Nvars );


	MPI_Bcast(w.data(), Nvars, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);


}


void setDefaultInitialValuesSimpleLumpHalf( double * w )
{
	// default initial values satisfy the equilibrium at the default simple lump
	w[0] = 20.0;
	w[1] = 50.0;
	w[2] = 11.0;


}

void ParamSimpleLumpHalf::importParams( DataContainer & prmData )
{
	// sets parameters to defaults
	setDefaultParams();

	string prmPrefix = "SimpleLumpHalf";

	//----------- 1. Lump models -----------//
	Csa = readPrmValue( "Csa", prmPrefix, prmData);
	Rsa = readPrmValue( "Rsa", prmPrefix, prmData);
	Rpv = readPrmValue( "Rpv", prmPrefix, prmData);
	Rpao = readPrmValue( "Rpao", prmPrefix, prmData);

	// fixed compartments
	Ppv_fixed = readPrmValue( "Ppv_fixed", prmPrefix, prmData);
	Psp_fixed = readPrmValue( "Psp_fixed", prmPrefix, prmData);


	//----------- 2. Varying elastance chambers -----------//
	// Left ventricle
	Emin_lv = readPrmValue( "Emin_lv", prmPrefix, prmData);
	Emax_lv = readPrmValue( "Emax_lv", prmPrefix, prmData);
	Ks_lv = readPrmValue( "Ks_lv", prmPrefix, prmData);
	V0_lv = readPrmValue( "V0_lv", prmPrefix, prmData);

	// Left atrium
	Emin_la = readPrmValue( "Emin_la", prmPrefix, prmData);
	Emax_la = readPrmValue( "Emax_la", prmPrefix, prmData);
	Ks_la = readPrmValue( "Ks_la", prmPrefix, prmData);
	V0_la = readPrmValue( "V0_la", prmPrefix, prmData);


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


	outputFile = readPrmString( "outputFile", prmPrefix, prmData);
	paramOutputFile = readPrmString( "paramOutputFile", prmPrefix, prmData);
	initialValuesFileName = readPrmString( "initialValuesFileName", prmPrefix, prmData);


}


