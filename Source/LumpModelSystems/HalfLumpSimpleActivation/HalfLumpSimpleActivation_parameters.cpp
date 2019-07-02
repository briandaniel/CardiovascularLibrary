/*
 * HalfLumpSimpleActivation_parameters.cpp
 *
 *  Created on: Nov 13, 2018
 *      Author: brian
 */



#include "HalfLumpSimpleActivation.hpp"



void ParamHalfSimple::paramPreloopComputations()
{
	dt = Tc/( (double) NtStepsCycle );
}



void ParamHalfSimple::importParams( DataContainer & prmData )
{

	setDefaultParams();

	string prmPrefix = "halfSimp";

	//----------- 1. Lump models -----------//
	Csa = readPrmValue( "Csa", prmPrefix, prmData);
	Csp = readPrmValue( "Csp", prmPrefix, prmData);

	Rsp = readPrmValue( "Rsp", prmPrefix, prmData);
	Rsa = readPrmValue( "Rsa", prmPrefix, prmData);

	// fixed compartments
	Ppv_fixed = readPrmValue( "Ppv_fixed", prmPrefix, prmData);
	Psv_fixed = readPrmValue( "Psv_fixed", prmPrefix, prmData);
	Rpv = readPrmValue( "Rpv", prmPrefix, prmData);

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

	//----------- 3. Time delay valves -----------//
	// rho is Converted from g/cm^3 to kg/cm^3 for consistency with kPa
	rho = readPrmValue( "rho", prmPrefix, prmData);

	// mitral valve
	Aeff_min_mv = readPrmValue( "Aeff_min_mv", prmPrefix, prmData);
	Aeff_max_mv = readPrmValue( "Aeff_max_mv", prmPrefix, prmData);
	leff_mv = readPrmValue( "leff_mv", prmPrefix, prmData);
	// opening rates converted from 1/( s * dyne/cm^2 ) to 1/(s*kPa)
	Kvo_mv = readPrmValue( "Kvo_mv", prmPrefix, prmData);
	Kvc_mv = readPrmValue( "Kvc_mv", prmPrefix, prmData);
	deltaP_open_mv = readPrmValue( "deltaP_open_mv", prmPrefix, prmData);
	deltaP_close_mv = readPrmValue( "deltaP_close_mv", prmPrefix, prmData);

	// aortic valve
	Aeff_min_aov = readPrmValue( "Aeff_min_aov", prmPrefix, prmData);
	Aeff_max_aov = readPrmValue( "Aeff_max_aov", prmPrefix, prmData);
	leff_aov = readPrmValue( "leff_aov", prmPrefix, prmData);
	// opening rates converted from 1/( s * dyne/cm^2 ) to 1/(s*kPa)
	Kvo_aov = readPrmValue( "Kvo_aov", prmPrefix, prmData);
	Kvc_aov = readPrmValue( "Kvc_aov", prmPrefix, prmData);
	deltaP_open_aov = readPrmValue( "deltaP_open_aov", prmPrefix, prmData);
	deltaP_close_aov = readPrmValue( "deltaP_close_aov", prmPrefix, prmData);

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
	Ta_ventricle = readPrmValue( "Ta_ventricle", prmPrefix, prmData);
	kt_ventricle = readPrmValue( "kt_ventricle", prmPrefix, prmData);

	// atrial activation
	Ts_atrium = readPrmValue( "Ts_atrium", prmPrefix, prmData);
	Ta_atrium = readPrmValue( "Ta_atrium", prmPrefix, prmData);
	kt_atrium = readPrmValue( "kt_atrium", prmPrefix, prmData);


	outputFile = readPrmString( "outputFile", prmPrefix, prmData);
	paramOutputFile = readPrmString( "paramOutputFile", prmPrefix, prmData);
	initialValuesFileName = readPrmString( "initialValuesFileName", prmPrefix, prmData);



}



void ParamHalfSimple::importInitialValues( double * w,  DataContainer & prmData )
{

	string prmPrefix = "initValsHalfSimp";

	// pull values from file
	readPrmVector( "w", prmPrefix, prmData,  w, Nvars );


	MPI_Bcast(w, Nvars, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);


}

void setDefaultInitialValuesHalfSimple( double * w )
{
	// default initial values satisfy the equilibrium at the default simple lump
	w[0] = 60.6685;
	w[1] = 88.0113;
	w[2] = 219.3750;
	w[3] = -0.1462;
	w[4] = 0.9999;
	w[5] = 0.0;
	w[6] = 11.7496;
	w[7] = 11.7177;

}

// sets parameters to defaults
void ParamHalfSimple::setDefaultParams(){

	//----------- 1. Lump models -----------//
	// Compliances (converted from ml/mmHg to ml/kPa)
	Csp = 100;
	Csa = 4;

	//  Hydraulic resistance (converted from  mmHg*s/ml to kPa*s/ml)
	Rsp = 0.15;
	Rsa = 0.01;

	// fixed values
	Ppv_fixed = 1.0;
	Psv_fixed = 1.0;
	Rpv = 0.001;


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


	//----------- 3. Time delay valves -----------//
	// rho is Converted from g/cm^3 to kg/cm^3 for consistency with kPa
	rho = 1.06e-3;

	// mitral valve
	Aeff_min_mv = 0.001;
	Aeff_max_mv = 7.7;
	leff_mv = 0.8;
	// opening rates converted from 1/( s * dyne/cm^2 ) to 1/(s*kPa)
	Kvo_mv = 0.02e4;
	Kvc_mv = 0.04e4;
	deltaP_open_mv = 0.01;
	deltaP_close_mv = 0.01;

	// aortic valve
	Aeff_min_aov = 0.001;
	Aeff_max_aov = 6.9;
	leff_aov = 1.5;
	// opening rates converted from 1/( s * dyne/cm^2 ) to 1/(s*kPa)
	Kvo_aov = 0.012e4;
	Kvc_aov = 0.012e4;
	deltaP_open_aov = 0.01;
	deltaP_close_aov = 0.01;


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

	//----------- 5. Activation model -----------//
	// ventricle activation
	Ta_ventricle = 0.4;
	Ts_ventricle = 0.1;
	kt_ventricle = -3;

	// atrial activation
	Ta_atrium = 0.25;
	Ts_atrium = -0.15;
	kt_atrium = -3;

	outputFile = "Output/HalfSimple_output.m";
	paramOutputFile = "Output/HalfSimple_params.m";
	initialValuesFileName = "Input/initialValuesHalfSimpleSystem.txt";

}


