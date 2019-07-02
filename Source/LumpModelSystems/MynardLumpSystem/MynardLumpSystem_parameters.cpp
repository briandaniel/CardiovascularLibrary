/*
 * MynardLumpSystem_parameters.cpp
 *
 *  Created on: Sep 14, 2018
 *      Author: brian
 */


#include "MynardLumpSystem.hpp"

void ParamMynard0D::paramPreloopComputations()
{

	Ts_ventricle = Tc+ventricleActivationOffset ;
	tau1_ventricle = tau1_multiplier_ventricle*Tc;
	tau2_ventricle = tau2_multiplier_ventricle*Tc;
	hillMaxVal_ventricle = computeHillMax( m1_ventricle, m2_ventricle, tau1_ventricle, tau2_ventricle, Tc );

	Ts_atrium = Tc+ventricleActivationOffset+atrialActivationOffset;
	tau1_atrium = tau1_multiplier_atrium*Tc;
	tau2_atrium = tau2_multiplier_atrium*Tc;
	hillMaxVal_atrium = computeHillMax( m1_atrium, m2_atrium, tau1_atrium, tau2_atrium, Tc );

}

void ParamMynard0D::importParams( DataContainer & prmData )
{

	setDefaultParams();

	string prmPrefix = "mynard0D";

	//----------- 1. Lump models -----------//

	Cpv = readPrmValue( "Cpv", prmPrefix, prmData);
	Cpp = readPrmValue( "Cpp", prmPrefix, prmData);
	Cpa = readPrmValue( "Cpa", prmPrefix, prmData);
	Csv = readPrmValue( "Csv", prmPrefix, prmData);
	Csp = readPrmValue( "Csp", prmPrefix, prmData);
	Csa = readPrmValue( "Csa", prmPrefix, prmData);

	Vu_pv = readPrmValue( "Vu_pv", prmPrefix, prmData);
	Vu_pp = readPrmValue( "Vu_pp", prmPrefix, prmData);
	Vu_pa = readPrmValue( "Vu_pa", prmPrefix, prmData);
	Vu_sv = readPrmValue( "Vu_sa", prmPrefix, prmData);
	Vu_sp = readPrmValue( "Vu_sp", prmPrefix, prmData);
	Vu_sa = readPrmValue( "Vu_sa", prmPrefix, prmData);

	Rpv = readPrmValue( "Rpv", prmPrefix, prmData);
	Rpp = readPrmValue( "Rpp", prmPrefix, prmData);
	Rpa = readPrmValue( "Rpa", prmPrefix, prmData);
	Rsv = readPrmValue( "Rsa", prmPrefix, prmData);
	Rsp = readPrmValue( "Rsp", prmPrefix, prmData);
	Rsa = readPrmValue( "Rsa", prmPrefix, prmData);


	//----------- 2. Varying elastance chambers -----------//
	// Left ventricle
	Emin_lv = readPrmValue( "Emin_lv", prmPrefix, prmData);
	Emax_lv = readPrmValue( "Emax_lv", prmPrefix, prmData);
	Ks_lv = readPrmValue( "Ks_lv", prmPrefix, prmData);
	V0_lv = readPrmValue( "V0_lv", prmPrefix, prmData);

	// 	Right ventricle
	Emin_rv = readPrmValue( "Emin_rv", prmPrefix, prmData);
	Emax_rv = readPrmValue( "Emax_rv", prmPrefix, prmData);
	Ks_rv = readPrmValue( "Ks_rv", prmPrefix, prmData);
	V0_rv = readPrmValue( "V0_rv", prmPrefix, prmData);

	// Left atrium
	Emin_la = readPrmValue( "Emin_la", prmPrefix, prmData);
	Emax_la = readPrmValue( "Emax_la", prmPrefix, prmData);
	Ks_la = readPrmValue( "Ks_la", prmPrefix, prmData);
	V0_la = readPrmValue( "V0_la", prmPrefix, prmData);

	// Right atrium
	Emin_ra = readPrmValue( "Emin_ra", prmPrefix, prmData);
	Emax_ra = readPrmValue( "Emax_ra", prmPrefix, prmData);
	Ks_ra = readPrmValue( "Ks_ra", prmPrefix, prmData);
	V0_ra = readPrmValue( "V0_ra", prmPrefix, prmData);


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

	// tricuspid valve
	Aeff_min_tri = readPrmValue( "Aeff_min_tri", prmPrefix, prmData);
	Aeff_max_tri = readPrmValue( "Aeff_max_tri", prmPrefix, prmData);
	leff_tri = readPrmValue( "leff_tri", prmPrefix, prmData);
	// opening rates converted from 1/( s * dyne/cm^2 ) to 1/(s*kPa)
	Kvo_tri = readPrmValue( "Kvo_tri", prmPrefix, prmData);
	Kvc_tri = readPrmValue( "Kvc_tri", prmPrefix, prmData);
	deltaP_open_tri = readPrmValue( "deltaP_open_tri", prmPrefix, prmData);
	deltaP_close_tri = readPrmValue( "deltaP_close_tri", prmPrefix, prmData);

	// pulmonary valve
	Aeff_min_pul = readPrmValue( "Aeff_min_pul", prmPrefix, prmData);
	Aeff_max_pul = readPrmValue( "Aeff_max_pul", prmPrefix, prmData);
	leff_pul = readPrmValue( "leff_pul", prmPrefix, prmData);
	// opening rates converted from 1/( s * dyne/cm^2 ) to 1/(s*kPa)
	Kvo_pul = readPrmValue( "Kvo_pul", prmPrefix, prmData);
	Kvc_pul = readPrmValue( "Kvc_pul", prmPrefix, prmData);
	deltaP_open_pul = readPrmValue( "deltaP_open_pul", prmPrefix, prmData);
	deltaP_close_pul = readPrmValue( "deltaP_close_pul", prmPrefix, prmData);


	//----------- 4. Others -----------//
	Vtotal = readPrmValue( "Vtotal", prmPrefix, prmData);
	Tc = readPrmValue( "Tc", prmPrefix, prmData);
	Ncyc = (int) readPrmValue( "Ncyc", prmPrefix, prmData);
	dt = readPrmValue( "dt", prmPrefix, prmData);

	// newton iteration parameters
	xMinDiff = readPrmValue( "xMinDiff", prmPrefix, prmData);
	fNewtMin = readPrmValue( "fNewtMin", prmPrefix, prmData);
	dX = readPrmValue( "dX", prmPrefix, prmData);
	maxNewtIter = readPrmValue( "maxNewtIter", prmPrefix, prmData);


	//----------- 5. Activation model -----------//

	// ventricle activation
	ventricleActivationOffset = readPrmValue( "ventricleActivationOffset", prmPrefix, prmData);
	m1_ventricle = readPrmValue( "m1_ventricle", prmPrefix, prmData);
	m2_ventricle = readPrmValue( "m2_ventricle", prmPrefix, prmData);
	tau1_multiplier_ventricle = readPrmValue( "tau1_multiplier_ventricle", prmPrefix, prmData);
	tau2_multiplier_ventricle = readPrmValue( "tau2_multiplier_ventricle", prmPrefix, prmData);

	// atrial activation
	atrialActivationOffset = readPrmValue( "atrialActivationOffset", prmPrefix, prmData);
	m1_atrium = readPrmValue( "m1_atrium", prmPrefix, prmData);
	m2_atrium = readPrmValue( "m2_atrium", prmPrefix, prmData);
	tau1_multiplier_atrium = readPrmValue( "tau1_multiplier_atrium", prmPrefix, prmData);
	tau2_multiplier_atrium = readPrmValue( "tau2_multiplier_atrium", prmPrefix, prmData);


	outputFile = readPrmString( "outputFile", prmPrefix, prmData);
	paramOutputFile = readPrmString( "paramOutputFile", prmPrefix, prmData);
	initialValuesFileName = readPrmString( "initialValuesFileName", prmPrefix, prmData);



}



void ParamMynard0D::importInitialValues( double * w,  DataContainer & prmData )
{

	string prmPrefix = "initValsMynard0D";

	// pull values from file
	readPrmVector( "w", prmPrefix, prmData,  w, Nvars );


	MPI_Bcast(w, Nvars, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);


}

void setDefaultInitialValuesMynard0D( double * w )
{
	// default initial values satisfy the equilibrium at the default simple lump
	w[0] = 59.306995;
	w[1] = 57.175841;
	w[2] = 90.807795;
	w[3] = 99.309165;
	w[4] = 234.426;
	w[5] = -0.13865024;
	w[6] = 149.42782;
	w[7] = -0.052826727;
	w[8] = 0.99997065;
	w[9] = 0;
	w[10] = 0.99968242;
	w[11] = 0;
	w[12] = 1.0471507;
	w[13] = 1.6617774;
	w[14] = 1.7190606;
	w[15] = 10.559152;
	w[16] = 10.588918;

}

// sets parameters to defaults
void ParamMynard0D::setDefaultParams(){

	//----------- 1. Lump models -----------//
	// Compliances (converted from ml/mmHg to ml/kPa)
	Cpv = 50;
	Cpp = 20;
	Cpa = 5.7;
	Csv = 300;
	Csp = 100;
	Csa = 4;

	// Zero pressure volumes (ml)
	Vu_pv =  120.0;
	Vu_pp =  123.0;
	Vu_pa =    0.0;
	Vu_sv = 2496.0;
	Vu_sp =  611.0;
	Vu_sa =    0.0;

	//  Hydraulic resistance (converted from  mmHg*s/ml to kPa*s/ml)
	Rpv = 0.001;
	Rpp = 0.01;
	Rpa = 0.004;
	Rsv = 0.00506;
	Rsp = 0.15;
	Rsa = 0.01;


	//----------- 2. Varying elastance chambers -----------//
	// Left ventricle
	Emin_lv = 0.0107;
	Emax_lv = 0.4;
	Ks_lv = 4e-9;
	V0_lv = 20;

	// 	Right ventricle
	Emin_rv = 0.0053;
	Emax_rv = 0.08;
	Ks_rv = 20e-9;
	V0_rv = 20;

	// Left atrium
	Emin_la = 0.0107;
	Emax_la = 0.0227;
	Ks_la = 10e-9;
	V0_la = 3;

	// Right atrium
	Emin_ra = 0.0053;
	Emax_ra = 0.02;
	Ks_ra = 10e-9;
	V0_ra = 17;


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

	// tricuspid valve
	Aeff_min_tri = 0.001;
	Aeff_max_tri = 8.0;
	leff_tri = 0.8;
	// opening rates converted from 1/( s * dyne/cm^2 ) to 1/(s*kPa)
	Kvo_tri = 0.03e4;
	Kvc_tri = 0.04e4;
	deltaP_open_tri = 0.01;
	deltaP_close_tri = 0.01;

	// pulmonary valve
	Aeff_min_pul = 0.001;
	Aeff_max_pul = 7.1;
	leff_pul = 1.0;
	// opening rates converted from 1/( s * dyne/cm^2 ) to 1/(s*kPa)
	Kvo_pul = 0.02e4;
	Kvc_pul = 0.02e4;
	deltaP_open_pul = 0.01;
	deltaP_close_pul = 0.01;

	//----------- 4. Others -----------//

	// total blood volume (ml)
	Vtotal = 5000;
	// Timing and computation
	Tc = 1.0;
	Ncyc = 5;
	dt = 0.001;

	// newton iteration parameters
	xMinDiff = 1e-3;
	fNewtMin = 1e-3;
	dX = 1e-7;
	maxNewtIter = 3;

	//----------- 5. Activation model -----------//
	// ventricle activation
	ventricleActivationOffset = 0.4;
	m1_ventricle = 1.32;
	m2_ventricle = 27.4;
	tau1_multiplier_ventricle = 0.269;
	tau2_multiplier_ventricle = 0.452;

	// atrial activation
	atrialActivationOffset = -0.15;
	m1_atrium = 1.32;
	m2_atrium = 13.1;
	tau1_multiplier_atrium = 0.11;
	tau2_multiplier_atrium = 0.18;

	outputFile = "Output/Mynard0DSystem_output.m";
	paramOutputFile = "Output/Mynard0DSystem_params.m";
	initialValuesFileName = "Input/initialValuesMynard0DSystem.txt";

}


