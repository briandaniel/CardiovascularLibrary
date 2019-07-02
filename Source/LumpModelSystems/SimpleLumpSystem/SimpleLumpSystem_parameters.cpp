/*
 * simpleLumpSystem_parameters.cpp
 *
 *  Created on: Sep 14, 2018
 *      Author: brian
 */


#include "SimpleLumpSystem.hpp"


void ParamSimpleLump::importParams( DataContainer & prmData )
{

	setDefaultParams();

	string prmPrefix = "simpleLump";

	//----------- 1. Lump models -----------//
	Ncyc = (int) readPrmValue( "Ncyc", prmPrefix, prmData);


	Cla = readPrmValue( "Cla", prmPrefix, prmData);
	Cra = readPrmValue( "Cra", prmPrefix, prmData);
	Cpv = readPrmValue( "Cpv", prmPrefix, prmData);
	Cpp = readPrmValue( "Cpp", prmPrefix, prmData);
	Cpa = readPrmValue( "Cpa", prmPrefix, prmData);
	Csv = readPrmValue( "Csv", prmPrefix, prmData);
	Csp = readPrmValue( "Csp", prmPrefix, prmData);
	Csa = readPrmValue( "Csa", prmPrefix, prmData);

	Vu_la = readPrmValue( "Vu_la", prmPrefix, prmData);
	Vu_ra = readPrmValue( "Vu_ra", prmPrefix, prmData);
	Vu_pv = readPrmValue( "Vu_pv", prmPrefix, prmData);
	Vu_pp = readPrmValue( "Vu_pp", prmPrefix, prmData);
	Vu_pa = readPrmValue( "Vu_pa", prmPrefix, prmData);
	Vu_sv = readPrmValue( "Vu_sa", prmPrefix, prmData);
	Vu_sp = readPrmValue( "Vu_sp", prmPrefix, prmData);
	Vu_sa = readPrmValue( "Vu_sa", prmPrefix, prmData);

	Rla = readPrmValue( "Rla", prmPrefix, prmData);
	Rra = readPrmValue( "Rra", prmPrefix, prmData);
	Rpv = readPrmValue( "Rpv", prmPrefix, prmData);
	Rpp = readPrmValue( "Rpp", prmPrefix, prmData);
	Rpa = readPrmValue( "Rpa", prmPrefix, prmData);
	Rsv = readPrmValue( "Rsa", prmPrefix, prmData);
	Rsp = readPrmValue( "Rsp", prmPrefix, prmData);
	Rsa = readPrmValue( "Rsa", prmPrefix, prmData);


	//----------- 2. Varying elastance chambers -----------//

	Vu_rv = readPrmValue( "Vu_rv", prmPrefix, prmData);
	P0_rv = readPrmValue( "P0_rv", prmPrefix, prmData);
	ke_rv = readPrmValue( "ke_rv", prmPrefix, prmData);
	Emax_rv = readPrmValue( "Emax_rv", prmPrefix, prmData);
	kr_rv = readPrmValue( "kr_rv", prmPrefix, prmData);

	Vu_lv = readPrmValue( "Vu_lv", prmPrefix, prmData);
	P0_lv = readPrmValue( "P0_lv", prmPrefix, prmData);
	ke_lv = readPrmValue( "ke_lv", prmPrefix, prmData);
	Emax_lv = readPrmValue( "Emax_lv", prmPrefix, prmData);
	kr_lv = readPrmValue( "kr_lv", prmPrefix, prmData);

	Vtotal = readPrmValue( "Vtotal", prmPrefix, prmData);
	Tc = readPrmValue( "Tc", prmPrefix, prmData);
	Ta = readPrmValue( "Ta", prmPrefix, prmData);
	Ncyc = (int) readPrmValue( "Ncyc", prmPrefix, prmData);
	dt = readPrmValue( "dt", prmPrefix, prmData);


	outputFile = readPrmString( "outputFile", prmPrefix, prmData);
	paramOutputFile = readPrmString( "paramOutputFile", prmPrefix, prmData);
	initialValuesFileName = readPrmString( "initialValuesFileName", prmPrefix, prmData);



}



void ParamSimpleLump::importInitialValues( double * w,  DataContainer & prmData )
{

	string prmPrefix = "initValuesSimpleLump";

	// pull values from file
	readPrmVector( "w", prmPrefix, prmData,  w, Nvars );


	MPI_Bcast(w, Nvars, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);


}


void setDefaultInitialValuesSimpleLump( double * w )
{
	// default initial values satisfy the equilibrium at the default simple lump
	w[0] = 0.6180;
	w[1] = 0.6746;
	w[2] = 1.0947;
	w[3] = 1.1347;
	w[4] = 0.5518;
	w[5] = 8.3615;
	w[6] = 8.3831;
	w[7] = 43.7406;
	w[8] = 60.6320;

}

// sets parameters to defaults
void ParamSimpleLump::setDefaultParams(){

	// Compliances (converted from ml/mmHg to ml/kPa)
	Cla = 230;
	Cpv = 50;
	Cpp = 20;
	Cpa = 5.7;
	Cra = 234.38;
	Csv = 300;
	Csp = 100;
	Csa = 4;

	// Zero pressure volumes (ml)
	Vu_la =   25.0;
	Vu_pv =  120.0;
	Vu_pp =  123.0;
	Vu_pa =    0.0;
	Vu_ra =   25.0;
	Vu_sv = 2496.0;
	Vu_sp =  611.0;
	Vu_sa =    0.0;

	//  Hydraulic resistance (converted from  mmHg*s/ml to kPa*s/ml)
	Rla = 0.000333;
	Rpv = 0.001;
	Rpp = 0.01;
	Rpa = 0.004;
	Rra = 0.000333;
	Rsv = 0.00506;
	Rsp = 0.15;
	Rsa = 0.01;

	// Right ventricle
	Vu_rv = 40.0;
	P0_rv = 0.2;
	ke_rv = 0.013;
	Emax_rv = 0.23;
	kr_rv = 0.0014;

	// Left ventricle
	Vu_lv = 16.77 ;
	P0_lv = 0.2;
	ke_lv = 0.014 ;
	Emax_lv = 0.4 ;
	kr_lv = 0.000375;

	// total blood volume (ml)
	Vtotal = 5000;

	// Timing and computation
	Tc = 1.0;
	Ta = 0.5;
	Ncyc = 5;
	dt = 0.001;


	outputFile = "Output/simpleLumpSystem_output.m";
	paramOutputFile = "Output/simpleLumpSystem_params.m";
	initialValuesFileName = "Input/initialValuesSimpleLumpSystem.txt";

}


