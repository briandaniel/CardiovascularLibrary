/*
 * HalfLumpSimpleActivation.hpp
 *
 *  Created on: Nov 13, 2018
 *      Author: brian
 */

#ifndef LUMPMODELSYSTEMS_HALFLUMPSIMPLEACTIVATION_HALFLUMPSIMPLEACTIVATION_HPP_
#define LUMPMODELSYSTEMS_HALFLUMPSIMPLEACTIVATION_HALFLUMPSIMPLEACTIVATION_HPP_



// Parallel commands
#define ROOT_ID 0	 // id of root processor

// other includes
#include <mpi.h>

// local includes
#include "../../ActivationFunctions/temporalActivationFunctions.hpp"
#include "../../LumpModels/LumpModel.hpp"
#include "../../VaryingElastanceModels/varyingElastance.hpp"
#include "../../ValveFunctions/ValveFunctions.hpp"

// Utility library includes
#include "UtilityFunctions/utilityFunctions.hpp"
#include "ImportExport/importExport.hpp"
#include "ParameterParser/ParameterParser.hpp"

// local param class definition
class ParamHalfSimple{

  public:

	//-------- 1. variables --------//
	// lump models
	double Csa, Csp; // capacitance
	double Rsa, Rsp;

	// Fixed pressure, also resistance for pulmonary vein system
	double Ppv_fixed, Psv_fixed;
	double Rpv;

	// varying elastance models
	double Emax_lv, Emin_lv, Ks_lv, V0_lv;
	double Emax_la, Emin_la, Ks_la, V0_la;

	// Activation model
	double Ta_ventricle, Ts_ventricle, kt_ventricle;
	double Ta_atrium, Ts_atrium, kt_atrium;

	// Valve models
	double rho; // same for all models
	double Aeff_max_mv, Aeff_min_mv, leff_mv, Kvo_mv, Kvc_mv, deltaP_open_mv, deltaP_close_mv;
	double Aeff_max_aov, Aeff_min_aov, leff_aov, Kvo_aov, Kvc_aov, deltaP_open_aov, deltaP_close_aov;

	// Timing and computation
	double Tc;
	int Ncyc, Nvars, NauxVars, NtStepsCycle;
	double dt;

	// newton iteration parameters
	double xMinDiff, fNewtMin, dX;
	int maxNewtIter;

	// Output file names
	string outputFile;
	string paramOutputFile;
	string initialValuesFileName;

	//-------- 2. functions --------//
	// simpleLumpSystem.cpp
	ParamHalfSimple() // sets default values
	{
		setDefaultParams();
		paramPreloopComputations();
		// the number of variables and output aux variables is fixed in this code:
		setNvars();
	}; // sets default values
	ParamHalfSimple(DataContainer & prmData ) // sets values according to paramFileName
	{
		importParams( prmData);
		paramPreloopComputations();
		// the number of variables and output aux variables is fixed in this code:
		setNvars();
	};
	~ParamHalfSimple(){}; // default destructor (does nothing)
	void setNvars(){
		Nvars = 8;
		NauxVars = 7;
	};

	// Mynard0DSystem_parameters.cpp
	void paramPreloopComputations();
	void setDefaultParams();
	void importParams( DataContainer & prmData );
	void importInitialValues( double * w,  DataContainer & prmData );

	// Mynard0DSystem_output.cpp
	void printParams(); // prints the params in matlab format

};

// Mynard0DSystem.cpp
void computeHalfSimpleTrapezoidIntegrationStep( double tn, double * w, ParamHalfSimple & lumpPrm, double * auxVars );
void computeHalfSimpleSystem( double * w0, ParamHalfSimple & lumpPrm );
void computeHalfSimpleSystemExport( double * w, double ** wStore, double ** auxStore, double * tStore, ParamHalfSimple & lumpPrm );

void computeHalfSimpleApproximateJacobian( double * GX, double * X, double dX, double ** J,
		double tnp1, double dt, double * wn, double * Fn, ParamHalfSimple & lumpPrm, double * auxVars );
void computeHalfSimpleGFunc( double tnp1, double dt, double * wn, double * wnp1_approx, double * Fn, ParamHalfSimple & lumpPrm, double * auxVars, double * G );
void computeHalfSimpleSystemTD( double t, double * w, ParamHalfSimple & lumpPrm, double * dwdt, double * auxVars );
void storeValuesHalfSimpleSystem( int stepIdx, double t, double * w, double * auxVars,
		double * tStore, double ** wStore, double ** auxStore, int Nvars, int NauxVars );

// Mynard0DSystem_parameters.cpp
void setDefaultInitialValuesHalfSimple( double * w );

// Mynard0DSystem_output.cpp
void outputModelResultHalfSimple( int Nsteps, ParamHalfSimple & lumpPrm, double * tStore, double ** wStore, double ** auxStore, double * wEnd );











#endif /* LUMPMODELSYSTEMS_HALFLUMPSIMPLEACTIVATION_HALFLUMPSIMPLEACTIVATION_HPP_ */
