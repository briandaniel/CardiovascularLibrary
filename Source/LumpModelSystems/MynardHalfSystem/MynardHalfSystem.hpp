/*
 * MynardHalfSystem.hpp
 *
 *  Created on: Oct 8, 2018
 *      Author: brian
 */

#ifndef LUMPMODELSYSTEMS_MYNARDHALFSYSTEM_MYNARDHALFSYSTEM_HPP_
#define LUMPMODELSYSTEMS_MYNARDHALFSYSTEM_MYNARDHALFSYSTEM_HPP_


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
class ParamMynardHalf0D {

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
	double m1_ventricle, m2_ventricle, tau1_ventricle, tau2_ventricle, Ts_ventricle, hillMaxVal_ventricle;
	double m1_atrium, m2_atrium, tau1_atrium, tau2_atrium, Ts_atrium, hillMaxVal_atrium;

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
	ParamMynardHalf0D() // sets default values
	{
		setDefaultParams();
		paramPreloopComputations();
		// the number of variables and output aux variables is fixed in this code:
		setNvars();
	}; // sets default values
	ParamMynardHalf0D(DataContainer & prmData) // sets values according to paramFileName
	{
		importParams( prmData );
		paramPreloopComputations();
		// the number of variables and output aux variables is fixed in this code:
		setNvars();
	};
	~ParamMynardHalf0D(){}; // default destructor (does nothing)
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
	void printParams( string outputFileName );
};

// Mynard0DSystem.cpp
void computeMynardHalfTrapezoidIntegrationStep( double tn, double * w, ParamMynardHalf0D & lumpPrm, double * auxVars );
void computeMynardHalf0DSystem( double * w0, ParamMynardHalf0D & lumpPrm );
void computeMynardHalf0DSystemExport( double * w, double ** wStore, double ** auxStore, double * tStore, ParamMynardHalf0D & lumpPrm );

void computeMynardHalfApproximateJacobian( double * GX, double * X, double dX, double ** J,
		double tnp1, double dt, double * wn, double * Fn, ParamMynardHalf0D & lumpPrm, double * auxVars );
void computeMynardHalfGFunc( double tnp1, double dt, double * wn, double * wnp1_approx, double * Fn, ParamMynardHalf0D & lumpPrm, double * auxVars, double * G );
void computeMynardHalf0DSystemTD( double t, double * w, ParamMynardHalf0D & lumpPrm, double * dwdt, double * auxVars );
void storeValuesMynardHalf0DSystem( int stepIdx, double t, double * w, double * auxVars,
		double * tStore, double ** wStore, double ** auxStore, int Nvars, int NauxVars );

// Mynard0DSystem_parameters.cpp
void setDefaultInitialValuesMynardHalf0D( double * w );

// Mynard0DSystem_output.cpp
void outputModelResultMynardHalf0D( int Nsteps, ParamMynardHalf0D & lumpPrm, double * tStore, double ** wStore, double ** auxStore, double * wEnd );









#endif /* LUMPMODELSYSTEMS_MYNARDHALFSYSTEM_MYNARDHALFSYSTEM_HPP_ */
