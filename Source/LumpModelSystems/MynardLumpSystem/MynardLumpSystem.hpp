/*
 * simpleLumpSystem.hpp
 *
 *  Created on: Sep 13, 2018
 *      Author: brian
 */

#ifndef LUMPMODELSYSTEMS_MYNARDLUMPSYSTEM_MYNARD0DSYSTEM_HPP_
#define LUMPMODELSYSTEMS_MYNARDLUMPSYSTEM_MYNARD0DSYSTEM_HPP_

// Parallel commands
#define ROOT_ID 0	 // id of root processor

// other includes
#include <mpi.h>

// local includes
#include "../../ActivationFunctions/temporalActivationFunctions.hpp"
#include "../../LumpModels/LumpModel.hpp"
#include "../../VaryingElastanceModels/varyingElastance.hpp"
#include "../../ValveFunctions/ValveFunctions.hpp"

#include "UtilityFunctions/utilityFunctions.hpp"
#include "ImportExport/importExport.hpp"
#include "ParameterParser/ParameterParser.hpp"

// local param class definition
class ParamMynard0D {

  public:

	//-------- 1. variables --------//
	// lump models
	double Cpv, Cpp, Cpa, Csv, Csp, Csa; // capacitance
	double Vu_pv, Vu_pp, Vu_pa, Vu_sv, Vu_sp, Vu_sa; // unstressed volumes
	double Rpv, Rpp, Rpa, Rsv, Rsp, Rsa;

	// varying elastance models
	double Emax_lv, Emin_lv, Ks_lv, V0_lv;
	double Emax_rv, Emin_rv, Ks_rv, V0_rv;
	double Emax_la, Emin_la, Ks_la, V0_la;
	double Emax_ra, Emin_ra, Ks_ra, V0_ra;

	// Activation model	tau1_ventricle_multiplier = 0.269;
	double atrialActivationOffset, ventricleActivationOffset;
	double tau1_multiplier_ventricle, tau2_multiplier_ventricle, tau1_multiplier_atrium, tau2_multiplier_atrium;
	double m1_ventricle, m2_ventricle, tau1_ventricle, tau2_ventricle, Ts_ventricle, hillMaxVal_ventricle;
	double m1_atrium, m2_atrium, tau1_atrium, tau2_atrium, Ts_atrium, hillMaxVal_atrium;

	// Valve models
	double rho; // same for all models
	double Aeff_max_mv, Aeff_min_mv, leff_mv, Kvo_mv, Kvc_mv, deltaP_open_mv, deltaP_close_mv;
	double Aeff_max_aov, Aeff_min_aov, leff_aov, Kvo_aov, Kvc_aov, deltaP_open_aov, deltaP_close_aov;
	double Aeff_max_tri, Aeff_min_tri, leff_tri, Kvo_tri, Kvc_tri, deltaP_open_tri, deltaP_close_tri;
	double Aeff_max_pul, Aeff_min_pul, leff_pul, Kvo_pul, Kvc_pul, deltaP_open_pul, deltaP_close_pul;

	// total blood volume
	double Vtotal;

	// Timing and computation
	double Tc;
	int Ncyc, Nvars, NauxVars;
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
	ParamMynard0D() // sets default values
	{
		setDefaultParams();
		paramPreloopComputations();
		// the number of variables and output aux variables is fixed in this code:
		setNvars();
	}; // sets default values
	ParamMynard0D(DataContainer & prmData) // sets values according to paramFileName
	{
		importParams( prmData );
		paramPreloopComputations();
		// the number of variables and output aux variables is fixed in this code:
		setNvars();
	};
	~ParamMynard0D(){}; // default destructor (does nothing)
	void setNvars(){
		Nvars = 17; // 4 + 2*4 + 5 = varying elastance variables + valve variables + lump variables
		NauxVars = 7;
	};

	// Mynard0DSystem_parameters.cpp
	void paramPreloopComputations();
	void setDefaultParams();
	void importParams( DataContainer & prmData );
	void importInitialValues( double * w, DataContainer & prmData );

	// Mynard0DSystem_output.cpp
	void printParams(); // prints the params in matlab format

};

// Mynard0DSystem.cpp
void computeMynardTrapezoidIntegrationStep( double tn, double * w, ParamMynard0D & lumpPrm, double * auxVars );
void computeMynard0DSystem( ParamMynard0D & lumpPrm );
void computeMynardApproximateJacobian( double * GX, double * X, double dX, double ** J,
		double tnp1, double dt, double * wn, double * Fn, ParamMynard0D & lumpPrm, double * auxVars );
void computeMynardGFunc( double tnp1, double dt, double * wn, double * wnp1_approx, double * Fn, ParamMynard0D & lumpPrm, double * auxVars, double * G );
void computeMynard0DSystemTD( double t, double * w, ParamMynard0D & lumpPrm, double * dwdt, double * auxVars );
void storeValuesMynard0DSystem( int stepIdx, double t, double * w, double * auxVars,
		double * tStore, double ** wStore, double ** auxStore, int Nvars, int NauxVars );

// Mynard0DSystem_parameters.cpp
void setDefaultInitialValuesMynard0D( double * w );

// Mynard0DSystem_output.cpp
void outputModelResultMynard0D( int Nsteps, ParamMynard0D & lumpPrm, double * tStore, double ** wStore, double ** auxStore, double * wEnd );








#endif /* LUMPMODELSYSTEMS_MYNARDLUMPSYSTEM_MYNARD0DSYSTEM_HPP_ */
