/*
 * SingleFiberMynardValveSystem.hpp
 *
 *  Created on: Dec 12, 2018
 *      Author: brian
 */

#ifndef LUMPMODELSYSTEMS_SINGLEFIBERMYNARDVALVESYSTEM_SINGLEFIBERMYNARDVALVESYSTEM_HPP_
#define LUMPMODELSYSTEMS_SINGLEFIBERMYNARDVALVESYSTEM_SINGLEFIBERMYNARDVALVESYSTEM_HPP_



// Parallel commands
#define ROOT_ID 0	 // id of root processor

// other includes
#include <mpi.h>
#include <cmath>

// nonlinear solver include
#include "NonlinearSolvers/NewtonsMethod.hpp"
#include "NonlinearSolvers/Nonl_Objective.hpp"

// local includes
#include "../../ActivationFunctions/temporalActivationFunctions.hpp"
#include "../../LumpModels/LumpModel.hpp"
#include "../../SingleFiberModel/SingleFiberModel.hpp"
#include "../../ValveFunctions/ValveFunctions.hpp"

#include "UtilityFunctions/utilityFunctions.hpp"
#include "ImportExport/importExport.hpp"
#include "ParameterParser/ParameterParser.hpp"


// local param class definition
class ParamSingleFiberMynardValve {

  public:

	//-------- 1. variables --------//

	// Fixed pressures
	double Ppv_fixed, Psp_fixed;

	// Flow resistance of pulmonary veins
	double Rpv;

	// lump models
	double Csa; // capacitance
	double Rsa;

	// varying elastance models
	double V0_lv, Vw_lv, Ta0_lv, Tp0_lv, cp_lv, cv_lv, v0_lv, ls0_lv, lsmax_lv, lsw_lv;
	double V0_la, Vw_la, Ta0_la, Tp0_la, cp_la, cv_la, v0_la, ls0_la, lsmax_la, lsw_la;

	// Activation model
	double m1_ventricle, m2_ventricle, tau1_ventricle, tau2_ventricle, Ts_ventricle, hillMaxVal_ventricle;
	double m1_atrium, m2_atrium, tau1_atrium, tau2_atrium, Ts_atrium, hillMaxVal_atrium;

	// Valve models
	double rho; // same for all models
	double Aeff_max_mv, Aeff_min_mv, leff_mv, Kvo_mv, Kvc_mv, deltaP_open_mv, deltaP_close_mv;
	double Aeff_max_aov, Aeff_min_aov, leff_aov, Kvo_aov, Kvc_aov, deltaP_open_aov, deltaP_close_aov;


	// Timing and computation
	double Tc, Ta;
	int Nvars, NauxVars;
	int Ncyc, NtStepsCycle;
	double dt;

	// newton iteration parameters
	double xMinDiff, fNewtMin, dX, maxNewtIter;

	// Output file names
	string outputFile;
	string paramOutputFile;
	string initialValuesFileName;

	//-------- 2. functions --------//
	// simpleLumpSystemHalf.cpp
	ParamSingleFiberMynardValve() // sets default values
	{
		setDefaultParams();
		// the number of variables and output aux variables is fixed in this code:
		setNvars();
		paramPreloopComputations();
	}; // sets default values
	ParamSingleFiberMynardValve(DataContainer & prmData) // sets values according to paramFileName
	{
		importParams( prmData);
		// the number of variables and output aux variables is fixed in this code:
		setNvars();
		paramPreloopComputations();
	};
	~ParamSingleFiberMynardValve(){}; // default destructor (does nothing)

	void setNvars(){
		Nvars = 9;
		NauxVars = 0;
	};

	// simpleLumpSystemHalf_parameters.cpp
	void setDefaultParams();
	void importParams( DataContainer & prmData );
	void importInitialValues( vector <double> & w,  DataContainer & prmData );
	void paramPreloopComputations();

	// simpleLumpSystemHalf_output.cpp
	void printParams(); // prints the params in matlab format
	void printParams( string outputFileName );

};



void computeSingleFiberMynardValveNewtonObj( vector<double> & wnp1, vector <double> & wn, vector <double > & Fn,
		double tn, double dt, ParamSingleFiberMynardValve lumpPrm, vector <double> & G );
void computeSingleFiberMynardValveFFunc( vector<double> & w, double t, ParamSingleFiberMynardValve lumpPrm, vector <double> & F,
		 double & Gplv, double & Gpla);

// Nonlinear solver objective
class LumpSingleFiberMynardResidual : public NonlObjective {

  private:
	int evals;

	vector <double> * wnPtr;
	vector <double> * FnPtr;
	double tn;
	double dt;
	ParamSingleFiberMynardValve * lumpPrmPtr;

  public:

	// main function (instantiated from virtual objective)
	void Feval(  vector <double> & X, vector <double> & F   )
	{
		evals++;

		// update model
		computeSingleFiberMynardValveNewtonObj( X, *wnPtr, *FnPtr, tn, dt, *lumpPrmPtr, F);

	}

	LumpSingleFiberMynardResidual(){
		evals = 0;
	}

	double getEvals(){	return evals; }

	void setPointers( double tIn, double dtIn, vector <double> & wn, vector <double> & Fn, ParamSingleFiberMynardValve & lumpPrm )
	{
		tn = tIn;
		dt = dtIn;
		wnPtr = &wn;
		FnPtr = &Fn;
		lumpPrmPtr = &lumpPrm;

	}

};

void computeSingleFiberMynardValve( vector <double> & w, ParamSingleFiberMynardValve & lumpPrm );

void computeSingleFiberMynardValveValues( vector <double> & w, ParamSingleFiberMynardValve & lumpPrm,
		vector<vector<double>> & wStore, vector<vector<double>> & auxStore, vector <double> & tStore );

void computeSingleFiberMynardValveTimeStep( double t, double dt, vector <double> & w,
		ParamSingleFiberMynardValve & lumpPrm );

void storeValuesSingleFiberMynardValve( vector<double> &w,  vector<double> &auxVars, double t, int stepIdx,
	ParamSingleFiberMynardValve & lumpPrm , vector<vector<double>> & wStore, vector<vector<double>> & auxStore, vector<double>& tStore);

void outputModelResultSingleFiberMynardValve( int Nsteps, vector<vector<double>> & wStore,
		vector<vector<double>> & auxStore, vector<double>& tStore, ParamSingleFiberMynardValve & lumpPrm );


#endif /* LUMPMODELSYSTEMS_SINGLEFIBERMYNARDVALVESYSTEM_SINGLEFIBERMYNARDVALVESYSTEM_HPP_ */
