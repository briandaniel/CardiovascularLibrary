/*
 * SingleFiberVarResSystem.hpp
 *
 *  Created on: Dec 11, 2018
 *      Author: brian
 */

#ifndef LUMPMODELSYSTEMS_SINGLEFIBERVARRESSYSTEM_SINGLEFIBERVARRESSYSTEM_HPP_
#define LUMPMODELSYSTEMS_SINGLEFIBERVARRESSYSTEM_SINGLEFIBERVARRESSYSTEM_HPP_


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
class ParamSingleFiberVarRes {

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
	double activeDurationAdjustment;

	// valve
	double Rmvo, Rmvc, Raovo, Raovc, beta;
	double Rpao;

	// Timing and computation
	double Tc, Ta;
	int Nvars, Nnewt, NauxVars;
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
	ParamSingleFiberVarRes() // sets default values
	{
		setDefaultParams();
		// the number of variables and output aux variables is fixed in this code:
		setNvars();
		paramPreloopComputations();
	}; // sets default values
	ParamSingleFiberVarRes(DataContainer & prmData) // sets values according to paramFileName
	{
		importParams( prmData);
		// the number of variables and output aux variables is fixed in this code:
		setNvars();
		paramPreloopComputations();
	};
	~ParamSingleFiberVarRes(){}; // default destructor (does nothing)

	void setNvars(){
		Nvars = 3;
		Nnewt = 3;
		NauxVars = 8;
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



void computeSingleFiberVarResNewtonObj( vector<double> &X, vector<double> &F,
		double t, vector <double> & w, ParamSingleFiberVarRes & lumpPrm );


// Nonlinear solver objective
class LumpSingleFiberResidual : public NonlObjective {

  private:
	int evals;

	double t;
	vector <double> * wPtr;
	 ParamSingleFiberVarRes * lumpPrmPtr;

  public:

	// main function (instantiated from virtual objective)
	void Feval(  vector <double> & X, vector <double> & F   )
	{
		evals++;

		// update model
		computeSingleFiberVarResNewtonObj( X, F, t, *wPtr, *lumpPrmPtr );

	}

	LumpSingleFiberResidual(){
		evals = 0;
	}

	double getEvals(){	return evals; }

	void setPointers( double tIn, vector <double> & wIn,  ParamSingleFiberVarRes & lumpPrm )
	{
		t = tIn;
		wPtr = &wIn;
		lumpPrmPtr = &lumpPrm;

	}

};



// other local functions
void computeSingleFiberVarRes( vector <double> & w0, vector <double> & X0, ParamSingleFiberVarRes & lumpPrm );
void computeSingleFiberVarResValues( vector <double> & w, vector <double> & X, ParamSingleFiberVarRes & lumpPrm,
		vector<vector<double>> & wStore, vector<vector<double>> & auxStore, vector <double> & tStore );
void computeSingleFiberVarResTD( double t, vector <double> & w, vector<double> & X,
		ParamSingleFiberVarRes & lumpPrm, vector <double> & dwdt, vector<double> & auxVars );
void computeSingleFiberVarResPressures( double t, vector <double> & w, vector<double> & X,
		ParamSingleFiberVarRes & lumpPrm, double & Plv, double & Pla, double & Ppao );
void computeSingleFiberVarResNewtonObj( vector<double> &X, vector<double> &F,
		double t, vector <double> & w, ParamSingleFiberVarRes & lumpPrm );
void storeValuesSingleFiberVarRes( vector<double> &w,  vector<double> &auxVars, double t, int stepIdx,
	ParamSingleFiberVarRes & lumpPrm , vector<vector<double>> & wStore, vector<vector<double>> & auxStore, vector<double>& tStore);


void outputModelResultSingleFiberVarRes(int Nsteps, vector<vector<double>> & wStore,
		vector<vector<double>> & auxStore, vector<double>& tStore, ParamSingleFiberVarRes & lumpPrm );

#endif /* LUMPMODELSYSTEMS_SINGLEFIBERVARRESSYSTEM_SINGLEFIBERVARRESSYSTEM_HPP_ */
