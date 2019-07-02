/*
 * simpleLumpSystem.hpp
 *
 *  Created on: Sep 13, 2018
 *      Author: brian
 */

#ifndef LUMPMODELSYSTEMS_SIMPLELUMPSYSTEMHALF_SIMPLELUMPSYSTEMHALF_HPP_
#define LUMPMODELSYSTEMS_SIMPLELUMPSYSTEMHALF_SIMPLELUMPSYSTEMHALF_HPP_

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
#include "../../VaryingElastanceModels/varyingElastance.hpp"
#include "../../ValveFunctions/ValveFunctions.hpp"

#include "UtilityFunctions/utilityFunctions.hpp"
#include "ImportExport/importExport.hpp"
#include "ParameterParser/ParameterParser.hpp"


// local param class definition
class ParamSimpleLumpHalf {

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
	double Emax_lv, Emin_lv, Ks_lv, V0_lv;
	double Emax_la, Emin_la, Ks_la, V0_la;

	// Activation model
	double m1_ventricle, m2_ventricle, tau1_ventricle, tau2_ventricle, Ts_ventricle, hillMaxVal_ventricle;
	double m1_atrium, m2_atrium, tau1_atrium, tau2_atrium, Ts_atrium, hillMaxVal_atrium;

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
	ParamSimpleLumpHalf() // sets default values
	{
		setDefaultParams();
		// the number of variables and output aux variables is fixed in this code:
		setNvars();
		paramPreloopComputations();
	}; // sets default values
	ParamSimpleLumpHalf(DataContainer & prmData) // sets values according to paramFileName
	{
		importParams( prmData);
		// the number of variables and output aux variables is fixed in this code:
		setNvars();
		paramPreloopComputations();
	};
	~ParamSimpleLumpHalf(){}; // default destructor (does nothing)

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



// local functions
void computeSimpleLumpSystemHalfNewtonObj( vector<double> &X, vector<double> &F,
		double t, vector <double> & w, ParamSimpleLumpHalf & lumpPrm );


// Nonlinear solver objective
class LumpResidual : public NonlObjective {

  private:
	int evals;

	double t;
	vector <double> * wPtr;
	ParamSimpleLumpHalf * lumpPrmPtr;

  public:

	// main function (instantiated from virtual objective)
	void Feval(  vector <double> & X, vector <double> & F   )
	{
		evals++;

		// update model
		computeSimpleLumpSystemHalfNewtonObj( X, F, t, *wPtr, *lumpPrmPtr );

	}

	LumpResidual(){
		evals = 0;
	}

	double getEvals(){	return evals; }

	void setPointers( double tIn, vector <double> & wIn, ParamSimpleLumpHalf & lumpPrm )
	{
		t = tIn;
		wPtr = &wIn;
		lumpPrmPtr = &lumpPrm;

	}

};


// simpleLumpSystem.cpp
void computeSimpleLumpSystemHalf( vector <double> & w0, vector <double> & X0, ParamSimpleLumpHalf & lumpPrm );

void computeSimpleLumpSystemHalfValues( vector <double> & w, vector <double> & X, ParamSimpleLumpHalf & lumpPrm,
		vector<vector<double>> & wStore, vector<vector<double>> & auxStore, vector <double> & tStore );

void computeSimpleLumpSystemHalfTD( double t, vector <double> & w, vector<double> & X,
		ParamSimpleLumpHalf & lumpPrm, vector <double> & dwdt, vector<double> & auxVars );

void computeSimpleLumpSystemPressures( double t, vector <double> & w, vector<double> & X,
		ParamSimpleLumpHalf & lumpPrm, double & Plv, double & Pla, double & Ppao );

void computeSimpleLumpSystemHalfNewtonObj( vector<double> &X, vector<double> &F,
		double t, vector <double> & w, ParamSimpleLumpHalf & lumpPrm );

void storeValuesSimpleLumpHalf( vector<double> &w,  vector<double> &auxVars, double t, int stepIdx,
		ParamSimpleLumpHalf & lumpPrm , vector<vector<double>> & wStore, vector<vector<double>> & auxStore, vector<double>& tStore);

// simpleLumpSystemHalf_parameters.cpp
void setDefaultInitialValuesSimpleLumpHalf( double * w );

// simpleLumpSystem_output.cpp
void outputModelResultSimpleLumpHalf(int Nsteps, vector<vector<double>> & wStore, vector<vector<double>> & auxStore, vector<double>& tStore, ParamSimpleLumpHalf & lumpPrm );




#endif /* LUMPMODELSYSTEMS_SIMPLELUMPSYSTEM_SIMPLELUMPSYSTEMHALF_HPP_ */
