/*
 * CDMLV_HalfSystem.hpp
 *
 *  Created on: Jan 29, 2019
 *      Author: brian
 */

#ifndef CDMLVMODELSYSTEMS_CDMLV_VRVALVES_HALFSYSTEM_CDMLV_HALFSYSTEM_HPP_
#define CDMLVMODELSYSTEMS_CDMLV_VRVALVES_HALFSYSTEM_CDMLV_HALFSYSTEM_HPP_

// Local includes
#include "lumpActiveModelParam.hpp"
#include "CDMLV_HalfSystem_Output.hpp"
#include "../../CDMLV/cdmlv.hpp"

// nonlinear solver include
#include "NonlinearSolvers/NewtonsMethod.hpp"
#include "NonlinearSolvers/Nonl_Objective.hpp"

// local functions
void computeLVSystemHalf( CDMLVModel & lv, ParamLumpActive & lumpPrm, string initialValuesFileName, string outFileName );
void computeLVSystemHalfVectorOutput( int Naux, vector<double> & X, vector<double> & q, vector <double> & w,
		CDMLVModel & lv, ParamLumpActive & lumpPrm, vector<double> & tStore, vector<vector<double>> & varStore, int & NtComputed );
void eulerStepCDMLVHalf( double dt, double t, vector<double> & q, vector<double> & w, CDMLVModel & lv, ParamLumpActive & lumpPrm,
		vector <double> & X, double & Plv, double & Ppao, double & Pla, double & Vlv, vector <double> & dqdt,  vector <double> & dwdt );
void evaluateModels_HalfSystem( double t, vector<double> & q, vector<double> & w, CDMLVModel & lv, ParamLumpActive & lumpPrm,
		vector <double> & X, double & Plv, double & Ppao, double & Pla, double & Vlv, vector <double> & dqdt,  vector <double> & dwdt );
void lumpTimeDerivatives_CDMLVHalf( double t, double Plv, double Ppao, double Pla, ParamLumpActive & lumpPrm, CDMLVModel & lv,
		vector <double> & w, vector <double> & dwdt );
void computeCDMLVHalf_TimeDerivatives_Pressures( double t, double Vla, double Psa, ParamLumpActive & lumpPrm, CDMLVModel & lv,
		vector <double> & X, vector <double> & dqdt, double & Plv, double & Ppao, double & Pla );
void computeNewtonFunction( vector<double> & X, double t, double Vla, double Psa, CDMLVModel & lv, ParamLumpActive & lumpPrm, vector<double> & F);
double computeTimeStepCDMLVHalf( double t, ParamLumpActive & lumpPrm );








// Nonlinear solver objective
class CDMLVHalfNewtObj : public NonlObjective {

  private:
	int evals;

	double t, Vla, Psa;
	CDMLVModel * lvPtr;
	ParamLumpActive * lumpPrmPtr;

  public:

	// main function (instantiated from virtual objective)
	void Feval(  vector <double> & X, vector <double> & F   )
	{
		evals++;

		// update model
		computeNewtonFunction(X, t, Vla, Psa, *lvPtr, *lumpPrmPtr, F);

	}

	CDMLVHalfNewtObj(){
		evals = 0;
	}

	double getEvals(){	return evals; }

	void setPointers( double tIn, double VlaIn, double PsaIn, CDMLVModel & lv, ParamLumpActive & lumpPrm )
	{
		t = tIn;
		Vla = VlaIn;
		Psa = PsaIn;
		lvPtr = &lv;
		lumpPrmPtr = &lumpPrm;

	}

};


#endif /* CDMLVMODELSYSTEMS_CDMLV_VRVALVES_HALFSYSTEM_CDMLV_HALFSYSTEM_HPP_ */
