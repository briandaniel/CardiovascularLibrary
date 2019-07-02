/*
 * GDMLVSystem_VRValves_Half.hpp
 *
 *  Created on: Dec 18, 2018
 *      Author: brian
 */

#ifndef GDMLVMODELSYSTEMS_GDMLVSYSTEMVARRESHALF_GDMLVSYSTEM_VRVALVES_HALF_HPP_
#define GDMLVMODELSYSTEMS_GDMLVSYSTEMVARRESHALF_GDMLVSYSTEM_VRVALVES_HALF_HPP_


// Local includes
#include "../../GDMLV/gdmlv.hpp"
#include "../../ValveFunctions/ValveFunctions.hpp"
#include "../../ActivationFunctions/temporalActivationFunctions.hpp"
#include "../../VaryingElastanceModels/varyingElastance.hpp"

#include "GDMLVSystem_VRValves_AltParams.hpp"
#include "GDMLVSystem_VRValves_Output.hpp"

// nonlinear solver include
#include "NonlinearSolvers/NewtonsMethod.hpp"
#include "NonlinearSolvers/Nonl_Objective.hpp"


// local functions
void computeNewtonFunction( vector<double> & X, double t, double Vla, double Psa, GdmLV & lv, ParamGdmlvAlt & altPrm, vector<double> & F);



// Nonlinear solver objective
class GdmLVHalfNewtObj : public NonlObjective {

  private:
	int evals;

	double t, Vla, Psa;
	GdmLV * lvPtr;
	ParamGdmlvAlt * altPrmPtr;

  public:

	// main function (instantiated from virtual objective)
	void Feval(  vector <double> & X, vector <double> & F   )
	{
		evals++;

		// update model
		computeNewtonFunction(X, t, Vla, Psa, *lvPtr, *altPrmPtr, F);

	}

	GdmLVHalfNewtObj(){
		evals = 0;
	}

	double getEvals(){	return evals; }

	void setPointers( double tIn, double VlaIn, double PsaIn, GdmLV & lv, ParamGdmlvAlt & altPrm )
	{
		t = tIn;
		Vla = VlaIn;
		Psa = PsaIn;
		lvPtr = &lv;
		altPrmPtr = &altPrm;

	}

};


// local functions
void computeLVSystemHalf( GdmLV & lv, ParamGdmlvAlt & altPrm, string initialValuesFileName, string outFileName );

void computeLVSystemHalfVectorOutput( int Naux, vector<double> & X, vector<double> & q, vector <double> & w,
		GdmLV & lv, ParamGdmlvAlt & altPrm, vector<double> & tStore, vector<vector<double>> & varStore, int & NtComputed );

void eulerStepGDMLVHalf( double dt, double t, vector<double> & q, vector<double> & w, GdmLV & lv, ParamGdmlvAlt & altPrm,
		vector <double> & X, double & Plv, double & Ppao, double & Pla, double & Vlv, vector <double> & dqdt, vector <double> & dwdt );

void evaluateModels_HalfSystem( double t, vector<double> & q, vector<double> & w, GdmLV & lv, ParamGdmlvAlt & altPrm,
		vector <double> & X, double & Plv, double & Ppao, double & Pla, double & Vlv, vector <double> & dqdt,  vector <double> & dwdt );

void computeGDMLVHalf_TimeDerivatives_Pressures( double t, double Vla, double Psa, ParamGdmlvAlt & altPrm, GdmLV & lv,
		vector <double> & X, vector <double> & dqdt, double & Plv, double & Ppao, double & Pla );

void lumpTimeDerivatives_GDMLVHalf( double t, double Plv, double Ppao, double Pla, ParamGdmlvAlt & altPrm, GdmLV & lv,
		vector <double> & w, vector <double> & dwdt );

double computeTimeStepGDMLVHalf( double t, ParamGdmLV & prm );







#endif /* GDMLVMODELSYSTEMS_GDMLVSYSTEMVARRESHALF_GDMLVSYSTEM_VRVALVES_HALF_HPP_ */
