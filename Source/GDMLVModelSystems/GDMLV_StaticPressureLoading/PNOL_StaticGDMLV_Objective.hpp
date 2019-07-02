/*
 * PNOL_StaticGDMLV_Objective.hpp
 *
 *  Created on: Nov 28, 2018
 *      Author: brian
 */

#ifndef GDMLVMODELSYSTEMS_GDMLV_STATICPRESSURELOADING_PNOL_STATICGDMLV_OBJECTIVE_HPP_
#define GDMLVMODELSYSTEMS_GDMLV_STATICPRESSURELOADING_PNOL_STATICGDMLV_OBJECTIVE_HPP_




#include <cmath>
#include "NonlinearSolvers/Nonl_Objective.hpp"

// Local includes
#include "../../GDMLV/gdmlv.hpp"
#include "../../ValveFunctions/ValveFunctions.hpp"
#include "../../ActivationFunctions/temporalActivationFunctions.hpp"
#include "../../VaryingElastanceModels/varyingElastance.hpp"


class StaticLVResidual : public NonlObjective {

  private:
	int evals;

	double Plv, At, Prv;

	GdmLV * lvPtr;

  public:

	// main function (instantiated from virtual objective)
	void Feval(  vector <double> & X, vector <double> & F   )
	{
		evals++;

		// update model
		lvPtr->computeStaticLV( X.data(), Plv, At, Prv );

		// compute residuals
		lvPtr->computeStaticF( F.data(), Plv );

	}

	StaticLVResidual(){
		evals = 0;
	}

	double getEvals(){	return evals; }

	void setPointers( GdmLV & lv, double PlvIn, double AtIn, double PrvIn )
	{
		lvPtr = &lv;
		Plv = PlvIn;
		At = AtIn;
		Prv = PrvIn;
	}

};


#endif /* GDMLVMODELSYSTEMS_GDMLV_STATICPRESSURELOADING_PNOL_STATICGDMLV_OBJECTIVE_HPP_ */
