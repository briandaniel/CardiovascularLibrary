/*
 * PNOL_CDMLV_StaticObjective.hpp
 *
 *  Created on: Jan 29, 2019
 *      Author: brian
 */

#ifndef CDMLVMODELSYSTEMS_CDMLV_STATICPRESSURELOADING_PNOL_CDMLV_STATICOBJECTIVE_HPP_
#define CDMLVMODELSYSTEMS_CDMLV_STATICPRESSURELOADING_PNOL_CDMLV_STATICOBJECTIVE_HPP_


#include <cmath>
#include "NonlinearSolvers/Nonl_Objective.hpp"

// Local includes
#include "../../CDMLV/cdmlv.hpp"

class StaticCDMLVResidual : public NonlObjective {

  private:
	int evals;

	double Plv, At ;

	CDMLVModel * lvPtr;

  public:

	// main function (instantiated from virtual objective)
	void Feval(  vector <double> & X, vector <double> & F   )
	{
		evals++;

		// update model
		// note: X := q
		lvPtr->computeStaticLV( X, At );

		// compute residuals
		lvPtr->computeStaticF( F, Plv );

	}

	StaticCDMLVResidual(){
		evals = 0;
	}

	double getEvals(){	return evals; }

	void setPointers( CDMLVModel & lv, double PlvIn, double AtIn )
	{
		lvPtr = &lv;
		Plv = PlvIn;
		At = AtIn;
	}

};




#endif /* CDMLVMODELSYSTEMS_CDMLV_STATICPRESSURELOADING_PNOL_CDMLV_STATICOBJECTIVE_HPP_ */
