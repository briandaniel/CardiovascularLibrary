/*
 * GDMLV_StaticLoading.hpp
 *
 *  Created on: Nov 28, 2018
 *      Author: brian
 */

#ifndef GDMLVMODELSYSTEMS_GDMLV_STATICPRESSURELOADING_GDMLV_STATICLOADING_HPP_
#define GDMLVMODELSYSTEMS_GDMLV_STATICPRESSURELOADING_GDMLV_STATICLOADING_HPP_

// Local includes
#include "PNOL_StaticGDMLV_Objective.hpp"
#include "../../GDMLV/gdmlv.hpp"
#include "../../ValveFunctions/ValveFunctions.hpp"
#include "../../ActivationFunctions/temporalActivationFunctions.hpp"
#include "../../VaryingElastanceModels/varyingElastance.hpp"

// nonlinear solver include
#include "NonlinearSolvers/NewtonsMethod.hpp"


// local functions
void computeStaticSolution( GdmLV & lv, double Plv, double At, double Prv, vector <double> & q, DataContainer & prmData );
void computeStaticSolutionContinuation( double deltaP, GdmLV & lv, double Plv, double At, double Prv, vector <double> & q, DataContainer & prmData );
void computeStaticSolutionContinuationRecomputes( double deltaP0, GdmLV & lv, double Plv, double At, double Prv, vector <double> & q, DataContainer & prmData );

#endif /* GDMLVMODELSYSTEMS_GDMLV_STATICPRESSURELOADING_GDMLV_STATICLOADING_HPP_ */
