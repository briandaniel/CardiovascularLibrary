/*
 * CDMLV_StaticLoading.hpp
 *
 *  Created on: Jan 29, 2019
 *      Author: brian
 */

#ifndef CDMLVMODELSYSTEMS_CDMLV_STATICPRESSURELOADING_CDMLV_STATICLOADING_HPP_
#define CDMLVMODELSYSTEMS_CDMLV_STATICPRESSURELOADING_CDMLV_STATICLOADING_HPP_


// Local includes
#include "PNOL_CDMLV_StaticObjective.hpp"
#include "../../CDMLV/cdmlv.hpp"

// nonlinear solver include
#include "NonlinearSolvers/NewtonsMethod.hpp"

// local functions
// computes the parameters q that satisfy the variational equations under static loading conditions
void computeStaticSolution( CDMLVModel & lv, double Plv, double At, vector <double> & q, string prmPrefix, DataContainer & prmData );

// computes the parameters q that satisfy the variational equations under static loading conditions
// computes incremental loading conditions from Plv = 0 up to Plv = Plv using steps deltaP
void computeStaticSolutionContinuation( double deltaP, CDMLVModel & lv, double Plv, double At, vector <double> & q, string prmPrefix, DataContainer & prmData );

// computes the parameters q that satisfy the variational equations under static loading conditions
// computes incremental loading conditions from Plv = 0 up to Plv = Plv using steps deltaP
// Recomputes with small increments if sltn fails
void computeStaticSolutionContinuationRecomputes( double deltaP0,  CDMLVModel & lv, double Plv, double At, vector <double> & q, string prmPrefix, DataContainer & prmData );

#endif /* CDMLVMODELSYSTEMS_CDMLV_STATICPRESSURELOADING_CDMLV_STATICLOADING_HPP_ */
