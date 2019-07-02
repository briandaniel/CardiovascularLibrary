/*
 * GDMLVSystem_VRValves.hpp
 *
 *  Created on: Nov 20, 2018
 *      Author: brian
 */

#ifndef GDMLVMODELSYSTEMS_GDMLVSYSTEMVARRESVALVES_GDMLVSYSTEM_VRVALVES_HPP_
#define GDMLVMODELSYSTEMS_GDMLVSYSTEMVARRESVALVES_GDMLVSYSTEM_VRVALVES_HPP_

// Local includes
#include "../../GDMLV/gdmlv.hpp"
#include "../../ValveFunctions/ValveFunctions.hpp"
#include "../../ActivationFunctions/temporalActivationFunctions.hpp"
#include "../../VaryingElastanceModels/varyingElastance.hpp"


// Local functions

int computeSolution( GdmLV & lv, int procID, int Nprocs );
void rk2Step( double * q, double * w, double t, double dt, GdmLV * lv,
					double &Plv, double &Ppao, double &Vlv, double *dqdt, double *dwdt );
void eulerStep( double * q, double * w, double t, double dt, GdmLV * lv,
					double &Plv, double &Ppao, double &Vlv, double *dqdt, double *dwdt );
void evaluateModels( double * q, double * w, double t, GdmLV * lv, double &Plv, double &Ppao, double &Vlv,  double *dqdt, double *dwdt );
bool computeLVTimeDerivatives( double * q, double t, double Pla, double Pao, double ** alpha, double * eta,	double * kappa,
	double * dVlvdq, double Plv_guess, double Ppao_guess, double *dqdt_guess, ParamGdmLV * prm, double *dqdt, double &Plv, double &Ppao);
void lvNewtonFunction( double * X, double ** alpha, double * eta,
		double * kappa, double * dVlvdq, double Pla, double Pao, ParamGdmLV * prm, double * H );
void lvNewtonJacobian( double * X, double ** alpha, double * eta,
		 	 	 	 	 double * dVlvdq, double Pla, double Pao, ParamGdmLV * prm, double ** J);
void analyticNewtonDerivatives( double Plv, double Ppao, double Pla, double Pao,
		ParamGdmLV * prm, double &dHNq1_dPlv, double &dHNq1_dPpao, double &dHNq2_dPlv, double &dHNq2_dPpao );
void lvNewtonFunction( double * X, double ** alpha, double * eta,
		double * kappa, double * dVlvdq, double Pla, double Pao, ParamGdmLV * prm, double * H );
double computeTimeStep( double t, ParamGdmLV * prm );
void pullW ( double * w, double & Pla, double & Psa, double & Ppv, double & Ppa, double & Ppp, double & Psp, double & Pra, double & Vrv );
void lumpTimeDerivatives( double t, double Plv, double Ppao, double Vlv, double* w, ParamGdmLV * prm,  double* dwdt );
double computePsv (  double Psa, double Psp, double Pra, double Ppa, double Ppp, double Ppv, double Pla, double Vrv, double Vlv, ParamGdmLV * prm );
void printInitialValues ( int procID, int rootID, string initFileName, double * q, double * w,
							double * dwdt, double * dqdt, double Plv, double Ppao, int Nq, int Nw );

#endif /* GDMLVMODELSYSTEMS_GDMLVSYSTEMVARRESVALVES_GDMLVSYSTEM_VRVALVES_HPP_ */













