/*
 * lumpModelParam.hpp
 *
 *  Created on: Jan 29, 2019
 *      Author: brian
 */

#ifndef CDMLVMODELSYSTEMS_CDMLV_VRVALVES_HALFSYSTEM_LUMPACTIVEMODELPARAM_HPP_
#define CDMLVMODELSYSTEMS_CDMLV_VRVALVES_HALFSYSTEM_LUMPACTIVEMODELPARAM_HPP_


// local includes
#include "../../ActivationFunctions/temporalActivationFunctions.hpp"
#include "../../LumpModels/LumpModel.hpp"
#include "../../SingleFiberModel/SingleFiberModel.hpp"
#include "../../ValveFunctions/ValveFunctions.hpp"

#include "UtilityFunctions/utilityFunctions.hpp"
#include "ImportExport/importExport.hpp"
#include "ParameterParser/ParameterParser.hpp"


// local param class definition
class ParamLumpActive {

  public:

	//-------- 1. variables --------//
	int Nw;

	// Fixed pressures
	double Ppv_fixed, Psp_fixed;

	// Flow resistance of pulmonary veins
	double Rpv;

	// lump models
	double Csa; // capacitance
	double Rsa;

	// varying elastance model
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
	int Ncyc;
	double dt;
	double dtMin, dtMax;

	// newton iteration parameters
	double xMinDiff, fNewtMin, dX, maxNewtIter;

	//-------- 2. functions --------//
	ParamLumpActive() // sets default values
	{
		activeDurationAdjustment = 0;
		Nw = 2;
		setDefaultParams();
		// the number of variables and output aux variables is fixed in this code:
		paramPreloopComputations(Tc);
	}; // sets default values
	ParamLumpActive( string prmPrefix, DataContainer & prmData) // sets values according to paramFileName
	{
		activeDurationAdjustment = 0;
		Nw = 2;
		importParams( prmPrefix, prmData);
		// the number of variables and output aux variables is fixed in this code:
		paramPreloopComputations(Tc);
	};
	~ParamLumpActive(){}; // default destructor (does nothing)

	// simpleLumpSystemHalf_parameters.cpp
	void setDefaultParams();
	void importParams( string prmPrefix, DataContainer & prmData );
	void paramPreloopComputations(double Tc);

};





#endif /* CDMLVMODELSYSTEMS_CDMLV_VRVALVES_HALFSYSTEM_LUMPACTIVEMODELPARAM_HPP_ */
