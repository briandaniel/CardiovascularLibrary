/*
 * simpleLumpSystem.hpp
 *
 *  Created on: Sep 13, 2018
 *      Author: brian
 */

#ifndef LUMPMODELSYSTEMS_SIMPLELUMPSYSTEM_SIMPLELUMPSYSTEM_HPP_
#define LUMPMODELSYSTEMS_SIMPLELUMPSYSTEM_SIMPLELUMPSYSTEM_HPP_

// Parallel commands
#define ROOT_ID 0	 // id of root processor

// other includes
#include <mpi.h>

// local includes
#include "../../ActivationFunctions/temporalActivationFunctions.hpp"
#include "../../LumpModels/LumpModel.hpp"
#include "../../VaryingElastanceModels/varyingElastance.hpp"

#include "UtilityFunctions/utilityFunctions.hpp"
#include "ImportExport/importExport.hpp"
#include "ParameterParser/ParameterParser.hpp"

// local param class definition
class ParamSimpleLump {

  public:

	//-------- 1. variables --------//
	// lump models
	double Cla, Cpv, Cpp, Cpa, Cra, Csv, Csp, Csa; // capacitance
	double Vu_la, Vu_pv, Vu_pp, Vu_pa, Vu_ra, Vu_sv, Vu_sp, Vu_sa; // unstressed volumes
	double Rla, Rpv, Rpp, Rpa, Rra, Rsv, Rsp, Rsa;

	// varying elastance models
	double P0_rv, ke_rv, Emax_rv, kr_rv, Vu_rv;
	double P0_lv, ke_lv, Emax_lv, kr_lv, Vu_lv;

	// total blood volume
	double Vtotal;

	// Timing and computation
	double Tc, Ta;
	int Ncyc, Nvars, NauxVars;
	double dt;

	// Output file names
	string outputFile;
	string paramOutputFile;
	string initialValuesFileName;

	//-------- 2. functions --------//
	// simpleLumpSystem.cpp
	ParamSimpleLump() // sets default values
	{
		setDefaultParams();
		// the number of variables and output aux variables is fixed in this code:
		setNvars();
	}; // sets default values
	ParamSimpleLump(DataContainer & prmData) // sets values according to paramFileName
	{
		importParams( prmData);
		// the number of variables and output aux variables is fixed in this code:
		setNvars();
	};
	~ParamSimpleLump(){}; // default destructor (does nothing)
	void setNvars(){
		Nvars = 9;
		NauxVars = 11;
	};

	// simpleLumpSystem_parameters.cpp
	void setDefaultParams();
	void importParams( DataContainer & prmData );
	void importInitialValues( double * w,  DataContainer & prmData );

	// simpleLumpSystem_output.cpp
	void printParams(); // prints the params in matlab format

};

// simpleLumpSystem.cpp
void computeSimpleLumpSystem( ParamSimpleLump & lumpPrm );
void computeSimpleLumpSystemTD( double t, double * w, ParamSimpleLump & lumpPrm, double * dwdt, double * auxVars );
void storeValuesSimpleLump( int stepIdx, double t, double * w, double * auxVars, double *	tStore, double * PlaStore, double * PpvStore,
		double * PppStore, double * PpaStore, double * PraStore, double * PspStore, double * PsaStore, double * VlvStore,
		double * VrvStore, double * PlvStore, double * PrvStore, double * Pmax_lvStore, double * Pmax_rvStore, double * PsvStore,
		double * FilvStore, double * FolvStore, double * FirvStore, double * ForvStore, double * At_lvStore, double * At_rvStore );

// simpleLumpSystem_parameters.cpp
void setDefaultInitialValuesSimpleLump( double * w );

// simpleLumpSystem_output.cpp
void outputModelResultSimpleLump( int Nsteps, ParamSimpleLump & lumpPrm, double * tStore, double * PlaStore, double * PpvStore,
		double * PppStore, double * PpaStore, double * PraStore, double * PspStore, double * PsaStore, double * VlvStore,
		double * VrvStore, double * PlvStore, double * PrvStore, double * Pmax_lvStore, double * Pmax_rvStore, double * PsvStore,
		double * FilvStore, double * FolvStore, double * FirvStore, double * ForvStore, double * At_lvStore, double * At_rvStore );

// local functions









#endif /* LUMPMODELSYSTEMS_SIMPLELUMPSYSTEM_SIMPLELUMPSYSTEM_HPP_ */
