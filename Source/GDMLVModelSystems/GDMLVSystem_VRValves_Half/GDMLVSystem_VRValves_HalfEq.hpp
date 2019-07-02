/*
 * GDMLVSystem_VRValves_HalfEq.hpp
 *
 *  Created on: Dec 21, 2018
 *      Author: brian
 */

#ifndef GDMLVMODELSYSTEMS_GDMLVSYSTEM_VRVALVES_HALF_GDMLVSYSTEM_VRVALVES_HALFEQ_HPP_
#define GDMLVMODELSYSTEMS_GDMLVSYSTEM_VRVALVES_HALF_GDMLVSYSTEM_VRVALVES_HALFEQ_HPP_


#include "GDMLVSystem_VRValves_Half.hpp"

// local functions
void computeEquilibriumCycleGdmLVAltRecompute( double volEqRequirement, double minEquilibriumCycles, double maxEquilibriumCycles, int maxRecomputes,
		vector <double> & q, vector <double> & w, vector <double> & X,	vector<double> & tStore, vector<double> & PlvStore,
		vector <vector <double> > & qStore,	GdmLV & lv, ParamGdmlvAlt & altPrm, DataContainer & prmData );

int computeEquilibriumCycleGdmLVAlt( double volEqRequirement, double minEquilibriumCycles, double maxEquilibriumCycles, vector <double> & q,
		vector <double> & w, vector <double> & X, vector<double> & tStore, vector<double> & PlvStore, vector <vector <double> > & qStore,
		GdmLV & lv, ParamGdmlvAlt & altPrm, DataContainer & prmData  );
#endif /* GDMLVMODELSYSTEMS_GDMLVSYSTEM_VRVALVES_HALF_GDMLVSYSTEM_VRVALVES_HALFEQ_HPP_ */
