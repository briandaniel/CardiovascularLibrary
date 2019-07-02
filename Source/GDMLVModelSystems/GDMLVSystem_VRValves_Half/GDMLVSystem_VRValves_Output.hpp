/*
 * GDMLVSystem_VRValves_Output.hpp
 *
 *  Created on: Dec 19, 2018
 *      Author: brian
 */

#ifndef GDMLVMODELSYSTEMS_GDMLVSYSTEM_VRVALVES_HALF_GDMLVSYSTEM_VRVALVES_OUTPUT_HPP_
#define GDMLVMODELSYSTEMS_GDMLVSYSTEM_VRVALVES_HALF_GDMLVSYSTEM_VRVALVES_OUTPUT_HPP_

// Local includes
#include "../../GDMLV/gdmlv.hpp"
#include "../../ValveFunctions/ValveFunctions.hpp"
#include "../../ActivationFunctions/temporalActivationFunctions.hpp"
#include "../../VaryingElastanceModels/varyingElastance.hpp"

#include "GDMLVSystem_VRValves_AltParams.hpp"


// local functions

void saveDataValuesGDMLVHalf( int stepIdx, double t, vector<double>&  q, vector<double> & w, vector<double>&  aux,
		vector<double> & tStore, vector<vector<double>> & varStore, ParamGdmlvAlt & altPrm, double Plv, double Vlv, double Pla, double Ppao  );

void exportDataValuesGDMLVHalf( string outFileName, int Ncomputed, int Nq, int Nw, int Naux, vector<double> & tStore, vector<vector<double>> & varStore );

void exportInitialValuesGDMLVHalf( string initialValuesFileName, vector<double> &X, vector<double> & q, vector<double> &w  );
void readInitialValuesGDMLVHalf( string initialValuesFileName, vector<double> &X, vector<double> & q, vector<double> &w  );


#endif /* GDMLVMODELSYSTEMS_GDMLVSYSTEM_VRVALVES_HALF_GDMLVSYSTEM_VRVALVES_OUTPUT_HPP_ */
