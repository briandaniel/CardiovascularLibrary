/*
 * CDMLV_HalfSystem_Output.hpp
 *
 *  Created on: Jan 29, 2019
 *      Author: brian
 */

#ifndef CDMLVMODELSYSTEMS_CDMLV_VRVALVES_HALFSYSTEM_CDMLV_HALFSYSTEM_OUTPUT_HPP_
#define CDMLVMODELSYSTEMS_CDMLV_VRVALVES_HALFSYSTEM_CDMLV_HALFSYSTEM_OUTPUT_HPP_

// Local includes
#include "lumpActiveModelParam.hpp"
#include "../../CDMLV/cdmlv.hpp"

// local functions
void saveDataValuesCDMLVHalf( int stepIdx, double t, vector<double>&  q, vector<double> & w, vector<double>&  aux,
		vector<double> & tStore, vector<vector<double>> & varStore, ParamLumpActive & lumpPrm, double Plv, double Vlv, double Pla, double Ppao  );
void exportDataValuesCDMLVHalf( string outFileName, int Ncomputed, int Nq, int Nw, int Naux, vector<double> & tStore, vector<vector<double>> & varStore );
void exportInitialValuesCDMLVHalf( string initialValuesFileName, vector<double> &X, vector<double> & q, vector<double> &w  );
void readInitialValuesCDMLVHalf( string initialValuesFileName, vector<double> &X, vector<double> & q, vector<double> &w  );


#endif /* CDMLVMODELSYSTEMS_CDMLV_VRVALVES_HALFSYSTEM_CDMLV_HALFSYSTEM_OUTPUT_HPP_ */
