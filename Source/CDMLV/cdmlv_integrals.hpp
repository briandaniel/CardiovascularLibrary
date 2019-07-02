/*
 * cdmlv_integrals.hpp
 *
 *  Created on: Jan 25, 2019
 *      Author: brian
 */

#ifndef CDMLV_CDMLV_INTEGRALS_HPP_
#define CDMLV_CDMLV_INTEGRALS_HPP_

#include "cdmlv_param.hpp"
#include "cdmlv_node.hpp"

// local functions
void computeRegularIntegrals( vector <CDMNode> & localNodes, int Nnodes, vector<vector<double>> & alpha,
				vector<double> & kappa, vector<double> & etae, vector<double> & dVregdq, double & Vreg, ParamCDMLV & prm );
void computeStaticRegularIntegrals( vector <CDMNode> & localNodes, int Nnodes, vector<double> & kappa,
		vector<double> & etae, double & Vreg, ParamCDMLV & prm );
void computeDynamicRegularIntegrals( vector <CDMNode> & localNodes, int Nnodes,
		vector<vector<double>> & alpha, vector<double> & dVregdq, ParamCDMLV & prm );
double computeReferenceWallVolume ( vector <CDMNode> & localNodes, int Nnodes );

#endif /* CDMLV_CDMLV_INTEGRALS_HPP_ */
