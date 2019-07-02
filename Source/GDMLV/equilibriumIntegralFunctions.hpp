/*
 * equilibriumIntegralFunctions.hpp
 *
 *  Created on: Sep 11, 2018
 *      Author: brian
 */

#ifndef GDMLV_EQUILIBRIUMINTEGRALFUNCTIONS_HPP_
#define GDMLV_EQUILIBRIUMINTEGRALFUNCTIONS_HPP_



// local includes
#include "gdmlv_node.hpp"

void computeProcessorIntegrals( vector <NodeGdmLV> & localNodes, int Nnodes, double ** alpha,
				double *kappa, double * etae, double *dVregdq, double & Vreg, ParamGdmLV * prm, double Prv );
void computeStaticProcessorIntegrals( vector <NodeGdmLV> & localNodes, int Nnodes,
				double *kappa, double * etae, double & Vreg, ParamGdmLV * prm, double Prv );
double computeReferenceWallVolume( vector <NodeGdmLV> & localNodes, int Nnodes );

#endif /* GDMLV_EQUILIBRIUMINTEGRALFUNCTIONS_HPP_ */
