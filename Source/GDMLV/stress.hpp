/*
 * stress.h
 *
 *  Created on: Sep 14, 2017
 *      Author: brian
 */

#ifndef STRESS_H_
#define STRESS_H_

#include "gdmlv_node.hpp"

// local functions
double tensionFunc( ParamGdmLV * prm, double lambda);
double activationFunc( double t, double eff_ed, ParamGdmLV * prm );

#endif /* STRESS_H_ */
