/*
 * LumpModel.hpp
 *
 *  Created on: Sep 13, 2018
 *      Author: brian
 */

#ifndef LUMPMODELS_LUMPMODEL_HPP_
#define LUMPMODELS_LUMPMODEL_HPP_



/* double evaluateLumpModel(...)
 *  Computes the rate of pressure change in a 0D linear lump parameter model
 *
 * The inputs are: P: chamber pressure, Pup: upstream chamber pressure, Pdown: downstream chamber pressure
 * The parameters are: C: chamber capacitance, R: chamber flow resistance, Rup: upstream chamber flow resistance
 * The output is: dPdt: the rate of change of the chamber pressure
 */
double evaluateLumpModel( double P, double Pup, double Pdown,
						  double C, double R, double Rup);

#endif /* LUMPMODELS_LUMPMODEL_HPP_ */
