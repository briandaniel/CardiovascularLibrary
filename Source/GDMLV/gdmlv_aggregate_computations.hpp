/*
 * gdmlv_aggregate_computations.hpp
 *
 *  Created on: Dec 7, 2018
 *      Author: brian
 */

#ifndef GDMLV_GDMLV_AGGREGATE_COMPUTATIONS_HPP_
#define GDMLV_GDMLV_AGGREGATE_COMPUTATIONS_HPP_


// model includes
#include "gdmlv.hpp"

// local functions
void computeAggregatesSeries( vector<vector<double>> & qStore, GdmLV & lv, vector <double> & VlvVec,
		vector <double> & shortAxisRadiusVec, vector <double> & longAxisLengthVec, vector <double> & meanTwistAngleVec );
void computeAggregates( vector<double> & q, GdmLV & lv, double & Vlv, double & shortAxisMeanRadius, double & longAxisLength, double & meanTwistAngle );
double computeLongAxisLength(  vector<double> & q, GdmLV & lv  );
// double computeShortAxisMeanRadius( GdmLV & lv );
double computeShortAxisMeanRadius( vector<double> & q, GdmLV & lv );
double computeMeanTwistAngle( GdmLV & lv  );


#endif /* GDMLV_GDMLV_AGGREGATE_COMPUTATIONS_HPP_ */
