/*
 * gdmlv_output.hpp
 *
 *  Created on: Nov 27, 2018
 *      Author: brian
 */

#ifndef GDMLV_GDMLV_OUTPUT_HPP_
#define GDMLV_GDMLV_OUTPUT_HPP_


#include "gdmlv.hpp"

#include "RigidMotion/rigidMotion.hpp"

// local definitions
void printStaticLVShapeMatlabFormat( double * q, string outputName, DataContainer & prmData);
void printStaticLVShapeMatlabFormatFixedRot( vector <double> & q, vector <double> & theta, vector <double> & shift, string outputName, DataContainer & prmData);
void printInitialLVShapeMatlabFormat( string outputName, DataContainer & prmData);
void printModelStateMatlabFormat( double * q, int Nmu, int Nnu, int Nphi, int NmuLid, int NphiLid, GdmLV & lv, RigidMotion & rm, string outputName );
void printModelValuesMatlabFormat( double * q, double At, GdmLV & lv, RigidMotion & rm, string outputName );
void printLVShapesTimeSeries( vector<vector<double>> & qStore, int Nsteps, int Nskip, string outputFolder, DataContainer & prmData);

#endif /* GDMLV_GDMLV_OUTPUT_HPP_ */


