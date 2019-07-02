/*
 * cdmlv_lid.hpp
 *
 *  Created on: Jan 25, 2019
 *      Author: brian
 */

#ifndef CDMLV_CDMLV_LID_HPP_
#define CDMLV_CDMLV_LID_HPP_



#include <iostream>
#include <string>
#include <vector>

// local headers
#include "cdmlv_param.hpp"

// Utility library headers
#include "ProlateSplines/prolateSplines.hpp"
#include "CubicSplines1D/cubicSplinesFunctions1D.hpp"
#include "../LVFunctions/LVFourierDeformation.hpp"
#include "../LVFunctions/LVInitial.hpp"
#include "ImportExport/importExport.hpp"



class LidCDMLV{

  public:

	double a;
	double nuUp0Min;
	double chmue;
	int Nq;
	int Nr, NmuLid, Nphi;

	// other doubles
	double dphi0;
	double Virreg1;
	double acube;
	double zMean;
	double Virreg2;

	// vectors that give the values at the endocardial boundary at the base
	vector<double> mu0Vec, nu0Vec, phi0Vec;
	vector<double> muVec, nuVec, phiVec;
	vector<double> nuUp0Vec;
	vector<double> dmuin0dnu0Vec, dmuin0dphi0Vec;

	vector<double> fc0Vec, dfc0dnu0Vec, dfc0dphi0Vec;
	vector<double> dmudphi0Vec, dnudphi0Vec, dphidphi0Vec;

	// second surface vectors at the boundary
	vector<double> rVec, zVec;


	// 2D arrays for surface 1
	vector<vector<double>> dS0;
	vector<vector<double>> mu0Surf1, nu0Surf1, phi0;
	vector<vector<double>> muSurf1, nuSurf1, phi;

	// 2D arrays for surface 2
	vector<vector<double>> zSurf1;
	vector<vector<double>> r, zLid, ur;
	vector<vector<double>> x0, y0, z0;
	vector<vector<double>> x, y, z;

	// lid displacement vector
	vector< vector<vector<double>> > uLid;

	// lid displacement vector derivatives
	vector< vector< vector<vector<double>> > > duLiddq;

	// lid normal vector
	vector< vector<vector<double>> >  vCross;

	// surface tangent vectors (required for integral)
	vector< vector<vector<double>> >  vur;
	vector< vector<vector<double>> > vphi0;

	// integration weight
	vector<vector<double>> cylWeight;

	// volume derivatives
	vector<double> dVirreg1dq;
	vector<double> dVirreg2dq;

	// functions

	// default constructor requires the size of the lid as an input (usually taken from prm)
	LidCDMLV(){
		Nphi = -1;
		NmuLid = -1;
	};


	// constructor
	void setupLid( ParamCDMLV & prm );

	// initial functions
	void setInitialPositions( ParamCDMLV & prm);

	// regular functions
	void exportVolumes( double & Virreg, vector<double> & dVirregdq );
	void computeLidSurfaceIntegrals( vector<double> & etal );
	void computeLid( vector<double> & q, ParamCDMLV & prm );
	void computeSurfaceNormals();
	void computeFirstSurfacePosition( vector<double> & q, ParamCDMLV & prm );
	void computeFirstIrregularVolume();
	void computeSecondSurfacePosition();
	void computeSecondIrregularVolume();
	void computeDisplacement();


};





#endif /* CDMLV_CDMLV_LID_HPP_ */
