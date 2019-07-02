/*
 * gdmlv_lid.hpp
 *
 *  Created on: Sep 11, 2018
 *      Author: brian
 */

#ifndef GDMLV_GDMLV_LID_HPP_
#define GDMLV_GDMLV_LID_HPP_


#include <iostream>
#include <string>
#include <vector>

// local headers
#include "gdmlv_params.hpp"

// Utility library headers
#include "ProlateSplines/prolateSplines.hpp"
#include "CubicSplines1D/cubicSplinesFunctions1D.hpp"
#include "../LVFunctions/LVBicubicDeformation.hpp"
#include "../LVFunctions/LVFourierDeformation.hpp"
#include "../LVFunctions/LVInitial.hpp"
#include "../LVFunctions/LVDeformation.hpp"
#include "ImportExport/importExport.hpp"



class LidGdmLV{

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
	double *mu0Vec, *nu0Vec, *phi0Vec;
	double *muVec, * nuVec, *phiVec;
	double *nuUp0Vec;
	double *dmuin0dnu0Vec , *dmuin0dphi0Vec ;


	double * fc0Vec, *dfc0dnu0Vec, *dfc0dphi0Vec;
	double * dmudphi0Vec, * dnudphi0Vec ,* dphidphi0Vec;

	// second surface vectors at the boundary
	double * rVec, * zVec;


	// 2D arrays for surface 1
	double ** dS0;
	double ** mu0Surf1, ** nu0Surf1, ** phi0;
	double ** muSurf1, **nuSurf1, **phi;


	// 2D arrays for surface 2
	double ** zSurf1;
	double ** r, **zLid, **ur;

	double ** x0, ** y0, **z0;
	double ** x, ** y, **z ;



	// lid displacement vector
	double *** uLid;

	// lid displacement vector derivatives
	double **** duLiddq;

	// lid normal vector
	double *** vCross;

	// surface tangent vectors (required for integral)
	double *** vur;
	double *** vphi0;

	// integration weight
	double ** cylWeight;


	// volume derivatives
	double * dVirreg1dq;
	double * dVirreg2dq;

	// functions

	// default constructor requires the size of the lid as an input (usually taken from prm)
	LidGdmLV(){
		Nphi = -1;
		NmuLid = -1;
	};


	// constructor
	void setLidSize( int NmuLidIn,  int NphiIn );
	void createLid( ParamGdmLV & prm );
	void createLid( int NphiIn, int NmuLidIn, ParamGdmLV & prm );

	// clean up
	void destroyLid();

	// initial functions
	void setInitialPositions( ParamGdmLV & prm);

	// regular functions
	void exportVolumes( double & Virreg, double * dVirregdq );
	void computeLidSurfaceIntegrals( double * etal );
	void computeLid( double * q, FourierDeformation & fourierDef, ParamGdmLV & prm );
	void computeSurfaceNormals();
	void computeFirstSurfacePosition( double * q, FourierDeformation & fourierDef, ParamGdmLV & prm );
	void computeFirstIrregularVolume();
	void computeSecondSurfacePosition();
	void computeSecondIrregularVolume();
	void computeDisplacement();


};



#endif /* GDMLV_GDMLV_LID_HPP_ */
