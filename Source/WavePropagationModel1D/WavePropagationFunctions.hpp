/*
 * WavePropagationFunctions.hpp
 *
 *  Created on: Jan 4, 2019
 *      Author: brian
 */

#ifndef WAVEPROPAGATIONMODEL1D_WAVEPROPAGATIONFUNCTIONS_HPP_
#define WAVEPROPAGATIONMODEL1D_WAVEPROPAGATIONFUNCTIONS_HPP_

#define PI 3.141592653589793

#include <math.h>
#include <vector>
#include "UtilityFunctions/utilityFunctions.hpp"

using namespace std;

class Vessel1D{

  public:

	// size / constants
	int Nx;
	double dx;
	double rho;
	double xi;
	int NxBndConst; // number of constant points on either side for computing the boundaries

	// variables
	vector<double> x;
	vector<double> A;
	vector<double> u;

	// spatially dependent constants
	vector<double> A0, alpha, mu;
	vector<double> R0;
	vector<double> A0Mid, alphaMid, muMid;


	// Vectors used in computations (do _not_ need initial values)
	vector <double> Amid;
	vector <double> umid;
	vector <double> F1;
	vector <double> F2;
	vector <double> S1;
	vector <double> S2;
	vector <double> F1mid;
	vector <double> F2mid;
	vector <double> S1mid;
	vector <double> S2mid;



	// Create constant 1D vessel
	void createConstantVessel( int NxIn, double xlength, double A0in, double nuIn, double xiIn, double rhoIn, double EIn, double muIn, double h0In  );


	// Set initial conditions
	void setZeroFlowInitialConditions() {
		for(int i = 0; i < A.size(); i++){ A[i] = A0[i];}
		for(int i = 0; i < u.size(); i++){ u[i] = 0;}
	}

	// local functions
	void update1DVesselLW( double dt );
	void computeOutgoingCharacteristics( double dt, double & wmInlet, double & wpOutlet );
	void computeBoundaries( double wpIn, double wmIn, double wpOut, double wmOut, double & Alb, double & ulb, double & Arb, double & urb );
	void updateBoundaries( double Alb, double ulb, double Arb, double urb  );

};

// local functions
void computeFlux( double A, double u, double alpha, double A0, double rho, double &F1, double &F2);
void computeSource( double A, double u, double mu, double rho, double xi, double &S1, double &S2 );

#endif /* WAVEPROPAGATIONMODEL1D_WAVEPROPAGATIONFUNCTIONS_HPP_ */
