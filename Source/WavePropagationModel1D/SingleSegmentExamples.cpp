/*
 * SingleSegmentExamples.cpp
 *
 *  Created on: Jan 7, 2019
 *      Author: brian
 */

#include "SingleSegmentExamples.hpp"



void computeFixedSingleSegment1DVessel( string outputFileName )
{

	if(getProcID()==ROOT_ID); cout << "Computing 1D vessel with constant properties." << endl;
	Vessel1D vessel;

	int Nt = 1000;
	int Nx = 100;

	double tMax,xMax,A0,nu,xi,rho,E,mu,h0;

	tMax = 1;
	xMax = 100;
	A0 = 4.5;
	nu = 0.5;
	xi = 8;
	rho = 1.06;
	E = 4e6;
	mu = 0.04;
	h0 = 0.16;

	double dt = tMax/(Nt-1);

	// Create constant parameter vessel
	vessel.createConstantVessel(Nx,xMax,A0,nu,xi,rho,E,mu,h0);

	// set zero flow initial conditions
	vessel.setZeroFlowInitialConditions();

	// main computation loop
	vector<vector<double>> A (Nt,vector<double>(Nx,0));
	vector<vector<double>> u (Nt,vector<double>(Nx,0));
	vector<double> t(Nt,0);

	copyVector( vessel.A, A[0] );
	copyVector( vessel.u, u[0] );
	t[0] = 0;


	for(int i = 0; i < Nt-1; i++ )
	{
		// Compute boundaries
		double wmInlet, wpOutlet;

		double wmOutlet = 0;
		vessel.computeOutgoingCharacteristics( dt, wmInlet, wpOutlet );

		// Compute an incomping wavepacket
		double uin = 1*sin(PI*t[i]/0.03)*exp(-pow(t[i]-0.3,2)/0.01);
		double wpInlet = 2*uin-wmInlet;

		double Alb, ulb, Arb, urb;
		vessel.computeBoundaries( wpInlet, wmInlet, wpOutlet, wmOutlet, Alb, ulb, Arb, urb );


		// Update domain
		vessel.update1DVesselLW(dt);

		// update boundaries
		vessel.updateBoundaries( Alb, ulb, Arb, urb );

		copyVector( vessel.A, A[i+1] );
		copyVector( vessel.u, u[i+1] );
		t[i+1] = t[i] + dt;


	}


	// output the result
	if( getProcID() == 0)
	{

		ofstream outFile( outputFileName.c_str() , ios::out );

		printMatlab2DArray( outFile, "A", A );
		printMatlab2DArray( outFile, "u", u );
		printMatlab1DArray( outFile, "t", t );
		printMatlab1DArray( outFile, "x", vessel.x );
		printMatlab1DArray( outFile, "A0", vessel.A0 );
		printMatlab1DArray( outFile, "alpha", vessel.alpha );

	}


}
