/*
 * WavePropagationFunctions.cpp
 *
 *  Created on: Jan 4, 2019
 *      Author: brian
 */

#include "WavePropagationFunctions.hpp"

// Lax-Wendroff finite-difference method to update the 1D wave propagation model*
// Note: boundary should be updated before this function is called
// variables --> 1. A: cross sectional area [cm^2], 2. u: propagation velocity [cm/s]
// dt: time step, dx: constant spatial resolution
void Vessel1D::update1DVesselLW( double dt )
{
	int Nx = A.size();


	// Compute flux/source at grid points
	for(int i = 0; i < Nx; i++)
	{
		computeFlux( A[i], u[i], alpha[i], A0[i], rho, F1[i], F2[i] );
		computeSource( A[i], u[i], mu[i], rho, xi, S1[i], S2[i]);
	}


	// Compute midpoint values
	double coef1 = dt/(2*dx);
	double coef2 = dt/4;
	for(int i = 0; i < Nx-1; i++)
	{
		Amid[i] = 0.5*( A[i+1] + A[i] ) - coef1*( F1[i+1] - F1[i] ) - coef2*( S1[i+1] + S1[i] );
		umid[i] = 0.5*( u[i+1] + u[i] ) - coef1*( F2[i+1] - F2[i] ) - coef2*( S2[i+1] + S2[i] );

		computeFlux( Amid[i], umid[i], alphaMid[i], A0Mid[i], rho, F1mid[i], F2mid[i] );
		computeSource( Amid[i], umid[i], muMid[i], rho, xi, S1mid[i], S2mid[i]);
	}


	// Update regular grid points (not boundary points)
	double coef3 = dt/dx;
	double coef4 = dt/2;
	for(int i = 1; i < Nx-1; i++)
	{
		A[i] = A[i] - coef3*( F1mid[i] - F1mid[i-1] ) - coef4 *( S1mid[i] + S1mid[i-1] );
		u[i] = u[i] - coef3*( F2mid[i] - F2mid[i-1] ) - coef4 *( S2mid[i] + S2mid[i-1] );
	}



}

// Computes flux terms
void computeFlux( double A, double u, double alpha, double A0, double rho, double &F1, double &F2)
{
	F1 = A*u;
	F2 = pow(u,2)/2 + (alpha/rho)*( sqrt( A/A0 ) - 1 );
}

// Compute source terms
void computeSource( double A, double u, double mu, double rho, double xi, double &S1, double &S2 )
{

	S1 = 0;
	S2 = xi*PI*(mu/rho)*(u/A);

}


// Should be computed before the main update
void Vessel1D::computeOutgoingCharacteristics( double dt, double & wmInlet, double & wpOutlet )
{

	vector<double> xShift (NxBndConst, 0);
	vector<double> wm_Inlet_prev (NxBndConst, 0);
	vector<double> wp_Outlet_prev (NxBndConst, 0);

	// Inlet boundary outgoing wave is wm
	for(int i = 0; i < NxBndConst; i++)
	{
		double lambda = u[i] - pow( A[i], 0.25 )*sqrt( R0[i] );
		double xloc = x[i];
		wm_Inlet_prev[i] = u[i] - 4*pow( A[i], 0.25 )*sqrt(R0[i]) + 4*pow(A0[i],0.25)*sqrt(R0[i]);
		xShift[i] = xloc + dt*lambda;
	}
	wmInlet = linearInterp( xShift, wm_Inlet_prev, x[0] );


	// Outgoing boundary outgoing wave is wp
	int j = 0;
	for(int i = Nx-NxBndConst; i < Nx; i++)
	{
		double lambda = u[i] + pow( A[i], 0.25 )*sqrt( R0[i] );
		double xloc = x[i];
		wp_Outlet_prev[j] = u[i] + 4*pow( A[i], 0.25 )*sqrt(R0[i]) - 4*pow(A0[i],0.25)*sqrt(R0[i]);
		xShift[j] = xloc + dt*lambda;
		j = j+1;

	}
	wpOutlet = linearInterp( xShift, wp_Outlet_prev, x[Nx-1] );

}





void Vessel1D::createConstantVessel( int NxIn, double xlength, double A0in, double nuIn, double xiIn, double rhoIn, double EIn, double muIn, double h0In  )
{
	// Number of boundary points used to compute characteristics, since this is a constant vessel the size here doesn't matter
	if(Nx < 10)
	{
		NxBndConst = Nx/2;
	}
	else
	{
		NxBndConst = 5;
	}

	// constant calculations
	Nx = NxIn;
	rho = rhoIn;
	xi = xiIn;
	dx = xlength/(Nx-1);
	double E = EIn;
	double h0 = h0In;
	double nu = nuIn;


	// Reset variable sizes
	x.resize(Nx);
	A.resize(Nx);
	u.resize(Nx);
	A0.resize(Nx);
	alpha.resize(Nx);
	mu.resize(Nx);
	F1.resize(Nx);
	F2.resize(Nx);
	S1.resize(Nx);
	S2.resize(Nx);
	R0.resize(Nx);

	Amid.resize(Nx-1);
	umid.resize(Nx-1);
	A0Mid.resize(Nx-1);
	alphaMid.resize(Nx-1);
	muMid.resize(Nx-1);
	F1mid.resize(Nx-1);
	F2mid.resize(Nx-1);
	S1mid.resize(Nx-1);
	S2mid.resize(Nx-1);


	for(int i = 0; i < Nx; i++ )
	{
		x[i] = i*dx;
		A0[i] = A0in;
		mu[i] = muIn;
		alpha[i] = sqrt( PI/A0[i] )*( E*h0 )/( 1.0 - pow(nu,2) );
		R0[i] = alpha[i]/(2*rho*sqrt(A0[i]));
	}

	for(int i = 0; i < Nx-1; i++ )
	{
		A0Mid[i] = A0in;
		muMid[i] = muIn;
		alphaMid[i] =  sqrt( PI/A0Mid[i] )*( E*h0 )/( 1.0 - pow(nu,2) );
	}


}


void Vessel1D::computeBoundaries( double wpIn, double wmIn, double wpOut, double wmOut,
		double & Alb, double & ulb, double & Arb, double & urb )
{
	Alb = pow( (wpIn - wmIn)/(8*sqrt(R0[0])) +  pow(A0[0],0.25), 4 );
	ulb = ( wpIn + wmIn )/2;

	Arb = pow( (wpOut - wmOut)/(8*sqrt(R0[Nx-1])) +  pow(A0[Nx-1],0.25), 4 );
	urb = ( wpOut + wmOut )/2;
}



void Vessel1D::updateBoundaries( double Alb, double ulb, double Arb, double urb  )
{
	A[0] = Alb;
	u[0] = ulb;

	A[Nx-1] = Arb;
	u[Nx-1] = urb;
}









