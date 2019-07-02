/*
 * GDMLVSystem_VRValves_HalfEq.cpp
 *
 *  Created on: Dec 21, 2018
 *      Author: brian
 */


#include "GDMLVSystem_VRValves_HalfEq.hpp"


void computeEquilibriumCycleGdmLVAltRecompute( double volEqRequirement, double minEquilibriumCycles, double maxEquilibriumCycles, int maxRecomputes,
		vector <double> & q, vector <double> & w, vector <double> & X,	vector<double> & tStore, vector<double> & PlvStore,
		vector <vector <double> > & qStore,	GdmLV & lv, ParamGdmlvAlt & altPrm, DataContainer & prmData )
{

	int	K = 0;
	lv.prm.Ncyc = 1; // only computes one cycle at each iteration

	double dtCycleStandard = lv.prm.dtMax;
	double dt = dtCycleStandard;
	int solutionSuccess = 0;

	while ( (!solutionSuccess && K < maxRecomputes) || K == 0)
	{
		if (K > 0)
		{
			dt = dtCycleStandard/pow(2,K);
			// preloop computations required to recompute time step
			cout << "Solution on procID = " << getProcID() << " failed for the " << K << " time, recomputing with dt = " << dt  << endl;
		}

		// compute the error
		lv.prm.dtMax = dt;
		lv.prm.dtMin = dt;

		solutionSuccess = computeEquilibriumCycleGdmLVAlt( volEqRequirement, minEquilibriumCycles, maxEquilibriumCycles,  q, w, X, tStore, PlvStore, qStore, lv, altPrm, prmData );

		K = K+1;
	}

	lv.prm.dtMax = dtCycleStandard;
	lv.prm.dtMin = dtCycleStandard;

}





int computeEquilibriumCycleGdmLVAlt( double volEqRequirement, double minEquilibriumCycles, double maxEquilibriumCycles, vector <double> & q,
		vector <double> & w, vector <double> & X, vector<double> & tStore, vector<double> & PlvStore, vector <vector <double> > & qStore,
		GdmLV & lv, ParamGdmlvAlt & altPrm, DataContainer & prmData  )
{
	int verbose = 0;
	int solutionSuccess = 1;

	vector<double> q0 (q.size(),0);
	vector<double> w0 (w.size(),0);
	vector<double> X0 (X.size(),0);

	for(int i = 0; i < q.size(); i++){ q0[i] = q[i]; }
	for(int i = 0; i < w.size(); i++){ w0[i] = w[i]; }
	for(int i = 0; i < X.size(); i++){ X0[i] = X[i]; }

	int Nsteps = ceil( (lv.prm.Ncyc*lv.prm.Tc)/lv.prm.dtMin ) + 2 ;

	int Naux = 11;
	vector<vector<double>> varStore( q.size()+ w.size() + Naux, vector<double>(Nsteps,0));

	qStore.resize(q.size());
	for(int i = 0; i < qStore.size();i++){qStore[i].resize(Nsteps);}
	tStore.resize(Nsteps);
	PlvStore.resize(Nsteps);

	// Get the number of processes
	int Nprocs,procID;
	MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	// local variables
	vector <double> VlvCurrent(Nsteps, 0);
	vector <double> VlvPrev (Nsteps, 0);
	vector <double> VlvAbsDiff (Nsteps, 0);

	// run cycles until equilibrium or timeout
	double eqError = volEqRequirement*2;
	int i = 0;

	while ( ( i < maxEquilibriumCycles && eqError > volEqRequirement ) || ( i < 3 ) || (i < minEquilibriumCycles) )
	{

		int NtComputed;
		computeLVSystemHalfVectorOutput( Naux, X, q, w, lv, altPrm, tStore, varStore, NtComputed );

		// compute error as volume change
		for(int k = 0; k < Nsteps; k++)
		{
			// aux = [ Vlv, Plv, Pla, Ppao, At_ventricle, At_atrium, Rmv, Raov, qmv, qaov, qinla ]
			VlvCurrent[k] = varStore[q.size()+ w.size() + 0][k];
			VlvAbsDiff[k] = fabs( VlvPrev[k] - VlvCurrent[k] )/Nsteps;
		}
		eqError = vectorNorm ( VlvAbsDiff.data(), Nsteps, 1 );

		if(getProcID() == 0 && verbose > 0)
		{
			cout << endl << "---------------------------------------------------------------" << endl;
			cout << "After cycle " << i << " out of maximum " <<  maxEquilibriumCycles << endl;
			cout << "VlvEqDiff = " << eqError << " with goal to drop below " << volEqRequirement << endl;
			cout << "Current time step is " << lv.prm.dtMax << endl;
			cout << "---------------------------------------------------------------" << endl << endl;
		}
		for(int k = 0; k < Nsteps; k++)
		{
			VlvPrev[k] = VlvCurrent[k];
		}



		// Exit solution loop if solver fails
		if( w[0] != w[0] )
		{

			solutionSuccess = 0;
			break;
		}

		i = i+1;
	}
	// cout << "Exiting solution loop after " << i << " cycles with error " << eqError << " while the tolerance is set to " << prmHSA.volumeEquilibrium  << endl;


	// if solver crashed reset initial values
	if( w[0] != w[0] )
	{

		for(int i = 0; i < q.size(); i++){ q[i] = q0[i]; }
		for(int i = 0; i < w.size(); i++){ w[i] = w0[i]; }
		for(int i = 0; i < X.size(); i++){ X[i] = X0[i]; }

	}

	for(int k = 0; k < Nsteps; k++)
	{
		// varStore = [qStore, wStore, auxStore]
		for(int i = 0; i < q.size(); i++)
		{
			qStore[i][k] = varStore[i][k];
		}
		// aux = [ Vlv, Plv, Pla, Ppao, At_ventricle, At_atrium, Rmv, Raov, qmv, qaov, qinla ]
		VlvCurrent[k] = varStore[q.size()+ w.size() + 0][k];
		PlvStore[k] = varStore[q.size()+ w.size() + 1][k];

	}


	return solutionSuccess;

}

