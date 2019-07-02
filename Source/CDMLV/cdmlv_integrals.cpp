/*
 * cdmlv_integrals.cpp
 *
 *  Created on: Jan 25, 2019
 *      Author: brian
 */

#include "cdmlv_integrals.hpp"


void computeRegularIntegrals( vector <CDMNode> & localNodes, int Nnodes, vector<vector<double>> & alpha,
				vector<double> & kappa, vector<double> & etae, vector<double> & dVregdq, double & Vreg, ParamCDMLV & prm )
{

	computeStaticRegularIntegrals( localNodes, Nnodes, kappa, etae, Vreg, prm );

	computeDynamicRegularIntegrals( localNodes, Nnodes, alpha, dVregdq, prm );


}


void computeStaticRegularIntegrals( vector <CDMNode> & localNodes, int Nnodes, vector<double> & kappa,
		vector<double> & etae, double & Vreg, ParamCDMLV & prm )
{


	vector<double> vtemp(3,0);
	vector<vector<double>> Atemp(3, vector<double>(3,0));

	// Set to zero before summation
	Vreg = 0.0;
	for(int i = 0; i < prm.Nq; i++)
	{
		etae[i] = 0;
		kappa[i] = 0;
	}

	// Sum values that were computed on the local nodes
	for (int k = 0; k < Nnodes; k++)
	{
		//compute kappa
		for(int i = 0; i < prm.Nq; i++)
		{
			for (int ii = 0; ii < 3; ii++){
				for(int jj = 0; jj < 3; jj++){
					Atemp[ii][jj] = localNodes[k].Se[ii][jj] + localNodes[k].Sa0[ii][jj];
				}
			}
			kappa[i] = kappa[i] + localNodes[k].IV0weight*matrixDoubleDot( localNodes[k].dEdq[i], Atemp );
		}

		// add forcing terms to kappa
		for(int i = 0; i < prm.Nq; i++)
		{
			// body force (opposite sign, other side of the equation)
			kappa[i] = kappa[i] - localNodes[k].IV0weight*dotProd( localNodes[k].bodyForcing, localNodes[k].dudq[i] );

			// sides surface integral (opposite sign, other side of the equation)
			// forcing term is zero everywhere except on the surfaces
			kappa[i] = kappa[i] - localNodes[k].IA0weight*localNodes[k].vCrossSideNorm*dotProd( localNodes[k].tractionForcingSide, localNodes[k].dudq[i] );

			// top/bottom surface integral (opposite sign, other side of the equation)
			// forcing term is zero everywhere except on the surfaces
			kappa[i] = kappa[i] - localNodes[k].IA0weight_top*localNodes[k].vCrossTopNorm*dotProd( localNodes[k].tractionForcingTop, localNodes[k].dudq[i] );
		}

		// compute etae (only on endocardial surface nodes)
		if( localNodes[k].endoSurfNode == 1 )
		{
			for(int i = 0; i < prm.Nq; i++)
			{

				etae[i] = etae[i] + localNodes[k].IA0weight* dotProd( localNodes[k].dudq[i], localNodes[k].vcrossSide );
			}
		}

		// Compute volume and volume derivatives
		if( localNodes[k].endoSurfNode == 1 )
		{
			Vreg = Vreg + localNodes[k].VI*localNodes[k].IA0weight;
		}

	}

}


void computeDynamicRegularIntegrals( vector <CDMNode> & localNodes, int Nnodes, vector<vector<double>> & alpha, vector<double> & dVregdq, ParamCDMLV & prm )
{

	vector<double> vtemp(3,0);
	vector<vector<double>> Atemp(3, vector<double>(3,0));

	// Set to zero before summation
	for(int i = 0; i < prm.Nq; i++)
	{
		dVregdq[i] = 0;
		for(int j = 0; j < prm.Nq; j++)
		{
			alpha[i][j] = 0.0;
		}
	}

	// Sum values that were computed on the local nodes
	for (int k = 0; k < Nnodes; k++)
	{
		// compute alpha
		for(int i = 0; i < prm.Nq; i++)
		{
			for(int j = 0; j < prm.Nq; j++)
			{

				for (int ii = 0; ii < 3; ii++){
					for(int jj = 0; jj < 3; jj++){
						Atemp[ii][jj] = localNodes[k].Svq[j][ii][jj] + localNodes[k].Saq[j][ii][jj];

					}
				}
				alpha[i][j] = alpha[i][j] + localNodes[k].IV0weight*matrixDoubleDot( localNodes[k].dEdq[i] , Atemp );
			}
		}

		// Compute volume and volume derivatives
		if( localNodes[k].endoSurfNode == 1 )
		{
			for(int i = 0; i < prm.Nq; i++)
			{
				dVregdq[i] = dVregdq[i] + localNodes[k].dVIdq[i]*localNodes[k].IA0weight;
			}
		}

	}
}


double computeReferenceWallVolume ( vector <CDMNode> & localNodes, int Nnodes )
{
	double wallVolumeLocal = 0;

	for (int k = 0; k < Nnodes; k++)
	{
		wallVolumeLocal = wallVolumeLocal + localNodes[k].IV0weight;
	}

	double wallVolume;
	MPI_Reduce( &wallVolumeLocal, &wallVolume, 1, MPI_DOUBLE, MPI_SUM, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast( &wallVolume, 1, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);

	return wallVolume;

}


