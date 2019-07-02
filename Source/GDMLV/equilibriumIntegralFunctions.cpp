/*
 * equilibriumIntegralFunctions.cpp
 *
 *  Created on: Sep 11, 2018
 *      Author: brian
 *
 *
 * copied from:
 * equilibriumFunctions.cpp
 *
 *  Created on: Sep 15, 2017
 *      Author: brian
 */


#include "equilibriumIntegralFunctions.hpp"



void computeProcessorIntegrals( vector <NodeGdmLV> & localNodes, int Nnodes, double ** alpha,
				double *kappa, double * etae, double *dVregdq, double & Vreg, ParamGdmLV * prm, double Prv )
{


	double * vtemp = new double [3];
	double ** Atemp = new double * [3];
	double * kappaRV = new double [prm->Nq];

	for(int i = 0; i < 3; i++){ Atemp[i] = new double [3]; }

	// Set to zero before summation
	Vreg = 0.0;
	for(int i = 0; i < prm->Nq; i++)
	{
		kappaRV[i] = 0;
		etae[i] = 0;
		dVregdq[i] = 0;
		kappa[i] = 0;
		for(int j = 0; j < prm->Nq; j++)
		{
			alpha[i][j] = 0.0;
		}
	}



	// Sum values that were computed on the local nodes
	for (int k = 0; k < Nnodes; k++)
	{

		// compute alpha
		for(int i = 0; i < prm->Nq; i++)
		{
			for(int j = 0; j < prm->Nq; j++)
			{

				for (int ii = 0; ii < 3; ii++){
					for(int jj = 0; jj < 3; jj++){
						Atemp[ii][jj] = localNodes[k].Svq[j][ii][jj] + localNodes[k].Saq[j][ii][jj];

					}
				}

				alpha[i][j] = alpha[i][j] + localNodes[k].IV0weight*matrixDoubleDot( localNodes[k].dEdq[i] , Atemp, 3, 3);
			}
		}



		// compute kappa
		for(int i = 0; i < prm->Nq; i++)
		{
			for (int ii = 0; ii < 3; ii++){
				for(int jj = 0; jj < 3; jj++){
					Atemp[ii][jj] = localNodes[k].Se[ii][jj] + localNodes[k].Sa0[ii][jj];
				}
			}

			kappa[i] = kappa[i] + localNodes[k].IV0weight*matrixDoubleDot( localNodes[k].dEdq[i], Atemp, 3, 3);
		}


		// add forcing terms to kappa
		for(int i = 0; i < prm->Nq; i++)
		{
			// body force (opposite sign, other side of the equation)
			kappa[i] = kappa[i] - localNodes[k].IV0weight*dotProd( localNodes[k].bodyForcing, localNodes[k].dudq[i], 3);

			// sides surface integral (opposite sign, other side of the equation)
			// forcing term is zero everywhere except on the surfaces
			kappa[i] = kappa[i] - localNodes[k].IA0weight*localNodes[k].vCrossSideNorm*dotProd( localNodes[k].tractionForcingSide, localNodes[k].dudq[i], 3 );

			// top/bottom surface integral (opposite sign, other side of the equation)
			// forcing term is zero everywhere except on the surfaces
			kappa[i] = kappa[i] - localNodes[k].IA0weight_top*localNodes[k].vCrossTopNorm*dotProd( localNodes[k].tractionForcingTop, localNodes[k].dudq[i], 3 );


		}


















		// Compute rv contribution at the current node, if required
		if (prm->useRVPressure == 1)
		{
			if( localNodes[k].rvSurfNode == 1 )
			{
				// reverse the direction of the vector for the RV surface
				for(int j = 0; j < 3; j++ )
					vtemp[j] = - localNodes[k].vcrossSide[j];

				for(int i = 0; i < prm->Nq; i++)
				{
					// this is only the local node contribution, it is summed inside kappa below
					double etarv_node = localNodes[k].IA0weight* dotProd( localNodes[k].dudq[i], vtemp, 3 );

					kappaRV[i] = kappaRV[i] - Prv*etarv_node;
				}
			}
		}
		else
		{
			for(int i = 0; i < prm->Nq; i++)
			{
				kappaRV[i] = 0.0;
			}
		}



		// compute etae (only on endocardial surface nodes)
		if( localNodes[k].endoSurfNode == 1 )
		{
			for(int i = 0; i < prm->Nq; i++)
			{
				/*
				double vNorm = sqrt( pow( localNodes[k].vcross[0], 2) + pow( localNodes[k].vcross[1], 2) + pow( localNodes[k].vcross[2], 2) );
				double fVal = sin(4*localNodes[k].phi);
				etae[i] = etae[i] + localNodes[k].IA0weight*vNorm;
				 */

				etae[i] = etae[i] + localNodes[k].IA0weight* dotProd( localNodes[k].dudq[i], localNodes[k].vcrossSide, 3 );



			}
		}


		// Compute volume and volume derivatives
		if( localNodes[k].endoSurfNode == 1 )
		{
			Vreg = Vreg + localNodes[k].VI*localNodes[k].IA0weight;
			for(int i = 0; i < prm->Nq; i++)
			{
				dVregdq[i] = dVregdq[i] + localNodes[k].dVIdq[i]*localNodes[k].IA0weight;
			}
		}

	}


	// add the contribution of the RV to the regular kappa
	for(int i = 0; i < prm->Nq; i++)
	{
		kappa[i] = kappa[i] + kappaRV[i];
	}

	// clean up
	for(int i = 0; i < 3; i++){ delete [] Atemp[i]; }
	delete [] Atemp;
	delete [] kappaRV;
	delete [] vtemp;
}






void computeStaticProcessorIntegrals( vector <NodeGdmLV> & localNodes, int Nnodes,
				double *kappa, double * etae, double & Vreg, ParamGdmLV * prm, double Prv )
{


	double * vtemp = new double [3];
	double ** Atemp = new double * [3];
	double * kappaRV = new double [prm->Nq];

	for(int i = 0; i < 3; i++){ Atemp[i] = new double [3]; }

	// Set to zero before summation
	Vreg = 0.0;
	for(int i = 0; i < prm->Nq; i++)
	{
		kappaRV[i] = 0;
		etae[i] = 0;
		kappa[i] = 0;
	}



	// Sum values that were computed on the local nodes
	for (int k = 0; k < Nnodes; k++)
	{


		// compute kappa
		for(int i = 0; i < prm->Nq; i++)
		{
			for (int ii = 0; ii < 3; ii++){
				for(int jj = 0; jj < 3; jj++){
					Atemp[ii][jj] = localNodes[k].Se[ii][jj] + localNodes[k].Sa0[ii][jj];
				}
			}

			kappa[i] = kappa[i] + localNodes[k].IV0weight*matrixDoubleDot( localNodes[k].dEdq[i], Atemp, 3, 3);
		}



		// add forcing terms to kappa
		for(int i = 0; i < prm->Nq; i++)
		{
			// body force (opposite sign, other side of the equation)
			kappa[i] = kappa[i] - localNodes[k].IV0weight*dotProd( localNodes[k].bodyForcing, localNodes[k].dudq[i], 3);

			// sides surface integral (opposite sign, other side of the equation)
			// forcing term is zero everywhere except on the surfaces
			kappa[i] = kappa[i] - localNodes[k].IA0weight*localNodes[k].vCrossSideNorm*dotProd( localNodes[k].tractionForcingSide, localNodes[k].dudq[i], 3 );

			// top/bottom surface integral (opposite sign, other side of the equation)
			// forcing term is zero everywhere except on the surfaces
			kappa[i] = kappa[i] - localNodes[k].IA0weight_top*localNodes[k].vCrossTopNorm*dotProd( localNodes[k].tractionForcingTop, localNodes[k].dudq[i], 3 );


		}
















		// Compute rv contribution at the current node, if required
		if (prm->useRVPressure == 1)
		{
			if( localNodes[k].rvSurfNode == 1 )
			{
				// reverse the direction of the vector for the RV surface
				for(int j = 0; j < 3; j++ )
					vtemp[j] = - localNodes[k].vcrossSide[j];

				for(int i = 0; i < prm->Nq; i++)
				{
					// this is only the local node contribution, it is summed inside kappa below
					double etarv_node = localNodes[k].IA0weight* dotProd( localNodes[k].dudq[i], vtemp, 3 );

					kappaRV[i] = kappaRV[i] - Prv*etarv_node;
				}
			}
		}



		// compute etae (only on endocardial surface nodes)
		if( localNodes[k].endoSurfNode == 1 )
		{
			for(int i = 0; i < prm->Nq; i++)
			{

				etae[i] = etae[i] + localNodes[k].IA0weight* dotProd( localNodes[k].dudq[i], localNodes[k].vcrossSide, 3 );
			}
		}

		// Compute volume and volume derivatives
		if( localNodes[k].endoSurfNode == 1 )
		{
			Vreg = Vreg + localNodes[k].VI*localNodes[k].IA0weight;
		}


	}


	// add the contribution of the RV to the regular kappa
	for(int i = 0; i < prm->Nq; i++)
	{
		kappa[i] = kappa[i] + kappaRV[i];
	}


	// clean up
	for(int i = 0; i < 3; i++){ delete [] Atemp[i]; }
	delete [] Atemp;
	delete [] kappaRV;
	delete [] vtemp;
}



double computeReferenceWallVolume ( vector <NodeGdmLV> & localNodes, int Nnodes )
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



