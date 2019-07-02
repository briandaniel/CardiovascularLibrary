/*
 * CDMLV_HalfSystem_Output.cpp
 *
 *  Created on: Jan 29, 2019
 *      Author: brian
 */

#include "CDMLV_HalfSystem_Output.hpp"



void saveDataValuesCDMLVHalf( int stepIdx, double t, vector<double>&  q, vector<double> & w, vector<double>&  aux,
		vector<double> & tStore, vector<vector<double>> & varStore, ParamLumpActive & lumpPrm, double Plv, double Vlv, double Pla, double Ppao  )
{



	// store values
	// compute aux variables
	double Rmv, Raov;
	valveResistances( Pla, Plv, Ppao, lumpPrm.beta, lumpPrm.Rmvc, lumpPrm.Rmvo,
					lumpPrm.Raovc, lumpPrm.Raovo, Rmv, Raov );

	double At_ventricle = twoHillActivation( t, lumpPrm.m1_ventricle, lumpPrm.m2_ventricle, lumpPrm.tau1_ventricle,
			lumpPrm.tau2_ventricle, lumpPrm.Tc, lumpPrm.Ts_ventricle, lumpPrm.hillMaxVal_ventricle );
	double At_atrium = twoHillActivation( t, lumpPrm.m1_atrium, lumpPrm.m2_atrium, lumpPrm.tau1_atrium,
			lumpPrm.tau2_atrium, lumpPrm.Tc, lumpPrm.Ts_atrium, lumpPrm.hillMaxVal_atrium );
	double qmv = (Pla - Plv)/Rmv;
	double qaov = (Plv - Ppao)/Raov;
	double qinla = ( lumpPrm.Ppv_fixed - Pla )/ lumpPrm.Rpv;

	// aux = [ Vlv, Plv, Pla, Ppao, At_ventricle, At_atrium, Rmv, Raov, qmv, qaov, qinla ]
	aux[0] = Vlv;
	aux[1] = Plv;
	aux[2] = Pla;
	aux[3] = Ppao;
	aux[4] = At_ventricle;
	aux[5] = At_atrium;
	aux[6] = Rmv;
	aux[7] = Raov;
	aux[8] = qmv;
	aux[9] = qaov;
	aux[10] = qinla;


	int i = 0;
	for(int k = 0; k < q.size(); k++)
	{
		varStore[i][stepIdx] = q[k];
		i++;
	}
	for(int k = 0; k < w.size(); k++)
	{
		varStore[i][stepIdx] = w[k];
		i++;
	}
	for(int k = 0; k < aux.size(); k++)
	{
		varStore[i][stepIdx] = aux[k];
		i++;
	}
	tStore[stepIdx] = t;

}


void exportDataValuesCDMLVHalf( string outFileName, int Ncomputed, int Nq, int Nw, int Naux, vector<double> & tStore, vector<vector<double>> & varStore )
{


	if( getProcID() == 0)
	{


		ofstream outFile( outFileName.c_str() , ios::out );

		int idx = 0;
		for(int i = 0; i < Nq; i++)
		{
			printMatlabArraySimple( outFile, "q"+to_string(i), varStore[idx].data(), Ncomputed); idx++;
		}

		// W = [ Vla, Psa ]
		printMatlabArraySimple( outFile, "Vla", varStore[idx].data(), Ncomputed); idx++;
		printMatlabArraySimple( outFile, "Psa", varStore[idx].data(), Ncomputed); idx++;

		// Aux vars
		printMatlabArraySimple( outFile, "Vlv", varStore[idx].data(), Ncomputed); idx++;
		printMatlabArraySimple( outFile, "Plv", varStore[idx].data(), Ncomputed); idx++;
		printMatlabArraySimple( outFile, "Pla", varStore[idx].data(), Ncomputed); idx++;
		printMatlabArraySimple( outFile, "Ppao", varStore[idx].data(), Ncomputed); idx++;
		printMatlabArraySimple( outFile, "At_ventricle", varStore[idx].data(), Ncomputed); idx++;
		printMatlabArraySimple( outFile, "At_atrium", varStore[idx].data(), Ncomputed); idx++;
		printMatlabArraySimple( outFile, "Rmv", varStore[idx].data(), Ncomputed); idx++;
		printMatlabArraySimple( outFile, "Raov", varStore[idx].data(), Ncomputed); idx++;
		printMatlabArraySimple( outFile, "qmv", varStore[idx].data(), Ncomputed); idx++;
		printMatlabArraySimple( outFile, "qaov", varStore[idx].data(), Ncomputed); idx++;
		printMatlabArraySimple( outFile, "qinla", varStore[idx].data(), Ncomputed); idx++;

		printMatlabArraySimple( outFile, "t", tStore.data(), Ncomputed);

		outFile.close();
	}


}



void exportInitialValuesCDMLVHalf( string initialValuesFileName, vector<double> &X, vector<double> & q, vector<double> &w  )
{

	if( getProcID() == 0)
	{

		ofstream outFile( initialValuesFileName.c_str() , ios::out );

		printMatlabArraySimple( outFile, "X", X.data(), X.size());
		printMatlabArraySimple( outFile, "q", q.data(), q.size());
		printMatlabArraySimple( outFile, "w", w.data(), w.size());

		outFile.close();
	}

}



void readInitialValuesCDMLVHalf( string initialValuesFileName, vector<double> &X, vector<double> & q, vector<double> &w  )
{

	if( getProcID() == 0)
	{

		ifstream inFile( initialValuesFileName.c_str() , ios::out );

		findArray( inFile, "X", X.size(), X.data(), 1 );
		findArray( inFile, "q", q.size(), q.data(), 1 );
		findArray( inFile, "w", w.size(), w.data(), 1 );

		inFile.close();

	}

    MPI_Bcast( X.data(), X.size(), MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD );
    MPI_Bcast( q.data(), q.size(), MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD );
    MPI_Bcast( w.data(), w.size(), MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD );


}





