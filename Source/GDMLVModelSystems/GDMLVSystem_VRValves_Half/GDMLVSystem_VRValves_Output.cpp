/*
 * GDMLVSystem_VRValves_Output.cpp
 *
 *  Created on: Dec 19, 2018
 *      Author: brian
 */



#include "GDMLVSystem_VRValves_Output.hpp"




void saveDataValuesGDMLVHalf( int stepIdx, double t, vector<double>&  q, vector<double> & w, vector<double>&  aux,
		vector<double> & tStore, vector<vector<double>> & varStore, ParamGdmlvAlt & altPrm, double Plv, double Vlv, double Pla, double Ppao  )
{



	// store values
	// compute aux variables
	double Rmv, Raov;
	valveResistances( Pla, Plv, Ppao, altPrm.beta, altPrm.Rmvc, altPrm.Rmvo,
					altPrm.Raovc, altPrm.Raovo, Rmv, Raov );

	double At_ventricle = twoHillActivation( t, altPrm.m1_ventricle, altPrm.m2_ventricle, altPrm.tau1_ventricle,
			altPrm.tau2_ventricle, altPrm.Tc, altPrm.Ts_ventricle, altPrm.hillMaxVal_ventricle );
	double At_atrium = twoHillActivation( t, altPrm.m1_atrium, altPrm.m2_atrium, altPrm.tau1_atrium,
			altPrm.tau2_atrium, altPrm.Tc, altPrm.Ts_atrium, altPrm.hillMaxVal_atrium );
	double qmv = (Pla - Plv)/Rmv;
	double qaov = (Plv - Ppao)/Raov;
	double qinla = ( altPrm.Ppv_fixed - Pla )/ altPrm.Rpv;

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


void exportDataValuesGDMLVHalf( string outFileName, int Ncomputed, int Nq, int Nw, int Naux,
		vector<double> & tStore, vector<vector<double>> & varStore )
{


	if( getProcID() == 0)
	{

		cout << "Printing output to file " << outFileName << endl;
		ofstream outFile( outFileName.c_str() , ios::out );

		int idx = 0;
		vector<vector<double>> qStore (  Nq, vector<double>(Ncomputed,0) );

		for(int i = 0; i < Nq; i++)
		{
			printMatlabArraySimple( outFile, "q"+to_string(i), varStore[idx].data(), Ncomputed);

			for(int j = 0; j < Ncomputed; j++)
			{
				qStore[i][j] = varStore[idx][j];
				//cout << j << "  " << qStore[i][j] << endl;
			}

			idx++;
		}

		printMatlab2DArray( outFile, "qStore",  qStore );


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



void exportInitialValuesGDMLVHalf( string initialValuesFileName, vector<double> &X, vector<double> & q, vector<double> &w  )
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



void readInitialValuesGDMLVHalf( string initialValuesFileName, vector<double> &X, vector<double> & q, vector<double> &w  )
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









