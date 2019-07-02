/*
 * cdmlv_node_fiber_directions.cpp
 *
 *  Created on: Jan 25, 2019
 *      Author: brian
 */


#include "cdmlv_node.hpp"


void CDMNode::computeFiberRotation( double * es, double * en, double * ef )
{


	double* emu = new double[3];
	double* enu = new double[3];
	double* ephi = new double[3];

    // 3. Compute prolate spheroid coordinate basis vectors
    emu[0] = a*chmu0*snu0*cphi0;
    emu[1] = a*chmu0*snu0*sphi0;
    emu[2] = a*shmu0*cnu0;

    enu[0] = a*shmu0*cnu0*cphi0;
    enu[1] = a*shmu0*cnu0*sphi0;
    enu[2] = -a*chmu0*snu0;

    ephi[0] = -a*shmu0*snu0*sphi0;
    ephi[1] = a*shmu0*snu0*cphi0;
    ephi[2] = 0;

	normalizeVector2norm( emu, 3 );
	normalizeVector2norm( enu, 3 );
	normalizeVector2norm( ephi, 3 );

	// 4. Compute the rotation matrix
	Q[0][0] = dotProd( emu, es, 3);
	Q[0][1] = dotProd( emu, en, 3);
	Q[0][2] = dotProd( emu, ef, 3);

	Q[1][0] = dotProd( enu, es, 3);
	Q[1][1] = dotProd( enu, en, 3);
	Q[1][2] = dotProd( enu, ef, 3);

	Q[2][0] = dotProd( ephi, es, 3);
	Q[2][1] = dotProd( ephi, en, 3);
	Q[2][2] = dotProd( ephi, ef, 3);

	matrixTranspose( Q, QT );

	// clean up
	delete[] emu;
	delete[] enu;
	delete[] ephi;

}



void CDMNode::computeFiberRotation( ParamCDMLV & prm )
{

	double* e1 = new double[3];
	double* e2 = new double[3];
	double* e3 = new double[3];
	double* e4 = new double[3];

	double* es = new double[3];
	double* en = new double[3];
	double* ef = new double[3];

	double* emu = new double[3];
	double* enu = new double[3];
	double* ephi = new double[3];


	computeFiberVectors( e1, e2, e3, e4, es, en, ef, prm );

    // 3. Compute prolate spheroid coordinate basis vectors
    emu[0] = a*chmu0*snu0*cphi0;
    emu[1] = a*chmu0*snu0*sphi0;
    emu[2] = a*shmu0*cnu0;

    enu[0] = a*shmu0*cnu0*cphi0;
    enu[1] = a*shmu0*cnu0*sphi0;
    enu[2] = -a*chmu0*snu0;

    ephi[0] = -a*shmu0*snu0*sphi0;
    ephi[1] = a*shmu0*snu0*cphi0;
    ephi[2] = 0;


	normalizeVector2norm( emu, 3 );
	normalizeVector2norm( enu, 3 );
	normalizeVector2norm( ephi, 3 );

	/*
	for(int i = 0; i < 3; i++)
	{
		ef[i] = ephi[i];
		es[i] = enu[i];
		en[i] = emu[i];
	}
	*/

	// 4. Compute the rotation matrix
	Q[0][0] = dotProd( emu, es, 3);
	Q[0][1] = dotProd( emu, en, 3);
	Q[0][2] = dotProd( emu, ef, 3);

	Q[1][0] = dotProd( enu, es, 3);
	Q[1][1] = dotProd( enu, en, 3);
	Q[1][2] = dotProd( enu, ef, 3);

	Q[2][0] = dotProd( ephi, es, 3);
	Q[2][1] = dotProd( ephi, en, 3);
	Q[2][2] = dotProd( ephi, ef, 3);

	matrixTranspose( Q, QT );


	delete[] e1;
	delete[] e2;
	delete[] e3;
	delete[] e4;
	delete[] es;
	delete[] en;
	delete[] ef;

	delete[] emu;
	delete[] enu;
	delete[] ephi;

}




void CDMNode::computeFiberVectors( double * e1, double * e2, double * e3, double * e4, double * es, double * en, double * ef, ParamCDMLV & prm )
{

	double dmu0dnu0, dmu0dphi0, dx0dnu0, dy0dnu0, dz0dnu0, dx0dphi0, dy0dphi0, dz0dphi0;
	double at, bt, ct, stheta, ctheta, dotProde13, w,  ut, ut0, G, psiTemp, signPsi, psi;
	double dotProde14, cpsi, spsi;

	double* crossProd = new double[3];
	double* crossProde13 = new double[3];
	double* crossProde14 = new double[3];

	double* emu = new double[3];
	double* enu = new double[3];
	double* ephi = new double[3];

	// 1. compute basis vectors
	dmu0dnu0 = (1-umu0)*dmuin0dnu0 + umu0*dmuout0dnu0;
	dmu0dphi0 = (1-umu0)*dmuin0dphi0 + umu0*dmuout0dphi0;

	dx0dnu0 = a*cphi0*( chmu0*dmu0dnu0*snu0 + cnu0*shmu0 );
	dy0dnu0 = a*sphi0*( chmu0*dmu0dnu0*snu0 + cnu0*shmu0 );
	dz0dnu0 = a*( shmu0*dmu0dnu0*cnu0 - chmu0*snu0 );

	dx0dphi0 = a*snu0*( chmu0*dmu0dphi0*cphi0 - sphi0*shmu0 );
	dy0dphi0 = a*snu0*( chmu0*dmu0dphi0*sphi0 + cphi0*shmu0 );
	dz0dphi0 = a*cnu0*shmu0*dmu0dphi0;

	e2[0] = dx0dnu0;
	e2[1] = dy0dnu0;
	e2[2] = dz0dnu0;

	e3[0] = dx0dphi0;
	e3[1] = dy0dphi0;
	e3[2] = dz0dphi0;


	// correct values at apex
	if ( fabs(nu0 - PI) < 1e-8 ){
		e3[0] = -sphi0;
		e3[1] = cphi0;
		e3[2] = 0;
	}

	normalizeVector2norm( e2, 3 );
	normalizeVector2norm( e3, 3 );


	crossProduct3( e2, e3, e1 );

	at = e3[2] - e1[2]*dotProd( e3, e1, 3 );
	crossProduct3( e1, e3, crossProd);
	bt = crossProd[2];
	ct = e1[2]*dotProd(e3,e1,3);

    ctheta = - ( at*ct - sqrt(pow(bt,2)*( pow(at,2) + pow(bt,2) - pow(ct,2) )) )/( pow(at,2)+pow(bt,2));
    stheta = - ( pow(bt,2)*pow(ct,2) + at*sqrt( pow(bt,2)*(pow(at,2) + pow(bt,2) - pow(ct,2) )))/( bt*(pow(at,2)+ pow(bt,2)) );

    dotProde13 = dotProd(e1,e3,3);
    crossProduct3( e1, e3, crossProde13 );

    // compute e4
    for (int k = 0; k < 3; k++)
    	e4[k] = e3[k]*ctheta + crossProde13[k]*stheta + e1[k]*dotProde13*(1-ctheta);

    // e4 should already be normalized but what the hell
    normalizeVector2norm( e4, 3 );

    // 2. Compute fiber angles and fiber coordinate basis vectors
    w = umu0;
    ut = PI - nu0;
    ut0 = PI - prm.nuTrans;
    G = (1 + cos((ut-ut0)*PI/ut0))/2*(ut < ut0) + (ut > ut0);
    psiTemp = prm.psi_in_b0 + (prm.psi_out_b0 - prm.psi_in_b0)*w;

    signPsi = psiTemp/fabs(psiTemp);
    if ( psiTemp == 0 )
    	signPsi = 1;

    psi = psiTemp + (1-G)*(PI/2*signPsi - psiTemp);
    cpsi = cos(psi);
    spsi = sin(psi);

    dotProde14 = dotProd(e1,e4,3);
    crossProduct3( e1, e4, crossProde14 );

    for(int k = 0; k < 3; k ++)
    {
    	es[k] = e1[k];
    	ef[k] = e4[k]*cpsi + crossProde14[k]*spsi + e1[k]*dotProde14*(1 - cpsi);
    }

    crossProduct3( ef, e1, en );


    // clean up
	delete[] crossProd;
	delete[] crossProde13;
	delete[] crossProde14;

	delete[] emu;
	delete[] enu;
	delete[] ephi;

}

















