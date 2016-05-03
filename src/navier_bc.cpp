#include <math.h>
#include <iostream>
#include <fstream>
#include "navier.h"

using namespace std;

void NavierStokesSolver::SetBCVelocity( double *br, double *bu,double *bv,double *bw )
{
	int    i,rid,iface,ic;
	double unormal,sav1n,sav2n,sav3n;
	double* initvalues = NULL;

	for( i=0; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;
		switch( regionMap[rid].type1 ){
		case(1):  // wall
			br[i]= Rn[ic];
			bu[i]= 0.;
			bv[i]= 0.;
			bw[i]= 0.;
			break;
		case(2):  // inlet
			initvalues = regionMap[rid].initvalues;
			br[i]= initvalues[4];
			bu[i]= initvalues[0];
			bv[i]= initvalues[1];
			bw[i]= initvalues[2];

			// if( DensityModel==1 ) br[i]= PressureReference/(Rcpcv*298.);
			break;
		case(3):  // outlet
			br[i]= Rn[ic];
			bu[i]= Un[ic];
			bv[i]= Vn[ic];
			bw[i]= Wn[ic];
			break;
		case(4):  // symmetric
			sav1n= Face[iface].n[0]/Face[iface].area;
			sav2n= Face[iface].n[1]/Face[iface].area;
			sav3n= Face[iface].n[2]/Face[iface].area;
			unormal = Un[ic]*sav1n + Vn[ic]*sav2n + Wn[ic]*sav3n;
			bu[i]= Un[ic] - unormal * sav1n;
			bv[i]= Vn[ic] - unormal * sav2n;
			bw[i]= Wn[ic] - unormal * sav3n;
			br[i]= Rn[ic];
			break;

		default:
			errorHandler->fatalRuntimeError("no such boundary type in velocity bc:",regionMap[rid].type1);
		}
	}

	// correct outlet boundary condition
	if( DensityModel==0 ){
	double massflowin=0., massflowout=0., rate, areaout=0.;
	for( i=0; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		if( regionMap[rid].type1==2 )
			massflowin  += br[i]*( bu[i]*Face[iface].n[0] +
					       bv[i]*Face[iface].n[1] +
			     		       bw[i]*Face[iface].n[2] );
		else if( regionMap[rid].type1==3 ){
			massflowout += br[i]*( bu[i]*Face[iface].n[0] +
					       bv[i]*Face[iface].n[1] +
					       bw[i]*Face[iface].n[2] );
			areaout += br[i]*Face[i].area;
		}
	}
	if( fabs(massflowout)>SMALL ){
		rate = - massflowin / massflowout;
		for( i=0; i<Nbnd; i++ )
		{
			rid   = Bnd[i].rid;
			iface = Bnd[i].face;
			ic    = Face[iface].cell1;
			if( regionMap[rid].type1==3 )
			{
				bu[i] *= rate ;
				bv[i] *= rate ;
				bw[i] *= rate ;

			}
		}
	}
	else
	{
		rate = (-massflowin - massflowout)/areaout;
		for( i=0; i<Nbnd; i++ )
		{
			rid   = Bnd[i].rid;
			iface = Bnd[i].face;
			ic    = Face[iface].cell1;
			if( regionMap[rid].type1==3 )
			{
				bu[i] += rate * Face[iface].n[0]/Face[iface].area;
				bv[i] += rate * Face[iface].n[1]/Face[iface].area;
				bw[i] += rate * Face[iface].n[2]/Face[iface].area;
			}
		}
	}
	}
}

void NavierStokesSolver::SetBCPressure(double*bp)
{
	int    i,rid,iface,ic;
	for( i=0; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;
		switch( regionMap[rid].type1 ){
		case(1):  // wall
			bp[i] = Pn[ic];
			break;
		case(2):  // inlet
		/*
			bp[i] = pin;
			break;
		*/
		case(3):  // outlet, back step
		/*
			bp[i] = pout;
			break;
		*/
		case(4):
			bp[i] = Pn[ic];
			break;
		default:
			cout<<"no such boundary type"<<i<<" "<<rid<<endl;
			exit(0);
		}
	}
}

void NavierStokesSolver::SetBCDeltaP(double*bp, double *dp)
{
	int    i,rid,iface,ic;
	for( i=0; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;
		switch( regionMap[rid].type1 ){
		case(1):  // wall
		case(2):  // inlet
		case(3):
		case(4):
			bp[i] = dp[ic];
			break;
		default:
			cout<<"no such boundary type"<<i<<" "<<rid<<endl;
			exit(0);
		}
	}
}

void NavierStokesSolver::SetBCTemperature( double *bt, double* diffCoef)
{
	/*
	if(Solve3DHeatConduction){
		HeatConductionSolver* ht = dynamic_cast<HeatConductionSolver*> (physicalModule["3dHeatConduction"]);
	}	
	*/
	//update btem
	//ht->coupledBoundCommunicationSolid2Fluid(Bnd,NCoupledBnd);//update coupled boundary flux
	int    i,rid,iface,ic;
	for( i=0; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;
		BdRegion& reg = regionMap[rid];
		switch( reg.type1 ){
		case(1):  // wall
			if(reg.type2==0){//given T
				//remain initial value;
			}else if(reg.type2==1){ //given flux
     		   	bt[i] = Tn[ic] - Bnd[i].q / (diffCoef[ic] *Face[iface].rlencos );// given flux // by CHENXUYI
     		   	//explicitly changed 2kind_bnd to 1kind_bnd
				//bt[i] is to cal gradient
				//flux remain initial
			}else if(reg.type2==2){//coupled  // deel with 2nd_boundary in fluid field
     		   	bt[i] = Tn[ic] - Bnd[i].q / (diffCoef[ic] *Face[iface].rlencos );// given flux // by CHENXUYI
				//update boudary flux
			}else{
				assert(false);	
			}
			break;
		case(2):  // inlet
			//remain initual value;
			break;
		case(3):
			bt[i]= Tn[ic];
			break;
		case(4):
			bt[i]= Tn[ic];
			break;
		default:
			cout<<"no such boundary type"<<i<<" "<<rid<<endl;
			exit(0);
		}
	}
}

void NavierStokesSolver::SetBCSpecies ( double **brs )
{
}

void NavierStokesSolver::SetBCKEpsilon(double *TESource,double *EDSource,double *ApTE,double *ApED, double *Prod )
{
	int    i,bndType,iface,ic;
	double Cmu25,vnor,vt1,vt2,vt3,vel,utau,yplus,tauw, eps,sav1n,sav2n,sav3n,vol;
	using namespace TurKEpsilonVar;



	for( i=0; i<Nbnd; i++ )
	{
		bndType   = regionMap[ Bnd[i].rid ].type1;

		iface = Bnd[i].face;
		ic    = Face[iface].cell1;

		if(bndType==1){
		 	 // wall
			BTE[i]= TE[ic];
			BED[i]= ED[ic];

			vol = Cell[ic].vol;
			// calculate the skin friction velocity
			sav1n = Face[iface].n[0]/(Face[iface].area+1.e-16);
			sav2n = Face[iface].n[1]/(Face[iface].area+1.e-16);
			sav3n = Face[iface].n[2]/(Face[iface].area+1.e-16);
			vnor= Un[ic]*sav1n + Vn[ic]*sav2n + Wn[ic]*sav3n;
			vt1 = Un[ic] - vnor*sav1n;
			vt2 = Vn[ic] - vnor*sav2n;
			vt3 = Wn[ic] - vnor*sav3n;
			vel = sqrt( vt1*vt1 + vt2*vt2 + vt3*vt3 );

			Cmu25 = pow(Cmu,0.25);
			utau  = Cmu25*sqrt( TE[ic] );
			yplus = Bnd[i].distance * utau * Rn[ic] / (VisLam[ic]+SMALL);
			tauw  = VisLam[ic] * vel/Bnd[i].distance;  // ?? which viscosity ??
			if( yplus<11. )  // inner
				utau = sqrt(tauw/Rn[ic]);
			else             // outer
				utau = Cmu25 * sqrt(TE[ic]);
			Prod[ic]= tauw * utau / (kappa*Bnd[i].distance);
			eps     = pow(Cmu,0.75)*pow(TE[ic],1.5)/(kappa*Bnd[i].distance);
			TESource[ic]  = Prod[ic] * vol ;
			ApTE    [ic] -= Rn[ic] * ED[ic] /TE[ic] * vol;
			ApTE    [ic] += Rn[ic] * eps    /TE[ic] * vol;

			ED      [ic] = eps;
			BED     [i ] = eps;
		}
	}

	for( i=0; i<Nbnd; i++ )
	{
		bndType   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;

		if(bndType==2){
			double *initvalues = regionMap[Bnd[i].rid].initvalues;
			BTE[i]= initvalues[6] ;  // turbulence intensity
			BED[i]= initvalues[7] ;
		}
	}
	
	for( i=0; i<Nbnd; i++ )
	{
		bndType   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;

		if(bndType==3){
			BTE[i]= TE[ic];
			BED[i]= ED[ic];
		}
	}

	for( i=0; i<Nbnd; i++ )
	{
		bndType   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;

		if(bndType==4){
			BTE[i]= TE[ic];
			BED[i]= ED[ic];
		}
	}
	
}
