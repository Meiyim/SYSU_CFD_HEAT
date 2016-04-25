#include <iostream>
#include <fstream>
#include <math.h>
#include "navier.h"
extern "C"{
#include "laspack/laspack.h"
}

using namespace std;

void HeadConductionSolver::UpdateEnergy( )
{
	int i,Iter;
	double *ESource,*ApE,IterRes,coef;

	Q_Constr(&As,   "matrixU",   Ncel, False, Rowws, Normal, True);
	ESource= new double[Ncel];
	ApE    = new double[Ncel];
	vec_init( ESource, Ncel, 0. );
	vec_init( ApE,     Ncel, 0. );
	// prepare the diffusion coefficient and source terms

	// part of viscous terms
	for(i=0;i!=Ncel;++i){
			ESource[i] = 0.0;	//set source
	}

	// boundary
	SetBCTemperature( BTem );
	// source terms, e.g., energy release, condensation/vaporization

	if( !IfSteady ){
	if(      TimeScheme==1 ){  // Euler forwards
		for( i=0; i<Ncel; i++ ){
			coef       = Rn[i]/dt * Cell[i].vol;
			ApE[i]    += coef;
			ESource[i]+= coef * Tnp[i];
		}
	}
	else if( TimeScheme==2 ){  // 2nd order BDF
		for( i=0; i<Ncel; i++ ){
			coef       = Rn[i]/dt * Cell[i].vol;
			ApE[i]    += 1.5*coef;
			ESource[i]+= coef * (2*Tnp[i]-0.5*Tnp2[i]);
		}
	}
	}

	// build matrix
	BuildScalarMatrix( 1, Tn,BTem,diffCoef,ESource,ApE );
	// Solve equations
	for( i=0; i<Ncel; i++ )
		xsol.Cmp[i+1]= Tn[i];
	SolveLinearEqu( CGIter, &As, &xsol, &bs, 500, SSORPrecond, 1.3, 1.e-8, &Iter, &IterRes );
	if( Iter>=500 && IterRes>1.e-8 ){
		cout<<"Energy cannot converge."<<Iter<<" "<<IterRes<<endl;
		exit(0);
	}
	for( i=0; i<Ncel; i++ )
		Tn[i] = xsol.Cmp[i+1];

	// clipping work
	delete [] ESource;
	delete [] ApE;

	Q_Destr ( &As );
}



void HeadConductionSolver::BuildScalarMatrix( int iSca, double *Phi,double *BPhi,double *DiffCoef, double *source, double *App )
{
	int i,j,iface, ip,in,ani[6],nj,rid,bnd;
	double app,apn[6],lambda,lambda2, Visc,dxc[3],
		dphidx,dphidy,dphidz, f,
		sav1,sav2,sav3,ViscAreaLen, 
		fcs=0., fde,fdi, dx[3],  sphi, pfi,pfh;

	Gradient ( Phi, BPhi,  dPhidX );

	for( i=0; i<Ncel; i++ )
	{

		app = App[i];
		nj  = 0 ;
		sphi= source[i];
		for( j=0;j<6;j++ ) apn[j] = 0.;

		for( j=0;j<Cell[i].nface;j++ )
		{
			iface  = Cell[i].face[j];
			ip     = Face[iface].cell1;
			in     = Face[iface].cell2;
			
			if( in<0 ) // boundary, i=ip naturally
			{
				bnd = Face[iface].bnd;
				rid = Bnd[bnd].rid;
				sav1    = Face[iface].n[0];
				sav2    = Face[iface].n[1];
				sav3    = Face[iface].n[2];
				// diffusion boundary
				Visc   = DiffCoef[i]; // This Should be changed using boundary condition. e.g., BVisTur[bnd]
				dphidx = dPhidX[i][0];
				dphidy = dPhidX[i][1];
				dphidz = dPhidX[i][2];
				vec_minus( dxc, Face[iface].x, Cell[i].x, 3 );
				switch( regionMap[rid].type1 ){
				case(1):
					if(regionMap[rid].type2==0){
						ViscAreaLen = Visc*Face[iface].rlencos;
						app += ViscAreaLen;
						fde  = Visc*( dphidx*sav1 + dphidy*sav2 + dphidz*sav3 );
						fdi  = ViscAreaLen*( dphidx*dxc[0] + dphidy*dxc[1] + dphidz*dxc[2] - BPhi[bnd] );
						break;
					}				
				case(4): //symmetric or given flux
					sphi += regionMap[rid].fixedValue * Face[iface].area;
					fde = 0.;
					fdi = 0.;
					break;
					// diffusion to implicit, only to central cell
				default:
					ErrorStop("no such rid");
				}
			}
			else // inner cell
			{
				// force i as face-left-cell, in as face-right-cell
				if( i != ip ){
					in  = ip;
					sav1    = -Face[iface].n[0];
					sav2    = -Face[iface].n[1];
					sav3    = -Face[iface].n[2];
					lambda  = 1.- Face[iface].lambda;
				}
				else{
					sav1    = Face[iface].n[0];
					sav2    = Face[iface].n[1];
					sav3    = Face[iface].n[2];
					lambda  = Face[iface].lambda;
				}
				
				lambda2= 1.- lambda;
				Visc   = lambda*DiffCoef[i]  + lambda2*DiffCoef[in];
				dphidx = lambda*dPhidX[i][0] + lambda2*dPhidX[in][0];
				dphidy = lambda*dPhidX[i][1] + lambda2*dPhidX[in][1];
				dphidz = lambda*dPhidX[i][2] + lambda2*dPhidX[in][2];


				// diffusion to implicit
				ViscAreaLen =  Visc*Face[iface].rlencos;
				app        +=  ViscAreaLen;
				apn[nj]    += -ViscAreaLen;
				ani[nj]     =  in;
				nj ++ ;

				// diffusion to source term. ( compressible or incompressible )
				fde = Visc*( dphidx*sav1 + dphidy*sav2 + dphidz*sav3 );
				vec_minus( dxc, Cell[in].x, Cell[i].x, 3 );
				fdi = ViscAreaLen*( dphidx*dxc[0]+dphidy*dxc[1]+dphidz*dxc[2] );
			}
			sphi += fde - fdi ;
		}

		// central cell coef is stored for later use
		// app   += Rn[i]/dt*Cell[i].vol;
		app   /= URF;

		Q_SetLen  (&As, i+1, nj + 1);
		// center cell
        		Q_SetEntry(&As, i+1, 0,  i+1, app);
		// off-diagonal
       	 	for( j=0; j<nj; j++ )
			Q_SetEntry(&As, i+1, j+1, ani[j]+1, apn[j]);

		// right hand side, including
		//   pressure, gravity and part of diffusion terms (explicit - implicit), 
		//   relaxation terms
		bs.Cmp[i+1] = sphi + (1.-URF)*app*Phi[i];
	}

}
