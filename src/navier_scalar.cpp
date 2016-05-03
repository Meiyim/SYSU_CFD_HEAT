#include <iostream>
#include <fstream>
#include <math.h>
#include "navier.h"
extern "C"{
#include "laspack/laspack.h"
}

using namespace std;

//-----------------------------------
// k-epsilon turbulence model
//------------------------------
void NavierStokesSolver::UpdateTurKEpsilon( )
{
	int i,Iter, rid,ic,iface;
	double *dudx,*dvdx,*dwdx,
		   s1,s2,s3,Dis, vol,coef,
		   fact, IterRes,
		   *Prod,*TESource,*EDSource,*VisTE,*VisED,*ApTE,*ApED,Pked,Dised;
	ofstream of;
	using namespace TurKEpsilonVar;

	Prod     = new double[Ncel];
	TESource = new double[Ncel];
	EDSource = new double[Ncel];
	VisTE    = new double[Ncel];
	VisED    = new double[Ncel];
	ApTE     = new double[Ncel];
	ApED     = new double[Ncel];
	vec_init( TESource,Ncel,0. );
	vec_init( EDSource,Ncel,0. );
	vec_init( ApTE,    Ncel,0. );
	vec_init( ApED,    Ncel,0. );
	// viscosity and source term
	for( i=0; i<Ncel; i++ )
	{
		vol = Cell[i].vol;
		// turbulent kinetic viscosity
        VisTE[i] = VisLam[i] + VisTur[i] / Sigma_k ;
        // turbulent dissipation viscosity
        VisED[i] = VisLam[i] + VisTur[i] / Sigma_e ;

		// source terms
		dudx = dUdX[i];
		dvdx = dVdX[i];
		dwdx = dWdX[i];
        // rate of production of turbulent energy (eq. 9.40)
        s1 = (dudx[0]+dudx[0])*dudx[0] + (dudx[1]+dvdx[0])*dudx[1] + (dudx[2]+dwdx[0])*dudx[2];
        s2 = (dvdx[0]+dudx[1])*dvdx[0] + (dvdx[1]+dvdx[1])*dvdx[1] + (dvdx[2]+dwdx[1])*dvdx[2];
        s3 = (dwdx[0]+dudx[2])*dwdx[0] + (dwdx[1]+dvdx[2])*dwdx[1] + (dwdx[2]+dwdx[2])*dwdx[2];

        Prod[i] = VisTur[i] * ( s1 + s2 + s3 ) ;
        // dissipation
        Dis     = Rn[i] * ED[i];
		// bouyancy production term: - Gi/(sigma_h,t rho) drho/dx
        // Pbouy =
		TESource[i] = Prod[i] * vol ;   //   - Dis  ;
		ApTE    [i] = Dis  /TE[i] * vol;        // add to diagonal element of linear system As

		//// ED production & dissipation. similarity to TE
		// fact = ED[i]/(TE[i]+SMALL);
		// Pked = Ceps1 * fact * Prod[i] * vol ;
		// Dised= Ceps2 * fact * Dis     * vol ;
		// EDSource[i] = Pked * vol;  //   - Dised;
		// ApED    [i] = Dised/ED[i] * vol;
	}
	// boundary for k, epsilon
	SetBCKEpsilon( TESource,EDSource,ApTE,ApED,Prod );
	// other source term. do not put the before the boundary, for boundary change Pk,epsilon mandatorily

	if( !IfSteady ){
	if(      TimeScheme==1 ){  // Euler forwards
		for( i=0; i<Ncel; i++ ){
			coef        = Rn[i]/dt*Cell[i].vol;
			ApTE[i]    += coef;
			TESource[i]+= coef * TEp[i];
			ApED[i]    += coef;
			EDSource[i]+= coef * EDp[i];
		}
	}
	else if( TimeScheme==2 ){  // 2nd order BDF
		for( i=0; i<Ncel; i++ ){
			coef        = Rn[i]/dt*Cell[i].vol;
			ApTE[i]    += 1.5*coef;
			TESource[i]+= coef * (2*TEp[i]-0.5*TEp2[i]);
			ApED[i]    += 1.5*coef;
			EDSource[i]+= coef * (2*EDp[i]-0.5*EDp2[i]);
		}
	}
	}


	//Solve turbulence kinetic energy
	Q_Constr(&As,   "matrixU",   Ncel, False, Rowws, Normal, True);
	// build matrix
	BuildScalarMatrix( TE,BTE,VisTE,TESource, ApTE ,NULL);
	// Solve equations
	for( i=0; i<Ncel; i++ )
		xsol.Cmp[i+1]= TE[i];
	SolveLinearEqu( GMRESIter, &As, &xsol, &bs, 500, SSORPrecond, 1.3, 1.e-8, &Iter, &IterRes );
	if( Iter>=500 && IterRes>1.e-8 ){
		cout<<"Energy cannot converge."<<Iter<<" "<<IterRes<<endl;
		exit(0);
	}
	for( i=0; i<Ncel; i++ ){
		TE[i] = CYCASMAX(1.e-9,xsol.Cmp[i+1]);
		/* if( TE[i]<0. ) {
			cout<<"error in te,"<<TE[i]<<endl;
			exit(0);
		} */
	}
	Q_Destr ( &As );


	for( i=0; i<Ncel; i++ )
	{
		// ED production & dissipation. similarity to TE
		fact = ED[i]/(TE[i]+SMALL)*Cell[i].vol;
		EDSource[i]  = Ceps1* fact * Prod[i];
		ApED    [i] += Ceps2* fact * Rn[i]  ;
	}


	// Solve turbulence dissipation rate
	Q_Constr(&As,   "matrixU",   Ncel, False, Rowws, Normal, True);
	// build matrix
	BuildScalarMatrix(ED,BED,VisED,EDSource,ApED,NULL );
	// ?? specieal treatment, force ED[i] on boundary cells to be BED[ib]  ??
	for( i=0; i<Nbnd; i++ ){
		rid   = Bnd[i].rid ;
		if( regionMap[rid].type1==1 ){
			iface = Bnd[i].face;
			ic    = Face[iface].cell1;
			Q_SetLen  (&As, ic+1, 1);
			Q_SetEntry(&As, ic+1, 0,  ic+1, 1.);
			bs.Cmp[ic+1] = ED[ic];
		}
	}
	// Solve equations
	for( i=0; i<Ncel; i++ )
		xsol.Cmp[i+1]= ED[i];
	SolveLinearEqu( GMRESIter, &As, &xsol, &bs, 500, SSORPrecond, 1.3, 1.e-8, &Iter, &IterRes );
	if( Iter>=500 && IterRes>1.e-8 ){
		cout<<"Energy cannot converge."<<Iter<<" "<<IterRes<<endl;
		exit(0);
	}
	for( i=0; i<Ncel; i++ ){
		ED[i] = CYCASMAX(1.e-12,xsol.Cmp[i+1]);
	}
	Q_Destr ( &As );
	

	// Calculate the turbulent viscosity
	for( i=0; i<Ncel; i++ )
    {
		if( ED[i]>1.e-12 )
			VisTur[i]= Cmu * Rn[i] * TE[i]*TE[i]/ (ED[i]+SMALL);
		else
			VisTur[i]= 0.;
    }

	delete []Prod;
	delete []TESource;
	delete []EDSource;
	delete []VisTE;
	delete []VisED;
	delete []ApTE;
	delete []ApED;


}

// a callback function to implement given flux bc for tempreture
// one shuold notice that gradient should solve before calling this function
void specialTreatmentForWallTempreture(NavierStokesSolver* solver, 
	  const int icell,const int jface, 
	  const double* DiffCoef,const double* BPhi,
	  double& app, double& fde, double& fdi, double& fcs ){
	int bnd = solver->Face[jface].bnd;
	int rid = solver->Bnd[bnd].rid;
	double sav1    = solver->Face[jface].n[0];
	double sav2    = solver->Face[jface].n[1];
	double sav3    = solver->Face[jface].n[2];
	double RUnormal= solver->RUFace[jface];

	// diffusion boundary
	double Visc   = DiffCoef[icell]; // This Should be changed using boundary condition. e.g., BVisTur[bnd]
	double dphidx = solver->dPhidX[icell][0];
	double dphidy = solver->dPhidX[icell][1];
	double dphidz = solver->dPhidX[icell][2];

	double dxc[3];
	vec_minus( dxc, solver->Face[jface].x, solver->Cell[icell].x, 3 );
	if(solver->regionMap[rid].type1 == 1 && solver->regionMap[rid].type2 == 1){//given flux;
		fde = solver->Bnd[bnd].q * solver->Face[jface].area;	
		fdi = 0.0;
	}else if(solver->regionMap[rid].type1 == 1 && solver->regionMap[rid].type2 == 2){//coupled the same as above
		fde = solver->Bnd[bnd].q * solver->Face[jface].area;
		fdi = 0.0;
	}else if(solver->regionMap[rid].type1==4){//symmetric heat
		fde = fdi = 0.0;
	}else{
		// diffusion to implicit, only to central cell
		double ViscAreaLen = Visc*solver->Face[jface].rlencos;
		app += ViscAreaLen;
		fde  = Visc*( dphidx*sav1 + dphidy*sav2 + dphidz*sav3 );
		fdi  = ViscAreaLen*( dphidx*dxc[0] + dphidy*dxc[1] + dphidz*dxc[2] - BPhi[bnd] );

	}

	// convection boundary
	double f=0;
	switch( solver->regionMap[rid].type1 ){
	case(1):     //---- Wall ----
				// convection to implicit, nothing
		fcs = 0.;
		break;
	case(2):     //---- Inlet ----
		// convection to implicit
		if( RUnormal>0 ){
			cout<<"reverse flow get out of inlet. stop!"<<endl;
			exit(0);
		}
		f    = CYCASMIN( RUnormal , 0.0 );
		app -= f;
		fcs  = f*BPhi[bnd];
		break;
	case(3):     //---- Outlet ----??????????? boundary equal inner cell, so ???
		if( RUnormal<0 ){
			cout<<"reverse flow get in of outlet. stop!"<<endl;
			// exit(0);
		}
		f    = CYCASMIN( RUnormal , 0.0 );
		app -= f;
		fcs  = f*BPhi[bnd];
		break;
	case(4):     //---- Symmetric ----
		// RUnormal = 0.
		fcs = 0.;
		break;
	default:
		cout<<"no this type of boundary! rid="<<rid<<endl;
		exit(0);
	}
}



//------------------------------------------------------------
// Energy transport equation : 
// p_(ro*T)/p_t + div(ro*V*T) = div(romda/cp*grad(T)) + QT/cp
//  note that cp is removed from time derivative
//-------------------------------------------------------
void NavierStokesSolver::UpdateEnergy( )
{
	int i,Iter;
	double *kcond,*ESource,*ApE,IterRes,coef;

	Q_Constr(&As,   "matrixU",   Ncel, False, Rowws, Normal, True);
	kcond  = new double[Ncel];
	ESource= new double[Ncel];
	ApE    = new double[Ncel];
	vec_init( ESource, Ncel, 0. );
	vec_init( ApE,     Ncel, 0. );
	// prepare the diffusion coefficient and source terms
	for( i=0; i<Ncel; i++ )
	{
		kcond[i] = cp*(VisLam[i]/prl + VisTur[i]/prte)/cp;
		kcond[i] *= 100;//debug: kcond seems too small

	}

	// part of viscous terms
	if( DensityModel==1 ){
	int c1,c2;
	double VisL,VisT,vmul,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,
		vlac,div,txx,tyy,tzz,txy,txz,tyz, vxg,vyg,vzg,btx,bty,btz,fvis,lambda,lambda2;
	for( i=0; i<Nfac; i++ )
	{
		c1 = Face[i].cell1;
		c2 = Face[i].cell2;
		if( c2<0 )
		{
			VisL = VisLam[c1];
			VisT = VisTur[c1];
			dudx = dUdX[c1][0];
			dudy = dUdX[c1][1];
			dudz = dUdX[c1][2];
			dvdx = dVdX[c1][0];
			dvdy = dVdX[c1][1];
			dvdz = dVdX[c1][2];
			dwdx = dWdX[c1][0];
			dwdy = dWdX[c1][1];
			dwdz = dWdX[c1][2];
			vxg  = Un[c1];
			vyg  = Vn[c1];
			vzg  = Wn[c1];
		}
		else
		{
			lambda = Face[i].lambda;
			lambda2= 1.-lambda;
			VisL = lambda*VisLam[c1]  + lambda2*VisLam[c2]  ;
			VisT = lambda*VisTur[c1]  + lambda2*VisTur[c2]  ;
			dudx = lambda*dUdX[c1][0] + lambda2*dUdX[c2][0] ;
			dudy = lambda*dUdX[c1][1] + lambda2*dUdX[c2][1] ;
			dudz = lambda*dUdX[c1][2] + lambda2*dUdX[c2][2] ;
			dvdx = lambda*dVdX[c1][0] + lambda2*dVdX[c2][0] ;
			dvdy = lambda*dVdX[c1][1] + lambda2*dVdX[c2][1] ;
			dvdz = lambda*dVdX[c1][2] + lambda2*dVdX[c2][2] ;
			dwdx = lambda*dWdX[c1][0] + lambda2*dWdX[c2][0] ;
			dwdy = lambda*dWdX[c1][1] + lambda2*dWdX[c2][1] ;
			dwdz = lambda*dWdX[c1][2] + lambda2*dWdX[c2][2] ;
			vxg  = lambda*Un[c1]      + lambda2*Un[c2];
			vyg  = lambda*Vn[c1]      + lambda2*Vn[c2];
			vzg  = lambda*Wn[c1]      + lambda2*Wn[c2];
		}

		vmul= VisL + VisT;
		vlac=  -2./3.*vmul;
		div= dudx+ dvdy+ dwdz;
		txx= 2.*vmul*dudx+ vlac*div ; //CXY: pressure should goes here ?
		tyy= 2.*vmul*dvdy+ vlac*div ;
		tzz= 2.*vmul*dwdz+ vlac*div ;
		txy= vmul*(dudy+dvdx);
		txz= vmul*(dudz+dwdx);
		tyz= vmul*(dvdz+dwdy);
		btx= vxg*txx+ vyg*txy+ vzg*txz ;
		bty= vxg*txy+ vyg*tyy+ vzg*tyz ;
		btz= vxg*txz+ vyg*tyz+ vzg*tzz ;
		fvis= ( btx*Face[i].n[0] + bty*Face[i].n[1] + btz*Face[i].n[2] ) / cp;//CXY: where is work done by pressure?
		ESource[c1] += fvis;
		if( c2>=0 )
		ESource[c2] -= fvis; //CXY: dissipation source term due to compressibility
                             //CX: why is this integral done on face? i supposed it should be done in volumn
	}
	}

	// boundary
	SetBCTemperature( BTem,kcond );
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
	BuildScalarMatrix( Tn,BTem,kcond,ESource,ApE,specialTreatmentForWallTempreture);
	// Solve equations
	for( i=0; i<Ncel; i++ )
		xsol.Cmp[i+1]= Tn[i];
	SolveLinearEqu( GMRESIter, &As, &xsol, &bs, 500, SSORPrecond, 1.3, 1.e-8, &Iter, &IterRes );
	if( Iter>=500 && IterRes>1.e-8 ){
		cout<<"Energy cannot converge."<<Iter<<" "<<IterRes<<endl;
		exit(0);
	}
	for( i=0; i<Ncel; i++ )
		Tn[i] = xsol.Cmp[i+1];

	// clipping work
	delete [] kcond;
	delete [] ESource;
	delete [] ApE;

	Q_Destr ( &As );
}

//------------------------------------------------------------
// species transport equation : 
// p_(ro*c)/p_t + div(ro*V*c) = div(D*grad(c)) + Qm
//    c is mass fraction, Qm could be injection,phase transition etc
//--------------------------------------------
void NavierStokesSolver::UpdateSpecies( )
{
	int i,is,Iter;
	double **DiffC, **ScSource,*ApS, IterRes,coef;

	ApS     = new double[Ncel]; 
	DiffC   = new_Array2D<double>(Nspecies,Ncel);
	ScSource= new_Array2D<double>(Nspecies,Ncel);
	init_Array2D( ScSource,Nspecies,Ncel, 0. );
	// prepare the diffusion coefficient and source terms
	for( is=0; is<Nspecies; is++){
		for( i=0; i<Ncel; i++ )
		{
			DiffC[is][i] = 1.e-5; // this should be calculated in Material class
			if( ! IfSteady )
				ScSource[is][i] += 0.;
		}
	}
	// boundary
	SetBCSpecies( BRS );
	// source terms, e.g., energy release

	for( is=0; is<Nspecies; is++ )
	{
		vec_init( ApS, Ncel, 0. );

		// transient time to source term and diagonal element
		if( !IfSteady ){
			if(      TimeScheme==1 ){  // Euler forwards
				for( i=0; i<Ncel; i++ ){
					coef            = Rn[i]/dt*Cell[i].vol;
					ApS[i]         += coef;
					ScSource[is][i]+= coef * RSnp[is][i];
				}
			}
			else if( TimeScheme==2 ){  // 2nd order BDF
				for( i=0; i<Ncel; i++ ){
					coef            = Rn[i]/dt*Cell[i].vol;
					ApS[i]         += 1.5*coef;
					ScSource[is][i]+= coef * (2*RSnp[is][i]-0.5*RSnp2[is][i]);
				}
			}
		}

		Q_Constr(&As,   "matrixU",   Ncel, False, Rowws, Normal, True);
		// build matrix
		BuildScalarMatrix( /*4+is,*/ RSn[is],BRS[is],DiffC[is],ScSource[is],ApS,NULL);
		// Solve equations
		for( i=0; i<Ncel; i++ ) 
			xsol.Cmp[i+1]= RSn[is][i];
		SolveLinearEqu( GMRESIter, &As, &xsol, &bs, 500, SSORPrecond, 1.3, 1.e-8, &Iter, &IterRes );
		if( Iter>=500 && IterRes>1.e-8 ){
			cout<<"Energy cannot converge."<<Iter<<" "<<IterRes<<endl;
			exit(0);
		}
		for( i=0; i<Ncel; i++ )
			RSn[is][i] = xsol.Cmp[i+1];

		Q_Destr ( &As );
	}
}

void NavierStokesSolver::BuildScalarMatrix(
	double *Phi,double *BPhi,double *DiffCoef, 
	double *source, double *App,
	CallBackType specialTreatmentOnBoundary = NULL)
	
{
	int i,j,iface, ip,in,ani[6],nj,rid,bnd;
	double app,apn[6],lambda,lambda2, Visc,dxc[3],
		dphidx,dphidy,dphidz, f,
		sav1,sav2,sav3,RUnormal,ViscAreaLen, 
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
			
			if(in < 0 && specialTreatmentOnBoundary != NULL){//call call back function to treet boundary
				specialTreatmentOnBoundary(this,i,iface,DiffCoef,BPhi,app,fde,fdi,fcs);

			}else if( in<0 && specialTreatmentOnBoundary == NULL) // boundary, i=ip naturally
			{
				bnd = Face[iface].bnd;
				rid = Bnd[bnd].rid;
				sav1    = Face[iface].n[0];
				sav2    = Face[iface].n[1];
				sav3    = Face[iface].n[2];
				RUnormal= RUFace[iface];

				// diffusion boundary
				Visc   = DiffCoef[i]; // This Should be changed using boundary condition. e.g., BVisTur[bnd]
				dphidx = dPhidX[i][0];
				dphidy = dPhidX[i][1];
				dphidz = dPhidX[i][2];
				vec_minus( dxc, Face[iface].x, Cell[i].x, 3 );
				switch( regionMap[rid].type1 ){
				case(1):
				case(2):
				case(3):
				case(4):
					/*fde = 0.;
					fdi = 0.;
					break;*/
					// diffusion to implicit, only to central cell
					ViscAreaLen = Visc*Face[iface].rlencos;
					app += ViscAreaLen;
					fde  = Visc*( dphidx*sav1 + dphidy*sav2 + dphidz*sav3 );
					fdi  = ViscAreaLen*( dphidx*dxc[0] + dphidy*dxc[1] + dphidz*dxc[2] - BPhi[bnd] );
					break;
				default:
					errorHandler->fatalRuntimeError("no such rid",rid);
				}
				
				// convection boundary
				switch( regionMap[rid].type1 ){
				case(1):     //---- Wall ----
					// convection to implicit, nothing
					fcs = 0.;
					break;
				case(2):     //---- Inlet ----
					// convection to implicit
					if( RUnormal>0 ){
						cout<<"reverse flow get out of inlet. stop!"<<endl;
						exit(0);
					}
					f    = CYCASMIN( RUnormal , 0.0 );
					app -= f;
					fcs  = f*BPhi[bnd];
					break;
				case(3):     //---- Outlet ----??????????? boundary equal inner cell, so ???
					if( RUnormal<0 ){
						cout<<"reverse flow get in of outlet. stop!"<<endl;
						// exit(0);
					}
					f    = CYCASMIN( RUnormal , 0.0 );
					app -= f;
					fcs  = f*BPhi[bnd];
					break;
				case(4):     //---- Symmetric ----
					// RUnormal = 0.
					fcs = 0.;
					break;
				default:
					cout<<"no this type of boundary! rid="<<rid<<endl;
					exit(0);
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
					RUnormal= -RUFace[iface];
					lambda  = 1.- Face[iface].lambda;
				}
				else{
					sav1    = Face[iface].n[0];
					sav2    = Face[iface].n[1];
					sav3    = Face[iface].n[2];
					RUnormal= RUFace[iface];
					lambda  = Face[iface].lambda;
				}
				
				lambda2= 1.- lambda;
				Visc   = lambda*DiffCoef[i]  + lambda2*DiffCoef[in];
				dphidx = lambda*dPhidX[i][0] + lambda2*dPhidX[in][0];
				dphidy = lambda*dPhidX[i][1] + lambda2*dPhidX[in][1];
				dphidz = lambda*dPhidX[i][2] + lambda2*dPhidX[in][2];

				// convection to implicit
				if( RUnormal<0. ){
					apn[nj] += RUnormal;
					app     -= RUnormal;
				}

				// diffusion to implicit
				ViscAreaLen =  Visc*Face[iface].rlencos;
				app        +=  ViscAreaLen;
				apn[nj]    += -ViscAreaLen;
				ani[nj]     =  in;
				nj ++ ;

				// convection to source term. (high order schemes)
				if( RUnormal>0. )
				{
					// low order interpolation, 1st upwind
					pfi= Phi[i];
					// high order interpolation, e.g., 2nd upwind
					vec_minus( dx, Face[iface].x, Cell[i].x, 3);
					pfh= Phi[i] + vec_dot( dPhidX[i], dx, 3 );
				}
				else
				{
					// low order interpolation, 1st upwind
					pfi= Phi[in];
					// high order interpolation, e.g., 2nd upwind
					vec_minus( dx, Face[iface].x, Cell[in].x, 3);
					pfh= Phi[in] + vec_dot( dPhidX[in], dx, 3 );
				}
				fcs = RUnormal*(pfh-pfi);
				
				// diffusion to source term. ( compressible or incompressible )
				fde = Visc*( dphidx*sav1 + dphidy*sav2 + dphidz*sav3 );
				vec_minus( dxc, Cell[in].x, Cell[i].x, 3 );
				fdi = ViscAreaLen*( dphidx*dxc[0]+dphidy*dxc[1]+dphidz*dxc[2] );
			}
			sphi += fde - fdi - fcs;
		}

		// central cell coef is stored for later use
		// app   += Rn[i]/dt*Cell[i].vol;
		app   /= URF[5];  // relaxation

		Q_SetLen  (&As, i+1, nj + 1);
		// center cell
        Q_SetEntry(&As, i+1, 0,  i+1, app);
		// off-diagonal
        for( j=0; j<nj; j++ )
			Q_SetEntry(&As, i+1, j+1, ani[j]+1, apn[j]);

		// right hand side, including
		//   pressure, gravity and part of diffusion terms (explicit - implicit), 
		//   relaxation terms
		bs.Cmp[i+1] = sphi + (1.-URF[5])*app*Phi[i];
	}

}
