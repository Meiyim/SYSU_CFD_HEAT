#include <iostream>
#include "BasicType.h"
using namespace std;

// Gradient calculation using Gaussian law
CycasSolver::CycasSolver():
    IfSteady(true),
    dt(0.0),
    TimeScheme(1),
    limiter(0),

    Nvrt(0),Ncel(0),Nfac(0),Nbnd(0),
    NCoupledBnd(0),
    Vert(NULL),
    Face(NULL),
    Cell(NULL),
    Bnd(NULL)
{
}

CycasSolver::~CycasSolver(){
    delete []Vert;
    delete []Face;
    delete []Cell;
    delete []Bnd;
}

int CycasSolver::Gradient( double *phi, double *Bphif, double **phigd )
{
    int    i,g, c1,c2;
    double lambda,pf;
    // using Gauss theorem
    for( i=0; i<Ncel; i++ )
        for(g=0;g<3;g++)
            phigd[i][g]= 0.;

    for( i=0; i<Nfac; i++ )
    {
        lambda = Face[i].lambda;
        c1     = Face[i].cell1;
        c2     = Face[i].cell2;

        if( c2>=0 ){
                pf = lambda*phi[c1] + (1.-lambda)*phi[c2];
	for( g=0;g<3;g++ ){
		phigd[c1][g] += pf * Face[i].n[g];
		phigd[c2][g] -= pf * Face[i].n[g];
	}
       }else{
            pf = Bphif[Face[i].bnd]; // how to add boundary condition ?
	    for( g=0;g<3;g++ ){
		phigd[c1][g] += pf * Face[i].n[g];
                   }

            }
    }
    //CHECK_ARRAY(Bphif,Nbnd);
    //CHECK_ARRAY(phigd[0],3*Ncel);
    for( i=0; i<Ncel; i++ ){
        for( g=0; g<3; g++ )
            phigd[i][g] /= Cell[i].vol;
     }

    if(      limiter==0 )
    {}
    else if( limiter==1 )
	Limiter_Barth( phi, phigd );
    else if( limiter==2 )
	Limiter_MLP  ( phi, phigd );
    else if( limiter==3 )
	Limiter_WENO ( phi, phigd );
    else
                errorHandler->fatalRuntimeError("no such limiter choice");

    return 0;
}






/////////////////////////////////////////////////
/// two kinds of limiters : MLP, WENO
///------------------------------------------

inline double Barth_fun( double us,double umin, double umax)
{
    double ff;
    if( us>1.0e-10 )
      ff= CYCASMIN(1., umax/us );
    else if( us< -1.0e-10 )
      ff= CYCASMIN(1., umin/us );
    else
      ff= 1.;
    return ff;
}
inline double smoothfun( double us,double umin,double umax, double epsilon2)
{
	double ff,umax2,umin2,us2;
    if( us>1.0e-10){
		umax2= umax * umax;
		us2  = us   * us;
		ff= 1./(CYCASSIGN(us)*(fabs(us)+1.e-12)) *
            ( (umax2+epsilon2)*us+ 2.*us2*umax )/( umax2+2.*us2+us*umax+ epsilon2);
	}
    else if( us< -1.0e-10 ){
		umin2= umin * umin;
		us2  = us   * us  ;
		ff= 1./(CYCASSIGN(us)*(fabs(us)+1.e-12)) *
            ( (umin2+epsilon2)*us+ 2.*us2*umin )/( umin2+2.*us2+us*umin+ epsilon2);
	}
    else
      ff= 1.;
    return ff;
}
int CycasSolver::Limiter_MLP(double UC[],double **GradU)
{
    int    i,j,iv;
    double *ficell, *UMax,*UMin, fi,umin,umax,dx,dy,dz,us;
    ficell = new double[Ncel];
    UMax   = new double[Nvrt];
    UMin   = new double[Nvrt];

    // vertex max and min value
    for( i=0; i<Nvrt; i++ )
    {
        UMax[i]= -1.e8;
        UMin[i]=  1.e8;
    }
    for( i=0; i<Ncel; i++ )
    {
        for(j=0; j<8; j++ ) // this may be a little wasteful, optimize later
        {
            iv = Cell[i].vertices[j];
            UMax[iv]= CYCASMAX( UMax[iv], UC[i] );
            UMin[iv]= CYCASMIN( UMin[iv], UC[i] );
        }
    }

    // MLP : find slope limit coef; limit gradients
    for( i=0; i<Ncel; i++ )
    {
        fi = 10.;
        for(j=0; j<8; j++ ) // this may be a little wasteful, 8 changed to more 
        {
            iv = Cell[i].vertices[j];
            umax= UMax[iv] - UC[i];
            umin= UMin[iv] - UC[i];

            dx = Vert[iv][0] - Cell[i].x[0];
            dy = Vert[iv][1] - Cell[i].x[1];
            dz = Vert[iv][2] - Cell[i].x[2];
            us = GradU[i][0]*dx + GradU[i][1]*dy + GradU[i][2]*dz ;
			// Barth
			fi = CYCASMIN( fi, Barth_fun(us,umin,umax) );
			// Venkatanishnan
			// fi = CYCASMIN( fi, smoothfun( us,umin,umax, 1./ Cell[i].vol**3 ) );
        }
        ficell[i]= fi;
    }


	//CHECK_ARRAY(UC,Ncel);
	//CHECK_ARRAY(ficell,Ncel);
	//CHECK_ARRAY(GradU[0],3*Ncel);
    for( i=0; i<Ncel; i++ )
    {
        GradU[i][0] = ficell[i] * GradU[i][0];
        GradU[i][1] = ficell[i] * GradU[i][1];
        GradU[i][2] = ficell[i] * GradU[i][2];
    }

	//CHECK_ARRAY(GradU[0],3*Ncel);
    delete [] ficell;
    delete [] UMax;
    delete [] UMin;
    return 1;
}

int CycasSolver::Limiter_Barth(double UC[],double **GradU)
{
    int    i,j,iv,in;
    double *ficell, fi,umin,umax,dx,dy,dz,us;
    ficell = new double[Ncel];

    // MLP : find slope limit coef; limit gradients
    for( i=0; i<Ncel; i++ )
    {
        fi = 10.;
		umax= UC[i];
        umin= UC[i];
		for(j=0; j<Cell[i].nface; j++ ) // this may be a little wasteful, optimize later
        {
			in = Cell[i].cell[j];
			if( in<0 )continue;
            umax= CYCASMAX( umax, UC[in] );
            umin= CYCASMIN( umin, UC[in] );
        }
		umax -= UC[i];
        umin -= UC[i];
		for(j=0; j<8; j++ ) // this may be a little wasteful, 8 changed to more 
		{
            iv = Cell[i].vertices[j];
            dx = Vert[iv][0] - Cell[i].x[0];
            dy = Vert[iv][1] - Cell[i].x[1];
            dz = Vert[iv][2] - Cell[i].x[2];
            us = GradU[i][0]*dx + GradU[i][1]*dy + GradU[i][2]*dz ;
			// Barth
			fi = CYCASMIN( fi, Barth_fun(us,umin,umax) );
        }
        ficell[i]= fi;
    }
    for( i=0; i<Ncel; i++ )
    {
        GradU[i][0] = ficell[i] * GradU[i][0];
        GradU[i][1] = ficell[i] * GradU[i][1];
        GradU[i][2] = ficell[i] * GradU[i][2];
    }

    delete [] ficell;
    return 1;
}

int CycasSolver::Limiter_WENO(double *UC,double **GradU)
{
    int i,j, ic,iface;
    double ss,wei, **gdTmp;

    gdTmp = new_Array2D<double>( Ncel, 3 );
	for( i=0; i<Ncel; i++ )
	{
		gdTmp[i][0]= 0.;
		gdTmp[i][1]= 0.;
		gdTmp[i][2]= 0.;
	}

    // WENO Limiter
    for( i=0; i<Ncel; i++ )
    {
        // non-linear weights
        ss = 0.;
        for( j=0; j<Cell[i].nface; j++ )
        {
            iface = Cell[i].face[j];
            ic    = Cell[i].cell[j];
            if( ic<0 ) continue;
            wei = 1./( GradU[ic][0]*GradU[ic][0] +
                       GradU[ic][1]*GradU[ic][1] +
                       GradU[ic][2]*GradU[ic][2] +1.e-16 );
            ss += wei;
			gdTmp[i][0] += wei*GradU[ic][0];
			gdTmp[i][1] += wei*GradU[ic][1];
			gdTmp[i][2] += wei*GradU[ic][2];
        }

		//// boundary cells, just quit
		//if(j<Cell[i].nface)
		//{
		//	gdTmp[i][0] = GradU[i][0];
		//	gdTmp[i][1] = GradU[i][1];
		//	gdTmp[i][2] = GradU[i][2];
		//	continue;
		//}

        wei     = 1./( GradU[i ][0]*GradU[i ][0] +
                       GradU[i ][1]*GradU[i ][1] +
                       GradU[i ][2]*GradU[i ][2] +1.e-16 );
        ss += wei;
		gdTmp[i][0] += wei*GradU[i][0];
		gdTmp[i][1] += wei*GradU[i][1];
		gdTmp[i][2] += wei*GradU[i][2];
		gdTmp[i][0] /= ss;
		gdTmp[i][1] /= ss;
		gdTmp[i][2] /= ss;
    }
    // substitution
    for( i=0; i<Ncel; i++ )
    {
        GradU[i][0]= gdTmp[i][0];
        GradU[i][1]= gdTmp[i][1];
        GradU[i][2]= gdTmp[i][2];
    }

    delete_Array2D( gdTmp, Ncel, 3 );
    return 1;
}
