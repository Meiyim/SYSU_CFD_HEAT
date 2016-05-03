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
    //CHECK_ARRAY(phi,Ncel);
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


//---------------GEOMETRY CONFIGURATION-------------//
/////////////////////////////////////////////////////
/// Other cell and face information
///--------------------------------------------------

double Tri3DArea(double x1[], double x2[], double x3[] )
{
    double a1,a2,a3, v1[3],v2[3];
    vec_minus(v1, x2,x1, 3 );
    vec_minus(v2, x3,x1, 3 );
    a1=  v1[1]*v2[2] - v1[2]*v2[1];
    a2= -v1[0]*v2[2] + v1[2]*v2[0];
    a3=  v1[0]*v2[1] - v1[1]*v2[0];
    return 0.5 * sqrt( a1*a1 + a2*a2 + a3*a3 );
}
double TetVolum(double x1[], double x2[], double x3[], double x4[])
{
    double v1[3], dx1[3],dx2[3],dx3[3];
    vec_minus( dx1, x2, x1, 3);
    vec_minus( dx2, x3, x1, 3);
    vec_minus( dx3, x4, x1, 3);
    vec_cross( v1,  dx1,dx2);
    return 1./6.*fabs( vec_dot( v1, dx3, 3 ) );
}


int CycasSolver::CellFaceInfo()
{
    int i,j, n1,n2,n3,n4, ic1,ic2,ic, iface,nc,jf,jc;
    double area1,area2, r1,r2, 
       xc1[3],xc2[3],xc3[3],xc4[3],xc5[3],xc6[3], dx[3],dx1[3],dx2[3], 
       xv1[3],xv2[3],xv3[3],xv4[3],xv5[3],xv6[3],xv7[3],xv8[3], 
       V1,V2,V3,V4,V5,V6;

    // cell ceneter estimation, not final ones
    for( i=0; i<Ncel; i++ )
    {
        for( j=0; j<3; j++ )
        Cell[i].x[j]=1./8.*( Vert[ Cell[i].vertices[0] ][j] + 
                             Vert[ Cell[i].vertices[1] ][j] + 
                             Vert[ Cell[i].vertices[2] ][j] + 
                             Vert[ Cell[i].vertices[3] ][j] + 
                             Vert[ Cell[i].vertices[4] ][j] + 
                             Vert[ Cell[i].vertices[5] ][j] + 
                             Vert[ Cell[i].vertices[6] ][j] +
                             Vert[ Cell[i].vertices[7] ][j] );
    }

    //-- face
    for( i=0; i<Nfac; i++ )
    {
        n1= Face[i].vertices[0];
        n2= Face[i].vertices[1];
        n3= Face[i].vertices[2];
        n4= Face[i].vertices[3];

        // x, area, n, lambda
        // face center, (gravity center???) I doubt that ???
        area1 = Tri3DArea( Vert[n1], Vert[n2], Vert[n4] );
        area2 = Tri3DArea( Vert[n2], Vert[n3], Vert[n4] );

        // face center is bary-center
        for(j=0;j<3;j++){
        xc1[j]= 1./3*( Vert[n1][j]+Vert[n2][j]+Vert[n4][j] );
        xc2[j]= 1./3*( Vert[n2][j]+Vert[n3][j]+Vert[n4][j] );
        }
        for(j=0;j<3;j++)
            Face[i].x[j] = ( area1*xc1[j] + area2*xc2[j] ) / (area1+area2);
    
        // face center is bary-center
        /*  for(j=0;j<3;j++)
        {
            if( n3==n4 ) // triangle
                Face[i].x[j]=1./3.*(Vert[n1][j]+Vert[n2][j]+Vert[n3][j] );
            else         // quadrilateral
                Face[i].x[j]=1./4.*(Vert[n1][j]+Vert[n2][j]+Vert[n3][j]+Vert[n4][j]);
        }  */

        // face normal vector and area
        vec_minus( dx1, Vert[n1], Vert[n3], 3 ); // dx1= vert[n1] - vert[n3]
        vec_minus( dx2, Vert[n2], Vert[n4], 3 ); // dx2= vert[n2] - vert[n4]
        vec_cross( Face[i].n, dx1, dx2 );     // face.n is normal area
        for( j=0;j<3;j++ )
            Face[i].n[j] = 0.5 * Face[i].n[j];
        ic = Face[i].cell1;
        vec_minus( dx, Face[i].x, Cell[ic].x, 3);
        if( vec_dot(Face[i].n, dx,3)<0. ){   // normal points from cell1 to cell2
            for( j=0;j<3;j++ )
                Face[i].n[j] =  -Face[i].n[j] ;
        }
        Face[i].area= vec_len( Face[i].n,3 );
        if( fabs( (area1+area2-Face[i].area)/Face[i].area )>1.e-4 ){
            cout<<i<<" area is not correct "<<area1+area2<<" "<<Face[i].area<<endl;
            cout<<Face[i].x[0]<<" "<<Face[i].x[1]<<" "<<Face[i].x[2]<<endl;
            cout<<Face[i].vertices[0]<<" "<<Face[i].vertices[1]<<" "
                <<Face[i].vertices[2]<<" "<<Face[i].vertices[3]<<endl;
            exit(0);
        }

        // face interpolation
        ic1= Face[i].cell1;
        ic2= Face[i].cell2;
        if( ic2>=0 ){
            vec_minus( dx1, Face[i].x, Cell[ic1].x, 3);
            vec_minus( dx2, Face[i].x, Cell[ic2].x, 3);
            r1 = vec_len( dx1,3 );
            r2 = vec_len( dx2,3 );
            Face[i].lambda = r2 / (r1+r2); // inverse-distance = (1/r1) / ( 1/r1 + 1/r2 )
        }
        else
            Face[i].lambda = 1.;
    }
    
    

    //-- cell
    for( i=0; i<Ncel; i++ )
    {
        for( j=0; j<3; j++ )
        {
            xv1[j] = Vert[ Cell[i].vertices[0] ] [j];
            xv2[j] = Vert[ Cell[i].vertices[1] ] [j];
            xv3[j] = Vert[ Cell[i].vertices[2] ] [j];
            xv4[j] = Vert[ Cell[i].vertices[3] ] [j];
            xv5[j] = Vert[ Cell[i].vertices[4] ] [j];
            xv6[j] = Vert[ Cell[i].vertices[5] ] [j];
            xv7[j] = Vert[ Cell[i].vertices[6] ] [j];
            xv8[j] = Vert[ Cell[i].vertices[7] ] [j];
        }

        // cell-gravity
        // Hexa is divided into 6 Tets
        V1= TetVolum(xv1,xv2,xv4,xv6);
        V2= TetVolum(xv1,xv5,xv4,xv6);
        V3= TetVolum(xv4,xv5,xv8,xv6);
        V4= TetVolum(xv2,xv3,xv4,xv7);
        V5= TetVolum(xv2,xv6,xv4,xv7);
        V6= TetVolum(xv8,xv4,xv6,xv7);
        Cell[i].vol= V1 + V2 + V3 + V4 + V5 + V6;

        for( j=0; j<3; j++ ){
            xc1[j]= 0.25*( xv1[j]+xv2[j]+xv4[j]+xv6[j] );
            xc2[j]= 0.25*( xv1[j]+xv5[j]+xv4[j]+xv6[j] );
            xc3[j]= 0.25*( xv4[j]+xv5[j]+xv8[j]+xv6[j] );
            xc4[j]= 0.25*( xv2[j]+xv3[j]+xv4[j]+xv7[j] );
            xc5[j]= 0.25*( xv2[j]+xv6[j]+xv4[j]+xv7[j] );
            xc6[j]= 0.25*( xv8[j]+xv4[j]+xv6[j]+xv7[j] );
        }
        double Vsum = V1+V2+V3+V4+V5+V6;
        for( j=0; j<3; j++ )
            Cell[i].x[j]=1./Vsum*( V1*xc1[j] + V2*xc2[j] +
                                   V3*xc3[j] + V4*xc4[j] +
                                   V5*xc5[j] + V6*xc6[j] );
        if( fabs((Cell[i].vol - (V1+V2+V3+V4+V5+V6))/Cell[i].vol) > 1.e-3 ){
            cout<< "error in cell volume calculation."<<i<<" "
                <<Cell[i].vol<<" "<<V1+V2+V3+V4+V5+V6<<endl;
            exit(0);
        }

        // neighbored cells
        for( j=0; j<Cell[i].nface; j++ )
        {
            iface = Cell[i].face[j];
            nc    = Face[iface].cell1;
            if( nc==i )
                nc= Face[iface].cell2;
            Cell[i].cell[j]= nc;
        }
    }

    
    //-- boundary
    for( i=0; i<Nbnd; i++ )
    {
        jf  = Bnd [i ].face  ;
        jc  = Face[jf].cell1 ;
        vec_minus( dx, Face[jf].x, Cell[jc].x,3 );
        Bnd[i].distance =  vec_dot(dx, Face[jf].n, 3) / Face[jf].area;
        if( Bnd[i].distance<0. )
        {
            cout<< "wall distance is negative." << i ;
            exit(0);
        }
    }

     //-- Face[].rlencos
    for( i=0; i<Nfac; i++ )
    {
        ic1= Face[i].cell1;
        ic2= Face[i].cell2;
        if( ic2<0 )
        {
            vec_minus( dx,Face[i].x,Cell[ic1].x,3 );
            // Face[i].rlencos= Face[i].area/vec_len(dx,3);
            Face[i].rlencos = Face[i].area / (vec_dot( dx,Face[i].n,3 )/Face[i].area);
        }
        else
        {
            vec_minus( dx,Cell[ic2].x,Cell[ic1].x,3 );
            // Face[i].rlencos= Face[i].area/vec_len(dx,3);
            Face[i].rlencos = Face[i].area / (vec_dot( dx,Face[i].n,3 )/Face[i].area);
        }

        // left and right auxiliary points Xpac, Xnac
        if( ic2>=0 ){
        double dx[3], ss;
        vec_minus( dx, Face[i].x, Cell[ic1].x, 3);
        ss = vec_dot(dx,Face[i].n,3)/ Face[i].area;
        for( j=0; j<3; j++ )
            dx[j] = ss * Face[i].n[j]/Face[i].area;
        vec_minus( Face[i].Xpac, Face[i].x, dx, 3);

        vec_minus( dx, Face[i].x, Cell[ic2].x, 3);
        ss = vec_dot(dx,Face[i].n,3)/ Face[i].area;
        for( j=0; j<3; j++ )
            dx[j] = ss * Face[i].n[j]/Face[i].area;
        vec_minus( Face[i].Xnac, Face[i].x, dx, 3);
        }
    }

    return 1;
}
