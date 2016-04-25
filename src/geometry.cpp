#include <iostream>
#include <string>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include "navier.h"
#include "tools.h"

using namespace std;

int HeadConductionSolver::ReadGridFile( )
{
    string   line;
    ifstream file;
    int      NumNodeInCell[8]={0,2,3,4,4,8,6,5},vertices[8],
		i, count,elem_type,ntags,p,tag[10];
		// bndType[10]={4,2,3,1,4,4,0}; // bndType change the tag in gmsh file to navier_bc types ( 4 types currently )

	file.open (GridFileName);
    if( ! file.is_open() )
	{
		cout<<"Grid file "<<GridFileName<<" cannot be opened. It does not exist or is busy!"<<endl;
		exit(0);
	}
    // Skip some lines
    file >> line;
    file >> line;
    file >> line;
    file >> line;
    file >> line;
    file >> line;
    // Read vertices
    file >> Nvrt;
    assert (Nvrt > 0);
    // realloc (Vert, Nvrt);
	Vert = new_Array2D<double>(Nvrt,3);
    for(i=0; i<Nvrt; i++){
        file >> count
             >> Vert[i][0]
             >> Vert[i][1]
             >> Vert[i][2];
    }
    file >> line;
    file >> line;
    // Read cells
    int n_elem;
    file >> n_elem;
    assert (n_elem > 0);
    
    Ncel = 0;
    Nbnd = 0;
    for(i=0; i<n_elem; ++i)
    {
        file >> count
             >> elem_type
             >> ntags;
		for( p=0; p<ntags; p++ )    // Dummy tags, maybe useful later
            file >> tag[p];
		for( p=0; p<NumNodeInCell[elem_type]; p++ ){
            file >> vertices[p];
			vertices[p] -- ;  // to agree with vertices numbering from 0 to Nvrt-1
		}

        // Some gmsh files have 2 and some have 3 tags
        assert( ntags==2 || ntags==3 );
  
        if( elem_type==2 || elem_type==3 ) // Boundary, Tri/Quad face
        {
            Bnd = (BoundaryData*) realloc(Bnd,(Nbnd+1)*sizeof(BoundaryData));
			Bnd[Nbnd].rid = tag[0];  // bndType[tag[0]];      // First tag is face type
            if(     elem_type==2 )
            {
                Bnd[Nbnd].vertices[0]= vertices[0];
                Bnd[Nbnd].vertices[1]= vertices[1];
                Bnd[Nbnd].vertices[2]= vertices[2];
                Bnd[Nbnd].vertices[3]= vertices[2];
            }
            else if( elem_type==3 )
            {
                Bnd[Nbnd].vertices[0]= vertices[0];
                Bnd[Nbnd].vertices[1]= vertices[1];
                Bnd[Nbnd].vertices[2]= vertices[2];
                Bnd[Nbnd].vertices[3]= vertices[3];
            }
			else
			{
				cout<<"boundary type is not correct."<<elem_type<<endl;
				abort();
			}
			Nbnd ++ ;
        }
        else if( elem_type==4 || elem_type==5 ||  // Tetrahedral/Hexa/Prism/Pyramid cell
                 elem_type==6 || elem_type==7 )
        {
            Cell = (CellData*) realloc(Cell,(Ncel+1)*sizeof(CellData));
            Cell[Ncel].rid = tag[0];
            // set all the cell types to 8-vertex hexa
            if(     elem_type==4 )   // tetra
            {
				// 0-3
				for( p=0; p<3; p++ )
				Cell[Ncel].vertices[p]= vertices[p];
				Cell[Ncel].vertices[3]= vertices[2];
				// 4-7
                for( p=4; p<8; p++ )
				Cell[Ncel].vertices[p]= vertices[3];
            }
			else if( elem_type==5 )  // hexa
			{
				for( p=0; p<8; p++ )
				Cell[Ncel].vertices[p]= vertices[p];
			}
            else if( elem_type==6 )  // prism
            {
				// 0-3
				for( p=0; p<3; p++ )
                Cell[Ncel].vertices[p]= vertices[p];
				Cell[Ncel].vertices[3]= vertices[2];
				// 4-7
                for( p=4; p<7; p++ )
                Cell[Ncel].vertices[p]= vertices[p-1];
                Cell[Ncel].vertices[7]= vertices[5];
            }
            else if( elem_type==7 )  // pyramid
            {
				// 0-3
				for( p=0; p<4; p++ )
				Cell[Ncel].vertices[p]= vertices[p];
				// 4-7
                for( p=4; p<8; p++ )
                Cell[Ncel].vertices[p]= vertices[4];
            }
			Ncel++ ;
        }
        else
        {
           cout << "Unknown element type !!!" << endl;
           abort();
        }
    }
    file.close ();

	OutputGrid();
    return 0;
}

// output grid to tecplot for validation
void HeadConductionSolver::OutputGrid()
{
	int i,j;
	ofstream of;
	of.open("grid.dat");
	of<<"variables="<<"\"x\","<<"\"y\","<<"\"z\""<<"\"rid\""<<endl;
	// output cells
	of<<"zone n="<<Nvrt<<", e="<<Ncel<<",VARLOCATION=([1-3]=NODAL,[4]=CELLCENTERED)"
    "DATAPACKING=BRICK, ZONETYPE = FEBRICK "<<endl;
    for( j=0; j<3; j++ ){
        for( i=0; i<Nvrt; i++ ){
            of<<Vert[i][j]<<" ";
            if( i%5==4 ) of<<endl;
        }
        of<<endl;
    }
    for(int i=0;i<Ncel;++i){
        of<<Cell[i].rid<<" ";
        if(i%5==4) of<<endl;
    }

    for( i=0; i<Ncel; i++ ){
        for( j=0;j<8;j++ )
            of<<Cell[i].vertices[j]+1<<" ";
        of<<endl;
    }

	of.close( );
}


int HeadConductionSolver::CreateFaces( )
{
    int i,j,id,n1,n2,n3,n4,n5,n6,n7,n8;
    int *NumNodeFace, **NodeFace;
    
    NumNodeFace = new int[Nvrt];
    NodeFace    = new_Array2D<int>( Nvrt,500 );
	for( i=0; i<Nvrt; i++ )
		NumNodeFace[i] = 0;
	for( i=0; i<Ncel; i++ )
		Cell[i].nface = 0;
	for( i=0; i<Nfac; i++ ){
		Face[i].bnd = -10;
		Face[i].lambda = 1.;
	}
    
    cout<<"begin faces construction ... ";
    Nfac = 0;
    // first boundary faces
    for( i=0; i<Nbnd; i++ )
    {
		Face = (FaceData*)realloc( Face, (Nfac+1)*sizeof(FaceData) );
        Bnd[i].face = Nfac;
        for( j=0;j<4;j++ )
            Face[Nfac].vertices[j] = Bnd[i].vertices[j];
        Face[Nfac].bnd     = i   ;  // bnd(i)%rid
        Face[Nfac].cell2   = -10 ;  // cell1 is not known, set in FindFace
        for( j=0;j<4;j++ )
        {
            id = Face[Nfac].vertices[j];
            NodeFace[id][NumNodeFace[id]++]= Nfac;
        }
		Nfac++;  // face accumulation is at last
    }
    // interior faces
    for( i=0;i<Ncel;i++ )
    {
        n1= Cell[i].vertices[0];
        n2= Cell[i].vertices[1];
        n3= Cell[i].vertices[2];
        n4= Cell[i].vertices[3];
        n5= Cell[i].vertices[4];
        n6= Cell[i].vertices[5];
        n7= Cell[i].vertices[6];
        n8= Cell[i].vertices[7];
        FindFace( i,n1,n2,n3,n4, Nfac, NumNodeFace,NodeFace );
        FindFace( i,n5,n6,n7,n8, Nfac, NumNodeFace,NodeFace );
        FindFace( i,n1,n2,n6,n5, Nfac, NumNodeFace,NodeFace );
        FindFace( i,n2,n3,n7,n6, Nfac, NumNodeFace,NodeFace );
        FindFace( i,n4,n3,n7,n8, Nfac, NumNodeFace,NodeFace );
        FindFace( i,n1,n4,n8,n5, Nfac, NumNodeFace,NodeFace );
    }
    cout<< "number of faces : "<< Nfac <<endl; // << maxval(NFaces,Nfac) << minval(NFaces,Nfac) ;
    delete [] NumNodeFace;
    delete_Array2D<int>( NodeFace, Nvrt, 500 );
    cout<<"finish constructing all faces"<<endl;
	return 1;
}

void HeadConductionSolver::FindFace( int ic, int n1, int n2, int n3, int n4, int &nf, 
        int *NumNodeface, int **NodeFace )
{
    int fid,k,irepeat,iface,ifn1,ifn2,ifn3,ifn4,id;
    irepeat = 0;
    if( n1==n2 || n1==n3 || n1==n4 ) irepeat= irepeat + 1;
    if( n2==n3 || n2==n4           ) irepeat= irepeat + 1;
    if( n3==n4                     ) irepeat= irepeat + 1;
    if( irepeat>1 ) return;
    // find 
    fid= -100;
    for( k=0; k<NumNodeface[n1]; k++ )
    {
        iface= NodeFace[n1][k];
        ifn1 = Face[iface].vertices[0];
        ifn2 = Face[iface].vertices[1];
        ifn3 = Face[iface].vertices[2];
        ifn4 = Face[iface].vertices[3];
        if( (ifn1==n1 || ifn1==n2 || ifn1==n3 || ifn1==n4 ) &&
            (ifn2==n1 || ifn2==n2 || ifn2==n3 || ifn2==n4 ) &&
            (ifn3==n1 || ifn3==n2 || ifn3==n3 || ifn3==n4 ) &&
            (ifn4==n1 || ifn4==n2 || ifn4==n3 || ifn4==n4 ) )
        {
            fid= iface;
            break;
        }
    }

    if( fid<0 )
    {
		Face = (FaceData*)realloc( Face, (Nfac+1)*sizeof(FaceData) );
        Face[nf].vertices[0]= n1;
        Face[nf].vertices[1]= n2;
        Face[nf].vertices[2]= n3;
        Face[nf].vertices[3]= n4;
        Face[nf].cell1      = ic;
        Face[nf].bnd        = -10;  // inner face
        Cell[ic].face[ Cell[ic].nface++ ] = nf;
        // add face to NodeFace
        for( k=0; k<4; k++ )
        {
            id = Face[nf].vertices[k];
            NodeFace[id][ NumNodeface[id]++ ] = nf;
        }
		nf = nf + 1; // face accumulation is at last
    }
    else
    {
        if( Face[fid].bnd >= 0 )    // boundary face
            Face[fid].cell1  = ic;
        else                       // inner face
            Face[fid].cell2  = ic;
        
        Cell[ic].face[ Cell[ic].nface++ ] = fid;
    }
}


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

int HeadConductionSolver::CellFaceInfo()
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

		// cell volume
        // way-1: |V|=Integ_V( div(x,0,0) ) = Integ_S( (x,0,0).n ).
        /* Cell[i].vol = 0.;
        for( j=0; j<Cell[i].nface; j++ )
        {
            int iface= Cell[i].face[j];
			if( Face[iface].cell1==i )
				Cell[i].vol += Face[iface].x[0] * Face[iface].n[0];
			else
				Cell[i].vol -= Face[iface].x[0] * Face[iface].n[0];
        } */
		// way-2:
		/* double ddp[3],dde[3],ddz[3];
		for( j=0;j<3;j++ ){
			ddp[j]=0.125*( -xv1[j]+xv2[j]+xv3[j]-xv4[j]-xv5[j]+xv6[j]+xv7[j]-xv8[j] );
			dde[j]=0.125*( -xv1[j]-xv2[j]+xv3[j]+xv4[j]-xv5[j]-xv6[j]+xv7[j]+xv8[j] );
			ddz[j]=0.125*( -xv1[j]-xv2[j]-xv3[j]-xv4[j]+xv5[j]+xv6[j]+xv7[j]+xv8[j] );
		}
		Cell[i].vol = ddp[0]*( dde[1]*ddz[2] - dde[2]*ddz[1] )
					 -ddp[1]*( dde[0]*ddz[2] - dde[2]*ddz[0] )
					 +ddp[2]*( dde[0]*ddz[1] - dde[1]*ddz[0] );
		Cell[i].vol *= 8.; */

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


	/* // output 
	ofstream of;
	of.open("face.dat");
	of<<Nfac<<endl;
	for( i=0; i<Nfac; i++ )
	{
		for( j=0; j<4; j++ )
			of<<Face[i].vertices[j]+1<<" ";
		of<<endl;
		of<<Face[i].n[0]<<" "<<Face[i].n[1]<<" "<<Face[i].n[2]<<" ";
		of<<Face[i].lambda<<" "<<Face[i].cell1+1<<" "<<Face[i].cell2+1;
		of<<endl;
	}
	of.close( );
	of.open("cell.dat");
	of<<Ncel<<endl;
	for( i=0; i<Ncel; i++ ){
		of<<Cell[i].vol<<" "<<Cell[i].x[0]<<" "<<Cell[i].x[1]<<" "<<Cell[i].x[2];
		of<<endl;
	}
	of.close(); */

    
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

int HeadConductionSolver::CheckAndAllocate()
{
	int i,c1,c2;
	// check if all boundaries are marked
	for( i=0; i<Nfac; i++ )
	{
		c1= Face[i].cell1;
		c2= Face[i].cell2;
		if( c2==-10 || c2>=0 ) continue;
		cout<<"error in face right hand side!"<<endl;
	}


	// allocate variables
	Tn = new double[Ncel];
    Rn = new double[Ncel];
    diffCoef = new double[Ncel];
	if( !IfSteady ){
	if( TimeScheme>=1 ){
    	Tnp = new double[Ncel];
        Rnp = new double[Ncel];
	}

	if( TimeScheme>=2 ){
    	Tnp2 = new double[Ncel];
        Rnp2 = new double[Ncel];
	}
	}

	dPhidX = new_Array2D<double>(Ncel,3);

	BTem= new double[Nbnd];

	// laspack working array
	V_Constr(&bs,   "rightU",    Ncel, Normal, True);
	V_Constr(&xsol, "rightU",    Ncel, Normal, True);

	cur_time = 0.;
	return 1;
}
