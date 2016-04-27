#include <iostream>
#include <string>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include "PreProcessor.h"
#include "HeatConductionSolver.h"
#include "tools.h"

using namespace std;

int PreProcessor::ReadGridFile( )
{
    string   line;
    ifstream file;
    int      NumNodeInCell[8]={0,2,3,4,4,8,6,5},vertices[8],
		i, count,elem_type,ntags,p,tag[10];
		// bndType[10]={4,2,3,1,4,4,0}; // bndType change the tag in gmsh file to navier_bc types ( 4 types currently )

	file.open (GridFileName);
    if( ! file.is_open() )
	{
        errorHandler->fatalRuntimeError("grid file not found",GridFileName);
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
    Bnd = new BoundaryData;//must allocate then realloc...
    Cell = new CellData;
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
void PreProcessor::OutputGrid()
{
    int i,j;
    ofstream of("grid.dat");
    of<<"variables="<<"\"x\","<<"\"y\","<<"\"z\""
        <<"\"rid\","<<endl;

    of<<"zone n="<<Nvrt<<", e="<<Ncel<<", VARLOCATION=([1-3]=NODAL,[4]=CELLCENTERED)"
        <<"DATAPACKING=BLOCK, ZONETYPE=FEBRICK"
        <<endl;

    for( j=0; j<3; j++ ){
        for( i=0; i<Nvrt; i++ ){
            of<<Vert[i][j]<<" ";
            if( i%5==4 ) of<<endl;
        }
        of<<endl;
    }
    for(int i=0; i<Ncel; i++ ){
        of<<Cell[i].rid<<"  ";
        if( i%5==0 ) of<<endl;
    }
    for( i=0; i<Ncel; i++ ){
        for( j=0;j<8;j++ )
            of<<Cell[i].vertices[j]+1<<" ";
        of<<endl;
    }

    of.close();
}


int PreProcessor::CreateFaces( )
{
    int i,j,id,n1,n2,n3,n4,n5,n6,n7,n8;
    int *NumNodeFace, **NodeFace;
    
    NumNodeFace = new int[Nvrt];
    NodeFace    = new_Array2D<int>( Nvrt,500 );
	for( i=0; i<Nvrt; i++ )
		NumNodeFace[i] = 0;
	for( i=0; i<Ncel; i++ )
		Cell[i].nface = 0;

    cout<<"begin faces construction ... ";
    Nfac = 0;
    Face = new FaceData;//alloc before realloc
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

void PreProcessor::FindFace( int ic, int n1, int n2, int n3, int n4, int &nf, 
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

int PreProcessor::CellFaceInfo()
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

//------this function is usded when solve3DHeatConduction---------
// input
// output
//----------------------------------------------------------------
int PreProcessor::findCoupledBoundary(int* new2Old,int coupledBoundId, int nc,//input
            CellData** retCells,FaceData** retFaces,BoundaryData** retBnd,double*** retVert,
            int& retnc, int& retnf, int& retnb,int& retncb, int& retnv)//output
{
        int nv=0,nf=0,nb=0;
        CellData* _cells = new CellData[nc];
        double **_verts = new_Array2D<double>(Nvrt,3);
        FaceData* _faces = new FaceData[Nfac];
        vector<BoundaryData> _bounds;
        _bounds.reserve(Nbnd);

        map<int,int> old2NewV;
        map<int,int> old2NewF;

        for(int i=0;i!=nc;++i){//create vert & cell
            _cells[i] = Cell[new2Old[i]]; 
            for(int j=0;j!=8;++j){
                int iv = _cells[i].vertices[j];
                if(old2NewV.find(iv)!=old2NewV.end()){
                    _cells[i].vertices[j] = old2NewV[iv];
                }else{
                    for(int k=0;k!=3;++k){
                        _verts[nv][k] = Vert[iv][k]; 
                    }
                    old2NewV.insert(make_pair(iv,nv));
                    nv++;
                }
            }
        }
        for(int i=0;i!=nc;++i){//create face
            for(int j=0;j!=_cells[i].nface;++j){
                int iface = _cells[i].face[j];
                if(old2NewF.find(iface)!=old2NewF.end()){
                    int newIface = old2NewF[iface];
                    _faces[newIface].cell2 = i;
                    _cells[i].face[j] = newIface;
                }else{
                    _faces[nf] = Face[iface];
                    _faces[nf].cell1 = i; 
                    _faces[nf].cell2 = -10; //to be set
                    _cells[i].face[j] = nf;
                    for(int k=0;k!=4;++k){
                        int iv = _faces[nf].vertices[k];
                        _faces[nf].vertices[k] = old2NewV[iv];
                    }
                    old2NewF.insert(make_pair(iface,nf));
                    nf++;
                }
            }
        }
        for(int i=0;i!=nc;++i){//config neighbouring
            for(int j=0;j!=_cells[i].nface;++j){
                int iface = _cells[i].face[j];//local
                if(_faces[iface].cell1 == i){
                    _cells[i].cell[j] = _faces[iface].cell2;
                }else if(_faces[iface].cell2 == i){
                    _cells[i].cell[j] = _faces[iface].cell1;
                }else {
                    assert(false); 
                }
            }
        }

        for(int i=0;i!=nf;++i){ //create coupled bound
            if(_faces[i].cell2<0 && _faces[i].bnd < 0){
                BoundaryData toPush;
                for(int k=0;k!=4;++k){
                    toPush.vertices[k] = _faces[i].vertices[k];
                }
                toPush.face = i;
                toPush.rid = coupledBoundId;
                _bounds.push_back(toPush);

                _faces[i].bnd = nb;
                nb++;
            }
        }

        int ncb = nb;
        for(int i=0;i!=nf;++i){
            if(_faces[i].cell2<0 && _faces[i].bnd > 0){
                BoundaryData toPush = Bnd[_faces[i].bnd];
                toPush.face = i;
                _bounds.push_back(toPush);
                _faces[i].bnd = nb;
                nb++;
            }
        }

        //compress size and return;
        retnc = nc;
        (*retCells) = _cells;
        _cells = NULL;

        retnv = nv;
        (*retVert) = new_Array2D<double>(nv,3);
        for(int i=0;i!=nv;++i){
            for(int j=0;j!=3;++j){
                (*retVert)[i][j] = _verts[i][j];
            }
        }
        delete_Array2D(_verts,Nvrt,3);

        retnf = nf;
        (*retFaces) = new FaceData[nf];
        for(int i=0;i!=nf;++i){
            (*retFaces)[i] = _faces[i];
        }
        delete []_faces;

        retnb = nb;
        retncb = ncb;
        (*retBnd) = new BoundaryData[nb]; 
        for(int i=0;i!=nb;++i){
            (*retBnd)[i] = _bounds[i];
        }

        return 0;
}

//------------- Build Solver Geometry and other param --------//
int PreProcessor::buildSolver(NavierStokesSolver* ns){
    if(!solve3DHeatConduction){
        //build NS_Solver
        ns->Vert = Vert;
        ns->Face = Face;
        ns->Cell = Cell;
        ns->Bnd  =  Bnd;
        ns->Nvrt = Nvrt;
        ns->Ncel = Ncel;
        ns->Nfac = Nfac;
        ns->Nbnd = Nbnd;

    }else{
        CycasSolver* sv = new CycasSolver;
        *sv = *ns;//copy
        HeatConductionSolver* ht = new HeatConductionSolver(*sv);//ht get all the basic parameter from ns
        delete sv;

        int coupledBoundId = -1;
        for(int i=0;i!=ns->regionMap.size()+1;++i){//insert a coupled bound region, a Wall in fact;
            if(ns->regionMap.find(i)==ns->regionMap.end()){
                BdRegion& reg = ns->regionMap[i];
                reg.name = "coupled boundary";
                reg.type1 = 1;
                reg.type2 = 2; //coupled
                coupledBoundId = i;
                break;
            }
        }

        int nc=0;
        int* new2Old = new int[Ncel];

        int ncS = 0;
        int* new2OldS = new int[Ncel];

        for(int i=0;i!=Ncel;++i){
            if(ns->regionMap.find(Cell[i].rid) == ns->regionMap.end()){
                errorHandler->fatalRuntimeError("no such rid for body, rid:",Cell[i].rid); 
            }
            BdRegion& reg = ns->regionMap[Cell[i].rid];
            if(reg.type1!=5){
                errorHandler->fatalRuntimeError("body should has rid 5"); 
            }
            if(reg.type2==0){//fluid
                new2Old[nc++] = i;
            }else{//solid
                new2OldS[ncS++] = i;
            }
        } 

        findCoupledBoundary(new2Old, coupledBoundId, nc,
                &ns->Cell, &ns->Face, &ns->Bnd, &ns->Vert,
                ns->Ncel,ns->Nfac,ns->Nbnd,ns->NCoupledBnd,ns->Nvrt
            );

        findCoupledBoundary(new2OldS, coupledBoundId, ncS,
                &ht->Cell, &ht->Face, &ht->Bnd, &ht->Vert,
                ht->Ncel,ht->Nfac,ht->Nbnd,ht->NCoupledBnd,ht->Nvrt
            );

        //advance config to ht;
        ht->URF = 1.0;  //seems useless

        ns->physicalModule["3dHeatConduction"] = ht;//apply module
        ht=NULL;

        delete []new2Old;
        delete []new2OldS;
    }
    return 0;
}

