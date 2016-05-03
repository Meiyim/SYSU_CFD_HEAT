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
    if( ! file.is_open() ){
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
    Face = new FaceData[Nbnd];//alloc before realloc
    // first boundary faces
    cout<<Nbnd<<endl;
    for( i=0; i<Nbnd; i++ )
    {
        Bnd[i].face = Nfac;
        for( j=0;j<4;j++ )
            Face[Nfac].vertices[j] = Bnd[i].vertices[j];
        Face[Nfac].bnd     = i   ;  // bnd(i)%rid
        Face[Nfac].cell1   = -10 ;  // cell1 is not known, set in FindFace
        Face[Nfac].cell2   = -10 ;  
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
        if( Face[fid].bnd >= 0 ){   // boundary face
               if( Face[fid].cell1 == -10 ){
                        Face[fid].cell1 = ic;
               }else if(Face[fid].cell2 == -10){
                        Face[fid].cell2 = ic;
               }else{
                    assert(false); 
               }
        }else{ //inner face
            Face[fid].cell2  = ic;
        }
        Cell[ic].face[ Cell[ic].nface++ ] = fid;
    }
}



//------this function is usded when solve3DHeatConduction---------
// input
// output
// not a optimized implementation, but, whatever, who the fxxk cares
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
                _cells[i].vertices[j] = nv;
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


    /*
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
    */



    for(int i=0;i!=nf;++i){
        if(_faces[i].cell2<0 && _faces[i].bnd >= 0){
            BoundaryData toPush = Bnd[_faces[i].bnd];
            toPush.face = i;
            _bounds.push_back(toPush);
            _faces[i].bnd = nb;
            nb++;
        }
    }

    int ncb = nb;
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
    printf("%d coupled boundary founded\n",nb - ncb);
    printf("%d  boundary  in total\n",nb);


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

        printf("found fluid cell %d, solid cell %d\n",nc,ncS);
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

