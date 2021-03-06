#ifndef _PREPROCESSOR_H_
#define _PREPROCESSOR_H_
#include "BasicType.h"
#include "navier.h"

class PreProcessor{
    int Nvrt,Ncel,Nfac,Nbnd;

    double       **Vert; // coordinate x,y,z
    FaceData     *Face;
    CellData     *Cell;
    BoundaryData *Bnd;  

    char   GridFileName[100];
    std::string paramFileName;

    //parameters
    bool solve3DHeatConduction;
public:
    PreProcessor():
        solve3DHeatConduction(false),
        paramFileName("param.in")
    {}
	//param
    void parseCommandLine(int argc, char** argv);
   void ReadParamFile   (NavierStokesSolver* ns);
	// Geometry
    int  ReadGridFile    ( );
	void OutputGrid      ( );
    int  CreateFaces     ( );
    void FindFace        ( int, int,int,int,int, int&, int*,int** );


    //build solver 
    int buildSolver(NavierStokesSolver* ns);
private:
    int findCoupledBoundary(int* new2Old,int coupledBoundId, int nc,//input
        CellData** retCells,FaceData** retFaces,BoundaryData** retBnd,double*** retVert, //ouput
        int& retnc, int& retnf, int& retnb, int& retncb, int& retnv);//output

};

#endif