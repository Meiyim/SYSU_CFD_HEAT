#pragma once
#include <stdlib.h>
#include "tools.h"
#include <string>
#include <map>
#include <assert.h>
extern "C"{ 
#include "laspack/laspack.h"
}
#include "terminalPrinter.h"

#define SMALL 1.e-16
// geometry, face & cell data
typedef struct 
{
    int    bnd;
    int    vertices[4];
    int    cell1,cell2;
    double x[3],n[3],area;   // face center and normal vector
    double lambda,   // lambda for left, (1.-lambda) for right
		rlencos;     // area/(|Xpn|*vect_cosangle(n,Xpn))
	double Xpac[3],Xnac[3]; // auxiliary points
}FaceData;
typedef struct
{
    int nface;
    int rid;   //same as boundary data
    int face[6], cell[6], vertices[8]; // maybe wasterful a bit. Make it dynamics to save memory
    //face[nface]: for the faces index
    //  all elements are treated as 8 nodes hexahdron,
    //  but with different number of faces, so judges will be used for avoid faces with one vertex
    double vol, x[3];
}CellData;
// connect boundary to faces, boundary to BdRegion
typedef struct
{
    int face;                // belongs to face...
    int vertices[4];         // the 4 vertices, to be done allocatable
    int rid;                 // region id as set in rtable
    double distance;         // normal distance from cell face
                             // center to cell center
    double  yplus;           // y+
    double  uplus;           // u+
    double  shear[3];        // shearstress components
    double  h;               // local heattransfer coef.
    double  q;               // local heat flux (in W/m2)
    double  T;               // local wall temperature
}BoundaryData;
// store the data set in one type of boundary
typedef struct
{
	std::string name;
	int  type1;  //=1: wall
				 //=4: symmetric
				 //=5: body
	// temperature
	int  type2; // for wall:(type1==1)
				// =0 : given temperature; 
	            // =1 : heat flux (0 for default adiabatic, others )
				// for body:(type1==5)
				// =0 : fluid field: u, v, w, ro, p, t ...
				// =1 : solid: ro, diffCoef

	double fixedValue;
	double fixedVector[3];
	double initvalues[10];

}BdRegion;


class HeadConductionSolver
{
public:
	HeadConductionSolver():outputCounter(0),printer(new TerminalPrinter){
	};
	~HeadConductionSolver();
    
	// this zone is added by CHENXUYI
	size_t outputCounter;
	TerminalPrinter* printer;
	// above is added by CHENXUYI
	char   GridFileName[100];
    map<int,BdRegion> regionMap; // record the region/bound info, Bnd[i].rid ---> regionMap

    // control param
    bool   IfSteady;
	int    TimeScheme;
	    // TurModel=0: laminar flow; =1: k-epsilon; =2: wait to implement
		// DensityModel=0: constant density; =1: perfect gas ro=p/(RT); =2: wait to implement
		// TimeScheme=1, Euler (default) ;  =2, BDF 2nd

    double* diffCoef;	//for each kind of material


	// relaxation factor
	double URF;
	int    limiter,MaxOuterStep;  // Maximum step at outer iteration

    // --- time evolution
	int    step,MaxStep;
    double dt, cur_time,total_time, Residual,ResidualSteady;

	// --- post process
	int    noutput,outputFormat;

    // geometry
    int          Nvrt, Ncel, Nfac, Nbnd;
    double       **Vert; // coordinate x,y,z
    FaceData     *Face;
    CellData     *Cell;
    BoundaryData *Bnd;

    // physical variables at cell center
    double   *Rn;
	double   *Tn;     // primitive vars
	double **dPhidX;
	// variables at face center
	// Boundary faces values
	double  *BTem;
	double Tin, Twall;
    
	// laspack library work array
	QMatrix As;                // As for non-sysmmetric, Ap for sysmmetric matrix
	Vector  bs,bu,bv,bw,bp, xsol; // right-hand-side vector

	// dual time unsteady simulation backup data, p = previous
	double   *Tnp, *Rnp,   // Euler
		    *Tnp2, *Rnp2;   // BDF 2nd

// Geometry
    int  ReadGridFile    ( );
	void OutputGrid      ( );
    int  CreateFaces     ( );
    void FindFace( int, int,int,int,int, int&, int*,int** );
    int  CellFaceInfo    ( );
	int  CheckAndAllocate( );

// Init flow field
    // read solver param, material, post, everything except 
    void InitSolverParam( );
	void ReadParamFile   ( );
    void InitFlowField  ( );
	


// Fluid calculation
    void NSSolve ( );

	// pressure


	// scalar. temperature, other passive variables
	// void ScalarTranport   ( double *Phi, double *BPhi, double *DiffCoef, double *source );
	void BuildScalarMatrix(int iSca,double *Phi, double *BPhi, double *DiffCoef, double *source, double *App);
	     // iSca= 1:temperature; =2: ke; =3: ed;  = 4-larger: species concentration
	void UpdateEnergy     ( );


	void SaveTransientOldData( );

    // gradients and divergence
	int  Gradient    ( double*, double*, double** );
    // int  Divergence  ( double*, double*, double*, double*, double*, double*, double* );
    int  Limiter_MLP  (double[],double **);
	int  Limiter_Barth(double[],double **);
    int  Limiter_WENO (double[],double **);

	// Boundary condition. one defact is arrays cost too much memory, especially 2D

	void SetBCTemperature( double *bt );

// Post process
    void OutputMoniter  ( );
    void Output2Tecplot ( );
	void Output2VTK     ( );
	void WriteBackupFile( );
	void ReadBackupFile ( );

private:

};

extern class HeadConductionSolver hcSolver;