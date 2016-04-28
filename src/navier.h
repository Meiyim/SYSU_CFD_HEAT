#pragma once
#include <stdlib.h>
#include "BasicType.h"
#include "tools.h"
extern "C"{ 
#include "laspack/laspack.h"
}
#include "terminalPrinter.h"

#define SMALL 1.e-16

class NavierStokesSolver: public CycasSolver
{
public:
	NavierStokesSolver();
	~NavierStokesSolver();
    
	size_t outputCounter;
	TerminalPrinter* printer;

	//Physical module
	map<std::string,CycasSolver* > physicalModule;

    // control param
    bool   SolveEnergy,SolveSpecies;
	int    TurModel,DensityModel, TimeScheme;
	    // TurModel=0: laminar flow; =1: k-epsilon; =2: wait to implement
		// DensityModel=0: constant density; =1: perfect gas ro=p/(RT); =2: wait to implement
		// TimeScheme=1, Euler (default) ;  =2, BDF 2nd
    double gama,ga1, cp,cv,prl,prte,Rcpcv, TempRef, gravity[3];
	double PressureReference;
	int    cellPressureRef;

	// --- numerical scheme
	// relaxation factor
	double URF[8];

    // --- time evolution
    double cur_time,total_time, Residual[10],ResidualSteady;
    int    step,MaxStep, MaxOuterStep;
    // --- post process
    int    noutput,outputformat;

    // physical variables at cell center
	double  *Rn, *Un, *Vn, *Wn,  *Pn, *Tn, *TE, *ED;     // primitive vars
	int     Nspecies;
	double  **RSn;                             // species density
	double  *VisLam, *VisTur;
	double  **dPdX, **dUdX,**dVdX,**dWdX, *Apr, **dPhidX;
	// variables at face center
	double  *RUFace;   // RUFace[Nfac]
	// Boundary faces values
	double  *BRo,*BU,*BV,*BW,*BPre, *BTem,**BRS, *BTE,*BED;
    
	// laspack library work array
	QMatrix As,Ap;                // As for non-sysmmetric, Ap for sysmmetric matrix
	Vector  bs,bu,bv,bw,bp, xsol; // right-hand-side vector

	// dual time unsteady simulation backup data, p = previous
	double  *Rnp, *Unp,  *Vnp,  *Wnp,  *Tnp,  *TEp,  *EDp, **RSnp,    // Euler
		    *Rnp2,*Unp2, *Vnp2, *Wnp2, *Tnp2, *TEp2, *EDp2,**RSnp2;   // BDF 2nd




// BasicSolver Interface
	virtual int  CheckAndAllocate( );
// Init flow field
    // read solver param, material, post, everything except 
    virtual void InitFlowField  ( );
// Post process
    virtual void Output2Tecplot (std::ofstream& ofile,int nvar);
// time scheme   
	virtual void SaveTransientOldData( );
    virtual int solve(){return 0;}; //NS_Solver dont need this
    virtual double getResidule() { return 0.0;};  //NS_Solver dont need this

// Fluid calculation
    void NSSolve ( );

    // velocity
    int  CalculateVelocity( );
	void BuildVelocityMatrix( );
	void CalRUFace ( );
	void CalRUFace2( );
	// pressure
    int  CalculatePressure( );
    void BuildPressureMatrix( );
	void CorrectRUFace2(double*);
	// scalar. temperature, other passive variables
	// void ScalarTranport   ( double *Phi, double *BPhi, double *DiffCoef, double *source );
	void BuildScalarMatrix(int iSca,double *Phi, double *BPhi, double *DiffCoef, double *source, double *App);
	     // iSca= 1:temperature; =2: ke; =3: ed;  = 4-larger: species concentration
	void UpdateEnergy     ( );
	void UpdateSpecies    ( );
	void UpdateTurKEpsilon( );


	// Boundary condition. one defact is arrays cost too much memory, especially 2D
	void SetBCVelocity( double*br,double*bu,double*bv,double*bw );
	void SetBCPressure( double*bp );
	void SetBCDeltaP  ( double*bp, double *dp );
	void SetBCTemperature( double *bt );
	void SetBCSpecies ( double **brs );
	void SetBCKEpsilon( double *TESource,double *EDSource,double *ApTE,double*ApED,double *Prod);

	//PostProcess
	bool shouldPostProcess(int,int,double);
	void Output();


private:

};

namespace TurKEpsilonVar
{
	const double 
	kappa    = 0.419 ,
	Cmu      = 0.09  ,
    Ceps1    = 1.44  ,
	Ceps2    = 1.92  ,
    LenSc    = 0.1   ,
    Sigma_k  = 1.0   ,    //! turbulence diff. coef. factors
    Sigma_e  = 1.3   ,    //! ie. turbulent Prandtl numbers
    Sigma_s  = 0.9   ;
}

