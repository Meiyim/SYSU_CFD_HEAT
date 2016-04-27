#pragma once

#include "BasicType.h"

class HeatConductionSolver: public CycasSolver
{
public:
	HeatConductionSolver();
	HeatConductionSolver(CycasSolver& ori);

	~HeatConductionSolver();
	// relaxation factor
	double URF;

    // physical variables at cell center
    double* diffCoef;	//for each kind of material
    double   *Rn;
	double   *Tn;     // primitive vars
	double **dPhidX;
	// dual time unsteady simulation backup data, p = previous
	double   *Tnp,    // Euler
		    *Tnp2 ;   // BDF 2nd
	// variables at face center
	// Boundary faces values
	double  *BTem;
    
	// laspack library work array
	QMatrix As;                // As for non-sysmmetric, Ap for sysmmetric matrix
	Vector  bs,xsol; // right-hand-side vector


	//Basic Solver interface
	virtual int  CheckAndAllocate( );
    // read solver param, material, post, everything except 
    virtual void InitFlowField  ( );
    //tempreture calculation called by nsSolver
    virtual int solve();
	// Post process output zone
    virtual void Output2Tecplot ( std::ofstream& ofile);
    // time scheme
    virtual void SaveTransientOldData();
    virtual double getResidule(); 

private:
	// scalar. temperature, other passive variables
	// void ScalarTranport   ( double *Phi, double *BPhi, double *DiffCoef, double *source );
	void BuildScalarMatrix(int iSca,double *Phi, double *BPhi, double *DiffCoef, double *source, double *App);
    // iSca= 1:temperature; =2: ke; =3: ed;  = 4-larger: species concentration
	void UpdateEnergy     ( );

	// Boundary condition. one defact is arrays cost too much memory, especially 2D
	void SetBCTemperature( double *bt );



};

extern class HeatConductionSolver hcSolver;