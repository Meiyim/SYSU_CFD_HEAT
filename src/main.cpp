#include <iostream>
#include "navier.h"

int startup( );
int simulate();
int postprocess();

int main( int argc, char *argv[] )
{
    startup ( );
    
    simulate( );
    
    postprocess( );
    return 0;
}

int startup( )
{
	// read input card, configure(original) & grid(gmsh)
	hcSolver.printer->printStarter();
	hcSolver.InitSolverParam( );
	// init constant
	
	// read mesh
	hcSolver.ReadGridFile( );
	// construct faces, compute geometry related info
	hcSolver.CreateFaces ( );
	hcSolver.CellFaceInfo( );
	hcSolver.CheckAndAllocate( );

	// init flow field
	hcSolver.InitFlowField();
	return 0;
}


int simulate( )
{
	// other modules

	// material property
	// GasProperties( );
	
	// Main CFD
	hcSolver.NSSolve();
	
	return 0;
}


int postprocess()
{
    return 0;
}



