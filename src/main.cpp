#include <iostream>
#include "BasicType.h"
#include "PreProcessor.h"
#include "navier.h"

ErrorHandler* errorHandler;

int main( int argc, char *argv[] )
{
	// read input card, configure(original) & grid(gmsh)
	errorHandler = new ErrorHandler;
	NavierStokesSolver* nsSolver = new NavierStokesSolver;
	PreProcessor* preProcess = new PreProcessor;
	nsSolver->printer->printStarter();
	
	// read mesh
	preProcess->ReadParamFile(nsSolver);
	preProcess->ReadGridFile();
	// construct faces, compute geometry related info
	preProcess->CreateFaces ( );
	preProcess->CellFaceInfo( );

	preProcess->buildSolver(nsSolver);

	nsSolver->CheckAndAllocate( );

	// init flow field
	nsSolver->InitFlowField();

	//------------ MAIN CFD SOLVE -----------//
	nsSolver->NSSolve();	

	//------------ Post Process -----------//
	printf("done \n");
	delete nsSolver;
	delete errorHandler;
	return 0;
}

