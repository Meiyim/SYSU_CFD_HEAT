#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include "navier.h"
#include "tools.h"
#include "terminalPrinter.h"

using namespace std;

class HeadConductionSolver hcSolver;

void HeadConductionSolver::NSSolve( )
{
	int iter;
	double tstart,tend;
	tstart = ttime( );
	if(IfSteady){
		dt = 1.e9;
		MaxOuterStep = MaxStep;
	}
   	 for( step=1; step<=MaxStep; step++ )
  	  {
  	  	for(iter=1;iter<=MaxOuterStep;++iter){
			if( (step-1)%10==0 )printer->printSectionHead(cur_time);
			cur_time += dt;
			if(!IfSteady){
				SaveTransientOldData( );
			}

			// outer iteration, drive residual to zero for every physical time step
			UpdateEnergy();
			printer->printStepStatus(step,iter,cur_time,dt,0.0);
			if( step%noutput==0 ){
				if( outputFormat==0 ) Output2Tecplot ( );  // exit(0);
				WriteBackupFile( );
			}
			
			if(false){//residual small enough
				break;
			}

		}
		if( cur_time >= total_time )break;
	}

   	 Output2Tecplot();
	printer->printEnding();
	
	//if( outputFormat==0 ) Output2Tecplot ( );
	//if( outputFormat==1 ) Output2VTK     ( );
	//WriteBackupFile( );
}

void HeadConductionSolver::InitSolverParam( )
{
	int i;
	// default values, can be reset in initflow
	MaxStep      = 10000 ;
	MaxOuterStep = 50  ;
	ResidualSteady= 1.e-6;

	//-- numerical scheme
	// relaxation factor
	URF = 0.6;  // u
	limiter = 0;

	total_time   = 0.    ;
	dt           = 1.0   ;    // no meaning in steady case
	TimeScheme   = 1     ;

	noutput  = 1;
	outputFormat = 0;  // 0 for tecplot; 1 for vtk
	// init parameters from param.in
	ReadParamFile( );

	// some parameters check work, e.g. 
}


void HeadConductionSolver::InitFlowField( )
{
	int i;

	//if( IfReadBackup ) 
		ReadBackupFile( );
	for( i=0; i<Ncel; i++ )
	{
		int rid = Cell[i].rid;
		assert(regionMap[rid].type1 == 5);//body region
		if( regionMap[rid].type2 == 0 ){//fluid

		}else if(regionMap[rid].type2 == 1){//solid
			double* initvalues;
			initvalues = regionMap[rid].initvalues;
			Tn[i] = initvalues[0];
			Rn[i] = initvalues[1];
			diffCoef[i] = initvalues[2];
		}

	}

	for(i=0;i!=Ncel;++i){//init this according to materials
		//diffCoef[i] = 3.0;
		//Rn[i] = 1.0;
	}


	// change grid boundary to actual boundary
}

void HeadConductionSolver::SaveTransientOldData( )
{
	int i,j;
	if(      TimeScheme==1 ){  // Euler forwards
		for(i=0;i<Ncel;i++)
		{
			Tnp [i]= Tn[i];
		}
	}
	else if( TimeScheme==2 ){  // 2nd order BDF
		for(i=0;i<Ncel;i++)
		{
			//?? this may go wrong if compiler does not execute from right to left ??
			Tnp2[i]= Tnp[i]= Tn[i];
		}
	}
	else
	{
		cout<<"No unsteady time advanced scheme? Are you kidding?"<<endl;
		exit(0);
	}
}

void OutArray2File(double arr[],int N, ofstream &of)
{
	for(int i=0; i<N; i++ ){
		of<<arr[i]<<"  ";
		if( i%5==0 ) of<<endl;
	}
}
void HeadConductionSolver::Output2Tecplot()
{
	int i,j;
	double *tmp=NULL,nvar;
	ofstream of;
	char tecTitle[256];
	sprintf(tecTitle,"tec/res%04d.dat",this->outputCounter++);
	of.open(tecTitle);
	of<<"variables="<<"\"x\","<<"\"y\","<<"\"z\""
		<<"\"t\","<<endl;

	nvar = 4;
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
	OutArray2File( Tn,Ncel,of );

	for( i=0; i<Ncel; i++ ){
		for( j=0;j<8;j++ )
			of<<Cell[i].vertices[j]+1<<" ";
		of<<endl;
	}

	of.close();
	delete []tmp;
}

void outputVTKScalar( char name[], double arr[],int N, ofstream &of)
{
}
void HeadConductionSolver::Output2VTK( )
{

}

void HeadConductionSolver::WriteBackupFile( )
{

}
void HeadConductionSolver::ReadBackupFile( )
{

}

void HeadConductionSolver::OutputMoniter( )
{
    
}

HeadConductionSolver::~HeadConductionSolver()
{
	// output the result before error
	// Output2Tecplot ( );

	cout<<"desturct object and free space"<<endl;
	// delete variables
	delete [] Tn;
	delete_Array2D( dPhidX,Ncel,3 );

	delete [] BTem;


    V_Destr ( &bs );
	V_Destr ( &xsol );
}

