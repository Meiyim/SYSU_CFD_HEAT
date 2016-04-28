#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include "navier.h"
#include "tools.h"
#include "terminalPrinter.h"

using namespace std;


void NavierStokesSolver::NSSolve( )
{
	int iter;
	double ResMax,ResMax0;
	//timespec tstart ={0,0};
	//timespec tend ={0,0};
	
	cur_time = 0.0;
	if(IfSteady){
		dt = CYCASHUGE_D;
		MaxOuterStep = MaxStep;
	}
	cur_time += dt;
	//MaxOuterStep=2;//test
   	for( step=1; step<=MaxStep; step++ ) //step : total time step
   	{
		if( !IfSteady ){
			SaveTransientOldData( );
		}
	// outer iteration, drive residual to zero for every physical time step
		for( iter=1; iter<MaxOuterStep; iter++ ){
			if( (iter-1)%10==0 )printer->printSectionHead();

			CalculateVelocity ( ); //interface communication U, V, W, Apr
			
			CalculatePressure ( ); //calculate deltaP and correct P,[R], U,V,W

			//------------   Physics Models   ------------//
			// scalar transportation
			//1. turbulence model
			if(TurModel==1) {
				UpdateTurKEpsilon( );
			}
			//2. energy couple
			if( SolveEnergy  ) {
				UpdateEnergy ( );
			}
			//3. species transport
			if( SolveSpecies ) UpdateSpecies( );//to be implemented
			//4. other physical models goes here
			//    e.g., condensation, combustion, wall heat conduction
			//------------   Record   ------------//
			if( shouldPostProcess(step,iter,cur_time) ){
				Output();	
			}
			/*-----------check if should break----------*/
			ResMax = vec_max( Residual,10);
			
			if( IfSteady ){
				//steady
				printer->printSteadyStatus(iter,ResMax);
				if( ResMax<ResidualSteady )break;
			}else{
				//unsteady
				if( iter == 1 ) ResMax0 = ResMax;
				if( ResMax<1.e-4 || ResMax0/(ResMax+1.e-16)>1000. ){
					printer->printStepStatus(step,iter,cur_time,dt,ResMax);
					break; // more reasonal to break : order drop 2
				}
			}

		}
		
		/*-----------record,tot file, restart file, etc.----------*/

		if(cur_time >= total_time){
			break;
		}else{//time advance
			cur_time+=dt;
		}	
	}

	//extra work before solve is complete 
	Output();

	printer->printEnding();

	return;
}


/*******************************************************/
//	determint if should backUp
//	currently the same frequency as post process
/*******************************************************/
bool NavierStokesSolver::shouldPostProcess(int timestep,int outiter, double now){
	if(IfSteady){
		return (outiter-1)%noutput == 0;
	}else{
		return false;
	}
}


void NavierStokesSolver::InitFlowField( )
{
	int i;

	for( i=0; i<Ncel; i++ )
	{
		int rid = Cell[i].rid;
		assert(regionMap[rid].type1==5&&regionMap[rid].type2==0);//fluid
		double* initvalues = regionMap[rid].initvalues;
		Un[i] = initvalues[0];
		Vn[i] = initvalues[1];
		Wn[i] = initvalues[2];
		Rn[i] = initvalues[3];
		Tn[i] = initvalues[4];
		if( DensityModel== 1 ) Rn[i]= (Pn[i]+PressureReference)/(Rcpcv*Tn[i]);

		VisLam[i]= initvalues[5]; // 0.6666667e-2;  // 1.458e-6 * pow(Tn[i],1.5) /(Tn[i]+110.5) ;
		VisTur[i]= 0.;
		if( TurModel==1 )
		{
			TE[i]    = initvalues[6];  // 1.e-4*(Un[i]*Un[i]+Vn[i]*Vn[i]+Wn[i]*Wn[i]);
			ED[i]    = initvalues[7];    // TurKEpsilonVar::Cmu * pow(TE[i],1.5) / 1.;
			VisTur[i]= Rn[i]*TurKEpsilonVar::Cmu * TE[i]*TE[i]/(ED[i]+SMALL);
		}
	}

	for( i=0;i<Nfac;i++ )
		RUFace[i] = 0.;

	//init physical module
	for(map<string,CycasSolver*>::iterator iter = physicalModule.begin();
		iter!=physicalModule.end(); ++iter){
		iter->second->InitFlowField();
	}
}

void NavierStokesSolver::SaveTransientOldData( )
{
	int i,j;
	if(      TimeScheme==1 ){  // Euler forwards
		for(i=0;i<Ncel;i++)
		{
			Rnp [i]= Rn[i];
			Unp [i]= Un[i];
			Vnp [i]= Vn[i];
			Wnp [i]= Wn[i];
			Tnp [i]= Tn[i];
			TEp [i]= TE[i];
			EDp [i]= ED[i];
			for( j=0; j<Nspecies; j++ )
				RSnp[i][j]= RSn[i][j];
			if( TurModel==1 )
			{
				TEp[i]= TE[i];
				EDp[i]= ED[i];
			}
		}
	}
	else if( TimeScheme==2 ){  // 2nd order BDF
		for(i=0;i<Ncel;i++)
		{
			//?? this may go wrong if compiler does not execute from right to left ??
			Rnp2[i]= Rnp[i]= Rn[i];
			Unp2[i]= Unp[i]= Un[i];
			Vnp2[i]= Vnp[i]= Vn[i];
			Wnp2[i]= Wnp[i]= Wn[i];
			Tnp2[i]= Tnp[i]= Tn[i];
			TEp2[i]= TEp[i]= TE[i];
			EDp2[i]= EDp[i]= ED[i];
			for( j=0; j<Nspecies; j++ )
				RSnp2[i][j]= RSnp[i][j]= RSn[i][j];
			if( TurModel==1 )
			{
				TEp2[i]= TEp[i]= TE[i];
				EDp2[i]= EDp[i]= ED[i];
			}
		}
	}
	else
	{
		errorHandler->fatalRuntimeError("No unsteady time advanced scheme? Are you kidding?");
	}

	for(map<string,CycasSolver*>::iterator iter = physicalModule.begin();
		iter!=physicalModule.end(); ++iter){
		iter->second->SaveTransientOldData();
	}
}



void NavierStokesSolver::Output(){
 	ofstream of;
 	int nvar;

	char tecTitle[256];
	sprintf(tecTitle,"tec/res%04lu.dat",this->outputCounter++);

	of.open(tecTitle);	
	of<<"variables="<<"\"x\","<<"\"y\","<<"\"z\""
		<<"\"p\","<<"\"u\","<<"\"v\","<<"\"w\","<<"\"ro\","<<"\"T\","
		<<"\"Mach\","<<"\"mu\"";
	nvar = 11;
	if( TurModel==1 ){
		of<<"\"te\""<<"\"ed\"";
		nvar = 13;
	}

	Output2Tecplot(of,nvar);

	// other module...
	for(map<string,CycasSolver*>::iterator iter = physicalModule.begin();
		iter!=physicalModule.end(); ++iter){
		iter->second->Output2Tecplot(of,nvar);
	}

	of.close();
}

void NavierStokesSolver::Output2Tecplot(std::ofstream& of,int nvar)
{
	int i,j;
	double *tmp;


	of<<endl;
	of<<"zone n="<<Nvrt<<", e="<<Ncel<<", VARLOCATION=([1-3]=NODAL,[4-"<<nvar<<"]=CELLCENTERED)"
		<<"DATAPACKING=BLOCK, ZONETYPE=FEBRICK"
		<<endl;
	for( j=0; j<3; j++ ){
		for( i=0; i<Nvrt; i++ ){
			of<<Vert[i][j]<<" ";
			if( i%5==4 ) of<<endl;
		}
	}
	OutArray2File( Pn,Ncel,of );
	OutArray2File( Un,Ncel,of );
	OutArray2File( Vn,Ncel,of );
	OutArray2File( Wn,Ncel,of );
	OutArray2File( Rn,Ncel,of );
	OutArray2File( Tn,Ncel,of );
	tmp = new double[Ncel];
	for( i=0; i<Ncel; i++ ){  // Mach / velocity magnitude
		if( DensityModel==1 )
			tmp[i]= sqrt( (Un[i]*Un[i]+Vn[i]*Vn[i]+Wn[i]*Wn[i])/(gama*(Pn[i]+PressureReference)/Rn[i]) );
		else
			tmp[i]= sqrt( (Un[i]*Un[i]+Vn[i]*Vn[i]+Wn[i]*Wn[i]) ); 
	}
	OutArray2File( tmp,Ncel,of );
	for( i=0; i<Ncel; i++ ){  // mu
		tmp[i]= VisLam[i] + VisTur[i];
	}
	OutArray2File( tmp,Ncel,of );
	if( TurModel==1 ){
		OutArray2File( TE,Ncel,of );
		OutArray2File( ED,Ncel,of );
	}

	for( i=0; i<Ncel; i++ ){
		for( j=0;j<8;j++ )
			of<<Cell[i].vertices[j]+1<<" ";
		of<<endl;
	}




	delete []tmp;
}


NavierStokesSolver::NavierStokesSolver():
	outputCounter(0),
	printer(new TerminalPrinter),
	TurModel(0),
	DensityModel(0),
	SolveEnergy(false),
	SolveSpecies(false),
	Nspecies(0),
	PressureReference(1.01325e5),
	cellPressureRef(0),
	
    MaxOuterStep(50),
    MaxStep(10000),
    cur_time(0.0),
    total_time(0.0),
    ResidualSteady(1.e-6),
    noutput(1),
    outputformat(0),

	gama(1.4),
	ga1(0.4),
	cp(1006.),
	cv(cp/gama),
	prl(0.72),
	prte(0.9),
	Rcpcv(cp-cv),
	TempRef(273.15),
	//all put NULL to avoid wild pointer
	Rn(NULL),Un(NULL),Vn(NULL),Wn(NULL),Tn(NULL),TE(NULL),ED(NULL),
	RSn(NULL),

	Pn(NULL),
	
	Rnp(NULL),Unp(NULL),Vnp(NULL),Wnp(NULL),Tnp(NULL),TEp(NULL),EDp(NULL),RSnp(NULL),
	Rnp2(NULL),Unp2(NULL),Vnp2(NULL),Wnp2(NULL),Tnp2(NULL),TEp2(NULL),EDp2(NULL),RSnp2(NULL),

	VisLam(NULL),VisTur(NULL),
	dPdX(NULL),dUdX(NULL),dVdX(NULL),dWdX(NULL),Apr(NULL),dPhidX(NULL),

	RUFace(NULL),

	BRo(NULL),BU(NULL),BV(NULL),BW(NULL),BPre(NULL),BTem(NULL),BRS(NULL),BTE(NULL),BED(NULL)
{
	for( int i=0; i<3; i++ )
		gravity[i] = 0.;
	//-- numerical scheme
	// relaxation factor
	URF[0]= 0.6;  // u
	URF[1]= 0.6;  // v
	URF[2]= 0.6;  // w
	URF[3]= 0.5;  // p
	URF[4]= 0.8;  // T
	URF[5]= 0.8;  // k
	URF[6]= 0.8;  // e
	URF[7]= 0.8;  // scalar
}


NavierStokesSolver::~NavierStokesSolver()
{
	// output the result before error
	// Output2Tecplot ( );

	cout<<"desturct object and free space"<<endl;
	// delete variables
	for(map<string,CycasSolver*>::iterator iter = physicalModule.begin();
		iter!=physicalModule.end(); ++iter){
			delete iter->second;
	}
   	delete [] Rn;
	delete [] Un;
	delete [] Vn;
	delete [] Wn;
	delete [] Pn;
	delete [] Tn;
	delete [] TE;
	delete [] ED;
	delete_Array2D(RSn,Nspecies,Ncel);

	delete [] Rnp;
	delete [] Unp;
	delete [] Vnp;
	delete [] Wnp;
	delete [] Tnp;
	delete [] TEp;
	delete [] EDp;
	delete_Array2D(RSnp,Nspecies,Ncel);
	delete [] Rnp2;
	delete [] Unp2;
	delete [] Vnp2;
	delete [] Wnp2;
	delete [] Tnp2;
	delete [] TEp2;
	delete [] EDp2;
	delete_Array2D(RSnp2,Nspecies,Ncel);

	delete [] VisLam;
	delete [] VisTur;
	delete_Array2D( dPdX,Ncel,3 );
	delete_Array2D( dUdX,Ncel,3 );
	delete_Array2D( dVdX,Ncel,3 );
	delete_Array2D( dWdX,Ncel,3 );
	delete [] Apr;
	delete_Array2D( dPhidX,Ncel,3 );

	delete [] BRo ;
	delete [] BU  ;
	delete [] BV  ;
	delete [] BW  ;
	delete [] BTem;
	delete [] BPre;


	delete_Array2D(BRS,Nspecies,Nbnd);
	delete [] BTE;
	delete [] BED;
	delete [] RUFace;

    V_Destr ( &bs );
	V_Destr ( &bu );
	V_Destr ( &bv );
	V_Destr ( &bw );
	V_Destr ( &bp );
	V_Destr ( &xsol );
}


int NavierStokesSolver::CheckAndAllocate()
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
   	Rn = new double[Ncel];
	Un = new double[Ncel];
	Vn = new double[Ncel];
	Wn = new double[Ncel];
	Pn = new double[Ncel];
	Tn = new double[Ncel];
	TE = new double[Ncel];
	ED = new double[Ncel];
	RSn= new_Array2D<double>(Nspecies,Ncel);

	if( !IfSteady ){
	if( TimeScheme>=1 ){
	Rnp = new double[Ncel];
	Unp = new double[Ncel];
	Vnp = new double[Ncel];
	Wnp = new double[Ncel];
	Tnp = new double[Ncel];
	TEp = new double[Ncel];
	EDp = new double[Ncel];
	RSnp= new_Array2D<double>(Nspecies,Ncel);
	}

	if( TimeScheme>=2 ){
	Rnp2 = new double[Ncel];
	Unp2 = new double[Ncel];
	Vnp2 = new double[Ncel];
	Wnp2 = new double[Ncel];
	Tnp2 = new double[Ncel];
	TEp2 = new double[Ncel];
	EDp2 = new double[Ncel];
	RSnp2= new_Array2D<double>(Nspecies,Ncel);
	}
	}

	VisLam = new double[Ncel];
	VisTur = new double[Ncel];
	dPdX   = new_Array2D<double>(Ncel,3);
	dUdX   = new_Array2D<double>(Ncel,3);
	dVdX   = new_Array2D<double>(Ncel,3);
	dWdX   = new_Array2D<double>(Ncel,3);
	Apr    = new double[Ncel];
	dPhidX = new_Array2D<double>(Ncel,3);

	BRo = new double[Nbnd];
	BU  = new double[Nbnd];
	BV  = new double[Nbnd];
	BW  = new double[Nbnd];
	BTem= new double[Nbnd];
	BPre= new double[Nbnd];
	BRS = new_Array2D<double>(Nspecies,Nbnd);
	BTE = new double[Nbnd];
	BED = new double[Nbnd];

	RUFace = new double[Nfac];

	
	// laspack working array
	V_Constr(&bs,   "rightU",    Ncel, Normal, True);
	V_Constr(&bu,   "rightU",    Ncel, Normal, True);
	V_Constr(&bv,   "rightU",    Ncel, Normal, True);
	V_Constr(&bw,   "rightU",    Ncel, Normal, True);
	V_Constr(&bp,   "rightU",    Ncel, Normal, True);
	V_Constr(&xsol, "rightU",    Ncel, Normal, True);

	cur_time = 0.;

	for(map<string,CycasSolver*>::iterator iter = physicalModule.begin();
		iter!=physicalModule.end(); ++iter){
		iter->second->CheckAndAllocate();
	}

	return 1;
}

