#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "PreProcessor.h"
#include "navier.h"

void PreProcessor::ReadParamFile(NavierStokesSolver* ns)
{
	int i,j, ikey;
	ifstream fin;
	char      *keyw[20];
	std::string str;

	for( i=0; i<20; i++ )
		keyw[i] = new char[30];

	fin.open(paramFileName.c_str());
	if( ! fin ){
		errorHandler->fatalRuntimeError("param file can not found",paramFileName.c_str());
	}
	while( getline( fin,str ) )
	{
		// skip comment line
		if( str[0]=='/' && str[1]=='/' ) continue;
		if( str.size()==1 )continue; //empty line with '\n'
		// clear keyw[] before putting value
		for( i=0; i<20; i++ )
			keyw[i][0] = '\0';
		// put commands+parameters into string array
		j=0; ikey=0;
		for( i=0; str[i]!='\0'; i++ )
		{
			if( str[i]!=',' )
			{
				keyw[ikey][j]= str[i];
				j++;
			}
			else
			{
				keyw[ikey][j]='\0';
				ikey++;
				j=0;
			}
		}
		keyw[ikey][j]='\0';
		
		// trim spaces in keywords
		for( i=0; i<=ikey; i++ )
		{
			keyw[i] = trimwhitespace( keyw[i] );
		}


		// analyze command and parameters
		//---- solver parameters
		if(      strcmp(keyw[0],"gridfile")==0 )
		{
			strcpy( GridFileName, keyw[1] );
		}
		else if( strcmp(keyw[0],"steady")==0 )
		{
			ns->IfSteady = true;
			ns->MaxStep  = atoi(keyw[1]);
			ns->ResidualSteady= atof(keyw[2]);
		}
		else if( strcmp(keyw[0],"transient")==0 )
		{
			ns->IfSteady = false;
			ns->dt = atof(keyw[1]);
			ns->total_time= atof(keyw[2]);
			if(      strcmp(keyw[3],"Euler")==0 )
				ns->TimeScheme = 1;
			else if( strcmp(keyw[3],"Dual")==0 )
				ns->TimeScheme = 2;
			else
			{
				errorHandler->fatalRuntimeError("no such unsteady time scheme!");
			}
			if( ikey>=4 ) ns->MaxOuterStep = atoi(keyw[4]);
		}
		else if( strcmp(keyw[0],"energy")==0 )
		{
			if(      strcmp(keyw[1],"on")==0 )
				ns->SolveEnergy= true;
			else if( strcmp(keyw[1],"off")==0 )
				ns->SolveEnergy= false;
			else
				errorHandler->fatalRuntimeError("energy key word is wrong");
		}
		else if( strcmp(keyw[0],"density")==0 )
		{
			ns->DensityModel = atoi(keyw[1]);
		}
		else if( strcmp(keyw[0],"turbulence")==0 )
		{
			if(      strcmp(keyw[1],"no")==0 )
				ns->TurModel= 0;
			else if( strcmp(keyw[1],"ke")==0 )
				ns->TurModel= 1;
			else
				errorHandler->fatalRuntimeError("turbulence key word is wrong");
		}
		else if( strcmp(keyw[0],"gravity")==0 )
		{
			ns->gravity[0] = atof( keyw[1] );
			ns->gravity[1] = atof( keyw[2] );
			ns->gravity[2] = atof( keyw[3] );
		}
		else if( strcmp(keyw[0],"PressureRef")==0 )
		{
			ns->PressureReference = atof(keyw[1]);
			if(strlen(keyw[2])>0) ns->cellPressureRef   = atoi(keyw[2])-1;
		}

		//--- numerical scheme
		else if( strcmp(keyw[0],"relaxation")==0 )
		{
			ns->URF[0] = ns->URF[1] = ns-> URF[2] = atof( keyw[1] );
			ns->URF[3] = atof( keyw[2] );
			ns->URF[4] = atof( keyw[3] );
		}
		else if( strcmp(keyw[0],"limiter")==0 )
		{
			if(      strcmp(keyw[1],"no")==0 )
				ns->limiter = 0;
			else if( strcmp(keyw[1],"Barth" )==0 )
				ns->limiter = 1;
			else if( strcmp(keyw[1],"MLP" )==0 )
				ns->limiter = 2;
			else if( strcmp(keyw[1],"WENO" )==0 )
				ns->limiter = 3;
			else
				errorHandler->fatalRuntimeError("error in limiter definition");
		}
		else if( strcmp(keyw[0],"initflow")==0 )
		{
			for( i=1; i<=ikey; i++ )
				ns->initvalues[i-1]= atof( keyw[i] );
		}


		//--- boundary condition
		else if( strcmp(keyw[0],"bound")==0 )
		{
			int bid= atoi( keyw[1] );
			if(ns->regionMap.find(bid)!=ns->regionMap.end()){
				errorHandler->fatalRuntimeError("repeating bid found",bid);
			}
			if(      strcmp(keyw[2],"Twall")==0 )
			{
				ns->regionMap[bid].name="Twall";
				ns->regionMap[bid].type1 = 1;
				ns->regionMap[bid].type2 = 0;//given T
				ns->regionMap[bid].fixedValue = atof(keyw[3]);
			}
			else if( strcmp(keyw[2],"Hwall")==0 ){
				ns->regionMap[bid].name="Hwall";
				ns->regionMap[bid].type1 = 1;
				ns->regionMap[bid].type2 = 1;//given flux
				ns->regionMap[bid].fixedValue = atof(keyw[3]);
			}
			else if( strcmp(keyw[2],"sym")==0 ){
				ns->regionMap[bid].name="sym";
				ns->regionMap[bid].type1 = 4;
				ns->regionMap[bid].type2 = 0;
				ns->regionMap[bid].fixedValue = 0.0;
			}
			else if( strcmp(keyw[2],"inlet")==0 ){
				ns->regionMap[bid].name="inlet";
				ns->regionMap[bid].type1 = 2;
				ns->regionMap[bid].type2 = 0;
				double* initvalues = ns->regionMap[bid].initvalues;
				initvalues[0] = atof(keyw[3]);//u
				initvalues[1] = atof(keyw[4]);//v
				initvalues[2] = atof(keyw[5]);//w
				initvalues[3] = atof(keyw[6]);//p
				initvalues[4] = atof(keyw[7]);//r
				initvalues[5] = atof(keyw[8]);//t
				if( ns->TurModel==1 ){
					initvalues[6] = atof(keyw[9]);//te
					initvalues[7] = atof(keyw[10]);//ed
				}
			}	
			else if( strcmp(keyw[2],"outlet")==0 ){
				ns->regionMap[bid].name="outlet";
				ns->regionMap[bid].type1 = 3;
				ns->regionMap[bid].type2 = 0;
				ns->regionMap[bid].fixedValue = atof(keyw[3]); //pout
			}
			else
			{
				errorHandler->fatalRuntimeError("unknown boundry type in bound",keyw[2]);
				exit(0);
			}
		}
		else if( strcmp(keyw[0],"volumn")==0 )
		{
			int bid= atoi( keyw[1] );
			if(ns->regionMap.find(bid)!=ns->regionMap.end()){
				cout<<"repeating bid found"<<endl;
			}
			if(      strcmp(keyw[2],"fluid")==0 )//fluid
			{
				ns->regionMap[bid].name= keyw[3];
				ns->regionMap[bid].type1 = 5;
				ns->regionMap[bid].type2 = 0;//fluid
				for(int i=4;i<=ikey;++i){
					ns->regionMap[bid].initvalues[i-4]= atof( keyw[i] );
				}
			}
			else if( strcmp(keyw[2],"solid")==0 ){// solid
				ns->regionMap[bid].name=keyw[3];
				ns->regionMap[bid].type1 = 5;
				ns->regionMap[bid].type2 = 1;//solid
				for(int i=4;i<=ikey;++i){
					ns->regionMap[bid].initvalues[i-4]= atof( keyw[i] );
				}
			}

			else
			{
				errorHandler->fatalRuntimeError("unknown boundry type in bound",keyw[2]);
			}
		}

		//--- post process 
		else if( strcmp(keyw[0],"output")==0 )
		{
			ns->noutput = atoi( keyw[1] );
			if(      strcmp(keyw[2],"vtk")==0 )
				ns->outputformat = 1;
			else if( strcmp(keyw[2],"tecplot")==0 )
				ns->outputformat = 0;
			else
				ns->outputformat = 0;
		}
		else if(strcmp(keyw[0],"3dheatconduction")==0){
			solve3DHeatConduction = true;	
			ns->Solve3DHeatConduction = true;
		}
		else
		{
			errorHandler->fatalRuntimeError("no such command: ",keyw[0]);
		}

	}

	fin.close( );
}

void PreProcessor::parseCommandLine(int argc, char** argv){
	for(int i=1;i<argc;++i){
		string command(argv[i]);
		if((++i)>argc){
			break;
		}
		string value(argv[i]);
		if(command=="-paramfile"){
			paramFileName = value;
		}
	}		
}
