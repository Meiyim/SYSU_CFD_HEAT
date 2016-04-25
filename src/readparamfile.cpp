#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "navier.h"

void HeadConductionSolver::ReadParamFile( )
{
	int i,j, ikey;
	ifstream fin;
	char      *keyw[20];
	std::string str;

	for( i=0; i<20; i++ )
		keyw[i] = new char[30];

	fin.open("param.in");
	if( ! fin.is_open( ) ){
		cout<<"param file does not exist!"<<endl;
		exit(0);
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
			IfSteady = true;
			MaxStep  = atoi(keyw[1]);
			ResidualSteady= atof(keyw[2]);
		}
		else if( strcmp(keyw[0],"transient")==0 )
		{
			IfSteady = false;
			dt = atof(keyw[1]);
			total_time= atof(keyw[2]);
			if(      strcmp(keyw[3],"Euler")==0 )
				TimeScheme = 1;
			else if( strcmp(keyw[3],"Dual")==0 )
				TimeScheme = 2;
			else
			{
				cout<<"no such unsteady time scheme!:"<<keyw[3]<<endl;
				exit(0);
			}
			if( ikey>=4 ) MaxOuterStep = atoi(keyw[4]);
		}


		//--- numerical scheme
		else if( strcmp(keyw[0],"relaxation")==0 )
		{
			URF = atof( keyw[1] );
		}
		else if( strcmp(keyw[0],"limiter")==0 )
		{
			if(      strcmp(keyw[1],"no")==0 )
				limiter = 0;
			else if( strcmp(keyw[1],"Barth" )==0 )
				limiter = 1;
			else if( strcmp(keyw[1],"MLP" )==0 )
				limiter = 2;
			else if( strcmp(keyw[1],"WENO" )==0 )
				limiter = 3;
			else
				ErrorStop("error in limiter definition");
		}

		//--- boundary condition
		else if( strcmp(keyw[0],"bound")==0 )
		{
			int bid= atoi( keyw[1] );
			if(regionMap.find(bid)!=regionMap.end()){
				cout<<"repeating bid found"<<endl;
			}
			if(      strcmp(keyw[2],"Twall")==0 )
			{
				regionMap[bid].name="Twall";
				regionMap[bid].type1 = 1;
				regionMap[bid].type2 = 0;//given T
				regionMap[bid].fixedValue = atof(keyw[3]);
			}
			else if( strcmp(keyw[2],"Hwall")==0 ){
				regionMap[bid].name="Hwall";
				regionMap[bid].type1 = 1;
				regionMap[bid].type2 = 1;//given flux
				regionMap[bid].fixedValue = atof(keyw[3]);
			}
			else if( strcmp(keyw[2],"sym")==0 ){
				regionMap[bid].name="sym";
				regionMap[bid].type1 = 4;
				regionMap[bid].type2 = 0;
				regionMap[bid].fixedValue = 0.0;
			}
			else
			{
				cout<<"unknown boundry type in bound"<<keyw[2]<<endl;
				exit(0);
			}
		}
		else if( strcmp(keyw[0],"volumn")==0 )
		{
			int bid= atoi( keyw[1] );
			if(regionMap.find(bid)!=regionMap.end()){
				cout<<"repeating bid found"<<endl;
			}
			if(      strcmp(keyw[2],"fluid")==0 )//fluid
			{
				regionMap[bid].name= keyw[3];
				regionMap[bid].type1 = 5;
				regionMap[bid].type2 = 0;//fluid
				for(int i=4;i<=ikey;++i){
					regionMap[bid].initvalues[i-4]= atof( keyw[i] );
				}
			}
			else if( strcmp(keyw[2],"solid")==0 ){// solid
				regionMap[bid].name=keyw[3];
				regionMap[bid].type1 = 5;
				regionMap[bid].type2 = 1;//solid
				for(int i=4;i<=ikey;++i){
					regionMap[bid].initvalues[i-4]= atof( keyw[i] );
				}
			}

			else
			{
				cout<<"unknown boundry type in bound"<<keyw[2]<<endl;
				exit(0);
			}
		}

		//--- post process 
		else if( strcmp(keyw[0],"output")==0 )
		{
			noutput = atoi( keyw[1] );
			if(      strcmp(keyw[2],"vtk")==0 )
				outputFormat = 1;
			else if( strcmp(keyw[2],"tecplot")==0 )
				outputFormat = 0;
			else
				outputFormat = 0;
		}
		else
		{
			cout<<"no such command, "<<keyw[0]<<endl;
			exit(0);
		}

	}

	fin.close( );
}
