#include <math.h>
#include <iostream>
#include <fstream>
#include "navier.h"

using namespace std;

void HeadConductionSolver::SetBCTemperature( double *bt )
{
	int    i,rid,iface,ic;
	for( i=0; i<Nbnd; i++ )
	{
		rid   = Bnd[i].rid;
		iface = Bnd[i].face;
		ic    = Face[iface].cell1;
		switch( regionMap[rid].type1 ){
		case(1):  // wall
			if(regionMap[rid].type2==0){
				bt[i] = regionMap[rid].fixedValue;//given T
			}else if(regionMap[rid].type2==1){
			   	//bt[i] = Tn[ic] - regionMap[rid].fixedValue / (diffCoef[ic] *Face[iface].rlencos );// given flux // by CHENXUYI
			   	bt[i] = 0.0;// bt will not be used if 2nd kind of boundary type
			}else{
				assert(false);
			}
			break;
		case(4):
			bt[i]=Tn[ic];//adiabatic method //this value wont be used
			break;
		default:
			cout<<"no such boundary type"<<i<<" "<<rid<<endl;
			exit(0);
		}
	}
}

