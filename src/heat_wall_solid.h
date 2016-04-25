#pragma once
#include <iostream>
#include "navier.h"

// two problems : 
// 1. wall heat transfer/condensation/vaporize, 
// 2. heat conduction in solids
class HeatTransfer   // : public HeadConductionSolver
{
protected:
    int    ;
    double ;
    FaceData  *face;
    CellData  *cell;
	
	const HeadConductionSolver &NS; // cannot modify everything of NS
    
public:
    double ConvectHeatTransferCoef(double Re, double Pr);
};

void HeatTransfer::Convect()
{
	NS.Un[i];
}


/***************************************************
1D heat conduct equation solving
   Cp * partial(T)/partial(t) = romda * partial^2(T)/partial(x)^2
   partial(T(t,0))  /partial(x)= -qL/romda,
   partial(T(t,Len))/partial(x)= -qR/romda,
   initial condition : T(0,x) is given through STem[]
return T(t+dt) using STem[].  n is the cell number in solving
***************************************************/
int  HeatConduct1D( int n, real STem[], real dt, real Len, real Cp, real romda, real qL, real qR );