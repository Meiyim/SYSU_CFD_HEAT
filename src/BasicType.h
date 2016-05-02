#ifndef _BASIC_TYPE_H_
#define _BASIC_TYPE_H_

#include <stdlib.h>
#include <fstream>
#include "tools.h"
#include <string>
#include <map>
#include <vector>
#include <limits>
#include <assert.h>
extern "C"{ 
    #include "laspack/laspack.h"
}

#define SMALL 1.e-16
#define CYCASHUGE_D std::numeric_limits<double>::max()
#define CYCASHUGE_I std::numeric_limits<int>::max()

class FaceData{
public:
    int    bnd;
    int    vertices[4];
    int    cell1,cell2;
    double x[3],n[3],area;   // face center and normal vector
    double lambda,   // lambda for left, (1.-lambda) for right
        rlencos;     // area/(|Xpn|*vect_cosangle(n,Xpn))
    double Xpac[3],Xnac[3]; // auxiliary points
};

class CellData{
public:
    int nface;
    int rid;   //same as boundary data
    int face[6], cell[6], vertices[8]; // maybe wasterful a bit. Make it dynamics to save memory
    //face[nface]: for the faces index
    //  all elements are treated as 8 nodes hexahdron,
    //  but with different number of faces, so judges will be used for avoid faces with one vertex
    double vol, x[3];
};
// connect boundary to faces, boundary to BdRegion
class BoundaryData{
public:
    BoundaryData():
        face(-1),
        rid(-1),
        distance(0.0),
        yplus(0.0),
        uplus(0.0),
        h(0.0),q(0.0)
    {
        vertices[0]= vertices[1]=vertices[2]=vertices[3] = -1;
        shear[0]=shear[1]=shear[2]= 0.0;
    }
    int face;                // belongs to face...
    int vertices[4];         // the 4 vertices, to be done allocatable
    int rid;                 // region id as set in rtable
    double distance;         // normal distance from cell face
                             // center to cell center
    double  yplus;           // y+
    double  uplus;           // u+
    double  shear[3];        // shearstress components
    double  h;               // local heattransfer coef.
    double  q;               // local heat flux (in W/m2)
};
// store the data set in one type of boundary
class BdRegion{
public:
    std::string name;
    int  type1;  //=1: wall
                 //=2: inlet
                 //=3: outlet
                 //=4: symmetric
                 //=5: body
    // temperature
    int  type2; // for wall:(type1==1)
                // =0 : given temperature; 
                // =1 : heat flux (0 for default adiabatic, others )
                // =2 : coupled

                // for body:(type1==5)
                // =0 : fluid field: u, v, w, ro, p, t ...
                // =1 : solid: ro, diffCoef

    double fixedValue;
    double initvalues[10];//u,v,w,p,ro,t,te,ed

};

class CycasSolver{
public:
    double Residule;
public: 
    //parameter 
    map<int,BdRegion> regionMap; // record the region/bound info, Bnd[i].rid ---> regionMap
    //time scheme
    bool   IfSteady;
    double dt;
    int    TimeScheme;
        // TurModel=0: laminar flow; =1: k-epsilon; =2: wait to implement
        // DensityModel=0: constant density; =1: perfect gas ro=p/(RT); =2: wait to implement
        // TimeScheme=1, Euler (default) ;  =2, BDF 2nd
    int    limiter;  // Maximum step at outer iteration


    //geometry
    int          Nvrt, Ncel, Nfac, Nbnd;
    int          NCoupledBnd;
    double       **Vert; // coordinate x,y,z
    FaceData     *Face;
    CellData     *Cell;
    BoundaryData *Bnd;  

    //method
public:
    CycasSolver();
    virtual ~CycasSolver();
    //solver interface
    virtual int solve(){return 0;};
    virtual void Output2Tecplot(std::ofstream& ofile,int nvar){};
    virtual int CheckAndAllocate() {return 0;};
    virtual void InitFlowField() {};
    virtual void SaveTransientOldData() {};
    virtual double getResidule() {return 0.0;};



    // Basic algorithm
    // gradients and divergence
    int  Gradient    ( double*, double*, double** );
    // int  Divergence  ( double*, double*, double*, double*, double*, double*, double* );
    int  Limiter_MLP  (double[],double **);
    int  Limiter_Barth(double[],double **);
    int  Limiter_WENO (double[],double **);

};

class ErrorHandler
{
public:
    void fatalRuntimeError(const char* msg){
        printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        printf("%s\n",msg);
        printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        exit(-1);
    }        
    void fatalRuntimeError(const char* msg,const char* detail){
        printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        printf("%s : %s\n",msg,detail);
        printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        exit(-1);
    }
    void fatalRuntimeError(const char* msg,int code){
        printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        printf("%s : %d\n",msg,code);
        printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        exit(-1);
    }
};

extern ErrorHandler* errorHandler;

#endif