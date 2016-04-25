#pragma once

class Material
{
public:
    Material(){ order=0; };
    ~Material(){ };
    // iso-volume   specific heat capacity
    double CalCv (double Tem);
	double CalVis(double Tem);
	double CalTherCoef(double Tem);
    // diffusion coefficients
    double CalBinDiff(int type, double Tem);

    // calculation species and their properties
    //   According to input-card, initialize these variables from the database
    char   name[50];
    int    order; // order<=5
    double CvCoef[5], VisCoef[5], ThermalCoef[5], *DijCoef[5];
};

extern class Material *MatSp; // should be an array stored the 

int MaterialInit( );
