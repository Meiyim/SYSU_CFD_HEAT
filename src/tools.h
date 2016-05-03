#ifndef  TOOLS_H
#define  TOOLS_H

#include <string>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <fstream>
extern "C"{
#include "laspack/laspack.h"
}
using namespace std;

#define CYCASMAX(x,y)  ((x)>(y)?(x):(y))
#define CYCASMIN(x,y)  ((x)<(y)?(x):(y))
#define CYCASSIGN(x)   ((x)>0?1:(-1))
// maybe SIGN still has problem, not for 0

// vector manipulation. should be defined as inline function for higher efficiency
void   vec_init  (double[], int, double );
void   vec_minus (double *x1, double *x2, double *x3, int n); // x1= x2 - x3
double vec_dot   (double[], double[], int);
double vec_len   (double[], int );
void   vec_cross (double[], double[], double[]); // only for C[3]= A[3] x B[3];
double vec_max   (double[], int );


void SolveLinearEqu( Vector* Func(QMatrix*, Vector*, Vector*, int,PrecondProcType, double),
			QMatrix *qa, Vector *x, Vector *b, int MaxIter, PrecondProcType PreCond, double omega,
					   double epsilon, int *Iter, double *IterRes);

void OutPlaceHolder2File(double placeHolder, int N, ofstream& of);
void OutArray2File(double arr[],int N, ofstream &of);
char *trimwhitespace(char *str);
double ttime (void);

template <typename T>
T** new_Array2D(int row, int col)
{
    int size = sizeof(T);
    int point_size = sizeof(T*);
    T **arr = (T **) malloc(point_size * row + size * row * col);
    if (arr != NULL)
    {   
        T *head = (T*)(arr + row); //made some change
        for (int i = 0; i < row; ++i)
        {
            arr[i] = &head[i*col]; //made some change
            for (int j = 0; j < col; ++j)
                new (&arr[i][j]) T;
        }
    }
    return (T**)arr;
}
template <typename T>
void delete_Array2D(T **arr, int row, int col)
{
    for (int i = 0; i < row; ++i)
        for (int j = 0; j < col; ++j)
            arr[i][j].~T();
    if (arr != NULL)
        free((void**)arr);
}
template <typename T>
void init_Array2D(T **arr, int row, int col, T val )
{
    for (int i = 0; i < row; ++i)
        for (int j = 0; j < col; ++j)
            arr[i][j] = val;
}


class Checker{
	string varname;
	vector<double> pool;
public:
	Checker(string n):varname(n){}
	
	void check(double v){
		pool.push_back(v);
	}

	void report(){
		double ret = 0.0;
		for(auto it=pool.begin();it!=pool.end();++it)
			ret += (*it);
		printf("%s:\t%e\tcount%lu\n",varname.c_str(),ret,pool.size());
	}
};

#define CHECK_ARRAY(arr,len) checkArray(arr,len,#arr)
#define CHECK_MEMBER_ARRAY(arr,mem,len) checkMemberArray(arr,len,mem,#arr #mem)

template<typename T>
void checkArray(T* arr, size_t len,const char* name){
	double ret = 0.0;
	for(int i=0;i!=len;++i){
		ret += (arr[i])*(arr[i]);
	}
	
	char temp[256];
	sprintf(temp,"%s: norm %30.9e\n",name,sqrt(ret));
	printf("%s",temp);
}
template<typename T>
void checkMemberArray(T* arr, size_t len, double T::* m_ptr,const char* name){
	double ret = 0.0;
	for(int i=0;i!=len;++i){
		ret += (arr[i].*m_ptr) * (arr[i].*m_ptr);
	}

	
	char temp[256];
	sprintf(temp,"%s: norm %30.9e\n",name,sqrt(ret));
	printf("%s",temp);
}


#endif
