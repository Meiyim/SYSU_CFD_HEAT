#ifndef TERMINAL_PRINTER_H
#define TERMINAL_PRINTER_H
#include <stdlib.h>
#include <fstream>
#include <string>
#include <deque>
#include <queue>
#include <vector>
using namespace std;


class TerminalPrinter{
public:
	void printStarter(){
		printf( "    1111111    111        111    1111111         111111          111111          \n");
		printf( "  11111111111  111        111  11111111111      11111111      11111111111        \n");
		printf( " 1111     1111 111        111 1111     1111    1111  1111    1111     1111       \n");
		printf( " 1111            111    111   1111            1111    1111   1111     1111       \n");
		printf( " 1111              111111     1111           1111      1111    111               \n");
		printf( " 1111               1111      1111           1111      1111      11111           \n");
		printf( " 1111               1111      1111           11111111111111         1111         \n");
		printf( " 1111               1111      1111           1111      1111  1111     1111       \n");
		printf( " 1111     1111      1111       1111    1111  1111      1111  1111     1111       \n");
		printf( "  11111111111       1111        1111111111   1111      1111   11111111111        \n");
	    printf( "    1111111         1111          111111     1111      1111      111111          \n");

	    printf("\n\n");
	    printf( " #    # # ####  ##### ###### #####  #    #  ##### ###### #    # #####  ######    \n");
	    printf( " #    # ##   # #        ##   #    # #    # #        ##   #    # #    # #         \n");
	    printf( " #    # #    #  ####    ##   #####  #    # #        ##   #    # #####  ######    \n");
	    printf( " #    # #    #      #   ##   #   #  #    # #        ##   #    # #   #  #         \n");
	    printf( "  ####  #    # #####    ##   #    #  ####   #####   ##    ####  #    # ######    \n");
		printf("                   *******************************************\n");
		printf("                     COMPUTATIONAL CODES FOR FLUID DYNAMICS\n");
		printf("                   *******************************************\n");
		printf("                   VERSION 2.0                     27 NOV 2015 \n");
		printf("                             EXECUTABLE ATTRIBUTES \n\n");
		printf("   All rights reserved. Unauthorized use, distribution or duplication is \n");
		printf("   prohibited. This product is subject to China Nuclear Power Technology \n");
		printf("   Research Institute and State Nuclear Power Software Development Center.\n");
		printf("--------------------------------------------------------------------------------\n\n");
	}
	void printStepStatus(int step,int piter ,double time,double dt,double res){
		printf("%15f\t%10d\t%10d\t%13.5f\t%13.5f\n",time,step,piter,dt,res);
	}
	void printSteadyStatus(int step,double res){
		printf("%15s\t%10d\t%10s\t%13s\t%15.10f\n","---",step,"---","---",res);
	}
	void printSectionHead(){
		printf("%15s\t%10s\t%10s\t%13s\t%15s\n","TIME","CAL STEPE","ITER","DELT","MAX RES");
	}
	void printEnding(){
		printf( " -------------------------------------------------------------------------------- \n\n");
		printf( "                      CPU REQUIREMENTS OF NUMERICAL SOLUTION\n\n");
		printf( " -------------------------------------------------------------------------------- \n");
		printf( "    The total problem used   : %10.1f seconds of CPU time\n", 0.0);
		printf( "    includes system time     : %10.1f seconds of CPU time\n", 0.0);
		printf( "    for an average of        : %10.1f  MICROSECONDS/CELL/CYCLE.\n", 0.0);
		printf( "    Total wall clock time of : %10.1f  seconds \n", 0.0);
	
	}
	~TerminalPrinter(){
	}
	TerminalPrinter(){
	}

};
#endif
