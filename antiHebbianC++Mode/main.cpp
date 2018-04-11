#include <iostream>
#include "antiHebbianC++Mode.h"
#include "RNGProcess.h"
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ppl.h>
using namespace std;
using namespace concurrency;

ofstream outfile1;
ofstream outfile2;
ofstream outfile3;
ofstream outfile4;

int main(void)
{
	double b = 0;
	int propor = 0;

	antiHebbian anti(b, propor);
	//b = 1;
	//propor = 40;
	anti.initial();
	anti.initialAdjacency();
	anti.tongji();
	double xx, xx1;
	xx = anti.cooperator1; xx1 = anti.defector1;


	cout << "cooperator:" << xx / (xx + xx1) << endl;
	cout << "defector:" << xx1 / (xx + xx1) << endl;
}