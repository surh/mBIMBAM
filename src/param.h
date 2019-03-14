#ifndef __PARAM_H__
#define __PARAM_H__

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "diploid.h"
#include "haploid.h"
#include "fpmath.h"
#include <iostream>
using namespace std;
		 
//differenet populations share theta while having different alpha and r;
class ModelParam
{                                    
private:
	double ** theta;  			// MxK, emission para;
	double *** alpha; 
	double ** r; 
public:
	ModelParam();
	~ModelParam(); 
	void Init(int, int, int); 
	inline void Settheta(int m, int k, double s) {theta[m][k] = s;}
	inline void Setalpha(int p, int m, int k, double s) {alpha[p][m][k] = s;}
	inline void Setr(int p, int m, double s) {r[p][m] = s;}
	//set is used to update parameters, one by one; 
	//while get functions return pointers;
	inline double ** Gettheta(void) {return theta;}            
	inline double ** Gettheta(int p) {return theta;}            
	
	inline double ** Getalpha(void) {return alpha[0];} 
	inline double * Getr(void) {return r[0];}
	//if don't specify which subpoplation, return the first one; 
	inline double ** Getalpha(int p) {return alpha[p];} 
	inline double * Getr(int p) {return r[p];}
};

#endif
