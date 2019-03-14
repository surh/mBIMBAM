#include <iostream>
#include <string.h>
#include <time.h>
#include "fpmath.h"
#include "param.h"

using namespace std;
ModelParam::ModelParam()
{
	;
}

ModelParam::~ModelParam(void)
{
	if (theta) {
		Free2DMatrix(theta); 
		theta = NULL;
	}
	if (alpha) {
		Free3DMatrix(alpha); 
		alpha = NULL;
	}
	if (r) {
		Free2DMatrix(r); 
		r = NULL; 
	}
}

void ModelParam::Init(int nLoci, int nK, int nSubPop)
{
	theta = Allocate2DMatrix(nLoci, nK);
	for (int m = 0; m < nLoci; m++)
		for (int k = 0; k < nK; k++)
		{
			theta[m][k] = gsl_ran_beta(gsl_r, 1.0, 1.0); 
			if(theta[m][k] > 0.999) theta[m][k] = 0.999; 
			else if(theta[m][k] < 0.001) theta[m][k] = 0.001; 
		}

	alpha = Allocate3DMatrix(nSubPop, nLoci, nK);
	r = Allocate2DMatrix(nSubPop, nLoci);
	double * dp = new double[nK]; 
	double * res = new double[nK]; 
	for (int k = 0; k < nK; k++)
		dp[k] = 1.0; // / nK; 
	for (int s = 0; s < nSubPop; s++)
	{
		for (int m = 0; m < nLoci; m++)
		{
			r[s][m] = 0.001; 
			gsl_ran_dirichlet(gsl_r, nK, dp, alpha[s][m]); 
		}
	}
	delete[] dp; 
	delete[] res; 
}
