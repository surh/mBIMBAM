#ifndef __DIPLOID_H__
#define __DIPLOID_H__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "fp.h"
#include "indiv.h"
#include "model.h"

using namespace std;

class DipInd : public Individual
{
private:
    double ** 	m_phi;  		// forward probability, for each invidual, this is an MxKxK matrix; 
	double ** 	m_beta; 		// backward probability.  
	double ** 	m_top; 			// MxK, the top block of C1 formula in appendix C;  
	double ** 	m_bot; 
	double **     m_prZ; 
	int ** 	    zpair;
	double ** 	phiz;
	double *     snp_dstr;       // snp distribution;   nx2; 
public:
	DipInd();
	~DipInd();

	double get_snpmgt(int m);  
	void get_snp_dstr(int m, double * r) {if(snp_dstr !=NULL) {r[0] = snp_dstr[2*m]; r[1] = snp_dstr[2*m+1];} else r[0] = -9;}
	void allocate_snp_mgt(int n) {mgt = (double *) Allocate1D(sizeof(double),  n); }
	void allocate_snp_dstr(int n) {snp_dstr = (double *) Allocate1D(sizeof(double),  n * 2); }
	void set_snp_dstr(int m, double p0, double p1) {snp_dstr[2*m] = p0; snp_dstr[2*m+1]=p1; }
	void set_snp_mgt(int m, double s) {mgt[m] = s;}
	double * get_snp_dstr(void) {return snp_dstr;}
	double * get_snp_mgt(void) {return mgt;}
	void norm_snp_dstr(int, int);
	void calc_snp_dstr(int nLoci, int nK, double ** theta); 
	
	//for cluster-based association; working with joint_impute to get cluster membership; 
	virtual double ** Gettopptr(void) {return m_top;}
	virtual double ** Getbotptr(void) {return m_bot;}
//	virtual double ** GetprZmptr(void) {return m_prZ;}       //m_phi = m_prZ; 
	virtual int ** Getzpair(void) {return zpair;}

	virtual void joint_imputation(class ModelParam *, int, short *, int, int, int, int *);
	virtual void CalcAll(int, int, class ModelParam *, int);	
	virtual void FreeMemAfterEM(void);
	
#if defined (IMPUTATION)
	virtual void MaskSNPs(int, int, int *);      
	virtual void ImputeMaskedCounting(int, int*, int*); 
	virtual double * Getmaf(void) {return maf;}
#endif 
};

#endif
