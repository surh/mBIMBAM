#ifndef __HAPLOID_H__
#define __HAPLOID_H__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "fp.h"
#include "indiv.h"
#include "model.h"

using namespace std;

class HapInd : public Individual
{
private:
	double ** m_phi; 						// forward MxK
	double ** m_beta; 					// backward MxK
	double ** m_top; 
	double ** m_bot; 
	double ** m_prZ; 
	double *  snp_dstr; 	             //nx2; 
public:
	HapInd();
	~HapInd();
	
	double get_snpmgt(int m) {if(mgt) return mgt[m]; else if(snp_dstr) return snp_dstr[2*m+1]; else return (double) (snpGT[m]-'0'); }; 
	void get_snp_dstr(int m, double * r) { r[0] = snp_dstr[2*m];r[1] = snp_dstr[2*m+1];}
	void allocate_snp_mgt(int n) {mgt = (double *) Allocate1D(sizeof(double),  n); }
	void allocate_snp_dstr(int n) {snp_dstr = (double *) Allocate1D(sizeof(double),  n * 2); }
	void set_snp_dstr(int m, double p0, double p1) {snp_dstr[2*m] = p0; snp_dstr[2*m+1]=p1; }
	void set_snp_mgt(int m, double s) {mgt[m] = s; }
	double * get_snp_dstr(void) {return snp_dstr;}
	double * get_snp_mgt(void) {return mgt;}
	void norm_snp_dstr(int, int);
	void calc_snp_dstr(int nLoci, int nK, double ** theta); 

	virtual inline int** Getzpair(void) {return NULL;}
	virtual double ** Gettopptr(void) {return m_top;}
	virtual double ** Getbotptr(void) {return m_bot;}
//	virtual double *** GetprZmptr(void) {return &m_prZ;}   
	
	virtual void CalcAll(int, int, class ModelParam *, int);
	virtual void FreeMemAfterEM(void); 
	virtual void joint_imputation(class ModelParam *, int, short *, int, int, int, int *);
	
#if defined (IMPUTATION)
	virtual void MaskSNPs(int, int, int *);      
	virtual void ImputeMaskedCounting(int, int*, int*); 
	virtual double * Getmaf(void) {return maf;}
#endif 
};

#endif 
