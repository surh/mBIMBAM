#ifndef __MODEL_H__                
#define __MODEL_H__

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>

#include "indiv.h" 
#include "diploid.h"
#include "haploid.h"
#include "fpmath.h"
#include "param.h"

#if defined (MPI_ENABLED)
#include "mpi.h"
#endif

using namespace std;

class SNPSummaryData
{
public:                               
	string rs; 
	char major; 
	char minor; 
	int ni;      	 //total number of individual;  
	double sg,sd, sgd, sg2, sd2, sy, syg, syd, sy2 ; 
}; 

class CaseCtrlData
{
public:
    int ni, np;         // ni = num of indiv, np = dim of beta;
    gsl_matrix * mx;    //design matrix, nxp; 
    gsl_vector * vy;    //phenotype; 
    gsl_vector * delta; //prior on beta ~ N(0, delta^2);
public:
    CaseCtrlData(int, int, double *); 
    ~CaseCtrlData(); 
};
		 
double func_f(const gsl_vector * beta, void * params);
void func_df(const gsl_vector * beta, void * params, gsl_vector * df);
void func_fdf(const gsl_vector * beta, void * params, double * f, gsl_vector * df);
double Hessian(const gsl_vector * beta, void * params);
double prDM1(CaseCtrlData * pD, gsl_vector * x);
//utility functions for case/ctrl calculations. 

double mylike_f(const gsl_vector * beta, void * params); 


class PosnBF 
{
public:
	vector<int> pos;
	double bf; 
	double var; 
	PosnBF(); 
	~PosnBF(); 
};

typedef vector< pair< char, char> > vSNP;
typedef map<string, vSNP> HASH;

class mmSNP 
{
public:
	char A; 
	char B; 
	double countA, countB, countQ; 
public:
	mmSNP() {A=B='N'; countA=countB=countQ=0;}
	~mmSNP(){;} 
	void assign(class mmSNP s) {A=s.A; B=s.B; countA=s.countA; countB=s.countB; countQ=s.countQ;}
};   //for merge SNPs; 

class cSNP
{
private:
	string rs; 
	int flag; 
	int chr; 
	long pos;     //the physical location; 
	long index;  //the position in genotype matrix; [0 .. nLoci-1];  
	char major; 
	char minor; 
	double maf; 
public:
	cSNP(){;} 
	~cSNP(){;} 
	cSNP(string, int, long, long, pair<char, char>, double); 

	string get_rs(void) {return rs;}
	int get_chr(void) {return chr;}
	long get_pos(void) {return pos;}
	long get_index(void) {return index;}
	int get_flag(void) {return flag;} 
	            
	void set_flag(void) {flag++;}
	void set_rs(string s) {rs = s;}
	void set_chr(int s) {chr = s;}
	void set_pos(long s) {pos = s;}
	void set_index(long s) {index = s;}
};

class cGENE
{
private:	
	string gn;
	int chr;          
	long beg, end;     //if beg == end == -1; then it's a non gene; 
	
public:
	cGENE(); 
	~cGENE();
	cGENE(string gn, int chr, long beg, long end); 
	void assign(class cGENE); 
	string get_name(void) {return gn;}
	int get_chr(void) {return chr;}
	long get_beg(void) {return beg;}
	long get_end(void) {return end;}

};

class ModelnData
{
public: 
	fstream fplog; // log;
	
private:
	int m_not_snp;              //if 1 then use snpmgt and ignore the range checking. note: must no missing data.  
	short ** snpIn01; 
	double ** snpmgt; 
	
	vector<class cSNP> vsnp; 
	vector<class cGENE> vgene; 

private:
    short 	nK; 				//number of clusters;  
    class 	ModelParam * pMP;   //nEMruns dimension. 
/***************************************************************************  
  above for model parameters                                                                    
***************************************************************************/ 
private:
#if defined (MPI_ENABLED)
    int     procID;
    int     nProc;  
    vector<int> * pvLoadDip;
	vector< vector<int> > vvLoad;    //vvLoad.at(procID).at(wp) = corresponding rows in snpInpr; 
	int * 	loadLoci;        	//for scan;
	int 	nSlice;          	//max_{all}{loadLoci[i+1]-loadLoci[i]}; 
	int *   loadPerm;        	//for multi-snp;
	int 	nMaxLoad;      		//max_{all processes} pvLoadDip.size(); 
    double *  pParamCollect; 
#endif      

private:
	int randSeed;               //user may specify random seed;
	int nEMRuns;                //how many different em runs;
	int nMaxSteps;              //maximum steps of each em run;
	int nWarmSteps; 
	int nImpute; 
	double total_likelihood; 
	double logLikelihood;      	//loglikelihood of pop. 
	string fnEM; 
	
	double ** sumiTop;
	double ** sumiBot;
	double *** sumiJmk;
    double ** sumiJm; 
/***************************************************************************  
  above for EM runs;
***************************************************************************/ 
private: 
	int nSubPop;  				//num. of subpopulations; 
	vector<int> vnSubHap;       //num. of haplotypes in subpop; 
	int nDip;    				//number of diploid individuals;
	int nHap;    				//number of haploid individuals; 
	int nIndiv; 
	int nCohort; 				//nIndiv = nCohort + nPanel; 
	int nPanel; 
    int nLoci;                  //number of snp loci, same for everyone;
	string hasPopStructure; 		
	vector<string> indivID;     // individual IDs; 
	class Individual ** pIndiv; 
	
	fstream logfile; 
	vector<string> vsRsnum; 
	vector<string> vGin;
	vector<string> vPin; 
	vector<string> vPos; 
	vector<string> vSSD; 
	string fnOutput;  			//output prefix; 
	string fnRsPos; 
	vector<int> vFileIndiv;       //num of diploids in genotype files; 
	vector<int> vFilePloid;     //ploid of input genotype files. 
	  
    string fnFILE;    //a dummy filename for flexible use; 
	string fnDOC; 
	
	vector<vector<char> > vvisQ; 
	
	map<string, pair<char, char> > mapRs2mm;  //rs to major, minor alleles.  
	map<string, double> mapRs2maf;              //rs to minor allele freq. 
	map<string, double> mapRs2var;            //rs to var = 2*maf*(1-maf). 
	map<string, long> mapRs2pos; 
	map<string, int> mapRs2chr; 
	map<string, class SNPSummaryData> mapRs2SNPSummary; 
	int nPH; 
	vector<vector<double> > vv_phval; 
	vector<vector<int> > vv_phdex;                                             
	vector<double> vsigma_a; 
	vector<double> vsigma_d; 
	vector<double> m_beta; 
	map<string, vector<double> > mapRs2beta; 
	
	int cc; //indicator for case_contrl; 

	int which_pheno; 
	vector<double>  m_pv;            //p-vale is [nLoci+1] array, the last number is region p-value; 
	vector<double>  m_sumbf; 

	string fnGene; 
	int nGeneFlank; 

	double m_current_maf;          //for bfmaf calculations. note this is dangerous!!!. 
private:
	int m_ignore_phase;            //ignore phase even when there is a '='; 
	int m_silence; 
	int m_allele_coding_mode; 
	double m_exclude_maf; 
	double m_exclude_miss; 
	int m_exclude_nopos; 
	vector<int> v_reduced_num_indiv; 
public: 
	inline void Setsilence(int s) {m_silence = s;}
	inline void Set_reduced_num_indiv(int s) {v_reduced_num_indiv.push_back(s);}
	inline void Setexcludemaf(double s) {m_exclude_maf = s;}
	inline void Setexcludemiss(double s) {m_exclude_miss = s;}
	inline void Setexcludenopos(int s) {m_exclude_nopos = s;}
	inline void Setallelecoding(int s){m_allele_coding_mode = s;} 
 	void print_progress_bar(int last, char* str, int p, int total);     
private:
#if defined (IMPUTATION)
	double percentMasked;       // nLoci * this = nMasked;
	vector<string> vMaskedSNPrs; 
	vector<string> vPARCrs; 
	vector<int> vPARCpos; 
	int  parcLoci; 
	int  nMasked; 
#endif	
	
	int  nMultiSnp; 
	int  nLevel; 
	vector<class PosnBF> vPerm;
	map<int, int> mrank2loc; 
	int m_num; 
	int m_df;     //degree of freedom; 
	int m_sortQ;  //if sort the snp according to bf value or not. 
	map<string, int> mPanelRsQ;   
public:
	ModelnData();
	~ModelnData();
	
#	if defined (MPI_ENABLED)
	inline void SetprocID(int s) {procID = s;}
	inline void SetnProc(int s) {nProc = s;}
	inline int  GetprocID(void) {return procID;}
	inline int  GetnProc(void) {return nProc;}
	void 		Loading(void); 
#	endif
    
	inline void set_ignore_phase(void) {m_ignore_phase = 1;}
	inline void Setmnotsnp(int s) {m_not_snp = s;}
	inline void SetsortQ(int s) {m_sortQ = s;}
	inline void Setdf(int s) {m_df = s;}
	inline void Setnum(int s) {m_num = s;}
	inline void Setsigmaa(double s) {vsigma_a.push_back(s);}
	inline void Setsigmad(double s) {vsigma_d.push_back(s);}

	inline void SetfnFILE(string s) {fnFILE.assign(s);}
	inline void SetfnDOC(string s) {fnDOC.assign(s);}
	inline int  GetnPanel(void) {return nPanel;}
	inline int  GetnCohort(void) {return nCohort;}
	inline void Setflank(int s) {nGeneFlank = s;}
	inline int  Getflank(void) {return nGeneFlank;}
	inline void SetfnGene(string s) {fnGene.assign(s);}
	inline string GetfnGene(void) {return fnGene;}            // for gene study. 
	inline void SetOnCaseCtrl(void) {cc = 1;}
	inline void SethasPopStructure(string s) {hasPopStructure.assign(s);}
	inline void SetvGin(string s) {vGin.push_back(s);}
	inline void SetvPin(string s) {vPin.push_back(s);}
	inline void SetvPos(string s) {vPos.push_back(s);}
	inline void SetvSSD(string s) {vSSD.push_back(s);}
	inline void SetfnOutput(string s) {fnOutput.assign(s);}
	inline string GetfnOutput(void) {return fnOutput;}
	inline void SetfnRsPos(string s) {fnRsPos.assign(s);}
	inline void SetnImpute(int s) {nImpute = s;}
	inline void SetnLevel(int s) {nLevel = s;}
	inline void SetnMultiSnp(int s) {nMultiSnp = s;}

	inline void SetnK(short s) {nK = s;}//only in CommandLine(...)
	inline void SetrandSeed(long int s){randSeed = s;}
	inline void SetnEMRuns(int s){nEMRuns = s;}
	inline void SetnWarmSteps(int s){nWarmSteps = s;}
	inline void SetnMaxSteps(int s){nMaxSteps = s;}
	inline void SetnDip(int s){nDip = s;}
	inline void SetnHap(int s){nHap = s;}
	inline void SetnLoci(int s){nLoci = s;}
	inline void SetnPH(int s){nPH = s;}	

	inline int GetnMaxSteps(void) {return nMaxSteps;}
	inline int GetnK(void) {return nK;}
	inline int GetnImpute(void) {return nImpute;}
	inline int GetnMultiSnp(void) {return nMultiSnp;}
	inline int GetnLevel(void) {return nLevel;}
	inline int Getnum_Gin(void) {return (int) vGin.size();}
	inline int Getnum_Pin(void) {return (int) vPin.size();}
	inline long int GetrandSeed(void){return randSeed;}
	inline int GetnLoci(void){return nLoci;}   
										//used by fileio.cpp to allocate memory for snpGT. 
	inline int GetnDip(void){return nDip;}
	inline int GetnHap(void){return nHap;}
	inline int GetnSubPop(void){return nSubPop;}
	inline int GetnEMRuns(void){return nEMRuns;}
	inline int GetnWarmSteps(void) {return nWarmSteps;}

#if defined (IMPUTATION)
	inline void Setmask(double s){percentMasked = s;}
	inline double Getmask(void) {return percentMasked;}; 
	void MaskSNP(void);					//randomly mask genotype. substitute 0,1's with ?'s.
	void MaskRand(void); 
	void ImputeMasked(void);             //missing genotype '?' imputation. 
#endif
	//interface
    void InitModelParam(void); 			//initialize parameters.
	void assign_population_label(void);			//structure the population even if only one subpop.

	void EM(int);     //expectation-maximization. if arg = 1, then only use panel. otherwise use all of them. 
	void SNP_Density(void); 
	
    void single_snp_importance_sampling(double *, double *);
	void single_snp_exact(int, double *); 
	void single_snp_mean_genotype(int, double *);
	
    double f_stat_mean(vector<double>* phval, vector<int>*, vector<int>* curindex, double ** prob_snp); 
	double f_stat(double *, double *, int); 
	double snp_lrt(vector<double>* phval, vector<int>* base, vector<int>* curindex, double ** prob_snp); 
	double snp_lrt(double * phval, double * genotype, int ni); 
	
	double cc_bf_core(double, double, int ni, int ns, class CaseCtrlData* pNul, class CaseCtrlData* pAlt); 
	double cc_bf(double, double, double *, double *, int);        //for eact 
	double cc_bf_gd(double, double, vector<double>*, vector<int>*, vector<int>*, double **);    //for genotype distribution
	double cc_bf_mgt(double, double, vector<double>*, vector<int>*, vector<int>*, double *);    //for mean genotype
	double cc_bf(double, double, vector<double>*, vector<int>*, short **, class PosnBF *);     //for multi-snp;
	double cc_bf(double, double, vector<double>*, vector<int>*, double **, class PosnBF *);     //for multi-snp;
	double cc_bf_mgt(double, double, vector<double>*, vector<int>*, double **, class PosnBF *);     //for multi-snp;
    double cc_bf(double, double, vector<double>&, vector<int>&, double *);     //for single snp; 
	//for case/ctrl bf calc. 
	double bf_core(double inv_va, double inv_vd, int ni, int ns, gsl_matrix * gXX, gsl_vector * gph);
	double calc_bf(double, double, double *, double *, int);     //for exact 
	double calc_bf_gd(double, double, vector<double>*, vector<int>*, vector<int>*, double **);  //for genotype distribution
	double calc_bf_mgt(double, double, vector<double>*, vector<int>*, vector<int>*, double *);  //for mean genotype
	double calc_bf(double, double, vector<double>*, vector<int>*, short **, class PosnBF *);   //for multi-snp;
	double calc_bf(double, double, vector<double>*, vector<int>*, double **, class PosnBF *);   //for multi-snp;
	double calc_bf_mgt(double, double, vector<double>*, vector<int>*, double **, class PosnBF *);   //for multi-snp;
	double calc_bf(double, double, vector<double>&, vector<int>&, double *);    //for single snp. 
	double calc_residual_logbf(double, double *, double *, int);     
	// variety of calc_bf's to adapt to different geno/phene input, main job is done in bf_core; 
	//for quatitative phenotype bf calc. 
	void single_snp_cluster(void); 
	void single_snp_pvalue(int); 
	void single_snp_analysis(int);     
	void multi_snp_analysis(int); 
	                                                    

	void read_gene(void); 
	int  read_rs2pos(string, long, long, long);  
	int  read_bimbam_genotype(int, long, long);
	int  read_bimbam_genotype_distribution(int, int, long, long);
	int  read_bimbam_phenotype(int);

	char merge_mmSNP(class mmSNP & ref, class mmSNP & s); 
	int compatibleQ(pair<char, char> p1, pair<char, char> p2); 
#if defined (POWER)
	void read_is_write();
	void read_gd_stat_write();
#endif          
		
	void write_genotype(int, int);
	void write_summary(double, double **); 
	void open_log(void); 
	void close_log(void); 

	void write_em_param(void); 
	int  read_em_param(string); 
	void clean_geno_pheno(void); 

	void read_snp_summary_data(string); 
	void get_snp_summary_data(int ns, class SNPSummaryData *);
	double calc_bf_of_snp_summary(double, double, class SNPSummaryData *);
	void combine_snp_summary_data(class SNPSummaryData *, class SNPSummaryData *); 
	void CombineStudy(void);
	void ProduceSummaryData(int, string); 
	void process_prior(void); 

	void write_for_gain(); 
	void fix_mean_genotype(); 
	void simulate_phenotype(int ncausal, double hh); 
	void for_barbara(); 

//multi-phenotype; 
private:
	int nn;   //n; 
	int np;   //p; 
	int nd;   //d; 
	//same notation as in paper. 

	double ** gt;    //nxp; 
	double ** ph;    //nxd; 
	//data; 

	int pm;
	double ** psi;   //dxd; 
	double sa;  //sigma_a;   K=diag(0, sa, sa, ...); 
	//prior; 

public:

   void logbf_summaries(gsl_matrix *, gsl_matrix * m1syy, gsl_vector * v1s11, int * zz, int pm, double sa, double * lbf);
   void logbf_rank1(double sa, double pi0, int pm, int * config, double * lbf); 
   void single_snp_multi_phenotype(int mode); 
   double compute_prior(int *, double); 
	
   double calc_bf_hpass(double sigma_a, double sigma_d, int ni, double* vph, vector<int>* cin, double ** mgt); 
	
   void hpass(void); 
};    

#endif

