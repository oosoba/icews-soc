#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h> // we use this to get process ID and help randomness

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>

#define BIG 1e6
#define BIGI 1e-6
#define EPSILON 1e-16

typedef struct {
	int n;	// number of states
	int m; 	// number of output symbols
	char *sym_list; // ordered list of symbols
	double	**a;	// transition probability a[i][j] -- prob that given you are in state i, you go to state j next
	double	**b;	// emission probability b[i][k] -- given you are in state i, you emit symbol k
	double	*pi;	// initial state distribution
	
	int *data; // array of data
	long t; // length of data
	
	double **alpha, **beta, **pstate;
	int *r_alpha, *r_beta, *r_pstate; // renormalization runs
	double log_l;
	double lambda2;
	int renorm;
	
	int **group;
	int *ngroup;
	int ngroups;
	
	int *proj;
	
	double hi_count, lo_count; // average fraction of first-symbol in high and low states
	int *hi_member; // states that are in the first-symbol high state
	int *lo_member; // states that are in the first-symbol high state
	int hi_n; // number of states
	double *pi2;
	double **pi_n;
	double *lambda_n;
	double *lambda_n_im;
	int n_modules;
	int *membership;
} HMM;

gsl_rng *rng();

HMM *new_hmm();
void delete_hmm(HMM *hmm);
HMM *copy_hmm(HMM *hmm);

int d_ind_compare(const void *elem1, const void *elem2);
int compare(const void *elem1, const void *elem2);

void read_obs(char *filename, HMM *hmm);
void read_hmm(char *filename, HMM *hmm);

void validate_hmm(HMM *hmm);
void init_hmm(HMM *hmm, unsigned long r_seed);
void print_hmm(HMM *hmm);
void write_hmm(char *filename, HMM *hmm);
void print_alpha(HMM *hmm, int num);

double forward_backward(HMM *hmm);
double forward_backward_fix(HMM *hmm);

void baum_welch(HMM *hmm, double **bnew);
void baum_welch_group(HMM *hmm);

HMM *fit(HMM *hmm, double tol);

void generate_string(HMM *hmm, long length, char *filename);

void init_tower(HMM *hmm, int pos, int len);

void pi(HMM *hmm);
int cmpfunc(const void *a, const void *b);
void pi2(HMM *hmm);
void pi_n(HMM *hmm, int n, double cut);

int *fb_state(HMM *hmm);
int *viterbi_state(HMM *hmm);
double *current(HMM *hmm);

void create_groups(HMM *hmm, int n);
void init_groups(HMM *hmm);
void project_states(HMM *hmm, int print_ans);
void simple_project(HMM *hmm);

void generate_string_pstart(HMM *hmm, long length, char *filename, double *pstart, int pos);

int hamming(int *one, int *two, int n);

void project_bw(HMM *a);
HMM *b_step(HMM *input);
HMM *a_step(HMM *input);
HMM *add(HMM *a, HMM *b, double factor);
HMM *multiply(HMM *a, double factor);

HMM *fa(HMM *a, double beta);
HMM *fb(HMM *a, double beta);