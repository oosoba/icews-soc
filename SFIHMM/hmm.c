#include "hmm.h"


int compare(const void *elem1, const void *elem2) {
	char f = *((char *)elem1);
	char s = *((char *)elem2);

	return (f > s) - (f < s);
}

int d_ind_compare(const void *elem1, const void *elem2) {
	double *f = *((double **)elem1);
	double *s = *((double **)elem2);

	return (f[1] > s[1]) - (f[1] < s[1]);
}

double simple_pow(int n) {
	double accum;
	int i;
	
	if (n == 0) {
		return 1.0;
	} else {
		accum=1.0;
		if (n > 0) {
			for(i=0;i<n;i++) {
				accum = accum*BIGI;
			}
		} else {
			for(i=0;i<-n;i++) {
				accum = accum*BIG;
			}			
		}
		return accum;
	}
}

void read_obs(char *filename, HMM *hmm) {
	FILE *f_in;
	char *obs;
	// char obs_sort[4]= { 'C', 'C', 'R', 'C'};
	char *obs_sort, running_counter;
	int pos;
	long i;
	
	f_in=fopen(filename, "r");
	fscanf(f_in, "%li\n", &(hmm->t));
	
	obs=(char *)malloc((hmm->t)*sizeof(char));
	obs_sort=(char *)malloc((hmm->t)*sizeof(char));
	
	for(i=0; i<(hmm->t); i++) {
		fscanf(f_in, "%c", &obs[i]);
		obs_sort[i]=obs[i];
	}
	fclose(f_in);
	
	qsort(obs_sort, hmm->t, sizeof(char), compare);
	
	
	if (hmm->sym_list == NULL) {
		
		running_counter=obs_sort[0];
		(hmm->m)=1;
		for(i=1; i<(hmm->t); i++) {
			if (obs_sort[i] != running_counter) {
				running_counter=obs_sort[i];
				hmm->m++;
			}
		}
		
		(hmm->sym_list)=(char *)malloc(sizeof(char)*hmm->m);
		pos=0;
		(hmm->sym_list)[pos]=obs_sort[0];
		pos=1;
		for(i=1; i<(hmm->t); i++) {
			if (obs_sort[i] != obs_sort[i-1]) {
				(hmm->sym_list)[pos]=obs_sort[i];
				pos++;
			}
		}
	}
	
	(hmm->data)=(int *)malloc((hmm->t)*sizeof(int));
	for(i=0; i<(hmm->t); i++) {
		pos=0;
		while(obs[i] != (hmm->sym_list)[pos]) {
			pos++;
		}
		(hmm->data)[i]=pos;
	}
	
	free(obs);
	free(obs_sort);
}

void init_tower(HMM *hmm, int pos, int len) { 
	// initializes the HMM to create a "tower" in position pos of length len
	// a tower is a primitive counter
	int i, j;
	const gsl_rng_type *T;
	gsl_rng *r;
	FILE *fn; 
	unsigned long r_seed;
	
	fn = fopen("/dev/urandom", "rb"); 		
	if (fread(&r_seed, sizeof(r_seed), 1, fn) != 1) 
		exit(-1); /* Failed! */

	gsl_rng_env_setup();

	T=gsl_rng_default;
	r=gsl_rng_alloc(T);
	gsl_rng_set(r, r_seed);

	hmm->n -= len;
	
	for(i=0;i<(hmm->n+len);i++) {
		for(j=hmm->n;j<(hmm->n+len);j++) {
			hmm->a[i][j]=0.0;
			hmm->a[j][i]=0.0;
		}
	}
	
	hmm->a[pos][hmm->n]=gsl_rng_uniform(r); // probability of entering the tower
	
	for(i=hmm->n;i<hmm->n+len;i++) {
		// up, stay, or out...
		if ((i+1) < hmm->n+len) {
			hmm->a[i][i+1]=gsl_rng_uniform(r);
		}
		if (i > hmm->n) {
			hmm->a[i][i-1]=gsl_rng_uniform(r);
		} else {
			hmm->a[i][pos]=gsl_rng_uniform(r);			
		}
		
		if (i == hmm->n+len-1) {
			hmm->a[i][i]=gsl_rng_uniform(r);
		}
		// hmm->a[i][pos]=gsl_rng_uniform(r);
		// hmm->a[i][i]=gsl_rng_uniform(r);
		
		// clean output symbols
		hmm->b[i]=(double *)malloc(hmm->m*sizeof(double));
		hmm->b[i][0]=1.0; //gsl_rng_uniform(r);
		for(j=1;j<hmm->m;j++) {
			hmm->b[i][j]=gsl_rng_uniform(r);
		}
	}
	
	hmm->n += len;
	
	gsl_rng_free(r);
	validate_hmm(hmm);
}

void init_groups(HMM *hmm) {
	int i, j, a, b;
	double denom;
	
	for(i=0;i<hmm->ngroups;i++) {
		for(j=0;j<hmm->ngroups;j++) {
			if (i != j) {
				
				for(a=0;a<hmm->ngroup[i];a++) {
					for(b=0;b<hmm->ngroup[j];b++) {
						
						hmm->a[hmm->group[i][a]][hmm->group[j][b]]=0;
						
					}
				}
				
			}
		}
	}
	
	for(i=0;i<hmm->n;i++) {
		denom=0;
		for(j=0;j<hmm->n;j++) {
			denom += hmm->a[i][j];
		}
		for(j=0;j<hmm->n;j++) {
			hmm->a[i][j] = hmm->a[i][j]/denom;
		}		
	}
}

void create_groups(HMM *hmm, int n) {
	int i, j;
	
	hmm->ngroups=n;
	hmm->ngroup=(int *)malloc(n*sizeof(int));
	hmm->group=(int **)malloc(n*sizeof(int *));
	
	for(i=0;i<hmm->ngroups;i++) {
		hmm->ngroup[i]=hmm->n/n;
		hmm->group[i]=(int *)malloc(hmm->ngroup[i]*sizeof(int));
		for(j=0;j<hmm->n/n;j++) {
			hmm->group[i][j]=(hmm->n/n)*i+j;
		}
	}
}

void project_states(HMM *hmm, int print_ans) {  // remember, alpha[t][state]
	int i, j, k, a, b, ai, bi, start, finish, nflips=0, original, fixed;
	double **current_state, tot, current, *best_prob, **flips=NULL, trans, *best_prob_curr, *best_prob_prev, top, bottom;
	double tot_alpha, tot_beta, denom_alpha, denom_beta;
		
	best_prob_prev=(double *)malloc(hmm->n*sizeof(double));
	best_prob_curr=(double *)malloc(hmm->n*sizeof(double));

	current_state=(double **)malloc(hmm->ngroups*sizeof(double *));
	for(i=0;i<hmm->ngroups;i++) {
		current_state[i]=(double *)malloc(2*sizeof(double *));
		current_state[i][0]=(double)i;
	}
	
	if (hmm->proj == NULL) {
		hmm->proj=(int *)malloc(hmm->t*sizeof(int));
	}
	
	for(i=0;i<hmm->t;i++) { // for each time-step
		
		for(j=0;j<hmm->ngroups;j++) { // initialize the probability of being in each state
			current_state[j][0]=(double)j;
			current_state[j][1]=0;
		}

		for(j=0;j<hmm->n;j++) {
			if (i > 0) {
				best_prob_prev[j]=best_prob_curr[j];	
			}
			best_prob_curr[j]=0;
		}
				
		tot=0;
		for(k=0;k<hmm->ngroups;k++) { // for each of the ngroups...
			current=0;
			for(a=0;a<hmm->ngroup[k];a++) {
				current += hmm->alpha[i][hmm->group[k][a]]*hmm->beta[i][hmm->group[k][a]];	// add in the probability of being in one of the states			
					
			}
			current_state[k][1] = current;
			tot += current;
		}
		
		for(k=0;k<hmm->ngroups;k++) {
			current_state[k][1] /= tot;
		}
				
		qsort(current_state, hmm->ngroups, sizeof(double *), d_ind_compare);
		hmm->proj[i]=(int)current_state[hmm->ngroups-1][0];

		if (print_ans) {
			printf("%i ", hmm->proj[i]);
		}
		
		if ((i > 0) && (hmm->proj[i] != hmm->proj[i-1])) {
			nflips++;
			flips=(double **)realloc(flips, nflips*sizeof(double *));
			flips[nflips-1]=(double *)malloc(sizeof(double)*2);
			flips[nflips-1][0]=nflips-1;
			flips[nflips-1][1]=current_state[hmm->ngroups-1][1];
		}
		
	}
	if (print_ans) {
		printf("\n");
	}
	
	qsort(flips, nflips, sizeof(double *), d_ind_compare);
	for(i=0;i<nflips-100000;i++) {
		start=flips[i][0];
		original=hmm->proj[start];
		fixed=hmm->proj[start-1];
		while(hmm->proj[start] == original) {
			hmm->proj[start]=fixed;
			start++;
		}
	}
	
	// NOW WE HAVE TO FIX THE ALPHA and BETA matrices
	for(i=0;i<hmm->t;i++) {
				
		for(j=0;j<hmm->ngroups;j++) {
			
			if (j != hmm->proj[i]) {
				
				for(k=0;k<hmm->ngroup[j];k++) {
					
					hmm->alpha[i][hmm->group[j][k]] *= 0;
					hmm->beta[i][hmm->group[j][k]] *= 0;	
					// do we renorm?
				}
				
			}	
		}
		
	}
	
}


double forward_backward(HMM *hmm) { // recomputes the alpha and beta matrices, returns current log-likelihoo
	int i, j, k, a;
	int t;
	double sum, full_sum, l;
	
	// make the alpha, beta, pstat matrices (and their renormalization factors)
	if (hmm->alpha == NULL) { 
		hmm->alpha = (double **)malloc(hmm->t*sizeof(double *));
		hmm->beta = (double **)malloc(hmm->t*sizeof(double *));
		
		for(i=0;i<hmm->t;i++) {
			hmm->alpha[i] = (double *)malloc(hmm->n*sizeof(double *));
			hmm->beta[i] = (double *)malloc(hmm->n*sizeof(double *));
		}
		
		hmm->r_alpha=(int *)malloc(hmm->t*sizeof(int));
		hmm->r_beta=(int *)malloc(hmm->t*sizeof(int));
		
		// initialize the initial states
		for(i=0;i<hmm->n;i++) { 
			hmm->alpha[0][i] = hmm->b[i][hmm->data[0]];
			hmm->beta[hmm->t-1][i] = 1.0;
		}
	}
		
	hmm->r_alpha[0]=0;
	hmm->r_beta[hmm->t-1]=0;
	
	// run through, doing all the alphas
	for(t=1;t<hmm->t;t++) { 
		full_sum=0.0;
		
		for(j=0;j<hmm->n;j++) { // for each j			
			sum = 0.0;
			for(i=0;i<hmm->n;i++) {
				sum += hmm->alpha[t-1][i]*hmm->a[i][j]*hmm->b[j][hmm->data[t]]; // *(i->j prob)*(emit data[t] in state j prob)
				// printf("%i to %i, emit %i at %i prob: %18.16lf\n", i, j, hmm->data[t], t, hmm->alpha[t-1][i]*hmm->a[i][j]*hmm->b[j][hmm->data[t]]);
				// printf("TERMS: %18.16lf %18.16lf %18.16lf\n", hmm->alpha[t-1][i], hmm->a[i][j], hmm->b[j][hmm->data[t]]);
			}				
			hmm->alpha[t][j]=sum;
			full_sum += sum;
		}
		
		hmm->r_alpha[t]=hmm->r_alpha[t-1];
		if (full_sum < BIGI) {
			hmm->r_alpha[t]++;
			for(j=0;j<hmm->n;j++) {
				hmm->alpha[t][j] = hmm->alpha[t][j]*BIG;
			}
		}
		
	}
	
	// run through, doing all the betas
	for(t=hmm->t-2;t>=0;t--) {
		full_sum=0.0;
		
		for(i=0;i<hmm->n;i++) { // for each i			
			sum = 0.0;
			for(j=0;j<hmm->n;j++) {
				sum += hmm->a[i][j]*hmm->b[j][hmm->data[t+1]]*hmm->beta[t+1][j]; // *(i->j prob)*(emit data[t] in state j prob)
			}				
			hmm->beta[t][i]=sum;
			full_sum += sum;			
		}
		
		hmm->r_beta[t]=hmm->r_beta[t+1];
		if (full_sum < BIGI) {
			hmm->r_beta[t]++;
			for(j=0;j<hmm->n;j++) {
				hmm->beta[t][j] *= BIG;
			}
		}

	}

	// for(i=0;i<hmm->t;i++) {
	// 	for(j=0;j<hmm->n;j++) {
	// 		if (hmm->alpha[i][j] < EPSILON) {
	// 			hmm->alpha[i][j]=EPSILON;
	// 		}
	// 		if (hmm->beta[i][j] < EPSILON) {
	// 			hmm->beta[i][j]=EPSILON;
	// 		}
	// 	}
	// }
	
	l = 0;
	for(i=0;i<hmm->n;i++) {
		l += hmm->alpha[0][i]*hmm->beta[0][i];
	}

	hmm->log_l = log(l) - log(BIG)*(hmm->r_alpha[0]+hmm->r_beta[0]);
	hmm->renorm = hmm->r_alpha[0]+hmm->r_beta[0];

	// if (fabs(hmm->log_l) > 1e200) {
	// 	printf("ERROR %lg %lg\n", hmm->alpha[0][0], hmm->beta[0][0]);
	// 	exit(-1);
	// }

	return hmm->log_l;
}

double forward_backward_fix(HMM *hmm) { // recomputes the alpha and beta matrices, returns current log-likelihoo
	int i, j, k, a;
	int t;
	double sum, full_sum, l, **set;
	
	// make the alpha, beta, pstat matrices (and their renormalization factors)
	if (hmm->alpha == NULL) { 
		hmm->alpha = (double **)malloc(hmm->t*sizeof(double *));
		hmm->beta = (double **)malloc(hmm->t*sizeof(double *));
		
		for(i=0;i<hmm->t;i++) {
			hmm->alpha[i] = (double *)malloc(hmm->n*sizeof(double *));
			hmm->beta[i] = (double *)malloc(hmm->n*sizeof(double *));
		}
		
		hmm->r_alpha=(int *)malloc(hmm->t*sizeof(int));
		hmm->r_beta=(int *)malloc(hmm->t*sizeof(int));
		
		// // initialize the initial states
		// set=(double **)malloc(hmm->n*sizeof(double *));
		// for(i=0;i<hmm->n;i++) { 
		// 	set[i]=(double *)malloc(2*sizeof(double));
		// 	set[i][0]=i;
		// 	set[i][1]=hmm->alpha[0][i];
		// }
		// qsort(set, hmm->n, sizeof(double *), d_ind_compare);
		// for(i=0;i<hmm->n;i++) {
		// 	hmm->alpha[i]=0;
		// }
		// hmm->alpha[(int)set[hmm->n-1][0]]=1;
		// 
		// // best
		// for(i=0;i<hmm->n;i++) { 
		// 	set[i][0]=i;
		// 	set[i][1]=hmm->beta[0][i];
		// }
		// qsort(set, hmm->n, sizeof(double *), d_ind_compare);
		// for(i=0;i<hmm->n;i++) {
		// 	hmm->beta[i]=0;
		// }
		// hmm->beta[(int)set[hmm->n-1][0]]=1;
		for(i=0;i<hmm->n;i++) { 
			hmm->alpha[0][i] = hmm->b[i][hmm->data[0]];
			hmm->beta[hmm->t-1][i] = 1.0;
		}
		
	}
		
	hmm->r_alpha[0]=0;
	hmm->r_beta[hmm->t-1]=0;
	
	// run through, doing all the alphas
	for(t=1;t<hmm->t;t++) { 
		full_sum=0.0;
		
		for(j=0;j<hmm->n;j++) { // for each j			
			sum = 0.0;
			for(i=0;i<hmm->n;i++) {
				sum += hmm->alpha[t-1][i]*hmm->a[i][j]*hmm->b[j][hmm->data[t]]; // *(i->j prob)*(emit data[t] in state j prob)
				// printf("%i to %i, emit %i at %i prob: %18.16lf\n", i, j, hmm->data[t], t, hmm->alpha[t-1][i]*hmm->a[i][j]*hmm->b[j][hmm->data[t]]);
				// printf("TERMS: %18.16lf %18.16lf %18.16lf\n", hmm->alpha[t-1][i], hmm->a[i][j], hmm->b[j][hmm->data[t]]);
			}				
			hmm->alpha[t][j]=sum;
			full_sum += sum;
		}
		
		hmm->r_alpha[t]=hmm->r_alpha[t-1];
		if (full_sum < BIGI) {
			hmm->r_alpha[t]++;
			for(j=0;j<hmm->n;j++) {
				hmm->alpha[t][j] = hmm->alpha[t][j]*BIG;
			}
		}
		
	}
	
	// run through, doing all the betas
	for(t=hmm->t-2;t>=0;t--) {
		full_sum=0.0;
		
		for(i=0;i<hmm->n;i++) { // for each i			
			sum = 0.0;
			for(j=0;j<hmm->n;j++) {
				sum += hmm->a[i][j]*hmm->b[j][hmm->data[t+1]]*hmm->beta[t+1][j]; // *(i->j prob)*(emit data[t] in state j prob)
			}				
			hmm->beta[t][i]=sum;
			full_sum += sum;			
		}
		
		hmm->r_beta[t]=hmm->r_beta[t+1];
		if (full_sum < BIGI) {
			hmm->r_beta[t]++;
			for(j=0;j<hmm->n;j++) {
				hmm->beta[t][j] *= BIG;
			}
		}

	}

	for(i=0;i<hmm->t;i++) {
		for(j=0;j<hmm->n;j++) {
			if (hmm->alpha[i][j] < 0) {
				hmm->alpha[i][j]=0;
			}
			if (hmm->beta[i][j] < 0) {
				hmm->beta[i][j]=0;
			}
		}
	}
	
	l = 0;
	for(i=0;i<hmm->n;i++) {
		l += hmm->alpha[0][i]*hmm->beta[0][i];
	}

	hmm->log_l = log(l) - log(BIG)*(hmm->r_alpha[0]+hmm->r_beta[0]);
	hmm->renorm = hmm->r_alpha[0]+hmm->r_beta[0];

	// if (fabs(hmm->log_l) > 1e200) {
	// 	printf("ERROR %lg %lg\n", hmm->alpha[0][0], hmm->beta[0][0]);
	// 	exit(-1);
	// }

	return hmm->log_l;
}

void baum_welch_group(HMM *hmm) {
	double num, sum, denom, term, **avg;
	int i, j, k, a, b, ip, jp;
	int t;
		
	avg=(double **)malloc(hmm->ngroups*sizeof(double *));
	for(i=0;i<hmm->ngroups;i++) {
		avg[i]=(double *)malloc(hmm->ngroups*sizeof(double));
	}
	for(i=0;i<hmm->t;i++) {
		for(j=0;j<hmm->n;j++) {
			if (hmm->alpha[i][j] < EPSILON) {
				hmm->alpha[i][j]=EPSILON;
			}
			if (hmm->beta[i][j] < EPSILON) {
				hmm->beta[i][j]=EPSILON;
			}
		}
	}
	
	for(a=0;a<hmm->ngroups;a++) {
		for(ip=0;ip<hmm->ngroup[a];ip++) {
			
			i=hmm->group[a][ip];
						
			denom = 0.0;
			for(t=0;t<hmm->t-1;t++) {
				term = hmm->alpha[t][i]*hmm->beta[t][i]*simple_pow(hmm->r_alpha[t]+hmm->r_beta[t]-hmm->renorm); //
				denom += term;
			}
		
			for(b=0;b<hmm->ngroups;b++) {
				
				if (a != b) {

					for(jp=0;jp<hmm->ngroup[b];jp++) {
						j=hmm->group[b][jp];

						num = 0.0;
						for(t=0;t<hmm->t-1;t++) {
							num += hmm->alpha[t][i]*hmm->b[j][hmm->data[t+1]]*hmm->beta[t+1][j]*simple_pow(hmm->r_alpha[t]+hmm->r_beta[t+1]-hmm->renorm);
							if (fabs(num) > 1e100) {
								printf("ERROR: %i %i %i %i %i %lg %lg %lg\n", t, i, hmm->r_alpha[t], hmm->r_beta[t+1], hmm->renorm, hmm->alpha[t][i], hmm->beta[t+1][j], simple_pow(hmm->r_alpha[t]+hmm->r_beta[t+1]-hmm->renorm));
								exit(-1);
							}
						}
						hmm->a[i][j] *= EPSILON + num/(EPSILON+denom);
						
					}
					
				}
				
			}
			
		}
	}
	
	for(a=0;a<hmm->ngroups;a++) {
		for(b=0;b<hmm->ngroups;b++) {
			if (a != b) {
				sum=0;
				for(ip=0;ip<hmm->ngroup[a];ip++) {
					for(jp=0;jp<hmm->ngroup[b];jp++) {
						sum += hmm->a[hmm->group[a][ip]][hmm->group[b][jp]];
					}
				}
				sum=sum/(double)(hmm->ngroup[a]*hmm->ngroup[b]);
				for(ip=0;ip<hmm->ngroup[a];ip++) {
					for(jp=0;jp<hmm->ngroup[b];jp++) {
						hmm->a[hmm->group[a][ip]][hmm->group[b][jp]]=sum;
					}
				}
			}
		}
	}

	
	for(i=0;i<hmm->n;i++) {
		denom=0;
		for(j=0;j<hmm->n;j++) {
			denom += hmm->a[i][j];
		}
		for(j=0;j<hmm->n;j++) {
			hmm->a[i][j] = hmm->a[i][j]/(1e-128+denom);
		}		
	}
				
	for(i=0;i<hmm->ngroups;i++) {
		free(avg[i]);
	}
	free(avg);
}

void baum_welch(HMM *hmm, double **bnew) {
	double num, denom, term;
	int i, j, k, a;
	long t;
	
	
	for(i=0;i<hmm->n;i++) {
		for(k=0;k<hmm->m;k++) {
			bnew[i][k]=0.0;
		}
	}
	
	for(i=0;i<hmm->n;i++) {
		
		denom = 0.0;
		for(t=0;t<hmm->t-1;t++) {
			term = hmm->alpha[t][i]*hmm->beta[t][i]*simple_pow(hmm->r_alpha[t]+hmm->r_beta[t]-hmm->renorm);
			denom += term;
			bnew[i][hmm->data[t]] += term;
		}
		for(k=0;k<hmm->m;k++) {
			bnew[i][k] = bnew[i][k]/denom;
		}
		
		for(j=0;j<hmm->n;j++) {
			num = 0.0;
			for(t=0;t<hmm->t-1;t++) {
				num += hmm->alpha[t][i]*hmm->b[j][hmm->data[t+1]]*hmm->beta[t+1][j]*simple_pow(hmm->r_alpha[t]+hmm->r_beta[t+1]-hmm->renorm);
			}
			hmm->a[i][j] *= num/denom;
		}
				
	}
	
	for(i=0;i<hmm->n;i++) {
		for(j=0;j<hmm->m;j++) {
			hmm->b[i][j]=bnew[i][j];
		}
	}
							
}

HMM *fit(HMM *hmm, double tol) {
	double log_l, new_log_l, **bnew;
	int i;
	
	log_l=forward_backward(hmm);
	validate_hmm(hmm);
	
	new_log_l=log_l+1.0;
	
	bnew=(double **)malloc(hmm->n*sizeof(double *));
	for(i=0;i<hmm->n;i++) {
		bnew[i]=(double *)malloc(hmm->m*sizeof(double));
	}
	
	while(fabs(log_l - new_log_l) > tol) {
		log_l=new_log_l;
		baum_welch(hmm, bnew);
		new_log_l=forward_backward(hmm);
	}
	
	for(i=0;i<hmm->n;i++) {
		free(bnew[i]);
	}
	free(bnew);
	
	return hmm;
}

void generate_string(HMM *hmm, long length, char *filename) { // starts with the actual first data point
	long i;
	int current;
	double *pstart;
	FILE *f_out, *fn;
	gsl_ran_discrete_t **a_samp, **b_samp, *start_samp;
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	unsigned long r_seed;
	
	if ((filename != NULL)) {
		f_out=fopen(filename, "w");
		fprintf(f_out, "%li\n", length);
	} else {
		printf("%li\n", length);
	}
		
	fn = fopen("/dev/urandom", "rb"); 		
	if (fread(&r_seed, sizeof(r_seed), 1, fn) != 1) 
		exit(-1); /* Failed! */
	
	T=gsl_rng_default;
	r=gsl_rng_alloc(T);
	gsl_rng_set(r, r_seed);
	
	// first, prepare some samplers
	a_samp=(gsl_ran_discrete_t **)malloc(hmm->n*sizeof(gsl_ran_discrete_t *));
	b_samp=(gsl_ran_discrete_t **)malloc(hmm->n*sizeof(gsl_ran_discrete_t *));
	for(i=0;i<hmm->n;i++) {
		a_samp[i]=gsl_ran_discrete_preproc(hmm->n, hmm->a[i]);
		b_samp[i]=gsl_ran_discrete_preproc(hmm->m, hmm->b[i]);
	}
	
	if (hmm->data != NULL) {
		pstart=(double *)malloc(hmm->n*sizeof(double));
		for(i=0;i<hmm->n;i++) {
			pstart[i]=hmm->b[i][hmm->data[0]];
		}
		start_samp=gsl_ran_discrete_preproc(hmm->n, pstart);

		// pick a start state
		current=gsl_ran_discrete(r, start_samp);		
	} else {
		current=0;
		for(i=1;i<(hmm->n)*100;i++) {
			current=gsl_ran_discrete(r, a_samp[current]);
		}			
	}
	
	if ((filename != NULL)) {
		fprintf(f_out, "%c", hmm->sym_list[(int)gsl_ran_discrete(r, b_samp[current])]);
	} else {
		printf("%c", hmm->sym_list[(int)gsl_ran_discrete(r, b_samp[current])]);
	}
	
	for(i=1;i<length;i++) {
		current=gsl_ran_discrete(r, a_samp[current]);
		if ((filename != NULL)) {
			fprintf(f_out, "%c", hmm->sym_list[(int)gsl_ran_discrete(r, b_samp[current])]);
		} else {
			printf("%c", hmm->sym_list[(int)gsl_ran_discrete(r, b_samp[current])]);
		}
	}	
	
	if ((filename != NULL)) {
		fprintf(f_out, "\n");
		fclose(f_out);
	} else {
		printf("\n");
	}
		
}

void generate_string_pstart(HMM *hmm, long length, char *filename, double *pstart, int pos) { // starts with the actual first data point
	long i;
	int current, sampled;
	FILE *f_out, *fn;
	gsl_ran_discrete_t **a_samp, **b_samp, *start_samp;
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	unsigned long r_seed;
		
	if ((filename != NULL)) {
		f_out=fopen(filename, "w");
		fprintf(f_out, "%li\n", length);
	}
		
	fn = fopen("/dev/urandom", "rb"); 		
	if (fread(&r_seed, sizeof(r_seed), 1, fn) != 1) 
		exit(-1); /* Failed! */
	
	T=gsl_rng_default;
	r=gsl_rng_alloc(T);
	gsl_rng_set(r, r_seed);
	
	// first, prepare some samplers
	a_samp=(gsl_ran_discrete_t **)malloc(hmm->n*sizeof(gsl_ran_discrete_t *));
	b_samp=(gsl_ran_discrete_t **)malloc(hmm->n*sizeof(gsl_ran_discrete_t *));
	for(i=0;i<hmm->n;i++) {
		a_samp[i]=gsl_ran_discrete_preproc(hmm->n, hmm->a[i]);
		b_samp[i]=gsl_ran_discrete_preproc(hmm->m, hmm->b[i]);
	}
	
	start_samp=gsl_ran_discrete_preproc(hmm->n, pstart);
	current=gsl_ran_discrete(r, start_samp);		
	
	sampled=(int)gsl_ran_discrete(r, b_samp[current]);
	if ((filename != NULL)) {
		fprintf(f_out, "%c", hmm->sym_list[sampled]);
	} else {
		// printf("%i ", current);
	}
	if (hmm->data != NULL) {
		hmm->data[pos]=sampled;
	}
	
	
	for(i=1;i<length;i++) {
		current=gsl_ran_discrete(r, a_samp[current]);
		sampled=(int)gsl_ran_discrete(r, b_samp[current]);
		if ((filename != NULL)) {
			fprintf(f_out, "%c", hmm->sym_list[sampled]);
		} else {
			// printf("%i ", current);
		}
		if (hmm->data != NULL) {
			hmm->data[pos+i]=sampled;
		}
	}	
	
	if ((filename != NULL)) {
		fprintf(f_out, "\n");
		fclose(f_out);
	} else {
		printf(" ");
	}
		
}

HMM *b_step(HMM *input) {
	HMM *dup;
	double logl;
	
	dup=copy_hmm(input);
	baum_welch_group(dup);
	// logl=forward_backward(dup);
	
	return dup;
}

HMM *a_step(HMM *input) {
	HMM *dup;
	
	dup=copy_hmm(input);
	project_states(dup, 0);
	
	return dup;
}

HMM *add(HMM *a, HMM *b, double factor) {
	HMM *dup;
	int t, i, j, delta_renorm;
	double beta_true, alpha_true, t_comb;
	
	dup=copy_hmm(a);
	
	for(t=0;t<a->t;t++) {
		for(i=0;i<a->n;i++) {
			
			// t_comb=log(a->alpha[t][i]) + factor*log(b->alpha[t][i]);
			// dup->alpha[t][i] = exp(t_comb);
			// 
			// t_comb=log(a->beta[t][i]) + factor*log(b->beta[t][i]);
			// dup->beta[t][i] = exp(t_comb);
			
			if (b->r_alpha[t] > a->r_alpha[t]) {
				dup->alpha[t][i] = a->alpha[t][i]*simple_pow(a->r_alpha[t]-b->r_alpha[t]) + factor*b->alpha[t][i];
				dup->r_alpha[t]=b->r_alpha[t]-a->renorm+b->renorm;				
			} else {
				dup->alpha[t][i] = a->alpha[t][i] + factor*b->alpha[t][i]*simple_pow(b->r_alpha[t]-a->r_alpha[t]) ;				
			}
			// 
			if (b->r_beta[t] > a->r_beta[t]) {
				dup->beta[t][i] = a->beta[t][i]*simple_pow(a->r_beta[t]-b->r_beta[t]) + factor*b->beta[t][i];
				dup->r_beta[t]=b->r_beta[t];
			} else {
				dup->beta[t][i] = a->beta[t][i] + factor*b->beta[t][i]*simple_pow(b->r_beta[t]-a->r_beta[t]) ;				
			}
			
		}
	}

	return dup;
}

HMM *multiply(HMM *a, double factor) {
	HMM *dup;
	int t, i, j;
	double t_comb;
	
	dup=copy_hmm(a);
	
	for(t=0;t<a->t;t++) {
		for(i=0;i<a->n;i++) {
			
			// t_comb=factor*log(EPSILON+a->alpha[t][i]);
			// dup->alpha[t][i] = exp(t_comb);
			// 
			// t_comb=factor*log(EPSILON+a->beta[t][i]);
			// dup->beta[t][i] = exp(t_comb);
			dup->alpha[t][i] = factor*a->alpha[t][i];
			dup->beta[t][i] = factor*a->beta[t][i];
		}
	}

	return dup;
	
	return dup;	
}

HMM *fa(HMM *a, double beta) {
	HMM *one, *two, *three, *four;
	
	one=copy_hmm(a);
	project_states(one, 0);
	
	two=add(one, a, -1.0);
	three=multiply(two, 1.0/beta);
	
	four=add(one, three, -1.0);

	delete_hmm(one);
	delete_hmm(two);
	delete_hmm(three);
	
	return four;
}

HMM *fb(HMM *a, double beta) {
	HMM *one, *two, *three, *four;
	double logl;
	int i;
	
	one=copy_hmm(a);
	project_bw(one);
	
	two=add(one, a, -1.0); // preserve the new AIJ
	three=multiply(two, 1.0/beta);
	four=add(one, three, 1.0);

	delete_hmm(one);
	delete_hmm(two);
	delete_hmm(three);

	return four;
}

void project_bw(HMM *a) {
	double old_logl=-1e64, logl;
	
	baum_welch_group(a);
	logl=forward_backward(a);	
	while(old_logl-logl > 1e-2) {
		old_logl=logl;
		baum_welch_group(a);
		logl=forward_backward(a);
	}
	
}

