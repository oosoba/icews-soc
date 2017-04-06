#include "hmm.h"

#define NORM 1e-2
#define NORM_INTERVENE 1e-12

gsl_rng *rng() {
	unsigned long s;
	FILE *fn;
	gsl_rng *r_global;
	
	fn = fopen("/dev/urandom", "rb"); 		
	if (fread(&s, sizeof(s), 1, fn) != 1) 
		exit(-1); /* Failed! */
	r_global=gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(r_global, s);
	fclose(fn);
	
	return r_global;
}

int hamming(int *one, int *two, int n) {
	int i, dist=0;
	
	for(i=0;i<n;i++) {
		if (one[i] != two[i]) {
			dist++;
		}
	}
	
	return dist;
}


HMM *new_hmm() {
	HMM *hmm;
	
	hmm=(HMM *)malloc(sizeof(HMM));
	hmm->a=NULL;
	hmm->b=NULL;
	hmm->alpha=NULL;
	hmm->data=NULL;
	hmm->sym_list=NULL;
	hmm->pi=NULL;
	hmm->proj=NULL;
	hmm->group=NULL;
	hmm->pi2=NULL;
	hmm->hi_member=NULL;
	hmm->lo_member=NULL;
	
	return hmm;
}

void delete_hmm(HMM *hmm) {
	int i, j;
	if (hmm != NULL) {
		for(i=0;i<hmm->n;i++) {
			free(hmm->a[i]);
			free(hmm->b[i]);
		}
			
		free(hmm->a);
		free(hmm->b);

		if (hmm->alpha != NULL) {
			for(i=0;i<hmm->t;i++) {
				free(hmm->alpha[i]);
				free(hmm->beta[i]);
			}
			free(hmm->alpha);
			free(hmm->beta);
			free(hmm->r_alpha);
			free(hmm->r_beta);
		}

		if (hmm->sym_list != NULL) {
			free(hmm->sym_list);
		}

		if (hmm->data != NULL) {
			free(hmm->data);
		}

		if (hmm->pi != NULL) {
			free(hmm->pi);
		}

		if (hmm->pi2 != NULL) {
			free(hmm->pi2);
		}

		if (hmm->hi_member != NULL) {
			free(hmm->hi_member);
		}

		if (hmm->lo_member != NULL) {
			free(hmm->lo_member);
		}

		if (hmm->proj != NULL) {
			free(hmm->proj);
		}

		if (hmm->group != NULL) {
			for(i=0;i<hmm->ngroups;i++) {
				free(hmm->group[i]);
			}
			free(hmm->group);
			free(hmm->ngroup);
		}
		
		free(hmm);
	}
}
// ubuntu 14.04 
HMM *copy_hmm(HMM *hmm) {
	int i, j;
	HMM *dup;
	
	dup=new_hmm();
	dup->n=hmm->n;
	dup->m=hmm->m;
	dup->t=hmm->t;
	
	dup->log_l=hmm->log_l;
	dup->lambda2=hmm->lambda2;

	dup->a=(double **)malloc(hmm->n*sizeof(double *));
	dup->b=(double **)malloc(hmm->n*sizeof(double *));
	
	for(i=0;i<hmm->n;i++) {
		dup->a[i]=(double *)malloc(hmm->n*sizeof(double));
		dup->b[i]=(double *)malloc(hmm->m*sizeof(double));
				
		for(j=0;j<hmm->n;j++) {
			dup->a[i][j]=hmm->a[i][j];
		}
		for(j=0;j<hmm->m;j++) {
			dup->b[i][j]=hmm->b[i][j];
		}
	}
	
	dup->sym_list=(char *)malloc(hmm->m*sizeof(char));
	for(i=0;i<hmm->m;i++) {
		dup->sym_list[i]=hmm->sym_list[i];
	}
	
	if (hmm->data != NULL) {
		dup->data=(int *)malloc(hmm->t*sizeof(int));
		for(i=0;i<hmm->t;i++) {
			dup->data[i]=hmm->data[i];
		}
	} else {
		dup->data=NULL;
	}
	
	if (hmm->pi != NULL) {
		dup->pi=(double *)malloc(hmm->n*sizeof(double));
		for(i=0;i<hmm->n;i++) {
			dup->pi[i]=hmm->pi[i];
		}
	} else {
		dup->pi=NULL;
	}
	
	if (hmm->alpha != NULL) {
		dup->r_alpha=(int *)malloc(hmm->t*sizeof(int));
		dup->r_beta=(int *)malloc(hmm->t*sizeof(int));
		dup->renorm=hmm->renorm;
		
		dup->alpha=(double **)malloc(hmm->t*sizeof(double *));
		dup->beta=(double **)malloc(hmm->t*sizeof(double *));
		
		for(i=0;i<hmm->t;i++) {
			dup->r_alpha[i]=hmm->r_alpha[i];
			dup->r_beta[i]=hmm->r_beta[i];

			dup->alpha[i]=(double *)malloc(hmm->n*sizeof(double));
			dup->beta[i]=(double *)malloc(hmm->n*sizeof(double));
			
			for(j=0;j<hmm->n;j++) {
				dup->alpha[i][j]=hmm->alpha[i][j];
				dup->beta[i][j]=hmm->beta[i][j];
			}
		}		
	} else {
		dup->alpha=NULL;
	}
	
	if (hmm->proj != NULL) {
		dup->proj=(int *)malloc(hmm->t*sizeof(int));
		for(i=0;i<hmm->t;i++) {
			dup->proj[i]=hmm->proj[i];
		}
	} else {
		dup->proj=NULL;
	}
	
	if (hmm->group != NULL) {
		dup->group=(int **)malloc(hmm->ngroups*sizeof(int *));
		dup->ngroup=(int *)malloc(hmm->ngroups*sizeof(int));
		dup->ngroups=hmm->ngroups;
		for(i=0;i<dup->ngroups;i++) {
			dup->group[i]=(int *)malloc(hmm->ngroup[i]*sizeof(int));
			dup->ngroup[i]=hmm->ngroup[i];
			for(j=0;j<hmm->ngroup[i];j++) {
				dup->group[i][j]=hmm->group[i][j];
			}
		}
	} else {
		dup->group=NULL;
	}
	
	return dup;
}

void init_hmm(HMM *hmm, unsigned long r_seed) {
	int i, j;
	const gsl_rng_type *T;
	double *dir_n, *dir_m; 
	gsl_rng *r;
	FILE *fn; 
	
	// printf("Starting init\n");
	
	// fn = fopen("/dev/urandom", "rb"); 		
	// if (fread(&r_seed, sizeof(r_seed), 1, fn) != 1) {
	// 	/* Failed!--use time instead; beware, could sync with other instances */	
	// 	printf("Warning: urandom read fail; using system clock\n");
	// 	r_seed=(unsigned long)time(NULL);
	// }
	// fclose(fn);
	// r_seed=(unsigned long)time(NULL);
	// printf("Init--got the seed %lu\n", r_seed);
	
	if (hmm->a == NULL) {
		hmm->a=(double **)malloc((hmm->n)*sizeof(double *));
		hmm->b=(double **)malloc((hmm->n)*sizeof(double *));
		for(i=0;i<hmm->n;i++) {
			hmm->a[i]=(double *)malloc((hmm->n)*sizeof(double));
			hmm->b[i]=(double *)malloc((hmm->m)*sizeof(double));
		}
	}
	// printf("Init--made the a and b\n");
	
	gsl_rng_env_setup();
	T=gsl_rng_default;
	r=gsl_rng_alloc(T);
	gsl_rng_set(r, r_seed);
	// printf("Init--made gsl rng set\n");

	dir_n=(double *)malloc((hmm->n)*sizeof(double));
	dir_m=(double *)malloc((hmm->m)*sizeof(double));
	for(i=0;i<hmm->n;i++) {
		dir_n[i]=1.0;
	}
	for(i=0;i<hmm->m;i++) {
		dir_m[i]=1.0;
	}

	for(i=0;i<hmm->n;i++) {
		gsl_ran_dirichlet(r, hmm->n, dir_n, hmm->a[i]);
		gsl_ran_dirichlet(r, hmm->m, dir_m, hmm->b[i]);
	}
	// printf("Init--dirichlet done\n");
	
	free(dir_n);
	free(dir_m);
	
	if (hmm->sym_list == NULL) {
		hmm->sym_list = (char *)malloc((hmm->m)*sizeof(char));
		for(i=0;i<hmm->m;i++) {
			hmm->sym_list[i]=(i+65);
		}
	}
	
	hmm->proj=NULL;
	hmm->group=NULL;
	
	// printf("Freeing rng r\n");
	gsl_rng_free(r);
}

void print_hmm(HMM *hmm) {
	int i,j;
	
	printf("log-l=%18.16lf\n", hmm->log_l);
	printf("N=%i\n", hmm->n);
	printf("M=%i\n", hmm->m);
	
	for(i=0;i<hmm->n;i++) {
		for(j=0;j<hmm->n;j++) {
			printf("%12.10lf ", hmm->a[i][j]);
		}
		printf("\n");
	}

	for(i=0;i<hmm->n;i++) {
		for(j=0;j<hmm->m;j++) {
			printf("%12.10lf ", hmm->b[i][j]);
		}
		printf("\n");
	}
	// printf("Output symbols: \n");
	// for(i=0;i<hmm->n;i++) {
	// 	printf("State %i: ", i);
	// 	for(j=0;j<hmm->m;j++) {
	// 		printf("%c (%18.16lf) ", hmm->sym_list[j], hmm->b[i][j]);
	// 	}
	// 	printf("\n");
	// }
	
	if (hmm->group != NULL) {
		printf("%i groups\n", hmm->ngroups);
		for(i=0;i<hmm->ngroups;i++) {
			printf("GROUP %i: ", i);
			for(j=0;j<hmm->ngroup[i];j++) {
				printf(" %i", hmm->group[i][j]);
			}
			printf("\n");
		}
	}
}

void write_hmm(char *filename, HMM *hmm) {
	int i, j;
	FILE *f_out;
	
	f_out=fopen(filename, "w");
	fprintf(f_out, "log-l=%18.16lf\n", hmm->log_l);
	fprintf(f_out, "N=%i\n", hmm->n);
	fprintf(f_out, "M=%i\n", hmm->m);
		
	for(i=0;i<hmm->n;i++) {
		for(j=0;j<hmm->n;j++) {
			fprintf(f_out, "%18.16lf ", hmm->a[i][j]);
		}
		fprintf(f_out, "\n");
	}
	
	for(i=0;i<hmm->n;i++) {
		for(j=0;j<hmm->m;j++) {
			fprintf(f_out, "%18.16lf ", hmm->b[i][j]);
		}
		fprintf(f_out, "\n");
	}
	for(i=0;i<hmm->m;i++) {
		fprintf(f_out, "%c", hmm->sym_list[i]);
		// printf("%i %c\n", hmm->m, hmm->sym_list[i]);
	}
	fprintf(f_out, "\n");
	if (hmm->pi != NULL) {
		for(i=0;i<hmm->n;i++) {
			fprintf(f_out, "%18.16lf ", hmm->pi[i]);
		}
		fprintf(f_out, "\n");
		fprintf(f_out, "lambda2=%18.16lf\n", hmm->lambda2);
	}
	fclose(f_out);	
}

void read_hmm(char *filename, HMM *hmm) {
	int i, j;
	char c;
	FILE *f_in;
	
	f_in=fopen(filename, "r");
	
	fscanf(f_in, "log-l=%lf\n", &(hmm->log_l));
	fscanf(f_in, "N=%i\n", &(hmm->n));
	fscanf(f_in, "M=%i\n", &(hmm->m));
	
	hmm->a=(double **)malloc(hmm->n*sizeof(double *));
	hmm->b=(double **)malloc(hmm->n*sizeof(double *)); // for each state, there are hmm->m symbol outputs
	
	for(i=0;i<hmm->n;i++) {
		hmm->a[i]=(double *)malloc(hmm->n*sizeof(double));
		for(j=0;j<hmm->n;j++) {
			fscanf(f_in, "%lf", &(hmm->a[i][j]));
		}
	}
	
	for(i=0;i<hmm->n;i++) {
		hmm->b[i]=(double *)malloc(hmm->m*sizeof(double));
		for(j=0;j<hmm->m;j++) {
			fscanf(f_in, "%lf", &(hmm->b[i][j]));
		}
	}
	fscanf(f_in, "%c", &c); // there's a space, and a carriage return to read
	fscanf(f_in, "%c", &c);
	
	hmm->sym_list = (char *)malloc((hmm->m)*sizeof(char));
	for(i=0;i<hmm->m;i++) {
		fscanf(f_in, "%c", &(hmm->sym_list[i]));
	}
	pi(hmm);
	
	fclose(f_in);	
}

void validate_hmm(HMM *hmm) {
	int i, j;
	double sum;
	
	for(i=0;i<hmm->n;i++) {
		sum=0.0;
		for(j=0;j<hmm->n;j++) {
			sum += hmm->a[i][j];
		}
		// if (fabs(sum - 1.0) > NORM) {
		// 	printf("TRANSITION MATRIX NOT NORMALIZED for %i\n", i);
		// 	exit(-1);
		// }
		if (fabs(sum - 1.0) > NORM_INTERVENE) {
			for(j=0;j<hmm->n;j++) {
				hmm->a[i][j] /= sum;
			}
		}		
	}

	for(i=0;i<hmm->n;i++) {
		sum=0.0;
		for(j=0;j<hmm->m;j++) {
			sum += hmm->b[i][j];
		}
		// if (fabs(sum - 1.0) > NORM) {
		// 	printf("EMISSION MATRIX NOT NORMALIZED for %i\n", i);
		// 	exit(-1);
		// }
		if (fabs(sum - 1.0) > NORM_INTERVENE) {
			for(j=0;j<hmm->m;j++) {
				hmm->b[i][j] /= sum;
			}
		}		
	}
	
	
}

void print_alpha(HMM *hmm, int num) {
	int i, j;
	
	for(i=0;i<num;i++) {
		printf("%i: ", i);
		for(j=0;j<hmm->n;j++) {
			printf("[%i] %lg (%lg) ",j,  hmm->alpha[i][j], hmm->beta[i][j]);
		}
		printf("\n");
	}
}





