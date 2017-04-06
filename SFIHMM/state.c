#include "hmm.h"

int *fb_state(HMM *hmm) {
	double logl, norm, trial, best;
	int *assign, i, j, jfound;
	
	logl=forward_backward(hmm);
	pi(hmm);
	for(i=0;i<hmm->n;i++) {
		hmm->alpha[0][i] *= hmm->pi[i];
	}
	
	assign=(int *)malloc(hmm->t*sizeof(int));

	for(i=0;i<hmm->t;i++) {
		best=-1e300;
		jfound=-1;
		for(j=0;j<hmm->n;j++) {
			trial=hmm->alpha[i][j]*hmm->beta[i][j]; // DON'T HAVE TO WORRY ABOUT NORMALIZATIONS IN r_alpha!
			if ((trial > best) || (jfound == -1)) {
				jfound=j;
				best=trial;
			}
		}
		assign[i]=jfound;
	}
	
	return assign;
}

double *current(HMM *hmm) {
	double logl, *out, norm;
	int *assign, i, j, jfound;
	
	logl=forward_backward(hmm);
	
	out=(double *)malloc(hmm->n*sizeof(double));
	norm=0;
	for(i=0;i<hmm->n;i++) {
		out[i]=hmm->alpha[hmm->t-1][i];
		norm += hmm->alpha[hmm->t-1][i];
	}
	for(i=0;i<hmm->n;i++) {
		out[i]=out[i]/norm;
	}
	
	return out;
}

int *viterbi_state(HMM *hmm) {
	int *assign, i, j, k, ifound, jfound;
	double logl, *past, best, trial, **t1, **t2;
	
	pi(hmm);
	
	assign=(int *)malloc(hmm->t*sizeof(int));
	
	t1=(double **)malloc(hmm->n*sizeof(double *));
	t2=(double **)malloc(hmm->n*sizeof(double *));
	for(i=0;i<hmm->n;i++) {
		t1[i]=(double *)malloc(hmm->t*sizeof(double));
		t2[i]=(double *)malloc(hmm->t*sizeof(double));
		
		if (hmm->b[i][hmm->data[0]] > 0) {
			t1[i][0]=log(hmm->pi[i]) + log(hmm->b[i][hmm->data[0]]);			
		} else {
			t1[i][0]=-1000*hmm->t;
		}
		t2[i][0]=0.0;
	}
	
	// MAKE LOGS HERE FOR EASE
	for(i=0;i<hmm->n;i++) {
		for(j=0;j<hmm->n;j++) {
			if (hmm->a[i][j] == 0.0) {
				hmm->a[i][j]=-1000*hmm->t;
			} else {
				hmm->a[i][j]=log(hmm->a[i][j]);
			}
		}
		for(j=0;j<hmm->m;j++) {
			if (hmm->b[i][j] == 0.0) {
				hmm->b[i][j]=-1000*hmm->t;
			} else {
				hmm->b[i][j]=log(hmm->b[i][j]);
			}
		}
	}
	
	// DO COMPUTATION OF BEST PATHS
	for(i=1;i<hmm->t;i++) {
		for(j=0;j<hmm->n;j++) {
			best=-1e200;
			ifound=-1;
			for(k=0;k<hmm->n;k++) {
				trial=t1[k][i-1] + hmm->a[k][j] + hmm->b[j][hmm->data[i]];
				if (trial > best) {
					best=trial;
					ifound=k;
				}
			}
			t1[j][i] = best;
			t2[j][i] = ifound;
		}
		// printf("\n");
	}
	
	// set final value
	best=-1e200;
	ifound=-1;
	for(k=0;k<hmm->n;k++) {
		trial=t1[k][hmm->t-1];
		if (trial > best) {
			best=trial;
			ifound=k;
		}
	}
	assign[hmm->t-1]=ifound;
	
	// NOW WORK BACKWARDS -- assignment of earlier is given by t2 of the next.
	for(i=hmm->t-1;i>0;i--) {
		assign[i-1]=t2[assign[i]][i];
	}
	
	// CLEANUP
	for(i=0;i<hmm->n;i++) {
		free(t1[i]);
		free(t2[i]);
	}
	free(t1);
	free(t2);
	// UNDO LOGS HERE FOR EASE
	for(i=0;i<hmm->n;i++) {
		for(j=0;j<hmm->n;j++) {
			if (hmm->a[i][j] == -1000*hmm->t) {
				hmm->a[i][j]=0;
			} else {
				hmm->a[i][j]=exp(hmm->a[i][j]);
			}
		}
		for(j=0;j<hmm->m;j++) {
			if (hmm->b[i][j] == -1000*hmm->t) {
				hmm->b[i][j]=0.0;
			} else {
				hmm->b[i][j]=exp(hmm->b[i][j]);
			}
		}
	}
	
	return assign;
}

void pi(HMM *hmm) {
	gsl_eigen_nonsymmv_workspace *ws;
	gsl_matrix *amat;
	gsl_vector_complex *eval;
	gsl_complex eval_i;
	gsl_vector_complex_view evec_i;
	gsl_matrix_complex *evec;
	gsl_complex z;
	int i, j;
	double norm=0, *ans, count_neg, count_pos;
	
	if (hmm->pi != NULL) {
		free(hmm->pi);
	}
	
	ws=gsl_eigen_nonsymmv_alloc(hmm->n);
	
	amat=gsl_matrix_alloc(hmm->n, hmm->n);
	eval=gsl_vector_complex_alloc(hmm->n);
	evec=gsl_matrix_complex_alloc(hmm->n, hmm->n);

	ans=(double *)malloc(hmm->n*sizeof(double));
	
	// set up the gsl matrix
	for(i=0;i<hmm->n;i++) {
		for(j=0;j<hmm->n;j++) {
			gsl_matrix_set(amat, j, i, hmm->a[i][j]);
		}
	}
	hmm->lambda2=-1.0;
	
	// do the computation!
	gsl_eigen_nonsymmv(amat, eval, evec, ws);
	
	for(i=0;i<hmm->n;i++) {
        eval_i=gsl_vector_complex_get(eval, i);
        evec_i=gsl_matrix_complex_column(evec, i);
		
		if (fabs(GSL_REAL(eval_i)-1) < 1e-12) {
	        for(j=0;j<hmm->n;j++) {
	            z=gsl_vector_complex_get(&evec_i.vector, j);
				ans[j]=sqrt(GSL_REAL(z)*GSL_REAL(z) + GSL_IMAG(z)*GSL_IMAG(z));
				norm += ans[j];
	        }
			for(j=0;j<hmm->n;j++) {
				ans[j]=ans[j]/norm;
			}
		} else {
			if (fabs(GSL_REAL(eval_i)) > hmm->lambda2) {
				hmm->lambda2=fabs(GSL_REAL(eval_i));
			}
		}
    }
    
	gsl_matrix_free(amat);
	gsl_vector_complex_free(eval);
	gsl_matrix_complex_free(evec);
	gsl_eigen_nonsymmv_free(ws);
	
	if (norm == 0) {
		printf("ERROR: UNIT EIGENVECTOR NOT FOUND\n");
	}
	hmm->pi=ans;
}

void pi2(HMM *hmm) {
	gsl_eigen_nonsymmv_workspace *ws;
	gsl_matrix *amat;
	gsl_vector_complex *eval;
	gsl_complex eval_i;
	gsl_vector_complex_view evec_i;
	gsl_matrix_complex *evec;
	gsl_complex z;
	int i, j, neg_i, pos_i, count_neg, *neg, *pos;
	double norm=0, *ans, a_neg, a_pos, p_neg, p_pos, *pi2;
	double **eigenval_list;
	
	if (hmm->pi == NULL) {
		pi(hmm);
	}
	if (hmm->pi2 != NULL) {
		free(hmm->pi2);
	}
	
	ws=gsl_eigen_nonsymmv_alloc(hmm->n);
	
	amat=gsl_matrix_alloc(hmm->n, hmm->n);
	eval=gsl_vector_complex_alloc(hmm->n);
	evec=gsl_matrix_complex_alloc(hmm->n, hmm->n);

	ans=(double *)malloc(hmm->n*sizeof(double));
	
	// set up the gsl matrix
	for(i=0;i<hmm->n;i++) {
		for(j=0;j<hmm->n;j++) {
			gsl_matrix_set(amat, j, i, hmm->a[i][j]);
		}
	}
	
	// do the computation!
	gsl_eigen_nonsymmv(amat, eval, evec, ws);
	
	// make a big array of eigenvalues
	eigenval_list=(double **)malloc(sizeof(double *)*hmm->n);
	for(i=0;i<hmm->n;i++) {
		eigenval_list[i]=(double *)malloc(sizeof(double)*2);
		eigenval_list[i][0]=(double)i;
		eigenval_list[i][1]=GSL_REAL(gsl_vector_complex_get(eval, i));
		// printf("%i %lf\n", (int)eigenval_list[i][0], eigenval_list[i][1]);
	}
	qsort(eigenval_list, hmm->n, sizeof(double *), d_ind_compare);
	// printf("Second highest eigenvalue is %i\n", (int)eigenval_list[hmm->n-2][0]);
	
	evec_i=gsl_matrix_complex_column(evec, (int)eigenval_list[hmm->n-2][0]);
	pi2=(double *)malloc(sizeof(double)*hmm->n);
	for(i=0;i<hmm->n;i++) {
		pi2[i]=GSL_REAL(gsl_vector_complex_get(&evec_i.vector, i));
		// printf("E2(%i) %lf\n", i, pi2[i]);
	}

	count_neg=0;
	for(i=0;i<hmm->n;i++) {
		if (pi2[i] < 0) {
			count_neg++;
		}
	}
	if ((count_neg == 0) || (count_neg == hmm->n)) {
		printf("Error: second eigenvector is all positive or all negative\n");
	}
	
	neg=(int *)malloc(sizeof(int)*count_neg);
	pos=(int *)malloc(sizeof(int)*(hmm->n-count_neg));
	neg_i=0;
	pos_i=0;
	a_neg=0;
	a_pos=0;
	p_neg=0;
	p_pos=0;
	for(i=0;i<hmm->n;i++) {
		if (pi2[i] < 0) {
			// printf("%i is negative\n", i);
			neg[neg_i]=i;
			neg_i++;
			a_neg += hmm->b[i][0]*hmm->pi[i];
			p_neg += hmm->pi[i];
		} else {
			// printf("%i is positive\n", i);
			pos[pos_i]=i;
			pos_i++;
			a_pos += hmm->b[i][0]*hmm->pi[i];
			p_pos += hmm->pi[i];
		}
	}
	// printf("Negative space: %lf\n", a_neg);
	// printf("Positive space: %lf\n", a_pos);
	
	if (a_neg/p_neg > a_pos/p_pos) {
		hmm->hi_n=count_neg;
		hmm->hi_member=neg;
		hmm->lo_member=pos;
		hmm->hi_count=a_neg/p_neg;
		hmm->lo_count=a_pos/p_pos;
	} else {
		hmm->hi_n=(hmm->n-count_neg);
		hmm->hi_member=pos;
		hmm->lo_member=neg;
		hmm->hi_count=a_pos/p_pos;
		hmm->lo_count=a_neg/p_neg;
	}
	for(i=0;i<hmm->n;i++) {
		free(eigenval_list[i]);
	}
	free(eigenval_list);
	
	gsl_matrix_free(amat);
	gsl_vector_complex_free(eval);
	gsl_matrix_complex_free(evec);
	gsl_eigen_nonsymmv_free(ws);
	
	hmm->pi2=pi2;
}

void pi_n(HMM *hmm, int n, double cut) {
	gsl_eigen_nonsymmv_workspace *ws;
	gsl_matrix *amat;
	gsl_vector_complex *eval;
	gsl_complex eval_i;
	gsl_vector_complex_view evec_i;
	gsl_matrix_complex *evec;
	gsl_complex z;
	int i, state, j, counter, num_modules, num_pos, neg_i, pos_i, count_neg, *neg, *pos;
	double norm=0, *ans, a_neg, a_pos, p_neg, p_pos, *pi2;
	double **eigenval_list;
	int *list, *uniq;
	
	if (hmm->pi_n == NULL) {
		hmm->pi_n=(double **)malloc(sizeof(double *)*n);
		hmm->lambda_n=(double *)malloc(sizeof(double)*n);
		hmm->lambda_n_im=(double *)malloc(sizeof(double)*n);
	}
	if (hmm->n < n) {
		printf("You've asked me to find more eigenvectors than the matrix has rows. Something bad is going to happen.\n");
	}
	
	ws=gsl_eigen_nonsymmv_alloc(hmm->n);
	
	amat=gsl_matrix_alloc(hmm->n, hmm->n);
	eval=gsl_vector_complex_alloc(hmm->n);
	evec=gsl_matrix_complex_alloc(hmm->n, hmm->n);

	ans=(double *)malloc(hmm->n*sizeof(double));
	
	// set up the gsl matrix
	for(i=0;i<hmm->n;i++) {
		for(j=0;j<hmm->n;j++) {
			gsl_matrix_set(amat, j, i, hmm->a[i][j]);
		}
	}
	
	// do the computation!
	gsl_eigen_nonsymmv(amat, eval, evec, ws);
	
	// make a big array of eigenvalues
	eigenval_list=(double **)malloc(sizeof(double *)*hmm->n);
	for(i=0;i<hmm->n;i++) {
		eigenval_list[i]=(double *)malloc(sizeof(double)*2);
		eigenval_list[i][0]=(double)i;
		eigenval_list[i][1]=GSL_REAL(gsl_vector_complex_get(eval, i));
		// printf("%i %lf\n", (int)eigenval_list[i][0], eigenval_list[i][1]);
	}
	qsort(eigenval_list, hmm->n, sizeof(double *), d_ind_compare);
	
	for(i=0;i<n;i++) {
		evec_i=gsl_matrix_complex_column(evec, (int)eigenval_list[hmm->n-i-1][0]);
		hmm->pi_n[i]=(double *)malloc(sizeof(double)*hmm->n);
		for(j=0;j<hmm->n;j++) {
			hmm->pi_n[i][j]=GSL_REAL(gsl_vector_complex_get(&evec_i.vector, j));
		}
		hmm->lambda_n[i]=(double)eigenval_list[hmm->n-i-1][1];
		hmm->lambda_n_im[i]=gsl_complex_arg((gsl_complex)gsl_vector_complex_get(eval, (int)eigenval_list[hmm->n-i-1][0]));
		hmm->lambda_n_im[i]=GSL_IMAG(gsl_vector_complex_get(&evec_i.vector, (int)eigenval_list[hmm->n-i-1][0]));
		// printf("%i eigenvector %lf %lf %lf\n", i, hmm->lambda_n[i], hmm->lambda_n_im[i], GSL_IMAG(gsl_vector_complex_get(eval, (int)eigenval_list[hmm->n-i-1][0])));
		// for(j=0;j<hmm->n;j++) {
		// 	printf("%lf ", hmm->pi_n[i][j]);
		// }
		// printf("\n");
	}

	num_pos=0;
	for(i=1;i<n;i++) {
		if ((hmm->lambda_n[i] > cut) && (fabs(hmm->lambda_n_im[i]) < 1e-4)) {
			num_pos++;
		}
	}
	// printf("%i eigenvectors found.\n", num_pos);
	hmm->membership=(int *)malloc(sizeof(int)*n);
	list=(int *)malloc(sizeof(int)*hmm->n);
	uniq=(int *)malloc(sizeof(int)*hmm->n);
	
	for(state=0;state<hmm->n;state++) {
		hmm->membership[state]=0;
		for(i=1;i<n;i++) {
			if ((hmm->lambda_n[i] > cut) && (fabs(hmm->lambda_n_im[i]) < M_PI/5)) { // want the module to have a relaxation time of at least 5 steps!
				if (hmm->pi_n[i][state] < 0) {
					hmm->membership[state]=2*hmm->membership[state];
				} else {
					hmm->membership[state]=2*hmm->membership[state]+1;
				}
			}
		}
	}
	
	for(i=0;i<hmm->n;i++) {
		// printf("%i is module %i\n", i, hmm->membership[i]);
		list[i]=hmm->membership[i];
	}
	qsort(list, hmm->n, sizeof(int), compare);
	// for(i=0;i<hmm->n;i++) {
	// 	printf("%i sorted\n", list[i]);
	// }
	
	num_modules=1;
	i=1;
	uniq[0]=list[0];
	while(i<hmm->n) {
		if (list[i] != uniq[num_modules-1]) {
			num_modules++;
			uniq[num_modules-1]=list[i];
		}
		i++;
	}
	// printf("%i modules found\n", num_modules);
	// for(i=0;i<num_modules;i++) {
	// 	printf("%i -> %i\n", uniq[i], i);
	// }

	for(i=0;i<hmm->n;i++) {
		for(j=0;j<num_modules;j++) {
			if (hmm->membership[i] == uniq[j]) {
				break;
			}
		}
		hmm->membership[i]=j;
	}
	// for(i=0;i<hmm->n;i++) {
	// 	printf("%i is module %i\n", i, hmm->membership[i]);
	// }
	hmm->n_modules=num_modules;
	
}

