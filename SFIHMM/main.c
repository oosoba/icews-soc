#include "hmm.h"

// specify number of states, then the input observations file, then (optionally) the initial conditions

// NOT PARALLEL
// hmm -g N_letters hmm.dat [filename] : generate an N_letter sequence realized from hmm.dat, using the stationary distribution, save to filename (optional)
// hmm -l obs.dat hmm.dat : give the likelihood of generating obs.dat given the machine hmm.dat
// hmm -v obs.dat hmm.dat : give the Viterbi maximum-likelihood path
// hmm -c obs.dat hmm.dat : give the current distribution over states
// hmm -s obs.dat hmm.dat : give the Forward-Backward estimate of state occupation
// hmm -p obs.dat hmm.dat tol : polish previous solution
// hmm -f N_iter N_states obs.dat : find an N_state HMM for obs.dat, write to obs.dat_OUT_[N_states]

int main (int argc, char *argv[]) {
	HMM  *hmm, *hbest;
	int i, in, iter, j, node=0, size, tower, isdet, *st, ham;
	int level, n_levels, start_states, end_states, nstates, best_states;
	long nobs;
	double t0, **message, bestlog_l, *logl_list, *logl_consolidate, tol, *curr, norm, logl_aic, trap_h, trap_l, n_h, n_l;
	char *filename_node, *command, *filename, current_state;
	FILE *f_out;
	FILE *fn; 
	unsigned long *r_seed,pid;
	

	t0=clock();

	if ((argc == 1) || (argv[1][0] != '-')) {

		printf("Greetings, Professor Falken. Please specify a command-line option.\n");
		
	} else {

		if (((int)argv[1][1] != 'f') & ((int)argv[1][1] != 'F')) {

			if (node == 0) {

				switch((int)argv[1][1]) {
					case 'g':

						sscanf(argv[2], "%li", &nobs); 
						hmm=new_hmm();
						read_hmm(argv[3], hmm);	

						pi(hmm);
						if (argc == 5) {
							generate_string(hmm, nobs, argv[4]);
						} else {
							generate_string(hmm, nobs, NULL);
						}					

						break;
					case 'l':
						hmm=new_hmm();
						read_hmm(argv[3], hmm);
						read_obs(argv[2], hmm);
						
						printf("log-l of HMM: %lf\n", forward_backward(hmm));

						break;
					case 'v':
						hmm=new_hmm();
						read_hmm(argv[3], hmm);
						read_obs(argv[2], hmm);
						logl_aic=forward_backward(hmm);
						printf("Log-l on observations: %lf\n", logl_aic);
						
						st=viterbi_state(hmm);
						for(i=0;i<hmm->t;i++) {
							printf("%i ", st[i]);
						}
						printf("\n");
						free(st);
						break;
					case 'c':
						hmm=new_hmm();
						read_obs(argv[2], hmm);
						read_hmm(argv[3], hmm);

						curr=current(hmm);
						for(i=0;i<hmm->n;i++) {
							printf("%18.16lf ", curr[i]);
						}
						printf("\n");
						free(curr);
						break;
					case 's':
						hmm=new_hmm();
						read_obs(argv[2], hmm);
						read_hmm(argv[3], hmm);
						
						st=fb_state(hmm);
						for(i=0;i<hmm->t;i++) {
							printf("%i ", st[i]);
						}
						printf("\n");
						free(st);
						break;
					case 'p':
						hmm=new_hmm();
						read_obs(argv[2], hmm);
						read_hmm(argv[3], hmm);
						sscanf(argv[4], "%lf", &tol); 
						
						printf("Reading in... %lf current log-l\n", forward_backward(hmm));
						hmm=fit(hmm, tol);
						printf("New is %lf log-l\n", forward_backward(hmm));
						forward_backward(hmm);
						filename_node=(char *)malloc(255*sizeof(char));
						sprintf(filename_node, "%s_polished", argv[3]);
						
						write_hmm(filename_node, hmm);
						break;
					case 'M':
						hmm=new_hmm();
						// sscanf(argv[2], "%i", &level); 
						// printf("%s %s %i\n", argv[3], argv[4], level);
						read_obs(argv[2], hmm);
						read_hmm(argv[3], hmm);
						if (argc == 5) {
							sscanf(argv[4], "%i", &n_levels); 
							pi_n(hmm, n_levels, 0.0);
						} else {
							pi_n(hmm, hmm->n, 0.0);
						}
						logl_aic=forward_backward(hmm);
						printf("Log-l on observations: %lf\n", logl_aic);
						
						
						printf("Eigenstructure reveals %i modules.\n", hmm->n_modules);
						for(i=0;i<hmm->n;i++) {
							printf("State %i is in module %i\n", i, hmm->membership[i]);
						}
						printf("Module reconstruction on %s:\n", argv[2]);
						
						current_state=-1;
						st=viterbi_state(hmm);
						for(i=0;i<hmm->t;i++) { // for each of the steps
							printf("M%i ", hmm->membership[st[i]]);
						}
						printf("\n");
						
						break;
					case 'C':
						hmm=new_hmm();
						// sscanf(argv[2], "%i", &level); 
						// printf("%s %s %i\n", argv[3], argv[4], level);
						logl_aic=forward_backward(hmm);
						printf("Log-l on observations: %lf\n", logl_aic);

						read_obs(argv[2], hmm);
						read_hmm(argv[3], hmm);
						pi_n(hmm, hmm->n, 0.95);

						printf("Cutting at 20 step relaxation time.\n");
						printf("Eigenstructure reveals %i modules.\n", hmm->n_modules);
						for(i=0;i<hmm->n;i++) {
							printf("State %i is in module %i\n", i, hmm->membership[i]);
						}
						printf("Module reconstruction on %s:\n", argv[2]);

						current_state=-1;
						st=viterbi_state(hmm);
						for(i=0;i<hmm->t;i++) { // for each of the steps
							printf("M%i ", hmm->membership[st[i]]);
						}
						printf("\n");

						break;
					case 'm':
						hmm=new_hmm();
						read_obs(argv[2], hmm);
						read_hmm(argv[3], hmm);
						pi2(hmm);
						logl_aic=forward_backward(hmm);
						printf("Log-l on observations: %lf\n", logl_aic);
						
						printf("P(%c) in H module: %lf\n", hmm->sym_list[0], hmm->hi_count);
						printf("P(%c) in L module: %lf\n", hmm->sym_list[0], hmm->lo_count);
						
						printf("States in H module: ");
						for(i=0;i<hmm->hi_n;i++) {
							printf("%i ", hmm->hi_member[i]);
						}
						printf("\n");
						printf("States in L module: ");
						for(i=0;i<(hmm->n-hmm->hi_n);i++) {
							printf("%i ", hmm->lo_member[i]);
						}
						printf("\nReconstruction:\n");
						
						
						trap_h=0;
						trap_l=0;
						n_h=0;
						n_l=0;
							
						current_state=-1;
						st=viterbi_state(hmm);
						for(i=0;i<hmm->t;i++) { // for each of the steps
							for(j=0;j<hmm->hi_n;j++) { // then count up from hi_n
								if (hmm->hi_member[j] == st[i]) {
									printf("H ");
									if (current_state != 0) {
										n_h++;
										current_state=0;
									}
									trap_h++;
									break;
								}
							}
							if (j == hmm->hi_n) {
								printf("L ");
								if (current_state != 1) {
									n_l++;
									current_state=1;
								}
								trap_l++;
							}
						}
						printf("\n");
						if (n_h > 0) {
							printf("Trapping time in H: %lf steps\n", (double)trap_h/(double)n_h);
						} else {
							printf("No residency in H\n");
						}
						if (n_l > 0) {
							printf("Trapping time in L: %lf steps\n", (double)trap_l/(double)n_l);
						} else {
							printf("No residency in L\n");
						}
						free(st);
						break;
				}
			}
			
		} else {	
			sscanf(argv[2], "%i", &iter);
			if (argc == 5) {
				sscanf(argv[3], "%i", &start_states); 
				filename=argv[4];
				end_states=start_states;
			} else {
				filename=argv[3];
				start_states=1;
				end_states=24;
			}
			logl_aic=-1e32;
			best_states=start_states;
			pid=(unsigned long)getpid();
			filename_node=(char *)malloc(512*sizeof(char));
			
			for(nstates=start_states;nstates<=end_states;nstates++) {
				logl_list=(double *)malloc(iter*sizeof(double));			

				hmm=new_hmm();
				hmm->n=nstates;
				read_obs(filename, hmm);
						
				sprintf(filename_node, "%s_OUT_%istates", filename, nstates);
				// FIND BEST Log-L
				hbest=NULL;		
				
				fn = fopen("/dev/urandom", "rb"); 	
				r_seed=(unsigned long *)malloc(sizeof(unsigned long)*iter)	;
				if (fread(r_seed, sizeof(unsigned long)*iter, 1, fn) != 1) {
					/* Failed!--use time instead; beware, could sync with other instances */	
					printf("Warning: urandom read fail\n");
					exit(-1);
				}
				fclose(fn);
				
				for(in=0;in<iter;in++) {
					// printf("Doing iteration %i\n", in);
					
					init_hmm(hmm, r_seed[in]+pid); // mix in the process ID just in case two versions got the same urandom
					// printf("Inited\n");
					hmm=fit(hmm, 1e-3);
					// printf("Fit\n");
					forward_backward(hmm);
					// printf("Final fit\n");
					
					logl_list[in]=hmm->log_l;
					if ((hbest == NULL) || (hbest->log_l < hmm->log_l)) {
						// printf("Found better one %lf.\n", hmm->log_l);
						delete_hmm(hbest);
						hbest=copy_hmm(hmm);
					}					
				}
				// printf("Finished iterations at stage %i\n", nstates);
				
				pi(hbest);
				write_hmm(filename_node, hbest);
				// printf("Written hbest.\n");
				
				if ((int)argv[1][1] == 'f') {
					// REFINEMENT STEP
					delete_hmm(hmm);
				
					hmm=new_hmm();
					read_hmm(filename_node, hmm);
					read_obs(filename, hmm);
					hmm=fit(hmm, 1e-6);
					forward_backward(hmm);
					// now skipping refinements
					write_hmm(filename_node, hmm);
					// printf("Done refinements.\n");					
				}
				
				if (argc != 5) {
					printf("%i states: AIC %lf\n", nstates, hmm->log_l-hmm->n*(hmm->n-1)-hmm->n*(hmm->m-1));
				}	
				if (logl_aic < hmm->log_l-hmm->n*(hmm->n-1)-hmm->n*(hmm->m-1)) {
					logl_aic=hmm->log_l-hmm->n*(hmm->n-1)-hmm->n*(hmm->m-1);
				} else {
					best_states=nstates-1;
					nstates=end_states+1; // terminate loop
				}
				
				// Delete_hmm(hmm);
				delete_hmm(hmm);
				delete_hmm(hbest);
				free(logl_list);
				free(r_seed);
				// printf("Written out.\n");
			}
			
			// printf("Opening best states...\n");
			sprintf(filename_node, "%s_OUT_%istates", filename, best_states); // now go and open up the "best_states" hmm (by AIC)
			
			hmm=new_hmm();
			read_hmm(filename_node, hmm);	
			pi(hmm);
			if (argc != 5) {
				printf("\nBest number of states: %i\n", best_states);
			} else {
				printf("Found %i state model.\n", best_states);
			}
			printf("log-likelihood: %lf\n", hmm->log_l);	
			printf("lambda2: %lf\n", hmm->lambda2);
			printf("tau (relaxation time): %lf steps\n", 1.0/(1.0-hmm->lambda2));			
			printf("decay time: %lf steps\n", -1.0/log(hmm->lambda2));
			printf("\n");			
			printf("Elapsed CPU time: %18.16lf seconds.\n", (double)(clock() - t0)/CLOCKS_PER_SEC);
			delete_hmm(hmm);
		}
		
	}
	exit(1);
}