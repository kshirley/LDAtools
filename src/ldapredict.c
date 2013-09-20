#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <R.h>

int *ivector();
int **imatrix();
double *dvector();
double **dmatrix();
void free_ivector();
void free_imatrix();
void free_dvector();
void free_dmatrix();
void sample();
void dsum();

void ldapredict(int *word_id_R, int *doc_id_R, int *T_R, int *n_chains_R, int *n_iter_R, int *topics_init_R, 
	 double *alpha_R, double *loglik_R, int *N_R, int *W_R, int *D_R, double *phi_R) {

  // Initialize arguments passed in:
  int n_chains = n_chains_R[0];
  int n_iter = n_iter_R[0];
  int T = T_R[0];
  int N = N_R[0];
  int W = W_R[0];
  int D = D_R[0];
  int i, k, t, w, d, g, j, x;
  double alpha = alpha_R[0];
  //double beta = beta_R[0];
	  
  // Read in the word_id and doc_id vectors:
  int *word_id, *doc_id;
  word_id = ivector(N);
  doc_id = ivector(N);
  for (i = 0; i < N; i++) {
    word_id[i] = word_id_R[i];
    doc_id[i] = doc_id_R[i];
  }

  // Internal objects:
  int **z, *r, **CDT, *docsum;
  double **phi, *s, Talpha = T * alpha;
  double *prob_vec, *prob_norm, *n1, *n2;
  double loglik = 0.0, ll = 0.0;
	//int *topicsum, **CTW;

  // Set up a the count variables:
  //CTW = imatrix(T, W);
  CDT = imatrix(D, T);
  docsum = ivector(D);
  //topicsum = ivector(T);

  // Set up paramters:
  z = imatrix(n_chains, N);
  //theta = dmatrix(D, T);
  phi = dmatrix(T, W);
	
  // placeholders for the dirichlet draws:
  prob_vec = dvector(T);
  prob_norm = dvector(T);
  n1 = dvector(T);
  n2 = dvector(T);
  s = dvector(1);
  r = ivector(1);
	
  // Initialize the latent topics from input:
  int **topics_init;
  topics_init = imatrix(n_chains, N);
  for (k = 0; k < n_chains; k++) {
    for (i = 0; i < N; i++) {
      topics_init[k][i] = topics_init_R[i*n_chains + k];
    }
  }
  //Rprintf("%d %d %d\n", topics_init[0][0], topics_init[1][0], topics_init[2][0]);
  //Rprintf("%d %d %d\n", topics_init[0][9], topics_init[1][9], topics_init[2][9]);

	// Initialize the topic-token distributions from input:
  for (t = 0; t < T; t++) {
    for (w = 0; w < W; w++) {
      phi[t][w] = phi_R[w*T + t];
    }
  }
	
  // start the MCMC:
  for (k = 0; k < n_chains; k++) {
    Rprintf("\nChain %d\n", k + 1);
    // initialize counts:
    for (t = 0; t < T; t++) {
      //for (w = 0; w < W; w++) {
      //  CTW[t][w] = 0;
      //}	
      for (d = 0; d < D; d++) {
        CDT[d][t] = 0;
      }
    }
    for (d = 0; d < D; d++) docsum[d] = 0;
    //for (t = 0; t < T; t++) topicsum[t] = 0;
    // fill up count matrices and fill initial values of z:
    for (i = 0; i < N; i++) {
      z[k][i] = topics_init[k][i];
      //CTW[z[k][i] - 1][word_id[i] - 1]++;
      CDT[doc_id[i] - 1][z[k][i] - 1]++;
      //topicsum[z[k][i] - 1]++;
      docsum[doc_id[i] - 1]++;
    }
    loglik_R[0*n_chains + k] = 0;
    // Start with iteration 1 in the Gibbs sampler:		
    for (g = 1; g < n_iter; g++) {
      if ((g + 1) % 10 == 0) Rprintf("%d:%d\n", k + 1, g + 1);
      // Loop through tokens:
      loglik = 0.0;
      for (i = 0; i < N; i++) {
        // decrement counts:
        //CTW[z[k][i] - 1][word_id[i] - 1]--;
        CDT[doc_id[i] - 1][z[k][i] - 1]--;
        //topicsum[z[k][i] - 1]--;
				docsum[doc_id[i] - 1]--;
				// loop through topics:
				ll = 0.0;
        for (t = 0; t < T; t++) {
					// compute probability of each topic:
					//n1[t] = CTW[t][word_id[i] - 1] + beta;
					//prob_vec[t] = (n1[t] * n2[t])/((topicsum[t] + Wbeta) * (docsum[doc_id[i] - 1] + Talpha));
					n2[t] = CDT[doc_id[i] - 1][t] + alpha;
					prob_vec[t] = phi[t][word_id[i] - 1] * n2[t]/(docsum[doc_id[i] - 1] + Talpha);					
					ll += prob_vec[t];
				}
        // normalize probability vector for topic draw:
				s[0] = 0.0;
				dsum(prob_vec, T, s);
				for (t = 0; t < T; t++) prob_norm[t] = prob_vec[t]/s[0];
				// draw the topic for this token:
				r[0] = 0;
				sample(1, T, prob_norm, r);
				z[k][i] = r[0];
				// increment counts:
				//CTW[z[k][i] - 1][word_id[i] - 1]++;
				CDT[doc_id[i] - 1][z[k][i] - 1]++;
        //topicsum[z[k][i] - 1]++;
				docsum[doc_id[i] - 1]++;
				loglik += log(ll);
      } // i=N
      loglik_R[g*n_chains + k] = loglik;
    } // g=G
  } // k=n_chains
	
  // Set up the file to store the last iteration of z:
  for (k = 0; k < n_chains; k++) {
    for (i = 0; i < N; i++) {
      topics_init_R[i*n_chains + k] = z[k][i];
    }
  }

  // Free space:
  free_ivector(word_id);
  free_ivector(doc_id);

  //free_imatrix(CTW, T);
  free_imatrix(CDT, D);
  free_ivector(docsum);
  //free_ivector(topicsum);
  
  free_imatrix(z, n_chains);
  //free_dmatrix(theta, D);
  free_dmatrix(phi, T);

  free_dvector(prob_vec);
  free_dvector(prob_norm);
  free_dvector(n1);
  free_dvector(n2);
  free_dvector(s);
  free_ivector(r);

  free_imatrix(topics_init, n_chains);
}

