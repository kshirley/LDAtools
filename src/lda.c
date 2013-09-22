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

void lda(int *word_id_R, int *doc_id_R, int *T_R, int *n_chains_R, int *n_iter_R, int *topics_init_R, 
	 double *alpha_R, double *beta_R, double *loglik_R, int *N_R, int *W_R, int *D_R) {

  // Initialize arguments passed in:
  int n_chains = n_chains_R[0];
  int n_iter = n_iter_R[0];
  int T = T_R[0];
  int N = N_R[0];
  int W = W_R[0];
  int D = D_R[0];
  int i, k, t, w, d, g, j, x;
  double alpha = alpha_R[0];
  double beta = beta_R[0];
	  
  // Read in the word_id and doc_id vectors:
  int *word_id, *doc_id;
  word_id = ivector(N);
  doc_id = ivector(N);
  for (i = 0; i < N; i++) {
    word_id[i] = word_id_R[i];
    doc_id[i] = doc_id_R[i];
  }

  // Internal objects:
  int **z, *r, **CTW, **CDT, *topicsum, *docsum;
  double **theta, **phi, *s, Wbeta = W * beta, Talpha = T * alpha;
  double *prob_vec, *prob_norm, *n1, *n2, *llprob;
  double loglik = 0.0, ll = 0.0;

  // Set up a the count variables:
  CTW = imatrix(T, W);
  CDT = imatrix(D, T);
  docsum = ivector(D);
  topicsum = ivector(T);

  // Set up paramters:
  z = imatrix(n_chains, N);
  theta = dmatrix(D, T);
  phi = dmatrix(T, W);
	
  // placeholders for the dirichlet draws:
  prob_vec = dvector(T);
  llprob = dvector(T);
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
			
  // start the MCMC:
  for (k = 0; k < n_chains; k++) {
    Rprintf("\nChain %d\n", k + 1);
    // initialize counts:
    for (t = 0; t < T; t++) {
      for (w = 0; w < W; w++) {
        CTW[t][w] = 0;
      }	
      for (d = 0; d < D; d++) {
        CDT[d][t] = 0;
      }
    }
    for (d = 0; d < D; d++) docsum[d] = 0;
    for (t = 0; t < T; t++) topicsum[t] = 0;
    // fill up count matrices and fill initial values of z:
    for (i = 0; i < N; i++) {
      z[k][i] = topics_init[k][i];
      CTW[z[k][i] - 1][word_id[i] - 1]++;
      CDT[doc_id[i] - 1][z[k][i] - 1]++;
      topicsum[z[k][i] - 1]++;
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
        CTW[z[k][i] - 1][word_id[i] - 1]--;
        CDT[doc_id[i] - 1][z[k][i] - 1]--;
        topicsum[z[k][i] - 1]--;
	docsum[doc_id[i] - 1]--;
	// loop through topics:
	ll = 0.0;
        for (t = 0; t < T; t++) {
	  // compute probability of each topic:
	  n1[t] = CTW[t][word_id[i] - 1] + beta;
	  n2[t] = CDT[doc_id[i] - 1][t] + alpha;
	  prob_vec[t] = (n1[t] * n2[t])/((topicsum[t] + Wbeta) * (docsum[doc_id[i] - 1] + Talpha));
	  if (z[k][i] - 1 == t) {
	    llprob[t] = ((n1[t] + 1) * (n2[t] + 1))/((topicsum[t] + 1 + Wbeta) * (docsum[doc_id[i] - 1] + 1 + Talpha));
	  } else {
	    llprob[t] = ((n1[t]) * (n2[t]))/((topicsum[t] + Wbeta) * (docsum[doc_id[i] - 1] + 1 + Talpha));
	  }
	  ll += llprob[t];
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
	CTW[z[k][i] - 1][word_id[i] - 1]++;
	CDT[doc_id[i] - 1][z[k][i] - 1]++;
        topicsum[z[k][i] - 1]++;
	docsum[doc_id[i] - 1]++;
	loglik += log(ll);
      } // i=N
      loglik_R[g*n_chains + k] = loglik;
    } // g=G
  } // k=K
	
  // Set up the file to store the last iteration of z:
  for (k = 0; k < n_chains; k++) {
    for (i = 0; i < N; i++) {
      topics_init_R[i*n_chains + k] = z[k][i];
    }
  }

  // Free space:
  free_ivector(word_id);
  free_ivector(doc_id);

  free_imatrix(CTW, T);
  free_imatrix(CDT, D);
  free_ivector(docsum);
  free_ivector(topicsum);
  
  free_imatrix(z, n_chains);
  free_dmatrix(theta, D);
  free_dmatrix(phi, T);

  free_dvector(prob_vec);
  free_dvector(llprob);
  free_dvector(prob_norm);
  free_dvector(n1);
  free_dvector(n2);
  free_dvector(s);
  free_ivector(r);

  free_imatrix(topics_init, n_chains);
}

// compiling instructions
// from ldatools/src/: R CMD shlib lda.c


// Allocate a vector of integers
int *ivector(int n){
  int *v;
  v = (int *)malloc(n * sizeof(int));
  return v;
}

// Initialize a matrix of integers
int **imatrix(int nrow, int ncol){
  int **mat;
  int i;
  mat = (int **)malloc(nrow * sizeof(int *));
  for (i=0; i<nrow; i++){
    mat[i] = (int *)malloc(ncol * sizeof(int));
    }
  return mat;
}

// Initialize a vector of doubles
double *dvector(int n){
  double *v;
  v = (double *)malloc(n * sizeof(double));
  return v;
}

// Initialize a matrix of doubles
double **dmatrix(int nrow, int ncol){
  double **mat;
  int i;
  mat = (double **)malloc(nrow * sizeof(double *));   
  for (i=0; i<nrow; i++){
    mat[i] = (double *)malloc(ncol * sizeof(double));
  }
  return mat;
}


// Free an ivector
void free_ivector(int *v){
  free(v);
}

// Free an imatrix
void free_imatrix(int **v, int nrow){
  int i;
  for (i=0; i<nrow; i++) free(v[i]);
  free(v);
}

// Free a dvector
void free_dvector(double *v){
  free(v);
}

// Free a dmatrix
void free_dmatrix(double **v, int nrow){
  int i;
  for (i=0; i<nrow; i++) free(v[i]);
  free(v);
}

// Replicate the Sample function from R, with replacement always:
void sample(int first, int last, double *prob, int *z){
  int i, j, l = last - first + 1;
  double x, s;
  x = drand48();
  s = 0.0;
  for (j = 0; j < l; j++){
    s += prob[j];
    if (x < s){
      z[0] = j+first;
      break;
    }
  }
}

// Sum of a double vector:
void dsum(double *x, int n, double *z){
  int i;
  for (i=0; i<n; i++) z[0] += x[i];
}


