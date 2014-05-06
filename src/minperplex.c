#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>

int *ivector();
void free_ivector();

void minperplex(int *token_id, int *N_R, int *W_R, int *D_R, double *alpha_R, double *beta_R, 
								double *logp, int *tokenfreq, int *docfreq, int *print_R) {
	
  // Initialize arguments passed in:
  int N = N_R[0];
  int W = W_R[0];
  int D = D_R[0];
  int d, i, w, seen;
	int print = print_R[0];
  double alpha = alpha_R[0];
  double beta = beta_R[0];
	double p, theta, phi;
  int *doctable;
	
	// Allocate and initialize doctable:
	doctable = ivector(W);
  for (w = 0; w < W; w++) doctable[w] = 0;

	seen = 0; // keep track of how many tokens have been seen
  // loop through documents
	for (d = 0; d < D; d++) {
		if ((d + 1) % print == 0) Rprintf("Document %d\n", d + 1);
    // increment the table of tokens for this document
  	for (i = 0; i < docfreq[d]; i++) doctable[token_id[seen + i] - 1]++;
		for (i = 0; i < docfreq[d]; i++) {
			p = 0.0;  // the probability of token i
			for (w = 0; w < W; w ++) {
				theta = (doctable[w] + alpha)/(docfreq[d] + alpha * W);
			  if (w == token_id[seen + i] - 1) {
				  phi = (tokenfreq[w] + beta)/(tokenfreq[w] + W * beta);
				} else {
				  phi = beta/(tokenfreq[w] + W * beta);				
				}
				p += theta * phi;
			}
			logp[seen + i] = log(p);
		}
    // decrement the token counts to set to zero for next document
  	for (i = 0; i < docfreq[d]; i++) {
	    doctable[token_id[seen + i] - 1]--;
		}		
    // set 'seen' to point to first token in the next document
		seen += docfreq[d];
  }
  // free the space allocated to doctable
  free_ivector(doctable);
}


