#include "mex.h"
#include "cokus.cpp"
#include <stdio.h>

//reference :Griffiths, T., & Steyvers, M. (2004).  Finding Scientific Topics. Proceedings of the National Academy of Sciences, 101 (suppl. 1), 5228-5235.

void CDSeqGibbsSampler( double ALPHA, double BETA, int totalGenes, int T, int samplesize, int NIter, int n, int *z, int *d, int *w, int *csGEPcount, int *SSPcount, int *ztot, int *order, double *probs, int OUTPUT)
{
  int wi,di,i,ii,j,topic, rp, temp, iter, wioffset, dioffset;
  double totprob, WBETA, r, max;

      for (i=0; i<n; i++)
      {
          wi = w[ i ];
          di = d[ i ];
          // pick a random topic 0..T-1
          topic = (int) ( (double) randomMT() * (double) T / (double) (4294967296.0 + 1.0) );
          z[ i ] = topic; // assign read to cell type
          csGEPcount[ wi*T + topic ]++; // increment csGEPcount count matrix
          SSPcount[ di*T + topic ]++; // increment SSPcount count matrix
          ztot[ topic ]++; // increment ztot matrix
      }
  
  
  for (i=0; i<n; i++) order[i]=i; // fill with increasing series
  for (i=0; i<(n-1); i++) {
      // pick a random integer between i and nw
      rp = i + (int) ((double) (n-i) * (double) randomMT() / (double) (4294967296.0 + 1.0));
      
      // switch contents on position i and position rp
      temp = order[rp];
      order[rp]=order[i];
      order[i]=temp;
  }
  
  WBETA = (double) (totalGenes*BETA);
  for (iter=0; iter<NIter; iter++) {
      if (OUTPUT >=1) {
          if ((iter % 10)==0) mexPrintf( "\t %d of %d MCMC iterations\n" , iter , NIter );
          if ((iter % 10)==0) mexEvalString("drawnow;");
      }
      for (ii = 0; ii < n; ii++) {
          i = order[ ii ]; // current read to assess
          
          wi  = w[i]; // current read index
          di  = d[i]; // current sample index  
          topic = z[i]; // current cell type assignment to read
          ztot[topic]--;  // substract this from counts
          
          wioffset = wi*T;
          dioffset = di*T;
          
          csGEPcount[wioffset+topic]--;
          SSPcount[dioffset+topic]--;
          
          totprob = (double) 0;
          for (j = 0; j < T; j++) {
              probs[j] = ((double) csGEPcount[ wioffset+j ] + (double) BETA)/( (double) ztot[j]+ (double) WBETA)*( (double) SSPcount[ dioffset+ j ] + (double) ALPHA);
              totprob += probs[j];
          }
          
          // sample a topic from the distribution
          r = (double) totprob * (double) randomMT() / (double) 4294967296.0;
          max = probs[0];
          topic = 0;
          while (r>max) {
              topic++;
              max += probs[topic];
          }
           
          z[i] = topic; // assign current word token i to topic j
          csGEPcount[wioffset + topic ]++; // and update counts
          SSPcount[dioffset + topic ]++;
          ztot[topic]++;         
      }
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  double *srwp, *srdp, *probs, *Z, *GID, *SID, *ZIN;
  double ALPHA,BETA;
  mwIndex *irwp, *jcwp, *irdp, *jcdp;
  int *z,*d,*w, *order, *csGEPcount, *SSPcount, *ztot;
  int totalGenes,T,samplesize,NIter,SEED,OUTPUT, nzmax, nzmaxwp, nzmaxdp, ntokens;
  int i,j,c,n,nt,wi,di;//, startcond;
  
  /* Check for proper number of arguments. */
  if (nrhs < 8) {
    mexErrMsgTxt("At least 8 input arguments required by Gibbs sampler");
  } else if (nlhs < 2) {
    mexErrMsgTxt("2 output arguments required by Gibbs sampler");
  }
  
  
  /* process the input arguments */
  if (mxIsDouble( prhs[ 0 ] ) != 1) mexErrMsgTxt("GID input vector must be a double precision matrix");
  if (mxIsDouble( prhs[ 1 ] ) != 1) mexErrMsgTxt("SID input vector must be a double precision matrix");
  
  // pointer to gene indices
  GID = mxGetPr( prhs[ 0 ] );
     
  // pointer to sample indices
  SID = mxGetPr( prhs[ 1 ] );
  
  // get the number of reads
  ntokens = mxGetM( prhs[ 0 ] ) * mxGetN( prhs[ 0 ] );
  
  if (ntokens == 0) mexErrMsgTxt("GID vector is empty"); 
  if (ntokens != ( mxGetM( prhs[ 1 ] ) * mxGetN( prhs[ 1 ] ))) mexErrMsgTxt("GID and SID vectors should have same number of entries");
  
  T    = (int) mxGetScalar(prhs[2]);
  if (T<=0) mexErrMsgTxt("Number of cell types must be greater than zero");
  
  NIter    = (int) mxGetScalar(prhs[3]);
  if (NIter<0) mexErrMsgTxt("Number of iterations must be positive");
  
  ALPHA = (double) mxGetScalar(prhs[4]);
  if (ALPHA<=0) mexErrMsgTxt("ALPHA must be greater than zero");
  
  BETA = (double) mxGetScalar(prhs[5]);
  if (BETA<=0) mexErrMsgTxt("BETA must be greater than zero");
  
  SEED = (int) mxGetScalar(prhs[6]);
  
  OUTPUT = (int) mxGetScalar(prhs[7]);
  
  // seeding
  seedMT( 1 + SEED * 2 ); // seeding only works on uneven numbers
  
   
  
  /* allocate memory */
  z  = (int *) mxCalloc( ntokens , sizeof( int ));
  
  d  = (int *) mxCalloc( ntokens , sizeof( int ));
  w  = (int *) mxCalloc( ntokens , sizeof( int ));
  order  = (int *) mxCalloc( ntokens , sizeof( int ));  
  ztot  = (int *) mxCalloc( T , sizeof( int ));
  probs  = (double *) mxCalloc( T , sizeof( double ));
  
  for (i=0; i<ntokens; i++) {
     w[ i ] = (int) GID[ i ] - 1;
     d[ i ] = (int) SID[ i ] - 1;
  }
  
  n = ntokens;
  
  totalGenes = 0;
  samplesize = 0;
  for (i=0; i<n; i++) {
     if (w[ i ] > totalGenes) totalGenes = w[ i ];
     if (d[ i ] > samplesize) samplesize = d[ i ];
  }
  totalGenes = totalGenes + 1;
  samplesize = samplesize + 1;
  
  csGEPcount  = (int *) mxCalloc( T*totalGenes , sizeof( int ));
  SSPcount  = (int *) mxCalloc( T*samplesize , sizeof( int ));
   
  
  /* run the model */
  CDSeqGibbsSampler( ALPHA, BETA, totalGenes, T, samplesize, NIter, n, z, d, w, csGEPcount, SSPcount, ztot, order, probs,OUTPUT );
  
  /* convert the full csGEPcount matrix into a sparse matrix */
  nzmaxwp = 0;
  for (i=0; i<totalGenes; i++) {
     for (j=0; j<T; j++)
         nzmaxwp += (int) ( *( csGEPcount + j + i*T )) > 0;
  }  
  
  // MAKE THE WP SPARSE MATRIX
  plhs[0] = mxCreateSparse( totalGenes,T,nzmaxwp,mxREAL);
  srwp  = mxGetPr(plhs[0]);
  irwp = mxGetIr(plhs[0]);
  jcwp = mxGetJc(plhs[0]);  
  n = 0;
  for (j=0; j<T; j++) {
      *( jcwp + j ) = n;
      for (i=0; i<totalGenes; i++) {
         c = (int) *( csGEPcount + i*T + j );
         if (c >0) {
             *( srwp + n ) = c;
             *( irwp + n ) = i;
             n++;
         }
      }    
  }  
  *( jcwp + T ) = n;    
   
  // MAKE THE DP SPARSE MATRIX
  nzmaxdp = 0;
  for (i=0; i<samplesize; i++) {
      for (j=0; j<T; j++)
          nzmaxdp += (int) ( *( SSPcount + j + i*T )) > 0;
  }  
  
  plhs[1] = mxCreateSparse( samplesize,T,nzmaxdp,mxREAL);
  srdp  = mxGetPr(plhs[1]);
  irdp = mxGetIr(plhs[1]);
  jcdp = mxGetJc(plhs[1]);
  n = 0;
  for (j=0; j<T; j++) {
      *( jcdp + j ) = n;
      for (i=0; i<samplesize; i++) {
          c = (int) *( SSPcount + i*T + j );
          if (c >0) {
              *( srdp + n ) = c;
              *( irdp + n ) = i;
              n++;
          }
      }
  }
  *( jcdp + T ) = n;
  
  plhs[ 2 ] = mxCreateDoubleMatrix( 1,ntokens , mxREAL );
  Z = mxGetPr( plhs[ 2 ] );
  for (i=0; i<ntokens; i++) Z[ i ] = (double) z[ i ] + 1;
}
