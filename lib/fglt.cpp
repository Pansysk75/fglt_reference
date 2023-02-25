/*!
  \file   fglt.cpp
  \brief  Fast Graphlet Transform source file

  \author Dimitris Floros
  \date   2020-08-18
*/

#ifdef HAVE_CONFIG_H
  #include "../config.h"
#endif

#ifdef HAVE_CILK_CILK_H
  #include <cilk/cilk.h>
  #include <cilk/cilk_api.h>
  #define FOR cilk_for
  #ifdef CILKSCALE
    #include <cilk/cilkscale.h>
  #endif
#else
  #define FOR for
#endif

#include "fglt.hpp"
#include <algorithm>

struct timeval tic(){
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv;
}
  
static double toc(struct timeval begin){
  struct timeval end;
  gettimeofday(&end, NULL);
  double stime = ((double) (end.tv_sec - begin.tv_sec) * 1000 ) +
    ((double) (end.tv_usec - begin.tv_usec) / 1000 );
  stime = stime / 1000;
  return(stime);
}



extern "C" int getWorkers(){
#ifdef HAVE_CILK_CILK_H
  return __cilkrts_get_nworkers();
#else
  return 1;
#endif
}

static void remove_neighbors
(
 int *isNgbh,
 mwIndex i,
 mwIndex *ii,
 mwIndex *jStart
){

  // --- remove neighbors
  for (mwIndex id_i = jStart[i]; id_i < jStart[i+1]; id_i++){

    // get the column (k)
    mwIndex k = ii[id_i];

    isNgbh[k] = 0;
  }
  
}


static void raw2net
(
 double        ** const f,
 double  const ** const d,
 mwIndex const          i
){
  f[ 0][i]  =  d[ 0][i];
  f[ 1][i]  =  d[ 1][i];
  f[ 2][i]  =  d[ 2][i] -  2 * d[ 4][i];
  f[ 3][i]  =  d[ 3][i] -      d[ 4][i];
  f[ 4][i]  =  d[ 4][i];
  f[ 5][i]  =  d[ 5][i] -  2 * d[ 9][i] -      d[10][i] -  2 * d[12][i] +  4 * d[13][i] +  2 * d[14][i] -  6 * d[15][i];
  f[ 6][i]  =  d[ 6][i] -      d[10][i] -  2 * d[11][i] -  2 * d[12][i] +  2 * d[13][i] +  4 * d[14][i] -  6 * d[15][i];
  f[ 7][i]  =  d[ 7][i] -      d[ 9][i] -      d[10][i] +  2 * d[13][i] +      d[14][i] -  3 * d[15][i];
  f[ 8][i]  =  d[ 8][i] -      d[11][i] +      d[14][i] -      d[15][i];
  f[ 9][i]  =  d[ 9][i] -  2 * d[13][i] +  3 * d[15][i];
  f[10][i]  =  d[10][i] -  2 * d[13][i] -  2 * d[14][i] +  6 * d[15][i];
  f[11][i]  =  d[11][i] -  2 * d[14][i] +  3 * d[15][i];
  f[12][i]  =  d[12][i] -      d[13][i] -      d[14][i] +  3 * d[15][i];
  f[13][i]  =  d[13][i] -  3 * d[15][i];
  f[14][i]  =  d[14][i] -  3 * d[15][i];
  f[15][i]  =  d[15][i];
}


static void compute_all_available
(
 double **f,
 mwIndex i
){

  f[0][i]   = 1;
  f[3][i]   = f[1][i] * ( f[1][i] - 1 ) * 0.5;
}


static void spmv_first_pass
(
 double *f2_i,
 double *f1,
 mwIndex i,
 mwIndex *jStart,
 mwIndex *ii
){

  // --- loop through every nonzero element A(i,k)
  for (mwIndex id_i = jStart[i]; id_i < jStart[i+1]; id_i++){

    // get the column (k)
    mwIndex k = ii[id_i];
      
    // --- matrix-vector products
    f2_i[0] += f1[k];

  }

  f2_i[0]  -= f1[i];
  
}

static void p2
(
 double *f4_i,
 double *c3,
 mwIndex i,
 mwIndex *jStart,
 mwIndex *ii,
 double *fl,
 int *pos,
 int *isNgbh,
 mwIndex *isUsed
){

  // setup the count of nonzero columns (j) visited for this row (i)
  mwIndex cnt = 0;

  // --- loop through every nonzero element A(i,k)
  for (mwIndex id_i = jStart[i]; id_i < jStart[i+1]; id_i++){

    // get the column (k)
    mwIndex k = ii[id_i];

    isNgbh[k] = id_i+1;
      
    // --- loop through all nonzero elemnts A(k,j)
    for (mwIndex id_k = jStart[k]; id_k < jStart[k+1]; id_k++){

      // get the column (j)
      mwIndex j = ii[id_k];

      if (i == j) continue;

      // if this column is not visited yet for this row (i), then set it
      if (!isUsed[j]) {
        fl[j]      = 0;  // initialize corresponding element
        isUsed[j]  = 1;  // set column as visited
        pos[cnt++] = j;  // add column position to list of visited
      }

      // increase count of A(i,j)
      fl[j]++;
        
    }

  }

  // --- perform reduction on [cnt] non-empty columns (j) 
  for (mwIndex l=0; l<cnt; l++) {

    // get next column number (j)
    mwIndex j = pos[l];


    if (isNgbh[j]) {
      c3[isNgbh[j]-1]  = fl[j];
        
      f4_i[0]  += fl[j];
    }
      
    // declare it non-used
    isUsed[j] = 0;
  }

  f4_i[0]  /= 2;
    
}


extern "C" int compute
(
 double ** const f,
 double ** const fn,
 mwIndex *ii,
 mwIndex *jStart,
 mwSize n,
 mwSize m,
 mwSize np
 ){

  struct timeval timer_all = tic();

#ifdef CILKSCALE
  wsp_t cs_start, cs_end;
  cs_start = wsp_getworkspan();
#endif



  FOR (mwSize i=0;i<n;i++) {
    // get degree of vertex (i)
    f[1][i] = jStart[i+1] - jStart[i];
  }
  
  
  // --- setup auxilliary vectors (size n)
  double *fl = (double *) calloc( n*np, sizeof(double) );
  int *pos = (int *) calloc( n*np, sizeof(int) );
  mwIndex *isUsed    = (mwIndex *) calloc( n*np, sizeof(mwIndex) );
  
  double *c3 = (double *) calloc( m, sizeof(double) );
  int *isNgbh = (int *) calloc( n*np, sizeof(int) );

  if ( fl == NULL || pos == NULL || isUsed == NULL || c3 == NULL || isNgbh == NULL ){
    printf( "Working memory allocation failed at auxilliary vectors, aborting...\n" );
    return 1;
  }
  
  // --- first pass
  FOR (mwIndex i=0; i<n;i++) {
#ifdef HAVE_CILK_CILK_H
    int ip = __cilkrts_get_worker_number();
#else
    int ip = 0;
#endif


    // d_4 d_10 d_12 d_14
    p2(  &f[4][i],
         c3, i, jStart, ii,
         &fl[ip*n], &pos[ip*n], &isNgbh[ip*n], &isUsed[ip*n] );

    
    // d_2 d_7
    spmv_first_pass( &f[2][i], f[1], i, jStart, ii );

    
    // d_3 d_6 d_8 d_11
    compute_all_available(f, i);
    remove_neighbors(&isNgbh[ip*n], i, ii, jStart);
    
  }

  free(fl);
  free(pos);
  free(isUsed);

  FOR (mwIndex i=0; i<n;i++){
    
    // transform to net
    raw2net( (double ** const) fn, (double const ** const) f, i );

    
  }
  
  free(isNgbh);
  free(c3);

#ifdef CILKSCALE
  cs_end = wsp_getworkspan();
#endif

  printf( "Total elapsed time: %.4f sec\n", toc( timer_all ) );

#ifdef CILKSCALE
  wsp_dump( wsp_sub( cs_end, cs_start ), "FGLT computation" );
#endif
  
  return 0;
  
}
