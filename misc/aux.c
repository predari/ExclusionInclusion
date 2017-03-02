
#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#include "aux.h"


/* Auxiliary routine for timing in miscoseconds */
long long timeval_diff(struct timeval *difference, struct timeval *end_time, struct timeval *start_time) {
  
  struct timeval temp_diff;
  if(difference == NULL) {
	difference = &temp_diff;
  }

  difference->tv_sec = end_time->tv_sec - start_time->tv_sec ;
  difference->tv_usec = end_time->tv_usec - start_time->tv_usec;
  
  while(difference->tv_usec < 0) {
    difference->tv_usec+= 1000000;
    difference->tv_sec -= 1;
  }

  return 1000000LL*difference->tv_sec+
                   difference->tv_usec;

} /* timeval_diff() */


/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda ) {
        lapack_int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )  
                        printf( " (%6.2f,%6.2f)", lapack_complex_double_real(a[i*lda+j]), lapack_complex_double_imag(a[i*lda+j]) );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a vector of integers */
void print_int_vector( char* desc, lapack_int n, lapack_int* a ) {
        lapack_int j;
        printf( "\n %s\n", desc );
        for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
        printf( "\n" );
}


 void prtdat(int nx, int ny, double *u1, char *fnam) {

   int ix;
   FILE *fp;
   // 
 fp = fopen(fnam, "w");
 //fprintf(fp, "%d %d\n", nx, ny);
   for (ix = 0; ix < nx*ny; ix++) {
     fprintf(fp, "%6.10f", u1[ix]);
     if (ix % nx == nx-1) 
       fprintf(fp, "\n");
     else
       fprintf(fp, " ");
     }
 fclose(fp);
 }
