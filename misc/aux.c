#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#include <assert.h>

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


void create_grcar( lapack_complex_double *a, lapack_int m, lapack_int lda ) {
  int i,j;
 /*  for (i = 0; i < lda* m ; i++ ) { */
 /*    if( i % (lda+1) == 0 || i % (lda+1) == 1 || i % (lda+1) == 2 || i % (lda+1) == 3 )  */
 /*      *(a+i) = 1 ; */
 /*    else if( i % (lda+1) == lda ) */
 /*      *(a+i) = -1 ; */
 /*    else *(a+i) = 0; */
 /* } */
   for (i = 0; i < m ; i++ ) {
       *(a+(i+1)*lda+i) = -1 ;
     for (j = 0; j < lda ; j++ ) {
       if( j == i || j == i + 1  || j == i + 2 || j == i + 3 )
	 *(a+i*lda+j) = 1 ;
       else if (*(a+i*lda+j)!= -1)
	 *(a+i*lda+j) = 0;
     }
   }
}


void scalar_matrix_mult(lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda, int s) {
        lapack_int i, j;

        for( i = 0; i < m; i++ ) {
	  for( j = 0; j < n; j++ ) {
	    //	    lapack_complex_double_real(a[i*lda+j]) = lapack_complex_double_real(a[i*lda+j]) * s;
	    a[i*lda+j] = a[i*lda+j] * s;
	  }
	}
}


/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda ) {
        lapack_int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )  
                        printf( " (%6.2f,%6.2f)", lapack_complex_double_real(a[i*lda+j]),
				lapack_complex_double_imag(a[i*lda+j]) );
                printf( "\n" );
        }
}

void save_array(const char * fname, int m, int n, double *ssv ) {
  FILE * fp;

  fp = fopen (fname, "w+");
   for (int i = 0; i < m * n ; i++) {
    if(!(i % n) && i != 0)
      fprintf(fp,"\n");
    fprintf(fp,"%.4f ",ssv[i]);
   }  
   /* fprintf(fp,"\n");   */
   fclose(fp);
}


void read_array(const char * fname, int m, int n, double *ssv ) {
  FILE * fp;
  char buf[n];
  assert(ssv);
  fp = fopen (fname, "r");
  assert(fp);
  int i = 0;
  while((fscanf(fp,"%lf",&ssv[i]))!=EOF)
    i++;
  fclose(fp);
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
