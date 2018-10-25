#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <lapacke.h>
#include <assert.h>

#include "aux.h"

#define M 8 // first diminsion of matrix
#define N 8 // second diminsion of matrix
#define GRID 4 // grid size
#define LDA N
#define LDU M
#define LDVT N
#define MAXNEPSILON 32 // unit32_t used for the activity array
#define NEPSILON 1


// following is the input values of the domain.
// in a first attempt, we copied them from GK paper.
// later on we should use one of the following methods: Gresh. Disk, FoV, matrix norms, QR
// see notes on which one is more appropriate for us.

#define YMAX      3.41
#define YMIN     -3.41
#define XMAX      3.27
#define XMIN      -0.91                 

#define min(a,b) ((a)>(b)?(b):(a))

int main()  {

  lapack_complex_double *a;
  lapack_int m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT, info;
  int svdPoints = 0;
  int gsize = GRID;
  
  /* Time variables */
  struct timeval earlier;
  struct timeval later;
  struct timeval interval;

  uint32_t * activity =  malloc((gsize*gsize)*sizeof( uint32_t));
  for (int i = 0; i < gsize*gsize ; i++)
    activity[i] = 0;

  a = malloc((lda*m)*sizeof(lapack_complex_double));
  create_grcar(a,m,lda);
  print_matrix("grcar",m,n,a,lda);	

  svdPoints = grid(m, n, gsize, NEPSILON, XMIN, XMAX, YMIN, YMAX, a, activity);

  printf("final activity: ");
  for (int i = 0; i < gsize * gsize ; i++) {
    if(!(i % gsize))
      printf("\n");
    printf("%d ",activity[i]);
  }  
  printf("\n");

  printf("The number of visited gridPoints is %d\n",svdPoints);
  printf("Gain percentage on point-wise calculations is %d out of 100 \n",
	 ((gsize*gsize-svdPoints)/gsize*gsize)*100);

  free(a);
  free (activity);
  return 0;
  
}



int grid(lapack_int m, lapack_int n, lapack_int gsize,
		 lapack_int nbepsilon, double x_min, double x_max, double y_min, double y_max,
		lapack_complex_double *a, uint32_t * activity ) {

  lapack_int lda = n;
  lapack_int ldu = m;
  lapack_int ldvt = n;
  lapack_int info = 0;
  int svdPoints = 0;
  
  /* Initialize grid */
  double stepx = (x_max-x_min)/gsize;
  double stepy=(y_max-y_min)/gsize;
  printf("Domain size is: X=[%f-%f]  Y=[%f-%f] (stepx=%f, stepy=%f)\n", x_min, x_max,
	 y_min, y_max, stepx, stepy);
  printf("Grid size is: X=%d  Y=%d \n", gsize, gsize);

  double *s;
  double *superb;
  lapack_complex_double *tmp;
  int i, j, iy;
  // todoo: include assert.h
  assert(nbepsilon != 0);
  double e[nbepsilon];
  assert(a);
  assert(activity);

    /* Time variables */
  struct timeval earlier;
  struct timeval later;
  struct timeval interval;

  for(int i = 0; i < nbepsilon; i++)
    e[i] = pow(0.1,(i+1));
  
  tmp = calloc((lda*m)*sizeof(lapack_complex_double),0);
  s = malloc(m*sizeof(double)); //  the singular value results are stores here
  superb = malloc(min(m,n)*sizeof(double));
  assert(tmp);

  if(gettimeofday(&earlier,NULL)) {
    perror("sixth gettimeofday()");
    exit(1);
  }
  
  for (int iy = 0; iy < gsize*gsize; iy++){
    for (int i = 0; i < lda*m ; i=i+(n+1))
      tmp[i]=a[i]+(x_min+(iy/gsize * stepx)+(y_min + (iy % gsize * stepy))*I); // a + z
    
    /* Compute SVD */
    info = LAPACKE_zgesvd( LAPACK_ROW_MAJOR, 'N', 'N',
			   m, n, a, lda, s,
			   NULL, ldu, NULL, ldvt, superb );


    *(activity+iy) = 1;
    svdPoints++;  
  }

  if(gettimeofday(&later,NULL)) {
    perror("sixth gettimeofday()");
    exit(1);
  }
  if( info > 0 ) {
    printf( "The algorithm computing SVD failed to converge.\n" );
    exit( 1 );
  }

  timeval_diff(&interval,&later,&earlier);
  printf(" (%ld seconds, %ld microseconds)\n", interval.tv_sec,(long) interval.tv_usec);

  printf("The singular values are stored here:\n");
  for (i = 0; i < m; i++)
    printf("%.2f ", s[i]);
  printf("\n");


  printf("activity (inside):");
  for (int i = 0; i < gsize * gsize ; i++) {
    if(!(i % gsize))
      printf("\n");
    printf("%d ",*(activity+i));
  }  
  printf("\n");
  
  free(tmp);
  free(s);
  free(superb);
  return svdPoints;
  
}


