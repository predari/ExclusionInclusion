#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <lapacke.h>

#include "aux.h"

#define M 8
#define N 8
#define GRID 50
#define LDA N
#define LDU M
#define LDVT N
#define MAXNEPSILON 32 // unit32_t used for the activity array
#define NEPSILON 1

#define YMAX      3.41
#define YMIN     -3.41
#define XMAX      3.27
#define XMIN      -0.91                 

/* COMPILE: gcc -o t test.c aux.c -lm -llapacke -lblas -L/usr/local/opt/lapack/lib (add lapack dir) -I/usr/local/Cellar/lapack/3.8.0/include (include path for lapacke) */
/* -L/usr/local/lib (link with blas) - done implicitly when adding dir for lapack*/

#define min(a,b) ((a)>(b)?(b):(a))

int main()  {
        /* not Locals */
  lapack_complex_double *a, *tmp;
  lapack_int m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT, info;
  int svdPoints = 0;
  
  /* Time variables */
  struct timeval earlier;
  struct timeval later;
  struct timeval interval;
  
  double  x_min, x_max,
    y_min, y_max,
    stepx, stepy;
  /*if (x_min==0.0)*/  x_min=XMIN;
  /*if (x_max==0.0)*/  x_max=XMAX;
  /*if (y_min==0.0)*/  y_min=YMIN;
  /*if (y_max==0.0)*/  y_max=YMAX;

  /* Initialize grid with random initial variables */
  stepx=(x_max-x_min)/(GRID);
  stepy=(y_max-y_min)/(GRID);
  printf("Domain size is: X=[%f-%f]  Y=[%f-%f] (stepx=%f, stepy=%f)\n", x_min, x_max, y_min, y_max, stepx, stepy);

  double *s;
  double *superb; 
  int i, j, iy;
  double e[MAXNEPSILON];
  uint32_t ** Activity;

  for(i = 0; i < NEPSILON; i++)
    e[i] = pow(0.1,(i+1));
  
  a = malloc((lda*m)*sizeof(lapack_complex_double));
  tmp = calloc((lda*m)*sizeof(lapack_complex_double),0);
  s = malloc(m*sizeof(double)); //  the svd result is here
  superb = malloc(min(m,n)*sizeof(double));
  if ((Activity = malloc(GRID*sizeof(uint32_t *))) == NULL) {
    fprintf(stderr,"Failed to malloc space for the data\n");
    exit(-1);
  }
  for (i=0;i<GRID;i++) {
    if ((Activity[i] = calloc(GRID*sizeof(uint32_t),0)) == NULL) {
      fprintf(stderr,"Failed to malloc space for the data\n");
      exit(-1);
    }
  }


  create_grcar(a,m,lda);
  print_matrix("grcar",m,n,a,lda);	

  if(gettimeofday(&earlier,NULL)) {
    perror("sixth gettimeofday()");
    exit(1);
  }
	
  for (iy = 0; iy < GRID*GRID; iy++){
    //  memcpy(tmp, a ,(lda*m)*sizeof(*tmp));
    // tmp is not right!!!
    for (i = 0; i < lda*m ; i=i+(n+1))
      tmp[i]=a[i]+(x_min+(iy/GRID * stepx)+(y_min + (iy%GRID * stepy))*I); // a + z
    
    /* Compute SVD */
    info = LAPACKE_zgesvd( LAPACK_ROW_MAJOR, 'N', 'N',
			   m, n, a, lda, s,
			   NULL, ldu, NULL, ldvt, superb );
    svdPoints++;
    
  }
  
  if(gettimeofday(&later,NULL)) {
    perror("sixth gettimeofday()");
    exit(1);
  }
  
  timeval_diff(&interval,&later,&earlier);
  printf(" (%ld seconds, %ld microseconds)\n", interval.tv_sec,(long) interval.tv_usec);
  

  if( info > 0 ) {
    printf( "The algorithm computing SVD failed to converge.\n" );
    exit( 1 );
  }
  printf("The singular values are stored here:\n");
  for (i = 0; i < m; i++)
    printf("%.2f ", s[i]);
  printf("\n");
  
  free(a);
  free(tmp);
  free(s);
  free(superb);
  return 0;
  
}
