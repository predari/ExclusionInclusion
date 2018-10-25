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
#define GRID 50 // grid size
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
  
  /* Initialize grid */
  stepx=(x_max-x_min)/(GRID);
  stepy=(y_max-y_min)/(GRID);
  printf("Domain size is: X=[%f-%f]  Y=[%f-%f] (stepx=%f, stepy=%f)\n", x_min, x_max, y_min, y_max, stepx, stepy);
  printf("Grid size is: X=%d  Y=%d \n", GRID,GRID);

  double *s;
  double *superb; 
  int i, j, iy;
  double e[MAXNEPSILON];
  uint32_t ** Activity;

  for(i = 0; i < NEPSILON; i++)
    e[i] = pow(0.1,(i+1));
  
  a = malloc((lda*m)*sizeof(lapack_complex_double));
  tmp = calloc((lda*m)*sizeof(lapack_complex_double),0);
  s = malloc(m*sizeof(double)); //  the singular value results are stores here
  superb = malloc(min(m,n)*sizeof(double));
  // allocate and initialize (with 0) the activity array
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

  printf("The number of visited gridPoints is %d\n",svdPoints);
  printf("Gain percentage on point-wise calculations is %d out of 100 \n",((GRID*GRID-svdPoints)/GRID*GRID)*100);
  
  free(a);
  free(tmp);
  free(s);
  free(superb);
  return 0;
  
}



uint32_t ** grid(lapack_int m, lapack_int n, lapack_int gsize,
		 lapack_int nbepsilon, double x_min, double x_max, double y_min, double y_max,
		 lapack_complex_double *a ) {

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
  uint32_t ** activity;
  assert(a);

    /* Time variables */
  struct timeval earlier;
  struct timeval later;
  struct timeval interval;

  for(int i = 0; i < nbepsilon; i++)
    e[i] = pow(0.1,(i+1));
  
  tmp = calloc((lda*m)*sizeof(lapack_complex_double),0);
  s = malloc(m*sizeof(double)); //  the singular value results are stores here
  superb = malloc(min(m,n)*sizeof(double));
  // allocate and initialize (with 0) the activity array
  if ((activity = malloc(gsize * sizeof(uint32_t *))) == NULL) {
    fprintf(stderr,"Failed to malloc space for the data\n");
    exit(-1);
  }
  for (int i = 0; i < gsize ; i++) {
    if ((activity[i] = calloc(gsize *sizeof(uint32_t),0)) == NULL) {
      fprintf(stderr,"Failed to malloc space for the data\n");
      exit(-1);
    }
  }
  
  assert(activity);
  assert(tmp);

  //create_grcar(a,m,lda);
  //print_matrix("grcar",m,n,a,lda);	

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

  return activity;
  
}
