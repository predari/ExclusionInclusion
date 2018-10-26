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
#define GRID 5 // grid size
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

struct diagnostics {
  struct timeval earlier;
  struct timeval later;
  int svdPoints;
  double *ssv;

};

int main()  {

  lapack_complex_double *a = NULL;
  lapack_int m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT;
  int svdPoints = 0;
  int gsize = GRID;
  
  /* Time variables */
  struct timeval interval;

  struct diagnostics * diag = NULL;
  uint32_t * activity =  malloc((gsize*gsize)*sizeof( uint32_t));
  for (int i = 0; i < gsize*gsize ; i++)
    activity[i] = 0;

  a = malloc((lda*m)*sizeof(lapack_complex_double));
  create_grcar(a,m,lda);
  print_matrix("grcar",m,n,a,lda);	


  diag = pseudospectra(m, n, gsize, NEPSILON, XMIN, XMAX, YMIN, YMAX, a, activity);
  assert(activity);
  assert(diag);

  printf("Diagnostics:\n");
  timeval_diff(&interval,&diag->later,&diag->earlier);
  printf("- Time: (%ld seconds, %ld microseconds)\n",
	 interval.tv_sec,(long) interval.tv_usec);
  assert(diag->ssv);
  printf("- Ssv table: ");
  for (int i = 0; i < gsize * gsize ; i++) {
    if(!(i % gsize))
      printf("\n");
    printf("%.2f ",diag->ssv[i]);
  }  
  printf("\n");  
  
  printf("- Activity table: ");
  for (int i = 0; i < gsize * gsize ; i++) {
    if(!(i % gsize))
      printf("\n");
    printf("%d ",activity[i]);
  }  
  printf("\n");

  printf(" - Number of visited gridPoints: %d\n",diag->svdPoints);
  printf(" - Gain percentage: %d \n",
	 ((gsize*gsize-diag->svdPoints)/gsize*gsize)*100);

  free(a);
  free (activity);
  if (diag)
    free(diag->ssv);
  free(diag);
  return 0;
  
}



struct diagnostics * pseudospectra(lapack_int m, lapack_int n, lapack_int gsize,
			  lapack_int nbepsilon, double x_min, double x_max,
			  double y_min, double y_max,
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
  double ssv;

  lapack_complex_double *acp = NULL;
  assert(nbepsilon != 0);
  double e[nbepsilon];
  assert(a);
  assert(activity);

  /* Diagnostic variables */
  struct diagnostics * diag = malloc(sizeof(struct diagnostics));
  diag->ssv = malloc((gsize*gsize)*sizeof(double));


  for(int i = 0; i < nbepsilon; i++)
    e[i] = pow(0.1,(i+1));
  
  acp =  malloc((lda*m)*sizeof(lapack_complex_double));


  if(gettimeofday(&diag->earlier,NULL)) {
    perror("sixth gettimeofday()");
    exit(1);
  }
  
  for (int iy = 0; iy < gsize*gsize; iy++){
    // recreating original matrix, since zgesvd writes on it (acp)
    create_grcar(acp, m,lda);
    for (int i = 0; i < lda*m ; i=i+(n+1))
      // is it correct ? before it was a + x, but I don't understand why now
      acp[i]=acp[i]+(x_min+(iy/gsize * stepx)+(y_min + (iy % gsize * stepy))*I); // a + z 

    ssv = grid(m,n,nbepsilon,e,acp,gsize,activity,iy);
    svdPoints++;
    diag->ssv[iy] = ssv;
  }

  if(gettimeofday(&diag->later,NULL)) {
    perror("sixth gettimeofday()");
    exit(1);
  }
  diag->svdPoints = svdPoints;

  printf("activity (inside):");
  for (int i = 0; i < gsize * gsize ; i++) {
    if(!(i % gsize))
      printf("\n");
    printf("%d ",*(activity+i));
  }  
  printf("\n");
  
  free(acp);
  return diag;
  
}


double grid(lapack_int m, lapack_int n,
	lapack_int nbepsilon, double * e,
	lapack_complex_double *a,
	lapack_int gsize, uint32_t * activity,
	int iy) {

  lapack_int lda = n;
  lapack_int ldu = m;
  lapack_int ldvt = n;
  lapack_int info = 0;

  double *s;
  double *superb;
  s = malloc(m*sizeof(double)); //  the singular value results are stores here
  superb = malloc(min(m,n)*sizeof(double));

  assert(a);
  assert(activity);
  printf("Gridpoint %d:\n",iy);
  printf("---------------------\n");
  print_matrix("acp",m,n,a,lda);  
  
  /* Compute SVD */
  // for some reason zgesvd overwrites a even though
  // I have N N in the 2nd 3rd arguments
  // so I have to clear it or something like that.
  
  info = LAPACKE_zgesvd( LAPACK_ROW_MAJOR, 'N', 'N',
			 m, n, a, lda, s,
			 NULL, ldu, NULL, ldvt, superb );
  if( info > 0 ) {
    printf( "- SVD convergence: failed \n" );
    exit( 1 );
  }
  *(activity+iy) = 1;
  printf("smallest singular value is [%.2f,%.2f,%.2f ... %.2f,%.2f, %.2f] \n",
	 s[0],s[1],s[2],s[m-3],s[m-2],s[m-1]);

  free(s);
  free(superb);
  return s[m-1];

}
