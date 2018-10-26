#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <lapacke.h>
#include <assert.h>

#include "aux.h"

#define M 20 // first diminsion of matrix
#define N 20 // second diminsion of matrix
#define GRID 10 // grid size
#define LDA N
#define LDU M
#define LDVT N
#define MAXNEPSILON 32 // unit32_t used for the activity array
#define NEPSILON 3


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

struct domain {
  double x_min;
  double x_max;
  double y_min;
  double y_max;
  double stepx;
  double stepy;
};


int main()  {
  // TODO: 
  /* USER INPUT PARAMETERS */
  /*   - NBEPSILON */
  /*   - m, n */
  /*   - GRID */
  /*   - a */
  /*   - XMIN, XMAX, YMIN, YMAX */
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
  scalar_matrix_mult (m, n, a, lda, -1);
  print_matrix("- grcar",m,n,a,lda);	

  /* Initialize grid */
  struct domain dm;
  dm.x_min = XMIN;
  dm.x_max = XMAX;
  dm.y_min = YMIN;
  dm.y_max = YMAX;
  dm.stepx = (dm.x_max - dm.x_min)/gsize;
  dm.stepy=(dm.y_max - dm.y_min)/gsize;
   
  diag = pseudospectra(m, n, gsize, NEPSILON, &dm, a, activity);
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
				   lapack_int nbepsilon, struct domain * dm,
				   lapack_complex_double *a, uint32_t * activity ) {

  lapack_int lda = n;
  lapack_int ldu = m;
  lapack_int ldvt = n;
  lapack_int info = 0;
  int svdPoints = 0;
  


  printf("Domain size is: X=[%f-%f]  Y=[%f-%f] (stepx=%f, stepy=%f)\n",
	 dm->x_min, dm->x_max,
	 dm->y_min, dm->y_max,
	 dm->stepx, dm->stepy);
  printf("Grid size is: X=%d  Y=%d \n", gsize, gsize);

  double ssv;

  lapack_complex_double *acp = NULL;
  assert(nbepsilon != 0);
  // todo: maybe e should be created in the main program
  double e[nbepsilon];
  assert(a);
  assert(activity);

  /* Diagnostic variables */
  struct diagnostics * diag = malloc(sizeof(struct diagnostics));
  diag->ssv = malloc((gsize*gsize)*sizeof(double));

  printf("e curves: ");
  for(int i = 0; i < nbepsilon; i++) {
    e[i] = pow(0.1,(i+1));
    printf(" %f ", e[i]);
  }
  printf("\n");
  acp =  malloc((lda*m)*sizeof(lapack_complex_double));


  if(gettimeofday(&diag->earlier,NULL)) {
    perror("sixth gettimeofday()");
    exit(1);
  }

  //printf("Pseudospectra with grid\n");
  printf("Pseudospectra with mog\n");  
  for (int iy = 0; iy < gsize*gsize; iy++){
    // recreating original matrix, since zgesvd writes on it (acp)
    create_grcar(acp, m,lda);
    scalar_matrix_mult (m, n, acp, lda, -1);
    //print_matrix("- grcar",m,n,a,lda);	

    for (int i = 0; i < lda*m ; i=i+(n+1))
      acp[i]=acp[i]+(dm->x_min+(iy/gsize * dm->stepx)+(dm->y_min + (iy % gsize * dm->stepy))*I); // -a + z 

    //ssv = grid(m,n,nbepsilon,e,acp,gsize,activity,iy);

    ssv = mog(m,n,nbepsilon,e,dm,acp,gsize,activity,iy);
    svdPoints++;
    diag->ssv[iy] = ssv;
  }

  if(gettimeofday(&diag->later,NULL)) {
    perror("sixth gettimeofday()");
    exit(1);
  }
  diag->svdPoints = svdPoints;

  /* printf("activity (inside):"); */
  /* for (int i = 0; i < gsize * gsize ; i++) { */
  /*   if(!(i % gsize)) */
  /*     printf("\n"); */
  /*   printf("%d ",*(activity+i)); */
  /* }   */
  /* printf("\n"); */
  
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
  s = malloc(m*sizeof(double)); // holds singular values
  superb = malloc(min(m,n)*sizeof(double));
  assert(a);
  assert(activity);
  
  /* printf("Gridpoint %d:\n",iy); */
  /* printf("---------------------\n"); */
  /* print_matrix("acp",m,n,a,lda);   */  
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
  printf("Point:%d --> smallest singular value is [%.2f,%.2f,%.2f ... %.2f,%.2f, %.2f] \n",
	 iy,s[0],s[1],s[2],s[m-3],s[m-2],s[m-1]);

  free(superb);
  return s[m-1];

}



double mog(lapack_int m, lapack_int n,
	   lapack_int nbepsilon, double * e,
	   struct domain * dm,
	   lapack_complex_double *a,
	   lapack_int gsize, uint32_t * activity,
	   int z) {

  lapack_int lda = n;
  lapack_int ldu = m;
  lapack_int ldvt = n;
  lapack_int info = 0;

  double *s;
  double *superb;
  s = malloc(m*sizeof(double)); // holds the singular values
  superb = malloc(min(m,n)*sizeof(double));
  double stepx = 0;
  double stepy = 0;
  int ir = 0;
  int jr = 0;
  int il = gsize;
  int jl = gsize;
  int start_i = 0;
  int start_j = 0;
  int end_i = gsize;
  int end_j = gsize;
  // todo: i indexes lines of matrix so [i-stepy, i, i+stepy] changes lines
  // todo: j indexes columns of matrix so [j-stepx, j, j+stepx] changes columns 
  
  int i = z/gsize;
  int j = z%gsize;
  assert(a);
  assert(activity);

  /* printf("Gridpoint %d:\n",z); */
  /* printf("---------------------\n"); */
  /* print_matrix("acp",m,n,a,lda);   */
  /* Compute SVD */
  // for some reason zgesvd overwrites a even though
  // I have N N in the 2nd 3rd arguments
  // so I have to clear it or something like that.
  //printf("Pseudospectra with mog\n");
  if(  *(activity+z) == 1) {
    printf("Point:%d is skipped!\n",z);
    return 0;
    // todo return a structure for that info. boolean skipped or not and the value of ssv when possible
  }
  info = LAPACKE_zgesvd( LAPACK_ROW_MAJOR, 'N', 'N',
			 m, n, a, lda, s,
			 NULL, ldu, NULL, ldvt, superb );
  if( info > 0 ) {
    printf( "- SVD convergence: failed \n" );
    exit( 1 );
  }
  *(activity+z) = 1;
  //printf("Point:%d --> smallest singular value is [%.2f,%.2f,%.2f ... %.2f,%.2f, %.2f] \n",
  //	 z,s[0],s[1],s[2],s[m-3],s[m-2],s[m-1]);
  if(s[m-1] > e[0]) { // first curve
    stepx = ceil((s[m-1] - e[0])/dm->stepx);
    stepy = ceil((s[m-1] - e[0])/dm->stepy);
    //printf("Disk c:%d s[m-1] - e[0]=%f and dm->stepx=%f ( %f , %f )\n", iy,s[m-1] - e[0], dm->stepx, stepx, stepy);
    start_i = i - stepy;
    if( start_i < 0 )
      start_i = 0;
    start_j = j - stepx;
    if( start_j < 0 )
      start_j = 0;
    end_i = i + stepy;
    if( end_i > gsize - 1 )
      end_i = gsize - 1;
    end_j = j + stepx;
    if( end_j > gsize - 1)
      end_j = gsize - 1;

    printf("array of z%d:(%d,%d) == start:(%d,%d) // end:(%d,%d) stepx: %f, stepy:%f \n",
	   z,i,j,start_i,start_j, end_i,end_j,stepx,stepy);
    for (int i = start_i ; i < end_i + 1; i++){
      for (int j = start_j ; j < end_j + 1; j++){
	//if (j < end_j) {
	  printf("point %d(%d,%d) in disk!\n",i*gsize + j, i, j);
	  *(activity+( i*gsize + j) ) = 1;
	  //}
	  //else continue;
      }
    }
  }
  
  free(superb);
  return s[m-1];

}
