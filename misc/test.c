#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
//#include <sys/time.h>
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


struct mog_status {
  int skip;
  double ssv;
};

struct disk {
  int si;
  int ei;
  int sj;
  int ej;
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

  /*     /\* Time variables *\/ */
  //struct timeval interval;

  struct diagnostics * diag = NULL;
  uint32_t * activity =  malloc((gsize*gsize)*sizeof( uint32_t));
  for (int i = 0; i < gsize*gsize ; i++)
    activity[i] = 0;

  a = malloc((lda*m)*sizeof(lapack_complex_double));
  create_grcar(a,m,lda);
  // todo: checkif input matrix is correct!!!!
  //print_matrix("grcar",m,n,a,lda);	
  scalar_matrix_mult (m, n, a, lda, -1);
  //print_matrix("- grcar",m,n,a,lda);	

  struct domain dm; 
  initDomain(&dm, gsize, XMIN, XMAX, YMIN, YMAX);
  diag = pseudospectra(m, n, gsize, NEPSILON, &dm, a, activity);
  assert(activity);
  assert(diag);
  printDiagnostics(m, n, &dm, gsize, diag, activity);
		     
  
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
  //lapack_int ldu = m;
  //lapack_int ldvt = n;
  //lapack_int info = 0;
  int svdPoints = 0;
  


  double ssv;
  struct mog_status * mg = NULL; //malloc(sizeof(struct mog_status));
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


  printf("Pseudospectra with mog\n");  
  for (int iy = 0; iy < gsize*gsize; iy++){
    // reload original matrix, since zgesvd writes on it (acp)
    create_grcar(acp, m,lda);
    scalar_matrix_mult (m, n, acp, lda, -1);
    //print_matrix("- grcar",m,n,a,lda);	

    for (int i = 0; i < lda*m ; i=i+(n+1))
      acp[i]=acp[i]+(dm->x_min+(iy/gsize * dm->stepx)+(dm->y_min + (iy % gsize * dm->stepy))*I); // -a + z 

    //ssv = grid(m,n,acp,gsize,activity,iy);
    mg = mog(m,n,nbepsilon,e,dm,acp,gsize,activity,iy);
    assert(mg);
    if(!(mg->skip))
      svdPoints++;
    diag->ssv[iy] = mg->ssv;
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
  if(mg)
    free(mg);
  free(acp);
  return diag;
  
}


double grid(lapack_int m, lapack_int n,
	    lapack_complex_double *a,
	    lapack_int gsize,
	    uint32_t * activity,
	    int iy) {

  lapack_int lda = n;
  lapack_int ldu = m;
  lapack_int ldvt = n;
  lapack_int info = 0;

  double *s;
  double *superb;
  s = malloc(m*sizeof(double));
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
  printf("Point:%d ssv:[%.2f,%.2f,%.2f ... %.2f,%.2f, %.2f]\n",
	 iy, s[0], s[1], s[2], s[m-3], s[m-2], s[m-1]);

  
  free(superb);
  return s[m-1];

}



struct mog_status * mog(lapack_int m, lapack_int n,
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
  assert(a);
  assert(activity);

  struct disk * dk = malloc(sizeof(struct disk));
  struct mog_status * mg = malloc(sizeof(struct mog_status));
  mg->skip = 0;
  mg->ssv = 0.0;
  
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
    mg->skip = 1;
    return mg;
  }
  info = LAPACKE_zgesvd( LAPACK_ROW_MAJOR, 'N', 'N',
			 m, n, a, lda, s,
			 NULL, ldu, NULL, ldvt, superb );
  if( info > 0 ) {
    printf( "- SVD convergence: failed \n" );
    exit( 1 );
  }
  *(activity+z) = 1;
  mg->ssv = s[m-1]; 
  if(s[m-1] > e[0]) {// first curve
    printf("Excluding disk (%d,%.2f)\n",z, s[m-1]-e[0]);
    locateDisk(s[m-1]-e[0], z, gsize, dm, dk);
    excludeDisk(gsize, dk, activity,1);
    //locateExcludeDisk(s[m-1]-e[0], z, gsize, dm, activity);
    printf("activity (inside):");
    for (int i = 0; i < gsize * gsize ; i++) {
      if(!(i % gsize))
        printf("\n");
      printf("%d ",*(activity+i));
    }
    printf("\n");

  }
  
  free(dk);
  free(superb);
  return mg;

}



struct mog_status * mmog(lapack_int m, lapack_int n,
			lapack_int nbepsilon, double * e,
			struct domain * dm,
			lapack_complex_double *a,
			lapack_int gsize,
			 // uint32_t ** activity,
			uint32_t * activity,
			int z) {
  
  lapack_int lda = n;
  lapack_int ldu = m;
  lapack_int ldvt = n;
  lapack_int info = 0;

  double *s;
  double *superb;
  s = malloc(m*sizeof(double));
  superb = malloc(min(m,n)*sizeof(double));
  assert(a);
  assert(activity);

  struct mog_status * mg = malloc(sizeof(struct mog_status));
  struct disk * dk = malloc(sizeof(struct disk));
  mg->skip = 0;
  mg->ssv = 0.0;
  
  for(int i = 0; i < nbepsilon; i++) {
  }
  if(*(activity+z) == 1) {
    printf("Point:%d is skipped!\n",z);
    mg->skip = 1;
    return mg;
  }
  info = LAPACKE_zgesvd( LAPACK_ROW_MAJOR, 'N', 'N',
			 m, n, a, lda, s,
			 NULL, ldu, NULL, ldvt, superb );
  if( info > 0 ) {
    printf( "- SVD convergence: failed \n" );
    exit( 1 );
  }
  *(activity+z) = 1;
  mg->ssv = s[m-1]; 
  if(s[m-1] > e[0]) {// first curve
    printf("Excluding disk (%d,%.2f)\n",z, s[m-1]-e[0]);
    locateDisk(s[m-1]-e[0], z, gsize, dm, dk);
    excludeDisk(gsize, dk, activity,1);
  }
  
  free(dk);
  free(superb);
  return mg;

}



void locateDisk(double radius, int center, int gsize, struct domain * dm, struct disk *dk) {

  double stepx = 0.0;
  double stepy = 0.0;
  int start_i = 0;
  int start_j = 0;
  int end_i = gsize;
  int end_j = gsize;
  // todo: i indexes lines of matrix so [i-stepy, i, i+stepy] changes lines
  // todo: j indexes columns of matrix so [j-stepx, j, j+stepx] changes columns 
  
  assert(dk && dm);  
  stepx = ceil(radius/dm->stepx);
  stepy = ceil(radius/dm->stepy);
  printf("gsize=%d\n",gsize);
  printf("Disk c:%d radius=%.2f and dm->stepx=%f ( %f , %f )\n", center,radius, dm->stepx, stepx, stepy);
  int i = center/gsize;
  int j = center%gsize;

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

  
    dk->si = start_i;
    dk->ei = end_i;
    dk->sj = start_j;
    dk->ej = end_j;
    printf("Disk start=(%d,%d), end=(%d,%d)\n",dk->si,dk->sj, dk->ei,dk->ej);


}


void excludeDisk(int gsize, struct disk * dk, uint32_t * activity, int value) {
  assert(dk && activity);
  for (int i = dk->si ; i < dk->ei + 1; i++){
    for (int j = dk->sj ; j < dk->ej + 1; j++){
      //if (j < end_j) {
      printf("point %d(%d,%d) in disk!\n",i*gsize + j, i, j);
      *(activity+( i*gsize + j) ) = value; // 1
      //}
      //else continue;
    }
  }
}



void locateExcludeDisk(double radius, int center, int gsize, struct domain * dm, uint32_t * activity) {

  double stepx = 0.0;
  double stepy = 0.0;
  int start_i = 0;
  int start_j = 0;
  int end_i = gsize;
  int end_j = gsize;
  // todo: i indexes lines of matrix so [i-stepy, i, i+stepy] changes lines
  // todo: j indexes columns of matrix so [j-stepx, j, j+stepx] changes columns 
  

  stepx = ceil(radius/dm->stepx);
  stepy = ceil(radius/dm->stepy);
  printf("Disk c:%d radius=%f and dm->stepx=%f ( %f , %f )\n", center,radius, dm->stepx, stepx, stepy);
  int i = center/gsize;
  int j = center%gsize;

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

    printf("Disk start=(%d,%d), end=(%d,%d)\n",start_i,start_j, end_i,end_j);
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


void initDomain(struct domain *dm, int gsize, int xmin, int xmax, int ymin, int ymax) {
  
  dm->x_min = xmin;
  dm->x_max = xmax;
  dm->y_min = ymin;
  dm->y_max = ymax;
  dm->stepx = (dm->x_max - dm->x_min)/gsize;
  dm->stepy=(dm->y_max - dm->y_min)/gsize;

}



void printDiagnostics(lapack_int m, lapack_int n,
		      struct domain *dm,
		      int gsize,
		      struct diagnostics * diag,
		      uint32_t * activity
		      ) {

    /* Time variables */
  struct timeval interval;
  
  printf("Input Diagnostics:\n");
  printf("- Matrix: grcar(%d,%d)\n",m,n);
  printf("- Gridsize: %d\n",gsize);
  printf("- Domain size is: X=[%f-%f]  Y=[%f-%f] (stepx=%f, stepy=%f)\n",
	 dm->x_min, dm->x_max,
	 dm->y_min, dm->y_max,
	 dm->stepx, dm->stepy);

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
  printf(" - Gain percentage: %d \%\n",
	 ((gsize*gsize-diag->svdPoints)/gsize*gsize));

}
