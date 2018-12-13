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
#define GRID 100 // grid size
#define LDA N
#define LDU M
#define LDVT N
#define MAXNEPSILON 32 // unit32_t used for the activity array
#define NEPSILON 3


// following is the input values of the domain.
// in a first attempt, we copied them from GK paper.
// later on we should use one of the following methods:
// Gresh. Disk, FoV, matrix norms, QR
// see notes on which one is more appropriate for us.

// fov -0.6368    2.9736   -3.1160    3.1160
//      -0.91     3.27     -3.41      3.41
#define YMAX    3.1160  
#define YMIN    -3.1160
#define XMAX     2.9736
#define XMIN    -0.6368                  

#define min(a,b) ((a)>(b)?(b):(a))

struct diagnostics {
  struct timeval earlier;
  struct timeval later;
  int svdPoints;
  double *ssv;
  double *pss;

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



int main(int argc, char * argv[])  {

  if(argc != 6) {
      perror("Wrong number of arguments. Usage: main method m n gsize nbepsilon");
      exit(1);
  }

  const char * meth = argv[1];
  int m_in =  atoi(argv[2]);
  int n_in =  atoi(argv[3]);
  int gsize = atoi(argv[4]);
  int nbepsilon = atoi(argv[5]);
  assert(gsize>0 && m_in>0 && n_in>0 && nbepsilon > 0);
  assert(meth);

  // TODO: 
  /* USER INPUT PARAMETERS */
  /*   - NBEPSILON */
  /*   - m, n */
  /*   - GRID size */
  /*   - method choice */
  /* LATER */
  /*   - a */
  /*   - XMIN, XMAX, YMIN, YMAX */  // check that min<max
  lapack_complex_double *a = NULL;
  lapack_int m = m_in, n = n_in;
  lapack_int lda = LDA, ldu = LDU, ldvt = LDVT;
  int svdPoints = 0;
  //int gsize = GRID;

  struct diagnostics * diag = NULL;
  uint32_t * activity =  malloc((gsize*gsize)*sizeof( uint32_t));
  for (int i = 0; i < gsize*gsize ; i++)
    activity[i] = 0;

  a = malloc((lda*m)*sizeof(lapack_complex_double));
  create_grcar(a,m,lda);
  //print_matrix("grcar",m,n,a,lda);	

  struct domain dm;
  initDomain(&dm, gsize, XMIN, XMAX, YMIN, YMAX);

  
  diag = pseudospectra(m, n, gsize, nbepsilon, meth, &dm, a, activity);
  assert(activity);
  assert(diag);
  printDiagnostics(m, n, &dm, gsize, diag, activity);
  printf("saving ssv in file\n");
  save_array("ssv.txt",gsize,gsize,diag->ssv);		     
  
  free(a);
  free (activity);
  if (diag)
    free(diag->ssv);
  free(diag);
  return 0;
  
}



struct diagnostics * pseudospectra(lapack_int m, lapack_int n, lapack_int gsize,
				   lapack_int nbepsilon, const char * meth, struct domain * dm,
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
  assert(dm);
  assert(meth);

  /* Diagnostic variables */
  struct diagnostics * diag = malloc(sizeof(struct diagnostics));
  diag->ssv = malloc((gsize*gsize)*sizeof(double));
  diag->pss = malloc((gsize*gsize)*sizeof(double));

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

  for (int iy = 0; iy < gsize*gsize; iy++){
    // reload original matrix, since zgesvd writes on it (acp)
    create_grcar(acp, m,lda);
    //print_matrix("grcar acp",m,n,acp,lda);
    for (int i = 0; i < lda*m ; i=i+(n+1))
      acp[i]=acp[i]-(dm->x_min+(iy/gsize * dm->stepx)+(dm->y_min + (iy % gsize * dm->stepy))*I);

    /* printf("gridpoint (%f %f) = ",dm->x_min+(iy/gsize * dm->stepx),(dm->y_min + (iy % gsize * dm->stepy))); */

    if(!strcmp(meth,"grid")) {
      ssv = grid(m,n,acp,gsize,activity,iy);
    //printf("%f\n",ssv);
    } else if(!strcmp(meth,"mog")) {
    mg = mog(m,n,nbepsilon,e,dm,acp,gsize,activity,iy);
    //printf("%f\n",mg->ssv);
    } else if(!strcmp(meth,"mmog")) {
    //mg = mmog(m,n,nbepsilon,e,dm,acp,gsize,activity,iy);
      printf("Sorry, mmog is not implemented yet!\n");
    }
    else {
      printf("Wrong method argument! Choose one of grid, mog or mmog.\n");
      exit(1);
    }
    if(mg != NULL) {
      assert(diag->ssv);
      diag->ssv[iy] = mg->ssv;
      if(!(mg->skip))
	svdPoints++;
      if (diag->ssv[iy] <= e[0])
	diag->pss[iy] = diag->ssv[iy];
      else
	diag->pss[iy] = 1;
    }
    else {
      assert(diag->ssv);
      diag->ssv[iy] = ssv;
      diag->pss[iy] = ssv; // not necessary for grid
    }
  }

  if(gettimeofday(&diag->later,NULL)) {
    perror("sixth gettimeofday()");
    exit(1);
  }
  // if mog or mmog have been used
  if (mg != NULL)
    diag->svdPoints = svdPoints;
  // if grid has been used
  else diag->svdPoints = gsize * gsize;

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
  

  if(  *(activity+z) == 1) {
    printf("Point:%d is skipped!\n",z);
    mg->skip = 1;
    mg->ssv = 0.11; // TODO: magic number. explain
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
  printf("Point %d is svded with value %.4f!\n",z, s[m-1]);
  
  if(s[m-1] > e[0]) {// first curve
    printf("Excluding disk (%d,%.4f)\n",z, s[m-1]-e[0]);
    //locateDisk(s[m-1]-e[0], z, gsize, dm, dk);
    //excludeDisk(gsize, dk, activity,1);
    locateExcludeDisk(s[m-1]-e[0], z, gsize, dm, activity);
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
  int evalue[nbepsilon];
  for(int i = 0; i < nbepsilon; i++) {
    evalue[i] = i + 1;
  }
  uint32_t ** eactivity;
  eactivity =  malloc((nbepsilon)*sizeof( uint32_t *));
  for(int j = 0; j < nbepsilon; j++)
    eactivity[j] = malloc((gsize*gsize)*sizeof( uint32_t));
  
  for(int j = 0; j < nbepsilon; j++) {
   for(int i = 0; i < gsize*gsize; i++) {
     eactivity[j][i] = 0;
   }
  } 

  /* for(int i = 0; i < gsize*gsize; i++) { */
  /*   for(int j = 0; j < nbepsilon; j++) { */
  /*     if(activity[i] > 0 && activity[i] <= evalue[j]) { */
  /* 	eactivity[j][i] = evalue[j]; */
  /* 	break; */
  /*     } */
  /*   } */
  /* } */

  struct mog_status * mg = malloc(sizeof(struct mog_status));
  struct disk * dk = malloc(sizeof(struct disk));
  mg->skip = 0;
  mg->ssv = 0.0;
  

  if(*(activity+z) == 1 || *(activity+z) == 2 || *(activity+z) == 3) {
    printf("Point:%d is skipped!\n",z);
    mg->skip = 1;
    if(*(activity+z) == 1)
      mg->ssv = 0.1;
    else if(*(activity+z) == 2)
      mg->ssv = 0.01;
    else if(*(activity+z) == 3)
      mg->ssv = 0.001;
    return mg;
  }
  info = LAPACKE_zgesvd( LAPACK_ROW_MAJOR, 'N', 'N',
			 m, n, a, lda, s,
			 NULL, ldu, NULL, ldvt, superb );
  if( info > 0 ) {
    printf( "- SVD convergence: failed \n" );
    exit( 1 );
  }
  for(int j = 0; j < nbepsilon; j++) {
    mg->ssv = s[m-1];
    if(s[m-1] > e[j]) {// first curve
      printf("excluding disk (%d,%.2f) (e[%d]))\n",z, s[m-1]-e[j],j);
      locateDisk(s[m-1]-e[j], z, gsize, dm, dk);
      excludeDisk(gsize, dk, eactivity[j], j + 1);
      /* printf("activity (before merge):"); */
      /* for (int i = 0; i < gsize * gsize ; i++) { */
      /* 	if(!(i % gsize)) */
      /* 	  printf("\n"); */
      /* 	printf("%d ", eactivity[j][i]); */
      /* } */
      /* printf("\n");  */
    }
  }
  mergeActivity(nbepsilon,eactivity,activity,gsize);
  
  printf("activity (after merge):");
  for (int i = 0; i < gsize * gsize ; i++) {
    if(!(i % gsize))
      printf("\n");
    printf("%d ",*(activity+i));
  }
  printf("\n");
  free(dk);
  free(superb);
  return mg;

}

void mergeActivity( int nbepsilon,uint32_t ** eactivity,uint32_t * activity, int gsize) {

  assert(activity && eactivity);
  for(int j = 0; j < gsize*gsize; j++) {
    if(eactivity[0][j] == 1) {
      activity[j] = 1;
      continue;
    }
    if(eactivity[0][j] == 0 && eactivity[1][j] == 2) {
      activity[j] = 2;
      continue;
    }
    if(eactivity[0][j] == 0 && eactivity[1][j] == 0 && eactivity[2][j] == 3) {
      activity[j] = 3;
      continue;
    }
  }
}


void locateDisk(double radius, int center, int gsize, struct domain * dm, struct disk *dk) {

  int stepx = 0;
  int stepy = 0;
  int start_i = 0;
  int start_j = 0;
  int end_i = gsize;
  int end_j = gsize;
  // todo: i indexes lines of matrix so [i-stepy, i, i+stepy] changes lines
  // todo: j indexes columns of matrix so [j-stepx, j, j+stepx] changes columns 
  
  assert(dk && dm);  
  stepx = floor(radius/dm->stepx);
  stepy = floor(radius/dm->stepy);

  printf("radius=%.4f dm->stepx=%.4f dm->stepy=%.4f\n",
	 radius,dm->stepx,dm->stepy );

  printf("stepx=%d stepy=%.d cstepx=%.4f cstepy=%.4f\n",
	 stepx, stepy,
	 (radius/dm->stepx),(radius/dm->stepy) );
    
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

      /* TODO: this is because we try to exclude discretized squares and not disks
       * so I have to make sure that the current point still belongs to the disk.
       * y = sqrt(x^2 + y^2)*/
      /* if(sqrt(pow(i*dm->stepx,2) + pow(j*dm->stepy,2)) <= radius) { */
      /* 	*(activity+( i*gsize + j) ) = value; // 1 exclude point */
      /* } */
      
      *(activity+( i*gsize + j) ) = value; // 1
           printf("point %d(%d,%d) in disk!\n",i*gsize + j, i, j);
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
  

  stepx = floor(radius/dm->stepx);
  stepy = floor(radius/dm->stepy);
  printf("Disk c:%d radius=%f and dm->stepx=%f ( %f , %f )\n", center,radius, dm->stepx, stepx, stepy);
  int zi = center/gsize;
  int zj = center%gsize;
  double realradius = 0.0;
    start_i = zi - stepy;
    if( start_i < 0 )
      start_i = 0;
    start_j = zj - stepx;
    if( start_j < 0 )
      start_j = 0;
    end_i = zi + stepy;
    if( end_i > gsize - 1 )
      end_i = gsize - 1;
    end_j = zj + stepx;
    if( end_j > gsize - 1)
      end_j = gsize - 1;
    
    printf("Disk start=(%d,%d), end=(%d,%d)\n",start_i,start_j, end_i,end_j);
    for (int i = start_i ; i < end_i + 1; i++){
      for (int j = start_j ; j < end_j + 1; j++){
	// TODO: does not work correctly yet
	realradius= sqrt(pow(abs(zi-i)*dm->stepy,2) + pow(abs(zj-j)*dm->stepx,2));
	printf("point (%d,%d) realradius=%f\n",i,j,realradius);
	if(realradius <= radius) {
	    printf("point %d(%d,%d) in disk!\n",i*gsize + j, i, j);
	    *(activity+( i*gsize + j) ) = 1; // value
	}
      }
    }
}


void initDomain(struct domain *dm, int gsize, double xmin, double xmax, double ymin, double ymax) {

  assert(dm);
  dm->x_min = xmin;
  dm->x_max = xmax;
  dm->y_min = ymin;
  dm->y_max = ymax;
  
  assert(dm->x_min < dm->x_max);
  assert(dm->y_min < dm->y_max);
 
  dm->stepx = (dm->x_max - dm->x_min)/(gsize-1);
  dm->stepy = (dm->y_max - dm->y_min)/(gsize-1);

}



void printDiagnostics(lapack_int m, lapack_int n,
		      struct domain *dm,
		      int gsize,
		      struct diagnostics * diag,
		      uint32_t * activity
		      ) {

    /* Time variables */
  struct timeval interval;
  int pps_count = 0;
  
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
    //    if(diag->ssv[i] <=0.1)
      printf("%.4f ",diag->ssv[i]);
      //    else
      //      printf("0.0000");
  }  
  printf("\n");  

  assert(diag->pss);
  printf("- Pss table:\n");
  for (int i = 0; i < gsize ; i++) {
    for (int j = 0; j <= gsize*(gsize-1) ; j=j+gsize) {
      printf("%.4f ",diag->pss[j+i]);
      if (diag->pss[j+i] <= 0.1)
	pps_count++;
    }
    printf("\n");
  }
  
  printf("- Activity table: ");
  for (int i = 0; i < gsize * gsize ; i++) {
    if(!(i % gsize))
      printf("\n");
    printf("%d ",activity[i]);
  }  
  printf("\n");

  printf(" - Number of visited gridPoints: %d\n",diag->svdPoints);
  printf(" - Number of points in the pseudospectra: %d\n",pps_count);
  printf(" - Gain percentage: %.2f\n", (double)( (gsize*gsize) - diag->svdPoints )/ (gsize*gsize) );
}
