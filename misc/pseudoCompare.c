#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
//#include <sys/time.h>
//#include <lapacke.h>
#include <assert.h>
#include <stdio.h>
#include <errno.h>

#include "aux.h"

// ps file1 file2 gsize #epsilon

#define min(a,b) ((a)>(b)?(b):(a))


int main(int argc, char * argv[])  {

  if(argc != 5) {
      perror("Wrong number of arguments");
      exit(1);
  }

  const char * f1 = argv[1];
  const char * f2 = argv[2];
  int gsize = atoi(argv[3]);
  int nbepsilon = atoi(argv[4]);
  assert(gsize>0);
  assert(f1 && f2);
  double ssv1[gsize*gsize];
  double ssv2[gsize*gsize];
  read_array(f1,gsize,gsize,ssv1);
  read_array(f2,gsize,gsize,ssv2);
  int same = 0;
  double e[nbepsilon];
  for(int i = 0; i < nbepsilon; i++) {
    e[i] = pow(0.1,(i+1));
  }
  
  for (int i = 0; i < gsize*gsize ; i++) {
    same = 0;
    if(ssv1[i] == ssv2[i])
      continue;
    else {
      for(int j = 0; j < nbepsilon; j++) {
	if(ssv1[i] <= e[j] && ssv2[i] <= e[j]) {
	  same = 1;
	  break;
	}
      }
      if (same == 0) {
	printf("point in position %d is not the same in those files \n",i);
	
      }
    }
  }
  return 0;
  
}
