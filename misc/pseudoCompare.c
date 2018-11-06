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


#define min(a,b) ((a)>(b)?(b):(a))


int main(int argc, int * argv[])  {

  if(argc != 4) {
      perror("Wrong number of arguments");
      exit(1);
  }

  const char * f1 = argv[1];
  const char * f2 = argv[2];
  int gsize = atoi(argv[3]);
  assrt(gsize>0);
  assert(f1 && f2);
  double ssv1[gsize*gsize];
  double ssv2[gsize*gsize];
  read_array(f1,gsize,gsize,ssv1);
  read_array(f2,gsize,gsize,ssv2);
  
  for (int i = 0; i < gsize*gsize ; i++) {
    
    //if(ssv1[i] == ssv2[i] || ssv1[i])
  }
    return 0;
  
}
