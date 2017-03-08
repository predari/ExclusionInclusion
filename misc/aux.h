#ifndef AUX_H
#define AUX_H

#include <lapacke.h>
#include <lapacke.h>


void create_grcar( lapack_complex_double *a,
                   lapack_int m,
                   lapack_int lda );

/* printing the content of the matrix*/
void print_matrix (char* desc, // description as a title
                   lapack_int m,
                   lapack_int n,
                   lapack_complex_double* a,
                   lapack_int lda
                   );

void print_int_vector (char* desc,
                       lapack_int n,
                       lapack_int* a
                       );

void prtdat (int nx,
             int ny,
             double * u1,
             char *fnam
            );

long long timeval_diff (struct timeval *difference,
                        struct timeval *end_time,
                        struct timeval *start_time
                        );


#endif