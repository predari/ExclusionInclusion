#ifndef AUX_H
#define AUX_H

#include <lapacke.h>
// todo: check if it complies with macos
#include <stdint.h>
#include <sys/time.h>

extern struct diagnostics diag;
extern struct domain dm;
extern struct mog_status mg;
extern struct disk dk;

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

/* multiple a matrix with a scalar a */
void scalar_matrix_mult (
                   lapack_int m,
                   lapack_int n,
                   lapack_complex_double* a,
                   lapack_int lda,
		   int s
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



struct diagnostics * pseudospectra(lapack_int m, lapack_int n, lapack_int gsize,
			  lapack_int nbepsilon, struct domain * dm,
			  lapack_complex_double *a, uint32_t * activity );

void save_array(const char * fname,
		int m, int n,
		double *ssv
		);

// todo: grid doesn't use all arguments for far
double grid(lapack_int m, lapack_int n,
	    lapack_complex_double *a,
	    lapack_int gsize,
	    uint32_t * activity,
	    int point
	    );


struct mog_status * mog(lapack_int m, lapack_int n,
	   lapack_int nbepsilon, double * e,
	   struct domain * dm,
	   lapack_complex_double *a,
	   lapack_int gsize, uint32_t * activity,
	   int iy);

struct mog_status * mmog(lapack_int m, lapack_int n,
			lapack_int nbepsilon, double * e,
			struct domain * dm,
			lapack_complex_double *a,
			lapack_int gsize,
			uint32_t * activity,
			 int z);

void mergeActivity( int nbepsilon,
		    uint32_t *eactivity[],
		    uint32_t * activity,
		    int gsize);

void excludeDisk(int gsize,
		 struct disk * dk,
		 uint32_t * activity,
		 int value);

void locateDisk(double radius,
		int center,
		int gsize,
		struct domain * dm,
		struct disk *dk);

void locateExcludeDisk( double radius,
			int center,
			int gsize,
			struct domain * dm,
			uint32_t * activity
			);

void initDomain(struct domain *dm,
		int gsize,
		double xmin, double xmax,
		double ymin, double ymax);

void printDiagnostics(lapack_int m, lapack_int n,
		      struct domain *dm,
		      int gsize,
		      struct diagnostics * diag,
		      uint32_t * activity
		      );


#endif
