#include <gsl/gsl_multifit.h>

extern "C" {

void c_multifit_linear ( const int n, const int p, const double xm[],
	         	 const double yv[], double cv[], double covm[],
                         double* chisq, int* ierr ) {

  gsl_matrix *x, *cov;
  gsl_vector *y, *c;
  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (n, p);

  x = gsl_matrix_alloc(n, p);
  y = gsl_vector_alloc(n);
  c = gsl_vector_alloc(p);
  cov = gsl_matrix_alloc(p,p);

  for (int i=0; i<n; i++) {
    gsl_vector_set(y, i, yv[i]);
    for (int j=0; j<p; j++) {
      gsl_matrix_set(x, i, j, xm[i*p+j]);
    }
  }
  *ierr = gsl_multifit_linear ( x, y, c, cov, chisq, work);

  for (int i=0; i<p; i++) {
    cv[i] = gsl_vector_get(c,i);
    for (int j = 0; j<p; j++) {
      covm[i*p+j] = gsl_matrix_get(cov,i,j);
    }
  }
  gsl_matrix_free(x);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(work);
}

}
