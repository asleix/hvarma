#include <complex.h>
#include <stdio.h>
#include <string.h>

/** Compute gradient as a the matrix from a linear map from the variables.
      - zcx: horizontal autocovariance.
      - cv: vertical autocovariance.
      - zcvx: horizontal-vertical crosscovariance.
      - mu, nu: prediction error weights.
      - p: arma filter order.
      - maxtau: maximum covariance lag.

      output:
      - mat: coefficient matrix, size (3*p+2, 3*p+2)
      - indep: independent term, size (3*p+2)
 **/

/** Compute autocovariance and crosscovariance of v and x1 + j*x2 **/
void covariances(const double *x1, const double *x2,  const double *v,
    double complex *zcx, double *cv, double complex *zcxv,
    int size, double maxtau)
{
    // Create complex horizontal signal
    double complex zx[size];
    for(int i = 0; i < size; ++i) zx[i] = x1[i] + I*x2[i];

    // Crosscovariances
    int shift = maxtau;
    for(int itau = 0; itau <= maxtau; ++itau) {
        for(int it = 0; it < size-itau; ++it) {
            zcxv[shift+itau] += zx[it+itau]*v[it];
        }
        zcxv[shift+itau] /= (double) size; // ML estimator
    }
    for(int itau = -maxtau; itau <= -1; ++itau) {
        for(int it = 0; it < size+itau; ++it) {
            zcxv[itau+shift] += zx[it]*v[it-itau];
        }
        zcxv[itau+shift] /= (double) size; // ML estimator
    }

    // Autocovariances
    for(int itau = 0; itau <= maxtau; ++itau) {  // forward lag
        for(int it = 0; it < size-itau; ++it) {
            cv[shift+itau] += v[it+itau]*v[it];
            zcx[shift+itau] += zx[it+itau]*conj(zx[it]);
        }
        cv[itau+shift] /= (double) size;
        zcx[itau+shift] /= (double) size;
    }
    for(int itau = -maxtau; itau <= -1; ++itau) { // backward lag
        cv[shift+itau] = cv[shift-itau];
        zcx[shift+itau] = conj(zcx[shift-itau]);
    }
}


/** Perform matrix coefficients loop.
    Here, mat is supposed to be dim (3p+3, 3p+3)  **/
void coefficients(const double complex *zcx, const double *cv,
    const double complex *zcxv, size_t size, double mat[size][size],
    double mu, double nu, int p, int maxtau)
{
    for(int i = 0; i <= p; ++i) {
        int ia = i;
        int ib1 = p+i+1;
        int ib2 = 2*p+i+2;
        for(int j = 0; j <= p; ++j) {
            int ja = j;
            int jb1 = p+j+1;
            int jb2 = 2*p+j+2;

            for(int tau = p; tau <= 2*maxtau; ++tau) {
                int iti=tau-i;
                int itj=tau-j;

                // Coefs. from vertical components
                if (ia != 0) {
                    if (ja != 0)
                        mat[ia][ja] += 2*(mu*creal(zcxv[iti]*conj(zcxv[itj]))
                                         + nu*creal(zcx[iti]*conj(zcx[itj])));
                    mat[ia][jb1] += -2*(mu*creal(zcxv[iti])*cv[itj]
                                        + nu*creal(zcx[iti]*conj(zcxv[itj])));
                    mat[ia][jb2] += -2*(mu*cimag(zcxv[iti])*cv[itj]
                                        + nu*cimag(zcx[iti]*conj(zcxv[itj])));
                }

                // Real part coefs. from horizontal comp.
                if(ja != 0)
                    mat[ib1][ja] += -2*(mu*cv[iti]*creal(zcxv[itj])
                                        + nu*creal(zcxv[iti]*conj(zcx[itj])));
                mat[ib1][jb1] += 2*(mu*cv[iti]*cv[itj]
                                    + nu*creal(zcxv[iti]*conj(zcxv[itj])));
                mat[ib1][jb2] += -2*nu*cimag(zcxv[iti]*conj(zcxv[itj]));

                // Imag part coefs. from horizontal comp.
                if(ja != 0)
                    mat[ib2][ja] += -2*(mu*cv[iti]*cimag(zcxv[itj])
                                        +nu*cimag(zcxv[iti]*conj(zcx[itj])));
                mat[ib2][jb1] += -2*nu*cimag(zcxv[iti]*conj(zcxv[itj]));
                mat[ib2][jb2] += 2*(mu*cv[iti]*cv[itj]
                                    + nu*creal(zcxv[iti]*conj(zcxv[itj])));
            }
        }
    }
}


/** Perform independent term  loop.
    Here, indep is supposed to be length (3p+3)  **/
void independent_term(const double complex *zcx, const double *cv,
    const double complex *zcxv, double *indep, double mu, double nu,
    int p, int maxtau)
{
    for(int i = 0; i <= p; ++i) {
        int ia=i;
        int ib1=i+p+1;
        int ib2=i+2*p+2;
        for(int tau = p; tau <= 2*maxtau; ++tau) {
            int iti=tau-i;
            int itj=tau;
            if(ia != 0)
                indep[ia] += -2*(mu*creal(zcxv[iti]*conj(zcxv[itj]))
                                + nu*creal(zcx[iti]*conj(zcx[itj])));
            indep[ib1] += 2*(mu*cv[iti]*creal(zcxv[itj])
                            + nu*creal(zcxv[iti]*conj(zcx[itj])));
            indep[ib2] += 2*(mu*cv[iti]*cimag(zcxv[itj])
                                +nu*cimag(zcxv[iti]*conj(zcx[itj])));
        }
    }
}


/** Function to compute gradient matrix. Calls previous
    functions to perform loops. **/
void gradient_matrix(const double complex *zcx, const double *cv,
    const double complex *zcxv, size_t size, double mat[size][size],
    double *indep, double mu, double nu, int p, int maxtau)
{
    // Fill coefficient matrix
    double cmat[size+1][size+1];
    memset( cmat, 0, (size+1)*(size+1)*sizeof(double) );
    coefficients(zcx, cv, zcxv, size+1, cmat, mu, nu, p, maxtau);
    // Fill independent term
    double cindep[size+1];
    memset( cindep, 0, (size+1)*sizeof(double) );
    independent_term(zcx, cv, zcxv, cindep, mu, nu, p, maxtau);

    // Fix index shift due to setting the first coef as constant
    for(int i = 0; i < (int)size; ++i) {
        indep[i] = cindep[i+1];
        for(int j = 0; j < (int)size; ++j) {
            mat[i][j] = cmat[i+1][j+1];
        }
    }
}

/** Main function to compute optimality conditions equations.
    Conists of finding variable values that set the gradient to 0. **/
void compute_equations(const double *x1, const double *x2, const double *v,
    double mu, double nu, size_t size, double mat[size][size], double *indep,
    size_t wsize, int p, int maxtau)
{

    // Compute covariances
    int cov_size = 2*maxtau+1;
    double cv[cov_size];
    memset( cv, 0, (cov_size)*sizeof(double) );
    double complex zcx[cov_size];
    memset( zcx, 0, (cov_size)*sizeof(double complex) );
    double complex zcxv[cov_size];
    memset( zcxv, 0, (cov_size)*sizeof(double complex) );

    covariances(x1, x2, v, zcx, cv, zcxv, wsize, maxtau);

    // Prediction weights
    if(mu == 0.) mu = 1./cv[maxtau];
    if(nu == 0.) nu = 1./creal(zcx[maxtau]);

    // Compute coef. matrix and indep. term
    gradient_matrix(zcx, cv, zcxv, size, mat, indep, mu, nu, p, maxtau);

}
