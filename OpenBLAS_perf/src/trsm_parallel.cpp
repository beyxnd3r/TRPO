#include "trsm.h"
#include <cblas.h>
#include <algorithm>

void trsm_double_parallel(
    Side side, Uplo uplo, Trans trans, Diag diag,
    int n, int m,
    double alpha,
    const double* A,
    double* B
) {
    const int BS = 12;

    if (!(side == Left && uplo == Lower && trans == NoTrans)) return;

    
    for (int i = 0; i < n * m; i++)
        B[i] *= alpha;

    for (int ii = 0; ii < n; ii += BS) {

        int ib = std::min(BS, n - ii);

        
        for (int i = 0; i < ib; i++) {
            int row = ii + i;

            for (int k = ii; k < row; k++) {
                double aik = A[row*n + k];

                for (int j = 0; j < m; j++) {
                    B[row*m + j] -= aik * B[k*m + j];
                }
            }

            if (diag == NonUnit) {
                double d = A[row*n + row];

                for (int j = 0; j < m; j++) {
                    B[row*m + j] /= d;
                }
            }
        }

        
        int next = ii + ib;
        if (next < n) {

            int rows = n - next;

            cblas_dgemm(
                CblasRowMajor,
                CblasNoTrans,
                CblasNoTrans,
                rows, m, ib,
                -1.0,
                &A[next*n + ii], n,
                &B[ii*m], m,
                1.0,
                &B[next*m], m
            );
        }
    }
}