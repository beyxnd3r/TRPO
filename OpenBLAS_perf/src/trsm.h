#pragma once

enum Side { Left, Right };
enum Uplo { Upper, Lower };
enum Trans { NoTrans, Transpose };
enum Diag { NonUnit, Unit };

void trsm_double_parallel(
    Side side, Uplo uplo, Trans trans, Diag diag,
    int n, int m,
    double alpha,
    const double* A,
    double* B
);