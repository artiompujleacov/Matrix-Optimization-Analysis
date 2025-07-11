#include "utils.h"
#include <cblas.h>
#include <stdlib.h>
#include <string.h>

double *my_solver(int N, double *A, double *B, double *x)
{
    double *C = (double *)malloc(N * N * sizeof(double));
    double *D = (double *)malloc(N * N * sizeof(double));
    double *tmp = (double *)malloc(N * sizeof(double));
    double *y = (double *)malloc(N * sizeof(double));
    double *result = (double *)malloc(N * sizeof(double));

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                N, N, N, 1.0, B, N, A, N, 0.0, C, N);

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                N, N, N, 1.0, C, N, A, N, 0.0, D, N);

    memcpy(tmp, x, N * sizeof(double));

    for (int i = 0; i < N; ++i)
    {
        cblas_dgemv(CblasRowMajor, CblasTrans, N, N, 1.0, C, N, tmp, 1, 0.0, y, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, C, N, y, 1, 0.0, tmp, 1);
    }

    cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, D, N, tmp, 1, 0.0, result, 1);

    free(C);
    free(D);
    free(tmp);
    free(y);

    return result;
}
