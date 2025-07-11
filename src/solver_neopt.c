/*
 * Tema 2 ASC
 * 2025 Spring
 */
#include "utils.h"
#include <stdlib.h>
#include <string.h>

/*
 * Implementare neoptimizată (fără BLAS, fără optimizări)
 */
double *my_solver(int N, double *A, double *B, double *x)
{
    double *At = (double *)malloc(N * N * sizeof(double));
    double *C = (double *)malloc(N * N * sizeof(double));
    double *Ct = (double *)malloc(N * N * sizeof(double));
    double *D = (double *)malloc(N * N * sizeof(double));
    double *y = (double *)malloc(N * sizeof(double));
    double *tmp = (double *)malloc(N * sizeof(double));
    double *result = (double *)malloc(N * sizeof(double));

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            At[j * N + i] = A[i * N + j];

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
        {
            C[i * N + j] = 0.0;
            for (int k = 0; k < N; ++k)
                C[i * N + j] += B[i * N + k] * At[k * N + j];
        }

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            Ct[j * N + i] = C[i * N + j];

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
        {
            D[i * N + j] = 0.0;
            for (int k = 0; k < N; ++k)
                D[i * N + j] += Ct[i * N + k] * A[k * N + j];
        }

    memcpy(tmp, x, N * sizeof(double));
    for (int iter = 0; iter < N; ++iter)
    {
        for (int i = 0; i < N; ++i)
        {
            y[i] = 0.0;
            for (int j = 0; j < N; ++j)
                y[i] += Ct[i * N + j] * tmp[j];
        }
        for (int i = 0; i < N; ++i)
        {
            tmp[i] = 0.0;
            for (int j = 0; j < N; ++j)
                tmp[i] += C[i * N + j] * y[j];
        }
    }

    for (int i = 0; i < N; ++i)
    {
        result[i] = 0.0;
        for (int j = 0; j < N; ++j)
            result[i] += D[i * N + j] * tmp[j];
    }

    free(At);
    free(C);
    free(Ct);
    free(D);
    free(y);
    free(tmp);
    return result;
}