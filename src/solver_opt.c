#include "utils.h"
#include <stdlib.h>

double *my_solver(int N, double *A, double *B, double *x)
{
    register int i, j, k;
    register double *At = (double *)malloc(N * N * sizeof(double));
    register double *C = (double *)calloc(N * N, sizeof(double));
    register double *Ct = (double *)malloc(N * N * sizeof(double));
    register double *D = (double *)calloc(N * N, sizeof(double));
    register double *v1 = (double *)malloc(N * sizeof(double));
    register double *v2 = (double *)malloc(N * sizeof(double));
    register double *res = (double *)malloc(N * sizeof(double));

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            At[j * N + i] = A[i * N + j];

    for (i = 0; i < N; i++)
    {
        register double *rowB = &B[i * N];
        register double *rowC = &C[i * N];
        for (k = 0; k < N; k++)
        {
            register double b = rowB[k];
            register double *rowAt = &At[k * N];
            for (j = 0; j <= N - 4; j += 4)
            {
                rowC[j] += b * rowAt[j];
                rowC[j + 1] += b * rowAt[j + 1];
                rowC[j + 2] += b * rowAt[j + 2];
                rowC[j + 3] += b * rowAt[j + 3];
            }
            for (; j < N; j++)
            {
                rowC[j] += b * rowAt[j];
            }
        }
    }

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            Ct[j * N + i] = C[i * N + j];

    for (i = 0; i < N; i++)
    {
        register double *rowCt = &Ct[i * N];
        register double *rowD = &D[i * N];
        for (k = 0; k < N; k++)
        {
            register double c = rowCt[k];
            register double *rowA = &A[k * N];
            for (j = 0; j <= N - 4; j += 4)
            {
                rowD[j] += c * rowA[j];
                rowD[j + 1] += c * rowA[j + 1];
                rowD[j + 2] += c * rowA[j + 2];
                rowD[j + 3] += c * rowA[j + 3];
            }
            for (; j < N; j++)
            {
                rowD[j] += c * rowA[j];
            }
        }
    }

    for (i = 0; i < N; i++)
        v1[i] = x[i];

    register double *cur = v1, *next = v2;
    for (int iter = 0; iter < N; iter++)
    {
        for (i = 0; i < N; i++)
        {
            register double sum = 0.0;
            register double *rowCt = &Ct[i * N];
            for (j = 0; j <= N - 4; j += 4)
            {
                sum += rowCt[j] * cur[j];
                sum += rowCt[j + 1] * cur[j + 1];
                sum += rowCt[j + 2] * cur[j + 2];
                sum += rowCt[j + 3] * cur[j + 3];
            }
            for (; j < N; j++)
            {
                sum += rowCt[j] * cur[j];
            }
            next[i] = sum;
        }

        register double *tmp = cur;
        cur = next;
        next = tmp;

        for (i = 0; i < N; i++)
        {
            register double sum = 0.0;
            register double *rowC = &C[i * N];
            for (j = 0; j <= N - 4; j += 4)
            {
                sum += rowC[j] * cur[j];
                sum += rowC[j + 1] * cur[j + 1];
                sum += rowC[j + 2] * cur[j + 2];
                sum += rowC[j + 3] * cur[j + 3];
            }
            for (; j < N; j++)
            {
                sum += rowC[j] * cur[j];
            }
            next[i] = sum;
        }

        tmp = cur;
        cur = next;
        next = tmp;
    }

    for (i = 0; i < N; i++)
    {
        register double sum = 0.0;
        register double *rowD = &D[i * N];
        int j = 0;
        for (; j <= N - 4; j += 4)
        {
            sum += rowD[j] * cur[j];
            sum += rowD[j + 1] * cur[j + 1];
            sum += rowD[j + 2] * cur[j + 2];
            sum += rowD[j + 3] * cur[j + 3];
        }
        for (; j < N; j++)
        {
            sum += rowD[j] * cur[j];
        }
        res[i] = sum;
    }

    free(At);
    free(C);
    free(Ct);
    free(D);
    free(v1);
    free(v2);

    return res;
}