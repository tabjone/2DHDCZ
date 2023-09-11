#include "solve_elliptic_equation.h"

void gauss_seidel(double **A, double *b, double *x, int N, int maxIterations, double tolerance) {
    double oldX[N];
    for (int i = 0; i < N; i++) {
        x[i] = 0.0; // initial guess, change this to last iteration of previous time step
    }

    for (int iter = 0; iter < maxIterations; iter++) {
        for (int i = 0; i < N; i++) {
            oldX[i] = x[i];
        }

        for (int i = 0; i < N; i++) {
            double sum1 = 0.0, sum2 = 0.0;
            for (int j = 0; j < i; j++) {
                sum1 += A[i][j] * x[j];
            }
            for (int j = i + 1; j < N; j++) {
                sum2 += A[i][j] * oldX[j];
            }

            x[i] = (b[i] - sum1 - sum2) / A[i][i];
        }

        double norm = 0.0;
        for (int i = 0; i < N; i++) {
            norm += (x[i] - oldX[i]) * (x[i] - oldX[i]);
        }
        if (sqrt(norm) < tolerance) {
            break;
        }
    }
}