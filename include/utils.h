#ifndef DGALERKIN_UTILS_H
#define DGALERKIN_UTILS_H

namespace lapack {
    // Compute the inverse of a square matrix A of size N*N.
    // Input can either be a array or a vector. Since c++11
    // vectors are also stored as a contiguous memory chunk.
    // Lapack assumes column major matrix, however as the
    // inverse of a transposed is the transposed of an inverse.
    // The matrix can either be column or row major.
    void inverse(double *A, int &N);

    // Computes the solution to a real system of linear equations
    //         A * X = B,
    // where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
    // On exit, the solution X is stored in B.
    void solve(double *A, double *B, int &N);

    // Normalize a vector/matrix with respect to Frobenius norm
    void normalize(double *A, int &N);

    // Dot product between 2 vectors
    double dot(double *A, double *b, int N);

    // Matrix/vector product:  y := alpha*A*x + y
    void linEq(double *A, double *X, double *Y, double &alpha, int &N);

    double minus(double *A, double *B, int N);

    double plus(double *A, double *B, int N);
}

namespace eigen {
    // Computes the solution to a real system of linear equations
    //         A * X = B,
    // where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
    // On exit, the solution X is stored in B.
    void solve(double *A, double *B, int &N);
}

#endif //DGALERKIN_UTILS_H
