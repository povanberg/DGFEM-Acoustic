#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <iomanip>

extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // Generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

    // Computes the solution to a real system of linear equations
    void dgesv_(int* N, int *NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);

    // Returns the value of the norm
    double dlange_(char* NORM, int* M, int *N, double* A, int* LDA, double* WORK);

    // Dot product
    double ddot_(int* N, double* DX, int* INCX, double* DY, int* INCY);

    // Matrix/vector product:  y := alpha*A*x + beta*y,
    void dgemv_(char& TRANS, int& M, int& N, double& a, double* A,
            int& LDA, double* X, int& INCX, double& beta, double* Y, int& INCY);
}


// Lapack with direct calling to fortran interface.
// => Migrate to Lapacke (conflict with Gmsh, segmentation fault)
namespace lapack {
    // Compute the inverse of a square matrix A of size N*N.
    // Input can either be a array or a vector. Since c++11
    // vectors are also stored as a contiguous memory chunk.
    // Lapack assumes column major matrix, however as the
    // inverse of a transposed is the transposed of an inverse.
    // The matrix can either be column or row major.
    void inverse(double *A, int &N) {
        int IPIV[N];
        int LWORK = N*N;
        double WORK[LWORK];
        int INFO;
        dgetrf_(&N,&N,A,&N,IPIV,&INFO);
        dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
    }

    // Computes the solution to a real system of linear equations
    //         A * X = B,
    // where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
    // On exit, the solution X is stored in B.
    void solve(double *A, double *B, int &N){
        int INFO;
        int IPIV[N];
        int NRHS = 1;
        dgesv_(&N, &NRHS, A, &N, IPIV, B, &N, &INFO);
    }

    // Normalize a vector/matrix with respect to Frobenius norm
    void normalize(double *A, int &N) {
        char NORM = 'F';
        int LDA = 1;
        double WORK[N];
        double norm;
        norm = dlange_(&NORM, &LDA, &N, A, &LDA, WORK);
        for(int i=0; i<N; i++){A[i] /= norm;};
    }

    // Matrix/vector product:  y := alpha*A*x + y
    void linEq(double *A, double *X, double *Y, double &alpha, int &N){
        char TRANS = 'T';
        int INC = 1;
        double beta=1;
        dgemv_(TRANS, N, N, alpha, A, N, X, INC, beta, Y, INC);
    }

    // Dot product between 2 vectors
    double dot(double *A, double *B, int N) {
        int INCX = 1;
        return ddot_(&N, A, &INCX, B, &INCX);
    }

    double minus(double *A, double *B, int N) {
        for(int i=0; i<N; ++i)
            A[i] -= B[i];
        //std::transform(A, A + N, B, A, std::minus<double>());
    }

    double plus(double *A, double *B, int N) {
        for(int i=0; i<N; ++i)
            A[i] += B[i];
        //std::transform(A, A + N, B, A, std::plus<double>());
    }
}

namespace eigen {

    // Computes the solution to a real system of linear equations
    //         A * X = B,
    // where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
    // On exit, the solution X is stored in B.
    void solve(double *A, double *B, int &N){
        Eigen::Map<Eigen::MatrixXd> A_eigen(A, N, N);
        Eigen::Map<Eigen::VectorXd> B_eigen(B, N);
        B_eigen = A_eigen.lu().solve(B_eigen);
    }

}

namespace display {
    template<typename Container>
    void print(const Container& cont, int row = 1, bool colMajor=false) {

        if(colMajor){
            for(int rowIt=0; rowIt<row; ++rowIt){
                int colIt = 0;
                for (auto const& x : cont) {
                    if(colIt%row == rowIt) {
                        std::cout << std::setprecision(4) << std::left << std::setw(10) << x << " ";
                    }
                    colIt++;
                }
                std::cout << std::endl;
            }
        }
        else {
            int colIt = 0;
            for (auto const& x : cont) {
                std::cout << std::setprecision(4) << std::left << std::setw(10) << x << " ";
                colIt++;
                if(colIt%row == 0)
                    std::cout << std::endl;
            }
        }

    }
}
