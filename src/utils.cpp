#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

#include "utils.h"

extern "C"
{
    // LU decomoposition of a general matrix
    void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);

    // Generate inverse of a matrix given its LU decomposition
    void dgetri_(int *N, double *A, int *lda, int *IPIV, double *WORK, int *lwork, int *INFO);

    // Computes the solution to a real system of linear equations
    void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);

    // Returns the value of the norm
    double dlange_(char *NORM, int *M, int *N, double *A, int *LDA, double *WORK);

    // Dot product
    double ddot_(int *N, double *DX, int *INCX, double *DY, int *INCY);

    // Matrix/vector product:  y := alpha*A*x + beta*y,
    void dgemv_(char &TRANS, int &M, int &N, double &a, double *A,
                int &LDA, double *X, int &INCX, double &beta, double *Y, int &INCY);
}

//! added by Sofiane KHELLADI in 11/03/2022 /////////////////////

namespace screen_display
{
    void write_string(std::string text, std::string color)
    {
#ifdef _WIN64
        std::cout << text << std::endl
                  << std::flush;
//#endif
#else
        //#ifdef _linux_
        std::cout << color << text << RESET << std::endl
                  << std::flush;
#endif // _linux_
    }

    void write_value(std::string name, double value, std::string unit, std::string color)
    {
#ifdef _WIN64
        if (fabs(value) >= 1.0e-3)
        {
            std::cout << std::setw(50) << std::left << name << std::right << std::fixed << std::setprecision(6) << std::setw(30) << value << std::right << " [" << unit << "]" << std::endl
                      << std::flush;
        }
        else
        {
            std::cout << std::setw(50) << std::left << name << std::right << std::scientific << std::setw(30) << value << std::right << " [" << unit << "]" << std::endl
                      << std::flush;
        }
//#endif
#else
        //#ifdef _linux_
        if (fabs(value) >= 1.0e-3)
        {
            std::cout << RESET << std::setw(50) << std::left << name << std::right << color << std::fixed << std::setprecision(6) << std::setw(30) << value << RESET << std::right << " [" << unit << "]" << std::endl
                      << std::flush;
        }
        else
        {
            std::cout << RESET << std::setw(50) << std::left << name << std::right << color << std::scientific << std::setw(30) << value << RESET << std::right << " [" << unit << "]" << std::endl
                      << std::flush;
        }
#endif
    }

    void write_value_r(std::string name, double value, std::string unit, std::string color, size_t precision)
    {
#ifdef _WIN64
        if (fabs(value) >= 1.0e-3)
        {
            std::cout << "\r" << std::setw(50) << std::left << name << std::right << std::fixed << std::setprecision(precision) << std::setw(30) << value << std::right << " [" << unit << "]";
        }
        else
        {
            std::cout << "\r" << std::setw(50) << std::left << name << std::right << std::scientific << std::setw(30) << value << std::right << " [" << unit << "]";
        }
//#endif
#else
        //#ifdef _linux_
        if (fabs(value) >= 1.0e-3)
        {
            std::cout << "\r" << RESET << std::setw(50) << std::left << name << std::right << color << std::fixed << std::setprecision(precision) << std::setw(30) << value << RESET << std::right << " [" << unit << "]";
        }
        else
        {
            std::cout << "\r" << RESET << std::setw(50) << std::left << name << std::right << color << std::scientific << std::setw(30) << value << RESET << std::right << " [" << unit << "]";
        }
#endif
    }
}
///////////////////////

// Lapack with direct calling to fortran interface.
// => Migrate to Lapacke (conflict with Gmsh, segmentation fault)
namespace lapack
{
    // Compute the inverse of a square matrix A of size N*N.
    // Input can either be a array or a vector. Since c++11
    // vectors are also stored as a contiguous memory chunk.
    // Lapack assumes column major matrix, however as the
    // inverse of a transposed is the transposed of an inverse.
    // The matrix can either be column or row major.
    void inverse(double *A, int &N)
    {
        int IPIV[N];
        int LWORK = N * N;
        double WORK[LWORK];
        int INFO;
        dgetrf_(&N, &N, A, &N, IPIV, &INFO);
        dgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO);
    }

    // Computes the solution to a real system of linear equations
    //         A * X = B,
    // where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
    // On exit, the solution X is stored in B.
    void solve(double *A, double *B, int &N)
    {
        int INFO;
        int IPIV[N];
        int NRHS = 1;
        dgesv_(&N, &NRHS, A, &N, IPIV, B, &N, &INFO);
    }

    // Normalize a vector/matrix with respect to Frobenius norm
    void normalize(double *A, int &N)
    {
        char NORM = 'F';
        int LDA = 1;
        double WORK[N];
        double norm;
        norm = dlange_(&NORM, &LDA, &N, A, &LDA, WORK);
        for (int i = 0; i < N; i++)
        {
            A[i] /= norm;
        };
    }

    // Matrix/vector product:  y := alpha*A*x + beta*y
    void linEq(double *A, double *X, double *Y, double &alpha, double beta, int &N)
    {
        char TRANS = 'T';
        int INC = 1;
        dgemv_(TRANS, N, N, alpha, A, N, X, INC, beta, Y, INC);
    }

    // Dot product between 2 vectors
    double dot(double *A, double *B, int N)
    {
        int INCX = 1;
        return ddot_(&N, A, &INCX, B, &INCX);
    }

    void minus(double *A, double *B, int N)
    {
        std::transform(A, A + N, B, A, std::minus<double>());
    }

    void plus(double *A, double *B, int N)
    {
        std::transform(A, A + N, B, A, std::plus<double>());
    }

    void plusTimes(double *A, double *B, double c, int N)
    {
        for (int i = 0; i < N; ++i)
            A[i] += B[i] * c;
    }
}

namespace eigen
{

    // Computes the solution to a real system of linear equations
    //         A * X = B,
    // where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
    // On exit, the solution X is stored in B.
    void solve(double *A, double *B, int &N)
    {
        Eigen::Map<Eigen::MatrixXd> A_eigen(A, N, N);
        Eigen::Map<Eigen::VectorXd> B_eigen(B, N);
        B_eigen = A_eigen.lu().solve(B_eigen);
    }

    void inverse(double *A, int &N)
    {
        Eigen::Map<Eigen::MatrixXd> A_eigen(A, N, N);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A_eigen);
        A_eigen = A_eigen.inverse();
    }

    // Normalize a vector/matrix with respect to Frobenius norm
    void normalize(double *A, int &N)
    {
        Eigen::Map<Eigen::VectorXd> A_eigen(A, N);
        A_eigen.normalize();
    }

    // Matrix/vector product:  y := alpha*A*x + beta*y
    void linEq(double *A, double *X, double *Y, double &alpha, double beta, int &N)
    {
        Eigen::Map<Eigen::VectorXd> X_eigen(X, N);
        Eigen::Map<Eigen::VectorXd> Y_eigen(Y, N);
        Eigen::Map<Eigen::MatrixXd> A_eigen(A, N, N);
        Y_eigen = beta * Y_eigen + alpha * A_eigen * X_eigen;
    }

    // Dot product between 2 vectors
    double dot(double *A, double *B, int N)
    {
        Eigen::Map<Eigen::VectorXd> A_eigen(A, N);
        Eigen::Map<Eigen::VectorXd> B_eigen(B, N);
        return A_eigen.dot(B_eigen);
    }

    void minus(double *A, double *B, int N)
    {
        Eigen::Map<Eigen::VectorXd> A_eigen(A, N);
        Eigen::Map<Eigen::VectorXd> B_eigen(B, N);
        A_eigen -= B_eigen;
    }

    void plus(double *A, double *B, int N)
    {
        Eigen::Map<Eigen::VectorXd> A_eigen(A, N);
        Eigen::Map<Eigen::VectorXd> B_eigen(B, N);
        A_eigen += B_eigen;
    }

    void plusTimes(double *A, double *B, double c, int N)
    {
        Eigen::Map<Eigen::VectorXd> A_eigen(A, N);
        Eigen::Map<Eigen::VectorXd> B_eigen(B, N);
        A_eigen += c * B_eigen;
    }

    void cross(double *A, double *B, double *OUT)
    {
        Eigen::Map<Eigen::Vector3d> A_eigen(A);
        Eigen::Map<Eigen::Vector3d> B_eigen(B);
        Eigen::Map<Eigen::Vector3d> OUT_eigen(OUT);
        OUT_eigen = A_eigen.cross(B_eigen);
    }
}
