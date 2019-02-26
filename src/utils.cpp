//
// Created by pierre-olivier on 26/02/19.
//

#include <iostream>
#include <vector>
#include "utils.h"

namespace lapack {

    extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
    }

    void inverse(double* A, int N)
    {
        int *IPIV = new int[N+1];
        int LWORK = N*N;
        double *WORK = new double[LWORK];
        int INFO;

        dgetrf_(&N,&N,A,&N,IPIV,&INFO);
        dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

        delete IPIV;
        delete WORK;
    }
}

namespace linealg {
    void printMatrix(double* A, int N) {
        for(int i =0; i<N; ++i){
            std::cout << "| ";
            for(int j =0; j<N; ++j){
                std::cout << A[i*N+j] << " ";
            }
            std::cout << "|" << std::endl;
        }
    }

    void printVector(std::vector<double> &A) {
        for (std::vector<double>::const_iterator i = A.begin(); i != A.end(); ++i)
            std::cout << *i << ' ';
        std::cout << std::endl;
    }
}
