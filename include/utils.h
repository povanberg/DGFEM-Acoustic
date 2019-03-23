#ifndef DGALERKIN_UTILS_H
#define DGALERKIN_UTILS_H

#include <iomanip>
#include <iostream>

namespace lapack {

    void inverse(double *A, int &N);

    void solve(double *A, double *B, int &N);

    void normalize(double *A, int &N);

    double dot(double *A, double *b, int N);

    void linEq(double *A, double *X, double *Y, double &alpha, double beta, int &N);

    void minus(double *A, double *B, int N);

    void plus(double *A, double *B, int N);

    void plusTimes(double *A, double *B, double c, int N);
}

namespace eigen {

    void inverse(double *A, int &N);

    void solve(double *A, double *B, int &N);

    void normalize(double *A, int &N);

    double dot(double *A, double *b, int N);

    void linEq(double *A, double *X, double *Y, double &alpha, double beta, int &N);

    void minus(double *A, double *B, int N);

    void plus(double *A, double *B, int N);

    void plusTimes(double *A, double *B, double c, int N);

    void cross(double *A, double *B, double *OUT);
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

#endif //DGALERKIN_UTILS_H
