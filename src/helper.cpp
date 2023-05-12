#define _USE_MATH_DEFINE
#include <cmath>
#include <iostream>

#include "eigen-3.4.0/Eigen/Dense"

#include "helper.h"

// calulates the squared norm of the difference vector of two vectors
double distance(Eigen::Vector3d a, Eigen::Vector3d b) {
    return (a - b).squaredNorm();
}

// convert a C array to a Eigen::Vector
Eigen::VectorXd arr2vec(const double arr[], int size) {
    Eigen::VectorXd converted(size);
    for (int i = 0; i < size; i++) {
        converted[i] = arr[i];
    }
    return converted;
}

// evaluate 0th order boys function
double boysf0(double r) {
    if (r >= 0.05) {
        // analytical definition of 0th order boys function
        return 0.5 * std::sqrt(M_PI / r) * std::erf(std::sqrt(r));
    } else {
        // use 5th order Taylor expansion if r is too small
        return  1.0 - (1.0 / 3.0)          * r
                    + (2.0e-1 / 3.0)       * r*r
                    - 4.761904761904761e-3 * r*r*r
                    + 1.763668430335097e-4 * r*r*r*r
                    - 4.008337341670675e-6 * r*r*r*r*r;
    }
}

// (n k) = n! / (k! * (n - k)!)
unsigned int binom(unsigned int n, unsigned int k) {
    return
      (        k> n  )? 0 :          // out of range
      (k==0 || k==n  )? 1 :          // edge
      (k==1 || k==n-1)? n :          // first
      (     k+k < n  )?              // recursive:
      (binom(n-1,k-1) * n)/k :       //  path to k=1   is faster
      (binom(n-1,k) * n)/(n-k);      //  path to k=n-1 is faster
}

// duplicates a half-vector for symmetric matrix of size nxn
Eigen::MatrixXd duplicate(const Eigen::VectorXd &half, int n) {
    Eigen::MatrixXd dupl = Eigen::MatrixXd::Zero(n * n, binom(n + 1, 2));
    
    // setup duplication matrix
    for (int j = 0; j < n; j++) {
        for (int i = j; i < n; i++) {
            Eigen::VectorXd u = Eigen::VectorXd::Zero(binom(n + 1, 2));
            u(j * n + i - ((j + 1) * j / 2)) = 1.0;
            
            Eigen::MatrixXd T = Eigen::MatrixXd::Zero(n, n);
            T(i, j) = 1.0;
            T(j, i) = 1.0;

            dupl += Eigen::VectorXd(
                Eigen::Map<Eigen::VectorXd>(
                    T.data(), T.cols() * T.rows()
                )
            ) * u.transpose();

        }
    }
    return dupl * half;
}

// unpack a half-vector into a full square matrix
Eigen::MatrixXd unpack(const Eigen::VectorXd &packed) {
    // duplicate the half-vector
    int resize = int(-0.5 + std::sqrt(0.25 + 2 * packed.size()));
    Eigen::VectorXd dupl = duplicate(packed, resize);

    // create new matrix with appropiate size
    Eigen::MatrixXd unpacked = Eigen::MatrixXd::Zero(resize, resize);

    // transform vector back into a matrix
    int ij = 0;
    for (int i = 0; i < resize; i++) {
        for (int j = 0; j < resize; j++) {
            unpacked(i, j) = dupl(ij);
            ij++;
        }
    }
    return unpacked;
}

// diagonalize square matrix, returning the diagonal matrix and eigenvectors
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> diag_nxn(const Eigen::MatrixXd &mat) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(mat);
    Eigen::VectorXd evals = solver.eigenvalues();
    Eigen::MatrixXd evecs = solver.eigenvectors();

    Eigen::MatrixXd evals_mat = Eigen::MatrixXd::Zero(evals.rows(), evals.rows());
    for (int i = 0; i < evals.rows(); i++) {
        evals_mat(i, i) = evals(i);
    }

    return std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>(evals_mat, evecs);
}