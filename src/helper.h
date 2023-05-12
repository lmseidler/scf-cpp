#include <vector>
#include <tuple>

#include "eigen-3.4.0/Eigen/Dense"

double distance(Eigen::Vector3d a, Eigen::Vector3d b);
Eigen::VectorXd arr2vec(const double *arr, int size);
double boysf0(double r);
unsigned int binom(unsigned int n, unsigned int k);
Eigen::MatrixXd duplicate(const Eigen::VectorXd &half);
Eigen::MatrixXd unpack(const Eigen::VectorXd &packed);
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> diag_nxn(const Eigen::MatrixXd &mat);