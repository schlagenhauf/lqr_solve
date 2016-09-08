#include <eigen3/Eigen/Core>
#include <iostream>

using Eigen::Matrix;
using Eigen::Dynamic;

typedef Matrix<double, Dynamic, Dynamic> Mtx;

void compGainMatrix(Mtx A, Mtx B, Mtx Q, Mtx R, Mtx N) {
  bool converged = false;
  Mtx P = Q;
  while (!converged) {
    P = A.transpose() * P * A - (A.transpose() * P * B + N) * (R + B.transpose()
        * P * B).inverse() * (B.transpose()



  }
}

int main() {

}
