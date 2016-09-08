/**
 * @file lqr_solve.cpp
 * @brief LQR solver for discrete time infinite horizon problems
 * @author Jonas Schlagenhauf
 * @version v0.1
 * @date 2016-09-08
 */

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>

using Eigen::Matrix;
using Eigen::Dynamic;

typedef Matrix<double, Dynamic, Dynamic> Mtx;

/**
 * @brief Computes the LQR gain matrix (usually denoted K) for a discrete time
 * infinite horizon problem.
 *
 * @param A State matrix of the underlying system
 * @param B Input matrix of the underlying system
 * @param Q Weight matrix penalizing the state
 * @param R Weight matrix penalizing the controls
 * @param N Weight matrix penalizing state / control pairs
 * @param K Pointer to the generated matrix (has to be a double/dynamic size
 * matrix!)
 * @param eps Delta between iterations that determines when convergence is
 * reached
 */
bool compGainMatrix(const Mtx &A, const Mtx &B, const Mtx &Q, const Mtx &R,
                    const Mtx &N, Mtx *K, double eps = 1e-15) {
  // check if dimensions are compatible
  if (A.rows() != A.cols() || B.rows() != A.rows() || Q.rows() != Q.cols() ||
      Q.rows() != A.rows() || R.rows() != R.cols() || R.rows() != B.cols() ||
      N.rows() != A.rows() || N.cols() != B.cols()) {
    std::cout << "One or more matrices have incompatible dimensions. Aborting."
              << std::endl;
    return false;
  }

  // precompute as much as possible
  Mtx B_T = B.transpose();
  Mtx Acal = A - B * R.inverse() * N.transpose();
  Mtx Acal_T = Acal.transpose();
  Mtx Qcal = Q - N * R.inverse() * N.transpose();

  // initialize P with Q
  Mtx P = Q;

  // iterate until P converges
  unsigned int numIterations = 0;
  Mtx Pold = P;
  while (true) {
    numIterations++;

    // compute new P
    P = Acal_T * P * Acal -
        Acal_T * P * B * (R + B_T * P * B).inverse() * B_T * P * Acal + Qcal;

    // update delta
    Mtx delta = P - Pold;
    if (fabs(delta.maxCoeff()) < eps) {
      std::cout << "Number of iterations until convergence: " << numIterations
                << std::endl;
      break;
    }
    Pold = P;
  }

  // compute K from P
  *K = (R + B_T * P * B).inverse() * (B_T * P * A + N.transpose());

  return true;
}

/**
 * @brief Little piece of code to test the solver
 *
 * @return 1 if something went wrong, 0 otherwise.
 */
int main() {
  Mtx A;
  Mtx B;
  Mtx Q;
  Mtx R;
  Mtx N;
  Mtx K;
  A.resize(4, 4);
  B.resize(4, 1);
  Q.resize(4, 4);
  R.resize(1, 1);
  N.resize(4, 1);

  A << 0.9904, 0.04772, 0.004251, 0.0007791, -0.3764, 0.9061, 0.167, 0.03211,
      0.002975, -0.004629, 0.9985, 0.04999, 0.1309, -0.1814, -0.06348, 0.9982;
  B << -0.00241, -0.09491, -9.478e-05, -0.0007852;
  Q << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
  R << 100;
  N << 0, 0, 0, 0;

  if (!compGainMatrix(A, B, Q, R, N, &K))
    return 1;

  Eigen::IOFormat fmt(4, 0, ", ", "\n", "[", "]");
  std::cout << "K: " << K.format(fmt) << std::endl;

  return 0;
}
