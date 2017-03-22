
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

using Eigen::MatrixXd;
typedef Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic> ArrayXXd;
typedef Eigen::Array<double,Eigen::Dynamic,1> ArrayXd;
int main(int argc, char const *argv[]) {
  Eigen::Matrix<double,2,2> Q;
  Q(0,0) = 2.0;
  Q(1,0) = 3.0;
  Q(0,1) = 4.0;
  Q(1,1) = 5.0;
  ArrayXXd A_a = Q.array();
  std::cout << A_a.inverse() << std::endl;
  /* code */
  return 0;
}
