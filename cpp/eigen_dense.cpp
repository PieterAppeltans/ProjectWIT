#include <iostream>
#include <Eigen/Dense>
#include <chrono>

<<<<<<< HEAD:code/test_eigen.cpp
// Compile command:
// g++ -std=c++14 test_eigen.cpp -O3 -DNDEBUG -o test.o -I/usr/local/include/eigen3
// Run command:
// ./test.o filename T(c) eta_u(%) eta_v(%)

// u = Cu = concentration O2, v = Cv = concentration CO2
=======
#include "mesh_reader_eigen.hpp"
#include "constants_dense.hpp"
#include "newton_raphson_dense.hpp"
>>>>>>> 660a1e9097d39d5cd40de087f10a291fbb31910e:cpp/eigen_dense.cpp

using Eigen::MatrixXd;
typedef Matrix<double, Dynamic, 1> VectorXd;

// Compile command:
// g++ -Wall -std=c++14 test_eigen.cpp -o test.o -I/usr/local/include/eigen3 && ./dense.o

// u = Cu = concentration O2, v = Cv = concentration CO2


class F
{
  /* Returns tuple containing two functors, one for the original expression and one for the Jacobian. */

  private:
    MatrixXd& Au_;
    MatrixXd& Av_;
    MatrixXd& B_;
    MatrixXd& C_;
    VectorXd& D_;

  public:
    F(MatrixXd& Au, MatrixXd& Av,MatrixXd& B,MatrixXd& C,VectorXd& D)
    :Au_(Au), Av_(Av),B_(B), C_(C),D_(D)
    {
    }

    VectorXd operator()(VectorXd& x)
    {
      int n = Au_.rows();
      VectorXd u = x.head(n);
      VectorXd v = x.tail(n);
      VectorXd func(2*n);
      func << Au_*u + B_*Ru(u,v) + hu*(C_*u - D_*uamb),Av_*v - B_*Rv(u,v) + hv*(C_*v - D_*vamb);
      return func;
    }
};


class J
{
  /* Returns tuple containing two functors, one for the original expression and one for the Jacobian. */

  private:
    MatrixXd& Au_;
    MatrixXd& Av_;
    MatrixXd& B_;
    MatrixXd& C_;
    VectorXd& D_;

  public:
    J(MatrixXd& Au, MatrixXd& Av,MatrixXd& B,MatrixXd& C,VectorXd& D)
    :Au_(Au), Av_(Av),B_(B), C_(C),D_(D)
    {
    }

    MatrixXd operator()(VectorXd& x)
    {
      int n = Au_.rows();
      VectorXd u = x.head(n);
      VectorXd v = x.tail(n);
      MatrixXd func(2*n,2*n);
      func << Au_ + B_*dRudu(u,v) + hu*C_,B_*dRudv(u,v),-B_*dRvdu(u,v),Av_ - (B_*dRvdv(u,v)) + hv*C_;
      return func;
    }
};


int main(int argc, char *argv[])
{
  setConstants(atof(argv[2])+273.15,atof(argv[3])/100,atof(argv[4])/100);
  auto t1 = std::chrono::high_resolution_clock::now();
  std::string file_name = argv[1];
  std::string location = "../triangle/"+ file_name +".1";
  MatrixXd vertices = mesh::read_vertices(location+".node");
  MatrixXd triangles = mesh::read_triangles(vertices,location+".ele");
  MatrixXi boundaries = mesh::read_boundaries(vertices,location+".poly");
  auto t2 = std::chrono::high_resolution_clock::now();

  std::cout << "Input data successfully read:" << std::endl;
  std::cout << "This took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " milliseconds" << std::endl;
  std::cout << "Number of vertices :" << vertices.rows() << std::endl;
  std::cout << "Number of triangles:" << triangles.rows() <<std::endl;
  t1 = std::chrono::high_resolution_clock::now();

  /* Calculate coefficient matrix B related to the respiration kinetics,
  part of right hand side in system of nonlinear equations */

  MatrixXd B = MatrixXd::Zero(vertices.rows(), vertices.rows());
  for (unsigned t = 0; t < triangles.rows(); ++t)
  {
    int a = triangles(t, 0);
    int b = triangles(t, 1);
    int c = triangles(t, 2);
    double area = triangles(t, 3);
    B(a, a) += area*(6.*vertices(a, 0) + 2.*vertices(b, 0) + 2.*vertices(c, 0));
    B(b, a) += area*(2.*vertices(a, 0) + 2.*vertices(b, 0) + vertices(c, 0));
    B(a,b) += area*(2.*vertices(a, 0) + 2.*vertices(b, 0) + vertices(c, 0));
    B(c, a) += area*(2.*vertices(a, 0) + vertices(b, 0) + 2.*vertices(c, 0));
    B(a,c) += area*(2.*vertices(a, 0) + vertices(b, 0) + 2.*vertices(c, 0));
    B(b, b) += area*(2.*vertices(a, 0) + 6.*vertices(b, 0) + 2.*vertices(c, 0));
    B(b, c) += area*(vertices(a, 0) + 2.*vertices(b, 0) + 2.*vertices(c, 0));
    B(c,b) += area*(vertices(a, 0) + 2.*vertices(b, 0) + 2.*vertices(c, 0));
    B(c, c) += area*(2.*vertices(a, 0) + 2.*vertices(b, 0) + 6.*vertices(c, 0));
  }
  B *= (1./60.);

  t2 = std::chrono::high_resolution_clock::now();
  std::cout << "B matrix successfully assembled" << std::endl;
  std::cout << "This took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " milliseconds" << std::endl;

  /* Calculate the first part of the stiffness matrix in the lefthand
  side of the the nonlinear system, A */

  t1 = std::chrono::high_resolution_clock::now();
  MatrixXd A_U = MatrixXd::Zero(vertices.rows(), vertices.rows());
  MatrixXd A_V = MatrixXd::Zero(vertices.rows(), vertices.rows());
  Eigen::Matrix<double,3,2> G;
  Eigen::Matrix<double,3,3> GGT_U;
  Eigen::Matrix<double,3,3> GGT_V;
  Eigen::Matrix<double,2,2> I_U;
  Eigen::Matrix<double,2,2> I_V;
  I_U(0,0) = DU_R;
  I_U(1,0) = 0;
  I_U(0,1) = 0;
  I_U(1,1) = DU_Z;
  I_V(0,0) = DV_R;
  I_V(1,1) = DV_Z;
  I_V(1,0) = 0;
  I_V(0,1) = 0;
  for (unsigned t = 0; t < triangles.rows(); ++t)
  {
    int a = triangles(t, 0);
    int b = triangles(t, 1);
    int c = triangles(t, 2);
    double area = triangles(t, 3);
    G(0, 0) = (vertices(b, 1) - vertices(c, 1));
    G(1, 0) = (vertices(c, 1) - vertices(a, 1));
    G(2, 0) = (vertices(a, 1) - vertices(b, 1));
    G(0, 1) = (vertices(c, 0) - vertices(b, 0));
    G(1, 1) = (vertices(a, 0) - vertices(c, 0));
    G(2, 1) = (vertices(b, 0) - vertices(a, 0));
    GGT_U = (1/(2*area))*((vertices(a, 0)+vertices(b, 0)+vertices(c, 0))/6)*(G*I_U*G.transpose());
    GGT_V = (1/(2*area))*((vertices(a, 0)+vertices(b, 0)+vertices(c, 0))/6)*(G*I_V*G.transpose());
    A_U(a, a) += GGT_U(0, 0);
    A_U(b, a) += GGT_U(1, 0);
    A_U(a, b) += GGT_U(1, 0);
    A_U(c, a) += GGT_U(2, 0);
    A_U(a, c) += GGT_U(2, 0);
    A_U(b, b) += GGT_U(1, 1);
    A_U(b, c) += GGT_U(1, 2);
    A_U(c, b) += GGT_U(1, 2);
    A_U(c, c) += GGT_U(2, 2);
    A_V(a, a) += GGT_V(0, 0);
    A_V(b, a) += GGT_V(1, 0);
    A_V(a, b) += GGT_V(1, 0);
    A_V(c, a) += GGT_V(2, 0);
    A_V(a, c) += GGT_V(2, 0);
    A_V(b, b) += GGT_V(1, 1);
    A_V(b, c) += GGT_V(1, 2);
    A_V(c, b) += GGT_V(1, 2);
    A_V(c, c) += GGT_V(2, 2);
  }

  t2 = std::chrono::high_resolution_clock::now();
  std::cout << "A matrices assembled." << std::endl;
  std::cout << "This took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " milliseconds" << std::endl;

  /* Calculate the second part of the stiffness matrix (C) and the
  second part of the righthand side (D) */

  MatrixXd C = MatrixXd::Zero(vertices.rows(), vertices.rows());
  VectorXd D = VectorXd::Zero(vertices.rows());
  for (unsigned b = 0; b < boundaries.rows(); ++b)
  {
    double len = sqrt(pow(vertices(boundaries(b, 0), 0) - vertices(boundaries(b, 1), 0), 2) +
      pow(vertices(boundaries(b, 0), 1) - vertices(boundaries(b, 1), 1), 2));
    C(boundaries(b, 0), boundaries(b, 0)) += len*(vertices(boundaries(b, 0), 0)/4 + vertices(boundaries(b, 1), 0)/12);
    C(boundaries(b, 0), boundaries(b, 1)) += len*(vertices(boundaries(b, 0), 0)/12 + vertices(boundaries(b, 1), 0)/12);
    C(boundaries(b, 1), boundaries(b, 0)) += len*(vertices(boundaries(b, 0), 0)/12 + vertices(boundaries(b, 1), 0)/12);
    C(boundaries(b, 1), boundaries(b, 1)) += len*(vertices(boundaries(b, 0), 0)/12 + vertices(boundaries(b, 1), 0)/4);
    D(boundaries(b,0)) += len*(vertices(boundaries(b,0),0)/3.+vertices(boundaries(b,1),0)/6.);
    D(boundaries(b,1)) += len*(vertices(boundaries(b,0),0)/6.+vertices(boundaries(b,1),0)/3.);
  }

  t2 = std::chrono::high_resolution_clock::now();

  std::cout << std::endl << std::endl << std::endl << std::endl << std::endl;
  std::cout << "C matrix and D vector assembled." << std::endl;
  std::cout << "This took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " milliseconds" << std::endl;

  /* Solve linearized problem with ambient concentrations as input_field
  to obtain starting concentration values, then solve the nonlinear
  system with Newton-Raphson */

  t1 = std::chrono::high_resolution_clock::now();

  t1 = std::chrono::high_resolution_clock::now();
  F F_funct(A_U, A_V, B, C, D);
  J J_funct(A_U, A_V, B, C, D);
  t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Functors are created" << std::endl;
  std::cout << "This took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " milliseconds" << std::endl;
  t1 = std::chrono::high_resolution_clock::now();

  VectorXd guess(vertices.rows()*2);
  VectorXd u_0 = (A_U+(Vmu/Kmu)*B+hu*C).colPivHouseholderQr().solve(hu*D*uamb);
  VectorXd v_0 = (A_V+hv*C).colPivHouseholderQr().solve(rq*(Vmu/Kmu)*B*u_0+hv*vamb*D);

  guess << u_0,v_0;
  t2 = std::chrono::high_resolution_clock::now();
  std::cout<< "Initial guess calculated" << std::endl;
  std::cout << "This took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " milliseconds" << std::endl;

  #ifdef DEBUG
  std::cout << "A_U =" << std::endl;
  std::cout << A_U << std::endl;
  std::cout << "A_V =" << std::endl;
  std::cout << A_V << std::endl;
  std::cout << "B =" << std::endl;
  std::cout << B << std::endl;
  std::cout << "C =" << std::endl;
  std::cout << C << std::endl;
  std::cout << "D =" << std::endl;
  std::cout << D << std::endl;
  std::cout << "u_0" << std::endl;
  std::cout << u_0 << std::endl;
  std::cout << "v_0" << std::endl;
  std::cout << v_0 <<std::endl;
  std::cout << "Initial function value = " << std::endl;
  std::cout << F_funct(guess) <<std::endl;
  std::cout << "Initial Jacobian = " << std::endl;
  std::cout << J_funct(guess);
  #endif

  newton_raphson(F_funct,J_funct,guess,pow(10,-17));

  /* Write out the result for python matplotlib code */

  int n = vertices.rows();
  VectorXd u = guess.head(n);
  VectorXd v = guess.tail(n);
  mesh::write_result(u,v);

  return 0;
}
