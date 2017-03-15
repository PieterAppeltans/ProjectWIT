#include "mesh_reader.hpp"
#include "newton_raphson.hpp"
#include "constants.hpp"

#include <tuple>
#include <iostream>
#include <math.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_blocked.hpp>
#include <chrono>

using namespace boost::numeric::ublas;

/* Configure boost-numeric-bindings in /usr/local to be able to use the compile statement on the next line. */
/* COMPILE WITH: g++ -Wall -std=c++14 -O3 -lstdc++ -o fem.o fem.cpp
mesh_reader.hpp needs those two flags for reading files */

class FJ {
  /* Returns tuple containing two functors, one for the original expression and one for the Jacobian. */

  private:
    matrix<double> Au_;
    matrix<double> Av_;
    matrix<double> B_;
    matrix<double> C_;
    vector<double> D_;

  public:
    FJ(matrix<double> Au, matrix<double> Av, matrix<double> B, matrix<double>C, vector<double>D) {
      Au_ = Au; Av_ = Av; B_ = B; C_ = C; D_ = D;
    }

    std::tuple<vector<double>, matrix<double>> operator()(vector<double> x) {
      /* x = [uT vT]T */
      int n = Au_.size1();
      vector<double> u = project(x,range(0,n));
      vector<double> v = project(x,range(n,2*n));
      vector<double> func(2*n);
      matrix<double> JAC(2*n, 2*n);
      project(func,range(0,n)) = prod(Au_,u) + prod(B_,Ru(u,v)) + hu*(prod(C_,u) - D_*uamb);
      project(func,range(n,2*n)) = -prod(Av_,v) + prod(B_,Rv(u,v)) - hv*(prod(C_,v) - D_*vamb);
      project(JAC,range(0,n),range(0, n)) = Au_ + prod(B_,dRudu(u,v)) + hu*C_;
      project(JAC,range(n,2*n),range(0,n)) = prod(B_,dRudv(u,v));
      project(JAC,range(n,2*n),range(0,n)) = prod(B_,dRvdu(u,v));
      project(JAC,range(n,2*n),range(n,2*n)) = -1*Av_ + prod(B_,dRvdv(u,v)) - hv*C_;
      return std::make_tuple(func, JAC);
    }
};

class F {
  /* Returns tuple containing two functors, one for the original expression and one for the Jacobian. */

  private:
    matrix<double> Au_;
    matrix<double> Av_;
    matrix<double> B_;
    matrix<double> C_;
    vector<double> D_;

  public:
    F(matrix<double> Au, matrix<double> Av, matrix<double> B, matrix<double> C, vector<double> D)
    :Au_(Au), Av_(Av),B_(B), C_(C),D_(D)
    {

    }

    vector<double> operator()(vector<double> x) {
      /* x = [uT vT]T */
      int n = Au_.size1();
      vector<double> u = project(x,range(0,n));
      vector<double> v = project(x,range(n,2*n));
      vector<double> func(2*n);
      project(func,range(0,n)) = prod(Au_,u) + prod(B_,Ru(u,v)) + hu*(prod(C_,u) - D_*uamb);
      project(func,range(n,2*n)) = prod(Av_,v) - prod(B_,Rv(u,v)) + hv*(prod(C_,v) - D_*vamb);
      return func;
    }
};

class J {
  /* Returns tuple containing two functors, one for the original expression and one for the Jacobian. */

  private:
    matrix<double> Au_;
    matrix<double> Av_;
    matrix<double> B_;
    matrix<double> C_;
    vector<double> D_;

  public:
    J(matrix<double> Au, matrix<double> Av, matrix<double> B, matrix<double> C, vector<double> D)
    :Au_(Au), Av_(Av),B_(B), C_(C),D_(D)
    {

    }

    matrix<double> operator()(vector<double> x) {
      /* x = [uT vT]T */
      int n = Au_.size1();
      vector<double> u = project(x,range(0,n));
      vector<double> v = project(x,range(n,2*n));
      matrix<double> JAC(2*n, 2*n);
      project(JAC,range(0,n),range(0, n)) = Au_ + prod(B_,dRudu(u,v)) + hu*C_;
      project(JAC,range(n,2*n),range(0,n)) = prod(B_,dRudv(u,v));
      project(JAC,range(n,2*n),range(0,n)) = prod(B_,dRvdu(u,v));
      project(JAC,range(n,2*n),range(n,2*n)) = -1*Av_ + prod(B_,dRvdv(u,v)) - hv*C_;
      return JAC;
    }
};



int main() {
  /* Read and store mesh information */
  auto t1 = std::chrono::high_resolution_clock::now();
  matrix<double> vertices = mesh::read_vertices();
  matrix<double> triangles = mesh::read_triangles(vertices);
  matrix<int> boundaries = mesh::read_boundaries(vertices);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Input data successfully read:" << std::endl;
  std::cout << "This took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " milliseconds" << std::endl;
  std::cout << "Number of vertices :" << vertices.size1() << std::endl;
  std::cout << "Number of triangles:" << triangles.size1() <<std::endl;
  /* Calculate righthand side vector integrals */
  t1 = std::chrono::high_resolution_clock::now();
  matrix<double> init_B(vertices.size1(), vertices.size1());
  symmetric_adaptor<matrix<double>, lower> B(init_B);
  B = set_zero(B);
  std::cout << "B BEFORE : " << B << std::endl;
  for (unsigned t = 0; t < triangles.size1(); ++t) {
    int a = triangles(t, 0);
    int b = triangles(t, 1);
    int c = triangles(t, 2);
    double area = triangles(t, 3);
    std::cout << std::endl;
    B(a, a) += area*(6.*vertices(a, 0) + 2.*vertices(b, 0) + 2.*vertices(c, 0));
    std::cout<< B(a,a) << std::endl;
    B(b, a) += area*(2.*vertices(a, 0) + 2.*vertices(b, 0) + vertices(c, 0));
    B(c, a) += area*(2.*vertices(a, 0) + vertices(b, 0) + 2.*vertices(c, 0));
    B(b, b) += area*(2.*vertices(a, 0) + 6.*vertices(b, 0) + 2.*vertices(c, 0));
    std::cout << "B(b,c): " << B(b,c) << "deel1: "<< area*(2.*vertices(a, 0) + 6.*vertices(b, 0) + 2.*vertices(c, 0)) << std::endl;
    B(b, c) += area*(vertices(a, 0) + 2.*vertices(b, 0) + 2.*vertices(c, 0));
    B(c, c) += area*(2.*vertices(a, 0) + 2.*vertices(b, 0) + 6.*vertices(c, 0));
    std::cout << B << std::endl;
  }
  std::cout << "before division: " << B << std::endl;
  B *= (1./60.);
  std::cout << "after divison: " << B << std::endl;
  t2 = std::chrono::high_resolution_clock::now();
  std::cout << "B matrix successfully assembled" << std::endl;
  std::cout << "This took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " milliseconds" << std::endl;

  /* Calculate first part of stiffness matrix A */
  t1 = std::chrono::high_resolution_clock::now();
  matrix<double> init_A_U(vertices.size1(), vertices.size1());
  matrix<double> init_A_V(vertices.size1(), vertices.size1());
  matrix<double> I_U(2,2);
  matrix<double> I_V(2,2);
  I_U(0,0) = DU_R;
  I_U(1,1) = DU_Z;
  I_V(0,0) = DV_R;
  I_V(1,1) = DV_Z;
  matrix<double> G(3, 2);
  matrix<double> GGT_U(3, 3);
  matrix<double> GGT_V(3, 3);
  matrix<double> temp(3, 3);
  symmetric_adaptor<matrix<double>, lower> A_U(init_A_U);
  symmetric_adaptor<matrix<double>, lower> A_V(init_A_V);
  A_U = set_zero(A_U);
  A_V = set_zero(A_V);
  for (unsigned t = 0; t < triangles.size1(); ++t) {
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
    temp = prod(G,I_U);
    GGT_U = (1/(2*area))*((vertices(a, 0)+vertices(b, 0)+vertices(c, 0))/6)*prod(temp, trans(G));
    temp = prod(G,I_V);
    GGT_V = (1/(2*area))*((vertices(a, 0)+vertices(b, 0)+vertices(c, 0))/6)*prod(temp, trans(G));
    A_U(a, a) += GGT_U(0, 0);
    A_U(b, a) += GGT_U(1, 0);
    A_U(c, a) += GGT_U(2, 0);
    A_U(b, b) += GGT_U(1, 1);
    A_U(b, c) += GGT_U(1, 2);
    A_U(c, c) += GGT_U(2, 2);
    A_V(a, a) += GGT_V(0, 0);
    A_V(b, a) += GGT_V(1, 0);
    A_V(c, a) += GGT_V(2, 0);
    A_V(b, b) += GGT_V(1, 1);
    A_V(b, c) += GGT_V(1, 2);
    A_V(c, c) += GGT_V(2, 2);
  }
  t2 = std::chrono::high_resolution_clock::now();
  std::cout << "A matrices assembled." << std::endl;
  std::cout << "This took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " milliseconds" << std::endl;
  /* Boundary condition integrals: second part of A and a constant vector term */
  t1 = std::chrono::high_resolution_clock::now();
  matrix<double> init_C(vertices.size1(), vertices.size1());
  symmetric_adaptor<matrix<double>, lower> C(init_C);
  C = set_zero(C);
  vector<double> D(vertices.size1());
  for (unsigned b = 0; b < boundaries.size1(); ++b) {
    double len = sqrt(pow(vertices(boundaries(b, 0), 0) - vertices(boundaries(b, 1), 0), 2) +
      pow(vertices(boundaries(b, 0), 1) - vertices(boundaries(b, 1), 1), 2));
    C(boundaries(b, 0), boundaries(b, 0)) += len*(vertices(boundaries(b, 0), 0)/4 + vertices(boundaries(b, 1), 0)/12);
    C(boundaries(b, 0), boundaries(b, 1)) += len*(vertices(boundaries(b, 0), 0)/12 + vertices(boundaries(b, 1), 0)/12);
    C(boundaries(b, 1), boundaries(b, 1)) += len*(vertices(boundaries(b, 0), 0)/12 + vertices(boundaries(b, 1), 0)/4);
    D(boundaries(b,0)) += len*(vertices(boundaries(b,0),0)/3.+vertices(boundaries(b,1),0)/6.);
    D(boundaries(b,1)) += len*(vertices(boundaries(b,0),0)/6.+vertices(boundaries(b,1),0)/3.);
  }
  std::cout << "A_U: " << A_U << std::endl;
  std::cout << "A_V: " << A_V << std::endl;
  std::cout << "B: " << B << std::endl;
  std::cout << "C: " << C << std::endl;
  std::cout << "D: " << D << std::endl;
  t2 = std::chrono::high_resolution_clock::now();
  std::cout << "C matrix and D vector assembled." << std::endl;
  std::cout << "This took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " milliseconds" << std::endl;
  /* Root finding for nonlinear system of equations */
  t1 = std::chrono::high_resolution_clock::now();
  F F_funct(A_U, A_V, B, C, D);
  J J_funct(A_U, A_V, B, C, D);
  t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Functors are created" << std::endl;
  std::cout << "This took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " milliseconds" << std::endl;
  t1 = std::chrono::high_resolution_clock::now();


  vector<double> guess = scalar_vector<double>(vertices.size1()*2,5.);
  matrix<double> inverse_u_0 (vertices.size1(),vertices.size1());
  InvertMatrix<double>(A_U+(Vmu/Kmu)*B+hu*C, inverse_u_0);
  vector<double> u_0 = prod(inverse_u_0, hu*D*uamb);
  matrix<double> inverse_v_0 (vertices.size1(),vertices.size1());
  InvertMatrix<double>(A_V+hv*C, inverse_v_0);
  vector<double> v_0 = prod(inverse_v_0, rq*(Vmu/Kmu)*prod(B,u_0)+hv*vamb*D);
  project(guess,range(0,vertices.size1())) = v_0;
  project(guess,range(vertices.size1(),2*vertices.size1())) = v_0;

  t2 = std::chrono::high_resolution_clock::now();
  std::cout<< "Initial guess calculated" << std::endl;
  std::cout << "This took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " milliseconds" << std::endl;
  newton_raphson(F_funct,J_funct,guess,pow(10,-7));
  std::cout << guess << std::endl;
  return 0;
}
