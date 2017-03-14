#include "mesh_reader.hpp"
#include "newton_raphson.hpp"
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


/* COMPILE WItH: g++ -Wall -std=c++14 -O3 -lstdc++ -o fem.o fem.cpp
mesh_reader.hpp needs those two flags for reading files */

/* Global variables for every storage scenario */
const double DU_R = 2.8*pow(10, -10);
const double DU_Z = 1.1*pow(10, -9);
const double DV_R = 2.32*pow(10, -9);
const double DV_Z = 6.97*pow(10, -9);
const double rq = 0.97;

/* Global variables for applicable storage scenario */
const double ETA_U = 0.208;
const double ETA_V = 0.0004;
const double T_CEL = 298.15;

const double vamb = (101300*ETA_U)/(8.314*T_CEL);
const double uamb = (101300*ETA_V)/(8.314*T_CEL);
const double Vmfv = 1.61*pow(10, -4)*exp((56700/8.314)*(1/293.15 - 1/T_CEL));
const double Vmu = 2.39*pow(10, -4)*exp((80200/8.314)*(1/293.15 - 1/T_CEL));
const double hu = 7*pow(10,-7);
const double hv = 7.5*pow(10,-7);
const double Kmu = 1.;
const double Kmv = 1.;

vector<double> Ru_simple(vector<double>u) {
  /* To calculate starting value of u for solving nonlinear system of equations */
  return (Vmu/0.4103)*u;
}

vector<double> Rv_simple(vector<double> u) {
  /* To calculate starting value of v for solving nonlinear system of equations */
  return rq*Ru_simple(u);
}

vector<double> Ru(vector<double> u, vector<double> v) {
  return element_div(Vmu*u,element_prod(scalar_vector<double>(u.size(),0.4103)+u,scalar_vector<double>(u.size(),1)+(v/27.2438)));
}

vector<double> Rv(vector<double> u, vector<double> v) {
  return 0.97*Ru(u,v)+ element_div(Vmfv*scalar_vector<double>(u.size(),1),scalar_vector<double>(u.size(),1)+(u/0.1149));
}

matrix<double> diagonalize(vector<double> v) {
  matrix<double> func(v.size(), v.size());
  for (unsigned l = 0; l < v.size(); ++l) {
    func(l, l) = v(l);
  }
  return func;
}

matrix<double> dRudu(vector<double> u, vector<double> v) {
  scalar_vector<double> one(u.size(),1.);
  vector<double> res_vec = element_div(Vmu*0.4103*one, element_prod(one+v/27.2438, element_prod(one*0.4103+u, one*0.4103+u)));
  return diagonalize(res_vec);
}

matrix<double> dRudv(vector<double> u, vector<double> v) {
  scalar_vector<double> one(u.size(),1.);
  vector<double> res_vec = element_div(-1*27.2438*Vmu*u, element_prod(0.1149*one+u, element_prod(27.2438*one+v, 27.2438*one+v)));
  return diagonalize(res_vec);
}

matrix<double> dRvdu(vector<double> u, vector<double> v) {
  scalar_vector<double> one(u.size(),1);
  vector<double> res_vec = element_div(Vmu*0.4103*one, element_prod(one+v/27.2438, element_prod(one*0.4103+u, one*0.4103+u)));
  return rq*dRudu(u, v) - diagonalize(element_div(0.1149*Vmfv*one, element_prod(0.1149*one+u, 0.1149*one+u)));
}

matrix<double> dRvdv(vector<double> u, vector<double> v) {
  return rq*dRudv(u, v);
}

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
    F(matrix<double> Au, matrix<double> Av, matrix<double> B, matrix<double>C, vector<double>D) {
      Au_ = Au; Av_ = Av; B_ = B; C_ = C; D_ = D;
    }

    vector<double> operator()(vector<double> x) {
      /* x = [uT vT]T */
      int n = Au_.size1();
      vector<double> u = project(x,range(0,n));
      vector<double> v = project(x,range(n,2*n));
      vector<double> func(2*n);
      project(func,range(0,n)) = prod(Au_,u) + prod(B_,Ru(u,v)) + hu*(prod(C_,u) - D_*uamb);
      project(func,range(n,2*n)) = -prod(Av_,v) + prod(B_,Rv(u,v)) - hv*(prod(C_,v) - D_*vamb);
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
    J(matrix<double> Au, matrix<double> Av, matrix<double> B, matrix<double>C, vector<double>D) {
      Au_ = Au; Av_ = Av; B_ = B; C_ = C; D_ = D;
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
  for (unsigned t = 0; t < triangles.size1(); ++t) {
    int a = triangles(t, 0);
    int b = triangles(t, 1);
    int c = triangles(t, 2);
    double area = triangles(t, 3);
    B(a, a) += area*(6*vertices(a, 0) + 2*vertices(b, 0) + 2*vertices(c, 0));
    B(b, a) += area*(2*vertices(a, 0) + 2*vertices(b, 0) + vertices(c, 0));
    B(c, a) += area*(2*vertices(a, 0) + vertices(b, 0) + 2*vertices(c, 0));
    B(b, b) += area*(2*vertices(a, 0) + 6*vertices(b, 0) + 2*vertices(c, 0));
    B(b, c) += area*(vertices(a, 0) + 2*vertices(b, 0) + 2*vertices(c, 0));
    B(c, c) += area*(2*vertices(a, 0) + 2*vertices(b, 0) + 6*vertices(c, 0));
  }
  B *= (1/60);
  t2 = std::chrono::high_resolution_clock::now();
  std::cout << "B matrix successfully assembled" << std::endl;
  std::cout << "This took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " milliseconds" << std::endl;

  /* Calculate first part of stiffness matrix A */
  t1 = std::chrono::high_resolution_clock::now();
  matrix<double> init_A(vertices.size1(), vertices.size1());
  matrix<double> I_U(2,2);
  matrix<double> I_V(2,2);
  I_U(0,0) = DU_R;
  I_U(1,1) = DU_Z;
  I_V(0,0) = DV_R;
  I_V(1,1) = DV_Z;
  matrix<double> G(3, 2);
  matrix<double> GGT_U(3, 3);
  matrix<double> GGT_V(3, 3);
  symmetric_adaptor<matrix<double>, lower> A_U(init_A);
  symmetric_adaptor<matrix<double>, lower> A_V(init_A);
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
    GGT_U = (1/(2*area))*((vertices(a, 0)+vertices(b, 0)+vertices(c, 0))/6)*block_prod<matrix<double>, 64>
      (block_prod<matrix<double>,64>(G,I_U), trans(G));
    GGT_V = (1/(2*area))*((vertices(a, 0)+vertices(b, 0)+vertices(c, 0))/6)*block_prod<matrix<double>, 64>
      (block_prod<matrix<double>,64>(G,I_V), trans(G));
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
  scalar_vector<double> one(vertices.size1(),1);
  vector<double> guess(vertices.size1()*2);
  project(guess,range(0,vertices.size1())) = Ru_simple(uamb*one);
  project(guess,range(vertices.size1(),2*vertices.size1())) = Rv_simple(uamb*one);
  t2 = std::chrono::high_resolution_clock::now();
  std::cout<< "Initial guess calculated" << std::endl;
  std::cout << "This took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " milliseconds" << std::endl;
  newton_raphson(F_funct,J_funct,guess,pow(10,-7));
  //std::cout << guess << std::endl;
  return 0;
}
