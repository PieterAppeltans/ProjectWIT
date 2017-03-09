#include "mesh_reader.hpp"

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
#include <boost/math/tools/roots.hpp>

using namespace boost::numeric::ublas;
using namespace boost::math::tools;
/* COMPILE WItH: g++ -Wall -std=c++14 -O3 -lstdc++ -o fem.o fem.cpp
mesh_reader.hpp needs those two flags for reading files */

/* Global variables for every storage scenario */
// DU/DV: sqrt of respectively radial and axial DU/DV squared
const double DU_R = 2.8*pow(10, -10);
const double DU_Z = 1.1*pow(10, -9);
const double DV_R = 2.32*pow(10, -9);
const double DV_Z = 6.97*pow(10, -9);
const double rq = 0.97;

/* Global variables for applicable storage scenario */
const double ETA_U = 0.208;
const double ETA_V = 0.0004;
const double T_CEL = 298.15;

const double Cvamb = (101300*ETA_U)/(8.314*T_CEL);
const double Cuamb = (101300*ETA_V)/(8.314*T_CEL);
const double Vmfv = 1.61*pow(10, -4)*exp((56700/8.314)*(1/293.15 - 1/T_CEL));
const double Vmu = 2.39*pow(10, -4)*exp((80200/8.314)*(1/293.15 - 1/T_CEL));
const double hu = 7*pow(10,-7);
const double hv = 7.5*pow(10,-7);
const double Kmu = 1.;
const double Kmv = 1.;

double Ru_simple(double u) {
  /* To calculate starting value of u for solving
  nonlinear system of equations */
  return Vmu*u/0.4103;
}

double Rv_simple(double u) {
  /* To calculate starting value of v for solving
  nonlinear system of equations */
  return rq*Ru_simple(u) + Vmfv;
}

class Ru {
  public:
    Ru() {
    }
    vector<double> operator()(vector<double> u, vector<double> v) const {
      return element_div(Vmu*u,element_prod(scalar_vector<double>(u.size(),0.4103)+u,scalar_vector<double>(u.size(),1)+(v/27.2438)));
    }

};

class Rv {
  public:
    Rv() { }
    vector<double> operator()(vector<double> u, vector<double> v) const {
      Ru Ru_funct;
      return 0.97*Ru_funct(u,v)+ element_div(Vmfv*scalar_vector<double>(u.size(),1),scalar_vector<double>(u.size(),1)+(u/0.1149));
    }

};

class F {
  private:
    matrix<double> Au_;
    matrix<double> Av_;
    matrix<double> B_;
    matrix<double> C_;
    vector<double> D_;
    int N_;

  public:
    F(matrix<double> Au,matrix<double> Av, matrix<double> B,matrix<double>C,vector<double>D,int N)
    {
      Au_ =Au;Av_ = Av; B_ = B;C_ = C;D_ =D;N_= N;
    };
    vector<double> operator()(vector<double> x){
      vector<double> result(2*N_);
      vector<double> temp(N_);
      vector<double> u = project(x,range(0,N_));
      vector<double> v = project(x,range(N_,2*N_));
      Ru Ru_funct;
      Rv Rv_funct;
      temp = prod(Au_, u)+ prod(B_,Ru_funct(u,v))+hu*(prod(C_,u)-D_*Cuamb);
      project(result,range(0,N_)) = temp;
      temp = -prod(Av_,v)+ prod(B_,Rv_funct(u,v))-hv*(prod(C_,v)-D_*Cvamb);
      project(result,range(N_,2*N_)) = temp;
      return result;
    };
};

matrix<double> dRudu(vector<double> u,vector<double> v){
  vector<double> t = element_div(scalar_vector<double>(u.size(),Vmu*Kmu),element_prod(element_prod(scalar_vector<double>(u.size(),1)+u,scalar_vector<double>(u.size(),1)+u),scalar_vector<double>(u.size(),1)+v/Kmv));
  return diagonal_matrix<double>(u.size(), t);
};
matrix<double> dRudv(vector<double> u,vector<double> v){
  vector<double> t = element_div(-Vmu*u,element_prod(scalar_vector<double>(u.size(),Kmu)+u,Kmv*element_prod(scalar_vector<double>(u.size(),1.),v/Kmv,scalar_vector<double>(u.size(),1.),v/Kmv)));
  return diagonal_matrix<double>(u.size(),t);
};
matrix<double> dRvdu(vector<double> u,vector<double> v){
  vector<double> dRudu = element_div(scalar_vector<double>(u.size(),Vmu*Kmu),element_prod(element_prod(scalar_vector<double>(u.size(),1)+u,scalar_vector<double>(u.size(),1)+u),scalar_vector<double>(u.size(),1)+v/Kmv));
  vector<double> t = rq*dRudu+Vmfv*element_div(scalar_vector<double>(u.size(),-1.),(Kmfu,element_prod(scalar_vector(u.size(),1)+u/Kmfu,scalar_vector(u.size(),1)+u/Kmfu)));
  return diagonal_matrix<double>(u.size(),t);
}
matrix<double> dRvdv(vector<double> u,vector<double> v){
  vector<double> t = rq*element_div(-Vmu*u,element_prod(scalar_vector<double>(u.size(),Kmu)+u,Kmv*element_prod(scalar_vector<double>(u.size(),1.),v/Kmv,scalar_vector<double>(u.size(),1.),v/Kmv)));
  return diagonal_matrix<double>(u.size(),t);
}
class J {
  private:
    matrix<double> Au_;
    matrix<double> Av_;
    matrix<double> B_;
    matrix<double> C_;
    vector<double> D_;
    int N_;

  public:
    J(matrix<double> Au,matrix<double> Av, matrix<double> B,matrix<double>C,vector<double>D,int N)
    {
      Au_ =Au; Av_ = Av; B_ = B; C_ = C; D_ =D; N_= N;
    };
    matrix<double> operator()(vector<double> u, vector<double> v) {
      matrix<double> result(2*N_,2*N_);
      matrix<double> temp(N_,N_);


      temp = Au_ + B_*dRudu(u,v)+hu*C_;
      project(result,range(0,N_),range(0,N_)) = temp;
      temp = B_*dRudv(u,v);
      project(result,range(0,N_),range(N_,2*N_)) = temp;
      temp = -B_*dRvdu(u,v);
      project(result,range(N_,2*N_),range(0,N_)) = temp;
      temp = A_v-B*dRvdv(u,v)+hv*C;
      project(result,range(N_,2*N_),range(N_,2*N_)) = temp;
      return result;
    };
};

int main() {
  /* Read and store mesh information */
  matrix<double> vertices = mesh::read_vertices();
  matrix<double> triangles = mesh::read_triangles(vertices);
  matrix<int> boundaries = mesh::read_boundaries(vertices);

  /* Calculate righthand side vector integrals */
  matrix<double> init_B(vertices.size1(), vertices.size1());
  vector<double> R(vertices.size1());
  vector<double> S(vertices.size1());
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

  /* Calculate first part of stiffness matrix A */
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
    GGT_U = (1/(2*area))*((vertices(a, 0)+vertices(b, 0)+vertices(c, 0))/6)*block_prod<matrix<double>, 64> (block_prod<matrix<double>,64>(G,I_U), trans(G));
    GGT_V = (1/(2*area))*((vertices(a, 0)+vertices(b, 0)+vertices(c, 0))/6)*block_prod<matrix<double>, 64> (block_prod<matrix<double>,64>(G,I_V), trans(G));
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

  /* Boundary condition integrals: second part of A and a constant vector term */
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
  F F_funct(A_U,A_V,B,C,D,vertices.size1());



  return 0;
}
