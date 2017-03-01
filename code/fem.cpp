#include "mesh_reader.hpp"

#include <iostream>
#include <math.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/operation.hpp>


using namespace boost::numeric::ublas;

/* COMPILE WItH: g++ -Wall -std=c++14 -O3 -lstdc++ -o fem.o fem.cpp
mesh_reader.hpp needs those two flags for reading files */

/* Global variables for every storage scenario */
// DU/DV: sqrt of respectively radial and axial DU/DV squared
const double DU_R = 2.8*pow(10, -10);
const double DU_Z = 1.1*pow(10, -9);
const double DV_R = 2.32*pow(10, -9);
const double DV_Z = 6.97*pow(10, -9);

/* Global variables for applicable storage scenario */
const double ETA_U = 0.208;
const double ETA_V = 0.0004;
const double T_CEL = 298.15;

double Cvamb(double const &t) {
  return (101300*ETA_U)/(8.314*t);
}

double Cuamb(double const &t) {
  return (101300*ETA_V)/(8.314*t);
}

double Vmfv(double const &t) {
  return 1.61*pow(10, -4)*exp((56700/8.314)*(1/293.15 - 1/t));
}

double Vmu(double const &t) {
  return 2.39*pow(10, -4)*exp((80200/8.314)*(1/293.15 - 1/t));
}

double Ru_simple(double &u, double &t) {
  /* To calculate starting value of u for solving
  nonlinear system of equations */
  return Vmu(t)*u/0.4103;
}

double Rv_simple(double u, double t) {
  /* To calculate starting value of v for solving
  nonlinear system of equations */
  return 0.97*Ru_simple(u, t) + Vmfv(t);
}

class Ru {

  private:
    double t;

  public:
    Ru(double ti) {
      t = ti;
    }
    double operator()(double u, double v) const {
      return (Vmu(t)*u)/((0.4103+u)*(1+(v/27.2438)));
    }

};

class Rv {

  private:
    double t;

  public:
    Rv(double ti) {
      t = ti;
    }
    double operator()(double u, double v) const {
      Ru Ru_funct(t);
      return 0.97*Ru_funct(u, v)+ Vmfv(t)/(1+(u/0.1149));
    }

};

int main() {
  /* Read and store mesh information */
  matrix<double> vertices = mesh::read_vertices();
  matrix<double> triangles = mesh::read_triangles(vertices);

  /* Calculate righthand side vector integrals */
  matrix<double> init_C(vertices.size1(), vertices.size1());
  vector<double> R(vertices.size1());
  vector<double> S(vertices.size1());
  symmetric_adaptor<matrix<double>, lower> C(init_C);
  for (unsigned t = 0; t < triangles.size1(); ++t) {
    int a = triangles(t, 0);
    int b = triangles(t, 1);
    int c = triangles(t, 2);
    double area = triangles(t, 3);
    C(a, a) += area*(6*vertices(a, 0) + 2*vertices(b, 0) + 2*vertices(c, 0));
    C(b, a) += area*(2*vertices(a, 0) + 2*vertices(b, 0) + vertices(c, 0));
    C(c, a) += area*(2*vertices(a, 0) + vertices(b, 0) + 2*vertices(c, 0));
    C(b, b) += area*(2*vertices(a, 0) + 6*vertices(b, 0) + 2*vertices(c, 0));
    C(b, c) += area*(vertices(a, 0) + 2*vertices(b, 0) + 2*vertices(c, 0));
    C(c, c) += area*(2*vertices(a, 0) + 2*vertices(b, 0) + 6*vertices(c, 0));
  }
  C *= (1/60);
  // multiplied by vector R nonlinear in u and v

  /* Calculate stiffness matrix A */
  matrix<double> init_A(vertices.size1(), vertices.size1());
  matrix<double> G(3, 2);
  matrix<double> GGT(3, 3);
  symmetric_adaptor<matrix<double>, lower> A(init_A);
  for (unsigned t = 0; t < triangles.size1(); ++t) {
    int a = triangles(t, 0);
    int b = triangles(t, 1);
    int c = triangles(t, 2);
    // double area = triangles(t, 3);
    G(0, 0) = (vertices(b, 1) - vertices(c, 1))*sqrt(DU_R);
    G(1, 0) = (vertices(c, 1) - vertices(a, 1))*sqrt(DU_R);
    G(2, 0) = (vertices(a, 1) - vertices(b, 1))*sqrt(DU_R);
    G(0, 1) = (vertices(c, 0) - vertices(b, 0))*sqrt(DU_Z);
    G(1, 1) = (vertices(a, 0) - vertices(c, 0))*sqrt(DU_Z);
    G(2, 1) = (vertices(b, 0) - vertices(a, 0))*sqrt(DU_Z);
    // GGT = (1/(2*area))*((vertices(a, 0)+vertices(b, 0)+vertices(c, 0))/6)*block_prod<matrix_type, 64> (G, trans(G));
    A(a, a) += G(0, 0);
    A(b, a) += G(1, 0);
    A(c, a) += G(2, 0);
    A(b, b) += G(1, 1);
    A(b, c) += G(1, 2);
    A(c, c) += G(2, 2);
  }

  /* Calculate Jacobian */

  return 0;
}
