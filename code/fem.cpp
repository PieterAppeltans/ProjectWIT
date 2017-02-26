#include "mesh_reader.hpp"

#include <iostream>
#include <math.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

using namespace boost::numeric::ublas;

/* COMPILE WItH: g++ -Wall -std=c++14 -O3 -lstdc++ -o fem.o fem.cpp
mesh_reader.hpp needs those two flags for reading files */

/* Global variables for every storage scenario */
// DU/DV: sqrt of respectively radial and axial DU/DV squared
const double DU = sqrt(pow(2.8*pow(10, -10), 2) + pow(1.1*pow(10, -9), 2));
const double DV = sqrt(pow(2.32*pow(10, -9), 2) + pow(6.97*pow(10, -9), 2));

/* Global variables for applicable storage scenario */
const double EtA_U = 0.208;
const double EtA_V = 0.0004;
const double t_CEL = 298.15;

int main() {
    /* Read and store mesh information */
    matrix<double> vertices = mesh::read_vertices();
    matrix<double> triangles = mesh::read_triangles(vertices);

    /* Calculate righthand side vector integrals */
    // Done for one triangle, index made a variable to be able to change it easily later
    int i = 0;
    int a = triangles(i, 0);
    int b = triangles(i, 1);
    int c = triangles(i, 2);
    matrix<double> init_C(3, 3);
    vector<double> R(3);
    symmetric_adaptor<matrix<double>, lower> C(init_C);
    C(0, 0) = 6*vertices(a, 0) + 2*vertices(b, 0) + 2*vertices(c, 0);
    C(1, 0) = 2*vertices(a, 0) + 2*vertices(b, 0) + vertices(c, 0);
    C(2, 0) = 2*vertices(a, 0) + vertices(b, 0) + 2*vertices(c, 0);
    C(1, 1) = 2*vertices(a, 0) + 6*vertices(b, 0) + 2*vertices(c, 0);
    C(1, 2) = vertices(a, 0) + 2*vertices(b, 0) + 2*vertices(c, 0);
    C(2, 2) = 2*vertices(a, 0) + 2*vertices(b, 0) + 6*vertices(c, 0);
    C *= (triangles(i,3)/60)*(1/DU);
    // then multiplicate with R column vector

    /* Calculate stiffness matrix */

    /* Calculate Jacobian */

    return 0;
}

double Cvamb(double &t) {
  return (101300*EtA_U)/(8.314*t);
}

double Cuamb(double &t) {
  return (101300*EtA_V)/(8.314*t);
}

double Vmfv(double &t) {
  return 1.61*pow(10, -4)*exp((56700/8.314)*(1/293.15 - 1/t));
}

double Vmu(double &t) {
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

double Ru(double &u, double &v, double &t) {
  return (Vmu(t)*u)/((0.4103+u)*(1+(v/27.2438)));
}

double Rv(double &u, double &v, double &t) {
  return 0.97*Ru(u, v, t)+ Vmfv(t)/(1+(u/0.1149));Z
}
