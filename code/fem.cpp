#include "mesh_reader.hpp"

#include <iostream>
#include <math.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

using namespace boost::numeric::ublas;

/* COMPILE WITH: g++ -Wall -std=c++14 -O3 -lstdc++ -o fem.o fem.cpp
mesh_reader.hpp needs those two flags for reading files */

int main(){
    /* Read and store mesh information */
    matrix<double> vertices = mesh::read_vertices();
    matrix<double> triangles = mesh::read_triangles(vertices);

    /* Define constants */
    // Du/Dv: sqrt of respectively radial and axial Du/Dv squared
    double Du = sqrt(pow(2.8*pow(10, -10), 2) + pow(1.1*pow(10, -9), 2));
    double Dv = sqrt(pow(2.32*pow(10, -9), 2) + pow(6.97*pow(10, -9), 2));

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
    C *= (triangles(i,3)/60)*(1/Du);

    /* Calculate stiffness matrix */

    /* Calculate Jacobian */

    return 0;
}
