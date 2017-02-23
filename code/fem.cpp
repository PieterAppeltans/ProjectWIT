#include "triangle.hpp"
#include "mesh_reader.hpp"
#include "node.hpp"

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

/* COMPILE WITH: g++ -Wall -std=c++14 -O3 -lstdc++ -o fem fem.cpp
mesh_reader.hpp needs those two flags for reading files */

int main(){
    /* Get node coordinates, boundary flag, triangle node id's, triangle areas */
    matrix<double> vertices = mesh::read_vertices();
    matrix<double> triangles = mesh::read_triangles(vertices);
    // std::cout << triangles << std::endl;

    /* Define constants */

    /* Calculate righthand side vector components */

    /* Calculate stiffness matrix */

    /* Calculate Jacobian */

    return 0;
}
