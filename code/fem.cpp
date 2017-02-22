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

void A(Triangle triangles[], int nb_triangles);
void  f();
void local(Triangle triangle, matrix<double> M);

int main(){
    matrix<double> triangles = mesh::read_triangles();
    matrix<double> vertices = mesh::read_vertices();
    std::cout << triangles << std::endl;
    std::cout << vertices << std::endl;
    return 0;
}

void f(){
  }
void A(Triangle triangles[],int nb_triangles, mapped_matrix<double> A){
  matrix<double> M (3,3);
  for (int i=0;i<nb_triangles;i++){
    local(triangles[i],M);
    for(int j= 0;j<3;j++){
      for(int k =0;k<3;k++){
        A(triangles[i].node(j).id(),triangles[i].node(k).id()) += M(j,k);
      }
    }
  }
}

void local(Triangle triangle, matrix<double> M){

}
