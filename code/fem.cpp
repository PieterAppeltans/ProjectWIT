#include "triangle.hpp"
#include "mesh_reader.hpp"
#include "node.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

/* COMPILE WITH: gcc fem.cpp -std=c++11 -lstdc++
mesh_reader.hpp needs those two flags for reading files */

void A(Triangle triangles[], int nb_triangles);
void  f();
void local(Triangle triangle, matrix<double> M);

int main(){
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
