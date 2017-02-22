#include "triangle.hpp"
#include "node.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
void A(Triangle triangles[],int nb_triangles);
void  f();
void local(Triangle triangle,boost::numeric::ublas::matrix<double> M);

int main(){
    return 0;
}

void f(){
  }
void A(Triangle triangles[],int nb_triangles,boost::numeric::ublas::mapped_matrix<double> A){
  boost::numeric::ublas::matrix<double> M (3,3);
  for (int i=0;i<nb_triangles;i++){
    local(triangles[i],M);
    for(int j= 0;j<3;j++){
      for(int k =0;k<3;k++){
        A(triangles[i].node(j).id(),triangles[i].node(k).id()) += M(j,k);
      }
    }
  }
}

void local(Triangle triangle,boost::numeric::ublas::matrix<double> M){

}
