#ifndef newton_raphson_hpp
#define newton_raphson_hpp 3141592

#include "constants.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>

using namespace boost::numeric::ublas;

template<class T>
bool InvertMatrix (const matrix<T>& input, matrix<T>& inverse) {
  typedef permutation_matrix<std::size_t> pmatrix;
  matrix<T> A(input);
  pmatrix pm(A.size1());
  int res = lu_factorize(A,pm);
  if( res != 0 ) return false;
  inverse.assign(identity_matrix<T>(A.size1()));
  lu_substitute(A, pm, inverse);
  return true;
}

template<typename U, typename V>void newton_raphson(U & F,V & J,vector<double>& x0,double tol){
  double res = norm_1(F(x0));
  //std::cout << "F(x0):" << F(x0) << std::endl;
  int size = x0.size();
  matrix<double> A;
  matrix<double> inverse(size,size);
  inverse = set_zero(inverse);
  std::cout << "Residu:" << res << std::endl;
  while(res > tol){
    A = J(x0);
    std::cout<< "A: " << A << std::endl;
    InvertMatrix(A, inverse);
    // //std::cout << "Inverse size " << inverse.size1() << " " << inverse.size2() << std::endl;
    // x0 = x0 - prod(inverse,F(x0));
    // res = norm_1(F(x0));
    // std::cout << "Residu: "<< res << std::endl;
    // std::cout << "Solution: "<< x0 << std::endl;
  }

}

#endif
