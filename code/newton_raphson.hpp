#ifndef newton_raphson_hpp
#define newton_raphson_hpp 3141592

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>

using namespace boost::numeric::ublas;
template<typename U, typename V>void newton_raphson(U & F,V & J,vector<double>& x0,double tol){
  double res = norm_1(F(x0));
  matrix<double> A;
  matrix<double> inverse;
  std::cout << res << std::endl;
  while(res > tol){
    A = J(x0);
    permutation_matrix<std::size_t> pm(A.size1());
    int res = lu_factorize(A,pm);
    if( res != 0 ) break;
   	// create identity matrix of "inverse"
   	inverse.assign(identity_matrix<double>(A.size1()));
   	// backsubstitute to get the inverse
   	lu_substitute(A, pm, inverse);
    x0 = x0 - prod(inverse,F(x0));
    res = norm_1(F(x0));
    std::cout << "Residu: "<< res << std::endl;
  }
}

#endif
