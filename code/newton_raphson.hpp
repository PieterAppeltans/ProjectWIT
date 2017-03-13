#ifndef newton_raphson_hpp
#define newton_raphson_hpp 3141592
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <iostream>
using namespace boost::numeric::ublas;
template<typename U, typename V>void newton_raphson(U const & F,V const & J,vector<double> const & x0,double tol){
  double res = norm_1(F(x0));
  while(res > tol){
    x0 = x0 - J(x0)^(-1)*F(x0);
    res = norm_1(F(x0));
    std::cout << "Residu: "<< res << std::endl;
  }
}
#endif
