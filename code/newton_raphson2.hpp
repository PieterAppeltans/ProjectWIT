#ifndef newton_raphson_hpp
#define newton_raphson_hpp 3141592


#include <iostream>

template<typename U, typename V>void newton_raphson(U & F,V & J,VectorXd& x0,double tol){
  int size = x0.rows();
  double res = tol+1;
  VectorXd x1;
  std::cout << "Residu:" << res << std::endl;
  while(res > tol){
    //x1= F(x0);
    MatrixXd A = J(x0);
    /*x1 = x0- J(x0).colPivHouseholderQr().solve(F(x0));
    res = (x0-x1).squaredNorm();
    x0 = x1;
    std::cout << "Residu: "<< res << std::endl;
    std::cout << "Solution: "<< x0 << std::endl;
    */
  }

}

#endif
