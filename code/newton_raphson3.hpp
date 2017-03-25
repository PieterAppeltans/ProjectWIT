#ifndef newton_raphson_hpp
#define newton_raphson_hpp 3141592


#include <iostream>

template<typename U, typename V>void newton_raphson(U& F,V& J,VectorXd& x0,double tol){
  int it = 0;
  int maxit = 1000;
  double res = tol+1;
  VectorXd x1,step;
  std::cout << "Residu:" << res << std::endl;
  SparseLU<SpMat> solver;
  while(res > tol && it <maxit){
    solver.compute(J(x0));
    x1 = x0- solver.solve(F(x0));
    res = (x1-x0).squaredNorm();
    x0 = x1;
    std::cout << "Residu: "<< res << std::endl;
    it += 1;
    }

}

#endif
