#ifndef newton_raphson_hpp
#define newton_raphson_hpp 3141592


#include <iostream>

void newton_raphson(VectorXd & F,SpMat & J,VectorXd& x0,double tol){
  int it = 0;
  int maxit = 1000;
  double res = tol+1;
  VectorXd x1;
  std::cout << "Residu:" << res << std::endl;
  while(res > tol && it <maxit){
  SparseLU<SpMat> solver;
    x1 = x0- solver.compute(J(x0)).solve(F(x0));
    res = (x1-x0).squaredNorm();
    x0 = x1;
    std::cout << "Residu: "<< res << std::endl;
    it += 1;
    }

}

#endif
