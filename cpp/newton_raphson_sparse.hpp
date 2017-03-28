#ifndef newton_raphson_hpp
#define newton_raphson_hpp 3141592


#include <iostream>
#include <Eigen/SparseCholesky>


template<typename U, typename V>
void newton_raphson(U& F,V& J,VectorXd& x0,double tol)
{
  int it = 0;
  int maxit = 1000;
  double res = tol+1;
  VectorXd x1,step;
  SparseLU<SpMat> solver;
  while(res > tol && it <maxit)
  {
    solver.analyzePattern(J(x0));
    solver.factorize(J(x0));
    x1 = x0- solver.solve(F(x0));
    res = (x1-x0).squaredNorm();
    x0 = x1;
    std::cout << "Residu: "<< res << std::endl;
    it += 1;
  }
}

template<typename U, typename V>
void quasi_newton_raphson(U& F,V& J,VectorXd& x0,double tol)
{
  int it = 0;
  int maxit = 1000;
  double res = tol+1;
  VectorXd x1,step;
  SparseLU<SpMat> solver;
  solver.analyzePattern(J(x0));
  solver.factorize(J(x0));
  while(res > tol && it <maxit)
  {
    x1 = x0- solver.solve(F(x0));
    res = (x1-x0).squaredNorm();
    x0 = x1;
    std::cout << "Residu: "<< res << std::endl;
    it += 1;
  }
}


#endif
