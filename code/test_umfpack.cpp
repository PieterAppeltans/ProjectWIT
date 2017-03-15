#include <iostream>
    #include <boost/numeric/bindings/traits/ublas_vector.hpp>
    #include <boost/numeric/bindings/traits/ublas_sparse.hpp>
    #include <boost/numeric/bindings/umfpack/umfpack.hpp>
    #include <boost/numeric/ublas/io.hpp>

    namespace ublas = boost::numeric::ublas;
    namespace umf = boost::numeric::bindings::umfpack;

    int main() {

      ublas::compressed_matrix<double, ublas::column_major, 0,
       ublas::unbounded_array<int>, ublas::unbounded_array<double> > A (5,5,12);
      ublas::vector<double> B (5), X (5);

      A(0,0) = 2.; A(0,1) = 3;
      A(1,0) = 3.; A(1,2) = 4.; A(1,4) = 6;
      A(2,1) = -1.; A(2,2) = -3.; A(2,3) = 2.;
      A(3,2) = 1.;
      A(4,1) = 4.; A(4,2) = 2.; A(4,4) = 1.;

      B(0) = 8.; B(1) = 45.; B(2) = -3.; B(3) = 3.; B(4) = 19.;

      umf::symbolic_type<double> Symbolic;
      umf::numeric_type<double> Numeric;

      umf::symbolic (A, Symbolic);
      umf::numeric (A, Symbolic, Numeric);
      umf::solve (A, X, B, Numeric);

      std::cout << X << std::endl;  // output: [5](1,2,3,4,5)
    } 
