#include <RcppParallel.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
using namespace RcppParallel;



/********************************************************************************************************/
/***********************************         PARALLEL COMPUTATION     ***********************************/
/********************************************************************************************************/
struct SquareRoot : public Worker
{
  // source matrix
  const RMatrix<double> input;

  // destination matrix
  RMatrix<double> output;

  // initialize with source and destination
  SquareRoot(const Rcpp::NumericMatrix input, Rcpp::NumericMatrix output)
    : input(input), output(output) {}

  // take the square root of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {
    std::transform(input.begin() + begin,
                   input.begin() + end,
                   output.begin() + begin,
                   ::sqrt);

  }
};


// [[Rcpp::export]]
Rcpp::NumericMatrix parallelMatrixSqrt(Rcpp::NumericMatrix x) {

  // allocate the output matrix
  Rcpp::NumericMatrix output(x.nrow(), x.ncol());

  // SquareRoot functor (pass input and output matrixes)
  SquareRoot squareRoot(x, output);

  // call parallelFor to do the work
  parallelFor(0, x.length(), squareRoot);

  // return the output matrix
  return output;
}










// generic function for kl_divergence
template <typename InputIterator1, typename InputIterator2>
inline double kl_divergence(InputIterator1 begin1, InputIterator1 end1,
                            InputIterator2 begin2) {

  // value to return
  double rval = 0;

  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;

  // for each input item
  while (it1 != end1) {

    // take the value and increment the iterator
    double d1 = *it1++;
    double d2 = *it2++;

    // accumulate if appropirate
    if (d1 > 0 && d2 > 0)
      rval += std::log(d1 / d2) * d1;
  }
  return rval;
}

// helper function for taking the average of two numbers
inline double average(double val1, double val2) {
  return (val1 + val2) / 2;
}
struct JsDistance : public Worker {

  // input matrix to read from
  const RMatrix<double> mat;

  // output matrix to write to
  RMatrix<double> rmat;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  JsDistance(const Rcpp::NumericMatrix mat, Rcpp::NumericMatrix rmat)
    : mat(mat), rmat(rmat) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for (std::size_t j = 0; j < i; j++) {

        // rows we will operate on
        RMatrix<double>::Row row1 = mat.row(i);
        RMatrix<double>::Row row2 = mat.row(j);

        // compute the average using std::tranform from the STL
        std::vector<double> avg(row1.length());
        std::transform(row1.begin(), row1.end(), // input range 1
                       row2.begin(),             // input range 2
                       avg.begin(),              // output range
                       average);                 // function to apply

        // calculate divergences
        double d1 = kl_divergence(row1.begin(), row1.end(), avg.begin());
        double d2 = kl_divergence(row2.begin(), row2.end(), avg.begin());

        // write to output matrix
        rmat(i,j) = sqrt(.5 * (d1 + d2));
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_parallel_js_distance(Rcpp::NumericMatrix mat) {

  // allocate the matrix we will return
  Rcpp::NumericMatrix rmat(mat.nrow(), mat.nrow());

  // create the worker
  JsDistance jsDistance(mat, rmat);

  // call it with parallelFor
  parallelFor(0, mat.nrow(), jsDistance);

  return rmat;
}
