#include <exception>
#include <Rcpp.h>

namespace {
  inline double cuber(const double x) {
    if (x <= 0) {
      return 0;
    } else {
      return x * x * x;
    }
  }
  
  inline double dCuber(const double x) {
    if (x <= 0) {
      return 0;
    } else {
      return 3 * x * x;
    }
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector
basis_vector(const Rcpp::NumericVector& knots,
	     const Rcpp::NumericVector& x) {
  if (knots.size() < 2) {
    throw std::runtime_error("Bad knots.");
  }
  Rcpp::NumericMatrix result(x.size(), knots.size());
  result(Rcpp::_, 0) = Rcpp::rep(1, x.size());
  result(Rcpp::_, 1) = x;

  for (R_xlen_t ind=0; ind < knots.size() - 2; ++ind) {
    const double last = *(knots.end() - 1);
    const double first = *(knots.begin());
    const double lam = (last - knots[ind + 1]) / (last - first);
    result(Rcpp::_, ind + 2) =
      Rcpp::sapply(x - knots[ind + 1], cuber)
      - lam * Rcpp::sapply(x - first, cuber)
      - (1 - lam) * Rcpp::sapply(x - last, cuber);
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix basis_matrix(const Rcpp::NumericMatrix& knots,
				 const Rcpp::NumericVector& x) {
  if (knots.ncol() < 2) {
    throw std::runtime_error("Bad knots.");
  }
  if (knots.nrow() != x.size()) {
    throw std::runtime_error("Mismatch between knots and points");
  }

  Rcpp::NumericMatrix result(x.size(), knots.ncol());

  result(Rcpp::_, 0) = Rcpp::rep(1, x.size());
  result(Rcpp::_, 1) = x;

  for (R_xlen_t ind_r=0; ind_r < result.nrow(); ++ind_r) {
    for (R_xlen_t ind=0; ind < knots.ncol() - 2;++ind) {
      const double last  = knots(ind_r, knots.ncol() - 1);
      const double first = knots(ind_r, 0);
      const double lam = (last - knots(ind_r, ind + 1)) / (last - first);
      result(ind_r, ind + 2) = cuber(x[ind_r] - knots(ind_r, ind + 1))
	- lam * cuber(x[ind_r] - first)
	- (1 - lam) * cuber(x[ind_r] - last);
    }
  }
  
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector
dbasis_vector(const Rcpp::NumericVector& knots,
	      const Rcpp::NumericVector& x) {
  if (knots.size() < 2) {
    throw std::runtime_error("Bad knots.");
  }
  Rcpp::NumericMatrix result(x.size(), knots.size());
  result(Rcpp::_, 0) = Rcpp::rep(0, x.size());
  result(Rcpp::_, 1) = Rcpp::rep(1, x.size());

  for (R_xlen_t ind=0; ind < knots.size() - 2; ++ind) {
    const double last  = *(knots.end() - 1);
    const double first = *(knots.begin());
    const double lam   = (last - knots[ind + 1]) / (last - first);
    result(Rcpp::_, ind + 2) =
      Rcpp::sapply(x - knots[ind + 1], dCuber)
      - lam * Rcpp::sapply(x - first, dCuber)
      - (1 - lam) * Rcpp::sapply(x - last, dCuber);
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix dbasis_matrix(const Rcpp::NumericMatrix& knots,
				  const Rcpp::NumericVector& x) {
  if (knots.ncol() < 2) {
    throw std::runtime_error("Bad knots.");
  }
  if (knots.nrow() != x.size()) {
    throw std::runtime_error("Mismatch between knots and points");
  }

  Rcpp::NumericMatrix result(x.size(), knots.ncol());

  result(Rcpp::_, 0) = Rcpp::rep(0, x.size());
  result(Rcpp::_, 1) = Rcpp::rep(1, x.size());

  for (R_xlen_t ind_r=0; ind_r < result.nrow(); ++ind_r) {
   for (R_xlen_t ind=0; ind < knots.ncol() - 2;++ind) {
      const double last  = knots(ind_r, knots.ncol() - 1);
      const double first = knots(ind_r, 0);
      const double lam   = (last - knots(ind_r, ind + 1)) / (last - first);

      result(ind_r, ind + 2) = dCuber(x[ind_r] - knots(ind_r, ind + 1))
	- lam * dCuber(x[ind_r] - first)
	- (1 - lam) * dCuber(x[ind_r] - last);
    }
  }
  
  return result;
}
