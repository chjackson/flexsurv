#include <Rcpp.h>

#include "distribution.h"
#include "rep_len.h"
#include "gengamma.h"
#include "mapply.h"

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector
dgengamma_work(const Rcpp::NumericVector& x,
	       const Rcpp::NumericVector& mu,
	       const Rcpp::NumericVector& sigma,
	       const Rcpp::NumericVector& Q,
	       const bool log) {
  if (x.size() == 0) { return x; }
  
  const R_xlen_t size
    = std::max(std::max(sigma.size(),
			Q.size()),
	       std::max(x.size(),
			mu.size()));
  return perhaps_exp(mapply(flexsurv::rep_len(x, size),
			    flexsurv::rep_len(mu, size),
			    flexsurv::rep_len(sigma, size),
			    flexsurv::rep_len(Q, size),
			    gengamma::density()),
		     log);
}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector
pgengamma_work(const Rcpp::NumericVector& q,
	       const Rcpp::NumericVector& mu,
	       const Rcpp::NumericVector& sigma,
	       const Rcpp::NumericVector& Q,
	       const bool lower_tail,
	       const bool give_log) {
  if (q.size() == 0) { return q; }
  const R_xlen_t size
    = std::max(std::max(sigma.size(),
			Q.size()),
	       std::max(q.size(),
			mu.size()));

  return mapply(flexsurv::rep_len(q, size),
		flexsurv::rep_len(mu, size),
		flexsurv::rep_len(sigma, size),
		flexsurv::rep_len(Q, size),
		gengamma::cdf(lower_tail, give_log));
}

// [[Rcpp::export(name="check.gengamma", rng=false)]]
Rcpp::LogicalVector check_gengamma(const Rcpp::NumericVector& mu,
				   const Rcpp::NumericVector& sigma,
				   const Rcpp::NumericVector& Q) {
  if ( (0==mu.size()) && (0==sigma.size()) && (0==Q.size()) ) {
    Rcpp::LogicalVector null_result(0);
    return null_result;
  }
  const R_xlen_t size = mu.size();
  return !Rcpp::mapply(mu,
		       flexsurv::rep_len(sigma, size),
		       flexsurv::rep_len(Q, size),
		       gengamma::bad);
}
