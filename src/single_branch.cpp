//
//  Copyright (c) 2025 Thijs Janzen
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include <cstdlib>    // std::getenv, std::atoi
#include <vector>
#include <chrono>
#include <string>
#include <utility>
#include <algorithm>
#include <memory>
#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "secsse_loglik.h"    // NOLINT [build/include_subdir]

template <typename ODE>
Rcpp::List calc_ll_single_branch(std::unique_ptr<ODE> od,
                                 const Rcpp::NumericVector& states,
                                 const Rcpp::NumericVector& forTime,
                                 const std::string& method,
                                 double atol,
                                 double rtol) {
  try {
    auto t0 = std::min(forTime[0], forTime[1]);
    auto t1 = std::max(forTime[0], forTime[1]);

    auto T0 = std::chrono::high_resolution_clock::now();

    auto states_out = std::vector<double>(states.begin(), states.end());

    auto workhorse = Integrator<ODE, odeintcpp::no_normalization>(
                              std::move(od), method, atol, rtol);

    workhorse(states_out, t0, t1);


    auto T1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> DT = (T1 - T0);


    return Rcpp::List::create(Rcpp::Named("states") = states_out,
                              Rcpp::Named("duration") = DT.count());
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (const char* msg) {
    Rcpp::Rcout << msg << std::endl;
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}

enum class string_code {
  interval2_NE,
  interval2_ES,
  interval2_EC,
  interval3_ES,
  interval3_NE,
  interval4
};

string_code hash_string(const std::string& s) {
  if (s == "interval2_NE") return string_code::interval2_NE;
  if (s == "interval2_ES") return string_code::interval2_ES;
  if (s == "interval2_EC") return string_code::interval2_EC;
  if (s == "interval3_ES") return string_code::interval3_ES;
  if (s == "interval3_NE") return string_code::interval3_NE;
  if (s == "interval4")    return string_code::interval4;
  
  return string_code::interval4;
}

// [[Rcpp::export]]
Rcpp::List cpp_solve(const double& lambda_c,
                     const double& lambda_a,
                     const double& mu,
                     const double& gamma,
                     const std::string& chosen_interval,
                     const std::string& inte_method,
                     const Rcpp::NumericVector& init_states,
                     const Rcpp::NumericVector& time,
                     double atol,
                     double rtol) {
  
  switch( hash_string(chosen_interval)) {
    case string_code::interval2_NE: 
      return calc_ll_single_branch(std::make_unique<loglik::interval2_NE>(lambda_c, lambda_a, mu, gamma), init_states, time, inte_method, atol, rtol);
    case string_code::interval2_ES:
      return calc_ll_single_branch(std::make_unique<loglik::interval2_ES>(lambda_c, lambda_a, mu, gamma), init_states, time, inte_method, atol, rtol);
    case string_code::interval2_EC:
      return calc_ll_single_branch(std::make_unique<loglik::interval2_EC>(lambda_c, lambda_a, mu, gamma), init_states, time, inte_method, atol, rtol);
    case string_code::interval3_ES:
      return calc_ll_single_branch(std::make_unique<loglik::interval3_ES>(lambda_c, lambda_a, mu, gamma), init_states, time, inte_method, atol, rtol);
    case string_code::interval3_NE:
      return calc_ll_single_branch(std::make_unique<loglik::interval3_NE>(lambda_c, lambda_a, mu, gamma), init_states, time, inte_method, atol, rtol);
    case string_code::interval4:
      return calc_ll_single_branch(std::make_unique<loglik::interval4   >(lambda_c, lambda_a, mu, gamma), init_states, time, inte_method, atol, rtol);
  }
  return NA_REAL;
}
