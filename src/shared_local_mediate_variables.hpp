#pragma once

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif //EIGEN_DONT_PARALLELIZE

#ifndef STRUCT_SHAREDLOCALMEDIATEVARIABLES_HPP
#define STRUCT_SHAREDLOCALMEDIATEVARIABLES_HPP

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

#include <string>
#include <array>
#include <vector>
#include "pseudo_dataframe.hpp"

struct SharedLocalMediateVariables {
  int n;
  int sims;
  int cat_0;
  int cat_1;
  std::string& treat();
  std::string& mediator();
  Eigen::Index& treat_i(bool post_flattening=false);
  Eigen::Index& mediator_i(bool post_flattening=false);
  Eigen::Index& interaction_term_i(bool post_flattening=false);
  std::vector<std::string>& terms();
  std::vector<std::string>& variables();
  
  bool isGlm_Y;
  bool FamilyY_is_binomial;
  bool Y_link_is_logit;
  bool interaction_term_present;
  
  std::array< std::array<int, 4>, 4 > tt_switch;
  Eigen::MatrixXd PredictM0;
  Eigen::MatrixXd PredictM1;
  Eigen::MatrixXd YModel;
  PseudoDataFrame y_data;
  
  Eigen::MatrixXd et1;
  Eigen::MatrixXd et2;
  Eigen::MatrixXd et3;
  Eigen::MatrixXd et4;
  SharedLocalMediateVariables();
  
  void initialize_from_environment(Rcpp::Environment & env);
  void store_result_diff(Eigen::MatrixXd &Pr1, Eigen::MatrixXd &Pr0, std::size_t e);
  void export_results(Rcpp::Environment & env);
};
#endif //STRUCT_SHAREDLOCALMEDIATEVARIABLES_HPP