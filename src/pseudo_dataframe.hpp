#pragma once

#ifndef PSEUDO_DATAFRAME_HPP
#define PSEUDO_DATAFRAME_HPP

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <string>

typedef std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd>> MatrixXd_Vector;

class PseudoDataFrame {
protected:
  
  std::string treat_name;
  std::string mediator_name;
  std::string interact_name;
  std::size_t num_rows;
  std::size_t treat_i;
  std::size_t mediator_i;
  std::size_t interact_i;
  bool interaction;
  MatrixXd_Vector data_elements;
  std::vector<std::string> element_names;
  std::vector<std::string> terms;
  std::vector<std::string> variables;
  std::vector<std::size_t> element_size;
  std::vector<bool> element_is_matrix;
  
public:
  PseudoDataFrame();
  void initialize_from_env(Rcpp::Environment &env, const std::string& source_variable_name, const std::string& treat_var_name, const std::string& mediator_var_name);
  std::size_t size() const;
  std::size_t rows() const;
  std::size_t cols() const;
  std::size_t column_size(std::size_t i) const;
  const std::string& treat_str() const;
  const std::string& mediator_str() const;
  const std::string& interact_str() const;
  Eigen::MatrixXd& operator[](std::size_t i);
  Eigen::MatrixXd& treat();
  Eigen::MatrixXd& mediator();
  void flatten_into_model_matrix(Eigen::MatrixXd& destmat);
  void to_model_matrix(Eigen::MatrixXd& destmat);
};

#endif //PSEUDO_DATAFRAME_HPP