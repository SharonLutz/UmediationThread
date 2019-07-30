// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif
#include <Rcpp.h>
#include <RcppEigen.h>
#include <exception>
#include <sstream>
#include  "mediation_helper.hpp"

// [[Rcpp::export]]
void mediate_helper(Rcpp::Environment &env){
  // Eigen::setNbThreads(1);
  try{
    MediationHelper h(env);
    // return;
    h();
  }catch(std::exception &e){
    std::stringstream ss;
    ss << "error in mediate_helper, exception caught: " << e.what() << std::endl;
    Rf_error(ss.str().c_str());
  }catch(...){
    Rf_error("error in mediate_helper: Unknown Exception type caught!");
  }
  // Eigen::setNbThreads(0);
}

// [[Rcpp::export]]
void mediate_helper_variable_exporter(Rcpp::Environment &env){
  // Eigen::setNbThreads(1);
  try{
    MediationHelper h(env, true);
    // return;
    h();
  }catch(std::exception &e){
    std::stringstream ss;
    ss << "error in mediate_helper_variable_exporter, exception caught: " << e.what() << std::endl;
    Rf_error(ss.str().c_str());
  }catch(...){
    Rf_error("error in mediate_helper_variable_exporter: Unknown Exception type caught!");
  }
  // Eigen::setNbThreads(0);
}

// [[Rcpp::export]]
void threaded_mediate_helper(Rcpp::Environment &env, long long int num_threads){
  
  // Eigen::setNbThreads(1);
  try{
    if(num_threads > 1){
      Eigen::initParallel();
      MediationHelper h(env, num_threads);
      // return;
      h();
    } else {
      mediate_helper(env);
    }
  }catch(std::exception &e){
    std::stringstream ss;
    ss << "error in threaded_mediate_helper, exception caught: " << e.what() << std::endl;
    Rf_error(ss.str().c_str());
  }catch(...){
    Rf_error("error in threaded_mediate_helper: Unknown Exception type caught!");
  }
  // Eigen::setNbThreads(0);
}
