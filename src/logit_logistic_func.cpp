#include "logit_logistic_func.hpp"

#include <cmath>
#include <ctgmath>

double logit(double p){
  return log( p / (1.0-p));
}

double logit2(double p){
  return log(p) - log(1.0-p);
}

double logit3(double p){
  return -(log((1.0/p)-1.0));
}

double logistic(double a){
  return  1.0 /(1.0 + exp(-a));
}

double logistic2(double a){
  return exp(a) / (exp(a) + 1);
}
