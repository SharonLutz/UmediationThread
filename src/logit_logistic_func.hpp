#pragma once

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif //EIGEN_DONT_PARALLELIZE

#ifndef FUNC_LOGIT_LOGISTIC_HPP
#define FUNC_LOGIT_LOGISTIC_HPP

double logit(double p);
double logit2(double p);
double logit3(double p);

double logistic(double a);
double logistic2(double a);

#endif //FUNC_LOGIT_LOGISTIC_HPP