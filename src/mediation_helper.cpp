// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif //EIGEN_DONT_PARALLELIZE
#include "mediation_helper.hpp"
#include "logit_logistic_func.hpp"
#include "pseudo_dataframe.hpp"
//#define MEDIATION_DEBUG_PRINTING

OuterLoopVars::OuterLoopVars(SharedLocalMediateVariables &sv, std::size_t e) :
  e(e),
  tt(sv.tt_switch.at(e)),
  Pr1(sv.n, sv.sims),
  Pr0(sv.n, sv.sims),
  cat_t(tt[0] ? sv.cat_1 : sv.cat_0 ),
  cat_t_ctrl(tt[0] ? sv.cat_0 : sv.cat_1 ),
  cat_c(tt[1] ? sv.cat_1 : sv.cat_0 ),
  cat_c_ctrl(tt[1] ? sv.cat_0 : sv.cat_1 ){}

MediationHelper::MediationHelper(Rcpp::Environment &env, long long int num_threads) : 
  num_threads(num_threads), 
  env(env),
  export_loop_vars(false),
  vars_exported(false),
  sv(), 
  threads(),
  promises(){
  sv.initialize_from_environment(env);
}

MediationHelper::MediationHelper(Rcpp::Environment &env, bool export_loop_vars) : 
  num_threads(1),
  env(env),
  export_loop_vars(export_loop_vars),
  vars_exported(false),
  sv(), 
  threads(),
  promises(){
  sv.initialize_from_environment(env);
}

void MediationHelper::inner_loop(OuterLoopVars &olv, std::size_t j){
#ifdef MEDIATION_DEBUG_PRINTING
  Rcpp::Rcout << "inner: " << j << " begin"<< std::endl;
#endif //MEDIATION_DEBUG_PRINTING
  PseudoDataFrame pred_data_t(sv.y_data), pred_data_c(sv.y_data);
  
#ifdef MEDIATION_DEBUG_PRINTING
  Rcpp::Rcout << "pred_data vars initialized with y_data"<< std::endl;
#endif //MEDIATION_DEBUG_PRINTING
  
  pred_data_t.treat().col(0).setConstant(olv.cat_t);
  pred_data_c.treat().col(0).setConstant(olv.cat_c);
  
#ifdef MEDIATION_DEBUG_PRINTING
  Rcpp::Rcout << "set category constants for pred.data matrices treatment cols"<< std::endl;
#endif //MEDIATION_DEBUG_PRINTING
  
  //PredictMt <- PredictM1[j,] * tt[3] + PredictM0[j,] * (1 - tt[3])
  //PredictMc <- PredictM1[j,] * tt[4] + PredictM0[j,] * (1 - tt[4])
  auto PredictMt = (sv.PredictM1.row(j) * olv.tt[2]) + (sv.PredictM0.row(j) * (1-olv.tt[2]));
  auto PredictMc = (sv.PredictM1.row(j) * olv.tt[3]) + (sv.PredictM0.row(j) * (1-olv.tt[3]));
  
#ifdef MEDIATION_DEBUG_PRINTING
  Rcpp::Rcout << "calculated PredictM_ values"<< std::endl;
#endif //MEDIATION_DEBUG_PRINTING
  //pred.data.t[,mediator] <- PredictMt
  //pred.data.c[,mediator] <- PredictMc
  
  pred_data_t.mediator().col(0) = PredictMt;
  pred_data_c.mediator().col(0) = PredictMc;
  
#ifdef MEDIATION_DEBUG_PRINTING
  Rcpp::Rcout << "stored predictM_ in pred.data vars"<< std::endl;
#endif //MEDIATION_DEBUG_PRINTING
  Eigen::MatrixXd ymat_t;
  Eigen::MatrixXd ymat_c;
  
  pred_data_t.to_model_matrix(ymat_t);
  pred_data_c.to_model_matrix(ymat_c);
  
#ifdef MEDIATION_DEBUG_PRINTING
  Rcpp::Rcout << "tranformed predictions -> model matrix"<< std::endl;
#endif //MEDIATION_DEBUG_PRINTING
  if(export_loop_vars && (!vars_exported)){
    Rcpp::Rcout << "attempting to export inner loop vars"<< std::endl;
    auto genv = Rcpp::Environment::global_env();
    // genv[["vanillaR.ymat.t"]] = ymat.t
    genv["rcpp.ymat.t"] = ymat_t;
    // genv[["vanillaR.ymat.c"]] = ymat.c
    genv["rcpp.ymat.c"] = ymat_c;
    // genv[["vanillaR.pred.data.t"]] = pred.data.t
    Eigen::MatrixXd tmp_t;
    Eigen::MatrixXd tmp_c;
    pred_data_t.flatten_into_model_matrix(tmp_t);
    pred_data_c.flatten_into_model_matrix(tmp_c);
    genv["rcpp.pred.data.t"] = tmp_t;
    // genv[["vanillaR.pred.data.c"]] = pred.data.c
    genv["rcpp.pred.data.c"] = tmp_c;
    // genv[["vanillaR.tt"]] = tt
    genv["rcpp.tt"] = olv.tt;
    // genv[["vanillaR.PredictMt"]] = PredictMt
    genv["rcpp.PredictMt"] = PredictMt;
    // genv[["vanillaR.PredictMc"]] = PredictMc
    genv["rcpp.PredictMc"] = PredictMc;
    // genv[["vanillaR.cat.t"]] = cat.t
    genv["rcpp.cat.t"] = olv.cat_t;
    // genv[["vanillaR.cat.c"]] = cat.c
    genv["rcpp.cat.c"] = olv.cat_c;
    // genv[["vanillaR.cat.t.ctrl"]] = cat.t.ctrl
    genv["rcpp.cat.t.ctrl"] = olv.cat_t_ctrl;
    // genv[["vanillaR.cat.c.ctrl"]] = cat.c.ctrl
    genv["rcpp.cat.c.ctrl"] = olv.cat_c_ctrl;
    vars_exported = true;
  }
  
  //Pr1[,j] <- t(as.matrix(YModel[j,])) %*% t(ymat.t)
  //Pr0[,j] <- t(as.matrix(YModel[j,])) %*% t(ymat.c)
  olv.Pr1.col(j) = (sv.YModel.row(j) * ymat_t.transpose()).row(0);
  olv.Pr0.col(j) = (sv.YModel.row(j) * ymat_c.transpose()).row(0);
#ifdef MEDIATION_DEBUG_PRINTING
  Rcpp::Rcout << "t(as.matrix(YModel[j,])) %*% t(ymat) calculated"<< std::endl;
  Rcpp::Rcout << "inner: " << j << " end"<< std::endl;
#endif //MEDIATION_DEBUG_PRINTING
}

void MediationHelper::outer_loop(std::size_t e) {
#ifdef MEDIATION_DEBUG_PRINTING
  Rcpp::Rcout<< "outer_loop begin: " << e << std::endl;
#endif //MEDIATION_DEBUG_PRINTING
  OuterLoopVars olv(sv,e);
#ifdef MEDIATION_DEBUG_PRINTING
  Rcpp::Rcout<< "outer_loop vars initialized: " << e << std::endl;
#endif //MEDIATION_DEBUG_PRINTING
  for(int j=0; j < sv.sims; ++j){
    inner_loop(olv,j);
  }
  
  // check for GLM, binomial, logit in model.y
  if(sv.FamilyY_is_binomial){
    
    if(export_loop_vars){
      auto genv = Rcpp::Environment::global_env();
      Rcpp::Rcout << "attempting to export outer loop vars"<< std::endl;
      std::stringstream ss1,ss2;
      ss1 << "rcpp.Pr1.pre_apply." << (e+1);
      ss2 << "rcpp.Pr0.pre_apply." << (e+1);
      genv[ss1.str()] = olv.Pr1;
      genv[ss2.str()] = olv.Pr0;
    }
    
    // TODO: if so, apply logistic(a) to Pr0 and Pr1
    olv.Pr1 = olv.Pr1.unaryExpr(std::ref(logistic));
    olv.Pr0 = olv.Pr0.unaryExpr(std::ref(logistic));
  }
  
  sv.store_result_diff(olv.Pr1, olv.Pr0, e);
  
#ifdef MEDIATION_DEBUG_PRINTING
  Rcpp::Rcout<< "outer_loop end: " << e << std::endl;
#endif //MEDIATION_DEBUG_PRINTING
  if(export_loop_vars){
    auto genv = Rcpp::Environment::global_env();
    Rcpp::Rcout << "attempting to export outer loop vars"<< std::endl;
    std::stringstream ss1,ss2;
    ss1 << "rcpp.Pr1." << (e+1);
    ss2 << "rcpp.Pr0." << (e+1);
    genv[ss1.str()] = olv.Pr1;
    genv[ss2.str()] = olv.Pr0;
  }
}

class protected_work_supplier{
protected:
  std::mutex m;
  std::size_t i;
  std::size_t limit;
public:
  protected_work_supplier(std::size_t limit) : i(0), limit(limit){}
  std::size_t operator()(){
    if(i == limit){
      return limit;
    } else {
      std::lock_guard<std::mutex> g(m);
      return i++;
    }
  }
  bool is_done(){
    std::lock_guard<std::mutex> g(m);
    return i==limit;
  }
};

void MediationHelper::outer_loop_with_threaded_inner_loop(std::size_t e){
  OuterLoopVars local_vars(sv,e);
  
  std::size_t num_sims = sv.sims;
  protected_work_supplier next_index_value(num_sims);
  
  auto inner_loop_thread = [&](std::size_t thread_index){
    try{
      std::size_t j;
      while(!next_index_value.is_done()){
        j = next_index_value();
        if(j >= num_sims){
          break;
        }
        this->inner_loop(local_vars, j);
      }
    } catch(...){
      // TODO, signal the error to the main thread somehow
    }
    promises[thread_index].set_value();
  };
  
  for(std::size_t i=0; i < num_threads; ++i){
    threads.emplace_back(inner_loop_thread, i);
  }
  
  // std::this_thread::sleep_for(std::chrono::milliseconds(10));
  
  for(std::size_t i=0;i<num_threads;++i){
    auto f = promises[i].get_future();
    f.wait();
  }
  
  for(std::size_t i=0;i<num_threads;++i){
    threads[i].join();
  }
  
  //TODO: if isGlm.y, apply inv.logit (logistic) to Pr1 and Pr0
  // sv.effects_tmp[e] = (Pr1 - Pr0);
  if(sv.FamilyY_is_binomial){
    local_vars.Pr1 = local_vars.Pr1.unaryExpr(std::ref(logistic));
    local_vars.Pr0 = local_vars.Pr0.unaryExpr(std::ref(logistic));
  }
  
  sv.store_result_diff(local_vars.Pr1, local_vars.Pr0, e);
}


void MediationHelper::launcher_n_2(){
  threads.clear();
  promises.clear();
  promises.resize(4);
  auto outer_loop_thread = [&](std::size_t e){
    try{
    this->outer_loop(e);
    this->promises[e].set_value();
    } catch (...){
      // TODO signal error to main thread somehow.
    }
  };
  //2 batches of 2, runs outer loop, auto-launching the inner loop
  for(std::size_t batch=0; batch < 2; ++batch){
    for(std::size_t i=0;i<2;++i){
      threads.emplace_back(outer_loop_thread, i + (batch*2));
    }
    // std::this_thread::sleep_for(std::chrono::milliseconds(10));
    for(std::size_t i=0;i<2;++i){
      auto f = promises[i + (batch*2)].get_future();
      f.wait();
    }
    for(std::size_t i=0;i<2;++i){
      threads[i + (batch*2)].join();
    }
  }
}

void MediationHelper::launcher_n_4(){
  threads.clear();
  promises.clear();
  promises.resize(4);
  auto outer_loop_thread = [&](std::size_t e){
    try {
    this->outer_loop(e);
    this->promises[e].set_value();
    } catch (...){
      // TODO signal error to main thread somehow
    }
  };
  
  for(std::size_t e=0;e<4;++e){
    threads.emplace_back(outer_loop_thread, e);
  }
  
  // std::this_thread::sleep_for(std::chrono::milliseconds(10));
  
  for(std::size_t e=0;e<4;++e){
    auto f = promises[e].get_future();
    f.wait();
  }
  for(std::size_t e=0;e<4;++e){
    threads[e].join();
  }
}

void MediationHelper::launcher_n_other(){
  threads.clear();
  promises.clear();
  for(std::size_t e=0;e<4;++e){
    promises.resize(num_threads);
    outer_loop_with_threaded_inner_loop(e);
    threads.clear();
    promises.clear();
  }
}

void MediationHelper::launcher_n_1(){
  threads.clear();
  promises.clear();
  for(std::size_t e=0;e<4;++e){
    outer_loop(e);
  }
}

void MediationHelper::operator()(){
  // Rcpp::Rcout << "running mediation with "<<num_threads<< " worker threads" << std::endl;
  switch(num_threads){
  case 1:
    launcher_n_1();
    break;
  case 2:
    launcher_n_2();
    break;
  case 4:
    launcher_n_4();
    break;
  default:
    launcher_n_other();
  break;
  }
  sv.export_results(env);
}