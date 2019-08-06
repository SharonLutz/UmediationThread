
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

#include "pseudo_dataframe.hpp"
#include <utility>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <iterator>
#include <Rcpp.h>
#include <RcppEigen.h>


PseudoDataFrame::PseudoDataFrame() :
  num_rows(0),
  treat_i(0),
  mediator_i(0),
  interact_i(0),
  interaction(false){}

std::size_t PseudoDataFrame::size() const{
  std::size_t elements_per_row=0;
  for(std::size_t i=0; i < this->cols(); ++i){
    elements_per_row += this->element_size.at(i);
  }
  return this->num_rows * elements_per_row;
}


std::size_t PseudoDataFrame::rows() const{
  return this->num_rows;
}

std::size_t PseudoDataFrame::cols() const{
  return this->element_names.size();
}

std::size_t PseudoDataFrame::column_size(std::size_t i) const{
  return this->element_size.at(i);
}

const std::string& PseudoDataFrame::treat_str() const{
  return this->treat_name;
}

const std::string& PseudoDataFrame::mediator_str() const{
  return this->mediator_name;
}

const std::string& PseudoDataFrame::interact_str() const{
  return this->interact_name;
}

Eigen::MatrixXd& PseudoDataFrame::operator[](std::size_t i){
  if(this->interaction && i == this->interact_i){
    Eigen::MatrixXd& treat_mat = this->data_elements[this->treat_i];
    Eigen::MatrixXd& med_mat = this->data_elements[this->mediator_i];
    this->data_elements.at(i) = treat_mat.array() * med_mat.array();
  } 
  return this->data_elements.at(i);
}

Eigen::MatrixXd& PseudoDataFrame::treat(){
  return this->data_elements.at(this->treat_i);
}

Eigen::MatrixXd& PseudoDataFrame::mediator(){
  return this->data_elements.at(this->mediator_i);
}


void PseudoDataFrame::flatten_into_model_matrix(Eigen::MatrixXd& destmat){
  std::size_t elements_per_row=0;
  for(std::size_t i=0; i < this->cols(); ++i){
    elements_per_row += this->element_size.at(i);
  }
  destmat.resize(this->num_rows,elements_per_row);
  destmat.setZero(this->num_rows,elements_per_row);
  std::size_t next_column=0;
  for(std::size_t i=0;i<this->cols();++i){
    Eigen::MatrixXd &element = this->operator[](i);
    std::size_t element_cols = element.cols();
    for(std::size_t j=0; j < element_cols ;++j){
      destmat.col(next_column) = element.col(j);
      ++next_column;
    }
  }
}

void PseudoDataFrame::to_model_matrix(Eigen::MatrixXd& destmat){
  this->flatten_into_model_matrix(destmat);
  destmat.col(0).setOnes();
}

void get_variable_and_term_names(Rcpp::RObject& y_data_obj, std::vector<std::string>&terms_and_vars, std::vector<std::string>& termvec, std::vector<std::string>& varvec){
  varvec.clear();
  termvec.clear();
  std::unordered_map<std::string, std::size_t> vars;
  std::size_t next_var_index=0;
  Rcpp::RObject terms_obj = y_data_obj.attr("terms");
  Rcpp::IntegerMatrix factor_matrix = terms_obj.attr("factors");
  Rcpp::GenericVector factor_dimnames = factor_matrix.attr("dimnames");
  Rcpp::CharacterVector variables_vec = factor_dimnames[0];
  Rcpp::CharacterVector terms_vec = factor_dimnames[1];
  varvec.reserve(variables_vec.size());
  termvec.reserve(terms_vec.size());
  for(auto& element : variables_vec){
    auto el_str = std::string(element);
    if(!vars.count(el_str)){
      vars[el_str] = next_var_index;
      varvec.push_back(el_str);
      ++next_var_index;
    }
  }
  for(auto& element : terms_vec){
    auto el_str = std::string(element);
    if(!vars.count(el_str)){
      vars[el_str] = next_var_index;
      termvec.push_back(el_str);
      ++next_var_index;
    }
  }
  terms_and_vars.clear();
  terms_and_vars.resize(vars.size(), "");
  for(auto& element : vars){
    terms_and_vars.at(element.second) = element.first;
  }
}

std::vector<std::pair<std::string, bool>> check_element_classes(Rcpp::RObject& y_data_obj){
  std::vector<std::pair<std::string, bool>> result;
  Rcpp::RObject terms_obj = y_data_obj.attr("terms");
  Rcpp::CharacterVector data_classes = terms_obj.attr("dataClasses");
  Rcpp::CharacterVector data_names = data_classes.names();
  for(long long int i =0 ; i < data_classes.size(); ++i){
    auto class_str = std::string(data_classes[i]);
    auto name_str = std::string(data_names[i]);
    bool is_a_matrix = class_str.find("matrix") != std::string::npos;
    result.push_back(std::make_pair(name_str, is_a_matrix));
  }
  return result;
}

void PseudoDataFrame::initialize_from_env(Rcpp::Environment &env, const std::string& source_variable_name, const std::string& treat_var_name, const std::string& mediator_var_name){
  this->treat_name= treat_var_name;
  this->mediator_name = mediator_var_name;
  
  std::string interact_str_1;
  std::string interact_str_2;
  
  {
    std::stringstream ss1,ss2;
    ss1 << this->treat_name;
    ss2 << this->mediator_name;
    ss1 << ':';
    ss2 << ':';
    ss1 << this->mediator_name;
    ss2 << this->treat_name;
    interact_str_1 = ss1.str();
    interact_str_2 = ss2.str();
  }
  {
    std::stringstream ss;
    ss << this->mediator_name << ":" << this->treat_name;
    this-> interact_name = ss.str();
  }
  
  bool found_treat=false;
  bool found_mediator=false;
  
  Rcpp::RObject R_y_data = env[source_variable_name];
  
  get_variable_and_term_names(R_y_data, this->element_names, this->terms, this->variables);
  Rcpp::CharacterVector row_names = R_y_data.attr("row.names");
  num_rows = row_names.length();
  auto matrix_checks = check_element_classes(R_y_data);
  this->element_is_matrix.resize(this->element_names.size(), false);
  this->element_size.resize(this->element_names.size(),1);
  for(std::size_t i = 0 ; i < matrix_checks.size(); ++i){
    if(this->element_names.at(i) != matrix_checks.at(i).first){
      std::stringstream ss;
      ss << "element name at position in data class list doesn't match variable name at position in element names: "<<i;
      Rf_error(ss.str().c_str());
    }
    this->element_is_matrix.at(i) = matrix_checks.at(i).second;
  }
  
  Rcpp::DataFrame R_y_data_df = R_y_data;
  
  for(long long int i =0; i < R_y_data_df.cols(); ++i){
    if(this->element_is_matrix.at(i)){
      //Rcpp::Rcout << "data element "<< i << " appears to be a matrix" << std::endl;
      Rcpp::NumericMatrix data = R_y_data_df[i];
      this->data_elements.emplace_back(Rcpp::as<Eigen::MatrixXd>(data));
      this->element_size.at(i) = data.cols();
    } else {
      //Rcpp::Rcout << "data element "<< i << " appears to be a vector" << std::endl;
      Rcpp::NumericVector data = R_y_data_df[i];
      this->data_elements.emplace_back(data.size(), 1);
      this->data_elements.back() = Rcpp::as<Eigen::VectorXd>(data).matrix();
      this->element_size.at(i) = 1;
    }
  }
  
  auto treat_iter = std::find(this->element_names.begin(), this->element_names.end(), this->treat_name);
  auto med_iter = std::find(this->element_names.begin(), this->element_names.end(), this->mediator_name);
  found_treat = treat_iter != this->element_names.end();
  found_mediator = med_iter != this->element_names.end();
  
  if(!(found_treat && found_mediator)){
    std::stringstream ss;
    ss << "error missing model variable(s):";
    if(!found_treat){
      ss << " treatment variable: " << this->treat_name << ';';
    }
    if(!found_mediator){
      ss << " mediator variable: " << this->mediator_name << ';';
    }
    Rf_error(ss.str().c_str());
  }
  
  this->treat_i = std::distance(this->element_names.begin(), treat_iter);
  this->mediator_i = std::distance(this->element_names.begin(), med_iter);
  
  auto interact1_iter = std::find(this->element_names.begin(), this->element_names.end(), interact_str_1);
  auto interact2_iter = std::find(this->element_names.begin(), this->element_names.end(), interact_str_2);
  bool interact_str_1_found = interact1_iter != this->element_names.end();
  bool interact_str_2_found = interact2_iter != this->element_names.end();
  if(interact_str_1_found ||interact_str_2_found){
    this->interaction = true;
    if(interact_str_1_found && interact_str_2_found){
      std::stringstream ss;
      ss << "error: both interaction strings found: "<< interact_str_1 << " , " << interact_str_2;
      Rf_error(ss.str().c_str());    
    }
    if(interact_str_1_found){
      this->interact_i = std::distance(this->element_names.begin(), interact1_iter);
    }
    if(interact_str_2_found){
      this->interact_i = std::distance(this->element_names.begin(), interact2_iter);
    }
  }
  
  if(this->interaction){
    Eigen::MatrixXd& treat_mat = this->data_elements[this->treat_i];
    Eigen::MatrixXd& med_mat = this->data_elements[this->mediator_i];
    if((treat_mat.rows() != med_mat.rows()) || (treat_mat.cols() != med_mat.cols())){
      std::stringstream ss;
      ss <<"dimension mismatch between treatment and mediator data elements";
      if(treat_mat.rows() != med_mat.rows()){
        ss <<"; ROWS: treat: "<<treat_mat.rows() << " != mediator: "<<med_mat.rows();
      }
      if(treat_mat.cols() != med_mat.cols()){
        ss <<"; COLS: treat: "<<treat_mat.cols() << " != mediator: "<<med_mat.cols();
      }
      this->element_is_matrix.at(this->interact_i) = this->element_is_matrix.at(this->treat_i) || this->element_is_matrix.at(this->mediator_i);
    }
    this->data_elements.emplace_back(treat_mat.rows(), treat_mat.cols());
  }
}