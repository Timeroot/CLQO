#pragma once

#include <vector>
#include <Eigen/Dense>
#include <utility>		// std::pair, std::get
#include <stdlib.h>     /* abs */
#include <algorithm>    // std::min, std::max
#include <tuple>

using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef unsigned int uint32_t;
typedef std::tuple<float,int,int> clause2;
typedef std::tuple<float,int,int,int> clause3;

//Represents a MAXQP problem
class Problem 
{ 
 public: 
  uint32_t nQP; //total variables in the QP
  
  MatrixXd& coeffs; //strictly uppertriangular
  
  //used to define a "shift" in the objective function. Doesn't
  //affect the search or what the ideal solution is, but is added
  //to the objective value
  float constantTerm;
  
  //Score a proposed solution
  float score(VectorXd sol);
  
  //Initialize a problem with 0 objective function
  Problem(uint32_t n);
  
  //Initialize a problem with given matrix and constant
  Problem(uint32_t n, MatrixXd& coeff, float constantTerm);
  
  //Initialize a problem from a MAX2SAT problem with per-clause and per-var weight
  //Requires one auxiliary variable to represent "true"
  static Problem* from2SAT(uint32_t n, std::vector<clause2>& clauses, std::vector<float>& literalWeights);
  
  //Initialize a problem from a MAX3SAT problem with per-clause and per-var weight
  //Requires one auxiliary variable to represnet "true", and one auxiliary
  //for each clause
  static Problem* from3SAT(uint32_t n, std::vector<clause3>& clause3s, std::vector<clause2>& clause2s, std::vector<float>& literalWeights);
  
  //Initialize a MAXQP problem from a max-clique instance
  static Problem* fromMaxClique(uint32_t n, bool** adjMat);
  
  //Initialize a MAXQP problem from a max-independent set instance
  //This is equivalent to fromMaxClique with a negated adjMat
  static Problem* fromIndSet(uint32_t n, bool** adjMat);
}; 