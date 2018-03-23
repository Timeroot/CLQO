#pragma once

#include <vector>
#include <Eigen/Dense>
#include <utility>		// std::pair, std::get
#include <stdlib.h>     /* abs */
#include <algorithm>    // std::min, std::max

using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef unsigned int uint32_t;
typedef std::pair<int,int> clause;
typedef std::pair<clause,float> wclause;

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
  
  //Initialize a problem from a 2SAT problem with per-clause and per-var weight
  static Problem* from2SAT(uint32_t n, std::vector<wclause>& clauses, std::vector<float>& literalWeights);
}; 