#pragma once

#include "Problem.hpp"
#include <glpk.h> //Linear programming toolkit
#include <Eigen/Sparse>

//Represents a constraint: (a*x[0] + b*x[1] + .. <= rightSide)
typedef struct {
	Eigen::SparseVector<float> coeffs;
	float rightSide;
} constraint;

//Represents a MAXQP solver that uses a linear relaxation, optionally
//with constraint learning and semidefiniteness.
class LPSolver 
{ 
  public: 
	//Problem it's trying to solve
	Problem* problem;
	
	//Current best score achieved, and the solution that does so 
	//TODO: add option for saving last K good solutions (so that they
	//can each be used in later branch/bounds, or local search improvements
	float lowerBound;
	VectorXd bestSol;
	
	//An upper bound on the score of this problem, as determined through
	//the exact score of the relaxation
	float upperBound;
	
	//Constructor: build solver for a given problem
	LPSolver(Problem* p);
	
	//Try to find successively better solutions
	//TODO: some kind of required bound on goodness?
	void solve();

  protected:
	//Convenience copy from Problem
	uint32_t nQP;
	
	//Number of variables in the _linear_ program
	//i.e. nQP-choose-2
	uint32_t nLP;
	
	//GLPK linear problem 
	glp_prob *lp;
	
	//Clauses generated so far (and possibly later removed)
	std::vector<constraint> active_clauses;
	std::vector<constraint> inactive_clauses;
	
	//Last solution to the linear program
	std::vector<float> currSol;
	
  private:
	
	//Given a variable vi and vj, get the index of vij
	uint32_t getLPVar(uint32_t x, uint32_t y);
	
	//given vij, find i and j
	std::pair<uint32_t,uint32_t> getQPVars(uint32_t v);
	
	//Find a minimal submatrix with at least one negative eigenvalue
	std::vector<uint32_t>& nonPSDcore();
	
	//Given an LP solution vector, extract a submatrix. Allocates a new matrix
	//for this purpose, must later be freed
	MatrixXd& getSubmatrix(std::vector<uint32_t> rows);
	
	MatrixXd& getMatrix();
	
	void roundToSol();
};

//TODO findConstraint, main solving loop