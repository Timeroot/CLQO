#include "Problem.hpp"
#include "LPSolver.hpp"

#include <iostream>
#include <random>
#include <chrono>
#include <cmath>

void setSampleProblem(uint32_t& variables, std::vector<clause3>& clause3s, std::vector<clause2>& clauses, std::vector<float>& literalWeights);

int main(){
	uint32_t variables;
	 //Each term is <(v1,v2), weight>, representing v1 OR v2.
	std::vector<clause3> clause3s;
	std::vector<clause2> clause2s;
	std::vector<float> literalWeights;
	setSampleProblem(variables, clause3s, clause2s, literalWeights);
	
	Problem* p = Problem::from3SAT(variables, clause3s, clause2s, literalWeights);
	LPSolver solver = LPSolver(p);
	solver.solve();
}

void setSampleProblem(uint32_t& variables, std::vector<clause3>& clause3s, std::vector<clause2>& clause2s, std::vector<float>& literalWeights){
	variables = 48;
	literalWeights.resize(variables);
	
	std::vector<int> realSol;
	std::default_random_engine generator(1234+std::chrono::high_resolution_clock::now().time_since_epoch().count()*0);
	std::uniform_int_distribution<int> truthGenerator(0,1);
	std::uniform_int_distribution<int> variableChoser(1,variables);
	std::cout << "True satisfaction: " << std::endl;
	for(uint32_t i=0;i<variables;i++){
		realSol.push_back(2*truthGenerator(generator) - 1);
		std::cout << realSol[i] << ", ";
	}
	std::cout  << std::endl;
	
	uint32_t bothTrueClause = 1*variables;
	uint32_t OneOneClause = 5*variables;
	uint32_t bothFalseClause = (variables <= 5 ? variables - 1 : 1*variables - 5);
	bool w1 = true;
	
	std::cout << "Should get at least " << (bothTrueClause + OneOneClause) << std::endl; 
	for(uint32_t i=0;i<bothTrueClause;i++){
		int varA = variableChoser(generator), varB;
		do varB = variableChoser(generator); while (varB == varA);
		varA *= realSol[varA-1];
		varB *= realSol[varB-1];
		//std::cout << "Clause " << varA << ", " << varB << std::endl;
		float weight = w1 ? 1 : (variables+variableChoser(generator))*variableChoser(generator)*0.1;
		clause2s.push_back({weight, varA, varB});
	}
	for(uint32_t i=0;i<OneOneClause;i++){
		int varA = variableChoser(generator), varB;
		do varB = variableChoser(generator); while (varB == varA);
		varA *= realSol[varA-1];
		varB *= -realSol[varB-1];
		//std::cout << "Clause " << varA << ", " << varB << std::endl;
		float weight = w1 ? 1 : (variables+variableChoser(generator))*variableChoser(generator)*0.1;
		clause2s.push_back({weight, varA, varB});
	}
	for(uint32_t i=0;i<bothFalseClause;i++){
		int varA = variableChoser(generator), varB;
		do varB = variableChoser(generator); while (varB == varA);
		varA *= -realSol[varA-1];
		varB *= -realSol[varB-1];
		//std::cout << "Clause " << varA << ", " << varB << std::endl;
		float weight = w1 ? 1 : (variables+variableChoser(generator))*variableChoser(generator)*0.1;
		clause2s.push_back({weight, varA, varB});
	}
	
	/*clause3s.push_back({1, 1, 2, 3});
	clause3s.push_back({1.5, -1, -4, 3});
	clause3s.push_back({1.6, -2, -3, 4});*/
}