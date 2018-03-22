#include <glpk.h>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> //for PSD tests
#include <Eigen/Cholesky> //for PSD tests
#include <utility>		// std::pair, std::get
#include <stdlib.h>     /* abs */
#include <algorithm>    // std::min, std::max
#include <stdio.h>
#include <stdexcept>
#include <string>
#include <iostream>
#include <random>
#include <chrono>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::Vector4d;

typedef unsigned int uint32_t;
typedef std::pair<int,int> clause;
typedef std::pair<clause,float> wclause;

#define PSD_EIGEN_TOL 0.00001
#define CONSTRAINT_FAIL_LIMIT 5
#define MAX_TRIES_ROUNDING 5

void setSampleProblem(uint32_t& variables, std::vector<wclause>& clauses, std::vector<float>& literalWeights){
	variables = 40;
	literalWeights.resize(variables);
	
	/*literalWeights[0] += 0.1;
	literalWeights[1] += 0.53;
	literalWeights[4] += 1.2;
	literalWeights[6] -= 3.7;
	literalWeights[7] += 0.53;
	literalWeights[9] += 1.2;
	literalWeights[11] -= 3.7;
	
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-1,-2), 1.3));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-3,1), 1.7));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-2,3), 0.6));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-4,5), 0.9));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-1,7), 0.4));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-7,2), 0.1));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-5,7), 0.4));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-6,5), 0.2));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-7,4), 0.3));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-4,3), 1.4));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-3,6), 0.8));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(4,6), 1.3));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-6,2), 1.8));
	
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-8,2), 0.1));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-5,8), 0.4));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-9,5), 0.2));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-7,9), 0.3));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-10,3), 1.4));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-3,10), 0.8));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-11,6), 1.3));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-5,11), 1.8));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-12,4), 1.7));
	clauses.push_back(std::pair<clause,float>(std::pair<int,int>(-5,12), 1.6));*/
	
	std::vector<int> realSol;
	std::default_random_engine generator(std::chrono::high_resolution_clock::now().time_since_epoch().count());
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
	uint32_t bothFalseClause = 1*variables - 5;
	for(uint32_t i=0;i<bothTrueClause;i++){
		int varA = variableChoser(generator), varB;
		do varB = variableChoser(generator); while (varB == varA);
		varA *= realSol[varA-1];
		varB *= realSol[varB-1];
		clauses.push_back(std::pair<clause,float>(std::pair<int,int>(varA,varB), 1.));
	}
	for(uint32_t i=0;i<OneOneClause;i++){
		int varA = variableChoser(generator), varB;
		do varB = variableChoser(generator); while (varB == varA);
		varA *= realSol[varA-1];
		varB *= -realSol[varB-1];
		clauses.push_back(std::pair<clause,float>(std::pair<int,int>(varA,varB), 1.));
	}
	for(uint32_t i=0;i<bothFalseClause;i++){
		int varA = variableChoser(generator), varB;
		do varB = variableChoser(generator); while (varB == varA);
		varA *= -realSol[varA-1];
		varB *= -realSol[varB-1];
		clauses.push_back(std::pair<clause,float>(std::pair<int,int>(varA,varB), 1.));
	}
}

unsigned long isqrt(unsigned long x)
{
	register unsigned long op, res, one;  
	
	op = x;
	res = 0;
	
	/* "one" starts at the highest power of four <= than the argument. */
	one = 1 << 30;  /* second-to-top bit set */
	while (one > op) one >>= 2;
	
	while (one != 0) {
		if (op >= res + one) {
			op -= res + one;
			res += one << 1;  // <-- faster than 2 * one
		}
		res >>= 1;
		one >>= 2;
	}
	return res;
}

//Given a variable vi and vj, get the index of vij
int getLPVar(int x, int y){
	if(y == x) throw std::runtime_error(std::string("Bad LP: ")+std::to_string(x)+", "+std::to_string(y));
	if(y > x) return getLPVar(y, x);
	return 1 + y + x*(x-1)/2;
}

//given vij, find i and j
std::pair<int,int> getQPVars(int v, int nLP){
	if(v > nLP) throw std::runtime_error(std::string("Bad vQP: ")+std::to_string(v)+", "+std::to_string(nLP));
	
	int x = floor(sqrt(2*v)+1./2);
	int y = v - 1 - x*(x-1)/2;
	
	if(y >= x) throw std::runtime_error(std::string("Bad QP: ")+std::to_string(x)+", "+std::to_string(y)+", "+std::to_string(v));
	
	return std::pair<int,int>(x,y);
}

MatrixXd& getSubmatrix(std::vector<int> rows, std::vector<float> currSol){
	
	MatrixXd& result = *new MatrixXd(rows.size(), rows.size());
	for(uint32_t i=0;i<rows.size();i++){
		result(i,i) = 1;
		
		for(uint32_t j=i+1;j<rows.size();j++)
			result(j,i) = result(i,j) = currSol[getLPVar(rows[i],rows[j])-1];
	}
	
	return result;
}

MatrixXd& getMatrix(std::vector<float> currSol, uint32_t nQP){
	MatrixXd& result = *new MatrixXd(nQP, nQP);
	for(uint32_t i=0;i<nQP;i++){
		result(i,i) = 1;
		
		for(uint32_t j=i+1;j<nQP;j++)
			result(j,i) = result(i,j) = currSol[getLPVar(i,j)-1];
	}
	
	return result;
}

//Given an assignment, find a minimal (but non necess. minimum) set of
//rows/columns that has a constraint that's it violating.
//If it's actually positive semidefinite, and there are not constraints
//that it's violating, return the empty list. (Congrats! This means you
//found the global optimum!) 
std::vector<int>& nonPSDcore(uint32_t nQP, uint32_t nLP, std::vector<float> currSol){
	//Two phases: one where we add rows hoping to make it not PSD,
	//and then later we remove rows to make it minimal.
	
	//Rows we're currently taking
	std::vector<int>& core = *new std::vector<int>();
	//Rows we have not taken
	std::vector<int> notInCore;
	for(uint32_t i=0;i<nQP;i++) notInCore.push_back(i);
	
	//put the rows in a random order
	std::random_shuffle(notInCore.begin(), notInCore.end());
	
	do{
		int newRow = notInCore.back();
		notInCore.pop_back();
		core.push_back(newRow);
		
		MatrixXd& subMat = getSubmatrix(core, currSol);
		//std::cout << subMat << std::endl;
		
		Eigen::SelfAdjointEigenSolver<MatrixXd> eig(subMat);
		auto evals = eig.eigenvalues();
		bool isPSD = true;
		for(int i=0; i<evals.size(); i++){
			if(evals[i] < -PSD_EIGEN_TOL){
				isPSD = false;
				break;
			}
		}
		delete &subMat;
		
		if(isPSD){
			if(notInCore.size() > 0)
				continue; //add more
			else {
				core.clear();
				return core; //nothing nonPSD, return empty vector
			}
		} else {
			//we found a non-PSD chunk!
			break;
		}
		
	}while(true);
	
	/*printf("nonPSD core found in vars \n");
	for(uint32_t i=0;i<core.size();i++){;
		printf("%d, ", core[i]);
	}
	printf("\n");*/
	
	int largeMatSize = core.size();
	for(int iter=0; iter<largeMatSize; iter++){
		int removedRow = core[0];
		//std::cout << "Pared down " << removedRow << std::endl;
		core.erase(core.begin());
		
		MatrixXd& subMat = getSubmatrix(core, currSol);
		//std::cout << subMat << std::endl;
		
		Eigen::SelfAdjointEigenSolver<MatrixXd> eig(subMat);
		auto evals = eig.eigenvalues();
		bool isPSD = true;
		for(int i=0; i<evals.size(); i++){
			if(evals[i] < -PSD_EIGEN_TOL){
				isPSD = false;
				break;
			}
		}
		delete &subMat;
		
		if(isPSD){
			//the row we removed was necessary for non-PSD-ness, add it back in
			core.push_back(removedRow);
			//std::cout << "kept." << std::endl;
		} else {
			//great, we didn't need it. do nothing and repeat.
		}
	}
	
	/*printf("nonPSD core found in vars \n");
	for(uint32_t i=0;i<core.size();i++){;
		printf("%d, ", core[i]);
	}
	printf("\n");*/
	return core;
}

void roundMat(glp_prob *lp, std::vector<float> currSol, MatrixXd& Amat, int nQP, int nLP, float constantTerm){
	
	//TODO call to an actual SDP solver
	
	MatrixXd& solMat = getMatrix(currSol, nQP);
	
	std::cout << "The matrix Sol:" << std::endl << solMat << std::endl;
	
	//for now --
	//Reduce 'magnitude' slightly to make PSD
	Eigen::SelfAdjointEigenSolver<MatrixXd> eig(solMat);
	auto evals = eig.eigenvalues();
	double minEigenvalue = evals[0] - 0.001;
	
	solMat *= -1. / (-1 + minEigenvalue);
	solMat += MatrixXd::Identity(nQP, nQP) * (minEigenvalue / (-1 + minEigenvalue));
	
	std::cout << "SDP-ified Sol:" << std::endl;
	
	std::cout << solMat << std::endl;
	Eigen::LLT<MatrixXd> lltOfA(solMat); // compute the Cholesky decomposition of A
	MatrixXd L = lltOfA.matrixL(); // retrieve factor L  in the decomposition
	// The previous two lines can also be written as "L = A.llt().matrixL()"
	std::cout << "The Cholesky factor L is" << std::endl << L << std::endl;
	
	delete &solMat;
	
	for(int tries=0; tries<MAX_TRIES_ROUNDING; tries++){
		//generate random dot vector
		VectorXd v(nQP);
		std::default_random_engine generator(std::chrono::high_resolution_clock::now().time_since_epoch().count());
		std::normal_distribution<double> distribution(0.0,1.0);
		for(int i=0;i<nQP;i++){
			v(i) = distribution(generator);
			if(i==0) v(i) = fabs(v(i)); //normalize solution to start with "1"
		}
		
		VectorXd resultVec = L * v;
		//std::cout << "Dot product vector" << resultVec.transpose() << std::endl;
		for(int i=0;i<nQP;i++)
			resultVec(i) = std::copysign(1.0,resultVec(i));
		
		/*std::cout << "Final assignment: " << std::endl;
		for(int i=0;i<nQP;i++)
			printf("%f, ", resultVec(i));
		puts("");*/
		
		printf("Score = %f\n", resultVec.transpose() * Amat * resultVec + constantTerm);
	}
}

//Returns a constraint that the submatrix violates, or 0-length
//if none is found.
std::vector<double>& findConstraint(MatrixXd& subMat){
	if(subMat.rows() == 1){
		Vector3d bestVec; 
		bestVec(0) = bestVec(1) = bestVec(2) = 1;
		float bestDot = 2; //maximum 'useful' score of 1
		for(int bits=0; bits<4; bits++){
			Vector3d test;
			test(0) = ((bits>>0) & 1)*2 - 1;
			test(1) = ((bits>>1) & 1)*2 - 1;
			test(2) = 1;
			float score = test.transpose() * subMat * test;
			if(score < bestDot){
				bestDot = score;
				bestVec = test;
			}
		}
		std::vector<double>& res = *new std::vector<double>(4);
		res[1] = bestVec(0)*bestVec(1);
		res[2] = bestVec(0)*bestVec(2);
		res[3] = bestVec(1)*bestVec(2);
		res[0] = -1;
		return res;
	}
	if(subMat.rows() == 1){
		Vector4d bestVec; 
		bestVec(0) = bestVec(1) = bestVec(2) = bestVec(3) = 1;
		float bestDot = 2; //maximum 'useful' score of 1
		for(int dir=0; dir<4; dir++){
				for(int bits=0; bits<4; bits++){
				Vector4d test;
				test(0) = ((bits>>0) & 1)*2 - 1;
				test(1) = ((bits>>1) & 1)*2 - 1;
				test(2) = 1;
				
				float temp = test(dir);
				test(dir) = 0;
				test(3) = temp;
				
				float score = test.transpose() * subMat * test;
				if(score < bestDot){
					bestDot = score;
					bestVec = test;
				}
			}
		}
		std::vector<double>& res = *new std::vector<double>(7);
		res[1] = bestVec(1)*bestVec(0);
		res[2] = bestVec(2)*bestVec(0);
		res[3] = bestVec(2)*bestVec(1);
		res[4] = bestVec(3)*bestVec(0);
		res[5] = bestVec(3)*bestVec(1);
		res[6] = bestVec(3)*bestVec(2);
		res[0] = -1;
		return res;
	}
	if(subMat.rows() == 1){
		VectorXd bestVec(5); 
		bestVec(0) = bestVec(1) = bestVec(2) = bestVec(3) = bestVec(4) = 1;
		float bestDot = 100; //maximum 'useful' score of 2
		for(int bits=0; bits<16; bits++){
			VectorXd test(5);
			test(0) = ((bits>>0) & 1)*2 - 1;
			test(1) = ((bits>>1) & 1)*2 - 1;
			test(2) = ((bits>>2) & 1)*2 - 1;
			test(3) = ((bits>>3) & 1)*2 - 1;
			test(4) = 1;
			
			float score = test.transpose() * subMat * test;
			if(score < bestDot){
				bestDot = score;
				bestVec = test;
			}
		}
		if(bestDot >= 2){
			//TODO the lifted 3 constraint
			return *new std::vector<double>;
		}
		
		std::vector<double>& res = *new std::vector<double>(11);
		res[1] = bestVec(1)*bestVec(0);
		res[2] = bestVec(2)*bestVec(0);
		res[3] = bestVec(2)*bestVec(1);
		res[4] = bestVec(3)*bestVec(0);
		res[5] = bestVec(3)*bestVec(1);
		res[6] = bestVec(3)*bestVec(2);
		res[7] = bestVec(4)*bestVec(0);
		res[8] = bestVec(4)*bestVec(1);
		res[9] = bestVec(4)*bestVec(2);
		res[10] = bestVec(4)*bestVec(3);
		res[0] = -2;
		return res;
	}
	printf("Unhandled size %d, returning no constraint.\n", (int)subMat.rows());
	return *new std::vector<double>;
}

void satToQP(uint32_t variables, std::vector<wclause> clauses, std::vector<float> literalWeights, uint32_t nQP, MatrixXd& Amat, float& constantTerm){
	for(uint32_t i=1;i<nQP;i++){
		//x with weight w corresponds to w/2 + (w/2)x
		float w = literalWeights[i-1];
		Amat(0, i) = w/2.;
		constantTerm += w/2.;
	}
	for(uint32_t c=0;c<clauses.size();c++){
		// (x OR y) with weight w on {-1,+1}
		//corresponds to (w/4)x + (w/4)y - (w/4)xy + 3/4w.
		//Negating a variable means negating appropriate term.
		float w = clauses[c].second;
		int xL = clauses[c].first.first;
		int yL = clauses[c].first.second;
		int xV = abs(xL), yV = abs(yL);
		Amat(0, xV) += (w/4)*(xL > 0 ? 1 : -1);
		Amat(0, yV) += (w/4)*(yL > 0 ? 1 : -1);
		Amat(std::min(xV,yV), std::max(xV,yV)) += -(w/4)*(xL > 0 ? 1 : -1)*(yL > 0 ? 1 : -1);
		constantTerm += 0.75*w;
	}
}

int main(void)
{
	//Literals are numbered [-variables,-1] and [1,variables].
	//A negative value indicates a negation clause.
	uint32_t variables;
	
	 //Each term is <(v1,v2), weight>, representing v1 OR v2.
	std::vector<wclause> clauses;
	std::vector<float> literalWeights;
	
	setSampleProblem(variables, clauses, literalWeights);
	
	//Number of variables in quadratic program
	uint32_t nQP = variables+1;
	
	//Upper diagonal form of the optimization matrix
	MatrixXd Amat(nQP,nQP);
	Amat.setZero();
	float constantTerm = 0;
	
	satToQP(variables, clauses, literalWeights, nQP, Amat, constantTerm);
	
	printf("Optimization A matrix: \n");
	std::cout << Amat << std::endl;
	
	//Number of variables in linear program
	uint32_t nLP = nQP*(nQP-1)/2;
	std::vector<float> currSol;
	currSol.resize(nLP);
	
	glp_prob *lp;
	lp = glp_create_prob();
	
	glp_set_prob_name(lp, "CLQO"); //Name the problem
	glp_set_obj_dir(lp, GLP_MAX);  //We maximize
	
	//each column corresponds to a variable we (linearly) optimize over
	glp_add_cols(lp, nLP);
	for(uint32_t i=1;i<=nLP;i++){
		//Each variable has DB (double bound) for [-1,+1]
		glp_set_col_bnds(lp, i, GLP_DB, -1., 1.);
		std::pair<int,int> xy = getQPVars(i, nLP);
		
		//Set objective weight
		glp_set_obj_coef(lp, i, Amat(xy.second, xy.first));
	}
	
	int succConstraintFail = 0;
	
solve:
	/* solve problem */
	int simplex_err = glp_simplex(lp, NULL);
	if(simplex_err != 0)
		printf("FAILED Error Code = %d\n", simplex_err);
	if(glp_get_status(lp) != GLP_OPT)
		printf("Simplex Optimality FAILED");
	
	for(uint32_t i=1;i<=nLP;i++){
		currSol[i-1] = glp_get_col_prim(lp, i);
		//printf("%f%s", currSol[i-1], i==nLP ? "\n" : ", ");
	}
	
findcore:
	std::vector<int> core = nonPSDcore(nQP, nLP, currSol);
	if(core.size() == 0){
		goto complete;
	} else {
		//We have a core, identify a new constraint
		MatrixXd coreMat = getSubmatrix(core, currSol);
		std::vector<double>& constraint = findConstraint(coreMat);
		if(constraint.size() == 0){
			succConstraintFail++;
			if(succConstraintFail == CONSTRAINT_FAIL_LIMIT)
				goto rounding;
			else
				goto findcore;
		} else {
			succConstraintFail=0;
			
			puts("Applying constraint");
			//apply the constraint
			//constraint consists of linear terms on the
			//(core.size() choose 2) variables, preceded by
			//a constant term. The variables are in the submatrix,
			//and they'll need to mapped back "up" to the full matrix:
			//turn ij into i and j, and then map i and j to the larger matrix,
			//and then map back down to ij.
			glp_add_rows(lp, 1);
			int rowNum = glp_get_num_rows(lp);
			//build indices to pass to GLPK's sparse representation
			int indices[constraint.size()]; //+1 for 1 indexing, -1 for ignoring constant
			for(uint32_t v=0; v<constraint.size()-1; v++){
				std::pair<int,int> ij = getQPVars(1+v, constraint.size()-1);
				int largeI = core[ij.first], largeJ = core[ij.second];
				int largeIJ = getLPVar(largeI, largeJ);
				indices[1+v] = largeIJ;
			}
			//printf("%d, %d, %d, %d, %d\n", constraint.size()-1, indices[0], indices[1], indices[2]);
			glp_set_mat_row(lp, rowNum, constraint.size()-1, indices, &(constraint[0]));
			//sum[ coeff[i]*x[i] ] >= coeff[-1]
			glp_set_row_bnds(lp, rowNum, GLP_LO, constraint[0], 0.0);
			goto solve;
		}
		
		
	rounding:
		puts("No more constraints detected. Rounding off.");
		roundMat(lp, currSol, Amat, nQP, nLP, constantTerm);
		goto cleanup;
	}
	
complete:
	puts("Global optimum found! Solution:");
	for(uint32_t i=1;i<nQP;i++){
		printf("%d, ", (int)lround(currSol[getLPVar(0,i)-1]));
	}
	printf("\n");
	
	{
		float score = constantTerm;
		printf("constant %f\n", constantTerm);
		for(int i=0;i<nQP;i++){
			for(int j=i+1;j<nQP;j++){
				score += (i==0 ? 1 : currSol[getLPVar(0,i)-1]) * currSol[getLPVar(0,j)-1] * Amat(i,j);
			}
		}
		printf("Score = %f\n", score);
	}
	
	goto cleanup;
	
cleanup:
	glp_delete_prob(lp);
	glp_free_env();
	return 0;
}