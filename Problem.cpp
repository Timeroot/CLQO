#include "Problem.hpp"

#include <iostream>

Problem::Problem(uint32_t n) : nQP(n), coeffs(*new MatrixXd(n,n)), constantTerm(0) {coeffs.setZero();}

Problem::Problem(uint32_t n, MatrixXd& coeff, float cT) : nQP(n), coeffs(coeff), constantTerm(cT) {}

Problem* Problem::from2SAT(uint32_t n, std::vector<clause2>& clauses, std::vector<float>& literalWeights){
	uint32_t nQP = n+1;
	
	Problem* res = new Problem(nQP);
	
	for(uint32_t i=1;i<nQP;i++){
		//x with weight w corresponds to w/2 + (w/2)x
		float w = literalWeights[i-1];
		res->coeffs(0, i) = w/2.;
		res->constantTerm += w/2.;
	}
	for(uint32_t c=0;c<clauses.size();c++){
		// (x OR y) with weight w on {-1,+1}
		//corresponds to (w/4)x + (w/4)y - (w/4)xy + 3/4w.
		//Negating a variable means negating appropriate term.
		float w = std::get<0>(clauses[c]);
		int xL = std::get<1>(clauses[c]);
		int yL = std::get<2>(clauses[c]);
		int xV = abs(xL), yV = abs(yL);
		res->coeffs(0, xV) += (w/4)*(xL > 0 ? 1 : -1);
		res->coeffs(0, yV) += (w/4)*(yL > 0 ? 1 : -1);
		res->coeffs(std::min(xV,yV), std::max(xV,yV)) += -(w/4)*(xL > 0 ? 1 : -1)*(yL > 0 ? 1 : -1);
		res->constantTerm += 0.75*w;
	}
	
	return res;
}

/* Thoughts on encoding logic gate constraints faster
x = AND(y,z)

x <= y
x <= z
3SAT: x|!y|!z

Allowed space:
x y z SCORE
0 0 0 1
0 0 1 1
0 1 0 1
0 1 1 0
1 1 1 1

Score func:
1 + x - (y*z)

In {-1,+1} version, each variable becomes (1+var)/2, so:

Allowed space:
x y z SCORE
- - - 1
- - + 1
- + - 1
- + + 0
+ + + 1

Score func:
5/4 + x/2 - y/4 - z/4 - (y*z)/4

But actually we have several options for how we want to score it, since we have 5 points to fix and 7 degrees of freedom. The following two expressions can be freely added or subtracted without affecting the true optimum:
1 + x - z - xz
1 + x - y - xy
or their sum + difference might be more attractive,
2 + 2x - y - z - xy - xz
z - y + xz - xy
If you wanted to minimize the L2 norm of the objective matrix (a reasonable goal, maybe), the optimum is to add -3/16 copies of (2x - y - z - xy - xz), so that the new objective is
7/8 + x/8 - y/16 - z/16 + (x*y)*3/16 + (x*z)*3/16 - (y*z)/4
If you wante the minimize the L1 norm of the objective, the optimum is -1/4, which gives a nice sparse objective:
3/4 + (x*y)/4 + (x*z)/4 - (y*z)/4
*/

Problem* Problem::from3SAT(uint32_t n, std::vector<clause3>& clause3s, std::vector<clause2>& clause2s, std::vector<float>& literalWeights){
	uint32_t m = clause3s.size(); //number of 3-clauses -> auxiliary variables
	
	//For each 3-clause, add a variable v
	//10 clauses: (x_i)(x_j)(x_k)(v)
	//       (~x_i|~x_j)(~x_j|~x_k)(~x_i|~x_k)
	//       (x_i|~v)(x_j|~v)(x_k|~v)
	for(uint32_t c=0; c<clause3s.size(); c++){
		auto c3 = clause3s[c];
		float w = std::get<0>(c3);
		
		int xi = std::get<1>(c3);
		int xj = std::get<2>(c3);
		int xk = std::get<3>(c3);
		int v = 1+n+c;
		
		
		literalWeights[abs(xi)] += w*(xi > 0 ? 1 : -1);
		literalWeights[abs(xj)] += w*(xj > 0 ? 1 : -1);
		literalWeights[abs(xk)] += w*(xk > 0 ? 1 : -1);
		literalWeights.push_back(w); //v
		
		clause2s.push_back({w, -xi,-xj});
		clause2s.push_back({w, -xj,-xk});
		clause2s.push_back({w, -xi,-xk});
		
		clause2s.push_back({w, xi,-v});
		clause2s.push_back({w, xj,-v});
		clause2s.push_back({w, xk,-v});
	}
	
	//Get problem as MAX2SAT instance
	Problem* res = from2SAT(n+m, clause2s, literalWeights);
	
	//Adjust constant -- the above reduction turns (0,w) into (6w,7w).
	for(uint32_t c=0; c<clause3s.size(); c++){
		auto c3 = clause3s[c];
		float w = std::get<0>(c3);
		res->constantTerm -= 6*w;
	}
	
	//cleanup modifications to literalWeights/clauses2
	literalWeights.erase(literalWeights.begin()+n, literalWeights.end());
	for(uint32_t c=0; c<clause3s.size(); c++){
		auto c3 = clause3s[c];
		float w = std::get<0>(c3);
		
		int xi = std::get<1>(c3);
		int xj = std::get<2>(c3);
		int xk = std::get<3>(c3);
		
		literalWeights[abs(xi)] -= w;
		literalWeights[abs(xj)] -= w;
		literalWeights[abs(xk)] -= w;
	}
	
	clause2s.erase(clause2s.end()-6*m, clause2s.end());
	
	return res;
}

//Turn MAX-CLIQUE into a MAXQP problem. Given an n-vertex graph,
//add an auxiliary variable v0 representing "true". Each other variable
//corresponds to a vertex. It has a weight of 1 (so that v0*vi has a
//coefficient of 1), representing a value of 1 of adding it to the clique.
//We don't want to add any two vertices in the 'clique' to not share an edge,
//so we have a penality of -k if they aren't connected in the adjacency matrix.
//Any value of k > 1 will suffice, to fully "discourage" adding a bad
//vertex. We take k=2 so that we get nice integral values along the way.
Problem* Problem::fromMaxClique(uint32_t n, bool** adjMat){
	uint32_t nQP = n+1;
	
	Problem* res = new Problem(nQP);
	
	for(uint32_t i=1;i<nQP;i++){
		res->coeffs(0, i) = 0.5;
		res->constantTerm += 0.5;
	}
	float k = 2;
	for(uint32_t i=1;i<nQP;i++){
		for(uint32_t j=i+1;j<nQP;j++){
			if(!adjMat[i][j]){
				//add a penalty
				//-k iff x==1 && y==1, or
				//-k/4(1 + x + y + xy)
				res->coeffs(0, i) -= k/4;
				res->coeffs(0, j) -= k/4;
				res->coeffs(i, j) -= k/4;
				res->constantTerm -= k/4;
			}
		}
	}
	
	return res;
}

Problem* Problem::fromIndSet(uint32_t n, bool** adjMat){
	uint32_t nQP = n+1;
	
	Problem* res = new Problem(nQP);
	
	for(uint32_t i=1;i<nQP;i++){
		res->coeffs(0, i) = 0.5;
		res->constantTerm += 0.5;
	}
	float k = 2;
	for(uint32_t i=1;i<nQP;i++){
		for(uint32_t j=i+1;j<nQP;j++){
			if(adjMat[i][j]){
				//add a penalty
				//-k iff x==1 && y==1, or
				//-k/4(1 + x + y + xy)
				res->coeffs(0, i) -= k/4;
				res->coeffs(0, j) -= k/4;
				res->coeffs(i, j) -= k/4;
				res->constantTerm -= k/4;
			}
		}
	}
	
	return res;
}

float Problem::score(VectorXd sol){
	return constantTerm + sol.transpose() * coeffs * sol;
}
