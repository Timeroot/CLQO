#include "Problem.hpp" 

Problem::Problem(uint32_t n) : nQP(n), coeffs(*new MatrixXd(n,n)), constantTerm(0) {coeffs.setZero();}

Problem::Problem(uint32_t n, MatrixXd& coeff, float cT) : nQP(n), coeffs(coeff), constantTerm(cT) {}

Problem* Problem::from2SAT(uint32_t n, std::vector<wclause>& clauses, std::vector<float>& literalWeights){
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
		float w = clauses[c].second;
		int xL = clauses[c].first.first;
		int yL = clauses[c].first.second;
		int xV = abs(xL), yV = abs(yL);
		res->coeffs(0, xV) += (w/4)*(xL > 0 ? 1 : -1);
		res->coeffs(0, yV) += (w/4)*(yL > 0 ? 1 : -1);
		res->coeffs(std::min(xV,yV), std::max(xV,yV)) += -(w/4)*(xL > 0 ? 1 : -1)*(yL > 0 ? 1 : -1);
		res->constantTerm += 0.75*w;
	}
	
	return res;
}

float Problem::score(VectorXd sol){
	return constantTerm + sol.transpose() * coeffs * sol;
}