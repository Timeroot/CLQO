#include "LPSolver.hpp"

#include <iostream>

using Eigen::Vector3d;
using Eigen::Vector4d;

std::vector<double>& findConstraint_3(MatrixXd& subMat);
std::vector<double>& findConstraint_4(MatrixXd& subMat);
std::vector<double>& findConstraint_5(MatrixXd& subMat);
std::vector<double>& findConstraint_6(MatrixXd& subMat);

//Debugging info in case a constraint isn't found when it should be
void fail_constraint(MatrixXd& subMat, const char* name, float bestDot);

//Returns a constraint that the submatrix violates, or 0-length
//if none is found.
std::vector<double>& LPSolver::findConstraint(MatrixXd& subMat){
	if(subMat.rows() == 3){
		return findConstraint_3(subMat);
	}
	if(subMat.rows() == 4){
		return findConstraint_4(subMat);
	}
	if(subMat.rows() == 5){
		return findConstraint_5(subMat);
	}
	/*if(subMat.rows() == 6){
		return findConstraint_6(subMat);
	}*/
	printf("Unhandled size %d, returning no constraint.\n", (int)subMat.rows());
	return *new std::vector<double>;
}

void fail_constraint(MatrixXd& subMat, const char* name, float bestDot){
	std::cout << "Constraint finding assertion error!" << std::endl;
	Eigen::SelfAdjointEigenSolver<MatrixXd> eig(subMat);
	auto evals = eig.eigenvalues();
	for(int i=0;i<evals.size();i++)
		std::cout << evals[i] << std::endl;
	std::cout << "Because" << std::endl;
	std::cout << subMat << std::endl;
	std::cout << "Constraint " << name << " ~ " << bestDot << std::endl;
	exit(3);
}

std::vector<double>& findConstraint_3(MatrixXd& subMat){
	Vector3d bestVec;
	bestVec(0) = bestVec(1) = bestVec(2) = 1;
	float bestDot = 0; //minimum 'useful' score of -1
	for(int bits=0; bits<4; bits++){
		Vector3d test;
		test(0) = ((bits>>0) & 1)*2 - 1;
		test(1) = ((bits>>1) & 1)*2 - 1;
		test(2) = 1;
		float score = test.transpose() * subMat * test;
		score = (score - 3)/2;
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
	
	if(bestDot >= -1)
		fail_constraint(subMat, "3", bestDot);
	return res;
}

std::vector<double>& findConstraint_4(MatrixXd& subMat){
	Vector4d bestVec; 
	bestVec(0) = bestVec(1) = bestVec(2) = bestVec(3) = 1;
	float bestDot = 0; //minimum 'useful' score of -1
	for(int dir=0; dir<4; dir++){
			for(int bits=0; bits<4; bits++){
			Vector4d test;
			test(0) = ((bits>>0) & 1)*2 - 1;
			test(1) = ((bits>>1) & 1)*2 - 1;
			test(2) = 1;
			test(3) = 0;
			
			float temp = test(dir);
			test(dir) = test(3);
			test(3) = temp;
			
			float score = test.transpose() * subMat * test;
			score = (score - 3)/2;
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
	
	if(bestDot >= -1)
		fail_constraint(subMat, "4", bestDot);
	return res;
}

std::vector<double>& findConstraint_5(MatrixXd& subMat){
	VectorXd bestVec(5); 
	bestVec(0) = bestVec(1) = bestVec(2) = bestVec(3) = bestVec(4) = 1;
	float bestDot = 0; //minimum 'useful' score of -2
	for(int dir=0; dir<5; dir++){
		for(int dir2=dir+1; dir2<5; dir2++){
			for(int bits=0; bits<4; bits++){
				VectorXd test(5);
				
				test(0) = ((bits>>0) & 1)*2 - 1;
				test(1) = ((bits>>1) & 1)*2 - 1;
				test(2) = 1;
				test(3) = 0;
				test(4) = 0;
				
				{float temp = test(dir);
				test(dir) = test(3);
				test(3) = temp;}
				
				
				{float temp = test(dir2);
				test(dir2) = test(4);
				test(4) = temp;}
				
				float score = test.transpose() * subMat * test;
				score = (score - 3)/2;
				if(score < bestDot){
					bestDot = score;
					bestVec = test;
				}
			}
		}
	}
	std::vector<double>& res = *new std::vector<double>(11);
	res[0] = -1;
	
	if(bestDot > -1){
		//It's option two: the non-lifted 3 constraint
		bestDot = 0; //minimum 'useful' score of -1
		
		for(int bits=0; bits<16; bits++){
			VectorXd test(5);
			test(0) = ((bits>>0) & 1)*2 - 1;
			test(1) = ((bits>>1) & 1)*2 - 1;
			test(2) = ((bits>>2) & 1)*2 - 1;
			test(3) = ((bits>>3) & 1)*2 - 1;
			test(4) = 1;
			
			float score = test.transpose() * subMat * test;
			score = (score - 5)/2; //this line shouldn't be uncommented ... but leaving it out
			//(a bug) makes stuff run faster, but definitely still correctly. Understand why!
			if(score < bestDot){
				bestDot = score;
				bestVec = test;
			}
		}
		res[0] = -2;
	
		if(bestDot >= -2)
			fail_constraint(subMat, "5", bestDot);
	}
	
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
	return res;
}

std::vector<double>& findConstraint_6(MatrixXd& subMat){
	VectorXd bestVec(6); 
	bestVec(0) = bestVec(1) = bestVec(2) = bestVec(3) = bestVec(4) = bestVec(5) = 1;
	float bestDot = -4; //minimum 'useful' score of -4
	
	for(int dir=0; dir<6; dir++){
		for(int bits=0; bits<32; bits++){
			VectorXd test(6);
			test(0) = ((bits>>0) & 1)*2 - 1;
			test(1) = ((bits>>1) & 1)*2 - 1;
			test(2) = ((bits>>2) & 1)*2 - 1;
			test(3) = ((bits>>3) & 1)*2 - 1;
			test(4) = ((bits>>4) & 1)*2 - 1;
			test(5) = 2;
			
			float temp = test(dir);
			test(dir) = test(5);
			test(5) = temp;
			
			float score = test.transpose() * subMat * test;
			score = (score - 9) / 2;
			if(score < bestDot){
				bestDot = score;
				bestVec = test;
			}
		}
	}
	
	if(bestDot < -4){
		//std::cout << "Constraint HP6!! :) ~ " << bestDot << std::endl;
	} else {
		//std::cout << "Constraint 6 NOT FOUND :( ~ " << bestDot << std::endl;
		return *new std::vector<double>(0);
	}
	
	std::vector<double>& res = *new std::vector<double>(1 + 6*5/2);
	
	res[0] = -4;
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
	res[11] = bestVec(5)*bestVec(0);
	res[12] = bestVec(5)*bestVec(1);
	res[13] = bestVec(5)*bestVec(2);
	res[14] = bestVec(5)*bestVec(3);
	res[15] = bestVec(5)*bestVec(4);
	return res;
}