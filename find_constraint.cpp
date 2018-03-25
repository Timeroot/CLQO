#include "LPSolver.hpp"

#include <iostream>

using Eigen::Vector3d;
using Eigen::Vector4d;

std::vector<double>& findConstraint_3(MatrixXd& subMat);
std::vector<double>& findConstraint_4(MatrixXd& subMat);
std::vector<double>& findConstraint_5(MatrixXd& subMat);

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
	printf("Unhandled size %d, returning no constraint.\n", (int)subMat.rows());
	return *new std::vector<double>;
}

std::vector<double>& findConstraint_3(MatrixXd& subMat){
	Vector3d bestVec; 
	bestVec(0) = bestVec(1) = bestVec(2) = 1;
	float bestDot = 100; //maximum 'useful' score of 1
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
	//std::cout << "Constraint HP3:3 ~ " << bestDot << std::endl;
	return res;
}

std::vector<double>& findConstraint_4(MatrixXd& subMat){
	Vector4d bestVec; 
	bestVec(0) = bestVec(1) = bestVec(2) = bestVec(3) = 1;
	float bestDot = 100; //maximum 'useful' score of 1
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
	//std::cout << "Constraint HP4:3 ~ " << bestDot << std::endl;
	return res;
}

std::vector<double>& findConstraint_5(MatrixXd& subMat){
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
	std::vector<double>& res = *new std::vector<double>(11);
	
	if(bestDot < 2){
		//std::cout << "Constraint HP5:5 ~ " << bestDot << std::endl;
		res[0] = -2;
	} else {
		//It's option two: the lifted 3 constraint
		//In practice this seems to occur only _very_ rarely!!
		bestDot = 2; //maximum 'useful' score of 1
		for(int dir=0; dir<5; dir++){
			for(int dir2=dir+1; dir2<5; dir2++){
				for(int bits=0; bits<4; bits++){
					VectorXd test(5);
					
					test(0) = ((bits>>0) & 1)*2 - 1;
					test(1) = ((bits>>1) & 1)*2 - 1;
					test(2) = ((bits>>2) & 1)*2 - 1;
					test(3) = 0;
					test(4) = 0;
					
					{float temp = test(dir);
					test(dir) = test(3);
					test(3) = temp;}
					
					
					{float temp = test(dir2);
					test(dir2) = test(4);
					test(4) = temp;}
					
					float score = test.transpose() * subMat * test;
					if(score < bestDot){
						bestDot = score;
						bestVec = test;
					}
				}
			}
		}
		//std::cout << "Constraint HP5:3 ~ " << bestDot  << std::endl;
		res[0] = -1;
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