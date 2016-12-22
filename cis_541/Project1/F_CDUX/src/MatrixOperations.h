#ifndef MATRIXOPERATIONS_H_
#define MATRIXOPERATIONS_H_

#include <iostream>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <sstream>

//Defining all data structures
class Matrix {
public:
	double A[4][4];

  static Matrix ComposeMatrices(const Matrix &M1, const Matrix &M2) {
  	Matrix rv;
  	for (int i = 0; i < 4; i++)
  		for (int j = 0; j < 4; j++) {
  			rv.A[i][j] = 0;
  			for (int k = 0; k < 4; k++)
  				rv.A[i][j] += M1.A[i][k] * M2.A[k][j];
  		}

  	return rv;
  }

  void TransformPoint(const double *ptIn, double *ptOut) {
  	ptOut[0] = ptIn[0] * A[0][0] + ptIn[1] * A[1][0] + ptIn[2] * A[2][0]
  			+ ptIn[3] * A[3][0];
  	ptOut[1] = ptIn[0] * A[0][1] + ptIn[1] * A[1][1] + ptIn[2] * A[2][1]
  			+ ptIn[3] * A[3][1];
  	ptOut[2] = ptIn[0] * A[0][2] + ptIn[1] * A[1][2] + ptIn[2] * A[2][2]
  			+ ptIn[3] * A[3][2];
  	ptOut[3] = ptIn[0] * A[0][3] + ptIn[1] * A[1][3] + ptIn[2] * A[2][3]
  			+ ptIn[3] * A[3][3];
  }

};

#endif
