#include <iostream>
#include <cmath>
#include <cuda.h>
#include <cuda_runtime_api.h>
#define DIM 127

__global__
void  calculateFieldValues(float*** fieldValues) {
	int x = blockIdx.x;
	int y = blockIdx.y;
	int z = blockIdx.z;
	float x_trans = x / 127.;
	float y_trans = y / 127.;
	float z_trans = z / 127.;
	float fieldValue = sqrtf(powf(0 - x_trans,2) 
				+ powf(0 - y_trans, 2) 
				+ powf(0 - z_trans, 2));
	fieldValues[x][y][z] = fieldValue;
}

int main(int argc, char* argv[]) {

	float*** fieldValues;
	cudaMalloc(fieldValues, DIM*sizeof(float**));
	for(int i = 0; i < DIM; i++) {
		cudaMalloc(fieldValues[i], DIM*sizeof(float*))
		for(int j = 0; j < DIM; j++) {
			cudaMalloc(fieldValues[i][j], DIM*sizeof(float));
		}
	}
	calculateFieldValues(fieldValues);
}
