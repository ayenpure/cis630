#include <cmath>
#include "VectorOperations.h"

using std::sqrt;

void vector_copy(double *destination, double *source) {
	destination[0] = source[0];
	destination[1] = source[1];
	destination[2] = source[2];
}

double vector_magnitude(double* quest_vec) {
	return sqrt( (quest_vec[0] * quest_vec[0])
						 + (quest_vec[1] * quest_vec[1])
						 + (quest_vec[2] * quest_vec[2]));
}

void normalize_vector(double* quest_normal) {
	double norm = vector_magnitude(quest_normal);
	quest_normal[0] = quest_normal[0] / norm;
	quest_normal[1] = quest_normal[1] / norm;
	quest_normal[2] = quest_normal[2] / norm;
}

double dot_product(double* vector_1, double* vector_2) {
	return (vector_1[0] * vector_2[0]) + (vector_1[1] * vector_2[1])
			+ (vector_1[2] * vector_2[2]);
}

void cross_product(double* vector_1, double* vector_2, double *cross_vec) {
	double product[3] = { (vector_1[1] * vector_2[2])
			- (vector_1[2] * vector_2[1]), (vector_1[2] * vector_2[0])
			- (vector_1[0] * vector_2[2]), (vector_1[0] * vector_2[1])
			- (vector_1[1] * vector_2[0]) };
	vector_copy(cross_vec, product);
}
