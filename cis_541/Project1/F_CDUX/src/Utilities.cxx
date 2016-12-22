#include <iostream>
#include <cmath>
#include <algorithm>
#include <string.h>
#include <vector>
#include <sstream>
#include "Utilities.h"

using std::cout;
using std::endl;
using std::min;
using std::max;
using std::abs;
using std::pow;
using std::tan;
using std::sin;

double ceil441(double f) {
	return ceil(f - 0.00001);
}

double floor441(double f) {
	return floor(f + 0.00001);
}

/**
 * This function interpolates the value at quest_point based on
 * values value_1 and value_2 at points point_1 and point_2 respectively
 */
double interpolate(double point_1, double point_2, double value_1,
		double value_2, double quest_point) {
	double diff = point_2 - point_1;
	double change_prop =
			(diff != 0) ? (quest_point - point_1) / (point_2 - point_1) : 0;
	double value_at_quest = value_1 + (change_prop * (value_2 - value_1));
	return value_at_quest;
}

double vector_magnitude(double* quest_vec) {
	return sqrt( (quest_vec[0] * quest_vec[0])
						 + (quest_vec[1] * quest_vec[1])
						 + (quest_vec[2] * quest_vec[2]));
}

void normalize_vector(double* quest_normal) {
	double norm = sqrt(
			(quest_normal[0] * quest_normal[0])
					+ (quest_normal[1] * quest_normal[1])
					+ (quest_normal[2] * quest_normal[2]));
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
	memcpy(cross_vec, product, 3 * sizeof(double));
}

void print_vector(double *print_vector) {
	cout << "{" << print_vector[0] << ", " << print_vector[1] << ", "
			<< print_vector[2] << "}" << endl;
}

double cot(double angle) {
	return (1 / tan(angle));
}

/**
 * This function interpolates over the vectors normal_1 and normal_2 to yield
 * the vector quest_normal for the point (quest_point[0],quest_point[1]) in space.
 */
void interpolate_vector(double point_1, double point_2, double* normal_1,
		double* normal_2, double quest_point, double* quest_normal) {
	double diff = point_2 - point_1;
	double proportion =
			(diff != 0) ? (quest_point - point_1) / (point_2 - point_1) : 0;
	double diff_vector[3] = { normal_2[0] - normal_1[0], normal_2[1]
			- normal_1[1], normal_2[2] - normal_1[2] };
	diff_vector[0] = proportion * diff_vector[0];
	diff_vector[1] = proportion * diff_vector[1];
	diff_vector[2] = proportion * diff_vector[2];
	quest_normal[0] = normal_1[0] + diff_vector[0];
	quest_normal[1] = normal_1[1] + diff_vector[1];
	quest_normal[2] = normal_1[2] + diff_vector[2];
	normalize_vector(quest_normal);
}
