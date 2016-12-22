#ifndef UTILITIES_H_
#define UTILITIES_H_

double ceil441(double f);

double floor441(double f);

double interpolate(double point_1, double point_2, double value_1,
		double value_2, double quest_point);

double vector_magnitude(double* quest_vec);

void normalize_vector(double* quest_normal);

double dot_product(double* vector_1, double* vector_2);

void cross_product(double* vector_1, double* vector_2, double *cross_vec);

void print_vector(double *print_vector);

double cot(double angle);

void interpolate_vector(double point_1, double point_2, double* normal_1,
		double* normal_2, double quest_point, double* quest_normal);

#endif
