#ifndef VECTOROPERATIONS_H_
#define VECTOROPERATIONS_H_

void vector_copy(double *destination, double *source);
double vector_magnitude(double* quest_vec);
void normalize_vector(double* quest_normal);
double dot_product(double* vector_1, double* vector_2);
void cross_product(double* vector_1, double* vector_2, double *cross_vec);

#endif
