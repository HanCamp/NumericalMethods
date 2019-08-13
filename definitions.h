//
//  definitions.h
//  NummetC
//
//  Created by Haniel Campos Alcantara Paulo on 8/12/19.
//  Copyright © 2019 Haniel Campos. All rights reserved.
//

#ifndef definitions_h
#define definitions_h

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>

typedef double (*univariate_function)(const double);
typedef void (*multivariate_function)(const double*, double*);
typedef void(*matrix_function)(const double*, double**);

double binomial(const double n, const double k);

struct Matrix {
  double** array;
  double** arrayTranspose;
  unsigned int nrows;
  unsigned int ncols;
};

struct Matrix* new_matrix(const unsigned int nrows, const unsigned int ncols);

struct Matrix* new_matrix_fill(const unsigned int nrows, const unsigned int ncols, const double fillVal);

struct Matrix* new_matrix_eye(const unsigned int nrows, const unsigned int ncols);

void print_vector(const double* v, const unsigned int dim, double precision);

double difference_norm(const double* v1, const double* v2, const unsigned int dim);

double inner_product(const double* v1, const double* v2, const unsigned int dim);

struct Matrix* matrix_product(const struct Matrix* m1, const struct Matrix* m2, struct Matrix* ptr);

struct DifferenceTable {
  double** table;
  unsigned long npoints;
};

struct DifferenceTable* new_difference_table(const double* val, const unsigned int npoints);

double forward_difference(const unsigned int i, const unsigned int degree, const double* vals, const unsigned int npoints);

double difference(const struct DifferenceTable* table, const unsigned int index, const unsigned int degree);

double central_difference(const struct DifferenceTable* table, const unsigned int index, const unsigned int degree);






#endif /* definitions_h */
