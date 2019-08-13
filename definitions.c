//
//  definititions.c
//  NummetC
//
//  Created by Haniel Campos Alcantara Paulo on 8/12/19.
//  Copyright © 2019 Haniel Campos. All rights reserved.
//

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "definitions.h"

double binomial(const double n, const double k) {
  return tgamma(n + 1) / (tgamma(k + 1) * tgamma(n - k + 1));
} 

struct Matrix* new_matrix(const unsigned int nrows, const unsigned int ncols) {
  struct Matrix* newMat = malloc(sizeof(struct Matrix));
  newMat->array = malloc(nrows * sizeof(double*));
  newMat->arrayTranspose = malloc(ncols * sizeof(double*));
  newMat->ncols = ncols;
  newMat->nrows = nrows;
  for (unsigned int i = 0; i < nrows; i++)
    (newMat->array)[i] = malloc(ncols * sizeof(double));
  for (unsigned int i = 0; i < ncols; i++)
    (newMat->arrayTranspose)[i] = malloc(nrows * sizeof(double));
  return newMat;
}

struct Matrix* new_matrix_fill(const unsigned int nrows, const unsigned int ncols, const double fillVal) {
  struct Matrix* newMat = malloc(sizeof(struct Matrix));
  newMat->array = malloc(nrows * sizeof(double*));
  newMat->arrayTranspose = malloc(ncols * sizeof(double*));
  newMat->ncols = ncols;
  newMat->nrows = nrows;
  for (unsigned int i = 0; i < nrows; i++) {
    (newMat->array)[i] = malloc(ncols * sizeof(double));
    for (unsigned int j = 0; j < ncols; j++)
      (newMat->array)[i][j] = fillVal;
  }
  for (unsigned int i = 0; i < ncols; i++) {
    (newMat->arrayTranspose)[i] = malloc(nrows * sizeof(double));
    for (unsigned int j = 0; j < nrows; j++)
      (newMat->array)[i][j] = fillVal;
  }
  return newMat;
}

struct Matrix* new_matrix_eye(const unsigned int nrows, const unsigned int ncols) {
  struct Matrix* newMat = malloc(sizeof(struct Matrix));
  newMat->array = malloc(nrows * sizeof(double*));
  newMat->arrayTranspose = malloc(ncols * sizeof(double*));
  newMat->ncols = ncols;
  newMat->nrows = nrows;
  for (unsigned int i = 0; i < nrows; i++) {
    (newMat->array)[i] = malloc(ncols * sizeof(double));
    for (unsigned int j = 0; j < ncols; j++) {
      if (i == j)
        (newMat->array)[i][j] = 1.0;
      else
        (newMat->array)[i][j] = 0.0;
    }
  }
  for (unsigned int i = 0; i < ncols; i++) {
    (newMat->arrayTranspose)[i] = malloc(nrows * sizeof(double));
    for (unsigned int j = 0; j < nrows; j++) {
      if (i == j)
        (newMat->array)[i][j] = 1.0;
      else
        (newMat->array)[i][j] = 0.0;
    }
  }
  return newMat;
}

void print_vector(const double* v, const unsigned int dim, double precision) {
  printf("[");
  for (unsigned int i = 0; i < dim; i++)
    printf(" %.*f", (int) -floor(log10(precision)), v[i]);
  printf(" ]\n");
}

double difference_norm(const double* v1, const double* v2, const unsigned int dim) {
  double val = 0;
  for (unsigned int i = 0; i < dim; i++)
    val += (v1[i] - v2[i]) * (v1[i] - v2[i]);
  return sqrt(val);
}

double inner_product(const double* v1, const double* v2, const unsigned int dim) {
  double val = 0;
  for (unsigned int i = 0; i < dim; i++)
    val += v1[i] * v2[i];
  return val;
}

struct Matrix* matrix_product(const struct Matrix* m1, const struct Matrix* m2, struct Matrix* ptr) {
  ptr->nrows = m1->nrows;
  ptr->ncols = m2->ncols;
  double** matrixArray = malloc(ptr->nrows * sizeof(double*));
  for (unsigned int i = 0; i < ptr->nrows; i++) {
    matrixArray[i] = malloc(ptr->ncols * sizeof(double));
    for (unsigned int j = 0; j < ptr->ncols; j++)
      matrixArray[i][j] = inner_product((m1->array)[i], (m2->arrayTranspose)[j], m1->ncols);
  }
  ptr->array = matrixArray;
  double** matrixArrayTranspose = malloc(ptr->ncols * sizeof(double*));
  for (unsigned int i = 0; i < ptr->ncols; i++) {
    matrixArrayTranspose[i] = malloc(ptr->nrows * sizeof(double));
    for (unsigned int j = 0; j < ptr->nrows; j++)
      matrixArrayTranspose[i][j] = ptr->array[j][i];
  }
  ptr->arrayTranspose = matrixArrayTranspose;
  return ptr;
}

struct DifferenceTable* new_difference_table(const double* val, const unsigned int npoints) {
  struct DifferenceTable* newTable = malloc(sizeof(struct DifferenceTable));
  newTable->npoints = npoints;
  newTable->table = malloc(npoints * sizeof(double*));
  for (int i = 0; i < npoints; i++) {
    (newTable->table)[i] = malloc((npoints - i - 1) * sizeof(double));
    for (int j = 1; j < npoints - i; j++)
      (newTable->table)[i][j] = forward_difference(i, j, val, npoints);
  }
  return newTable;
}

double forward_difference(const unsigned int i, const unsigned int degree, const double* vals, const unsigned int npoints) {
  if (i + degree > npoints - 1 || degree < 1) {
    fprintf(stderr, "ERROR: Invalid difference operator degree\n");
    exit(1);
  }
  if (degree == 1)
    return vals[i + 1] - vals[i];
  else
    return forward_difference(i + 1, degree - 1, vals, npoints) - forward_difference(i, degree - 1, vals, npoints);
}

double difference(const struct DifferenceTable* table, const unsigned int index, const unsigned int degree) {
  return (table->table)[index][degree - 1];
}

double central_difference(const struct DifferenceTable* table, const unsigned int index, const unsigned int degree) {
  if (degree % 2) {
    fprintf(stderr, "Central difference is undefined for given degree");
    exit(1);
  }
  return difference(table, index - degree / 2, degree);
}
