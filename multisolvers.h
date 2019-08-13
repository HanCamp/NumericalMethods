//
//  multisolvers.h
//  NummetC
//
//  Created by Haniel Campos Alcantara Paulo on 8/12/19.
//  Copyright © 2019 Haniel Campos. All rights reserved.
//

#ifndef multisolvers_h
#define multisolvers_h

#include <stddef.h>
#include <stdbool.h>
#include "definitions.h"

// IMPLEMENTATION OF MULTIVARIATE LINEAR ITERATION
// CONVERGENCE: LINEAR
void linear_iteration_multi(const multivariate_function f, double* x0, double* tmp, const unsigned int dimension, const unsigned int max_iter, const double precision, const bool verbose);

// IMPLEMENTATION OF MULTIVARIATE AITKEN'S Δ SQUARED PROCESS
// CONVERGENCE: LINEAR
void aitkens_delta_multi(const multivariate_function f, double* x0, double* tmp1, double* tmp2, const unsigned int dimension, const unsigned int max_iter, const double precision, const bool verbose);

// IMPLEMENTATION OF MULTIVARIATE NEWTON'S METHOD
// CONVERGENCE: QUADRATIC
void newton_multi(const multivariate_function f, const matrix_function J, double* x0, double* tmp, const unsigned int dimension, const unsigned int max_iter, const double precision, const bool verbose);

#endif /* multisolvers_h */
