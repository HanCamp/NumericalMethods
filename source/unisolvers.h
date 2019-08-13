//
//  unisolvers.h
//  NumericalMethods
//
//  Created by Haniel Campos Alcantara Paulo on 8/11/19.
//  Copyright © 2019 Haniel Campos. All rights reserved.
//

#ifndef unisolvers_h
#define unisolvers_h

#include <stddef.h>
#include <stdbool.h>
#include "definitions.h"

// IMPLEMENTATIONS OF THE BISECTION METHOD
// CONVERGENCE: LINEAR (GUARANTEED)
double bisection(const univariate_function f, const double x0, const double x1, const unsigned int max_iter, const double precision, const bool verbose);

// IMPLEMENTATIONS OF UNIVARIATE LINEAR ITERATION
// CONVERGENCE: LINEAR
double linear_iteration(const univariate_function f, const double x0, const unsigned int max_iter, const double precision, const bool verbose);

// IMPLEMENTATIONS OF UNIVARIATE AITKEN'S Δ SQUARED PROCESS
// CONVERGENCE: LINEAR BUT ALWAYS BETTER THAN LINEAR ITERATION
double aitkens_delta(const univariate_function f, const double x0, const unsigned int max_iter, const double precision, const bool verbose);

// IMPLEMENTATIONS OF UNIVARIATE NEWTON'S METHOD
// CONVERGENCE: QUADRATIC
double newton(const univariate_function f, const univariate_function fp, const double x0, const unsigned int max_iter, const double precision, const bool verbose);

// IMPLEMENTATIONS OF SECANT METHOD
// CONVERGENCE: SUPERLINEAR BUT NOT QUADRATIC
double secant(const univariate_function f, const double x1, const double x0, const unsigned int max_iter, const double precision, const bool verbose);

// IMPLEMENTATION OF FALSE POSITION METHOD
// CONVERGENCE: SUPERLINEAR BUT NOT QUADRATIC (GUARANTEED)

double false_position(const univariate_function f, const double x1, const double x0, const unsigned int max_iter, const double precision, const bool verbose);

// IMPLEMENTATION OF MÜLLER'S METHOD
double muller(const univariate_function f, const double x2, const double x1, const double x0, const unsigned int max_iter, const double precision, const bool verbose);


#endif /* unisolvers_h */
