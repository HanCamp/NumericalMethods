//
//  unisolvers.c
//  NummetC
//
//  Created by Haniel Campos Alcantara Paulo on 8/11/19.
//  Copyright © 2019 Haniel Campos. All rights reserved.
//

#include <math.h>
#include <fenv.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include "unisolvers.h"
#include "definitions.h"

// BISECTION METHOD IMPLEMENTATION
double bisection(const univariate_function f, const double x0, const double x1, const unsigned int max_iter, const double precision, const bool verbose) {
  static unsigned int count = 1;
  if (!(f(x0) * f(x1) < 0)) {
    fprintf(stderr, "ERROR: Invalid arguments given, approximations must result in function values with opposite signs\n");
    exit(1);  
  }
  if (max_iter && count > max_iter) {
    fprintf(stderr, "Bisection method was unable to converge in %i iterations\n", max_iter);
    printf("Last value: %.*f\n", (int) -floor(log10(precision)), x0); 
    return x0;
  }
  if (fetestexcept(FE_INVALID)) {
    fprintf(stderr, "ERROR: Invalid argument detected (check for complex or outside domain arguments)\n");
    exit(1);
  }
  if (fetestexcept(FE_OVERFLOW) || fetestexcept(FE_UNDERFLOW)){
    fprintf(stderr, "ERROR: Floating point operations have gone outside of representable range\n");
    fprintf(stderr, "(check if process might be diverging)\n");
    exit(1);
  }
  if (fetestexcept(FE_DIVBYZERO)){
    fprintf(stderr, "ERROR: Division by zero detected\n");
    exit(1);
  }
  const double x2 = (x0 + x1) / 2;
  if (fabs(x0 - x2) < precision || fabs(x1 - x2) < precision) {
    if (verbose)
      printf("Bisection method converged to %.*f\n", (int) -floor(log10(precision)), x2);
    return x2;
  }
  if (verbose)
    fprintf(stdout, "Iteration #%i : %.*f\n", count, (int) -floor(log10(precision)), x2);
  count++;
  if (f(x2) * f(x0) < 0)
    return bisection(f, x2, x0, max_iter, precision, verbose);
  else
    return bisection(f, x2, x1, max_iter, precision, verbose);
}

// UNIVARIATE LINEAR ITERATION IMPLEMENTATION
double linear_iteration(const univariate_function f, const double x0, const unsigned int max_iter, const double precision, const bool verbose) {
  static unsigned int count = 1;
  if (max_iter && count > max_iter) {
    fprintf(stderr, "Linear iteration was unable to converge in %i iterations\n", max_iter);
    printf("Last value: %.*f\n", (int) -floor(log10(precision)), x0);
    return x0;
  }
  if (fetestexcept(FE_INVALID)) {
    fprintf(stderr, "ERROR: Invalid argument detected (check for complex or outside domain arguments)\n");
    exit(1);
  }
  if (fetestexcept(FE_OVERFLOW) || fetestexcept(FE_UNDERFLOW)){
    fprintf(stderr, "ERROR: Floating point operations have gone outside of representable range\n");
    fprintf(stderr, "(check if process might be diverging)\n");
    exit(1);
  }
  if (fetestexcept(FE_DIVBYZERO)){
    fprintf(stderr, "ERROR: Division by zero detected\n");
    exit(1);
  }
  const double x1 = f(x0);
  if (fabs(x0 - x1) < precision) {
    if (verbose)
      printf("Linear iteration converged to %.*f\n", (int) -floor(log10(precision)), x1);
    return x1;
  }
  if (verbose)
    printf("Iteration #%i\t : %.*f\n", count, (int) -floor(log10(precision)), x1);
  count++;
  return linear_iteration(f, x1, max_iter, precision, verbose);
}


// UNIVARIATE AITKEN'S Δ SQUARED PROCESS IMPLEMENTATION
double aitkens_delta(const univariate_function f, const double x0, const unsigned int max_iter, const double precision, const bool verbose) {
  static unsigned int count = 1;
  if (fetestexcept(FE_INVALID)) {
    fprintf(stderr, "ERROR: Invalid argument detected (check for complex or outside domain arguments)\n");
    exit(1);
  }
  if (fetestexcept(FE_OVERFLOW) || fetestexcept(FE_UNDERFLOW)){
    fprintf(stderr, "ERROR: Floating point operations have gone outside of representable range\n");
    fprintf(stderr, "(check if process might be diverging)\n");
    exit(1);
  }
  if (fetestexcept(FE_DIVBYZERO)){
    fprintf(stderr, "ERROR: Division by zero detected\n");
    exit(1);
  }
  if (max_iter && count > max_iter) {
    fprintf(stderr, "Aitken's Δ squared process wasn't able to converge in %i iterations.\n", count);
    printf("Last value : %.*f\n", (int) -floor(log10(precision)),  x0);
    return x0;
  }
  const double x1 = f(x0);
  if (fabs(x0 - x1) < precision) {
    if (verbose)
      printf("Aitken's Δ squared process converged to %.*f\n", (int) -floor(log10(precision)), x1);
    return x1;
  }
  if (verbose)
    printf("Iteration #%i\t : %.*f\n", count, (int) -floor(log10(precision)), x1);
  count++;
  const double x2 = f(x1);
  if (fabs(x1 - x2) < precision) {
    if (verbose)
      printf("Aitken's Δ squared process converged to %.*f\n", (int) -floor(log10(precision)), x2);
    return x2;
  }
  if (verbose)
    printf("Iteration #%i\t : %.*f\n", count, (int) -floor(log10(precision)), x2);
  count++;
  const double xcorr = (x0 * x2 - x1 * x1) / (x0 + x2 - 2 * x1);
  if (fabs(x2 - xcorr) < precision) {
    if (verbose)
      printf("Aitken's Δ squared process converged to %.*f\n", (int) -floor(log10(precision)), xcorr);
    return xcorr;
  }
  if (verbose)
    printf("Iteration #%i\t : %.*f\n", count, (int) -floor(log10(precision)), xcorr);
  count++;
  return aitkens_delta(f, xcorr, max_iter, precision, verbose);
}

// UNIVARIATE NEWTON'S METHOD IMPLEMENTATION
double newton(const univariate_function f, const univariate_function fp, const double x0, const unsigned int max_iter, const double precision, const bool verbose) {
  static unsigned int count = 1;
  if (fetestexcept(FE_INVALID)) {
    fprintf(stderr, "ERROR: Invalid argument detected (check for complex or outside domain arguments)\n");
    exit(1);
  }
  if (fetestexcept(FE_OVERFLOW) || fetestexcept(FE_UNDERFLOW)){
    fprintf(stderr, "ERROR: Floating point operations have gone outside of representable range\n");
    fprintf(stderr, "(check if process might be diverging)\n");
    exit(1);
  }
  if (fetestexcept(FE_DIVBYZERO)){
    fprintf(stderr, "ERROR: Division by zero detected\n");
    exit(1);
  }
  if (max_iter && count > max_iter) {
    fprintf(stderr, "Newton's method wasn't able to converge in %i iterations.\n", count);
    printf("Last value : %.*f", (int) -floor(log10(precision)),  x0);
    return x0;
  }
  const double x1 = x0 - f(x0) / fp(x0);
  if (fabs(x0 - x1) < precision) {
    if (verbose)
      printf("Newton's method converged to %.*f\n", (int) -floor(log10(precision)), x1);
    return x1;
  }
  if (verbose)
    printf("Iteration #%i\t : %.*f\n", count, (int) -floor(log10(precision)), x1);
  count++;
  return newton(f, fp, x1, max_iter, precision, verbose);
}

// SECANT METHOD IMPLEMENTATION
double secant(const univariate_function f, const double x1, const double x0, const unsigned int max_iter, const double precision, const bool verbose) {
  static unsigned int count = 1;
  if (fetestexcept(FE_INVALID)) {
    fprintf(stderr, "ERROR: Invalid argument detected (check for complex or outside domain arguments)\n");
    exit(1);
  }
  if (fetestexcept(FE_OVERFLOW) || fetestexcept(FE_UNDERFLOW)){
    fprintf(stderr, "ERROR: Floating point operations have gone outside of representable range\n");
    fprintf(stderr, "(check if process might be diverging)\n");
    exit(1);
  }
  if (fetestexcept(FE_DIVBYZERO)){
    fprintf(stderr, "ERROR: Division by zero detected\n");
    exit(1);
  }
  if (max_iter && count > max_iter) {
    fprintf(stderr, "Secant method wasn't able to converge in %i iterations.\n", count);
    printf("Last value : %.*f", (int) -floor(log10(precision)),  x0);
    return x0;
  }
  const double x2 = (x0 * f(x1) - x1 * f(x0)) / (f(x1) - f(x0));
  if (fabs(x0 - x1) < precision) {
    if (verbose)
      printf("Secant method converged to %.*f\n", (int) -floor(log10(precision)), x1);
    return x1;
  }
  if (verbose) {
    printf("Iteration #%i\t : %.*f\n", count, (int) -floor(log10(precision)), x1);
  }
  count++;
  return secant(f, x2, x1, max_iter, precision, verbose);
}

// FALSE POSITION METHOD (REGULA FALSI) IMPLEMENTATION
double false_position(const univariate_function f, const double x1, const double x0, const unsigned int max_iter, const double precision, const bool verbose) {
  static unsigned int count = 1;
  if (fetestexcept(FE_INVALID)) {
    fprintf(stderr, "ERROR: Invalid argument detected (check for complex or outside domain arguments)\n");
    exit(1);
  }
  if (fetestexcept(FE_OVERFLOW) || fetestexcept(FE_UNDERFLOW)){
    fprintf(stderr, "ERROR: Floating point operations have gone outside of representable range\n");
    fprintf(stderr, "(check if process might be diverging)\n");
    exit(1);
  }
  if (fetestexcept(FE_DIVBYZERO)){
    fprintf(stderr, "ERROR: Division by zero detected\n");
    exit(1);
  }
  if (max_iter && count > max_iter) {
    fprintf(stderr, "False position method wasn't able to converge in %i iterations.\n", count);
    printf("Last value : %.*f", (int) -floor(log10(precision)),  x0);
    return x0;
  }
  const double x2 = (x0 * f(x1) - x1 * f(x0)) / (f(x1) - f(x0));
  if (fabs(x2 - x1) < precision || fabs(x2 - x0) < precision) {
    if (verbose)
      printf("False position method converged to %.*f\n", (int) -floor(log10(precision)), x2);
    return x2;
  }
  if (verbose) {
    printf("Iteration #%i\t : %.*f\n", count, (int) -floor(log10(precision)), x2);
  }
  count++;
  if (f(x2) * f(x0) < 1)
    return false_position(f, x2, x0, max_iter, precision, verbose);
  else
    return false_position(f, x1, x2, max_iter, precision, verbose);
}

// MÜLLER'S PROCESS IMPLEMENTATION
double muller(const univariate_function f, const double x2, const double x1, const double x0, const unsigned int max_iter, const double precision, const bool verbose) {
  static unsigned int count = 1;
  if (fetestexcept(FE_INVALID)) {
    fprintf(stderr, "ERROR: Invalid argument detected (check for complex or outside domain arguments)\n");
    exit(1);
  }
  if (fetestexcept(FE_OVERFLOW) || fetestexcept(FE_UNDERFLOW)){
    fprintf(stderr, "ERROR: Floating point operations have gone outside of representable range\n");
    fprintf(stderr, "(check if process might be diverging)\n");
    exit(1);
  }
  if (fetestexcept(FE_DIVBYZERO)){
    fprintf(stderr, "ERROR: Division by zero detected\n");
    exit(1);
  }
  if (max_iter && count > max_iter) {
    fprintf(stderr, "Müller's process wasn't able to converge in %i iterations.\n", count);
    printf("Last value : %.*f", (int) -floor(log10(precision)),  x0);
    return x0;
  }
  const double lambda0 = (x2 - x1) / (x1 - x0);
  const double delta = 1 + lambda0;
  const double g = f(x0) * lambda0 * lambda0 - f(x1) * delta * delta + f(x2) * (lambda0 + delta);
  const double lambda1 = (g > 0) ? (-2 * f(x2) * delta) / (g + sqrt(g * g - 4 * f(x2) * delta * lambda0 * (f(x0)  * lambda0 - f(x1) * delta + f(x2)))) : (-2 * f(x2) * delta) / (g - sqrt(g * g - 4 * f(x2) * delta * lambda0 * (f(x0)  * lambda0 - f(x1) * delta + f(x2))));
  const double x3 = x2 + lambda1 * (x2 - x1);
  if (fabs(x2 - x3) < precision) {
    if (verbose) {
      printf("Müller's process converged to %.*f\n", (int) -floor(log10(precision)), x3);
    }
    return x3;
  }
  if (verbose) {
    printf("Iteration #%i\t : %.*f\n", count, (int) -floor(log10(precision)), x3);
  }
  count++;
  return muller(f, x3, x2, x1, max_iter, precision, verbose);
}

