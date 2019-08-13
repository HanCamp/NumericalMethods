//
//  multisolvers.c
//  NummetC
//
//  Created by Haniel Campos Alcantara Paulo on 8/12/19.
//  Copyright © 2019 Haniel Campos. All rights reserved.
//

#include <math.h>
#include <string.h>
#include <fenv.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include "multisolvers.h"
#include "definitions.h"

void linear_iteration_multi(const multivariate_function f, double* x0, double* tmp, const unsigned int dimension, const unsigned int max_iter, const double precision, const bool verbose) {
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
    fprintf(stderr, "Linear iteration wasn't able to converge in %i iterations.", max_iter);
    if (verbose) {
      printf("Last value : ");
      print_vector(x0, dimension, precision);
    }
  } else {
    f(x0, tmp);
    if (difference_norm(x0, tmp, dimension) < precision) {
      if (verbose) {
        printf("Linear iteration converged to ");
        print_vector(tmp, dimension, precision);
      }
      memmove(x0, tmp, dimension * sizeof(double));
    } else {
      if (verbose) {
        printf("Iteration #%i\t : ", count);
        print_vector(tmp, dimension, precision);
      }
      memmove(x0, tmp, dimension * sizeof(double));
      count++;
      linear_iteration_multi(f, x0, tmp, dimension, max_iter, precision, verbose);
    }
  }
}

void aitkens_delta_multi(const multivariate_function f, double* x0, double* tmp1, double* tmp2, const unsigned int dimension, const unsigned int max_iter, const double precision, const bool verbose) {
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
    fprintf(stderr, "Aitken's Δ squared process wasn't able to converge in %i iterations.", max_iter);
    if (verbose) {
      printf("Last value : ");
      print_vector(x0, dimension, precision);
    }
  } else {
    f(x0, tmp1);
    if (difference_norm(x0, tmp1, dimension) < precision) {
      if (verbose) {
        printf("Aitken's Δ squared process converged to ");
        print_vector(tmp1, dimension, precision);
      }
      memmove(x0, tmp1, dimension * sizeof(double));
    } else {
      if (verbose) {
        printf("Iteration #%i\t : ", count);
        print_vector(tmp1, dimension, precision);
      }
      count++;
      f(tmp1, tmp2);
      if (difference_norm(tmp1, tmp2, dimension) < precision) {
        if (verbose) {
          printf("Aitken's Δ squared process converged to ");
          print_vector(tmp2, dimension, precision);
        }
        memmove(x0, tmp2, dimension * sizeof(double));
      } else {
        if (verbose) {
          printf("Iteration #%i\t : ", count);
          print_vector(tmp2, dimension, precision);
        }
        count++;
        for (unsigned int i = 0; i < dimension; i++)
          x0[i] = (x0[i] * tmp2[i] - tmp1[i] * tmp1[i]) / (x0[i] + tmp2[i] - 2 * tmp1[i]);
        if (difference_norm(tmp2, x0, dimension) < precision) {
          if (verbose) {
            printf("Aitken's Δ squared process converged to ");
            print_vector(x0, dimension, precision);
          }
        } else {
          aitkens_delta_multi(f, x0, tmp1, tmp2, dimension, max_iter, precision, verbose);
        }
      }
    }
  }
}
