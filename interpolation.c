//
//  interpolation.c
//  NummetC
//
//  Created by Haniel Campos Alcantara Paulo on 8/13/19.
//  Copyright © 2019 Haniel Campos. All rights reserved.
//

#include <math.h>
#include <string.h>
#include <fenv.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include "definitions.h"
#include "interpolation.h"

double lagrange(const double xinput, const double* xval, const double* fval, const unsigned int npoints) {
  double pval = 0;
  for (unsigned int i = 0; i < npoints; i++) {
    double L_i = 1;
    for (unsigned int j = 0; j < npoints; j++) {
      if (j != i)
        L_i *= (xinput - xval[j]) / (xval[i] - xval[j]);
    }
    pval += fval[i] * L_i;
  }
  return pval;
}

double aitken_helper(const double xinput, const double* xval, const double* fval, const unsigned int count, const unsigned int npoints)  {
  if (npoints - count > 1) {
    double polyVals[(npoints - count - 1) * sizeof(double)];
    for (unsigned int i = 0; i < npoints - count - 1; i++)
      polyVals[i] = ((xval[count + i + 1] - xinput) * fval[0] - (xval[count] - xinput) * fval[i + 1]) / (xval[count + i + 1] - xval[count]);
    return aitken_helper(xinput, xval, polyVals, count + 1, npoints);
  } else {
    return fval[0];
  }
  return 0;
}

double aitken(const double xinput, const double* xval, const double* fval, const unsigned int npoints) {
  return aitken_helper(xinput, xval, fval, 0, npoints);
}

double newton_ascending2(const double xinput, const double* xval, const double* fval, const unsigned int degree, const double h_width, const  unsigned int npoints) {
  if (degree > npoints - 1) {
    fprintf(stderr, "ERROR:Invalid degree provided\n");
    exit(1);
  }
  struct DifferenceTable* differenceTable = new_difference_table(fval, npoints);
  double s = (xinput - xval[0]) / h_width;
  double polyVal = fval[0];
  for (unsigned int i = 1; i <= degree; i++)
    polyVal += difference(differenceTable, 0, i) * binomial(s, i);
  free(differenceTable);
  return polyVal;
}

double newton_ascending(const struct DifferenceTable* differenceTable, const double xinput, const double* xval, const double* fval, const unsigned int degree, const double h_width, const  unsigned int npoints) {
  if (degree > npoints - 1) {
    fprintf(stderr, "ERROR:Invalid degree provided\n");
    exit(1);
  }
  double s = (xinput - xval[0]) / h_width;
  double polyVal = fval[0];
  for (unsigned int i = 1; i <= degree; i++)
    polyVal += difference(differenceTable, 0, i) * binomial(s, i);
  return polyVal;
}

double newton_descending2(const double xinput, const unsigned int index, const double* xval, const double* fval, const unsigned int degree, const double h_width, const unsigned int npoints) {
  if (degree > npoints - 1) {
    fprintf(stderr, "ERROR: Invalid degree provided\n");
    exit(1);
  }
  struct DifferenceTable* differenceTable = new_difference_table(fval, npoints);
  double s = (xinput - xval[index]) / h_width;
  double polyVal = fval[index];
  for (unsigned int i = 1; i <= degree; i++)
    polyVal += difference(differenceTable, index - i, i) * binomial(s + i - 1, i);
  free(differenceTable);
  return polyVal;
}

double newton_descending(const struct DifferenceTable* differenceTable, const double xinput, const unsigned int index, const double* xval, const double* fval, const unsigned int degree, const double h_width, const unsigned int npoints) {
  if (degree - index > npoints - 1) {
    fprintf(stderr, "ERROR: Invalid degree provided\n");
    exit(1);
  }
  double s = (xinput - xval[index]) / h_width;
  double polyVal = fval[index];
  for (unsigned int i = 1; i <= degree; i++)
    polyVal += difference(differenceTable, index - i, i) * binomial(s + i - 1, i);
  return polyVal;
}

double stirling2(const double xinput, const unsigned int index, const double* xval, const double* fval, const unsigned int degree, const double h_width, const unsigned int npoints) {
  if (degree > npoints - 1) {
    fprintf(stderr, "ERROR: Invalid degree provided\n");
    exit(1);
  }
  struct DifferenceTable* differenceTable = new_difference_table(fval, npoints);
  double s = (xinput - xval[index]) / h_width;
  double polyVal = fval[index];
  for (unsigned int i = 1; i <= degree; i++) {
    if (!(i % 2)) {
      polyVal += 0.5 * (binomial(s + i / 2 - 1, i) + binomial(s + i / 2, i)) * difference(differenceTable, index - i / 2, i);
    } else {
      polyVal += 0.5 * binomial(s + (i - 1) / 2, i) * (difference(differenceTable, index - (i - 1) / 2, i) + difference(differenceTable, index - 1 - (i - 1) / 2, i));
    }
  }
  free(differenceTable);
  return polyVal;
}

double stirling(const struct DifferenceTable* differenceTable, const double xinput, const unsigned int index, const double* xval, const double* fval, const unsigned int degree, const double h_width, const unsigned int npoints) {
  if (degree + index > npoints - 1) {
    fprintf(stderr, "ERROR: Invalid degree provided\n");
    exit(1);
  }
  double s = (xinput - xval[index]) / h_width;
  double polyVal = fval[index];
  for (unsigned int i = 1; i <= degree; i++) {
    if (!(i % 2)) {
      polyVal += 0.5 * (binomial(s + i / 2 - 1, i) + binomial(s + i / 2, i)) * difference(differenceTable, index - i / 2, i);
    } else {
      polyVal += 0.5 * binomial(s + (i - 1) / 2, i) * (difference(differenceTable, index - (i - 1) / 2, i) + difference(differenceTable, index - 1 - (i - 1) / 2, i));
    }
  }
  return polyVal;
}

double everett2(const double xinput, const unsigned int index, const double* xval, const double* fval, const unsigned int degree, const double h_width, const unsigned int npoints) {
  if (degree + index > npoints - 1) {
    fprintf(stderr, "Invalid degree provided\n");
    exit(1);
  }
  struct DifferenceTable* differenceTable = new_difference_table(fval, npoints);
  double s = (xinput - xval[index]) / h_width;
  double p = 1 - s;
  double polyVal = binomial(p, 1) * fval[index] + binomial(s, 1) * fval[index + 1];
  for (unsigned int i = 1; i <= degree / 2; i++) {
    polyVal += binomial(p + i, 2 * i + 1) * central_difference(differenceTable, index, 2 * i);
    polyVal += binomial(s + i, 2 * i + 1) * central_difference(differenceTable, index + 1, 2 * i);
  }
  free(differenceTable);
  return polyVal;
}

double everett(const struct DifferenceTable* differenceTable, const double xinput, const unsigned int index, const double* xval, const double* fval, const unsigned int degree, const double h_width, const unsigned int npoints) {
  if (degree > npoints - 1) {
    fprintf(stderr, "Invalid degree provided\n");
    exit(1);
  }
  double s = (xinput - xval[index]) / h_width;
  double p = 1 - s;
  double polyVal = binomial(p, 1) * fval[index] + binomial(s, 1) * fval[index + 1];
  for (unsigned int i = 1; i <= degree / 2; i++) {
    polyVal += binomial(p + i, 2 * i + 1) * central_difference(differenceTable, index, 2 * i);
    polyVal += binomial(s + i, 2 * i + 1) * central_difference(differenceTable, index + 1, 2 * i);
  }
  return polyVal;
}
