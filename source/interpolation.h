//
//  interpolation.h
//  NummetC
//
//  Created by Haniel Campos Alcantara Paulo on 8/13/19.
//  Copyright © 2019 Haniel Campos. All rights reserved.
//

#ifndef interpolation_h
#define interpolation_h

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fenv.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include "definitions.h"

double lagrange(const double xinput, const double* xval, const double* fval, const unsigned int npoints);

double aitken(const double xinput, const double* xval, const double* fval, const unsigned int npoints);

double newton_ascending(const struct DifferenceTable* table, const double xinput, const double* xval, const double* fval, const unsigned int degree, const double h_width, const  unsigned int npoints);

double newton_ascending2(const double xinput, const double* xval, const double* fval, const unsigned int degree, const double h_width, const  unsigned int npoints);

double newton_descending(const struct DifferenceTable* table, const double xinput, const unsigned int index, const double* xval, const double* fval, const unsigned int degree, const double h_width, const unsigned int npoints);

double newton_descending2(const double xinput, const unsigned int index, const double* xval, const double* fval, const unsigned int degree, const double h_width, const unsigned int npoints);

double stirling(const struct DifferenceTable* differenceTable, const double xinput, const unsigned int index, const double* xval, const double* fval, const unsigned int degree, const double h_width, const unsigned int npoints);

double stirling2(const double xinput, const unsigned int index, const double* xval, const double* fval, const unsigned int degree, const double h_width, const unsigned int npoints);

double everett(const struct DifferenceTable* differenceTable, const double xinput, const unsigned int index, const double* xval, const double* fval, const unsigned int degree, const double h_width, const unsigned int npoints);

double everett2(const double xinput, const unsigned int index, const double* xval, const double* fval, const unsigned int degree, const double h_width, const unsigned int npoints);

#endif /* interpolation_h */
