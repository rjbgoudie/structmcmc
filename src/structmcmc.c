/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1999-2012   The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

/* Integer vector tabulation */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <R_ext/Applic.h>

/* read-only, apart from 'ans' */
void R_tabulate(int *x, int *n, int *nbin, int *ans)
{
  int i;
  if(*n < 1) return;
  for(i = 0 ; i < *n ; i++)
    if(x[i] != R_NaInt && x[i] > 0 && x[i] <= *nbin)
      ans[x[i] - 1]++;
}


/*
   1-bit boolean vectors for R
   first bit is stored in lowest (rightmost) bit of forst word
   remember that rightshifting is dangerous because we use the sign position
   Copyright 2008 Jens Oehlschl‰gel
*/

/*
 & bitwise and
 | bitwise or
 ^ bitwise xor
 ~ bitwise not
*/

void bit_or(int *b1, int *b2, int *ret, int n){
  register int i;
  for (i=0; i<n; i++){
    ret[i] = b1[i] | b2[i];
  }
}

SEXP R_bit_or(SEXP b1_, SEXP b2_, SEXP ret_){
  SEXP ret_copy_;
  PROTECT(ret_copy_ = duplicate(ret_));
  int *b1 = INTEGER(b1_);
  int *b2 = INTEGER(b2_);
  int *ret = INTEGER(ret_copy_);
  int n = LENGTH(b1_);
  bit_or(b1, b2, ret, n);
  UNPROTECT(1);
  return(ret_copy_);
}

SEXP R_bit_or_list(SEXP rowsNotAllowed_, SEXP ret_){
  SEXP ret_copy_;
  PROTECT(ret_copy_ = duplicate(ret_));
  int *ret = INTEGER(ret_copy_);
  int n = LENGTH(ret_copy_);
  int m = length(rowsNotAllowed_);

  int i;
  int j;

  // loop over indicators
  for (i = 0; i < n; i++){
    SEXP list_item = VECTOR_ELT(rowsNotAllowed_, 0);
    int now = INTEGER(list_item)[i];

    // loop over list items of rowsNotAllowed_
    for (j = 1; j < m; j++){
      SEXP list_item = VECTOR_ELT(rowsNotAllowed_, j);
      int now_this = INTEGER(list_item)[i];

      now = now | now_this;
    }
    ret[i] = now;
  }
  UNPROTECT(1);
  return(ret_copy_);
}




SEXP R_bit_and_list(SEXP rowsNotAllowed_, SEXP ret_){
  SEXP ret_copy_;
  PROTECT(ret_copy_ = duplicate(ret_));
  int *ret = INTEGER(ret_copy_);
  int n = LENGTH(ret_copy_);
  int m = length(rowsNotAllowed_);

  int i;
  int j;

  // loop over indicators
  for (i = 0; i < n; i++){
    SEXP list_item = VECTOR_ELT(rowsNotAllowed_, 0);
    int now = INTEGER(list_item)[i];

    // loop over list items of rowsNotAllowed_
    for (j = 1; j < m; j++){
      SEXP list_item = VECTOR_ELT(rowsNotAllowed_, j);
      int now_this = INTEGER(list_item)[i];

      now = now & now_this;
    }
    ret[i] = now;
  }
  UNPROTECT(1);
  return(ret_copy_);
}

void bit_setdiff(int *b1, int *b2, int *ret, int n){
  register int i;
  for (i=0; i<n; i++){
    ret[i] = b1[i] & (~b2[i]);
  }
}

SEXP R_bit_setdiff(SEXP b1_, SEXP b2_, SEXP ret_){
  SEXP ret_copy_;
  PROTECT(ret_copy_ = duplicate(ret_));
  int *b1 = INTEGER(b1_);
  int *b2 = INTEGER(b2_);
  int *ret = INTEGER(ret_copy_);
  int n = LENGTH(b1_);
  bit_setdiff(b1, b2, ret, n);
  UNPROTECT(1);
  return(ret_copy_);
}



#include <Rdefines.h>
#include <Rmath.h>



// Based upon code by Henrik Bengtsson from matrixStats
// Artistic-2.0 License

/*
 logSumExp_double(x):

  1. Scans for the maximum value of x=(x[0], x[1], ..., x[n-1])
  2. Computes result from 'x'.

  NOTE: The above sweeps the "contiguous" 'x' vector twice.
*/
double logSumExp_double(double *x, int n) {
  int ii, iMax;
  double xii, xMax, sum;

  /* Quick return? */
  if (n == 0) {
    return(R_NegInf);
  } else if (n == 1) {
    return(x[0]);
  }

  /* Find the maximum value */
  iMax = 0;
  xMax = x[0];
  for (ii=1; ii < n; ii++) {
    /* Get the ii:th value */
    xii = x[ii];

    if (ISNAN(xii)) {
      return(R_NaReal);
    }

    if (xii > xMax) {
      iMax = ii;
      xMax = xii;
    }

    if (ii % 1000000 == 0) {
      R_CheckUserInterrupt();
    }
  } /* for (ii ...) */

  /* Sum differences */
  sum = 0.0;
  for (ii=0; ii < n; ii++) {
    if (ii == iMax) {
      continue;
    }

    /* Get the ii:th value */
    xii = x[ii];

    if (ISNAN(xii)) {
      return(R_NaReal);
    } else {
      sum += exp(xii - xMax);
    }

    if (ii % 1000000 == 0) {
      R_CheckUserInterrupt();
    }
  } /* for (ii ...) */

  sum = xMax + log1p(sum);

  return(sum);
} /* logSumExp_double() */


SEXP logSumExp(SEXP lx) {
  double *x;
  int n;

  /* Argument 'lx': */
  if (!isReal(lx)) {
    error("Argument 'lx' must be a numeric vector.");
  }

  /* Get the values */
  x = REAL(lx);
  n = length(lx);

  return(Rf_ScalarReal(logSumExp_double(x, n)));
} /* logSumExp() */
