/* 
 * Functions to compute dissimilarity between two matrices
 * X and Y.
 *
 * Based on code from vegdist by Jari Oksanen:
 *
 * (C) 2001-2009, Jari Oksanen
 * (C) 2009-2012 Gavin L. Simpson
 *
 * Licene: GPL 2
 */

/* Standard R headers */

#include <R.h>
#include <Rmath.h>
#include <math.h>

/* Indices */
/* Note we don't actually call all of these via xx_distance
 * Some are called via direct methods, but we include the 
 * indices here to allow the pattern matching to work 
 * correctly
 */
#define EUCLIDEAN 1
#define SQEUCLIDEAN 2
#define CHORD 3
#define SQCHORD 4
#define BRAY 5
#define CHISQUARE 6
#define SQCHISQUARE 7
#define INFORMATION 8
#define CHIDISTANCE 9
#define MANHATTAN 10
#define KENDALL 11
#define GOWER 12
#define ALTGOWER 13
#define MIXED 14

/* Distance functions */

/* Euclidean distance */
double c_xx_euclidean(double *x, int nr, int nc, int i1, int i2)
{
    double dist, dev;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dev = x[i1] - x[i2];
	    dist += dev*dev;
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if (count == 0) return NA_REAL;
    return sqrt(dist);
}

/* Squared Euclidean distance */
double c_xx_sq_euclidean(double *x, int nr, int nc, int i1, int i2)
{
    double dist, dev;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dev = x[i1] - x[i2];
	    dist += dev*dev;
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if (count == 0) return NA_REAL;
    return dist;
}

/* Chord distance */
double c_xx_chord(double *x, int nr, int nc, int i1, int i2)
{
    double dist, dev;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dev = sqrt(x[i1]) - sqrt(x[i2]);
	    dist += dev*dev;
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if (count == 0) return NA_REAL;
    return sqrt(dist);
}

/* Squared Chord distance */
double c_xx_sq_chord(double *x, int nr, int nc, int i1, int i2)
{
    double dist, dev;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dev = sqrt(x[i1]) - sqrt(x[i2]);
	    dist += dev*dev;
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if (count == 0) return NA_REAL;
    return dist;
}

/*  Bray-Curtis */
double c_xx_bray(double *x, int nr, int nc, int i1, int i2)
{
    double dist, total;
    int count, j;
    
    total = 0.0;
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dist += fabs(x[i1] - x[i2]);
	    total += x[i1] + x[i2];
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if (count==0) return NA_REAL;
    dist /= total;
    return dist;
}

/*  chi square */
double c_xx_chi_square(double *x, int nr,  int nc, int i1, int i2)
{
    double dist, dev, nomin, denomin;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    if(x[i1] != 0 || x[i2] != 0) {
		dev = x[i1] - x[i2];
		nomin = dev*dev;
		denomin = x[i1] + x[i2];
		dist += nomin / denomin;
		count++;
	    }
	}
	i1 += nr;
	i2 += nr;
    }
    if (count==0) return NA_REAL;
    return sqrt(dist);
}

/*  square chi square */
double c_xx_sq_chi_square(double *x, int nr, int nc, int i1, int i2)
{
    double dist, dev, nomin, denomin;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    if(x[i1] != 0 || x[i2] != 0) {
		dev = x[i1] - x[i2];
		nomin = dev*dev;
		denomin = x[i1] + x[i2];
		dist += nomin / denomin;
		count++;
	    }
	}
	i1 += nr;
	i2 += nr;
    }
    if (count==0) return NA_REAL;
    return dist;
}

/* information statistic */
double c_xx_information(double *x, int nr, int nc, int i1, int i2)
{
    double dist, XY, A, B, Adist, Bdist;
    int count, j;
    
    count = 0;
    dist = 0.0;
    Adist = 0.0;
    Bdist = 0.0;
    A = 0.0;
    B = 0.0;
    for(j=0; j<nc; j++) {
	if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    XY = x[i1] + x[i2];
	    A += x[i1] * (log((2 * x[i1]) / XY)/log(2));
	    B += x[i2] * (log((2 * x[i2]) / XY)/log(2));
	    if(R_FINITE(A)) {
		Adist += A;
	    }
	    if(R_FINITE(B)) {
		Bdist += B;
	    }
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if(count==0) return NA_REAL;
    dist = A + B;
    return dist;
}

/*  chi square distance*/
/* currently not correct */
double c_xx_chi_distance(double *x, int nr, int nc, int i1, int i2)
{
    double dist, dev, nomin;
    int count, j;
    
    count = 0;
    dist = 0.0;
    nomin = 0.0;

    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dev = x[i1] - x[i2];
	    nomin += dev*dev;
	}
	i1 += nr;
	i2 += nr;
    }
    if (count==0) return NA_REAL;
    return dist;
}

/* Manhattan metric */
double c_xx_manhattan(double *x, int nr, int nc, int i1, int i2)
{
    double dist;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dist += fabs( x[i1] - x[i2] );
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if (count == 0) return NA_REAL;
    return dist;
}

/* Gower's distance */
/* Needs to be preprocessed by dividing by Maxi - Mini
   in the R wrapper */
double c_xx_gower(double *x, int nr, int nc, int i1, int i2)
{
    double dist;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dist += fabs( x[i1] - x[i2] );
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if (count == 0) return NA_REAL;
    return dist;
}

/* Alternative Gower's distance */
/* Needs to be preprocessed by dividing by Maxi - Mini
   in the R wrapper */
double c_xx_alt_gower(double *x, int nr, int nc, int i1, int i2)
{
    double dist;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dist += fabs( x[i1] - x[i2] );
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if (count == 0) return NA_REAL;
    return sqrt(2 * dist);
}

/* Driver */

static double (*distfun)(double*, int, int, int, int);

void c_xx_distance(double *x, int *nr, int *nc, double *d, 
		 int *diag, int *method)
{
    int dc, i, j, ij;
    switch(*method) {
    case EUCLIDEAN:
	distfun = c_xx_euclidean;
	break;
    case SQEUCLIDEAN:
	distfun = c_xx_sq_euclidean;
	break;
    case CHORD:
	distfun = c_xx_chord;
	break;
    case SQCHORD:
	distfun = c_xx_sq_chord;
	break;
    case BRAY:
	distfun = c_xx_bray;
	break;
    case CHISQUARE:
	distfun = c_xx_chi_square;
	break;
    case SQCHISQUARE:
	distfun = c_xx_sq_chi_square;
	break;
    case INFORMATION:
	distfun = c_xx_information;
	break;
    case MANHATTAN:
	distfun = c_xx_manhattan;
	break;
    case GOWER:
	distfun = c_xx_gower;
	break;
    case ALTGOWER:
	distfun = c_xx_alt_gower;
	break;
    default:
	error("Unknown distance in the internal C function");
    }
    
    dc = (*diag) ? 0 : 1;
    ij = 0;
    for (j=0; j <= *nr; j++)
	for (i=j+dc; i < *nr; i++) {
	    d[ij++] = distfun(x, *nr, *nc, i, j);
	}
}

/* .Call interface */

/* Based on code from distance.c in R Sources*/

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998-2012  The R Core Team
 *  Copyright (C) 2002, 2004  The R Foundation
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

#include <Rinternals.h>

SEXP Cdistxx(SEXP x, SEXP smethod)
{
    SEXP ans;
    int nr = nrows(x), nc = ncols(x), method = asInteger(smethod);
    int N, diag = 0;
    N = (double)nr * (nr-1)/2; /* avoid overflow for N ~ 50,000 */
    PROTECT(ans = allocVector(REALSXP, N));
    c_xx_distance(REAL(x), &nr, &nc, REAL(ans), &diag, &method);
    UNPROTECT(1);
    return ans;
}

/* Kendall's coefficient -------------------------------------------
 *
 * Should be called separately from the underlying R code,
 * not via xy_distance.
 *
 * maxi: the max abundance for each species
 *
*/

/* inner-most loop */
double c_xx_KENDALL(double *x, int nr, int nc, int i1, int i2,
		      double *maxi)
{
  double dist, dev;
  int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dev = (x[i1] >= x[i2]) ? x[i2] : x[i1];
	    dist += maxi[j] - dev;
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if (count == 0) return NA_REAL;
    return dist;
}

/* driver, called by .call wrapper Ckendallxx */
void c_xx_kendall(double *x, int *nr, int *nc, double *d, 
		    int *diag, double *maxi)
{
int dc, i, j, ij;

dc = (*diag) ? 0 : 1;
ij = 0;
    for(j=0; j <= *nr; j++) {
	for(i=j+dc; i < *nr; i++) {
	    d[ij++] = c_xx_KENDALL(x, *nr, *nc, i, j, maxi);
	}
    }
}

/* .call wrapper  */
SEXP Ckendallxx(SEXP x, SEXP smaxi)
{
  SEXP ans;
  int nr = nrows(x), nc = ncols(x);
  int N, diag = 0;
  double maxi = asReal(smaxi);
  N = (double)nr * (nr-1)/2; /* avoid overflow for N ~ 50,000 */
  PROTECT(ans = allocVector(REALSXP, N));
  c_xx_kendall(REAL(x), &nr, &nc, REAL(ans), &diag, &maxi);
  UNPROTECT(1);
  return ans;
}

/* Chi square distance ------------------------------------------------
 *
 * Should be called separately from the underlying R code,
 * not via xy_distance.
 *
 * csum: species sums
 *
 */

double c_xx_CHISQ_DIST(double *x, int nr, int nc, int i1, int i2,
			 double *csum, double ccsum)
{
    double dist, dev, denom, nomin;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    dev = x[i1] - x[i2];
	    nomin = dev*dev;
	    denom = csum[j] / ccsum;
	    dist += nomin / denom;
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if (count == 0) return NA_REAL;
    return sqrt(dist);
}

void c_xx_chisq_dist(double *x, int *nr, int *nc, double *d, 
		       int *diag,  double *csum)
{
  int dc, i, j, k, ij;
  double ccsum;
  
  ccsum = 0.0;
  
  ij = 0;
  dc = (*diag) ? 0 : 1;
  for(k=0; k < *nc; k++) {
    ccsum += csum[k];
  }
  
  for(j=0; j < *nr; j++) {
    for(i=j+dc; i < *nr; i++) {
      d[ij++] = c_xx_CHISQ_DIST(x, *nr, *nc, i, j, 
				csum, ccsum);
    }
  }
}

/* .call wrapper  */
SEXP Cchisqdistxx(SEXP x, SEXP scsum)
{
  SEXP ans;
  int nr = nrows(x), nc = ncols(x);
  int N, diag = 0;
  double csum = asReal(scsum);
  N = (double)nr * (nr-1)/2; /* avoid overflow for N ~ 50,000 */
  PROTECT(ans = allocVector(REALSXP, N));
  c_xx_chisq_dist(REAL(x), &nr, &nc, REAL(ans), &diag, &csum);
  UNPROTECT(1);
  return ans;
}


/* Gower's coefficient for mixed data ---------------------------------
 *
 * Should be called separately from the underlying R code,
 * not via xy_distance.
 *
 * vtype  : variable type
 *          1 == Symmetric Binary
 *          2 == Asymmetric Binary
 *          3 == Nominal (class/factor)
 *          4 == Ordinal (ordered factor)
 *          5 == Quantitative
 * weights: variable weights
 * R      : variable range (max - min)
 *
 */
double c_xx_MIXED(double *x, int nr, int nc, int i1, int i2, 
		  int *vtype, double *weights, double *R, 
		  double wsum)
{
  double dist, dev;
  int count, j;
  
  count = 0;
  dist = 0.0;
    wsum = 0.0;
    
    for (j=0; j<nc; j++) {
      if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	if(vtype[j] == 1) {
	  dev = (x[i1] == x[i2]) ? 1 : 0;
	  dist += dev * weights[j];
	}
	if(vtype[j] == 2) { // Asymmetric binary
	  /*dev = (x[i1] == x[i2]) ? 1 : 0;
	    dist += dev * weights[j]; */
	  if((x[i1] != 0) || (x[i2] != 0)) {
	    // both x1 and x2 not zero for this variables
	    dev = (x[i1] == x[i2]) ? 1 : 0;
	    dist += dev * weights[j];
	  } else {
	    /* set jth current weight to zero and do not
	       increment dist as ignoring double zero
	       We actually subtract the weight as it gets added
	       later on.
	    */
	    wsum -= weights[j];
	  }
	}
	if(vtype[j] == 3) { // Nominal
	  dev = (x[i1] == x[i2]) ? 1 : 0;
	  dist += dev * weights[j];
	}
	if(vtype[j] == 4) { // Ordinal
	  /* ordinal data currently handled as Nominal */
	  dev = (x[i1] == x[i2]) ? 1 : 0;
	  dist += dev * weights[j];
	  break;
	}
	if(vtype[j] == 5) {
	  dev = 1 - (fabs(x[i1] - x[i2]) / R[j]);
	  dist += dev * weights[j];
	}
	count++;
	// only summing weights for non-NA comparisons
	wsum += weights[j];
      }
      i1 += nr;
      i2 += nr;
    }
    if (count == 0) return NA_REAL;
    return 1 - (dist / wsum);
}

void c_xx_mixed(double *x, int *nr, int *nc, double *d, 
		int *diag, int *vtype, double *weights, double *R)
{
  int dc, i, j, k, ij;
  double wsum;
  
  wsum = 0.0;
  
  ij = 0;
  
  dc = (*diag) ? 0 : 1;
  
  for(k=0; k <*nc; k++) {
    wsum += weights[k];
  }
  
  for(j=0; j < *nr; j++) {
    for(i=j+dc; i < *nr; i++) {
      d[ij++] = c_xx_MIXED(x, *nr, *nc, i, j, vtype,
			   weights, R, wsum);
    }
  }
}

/* .call wrapper  */
SEXP Cmixedxx(SEXP x, SEXP svtype, SEXP sweights, SEXP sR)
{
  SEXP ans;
  int nr = nrows(x), nc = ncols(x), vtype = asInteger(svtype);
  int N, diag = 0;
  double weights = asReal(sweights), R = asReal(sR);
  N = (double)nr * (nr-1)/2; /* avoid overflow for N ~ 50,000 */
  PROTECT(ans = allocVector(REALSXP, N));
  c_xx_mixed(REAL(x), &nr, &nc, REAL(ans), &diag, &vtype, &weights, &R);
  UNPROTECT(1);
  return ans;
}
