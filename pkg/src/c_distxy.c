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
/* Note we don't actually call all of these via xy_distance
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
double c_xy_euclidean(double *x, double *y, int nr1, int nr2, 
		    int nc, int i1, int i2)
{
    double dist, dev;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(y[i2])) {
	    dev = x[i1] - y[i2];
	    dist += dev*dev;
	    count++;
	}
	i1 += nr1;
	i2 += nr2;
    }
    if (count == 0) return NA_REAL;
    return sqrt(dist);
}

/* Squared Euclidean distance */
double c_xy_sq_euclidean(double *x, double *y, int nr1, int nr2, 
		       int nc, int i1, int i2)
{
    double dist, dev;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(y[i2])) {
	    dev = x[i1] - y[i2];
	    dist += dev*dev;
	    count++;
	}
	i1 += nr1;
	i2 += nr2;
    }
    if (count == 0) return NA_REAL;
    return dist;
}

/* Chord distance */
double c_xy_chord(double *x, double *y, int nr1, int nr2, 
		int nc, int i1, int i2)
{
    double dist, dev;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(y[i2])) {
	    dev = sqrt(x[i1]) - sqrt(y[i2]);
	    dist += dev*dev;
	    count++;
	}
	i1 += nr1;
	i2 += nr2;
    }
    if (count == 0) return NA_REAL;
    return sqrt(dist);
}

/* Squared Chord distance */
double c_xy_sq_chord(double *x, double *y, int nr1, int nr2, 
		   int nc, int i1, int i2)
{
    double dist, dev;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(y[i2])) {
	    dev = sqrt(x[i1]) - sqrt(y[i2]);
	    dist += dev*dev;
	    count++;
	}
	i1 += nr1;
	i2 += nr2;
    }
    if (count == 0) return NA_REAL;
    return dist;
}

/*  Bray-Curtis */
double c_xy_bray(double *x, double *y, int nr1, int nr2, 
       	       int nc, int i1, int i2)
{
    double dist, total;
    int count, j;
    
    total = 0.0;
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(y[i2])) {
	    dist += fabs(x[i1] - y[i2]);
	    total += x[i1] + y[i2];
	    count++;
	}
	i1 += nr1;
	i2 += nr2;
    }
    if (count==0) return NA_REAL;
    dist /= total;
    return dist;
}

/*  chi square */
double c_xy_chi_square(double *x, double *y, int nr1, int nr2, 
		     int nc, int i1, int i2)
{
    double dist, dev, nomin, denomin;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(y[i2])) {
	    if(x[i1] != 0 || y[i2] != 0) {
		dev = x[i1] - y[i2];
		nomin = dev*dev;
		denomin = x[i1] + y[i2];
		dist += nomin / denomin;
		count++;
	    }
	}
	i1 += nr1;
	i2 += nr2;
    }
    if (count==0) return NA_REAL;
    return sqrt(dist);
}

/*  square chi square */
double c_xy_sq_chi_square(double *x, double *y, int nr1, int nr2, 
			int nc, int i1, int i2)
{
    double dist, dev, nomin, denomin;
    int count, j;
	
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(y[i2])) {
	    if(x[i1] != 0 || y[i2] != 0) {
		dev = x[i1] - y[i2];
		nomin = dev*dev;
		denomin = x[i1] + y[i2];
		dist += nomin / denomin;
		count++;
	    }
	}
	i1 += nr1;
	i2 += nr2;
    }
    if (count==0) return NA_REAL;
    return dist;
}

/* information statistic */
double c_xy_information(double *x, double *y, int nr1, int nr2, 
		      int nc, int i1, int i2)
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
	if(R_FINITE(x[i1]) && R_FINITE(y[i2])) {
	    XY = x[i1] + y[i2];
	    A += x[i1] * (log((2 * x[i1]) / XY)/log(2));
	    B += y[i2] * (log((2 * y[i2]) / XY)/log(2));
	    if(R_FINITE(A)) {
		Adist += A;
	    }
	    if(R_FINITE(B)) {
		Bdist += B;
	    }
	    count++;
	}
	i1 += nr1;
	i2 += nr2;
    }
    if(count==0) return NA_REAL;
    dist = A + B;
    return dist;
}

/* Manhattan metric */
double c_xy_manhattan(double *x, double *y, int nr1, int nr2, 
		    int nc, int i1, int i2)
{
    double dist;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(y[i2])) {
	    dist += fabs( x[i1] - y[i2] );
	    count++;
	}
	i1 += nr1;
	i2 += nr2;
    }
    if (count == 0) return NA_REAL;
    return dist;
}

/* Gower's distance */
/* Needs to be preprocessed by dividing by Maxi - Mini
   in the R wrapper */
double c_xy_gower(double *x, double *y, int nr1, int nr2, 
		int nc, int i1, int i2)
{
    double dist;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(y[i2])) {
	    dist += fabs( x[i1] - y[i2] );
	    count++;
	}
	i1 += nr1;
	i2 += nr2;
    }
    if (count == 0) return NA_REAL;
    return dist;
}

/* Alternative Gower's distance */
/* Needs to be preprocessed by dividing by Maxi - Mini
   in the R wrapper */
double c_xy_alt_gower(double *x, double *y, int nr1, int nr2, 
		    int nc, int i1, int i2)
{
    double dist;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(y[i2])) {
	    dist += fabs( x[i1] - y[i2] );
	    count++;
	}
	i1 += nr1;
	i2 += nr2;
    }
    if (count == 0) return NA_REAL;
    return sqrt(2 * dist);
}

/* Driver */

void c_xy_distance(double *x, double *y, int *nrx, int *nry,
		 int *nc, double *d, int *method)
{
    int i, j;
    size_t  ij;  /* can exceed 2^31 - 1 */
    double (*distfun)(double*, double*, int, int, int, int, int) = NULL;
    switch(*method) {
    case EUCLIDEAN:
	distfun = c_xy_euclidean;
	break;
    case SQEUCLIDEAN:
	distfun = c_xy_sq_euclidean;
	break;
    case CHORD:
	distfun = c_xy_chord;
	break;
    case SQCHORD:
	distfun = c_xy_sq_chord;
	break;
    case BRAY:
	distfun = c_xy_bray;
	break;
    case CHISQUARE:
	distfun = c_xy_chi_square;
	break;
    case SQCHISQUARE:
	distfun = c_xy_sq_chi_square;
	break;
    case INFORMATION:
	distfun = c_xy_information;
	break;
    case MANHATTAN:
	distfun = c_xy_manhattan;
	break;
    case GOWER:
	distfun = c_xy_gower;
	break;
    case ALTGOWER:
	distfun = c_xy_alt_gower;
	break;
    default:
	    error("Unknown distance in the internal C function");
    }
    
    ij = 0;
    for (j=0; j < *nrx; j++)
	for (i=0; i < *nry; i++) {
	    d[ij++] = distfun(x, y, *nrx, *nry, *nc, j, i);
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

SEXP Cdistxy(SEXP x, SEXP y, SEXP smethod)
{
    SEXP ans;
    int nrx = nrows(x), nry = nrows(y), nc = ncols(x), method = asInteger(smethod);
    int N;
    N = (double)nrx * (double)nry; /* avoid overflow for N ~ 50,000 */
    PROTECT(ans = allocVector(REALSXP, N));
    c_xy_distance(REAL(x), REAL(y), &nrx, &nry, &nc, REAL(ans), &method);
    UNPROTECT(1);
    return ans;
}
