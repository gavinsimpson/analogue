/* 
 * Functions to compute dissimilarity for a single matrix
 *
 * Based on code from vegdist by Jari Oksanen:
 *
 * (C) 2001-2009, Jari Oksanen
 * (C) 2009 Gavin L. Simpson
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
double xx_euclidean(double *x, int nr, int nc, int i1, int i2)
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
double xx_sq_euclidean(double *x, int nr, int nc, int i1, int i2)
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
double xx_chord(double *x, int nr, int nc, int i1, int i2)
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
double xx_sq_chord(double *x, int nr, int nc, int i1, int i2)
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
double xx_bray(double *x, int nr, int nc, int i1, int i2)
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
double xx_chi_square(double *x, int nr,  int nc, int i1, int i2)
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
double xx_sq_chi_square(double *x, int nr, int nc, int i1, int i2)
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
double xx_information(double *x, int nr, int nc, int i1, int i2)
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
double xx_chi_distance(double *x, int nr, int nc, int i1, int i2)
{
    double dist, dev, nomin;
    int count, j;
    
    count = 0;
    dist = 0.0;
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
double xx_manhattan(double *x, int nr, int nc, int i1, int i2)
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
double xx_gower(double *x, int nr, int nc, int i1, int i2)
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
double xx_alt_gower(double *x, int nr, int nc, int i1, int i2)
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

void xx_distance(double *x, int *nr, int *nc, double *d, 
		 int *diag, int *method)
{
    int dc, i, j, ij;
    switch(*method) {
    case EUCLIDEAN:
	distfun = xx_euclidean;
	break;
    case SQEUCLIDEAN:
	distfun = xx_sq_euclidean;
	break;
    case CHORD:
	distfun = xx_chord;
	break;
    case SQCHORD:
	distfun = xx_sq_chord;
	break;
    case BRAY:
	distfun = xx_bray;
	break;
    case CHISQUARE:
	distfun = xx_chi_square;
	break;
    case SQCHISQUARE:
	distfun = xx_sq_chi_square;
	break;
    case INFORMATION:
	distfun = xx_information;
	break;
    case MANHATTAN:
	distfun = xx_manhattan;
	break;
    case GOWER:
	distfun = xx_gower;
	break;
    case ALTGOWER:
	distfun = xx_alt_gower;
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

/*
 * Kendall's coefficient
 *
 * Should be called separately from the underlying R code,
 * not via xy_distance.
 *
 * maxi: the max abundance for each species
 *
 */
double xx_KENDALL(double *x, int nr, int nc, int i1, int i2,
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

void xx_kendall(double *x, int *nr, int *nc, double *d, 
		double *maxi)
{
    int i, j, ij;
    double dist;

    dist = 0.0;
    
    ij = 0;
    for(j=0; j <= *nr; j++) {
	for(i=0; i < *nr; i++) {
	    d[ij++] = xx_KENDALL(x, *nr, *nc, i, j, maxi); 
	}
    }
}

/*
 * Chi square distance
 *
 * Should be called separately from the underlying R code,
 * not via xy_distance.
 *
 * csum: species sums
 *
 */
double xx_CHISQ_DIST(double *x, int nr, int nc, int i1, int i2,
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

void xx_chisq_dist(double *x, int *nr, int *nc, double *d, 
		   double *csum)
{
    int i, j, k, ij;
    double ccsum;

    ccsum = 0.0;
    
    ij = 0;

    for(k=0; k < *nc; k++) {
	ccsum += csum[k];
    }

    for(j=0; j < *nr; j++) {
	for(i=0; i < *nr; i++) {
	    d[ij++] = xx_CHISQ_DIST(x, *nr, *nc, i, j, 
				    csum, ccsum);
	}
    }
}

/*
 * Gower's coefficient for mixed data
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
double xx_MIXED(double *x, int nr, int nc, int i1, int i2, 
		int *vtype, double *weights, double *R, 
		double wsum)
{
    double dist, dev;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    if(vtype[j] == 1) {
		dev = (x[i1] == x[i2]) ? 1 : 0;
		dist += dev * weights[j];
	    }
	    if(vtype[j] == 2) {
		dev = (x[i1] == x[i2]) ?   : 0;
		dist += dev * weights[j];
	    }
	    if(vtype[j] == 3) {
		dev = (x[i1] == x[i2]) ? 1 : 0;
		dist += dev * weights[j];
	    }
	    if(vtype[j] == 4) {
		/* ordinal data current not handled 
		 * so don't call this yet
		 */
	    }
	    if(vtype[j] == 5) {
		dev = 1 - (fabs(x[i1] - x[i2]) / R[j]);
		dist += dev * weights[j];
	    }
	    count++;
	}
	i1 += nr;
	i2 += nr;
    }
    if (count == 0) return NA_REAL;
    return 1 - (dist / wsum);
}

void xx_mixed(double *x, int *nr, int *nc, double *d, 
	      int *vtype, double *weights, double *R)
{
    int i, j, k, ij;
    double wsum;

    wsum = 0.0;
    
    ij = 0;

    for(k=0; k <*nc; k++) {
	wsum += weights[k];
    }

    for(j=0; j < *nr; j++) {
	for(i=0; i < *nr; i++) {
	    d[ij++] = xx_MIXED(x, *nr, *nc, i, j, vtype,
			       weights, R, wsum);
	}
    }
}
