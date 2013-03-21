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
double xy_euclidean(double *x, double *y, int nr1, int nr2, 
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
double xy_sq_euclidean(double *x, double *y, int nr1, int nr2, 
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
double xy_chord(double *x, double *y, int nr1, int nr2, 
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
double xy_sq_chord(double *x, double *y, int nr1, int nr2, 
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
double xy_bray(double *x, double *y, int nr1, int nr2, 
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
double xy_chi_square(double *x, double *y, int nr1, int nr2, 
		     int nc, int i1, int i2)
{
    double dist, dev, nomin, denomin;
    int count, j;
    
    count = 0;
    dist = 0.0;
    dev = 0.0;

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
double xy_sq_chi_square(double *x, double *y, int nr1, int nr2, 
			int nc, int i1, int i2)
{
    double dist, dev, nomin, denomin;
    int count, j;
	
    count = 0;
    dist = 0.0;
    dev = 0.0;

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
double xy_information(double *x, double *y, int nr1, int nr2, 
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
double xy_manhattan(double *x, double *y, int nr1, int nr2, 
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
double xy_gower(double *x, double *y, int nr1, int nr2, 
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
double xy_alt_gower(double *x, double *y, int nr1, int nr2, 
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

static double (*distfun)(double*, double*, int, int, int, 
			 int, int);

void xy_distance(double *x, double *y, int *nr1, int *nr2,
		 int *nc, double *d, int *method)
{
    int i, j, ij;
    switch(*method) {
    case EUCLIDEAN:
	distfun = xy_euclidean;
	break;
    case SQEUCLIDEAN:
	distfun = xy_sq_euclidean;
	break;
    case CHORD:
	distfun = xy_chord;
	break;
    case SQCHORD:
	distfun = xy_sq_chord;
	break;
    case BRAY:
	distfun = xy_bray;
	break;
    case CHISQUARE:
	distfun = xy_chi_square;
	break;
    case SQCHISQUARE:
	distfun = xy_sq_chi_square;
	break;
    case INFORMATION:
	distfun = xy_information;
	break;
    case MANHATTAN:
	distfun = xy_manhattan;
	break;
    case GOWER:
	distfun = xy_gower;
	break;
    case ALTGOWER:
	distfun = xy_alt_gower;
	break;
    default:
	error("Unknown distance in the internal C function");
    }
    
    ij = 0;
    for (j=0; j < *nr1; j++)
	for (i=0; i < *nr2; i++) {
	    d[ij++] = distfun(x, y, *nr1, *nr2, *nc, j, i);
	}
}

/* 
 * These functions are called directly as they don't fit the
 * nice, ordered manner of the coefficients above
 *
 */

/*
 * Kendall's coefficient
 *
 * Should be called separately from the underlying R code,
 * not via xy_distance.
 *
 * maxi: the max abundance for each species
 *
 */
double xy_KENDALL(double *x, double *y, int nr1, int nr2, 
		  int nc, int i1, int i2, double *maxi)
{
    double dist, dev;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(y[i2])) {
	    dev = (x[i1] >= y[i2]) ? y[i2] : x[i1];
	    dist += maxi[j] - dev;
	    count++;
	}
	i1 += nr1;
	i2 += nr2;
    }
    if (count == 0) return NA_REAL;
    return dist;
}

void xy_kendall(double *x, double *y, int *nr1, int *nr2,
		int *nc, double *d, double *maxi)
{
    int i, j, ij;
    
    ij = 0;
    for(j=0; j < *nr1; j++) {
	for(i=0; i < *nr2; i++) {
	    d[ij++] = xy_KENDALL(x, y, *nr1, *nr2, *nc, j, i,
				 maxi);
	}
    }
}

/*
 * Chi square distance
 *
 * Should be called separately from the underlying R code,
 * not via xy_distance.
 *
 * Needs to be pre processed in R
 *
 * csum: species sums
 *
 */
double xy_CHISQ_DIST(double *x, double *y, int nr1, int nr2, 
		     int nc, int i1, int i2, double *csum,
		     double ccsum)
{
    double dist, dev, denom, nomin;
    int count, j;
    
    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(y[i2])) {
	    dev = x[i1] - y[i2];
	    nomin = dev*dev;
	    denom = csum[j] / ccsum;
	    dist += nomin / denom;
	    count++;
	}
	i1 += nr1;
	i2 += nr2;
    }
    if (count == 0) return NA_REAL;
    return sqrt(dist);
}

void xy_chisq_dist(double *x, double *y, int *nr1, int *nr2,
		   int *nc, double *d, double *csum)
{
    int i, j, k, ij;
    double ccsum;

    ccsum = 0.0;
    
    ij = 0;

    for(k=0; k < *nc; k++) {
	ccsum += csum[k];
    }

    for(j=0; j < *nr1; j++) {
	for(i=0; i < *nr2; i++) {
	    d[ij++] = xy_CHISQ_DIST(x, y, *nr1, *nr2, *nc, 
				    j, i, csum, ccsum);
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
double xy_MIXED_new(double *x, double *y, int nr1, int nr2, 
		    int nc, int i1, int i2, int *vtype,
		    double *weights, double *R, int tmin, int tmax, int ord)
{
    double dist, dev, wsum; //, *curweights;
    int count, j;
    
    count = 0;
    dist = 0.0;
    wsum = 0.0;
    //curweights = weights; /* current weights */

    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(y[i2])) {
	    if(vtype[j] == 1) { // Symmetric binary
		dev = (x[i1] == y[i2]) ? 1 : 0;
		dist += dev * weights[j];
	    }
	    if(vtype[j] == 2) { // Asymmetric binary
		if((x[i1] != 0) || (y[i2] != 0)) {
		    // both x and y not zero for this variables
		    dev = (x[i1] == y[i2]) ? 1 : 0;
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
		dev = (x[i1] == y[i2]) ? 1 : 0;
		dist += dev * weights[j];
	    }
	    if(vtype[j] == 4) { // Ordinal
		/* ordinal data current not handled 
		 * so don't call this yet
		 */
		switch(ord) {
		case 1: { // classic gower as per nominal
		    dev = (x[i1] == y[i2]) ? 1 : 0;
		    dist += dev * weights[j];
		    break;
		}
		case 2: { // podanis rank version
		    if(x[i1] == y[i2]) {
			dev = 1;
		    } else {
			dev = (fabs(x[i1] - y[i2])) / 
			    (R[j] - (tmax - 1)/2 - (tmin - 1)/2);
		    }
		    break;
		}
		case 3: { // podanis metric version treat as Quantitative??
		    dev = 1 - (fabs(x[i1] - y[i2]) / R[j]);
		    dist += dev * weights[j];
		    break;
		}
		default: {
		    dist += 0;
		    break;
		}
		}
	    }
	    if(vtype[j] == 5) { // Quantitative
		dev = 1 - (fabs(x[i1] - y[i2]) / R[j]);
		dist += dev * weights[j];
	    }
	    count++;
	    // only summing weights for non-NA comparisons
	    wsum += weights[j];
	}
	i1 += nr1;
	i2 += nr2;
    }
    if (count == 0) return NA_REAL;
    return 1 - (dist / wsum);
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
double xy_MIXED(double *x, double *y, int nr1, int nr2, 
		int nc, int i1, int i2, int *vtype,
		double *weights, double *R, double wsum)
{
    double dist, dev;
    int count, j;
    
    count = 0;
    dist = 0.0;
    wsum = 0.0;
    //curweights = weights; /* current weights */
    
    for (j=0; j<nc; j++) {
	if (R_FINITE(x[i1]) && R_FINITE(y[i2])) {
	    if(vtype[j] == 1) { // Symmetric binary
		dev = (x[i1] == y[i2]) ? 1 : 0;
		dist += dev * weights[j];
	    }
	    if(vtype[j] == 2) { // Asymmetric binary
		if((x[i1] != 0) || (y[i2] != 0)) {
		    // both x and y not zero for this variables
		    dev = (x[i1] == y[i2]) ? 1 : 0;
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
		dev = (x[i1] == y[i2]) ? 1 : 0;
		dist += dev * weights[j];
	    }
	    if(vtype[j] == 4) { // Ordinal
		dev = (x[i1] == y[i2]) ? 1 : 0;
		dist += dev * weights[j];
		break;
		/* ordinal data current not handled 
		 * so don't call this yet
		 */
		/* switch(ord) {
		   case 1: { // classic gower as per nominal
		   dev = (x[i1] == y[i2]) ? 1 : 0;
		   dist += dev * weights[j];
		   break;
		   }
		   case 2: { // podanis rank version
		   if(x[i1] == y[i2]) {
		   dev = 1;
		   } else {
		   dev = (fabs(x[i1] - y[i2])) / 
		   (R[j] - (tmax - 1)/2 - (tmin - 1)/2);
		   }
		   break;
		   }
		   case 3: { // podanis metric version treat as Quantitative??
		   dev = 1 - (fabs(x[i1] - y[i2]) / R[j]);
		   dist += dev * weights[j];
		   break;
		   }
		   default: {
		   dist += 0;
		   break;
		   }
		   }*/
	    }
	    if(vtype[j] == 5) { // Quantitative
		dev = 1 - (fabs(x[i1] - y[i2]) / R[j]);
		dist += dev * weights[j];
	    }
	    count++;
	    // only summing weights for non-NA comparisons
	    wsum += weights[j];
	}
	i1 += nr1;
	i2 += nr2;
    }
    if (count == 0) return NA_REAL;
    return 1 - (dist / wsum);
}

double xy_calcTI(double *x, double *y, int nr1, int nr2, int nc, int i1, int i2)
{
    int k;
    double ti;

    ti = 0.0;

    for (k=0; k<nc; k++) {
	ti += (x[i1] == y[i2]) ? 1.0 : 0.0;
	i1 += nr1;
	i2 += nr2;
    }
    return ti;
}

void xy_mixed_new(double *x, double *y, int *nr1, int *nr2,
		  int *nc, double *d, int *vtype, double *weights, 
		  double *R, int tmin, int tmax, int *podani, 
		  double *tmatx, double *tmaty)
{
    int i, j, ij, ord = *podani;
    //double tmat = *d;
    
    ij = 0;

    for(j=0; j < *nr1; j++) {
	for(i=0; i < *nr2; i++) {
	    tmatx[ij] = xy_calcTI(x, y, *nr1, *nr2, *nc, j, i);
	    //tmaty[ij++] = xy_calcTI(
	}
    }

    ij = 0;

    for(j=0; j < *nr1; j++) {
	for(i=0; i < *nr2; i++) {
	    d[ij++] = xy_MIXED_new(x, y, *nr1, *nr2, *nc, j, i,
				   vtype, weights, R, tmin, tmax, ord);
	}
    }
}

void xy_mixed(double *x, double *y, int *nr1, int *nr2,
	      int *nc, double *d, int *vtype, double *weights, 
	      double *R)
{
    int i, j, k, ij;
    double wsum;
    
    wsum = 0.0;
    
    ij = 0;
    
    for(k=0; k <*nc; k++) {
	wsum += weights[k];
    }
    
    for(j=0; j < *nr1; j++) {
	for(i=0; i < *nr2; i++) {
	    d[ij++] = xy_MIXED(x, y, *nr1, *nr2, *nc, j, i,
			       vtype, weights, R, wsum);
	}
    }
}
