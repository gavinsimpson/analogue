/* Standard R headers */

#include <R.h>
#include <Rmath.h>
#include <math.h>

void WA(double *x, double *w, int *nr, int *nc, double *stat)
{
    int i, j, ij;
    double sumW, sumWX;

    for (i=0; i<*nc; i++) {
	sumW = 0.0;
	sumWX = 0.0;
	for(j=0; j<*nr; j++) {
	    ij = (i * *nr) + j;
	    sumW += w[ij];
	    sumWX += w[ij] * x[j];
	}
	stat[i] = sumWX / sumW;
    }
}

void WTOL(double *x, double *w, double *opt, int *nr, int *nc, double *stat)
{
    int i, j, ij;
    /* sumW == sum of weights
       sumWX == sum of weights * Env
       eXbar == Env - Optima
    */
    double sumW, sumWX, eXbar;

    for (i=0; i<*nc; i++) {
	sumW = 0.0;
	sumWX = 0.0;
	for(j=0; j<*nr; j++) {
	    ij = (i * *nr) + j;
	    eXbar = x[j] - opt[i];
	    eXbar = eXbar * eXbar;
	    sumW += w[ij];
	    sumWX += w[ij] * eXbar;
	}
	stat[i] = sqrt(sumWX / sumW);
    }
}

void WATpred(double *spp, double *opt, double *tol2, 
	     int *nr, int *nc, int *want, double *stat)
{
    int i, j, ij;
    /* spp == species data
       opt == species optima
       tol == species tolerances^2
       nr, nc == n rows, n cols
    */
    
    double nomin, denom;

    for (i=0; i<*nr; i++) {
	nomin = 0.0;
	denom= 0.0;
	for(j=0; j<*nc; j++) {
	    ij = i + (j * *nr);
	    nomin += (spp[ij] * opt[j]) / tol2[j];
	    if(*want == 0)
		denom += spp[ij] / tol2[j];
	}
	stat[i] = nomin / denom;
    }
}
