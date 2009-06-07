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
