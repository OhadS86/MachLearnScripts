/*
 * Matlab usage: x = sght(v, ind, s1, s2)
 */
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <cmath>
#include <mex.h>

using namespace std;

int n, g;
int flmt, glmt;
double *v, *x;
int *vidx;
int *gidx;

// 1-based
double ***d;
int ***path;
double **selection_value;

bool cmp (int x, int y) {
    return fabs(v[x])  > fabs(v[y]) + 1e-9;
}

void dp_preprocess () {
	int i, j, k;

	d = (double ***)malloc((g + 1) * sizeof(double **));
	path = (int ***)malloc((g + 1) * sizeof(int **));
	for (i = 0; i <= g; ++i) {
		d[i] = (double **)malloc((glmt + 1) * sizeof(double *));
		path[i] = (int **)malloc((glmt + 1) * sizeof(int *));
		for (j = 0; j <= glmt; ++j) {
			d[i][j] = (double *)malloc((flmt + 1) * sizeof(double));
			path[i][j] = (int *)malloc((flmt + 1) * sizeof(int));
			for (k = 0; k <= flmt; ++k) {
				d[i][j][k] = .0;
				path[i][j][k] = 0;
			}
		}
	}
	int up = 0;
	for (i = 1; i <= g; ++i)
		up = max(up, gidx[i] - gidx[i-1]);

	selection_value = (double **)malloc((g + 1) * sizeof(double *));
	for (i = 0; i <= g; ++i)
		selection_value[i] = (double *)malloc((up + 1) * sizeof(double));

	vidx = (int *)malloc(n * sizeof(int));
	for (i = 0; i < n; ++i) vidx[i] = i;

	for (i = 1; i <= g; ++i) {
		sort(vidx + gidx[i-1], vidx + gidx[i], cmp);
		selection_value[i][0] = .0;
		for (j = 1; j <= gidx[i] - gidx[i-1]; ++j)
			selection_value[i][j] = selection_value[i][j-1] +
				v[vidx[gidx[i-1]+j-1]] * v[vidx[gidx[i-1]+j-1]];
	}
}

void dp () {
	int i, j, k, t;

	dp_preprocess ();

	for (i = 1; i <= g; ++i) {
		for (j = 1; j <= glmt; ++j)
			for (k = 1; k <= flmt; ++k) {
				d[i][j][k] = d[i-1][j][k];
				int u = min(k, gidx[i] - gidx[i-1]); // min(|G_i|, k)
				int max_idx = 0;
				double best = d[i][j][k], v;
				for (t = 1; t <= u; ++t) {
					if ((v = d[i-1][j-1][k-t] + selection_value[i][t]) > best) {
						best = v;
						max_idx = t;
					}
				}
				d[i][j][k] = best;
				path[i][j][k] = max_idx;
			}
	}
}

void calc_sol () {
	int i, j, k;
	memset (x, 0, n * sizeof(double));

	for (i = g, j = glmt, k = flmt; i >= 1; --i) {
		int num = path[i][j][k];
		k -= num;
		if (num) --j;
		for (int t = 1; t <= num; ++t)
			x[vidx[gidx[i-1]+t-1]] = v[vidx[gidx[i-1]+t-1]];
	}
}

void destruct () {
	int i, j;

//free(v); v = NULL;
	free(vidx); vidx = NULL;
	free(gidx); gidx = NULL;
	for (i = 0; i <= g; ++i) {
		for (j = 0; j <= glmt; ++j) {
			free(d[i][j]); d[i][j] = NULL;
			free(path[i][j]); path[i][j] = NULL;
		}
		free(d[i]); d[i] = NULL;
		free(path[i]); path[i] = NULL;
	}
	free(d); d = NULL;
	free(path); path = NULL;
    
	for (i = 0; i <= g; ++i) {
		free(selection_value[i]);
		selection_value[i] = NULL;
	}
	free(selection_value);
	selection_value = NULL;
}

void mexFunction (int nlhs, mxArray* plhs[],
				  int nrhs, const mxArray* prhs[])
{
	v = mxGetPr(prhs[0]);
	double* gidx_double = mxGetPr(prhs[1]);
    
	n = mxGetNumberOfElements(prhs[0]);
	g = mxGetNumberOfElements(prhs[1]) - 1;
    
	gidx = (int *)malloc((g + 1) * sizeof(int));
	for (int i = 0; i <= g; ++i) gidx[i] = (int)gidx_double[i];
    
	flmt = mxGetScalar(prhs[2]);
	glmt = mxGetScalar(prhs[3]);
    
	double eps = 1e-7;
    
	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
	x = mxGetPr(plhs[0]);
    
	dp ();
	calc_sol ();
//mexPrintf ("sglt: Optimal solution is %.5lf\n", d[g][glmt][flmt]);
    
	destruct ();
}

