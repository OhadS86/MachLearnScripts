/*function x = sglp(v, s1, s2, ind);

 Solve the following group projection:

    \min_{x}      1/2 \| x - v \|_2^2
    subject to    \|x\|_1 <= s1
 *                \sum \|x_{G_i}\|_2 <= radius,

    where the group information is contained in the "ind" array
    Require SLEP package, the "eplb" function.
*/

#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include "SLEP_4.0/SLEP/CFiles/q1/epph.h"

#define MAX(x, y) ((x) > (y)? (x) : (y))

double l2Norm(double *x, int n)
{
	int i;
	double ans = .0;
	for(i = 0; i < n; i++)
		ans += x[i] * x[i];
	return sqrt(ans);
}

double l1Norm(double *x, int n)
{
	int i;
	double ans = .0;
	for(i = 0; i < n; i++)
		ans += fabs(x[i]);
	return ans;
}

double grpNorm(double *x, double *ind, int grpNum)
{
	int be, en, i, j;
	double sum, ans = .0;
	for(i = 0; i < grpNum; i++)
	{
		sum = .0;
		be = (int)ind[i], en = (int)ind[i+1];
		for(j = be; j < en; j++)
			sum += x[j] * x[j];
		ans += sqrt(sum);
	}
	return ans;
}

double calFunVal(double *x, double *v, int n)
{
	int i;
	double ans = .0;
	for(i = 0; i < n; i++)
		ans += (x[i] - v[i]) * (x[i] - v[i]);
	return 0.5 * ans;
}

void soft_thresholding(double *ans, double *v, int n, double z)
{
	int i;
	for(i = 0; i < n; i++)
		ans[i] = MAX(0, v[i] - z) - MAX(0, -v[i] - z);
}

int cmp(double x, double y)
{
	return (fabs(y+1) < 1e-15 || y-x > 1e-15);
}
	
/* epglb(x_c2, v, n, s2, ind, grpNum, tmpW, tmpU); */
void epglb(double *x, double *v, int n, double radius, double *ind, int grpNum,
		double *w, double *u)
{
	int i, j, be, en, steps;
	double lambda, sum, lambda0;
	
	for(i = 0; i < grpNum; i++)
	{
		sum = .0;
		be = (int)ind[i], en = (int)ind[i+1];
		for(j = be; j < en; j++)
			sum += v[j] * v[j];
		w[i] = sqrt(sum);
	}
	
	lambda = lambda0 = .0;
	
	eplb(u, &lambda, &steps, w, grpNum, radius, lambda0);
	for(i = 0, j = 0; j < grpNum; j++)
	{
		be = (int)ind[j], en = (int)ind[j+1];
      
		if(fabs(w[j]) > 1e-15)
			for(i = be; i < en; i++)
				x[i] = u[j] * v[i] / w[j];
		else
			for(i = be; i < en; i++)
				x[i] = .0;
	} 
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *v, *ans, s1, s2;
	double *ind;
	double *x;
	double *u, *tmpV, *tmpU, *tmpW;
	double up, mid, low, h_s1, lambda0, lambda;
	int n, grpNum;
	int i, j, k;
	
	if(nlhs != 1)
		mexErrMsgTxt("Error: Number of output parameter does not match.");
	
	if(nrhs !=4)
		mexErrMsgTxt("Error: Number of input parameter does not match.");
	
	/*
	 *  Input Parameter handle
	 */
	v = mxGetPr(prhs[0]);
	s1 = mxGetScalar(prhs[1]);
	s2 = mxGetScalar(prhs[2]);
	ind = mxGetPr(prhs[3]);
	
	
	n = (int)mxGetNumberOfElements(prhs[0]);
	grpNum = (int)mxGetNumberOfElements(prhs[3])  - 1;
	
	/*
	 *	Output Parameter handle
	 */
	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
	ans = mxGetPr(plhs[0]);
	
	/* Temporary variables. */
	lambda = lambda0 = .0;
	
	x = (double *)malloc(n * sizeof(double));
	tmpV = (double *)malloc(n * sizeof(double));

	tmpU = (double *)malloc(grpNum * sizeof(double));
	tmpW = (double *)malloc(grpNum * sizeof(double));
	
	if(l1Norm(v, n) <= s1 + 1e-12 && grpNorm(v, ind, grpNum) <= s2 + 1e-12)
	{
		for(i = 0; i < n; i++)
			ans[i] = v[i];
		
		free(x);
		free(tmpV);
		free(tmpU);
		free(tmpW);
		return;
	}
	
	lambda0 = .0;
	
	eplb(x, &lambda, &k, v, n, s1, lambda0);
	
	if(grpNorm(x, ind, grpNum) < s2 - 1e-12)
	{
		for(i = 0; i < n; i++)
			ans[i] = x[i];
		
		free(x);
		free(tmpV);
		free(tmpU);
		free(tmpW);
		return;
	}
	
	epglb(x, v, n, s2, ind, grpNum, tmpW, tmpU);
	
	if(l1Norm(x, n) < s1 - 1e-12)
	{
		for(i = 0; i < n; i++)
			ans[i] = x[i];
		
		free(x);
		free(tmpV);
		free(tmpU);
		free(tmpW);
		return;
	}
	up = 1e15, low = 1e-12;
	while(up - low > 1e-7)
	{
		if(low > 1e6 && up - low < 1e-5)
			break;
		
		mid = (up + low) / 2;
		soft_thresholding(tmpV, v, n, mid);
		
		if(grpNorm(tmpV, ind, grpNum) < s2 - 1e-15)
			up = mid;
		else{
			epglb(x, tmpV, n, s2, ind, grpNum, tmpW, tmpU);
			h_s1 = l1Norm(x, n);
			if(h_s1 < s1 + 1e-15)
				up = mid;
			else
				low = mid;
		}
	}
	soft_thresholding(tmpV, v, n, up);
	epglb(x, tmpV, n, s2, ind, grpNum, tmpW, tmpU);
	
	for(i = 0; i < n; i++)
		ans[i] = x[i];
	
	free(x);
	free(tmpV);
	free(tmpU);
	free(tmpW);
}
