/*function [ x, lambda ] = glLeastC( v, radius, ind, w, u )

 Solve the following group projection:

    \min_{x}      \| x - v \|_2
    subject to    \sum \|x_{G_i}\|_2 <= radius,

    where the group information is contained in the "ind" array
    Require SLEP package, the "eplb" function.
*/

#include <math.h>
#include "mex.h"
#include "SLEP_4.0/SLEP/CFiles/q1/epph.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	double *v, *ind, *w, *u, *x, *lambda;
	double radius, lambda0, sum;

	int n, grpNum;
	int i, j, k, steps, be, en;

	if(nlhs != 2)
		mexErrMsgTxt("Error: Number of output parameter does not match.");
	
	if(nrhs !=5)
		mexErrMsgTxt("Error: Number of input parameter does not match.");
    
	n = mxGetNumberOfElements(prhs[0]);
	grpNum = mxGetNumberOfElements(prhs[2])  - 1;
	lambda0 = .0;
    
	v = mxGetPr(prhs[0]);
	radius = mxGetScalar(prhs[1]);
	ind = mxGetPr(prhs[2]);
	w = mxGetPr(prhs[3]);
	u = mxGetPr(prhs[4]);
	

	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	x = mxGetPr(plhs[0]);
	lambda = mxGetPr(plhs[1]);
    
	for(i = 0; i < grpNum; i++)
	{
		sum = .0;
		be = (int)ind[i], en = (int)ind[i+1];
		for(j = be; j < en; j++)
			sum += v[j] * v[j];
		w[i] = sqrt(sum);
	}

	eplb(u, lambda, &steps, w, grpNum, radius, lambda0);

	for(i = 0, j = 0; j < grpNum; j++)
	{
		be = (int)ind[j], en = (int)ind[j+1];
        
		if(fabs(w[j]) > 1e-15)
		for(i = be; i < en; i++)
			x[i] = u[j] * v[i] / w[j];
		else
			for(i = be; i < en; i++)
				x[i] = 0;
	} 
}

