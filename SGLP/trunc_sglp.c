/*function x = sglp(v, s1, s2, ind, supp1, suppG);

 Solve the following truncated group projection:

		\min_{x}      \| x - v \|_2
		subject to    \| x_A \|_1 <= s1
					  \sum \|x_B_{G_i}\|_2 <= s2
    
    x_A: a subsequence of x, indexed by supp1
    x_B: a subsequence of x, indexed by suppG. suppG \subset supp1
    Group information is contained in the "ind" array
    where the group information is contained in the "ind" array
    
 * Require SLEP package, the "eplb" function.
*/
#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include "SLEP_4.0/SLEP/CFiles/q1/epph.h"

#define MAX(x, y) ((x) > (y)? (x) : (y))

double l1Norm(double *x, double *supp1, int sizeSupp1)
{
	int i;
	double ans = .0;
	
	for(i = 0; i < sizeSupp1; i++)
		ans += fabs(x[(int)supp1[i]-1]);
	return ans;
}

double grpNorm(double *x, double *ind, double *suppG, int sizeSuppG)
{
	int be, en, i, j, k;
	double sum, ans = .0;

	for(i = 0; i < sizeSuppG; i++)
	{
		k = (int)suppG[i];
		sum = .0;
		be = (int)ind[k-1], en = (int)ind[k];
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

// epglb(x_c2, v, n, s2, ind, grpNum, tmpW, tmpU);
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
		{
			double fct = u[j] / w[j];
            for(i = be; i < en; i++)
                x[i] = fct * v[i];
		}
        else
            for(i = be; i < en; i++)
                x[i] = .0;
    } 
	
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *v, *ans, s1, s2;
	double *ind, *supp1, *suppG;
	double *x;
	double *u, *tmpV, *tmpU, *tmpW, *v_1, *v_g;
	double up, mid, low, h_s1, lambda0, lambda;
	int n, grpNum, sizeSupp1, sizeSuppG;
	int i, j, k, t, be, en;
	
	if(nlhs != 1)
		mexErrMsgTxt("Error: Number of output parameter does not match.");
	
	if(nrhs !=6)
		mexErrMsgTxt("Error: Number of input parameter does not match.");
	
	/*
	 *  Input Parameter handle
	 */
	v = mxGetPr(prhs[0]);
	s1 = mxGetScalar(prhs[1]);
	s2 = mxGetScalar(prhs[2]);
	ind = mxGetPr(prhs[3]);
	supp1 = mxGetPr(prhs[4]);
	suppG = mxGetPr(prhs[5]);
	
	
	n = (int)mxGetNumberOfElements(prhs[0]);
	grpNum = (int)mxGetNumberOfElements(prhs[3])  - 1;
	sizeSupp1 = (int)mxGetNumberOfElements(prhs[4]);
	sizeSuppG = (int)mxGetNumberOfElements(prhs[5]);
	
	
	/*
	 *	Output Parameter handle
	 */
	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
	ans = mxGetPr(plhs[0]);
	
	// Temporary variables.
	lambda = lambda0 = .0;
	
	x = (double *)malloc(n * sizeof(double));
	tmpV = (double *)malloc(n * sizeof(double));
	v_1 = (double *)malloc(n * sizeof(double));
	v_g = (double *)malloc(n * sizeof(double));
	
	
	tmpU = (double *)malloc(grpNum * sizeof(double)); 
	tmpW = (double *)malloc(grpNum * sizeof(double)); 
	
	if(l1Norm(v, supp1, sizeSupp1) <= s1 + 1e-12 && grpNorm(v, ind, suppG, sizeSuppG) <= s2 + 1e-12)
	{
		for(i = 0; i < n; i++)
			ans[i] = v[i];
		
		free(x);
		free(tmpV);
		free(v_1);
		free(v_g);
		free(tmpU);
		free(tmpW);
	
		return;
	}
	
	lambda0 = .0;
	
	
	
/*
 *	Step 1: Check if L1 ball constraint is active
 */	
	
	// Preprocessing before projection onto L_1 ball
	for(i = 0, j = 0; i < n; i++)
	{
		if(j < sizeSupp1 && i+1 == (int)supp1[j]) {
			v_1[i] = v[i]; 
			j++;
		}
		else{
			if(j == sizeSupp1) 
				while(i < n) v_1[i++] = .0;
			else{
				while(i+1 != (int)supp1[j]) v_1[i++] = .0;
				i--;
			}
		}
	}
	
	eplb(x, &lambda, &k, v_1, n, s1, lambda0); 
	
	// Postprocessing after projection onto L_1 ball
	for(i = 0, j = 0; i < n; i++)
	{
		if(j < sizeSupp1 && i+1 == (int)supp1[j]) j++;
		else{
			if(j == sizeSupp1)
				while(i < n) 
				{
					x[i] = v[i];
					i++;
				}
			else{
				while(i+1 != (int)supp1[j]) {x[i] = v[i]; i++;}
				i--;
			}
		}
	}
	if(grpNorm(x, ind, suppG, sizeSuppG) < s2 - 1e-12)
	{
		for(i = 0; i < n; i++)
			ans[i] = x[i];
		free(x);
		free(tmpV);
		free(v_1);
		free(v_g);
		free(tmpU);
		free(tmpW);
		return;
	}
	
/*
 *	Step 2: Check if L_G ball constraint is active
 */	
	
	// Preprocessing before projection onto L_G ball
	for(i = 0, j = 0; i < grpNum; i++)
	{
		if(j < sizeSuppG && i + 1 == (int)suppG[j]){
			be = (int)ind[i], en = (int)ind[i+1];
			for(k = be; k < en; k++)
				v_g[k] = v[k];
			j++;
		}
		else{
			be = (int)ind[i], en = (int)ind[i+1];
			for(k = be; k < en; k++)
				v_g[k] = .0;
		}
	}
	epglb(x, v_g, n, s2, ind, grpNum, tmpW, tmpU); 
	
	// Postprocessing after projection onto L_G ball
	for(i = 0, j = 0; i < grpNum; i++)
	{
		if(j < sizeSuppG && i + 1 == (int)suppG[j]) j++;
		else{
			be = (int)ind[i], en = (int)ind[i+1];
			for(k = be; k < en; k++)
				x[k] = v[k];
		}
	}
	
	if(l1Norm(x, supp1, sizeSupp1) < s1 - 1e-12)
	{
		for(i = 0; i < n; i++)
			ans[i] = x[i];
		free(x);
		free(tmpV);
		free(v_1);
		free(v_g);
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
		
		if(grpNorm(tmpV, ind, suppG, sizeSuppG) < s2 - 1e-15)
			up = mid;
		else{
			// Preprocessing before projection onto L_G ball
			for(i = 0, j = 0 ; i < grpNum; i++)
			{
				if(j < sizeSuppG && i + 1 == (int)suppG[j]){
					be = (int)ind[i], en = (int)ind[i+1];
					for(k = be; k < en; k++)
						v_g[k] = tmpV[k];
					j++;
				}
				else{
					be = (int)ind[i], en = (int)ind[i+1];
					for(k = be; k < en; k++)
						v_g[k] = .0;
				}
			}
	
			epglb(x, v_g, n, s2, ind, grpNum, tmpW, tmpU);
			
			// Postprocessing after projection onto L_G ball
			for(i = 0, j = 0; i < grpNum; i++)
			{
				if(j < sizeSuppG && i + 1 == (int)suppG[j]) j++;
				else{
					be = (int)ind[i], en = (int)ind[i+1];
					for(k = be; k < en; k++)
						x[k] = tmpV[k];
				}
			}// Notice that the x here is not really the \hat{x}, however
			// it doesn't matter for the following calculation.
			
/*
 *	Calculate \hat{s_1}.
 */
			h_s1 = l1Norm(x, supp1, sizeSupp1);		
			if(h_s1 < s1 + 1e-15)
				up = mid;
			else
				low = mid;
		}
	}
	soft_thresholding(tmpV, v, n, up);
	
	// Preprocessing before projection onto L_G ball
	for(i = 0, j = 0; i < grpNum; i++)
	{
		if(j < sizeSuppG && i + 1 == (int)suppG[j]){
			be = (int)ind[i], en = (int)ind[i+1];
			for(k = be; k < en; k++)
				v_g[k] = tmpV[k];
			j++;
		}
		else{
			be = (int)ind[i], en = (int)ind[i+1];
			for(k = be; k < en; k++)
				v_g[k] = 0;
		}
	}
	epglb(x, v_g, n, s2, ind, grpNum, tmpW, tmpU);
	for(i = 0, j = 0, t = 0; i < grpNum; i++)
	{
		if(j < sizeSuppG && i + 1 == (int)suppG[j]) j++;
		else{
			be = (int)ind[i], en = (int)ind[i+1];
			for(k = be; k < en; k++){
				while(t < sizeSupp1 && (int)supp1[t] < k+1) t++;
				if(t < sizeSupp1 && (int)supp1[t] == k + 1)
					x[k] = tmpV[k];
				else x[k] = v[k];
			}
		}
	}
	for(i = 0; i < n; i++)
		ans[i] = x[i];
	
	free(x);
	free(tmpV);
	free(v_1);
	free(v_g);
	free(tmpU);
	free(tmpW);
}