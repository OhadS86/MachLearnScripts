## Sparse Group Lasso Projection (SGLP) problem

	sglp.c
	sglp_admm.m
	sglp_dykstra.m
	glLeast.c

## Constrained Sparse Group Lasso Problem (CSGLP)

	sglLeastC.m
	sglLeastC_admm.m

## DC programming for truncated norm

	trunc_sglp.c
	trunc_sglLeastC.m
	trunc_DC_sglLeastC.m

## Usage

The SLEP package is required. You may download the latest version <a href="http://www.public.asu.edu/~jye02/Software/SLEP/">here</a>.
Use mex to compile the c files:

	mex sglp.c/glLeast.c/trunc_sglp.c
 


 
	
	
	
