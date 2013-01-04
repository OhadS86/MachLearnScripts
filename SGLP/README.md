# Sparse Group Lasso Projection (SGLP) problem
	sglp.c
	sglp_admm.m
	sglp_dykstra.m
	glLeast.c

# Constrained Sparse Group Lasso Problem (CSGLP)
	sglLeastC.m
	sglLeastC_admm.m

# DC programming for truncated norm
	trunc_sglp.c
	trunc_sglLeastC.m
	trunc_DC_sglLeastC.m

# Usage
The SLEP package is required. You may download the latest version at http://www.public.asu.edu/~jye02/Software/SLEP/.
Compile the c files using ordinary mex tool:

	mex sglp.c/glLeast.c/trunc_sglp.c
 


 
	
	
	
