iSFS
============

*incomplete Source and Feature Selection (iSFS) model*

> iSFS is a feature learning model for multi-modal block-wise missing data. 
Linear models on both of the feature-level and source-level are learned simultaneously in a joint formulation.

## Dependency
 - [SLEP](http://www.public.asu.edu/~jye02/Software/SLEP/)
 - [Random Forest](http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm)

Both of them are included in `Tools`

## Setup
 - Compile Random Forest and SLEP according to their manual 
 - Initialize iSFS by running `Init.m` under `Tools/iMAD/Function`
 - Put your (block-wise missing) data under folder `Data/`
 - Implement you own `processData` under folder `Tools/iMAD/Function/processData.m`

## Usage
Run `Main.m` under `Modules/`



