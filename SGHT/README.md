Sparse Group Hard Thresholding (SGHT)

> Sparse Group Feature Selection (SGFS) via SGHT.


## Introduction
The `.m` files are the matlab interfaces for solving the Sparse Group Feature Selectio problem while the `sght_*.cpp` files contain the key proximal parts, i.e., the **Sparse Group Hard Thresholding (SGHT)** problem. Currently there are three versions available:

 - `sght.cpp`: **recommended**, best readability, regular DP format
 - `sght_external.cpp`: use external memory, complicated indices conversion
 - `sght_persisten.cpp`: **unstable and deprecated**, persistent the DP table, complicated indices conversion, needs `clear sght_persistent` properly

## Functions

| Tables              | FISTA           | ISTA  | Barzilai-Borwein | Const | Lipschiz | Sufficient Decrease|
| --------------------|:---------------:|:-----:|:----------------:|:-----:|:--------:|:------------------:|
| `sghtFISTA.m`       | Y               |       |Y                 |       |Y         |                    |
| `sghtFISTAConst.m`  | Y               |       |                  |Y      |Y         |                    |
| `sghtISTA.m`        |                 |Y      |Y                 |       |          |Y                   |   
| `sghtISTAConst.m`   |                 |Y      |                  |Y      |          |Y                   |
| `sghtISTAWolfe.m`   |                 |Y      |Y                 |       |Y         |                    |   


 - framework: FISTA/ISTA:
 - Line search initialization: Barzilai-Borwein/const
 - Line search criterion: Lipschiz/sufficient decrease

## Usage

See `test_sght.m` for details of calling the functions in matlab. Make sure do `mex sght.cpp` before that.
