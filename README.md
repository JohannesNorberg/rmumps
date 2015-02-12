# RMUMPS - R Interface for MUMPS

Written by Mikko Orispaa (mikko.orispaa@gmail.com)

Copyright (c) 2015, Mikko Orispaa

Licensed under MIT License. See file LICENSE for details.

**MUMPS stands for a MUltifrontal Massively Parallel sparse direct Solver.**

## History

2015-02-04
- Development started

2015-02-09
- Ability to solve unsymmetric and symmetric systems
- Ability to have (dense) matrix right hand side
- tagged as version 0.2.

2015-02-10
- Ability to calculate inverse matrix diagonal
- ICNTL(7) fixed for 0 for unsymmetric matrices. Default value 7 hangs the calculation
  (MATIS problem)
- tagged as 0.2-1

2015-02-12
- Ability to calculate arbitrary inverse matrix elements.




## INSTALLATION

You need MUMPS libraries installed to use rmumps package.

Write instructions for installing MUMPS here...
