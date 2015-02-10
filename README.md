# RMUMPS - R Interface for MUMPS

Written by Mikko Orispaa (mikko.orispaa@gmail.com)
(c) 2015, Mikko Orispaa

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



## Roadmap (Delete before release)

### Version 0.1
+ basic function
+ able to solve symmetric sparse equations with ONE (dense) RHS

### Version 0.2
+ Able to solve systems with multiple (matrix) RHS

### Version 0.4
- Able to calculate elements of A^{-1}
- Use of MUMPS parameters
+ Able to solve general sparse systems

### Version 0.6
- able to solve systems with sparse matrix RHS

### Version 1.0
- documentation added
- INSTALL file with installation instructions written
