# RMUMPS - R Interface for MUMPS

Written by Mikko Orispaa (mikko.orispaa@gmail.com)

Copyright (c) 2015, Mikko Orispaa

Licensed under MIT License. See file LICENSE for details.


## MUMPS (MUltifrontal Massively Parallel sparse direct Solver)

MUMPS is a public domain parallel sparse solver based on public domain software eveloped during the Esprit IV European project PARASOL (1996-1999).

For details, see: http://mumps.enseeiht.fr


## INSTALLATION OF MUMPS

You need MUMPS libraries installed to use rmumps package.

(Note that you also need MPI implementation such as MPICH with mpicc compiler installed in your system.)

### Mac OS X

The easiest way to install MUMPS is to use MacPorts package manager

´´´
sudo port install mumps
´´´

as it will also install all required dependencies.

It is also possible to download the source code from MUMPS webpage (http://mumps.enseeiht.fr) and compile it yourself. In this case you have to handle the dependencies yourself.

### Linux

MUMPS is available at least in Ubuntu repositories:
´´´
sudo apt-get install libmumps-4.10.0 
´´´
(or whatever the current version number is).

It is also available for CentOS & yum:
´´´
sudo yum install MUMPS-openmpi
´´´

It is also possible to download the source code from MUMPS webpage (http://mumps.enseeiht.fr) and compile it yourself. In this case you have to handle the dependencies yourself.

## INSTALLATION OF RMUMPS

The rmumps package assumes that the MUMPS libraries are either in /opt/local/lib or somewhere else where mpicc compiler can find them. If this is not the case, you need to modify the file src/Makefile by changing the variables IDIR and LDIR to point to the installation path of MUMPS.


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
