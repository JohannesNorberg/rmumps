# RMUMPS - R Interface for MUMPS

Written by Mikko Orispaa (mikko.orispaa@gmail.com)

Copyright 2017 Mikko Orispää <mikko.orispaa@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.


## RMUMPS

RMUMPS is a simple non-sophisticated R interface for MUMPS library. It uses binary files to transmit data between R and a external MUMPS driver program which is included in the R package. To install rmumps, MUMPS libraries and a working MPI installation (with mpicc and mpirun) are required.


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

2015-02-24
- Fixed a bug that prevented using sparse matrices with over 2^28 non-zero element (writeBin can
  not handle data chunks larger than 2^31 bytes)
