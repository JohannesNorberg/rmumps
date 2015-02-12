MPICC=mpicc
IDIR=/opt/local/include
LDIR=/opt/local/lib
LIBS=-ldmumps
CFLAGS=-O3 -I$(IDIR) -L$(LDIR)




all: build


build: mumpsdrv.c
	$(MPICC) $(CFLAGS) mumpsdrv.c -o mumpsdrv $(LIBS)
