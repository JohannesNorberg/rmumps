MPICC=mpicc
IDIR=/opt/local/include
LDIR=/opt/local/lib
LIBS=-ldmumps -lmumps_common
CFLAGS=-O3 -I$(IDIR) -L$(LDIR)




all: build


build: mumpsdrv.c
	$(MPICC) mumpsdrv.c -o mumpsdrv $(LIBS) $(CFLAGS)
