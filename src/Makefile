MPICC=mpicc
CFLAGS=-O3
LIBS=-ldmumps -lmumps_common -lscalapack-openmpi -lopenblas -lpthread

all: mumpsdrv

mumpsdrv: mumpsdrv.c
    $(MPICC) $(CFLAGS) mumpsdrv.c -o mumpsdrv $(LIBS)

clean:
    rm -f mumpsdrv
