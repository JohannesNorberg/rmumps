MPICC=mpicc
CFLAGS=-O3
LIBS=-ldmumps -lmumps_common -lscalapack-openmpi -lopenblas -lpthread

all: mumpsdrv

mumpsdrv: mumpsdrv.c
	$(MPICC) $(CFLAGS) mumpsdrv.c -o mumpsdrv $(LIBS)

clean:
	rm -f mumpsdrv

# Mac
# MPICC=mpicc
# IDIR=/opt/local/include
# LDIR=/opt/local/lib
# LIBS=-ldmumps -lmumps_common
# CFLAGS=-O3 -I$(IDIR) -L$(LDIR)
# 
# 
# 
# 
# all: build
# 
# 
# build: mumpsdrv.c
#   $(MPICC) mumpsdrv.c -o mumpsdrv $(LIBS) $(CFLAGS)
