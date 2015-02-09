all: build


build: mumpsdrv.c
	mpicc -O3 mumpsdrv.c -o mumpsdrv -I /opt/local/include -L /opt/local/lib -ldmumps
