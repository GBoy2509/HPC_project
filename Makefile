# Makefile for compiling with .c .h
# To compile and link write simply "make"
# make: normal compiling
# make clean: remove object and module files along with error logs

#List of files needed to compile
OBJ = heat_trans.o
MYDIR = ./
# FCOMPILER = icc -c -std=c99 -O3 -I/share/base/hdf5/1.10.4-intel-18.4/include -L/share/base/hdf5/1.10.4-intel-18.4/lib -lhdf5

FCOMPILER = icc -c -std=c99 -O3 -lhdf5

#Compile object files
heat_trans:	$(OBJ)
		icc -O3 -o heat_trans $(OBJ) -lm -lhdf5

heat_trans.o:	$(MYDIR)heat_trans.c global.h
		$(FCOMPILER) $(MYDIR)heat_trans.c

#eulersch.o:	$(MYDIR)eulersch.c eulersch.h global.h
#		$(FCOMPILER) $(MYDIR)eulersch.c
		
#Remove Old Object files, useful when recompiling
.PHONY: clean
clean:
	rm -rf *.o *.h5 *.err *.out heat_trans

