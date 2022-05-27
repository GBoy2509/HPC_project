# Makefile for compiling with .c .h
# To compile and link write simply "make"
# make: normal compiling
# make clean: remove object and module files along with error logs

#List of files needed to compile
OBJ = heat_trans.o
MYDIR = ./
FCOMPILER = gcc -c -std=c99 -O3

#Compile object files
heat_trans:	$(OBJ)
		gcc -O3 -o heat_trans $(OBJ) -lm

heat_trans.o:	$(MYDIR)heat_trans.c global.h
		$(FCOMPILER) $(MYDIR)heat_trans.c

#eulersch.o:	$(MYDIR)eulersch.c global.h eulersch.h
#		$(FCOMPILER) $(MYDIR)eulersch.c
		
#Remove Old Object files, useful when recompiling
.PHONY: clean
clean:
	rm -rf *.o heat_trans

