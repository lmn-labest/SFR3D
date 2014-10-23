#!/bin/sh
PATH_INCLUDE="include"    
PATH_LIB="lib"    
PRENAME=mfvCell
CC=gcc
OPENMP=no
OT=O3
DEBUG=no

#------------------gerando o nome do excutavel-------------
ifeq ($(CC),icc)
  ifeq ($(OPENMP),yes)
    COMPILER_NAME=intel_omp_$(OT)
  else
    COMPILER_NAME=intel_$(OT)
  endif  
endif

ifeq ($(CC),gcc)
  ifeq ($(OPENMP),yes)
    COMPILER_NAME=gnu_omp_$(OT)
  else
    COMPILER_NAME=gnu_$(OT)
  endif  
endif
  
NAME+=$(PRENAME)_$(COMPILER_NAME)
#-------------------Fontes--------------------------------
fontes = \
src/Adjcency.c\
src/CellLoop.c\
src/CellLib.c\
src/Csr.c\
src/Datastruct.c\
src/Debug.c\
src/File.c\
src/Graph.c\
src/HccaBlas.c\
src/Main.c\
src/Memoria.c\
src/Numeq.c\
src/Rcm.c\
src/ReadFile.c\
src/Reord.c\
src/Vtk.c\
src/WriteVtk.c\
src/WriteMtx.c
#-------------------Flags necessarios--------------------------------
NFLAGS=-I$(PATH_INCLUDE) -L$(PATH_LIB) -D_MMIO_  -D_DEBUG_GEOM \
        -fopt-info-optimized-missed=logOpt.txt 
LDFLAGS=-lmmio -lmetis-x86_64 
#--------------------compiladores------------------------------------
# intel ifort
ifeq ($(CC),icc)
  LDFLAGS += 
  OFLAGS  +=  
  ifeq ($(OPENMP),yes)
    OFLAGS  += -openmp
  endif
endif
# gnu gcc
ifeq ($(CC),gcc)
  LDFLAGS +=  -lm 
  OFLAGS  +=  -Wall -ansi -std=c99 -pedantic-errors 
  ifeq ($(OPENMP),yes)
    OFLAGS  += -fopenmp
  endif
endif
#--------------------------------------------------------------------
#---------------------------Debug------------------------------------
ifeq ($(DEBUG),yes)
  OFLAGS += -g -$(OT)	
else
  OFLAGS += -$(OT) 
endif
#--------------------------------------------------------------------
CFLAGS= $(NFLAGS) $(OFLAGS) 


.SUFFIXES: 
.SUFFIXES: .c .h .o
OBJS= $(fontes:%.c=%.o)
build:	$(OBJS) 
	ls bin || mkdir -p bin
	$(CC) $(CFLAGS) $(OBJS) -o bin/$(NAME)  $(LDFLAGS)

tags:
	ctags -R src/*.c include/*.h

.PHONY: cleantags
cleantags:
	@rm -fv tags
	
.PHONY: clean
clean:  
	@rm -fv src/*.o
	@rm -fv ../src/*.o
	@rm -fv bin/$(NAME)

.PHONY: cleanall
cleanall:  
	@rm -fv tags
	@rm -fv src/*.o
	@rm -fv ../src/*.o
	@rm -fv bin/$(NAME)


.PHONY: help
help:
	@echo "Autor :$(AUTHOR)                              "
	@echo "Makefile para prepar para sitemas linux.      "
	@echo -e "\E[7;32mOpcoes:\E[1;0m                      "
	@echo "build         - compila o prepar              "
	@echo "build_modules - gera os modulos               "
	@echo "tags          - gera os tags                  "
	@echo "cleantags     - limpa os tags                 "
	@echo "clean         - limpa os obj, bin e mod       "
	@echo "cleaall       - limpa tudo obj,bin,mod e tags "

# DO NOT DELETE

# DO NOT DELETE
