FC = g++
LD = g++
#FFLAGS = -g -O0
#FFLAGS = -O3
FFLAGS = -Ofast
LDFLAGS += -lgfortran -I #-l
#FFLAGS += -fallow-argument-mismatch
#-----------------------------------------------------------------#
SRCDIR = ../src
OBJDIR = ../obj
RUNDIR = ../run
EXTDIR = /home/numazawa/Solver
TARGET = ./a.out
OBJ = main.o vis.o vis_node.o amr.o smac.o iflow.o
#MKL--------------------------------------------------------------#
BLASLAPACK   = $(EXTDIR)/AGMG/lapack-3.4.2/liblapack.a
BLASLAPACK  += $(EXTDIR)/AGMG/lapack-3.4.2/librefblas.a
MKLlib     = $(BLASLAPACK)
#AGMG-GNU-----------------------------------------------------#
AGMGlib    +=   $(EXTDIR)/AGMG/AGMG_3.2.4-aca/SRC_seq/agmg_s.a
AGMGobj     =   $(EXTDIR)/AGMG/AGMG_3.2.4-aca/SRC_seq/dagmg.o
AGMGobj    +=   $(EXTDIR)/AGMG/AGMG_3.2.4-aca/SRC_seq/dagmg_mumps.o
LDFLAGS    +=   $(AGMGlib) $(MKLlib) 
#AMGS-GNU-----------------------------------------------------#
AMGSlib     =   $(EXTDIR)/AMGS/amgs-s-1.1/amgs_s.a
LDFLAGS    +=   $(AMGSlib) -fopenmp
#-----------------------------------------------------------------#

all: $(OBJ)  
	$(LD) -o $(TARGET) $(AGMGobj) $(OBJ) $(LDFLAGS) 
main.o: $(SRCDIR)/main.cpp  
	$(FC) $(CPPFLAGS) $(FFLAGS) -c  $(SRCDIR)/main.cpp 
vis.o: $(SRCDIR)/vis.cpp 
	$(FC) $(CPPFLAGS) $(FFLAGS) -c  $(SRCDIR)/vis.cpp 
vis_node.o: $(SRCDIR)/vis_node.cpp 
	$(FC) $(CPPFLAGS) $(FFLAGS) -c  $(SRCDIR)/vis_node.cpp 
amr.o: $(SRCDIR)/amr.cpp 
	$(FC) $(CPPFLAGS) $(FFLAGS) -c  $(SRCDIR)/amr.cpp 
smac.o: $(SRCDIR)/smac.cpp 
	$(FC) $(CPPFLAGS) $(FFLAGS) -c  $(SRCDIR)/smac.cpp  
iflow.o: $(SRCDIR)/iflow.cpp 
	$(FC) $(CPPFLAGS) $(FFLAGS) -c  $(SRCDIR)/iflow.cpp 
clean : 
	rm  -f $(OBJ) *.mod $(TARGET)
