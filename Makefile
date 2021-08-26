# Fortran compiler
FC = /Users/huangruizhe/opt/usr/local/bin/mpif90
# FC = gfortran

# Quick non-optimized compile
# FFLAGS = -g -fno-silent -ftrapping-math -ftrapv -Wimplicit -Wall

# Non-optimized compile with profiling
# FFLAGS = -pg -fno-silent -ftrapping-math -ftrapv -Wimplicit -Wall

# Optimized compile
# FFLAGS = -O3 -mcmodel=large -shared_intel
FFLAGS = -O3 -Wall


# Optimized compile with profiling
# FFLAGS = -pg -O3 -fno-silent -Wimplicit -Wall

# Source files
SRCDIR = src
OBJDIR = obj
BINDIR = bin

EXE=flare

TARGET = $(BINDIR)/$(EXE)

VPATH = $(SRCDIR):$(OBJDIR)

OPTIONS1 = -fno-range-check -fcheck=all
OPTIONS2 = -J $(OBJDIR)

SRC = integrate.F90 func.F90 pdf.F90 flare.F90 # interpLamFlame.f

# Include files
# INC = integrate.inc data.inc

# Libraries: need separate compilation in their own subdirectories
# FFTs; random number generators
# MPI libraries
LIBS = -lmpi
#LIBDIRS = -L$$HOME/lib64 -Wl, -rpath, $$HOME/lib
# Object files are the same as the source files but with .o rather than .f
OBJS=$(SRC:.F90=.o)

# Default compilation of .f files to .o
# .f.o: $(FC) $(INC) -c $<
.F90.o: $(FC)  -c $<

%.o:%.F90
	@mkdir -p $(OBJDIR)
	$(FC) $(FFLAGS) $(OPTIONS1) $(OPTIONS2) -c -o $(OBJDIR)/$@  $<

# Full make of riemannCZ
all: flare

flare: $(OBJS)
	@mkdir -p $(BINDIR)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJDIR)/*.o $(LIBS)

# Add -pg for profiling
# riemann: $(OBJS)
#	$(FC) $(OBJS) -pg -o riemann51ZptsSerial $(LIBS)

# Command to clean up after
clean:
	rm -fv $(OBJDIR)/*.o $(OBJDIR)/*.mod $(TARGET)
