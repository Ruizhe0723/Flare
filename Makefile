# Fortran compiler
FC = mpif90
#FC = ifort

# Quick non-optimized compile
# FFLAGS = -g -fno-silent -ftrapping-math -ftrapv -Wimplicit -Wall

# Non-optimized compile with profiling
# FFLAGS = -pg -fno-silent -ftrapping-math -ftrapv -Wimplicit -Wall

# Optimized compile
# FFLAGS = -O3 -mcmodel=large -shared_intel
FFLAGS = -O3 -mcmodel=large -Wall


# Optimized compile with profiling
# FFLAGS = -pg -O3 -fno-silent -Wimplicit -Wall

# Source files
SRCDIR = src
OBJDIR = obj

VPATH = $(SRCDIR):$(OBJDIR)

OPTIONS1 = -fno-range-check -fcheck=all
OPTIONS2 = -J $(OBJDIR)

SRC =   riemannCZ.f pdfsSrcs.f func.f #interpLamFlame.f

OBJ = $(SRC:.f=.o)

# Include files
INC = integrate.inc data.inc

# Libraries: need separate compilation in their own subdirectories
# FFTs; random number generators
# MPI libraries
LIBS = -lmpi
#LIBDIRS = -L$$HOME/lib64 -Wl, -rpath, $$HOME/lib
# Object files are the same as the source files but with .o rather than .f
OBJS=$(SRC:.f=.o)

# Default compilation of .f files to .o
.f.o: $(FC) $(INC) -c $<
#.f.o: $(FC)  -c $<

%.o:%.f
	@mkdir -p $(OBJDIR)
	$(FC) $(FFLAGS) $(OPTIONS1) $(OPTIONS2) -c -o $(OBJDIR)/$@  $<

# Full make of riemannCZ
all: riemann

riemann: $(OBJS)
	$(FC) $(FFLAGS) -o riemann $(OBJDIR)/*.o $(LIBS)

# Add -pg for profiling
# riemann: $(OBJS)
#	$(FC) $(OBJS) -pg -o riemann51ZptsSerial $(LIBS)

# Command to clean up after
clean:
	rm -f *.o riemann fort.98 log fort.* *~ unit* $(OBJDIR)/*.o
