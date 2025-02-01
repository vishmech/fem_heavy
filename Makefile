# Compiler and flags
FC = gfortran
FLAGS = -fbounds-check
LIBS = -L/usr/local/pgplot -L/usr/lib/X11 -lpgplot -lX11

# Source files, objects, and modules
SRCS = types.f90 fedata.f90 plane42.f90 link1.f90 processor.f90 numeth.f90 fea.f90 main.f90
OBJS = $(SRCS:.f90=.o)
EXEC = FEA_RESULTS

# Default target: Build the executable
all: $(EXEC)

# Link object files to create the executable
$(EXEC): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBS) $(FLAGS)

# Rule to compile .f90 files into .o files and generate .mod files
%.o: %.f90
	$(FC) -c $< $(FLAGS)

# Clean up object files, modules, and executable
clean:
	rm -f $(OBJS) *.mod $(EXEC)

# Rebuild everything from scratch
rebuild: clean all
