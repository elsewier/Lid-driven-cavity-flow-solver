# Path to your LAPACK installation source directory
TOPSRCDIR = /home/ali/Downloads/lapack-3.11.0
include $(TOPSRCDIR)/make.inc

# Compiler and flags
FC = gfortran
FFLAGS = -g -Wall -O2 $(OPTS) 

# --- SOURCE FILES ---
MOD_SRCS = mod.f90 grid.f90 solver.f90 
MAIN_SRC = main.f90
CHECKGRID_SRC = check_grid.f90
CHECKDERIV_SRC = check_derivatives.f90 # Added new derivative checker source

# --- OBJECT FILES ---
MOD_OBJS = $(MOD_SRCS:.f90=.o)
SOLVER_MAIN_OBJ = $(MAIN_SRC:.f90=.o)
CHECKGRID_OBJ = $(CHECKGRID_SRC:.f90=.o)
CHECKDERIV_OBJ = $(CHECKDERIV_SRC:.f90=.o) # Added new object file

# --- EXECUTABLE NAMES ---
SOLVER_EXE = cavity_solver
CHECKGRID_EXE = check_grid
CHECKDERIV_EXE = check_derivatives # Added new executable name

# Default target: build the main solver
.PHONY: all
all: $(SOLVER_EXE)

# Rule to build the main solver executable
$(SOLVER_EXE): $(MOD_OBJS) $(SOLVER_MAIN_OBJ)
	$(FC) $(FFLAGS) $^ -o $@ $(LAPACKLIB) $(BLASLIB)

# Rule to build the grid checking utility
.PHONY: check_grid
check_grid: $(CHECKGRID_EXE)
$(CHECKGRID_EXE): $(MOD_OBJS) $(CHECKGRID_OBJ)
	$(FC) $(FFLAGS) $^ -o $@

# --- NEW RULE: Build the derivative checking utility ---
# This requires LAPACK and BLAS for the linear solve (DGESV)
.PHONY: check_derivatives
check_derivatives: $(CHECKDERIV_EXE)
$(CHECKDERIV_EXE): $(MOD_OBJS) $(CHECKDERIV_OBJ)
	$(FC) $(FFLAGS) $^ -o $@ $(LAPACKLIB) $(BLASLIB)

# Generic rule to compile .f90 source files to .o object files
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Rule to clean up build files
.PHONY: clean
clean:
	rm -f *.o *.mod $(SOLVER_EXE) $(CHECKGRID_EXE) $(CHECKDERIV_EXE)
