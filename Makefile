TOPSRCDIR = /home/ali/Downloads/lapack-3.11.0
include $(TOPSRCDIR)/make.inc

FC = gfortran
FFLAGS = -g -Wall -O2 $(OPTS) 

MOD_SRCS = mod.f90 grid.f90 solver.f90
MAIN_SRC = main.f90
CHECKGRID_SRC = check_grid.f90

MOD_OBJS = $(MOD_SRCS:.f90=.o)
SOLVER_MAIN_OBJ = $(MAIN_SRC:.f90=.o)
CHECKGRID_OBJ = $(CHECKGRID_SRC:.f90=.o)

SOLVER_EXE = cavity_solver
CHECKGRID_EXE = check_grid

.PHONY: all
all: $(SOLVER_EXE)

$(SOLVER_EXE): $(MOD_OBJS) $(SOLVER_MAIN_OBJ)
	$(FC) $(FFLAGS) $^ -o $@ $(LAPACKLIB) $(BLASLIB)

.PHONY: check_grid
check_grid: $(CHECKGRID_EXE)
$(CHECKGRID_EXE): $(MOD_OBJS) $(CHECKGRID_OBJ)
	$(FC) $(FFLAGS) $^ -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o *.mod $(SOLVER_EXE) $(CHECKGRID_EXE) 
