
FC := gfortran

FFLAGS := -g -Wall -O2 


MOD_SRCS := mod.f90 grid.f90 solver.f90

MAIN_SRC := main.f90 

# object files
MOD_OBJS := $(MOD_SRCS:.f90=.o)
SOLVER_MAIN_OBJ := $(MAIN_SRC:.f90=.o)


# exec name 
SOLVER_EXE := cavity_solver


.PHONY: all 
all: $(SOLVER_EXE)

$(SOLVER_EXE): $(MOD_OBJS) $(SOLVER_MAIN_OBJ)
	$(FC) $(FFLAGS) $^ -o $@


%.o: %.f90 
	$(FC) $(FFLAGS) -c $< -o $@

.PHONY: clean
clean: 
	rm -f *.o *.mod $(SOLVER_EXE)



