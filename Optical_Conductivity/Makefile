FC = mpifort

SRC_DIR = src
OBJ_DIR = build

BUILD ?= release

FFLAGS_COMMON = -std=f2008 -mcpu=apple-m2 -ftree-vectorize -funroll-loops -fno-math-errno
LDFLAGS = -framework Accelerate

ifeq ($(BUILD),debug)
  FFLAGS = -O0 -g -Wall -Wextra -fcheck=all -fbacktrace -fbounds-check $(FFLAGS_COMMON)
else
  FFLAGS = -Ofast $(FFLAGS_COMMON)
  # FFLAGS  += -flto
  # LDFLAGS += -flto
endif

# --- source files (WITHOUT path)
SRCS = math_utils.f90 mpi_utils.f90 qc_data.f90 integration_vec.f90 \
       write_stuff.f90 qc_farm.f90 qc_cond.f90 qc_mod.f90 qc_main.f90

# --- derived object list
OBJS = $(SRCS:%.f90=$(OBJ_DIR)/%.o)

all: qc_project

qc_project: $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

# --- pattern rule: src/*.f90 -> build/*.o
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# --- ensure build dir exists
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	rm -rf $(OBJ_DIR) *.mod qc_project
