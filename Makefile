FC = mpif90
FFLAGS ?= -O2 -cpp -ffree-form
TARGET ?= gtc_local_num_mode13
MPICC ?= mpicc
MPICXX ?= mpicxx
MPI_CC_BACKEND ?= clang
MPI_CXX_BACKEND ?= clang++
GTC_CFLAGS ?= -O3 -march=native -std=c11 -Wall -Wextra
GTC_CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra
GTC_CLIBS ?=
C_TARGET ?= gtc_c
POCKETFFT ?= 1
OPENMP ?= 1
GTC_OMPFLAGS ?= -fopenmp
GTC_OMPLIBS ?= -lomp

ifeq ($(OPENMP),1)
GTC_CFLAGS += $(GTC_OMPFLAGS) -DGTC_USE_OPENMP
GTC_CLIBS += $(GTC_OMPLIBS)
endif

SRC = module.F90 setup.F90 ran_num_gen.F90 set_random_values.f90 \
      function.F90 load.F90 restart.F90 diagnosis.F90 snapshot.F90 \
      chargei.F90 poisson.F90 smooth.F90 field.F90 pushi.F90 shifti.F90 \
      fft_gl.F90 tracking.F90 main.F90

C_SRC = c/module.c c/setup.c c/ran_num_gen.c c/set_random_values.c \
        c/function.c c/load.c c/restart.c c/diagnosis.c c/snapshot.c \
        c/chargei.c c/poisson.c c/smooth.c c/field.c c/pushi.c c/shifti.c \
        c/fft_gl.c c/main.c
ifeq ($(POCKETFFT),1)
GTC_CFLAGS += -DGTC_USE_POCKETFFT
CXX_SRC = c/fft_pocketfft.cpp
else
CXX_SRC =
endif
C_OBJ = $(C_SRC:.c=.o)
CXX_OBJ = $(CXX_SRC:.cpp=.o)

OBJ = $(SRC:.F90=.o)
OBJ := $(OBJ:.f90=.o)

.PHONY: all c clean

all: $(TARGET)

c: $(C_TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

$(C_TARGET): $(C_OBJ) $(CXX_OBJ)
ifeq ($(POCKETFFT),1)
	MPICH_CXX=$(MPI_CXX_BACKEND) $(MPICXX) $(GTC_CXXFLAGS) -o $@ $(C_OBJ) $(CXX_OBJ) -lm $(GTC_CLIBS)
else
	MPICH_CC=$(MPI_CC_BACKEND) $(MPICC) $(GTC_CFLAGS) -o $@ $(C_OBJ) -lm $(GTC_CLIBS)
endif

c/%.o: c/%.c c/gtc.h
	MPICH_CC=$(MPI_CC_BACKEND) $(MPICC) $(GTC_CFLAGS) -c -o $@ $<

c/%.o: c/%.cpp c/gtc.h c/pocketfft_hdronly.h
	MPICH_CXX=$(MPI_CXX_BACKEND) $(MPICXX) $(GTC_CXXFLAGS) -c -o $@ $<

%.o: %.F90 module.o
	$(FC) $(FFLAGS) -c $<

%.o: %.f90 module.o
	$(FC) $(FFLAGS) -c $<

module.o: module.F90
	$(FC) $(FFLAGS) -c $<

setup.o: setup.F90 module.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(TARGET) $(C_TARGET) $(OBJ) $(C_OBJ) $(CXX_OBJ) *.mod
