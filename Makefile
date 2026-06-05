FC = mpif90
comma := ,
FFLAGS ?= -O2 -cpp -ffree-form
TARGET ?= gtc_local_num_mode13
MPICC ?= mpicc
MPICXX ?= mpicxx
LLVM_PREFIX ?= $(if $(wildcard /opt/anaconda3/bin/clang),/opt/anaconda3,)
MPI_CC_BACKEND ?= $(if $(LLVM_PREFIX),$(LLVM_PREFIX)/bin/clang,clang)
MPI_CXX_BACKEND ?= $(if $(wildcard $(LLVM_PREFIX)/bin/clang++),$(LLVM_PREFIX)/bin/clang++,$(if $(wildcard $(LLVM_PREFIX)/bin/clang++-20),$(LLVM_PREFIX)/bin/clang++-20,clang++))
GTC_CFLAGS ?= -O3 -march=native -std=c11 -Wall -Wextra
GTC_CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra
GTC_OBJCFLAGS ?= -fobjc-arc
GTC_CLIBS ?=
C_TARGET ?= gtc_c
POCKETFFT ?= 1
OPENMP ?= 1
METAL ?= 0
OPENMP_PREFIX ?= $(LLVM_PREFIX)
GTC_OMPFLAGS ?= -fopenmp $(if $(wildcard $(OPENMP_PREFIX)/include/omp.h),-I$(OPENMP_PREFIX)/include,)
GTC_OMPLIBS ?= $(if $(wildcard $(OPENMP_PREFIX)/lib/libomp.dylib),-L$(OPENMP_PREFIX)/lib -lomp -Wl$(comma)-rpath$(comma)$(OPENMP_PREFIX)/lib,-lomp)
GTC_METALLIBS ?= -framework Foundation -framework Metal

ifeq ($(OPENMP),1)
GTC_CFLAGS += $(GTC_OMPFLAGS) -DGTC_USE_OPENMP
GTC_CLIBS += $(GTC_OMPLIBS)
endif

ifeq ($(METAL),1)
GTC_CFLAGS += -DGTC_USE_METAL
GTC_CLIBS += $(GTC_METALLIBS)
OBJC_SRC = c/metal_gpu.m
else
OBJC_SRC =
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

C_BUILD_DIR ?= build/c-openmp$(OPENMP)-metal$(METAL)-pocketfft$(POCKETFFT)
ifeq ($(METAL),1)
C_OBJ = $(patsubst c/%.c,$(C_BUILD_DIR)/%.o,$(C_SRC))
CXX_OBJ = $(patsubst c/%.cpp,$(C_BUILD_DIR)/%.o,$(CXX_SRC))
OBJC_OBJ = $(patsubst c/%.m,$(C_BUILD_DIR)/%.o,$(OBJC_SRC))
else
C_OBJ = $(C_SRC:.c=.o)
CXX_OBJ = $(CXX_SRC:.cpp=.o)
OBJC_OBJ =
endif
CLEAN_C_TARGETS = $(sort $(C_TARGET) gtc_c gtc_c_metal)

OBJ = $(SRC:.F90=.o)
OBJ := $(OBJ:.f90=.o)

.PHONY: all c clean

all: $(TARGET)

c: $(C_TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

$(C_TARGET): $(C_OBJ) $(CXX_OBJ) $(OBJC_OBJ)
ifeq ($(POCKETFFT),1)
	MPICH_CXX=$(MPI_CXX_BACKEND) $(MPICXX) $(GTC_CXXFLAGS) -o $@ $(C_OBJ) $(CXX_OBJ) $(OBJC_OBJ) -lm $(GTC_CLIBS)
else
	MPICH_CC=$(MPI_CC_BACKEND) $(MPICC) $(GTC_CFLAGS) -o $@ $(C_OBJ) $(OBJC_OBJ) -lm $(GTC_CLIBS)
endif

c/%.o: c/%.c c/gtc.h
	MPICH_CC=$(MPI_CC_BACKEND) $(MPICC) $(GTC_CFLAGS) -c -o $@ $<

c/%.o: c/%.cpp c/gtc.h c/pocketfft_hdronly.h
	MPICH_CXX=$(MPI_CXX_BACKEND) $(MPICXX) $(GTC_CXXFLAGS) -c -o $@ $<

$(C_BUILD_DIR):
	mkdir -p $@

$(C_BUILD_DIR)/%.o: c/%.c c/gtc.h | $(C_BUILD_DIR)
	MPICH_CC=$(MPI_CC_BACKEND) $(MPICC) $(GTC_CFLAGS) -c -o $@ $<

$(C_BUILD_DIR)/%.o: c/%.cpp c/gtc.h c/pocketfft_hdronly.h | $(C_BUILD_DIR)
	MPICH_CXX=$(MPI_CXX_BACKEND) $(MPICXX) $(GTC_CXXFLAGS) -c -o $@ $<

$(C_BUILD_DIR)/%.o: c/%.m c/gtc.h | $(C_BUILD_DIR)
	MPICH_CC=$(MPI_CC_BACKEND) $(MPICC) $(GTC_CFLAGS) $(GTC_OBJCFLAGS) -c -o $@ $<

%.o: %.F90 module.o
	$(FC) $(FFLAGS) -c $<

%.o: %.f90 module.o
	$(FC) $(FFLAGS) -c $<

module.o: module.F90
	$(FC) $(FFLAGS) -c $<

setup.o: setup.F90 module.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(TARGET) $(CLEAN_C_TARGETS) $(OBJ) $(C_OBJ) $(CXX_OBJ) c/*.o *.mod
	rm -rf build
