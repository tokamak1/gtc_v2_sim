FC = mpif90
FFLAGS ?= -O2 -cpp -ffree-form
TARGET ?= gtc_local_num_mode13

SRC = module.F90 setup.F90 ran_num_gen.F90 set_random_values.f90 \
      function.F90 load.F90 restart.F90 diagnosis.F90 snapshot.F90 \
      chargei.F90 poisson.F90 smooth.F90 field.F90 pushi.F90 shifti.F90 \
      fft_gl.F90 tracking.F90 main.F90

OBJ = $(SRC:.F90=.o)
OBJ := $(OBJ:.f90=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

%.o: %.F90 module.o
	$(FC) $(FFLAGS) -c $<

%.o: %.f90 module.o
	$(FC) $(FFLAGS) -c $<

module.o: module.F90
	$(FC) $(FFLAGS) -c $<

setup.o: setup.F90 module.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(TARGET) $(OBJ) *.mod
