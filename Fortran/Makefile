SRC=MD.f control.f util.f
OBJ=$(SRC:.f=.o)

FFLAGS= -g
FC=pgf90

all: MD

MD: $(OBJ)
	$(FC) $(FFLAGS)  -o $@  $(OBJ)

clean:
	rm -f MD $(OBJ)

$(OBJ): coord.inc Makefile

MD.o: MD.f
	$(FC) $(FFLAGS) -c MD.f

control.o: control.f
	$(FC) $(FFLAGS) -c control.f

util.o: util.f
	$(FC) $(FFLAGS) -c util.f




