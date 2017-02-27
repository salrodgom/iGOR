LINKFLAGS_FOR=-O3 -march=native
COMP_FOR = gfortran
all:
	${COMP_FOR} ${LINKFLAGS_FOR} -c random.f95
	${COMP_FOR} ${LINKFLAGS_FOR} -c distance.f95
	${COMP_FOR} ${LINKFLAGS_FOR} -c cations.f95
	${COMP_FOR} ${LINKFLAGS_FOR} -c output.f95
	${COMP_FOR} ${LINKFLAGS_FOR} -c cell.f95
	${COMP_FOR} ${LINKFLAGS_FOR} -c volume.f95
	${COMP_FOR} ${LINKFLAGS_FOR} -c Si_Al.f95
	${COMP_FOR} -o si_al Si_Al.o cell.o distance.o cations.o random.o output.o volume.o
clean:
	rm *.o
