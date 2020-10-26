PROG  = soc_td
SRCS = soc_td.f90 init_prep.f90 data_input.f90 prec_mod.f90 soc_mo.f90 dip_mo.f90 basis_match.f90 norm_mo.f90 soc_state.f90 trans_dipole.f90

OBJS = soc_td.o init_prep.o data_input.o prec_mod.o soc_mo.o dip_mo.o basis_match.o norm_mo.o soc_state.o trans_dipole.o

LIBDIR = /opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64
INCLUDEDIR = 
LIBS = -Wl,--start-group $(LIBDIR)/libmkl_intel_ilp64.a\
       $(LIBDIR)/libmkl_core.a $(LIBDIR)/libmkl_sequential.a\
       -Wl,--end-group\
       -lpthread -lm
LIBS =

F90 =ifort
F90FLAGS = -O -heap-arrays 1024

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90
.f90.o:
	$(F90) $(F90FLAGS) -c $<

#soc_td.o: prec_mod.o data_input.o init_prep.o
soc_td.o:  init_prep.o soc_mo.o dip_mo.o basis_match.o norm_mo.o soc_state.o trans_dipole.o
data_input.o: prec_mod.o
init_prep.o: data_input.o prec_mod.o
soc_mo.o: prec_mod.o data_input.o init_prep.o 
dip_mo.o: prec_mod.o data_input.o init_prep.o 
basis_match.o: prec_mod.o data_input.o init_prep.o 
norm_mo.o: prec_mod.o data_input.o init_prep.o 
soc_state.o: prec_mod.o data_input.o init_prep.o soc_mo.o
trans_dipole.o: prec_mod.o data_input.o init_prep.o dip_mo.o


