#!/bin/sh
#
# **********************
# *** MolSOC program ***
# **********************
#
# Creation: S. G. Chiodo
#
# *** Alpha OSF1 Fortran Compiler; f90 ***
#
#F90 = f90
#FLAG = -r8 -fast -check nopower
#LINK = -L../lib -ldmy
#LINK = -lcxml
#INCLUDE = ../include
#
# *** Intel Compiler;  efc ***
#
#F90 = ifort
#FLAG = -r8 
#LINK = -L../lib -ldmy
#LINK = -lmkl -lpthread -lguide
#INCLUDE = ../include
#
# *** PGI Fortran Compiler; pgf90 ***
#
#F90 = pgf90
#FLAG = -r8 -fast
#LINK = -L../lib -ldmy
#INCLUDE = ../include
#
F90 = gfortran
FLAG = -fdefault-real-8 
LINK = -L../lib -ldmy 
INCLUDE = ../include
#
PROGNAME = molsoc0.1.exe
#
SHELL = /bin/sh
CP = /bin/cp
MV = /bin/mv
RM = /bin/rm -f

OBJECTS = \
	altmos.o atomdata.o atomdis.o bldimat.o bldimso.o bldimom.o \
	bldover.o bldpfs0.o bldpfs2.o bldsoc.o bldsoc0.o bldsoc2.o \
	chkmult.o clebsch.o cputime.o copyr.o crftable.o demon.o \
	discone.o disczero.o errors.o fileio.o finele.o \
	finnel.o finkey.o finstr.o gammaf.o \
	gaussian.o genop.o getdate.o iniocc.o initialize.o ioscr.o \
	mathcon.o matrelem.o molsoc.o moment.o momesym.o \
	moocc.o mos.o motrsc.o muldens.o multipole.o norgtf.o \
	norsto.o orbord.o orbsym.o ortho.o ortsvdlap.o overabu.o \
	overfs.o overlap.o pointer.o readbas.o readinp.o readmos.o readshgf.o \
	rr1socin.o rr2socin.o rrmomint.o rrovint.o scangm.o scangmf.o \
	scangmt.o sdist.o shellop.o socone.o \
	socfst0.o socfst2.o socint.o sozeff.o \
	traatb.o transc.o turbomole.o turbord.o \
	kssoc.o ucc.o vecfso.o vecfso2.o vprod.f writout.o wrthead.o xyzcomp.o \
	rr21.o rr22.o rr23.o rr24.o

all:	s cleano utils

s:	$(OBJECTS)
	$(F90) $(FLAG) -o $(PROGNAME) $(OBJECTS) $(LINK)

clean:
	${RM} ${PROGNAME} *.o *~ core

.f.o:
	$(F90) $(FLAG) -I$(INCLUDE) -c $*.f

cleano:
	${RM} *.o *~ core

utils:
	${MV} molsoc0.1.exe ../bin
