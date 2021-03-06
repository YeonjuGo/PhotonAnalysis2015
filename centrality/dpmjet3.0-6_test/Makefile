########################################################################
#                                                                      #
#                      makefile for DPMJET 3.0                         #
#                                                                      #
#                      (S. Roesler 21-10-2008)                         #
#                                                                      #
########################################################################

# without FLUKA evaporation
USER   = user3.0-6
DPMJET = dpmjet3.0-6

# with FLUKA evaporation
# USER   = user3.0-5F
# DPMJET = dpmjet3.0-5F

PHOJET = phojet1.12-35c4
PYTHIA = pythia6115dpm3v1

EXE = $(DPMJET)

######################### compiler - options ###########################

#  DEC (Alfa RISC)
#  ---------------
# OPT  = -C -v -fpe4 -check underflow -O2
# FC = f77
# LD = f77

#  HP (series 700)
#  ---------------
# OPT = -K +T +E1 +ppu -O
# FC = f77
# LD = f77

#  Linux g77
#  ---------
# OPT = -C -Wall -malign-double -mpentium -O
#OPT = -c -O3 -g -mcpu=pentiumpro -fexpensive-optimizations -funroll-loops -fstrength-reduce -Wall -fno-automatic -fno-silent -m32 -ffortran-bounds-check -v
#OPT = -c -C -std=legacy -O  -O3 -g -fexpensive-optimizations -funroll-loops -fno-automatic -fbounds-check -v -fPIC
#FC = g77
#LD = g77
#TRAP = fpe.o

#  Linux gfortran
#  --------------
# OPT = -c -C -std=legacy -O  -O3 -g -fexpensive-optimizations -funroll-loops -fno-automatic -m32 -fbounds-check -v
 OPT = -c -C -std=legacy -O  -O3 -g -fexpensive-optimizations -funroll-loops -fno-automatic -fbounds-check -v -fPIC
 FC = gfortran
 LD = gfortran
 TRAP = fpe.o

#  IBM (RISC)
#  ----------
# OPT = -C -v -qextchk -O2
# FC = xlf
# LD = xlf

############################# libraries ################################

# LIB = -L/cern/pro/lib -lpdflib -lmathlib -lpacklib -lkernlib
LIB = -L $(FLUPRO) -lflukahp
INC = -I$(FLUPRO)/flukapro

########################################################################

all:  $(EXE)

clean:
	rm -f $(EXE) $(USER).o $(DPMJET).o $(PHOJET).o $(PYTHIA).o \
              *.prj *.chk core

$(EXE): $(USER).o $(DPMJET).o $(PHOJET).o $(PYTHIA).o Makefile

	$(LD) -o $(EXE) $(USER).o $(DPMJET).o $(PHOJET).o $(PYTHIA).o \
                        $(LIB) $(TRAP)

.f.o:
	$(FC) -c $(OPT) $(INC) $<
