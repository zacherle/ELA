PDIR:=$(HOME)/lab/ela
BIN_DIR:=.
OBJ_DIR:=$(PDIR)/obj
LIBSRC_DIR:=$(PDIR)/libsrc

FC:=gfortran
CC:=g++

CFLAGS=-g -O0 -Wall
FFLAGS=-g -O0 -Wall  # -ffunction-sections -Wl,--gc-sections
LDFLAGS = -L $(PDIR)/lib -L $(HOME)/prj/datetime-fortran/lib
CPPFLAGS= -I $(PDIR)/lib -I $(HOME)/prj/datetime-fortran/include # -I $(HOME)/prj/minpack/src
LIBS= -ldatetime -llapack -lblas
LIB=$(LIBSRC_DIR)/libela.a

.PHONY: all clean

APP_NAMES:=ela
APPS:=$(APP_NAMES:%=$(BIN_DIR)/%)

#LDFLAGS += -fcheck=bounds,pointer

all: $(APPS)

clean:
	rm *.o $(APPS)

$(BIN_DIR)/ela: ela_main.o ela_d.o loc_hypo.o lmder2.o  $(LIB)
	$(CC)  $(LDFLAGS) $^ -o $@ -lboost_program_options  -lgfortran -Wl,-Bstatic -lminpack -Wl,-Bdynamic $(LIBS)

lmder2.o: lmder2.for
	$(FC) $(FFLAGS) -c $<

#% %.o:
#	$(CC)  $(LDFLAGS) $^ -o $@ $(LIBS)

%.o: %.f90
	$(FC) -cpp $(FFLAGS)  -c $< -o $@ $(CPPFLAGS) -I$(LIBSRC_DIR)

%.o: %.cpp
	$(CC) -cpp $(CFLAGS)  -c $< -o $@ $(CPPFLAGS) -I$(LIBSRC_DIR)

# makedepf90 *.f90 *.for *.cpp > ela_depend.lst
-include ela_depend.lst
