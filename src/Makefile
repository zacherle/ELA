PDIR:=$(HOME)/lab/ela
BIN_DIR:=.
OBJ_DIR:=$(PDIR)/obj
LIBSRC_DIR:=$(PDIR)/libsrc

FC:=gfortran
CC:=g++

CFLAGS=-g -O0 -Wall
FFLAGS=-g -O0 -Wall  # -ffunction-sections -Wl,--gc-sections
LDFLAGS = -L $(PDIR)/lib 
CPPFLAGS= -I $(PDIR)/include
LIBS= -llapack -lblas
LIB=$(LIBSRC_DIR)/libela.a

.PHONY: all clean

APP_NAMES:=ela
APPS:=$(APP_NAMES:%=$(BIN_DIR)/%)

#LDFLAGS += -fcheck=bounds,pointer

all: $(APPS)

clean:
	rm *.o $(APPS)

$(BIN_DIR)/ela: ela_main.o ela_d.o loc_hypo.o $(LIB)
	$(CC)  $(LDFLAGS) $^ -o $@ -Wl,-Bstatic -lboost_program_options -lminpack -Wl,-Bdynamic -lgfortran $(LIBS)

#lmder2.o: lmder2.for
#	$(FC) $(FFLAGS) -c $<

#% %.o:
#	$(CC)  $(LDFLAGS) $^ -o $@ $(LIBS)

#%.o: %.f90
#	$(FC) -cpp $(FFLAGS)  -c $< -o $@ $(CPPFLAGS) -I$(LIBSRC_DIR)

%.o: %.cpp
	$(CC) -cpp $(CFLAGS)  -c $< -o $@ $(CPPFLAGS) -I$(LIBSRC_DIR)


ela_d.o: ../include/param.h ../include/calc_covar.h ../include/arrivals.h
ela_d.o: ../include/hypfile.h ../include/stations.h ../include/fw_tdx.h
ela_d.o: ../include/gather.h
ela_main.o: ../include/stations.h ../include/layers.h ../include/hypfile.h
ela_main.o: ../include/arrivals.h
loc_hypo.o: ../include/fw_tdx.h
