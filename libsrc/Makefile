LIB:=libela.a

FC:=gfortran
CC:=gcc
CPP:=g++
AR:=ar


CFLAGS=-g -O0 -Wall
FFLAGS=-g -O0 -Wall # -ffunction-sections -Wl,--gc-sections
CPPFLAGS= -I../include
#CPPFLAGS= -I$(HOME)/lab/hypoela/lib -I$(HOME)/prj/datetime-fortran/include -I../include
    # -I $(HOME)/prj/minpack/src

.PHONY: all clean

all: $(LIB)

clean:
	rm *.o *.mod $(LIB)

#c_f_elalib.f90 gap.c
SRC:= hypfile.cpp \
    arrivals.cpp \
    stations.cpp \
    magni.cpp \
    wgsjtsk.cpp \
    rt1Dlayers.f90 \
    calc_covar.f90 \
    layers.cpp \
    gather.cpp \
    fw_tdx.cpp \
    output_list.cpp \
    twopoint.cpp \
    hy3file.c

OBJ_2 = $(SRC:.c=.o)
OBJ_1 = $(OBJ_2:.cpp=.o)
OBJ   = $(OBJ_1:.f90=.o)

do:
	echo $(OBJ)

%.o: %.f90
	$(FC) -cpp $(FFLAGS)  -c $< -o $@ $(CPPFLAGS)

%.o: %.c
	$(CC) -cpp $(CFLAGS)  -c $< -o $@ $(CPPFLAGS)

%.o: %.cpp
	$(CPP) -cpp $(CFLAGS)  -c $< -o $@ $(CPPFLAGS)

#makedepf90 *.f90 > libela_depend_1.lst
#-include libela_depend.lst

#xy2fl.o: xy2fl.for
#	$(FC) $(FFLAGS) -c -ffixed-line-length-170 -fd-lines-as-comments xy2fl.for

$(LIB): $(OBJ)
	$(AR) rcs $(LIB) $(OBJ)

#

arrivals.o: ../include/hypfile.h ../include/datetime.hpp
arrivals.o: ../include/stations.h ../include/arrivals.h
fw_tdx.o: ../include/gather.h
gather.o: ../include/arrivals.h ../include/hypfile.h ../include/stations.h
gather.o: ../include/twopoint.h ../include/gather.h
hypfile.o: ../include/hypfile.h
layers.o: ../include/layers.h
layers1.o: ../include/layers.h
magni.o: ../include/arrivals.h ../include/hypfile.h ../include/stations.h
output_list.o: ../include/version.h ../include/arrivals.h
output_list.o: ../include/hypfile.h ../include/stations.h
output_list.o: ../include/const_raddeg.h ../include/gather.h
output_list.o: ../include/wgsjtsk.h ../include/hy3file.h /usr/include/time.h
output_list.o: ../include/magni.h
stations.o: ../include/stations.h
twopoint.o: ../include/rt1Dlayers.h ../include/layers.h
rt1Dlayers.o : rt1Dlayers.f90 
calc_covar.o : calc_covar.f90 
