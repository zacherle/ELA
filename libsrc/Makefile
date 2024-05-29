LIB:=libela.a

FC:=gfortran
CC:=gcc
AR:=ar


CFLAGS=-g -O0 -Wall
FFLAGS=-g -O0 -Wall # -ffunction-sections -Wl,--gc-sections
CPPFLAGS= -I$(HOME)/lab/hypoela/lib -I$(HOME)/prj/datetime-fortran/include -I../include
    # -I $(HOME)/prj/minpack/src

.PHONY: all clean

all: $(LIB)

clean:
	rm *.o *.mod $(LIB)

SRC:= c_f_elalib.f90 \
    version.f90 \
    hypfile.f90 \
    arrivals.f90 \
    stations.f90 \
    const_raddeg.f90 \
    magni.f90 \
    xy2fl.for \
    gap.c \
    rt1Dlayers.f90 \
    calc_covar.f90 \
    layers.f90 \
    gather.f90 \
    fw_tdx.f90 \
    output_list.f90 \
    twopoint.f90 \



OBJ_2 = $(SRC:.c=.o)
OBJ_1 = $(OBJ_2:.for=.o)
OBJ   = $(OBJ_1:.f90=.o)

do:
	echo $(OBJ)

%.o: %.f90
	$(FC) -cpp $(FFLAGS)  -c $< -o $@ $(CPPFLAGS)

%.o: %.c
	$(CC) -cpp $(CFLAGS)  -c $< -o $@ $(CPPFLAGS)


#makedepf90 *.f90 > libela_depend_1.lst
-include libela_depend.lst

xy2fl.o: xy2fl.for
	$(FC) $(FFLAGS) -c -ffixed-line-length-170 -fd-lines-as-comments xy2fl.for

$(LIB): $(OBJ)
	$(AR) rcs $(LIB) $(OBJ)
