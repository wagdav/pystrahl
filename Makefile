NCDFLIB  = $(HOME)/strahl/source/libnetcdf/lib
NCDFINC  = $(HOME)/strahl/source/libnetcdf/include

LIBFLAGS = -L$(NCDFLIB) -I$(NCDFINC)

#FC  =  f77 
FC  =  ifort
FFLAGS  = $(LIBFLAGS)
FOPTIMFLAGS=-g



STRAHL_PATH=./strahl

OBJstrahl = 	${STRAHL_PATH}/time_steps.o\
		${STRAHL_PATH}/read_parameter.o\
	    	${STRAHL_PATH}/neutrals.o\
	     	${STRAHL_PATH}/save_param.o\
		${STRAHL_PATH}/atomic_data.o\
	       	${STRAHL_PATH}/plad.o\
	       	${STRAHL_PATH}/read_atomdat.o\
		${STRAHL_PATH}/emissiv.o \
		${STRAHL_PATH}/math.o\
		${STRAHL_PATH}/num_recip.o\
	       	${STRAHL_PATH}/strahl_util.o\
	       	${STRAHL_PATH}/netcdfutil.o\
		${STRAHL_PATH}/impden.o\
	       	${STRAHL_PATH}/neo_trans_fsm.o\
	       	${STRAHL_PATH}/neo_trans_atp.o\
	       	${STRAHL_PATH}/neo_trans.o\
	       	${STRAHL_PATH}/neoart.o\
		${STRAHL_PATH}/saw_mix.o\
	       	${STRAHL_PATH}/read_grid.o\
		${STRAHL_PATH}/radiation_diag.o\
	       	${STRAHL_PATH}/radiation_total.o\
	       	${STRAHL_PATH}/save_dist.o\
		${STRAHL_PATH}/strahl.o


pystrahl:
	f2py -c wrap_strahl.pyf --fcompiler=intelem -L${STRAHL_PATH}\
	    ${LIBFLAGS} -lnetcdf ${OBJstrahl} 
# NetCDF and STRAHL must be compiled with the -fPIC option
#
# for netcdf
# ./configure --prefix=/home/dwagner/strahl/source/libnetcdf/ FCFLAGS=-fPIC CFLAGS=-fPIC FC=ifort
