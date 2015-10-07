
# --Choose your system

SYSTYPE = "Basic"

HDF5PATH = /usr/local/hdf5-1.8.7-linux-x86_64-static/
CFITSIOPATH = /data2/aaorsi/cfits
ifeq ($(SYSTYPE), "Basic")
	CC        = $(HDF5PATH)/bin/h5cc
	OPTIMIZE  = -g
	CFLAGS    = $(OPTIMIZE) -I$(HDF5PATH)/include -I$(CFITSIOPATH)/include # -I/usr/local/include # -DDEBUG
	LIBS      = -lm -L$(HDF5PATH)/lib -L$(CFITSIOPATH)/lib -lnsl -lcfitsio   # -lgsl -lgslcblas 
endif

#---------------------------------------------------------

EXEC   = get_mags 

OBJS   = main.o filters.o read_saghdf5.o igm_attenuation.o write_hdf5.o write_fits.o

INCL   = proto.h allvars.h $(HDF5PATH)/include/hdf5.h $(HDF5PATH)/include/hdf5_hl.h $(CFITSIOPATH)/include/fitsio.h Makefile

#---------------------------------------------------------

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o  $(EXEC)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC) core

#---------------------------------------------------------

