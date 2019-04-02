MAKE = make

F90 = ftn
CXX = CC
CC = cc

CPP_OPTS = -D_MPI_ 

F90_OPTS = -O3 -cpp $(CPP_OPTS) -qopenmp -Wl -z muldefs -funroll-loops -check bound

F90_LINK_OPTS = $(F90_OPTS)

LAPACK_LIB = -mkl=cluster


#### ========== compile with libXC support ==========
###CPP_OPTS := $(CPP_OPTS) -D_LIBXC_
###XC_LIB = ../libxc.a

#### ========== compile with NFFT support ==========
###CPP_OPTS := $(CPP_OPTS) -D_NFFT_
###NFFT_INC = -I$(HOME)/local/include
###NFFT_LIB = $(HOME)/local/lib/libnfft3.a $(HOME)/local/lib/libfftw3.a

# === compile with HDF5 library support ===
CPP_OPTS := $(CPP_OPTS) -D_HDF5_
###HDF5_LIB = /global/homes/z/zhang714/install_excitingplus_cori_2017Jan04/libs/hdf5-1.8.12/build/lib/libhdf5hl_fortran.a /global/homes/z/zhang714/install_excitingplus_cori_2017Jan04/libs/hdf5-1.8.12/build/lib/libhdf5_hl.a /global/homes/z/zhang714/install_excitingplus_cori_2017Jan04/libs/hdf5-1.8.12/build/lib/libhdf5_fortran.a /global/homes/z/zhang714/install_excitingplus_cori_2017Jan04/libs/hdf5-1.8.12/build/lib/libhdf5.a -lz -lm

# ========== compile with SIRIUS library support ==========
CPP_OPTS := $(CPP_OPTS) -D_SIRIUS_
SIRIUS_LIB = $(HOME)/install_sirius600_cori/SIRIUS-develop/src/libsirius.a $(HOME)/install_sirius600_cori/SIRIUS-develop/src/libsirius_f.a

# ========== collect all libraries under one name ==========
LIBS = $(LAPACK_LIB) $(HDF5_LIB) $(SIRIUS_LIB)


