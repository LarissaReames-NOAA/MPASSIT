
if [[ $- != *i* ]]; then
    source /usr/share/Modules/init/bash
fi

module purge
module load hpcx-mt-ompi-gcc
module load gcc/13.2.0

# ZLIB:  ./configure --prefix=/scratch/ywang/tools/gnu
# HDF5:  CC=mpicc ./configure --with-zlib=/scratch/ywang/tools/gnu --prefix=/scratch/ywang/tools/gnu --enable-parallel
# netcdf-C: CC=mpicc CPPFLAGS='-I/scratch/ywang/tools/gnu/include' LDFLAGS='-L/scratch/ywang/tools/gnu/lib'
#           ./configure --prefix=/scratch/ywang/tools/gnu --disable-libxml2 --disable-shared --enable-parallel-tests
# netcdf-fortran:
#       export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/scratch/ywang/tools/gnu/lib
#       CC=mpicc CPPFLAGS='-I/scratch/ywang/tools/gnu/include' LDFLAGS='-L/scratch/ywang/tools/gnu/lib -lnetcdf -lhdf5_hl -lhdf5 -lm -lz -ldl -lbz2 -lzstd -lcurl -lstdc++'
#       LIBS='-L/scratch/ywang/tools/gnu/lib -lnetcdf -lhdf5_hl -lhdf5 -lm -lz -ldl -lbz2 -lzstd -lcurl -lstdc++'
#       ./configure --prefix=/scratch/ywang/tools/gnu --disable-shared --enable-parallel-tests


export CC=/scratch/software/gcc/gcc-13.2.0/bin/gcc
export FC=/scratch/software/gcc/gcc-13.2.0/bin/gfortran
export F77=/scratch/software/gcc/gcc-13.2.0/bin/gfortran
export F90=/scratch/software/gcc/gcc-13.2.0/bin/gfortran
export CXX=g++
export MPICC=mpicc
export MPIF77=mpif90
export MPIF90=mpif90

export ZLIB=/scratch/ywang/tools/gnu
export H5DIR=/scratch/ywang/tools/gnu
export NETCDF=/scratch/ywang/tools/gnu
export PNETCDF=/scratch/ywang/tools/pnetcdf-1.12.3
export ESMFMKFILE=/scratch/ywang/tools/gnu/esmf-8.6.0/lib/libO/Linux.gfortran.64.openmpi.default/esmf.mk

#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/scratch/software/intel/grib2/lib:/usr/lib64
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/scratch/ywang/MPAS/tools/lib
export LD_LIBRARY_PATH=${NETCDF}/lib:/lib64:${LD_LIBRARY_PATH}

export PATH=./:/bin:${NETCDF}/bin:/scratch/ywang/MPAS/tools/bin:$PATH

gfortran --version
if [[ $- == *i* ]]; then
    appendpath /scratch/ywang/MPAS/tools/bin /scratch/ywang/MPAS/mpas_runscripts/scripts
    #appendlibrary /scratch/ywang/MPAS/tools/lib

    function mpastime {
       #rand=$RANDOM
       tmpfile=$(mktemp -t timestep_XXX)
       grep -E 'Timing for integration step:|Begin timestep' $1 > ${tmpfile}
       paste -d " " - - < ${tmpfile}
       rm -rf ${tmpfile}
    }
fi
