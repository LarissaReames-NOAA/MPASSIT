
if [[ $- != *i* ]]; then
    source /usr/share/Modules/init/bash
fi

module purge

module load hpcx-mt-ompi-gcc
module load gcc/13.2.0

export GNU_LIBPATH=/scratch/ywang/tools/gnu

export NETCDF=${GNU_LIBPATH}
export PNETCDF=${GNU_LIBPATH}

export LD_LIBRARY_PATH=/scratch/software/gcc/gcc-13.2.0/lib64:${LD_LIBRARY_PATH}:/scratch/ywang/MPAS/tools/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${GNU_LIBPATH}/lib:${GNU_LIBPATH}/lib64

export PATH=.:$PATH:${GNU_LIBPATH}/bin:/scratch/ywang/MPAS/tools/bin
