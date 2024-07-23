help([[
This module loads libraries for building the UFS SRW App on
the NSSL machine Vecna using Intel
]])

whatis([===[Loads libraries needed for building the UFS SRW App on Vecna ]===])

load(pathJoin("compiler", "latest"))
load(pathJoin("mpi", "latest"))

prepend_path("MODULEPATH","/scratch/ywang/tools/intel/hpc-stack/modulefiles/stack")
load(pathJoin("hpc", os.getenv("hpc_ver") or "1.2.0"))
load(pathJoin("hpc-intel", os.getenv("hpc_intel_ver") or "2021.10.0"))
load("hpc-openmpi/hpc-x-intel-classic")

load("esmf")
load("netcdf")

setenv("MPI_HOME","/scratch/software/hpc-x-intel-classic/hpcx-ompi")

setenv("CMAKE_C_COMPILER","mpicc")
setenv("CMAKE_CXX_COMPILER","mpicxx")
setenv("CMAKE_Fortran_COMPILER","mpifort")
setenv("CMAKE_Platform","vecna.intel")

