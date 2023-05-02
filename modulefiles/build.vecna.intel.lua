help([[
This module loads libraries for building the UFS SRW App on
the NSSL machine Vecna using Intel
]])

whatis([===[Loads libraries needed for building the UFS SRW App on Vecna ]===])

load(pathJoin("compiler", "latest"))
load(pathJoin("mpi", "latest"))

prepend_path("MODULEPATH","/scratch/ywang/tools/hpc-stack/modulefiles/stack")
load(pathJoin("hpc", os.getenv("hpc_ver") or "1.2.0"))
load(pathJoin("hpc-intel", os.getenv("hpc_intel_ver") or "2021.8.0"))
load(pathJoin("hpc-impi", os.getenv("hpc_impi_ver") or "2021.8.0"))

load("esmf")
load("netcdf")

setenv("CMAKE_C_COMPILER","mpiicc")
setenv("CMAKE_CXX_COMPILER","mpiicpc")
setenv("CMAKE_Fortran_COMPILER","mpiifort")
setenv("CMAKE_Platform","vecna.intel")

