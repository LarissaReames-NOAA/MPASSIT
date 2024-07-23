help([[
This module loads libraries for building MPASSIT on NOAA's Jet
]])

whatis([===[Loads libraries needed for building MPASSIT on Jet ]===])

prepend_path("MODULEPATH","/lfs4/NAGAPE/hpc-wof1/ywang/MPAS/mpas_scripts/modules/")

-- Set up path for CMake
prepend_path("PATH","/apps/cmake/3.28.1/bin")

--Load intel compiler (gnu needed for this version of intel) and intel MPI
load("gnu/13.2.0")
load("intel/2023.2.0")
load("impi/2023.2.0")

--Set up ld library paths that are required for netcdf (this one needs hdf5 which needs szip) 
append_path("LD_LIBRARY_PATH","/apps/szip/2.1/lib")
prepend_path("LD_LIBRARY_PATH","/apps/hdf5/1.10.5/intel_2023.2.0-impi/lib")

--Export environmental variable NETCDF pointing to NetCDF library
setenv("NETCDF","/apps/netcdf/4.7.0/intel_2023.2.0-impi")

--Export environmental variable pointing to ESMF file esmf.mk
setenv("ESMFMKFILE", "/lfs4/NAGAPE/hpc-wof1/ywang/tools/esmf-8.6.0/lib/libO/Linux.intel.64.intelmpi.default/esmf.mk")
