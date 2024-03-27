help([[
Load environment to compile MPASSIT for MPAS-WoFS on Jet
]])

whatis("Description: MPASSIT build environment")

cmake_ver=os.getenv("cmake_ver") or "3.28.1"
load(pathJoin("cmake", cmake_ver))

load(pathJoin("gnu","9.2.0"))

load(pathJoin("intel","2023.2.0"))
load(pathJoin("impi","2023.2.0"))

load(pathJoin("pnetcdf","1.12.3"))
load("szip")
load(pathJoin("hdf5parallel","1.10.5"))
load(pathJoin("netcdf-hdf5parallel","4.7.0"))

setenv("PNETCDF","/apps/pnetcdf/1.12.3/intel_2023.2.0-impi")

prepend_path("LD_LIBRARY_PATH", "/lfs4/NAGAPE/wof/miniconda3_RL/lib")   -- contains the grib2 libraries that are needed for the WPS compile
prepend_path("CPATH",           "/usr/include/tirpc")                   -- only be important to the WRF build

-- # the Jasper environment settings for WPS
setenv("JASPERLIB", "/lfs4/NAGAPE/wof/miniconda3_RL/lib")
setenv("JASPERINC", "/lfs4/NAGAPE/wof/miniconda3_RL/lib/include/jasper")

-- # ESMF V8.6 for MPASSIT
setenv("ESMFMKFILE", "/lfs4/NAGAPE/hpc-wof1/ywang/tools/esmf-8.6.0/lib/libO/Linux.intel.64.intelmpi.default/esmf.mk")

setenv("CMAKE_C_COMPILER", "mpiicc")
setenv("CMAKE_CXX_COMPILER", "mpiicpc")
setenv("CMAKE_Fortran_COMPILER", "mpiifort")
