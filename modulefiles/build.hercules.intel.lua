help([[
This module loads libraries for building the UFS SRW App on
the MSU machine Hercules using intel-oneapi-compilers/2022.2.1
]])

whatis([===[Loads libraries needed for building the UFS SRW App on Orion ]===])

load("contrib")
load("noaatools")

load(pathJoin("cmake", os.getenv("cmake_ver") or "3.26.3"))

--prepend_path("MODULEPATH","/work/noaa/epic/role-epic/contrib/hercules/hpc-stack/intel-2022.2.1/modulefiles/stack")
prepend_path("MODULEPATH","/work2/noaa/wof/ywang/tools/hpc-stack/modulefiles/stack")
load(pathJoin("hpc", os.getenv("hpc_ver") or "1.2.0"))
load(pathJoin("hpc-intel-oneapi-compilers", os.getenv("hpc_intel_ver") or "2022.2.1"))
load(pathJoin("hpc-intel-oneapi-mpi", os.getenv("hpc_mpi_ver") or "2021.7.1"))

--load("srw_common")
load_any("jasper/2.0.25","jasper/2.0.32")
load_any("zlib/1.2.11","zlib/1.2.13")
load("libpng/1.6.37")

load_any("netcdf/4.9.2", "netcdf-c/4.9.2")
load_any("netcdf/4.9.2", "netcdf-fortran/4.6.0")
--load_any("pio/2.5.10", "parallelio/2.5.9")
load("esmf/8.4.2")
load("fms/2023.01")

load("bacio/2.4.1")
load("g2/3.4.5")
load("crtm/2.4.0")
load("g2tmpl/1.10.2")
load("ip/3.3.3")
load("sp/2.3.3")
load("w3emc/2.9.2")

load_any("gftl-shared/v1.5.0", "gftl-shared/1.5.0")
load_any("yafyaml/v0.5.1", "yafyaml/0.5.1")
load("mapl/2.35.2-esmf-8.4.2")

load("nemsio/2.5.4")
load("sfcio/1.4.1")
load("sigio/2.3.2")
load("w3nco/2.4.1")
load_any("wrf_io/1.2.0","wrf-io/1.2.0")

load("pnetcdf")
load("wgrib2/2.0.8")
--

load(pathJoin("nccmp", os.getenv("nccmp_ver") or "1.9.1.0"))
load(pathJoin("nco", os.getenv("nco_ver") or "5.0.6"))

setenv("CFLAGS","-diag-disable=10441")
setenv("FFLAGS","-diag-disable=10441")

setenv("CMAKE_C_COMPILER","mpiicc")
setenv("CMAKE_CXX_COMPILER","mpiicpc")
setenv("CMAKE_Fortran_COMPILER","mpiifort")
setenv("CMAKE_Platform","hercules.intel")

