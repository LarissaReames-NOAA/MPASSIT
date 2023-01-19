 module write_data

 use program_setup, only   : output_file
 use datetime_module, only: datetime, timedelta, clock
 private

 public :: write_to_file

 contains

!> Write data on the target grid
!!
!!
!! @param[in] localpet  ESMF local persistent execution thread
!! @author Larissa Reames CIWRO/NOAA/NSSL

 subroutine write_to_file(localpet)

 use esmf
 use netcdf
 use mpi

 use program_setup, only           : interp_diag, interp_hist, &
                                     wrf_mod_vars

 use model_grid, only              : i_target, j_target, &
                                     ip1_target, jp1_target, &
                                     nz_input, nzp1_input, &
                                     nsoil_input, &
                                     dx, start_time, &
                                     config_dt, &
                                     strlen, valid_time, &
                                     lsm_scheme, mp_scheme, &
                                     conv_scheme, &
                                     cen_lat, cen_lon, truelat1, &
                                     truelat2, moad_cen_lat, &
                                     stand_lon, pole_lat, pole_lon, &
                                     map_proj, map_proj_char, &
                                     longitude_target_grid, &
                                     latitude_target_grid, &
                                     latitude_u,longitude_u, &
                                     latitude_v, longitude_v, &
                                     zs_target_grid, &
                                     hgt_target_grid, &
                                     target_diag_bundle, &
                                     n_diag_fields, &
                                     target_hist_bundle_2d_patch, &
                                     target_hist_bundle_2d_cons, &
                                     target_hist_bundle_2d_nstd, &
                                     target_hist_bundle_3d_nz, &  
                                     target_hist_bundle_3d_nzp1, &
                                     target_hist_bundle_soil, &
                                     n_hist_fields_2d_patch, &
                                     n_hist_fields_2d_cons, &
                                     n_hist_fields_2d_nstd, &
                                     n_hist_fields_3d_nz, &
                                     n_hist_fields_3d_nzp1, &
                                     n_hist_fields_soil, &
                                     target_diag_units, &
                                     target_hist_units_2d_cons, &
                                     target_hist_units_2d_nstd, &
                                     target_hist_units_2d_patch, &
                                     target_hist_units_3d_nzp1, &
                                     target_hist_units_3d_nz, &
                                     target_hist_units_soil, &
                                     target_diag_longname, &
                                     target_hist_longname_2d_cons, &
                                     target_hist_longname_2d_nstd, &
                                     target_hist_longname_2d_patch, &
                                     target_hist_longname_3d_nzp1, &
                                     target_hist_longname_3d_nz, &
                                     target_hist_longname_soil

 implicit none

 integer, intent(in)              :: localpet

 character(len=128)               :: outfile
 character(len=50)                :: varname
 character(len=20)                 :: tempstr(1,19)

 integer, parameter               :: Datestrlen=19
 integer                          :: error, ncid, n, rc, i, j, k, m
 integer                          :: header_buffer_val = 16384
 integer                          :: dim_time, dim_lon, dim_lat, dim_z, dim_zp1, dim_soil
 integer                          :: dim_lonp, dim_latp, dim_str, dim_lon_stag, dim_lat_stag
 integer                          :: id_lat, id_lon, id_z, id_zs, id_times, id_xtime, id_itime
 integer                          :: id_latu, id_latv, id_lonu, id_lonv, id_ph, id_mu, id_hgt, id_ptop
 integer                          :: n2d, n3d, ndims
 integer                          :: sy,sm,sd,sh,smi,ss,  vy,vm,vd,vh,vmi,vs
 integer, allocatable             :: id_vars2(:), id_vars3_nz(:), id_vars3_nzp1(:), &
                                     id_vars_soil(:)

 real(esmf_kind_r8), allocatable  :: dum2d(:,:), dum2dt(:,:,:), &
                                     dum2dtu(:,:,:), dum2dtv(:,:,:), &
                                     dum3d(:,:,:), dum3dt(:,:,:,:), &
                                     dum3dp1(:,:,:), dum3dp1t(:,:,:,:), &
                                     dumsoil(:,:,:), dumsoilt(:,:,:,:), &
                                     dumsmall(:,:), dum3dtmp(:,:,:), dum1d(:)
 
 type(esmf_field), allocatable    :: fields(:), field_write_2d(:), field_extra3(:)
 type (timedelta)                 :: xtime_dt

 n2d = n_diag_fields + n_hist_fields_2d_patch + n_hist_fields_2d_nstd + n_hist_fields_2d_cons
 allocate(field_write_2d(n2d),id_vars2(n2d))
 allocate(field_extra3(n2d))  !allocate large incase all diag fields are 3d
 allocate(id_vars3_nz(n_hist_fields_3d_nz+1))
 allocate(id_vars3_nzp1(n_hist_fields_3d_nzp1))
 allocate(id_vars_soil(n_hist_fields_soil))
 
 if (localpet ==0) then
   allocate(dumsmall(nsoil_input,1))
   allocate(dum2d(i_target,j_target))
   allocate(dum2dt(i_target,j_target,1))
   allocate(dum2dtu(ip1_target,j_target,1))
   allocate(dum2dtv(i_target,jp1_target,1))
   allocate(dum3d(i_target,j_target,nz_input))
   allocate(dum3dt(i_target,j_target,nz_input,1))
   allocate(dum3dp1(i_target,j_target,nzp1_input))
   allocate(dum3dp1t(i_target,j_target,nzp1_input,1))
   allocate(dumsoil(i_target,j_target,nsoil_input))
   allocate(dumsoilt(i_target,j_target,nsoil_input,1))
   allocate(dum1d(1))
 else
   allocate(dumsmall(0,0))
   allocate(dum2d(0,0))
   allocate(dum2dt(0,0,0))
   allocate(dum3d(0,0,0))
   allocate(dum3dt(0,0,0,0))
   allocate(dum3dp1(0,0,0))
   allocate(dum3dp1t(0,0,0,0))
   allocate(dumsoil(0,0,0))
   allocate(dumsoilt(0,0,0,0))
   allocate(dum1d(0))
 endif
 
if (localpet == 0) then

!--- open the file
   error = nf90_create(output_file, NF90_NETCDF4, ncid)
   call netcdf_err(error, 'CREATING FILE '//trim(output_file) )

!--- define dimension
   error = nf90_def_dim(ncid, 'Time', NF90_UNLIMITED , dim_time)
   call netcdf_err(error, 'DEFINING Time DIMENSION' )
   error = nf90_def_dim(ncid, 'west_east', i_target, dim_lon)
   call netcdf_err(error, 'DEFINING LON DIMENSION' )
   error = nf90_def_dim(ncid, 'west_east_stag', i_target+1, dim_lon_stag)
   call netcdf_err(error, 'DEFINING STAGGERED LON DIMENSION' )
   error = nf90_def_dim(ncid, 'south_north', j_target, dim_lat)
   call netcdf_err(error, 'DEFINING LAT DIMENSION' )
   error = nf90_def_dim(ncid, 'south_north_stag', j_target+1, dim_lat_stag)
   call netcdf_err(error, 'DEFINING STAGGERED LAT DIMENSION' )
   error = nf90_def_dim(ncid, 'bottom_top', nz_input, dim_z)
   call netcdf_err(error, 'DEFINING VERTICAL DIMENSION' )
   error = nf90_def_dim(ncid, 'bottom_top_stag', nzp1_input, dim_zp1)
   call netcdf_err(error, 'DEFINING VERTICALP1 DIMENSION' )
   error = nf90_def_dim(ncid, 'soil_layers_stag', nsoil_input, dim_soil)
   call netcdf_err(error, 'DEFINING VERTICALP1 DIMENSION' )
   error = nf90_def_dim(ncid, 'StrLen', Datestrlen, dim_str)
   call netcdf_err(error, 'DEFINING STRLEN DIMENSION' )
   
 !--- define global attributes
   error = nf90_put_att(ncid, NF90_GLOBAL, 'WEST-EAST_GRID_DIMENSION', i_target+1)
   call netcdf_err(error, 'DEFINING WEST-EAST GRID DIMENSION GLOBAL ATTRIBUTE')
   
   error = nf90_put_att(ncid, NF90_GLOBAL, 'SOUTH-NORTH_GRID_DIMENSION', j_target+1)
   call netcdf_err(error, 'DEFINING NORTH-SOUTH GRID DIMENSION GLOBAL ATTRIBUTE')
   
   error = nf90_put_att(ncid, NF90_GLOBAL, 'BOTTOM-TOP_GRID_DIMENSION', nz_input+1)
   call netcdf_err(error, 'DEFINING BOTTOM-TOP GRID DIMENSION GLOBAL ATTRIBUTE')
   
   error = nf90_put_att(ncid, NF90_GLOBAL, 'SIMULATION_START_DATE', start_time)
   call netcdf_err(error, 'DEFINING SUMLATION START DATE GLOBAL ATTRIBUTE')
   
   error = nf90_put_att(ncid, NF90_GLOBAL, 'START_DATE', start_time)
   call netcdf_err(error, 'DEFINING START DATE GLOBAL ATTRIBUTE')
   
   error = nf90_put_att(ncid, NF90_GLOBAL, 'DX', dx)
   call netcdf_err(error, 'DEFINING DX GLOBAL ATTRIBUTE')
   
   error = nf90_put_att(ncid, NF90_GLOBAL, 'SF_SURFACE_PHYSICS', lsm_scheme)
   call netcdf_err(error, 'DEFINING SF SURFACE PHYSICS GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'MP_PHYSICS', mp_scheme)
   call netcdf_err(error, 'DEFINING MP PHYSICS GLOBAL ATTRIBUTE')
   
   error = nf90_put_att(ncid, NF90_GLOBAL, 'CU_PHYSICS', conv_scheme)
   call netcdf_err(error, 'DEFINING CU PHYSICS GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'CEN_LAT', cen_lat)
   call netcdf_err(error, 'DEFINING CEN_LAT GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'CEN_LON', cen_lon)
   call netcdf_err(error, 'DEFINING CEN_LON GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'TRUELAT1', truelat1)
   call netcdf_err(error, 'DEFINING TRUELAT1 GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'TRUELAT2', truelat2)
   call netcdf_err(error, 'DEFINING TRUELAT2 GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'MOAD_CEN_LAT', cen_lat)
   call netcdf_err(error, 'DEFINING MOAD_CEN_LAT GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'STAND_LON', stand_lon)
   call netcdf_err(error, 'DEFINING STAND_LON GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'POLE_LAT', pole_lat)
   call netcdf_err(error, 'DEFINING POLELAT GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'POLE_LON', pole_lon)
   call netcdf_err(error, 'DEFINING POLE_LON GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'POL_ELAT', pole_lat)
   call netcdf_err(error, 'DEFINING POLE_LAT GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'MAP_PROJ', map_proj)
   call netcdf_err(error, 'DEFINING MAP_PROJ GLOBAL ATTRIBUTE')

   error = nf90_put_att(ncid, NF90_GLOBAL, 'MAP_PROJ_CHAR', map_proj_char)
   call netcdf_err(error, 'DEFINING MAP_PROJ_CHAR GLOBAL ATTRIBUTE')

!--- define fields

   error = nf90_def_var(ncid, 'XLONG', NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_lon)
   call netcdf_err(error, 'DEFINING GEOLON FIELD' )
   error = nf90_put_att(ncid, id_lon, "description", "LONGITUDE, WEST IS NEGATIVE")
   call netcdf_err(error, 'DEFINING GEOLON NAME' )
   error = nf90_put_att(ncid, id_lon, "units", "degree_east")
   call netcdf_err(error, 'DEFINING GEOLON UNITS' )
   error = nf90_put_att(ncid, id_lon, "MemoryOrder", "XY " )
    call netcdf_err(error, 'DEFINING MEMORYORDER' )
    error = nf90_put_att(ncid, id_lon, "coordinates", "XLONG XLAT")
    call netcdf_err(error, 'DEFINING COORD' )
   error = nf90_put_att(ncid, id_lon, "stagger", "")
   call netcdf_err(error, 'DEFINING STAGGER' )
   error = nf90_put_att(ncid, id_lon, "FieldType", 104)
   call netcdf_err(error, 'DEFINING FieldType' )
   
   error = nf90_def_var(ncid, 'XLONG_U', NF90_FLOAT, (/dim_lon_stag, dim_lat, dim_time/), id_lonu)
   call netcdf_err(error, 'DEFINING GEOLON FIELD' )
   error = nf90_put_att(ncid, id_lonu, "description", "LONGITUDE, WEST IS NEGATIVE")
   call netcdf_err(error, 'DEFINING GEOLON NAME' )
   error = nf90_put_att(ncid, id_lonu, "units", "degree_east")
   call netcdf_err(error, 'DEFINING GEOLON UNITS' )
   error = nf90_put_att(ncid, id_lonu, "MemoryOrder", "XY " )
    call netcdf_err(error, 'DEFINING MEMORYORDER' )
    error = nf90_put_att(ncid, id_lonu, "coordinates", "XLONG_U XLAT_U")
    call netcdf_err(error, 'DEFINING COORD' )
   error = nf90_put_att(ncid, id_lonu, "stagger", "X")
   call netcdf_err(error, 'DEFINING STAGGER' )
   error = nf90_put_att(ncid, id_lonu, "FieldType", 104)
   call netcdf_err(error, 'DEFINING FieldType' )
   
   error = nf90_def_var(ncid, 'XLONG_V', NF90_FLOAT, (/dim_lon, dim_lat_stag, dim_time/), id_lonv)
   call netcdf_err(error, 'DEFINING GEOLON FIELD' )
   error = nf90_put_att(ncid, id_lonv, "description", "LONGITUDE, WEST IS NEGATIVE")
   call netcdf_err(error, 'DEFINING GEOLON NAME' )
   error = nf90_put_att(ncid, id_lonv, "units", "degree_east")
   call netcdf_err(error, 'DEFINING GEOLON UNITS' )
   error = nf90_put_att(ncid, id_lonv, "MemoryOrder", "XY " )
    call netcdf_err(error, 'DEFINING MEMORYORDER' )
    error = nf90_put_att(ncid, id_lonv, "coordinates", "XLONG_V XLAT_V")
    call netcdf_err(error, 'DEFINING COORD' )
   error = nf90_put_att(ncid, id_lonv, "stagger", "Y")
   call netcdf_err(error, 'DEFINING STAGGER' )
   error = nf90_put_att(ncid, id_lonv, "FieldType", 104)
   call netcdf_err(error, 'DEFINING FieldType' )

   error = nf90_def_var(ncid, 'XLAT', NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_lat)
   call netcdf_err(error, 'DEFINING GEOLAT FIELD' )
   error = nf90_put_att(ncid, id_lat, "description", "LATITUDE, SOUTH IS NEGATIVE")
   call netcdf_err(error, 'DEFINING GEOLAT NAME' )
   error = nf90_put_att(ncid, id_lat, "units", "degree_north")
   call netcdf_err(error, 'DEFINING GEOLAT UNITS' )
   error = nf90_put_att(ncid, id_lat, "MemoryOrder", "XY " )
    call netcdf_err(error, 'DEFINING MEMORYORDER' )
    error = nf90_put_att(ncid, id_lat, "coordinates", "XLONG XLAT")
    call netcdf_err(error, 'DEFINING COORD' )
   error = nf90_put_att(ncid, id_lat, "stagger", "")
   call netcdf_err(error, 'DEFINING STAGGER' )
   error = nf90_put_att(ncid, id_lat, "FieldType", 104)
   call netcdf_err(error, 'DEFINING FieldType' )
   
   error = nf90_def_var(ncid, 'XLAT_U', NF90_FLOAT, (/dim_lon_stag,dim_lat, dim_time/), id_latu)
   call netcdf_err(error, 'DEFINING GEOLAT FIELD' )
   error = nf90_put_att(ncid, id_latu, "description", "LATITUDE, SOUTH IS NEGATIVE")
   call netcdf_err(error, 'DEFINING GEOLAT NAME' )
   error = nf90_put_att(ncid, id_latu, "units", "degree_north")
   call netcdf_err(error, 'DEFINING GEOLAT UNITS' )
   error = nf90_put_att(ncid, id_latu, "MemoryOrder", "XY " )
    call netcdf_err(error, 'DEFINING MEMORYORDER' )
    error = nf90_put_att(ncid, id_latu, "coordinates", "XLONG_U XLAT_U")
    call netcdf_err(error, 'DEFINING COORD' )
   error = nf90_put_att(ncid, id_latu, "stagger", "X")
   call netcdf_err(error, 'DEFINING STAGGER' )
   error = nf90_put_att(ncid, id_latu, "FieldType", 104)
   call netcdf_err(error, 'DEFINING FieldType' )
   
   error = nf90_def_var(ncid, 'XLAT_V', NF90_FLOAT, (/dim_lon,dim_lat_stag, dim_time/), id_latv)
   call netcdf_err(error, 'DEFINING GEOLAT FIELD' )
   error = nf90_put_att(ncid, id_latv, "description", "LATITUDE, SOUTH IS NEGATIVE")
   call netcdf_err(error, 'DEFINING GEOLAT NAME' )
   error = nf90_put_att(ncid, id_latv, "units", "degree_north")
   call netcdf_err(error, 'DEFINING GEOLAT UNITS' )
   error = nf90_put_att(ncid, id_latv, "MemoryOrder", "XY " )
    call netcdf_err(error, 'DEFINING MEMORYORDER' )
    error = nf90_put_att(ncid, id_latv, "coordinates", "XLONG_V XLAT_V")
    call netcdf_err(error, 'DEFINING COORD' )
   error = nf90_put_att(ncid, id_latv, "stagger", "Y")
   call netcdf_err(error, 'DEFINING STAGGER' )
   error = nf90_put_att(ncid, id_latv, "FieldType", 104)
   call netcdf_err(error, 'DEFINING FieldType' )
   
   error = nf90_def_var(ncid, 'Z_C', NF90_FLOAT, (/dim_lon,dim_lat,dim_zp1, dim_time/), id_z)
   call netcdf_err(error, 'DEFINING Z_C FIELD' )
   error = nf90_put_att(ncid, id_z, "description", "Layer center height above mean sea level")
   call netcdf_err(error, 'DEFINING Z_C NAME' )
   error = nf90_put_att(ncid, id_z, "units", "m AMSL")
   call netcdf_err(error, 'DEFINING Z_C UNITS' )
   error = nf90_put_att(ncid, id_z, "MemoryOrder", "XYZ ")
    call netcdf_err(error, 'DEFINING MEMORYORDER' )
    error = nf90_put_att(ncid, id_z, "coordinates", "XLAT XLONG Z_C")
    call netcdf_err(error, 'DEFINING COORD' )
   error = nf90_put_att(ncid, id_z, "stagger", "")
   call netcdf_err(error, 'DEFINING STAGGER' )
   error = nf90_put_att(ncid, id_z, "FieldType", 104)
   call netcdf_err(error, 'DEFINING FieldType' )
   
   error = nf90_def_var(ncid, 'ZS', NF90_FLOAT, (/dim_soil, dim_time/), id_zs)
   call netcdf_err(error, 'DEFINING ZS FIELD' )
   error = nf90_put_att(ncid, id_zs, "description", "DEPTHS OF CENTERS OF SOIL LAYERS")
   call netcdf_err(error, 'DEFINING ZS NAME' )
   error = nf90_put_att(ncid, id_zs, "units", "m")
   call netcdf_err(error, 'DEFINING ZS UNITS' )
   error = nf90_put_att(ncid, id_zs, "MemoryOrder", "X")
    call netcdf_err(error, 'DEFINING MEMORYORDER' )
    error = nf90_put_att(ncid, id_zs, "coordinates", "ZS XTIME")
    call netcdf_err(error, 'DEFINING COORD' )
   error = nf90_put_att(ncid, id_zs, "stagger", "")
   call netcdf_err(error, 'DEFINING STAGGER' )
   error = nf90_put_att(ncid, id_zs, "FieldType", 104)
   call netcdf_err(error, 'DEFINING FieldType' )

   error = nf90_def_var(ncid, 'HGT', NF90_FLOAT, (/dim_lon,dim_lat,dim_time/), id_hgt)
   call netcdf_err(error, 'DEFINING HGT FIELD' )
   error = nf90_put_att(ncid, id_hgt, "description", "TERRAIN HEIGHT ")
   call netcdf_err(error, 'DEFINING HGT NAME' )
   error = nf90_put_att(ncid, id_hgt, "units", "m AMSL")
   call netcdf_err(error, 'DEFINING HGT UNITS' )
   error = nf90_put_att(ncid, id_hgt, "MemoryOrder", "XY ")
    call netcdf_err(error, 'DEFINING MEMORYORDER' )
    error = nf90_put_att(ncid, id_hgt, "coordinates", "XLAT XLONG ")
    call netcdf_err(error, 'DEFINING COORD' )
   error = nf90_put_att(ncid, id_hgt, "stagger", "")
   call netcdf_err(error, 'DEFINING STAGGER' )
   error = nf90_put_att(ncid, id_hgt, "FieldType", 104)
   call netcdf_err(error, 'DEFINING FieldType' )

   error = nf90_def_var(ncid, 'Times', NF90_CHAR, (/dim_str, dim_time/), id_times)
   call netcdf_err(error, 'DEFINING Times FIELD' )
   error = nf90_put_att(ncid, id_times, "description", "Times")
   call netcdf_err(error, 'DEFINING Times NAME' )
   error = nf90_put_att(ncid, id_times, "units", "m")
   call netcdf_err(error, 'DEFINING Times UNITS' )
   error = nf90_put_att(ncid, id_times, "coordinates", "Time")
   call netcdf_err(error, 'DEFINING Times COORD' )
   error = nf90_put_att(ncid, id_times, "stagger", "")
   call netcdf_err(error, 'DEFINING STAGGER' )
   error = nf90_put_att(ncid, id_times, "FieldType", 104)
   call netcdf_err(error, 'DEFINING FieldType' )

   error = nf90_def_var(ncid, 'ITIMESTEP', NF90_INT, (/dim_time/),id_itime)
   call netcdf_err(error, 'DEFINING ITIMESTEP FIELD' )
   error = nf90_put_att(ncid, id_itime, "description", "")
   call netcdf_err(error, 'DEFINING ITIMESTEP NAME' )
   error = nf90_put_att(ncid, id_itime, "units", "")
   call netcdf_err(error, 'DEFINING ITIMESTEP UNITS' )
   error = nf90_put_att(ncid, id_itime, "stagger", "")
   call netcdf_err(error, 'DEFINING STAGGER' )
   error = nf90_put_att(ncid, id_itime, "FieldType", 106)
   call netcdf_err(error, 'DEFINING FieldType' )
    error = nf90_put_att(ncid, id_itime, "MemoryOrder", "O ")
   call netcdf_err(error, 'DEFINING MemoryOrder' )

   error = nf90_def_var(ncid, 'XTIME', NF90_FLOAT, (/dim_time/),id_xtime)
   call netcdf_err(error, 'DEFINING XTIME FIELD' )
   error = nf90_put_att(ncid, id_xtime, "description", "minutes since "//start_time)
   call netcdf_err(error, 'DEFINING XTIME NAME' )
   error = nf90_put_att(ncid, id_xtime, "units", "minutes since "//start_time)
   call netcdf_err(error, 'DEFINING XTIME UNITS' )
   error = nf90_put_att(ncid, id_xtime, "stagger", "")
   call netcdf_err(error, 'DEFINING STAGGER' )
   error = nf90_put_att(ncid, id_xtime, "FieldType", 104)
   call netcdf_err(error, 'DEFINING FieldType' )
    error = nf90_put_att(ncid, id_xtime, "MemoryOrder", "O ")
   call netcdf_err(error, 'DEFINING MemoryOrder' )

 endif
   
   k = 0
   m = 0
   if (interp_diag .and. n_diag_fields>0) then
        allocate(fields(n_diag_fields))
        call ESMF_FieldBundleGet(target_diag_bundle, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=error) 
       if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
              call error_handler("IN FieldBundleGet", error)
        do i = 1, n_diag_fields  
            call ESMF_FieldGet(fields(i),name=varname,rc=error)
            if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
                  call error_handler("IN FieldGet", error)
             call ESMF_FieldGet(fields(i), dimCount=ndims, rc=rc)
                if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
                  call error_handler("IN FieldGet", error) 
             
             if (ndims==2) then
               k = k+1
               field_write_2d(k) = fields(i) 
               if (localpet==0) then
                    print*,"- DEFINE 2d diag ON FILE TARGET GRID ", varname
                    error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_vars2(k))
                    call netcdf_err(error, 'DEFINING VAR' )
                    error = nf90_put_att(ncid, id_vars2(k), "MemoryOrder", "XY " )
                    call netcdf_err(error, 'DEFINING MEMORYORDER' )
                    error = nf90_put_att(ncid, id_vars2(k), "coordinates", "XLONG XLAT XTIME")
                    call netcdf_err(error, 'DEFINING COORD' )
                    error = nf90_put_att(ncid, id_vars2(k), "units", target_diag_units(i))
                    call netcdf_err(error, 'DEFINING UNITS' )
                    error = nf90_put_att(ncid, id_vars2(k), "description", target_diag_longname(i))
                    call netcdf_err(error, 'DEFINING LONG_NAME' )
                    error = nf90_put_att(ncid, id_vars2(k), "stagger", "")
                    call netcdf_err(error, 'DEFINING STAGGER' )
                    error = nf90_put_att(ncid, id_vars2(k), "FieldType", 104)
                    call netcdf_err(error, 'DEFINING FieldType' )
              endif
            else
              m = m+1
              field_extra3(m) = fields(i)
              if (localpet==0) then  
                    print*,"- DEFINE 3d diag ON FILE TARGET GRID ", varname 
                    error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat,dim_z, dim_time/), id_vars3_nz(m))
                    call netcdf_err(error, 'DEFINING VAR' )
                    error = nf90_put_att(ncid, id_vars3_nz(m), "MemoryOrder", "XYZ ")
                    call netcdf_err(error, 'DEFINING MEMORYORDER' )
                    error = nf90_put_att(ncid, id_vars3_nz(m), "coordinates", "XLONG XLAT XTIME")
                    call netcdf_err(error, 'DEFINING COORD' )
                    error = nf90_put_att(ncid, id_vars3_nz(m), "units", target_diag_units(i))
                    call netcdf_err(error, 'DEFINING UNITS' )
                    error = nf90_put_att(ncid, id_vars3_nz(m), "description", target_diag_longname(i))
                    call netcdf_err(error, 'DEFINING LONG_NAME' )
                    error = nf90_put_att(ncid, id_vars3_nz(m), "stagger", "")
                    call netcdf_err(error, 'DEFINING STAGGER' )
                    error = nf90_put_att(ncid, id_vars3_nz(m), "FieldType", 104)
                    call netcdf_err(error, 'DEFINING FieldType' )
              endif
            endif
       enddo
       deallocate(fields)
       !n2d = k
       n3d = m
   endif
   
   if (interp_hist) then
    if(n_hist_fields_2d_cons>0) then
        allocate(fields(n_hist_fields_2d_cons))
        call ESMF_FieldBundleGet(target_hist_bundle_2d_cons, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=error) 
       if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
              call error_handler("IN FieldBundleGet", error)
       do i = 1, n_hist_fields_2d_cons 
            k = k+1
            call ESMF_FieldGet(fields(i),name=varname,rc=error)
            if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
                  call error_handler("IN FieldGet", error)
      
            if (localpet==0) print*,"- DEFINE ON FILE TARGET GRID ", varname
            field_write_2d(k) = fields(i)
            if (localpet==0) then
                error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_vars2(k))
                call netcdf_err(error, 'DEFINING VAR' )
                error = nf90_put_att(ncid, id_vars2(k), "MemoryOrder", "XY " )
                call netcdf_err(error, 'DEFINING MEMORYORDER' )
                error = nf90_put_att(ncid, id_vars2(k), "coordinates", "XLONG XLAT XTIME")
                call netcdf_err(error, 'DEFINING COORD' )
                error = nf90_put_att(ncid, id_vars2(k), "units", target_hist_units_2d_cons(i))
                call netcdf_err(error, 'DEFINING UNITS' )
                error = nf90_put_att(ncid, id_vars2(k), "description", target_hist_longname_2d_cons(i))
                call netcdf_err(error, 'DEFINING LONG_NAME' )   
                error = nf90_put_att(ncid, id_vars2(k), "stagger", "")
                call netcdf_err(error, 'DEFINING STAGGE' )
                error = nf90_put_att(ncid, id_vars2(k), "FieldType", 104)
                call netcdf_err(error, 'DEFINING FieldType' )
            endif
       enddo
       deallocate(fields)
    endif
    
    if(n_hist_fields_2d_patch>0) then
        allocate(fields(n_hist_fields_2d_patch))
        call ESMF_FieldBundleGet(target_hist_bundle_2d_patch, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=error) 
       if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
              call error_handler("IN FieldBundleGet", error)
       do i = 1, n_hist_fields_2d_patch 
            k = k+1
            call ESMF_FieldGet(fields(i),name=varname,rc=error)
            if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
                  call error_handler("IN FieldGet", error)
      
            if (localpet==0) print*,"- DEFINE ON FILE TARGET GRID ", varname
            field_write_2d(k) = fields(i)
            if (localpet==0) then
                error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_vars2(k))
                call netcdf_err(error, 'DEFINING VAR' )
                error = nf90_put_att(ncid, id_vars2(k), "MemoryOrder", "XY " )
                call netcdf_err(error, 'DEFINING MEMORYORDER' )
                error = nf90_put_att(ncid, id_vars2(k), "coordinates", "XLONG XLAT XTIME")
                call netcdf_err(error, 'DEFINING COORD' )   
                error = nf90_put_att(ncid, id_vars2(k), "units", target_hist_units_2d_patch(i))
                call netcdf_err(error, 'DEFINING UNITS' )
                error = nf90_put_att(ncid, id_vars2(k), "description", target_hist_longname_2d_patch(i))
                call netcdf_err(error, 'DEFINING LONG_NAME' )   
                error = nf90_put_att(ncid, id_vars2(k), "stagger", "")
                call netcdf_err(error, 'DEFINING STAGGE' )
                error = nf90_put_att(ncid, id_vars2(k), "FieldType", 104)
                call netcdf_err(error, 'DEFINING FieldType' )
            endif
       enddo
       deallocate(fields)
    endif
    
    if(n_hist_fields_2d_nstd>0) then
        allocate(fields(n_hist_fields_2d_nstd))
        call ESMF_FieldBundleGet(target_hist_bundle_2d_nstd, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=error) 
       if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
              call error_handler("IN FieldBundleGet", error)
       do i = 1, n_hist_fields_2d_nstd
            k = k+1
            call ESMF_FieldGet(fields(i),name=varname,rc=error)
            if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
                  call error_handler("IN FieldGet", error)
            field_write_2d(k) = fields(i)
            if (localpet==0) then
                print*,"- DEFINE ON FILE TARGET GRID ", varname
                error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_vars2(k))
                call netcdf_err(error, 'DEFINING VAR' )
                error = nf90_put_att(ncid, id_vars2(k), "MemoryOrder", "XY " )
                call netcdf_err(error, 'DEFINING MEMORYORDER' )
                error = nf90_put_att(ncid, id_vars2(k), "coordinates", "XLONG XLAT XTIME")
                call netcdf_err(error, 'DEFINING COORD' )
                error = nf90_put_att(ncid, id_vars2(k), "units", target_hist_units_2d_nstd(i))
                call netcdf_err(error, 'DEFINING UNITS' )
                error = nf90_put_att(ncid, id_vars2(k), "description", target_hist_longname_2d_nstd(i))
                call netcdf_err(error, 'DEFINING LONG_NAME' )   
                error = nf90_put_att(ncid, id_vars2(k), "stagger", "")
                call netcdf_err(error, 'DEFINING STAGGER' )
                error = nf90_put_att(ncid, id_vars2(k), "FieldType", 104)
                call netcdf_err(error, 'DEFINING FieldType' )
            endif   
       enddo
       deallocate(fields)
    endif
    
    if (n_hist_fields_soil>0) then
        allocate(fields(n_hist_fields_soil))
        call ESMF_FieldBundleGet(target_hist_bundle_soil, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=error) 
       if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
              call error_handler("IN FieldBundleGet", error)
        do i = 1, n_hist_fields_soil 
            call ESMF_FieldGet(fields(i),name=varname,rc=error)
            if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
                  call error_handler("IN FieldGet", error)
            if (localpet==0) then
                print*,"- DEFINE ON FILE TARGET GRID ", varname
                error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat, dim_soil, dim_time/), id_vars_soil(i))
                call netcdf_err(error, 'DEFINING VAR' )
                error = nf90_put_att(ncid, id_vars_soil(i), "MemoryOrder", "XYZ ")
                call netcdf_err(error, 'DEFINING MEMORYORDER' )
                error = nf90_put_att(ncid, id_vars_soil(i), "coordinates", "XLONG XLAT XTIME")
                call netcdf_err(error, 'DEFINING COORD' )
                error = nf90_put_att(ncid, id_vars_soil(i), "units", target_hist_units_soil(i))
                call netcdf_err(error, 'DEFINING UNITS' )
                error = nf90_put_att(ncid, id_vars_soil(i), "description", target_hist_longname_soil(i))
                call netcdf_err(error, 'DEFINING LONG_NAME' )   
                error = nf90_put_att(ncid, id_vars_soil(i), "stagger", "")
                call netcdf_err(error, 'DEFINING STAGGER' )
                error = nf90_put_att(ncid, id_vars_soil(i), "FieldType", 104)
                call netcdf_err(error, 'DEFINING FieldType' )
            endif   
        enddo
        deallocate(fields)
    endif
    n2d = k
    if (n_hist_fields_3d_nz>0) then
        allocate(fields(n_hist_fields_3d_nz))
        call ESMF_FieldBundleGet(target_hist_bundle_3d_nz, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=error) 
       if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
              call error_handler("IN FieldBundleGet", error)
        
        do i = 1, n_hist_fields_3d_nz  
            call ESMF_FieldGet(fields(i),name=varname,rc=error)
            if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
                  call error_handler("IN FieldGet", error)
            if (localpet==0) then
                if (wrf_mod_vars .and. trim(varname)=='U') then
                    print*,"- DEFINE ON FILE STAGGERED TARGET GRID ", varname
                    error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon_stag,dim_lat,dim_z, dim_time/), id_vars3_nz(i+n3d))
                    call netcdf_err(error, 'DEFINING VAR' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d), "MemoryOrder", "XYZ ")
                    call netcdf_err(error, 'DEFINING MEMORYORDER' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d), "coordinates", "XLONG_U XLAT_U XTIME")
                    call netcdf_err(error, 'DEFINING COORD' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d), "units", target_hist_units_3d_nz(i))
                    call netcdf_err(error, 'DEFINING UNITS' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d), "description", target_hist_longname_3d_nz(i))
                    call netcdf_err(error, 'DEFINING LONG_NAME' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d), "stagger", "X")
                    call netcdf_err(error, 'DEFINING STAGGER' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d),"FieldType", 104)
                    call netcdf_err(error, 'DEFINING FieldType' )
                elseif (wrf_mod_vars .and. trim(varname)=='V') then
                    print*,"- DEFINE ON FILE STAGGERED TARGET GRID ", varname
                    error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat_stag,dim_z, dim_time/), id_vars3_nz(i+n3d))
                    call netcdf_err(error, 'DEFINING VAR' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d), "MemoryOrder", "XYZ ")
                    call netcdf_err(error, 'DEFINING MEMORYORDER' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d), "coordinates", "XLONG_V XLAT_V XTIME")
                    call netcdf_err(error, 'DEFINING COORD' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d), "units", target_hist_units_3d_nz(i))
                    call netcdf_err(error, 'DEFINING UNITS' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d), "description", target_hist_longname_3d_nz(i))
                    call netcdf_err(error, 'DEFINING LONG_NAME' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d), "stagger", "Y")
                    call netcdf_err(error, 'DEFINING STAGGER' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d),"FieldType", 104)
                    call netcdf_err(error, 'DEFINING FieldType' )
                else                
                    print*,"- DEFINE ON FILE TARGET GRID ", varname
                    error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat,dim_z, dim_time/), id_vars3_nz(i+n3d))
                    call netcdf_err(error, 'DEFINING VAR' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d), "MemoryOrder", "XYZ ")
                    call netcdf_err(error, 'DEFINING MEMORYORDER' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d), "coordinates", "XLONG XLAT XTIME")
                    call netcdf_err(error, 'DEFINING COORD' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d), "units", target_hist_units_3d_nz(i))
                    call netcdf_err(error, 'DEFINING UNITS' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d), "description", target_hist_longname_3d_nz(i))
                    call netcdf_err(error, 'DEFINING LONG_NAME' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d), "stagger", "")
                    call netcdf_err(error, 'DEFINING STAGGER' )
                    error = nf90_put_att(ncid, id_vars3_nz(i+n3d),"FieldType", 104)
                    call netcdf_err(error, 'DEFINING FieldType' )
                endif

                if (wrf_mod_vars .and. trim(varname)=='MUB') then
                    print*,"- DEFINE ON FILE STAGGERED TARGET GRID MU"
                    error = nf90_def_var(ncid, 'MU', NF90_FLOAT, (/dim_lon,dim_lat,dim_z, dim_time/), id_mu)
                    call netcdf_err(error, 'DEFINING VAR' )
                    error = nf90_put_att(ncid, id_mu, "MemoryOrder", "XYZ ")
                    call netcdf_err(error, 'DEFINING MEMORYORDER' )
                    error = nf90_put_att(ncid, id_mu, "coordinates", "XLONG XLAT XTIME")
                    call netcdf_err(error, 'DEFINING COORD')
                    error = nf90_put_att(ncid, id_mu, "units", target_hist_units_3d_nz(i))
                    call netcdf_err(error, 'DEFINING UNITS')
                    error = nf90_put_att(ncid, id_mu, "description", 'Perturbation '//target_hist_longname_3d_nz(i))
                    call netcdf_err(error, 'DEFINING LONG_NAME' )
                    error = nf90_put_att(ncid, id_mu, "stagger", "")
                    call netcdf_err(error, 'DEFINING STAGGER' )
                    error = nf90_put_att(ncid, id_mu,"FieldType", 104)
                    call netcdf_err(error, 'DEFINING FieldType' )
                 endif        

                 if (wrf_mod_vars .and. trim(varname)=='P_HYD') then
                    print*,"- DEFINE ON FILE STAGGERED TARGET GRID P_TOP"
                    error = nf90_def_var(ncid, 'P_TOP', NF90_FLOAT, (/dim_time/), id_ptop)
                    call netcdf_err(error, 'DEFINING VAR' )
                    error = nf90_put_att(ncid, id_ptop, "MemoryOrder", "0 ")
                    call netcdf_err(error, 'DEFINING MEMORYORDER' )
                    error = nf90_put_att(ncid, id_ptop, "units", target_hist_units_3d_nz(i))
                    call netcdf_err(error, 'DEFINING UNITS')
                    error = nf90_put_att(ncid, id_ptop, "description", 'PRESSURE TOP OF THE MODEL')
                    call netcdf_err(error, 'DEFINING LONG_NAME' )
                    error = nf90_put_att(ncid, id_ptop, "stagger", "")
                    call netcdf_err(error, 'DEFINING STAGGER' )
                    error = nf90_put_att(ncid, id_ptop,"FieldType", 104)
                    call netcdf_err(error, 'DEFINING FieldType' )
                 endif 
            endif   
        enddo
        deallocate(fields)
    endif
    
    if (n_hist_fields_3d_nzp1>0) then
        allocate(fields(n_hist_fields_3d_nzp1))
        call ESMF_FieldBundleGet(target_hist_bundle_3d_nzp1, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=error) 
       if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
              call error_handler("IN FieldBundleGet", error)
        do i = 1, n_hist_fields_3d_nzp1
            call ESMF_FieldGet(fields(i),name=varname,rc=error)
            if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
                  call error_handler("IN FieldGet", error)
            if (localpet==0) then
                print*,"- DEFINE ON FILE TARGET GRID ", varname
                error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat,dim_zp1, dim_time/), id_vars3_nzp1(i))
                call netcdf_err(error, 'DEFINING VAR' )
                error = nf90_put_att(ncid, id_vars3_nzp1(i), "MemoryOrder", "XYZ ")
                call netcdf_err(error, 'DEFINING MEMORYORDER' )
                error = nf90_put_att(ncid, id_vars3_nzp1(i), "coordinates", "XLONG XLAT XTIME")
                call netcdf_err(error, 'DEFINING COORD' )
                if (wrf_mod_vars .and. trim(varname)=='PHB') then
                    error = nf90_put_att(ncid, id_vars3_nzp1(i), "units", "gpm")
                    call netcdf_err(error, 'DEFINING UNITS' )
                    error = nf90_put_att(ncid, id_vars3_nzp1(i), "description", "Base Geopotential Height")
                    call netcdf_err(error, 'DEFINING LONG_NAME' )   
                else
                    error = nf90_put_att(ncid, id_vars3_nzp1(i), "units", target_hist_units_3d_nzp1(i))
                    call netcdf_err(error, 'DEFINING UNITS' )
                    error = nf90_put_att(ncid, id_vars3_nzp1(i), "description", target_hist_longname_3d_nzp1(i))
                    call netcdf_err(error, 'DEFINING LONG_NAME' )   
                endif   
                error = nf90_put_att(ncid, id_vars3_nzp1(i),"stagger", "Z")
                call netcdf_err(error, 'DEFINING STAGGER' )
                error = nf90_put_att(ncid, id_vars3_nzp1(i),"FieldType", 104)
                call netcdf_err(error, 'DEFINING FieldType' )   
            
            
                if (wrf_mod_vars .and. trim(varname)=='PHB') then
                    print*,"- DEFINE ON FILE STAGGERED TARGET GRID MU"
                    error = nf90_def_var(ncid, 'PH', NF90_FLOAT, (/dim_lon,dim_lat,dim_zp1, dim_time/), id_ph)
                    call netcdf_err(error, 'DEFINING VAR' )
                    error = nf90_put_att(ncid, id_ph, "MemoryOrder", "XYZ ")
                    call netcdf_err(error, 'DEFINING MEMORYORDER' )
                    error = nf90_put_att(ncid, id_ph, "coordinates", "XLONG XLAT XTIME")
                    call netcdf_err(error, 'DEFINING COORD')
                    error = nf90_put_att(ncid, id_ph, "units", "gpm")
                    call netcdf_err(error, 'DEFINING UNITS')
                    error = nf90_put_att(ncid, id_ph, "description", 'Perturbation Geopotential Height')
                    call netcdf_err(error, 'DEFINING LONG_NAME' )
                    error = nf90_put_att(ncid, id_ph, "stagger", "Z")
                    call netcdf_err(error, 'DEFINING STAGGER' )
                    error = nf90_put_att(ncid, id_ph,"FieldType", 104)
                    call netcdf_err(error, 'DEFINING FieldType' )
                endif
            endif 
        enddo
        deallocate(fields)
    endif
   endif !write hist    
  
 if (localpet==0) then 
   error = nf90_enddef(ncid, header_buffer_val,4,0,4)
   call netcdf_err(error, 'DEFINING HEADER' )
 endif

!--- write fields 

!  longitude

   if (localpet==0) print*,"- CALL FieldGather FOR TARGET GRID LONGITUDE"
   call ESMF_FieldGather(longitude_target_grid, dum2d, rootPet=0, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", error)


 if (localpet ==0) then
   dum2dt(:,:,1) = dum2d
   error = nf90_put_var( ncid, id_lon, dum2dt, count=(/i_target,j_target,1/))
   call netcdf_err(error, 'WRITING LONGITUDE RECORD' )
 endif

!  latitude

   if (localpet==0) print*,"- CALL FieldGather FOR TARGET GRID LATITUDE"
   call ESMF_FieldGather(latitude_target_grid, dum2d, rootPet=0, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", error)


 if (localpet ==0) then
   dum2dt(:,:,1) = dum2d
   error = nf90_put_var( ncid, id_lat, dum2dt,count=(/i_target,j_target,1/))
   call netcdf_err(error, 'WRITING LATITUDE RECORD' )
 endif

!  longitude on u grid

 if (localpet ==0) then
   dum2dtu(:,:,1) = longitude_u
   error = nf90_put_var( ncid, id_lonu, dum2dtu, count=(/ip1_target,j_target,1/))
   call netcdf_err(error, 'WRITING XLONG_U RECORD' )
 endif

!  latitude on u grid

 if (localpet ==0) then
   dum2dtu(:,:,1) = latitude_u
   error = nf90_put_var( ncid, id_latu, dum2dtu,count=(/ip1_target,j_target,1/))
   call netcdf_err(error, 'WRITING XLAT_U RECORD' )
 endif

!  longitude on v grid

 if (localpet ==0) then
   dum2dtv(:,:,1) = longitude_v
   error = nf90_put_var( ncid, id_lonv, dum2dtv, count=(/i_target,jp1_target,1/))
   call netcdf_err(error, 'WRITING XLONG_U RECORD' )
 endif

!  latitude on v grid

 if (localpet ==0) then
   dum2dtv(:,:,1) = latitude_v
   error = nf90_put_var( ncid, id_latv, dum2dtv,count=(/i_target,jp1_target,1/))
   call netcdf_err(error, 'WRITING XLAT_U RECORD' )
 endif
 deallocate(longitude_u,latitude_u,longitude_v,latitude_v) 
!  z_s

 if (localpet==0) print*,"- WRITE TO FILE TARGET GRID Z_S"
 if (localpet ==0) then
   error = nf90_put_var( ncid, id_zs, zs_target_grid,count=(/nsoil_input,1/))
   call netcdf_err(error, 'WRITING ZS RECORD' )
 endif

!  hgt

 if (localpet==0)  print*,"- CALL FieldGather FOR TARGET GRID LATITUDE"
   call ESMF_FieldGather(hgt_target_grid, dum2d, rootPet=0, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGather", error)

 if (localpet==0) print*,"- WRITE TO FILE TARGET GRID HGT"
 if (localpet ==0) then
   dum2dt(:,:,1) = dum2d
   error = nf90_put_var( ncid, id_hgt, dum2dt,count=(/i_target,j_target,1/))
   call netcdf_err(error, 'WRITING HGT RECORD' )
 endif

!  times

 if (localpet==0)  print*,"- WRITE TO FILE TARGET GRID Times"
 if (localpet ==0) then
   tempstr(1,:) = valid_time(1,1:Datestrlen)
   error = nf90_put_var( ncid, id_times, tempstr, start = (/1,1/), count=(/Datestrlen,1/))
   call netcdf_err(error, 'WRITING TIMES RECORD' )
 endif

!  xtime

 if (localpet==0)  print*,"- WRITE TO FILE TARGET GRID ITIMESTEP"
 if (localpet ==0) then
   sy = substr(start_time,1,4)
   sm = substr(start_time,6,7)
   sd = substr(start_time,9,10)
   sh = substr(start_time,12,13)
   smi = substr(start_time,15,16)
   ss = substr(start_time,18,19)

   vy = substr(valid_time(1,1),1,4)
   vm = substr(valid_time(1,1),6,7)
   vd = substr(valid_time(1,1),9,10)
   vh = substr(valid_time(1,1),12,13)
   vmi = substr(valid_time(1,1),15,16)
   vs = substr(valid_time(1,1),18,19)
   xtime_dt = datetime(sy,sm,sd,sh,smi,ss) - datetime(vy,vm,vd,vh,vmi,vs)

   error = nf90_put_var( ncid, id_xtime, (/xtime_dt%total_seconds()/60.0/), count=(/1/))
   call netcdf_err(error, 'WRITING XTIME RECORD' )
 endif

 !  itimestep

 if (localpet==0) print*,"- WRITE TO FILE TARGET GRID ITIMESTEP"
 if (localpet ==0) then
   error = nf90_put_var( ncid, id_itime, (/int(xtime_dt%total_seconds()/config_dt)/), count=(/1/))
   call netcdf_err(error, 'WRITING ITIMESTEP RECORD' )
 endif
 
 deallocate(dumsmall)
 
 !  2d fields

 do i = 1, n2d
    
    call ESMF_FieldGet(field_write_2d(i), name=varname, rc=error)
    if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGet", error)
      
    if (localpet==0) print*,"- CALL FieldGather FOR TARGET GRID ", trim(varname) 
    call ESMF_FieldGather(field_write_2d(i), dum2d, rootPet=0, rc=error)
    if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", error)

    if (localpet==0) then
        print*, "- WRITE TO FILE ", trim(varname)
        dum2dt(:,:,1) = dum2d
        error = nf90_put_var( ncid, id_vars2(i), dum2dt,count=(/i_target,j_target,1/))
        call netcdf_err(error, 'WRITING RECORD' )
    endif
 enddo
 deallocate(field_write_2d)
 deallocate(dum2d, dum2dt)

 !    3d fields from diaglist
 if (n3d>0) then
 do i = 1, n3d
    call ESMF_FieldGet(field_extra3(i), name=varname, rc=error)
    if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", error)

    if (localpet==0) print*,"- CALL FieldGather FOR TARGET GRID ", trim(varname)
    call ESMF_FieldGather(field_extra3(i), dum3d, rootPet=0, rc=error)
    if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGather", error)

    if (localpet==0) then
        print*, trim(varname), minval(dum3d), maxval(dum3d)
        dum3dt(:,:,:,1) = dum3d
        error = nf90_put_var( ncid, id_vars3_nz(i), dum3dt,count=(/i_target,j_target,nz_input,1/))
        call netcdf_err(error, 'WRITING RECORD' )
    endif
 enddo
 endif
 ! 3d soil fields
 
 if (interp_hist .and. n_hist_fields_soil>0) then
     allocate(fields(n_hist_fields_soil))
     call ESMF_FieldBundleGet(target_hist_bundle_soil, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=error) 
       if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
              call error_handler("IN FieldBundleGet", error)
              
     do i = 1, n_hist_fields_soil  
        call ESMF_FieldGet(fields(i), name=varname, rc=error)
        if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldGet", error)
      
        if (localpet==0) print*,"- CALL FieldGather FOR TARGET GRID ", trim(varname) 
        call ESMF_FieldGather(fields(i), dumsoil, rootPet=0, rc=error)
        if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldGather", error)


        if (localpet==0) then
            print*, trim(varname), minval(dumsoil), maxval(dumsoil)
            dumsoilt(:,:,:,1) = dumsoil
            error = nf90_put_var( ncid, id_vars_soil(i), dumsoilt,count=(/i_target,j_target,nsoil_input,1/))
            call netcdf_err(error, 'WRITING RECORD' )
        endif
     enddo
     deallocate(fields)
 endif
 deallocate(dumsoil,dumsoilt)
 
 ! 3d nz fields
 
 if (interp_hist .and. n_hist_fields_3d_nz>0) then
     allocate(fields(n_hist_fields_3d_nz))
     call ESMF_FieldBundleGet(target_hist_bundle_3d_nz, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=error) 
       if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
              call error_handler("IN FieldBundleGet", error)
              
     do n = 1, n_hist_fields_3d_nz  
        call ESMF_FieldGet(fields(n), name=varname, rc=error)
        if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldGet", error)
        call ESMF_FieldGather(fields(n), dum3d, rootPet=0, rc=error)
         if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
             call error_handler("IN FieldGather", error)
        if (localpet==0) then
            if (wrf_mod_vars .and. trim(varname)=='U') then
                allocate(dum3dtmp(i_target+1,j_target,nz_input))
                call create_stagger(dum3d,nz_input, j_target, i_target, &
                                    nz_input, j_target, i_target+1, 1,dum3dtmp)
                error = nf90_put_var( ncid, id_vars3_nz(n+n3d), dum3dtmp,count=(/i_target+1,j_target,nz_input,1/))
                call netcdf_err(error, 'WRITING RECORD' )
                print*, trim(varname), minval(dum3dtmp), maxval(dum3dtmp)
                deallocate(dum3dtmp)
            elseif (wrf_mod_vars .and. trim(varname)=='V') then
                allocate(dum3dtmp(i_target,j_target+1,nz_input))
                call create_stagger(dum3d,nz_input, j_target, i_target,  &
                                    nz_input, j_target+1, i_target, 2,dum3dtmp)
                print*, trim(varname), minval(dum3dtmp), maxval(dum3dtmp)
                error = nf90_put_var( ncid, id_vars3_nz(n+n3d), dum3dtmp, &
                                        count=(/i_target,j_target+1,nz_input,1/))
                call netcdf_err(error, 'WRITING RECORD' )
                deallocate(dum3dtmp)
            else
                if (wrf_mod_vars .and. trim(varname)=='T') then
                    do i = 1,i_target
                    do j = 1,j_target
                      if (dum3d(i,j,1) < 10.0_esmf_kind_r8) continue
                      dum3dt(i,j,:,1) = dum3d(i,j,:) - 300.0
                    enddo
                    enddo
                else    
                    dum3dt(:,:,:,1) = dum3d(:,:,:)
                endif
                print*, trim(varname), minval(dum3dt), maxval(dum3dt)
                error = nf90_put_var( ncid, id_vars3_nz(n+n3d), dum3dt, &
                                        count=(/i_target,j_target,nz_input,1/))
                call netcdf_err(error, 'WRITING RECORD' )
                
                if (wrf_mod_vars .and. trim(varname)=='MUB') then
                    dum3dt(:,:,:,1) = 0.0_esmf_kind_r8
                    print*, 'MU', minval(dum3d), maxval(dum3d)
                    error = nf90_put_var( ncid, id_mu, dum3dt, &
                                        count=(/i_target,j_target,nz_input,1/))
                    call netcdf_err(error, 'WRITING RECORD' )
                endif

                if (wrf_mod_vars .and. trim(varname)=='P_HYD') then
                    dum1d(1) = maxval(dum3d)
                    print*, "initial p_top = ", dum1d(1)
                    do i = 1,i_target
                    do j = 1,j_target
                      if (dum3d(i,j,nz_input) >= 10.0_esmf_kind_r8) then
                          dum1d(1) = min(dum3d(i,j,nz_input)* 0.80_esmf_kind_r8, dum1d(1))
                      endif
                    enddo
                    enddo
                    print*, "P_TOP", dum1d(1)
                    error = nf90_put_var( ncid, id_ptop, dum1d, &
                                        count=(/1/))
                    call netcdf_err(error, 'WRITING RECORD' )
                endif
            endif
        endif
     enddo
     deallocate(fields)
 endif
 
  ! 3d nzp1 fields
 
 if (interp_hist .and. n_hist_fields_3d_nzp1>0) then
     allocate(fields(n_hist_fields_3d_nzp1))
     call ESMF_FieldBundleGet(target_hist_bundle_3d_nzp1, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=error) 
       if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
              call error_handler("IN FieldBundleGet", error)
              
     do i = 1, n_hist_fields_3d_nzp1  
        call ESMF_FieldGet(fields(i), name=varname, rc=error)
        if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldGet", error)
      
        if (localpet==0) print*,"- CALL FieldGather FOR TARGET GRID ", trim(varname) 
        call ESMF_FieldGather(fields(i), dum3dp1, rootPet=0, rc=error)
        if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldGather", error)
          
        if (localpet==0) then  
            if(trim(varname)=='PHB') then
                do n = 1,i_target
                do j = 1,j_target
                do k = 2,nzp1_input
                    dum3dt(n,j,k-1,1) = 0.5 * ( dum3dp1(n,j,k) + dum3dp1(n,j,k-1) )
                enddo
                enddo
                enddo
            
                error = nf90_put_var( ncid, id_z, dum3dt, count=(/i_target,j_target,nz_input,1/))
                call netcdf_err(error, 'WRITING RECORD' )
                
                dum3dp1 = dum3dp1 * 9.81
            endif
         
            print*, trim(varname), minval(dum3dp1), maxval(dum3dp1)
            dum3dp1t(:,:,:,1) = dum3dp1
            error = nf90_put_var( ncid, id_vars3_nzp1(i), dum3dp1t, &
                                        count=(/i_target,j_target,nz_input+1,1/))
            call netcdf_err(error, 'WRITING RECORD' )
            
            if (wrf_mod_vars .and. trim(varname)=='PHB') then
                    dum3dt(:,:,:,1) = 0.0_esmf_kind_r8
                    error = nf90_put_var( ncid, id_ph, dum3dt, &
                                        count=(/i_target,j_target,nz_input+1,1/))
                    call netcdf_err(error, 'WRITING RECORD' )
            endif


        endif
     enddo
     deallocate(fields)
 endif
 

 deallocate(dum3d, dum3dp1, dum3dt, dum3dp1t, dum1d)
 deallocate(id_vars2, id_vars3_nz, id_vars3_nzp1, id_vars_soil)
 deallocate(target_hist_longname_2d_cons, target_hist_longname_2d_nstd)
 deallocate(target_hist_longname_2d_patch, target_hist_longname_3d_nz)
 deallocate(target_hist_longname_3d_nzp1, target_diag_longname)
 deallocate(target_hist_units_2d_cons, target_hist_units_2d_nstd)
 deallocate(target_hist_units_2d_patch, target_hist_units_3d_nz)
 deallocate(target_hist_units_3d_nzp1, target_diag_units)

 if (localpet == 0) error = nf90_close(ncid)

 end subroutine write_to_file
 
 subroutine create_stagger(array_in,nz_in, ny_in, nx_in,  &
                            nz_out, ny_out, nx_out, stag_dim, &
                            array_out)
    use esmf
    
    implicit none
    integer, intent(IN)               :: nz_in, ny_in, nx_in, nz_out, ny_out, nx_out, stag_dim
    real(esmf_kind_r8), intent(IN)    :: array_in(nx_in,ny_in,nz_in)
    real(esmf_kind_r8), intent(INOUT) :: array_out(nx_out,ny_out,nz_out)
    
    integer                           :: i, j, k
    do k = 1, nz_out
        if (stag_dim == 1) then
            do j = 1,ny_out
                array_out(1,j,k) = array_in(1,j,k)
                array_out(2,j,k) = array_in(1,j,k)
                do i = 3,nx_out
                    array_out(i,j,k) = 2.0*array_in(i-1,j,k)-array_out(i-1,j,k)
                enddo
            enddo
        elseif(stag_dim == 2) then
            do i = 1,nx_out
                array_out(i,1,k) = array_in(i,1,k)
                array_out(i,2,k) = array_in(i,1,k)
                do j = 3,ny_out
                    array_out(i,j,k) = 2.0*array_in(i,j-1,k)-array_out(i,j-1,k)
                enddo
            enddo
        endif
    enddo
 end subroutine create_stagger

 elemental function substr( s, a, b ) result( val )
        character(*), intent(in) :: s
        integer,      intent(in) :: a, b
        character(len(s)) :: res
        integer           :: val

        res = s( a : b )

        read(res,*) val
       
    endfunction
 
 end module write_data
