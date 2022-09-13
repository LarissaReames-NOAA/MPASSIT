 module write_data

 use program_setup, only   : output_file

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

 use program_setup, only           : interp_diag, interp_hist

 use model_grid, only              : i_target, j_target, &
                                     ip1_target, jp1_target, &
                                     nz_input, nzp1_input, &
                                     nsoil_input, &
                                     dx, start_time, &
                                     lsm_scheme, mp_scheme, &
                                     conv_scheme, &
                                     longitude_target_grid, &
                                     latitude_target_grid, &
                                     zs_target_grid, &
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


 integer                          :: error, ncid, n, rc, i, j, k
 integer                          :: header_buffer_val = 16384
 integer                          :: dim_time, dim_lon, dim_lat, dim_z, dim_zp1, dim_soil
 integer                          :: dim_lonp, dim_latp
 integer                          :: id_lat, id_lon, id_z, id_zs
 integer                          :: n2d
 integer, allocatable             :: id_vars2(:), id_vars3_nz(:), id_vars3_nzp1(:), &
 									 id_vars_soil(:)

 real(esmf_kind_r8), allocatable  :: dum2d(:,:), dum2dt(:,:,:), &
 									 dum3d(:,:,:), dum3dt(:,:,:,:), &
 									 dum3dp1(:,:,:), dum3dp1t(:,:,:,:), &
  									 dumsoil(:,:,:), dumsoilt(:,:,:,:), &
  									 dumsmall(:,:)
 
 type(esmf_field), allocatable    :: fields(:), field_write_2d(:)

 n2d = n_diag_fields + n_hist_fields_2d_patch + n_hist_fields_2d_nstd + n_hist_fields_2d_cons
 allocate(field_write_2d(n2d),id_vars2(n2d))
 allocate(id_vars3_nz(n_hist_fields_3d_nz))
 allocate(id_vars3_nzp1(n_hist_fields_3d_nzp1))
 allocate(id_vars_soil(n_hist_fields_soil))
 
 if (localpet ==0) then
   allocate(dumsmall(nsoil_input,1))
   allocate(dum2d(i_target,j_target))
   allocate(dum2dt(i_target,j_target,1))
   allocate(dum3d(i_target,j_target,nz_input))
   allocate(dum3dt(i_target,j_target,nz_input,1))
   allocate(dum3dp1(i_target,j_target,nzp1_input))
   allocate(dum3dp1t(i_target,j_target,nzp1_input,1))
   allocate(dumsoil(i_target,j_target,nsoil_input))
   allocate(dumsoilt(i_target,j_target,nsoil_input,1))
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
 endif
 
if (localpet == 0) then

!--- open the file
   error = nf90_create(output_file, NF90_NETCDF4, ncid)
   call netcdf_err(error, 'CREATING FILE '//trim(output_file) )

!--- define dimension
   error = nf90_def_dim(ncid, 'Time', 1 , dim_time)
   call netcdf_err(error, 'DEFINING Time DIMENSION' )
   error = nf90_def_dim(ncid, 'west_east', i_target, dim_lon)
   call netcdf_err(error, 'DEFINING LON DIMENSION' )
   error = nf90_def_dim(ncid, 'south_north', j_target, dim_lat)
   call netcdf_err(error, 'DEFINING LAT DIMENSION' )
   error = nf90_def_dim(ncid, 'bottom_top', nz_input, dim_z)
   call netcdf_err(error, 'DEFINING VERTICAL DIMENSION' )
   error = nf90_def_dim(ncid, 'bottom_top_stag', nzp1_input, dim_zp1)
   call netcdf_err(error, 'DEFINING VERTICALP1 DIMENSION' )
   error = nf90_def_dim(ncid, 'soil_layers_stag', nsoil_input, dim_soil)
   call netcdf_err(error, 'DEFINING VERTICALP1 DIMENSION' )
   
 !--- define global attributes
   error = nf90_put_att(ncid, NF90_GLOBAL, 'WEST-EAST_GRID_DIMENSION', i_target)
   call netcdf_err(error, 'DEFINING WEST-EAST GRID DIMENSION GLOBAL ATTRIBUTE')
   
   error = nf90_put_att(ncid, NF90_GLOBAL, 'NORTH-SOUTH_GRID_DIMENSION', j_target)
   call netcdf_err(error, 'DEFINING NORTH-SOUTH GRID DIMENSION GLOBAL ATTRIBUTE')
   
   error = nf90_put_att(ncid, NF90_GLOBAL, 'BOTTOM-TOP_GRID_DIMENSION', nz_input)
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

!--- define fields

   error = nf90_def_var(ncid, 'XLONG', NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_lon)
   call netcdf_err(error, 'DEFINING GEOLON FIELD' )
   error = nf90_put_att(ncid, id_lon, "long_name", "Longitude")
   call netcdf_err(error, 'DEFINING GEOLON NAME' )
   error = nf90_put_att(ncid, id_lon, "units", "degrees_east")
   call netcdf_err(error, 'DEFINING GEOLON UNITS' )
   error = nf90_put_att(ncid, id_lon, "MemoryOrder", "XY")
	call netcdf_err(error, 'DEFINING MEMORYORDER' )
	error = nf90_put_att(ncid, id_lon, "coordinates", "XLONG XLAT")
	call netcdf_err(error, 'DEFINING COORD' )

   error = nf90_def_var(ncid, 'XLAT', NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_lat)
   call netcdf_err(error, 'DEFINING GEOLAT FIELD' )
   error = nf90_put_att(ncid, id_lat, "long_name", "Latitude")
   call netcdf_err(error, 'DEFINING GEOLAT NAME' )
   error = nf90_put_att(ncid, id_lat, "units", "degrees_north")
   call netcdf_err(error, 'DEFINING GEOLAT UNITS' )
   error = nf90_put_att(ncid, id_lat, "MemoryOrder", "XY")
	call netcdf_err(error, 'DEFINING MEMORYORDER' )
	error = nf90_put_att(ncid, id_lat, "coordinates", "XLONG XLAT")
	call netcdf_err(error, 'DEFINING COORD' )
   
   error = nf90_def_var(ncid, 'Z_C', NF90_FLOAT, (/dim_lon,dim_lat,dim_zp1, dim_time/), id_z)
   call netcdf_err(error, 'DEFINING Z_C FIELD' )
   error = nf90_put_att(ncid, id_z, "long_name", "Layer center height above mean sea level")
   call netcdf_err(error, 'DEFINING Z_C NAME' )
   error = nf90_put_att(ncid, id_z, "units", "m AMSL")
   call netcdf_err(error, 'DEFINING Z_C UNITS' )
   error = nf90_put_att(ncid, id_z, "MemoryOrder", "XYZ")
	call netcdf_err(error, 'DEFINING MEMORYORDER' )
	error = nf90_put_att(ncid, id_z, "coordinates", "XLONG XLAT Z_C")
	call netcdf_err(error, 'DEFINING COORD' )
   
   error = nf90_def_var(ncid, 'ZS', NF90_FLOAT, (/dim_soil, dim_time/), id_zs)
   call netcdf_err(error, 'DEFINING ZS FIELD' )
   error = nf90_put_att(ncid, id_zs, "long_name", "DEPTHS OF CENTERS OF SOIL LAYERS")
   call netcdf_err(error, 'DEFINING ZS NAME' )
   error = nf90_put_att(ncid, id_zs, "units", "m")
   call netcdf_err(error, 'DEFINING ZS UNITS' )
   error = nf90_put_att(ncid, id_zs, "MemoryOrder", "X")
	call netcdf_err(error, 'DEFINING MEMORYORDER' )
	error = nf90_put_att(ncid, id_zs, "coordinates", "ZS")
	call netcdf_err(error, 'DEFINING COORD' )
 endif
   
   k = 0
   if (interp_diag .and. n_diag_fields>0) then
   		allocate(fields(n_diag_fields))
   		call ESMF_FieldBundleGet(target_diag_bundle, fieldList=fields, & 
 						  itemorderflag=ESMF_ITEMORDER_ADDORDER, &
 						  rc=error) 
	   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
			  call error_handler("IN FieldBundleGet", error)
		do i = 1, n_diag_fields  
			k = k+1
			call ESMF_FieldGet(fields(i),name=varname,rc=error)
			if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
				  call error_handler("IN FieldGet", error)
	  
			print*,"- DEFINE ON FILE TARGET GRID ", varname
			field_write_2d(k) = fields(i)
			if (localpet==0) then
				error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_vars2(k))
				call netcdf_err(error, 'DEFINING VAR' )
				error = nf90_put_att(ncid, id_vars2(k), "MemoryOrder", "XY")
				call netcdf_err(error, 'DEFINING MEMORYORDER' )
				error = nf90_put_att(ncid, id_vars2(k), "coordinates", "XLONG XLAT XTIME")
				call netcdf_err(error, 'DEFINING COORD' )
				error = nf90_put_att(ncid, id_vars2(k), "units", target_diag_units(i))
				call netcdf_err(error, 'DEFINING UNITS' )
				error = nf90_put_att(ncid, id_vars2(k), "long_name", target_diag_longname(i))
				call netcdf_err(error, 'DEFINING LONG_NAME' )			
			endif
	   enddo
	   deallocate(fields)
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
	  
			print*,"- DEFINE ON FILE TARGET GRID ", varname
			field_write_2d(k) = fields(i)
			if (localpet==0) then
				error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_vars2(k))
				call netcdf_err(error, 'DEFINING VAR' )
				error = nf90_put_att(ncid, id_vars2(k), "MemoryOrder", "XY")
				call netcdf_err(error, 'DEFINING MEMORYORDER' )
				error = nf90_put_att(ncid, id_vars2(k), "coordinates", "XLONG XLAT XTIME")
				call netcdf_err(error, 'DEFINING COORD' )
				error = nf90_put_att(ncid, id_vars2(k), "units", target_hist_units_2d_cons(i))
				call netcdf_err(error, 'DEFINING UNITS' )
				error = nf90_put_att(ncid, id_vars2(k), "long_name", target_hist_longname_2d_cons(i))
				call netcdf_err(error, 'DEFINING LONG_NAME' )		
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
	  
			print*,"- DEFINE ON FILE TARGET GRID ", varname
			field_write_2d(k) = fields(i)
			if (localpet==0) then
				error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_vars2(k))
				call netcdf_err(error, 'DEFINING VAR' )
				error = nf90_put_att(ncid, id_vars2(k), "MemoryOrder", "XY")
				call netcdf_err(error, 'DEFINING MEMORYORDER' )
				error = nf90_put_att(ncid, id_vars2(k), "coordinates", "XLONG XLAT XTIME")
				call netcdf_err(error, 'DEFINING COORD' )	
				error = nf90_put_att(ncid, id_vars2(k), "units", target_hist_units_2d_patch(i))
				call netcdf_err(error, 'DEFINING UNITS' )
				error = nf90_put_att(ncid, id_vars2(k), "long_name", target_hist_longname_2d_patch(i))
				call netcdf_err(error, 'DEFINING LONG_NAME' )	
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
				error = nf90_put_att(ncid, id_vars2(k), "MemoryOrder", "XY")
				call netcdf_err(error, 'DEFINING MEMORYORDER' )
				error = nf90_put_att(ncid, id_vars2(k), "coordinates", "XLONG XLAT XTIME")
				call netcdf_err(error, 'DEFINING COORD' )
				error = nf90_put_att(ncid, id_vars2(k), "units", target_hist_units_2d_nstd(i))
				call netcdf_err(error, 'DEFINING UNITS' )
				error = nf90_put_att(ncid, id_vars2(k), "long_name", target_hist_longname_2d_nstd(i))
				call netcdf_err(error, 'DEFINING LONG_NAME' )	
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
				error = nf90_put_att(ncid, id_vars_soil(i), "MemoryOrder", "XYZ")
				call netcdf_err(error, 'DEFINING MEMORYORDER' )
				error = nf90_put_att(ncid, id_vars_soil(i), "coordinates", "XLONG XLAT ZS TIME")
				call netcdf_err(error, 'DEFINING COORD' )
				error = nf90_put_att(ncid, id_vars_soil(i), "units", target_hist_units_soil(i))
				call netcdf_err(error, 'DEFINING UNITS' )
				error = nf90_put_att(ncid, id_vars_soil(i), "long_name", target_hist_longname_soil(i))
				call netcdf_err(error, 'DEFINING LONG_NAME' )	
			endif	
	   	enddo
	   	deallocate(fields)
 	endif
   
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
				print*,"- DEFINE ON FILE TARGET GRID ", varname
				error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat,dim_z, dim_time/), id_vars3_nz(i))
				call netcdf_err(error, 'DEFINING VAR' )
				error = nf90_put_att(ncid, id_vars3_nz(i), "MemoryOrder", "XYZ")
				call netcdf_err(error, 'DEFINING MEMORYORDER' )
				error = nf90_put_att(ncid, id_vars3_nz(i), "coordinates", "XLONG XLAT Z_C TIME")
				call netcdf_err(error, 'DEFINING COORD' )
				error = nf90_put_att(ncid, id_vars3_nz(i), "units", target_hist_units_3d_nz(i))
				call netcdf_err(error, 'DEFINING UNITS' )
				error = nf90_put_att(ncid, id_vars3_nz(i), "long_name", target_hist_longname_3d_nz(i))
				call netcdf_err(error, 'DEFINING LONG_NAME' )	
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
				error = nf90_put_att(ncid, id_vars3_nzp1(i), "MemoryOrder", "XYZ")
				call netcdf_err(error, 'DEFINING MEMORYORDER' )
				error = nf90_put_att(ncid, id_vars3_nzp1(i), "coordinates", "XLONG XLAT Z TIME")
				call netcdf_err(error, 'DEFINING COORD' )
				error = nf90_put_att(ncid, id_vars3_nzp1(i), "units", target_hist_units_3d_nzp1(i))
				call netcdf_err(error, 'DEFINING UNITS' )
				error = nf90_put_att(ncid, id_vars3_nzp1(i), "long_name", target_hist_longname_3d_nzp1(i))
				call netcdf_err(error, 'DEFINING LONG_NAME' )		
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

   print*,"- CALL FieldGather FOR TARGET GRID LONGITUDE"
   call ESMF_FieldGather(longitude_target_grid, dum2d, rootPet=0, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", error)


 if (localpet ==0) then
   dum2dt(:,:,1) = dum2d
   error = nf90_put_var( ncid, id_lon, dum2dt)
   call netcdf_err(error, 'WRITING LONGITUDE RECORD' )
 endif

!  latitude


   print*,"- CALL FieldGather FOR TARGET GRID LATITUDE"
   call ESMF_FieldGather(latitude_target_grid, dum2d, rootPet=0, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", error)


 if (localpet ==0) then
   dum2dt(:,:,1) = dum2d
   error = nf90_put_var( ncid, id_lat, dum2dt)
   call netcdf_err(error, 'WRITING LATITUDE RECORD' )
 endif
 
!  z_s

 print*,"- WRITE TO FILE TARGET GRID Z_S"
 if (localpet ==0) then
   error = nf90_put_var( ncid, id_zs, zs_target_grid)
   call netcdf_err(error, 'WRITING ZS RECORD' )
 endif
 
 deallocate(dumsmall)
 
 !  2d fields

 do i = 1, n2d
	call ESMF_FieldGet(field_write_2d(i), name=varname, rc=error)
	if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
	  call error_handler("IN FieldGet", error)
	  
	print*,"- CALL FieldGather FOR TARGET GRID ", trim(varname) 
	call ESMF_FieldGather(field_write_2d(i), dum2d, rootPet=0, rc=error)
	if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
	  call error_handler("IN FieldGather", error)


	if (localpet==0) then
		print*, trim(varname), minval(dum2d), maxval(dum2d)
		dum2dt(:,:,1) = dum2d
		error = nf90_put_var( ncid, id_vars2(i), dum2dt)
		call netcdf_err(error, 'WRITING RECORD' )
	endif
 enddo
 deallocate(field_write_2d)
 deallocate(dum2d, dum2dt)
 
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
	  
		print*,"- CALL FieldGather FOR TARGET GRID ", trim(varname) 
		call ESMF_FieldGather(fields(i), dumsoil, rootPet=0, rc=error)
		if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
		  call error_handler("IN FieldGather", error)


		if (localpet==0) then
			print*, trim(varname), minval(dumsoil), maxval(dumsoil)
			dumsoilt(:,:,:,1) = dumsoil
			error = nf90_put_var( ncid, id_vars_soil(i), dumsoilt)
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
			  
	 do i = 1, n_hist_fields_3d_nz  
		call ESMF_FieldGet(fields(i), name=varname, rc=error)
		if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
		  call error_handler("IN FieldGet", error)
	  
		print*,"- CALL FieldGather FOR TARGET GRID ", trim(varname) 
		call ESMF_FieldGather(fields(i), dum3d, rootPet=0, rc=error)
		if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
		  call error_handler("IN FieldGather", error)


		if (localpet==0) then
			print*, trim(varname), minval(dum3d), maxval(dum3d)
			dum3dt(:,:,:,1) = dum3d
			error = nf90_put_var( ncid, id_vars3_nz(i), dum3dt)
			call netcdf_err(error, 'WRITING RECORD' )
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
	  
		print*,"- CALL FieldGather FOR TARGET GRID ", trim(varname) 
		call ESMF_FieldGather(fields(i), dum3dp1, rootPet=0, rc=error)
		if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
		  call error_handler("IN FieldGather", error)
		  
		if (localpet==0) then  
			if(trim(varname)=='Z') then
				do n = 1,i_target
				do j = 1,j_target
				do k = 2,nzp1_input
					dum3dt(n,j,k-1,1) = 0.5 * ( dum3dp1(n,j,k) + dum3dp1(n,j,k-1) )
				enddo
				enddo
				enddo
			
				error = nf90_put_var( ncid, id_z, dum3dt)
				call netcdf_err(error, 'WRITING RECORD' )
			endif
         
			print*, trim(varname), minval(dum3dp1), maxval(dum3dp1)
			dum3dp1t(:,:,:,1) = dum3dp1
			error = nf90_put_var( ncid, id_vars3_nzp1(i), dum3dp1t)
			call netcdf_err(error, 'WRITING RECORD' )

		endif
	 enddo
	 deallocate(fields)
 endif
 

 deallocate(dum3d, dum3dp1, dum3dt, dum3dp1t)
 deallocate(id_vars2, id_vars3_nz, id_vars3_nzp1, id_vars_soil)
 deallocate(target_hist_longname_2d_cons, target_hist_longname_2d_nstd)
 deallocate(target_hist_longname_2d_patch, target_hist_longname_3d_nz)
 deallocate(target_hist_longname_3d_nzp1, target_diag_longname)
 deallocate(target_hist_units_2d_cons, target_hist_units_2d_nstd)
 deallocate(target_hist_units_2d_patch, target_hist_units_3d_nz)
 deallocate(target_hist_units_3d_nzp1, target_diag_units)

 if (localpet == 0) error = nf90_close(ncid)

 end subroutine write_to_file


 end module write_data
