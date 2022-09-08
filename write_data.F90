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

	 use program_setup, only        : data_to_interp
	 
	 implicit none
	 
	 integer, intent(in)            :: localpet
	 
	 if (data_to_interp=='diag') then
	 	call write_diag_data(localpet)
	 elseif (data_to_interp=='init') then
	 	call write_init_data(localpet)
	 endif
	 
 end subroutine write_to_file

!> Write diag data on the target grid
!!
!!
!! @param[in] localpet  ESMF local persistent execution thread
!! @author Larissa Reames CIWRO/NOAA/NSSL

 subroutine write_diag_data(localpet)

 use esmf
 use netcdf

 use model_grid, only              : i_target, j_target, &
                                     ip1_target, jp1_target, &
                                     longitude_target_grid, &
                                     latitude_target_grid, &
                                     target_diag_bundle, &
                                     n_diag_fields

 implicit none

 integer, intent(in)              :: localpet

 character(len=128)               :: outfile
 character(len=50)                :: varname

 integer                          :: error, ncid, n, rc, i
 integer                          :: header_buffer_val = 16384
 integer                          :: dim_lon, dim_lat
 integer                          :: dim_lonp, dim_latp
 integer                          :: id_lat, id_lon
 integer			              :: id_vars(n_diag_fields)

 real(esmf_kind_r8), allocatable  :: data_one_tile(:,:)
  real(esmf_kind_r8), allocatable  :: dum2d(:,:)
 
 type(esmf_field)                 :: fields(n_diag_fields)


 if (localpet ==0) then
   allocate(dum2d(i_target,j_target))
 else
   allocate(dum2d(0,0))
 endif
 
 call ESMF_FieldBundleGet(target_diag_bundle, fieldList=fields, & 
 						  itemorderflag=ESMF_ITEMORDER_ADDORDER, &
 						  rc=error) 
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
		  call error_handler("IN FieldBundleGet", error)

if (localpet == 0) then

!--- open the file
   error = nf90_create(output_file, NF90_NETCDF4, ncid)
   call netcdf_err(error, 'CREATING FILE '//trim(output_file) )

!--- define dimension
   error = nf90_def_dim(ncid, 'west_east', i_target, dim_lon)
   call netcdf_err(error, 'DEFINING LON DIMENSION' )
   error = nf90_def_dim(ncid, 'south_north', j_target, dim_lat)
   call netcdf_err(error, 'DEFINING LAT DIMENSION' )


!--- define fields
   error = nf90_def_var(ncid, 'XLONG', NF90_FLOAT, (/dim_lon,dim_lat/), id_lon)
   call netcdf_err(error, 'DEFINING GEOLON FIELD' )
   error = nf90_put_att(ncid, id_lon, "long_name", "Longitude")
   call netcdf_err(error, 'DEFINING GEOLON NAME' )
   error = nf90_put_att(ncid, id_lon, "units", "degrees_east")
   call netcdf_err(error, 'DEFINING GEOLON UNITS' )

   error = nf90_def_var(ncid, 'XLAT', NF90_FLOAT, (/dim_lon,dim_lat/), id_lat)
   call netcdf_err(error, 'DEFINING GEOLAT FIELD' )
   error = nf90_put_att(ncid, id_lat, "long_name", "Latitude")
   call netcdf_err(error, 'DEFINING GEOLAT NAME' )
   error = nf90_put_att(ncid, id_lat, "units", "degrees_north")
   call netcdf_err(error, 'DEFINING GEOLAT UNITS' )
   
   do i = 1, n_diag_fields  
	call ESMF_FieldGet(fields(i),name=varname,rc=error)
	if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
		  call error_handler("IN FieldGet", error)
	  
	print*,"- DEFINE ON FILE TARGET GRID ", varname
	error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat/), id_vars(i))
	call netcdf_err(error, 'DEFINING VAR' )
	error = nf90_put_att(ncid, id_vars(i), "coordinates", "xlong xlat")
	call netcdf_err(error, 'DEFINING COORD' )	
   enddo
   
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
   error = nf90_put_var( ncid, id_lon, dum2d)
   call netcdf_err(error, 'WRITING LONGITUDE RECORD' )
 endif

!  latitude


   print*,"- CALL FieldGather FOR TARGET GRID LATITUDE"
   call ESMF_FieldGather(latitude_target_grid, dum2d, rootPet=0, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", error)


 if (localpet ==0) then
   error = nf90_put_var( ncid, id_lat, dum2d)
   call netcdf_err(error, 'WRITING LATITUDE RECORD' )
 endif

!   diag fields

 do i = 1, n_diag_fields  
	call ESMF_FieldGet(fields(i), name=varname, rc=error)
	if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
	  call error_handler("IN FieldGet", error)
	  
	print*,"- CALL FieldGather FOR TARGET GRID ", trim(varname) 
	call ESMF_FieldGather(fields(i), dum2d, rootPet=0, rc=error)
	if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
	  call error_handler("IN FieldGather", error)


	if (localpet==0) then
		print*, trim(varname), minval(dum2d), maxval(dum2d)
		error = nf90_put_var( ncid, id_vars(i), dum2d)
		call netcdf_err(error, 'WRITING RECORD' )
	endif
 enddo

 deallocate(dum2d)

if (localpet == 0) error = nf90_close(ncid)

 end subroutine write_diag_data
 
 !> Write init data on the target grid
!!
!!
!! @param[in] localpet  ESMF local persistent execution thread
!! @author Larissa Reames CIWRO/NOAA/NSSL

 subroutine write_init_data(localpet)

 use esmf
 use netcdf

 use model_grid, only              : i_target, j_target, &
                                     ip1_target, jp1_target, &
                                     nz_input, nzp1_input, &
                                     longitude_target_grid, &
                                     latitude_target_grid, &
                                     zgrid_target_grid, &
                                     zgrid_input_grid, &
                                     target_init_bundle_2d, &
                                     target_init_bundle_3d, &
                                     n_init_fields_2d, &
                                     n_init_fields_3d


 implicit none

 integer, intent(in)              :: localpet

 character(len=128)               :: outfile
 character(len=50)                :: varname

 integer                          :: error, ncid, n, rc, i
 integer                          :: header_buffer_val = 16384
 integer                          :: dim_lon, dim_lat, dim_z, dim_zp1
 integer                          :: dim_lonp, dim_latp
 integer                          :: id_lat, id_lon, id_z
 integer			              :: id_vars2(n_init_fields_2d), &
 									 id_vars3(n_init_fields_3d)

 real(esmf_kind_r8), allocatable  :: dum2d(:,:), dum3d(:,:,:)
 
 type(esmf_field)                 :: fields2(n_init_fields_2d), &
 									 fields3(n_init_fields_3d)


 if (localpet ==0) then
   allocate(dum2d(i_target,j_target))
   allocate(dum3d(i_target,j_target,nzp1_input))
 else
   allocate(dum2d(0,0))
   allocate(dum3d(0,0,0))
 endif
 
 call ESMF_FieldBundleGet(target_init_bundle_2d, fieldList=fields2, & 
 						  itemorderflag=ESMF_ITEMORDER_ADDORDER, &
 						  rc=error) 
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
		  call error_handler("IN FieldBundleGet", error)
		  
 call ESMF_FieldBundleGet(target_init_bundle_3d, fieldList=fields3, & 
 						  itemorderflag=ESMF_ITEMORDER_ADDORDER, &
 						  rc=error) 
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
		  call error_handler("IN FieldBundleGet", error)

if (localpet == 0) then

!--- open the file
   error = nf90_create(output_file, NF90_NETCDF4, ncid)
   call netcdf_err(error, 'CREATING FILE '//trim(output_file) )

!--- define dimension
   error = nf90_def_dim(ncid, 'west_east', i_target, dim_lon)
   call netcdf_err(error, 'DEFINING LON DIMENSION' )
   error = nf90_def_dim(ncid, 'south_north', j_target, dim_lat)
   call netcdf_err(error, 'DEFINING LAT DIMENSION' )
   error = nf90_def_dim(ncid, 'bottom_top', nz_input, dim_z)
   call netcdf_err(error, 'DEFINING VERTICAL DIMENSION' )
   error = nf90_def_dim(ncid, 'bottom_top_stag', nzp1_input, dim_zp1)
   call netcdf_err(error, 'DEFINING VERTICALP1 DIMENSION' )


!--- define fields
   error = nf90_def_var(ncid, 'XLONG', NF90_FLOAT, (/dim_lon,dim_lat/), id_lon)
   call netcdf_err(error, 'DEFINING GEOLON FIELD' )
   error = nf90_put_att(ncid, id_lon, "long_name", "Longitude")
   call netcdf_err(error, 'DEFINING GEOLON NAME' )
   error = nf90_put_att(ncid, id_lon, "units", "degrees_east")
   call netcdf_err(error, 'DEFINING GEOLON UNITS' )

   error = nf90_def_var(ncid, 'XLAT', NF90_FLOAT, (/dim_lon,dim_lat/), id_lat)
   call netcdf_err(error, 'DEFINING GEOLAT FIELD' )
   error = nf90_put_att(ncid, id_lat, "long_name", "Latitude")
   call netcdf_err(error, 'DEFINING GEOLAT NAME' )
   error = nf90_put_att(ncid, id_lat, "units", "degrees_north")
   call netcdf_err(error, 'DEFINING GEOLAT UNITS' )
   
   error = nf90_def_var(ncid, 'Z', NF90_FLOAT, (/dim_lon,dim_lat,dim_zp1/), id_z)
   call netcdf_err(error, 'DEFINING ZGRID FIELD' )
   error = nf90_put_att(ncid, id_z, "long_name", "Height above mean sea level")
   call netcdf_err(error, 'DEFINING ZGRID NAME' )
   error = nf90_put_att(ncid, id_z, "units", "m AMSL")
   call netcdf_err(error, 'DEFINING zgrid UNITS' )
   
   do i = 1, n_init_fields_2d  
	call ESMF_FieldGet(fields2(i),name=varname,rc=error)
	if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
		  call error_handler("IN FieldGet", error)
	  
	print*,"- DEFINE ON FILE TARGET GRID ", varname
	error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat/), id_vars2(i))
	call netcdf_err(error, 'DEFINING VAR' )
	error = nf90_put_att(ncid, id_vars2(i), "coordinates", "xlong xlat")
	call netcdf_err(error, 'DEFINING COORD' )	
   enddo
   
   do i = 1, n_init_fields_3d  
	call ESMF_FieldGet(fields3(i),name=varname,rc=error)
	if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
		  call error_handler("IN FieldGet", error)
	  
	print*,"- DEFINE ON FILE TARGET GRID ", varname
	error = nf90_def_var(ncid, varname, NF90_FLOAT, (/dim_lon,dim_lat,dim_z/), id_vars3(i))
	call netcdf_err(error, 'DEFINING VAR' )
	error = nf90_put_att(ncid, id_vars3(i), "coordinates", "xlong xlat z")
	call netcdf_err(error, 'DEFINING COORD' )	
   enddo
   
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
   error = nf90_put_var( ncid, id_lon, dum2d)
   call netcdf_err(error, 'WRITING LONGITUDE RECORD' )
 endif

!  latitude


   print*,"- CALL FieldGather FOR TARGET GRID LATITUDE"
   call ESMF_FieldGather(latitude_target_grid, dum2d, rootPet=0, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", error)


 if (localpet ==0) then
   error = nf90_put_var( ncid, id_lat, dum2d)
   call netcdf_err(error, 'WRITING LATITUDE RECORD' )
 endif
 
!  z
   print*,"- CALL FieldGather FOR TARGET GRID ZGRID"
   call ESMF_FieldGather(zgrid_target_grid, dum3d, rootPet=0, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGather", error)


 if (localpet ==0) then
   error = nf90_put_var( ncid, id_z, dum3d)
   call netcdf_err(error, 'WRITING ZGRID RECORD' )
 endif

!   2d init fields

 do i = 1, n_init_fields_2d  
	call ESMF_FieldGet(fields2(i), name=varname, rc=error)
	if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
	  call error_handler("IN FieldGet", error)
	  
	print*,"- CALL FieldGather FOR TARGET GRID ", trim(varname) 
	call ESMF_FieldGather(fields2(i), dum2d, rootPet=0, rc=error)
	if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
	  call error_handler("IN FieldGather", error)


	if (localpet==0) then
		print*, trim(varname), minval(dum2d), maxval(dum2d)
		error = nf90_put_var( ncid, id_vars2(i), dum2d)
		call netcdf_err(error, 'WRITING RECORD' )
	endif
 enddo

 deallocate(dum2d)
 
 ! 3d init fields
 
 do i = 1, n_init_fields_3d  
	call ESMF_FieldGet(fields3(i), name=varname, rc=error)
	if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
	  call error_handler("IN FieldGet", error)
	  
	print*,"- CALL FieldGather FOR TARGET GRID ", trim(varname) 
	call ESMF_FieldGather(fields3(i), dum3d(:,:,1:nz_input), rootPet=0, rc=error)
	if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
	  call error_handler("IN FieldGather", error)


	if (localpet==0) then
		print*, trim(varname), minval(dum3d), maxval(dum3d)
		error = nf90_put_var( ncid, id_vars3(i), dum3d(:,:,1:nz_input))
		call netcdf_err(error, 'WRITING RECORD' )
	endif
 enddo

 deallocate(dum3d)

if (localpet == 0) error = nf90_close(ncid)

 end subroutine write_init_data


 end module write_data
