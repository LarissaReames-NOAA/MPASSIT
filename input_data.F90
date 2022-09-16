!> @file
!! @brief Read atmospheric and surface data from MPAS diag, hist, or forecast files.
!! @author  Larissa Reames CIWRO/NOAA/NSSL 

!> Read atmospheric and/or surface input MPAS grid.
!! Supported formats include INITIalization file, forecast diag file (only 2d diagnostic
!! fields), and full 3d forecast files.
!!
!! Public variables are defined below: "input" indicates field
!! associated with the input grid.
!!
!! @author  Larissa Reames CIWRO/NOAA/NSSL 

 module input_data

 use esmf
 use netcdf

 use program_setup, only          : hist_file_input_grid, &
                                    diag_file_input_grid, &
                                    grid_file_input_grid, &
                                    interp_diag, interp_hist

 use model_grid, only             : input_grid,        &
                                    nCells_input, nVert_input,  &
                                    nz_input, nzp1_input, &
                                    nsoil_input, strlen, &
                                    valid_time, config_dt, &
                                    start_time, lsm_scheme, &
                                    mp_scheme, conv_scheme, &
                                    cell_latitude_input_grid, &
                                    cell_longitude_input_grid, &
                                    zgrid_input_grid, &
                                    input_diag_bundle, &
                                    target_diag_names, &
                                    input_hist_bundle_2d_cons, &
                                    input_hist_bundle_2d_patch, &
                                    input_hist_bundle_2d_nstd, &
                                    input_hist_bundle_3d_nz, &
                                    input_hist_bundle_3d_nzp1, &
                                    input_hist_bundle_soil, &
                                    target_hist_names_2d_cons, &
                                    target_hist_names_2d_nstd, &
                                    target_hist_names_2d_patch, &
                                    target_hist_names_3d_nz, &
                                    target_hist_names_3d_nzp1, &
                                    target_hist_names_soil, &
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
                                    target_hist_longname_soil, &
                                    n_diag_fields, &
                                    n_hist_fields_2d_cons, &
                                    n_hist_fields_2d_nstd, &
                                    n_hist_fields_2d_patch, &
                                    n_hist_fields_3d_nz, &
                                    n_hist_fields_3d_nzp1, &
                                    n_hist_fields_soil, &
                                    elemIDs, nCellsPerPET, &
                                    nodeIDs

 implicit none

 private
 public :: read_input_data
                                         
 contains

!> Read input grid data driver.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author Larissa Reames CIWRO/NOAA/NSSL   
 subroutine read_input_data(localpet)

 implicit none

 include 'mpif.h'

 integer, intent(in)             :: localpet
 
 if (interp_diag) then
 	call read_input_diag_data(localpet)
 endif
 
 if (interp_hist) then
 	call read_input_hist_data(localpet)
 endif
 
 if (.not. interp_diag .and. .not. interp_hist) then
 	call error_handler(" SET INTERP_DIAG AND/OR INTERP_HIST TO TRUE TO OBTAIN OUTPUT", -1)
 endif
 
 end subroutine read_input_data
 
!> Read input grid diag data.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author Larissa Reames CIWRO/NOAA/NSSL
 subroutine read_input_diag_data(localpet)

 character(len=500)              :: the_file
 character(len=50)               :: vname
 
 integer, intent(in)             :: localpet
 integer                         :: error, ncid, rc
 integer                         :: id_dim
 integer                         :: id_var, i, j, nodes
 
 type(esmf_field),allocatable    :: fields(:)
 
 real(esmf_kind_r8), allocatable :: dummy(:)

 real(esmf_kind_r8), pointer     :: varptr(:)
 
 call init_input_diag_fields()

 print*,"- READ INPUT DIAG DATA."


 the_file = trim(diag_file_input_grid)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening: '//trim(the_file) )


!---------------------------------------------------------------------------
! Initialize esmf atmospheric fields.
!---------------------------------------------------------------------------
 allocate(fields(n_diag_fields))
 allocate(target_diag_units(n_diag_fields))
 allocate(target_diag_longname(n_diag_fields))
 call ESMF_FieldBundleGet(input_diag_bundle, fieldList=fields, & 
 						  itemorderflag=ESMF_ITEMORDER_ADDORDER, &
 						  rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
	 call error_handler("IN FieldBundleGet", rc)
	 
 call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
	 call error_handler("IN MeshGet", rc)
 
 allocate(dummy(nCells_input))
 
 do i = 1,n_diag_fields

 	call ESMF_FieldGet(fields(i), name=vname, rc=rc)
 	if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
	 call error_handler("IN FieldGet", rc)
 	
 	call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc)
 	if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
	 call error_handler("IN FieldGet", rc)
 	
	print*,"- READ ", trim(vname)
	error=nf90_inq_varid(ncid, trim(vname), id_var)
	call netcdf_err(error, 'reading field id' )
	error=nf90_get_var(ncid, id_var, dummy)
	call netcdf_err(error, 'reading field' )
	error=nf90_get_att(ncid,id_var,'units',target_diag_units(i))
	call netcdf_err(error, 'reading field units' )
	error=nf90_get_att(ncid,id_var,'long_name',target_diag_longname(i))
	call netcdf_err(error, 'reading field long_name' )
		
	print*,"- SET ON MESH ", trim(vname)
	do j = 1, nCellsPerPET
		varptr(j) = dummy(elemIDs(j))
	enddo
	

	
	print*, localpet, minval(varptr), maxval(varptr)	
	nullify(varptr)
 enddo

 error = nf90_close(ncid)


 deallocate( dummy)

 end subroutine read_input_diag_data
 
 subroutine init_input_diag_fields
 
 	implicit none
 	
 	integer                       :: i, rc
 	type(esmf_field),allocatable  :: diag_fields(:)
 	character(50), allocatable 	  :: input_diag_names(:)
 	character(50)                 :: fname
 	
 	fname = 'diaglist'
 	
 	call read_varlist(fname,n_diag_fields,input_diag_names, target_diag_names)
 	
 	allocate(diag_fields(n_diag_fields))
 	print*,"- INITIALIZE INPUT DIAG FIELDS."
 	do i = 1, n_diag_fields
 	
 		print*, "- INIT FIELD ", input_diag_names(i)	
		diag_fields(i) = ESMF_FieldCreate(input_grid, & 
							typekind=ESMF_TYPEKIND_R8, &
							meshloc=ESMF_MESHLOC_ELEMENT, &
							name=input_diag_names(i), rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
		call error_handler("IN FieldCreate", rc)
    enddo
    
    
 	input_diag_bundle = ESMF_FieldBundleCreate(fieldList=diag_fields, & 
 									name="input diag data", rc=rc)
 	if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldBundleCreate", rc)
 
 	deallocate(diag_fields, input_diag_names)
 end subroutine init_input_diag_fields
 
!> Read input grid hist data.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author Larissa Reames CIWRO/NOAA/NSSL
 subroutine read_input_hist_data(localpet)

 character(len=500)              :: the_file
 character(len=50)               :: vname, att_text
 
 integer, intent(in)             :: localpet
 integer                         :: error, ncid, rc
 integer                         :: id_dim
 integer                         :: id_var, i, j, nodes

 type(esmf_field),allocatable    :: fields(:)
 
 real(esmf_kind_r8), allocatable :: dummy2(:,:), dummy3(:,:,:)

 real(esmf_kind_r8), pointer     :: varptr(:), varptr2(:,:)

 
 call init_input_hist_fields()

 print*,"- READ INPUT HIST DATA."


 the_file = trim(hist_file_input_grid)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening: '//trim(the_file) )
 
!---------------------------------------------------------------------------
! Read global attributes for use when creating output file
!---------------------------------------------------------------------------
 
 print*,'- READ GLOBAL ATTRIBUTE LSM SCHEME'
 error = nf90_get_att(ncid,NF90_GLOBAL,'config_lsm_scheme',att_text)
 call netcdf_err(error, 'reading config_lsm_scheme')
 if (trim(att_text) == 'noah') then
 	lsm_scheme = 2
 elseif (trim(att_text) == 'ruc') then
 	lsm_scheme = 3
 endif
 
 print*,'- READ GLOBAL ATTRIBUTE START TIME'
 error = nf90_get_att(ncid,NF90_GLOBAL,'config_start_time',start_time)
 call netcdf_err(error, 'reading config_start_time')
 
 print*,'- READ GLOBAL ATTRIBUTE LMP SCHEME'
 error = nf90_get_att(ncid,NF90_GLOBAL,'config_microp_scheme',att_text)
 call netcdf_err(error, 'reading config_microp_scheme')
  if (trim(att_text) == 'mp_thompson') then
 	mp_scheme = 8
 elseif (trim(att_text) == 'nssl') then
 	mp_scheme = 17
 endif
 
 
 print*,'- READ GLOBAL ATTRIBUTE CONVECTION SCHEME'
 error = nf90_get_att(ncid,NF90_GLOBAL,'config_convection_scheme',att_text)
 call netcdf_err(error, 'reading config_conv_scheme')
 
  if (trim(att_text) == 'cu_ntiedke') then
 	conv_scheme = 16
 elseif (trim(att_text) == 'cu_kainfritsch') then
 	conv_scheme = 1
 elseif (trim(att_text) == 'cu_grellfreitas') then
 	conv_scheme = 3
 endif

 print*,'- READ GLOBAL ATTRIBUTE CONFIG_DT'
 error = nf90_get_att(ncid,NF90_GLOBAL,'config_dt',config_dt)
 call netcdf_err(error, 'reading config_dt')

 vname = 'xtime'
 print*, '- READ TIMES VARIABLE'
 error = nf90_inq_dimid(ncid,'StrLen',id_var)
 call netcdf_err(error, 'reading strlen dim id')
 error = nf90_inquire_dimension(ncid,id_var,len=strlen)
 allocate(valid_time(1,strlen))
 error = nf90_inq_varid(ncid,vname,id_var)
 call netcdf_err(error, 'reading xtime id')
 error = nf90_get_var(ncid, id_var, valid_time)
 call netcdf_err(error, 'getting xtime')

!---------------------------------------------------------------------------
! Initialize 2d esmf atmospheric fields for bilinear/patch interpolation
!---------------------------------------------------------------------------

 if (n_hist_fields_2d_patch > 0) then
 	allocate(fields(n_hist_fields_2d_patch))
 	allocate(target_hist_units_2d_patch(n_hist_fields_2d_patch))
    allocate(target_hist_longname_2d_patch(n_hist_fields_2d_patch))
	 call ESMF_FieldBundleGet(input_hist_bundle_2d_patch, fieldList=fields, & 
							  itemorderflag=ESMF_ITEMORDER_ADDORDER, &
							  rc=rc)
	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldBundleGet", rc)
	 
	 call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN MeshGet", rc)
 
	 allocate(dummy2(nCells_input,1))
 
	 do i = 1,n_hist_fields_2d_patch

		call ESMF_FieldGet(fields(i), name=vname, rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldGet", rc)
	
		call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldGet", rc)
	
		print*,"- READ ", trim(vname)
		error=nf90_inq_varid(ncid, trim(vname), id_var)
		call netcdf_err(error, 'reading field id' )
		error=nf90_get_var(ncid, id_var, dummy2)
		call netcdf_err(error, 'reading field' )
		error=nf90_get_att(ncid,id_var,'units',target_hist_units_2d_patch(i))
		call netcdf_err(error, 'reading field units' )
		error=nf90_get_att(ncid,id_var,'long_name',target_hist_longname_2d_patch(i))
		call netcdf_err(error, 'reading field long_name' )
		
		print*,"- SET ON MESH ", trim(vname)
		do j = 1, nCellsPerPET
			varptr(j) = dummy2(elemIDs(j),1)
		enddo

		print*, localpet, minval(varptr), maxval(varptr)	
		nullify(varptr)
	 enddo
	 deallocate(dummy2)
	 deallocate(fields)
 endif
 
!---------------------------------------------------------------------------
! Initialize 2d esmf atmospheric fields for conservative interpolation
!---------------------------------------------------------------------------

 if (n_hist_fields_2d_cons > 0) then
 	allocate(fields(n_hist_fields_2d_cons))
 	allocate(target_hist_units_2d_cons(n_hist_fields_2d_cons))
    allocate(target_hist_longname_2d_cons(n_hist_fields_2d_cons))
	 call ESMF_FieldBundleGet(input_hist_bundle_2d_cons, fieldList=fields, & 
							  itemorderflag=ESMF_ITEMORDER_ADDORDER, &
							  rc=rc)
	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldBundleGet", rc)
	 
	 call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN MeshGet", rc)
 
	 allocate(dummy2(nCells_input,1))
 
	 do i = 1,n_hist_fields_2d_cons

		call ESMF_FieldGet(fields(i), name=vname, rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldGet", rc)
	
		call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldGet", rc)
	
		print*,"- READ ", trim(vname)
		error=nf90_inq_varid(ncid, trim(vname), id_var)
		call netcdf_err(error, 'reading field id' )
		error=nf90_get_var(ncid, id_var, dummy2)
		call netcdf_err(error, 'reading field' )
		error=nf90_get_att(ncid,id_var,'units',target_hist_units_2d_cons(i))
		call netcdf_err(error, 'reading field units' )
		error=nf90_get_att(ncid,id_var,'long_name',target_hist_longname_2d_cons(i))
		call netcdf_err(error, 'reading field long_name' )
		
		print*,"- SET ON MESH ", trim(vname)
		do j = 1, nCellsPerPET
			varptr(j) = dummy2(elemIDs(j),1)
		enddo
	

	
		print*, vname, minval(varptr), maxval(varptr)	
		nullify(varptr)
	 enddo
	 deallocate(dummy2)
	 deallocate(fields)
 endif
 
!---------------------------------------------------------------------------------
! Initialize 2d esmf atmospheric fields for nearest source to dest interpolation
!---------------------------------------------------------------------------------

 if (n_hist_fields_2d_nstd > 0) then
 	allocate(fields(n_hist_fields_2d_nstd))
 	allocate(target_hist_units_2d_nstd(n_hist_fields_2d_nstd))
    allocate(target_hist_longname_2d_nstd(n_hist_fields_2d_nstd))
	 call ESMF_FieldBundleGet(input_hist_bundle_2d_nstd, fieldList=fields, & 
							  itemorderflag=ESMF_ITEMORDER_ADDORDER, &
							  rc=rc)
	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldBundleGet", rc)
	 
	 call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN MeshGet", rc)
 
	 allocate(dummy2(nCells_input,1))
 
	 do i = 1,n_hist_fields_2d_nstd

		call ESMF_FieldGet(fields(i), name=vname, rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldGet", rc)
	
		call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldGet", rc)
	
		print*,"- READ ", trim(vname)
		error=nf90_inq_varid(ncid, trim(vname), id_var)
		call netcdf_err(error, 'reading field id' )
		error=nf90_get_var(ncid, id_var, dummy2)
		call netcdf_err(error, 'reading field' )
		error=nf90_get_att(ncid,id_var,'units',target_hist_units_2d_nstd(i))
		call netcdf_err(error, 'reading field units' )
		error=nf90_get_att(ncid,id_var,'long_name',target_hist_longname_2d_nstd(i))
		call netcdf_err(error, 'reading field long_name' )
		
		print*,"- SET ON MESH ", trim(vname)
		do j = 1, nCellsPerPET
			varptr(j) = dummy2(elemIDs(j),1)
		enddo
	

	
		print*, vname, minval(varptr), maxval(varptr)	
		nullify(varptr)
	 enddo
	 deallocate(dummy2)
	 deallocate(fields)
 endif
 
!---------------------------------------------------------------------------
! Initialize 3d esmf soil fields for bilinear/patch interpolation
!---------------------------------------------------------------------------

 if (n_hist_fields_soil > 0) then
 	allocate(fields(n_hist_fields_soil))
 	allocate(target_hist_units_soil(n_hist_fields_soil))
    allocate(target_hist_longname_soil(n_hist_fields_soil))
	 call ESMF_FieldBundleGet(input_hist_bundle_soil, fieldList=fields, & 
							  itemorderflag=ESMF_ITEMORDER_ADDORDER, &
							  rc=rc)
	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldBundleGet", rc)
	 
	 call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN MeshGet", rc)
 
	 allocate(dummy3(nsoil_input,nCells_input,1))
 
	 do i = 1,n_hist_fields_soil

		call ESMF_FieldGet(fields(i), name=vname, rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldGet", rc)
	
		call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldGet", rc)
	
		print*,"- READ ", trim(vname)
		error=nf90_inq_varid(ncid, trim(vname), id_var)
		call netcdf_err(error, 'reading field id' )
		error=nf90_get_var(ncid, id_var, dummy3)
		call netcdf_err(error, 'reading field' )
		error=nf90_get_att(ncid,id_var,'units',target_hist_units_soil(i))
		call netcdf_err(error, 'reading field units' )
		error=nf90_get_att(ncid,id_var,'long_name',target_hist_longname_soil(i))
		call netcdf_err(error, 'reading field long_name' )
		
		print*,"- SET ON MESH ", trim(vname)
		do j = 1, nCellsPerPET
			varptr2(j,:) = dummy3(:,elemIDs(j),1)
		enddo

		print*, vname, minval(varptr2), maxval(varptr2)
		nullify(varptr2)
	 enddo
	 deallocate(dummy3)
	 deallocate(fields)
 endif

 
!---------------------------------------------------------------------------
! Initialize 3d esmf atmospheric fields with nVertLevels vertical dimension
!---------------------------------------------------------------------------

 if (n_hist_fields_3d_nz > 0 ) then
 	 allocate(fields(n_hist_fields_3d_nz))
 	 allocate(target_hist_units_3d_nz(n_hist_fields_3d_nz))
    allocate(target_hist_longname_3d_nz(n_hist_fields_3d_nz))
	 call ESMF_FieldBundleGet(input_hist_bundle_3d_nz, fieldList=fields, & 
							  itemorderflag=ESMF_ITEMORDER_ADDORDER, &
							  rc=rc)
	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldBundleGet", rc)
	 
	 call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN MeshGet", rc)
 
	 allocate(dummy3(nz_input,nCells_input,1))
 
	 do i = 1,n_hist_fields_3d_nz

		call ESMF_FieldGet(fields(i), name=vname, rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldGet", rc)
	
		call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldGet", rc)
	
		print*,"- READ ", trim(vname)
		error=nf90_inq_varid(ncid, trim(vname), id_var)
		call netcdf_err(error, 'reading field id' )
		error=nf90_get_var(ncid, id_var, dummy3)
		call netcdf_err(error, 'reading field' )
		error=nf90_get_att(ncid,id_var,'units',target_hist_units_3d_nz(i))
		call netcdf_err(error, 'reading field units' )
		error=nf90_get_att(ncid,id_var,'long_name',target_hist_longname_3d_nz(i))
		call netcdf_err(error, 'reading field long_name' )
		
		print*,"- SET ON MESH ", trim(vname)
		do j = 1, nCellsPerPET
			varptr2(j,:) = dummy3(:,elemIDs(j),1)
		enddo
	
		print*, vname, minval(varptr2), maxval(varptr2)
		nullify(varptr2)
	 enddo
	 deallocate(dummy3)
	 deallocate(fields)
 endif

!-------------------------------------------------------------------------------
! Initialize 3d esmf atmospheric fields with nVertLevels+1 vertical dimension
!-------------------------------------------------------------------------------
 
 if (n_hist_fields_3d_nzp1 > 0 ) then
 	 allocate(fields(n_hist_fields_3d_nzp1))
 	 allocate(target_hist_units_3d_nzp1(n_hist_fields_3d_nzp1))
     allocate(target_hist_longname_3d_nzp1(n_hist_fields_3d_nzp1))
	 call ESMF_FieldBundleGet(input_hist_bundle_3d_nzp1, fieldList=fields, & 
							  itemorderflag=ESMF_ITEMORDER_ADDORDER, &
							  rc=rc)
	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldBundleGet", rc)
	 
	 call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN MeshGet", rc)
 
	 allocate(dummy3(nzp1_input,nCells_input,1))
 
	 do i = 1,n_hist_fields_3d_nzp1

		call ESMF_FieldGet(fields(i), name=vname, rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldGet", rc)
	
		call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
		 call error_handler("IN FieldGet", rc)
	
		print*,"- READ ", trim(vname)
		error=nf90_inq_varid(ncid, trim(vname), id_var)
		call netcdf_err(error, 'reading field id' )
		error=nf90_get_var(ncid, id_var, dummy3)
		call netcdf_err(error, 'reading field' )
		error=nf90_get_att(ncid,id_var,'units',target_hist_units_3d_nzp1(i))
		call netcdf_err(error, 'reading field units' )
		error=nf90_get_att(ncid,id_var,'long_name',target_hist_longname_3d_nzp1(i))
		call netcdf_err(error, 'reading field long_name' )
		
		print*,"- SET ON MESH ", trim(vname)
		do j = 1, nCellsPerPET
			varptr2(j,:) = dummy3(:,elemIDs(j),1)
		enddo
		
		print*, vname, minval(varptr2), maxval(varptr2)
		nullify(varptr2)
	 enddo
	 deallocate(dummy3)
	 deallocate(fields)
 endif
 

 error = nf90_close(ncid)


 end subroutine read_input_hist_data
 
 subroutine init_input_hist_fields
 
 	implicit none
 	
 	integer                       :: i, j, k, n, rc 
 	integer                       :: n_hist_fields_2d, n_hist_fields_3d
 	!type(esmf_field), allocatable :: hist_fields_2d_cons(:), &
 !								   hist_fields_2d_nstd(:), &
 !								   hist_fields_2d_patch(:)
 	character(50), allocatable    :: input_hist_names_2d(:), &
 									 input_hist_names_2d_cons(:), &
 									 input_hist_names_2d_nstd(:), &
 									 input_hist_names_2d_patch(:), &
 									 input_hist_names_3d(:), &
 								     input_hist_names_3d_nz(:), &
 								     input_hist_names_soil(:), &
 								     input_hist_names_3d_nzp1(:), &
 								     target_hist_names_2d(:), &
 								     target_hist_names_3d(:)
 	type(esmf_field), allocatable :: fields(:)
 	character(50)                 :: cons_vars(2), nstd_vars(4), nzp1_vars(2)
 	character(50)                 :: fname
 	
 	cons_vars = (/'snow','snowh'/)
 	nstd_vars = (/'ivgtyp','isltyp', 'xland','landmask'/)
 	nzp1_vars = (/'zgrid','w'/)
 	n_hist_fields_2d_cons = 0
 	n_hist_fields_2d_nstd = 0
 	n_hist_fields_2d_patch = 0
 	n_hist_fields_3d_nz = 0
 	n_hist_fields_3d_nzp1 = 0
 	
 	fname = 'histlist_2d'
 	call read_varlist(fname,n_hist_fields_2d,input_hist_names_2d, target_hist_names_2d)
 	fname = 'histlist_3d'
 	call read_varlist(fname,n_hist_fields_3d,input_hist_names_3d, target_hist_names_3d)
 	fname = 'histlist_soil'
 	call read_varlist(fname,n_hist_fields_soil,input_hist_names_soil, target_hist_names_soil)
 	
 	do i = 1, n_hist_fields_2d
 		if (any(cons_vars == input_hist_names_2d(i))) then
 			n_hist_fields_2d_cons = n_hist_fields_2d_cons + 1 
 		elseif (any(nstd_vars == input_hist_names_2d(i))) then
 			n_hist_fields_2d_nstd = n_hist_fields_2d_nstd + 1 
 		else
 		    n_hist_fields_2d_patch = n_hist_fields_2d_patch + 1
 		endif
 	enddo
 	
 	allocate(input_hist_names_2d_cons(n_hist_fields_2d_cons))
 	allocate(input_hist_names_2d_nstd(n_hist_fields_2d_nstd))
 	allocate(input_hist_names_2d_patch(n_hist_fields_2d_patch))
 	allocate(target_hist_names_2d_cons(n_hist_fields_2d_cons))
 	allocate(target_hist_names_2d_nstd(n_hist_fields_2d_nstd))
 	allocate(target_hist_names_2d_patch(n_hist_fields_2d_patch))

 	
 	j = 0
 	k = 0
 	n = 0
 	do i = 1, n_hist_fields_2d
 		if (any(cons_vars == input_hist_names_2d(i))) then
 			j = j+1
 			input_hist_names_2d_cons(j) = input_hist_names_2d(i)
 			target_hist_names_2d_cons(j) = target_hist_names_2d(i)
 			
 		elseif (any(nstd_vars == input_hist_names_2d(i))) then
 			k = k+1
 			input_hist_names_2d_nstd(k) = input_hist_names_2d(i)
 			target_hist_names_2d_nstd(k) = target_hist_names_2d(i)
 		else 
 			n = n+1
 			input_hist_names_2d_patch(n) = input_hist_names_2d(i)
 			target_hist_names_2d_patch(n) = target_hist_names_2d(i)
 		endif
 	enddo
 	
 	do i = 1, n_hist_fields_3d
 		if (any(nzp1_vars == input_hist_names_3d(i))) then
 			n_hist_fields_3d_nzp1 = n_hist_fields_3d_nzp1 + 1
 		else
 		    n_hist_fields_3d_nz = n_hist_fields_3d_nz + 1
 		endif
 	enddo
 	
 	allocate(input_hist_names_3d_nz(n_hist_fields_3d_nz))
 	allocate(input_hist_names_3d_nzp1(n_hist_fields_3d_nzp1))
 	allocate(target_hist_names_3d_nz(n_hist_fields_3d_nz))
 	allocate(target_hist_names_3d_nzp1(n_hist_fields_3d_nzp1))
	
 	j = 0
 	k = 0
 	do i = 1, n_hist_fields_3d
 		if (any(nzp1_vars == input_hist_names_3d(i))) then
 			j = j+1
 			input_hist_names_3d_nzp1(j) = input_hist_names_3d(i)
 			target_hist_names_3d_nzp1(j) = target_hist_names_3d(i)
 			
 		else
 			k = k+1 
 			input_hist_names_3d_nz(k) = input_hist_names_3d(i)
 			target_hist_names_3d_nz(k) = target_hist_names_3d(i)

 		endif
 	enddo

 	
 	print*,"- INITIALIZE INPUT HIST FIELDS."
 	if (n_hist_fields_2d_cons > 0) then	
 		allocate(fields(n_hist_fields_2d_cons))
		do i = 1, n_hist_fields_2d_cons
	
			print*, "- INIT FIELD ", input_hist_names_2d_cons(i)

			fields(i) = ESMF_FieldCreate(input_grid, & 
								typekind=ESMF_TYPEKIND_R8, &
								meshloc=ESMF_MESHLOC_ELEMENT, &
								name=input_hist_names_2d_cons(i), rc=rc)
			if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
			call error_handler("IN FieldCreate", rc)
		enddo
	
	
		input_hist_bundle_2d_cons = ESMF_FieldBundleCreate(fieldList=fields, & 
										name="input hist 2d data cons", rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
			call error_handler("IN FieldBundleCreate", rc)
		deallocate(fields)
    endif
    
    if (n_hist_fields_2d_nstd > 0) then	
    	allocate(fields(n_hist_fields_2d_nstd))
		do i = 1, n_hist_fields_2d_nstd
	
			print*, "- INIT FIELD ", input_hist_names_2d_nstd(i)

			fields(i) = ESMF_FieldCreate(input_grid, & 
								typekind=ESMF_TYPEKIND_R8, &
								meshloc=ESMF_MESHLOC_ELEMENT, &
								name=input_hist_names_2d_nstd(i), rc=rc)
			if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
			call error_handler("IN FieldCreate", rc)
		enddo
	
	
		input_hist_bundle_2d_nstd = ESMF_FieldBundleCreate(fieldList=fields, & 
										name="input hist 2d data nstd", rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
			call error_handler("IN FieldBundleCreate", rc)
			
		deallocate(fields)
	endif
    	
    if (n_hist_fields_2d_patch > 0) then	
    	allocate(fields(n_hist_fields_2d_patch))
		do i = 1, n_hist_fields_2d_patch
	
			print*, "- INIT FIELD ", input_hist_names_2d_patch(i)

			fields(i) = ESMF_FieldCreate(input_grid, & 
								typekind=ESMF_TYPEKIND_R8, &
								meshloc=ESMF_MESHLOC_ELEMENT, &
								name=input_hist_names_2d_patch(i), rc=rc)
			if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
			call error_handler("IN FieldCreate", rc)
		enddo
	
	
		input_hist_bundle_2d_patch = ESMF_FieldBundleCreate(fieldList=fields, & 
										name="input hist 2d data patch", rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
			call error_handler("IN FieldBundleCreate", rc)
			
		deallocate(fields)
	endif
	
	if (n_hist_fields_soil > 0) then
		allocate(fields(n_hist_fields_soil))
		do i = 1, n_hist_fields_soil
	
			print*, "- INIT FIELD ", input_hist_names_soil(i)

			fields(i) = ESMF_FieldCreate(input_grid, & 
								typekind=ESMF_TYPEKIND_R8, &
								meshloc=ESMF_MESHLOC_ELEMENT, &
								name=input_hist_names_soil(i), & 
								ungriddedLBound=(/1/), &
								ungriddedUBound=(/nsoil_input/), rc=rc)
			if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
			call error_handler("IN FieldCreate", rc)
		enddo
	
		input_hist_bundle_soil = ESMF_FieldBundleCreate(fieldList=fields, & 
										name="input hist soil data", rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
			call error_handler("IN FieldBundleCreate", rc)
			
		deallocate(fields)
	endif

	if (n_hist_fields_3d_nz > 0) then
		allocate(fields(n_hist_fields_3d_nz))
		do i = 1, n_hist_fields_3d_nz
	
			print*, "- INIT FIELD ", input_hist_names_3d_nz(i)

			fields(i) = ESMF_FieldCreate(input_grid, & 
								typekind=ESMF_TYPEKIND_R8, &
								meshloc=ESMF_MESHLOC_ELEMENT, &
								name=input_hist_names_3d_nz(i), & 
								ungriddedLBound=(/1/), &
								ungriddedUBound=(/nz_input/), rc=rc)
			if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
			call error_handler("IN FieldCreate", rc)
		enddo
	
		input_hist_bundle_3d_nz = ESMF_FieldBundleCreate(fieldList=fields, & 
										name="input hist 3d nz data", rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
			call error_handler("IN FieldBundleCreate", rc)
			
		deallocate(fields)
	endif
    	
    if (n_hist_fields_3d_nzp1 > 0) then
    	allocate(fields(n_hist_fields_3d_nzp1))
		do i = 1, n_hist_fields_3d_nzp1
	
			print*, "- INIT FIELD ", input_hist_names_3d_nzp1(i)

			fields(i) = ESMF_FieldCreate(input_grid, & 
								typekind=ESMF_TYPEKIND_R8, &
								meshloc=ESMF_MESHLOC_ELEMENT, &
								name=input_hist_names_3d_nzp1(i), & 
								ungriddedLBound=(/1/), &
								ungriddedUBound=(/nzp1_input/), rc=rc)
			if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
			call error_handler("IN FieldCreate", rc)
		enddo
	
		input_hist_bundle_3d_nzp1 = ESMF_FieldBundleCreate(fieldList=fields, & 
										name="input hist 3d nzp1 data", rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
			call error_handler("IN FieldBundleCreate", rc)
			
		deallocate(fields)
	endif
	
 	deallocate(input_hist_names_2d, target_hist_names_3d)
 	deallocate(input_hist_names_2d_cons)
 	deallocate(input_hist_names_2d_nstd)
 	deallocate(input_hist_names_2d_patch)
 	deallocate(input_hist_names_3d_nz)
 	deallocate(input_hist_names_3d_nzp1)
 	deallocate(input_hist_names_soil)
 	
 
 end subroutine init_input_hist_fields

subroutine read_varlist(file,nfields,field_names,field_names_target)

	implicit none
	
   character(50), INTENT(IN)				:: file
   integer, INTENT(OUT)      				:: nfields
   character(50), INTENT(OUT), ALLOCATABLE  :: field_names(:), field_names_target(:)
   
   integer :: k, istat
   character(200) :: line
	
   open(14, file=trim(file), form='formatted', iostat=istat)
   if (istat /= 0) then
     call error_handler("OPENING VARLIST FILE", istat)
   endif

   nfields = 0

   !Loop over lines of file to count the number of variables
   do
     read(14, '(A)', iostat=istat) line 
     if (istat/=0) exit
     if ( trim(line) .eq. '' ) cycle
     nfields = nfields+1
   enddo
   
   print*, "READING ", nfields, " FIELDS ACCORDING TO ", trim(file)
   if ( nfields == 0) call error_handler("VARLIST FILE IS EMPTY.", -1)

   allocate(field_names(nfields))
   allocate(field_names_target(nfields))

   rewind(14)
    do k = 1,nfields
      read(14, *, iostat=istat) field_names(k), field_names_target(k)
     if (istat /= 0) call error_handler("READING VARLIST FILE", istat)
    
    enddo
   close(14)

end subroutine read_varlist
 end module input_data
