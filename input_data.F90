!> @file
!! @brief Read atmospheric and surface data from MPAS diag, init, or forecast files.
!! @author  Larissa Reames CIWRO/NOAA/NSSL 

!> Read atmospheric and/or surface input MPAS grid.
!! Supported formats include initialization file, forecast diag file (only 2d diagnostic
!! fields), and full 3d forecast files.
!!
!! Public variables are defined below: "input" indicates field
!! associated with the input grid.
!!
!! @author  Larissa Reames CIWRO/NOAA/NSSL 

 module input_data

 use esmf
 use netcdf

 use program_setup, only          : init_file_input_grid, &
                                    diag_file_input_grid, &
                                    fcst_file_input_grid, &
                                    data_to_interp

 use model_grid, only             : input_grid,        &
                                    nCells_input, nVert_input,  &
                                    nz_input, nzp1_input, &
                                    cell_latitude_input_grid, &
                                    cell_longitude_input_grid, &
                                    zgrid_input_grid, &
                                    input_diag_bundle, &
                                    input_init_bundle_2d, &
                                    input_init_bundle_3d, &
                                    n_diag_fields, &
                                    n_init_fields_2d, &
                                    n_init_fields_3d, &
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
 
 if (data_to_interp == 'diag') then
 	call read_input_diag_data(localpet)
 elseif (data_to_interp == 'init') then
 	call read_input_init_data(localpet)
 else
 	call error_handler(" ONLY INTERPOLATION OF DIAG OR INIT FILE IS CURRENTLY SUPPORTED")
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
 
 type(esmf_field)                :: fields(n_diag_fields)
 
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
 	
 	integer                     :: i, rc
 	type(esmf_field)            :: diag_fields(n_diag_fields)
 	character(50)			:: input_diag_names(n_diag_fields)
 	
 	input_diag_names = (/'u10','v10','q2','th2m'/)
 	
 	print*,"- INITIALIZE INPUT DIAG FIELDS."
 	do i = 1, n_diag_fields
 	
 		print*, "- INIT FIELD ", input_diag_names(i)
 		
 		if (input_diag_names(i) == 'vorticity_500hPa') then
 			diag_fields(i) = ESMF_FieldCreate(input_grid, & 
 							typekind=ESMF_TYPEKIND_R8, &
                            meshloc=ESMF_MESHLOC_NODE, &
							name=input_diag_names(i), rc=rc)
			if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
			call error_handler("IN FieldCreate", rc)
 		else
			diag_fields(i) = ESMF_FieldCreate(input_grid, & 
								typekind=ESMF_TYPEKIND_R8, &
								meshloc=ESMF_MESHLOC_ELEMENT, &
								name=input_diag_names(i), rc=rc)
			if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
			call error_handler("IN FieldCreate", rc)
		endif
    enddo
    
    
 	input_diag_bundle = ESMF_FieldBundleCreate(fieldList=diag_fields, & 
 									name="input diag data", rc=rc)
 	if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldBundleCreate", rc)
 
 
 end subroutine init_input_diag_fields
 
!> Read input grid init data.
!!
!! @param[in] localpet  ESMF local persistent execution thread 
!! @author Larissa Reames CIWRO/NOAA/NSSL
 subroutine read_input_init_data(localpet)

 character(len=500)              :: the_file
 character(len=50)               :: vname
 
 integer, intent(in)             :: localpet
 integer                         :: error, ncid, rc
 integer                         :: id_dim
 integer                         :: id_var, i, j, nodes
 
 type(esmf_field)                :: fields(n_diag_fields)
 
 real(esmf_kind_r8), allocatable :: dummy2(:,:), dummy3(:,:,:)

 real(esmf_kind_r8), pointer     :: varptr(:), varptr2(:,:)
 
 call init_input_init_fields()

 print*,"- READ INPUT INIT DATA."

 the_file = trim(init_file_input_grid)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening: '//trim(the_file) )


!---------------------------------------------------------------------------
! Initialize 2d esmf atmospheric fields.
!---------------------------------------------------------------------------

 call ESMF_FieldBundleGet(input_init_bundle_2d, fieldList=fields, & 
 						  itemorderflag=ESMF_ITEMORDER_ADDORDER, &
 						  rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
	 call error_handler("IN FieldBundleGet", rc)
	 
 call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
	 call error_handler("IN MeshGet", rc)
 
 allocate(dummy2(nCells_input,1))
 
 do i = 1,n_init_fields_2d

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
		
	print*,"- SET ON MESH ", trim(vname)
	do j = 1, nCellsPerPET
		varptr(j) = dummy2(elemIDs(j),1)
	enddo
	

	
	print*, localpet, minval(varptr), maxval(varptr)	
	nullify(varptr)
 enddo
 deallocate(dummy2)
!---------------------------------------------------------------------------
! Initialize 3d esmf atmospheric fields.
!---------------------------------------------------------------------------

 call ESMF_FieldBundleGet(input_init_bundle_3d, fieldList=fields, & 
 						  itemorderflag=ESMF_ITEMORDER_ADDORDER, &
 						  rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
	 call error_handler("IN FieldBundleGet", rc)
	 
 call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
	 call error_handler("IN MeshGet", rc)
 
 allocate(dummy3(nz_input,nCells_input,1))
 
 do i = 1,n_init_fields_3d

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
		
	print*,"- SET ON MESH ", trim(vname)
	do j = 1, nCellsPerPET
		varptr2(j,:) = dummy3(:,elemIDs(j),1)
	enddo
	
	print*, localpet, minval(varptr2), maxval(varptr2)	
	nullify(varptr2)
 enddo
 deallocate(dummy3)
 allocate(dummy2(nzp1_input,nCells_input))

 	
 call ESMF_FieldGet(zgrid_input_grid, farrayPtr=varptr2, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
  call error_handler("IN FieldGet", rc)

 print*,"- READ INPUT ZGRID"
 error=nf90_inq_varid(ncid, 'zgrid', id_var)
 call netcdf_err(error, 'reading field id' )
 error=nf90_get_var(ncid, id_var, dummy2)
 call netcdf_err(error, 'reading field' )
	
 print*,"- SET ON MESH ZGRID"
 do j = 1, nCellsPerPET
 	varptr2(j,:) = dummy2(:,elemIDs(j))
 enddo
 
 deallocate(dummy2)
 print*, localpet, minval(varptr2), maxval(varptr2)	
 nullify(varptr2)

 error = nf90_close(ncid)


 end subroutine read_input_init_data
 
 subroutine init_input_init_fields
 
 	implicit none
 	
 	integer                     :: i, rc
 	type(esmf_field)            :: init_fields_2d(n_init_fields_2d), &
 								   init_fields_3d(n_init_fields_3d), &
 								   zgrid
 	character(50)			    :: input_init_names_2d(n_init_fields_2d), &
 								   input_init_names_3d(n_init_fields_3d)
 	
 	input_init_names_2d = (/'xland','u10'/)
 	input_init_names_3d = (/'qv','qc'/)
 	
 	print*,"- INITIALIZE INPUT DIAG FIELDS."
 	do i = 1, n_init_fields_2d
 	
 		print*, "- INIT FIELD ", input_init_names_2d(i)

		init_fields_2d(i) = ESMF_FieldCreate(input_grid, & 
							typekind=ESMF_TYPEKIND_R8, &
							meshloc=ESMF_MESHLOC_ELEMENT, &
							name=input_init_names_2d(i), rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
		call error_handler("IN FieldCreate", rc)
    enddo
    
    
 	input_init_bundle_2d = ESMF_FieldBundleCreate(fieldList=init_fields_2d, & 
 									name="input init 2d data", rc=rc)
 	if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldBundleCreate", rc)
    	

	do i = 1, n_init_fields_3d
 	
 		print*, "- INIT FIELD ", input_init_names_3d(i)

		init_fields_3d(i) = ESMF_FieldCreate(input_grid, & 
							typekind=ESMF_TYPEKIND_R8, &
							meshloc=ESMF_MESHLOC_ELEMENT, &
							name=input_init_names_3d(i), & 
							ungriddedLBound=(/1/), &
                            ungriddedUBound=(/nz_input/), rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
		call error_handler("IN FieldCreate", rc)
    enddo
    
    
 	input_init_bundle_3d = ESMF_FieldBundleCreate(fieldList=init_fields_3d, & 
 									name="input init 3d data", rc=rc)
 	if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldBundleCreate", rc)
    	
    zgrid_input_grid = ESMF_FieldCreate(input_grid, & 
							typekind=ESMF_TYPEKIND_R8, &
							meshloc=ESMF_MESHLOC_ELEMENT, &
							name='zgrid',  &
							ungriddedLBound=(/1/), &
                            ungriddedUBound=(/nzp1_input/), rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
		call error_handler("IN FieldCreate", rc)
 
 
 end subroutine init_input_init_fields


 end module input_data
