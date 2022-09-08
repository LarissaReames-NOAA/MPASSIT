!> @file
!! @brief Interpolate atmospheric and surface data from MPAS diag, init, or forecast files
!! to a target grid.
!! @author  Larissa Reames CIWRO/NOAA/NSSL 

!> Public variables are defined below: "target" indicates field
!! associated with the input grid.
!!
!! @author  Larissa Reames CIWRO/NOAA/NSSL 

 module interp

 use esmf
 use netcdf

 use program_setup, only          : init_file_input_grid, &
                                    diag_file_input_grid, &
                                    fcst_file_input_grid, &
                                    data_to_interp

 use model_grid, only             : input_grid, target_grid, &
                                    nCells_input, nVert_input,  &
                                    nz_input, nzp1_input, &
                                    cell_latitude_input_grid, &
                                    cell_longitude_input_grid, &
                                    zgrid_input_grid, &
                                    zgrid_target_grid, &
                                    target_diag_bundle, & 
 									input_diag_bundle, &
 									input_init_bundle_2d, &
 									input_init_bundle_3d, &
 									target_init_bundle_2d, &
 									target_init_bundle_3d, &
 									n_diag_fields, nCellsPerPET, &
 									n_init_fields_2d, &
 									n_init_fields_3d

 implicit none

 private
 
 public :: interp_data
 
 contains
 
 subroutine interp_data(localpet)
 
	 implicit none
	 
	 integer, intent(in)   :: localpet
 
	 if (data_to_interp == 'diag') then
		call interp_diag_data(localpet)
	 elseif (data_to_interp == 'init') then
		call interp_init_data(localpet)
	 endif
 
 end subroutine interp_data
 
 subroutine interp_diag_data(localpet)
 
 	implicit none
 	
 	integer, intent(in)              :: localpet
 	integer 						 :: i, j, k,  rc
 	integer                          :: isrctermprocessing, nnodes
 	integer, allocatable             :: nodecoords(:)
 	type(ESMF_RegridMethod_Flag)     :: method
 	type(esmf_routehandle)           :: rh_diag
 	real(esmf_kind_r8), pointer :: temp2(:,:)
 	real(esmf_kind_r8),allocatable :: temp1(:), temp1lat(:), temp1lon(:)
 	
 	
 	call init_target_diag_fields
 	
 	isrctermprocessing = 1
 	method = ESMF_REGRIDMETHOD_PATCH
 	
 	print*,"- CREATE DIAG BUNDLE REGRID ROUTEHANDLE"
 	
 	call ESMF_FieldBundleRegridStore(input_diag_bundle, target_diag_bundle, &
 	                                 regridmethod=method, &
 	                                 routehandle=rh_diag, &
 	                                 srcTermProcessing=isrctermprocessing, &
 	                                 extrapMethod=ESMF_EXTRAPMETHOD_NONE, &
 	                                 rc=rc)
 	                                 
 	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldBundleRegridStore", rc)
 	  
 	  
 	 print*,"- REGRID DIAG FIELDS "                                
 	 call ESMF_FieldBundleRegrid(input_diag_bundle, target_diag_bundle, rh_diag, rc=rc)
 	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldBundleRegrid", rc)
    	
  end subroutine interp_diag_data
 		 		
  subroutine init_target_diag_fields
 
 	implicit none
 	
 	integer                     :: i, rc
 	character(50)               :: target_diag_names(n_diag_fields)
 	type(esmf_field)            :: diag_fields(n_diag_fields)
 	
 	target_diag_names = (/'U10','V10','Q2','TH2'/)
 	print*,"- INITIALIZE TARGET DIAG FIELDS."
 	do i = 1, n_diag_fields
 		print*, "- INIT FIELD ", target_diag_names(i)
 		diag_fields(i) = ESMF_FieldCreate(target_grid, & 
 							typekind=ESMF_TYPEKIND_R8, &
                            staggerloc=ESMF_STAGGERLOC_CENTER, &
							name=target_diag_names(i), rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldCreate", rc)
    enddo
    
    
 	target_diag_bundle = ESMF_FieldBundleCreate(fieldList=diag_fields, & 
 									name="target diag data", rc=rc)
 	if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldBundleCreate", rc)
 
 
 end subroutine init_target_diag_fields
 
  subroutine interp_init_data(localpet)
 
 	implicit none
 	
 	integer, intent(in)              :: localpet
 	integer 						 :: i, j, k,  rc
 	integer                          :: isrctermprocessing, nnodes
 	integer, allocatable              :: nodecoords(:)
 	type(esmf_routehandle)           :: rh_init
 	real(esmf_kind_r8), pointer :: temp2(:,:)
 	real(esmf_kind_r8),allocatable :: temp1(:), temp1lat(:), temp1lon(:)
 	
 	
 	call init_target_init_fields
 	
 	isrctermprocessing = 1
 	
 	print*,"- CREATE INIT BUNDLE REGRID ROUTEHANDLE"
 	
 	call ESMF_FieldBundleRegridStore(input_init_bundle_2d, target_init_bundle_2d, &
 	                                 regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
 	                                 polemethod=ESMF_POLEMETHOD_ALLAVG, &
 	                                 routehandle=rh_init, &
 	                                 srcTermProcessing=isrctermprocessing, &
 	                                 extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD, &
 	                                 rc=rc)
 	                                 
 	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldBundleRegridStore", rc)
 	  
 	  
 	 print*,"- REGRID INIT FIELDS "                                
 	 call ESMF_FieldBundleRegrid(input_init_bundle_2d, target_init_bundle_2d, rh_init, rc=rc)
 	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldBundleRegrid", rc)
    
     call ESMF_FieldBundleRegrid(input_init_bundle_3d, target_init_bundle_3d, rh_init, rc=rc)
 	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldBundleRegrid", rc)
    
    print*,"- CALL FieldRegridRelease."
    call ESMF_FieldRegridRelease(routehandle=rh_init, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegridRelease", rc)
    	
    call ESMF_FieldRegridStore(zgrid_input_grid, zgrid_target_grid, &
							 regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
							 polemethod=ESMF_POLEMETHOD_ALLAVG, &
							 routehandle=rh_init, &
							 srcTermProcessing=isrctermprocessing, &
							 extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD, &
							 rc=rc)
 	                                 
 	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldRegridStore", rc)
    	
    call ESMF_FieldRegrid(zgrid_input_grid, zgrid_target_grid, rh_init, rc=rc)
 	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldRegrid", rc)
    	
 print*,"- CALL FieldRegridRelease."
 call ESMF_FieldRegridRelease(routehandle=rh_init, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegridRelease", rc)
    
    	
  end subroutine interp_init_data
 		 		
  subroutine init_target_init_fields
 
 	implicit none
 	
 	integer                     :: i, rc
 	character(50)               :: target_init_names_2d(n_init_fields_2d), &
 								   target_init_names_3d(n_init_fields_3d)
 	type(esmf_field)            :: init_fields_2d(n_init_fields_2d), &
 	                               init_fields_3d(n_init_fields_3d)
 	
 	target_init_names_2d = (/'XLAND','U10'/)
 	target_init_names_3d = (/'QVAPOR','QCLOUD'/)
 	
 	print*,"- INITIALIZE TARGET DIAG FIELDS."
 	do i = 1, n_init_fields_2d
 		print*, "- INIT FIELD ", target_init_names_2d(i)
 		init_fields_2d(i) = ESMF_FieldCreate(target_grid, & 
 							typekind=ESMF_TYPEKIND_R8, &
                            staggerloc=ESMF_STAGGERLOC_CENTER, &
							name=target_init_names_2d(i), rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldCreate", rc)
    enddo
    
    
 	target_init_bundle_2d = ESMF_FieldBundleCreate(fieldList=init_fields_2d, & 
 									name="target init 2d data", rc=rc)
 	if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldBundleCreate", rc)
    	
    do i = 1, n_init_fields_3d
 		print*, "- INIT FIELD ", target_init_names_3d(i)
 		init_fields_3d(i) = ESMF_FieldCreate(target_grid, & 
 							typekind=ESMF_TYPEKIND_R8, &
                            staggerloc=ESMF_STAGGERLOC_CENTER, &
							name=target_init_names_3d(i), &
							ungriddedLBound=(/1/), &
                            ungriddedUBound=(/nz_input/), rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldCreate", rc)
    enddo
    
    
 	target_init_bundle_3d = ESMF_FieldBundleCreate(fieldList=init_fields_3d, & 
 									name="target init 3d data", rc=rc)
 	if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldBundleCreate", rc)
    	
    zgrid_target_grid = ESMF_FieldCreate(target_grid, & 
 							typekind=ESMF_TYPEKIND_R8, &
                            staggerloc=ESMF_STAGGERLOC_CENTER, &
							name='Z', &
							ungriddedLBound=(/1/), &
                            ungriddedUBound=(/nzp1_input/), rc=rc)
		if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldCreate", rc)
 
 
 end subroutine init_target_init_fields
end module interp
 
 