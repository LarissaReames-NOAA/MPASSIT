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
                                    target_diag_bundle, & 
 									input_diag_bundle, &
 									n_diag_fields, nCellsPerPET

 implicit none

 private
 
 
 public :: interp_data
 
 contains
 
 subroutine interp_data(localpet)
 
	 implicit none
	 
	 integer, intent(in)   :: localpet
 
	 if (data_to_interp == 'diag') then
		call interp_diag_data(localpet)
	 endif
 
 end subroutine interp_data
 
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
 
 subroutine interp_diag_data(localpet)
 
 	implicit none
 	
 	integer, intent(in)              :: localpet
 	integer 						 :: i, j, k,  rc
 	integer                          :: isrctermprocessing, nnodes
 	integer, allocatable              :: nodecoords(:)
 	type(esmf_routehandle)           :: rh_diag
 	real(esmf_kind_r8), pointer :: temp2(:,:)
 	real(esmf_kind_r8),allocatable :: temp1(:), temp1lat(:), temp1lon(:)
 	
 	
 	call init_target_diag_fields
 	
 	isrctermprocessing = 1
 	
 	
 	print*,"- CREATE DIAG BUNDLE REGRID ROUTEHANDLE"
 	
 	call ESMF_FieldBundleRegridStore(input_diag_bundle, target_diag_bundle, &
 	                                 regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
 	                                 polemethod=ESMF_POLEMETHOD_ALLAVG, &
 	                                 routehandle=rh_diag, &
 	                                 srcTermProcessing=isrctermprocessing, &
 	                                 extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD, &
 	                                 rc=rc)
 	                                 
 	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldBundleRegridStore", rc)
 	  
 	  
 	 print*,"- REGRID DIAG FIELDS "                                
 	 call ESMF_FieldBundleRegrid(input_diag_bundle, target_diag_bundle, rh_diag, rc=rc)
 	 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    	call error_handler("IN FieldBundleRegrid", rc)
    	
  end subroutine interp_diag_data
 		
end module interp
 
 