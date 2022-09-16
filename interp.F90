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

 use program_setup, only          : hist_file_input_grid, &
                                    diag_file_input_grid, &
                                    grid_file_input_grid, &
                                    interp_diag, interp_hist

 use model_grid, only             : input_grid, target_grid, &
                                    nCells_input, nVert_input,  &
                                    nz_input, nzp1_input, &
                                    nsoil_input, &
                                    cell_latitude_input_grid, &
                                    cell_longitude_input_grid, &
                                    zgrid_input_grid, &
                                    zgrid_target_grid, &
                                    target_diag_bundle, & 
                                    target_diag_names, &
                                    input_diag_bundle, &
                                    input_hist_bundle_2d_cons, &
                                    input_hist_bundle_2d_nstd, &
                                    input_hist_bundle_2d_patch, &
                                    input_hist_bundle_3d_nz, &
                                    input_hist_bundle_3d_nzp1, &
                                    input_hist_bundle_soil, &
                                    target_hist_bundle_2d_patch, &
                                    target_hist_bundle_2d_cons, &
                                    target_hist_bundle_2d_nstd, &
                                    target_hist_bundle_3d_nz, &  
                                    target_hist_bundle_3d_nzp1, &
                                    target_hist_bundle_soil, &
                                    target_diag_names, &
                                    target_hist_names_2d_cons, &
                                    target_hist_names_2d_nstd, &
                                    target_hist_names_2d_patch, &
                                    target_hist_names_3d_nz, &
                                    target_hist_names_3d_nzp1, &
                                    target_hist_names_soil, &
                                    n_diag_fields, nCellsPerPET, &
                                    n_hist_fields_2d_cons, &
                                    n_hist_fields_2d_nstd, &
                                    n_hist_fields_2d_patch, &
                                    n_hist_fields_3d_nz, &
                                    n_hist_fields_3d_nzp1, &
                                    n_hist_fields_soil

 implicit none

 private
 
 public :: interp_data
 
 type(esmf_routehandle)           :: rh_patch
 
 contains
 
 subroutine interp_data(localpet)
 
     implicit none
     
     integer, intent(in)   :: localpet
 
     if (interp_diag) then
        call interp_diag_data(localpet)
     endif
     if (interp_hist) then
        call interp_hist_data(localpet)
     endif
 
 end subroutine interp_data
 
 subroutine interp_diag_data(localpet)
 
    implicit none
    
    integer, intent(in)              :: localpet
    integer                          :: i, j, k,  rc
    integer                          :: isrctermprocessing, nnodes
    integer, allocatable             :: nodecoords(:)
    type(ESMF_RegridMethod_Flag)     :: method
    real(esmf_kind_r8), pointer :: temp2(:,:)
    real(esmf_kind_r8),allocatable :: temp1(:), temp1lat(:), temp1lon(:)
    
    
    call init_target_diag_fields
    
    isrctermprocessing = 1
    method = ESMF_REGRIDMETHOD_BILINEAR
    
    print*,"- CREATE DIAG BUNDLE REGRID ROUTEHANDLE"
    
    call ESMF_FieldBundleRegridStore(input_diag_bundle, target_diag_bundle, &
                                     regridmethod=method, &
                                     routehandle=rh_patch, &
                                     srcTermProcessing=isrctermprocessing, &
                                     rc=rc)
                                     
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleRegridStore", rc)
      
      
     print*,"- REGRID DIAG FIELDS "                                
     call ESMF_FieldBundleRegrid(input_diag_bundle, target_diag_bundle, rh_patch, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleRegrid", rc)
        
  end subroutine interp_diag_data
                
  subroutine init_target_diag_fields
 
    implicit none
    
    integer                     :: i, rc
    type(esmf_field)            :: diag_fields(n_diag_fields)
    
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
 
  subroutine interp_hist_data(localpet)
 
    implicit none
    
    integer, intent(in)              :: localpet
    integer                          :: rc, nfields
    integer                          :: isrctermprocessing
    type(ESMF_RegridMethod_Flag)     :: method
    type(ESMF_RouteHandle)           :: rh_cons, rh_nstd
    type(ESMF_field)                 :: fields(6)
    real(esmf_kind_r8), pointer      :: varptr2(:,:), varptr3(:,:,:)
    
    
    call init_target_hist_fields
    
    isrctermprocessing = 1
    
    !if (.not. interp_diag) then
        method = ESMF_REGRIDMETHOD_BILINEAR
        print*,"- CREATE HIST BUNDLE PATCH REGRID ROUTEHANDLE"
    
        call ESMF_FieldBundleRegridStore(input_hist_bundle_2d_patch, target_hist_bundle_2d_patch, &
                                         regridmethod=method, &
                                         routehandle=rh_patch, &
                                         srcTermProcessing=isrctermprocessing, &
                                         rc=rc)
                                     
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleRegridStore", rc)
    !endif
     
     
      
     print*,"- PATCH REGRID INIT FIELDS "                                
     call ESMF_FieldBundleRegrid(input_hist_bundle_2d_patch, target_hist_bundle_2d_patch, rh_patch, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleRegrid", rc)

     call ESMF_FieldBundleRegridStore(input_hist_bundle_3d_nz, target_hist_bundle_3d_nz, &
                                         regridmethod=method, &
                                         routehandle=rh_patch, &
                                         srcTermProcessing=isrctermprocessing, &
                                         rc=rc)

         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
            call error_handler("IN FieldBundleRegridStore", rc)
    
     call ESMF_FieldBundleRegrid(input_hist_bundle_3d_nz, target_hist_bundle_3d_nz, rh_patch, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleRegrid", rc)

    print*,"- CREATE HIST BUNDLE PATCH REGRID ROUTEHANDLE"
    
    call ESMF_FieldBundleRegridStore(input_hist_bundle_3d_nzp1, target_hist_bundle_3d_nzp1, &
                                         regridmethod=method, &
                                         routehandle=rh_patch, &
                                         srcTermProcessing=isrctermprocessing, &
                                         rc=rc)
                                     
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleRegridStore", rc)
        
    call ESMF_FieldBundleRegrid(input_hist_bundle_3d_nzp1, target_hist_bundle_3d_nzp1, rh_patch, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleRegrid", rc)
        
    print*,"- CREATE HIST BUNDLE CONSERVATIVE REGRID ROUTEHANDLE"
    method = ESMF_REGRIDMETHOD_CONSERVE
    call ESMF_FieldBundleRegridStore(input_hist_bundle_2d_cons, target_hist_bundle_2d_cons, &
                                         regridmethod=method, &
                                         routehandle=rh_cons, &
                                         srcTermProcessing=isrctermprocessing, &
                                         rc=rc)
                                     
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleRegridStore", rc)
        
    call ESMF_FieldBundleRegrid(input_hist_bundle_2d_cons, target_hist_bundle_2d_cons, rh_cons, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleRegrid", rc)
    
    print*,"- CREATE HIST BUNDLE NSTD REGRID ROUTEHANDLE"
    method = ESMF_REGRIDMETHOD_NEAREST_STOD
    call ESMF_FieldBundleRegridStore(input_hist_bundle_2d_nstd, target_hist_bundle_2d_nstd, &
                                         regridmethod=method, &
                                         routehandle=rh_nstd, &
                                         srcTermProcessing=isrctermprocessing, &
                                         rc=rc)
                                     
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleRegridStore", rc)
        
    call ESMF_FieldBundleRegrid(input_hist_bundle_2d_nstd, target_hist_bundle_2d_nstd, rh_nstd, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleRegrid", rc)
    
    if (n_hist_fields_soil>0) then  
        call ESMF_FieldBundleRegridStore(input_hist_bundle_soil, target_hist_bundle_soil, &
                                         regridmethod=method, &
                                         routehandle=rh_nstd, &
                                         srcTermProcessing=isrctermprocessing, &
                                         rc=rc)
        call ESMF_FieldBundleRegrid(input_hist_bundle_soil, target_hist_bundle_soil, rh_nstd, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleRegrid", rc)
    endif
        
    print*,"- CALL FieldRegridRelease."
    call ESMF_FieldBundleRegridRelease(routehandle=rh_patch, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegridRelease", rc)
    
    print*,"- CALL FieldRegridRelease."
    call ESMF_FieldBundleRegridRelease(routehandle=rh_cons, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegridRelease", rc)
      
      print*,"- CALL FieldRegridRelease."
    call ESMF_FieldBundleRegridRelease(routehandle=rh_nstd, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegridRelease", rc)
  end subroutine interp_hist_data
                
  subroutine init_target_hist_fields
 
    implicit none
    
    integer                      :: i, rc
    character(50), allocatable   :: field_names(:)
    type(esmf_field),allocatable :: fields(:)
    
    
    print*,"- INITIALIZE TARGET 2D CONS HIST FIELDS."
    
    if (n_hist_fields_2d_cons > 0) then
        allocate(fields(n_hist_fields_2d_cons))
        do i = 1, n_hist_fields_2d_cons
        
            print*, "- INIT FIELD ", target_hist_names_2d_cons(i)
            fields(i) = ESMF_FieldCreate(target_grid, & 
                                typekind=ESMF_TYPEKIND_R8, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, &
                                name=target_hist_names_2d_cons(i), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo
    
    
        target_hist_bundle_2d_cons = ESMF_FieldBundleCreate(fieldList=fields, & 
                                        name="target init 2d cons data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)
        deallocate(fields)
    endif
    
    print*,"- INITIALIZE TARGET 2D NSTD HIST FIELDS."
    
    if (n_hist_fields_2d_nstd > 0) then
        allocate(fields(n_hist_fields_2d_nstd))
        do i = 1, n_hist_fields_2d_nstd
        
            print*, "- INIT FIELD ", target_hist_names_2d_nstd(i)
            fields(i) = ESMF_FieldCreate(target_grid, & 
                                typekind=ESMF_TYPEKIND_R8, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, &
                                name=target_hist_names_2d_nstd(i), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo
    
    
        target_hist_bundle_2d_nstd = ESMF_FieldBundleCreate(fieldList=fields, & 
                                        name="target init 2d nstd data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)
        deallocate(fields)
    endif
    
    print*,"- INITIALIZE TARGET 2D PATCH HIST FIELDS."
    
    if (n_hist_fields_2d_patch > 0) then
        allocate(fields(n_hist_fields_2d_patch))
        do i = 1, n_hist_fields_2d_patch
        
            print*, "- INIT FIELD ", target_hist_names_2d_patch(i)
            fields(i) = ESMF_FieldCreate(target_grid, & 
                                typekind=ESMF_TYPEKIND_R8, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, &
                                name=target_hist_names_2d_patch(i), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo
    
    
        target_hist_bundle_2d_patch = ESMF_FieldBundleCreate(fieldList=fields, & 
                                        name="target init 2d patch data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)
        deallocate(fields)
    endif
    
    if (n_hist_fields_soil>0) then
        allocate(fields(n_hist_fields_soil))    
        do i = 1, n_hist_fields_soil
            print*, "- INIT FIELD ", target_hist_names_soil(i)
            fields(i) = ESMF_FieldCreate(target_grid, & 
                                typekind=ESMF_TYPEKIND_R8, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, &
                                name=target_hist_names_soil(i), &
                                ungriddedLBound=(/1/), &
                                ungriddedUBound=(/nsoil_input/), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo
    
    
        target_hist_bundle_soil= ESMF_FieldBundleCreate(fieldList=fields, & 
                                        name="target init soil data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)
        deallocate(fields)
    endif
    
    if (n_hist_fields_3d_nz>0) then
        allocate(fields(n_hist_fields_3d_nz))   
        do i = 1, n_hist_fields_3d_nz
            print*, "- INIT FIELD ", target_hist_names_3d_nz(i)
            fields(i) = ESMF_FieldCreate(target_grid, & 
                                typekind=ESMF_TYPEKIND_R8, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, &
                                name=target_hist_names_3d_nz(i), &
                                ungriddedLBound=(/1/), &
                                ungriddedUBound=(/nz_input/), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo
    
    
        target_hist_bundle_3d_nz = ESMF_FieldBundleCreate(fieldList=fields, & 
                                        name="target init 3d nz data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)
        deallocate(fields)
    endif
    
    if (n_hist_fields_3d_nzp1>0) then
        allocate(fields(n_hist_fields_3d_nzp1)) 
        do i = 1, n_hist_fields_3d_nzp1
            print*, "- INIT FIELD ", target_hist_names_3d_nzp1(i)
            fields(i) = ESMF_FieldCreate(target_grid, & 
                                typekind=ESMF_TYPEKIND_R8, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, &
                                name=target_hist_names_3d_nzp1(i), &
                                ungriddedLBound=(/1/), &
                                ungriddedUBound=(/nzp1_input/), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo
    
    
        target_hist_bundle_3d_nzp1 = ESMF_FieldBundleCreate(fieldList=fields, & 
                                        name="target init 3d nz data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)
        deallocate(fields)
    endif
    deallocate(target_hist_names_2d_cons, target_hist_names_2d_nstd)
    deallocate(target_hist_names_2d_patch, target_hist_names_3d_nz)
    deallocate(target_hist_names_3d_nzp1)
 
 end subroutine init_target_hist_fields
end module interp
 
 
