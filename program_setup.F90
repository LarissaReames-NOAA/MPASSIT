!> @file
!! @brief Set up program execution.
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD

!> This module contains code to read the setup namelist file.
!!
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD
 module program_setup

 use esmf
 use ESMF_LogPublicMod
 use constants_module
 use utils_mod

 implicit none

 private

 type(ESMF_LogKind_Flag), public :: LogType 
 
 ! Namelist variables
 character(len=500), public      :: grid_file_input_grid = "NULL" !< Full path of MPAS file containing grid information
 character(len=500), public      :: diag_file_input_grid = "NULL" !< Full path of input diagnostic MPAS data
 character(len=500), public      :: hist_file_input_grid = "NULL" !< Full path of input history MPAS data
 character(len=500), public      :: file_target_grid = "NULL" !<Full path of file containing target 
 															  !<grid information for target_grid_type='file'
 character(len=500), public      :: output_file = "NULL" !< Full path of output file
 logical, public                 :: interp_diag = .false. !< Read data from diag file?
 logical, public                 :: interp_hist = .false. !< Read data from hist file?
 logical, public                 :: wrf_mod_vars = .false. !< Whether to modify variable values/dimensions 
                                                           !< to conform to WRF format. Set to true for
                                                           !< UPP-compatible output
 character(len=500), public      :: target_grid_type	  !< Grid type to interpolate data to
 														  !< Valid options: 'file', 'lcc','ll','ll_global'														  
 !! These entries are only valid for target_grid_type = 'lcc','ll','ll_global'
 integer, public 				 :: i_target			  !< # staggered east-west grid points in target grid
 integer, public				 :: j_target     		  !< # staggered north-south grid points in target grid									  
 real, public					 :: truelat1 = NAN     	  !< First true latitude (all projections)
 real, public					 :: truelat2 = NAN     	  !< Second true latitude (LCC only)
 real, public					 :: stand_lon = NAN		  !< Longitude parallel to y-axis (-180->180E)
 real, public					 :: dx = NAN			  !< Grid cell east-west dimension(meters)
 real, public					 :: dy = NAN			  !< Grid cell north-south dimension(meters)
 real, public					 :: ref_lat				  !< Latitude of reference point
 real, public					 :: ref_lon   	 	      !< Longitude of reference point
 real, public					 :: ref_x				  !< Grid-relative e-w index of reference point
 														  !< Defaults to grid center (nx/2)
 real, public					 :: ref_y		  		  !< Grid-relative n-s index of reference point
 														  !< Defaults to grid center (ny/2)	
  
 
 !These aren't namelist variables but they're created directly from them
 real, public					 :: dxkm				  !< grid-cell east-west dimension (meters)
 real, public					 :: dykm 				  !< grid-cell north-south dimension (meters)
 real, public					 :: dlondeg				  !< grid-cell east-west dimension (deg)
 real, public					 :: dlatdeg				  !< grid-cell north-south dimension (deg)
 real, public					 :: known_lat			  !< Latitude of reference point
 real, public					 :: known_lon   	      !< Longitude of reference point
 real, public					 :: known_x				  !< Grid-relative e-w index of reference point
 real, public					 :: known_y		  		  !< Grid-relative n-s index of reference point	
 real, public					 :: proj_code			  !< Integer code corresponding to the requested
 														  ! target grid projection type
  
 public :: read_setup_namelist

 contains

!> Reads program configuration namelist.
!!
!! @param filename the name of the configuration file (defaults to
!! ./fort.41).
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD
 subroutine read_setup_namelist(unum, filename)
 implicit none

 character(len=*), intent(in), optional :: filename
 integer, intent(in), optional          :: unum
 character(:), allocatable :: filename_to_use
 integer                   :: unit_to_use
 logical                   :: esmf_log

 integer                     :: is, ie, ierr
 
 !Namelist variables that are used to create global variables
 real						:: dx,dy


 namelist /config/ grid_file_input_grid, diag_file_input_grid, hist_file_input_grid, &
            file_target_grid, output_file, interp_diag, interp_hist, &
            wrf_mod_vars, esmf_log,target_grid_type,nx,ny,dx,dy,ref_lat,ref_lon,ref_x,ref_y,&
            truelat1,truelat2,stand_lon

  ref_x = NAN
  ref_y = NAN
  ref_lat = NAN
  ref_lon = NAN
  dx = NAN
  dy = NAN 
  pole_lat = 90.0
  pole_lon = 0.0                   

 !print*,"- READ SETUP NAMELIST"

 if (present(filename)) then
    filename_to_use = filename
 else
    filename_to_use = "./fort.41"
 endif

 if (present(unum)) then
     unit_to_use = unum
 else
     unit_to_use = 41
 endif

 open(unit_to_use, file=filename_to_use, iostat=ierr)
 if (ierr /= 0) call error_handler("OPENING SETUP NAMELIST.", ierr)
 read(unit_to_use, nml=config, iostat=ierr)
 if (ierr /= 0) call error_handler("READING SETUP NAMELIST.", ierr)
 close (unit_to_use)
 
 if (esmf_log) then
   LogType = ESMF_LOGKIND_MULTI_ON_ERROR
 else
   LogType = ESMF_LOGKIND_NONE
 endif 
 
 dxkm = dx
 dykm = dy
 
 known_lat = ref_lat
 known_lon = ref_lon
 known_x = ref_x
 known_y = ref_y
 
 if (trim(target_grid_type)=='ll' .or. trim(target_grid_type)=='ll_global') then

	 ! If no dx,dy specified, assume global grid
	 if (dx == NAN .and. dy == NAN) then
	    if (trim(target_grid_type) .ne. 'll_global') then
	    	call error_handler("For lat-lon projection, if dx/dy are not specified "// &
	    	"a global grid is assumed. Please set dx/dy or change target_grid_type to "// &
	    	"'ll_global'", ERROR)
	    endif
		dlondeg = 360. / (nx)   ! Here, we really do not want e_we-s_we+1
		dlatdeg = 180. / (ny)   ! Here, we really do not want e_we-s_we+1
		known_x = 1.
		known_y = 1.
		known_lon = stand_lon + dlondeg/2.
		known_lat = -90. + dlatdeg/2.
		dxkm = EARTH_RADIUS_M * PI * 2.0 / (e_we(1)-s_we(1))
		dykm = EARTH_RADIUS_M * PI       / (e_sn(1)-s_sn(1))

	 ! If dx,dy specified, however, assume regional grid
	 else
		dlatdeg = dy
		dlondeg = dx
		dxkm = dlondeg * EARTH_RADIUS_M * PI * 2.0 / 360.0
		dykm = dlatdeg * EARTH_RADIUS_M * PI * 2.0 / 360.0
		if (known_lat == NAN .or. known_lon == NAN) then
		   call error_handler('For lat-lon projection, if dx/dy are specified, '// &
					'a regional domain is assumed, and a ref_lat,ref_lon must also be specified',ERROR)
		end if
	 end if
 end if
 
 ! Manually set truelat2 = truelat1 if truelat2 not specified for Lambert
  if (trim(target_grid_type) == 'lcc' .and. truelat2 == NAN) then
	 if (truelat1 == NAN) call error_handler("No TRUELAT1 specified for Lambert conformal projection.",ERROR)) 
	 truelat2 = truelat1
  end if
  
  ! If the user hasn't supplied a known_x and known_y, assume the center of domain 1
  if (known_x == NAN .and. known_y == NAN) then
	known_x = ixdim(1) / 2.
	known_y = jydim(1) / 2.
  else if (known_x == NAN .or. known_y == NAN) then
	call error_handler('In namelist, neither or both of ref_x, ref_y must be specified.',ERROR)
  end if 

 return
 
 end subroutine read_setup_namelist

 end module program_setup
