!> @file
!! @brief Set up program execution.
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD

!> This module contains code to read the setup namelist file.
!!
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD
 module program_setup

 use esmf
 use ESMF_LogPublicMod

 implicit none

 private

 type(ESMF_LogKind_Flag), public :: LogType 
 character(len=500), public      :: grid_file_input_grid = "NULL" !< Full path of MPAS file containing grid information
 character(len=500), public      :: diag_file_input_grid = "NULL" !< Full path of input diagnostic MPAS data
 character(len=500), public      :: hist_file_input_grid = "NULL" !< Full path of input history MPAS data
 character(len=500), public      :: file_target_grid = "NULL" !<Full path of file containing target grid information
 character(len=500), public      :: output_file = "NULL" !< Full path of output file
 logical, public                 :: interp_diag = .false. !< Read data from diag file?
 logical, public                 :: interp_hist = .false. !< Read data from hist file?
 logical, public                 :: wrf_mod_vars = .false. !< Whether to modify variable values/dimensions 
                                                           !< to conform to WRF format. Set to true for
                                                           !< UPP-compatible output
                                                           

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


 namelist /config/ grid_file_input_grid, diag_file_input_grid, hist_file_input_grid, &
            file_target_grid, output_file, interp_diag, interp_hist, &
                        wrf_mod_vars, esmf_log

 print*,"- READ SETUP NAMELIST"

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

 return
 
 end subroutine read_setup_namelist

 end module program_setup
