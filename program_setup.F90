!> @file
!! @brief Set up program execution.
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD

!> This module contains code to read the setup namelist file.
!!
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD
 module program_setup

 implicit none

 private
 
 character(len=500), public      :: grid_file_input_grid = "NULL" !< Full path of MPAS file containing grid information
 character(len=500), public      :: diag_file_input_grid = "NULL" !< Full path of input diagnostic MPAS data
 character(len=500), public      :: hist_file_input_grid = "NULL" !< Full path of input history MPAS data
 character(len=500), public      :: file_target_grid = "NULL" !<Full path of file containing target grid information
 character(len=500), public      :: output_file = "NULL" !< Full path of output file
 logical, public      			 :: interp_diag = .false. !< Read data from diag file?
 logical, public        		 :: interp_hist = .false. !< Read data from hist file?

 public :: read_setup_namelist


 contains

!> Reads program configuration namelist.
!!
!! @param filename the name of the configuration file (defaults to
!! ./fort.41).
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD
 subroutine read_setup_namelist(filename)
 implicit none

 character(len=*), intent(in), optional :: filename
 character(:), allocatable :: filename_to_use
 

 integer                     :: is, ie, ierr


 namelist /config/ grid_file_input_grid, diag_file_input_grid, hist_file_input_grid, &
 			file_target_grid, output_file, interp_diag, interp_hist

 print*,"- READ SETUP NAMELIST"

 if (present(filename)) then
    filename_to_use = filename
 else
    filename_to_use = "./fort.41"
 endif

 open(41, file=filename_to_use, iostat=ierr)
 if (ierr /= 0) call error_handler("OPENING SETUP NAMELIST.", ierr)
 read(41, nml=config, iostat=ierr)
 if (ierr /= 0) call error_handler("READING SETUP NAMELIST.", ierr)
 close (41)
 
 
 return
 
 end subroutine read_setup_namelist

 end module program_setup
