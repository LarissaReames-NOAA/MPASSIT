!> @file
!! @brief Contains utility routines.
!!

!> General error handler.
!!
!! @param[in] string  error message
!! @param[in] rc      error status code
 subroutine error_handler(string, rc)

 use mpi

 implicit none

 character(len=*), intent(in)    :: string
 
 integer,          intent(in)    :: rc

 integer :: ierr

 print*,"- FATAL ERROR: ", string
 print*,"- IOSTAT IS: ", rc
 call mpi_abort(mpi_comm_world, 999, ierr)

 end subroutine error_handler

!> Error handler for netcdf
!!
!! @param[in] err     error status code
!! @param[in] string  error message
 subroutine netcdf_err( err, string )

 use mpi
 use netcdf

 implicit none
 integer, intent(in) :: err
 character(len=*), intent(in) :: string
 character(len=256) :: errmsg
 integer :: iret

 if( err.EQ.NF90_NOERR )return
 errmsg = NF90_STRERROR(err)
 print*,''
 print*,'FATAL ERROR: ', trim(string), ': ', trim(errmsg)
 print*,'STOP.'
 call mpi_abort(mpi_comm_world, 999, iret)

 return
 end subroutine netcdf_err
 

