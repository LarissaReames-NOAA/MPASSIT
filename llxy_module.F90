!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE LLXY_MODULE
!
! This module handles transformations between model grid coordinates and 
!   latitude-longitude coordinates. The actual transformations are done through
!   the map_utils module. 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module llxy_module

   use program_setup
   use map_utils_mod
   use utils_mod
   use misc_definitions_module
 
   ! Parameters
   integer, parameter :: MAX_SOURCE_LEVELS = 20
 
   ! Variables
   integer :: current_nest_number
   integer :: SOURCE_PROJ = 0
   ! The following arrays hold values for all available domains 
   ! NOTE: The entries in the arrays for "domain 0" are used for projection
   !       information of user-specified source data
   type (proj_info),public :: proj_stack
 
   ! The projection and domain that we have computed constants for
   integer :: computed_proj = INVALID
   integer :: computed_domain = INVALID
 
   contains
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! Name: push_source_projection
   !
   ! Purpose: 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   subroutine push_source_projection(iprojection, user_stand_lon, user_truelat1, user_truelat2, &
                                  user_dxkm, user_dykm, user_dlat, user_dlon, user_known_x, &
                                  user_known_y, user_known_lat, user_known_lon, &
                                  user_pole_lat, user_pole_lon, &   
                                  user_centerlat, user_centerlon, &   
                                  user_centeri, user_centerj, &   
                                  earth_radius)
 
      implicit none
  
      ! Arguments
      integer, intent(in) :: iprojection
      real, intent(in) :: user_stand_lon, user_truelat1, user_truelat2, user_dxkm, user_dykm, &
                          user_dlat, user_dlon, &
                          user_known_x, user_known_y, user_known_lat, user_known_lon
      real, intent(in), optional :: earth_radius
      real, intent(in), optional :: user_centerlon, user_centerlat, user_pole_lat, user_pole_lon 
      real, intent(in), optional :: user_centerj, user_centeri   
      
  
      call map_init(proj_stack)

      if (iprojection == PROJ_LATLON) then
         call map_set(iprojection, proj_stack, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      nxmax=nint(360.0 / user_dlon), &
                      latinc=user_dlat, &
                      loninc=user_dlon, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_MERC) then
         call map_set(iprojection, proj_stack, &
                      truelat1=user_truelat1, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      dx=user_dxkm, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_CYL) then
        call error_handler('Should not have PROJ_CYL as projection for ' &
                          //'target data in push_source_projection().',ERROR_CODE)
  
      else if (iprojection == PROJ_CASSINI) then
 
         call map_set(iprojection, proj_stack, &
                      latinc=user_dlat, &
                      loninc=user_dlon, &
                      stdlon=user_stand_lon, &  
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      lat0=user_pole_lat, &
                      lon0=user_pole_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_LC) then
         call map_set(iprojection, proj_stack, &
                      truelat1=user_truelat1, &
                      truelat2=user_truelat2, &
                      stdlon=user_stand_lon, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      dx=user_dxkm, &
                      r_earth=earth_radius)

      else if (iprojection == PROJ_ALBERS_NAD83) then
         call map_set(iprojection, proj_stack, &
                      truelat1=user_truelat1, &
                      truelat2=user_truelat2, &
                      stdlon=user_stand_lon, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      dx=user_dxkm, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_PS) then
         call map_set(iprojection, proj_stack, &
                      truelat1=user_truelat1, &
                      stdlon=user_stand_lon, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      dx=user_dxkm, &
                      r_earth=earth_radius)

      else if (iprojection == PROJ_PS_WGS84) then
         call map_set(iprojection, proj_stack, &
                      truelat1=user_truelat1, &
                      stdlon=user_stand_lon, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      dx=user_dxkm, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_GAUSS) then
         call map_set(iprojection, proj_stack, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      nxmax=nint(360.0 / user_dlon), &
                      nlat=nint(user_dlat), &
                      loninc=user_dlon, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_ROTLL) then
  ! BUG: Implement this projection.
  
      end if
     
   end subroutine push_source_projection
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! Name: lltoxy
   !
   ! Purpose:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   subroutine xytoll(x, y, xlat, xlon, stagger, comp_ll)
 
      implicit none
  
      ! Arguments
      integer, intent(in) :: stagger
      real, intent(in) :: x, y
      real, intent(out) :: xlat, xlon
      logical, optional, intent(in) :: comp_ll

      ! Local variables
      real :: rx, ry
      logical :: save_comp_ll
  
      ! Account for grid staggering; we cannot modify x and y, so modify local
      !   copies of them
      if (stagger == U) then
         rx = x - 0.5
         ry = y
      else if (stagger == V) then
         rx = x
         ry = y - 0.5
      else if (stagger == HH) then
         proj_stack%stagger = HH
         rx = x
         ry = y
      else if (stagger == VV) then
         proj_stack%stagger = VV
         rx = x
         ry = y
      else if (stagger == CORNER) then
         proj_stack%stagger = CORNER
         rx = x - 0.5
         ry = y - 0.5
      else
         rx = x
         ry = y
      end if

      if (present(comp_ll)) then
         save_comp_ll = proj_stack%comp_ll
         proj_stack%comp_ll = comp_ll
      end if
  
      call ij_to_latlon(proj_stack, rx, ry, xlat, xlon)

      if (present(comp_ll)) then
         proj_stack%comp_ll = save_comp_ll
      end if
 
   end subroutine xytoll

end module llxy_module
