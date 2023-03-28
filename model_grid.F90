!> @file
!! @brief Specify input and target model grids.
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD

!> Sets up the ESMF grid objects for the input data grid and target grid.
!!
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD
 module model_grid

 use esmf
 use ESMF_LogPublicMod

 implicit none

 private

 integer, public                        :: nCells_input
                                           !< cells on input grid
 integer, public                        :: nVert_input
                                           !< vertices on input grid
 integer, public                        :: nz_input
                                           !< number of input grid atm layers
 integer, public                        :: nzp1_input
                                           !< number of input grid atm layer interfaces
 integer, public                        :: nsoil_input
                                           !< number of input soil levels                                           
 real, public                           :: dx
                                           !< grid size (m) of target grid
 character(50), public                  :: start_time
                                           !< simulation start time input data
 character(20), public,allocatable      :: valid_time(:,:)
                                           !< valid forecast time input data
 integer, public                        :: strlen
                                           !< StrLen on input file
 real, public                           :: config_dt
                                           !< timestep (seconds) of input data
 integer, public                        :: lsm_scheme
                                           !< land surface scheme input data
 integer, public                        :: mp_scheme
                                           !< microphysics scheme input data
 integer, public                        :: conv_scheme
                                           !< convection scheme input data 
 real, public                           :: cen_lat
                                           !< target grid projection center latitude
 real, public                           :: cen_lon
                                           !< target grid projection center longitude
 real, public                           :: truelat1
                                           !< target grid projection 1st true latitude
 real, public                           :: truelat2
                                           !< target grid projection 2nd true latitude
 real, public                           :: moad_cen_lat
                                           !< target grid moad center latitude
 real, public                           :: stand_lon
                                           !< target grid projection standard longitude
 real, public                           :: pole_lat
                                           !< target grid projection pole latitude
 real, public                           :: pole_lon
                                           !< target grid projection pole longitude
 real,allocatable,public                :: latitude_u(:,:)
                                           !< target grid latitude on the
                                           !< staggered u grid
 real,allocatable,public                :: latitude_v(:,:)
                                           !< target grid latitude on the
                                           !< staggered v grid
 real,allocatable,public                :: longitude_u(:,:)
                                           !< target grid longitude on the
                                           !< staggered u grid
 real,allocatable,public                :: longitude_v(:,:)
                                           !< target grid longitude on the
                                           !< staggered v grid
 integer, public                        :: map_proj
                                           !< target grid map projection integer label
 character(50), public                  :: map_proj_char
                                           !< target grid map projection character string
 integer, public                        :: i_target
                                           !< i dimension of each global tile, 
                                           !! or of a nest, target grid.
 integer, public                        :: j_target
                                           !< j dimension of each global tile,
                                           !! or of a nest, target grid.
 integer, public                        :: ip1_target
                                           !< ip1_target plus 1
 integer, public                        :: jp1_target
                                           !< jp1_target plus 1
                                           
 integer, allocatable, public           :: elemIDs(:)
                                           !< IDs of the elements on present PET   
 integer, allocatable, public           :: nodeIDs(:)
                                            !< IDs of the nodes on present PET  
 integer, public                        :: nCellsPerPET
                                            !< Number of cells on this PET
                                            
 type(esmf_mesh),  public               :: input_grid
                                           !< input grid esmf grid object
 type(esmf_grid),  public               :: target_grid
                                           !< target grid esmf grid object.

 type(esmf_field),  public              :: cell_latitude_input_grid
                                           !< latitude of grid center, input grid
 type(esmf_field),  public              :: cell_longitude_input_grid
                                           !< longitude of grid center, input grid
                                           
 type(esmf_field),  public              :: node_latitude_input_grid
                                           !< latitude of grid center, input grid
 type(esmf_field),  public              :: node_longitude_input_grid
                                           !< longitude of grid center, input grid

  type(esmf_field),  public              :: zgrid_input_grid
                                           !< esmf field to hold level height on input grid                                            
  type(esmf_field),  public              :: zgrid_target_grid
                                           !< esmf field to hold level height on target grid 

 type(esmf_field),  public              :: latitude_target_grid
                                           !< latitude of grid center, target grid
 type(esmf_field),  public              :: longitude_target_grid
                                           !< longitude of grid center, target grid
 real(esmf_kind_r8), allocatable , public :: zs_target_grid(:,:)
                                          !< soil center depth, target grid
 type(esmf_field), public               :: hgt_input_grid, hgt_target_grid
                                          !< surface elevation, target grid
                                           
 integer, public                       :: n_diag_fields
                                          !< number of fields read from the diag file
 type(esmf_fieldbundle), public        :: input_diag_bundle    
                                          !< bundle to hold input diag fields
 type(esmf_fieldbundle), public        :: target_diag_bundle
                                          !< bundle to hold target diag fields
                                          
 integer, public                       :: n_hist_fields_2d_patch
                                          !< number of 2d fields read from the hist file
                                          !< to use with patch regridding
 integer, public                       :: n_hist_fields_2d_cons
                                          !< number of 2d fields read from the hist file
                                          !< to use with conservative regridding
 integer, public                       :: n_hist_fields_2d_nstd
                                          !< number of 2d fields read from the hist file
                                          !< to use with nearest source to destination 
                                          !< regridding
 integer, public                       :: n_hist_fields_3d_nz
                                          !< number of 3d fields read from the hist file
                                          !< with vertical dimension nVertLevels
 integer, public                       :: n_hist_fields_3d_nzp1
                                          !< number of 3d fields read from the hist file
                                          !< with vertical dimension nVertLevelsp1
 integer, public                      :: n_hist_fields_soil
                                          !< number of soil fields read from the hist file

 integer, public                       :: diag_out_interval
                                         !< output_interval from diag file
 character(50), allocatable, public    :: target_diag_names(:), &
                                          target_hist_names_2d_cons(:), &
                                          target_hist_names_2d_nstd(:), &
                                          target_hist_names_2d_patch(:), &
                                          target_hist_names_3d_nzp1(:), &
                                          target_hist_names_3d_nz(:), &
                                          target_hist_names_soil(:)
                                          !< Arrays to hold target field names
 character(50), allocatable, public    :: target_diag_units(:), &
                                          target_hist_units_2d_cons(:), &
                                          target_hist_units_2d_nstd(:), &
                                          target_hist_units_2d_patch(:), &
                                          target_hist_units_3d_nzp1(:), &
                                          target_hist_units_3d_nz(:), &
                                          target_hist_units_soil(:)
                                          !< Arrays to hold target field units
 character(200), allocatable, public   :: target_diag_longname(:), &
                                          target_hist_longname_2d_cons(:), &
                                          target_hist_longname_2d_nstd(:), &
                                          target_hist_longname_2d_patch(:), &
                                          target_hist_longname_3d_nzp1(:), &
                                          target_hist_longname_3d_nz(:), &
                                          target_hist_longname_soil(:)
                                          !< Arrays to hold target field longname                                         
 type(esmf_fieldbundle), public        :: input_hist_bundle_2d_patch, &
                                          input_hist_bundle_2d_cons, &
                                          input_hist_bundle_2d_nstd, &
                                          input_hist_bundle_3d_nz, &  
                                          input_hist_bundle_3d_nzp1, &
                                          input_hist_bundle_soil
                                          !< bundles to hold input hist fields
 type(esmf_fieldbundle), public        :: target_hist_bundle_2d_patch, &
                                          target_hist_bundle_2d_cons, &
                                          target_hist_bundle_2d_nstd, &
                                          target_hist_bundle_3d_nz, &  
                                          target_hist_bundle_3d_nzp1, &
                                          target_hist_bundle_soil
                                          !< bundles to hold target hist fields

 public :: define_target_grid
 public :: define_input_grid
 public :: cleanup_input_target_grid_data

 contains

!> Define input grid object for MPAS input data.
!!
!! @param [in] localpet ESMF local persistent execution thread 
!! @param [in] npets  Number of persistent execution threads
!! @author  Larissa Reames CIWRO/NOAA/NSSL/FRDD

 subroutine define_input_grid(localpet,npets)

 use mpi
 use netcdf
 use program_setup, only       : grid_file_input_grid
 implicit none

 character(len=500)           :: the_file, dimname

 integer, intent(in)          :: localpet, npets

 integer                      :: error, i, j, k, rc, n, lmi(1), lma(1)

 integer                               :: ncid,id_var, id_dim, dimsize, nVertThis
 integer                               :: nCells, nVertices, maxEdges, dimids(2)
 integer(esmf_kind_i8)                 :: cell_start, cell_end, temp(1)
 integer, allocatable                  :: elemTypes2(:), vertOnCell(:,:), &
                                          nodesPET(:), nodeIDs_temp(:), &
                                          elementConn_temp(:), elementConn(:)                             
 real(esmf_kind_r8), allocatable       :: latCell(:), lonCell(:), &
                                          latVert(:), lonVert(:), &
                                          nodeCoords(:), &
                                          nodeCoords_temp(:), &
                                          elemCoords(:), dummy(:)
 real(esmf_kind_r8), pointer           :: data_1d(:), data_1d2(:)
 real(esmf_kind_r8), parameter         :: PI=4.D0*DATAN(1.D0)


 the_file = grid_file_input_grid


 if (localpet==0) print*,'- OPEN MPAS INPUT FILE: ',trim(the_file)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 if (error /=0) call error_handler("OPENING MPAS INPUT FILE",error)

 if (localpet==0) print*,'- READ nCells'
 error = nf90_inq_dimid(ncid,'nCells', id_dim)
 call netcdf_err(error, 'reading nCells id')
 
 error=nf90_inquire_dimension(ncid,id_dim,len=nCells)
 call netcdf_err(error, 'reading nCells')
 
 nCells_input = nCells
 
  if (localpet==0) print*,'- READ nVertices'
 error = nf90_inq_dimid(ncid,'nVertices',id_dim)
 call netcdf_err(error, 'reading nVertices id')
 
 error=nf90_inquire_dimension(ncid,id_dim,len=nVertices)
 call netcdf_err(error, 'reading nVertices')
 
 nVert_input = nVertices
 
  if (localpet==0) print*,'- READ nVertLevels'
 error = nf90_inq_dimid(ncid,'nVertLevels',id_dim)
 call netcdf_err(error, 'reading nVertLevels id')
 
 error=nf90_inquire_dimension(ncid,id_dim,len=nz_input)
 call netcdf_err(error, 'reading nVertLevels')
 
  if (localpet==0) print*,'- READ nVertLevelsP1'
 error = nf90_inq_dimid(ncid,'nVertLevelsP1',id_dim)
 call netcdf_err(error, 'reading nVertLevelsP1 id')
 
 error=nf90_inquire_dimension(ncid,id_dim,len=nzp1_input)
 call netcdf_err(error, 'reading nVertLevelsP1')
 
  if (localpet==0) print*,'- READ maxEdges'
 error = nf90_inq_dimid(ncid,'maxEdges',id_dim)
 call netcdf_err(error, 'reading maxEdges id')
 
 error=nf90_inquire_dimension(ncid,id_dim,len=maxEdges)
 call netcdf_err(error, 'reading maxEdges')
 
  if (localpet==0) print*,'- READ nSoilLevels'
 error = nf90_inq_dimid(ncid,'nSoilLevels',id_dim)
 call netcdf_err(error, 'reading nSoilLevels id')
 
 error=nf90_inquire_dimension(ncid,id_dim,len=nsoil_input)
 call netcdf_err(error, 'reading nSoilLevels')
 
 
 allocate(latCell(nCells))
 allocate(lonCell(nCells))
 
 allocate(latVert(nVertices))
 allocate(lonVert(nVertices))
 
 allocate(vertOnCell(maxEdges,nCells))
 allocate(zs_target_grid(nsoil_input,1))
 
 allocate(dummy(nCells))
 ! GET CELL CENTER LAT/LON
 if (localpet==0) print*,'- READ LONCELL ID'
 error=nf90_inq_varid(ncid, 'lonCell', id_var)
 call netcdf_err(error, 'reading lonCell id')

 if (localpet==0) print*,'- READ LONCELL'
 error=nf90_get_var(ncid, id_var, lonCell)
 call netcdf_err(error, 'reading lonCell')
 
 if (localpet==0) print*,'- READ LATCELL ID'
 error=nf90_inq_varid(ncid, 'latCell', id_var)
 call netcdf_err(error, 'reading latCell id')

 if (localpet==0) print*,'- READ LATCELL'
 error=nf90_get_var(ncid, id_var, latCell)
 call netcdf_err(error, 'reading latCell')
 

 ! GET VERTEX LAT/LON
 if (localpet==0) print*,'- READ LONVERTEX ID'
 error=nf90_inq_varid(ncid, 'lonVertex', id_var)
 call netcdf_err(error, 'reading lonVertex id')

 if (localpet==0) print*,'- READ LONVERTEX'
 error=nf90_get_var(ncid, id_var, lonVert)
 call netcdf_err(error, 'reading lonVertex')
 
 if (localpet==0) print*,'- READ LATVERTEX ID'
 error=nf90_inq_varid(ncid, 'latVertex', id_var)
 call netcdf_err(error, 'reading latVertex id')

 if (localpet==0) print*,'- READ LATVERTEX'
 error=nf90_get_var(ncid, id_var, latVert)
 call netcdf_err(error, 'reading latVertex')
 
  ! SOIL CENTER DEPTHS
 if (localpet==0) print*,'- READ ZS ID'
 error=nf90_inq_varid(ncid, 'zs', id_var)
 call netcdf_err(error, 'reading zs id')

 if (localpet==0) print*,'- READ ZS'
 error=nf90_get_var(ncid, id_var, zs_target_grid)
 call netcdf_err(error, 'reading ZS')

 if (localpet==0) print*,'- READ HGT'
 error=nf90_inq_varid(ncid, 'ter', id_var)
 call netcdf_err(error, 'reading ter id')
 error=nf90_get_var(ncid, id_var, dummy)
 call netcdf_err(error, 'reading ter')

 if (localpet==0) print*,"- NUMBER OF CELLS ON INPUT GRID ", nCells_input
 if (localpet==0) print*,"- NUMBER OF NODES ON INPUT GRID ", nVert_input

!-----------------------------------------------------------------------
! Create ESMF grid object for the model grid.
!-----------------------------------------------------------------------
nCellsPerPET = ceiling(real(nCells)/real(npets))

 if (localpet==0) print*,'- READ verticesOnCell ID'
 error=nf90_inq_varid(ncid, 'verticesOnCell', id_var)
 call netcdf_err(error, 'reading verticesOnCell id')

 if (localpet==0) print*,'- READ verticesOnCell'
 error=nf90_get_var(ncid, id_var, vertOnCell)
 call netcdf_err(error, 'reading verticesOnCell')

 error = nf90_close(ncid)

 cell_start = localpet*nCellsPerPET+1
 cell_end = min(localpet*nCellsPerPET+nCellsPerPET,nCells)
 nCellsPerPET = cell_end - cell_start + 1

 ! Allocate and fill element corner coordinate array. 
 allocate(nodeCoords_temp(2*maxEdges*nCellsPerPET))
 allocate(nodeIDs_temp(maxEdges*nCellsPerPET))
 allocate(elementConn_temp(maxEdges*nCellsPerPET))
 allocate(elemCoords(2*nCellsPerPET))
 allocate(elemTypes2(nCellsPerPET))
 allocate(elemIDs(nCellsPerPET))
 
 nVertThis = 0
 k = 0
 do i = cell_start,cell_end
    j = i - cell_start + 1
    elemIDs(j) = i
    elemTypes2(j) = 0
    elemCoords(2*j-1) = lonCell(i)*180.0_esmf_kind_r8/PI
    if (elemCoords(2*j-1) > 180.0_esmf_kind_r8) then
        elemCoords(2*j-1) = elemCoords(2*j-1) - 360.0_esmf_kind_r8
    endif
    elemCoords(2*j) = latCell(i)*180.0_esmf_kind_r8/PI
    do n = 1,maxEdges
        if (vertOnCell(n,i)>0) then 
            nVertThis = nVertThis + 1
            
            ! Make sure we don't duplicate nodeIDs or nodeCoords on any PET
            if (.not. any(nodeIDs_temp(1:k)==vertOnCell(n,i))) then
                k = k+1
                nodeCoords_temp(2*k-1) = &
                    lonVert(vertOnCell(n,i)) *180.0_esmf_kind_r8/PI
                if (nodeCoords_temp(2*k-1) > 180.0_esmf_kind_r8) then
                    nodeCoords_temp(2*k-1)  =  &
                        nodeCoords_temp(2*k-1) - 360.0_esmf_kind_r8
                endif
                nodeCoords_temp(2*k) = & 
                    latVert(vertOnCell(n,i))*180.0_esmf_kind_r8/PI

                nodeIDs_temp(k) = vertOnCell(n,i)
            endif
            
            elemTypes2(j) = elemTypes2(j) + 1
            
            !This will have duplicates by defhistion
            temp = FINDLOC(nodeIDS_temp, vertOnCell(n,i))
            elementConn_temp(nVertThis) = temp(1)
            
            
        endif
        
    enddo
 enddo
 allocate(nodeCoords(2*k), nodeIDs(k), elementConn(nVertThis))
 nodeCoords = nodeCoords_temp(1:k*2)
 nodeIDs = nodeIDs_temp(1:k)
 elementConn = elementConn_temp(1:nVertThis)
 
 input_grid = ESMF_MeshCreate(parametricDim=2, & 
                     spatialDim=2, &
                     nodeIDs= nodeIDs, &
                     nodeCoords = nodeCoords, &
                     elementIDs = elemIDs, &
                     elementTypes=elemTypes2, & 
                     elementConn = elementConn, &
                     elementCoords=elemCoords, &
                     coordSys=ESMF_COORDSYS_SPH_DEG, & 
                     rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN MeshCreate", rc)
                     
    
 ! After the creation we are through with the arrays, so they may be deallocated.
 deallocate(elemCoords,elemTypes2)
 deallocate(nodeCoords_temp, nodeCoords)
 deallocate(elementConn_temp, elementConn)
 deallocate(nodeIDs_temp)


 !-----------------------------------------------------------------------
 ! Create lat/lon arrays on input grid
 !-----------------------------------------------------------------------
   
 if (localpet==0) print*,"- CALL FieldCreate FOR INPUT GRID CELL LATITUDE."
 cell_latitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   meshloc=ESMF_MESHLOC_ELEMENT, &
                                   name="input_grid_cell_latitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*,"- CALL FieldCreate FOR INPUT GRID CELL LONGITUDE."
 cell_longitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   meshloc=ESMF_MESHLOC_ELEMENT, &
                                   name="input_grid_cell_longitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN FieldCreate", error)
    
 call ESMF_FieldGet(cell_latitude_input_grid, farrayPtr=data_1d,rc = rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN FieldGet", error)

 call ESMF_FieldGet(cell_longitude_input_grid,farrayPtr=data_1d2, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN FieldGet", error)

 do i = 1, nCellsPerPET
    data_1d(i) = latVert(elemIDs(i))
    data_1d2(i) = lonVert(elemIDs(i))
 enddo

 nullify(data_1d,data_1d2)

 if (localpet==0) print*,"- CALL FieldCreate FOR INPUT GRID HGT."
 hgt_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   meshloc=ESMF_MESHLOC_ELEMENT, &
                                   name="input_grid_hgt", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))&
    call error_handler("IN FieldCreate", error)

 call ESMF_FieldGet(hgt_input_grid, farrayPtr=data_1d,rc = rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))&
    call error_handler("IN FieldGet", error)

 do i = 1, nCellsPerPET
    data_1d(i) = dummy(elemIDs(i))
 enddo
 nullify(data_1d)
 deallocate(lonCell, latCell, lonVert, latVert,vertOnCell,dummy)


 end subroutine define_input_grid

!> Setup the esmf grid object for the target grid.
!!
!! @param [in] localpet ESMF local persistent execution thread 
!! @param [in] npets Number of persistent execution threads
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD   
 subroutine define_target_grid(localpet, npets)

 use netcdf
 use program_setup, only       : file_target_grid

 implicit none

 character(len=500)           :: the_file

 integer, intent(in)          :: localpet, npets

 integer                      :: error, extra, i, j, clb(2), cub(2) 


 real(esmf_kind_r8), allocatable       :: latitude(:,:), longitude(:,:), &
                                          dum2d(:,:)
 integer                               :: ncid,id_var, id_dim 
 real(esmf_kind_r8), pointer           :: lat_src_ptr(:,:), lon_src_ptr(:,:)


 the_file = file_target_grid
    
 if (localpet==0) print*,'- OPEN WRF INPUT FILE: ',trim(the_file)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 if (error /=0) call error_handler("OPENING WRF INPUT FILE",error)

 if (localpet==0) print*,'- READ WEST_EAST ID'
 error=nf90_inq_dimid(ncid, 'west_east', id_dim)
 call netcdf_err(error, 'reading west_east id')

 if (localpet==0) print*,'- READ WEST_EAST'
 error=nf90_inquire_dimension(ncid,id_dim,len=i_target)
 call netcdf_err(error, 'reading west_east')

 if (localpet==0) print*,'- READ SOUTH_NORTH ID'
 error=nf90_inq_dimid(ncid, 'south_north', id_dim)
 call netcdf_err(error, 'reading south_north id')

 if (localpet==0) print*,'- READ SOUTH_NORTH'
 error=nf90_inquire_dimension(ncid,id_dim,len=j_target)
 call netcdf_err(error, 'reading south_north')

 ip1_target = i_target + 1
 jp1_target = j_target + 1 
 allocate(latitude(i_target,j_target))
 allocate(longitude(i_target,j_target))
 allocate(latitude_u(ip1_target,j_target))
 allocate(longitude_u(ip1_target,j_target))
 allocate(latitude_v(i_target,jp1_target))
 allocate(longitude_v(i_target,jp1_target))
 allocate(dum2d(i_target,j_target))
 
 if (localpet==0) print*,'- READ LONGITUDE ID'
 error=nf90_inq_varid(ncid, 'XLONG', id_var)
 if (error .ne. 0) then
   error=nf90_inq_varid(ncid, 'XLONG_M', id_var)
   call netcdf_err(error, 'reading longitude id')
 endif

 if (localpet==0) print*,'- READ LONGITUDE'
 error=nf90_get_var(ncid, id_var, longitude)
 call netcdf_err(error, 'reading longitude')
 
 if (localpet==0) print*,'- READ LATITUDE ID'
 error=nf90_inq_varid(ncid, 'XLAT', id_var)
 if (error .ne. 0) then
   error=nf90_inq_varid(ncid, 'XLAT_M', id_var)
   call netcdf_err(error, 'reading latitude id')
 endif

 if (localpet==0) print*,'- READ LATITUDE'
 error=nf90_get_var(ncid, id_var, latitude)
 call netcdf_err(error, 'reading latitude')

 if (localpet==0) print*,'- READ LONGITUDE_U ID'
 error=nf90_inq_varid(ncid, 'XLONG_U', id_var)
 call netcdf_err(error, 'reading xlong_u id')

 if (localpet==0) print*,'- READ LONGITUDE_U'
 error=nf90_get_var(ncid, id_var, longitude_u)
 call netcdf_err(error, 'reading xlong_u')

 if (localpet==0) print*,'- READ LATITUDE_U ID'
 error=nf90_inq_varid(ncid, 'XLAT_U', id_var)
 call netcdf_err(error, 'reading xlat_u id')
 

 if (localpet==0) print*,'- READ LATITUDE_U'
 error=nf90_get_var(ncid, id_var, latitude_u)
 call netcdf_err(error, 'reading xlat_u')

 if (localpet==0) print*,'- READ LONGITUDE_V ID'
 error=nf90_inq_varid(ncid, 'XLONG_V', id_var)
 call netcdf_err(error, 'reading xlong_v id')

 if (localpet==0) print*,'- READ LONGITUDE_v'
 error=nf90_get_var(ncid, id_var, longitude_v)
 call netcdf_err(error, 'reading xlong_v')

 if (localpet==0) print*,'- READ LATITUDE_v ID'
 error=nf90_inq_varid(ncid, 'XLAT_V', id_var)
 call netcdf_err(error, 'reading xlat_v id')


 if (localpet==0) print*,'- READ LATITUDE_v'
 error=nf90_get_var(ncid, id_var, latitude_v)
 call netcdf_err(error, 'reading xlat_v')

  if (localpet==0) print*,'- READ HGT ID'
 error=nf90_inq_varid(ncid, 'HGT', id_var)
 if (error .ne. 0) then
   error=nf90_inq_varid(ncid, 'HGT_M', id_var)
   call netcdf_err(error, 'reading hgt id')
 endif

 if (localpet==0) print*,'- READ HGT'
 error=nf90_get_var(ncid, id_var, dum2d)
 call netcdf_err(error, 'reading hgt')
    
 if (localpet==0) print*,'- READ GLOBAL ATTRIBUTE DX'
 error = nf90_get_att(ncid,NF90_GLOBAL,'DX',dx)
 call netcdf_err(error, 'reading dx')

 if (localpet==0) print*,"- I/J DIMENSIONS OF THE TARGET GRID TILES ", i_target, j_target

 ip1_target = i_target + 1
 jp1_target = j_target + 1

 if (localpet==0) print*, '- READING GLOBAL ATTRIBUTES'
 error = nf90_get_att(ncid, NF90_GLOBAL, 'CEN_LAT', cen_lat)
 call netcdf_err(error, 'GETTING CEN_LAT GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'CEN_LON', cen_lon)
 call netcdf_err(error, 'GETTING CEN_LON GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'TRUELAT1', truelat1)
 call netcdf_err(error, 'GETTING TRUELAT1 GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'TRUELAT2', truelat2)
 call netcdf_err(error, 'GETTING TRUELAT2 GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'MOAD_CEN_LAT', cen_lat)
 call netcdf_err(error, 'GETTING MOAD_CEN_LAT GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'STAND_LON', stand_lon)
 call netcdf_err(error, 'GETTING STAND_LON GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'POLE_LAT', pole_lat)
 call netcdf_err(error, 'GETTING POLELAT GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'POLE_LON', pole_lon)
 call netcdf_err(error, 'GETTING POLE_LON GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'MAP_PROJ', map_proj)
 call netcdf_err(error, 'GETTING MAP_PROJ GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'MAP_PROJ_CHAR', map_proj_char)
 if (error .ne. 0) then
   if (map_proj == 1) then 
     map_proj_char = "Lambert Conformal"
   else
     map_proj_char = "Lat/Lon"
   endif
 endif
 
 error = nf90_close(ncid)

!-----------------------------------------------------------------------
! Create ESMF grid object for the model grid.
!-----------------------------------------------------------------------

 if (localpet==0) print*,"- CALL GridCreateNoPeriDim FOR TARGET MODEL GRID"
 target_grid = ESMF_GridCreateNoPeriDim(maxIndex=(/i_target,j_target/), & 
                                       indexflag=ESMF_INDEX_GLOBAL, &
                                       rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN GridCreateNoPeriDim", error)


!-----------------------------------------------------------------------
! Read the mask and lat/lons.
!-----------------------------------------------------------------------

 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE."
 latitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_latitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE."
 longitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_longitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID HGT."
 hgt_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_hgt", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))&
    call error_handler("IN FieldCreate", error)
    
 if (localpet==0) print*,"- CALL FieldScatter FOR TARGET GRID LATITUDE. "
 call ESMF_FieldScatter(latitude_target_grid, real(latitude,esmf_kind_r8), rootpet=0, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
   call error_handler("IN FieldScatter", error)
   
 if (localpet==0) print*,"- CALL FieldScatter FOR TARGET GRID LONGITUDE."
 call ESMF_FieldScatter(longitude_target_grid, real(longitude,esmf_kind_r8), rootpet=0, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
   call error_handler("IN FieldScatter", error)

  if (localpet==0) print*,"- CALL FieldScatter FOR TARGET GRID HGT."
 call ESMF_FieldScatter(hgt_target_grid, real(dum2d,esmf_kind_r8), rootpet=0, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__))&
   call error_handler("IN FieldScatter", error)
   
 if (localpet==0) print*,"- CALL GridAddCoord FOR INPUT GRID."
 call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN GridAddCoord", error)
   
 if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT GRID X-COORD."
   nullify(lon_src_ptr)
   call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_CENTER, &
                          coordDim=1, &
                          farrayPtr=lon_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call error_handler("IN GridGetCoord", error)

   if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
   nullify(lat_src_ptr)
   call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_CENTER, &
                          coordDim=2, &
                          computationalLBound=clb, &
                          computationalUBound=cub, &
                          farrayPtr=lat_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call error_handler("IN GridGetCoord", error)
   
    do j = clb(2),cub(2)
      do i = clb(1), cub(1)
        lon_src_ptr(i,j)=real(longitude(i,j),esmf_kind_r8)
        lat_src_ptr(i,j)=real(latitude(i,j),esmf_kind_r8)
      enddo
    enddo
                
  nullify(lon_src_ptr)
  nullify(lat_src_ptr)

 if (localpet==0) print*,"- CALL GridAddCoord FOR TARGET GRID."
 call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_CORNER, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", error)

 if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT GRID X-COORD."

 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=1, &
                        farrayPtr=lon_src_ptr, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord", error)

 if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."

 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_src_ptr, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord", error)
    
 call get_cell_corners(latitude, longitude, lat_src_ptr, lon_src_ptr, dx, clb, cub)

 nullify(lon_src_ptr)
 nullify(lat_src_ptr)
 deallocate(longitude)
 deallocate(latitude)
 deallocate(dum2d)

 end subroutine define_target_grid

 !> For grids with equal cell sizes (e.g., lambert conformal), get lat and on of the grid 
!! cell corners 
!!
!! @param [in]  latitude  2d array of grid latitude
!! @param [in]  longitude 2d array of grid longitude
!! @param [in]  dx grid cell size in meters
!! @param [in]  clb  computational lower bound as returned from ESMF_GridGetCoord
!! @param [in]  cub  computational upper bound as returned from ESMF_GridGetCoord
!! @param [out] latitude_sw  returned latitude at sw corner of grid cells
!! @param [out]elongitude_sw  returned longitude at sw corner of grid cells
!! @author Larissa Reames CIWRO/NSSL/FRDD
 
  subroutine get_cell_corners( latitude, longitude, latitude_sw, longitude_sw, dx,clb,cub)
  implicit none

  real(esmf_kind_r8), intent(in)    :: latitude(i_target,j_target)
  real(esmf_kind_r8), intent(inout), pointer   :: latitude_sw(:,:)
  real(esmf_kind_r8), intent(in)    :: longitude(i_target, j_target)
  real(esmf_kind_r8), intent(inout), pointer   :: longitude_sw(:,:)
  real(esmf_kind_r8), intent(in)    :: dx !grid cell side size (m)

  integer, intent(in) :: clb(2), cub(2)

  real(esmf_kind_r8)                :: lat1, lon1, lat2, lon2, d, brng


  real(esmf_kind_r8), parameter    :: pi = 3.14159265359
  real(esmf_kind_r8), parameter    :: R =  6370000.0
  real(esmf_kind_r8), parameter    :: bearingInDegrees = 135.0

  integer                           :: i, j, e

  d = sqrt((dx**2.0_esmf_kind_r8)/2.0_esmf_kind_r8)

  do j = clb(2),cub(2)
   do i = clb(1), cub(1)
         if (j == jp1_target .and. i == ip1_target) then
       lat1 = latitude(i_target,j_target)  * ( pi / 180.0_esmf_kind_r8 )
       lon1 = longitude(i_target,j_target) * ( pi / 180.0_esmf_kind_r8 )
             brng = 315.0_esmf_kind_r8 * pi / 180.0_esmf_kind_r8
             lat2 = asin( sin( lat1 ) * cos( d / R ) + cos( lat1 ) * sin( d / R ) * cos( brng ) );
             lon2= lon1 + atan2( sin( brng ) * sin( d / R ) * cos( lat1 ), cos( d / R ) - sin( lat1 ) * sin( lat2 ) );
             latitude_sw(ip1_target,jp1_target) = lat2 * 180.0_esmf_kind_r8 / pi
             longitude_sw(ip1_target,jp1_target) = lon2 * 180.0_esmf_kind_r8 / pi
             cycle
         endif
         
     if (i == ip1_target) then
       brng = 225.0_esmf_kind_r8 * pi / 180.0_esmf_kind_r8
       lat1 = latitude(i_target,j)  * ( pi / 180.0_esmf_kind_r8 )
       lon1 = longitude(i_target,j) * ( pi / 180.0_esmf_kind_r8 )
       lat2 = asin( sin( lat1 ) * cos( d / R ) + cos( lat1 ) * sin( d / R ) * cos( brng ) );
       lon2= lon1 + atan2( sin( brng ) * sin( d / R ) * cos( lat1 ), cos( d / R ) - sin( lat1 ) * sin( lat2 ) );
       latitude_sw(ip1_target,j) = lat2 * 180.0_esmf_kind_r8 / pi
       longitude_sw(ip1_target,j) = lon2 * 180.0_esmf_kind_r8 / pi
       cycle
     endif

     if (j == jp1_target) then
       brng = 45.0_esmf_kind_r8 * pi / 180.0_esmf_kind_r8
       lat1 = latitude(i,j_target)  * ( pi / 180.0_esmf_kind_r8 )
       lon1 = longitude(i,j_target) * ( pi / 180.0_esmf_kind_r8 )
       lat2 = asin( sin( lat1 ) * cos( d / R ) + cos( lat1 ) * sin( d / R ) * cos( brng ) );
       lon2= lon1 + atan2( sin( brng ) * sin( d / R ) * cos( lat1 ), cos( d / R ) - sin( lat1 ) * sin( lat2 ) );
       latitude_sw(i,jp1_target) = lat2 * 180.0_esmf_kind_r8 / pi
       longitude_sw(i,jp1_target) = lon2 * 180.0_esmf_kind_r8 / pi
       cycle
     endif

     lat1 = latitude(i,j)  * ( pi / 180.0_esmf_kind_r8 )
     lon1 = longitude(i,j) * ( pi / 180.0_esmf_kind_r8 )
     
     brng = bearingInDegrees * ( pi / 180.0_esmf_kind_r8 );
     lat2 = asin( sin( lat1 ) * cos( d / R ) + cos( lat1 ) * sin( d / R ) * cos( brng ) );
     lon2= lon1 + atan2( sin( brng ) * sin( d / R ) * cos( lat1 ), cos( d / R ) - sin( lat1 ) * sin( lat2 ) );

     latitude_sw(i,j) = lat2 * 180.0_esmf_kind_r8 / pi
     longitude_sw(i,j) = lon2 * 180.0_esmf_kind_r8 / pi
     
   enddo
 enddo

 end subroutine get_cell_corners
 
 subroutine cleanup_input_target_grid_data(localpet)
 
 use program_setup, only    : interp_diag, interp_hist
 

 implicit none
 integer, intent(in)              :: localpet
 integer                          :: rc, i
 type(esmf_field), allocatable    :: fields(:)

 if (localpet==0) print*,"- DESTROY MODEL DATA."
 
 call ESMF_FieldDestroy(node_latitude_input_grid,rc=rc)
 call ESMF_FieldDestroy(node_longitude_input_grid,rc=rc)
 
 call ESMF_FieldDestroy(cell_latitude_input_grid,rc=rc)
 call ESMF_FieldDestroy(cell_longitude_input_grid,rc=rc)
 
 call ESMF_FieldDestroy(latitude_target_grid, rc=rc)
 call ESMF_FieldDestroy(longitude_target_grid, rc=rc)

 if (interp_diag) then
    allocate(fields(n_diag_fields))
    call ESMF_FieldBundleGet(input_diag_bundle, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc) 
    do i = 1, n_diag_fields
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo
    
    call ESMF_FieldBundleDestroy(input_diag_bundle)
    
    call ESMF_FieldBundleGet(target_diag_bundle, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc) 
    do i = 1, n_diag_fields
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo
    
    call ESMF_FieldBundleDestroy(target_diag_bundle)
    deallocate(fields)
 endif
 
 if (n_hist_fields_2d_cons>0) then
    allocate(fields(n_hist_fields_2d_cons))
    call ESMF_FieldBundleGet(input_hist_bundle_2d_cons, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc) 
    do i = 1, n_hist_fields_2d_cons
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo
    
    call ESMF_FieldBundleDestroy(input_hist_bundle_2d_cons)
    
    call ESMF_FieldBundleGet(target_hist_bundle_2d_cons, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc) 
    do i = 1, n_hist_fields_2d_cons
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo
    
    call ESMF_FieldBundleDestroy(target_hist_bundle_2d_cons)
    deallocate(fields)
 endif
 
  if (n_hist_fields_2d_nstd>0) then
    allocate(fields(n_hist_fields_2d_nstd))
    call ESMF_FieldBundleGet(input_hist_bundle_2d_nstd, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc) 
    do i = 1, n_hist_fields_2d_nstd
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo
    
    call ESMF_FieldBundleDestroy(input_hist_bundle_2d_nstd)
    
    call ESMF_FieldBundleGet(target_hist_bundle_2d_nstd, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc) 
    do i = 1, n_hist_fields_2d_nstd
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo
    
    call ESMF_FieldBundleDestroy(target_hist_bundle_2d_nstd)
    deallocate(fields)
 endif
 
 if (n_hist_fields_2d_patch>0) then
    allocate(fields(n_hist_fields_2d_patch))
    call ESMF_FieldBundleGet(input_hist_bundle_2d_patch, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc) 
    do i = 1, n_hist_fields_2d_patch
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo
    
    call ESMF_FieldBundleDestroy(input_hist_bundle_2d_patch)
    
    call ESMF_FieldBundleGet(target_hist_bundle_2d_patch, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc) 
    do i = 1, n_hist_fields_2d_patch
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo
    
    call ESMF_FieldBundleDestroy(target_hist_bundle_2d_patch)
    deallocate(fields)
 endif
 
  if (n_hist_fields_3d_nz>0) then
    allocate(fields(n_hist_fields_3d_nz))
    call ESMF_FieldBundleGet(input_hist_bundle_3d_nz, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc) 
    do i = 1, n_hist_fields_3d_nz
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo
    
    call ESMF_FieldBundleDestroy(input_hist_bundle_3d_nz)
    
    call ESMF_FieldBundleGet(target_hist_bundle_3d_nz, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc) 
    do i = 1, n_hist_fields_3d_nz
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo
    
    call ESMF_FieldBundleDestroy(target_hist_bundle_3d_nz)
    deallocate(fields)
 endif
 
 if (n_hist_fields_3d_nzp1>0) then
    allocate(fields(n_hist_fields_3d_nzp1))
    call ESMF_FieldBundleGet(input_hist_bundle_3d_nzp1, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc) 
    do i = 1, n_hist_fields_3d_nzp1
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo
    
    call ESMF_FieldBundleDestroy(input_hist_bundle_3d_nzp1)
    
    call ESMF_FieldBundleGet(target_hist_bundle_3d_nzp1, fieldList=fields, & 
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc) 
    do i = 1, n_hist_fields_3d_nzp1
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo
    
    call ESMF_FieldBundleDestroy(target_hist_bundle_3d_nzp1)
    deallocate(fields)
 endif
    

 call ESMF_MeshDestroy(input_grid, rc=rc)
 
 call ESMF_GridDestroy(target_grid, rc=rc)

 end subroutine cleanup_input_target_grid_data
 
subroutine unique_sort(val,nvals, final)
    implicit none
    integer, intent(in) :: nvals
    integer, intent(in) :: val(nvals)
    integer, intent(inout), allocatable :: final(:)
    integer :: i = 0, min_val, max_val
    integer, dimension(nvals) :: unique

    min_val = minval(val)-1
    max_val = maxval(val)
    do while (min_val<max_val)
        i = i+1
        min_val = minval(val, mask=val>min_val)
        unique(i) = min_val
    enddo
    allocate(final(i), source=unique(1:i))   !<-- Or, just use unique(1:i) 
end subroutine unique_sort

 end module model_grid
