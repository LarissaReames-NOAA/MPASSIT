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


 type(esmf_field),  public              :: latitude_target_grid
                                           !< latitude of grid center, target grid
 type(esmf_field),  public              :: longitude_target_grid
                                           !< longitude of grid center, target grid
                                           
 integer, parameter, public            :: n_diag_fields=4
 type(esmf_fieldbundle), public        :: input_diag_bundle     
 type(esmf_fieldbundle), public        :: target_diag_bundle                                      
 

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
 use program_setup, only       : init_file_input_grid
 implicit none

 character(len=500)           :: the_file, dimname

 integer, intent(in)          :: localpet, npets

 integer                      :: error, i, j, k, rc, n, lmi(1), lma(1)

 integer                               :: ncid,id_var, id_dim, dimsize, nVertThis
 integer                               :: nCells, nVertices, maxEdges, dimids(2)
 integer(esmf_kind_i8)                 :: cell_start, cell_end, temp(1)
 integer, allocatable				   :: elemTypes2(:), vertOnCell(:,:), &
 										  nodesPET(:), nodeIDs_temp(:), &
 										  elementConn_temp(:), elementConn(:)							  
 real(esmf_kind_r8), allocatable       :: latCell(:), lonCell(:), &
 										  latVert(:), lonVert(:), &
 										  nodeCoords(:), &
 										  nodeCoords_temp(:), &
 										  elemCoords(:)
 real(esmf_kind_r8), pointer           :: data_1d(:), data_1d2(:)
 real(esmf_kind_r8), parameter         :: PI=4.D0*DATAN(1.D0)


 the_file = init_file_input_grid


 print*,'- OPEN MPAS INPUT FILE: ',trim(the_file)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 if (error /=0) call error_handler("OPENING MPAS INPUT FILE",error)

 print*,'- READ nCells'
 error = nf90_inq_dimid(ncid,'nCells', id_dim)
 call netcdf_err(error, 'reading nCells id')
 
 error=nf90_inquire_dimension(ncid,id_dim,len=nCells)
 call netcdf_err(error, 'reading nCells')
 
 nCells_input = nCells
 
 print*,'- READ nVertices'
 error = nf90_inq_dimid(ncid,'nVertices',id_dim)
 call netcdf_err(error, 'reading nVertices id')
 
error=nf90_inquire_dimension(ncid,id_dim,len=nVertices)
 call netcdf_err(error, 'reading nVertices')
 
 nVert_input = nVertices
 
 print*,'- READ nVertLevels'
 error = nf90_inq_dimid(ncid,'nVertLevels',id_dim)
 call netcdf_err(error, 'reading nVertLevels id')
 
 error=nf90_inquire_dimension(ncid,id_dim,len=nz_input)
 call netcdf_err(error, 'reading nVertLevels')
 
  print*,'- READ nVertLevelsP1'
 error = nf90_inq_dimid(ncid,'nVertLevelsP1',id_dim)
 call netcdf_err(error, 'reading nVertLevelsP1 id')
 
 error=nf90_inquire_dimension(ncid,id_dim,len=nzp1_input)
 call netcdf_err(error, 'reading nVertLevelsP1')
 
 print*,'- READ maxEdges'
 error = nf90_inq_dimid(ncid,'maxEdges',id_dim)
 call netcdf_err(error, 'reading maxEdges id')
 
 error=nf90_inquire_dimension(ncid,id_dim,len=maxEdges)
 call netcdf_err(error, 'reading maxEdges')

 
 allocate(latCell(nCells))
 allocate(lonCell(nCells))
 
 allocate(latVert(nVertices))
 allocate(lonVert(nVertices))
 
 allocate(vertOnCell(maxEdges,nCells))
 
 
 ! GET CELL CENTER LAT/LON
 print*,'- READ LONCELL ID'
 error=nf90_inq_varid(ncid, 'lonCell', id_var)
 call netcdf_err(error, 'reading lonCell id')

 print*,'- READ LONCELL'
 error=nf90_get_var(ncid, id_var, lonCell)
 call netcdf_err(error, 'reading lonCell')
 
 print*,'- READ LATCELL ID'
 error=nf90_inq_varid(ncid, 'latCell', id_var)
 call netcdf_err(error, 'reading latCell id')

 print*,'- READ LATCELL'
 error=nf90_get_var(ncid, id_var, latCell)
 call netcdf_err(error, 'reading latCell')
 

 ! GET VERTEX LAT/LON
 print*,'- READ LONVERTEX ID'
 error=nf90_inq_varid(ncid, 'lonVertex', id_var)
 call netcdf_err(error, 'reading lonVertex id')

 print*,'- READ LONVERTEX'
 error=nf90_get_var(ncid, id_var, lonVert)
 call netcdf_err(error, 'reading lonVertex')
 
 print*,'- READ LATVERTEX ID'
 error=nf90_inq_varid(ncid, 'latVertex', id_var)
 call netcdf_err(error, 'reading latVertex id')

 print*,'- READ LATVERTEX'
 error=nf90_get_var(ncid, id_var, latVert)
 call netcdf_err(error, 'reading latVertex')


 print*,"- NUMBER OF CELLS ON INPUT GRID ", nCells_input
 print*,"- NUMBER OF NODES ON INPUT GRID ", nVert_input

!-----------------------------------------------------------------------
! Create ESMF grid object for the model grid.
!-----------------------------------------------------------------------
nCellsPerPET = ceiling(real(nCells)/real(npets))

 print*,'- READ verticesOnCell ID'
 error=nf90_inq_varid(ncid, 'verticesOnCell', id_var)
 call netcdf_err(error, 'reading verticesOnCell id')

 print*,'- READ verticesOnCell'
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
 			
 			!This will have duplicates by definition
 			temp = FINDLOC(nodeIDS_temp, vertOnCell(n,i))
 			elementConn_temp(nVertThis) = temp(1)
 			
 			
		endif
		
	enddo
 enddo
 allocate(nodeCoords(2*k), nodeIDs(k), elementConn(nVertThis))
 nodeCoords = nodeCoords_temp(1:k*2)
 nodeIDs = nodeIDs_temp(1:k)
 elementConn = elementConn_temp(1:nVertThis)
 print*, elementConn(1:elemTypes2(1)), nodeCoords(1:elemTypes2(1)*2)
 
 
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
   
 print*,"- CALL FieldCreate FOR INPUT GRID CELL LATITUDE."
 cell_latitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   meshloc=ESMF_MESHLOC_ELEMENT, &
                                   name="input_grid_cell_latitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR INPUT GRID CELL LONGITUDE."
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

 print*, localpet, maxval(elemIDs), minval(elemIDs)
 do i = 1, nCellsPerPET
 	data_1d(i) = latVert(elemIDs(i))
    data_1d2(i) = lonVert(elemIDs(i))
 enddo

 nullify(data_1d,data_1d2)
 deallocate(lonCell, latCell, lonVert, latVert,vertOnCell)


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


 real(esmf_kind_r8), allocatable       :: latitude(:,:), longitude(:,:)
 integer                               :: ncid,id_var, id_dim 
 real(esmf_kind_r8), pointer           :: lat_src_ptr(:,:), lon_src_ptr(:,:)
 real(esmf_kind_r8)                    :: dx


 the_file = file_target_grid
    
 print*,'- OPEN WRF INPUT FILE: ',trim(the_file)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 if (error /=0) call error_handler("OPENING WRF INPUT FILE",error)

 print*,'- READ WEST_EAST ID'
 error=nf90_inq_dimid(ncid, 'west_east', id_dim)
 call netcdf_err(error, 'reading west_east id')

 print*,'- READ WEST_EAST'
 error=nf90_inquire_dimension(ncid,id_dim,len=i_target)
 call netcdf_err(error, 'reading west_east')

 print*,'- READ SOUTH_NORTH ID'
 error=nf90_inq_dimid(ncid, 'south_north', id_dim)
 call netcdf_err(error, 'reading south_north id')

 print*,'- READ SOUTH_NORTH'
 error=nf90_inquire_dimension(ncid,id_dim,len=j_target)
 call netcdf_err(error, 'reading south_north')
 
 allocate(latitude(i_target,j_target))
 allocate(longitude(i_target,j_target))
 
 print*,'- READ LONGITUDE ID'
 error=nf90_inq_varid(ncid, 'XLONG', id_var)
 call netcdf_err(error, 'reading longitude id')

 print*,'- READ LONGITUDE'
 error=nf90_get_var(ncid, id_var, longitude)
 call netcdf_err(error, 'reading longitude')
 
 print*,'- READ LATITUDE ID'
 error=nf90_inq_varid(ncid, 'XLAT', id_var)
 call netcdf_err(error, 'reading latitude id')

 print*,'- READ LATITUDE'
 error=nf90_get_var(ncid, id_var, latitude)
 call netcdf_err(error, 'reading latitude')
    
 print*,'- READ GLOBAL ATTRIBUTE DX'
 error = nf90_get_att(ncid,NF90_GLOBAL,'DX',dx)
 call netcdf_err(error, 'reading dx')
 error = nf90_close(ncid)
 print*,"- I/J DIMENSIONS OF THE TARGET GRID TILES ", i_target, j_target

 ip1_target = i_target + 1
 jp1_target = j_target + 1

!-----------------------------------------------------------------------
! Create ESMF grid object for the model grid.
!-----------------------------------------------------------------------

 print*,"- CALL GridCreateNoPeriDim FOR TARGET MODEL GRID"
 target_grid = ESMF_GridCreateNoPeriDim(maxIndex=(/i_target,j_target/), & 
                                       indexflag=ESMF_INDEX_GLOBAL, &
                                       rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN GridCreateNoPeriDim", error)


!-----------------------------------------------------------------------
! Read the mask and lat/lons.
!-----------------------------------------------------------------------

 print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE."
 latitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_latitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN FieldCreate", error)

 print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE."
 longitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_longitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN FieldCreate", error)
    
 print*,"- CALL FieldScatter FOR TARGET GRID LATITUDE. "
 call ESMF_FieldScatter(latitude_target_grid, real(latitude,esmf_kind_r8), rootpet=0, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
   call error_handler("IN FieldScatter", error)
   
 print*,"- CALL FieldScatter FOR TARGET GRID LONGITUDE."
 call ESMF_FieldScatter(longitude_target_grid, real(longitude,esmf_kind_r8), rootpet=0, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
   call error_handler("IN FieldScatter", error)
   
 print*,"- CALL GridAddCoord FOR INPUT GRID."
 call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
    call error_handler("IN GridAddCoord", error)
   
 print*,"- CALL GridGetCoord FOR INPUT GRID X-COORD."
   nullify(lon_src_ptr)
   call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_CENTER, &
                          coordDim=1, &
                          farrayPtr=lon_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call error_handler("IN GridGetCoord", error)

   print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
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
  print*, localpet, minval(lon_src_ptr), maxval(lon_src_ptr), &
 			minval(lat_src_ptr), maxval(lat_src_ptr)			
  nullify(lon_src_ptr)
  nullify(lat_src_ptr)

 deallocate(longitude)
 deallocate(latitude)


 end subroutine define_target_grid
 
 subroutine cleanup_input_target_grid_data
 
 use program_setup, only    : data_to_interp

 implicit none

 integer                          :: rc, i
 type(esmf_field), allocatable    :: fields(:)

 print*,"- DESTROY MODEL DATA."
 
 call ESMF_FieldDestroy(node_latitude_input_grid,rc=rc)
 call ESMF_FieldDestroy(node_longitude_input_grid,rc=rc)
 
 call ESMF_FieldDestroy(cell_latitude_input_grid,rc=rc)
 call ESMF_FieldDestroy(cell_longitude_input_grid,rc=rc)
 
 call ESMF_FieldDestroy(latitude_target_grid, rc=rc)
 call ESMF_FieldDestroy(longitude_target_grid, rc=rc)

 if (data_to_interp == 'diag') then
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
