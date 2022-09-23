This program interpolates fields from an MPAS mesh (the "source" mesh)
to a regular WRF grid (the "target" grid). Integer fields use
nearest neighbor interpolation, snow-related fields use conservative
regridding, and all other fields use bilinear interpolation.

The main point of this program was to be able to interpolate fields from
an MPAS mesh to a regular grid for data visualization and post-processing
purposes. 



Compiling instructions: 

1) Edit build.sh script.  Needed libraries: NETCDF, ESMF, MPI.
  Note: If working on a RDHPCS machine (Hera, Jet, etc.) no edits should
  be necessary 

2) run build.sh


Run instructions:

1) Copy all variable list files from ./parm in to run directory. Each file is formatted like:

rainc                           RAINC
rainnc                          RAINNC
refl10cm_max                    REFL10CM_MAX
[...]

with the variable name in the MPAS file in the first column, and the desired output
name in the second file. Variables are separated in to four categories: 

a) Variables from the diag file (diaglist)

b) 2-dimensional variables from the history file (histlist_2d)

c) 3-dimensional variables from the history file (histlist_3d)

d) 3-dimensional soil variables from the history file (histlist_soil)

The user can edit these to their specifications, if they so wish.

2) Create a namelist file of any name in run directory (e.g., namelist.input)

Input/namelist:

&config
  grid_file_input_grid="/scratch/wof/mpas/init.nc"
  hist_file_input_grid="/scratch/wof/mpas/init.nc"
  diag_file_input_grid="/scratch/wof/mpas/diag.2019-05-18_00.00.00.nc"
  file_target_grid="/scratch/wicker/27April2011/ICs/wrfinput_d01"
  output_file="/scratch/larissa.reames/out_hist_diag.nc"
  interp_diag=.true.
  interp_hist=.true.
  esmf_log=.false.
/

grid_file_input_grid : Full path of MPAS file containing grid information

hist_file_input_grid : Full path of input history MPAS data

diag_file_input_grid : Full path of input diag MPAS data

file_target_grid : Full path of WRF file containing target grid information, 
                   supported types: wrfout, wrfinput, geo_em

output_file : Full path of output file

interp_diag : Whether to interpolate fields from the diag file (T/F)

interp_hist : Whether to interpolate fields from the hist file (T/F)

wrf_mod_vars : Whether to modify output variables to conform to WRF shapes (i.e. staggered winds) (default:.false.)

esmf_log    : Whether to output ESMF files (PET) if an ESMF error is encounters (default:.false.)


3) Submit to job queue (most grids won't need more then 3 nodes/72 processes)

 srun bin/mpassit namelist.input  [or whatever your chosen namelist is named]
