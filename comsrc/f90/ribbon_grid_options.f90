module ribbon_grid_options

  save

  integer :: rg_grid_opt
  integer :: rg_block_av
  integer :: rg_min_cells

  real :: rg_max_r_sep
  real :: rg_max_s_sep

  ! Minimum and maximum S values for building a reduced grid on a subset of
  ! the intersection data.
  real :: rg_int_win_mins,rg_int_win_maxs
  
  ! Boundary points for building an expanded range grid. These are used as the 
  ! corner points of a box which connect to the embedded intersection data for 
  ! grid generation. 
  ! integer :: rg_box_opt 
  real :: rg_minr,rg_maxr,rg_mins,rg_maxs
  
  ! cutoff length limiting range of ring generation
  real :: lcutoff

  ! factor used to control spacing of cells between "surfaces" 
  ! Used as an exponential factor
  integer :: cell_spacing_option
  real :: cell_spacing_factor

  ! Option defining the input data file format
  integer :: ribbon_input_format_opt

  character*256 :: rg_castem_data

end module ribbon_grid_options
