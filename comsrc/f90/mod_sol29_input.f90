module mod_sol29_input
  implicit none
  
contains 

  subroutine sol29_initialize_unstructured_input
    use mod_solcommon
    implicit none
  
	! Explanations down below.
    vr_gamma_loc     = 0.0
    vr_gamma_scale   = 0.0
    vr_gamma_a       = 0.0
    vr_gamma_c       = 0.0
    timestep         = 0.0
    blob_ne          = 0.0
    blob_te          = 0.0
    tau_ne           = 0.0
    tau_te           = 0.0
    seed_targ_te     = 0.0
    blob_radius      = 0.0
    nblobs           = 0
    niterations      = 0
    runeir29         = 0
    load_divimp      = 0
    load_divimp_path = ''
    blob_length      = 0.0
    blob_freq        = 0.0
    vr_offset        = 0.0
    frac_holes       = 0.0
    hole_ne          = -1   ! -1 = default to blob value
    hole_te          = -1
    hole_tau_ne      = -1
    hole_tau_te      = -1
    tau_rad_start    = 0.0
    blob_vr_type     = 0
    vr_gauss_loc     = 0.0
    vr_gauss_scale   = 0.0
  
  end subroutine sol29_initialize_unstructured_input

  subroutine sol29_unstructured_input(tag, line)
    use mod_solcommon
    implicit none
    
    character*(*) :: tag, line
    
    ! Parameters for the generalized gamma function that defines the
	! radial velocity distribution of blobs crossing the separatrix.
	! It is difficult to assign intuitive meaning to each parameter
	! since it is a bit complicated of a distribution, but it is one
	! that seems to work well. The distribution can be seen here:
	! https://docs.scipy.org/doc/scipy/reference/generated/scipy.
	! stats.gengamma.html
    if (tag(1:3).eq.'X01') then
      call readr(line, vr_gamma_loc, 0.0, sol22_machhi, &
        'SOL29: vr_gamma_loc') 
    elseif (tag(1:3).eq.'X02') then
      call readr(line, vr_gamma_scale, 0.0, sol22_machhi, &
        'SOL29: vr_gamma_scale')      
    elseif (tag(1:3).eq.'X03') then
      call readr(line, vr_gamma_a, 0.0, sol22_machhi, &
        'SOL29: vr_gamma_a')   
    elseif (tag(1:3).eq.'X04') then
      call readr(line, vr_gamma_c, 0.0, sol22_machhi, &
        'SOL29: vr_gamma_c')
    
    ! Analogous to qtim, just using a different variable.
    elseif (tag(1:3).eq.'X05') then
      call readr(line, timestep, 0.0, sol22_machhi, &
        'SOL29: timestep')
    
    ! The starting ne, Te of blobs crossing the separatrix. The values
    ! exponentially decay in time according to the characteristic times
    ! specified (every tau_ne seconds the density decreases by 1/e). 
    elseif (tag(1:3).eq.'X06') then
      call readr(line, blob_ne, 0.0, sol22_machhi, &
        'SOL29: blob_ne')   
    elseif (tag(1:3).eq.'X07') then
      call readr(line, blob_te, 0.0, sol22_machhi, &
        'SOL29: blob_te') 
    elseif (tag(1:3).eq.'X08') then
      call readr(line, tau_ne, 0.0, sol22_machhi, &
        'SOL29: tau_ne')  
    elseif (tag(1:3).eq.'X09') then
      call readr(line, tau_te, 0.0, sol22_machhi, &
        'SOL29: tau_te')
        
    ! Starting target Te to calculate the parallel flows on the first 
    ! iteration.
    elseif (tag(1:3).eq.'X10') then
      call readr(line, seed_targ_te, 0.0, sol22_machhi, &
        'SOL29: seed_targ_te')
        
    ! Width of the blobs.
    elseif (tag(1:3).eq.'X11') then
      call readr(line, blob_radius, 0.0, sol22_machhi, &
        'SOL29: blob_radius')
        
    ! How many blobs to launch and how many iterations to run the 
    ! simulation for. The higher the better here, it will only give 
    ! better statistics and won't change the fundamental results.
    elseif (tag(1:3).eq.'X12') then
      call readi(line, nblobs, 0, sol22_machhi, &
        'SOL29: nblobs')  
    elseif (tag(1:3).eq.'X13') then
      call readi(line, niterations, 0, sol22_machhi, &
        'SOL29: niterations')
        
    ! Whether or not to iterate with EIRENE for the contribution to the
    ! plasma density from the neutrals.
    elseif (tag(1:3).eq.'X14') then
      call readi(line, runeir29, 0, 1, &
        'SOL29: run EIRENE07')
        
    ! Load a previous DIVIMP run to get the contribution to the plasma
    ! density from the impurities.
    elseif (tag(1:3).eq.'X15') then
      call readi(line, load_divimp, 0, 1, &
        'SOL29: load previous DIVIMP run')
    elseif (tag(1:3).eq.'X16') then
      call readc(line, load_divimp_path, &
        'SOL29: path the previous DIVIMP run')
        
    ! Length of blobs. Assumed constant.
    elseif (tag(1:3).eq.'X17') then
      call readr(line, blob_length, 0.0, sol22_machhi, &
        'SOL29: blob_length')
    
    ! Measured blob frequency. Assumed constant everywhere.
    elseif (tag(1:3).eq.'X18') then
      call readr(line, blob_freq, 0.0, sol22_machhi, &
        'SOL29: blob_freq')
        
    ! Constant offset of the radial velocity. When a value for the blob
    ! is chosen, this value is then added onto it.
    elseif (tag(1:3).eq.'X19') then
      call readr(line, vr_offset, -sol22_machhi, sol22_machhi, &
        'SOL29: vr_offset')
        
    ! How many holes (inward moving blobs) to follow as a fraction of
    ! nblobs.
    elseif (tag(1:3).eq.'X20') then
      call readr(line, frac_holes, 0.0, sol22_machhi, &
        'SOL29: frac_holes')
    
    ! Similar to the above, just for the holes.
    elseif (tag(1:3).eq.'X21') then
      call readr(line, hole_ne, -1.0, sol22_machhi, &
        'SOL29: hole_ne')   
    elseif (tag(1:3).eq.'X22') then
      call readr(line, hole_te, -1.0, sol22_machhi, &
        'SOL29: hole_te') 
    elseif (tag(1:3).eq.'X23') then
      call readr(line, hole_tau_ne, -1.0, sol22_machhi, &
        'SOL29: hole_tau_ne')  
    elseif (tag(1:3).eq.'X24') then
      call readr(line, hole_tau_te, -1.0, sol22_machhi, &
        'SOL29: hole_tau_te')
        
    ! The starting radial location of the  blobs/holes are generated 
    ! using an exponential that decays away from the separatrix 
    ! characterized by tau_rad_start. Max probability at the separatrix,
    ! and then decays from there outwards. This value is mapped to the 
    ! OMP.
    elseif (tag(1:3).eq.'X25') then
      call readr(line, tau_rad_start, 0.0, sol22_machhi, &
        'SOL29: tau_rad_start')
        
    ! What type of radial velocity distribution to use for the blobs and
    ! holes.
    ! 0 = Gaussian distribution
    ! 1 = Generalized Gamma distribution (see X01-X04). 
    elseif (tag(1:3).eq.'X26') then
      call readi(line, blob_vr_type, 0, 1, &
        'SOL29: blob_vr_type')
        
    ! The characteristics of the Gaussian distribution for the blob/hole
    ! vr distribution (when blob_vr_type = 0). 
    elseif (tag(1:3).eq.'X27') then
      call readr(line, vr_gauss_loc, 0.0, sol22_machhi, &
        'SOL29: vr_gauss_loc')
    elseif (tag(1:3).eq.'X28') then
      call readr(line, vr_gauss_scale, 0.0, sol22_machhi, &
        'SOL29: vr_gauss_loc')
        
        
    ! Assign the blob values to the hole parameters if -1.
    if (hole_ne.lt.0.0) then
      hole_ne = blob_ne
    endif
    if (hole_te.lt.0.0) then
      hole_te = blob_te
    endif
    if (hole_tau_ne.lt.0.0) then
      hole_tau_ne = tau_ne
    endif
    if (hole_tau_te.lt.0.0) then
      hole_tau_te = tau_te
    endif 
     
        
     
        
    endif
  
  end subroutine sol29_unstructured_input
  
end module mod_sol29_input

