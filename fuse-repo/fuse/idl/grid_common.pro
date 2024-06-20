;
; ======================================================================
;
PRO grid_SetupParameters
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state

  param = {  $
    FORWARD       :  1.0D ,  $
    BACKWARD      : -1.0D ,  $
    CORE          :  1    ,  $
    SOL           :  2    ,  $
    SOL_LFS       :  3    ,  $
    SOL_HFS       :  4    ,  $
    PFZ           :  5    ,  $
    PFZ_SECONDARY :  6    ,  $
    HFS           :  1    ,  $
    LFS           :  2    ,  $
;    LOWER_NULL    :  1    ,  $
;    UPPER_NULL    :  2    ,  $

    TRUE          :  1    ,  $
    FALSE         :  0    ,  $

    LSN           :  1    ,  $
    USN           :  2	  ,  $
    CDN           :  3	  ,  $
    LDN           :  4	  ,  $
    UDN           :  5	  ,  $
    LIMITED       :  6	  ,  $

    EXPONENTIAL   :  1    ,  $ 
    TANH          :  2       $ 
    }

  ctr_min_length = MAKE_ARRAY(7,VALUE=0.10D)

  option = {  $

    ctr_break_distance    : 0.10D         ,  $  ; distance threshold between consecutive points on an IDL PSI contour which indicates a break between contour segments
    ctr_2nd_xpt_proximity : 0.10D         ,  $  ; distance between a tangency point and the 2nd x-point below which the slow PSI scan is triggered
    ctr_2nd_xpt_frac      : 0.05D         ,  $  ; fraction of distance between 1st and 2nd separtrices (in PSI) when the slow PSI scan is triggered
    ctr_2nd_xpt_slow      : 0.05D         ,  $  ; fraction of the PSI scanning speed taken as the 2nd x-point is approached
    ctr_2nd_xpt_parallel  : 0.05D         ,  $  ; distance between a tangency point and the 2nd x-point below which a parallel contour boundary marker is used (rather than perpendicular to the wall) 
    ctr_minimum_length    : ctr_min_length,  $  ; minimum length of the contour, which stops outward expansion of the grid (after a certain number of counts)
    ctr_minimum_count     : 3             ,  $  ; number of times the contour length is below the minimum length before the scan is stopped

    rad_res_core_n   : 51        ,  $  ; radial resolution at the separatrix, core
    rad_res_core_min : 0.002D    ,  $  ; 
    rad_res_core_max : 0.100D    ,  $  ; 
    rad_res_core_exp : 2.000D    ,  $  ; 
    rad_res_core_opt : param.TANH,  $  ; 

    rad_res_sol_n   : 101        ,  $  ; radial resolution at the separatrix, core
    rad_res_sol_min : 0.001D     ,  $  ; 
    rad_res_sol_max : 0.020D     ,  $  ; 
    rad_res_sol_exp : 1.000D     ,  $  ; 
    rad_res_sol_opt : param.TANH ,  $  ; 

    pol_res_min           : 0.001D ,  $  ; poloidal resoluiton at the separatrix targets
    pol_res_max           : 0.050D ,  $  ; maxiumum poloidal resolution along a separatrix contour
    pol_res_opt           : 1      ,  $  ; poloidal distribution of cells: 1-TANH
    pol_res_exp           : 1.000D ,  $  ; parameter for distribution

    pol_2nd_roi           : 0.03D  ,  $  ; region of interest for poloidal refinement near the secondary x-point, as a % of poloidal field line length
    pol_2nd_niter         : 2      ,  $  ; number of iterations when applying poloidal refinement at the secondary x-point

    tar_threshold_size    : 0.05D  ,  $  ; maximum length of a target segment
    tar_threshold_angle   : 10.0D  ,  $  ; angle between contour and wall below which the target size is increased (important near tangency)

    ;                        CORE,   SOL, SOL_LFS, SOL_HFS,  PFZ, PFZ_SECONDARY
    ctr_boundary          : [0.0D, 0.0D0,    0.0D,    0.0D, 0.0D,          0.0D],  $ ; PSI boundaries for the grid
    

    debug                 : 0           ,  $
    xrange                : [0.0D,0.0D] ,  $
    yrange                : [0.0D,0.0D] ,  $
    
    dummy : 0  $
    }


  state = {  $

    version      : 1.0,  $

    iteration    :  1 ,  $  current grid generation integration (can be > 1 if refining a malformed grid)
    section_i1   : -1 ,  $  poloidal segment list
    section_i2   : -1 ,  $

    refinement_rad :  0 ,  $  automated radial refinement parameters
    ref_rad_x1     : -1 ,  $
    ref_rad_y1     : -1 ,  $
    ref_rad_x2     : -1 ,  $
    ref_rad_y2     : -1 ,  $
    ref_rad_psi1   : -1 ,  $
    ref_rad_psi2   : -1 ,  $ 

    refinement_pol :  0 ,  $  automated poloidal refinement parameters
    ref_pol_i1     : -1 ,  $
    ref_pol_i2     : -1 ,  $
    ref_pol_pos    : -1 ,  $

    diverted       :  1 ,  $    

    geometry       : -1 ,  $  description of the magnetic geometry, i.e. LSN, CDN, etc. 

    dummy : 0  $
    }


  ; ------------------------------------------------------------------
END
