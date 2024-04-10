; ======================================================================
;
; GRID
;
; This program will generate non-orthogonal extended grids for use with 
; OSM-EIRENE-DIVIMP.  Discontinuous targets ("broken targets") of 
; arbitrary complexity can be represented.  
;
; ======================================================================
;
; WARNING: This is pre-alpha code and so is a complete mess -- enter at 
; your own risk.  The noise to your screen is extensive and there is 
; currently no way to turn it off, i.e. removing the /debug keyword
; won't help you very much.
; 
; ======================================================================
;
; Instructions:
; -------------
;
; It is strongly recommended that only lines in the "USER SECTION" are 
; modified.
;
; 1. Compile the code:
;
;    IDL> @grid_make
;
; 1. Set the radial boundaries of the grid:
; 
;    For the appropriate <machine> case block in the user section, edit
;    user_finish and option.ctr_boundary to set the outer radial boundary
;    for the grid.  Check your settings with:
; 
;    IDL> a=grid_Run(/<machine>,shot='<shot>',equ='<.equ file>',/debug,/preview)
; 
;    where <shot> is the name of the directory where the .equ file is 
;    stored.
;
;    You can try to see where the radial boundaries of the grid are from the 
;    coloured contours -- it's confusing now, I admit.  Will clean it up "soon".
;
; 2. Generate the grid:
;
;    IDL> a=grid_Run(/<machine>,shot='<shot>',equ='<.equ file>',/debug,/all,/save,/check)
;
;    When the grid is finished you'll see:
;
;    IDL> FNAME: grid_data/<name of grid>
;
; Example:
; --------
;
; For the C-Mod example case included with FUSE distribution:
;
;   IDL> a=grid_Run(/cmod,shot='1100303017',equ='cmod_1100303017_01380.x16.equ',/debug,/preview)
;   IDL> a=grid_Run(/cmod,shot='1100303017',equ='cmod_1100303017_01380.x16.equ',/debug,/save,/all,/check)
;

;PRO grid_Test
;  print,'hello'
;END

FUNCTION grid_Run, machine  = machine  , $
                   iter     = iter     , $
                   cfetr    = cfetr    , $
                   mast     = mast     , $
                   cmod     = cmod     , $
                   west     = west     , $
                   jet      = jet      , $
                   aug      = aug      , $
                   d3d      = d3d      , $
                   east     = east     , $
                   shot     = shot     , $
                   equ      = equ      , $
                   debug    = debug    , $
                   save     = save     , $
                   step     = step     , $
                   preview  = preview  , $
                   check    = check    , $
                   all      = all      , $
                   limiter  = limiter  , $
                   boundary = boundary , $
                   mswin  = mswin  , $
                   settings = settings
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  IF (KEYWORD_SET(iter )) THEN machine = 'iter' 
  IF (KEYWORD_SET(cfetr)) THEN machine = 'cfetr' 
  IF (KEYWORD_SET(mast )) THEN machine = 'mast' 
  IF (KEYWORD_SET(cmod )) THEN machine = 'cmod' 
  IF (KEYWORD_SET(west )) THEN machine = 'west' 
  IF (KEYWORD_SET(jet  )) THEN machine = 'jet'  
  IF (KEYWORD_SET(aug  )) THEN machine = 'aug'  
  IF (KEYWORD_SET(d3d  )) THEN machine = 'd3d'  
  IF (KEYWORD_SET(east )) THEN machine = 'east' 

  IF (NOT KEYWORD_SET(step )) THEN step  = 1
  IF (NOT KEYWORD_SET(debug)) THEN debug = 0


  grid_SetupParameters


  IF (KEYWORD_SET(shot)) THEN BEGIN
;
;   --------------------------------------------------------------------
;   USER SECTION  
;   ----------------------------------------------------------------------
;
    CASE machine OF
      ; ----------------------------------------------------------------
      'cfetr': BEGIN

        user_step   = -0.003D  ; the psin step size when scanning for wall contact points ("tangency points")
        user_finish =  3.000D   ; PSIn "base" for setting the radial boundaries of the grid; see /preview
    
        IF (NOT KEYWORD_SET(boundary)) THEN  $
                        ; DN max   ; DN                      ; DN thin    ;  SN        ; default 
          boundary = [ -0.40D,  $  ; -0.40D,  $  ;-0.40D,  $  ; -0.40D,  $ ; -0.40D,  $  ; -0.40D,  $  ; CORE PSIn boundary
                        0.50D,  $  ;  0.50D,  $  ; 0.50D,  $  ;  0.50D,  $ ;  0.20D,  $  ;  1.00D,  $  ; SOL      
                        2.25D,  $  ;  1.75D,  $  ; 1.75D,  $  ;  0.70D,  $ ;  1.00D,  $  ;  1.00D,  $  ; SOL, LFS            
                        0.85D,  $  ;  0.65D,  $  ; 0.65D,  $  ;  0.30D,  $ ;  0.90D,  $  ;  0.90D,  $  ; SOL, HFS            
                       -0.07D,  $  ; -0.07D,  $  ;-0.09D,  $  ; -0.05D,  $ ; -0.05D,  $  ; -0.14D,  $  ; PFZ                 .7
                       -0.045D]    ; -0.01D]     ;-0.025D]    ; -0.01D]    ; -0.25D]     ; -0.25D]     ; PFZ, SECONDARY      

        option.ctr_boundary = boundary * user_finish
    
        option.pol_res_min = 0.010D  ; sets the poloidal length of the cells at the target (m)
        option.pol_res_max = 0.500D  ; sets the poloidal length of the cells upstream (m)
        option.pol_res_exp = 1.0D    ; ignore

        option.rad_res_core_n   = 51      ; ignore 
        option.rad_res_core_min = 0.0005D ; 0.001D  ; radial size of the core ring closest to the separatrix (m)
        option.rad_res_core_max = 0.020D  ; radial size of the core rings far from the separatrix (m)
        option.rad_res_core_exp = 0.500D  ; ignore
        option.rad_res_core_opt = 2       ; ignore

        option.rad_res_sol_n   = 101      ; ignore
        option.rad_res_sol_min = 0.0005D ; 0.001D  ; radial size of the SOL ring closest to the separatrix (m)
        option.rad_res_sol_max = 0.060D ; 0.010D ; 0.200D ; 0.050D   ; radial size of the SOL rings far from the separatrix (m)
        option.rad_res_sol_exp = 1.200D   ; ignore
        option.rad_res_sol_opt = 2        ; ignore

        option.tar_threshold_size = 0.50D ; 0.25D ; maximum length of a target segment

        ; Low res
;        option.pol_res_min = 0.10D0 ; 0.100D  
;        option.pol_res_max = 0.500D  
   
;        option.rad_res_core_min = 0.01D 
;        option.rad_res_core_max = 0.20D 

;        option.rad_res_sol_min = 0.01D0 ; 0.010D   
;        option.rad_res_sol_max = 0.20D0 ; 0.05D0 ; 0.1D0 ; 0.200D  
   
;        option.tar_threshold_size = 10.00D ; 1.00D 
    
        xrange = [ 4.5,10.5]  ; sets the extent of the graphics window
        yrange = [-6.0, 6.5]   

        xpoint_zone = [5.0,8.0,-4.5,6.0]

        change_sign = 1

        IF (KEYWORD_SET(mswin)) THEN BEGIN
          wall_path = mswin + '\shots\' + machine + '\default\'
        ENDIF ELSE BEGIN
          wall_path = '$FUSEHOME/shots/' + machine + '/default/'                 
        ENDELSE
        wall_file = 'vessel_wall.dat'                           ; no pump plenum entrance
        END
      ; ----------------------------------------------------------------
      'cmod': BEGIN

        user_step   = -0.0001D  ; the psin step size when scanning for wall contact points ("tangency points")
        user_finish =  0.025D   ; PSIn "base" for setting the radial boundaries of the grid; see /preview
    
        IF (NOT KEYWORD_SET(boundary)) THEN  $
                                   ; g1101103020.01400.x16.equ     ; users guide
          boundary = [ -0.50D,  $  ;-0.50D,                     $  ; -0.50D,  $  ; CORE PSIn boundary
                       0.20D, $ ; 1.00D,  $  ; 0.20D,                     $  ;  0.20D,  $  ; SOL      
                        1.00D,  $  ; 0.90D,                     $  ;  1.00D,  $  ; SOL, LFS     
                        1.00D,  $  ; 0.25D,                     $  ;  1.00D,  $  ; SOL, HFS
                       -0.10D,  $  ;-0.10D,                     $  ; -0.10D,  $  ; PFZ 
                       -0.10D]     ;-0.050]                        ; -0.10D]     ; PFZ, SECONDARY

        option.ctr_boundary = boundary * user_finish
    
        option.pol_res_min = 0.001D  ; sets the poloidal length of the cells at the target (m)
        option.pol_res_max = 0.050D  ; sets the poloidal length of the cells upstream (m)
        option.pol_res_exp = 1.0D    ; ignore
    
        option.rad_res_core_n   = 51      ; ignore 
        option.rad_res_core_min = 0.001D  ; radial size of the core ring closest to the separatrix (m)
        option.rad_res_core_max = 0.050D  ; radial size of the core rings far from the separatrix (m)
        option.rad_res_core_exp = 1.000D  ; ignore
        option.rad_res_core_opt = 2       ; ignore
    
        option.rad_res_sol_n   = 101      ; ignore
        option.rad_res_sol_min = 0.001D   ; radial size of the core ring closest to the separatrix (m)
        option.rad_res_sol_max = 0.050D   ; radial size of the core rings far from the separatrix (m)
        option.rad_res_sol_exp = 1.000D   ; ignore
        option.rad_res_sol_opt = 2        ; ignore
    
        xrange = [ 0.3,1.2]  ; sets the extent of the graphics window
        yrange = [-0.7,0.7]
;        xrange = [ 0.40, 0.60]   ; inner CC and nose
;        yrange = [-0.40,-0.10]

        xpoint_zone = [0.52,0.75,-0.45,0.45]

        change_sign = 1  

        IF (KEYWORD_SET(mswin)) THEN BEGIN
          wall_path = mswin + '\shots\' + machine + '\default\'
        ENDIF ELSE BEGIN
          wall_path = '$FUSEHOME/shots/' + machine + '/default/'
        ENDELSE
        wall_file = 'vessel_wall4.dat'                           ; no pump plenum entrance
        END
      ; ----------------------------------------------------------------
      'east': BEGIN
        user_step   = -0.0001D  ; the psin step size when scanning for wall contact points ("tangency points")
        user_finish =  0.100D   ; PSIn "base" for setting the radial boundaries of the grid; see /preview
    
        IF (NOT KEYWORD_SET(boundary)) THEN  $
                       ; tests       g041587.006060.x8.equ
          boundary = [  -0.10D,  $  ; -0.40D,  $  ; -0.40D,  $  ; CORE PSIn boundary
                         0.10D,  $  ;  1.00D,  $  ;  1.00D,  $  ; SOL      
                         0.10D,  $  ;  0.50D,  $  ;  0.30D,  $  ; SOL, LFS     
                         0.10D,  $  ;  0.10D,  $  ;  0.13D,  $  ; SOL, HFS
                        -0.02D,  $  ; -0.10D,  $  ; -0.10D,  $  ; PFZ 
                        -0.01D]	    ; -0.05D]	  ; -0.05D]     ; PFZ, SECONDARY

        option.ctr_boundary = boundary * user_finish
    
        option.pol_res_min = 0.002D  ; sets the poloidal length of the cells at the target
        ; tests
        option.pol_res_max = 0.050D  ; sets the poloidal length of the cells upstream
;        option.pol_res_max = 0.100D  ; sets the poloidal length of the cells upstream
        option.pol_res_exp = 1.0D    ; ignore
    
        option.rad_res_core_n   = 51      ; ignore 
        option.rad_res_core_min = 0.002D  ; radial size of the core ring closest to the separatrix 
        option.rad_res_core_max = 0.050D  ; radial size of the core rings far from the separatrix
        option.rad_res_core_exp = 1.000D  ; ignore
        option.rad_res_core_opt = 2       ; ignore
    
        option.rad_res_sol_n   = 101      ; ignore
        option.rad_res_sol_min = 0.001D   ; radial size of the core ring closest to the separatrix 
        option.rad_res_sol_max = 0.050D   ; radial size of the core rings far from the separatrix
        option.rad_res_sol_exp = 1.000D   ; ignore
        option.rad_res_sol_opt = 2        ; ignore

        option.tar_threshold_size = 0.10D

        xrange      = [ 1.3,2.8]  ; sets the extent of the graphics window
        yrange      = [-1.3,1.3]

        xpoint_zone = [1.50,2.10,-1.10,1.10]    

        psi_zone = [0.15,100.0,-100.0,100.0]

                      ; tests
        change_sign = 0        ; 1

        wall_path = '~/fuse/shots/' + machine + '/default/'
        wall_file = 'main_wall.dat'                           

        ; tests
        wall_path = '$FUSEHOME/shots/' + machine + '/' + shot + '/'
        wall_file = 'main_wall_file_update.dat'
        END
      ; ----------------------------------------------------------------
      'd3d': BEGIN
        user_step   = -0.0001D  ; the psin step size when scanning for wall contact points ("tangency points")
        user_finish =  0.1000D   ; PSIn "base" for setting the radial boundaries of the grid; see /preview

        IF (NOT KEYWORD_SET(boundary)) THEN  $
                                   ; g149094.01990      ; g149094.01990         ; reference
                                   ; USN with iw        ; USN with iw 
                                   ; intersection       ; intersection            
          option.ctr_boundary = [  -0.75D*user_finish,  $  ;-1.80D*user_finish,  $  ;-1.80D*user_finish,  $  ; CORE PSIn boundary
                                    0.15D*user_finish,  $  ; 0.25D*user_finish,  $  ; 1.00D*user_finish,  $  ; SOL      
                                    1.80D*user_finish,  $  ; 1.80D*user_finish,  $  ; 1.80D*user_finish,  $  ; SOL, LFS     
                                    0.50D*user_finish,  $  ; 0.50D*user_finish,  $  ; 0.50D*user_finish,  $  ; SOL, HFS
                                   -0.10D*user_finish,  $  ;-0.10D*user_finish,  $  ;-0.18D*user_finish,  $  ; PFZ 
                                   -0.30D*user_finish]     ;-0.30D*user_finish]     ;-0.30D*user_finish]     ; PFZ, SECONDARY
    
        option.pol_res_min = 0.002D  ; sets the poloidal length of the cells at the target
        option.pol_res_max = 0.100D  ; sets the poloidal length of the cells upstream
        option.pol_res_exp = 1.0D    ; ignore
    
        option.rad_res_core_n   = 51      ; ignore 
        option.rad_res_core_min = 0.002D  ; radial size of the core ring closest to the separatrix 
        option.rad_res_core_max = 0.100D  ; radial size of the core rings far from the separatrix
        option.rad_res_core_exp = 2.000D  ; ignore
        option.rad_res_core_opt = 2       ; ignore
    
        option.rad_res_sol_n   = 101      ; ignore
        option.rad_res_sol_min = 0.001D   ; radial size of the core ring closest to the separatrix 
        option.rad_res_sol_max = 0.002D   ; radial size of the core rings far from the separatrix
        option.rad_res_sol_exp = 1.000D   ; ignore
        option.rad_res_sol_opt = 2        ; ignore

        xrange      = [ 0.8,2.5]  ; sets the extent of the graphics window
        yrange      = [-1.5,1.5]

        xpoint_zone = [1.1,1.9,-1.38,1.36]    

        change_sign = 1

        wall_path = '$FUSEHOME/shots/' + machine + '/default/'
        wall_file = 'wall_d3d_shelf_upper_closed.dat'                           ; no pump plenum entrance
        END
      ; ----------------------------------------------------------------
      'aug': BEGIN
        user_step   = -0.00001D  ; the psin step size when scanning for wall contact points ("tangency points")
        user_finish =  0.1000D   ; PSIn "base" for setting the radial boundaries of the grid; see /preview

        IF (NOT KEYWORD_SET(boundary)) THEN  $
                      ; 27100 @ 2.4 s          ; reference                                   
          boundary = [ -0.50D ,  $  ; -0.50D,  $  ; CORE PSIn boundary
                        1.00D ,  $  ;  1.00D,  $  ; SOL      
                        0.80D ,  $  ;  0.20D,  $  ; SOL, LFS     
                        0.50D ,  $  ;  0.15D,  $  ; SOL, HFS
                       -0.11D ,  $  ; -0.07D,  $  ; PFZ 
                       -0.002D]     ; -0.02D]     ; PFZ, SECONDARY
    
        option.ctr_boundary = boundary * user_finish

        ; 27100 @ 2.4 s
        option.ctr_minimum_length[param.PFZ          ] = 0.001D 
        option.ctr_minimum_length[param.PFZ_SECONDARY] = 0.001D 

        option.pol_res_min = 0.002D  ; sets the poloidal length of the cells at the target
        option.pol_res_max = 0.100D  ; sets the poloidal length of the cells upstream
        option.pol_res_exp = 1.0D    ; ignore
    
        option.rad_res_core_n   = 51      ; ignore 
        option.rad_res_core_min = 0.002D  ; radial size of the core ring closest to the separatrix 
        option.rad_res_core_max = 0.100D  ; radial size of the core rings far from the separatrix
        option.rad_res_core_exp = 2.000D  ; ignore
        option.rad_res_core_opt = 2       ; ignore
    
        option.rad_res_sol_n   = 101      ; ignore
        option.rad_res_sol_min = 0.001D   ; radial size of the core ring closest to the separatrix 
        option.rad_res_sol_max = 0.050D   ; radial size of the core rings far from the separatrix
        option.rad_res_sol_exp = 1.000D   ; ignore
        option.rad_res_sol_opt = 2        ; ignore

        xrange      = [ 0.9, 2.4]  ; full vessel
        yrange      = [-1.4, 1.4]
;        xrange      = [ 1.2, 1.7]  ; lower divertor
;        yrange      = [-1.3,-0.8]
;        xrange      = [ 1.2, 1.7]  ; upper divertor
;        yrange      = [ 0.9, 1.3]

        xpoint_zone = [1.30,1.80,-1.20,1.35]

        change_sign = 0

        wall_path = '$FUSEHOME/shots/' + machine + '/default/'
;        wall_file = 'vessel_wall.dat'
;        wall_file = 'vessel_wall_blocky_3.dat'   
;        wall_file = 'vessel_wall_solps_blocky.dat'   
        wall_file = 'vessel_wall_solps.dat'   
        END                                      
     ; ----------------------------------------------------------------
      'mast': BEGIN
        user_step   = -0.0005D  ; the psin step size when scanning for wall contact points ("tangency points")
        user_finish =  0.030D   ; PSIn "base" for setting the radial boundaries of the grid; see /preview
    
        IF (NOT KEYWORD_SET(boundary)) THEN  $

                                 ; 29125                 ; reference   
          option.ctr_boundary = [-0.30D*user_finish,  $  ;-0.30D*user_finish,  $  ; -0.30D*user_finish,  $  ; CORE PSIn boundary
                                  1.00D*user_finish,  $  ; 1.00D*user_finish,  $  ;  1.00D*user_finish,  $  ; SOL      
                                  1.30D*user_finish,  $  ; 0.80D*user_finish,  $  ;  0.80D*user_finish,  $  ; SOL, LFS     
                                  0.12D*user_finish,  $  ; 0.20D*user_finish,  $  ;  0.20D*user_finish,  $  ; SOL, HFS
                                 -0.15D*user_finish,  $  ;-0.30D*user_finish,  $  ; -0.30D*user_finish,  $  ; PFZ 
                                 -0.50D*user_finish]	 ;-0.30D*user_finish]	  ; -0.30D*user_finish]	  ; PFZ, SECONDARY
    
        option.pol_res_min = 0.001D  ; sets the poloidal length of the cells at the target
        option.pol_res_max = 0.050D  ; sets the poloidal length of the cells upstream
        option.pol_res_exp = 1.0D    ; ignore
    
        option.rad_res_core_n   = 51      ; ignore 
        option.rad_res_core_min = 0.002D  ; radial size of the core ring closest to the separatrix 
        option.rad_res_core_max = 0.100D  ; radial size of the core rings far from the separatrix
        option.rad_res_core_exp = 1.000D  ; ignore
        option.rad_res_core_opt = 2       ; ignore
    
        option.rad_res_sol_n   = 101      ; ignore
        option.rad_res_sol_min = 0.001D   ; radial size of the core ring closest to the separatrix 
        option.rad_res_sol_max = 0.020D   ; radial size of the core rings far from the separatrix
        option.rad_res_sol_exp = 1.000D   ; ignore
        option.rad_res_sol_opt = 2        ; ignore

        xrange      = [ 0.05, 2.2]  ; sets the extent of the graphics window
        yrange      = [-2.50, 2.5]

        xpoint_zone = [0.4,1.0,-1.5,1.5]

        psi_zone = [0.15,100.0,-100.0,100.0]

        change_sign = 1

        wall_path = '$FUSEHOME/shots/' + machine + '/default/'
        wall_file = 'main_wall_extended_outer_target.dat'
        END
     ; ----------------------------------------------------------------
      'west': BEGIN
        user_step   = -0.0001D  ; the psin step size when scanning for wall contact points ("tangency points")
        user_finish =  0.100D   ; PSIn "base" for setting the radial boundaries of the grid; see /preview
    
        IF (NOT KEYWORD_SET(boundary)) THEN  $

                                 ; reference   
          option.ctr_boundary = [-0.50D*user_finish,  $  ; CORE PSIn boundary
                                  0.60D*user_finish,  $  ; SOL      
                                  0.30D*user_finish,  $  ; SOL, LFS     
                                  0.50D*user_finish,  $  ; SOL, HFS
                                 -0.05D*user_finish,  $  ; PFZ 
                                 -0.02D*user_finish]	 ; PFZ, SECONDARY
    
        option.pol_res_min = 0.001D  ; sets the poloidal length of the cells at the target
        option.pol_res_max = 0.050D  ; sets the poloidal length of the cells upstream
        option.pol_res_exp = 1.0D    ; ignore
    
        option.rad_res_core_n   = 51      ; ignore 
        option.rad_res_core_min = 0.002D  ; radial size of the core ring closest to the separatrix 
        option.rad_res_core_max = 0.100D  ; radial size of the core rings far from the separatrix
        option.rad_res_core_exp = 1.000D  ; ignore
        option.rad_res_core_opt = 2       ; ignore
    
        option.rad_res_sol_n   = 101      ; ignore
        option.rad_res_sol_min = 0.001D   ; radial size of the core ring closest to the separatrix 
        option.rad_res_sol_max = 0.020D   ; radial size of the core rings far from the separatrix
        option.rad_res_sol_exp = 1.000D   ; ignore
        option.rad_res_sol_opt = 2        ; ignore

        xrange = [ 1.4,3.3]  ; sets the extent of the graphics window
        yrange = [-1.0,1.0]

;        xrange = [2.2 ,2.4 ]  ; sets the extent of the graphics window
;        yrange = [0.64,0.67]

        xpoint_zone = [2.15,3.0,-0.73,0.75]

        psi_zone = [1.7,3.1,-0.8,0.8]

        change_sign = 0

        wall_path = '$FUSEHOME/shots/' + machine + '/psi/'
        wall_file = 'soledge2D.wall_segments_simple_inverted'
        END
     ; ----------------------------------------------------------------
      'jet': BEGIN
        user_step   = -0.0001D ; -0.001D  ; the psin step size when scanning for wall contact points ("tangency points")
        user_finish =  0.300D  ; PSIn "base" for setting the radial boundaries of the grid; see /preview
    
        IF (NOT KEYWORD_SET(boundary)) THEN  $
                     ; 83559     ; limiter      ;  83559    ;
          boundary = [-0.60D , $ ; -3.10D ,  $  ; -0.60D , $; -0.60D ,  $  ; -0.60D ,  $  ; CORE PSIn boundary
                       0.55D , $ ;  0.50D ,  $  ;  0.55D , $;  0.58D ,  $  ;  0.15D ,  $  ; SOL      
                       1.00D , $ ;  1.00D ,  $  ;  1.00D , $;  1.00D ,  $  ;  1.00D ,  $  ; SOL, LFS     
                       1.00D , $ ;  1.00D ,  $  ;  1.00D , $;  1.00D ,  $  ;  1.00D ,  $  ; SOL, HFS
                      -0.08D , $ ; -0.04D ,  $  ; -0.08D , $; -0.04D ,  $  ; -0.04D ,  $  ; PFZ 
                      -0.30D ]	 ; -0.30D ]     ; -0.30D ]  ;      ; -0.30D ]	   ; -0.30D ]  ; PFZ, SECONDARY

        option.ctr_boundary = boundary * user_finish
    
        option.ctr_minimum_length = 0.1D ; 0.3D 

        option.tar_threshold_size = 0.10D
 

        ; default
        option.pol_res_min = 0.001D  ; sets the poloidal length of the cells at the target
        option.pol_res_max = 0.050D  ; sets the poloidal length of the cells upstream

        ; limiter grid
        option.pol_res_min = 0.010D  ; sets the poloidal length of the cells at the target
        option.pol_res_max = 0.100D  ; sets the poloidal length of the cells upstream

        ; 83559
        option.pol_res_min = 0.005D  ; sets the poloidal length of the cells at the target
        option.pol_res_max = 0.100D  ; sets the poloidal length of the cells upstream

        option.pol_res_exp = 1.0D    ; ignore
    
        option.rad_res_core_n   = 51      ; ignore 
        option.rad_res_core_min = 0.002D  ; radial size of the core ring closest to the separatrix 
        option.rad_res_core_max = 0.100D  ; radial size of the core rings far from the separatrix
        option.rad_res_core_exp = 1.000D  ; ignore
        option.rad_res_core_opt = 2       ; ignore
    
        option.rad_res_sol_n   = 101      ; ignore
        option.rad_res_sol_min = 0.001D   ; radial size of the core ring closest to the separatrix 
        option.rad_res_sol_max = 0.020D   ; radial size of the core rings far from the separatrix
        option.rad_res_sol_exp = 1.000D   ; ignore
        option.rad_res_sol_opt = 2        ; ignore

        xrange      = [ 1.6, 4.1]  ; sets the extent of the graphics window
        yrange      = [-2.0, 2.3]
;       xrange      = [ 3.885, 3.895]  ; outer wall limiter
;       yrange      = [-0.0, 0.5]

;        xpoint_zone = [1.5,4.0,-2.0,3.0]
        xpoint_zone = [2.4,3.1,-1.7,2.0]

        change_sign = 1

        wall_path = '$FUSEHOME/shots/' + machine + '/default/'
;        wall_file = 'vessel_wall_vertical_lfs.txt'
;        wall_file = 'vessel_wall_blocky.txt'
        wall_file = 'vessel_wall_blocky2.txt'  ; for the limter grid

;        Used to make grid_jet_83559_52000 (I think...)
        ; 83559:
        wall_path = '$FUSEHOME/shots/' + machine + '/' + shot + '/'
        wall_file = 'wall_jet_blocky_adjusted_plate.dat'

        END
      ; ----------------------------------------------------------------
      ELSE: BEGIN
        PRINT,'ERROR grid_Input: Unknown MACHINE'
        PRINT,'  MACHINE = ',machine
        RETURN, -1
        END
    ENDCASE

    IF (KEYWORD_SET(mswin)) THEN BEGIN
      equ_path  = mswin + '\shots\' + machine + '\' + shot + '\'
    ENDIF ELSE BEGIN
      equ_path  = '$FUSEHOME/shots/' + machine + '/' + shot + '/'    
    ENDELSE
    equ_file  = equ

    machine = machine + '_' + shot  ; not great, fix...

    IF (KEYWORD_SET(settings)) THEN BEGIN

      values = {  $
        version  : 1.0     ,  $
        boundary : boundary   } 

      RETURN, values
    ENDIF

  ENDIF ELSE BEGIN

;   --------------------------------------------------------------------
;   DEVELOPMENT SECTION
;   -------------------------------------------------------------------- 

    CASE machine OF
      'iter': shot = iter
      'mast': shot = mast
      'cmod': shot = cmod
      'west': shot = west
      'jet' : shot = jet 
      'aug' : shot = aug 
      'd3d' : shot = d3d 
      'east': shot = east
      ELSE: BEGIN
        PRINT,'ERROR grid_Input: MACHINE not set'
        RETURN,-1
        END
    ENDCASE

    machine = machine + "_" + STRTRIM(STRING(shot),2)

    ; Limit the region (or not in the initialisation) of the applicability 
    ; of the PSI matrix:
    psi_zone = [0.0,100.0,-100.0,100.0]

;   orange - o-point
;   purple - primary x-point
;   green  - secondary x-point
  
    CASE machine OF
      ; ----------------------------------------------------------------
      'iter_4': BEGIN
        user_step   = -0.003D ; -0.003D
        user_finish =  3.0D ; 3.0D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
$          option.ctr_boundary = [-0.5D*user_finish, user_finish, user_finish, user_finish, -0.18D*user_finish, -0.18D*user_finish]
$          option.ctr_boundary = [-0.5D*user_finish, user_finish, user_finish, 0.59D*user_finish, -0.05*user_finish, -0.25D*user_finish]  ; b?
          option.ctr_boundary = [-0.5D*user_finish, user_finish, user_finish, 0.65D*user_finish, -0.12*user_finish, -0.25D*user_finish]  ; c, with dome
    
        option.ctr_minimum_length            = 0.5D ; b, c
        option.ctr_minimum_length[param.PFZ] = 0.3D ; c
;        option.ctr_minimum_length = 0.75D ; 0.5D ; 0.3D

        option.tar_threshold_size  = 0.10D ; c  maximum length of a target segment
;        option.tar_threshold_angle = 1.0D  ; c  angle between contour and wall below which the target size is increased (important near tangency)

        option.pol_res_min      = 0.010D  ; 0.100D
        option.pol_res_max      = 0.200D	; 0.500D
        option.pol_res_exp      = 1.0D	; 1.0D
          				  
        option.rad_res_core_n   = 51   	; 51   
        option.rad_res_core_min = 0.001D 	; 0.100D 
        option.rad_res_core_max = 0.020D 	; 0.100D 
        option.rad_res_core_exp = 0.500D 	; 2.000D 
        option.rad_res_core_opt = 2      	; 2      

        ; a, b          				  
        ;option.rad_res_sol_n    = 101   	; 101   
        ;option.rad_res_sol_min  = 0.001D	; 0.100D
        ;option.rad_res_sol_max  = 0.020D  ; 0.300D  ; 0.020D  ; 0.100D
        ;option.rad_res_sol_exp  = 1.200D	; 1.000D
        ;option.rad_res_sol_opt  = 2     	; 2     

        ; c
        option.rad_res_sol_n    = 101   	; 101   
        option.rad_res_sol_min  = 0.001D	; 0.100D
        option.rad_res_sol_max  = 0.022D  ; 0.300D  ; 0.020D  ; 0.100D
        option.rad_res_sol_exp  = 1.200D	; 1.000D
        option.rad_res_sol_opt  = 2     	; 2     
    
        xrange      = [ 3.50,8.75]     
        yrange      = [-5.00,5.00]
;        xrange      = [ 4.50, 5.50]  ; dome     
;        yrange      = [-4.50,-3.50]
;        xrange      = [ 4.90, 5.10]  ; dome     
;        yrange      = [-3.80,-3.60]
        xpoint_zone = [4.0,7.0,-3.8,5.5]
        wall_path = '/home/ITER/lisgos/divimp/shots/iter/1514/'
;        wall_file = 'psi_wall_simple2.dat' ; b
        wall_file = 'psi_wall_simple3_dome_blocky2.dat' ; a,b,c
;        wall_path = '/home/ITER/lisgos/fuse_data/test/shots/grid/'
;        wall_file = 'wall_iter.dat
        change_sign =  0
        equ_path = '/home/ITER/lisgos/divimp/shots/iter/equilibria/'
        equ_file = 'feat_001.x8.equ'
        END
      ; ----------------------------------------------------------------
      'iter_10': BEGIN

; a=grid_Run(iter=10,/debug,/save,/all,/check)

        user_step   = -0.003D ; -0.003D; -0.003D
        user_finish =  3.0D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
;         option.ctr_boundary = [-0.3D*user_finish, 0.3D*user_finish, user_finish, user_finish, -0.25D * user_finish, -0.25D*user_finish]
        ; test  option.ctr_boundary = [-0.3D*user_finish, user_finish, user_finish, user_finish, -0.25D * user_finish, -0.25D*user_finish]
        ;   option.ctr_boundary = [-0.4D*user_finish, user_finish, user_finish, 0.75D*user_finish, -0.0580D * user_finish, -0.25D*user_finish]
        ; 10b:  option.ctr_boundary = [-0.4D*user_finish, user_finish, user_finish, 0.44D*user_finish, -0.0580D * user_finish, -0.25D*user_finish]
        ; 10c:  option.ctr_boundary = [-0.4D*user_finish, user_finish, user_finish, 0.90D*user_finish, -0.0580D * user_finish, -0.25D*user_finish]
        ; 10d: 
;        option.ctr_boundary = [-0.4D*user_finish, user_finish, user_finish, 0.90D*user_finish, -0.1400D * user_finish, -0.25D*user_finish]  ; lid
        option.ctr_boundary = [-0.4D*user_finish, user_finish, user_finish, 0.90D*user_finish, -0.2000D * user_finish, -0.25D*user_finish] ; lid, no dome
        ; 10e
        ; option.ctr_boundary = [-3.94D*user_finish, user_finish, user_finish, 0.90D*user_finish, -0.1400D * user_finish, -0.25D*user_finish]
    
        option.ctr_minimum_length = 0.25D ; 0.25D
    
                                         ;  10d?    ; test
        option.pol_res_min      = 0.010D ;  0.010D  ; 0.100D
        option.pol_res_max      = 0.500D ;  0.200D  ; 0.500D
        option.pol_res_exp      = 1.0D	 ;   1.0D   ; 1.0D  
          				        
        option.rad_res_core_n   = 51   	; 51    
        option.rad_res_core_min = 0.001D 	; 0.100D
        ; option.rad_res_core_max = 0.200D 	; 10 e
        option.rad_res_core_max = 0.020D 	; 0.100D
        option.rad_res_core_exp = 0.500D 	; 2.000D
        option.rad_res_core_opt = 2      	; 2     
          				        
        option.rad_res_sol_n    = 101   	; 101   
        option.rad_res_sol_min  = 0.001D	; 0.100D
        option.rad_res_sol_max  = 0.050D    ; 0.10D ; f
;        option.rad_res_sol_max  = 0.020D        ; 0.300D   ; 0.020D  ; 0.100D
        option.rad_res_sol_exp  = 1.200D	; 1.000D
        option.rad_res_sol_opt  = 2     	; 2     
    
        option.tar_threshold_size = 0.25D ; f
;        option.tar_threshold_size = 0.10D
    
        xrange      = [ 3.5,8.75]     
        yrange      = [-5.0,5.0]
        ;xrange      = [ 5.0, 6.5]     ; lower, outer target
        ;yrange      = [-4.5,-3.0]
        ;xrange      = [ 4.5, 5.5]     ; lower x-point
        ;yrange      = [-4.0,-2.5]
    
        xpoint_zone = [4.0,7.0,-3.8,5.5]
        wall_path = '/home/ITER/lisgos/divimp/shots/iter/1514/'
;        wall_file = 'psi_wall_simple.dat'
;        wall_file = 'psi_wall_simple2.dat'
;        wall_file = 'psi_wall_simple3.dat' ; 10c 
;        wall_file = 'psi_wall_simple3_dome_blocky.dat'
        wall_file = 'psi_wall_simple3_dome_blocky2_lid_no_dome.dat'
        change_sign =  0
        equ_path = '/home/ITER/lisgos/divimp/shots/iter/1514/'
;        equ_file = 'Baseline2008-li0.70.x4.equ'
        equ_file = 'Baseline2008-li0.70.x8.equ'
        END
      ; ----------------------------------------------------------------    
      'mast_1': BEGIN
        user_step   = -0.0005D
        user_finish =  0.03D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.3D*user_finish, user_finish, 0.8D*user_finish, 0.2D*user_finish, -0.3D*user_finish, -0.3D*user_finish]
    
        option.pol_res_min = 0.050D  ;  ok        bad
        option.pol_res_max = 0.200D  ;  0.050D  ; 0.100D 
        option.pol_res_exp = 1.0D
    
        option.rad_res_core_n   = 51   
        option.rad_res_core_min = 0.050D 
        option.rad_res_core_max = 0.100D 
        option.rad_res_core_exp = 1.000D 
        option.rad_res_core_opt = 2      
    
        option.rad_res_sol_n   = 101   
        option.rad_res_sol_min = 0.050D
        option.rad_res_sol_max = 0.100D   ; 0.100D
        option.rad_res_sol_exp = 1.000D
        option.rad_res_sol_opt = 2     
    
        xrange      = [ 0.05, 2.2]
        yrange      = [-2.5 , 2.5]
;        xrange      = [ 1.35, 1.4] ; outer midplane - zoom
;        yrange      = [-0.7 , 0.7]
;        xrange      = [ 0.50, 1.4] ; outer midplane
;        yrange      = [-1.5 , 1.5]
;        xrange      = [ 0.2 , 1.1] ; lower x-point
;        yrange      = [-1.5 ,-0.5]
;        xrange      = [ 0.2 , 1.1] ; upper x-point
;        yrange      = [ 0.5 , 1.5]
;        xrange      = [ 0.4 , 1.6] ; upper, outer target
;        yrange      = [ 0.9 , 1.9]
;        xrange      = [ 0.2, 0.30] ; inner midplane
;        yrange      = [-0.2 , 0.2]
        xpoint_zone = [0.4,1.0,-1.5,1.5]
        wall_path = '/home/ITER/lisgos/divimp/shots/mast/default/'
        wall_file = 'main_wall.dat'
        change_sign = 0
        ; trouble with multiple target intersections when far out
        equ_path = '~/fuse_data/mast/shots/24867/'
        equ_file = '24867_335.equ'
        ; Fails on the inner wall because of limited efit data
        ;equ_path = '~/divimp/shots/mast/24860/'
        ;equ_file = 'carre.24860_240.equ'
        ; works, but get into trouble if going to far out(>0.1) as multiple wall interesections appear
        ; can perhaps have a selection rule based on the appearance of an odd number of intersections, i.e.
        ; 3 rather than 4...
        equ_path = '/home/ITER/lisgos/divimp/shots/mast/13018_250/'
        equ_file = 'sonnet_13018_250.equ'
        END
      ; ----------------------------------------------------------------    
      'mast_2': BEGIN
        user_step   = -0.0001D
        user_finish =  0.01D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.7D*user_finish, user_finish, 1.0D*user_finish, 0.15D*user_finish, -0.6D*user_finish, -0.6D*user_finish]
    
        xrange      = [ 0.05, 2.2]
        yrange      = [-2.5 , 2.5]
    
        xpoint_zone = [0.4,1.0,-1.5,1.5]
        psi_zone    = [0.15,100.0,-100.0,100.0]
        wall_path = '/home/ITER/lisgos/divimp/shots/mast/default/'
        wall_file = 'main_wall.dat'
        change_sign = 1
        equ_path = '/home/ITER/lisgos/divimp/shots/mast/26798/'
        equ_file = 'mast_26798_220.x32.equ'
        END
      ; ----------------------------------------------------------------    
      'mast_3': BEGIN
        user_step   = -0.0001D
        user_finish =  0.01D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.7D*user_finish, user_finish, 1.0D*user_finish, 0.15D*user_finish, -0.6D*user_finish, -0.6D*user_finish]
    
        xrange      = [ 0.05, 2.2]
        yrange      = [-2.5 , 2.5]
    
        xpoint_zone = [0.4,1.0,-1.5,1.5]
        psi_zone    = [0.15,100.0,-100.0,100.0]
        wall_path = '/home/ITER/lisgos/divimp/shots/mast/default/'
        wall_file = 'main_wall.dat'
        change_sign = 1
        equ_path = '/home/ITER/lisgos/divimp/shots/mast/27741/'
        equ_file = 'mast_27741.220.equ'
        END
      ; ----------------------------------------------------------------
      'cmod_1': BEGIN
        user_step   = -0.0005D0
        user_finish =  0.025D
    
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [ -3.0D*user_finish, user_finish, user_finish, user_finish, -0.1D*user_finish, -0.1D*user_finish]
    
        option.pol_res_min = 0.001D
        option.pol_res_max = 0.010D
        option.pol_res_exp = 1.0D
    
        option.rad_res_core_n   = 51   
        option.rad_res_core_min = 0.050D 
        option.rad_res_core_max = 0.100D 
        option.rad_res_core_exp = 1.000D 
        option.rad_res_core_opt = 2      
    
        option.rad_res_sol_n   = 101   
        option.rad_res_sol_min = 0.010D
        option.rad_res_sol_max = 0.200D   ; 0.100D
        option.rad_res_sol_exp = 1.000D
        option.rad_res_sol_opt = 2     
    
        xrange      = [0.3,1.2]
        yrange      = [-0.7,0.7]
;        xrange      = [0.4,0.5]   ; inner midplane
;        yrange      = [-0.1,0.1]
;        xrange      = [0.5,0.8]   ; upper x-point
;        yrange      = [0.2,0.7]
;        xrange      = [0.69 ,0.71 ]   ; upper x-point
;        yrange      = [0.405,0.415]
        xpoint_zone = [0.52,0.75,-0.45,0.45]
        wall_path = '/home/ITER/lisgos/divimp/shots/cmod/1100303017_01380/'
        wall_file = 'vessel_wall4.dat'  ; no pump plenum entrance
;        wall_file = 'vessel_wall.dat'  
        change_sign = 1
        ; works... but need to re-address the problem of "funny business" in the param.PFZ, see vessel_wall4.dat file
        equ_path  = '/home/ITER/lisgos/divimp/shots/cmod/1100303017_01380/'
        equ_file  = 'cmod.1100303017.01380.x16.equ'
        END
      ; ----------------------------------------------------------------
      'cmod_2': BEGIN
        user_step   = -0.0005D0
        user_finish =  0.025D
    
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [ -3.0D*user_finish, user_finish, user_finish, user_finish, -0.1D*user_finish, -0.1D*user_finish]
    
        option.pol_res_min = 0.001D
        option.pol_res_max = 0.010D
        option.pol_res_exp = 1.0D
    
        option.rad_res_core_n   = 51   
        option.rad_res_core_min = 0.050D 
        option.rad_res_core_max = 0.100D 
        option.rad_res_core_exp = 1.000D 
        option.rad_res_core_opt = 2      
    
        option.rad_res_sol_n   = 101   
        option.rad_res_sol_min = 0.010D
        option.rad_res_sol_max = 0.200D   ; 0.100D
        option.rad_res_sol_exp = 1.000D
        option.rad_res_sol_opt = 2     
    
        xrange      = [0.3,1.2]
        yrange      = [-0.7,0.7]
;        xrange      = [0.4,0.5]   ; inner midplane
;        yrange      = [-0.1,0.1]
;        xrange      = [0.5,0.8]   ; upper x-point
;        yrange      = [0.2,0.7]
;        xrange      = [0.69 ,0.71 ]   ; upper x-point
;        yrange      = [0.405,0.415]
        xpoint_zone = [0.52,0.75,-0.45,0.45]
        wall_path = '/home/ITER/lisgos/divimp/shots/cmod/1100303017_01380/'
        wall_file = 'vessel_wall4.dat'  ; no pump plenum entrance
;        wall_file = 'vessel_wall5.dat'  ; inner wall further in to test limiter grids
;        wall_file = 'vessel_wall.dat'   ; full wall, with plenum entrance
        change_sign = 1
        ; works... but need to re-address the problem of "funny business" in the param.PFZ, see vessel_wall4.dat file
        equ_path  = '/home/ITER/lisgos/divimp/shots/cmod/1100303017_01380/'
        equ_file  = 'cmod.1100303017.01380.x4.equ'
        END
      ; ----------------------------------------------------------------
      'west_1': BEGIN
        user_step   = -0.01D
        user_finish =  0.80D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.1D*user_finish, user_finish, user_finish, user_finish, -0.3D*user_finish, -0.3D*user_finish]
        xrange      = [1.4,3.3]
        yrange      = [-1.0,1.0]
        xpoint_zone = [2.15,2.7,-0.73,0.75]
        wall_path = '/home/ITER/lisgos/divimp/shots/ts2/upgrade/'
        wall_file = 'vessel_wall.dat'
        change_sign = 0
        ; dies on corner - might be same as D3D problem
         equ_path = '/home/ITER/lisgos/divimp/shots/ts2/600kA_LN_1cm/'
         equ_file = 'TS2_600kA_LN.x2.cr.equ'
        ; gets pissed off by the second x-point that's just inside the vessel
        ; equ_path = '/home/ITER/lisgos/divimp/shots/ts2/1MA_LN_1cm/'
        ; equ_file = 'TS2_1MA_LN.x2.equ'
        ; pissed off by CDNness... to be expected 
        ;equ_file = '/home/ITER/lisgos/divimp/shots/ts2/600kA_CDN/TS2_600kA.x2.equ'
        ;equ_file = '/home/ITER/lisgos/divimp/shots/ts2/600kA_CDN/TS2_600kA.x2.equ'
        END
      ; ----------------------------------------------------------------
      'jet_1': BEGIN
        user_step   = -0.001D ; -0.0001D ; -0.001D
        user_finish =  0.3D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.6D*user_finish, user_finish, user_finish, user_finish, -0.04D*user_finish, -0.3D*user_finish]
    
        option.ctr_minimum_length = 0.1D ; 0.3D 
    
        option.tar_threshold_size = 0.10D
        option.pol_res_min = 0.001D
        option.pol_res_max = 0.05D ; 0.250D
        option.pol_res_exp = 1.0D
    
;        option.rad_res_sol_n   = 101   
;        option.rad_res_sol_min = 0.100D
;        option.rad_res_sol_max = 0.100D
;        option.rad_res_sol_exp = 1.000D
;        option.rad_res_sol_opt = 2     
    
        xrange      = [ 1.6, 4.1]  ; full
        yrange      = [-2.0, 2.3]
;        xrange      = [ 2.3, 3.0]  ; divertor
;        yrange      = [-1.8,-1.2]
;        xrange      = [ 2.45, 2.65]  ; divertor
;        yrange      = [-1.8,-1.5]
    
;        xrange = [ 1.8037636     ,  1.8170778]
;        yrange = [ -0.10697396    ,  0.55873864]
    
;       xrange = [ 1.8037636     ,  1.82]
;       yrange = [ 0.55, 0.58]
    
;       xrange = [ 1.8037636     ,  2.0000000]
;       yrange = [ -1.10697396    ,  0.55873864]
    
;        xres = 0.05
;        yres = 0.2
;        xrange = [ 3.84 -xres, 3.84 +xres]
;        yrange = [-0.180-yres,-0.180+yres]
;        yrange = [-0.255-yres,-0.255+yres]
    
;        xrange      = [ 2.3, 2.9]  ; x-point
;        yrange      = [-1.8,-1.3]
    
        xpoint_zone = [2.4,3.1,-1.8,2.0]
        wall_path = '/home/ITER/lisgos/divimp/shots/jet/default/'
;        wall_file = 'vessel_wall.dat'
        wall_file = 'vessel_wall_blocky.txt'
        change_sign = 1
        equ_path = '/home/ITER/lisgos/divimp/shots/jet/68124_49000/'
        equ_file = 'JET_68124_49000.x4.equ' ; LSND
        END
      ; ----------------------------------------------------------------
      'jet_2': BEGIN
        user_step   = -0.001D
        user_finish =  0.2D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.5D*user_finish, user_finish, user_finish, user_finish, -0.08D*user_finish, -0.1D*user_finish]
    
        option.ctr_minimum_length = 0.1D
        option.tar_threshold_size = 0.15D  ; 0.50D   ;  *** trying to get this smaller ***
        option.pol_res_min = 0.001D ; 0.001D
        option.pol_res_max = 0.10D ; 0.250D
        option.pol_res_exp = 1.0D
    
;        option.rad_res_sol_n   = 101   
        option.rad_res_sol_min = 0.0005D                                          ; 0.001 -- these are safe
        option.rad_res_sol_max = 0.010D    ; *** still playihg around here **     ; 0.05    GO BACK TO THESE, lack of /save in call from IDL screwed me up...
        option.rad_res_sol_exp = 1.000D
;        option.rad_res_sol_opt = 2     
    
;        option.rad_res_core_n   = 51           ; radial resolution at the separatrix, core
        option.rad_res_core_min = 0.0005D       ; 
;        option.rad_res_core_max = 0.500D        ; 
;        option.rad_res_core_exp = 1.000D        ; 
;        option.rad_res_core_opt = param.TANH   ; 
    
        xrange      = [ 1.6, 4.1]  ; full
        yrange      = [-2.0, 2.3]
;        xrange      = [ 2.3, 3.0]  ; divertor
;        yrange      = [-1.8,-1.2]
;        xrange      = [ 1.6, 4.1]  ; lower half
;        yrange      = [-1.8,-0.5]
        xrange      = [ 2.0, 2.7]  ; crown
        yrange      = [ 1.6, 2.1]
;        xrange      = [ 2.0, 2.2]  ; local zoom
;        yrange      = [ 1.6, 1.8]
        xpoint_zone = [2.4,3.1,-1.8,2.0]
        wall_path = '/home/ITER/lisgos/divimp/shots/jet/default/'
        wall_file = 'vessel_wall_blocky_adjusted_plate2.txt'
        change_sign = 1
        equ_path = '/home/ITER/lisgos/divimp/shots/jet/80295/'
        equ_file = 'jet_80295_18600.16x.equ'
        END
      ; ----------------------------------------------------------------
      'aug_1': BEGIN
        user_step   = -0.0005D
        user_finish =  0.1D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.1D*user_finish, user_finish, user_finish, user_finish, -0.3D*user_finish, -0.3D*user_finish]
        xrange      = [ 0.9, 2.4]  ; full
        yrange      = [-1.4, 1.4]
;        xrange      = [ 1.2, 1.7]  ; lower divertor
;        yrange      = [-1.3,-0.8]
        xpoint_zone = [ 1.30,1.75,-1.20,1.35]
        wall_path = '/home/ITER/lisgos/divimp/shots/aug/default/'
        wall_file = 'vessel_wall.dat'
        change_sign = 0
        equ_path  = '/home/ITER/lisgos/divimp/shots/aug/22575_3000/' 
        equ_file  = '22575_72x18.sonnet.equ'
        END
    
      ; ----------------------------------------------------------------
      'aug_2': BEGIN
        user_step   = -0.0005D
        user_finish =  0.05D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
                               ; core               SOL          LFS              HFS          PFZ                secondary PFZ
          option.ctr_boundary = [-0.5D*user_finish, user_finish, 0.9*user_finish, user_finish, -0.1D*user_finish, -0.05D*user_finish]
        xrange      = [ 0.9, 2.4]  ; full
        yrange      = [-1.4, 1.4]
;        xrange      = [ 1.2, 1.7]  ; lower divertor
;        yrange      = [-1.3,-0.8]
        xpoint_zone = [ 1.30,1.75,-1.20,1.35]
        wall_path = '/home/ITER/lisgos/divimp/shots/aug/default/'
        wall_file = 'vessel_wall.dat'
        change_sign = 0
        equ_path  = '/home/ITER/lisgos/fuse/src/equtrn/'
        equ_file  = '27332.3900.AUGD.EQH.00.x2.equ'
        END
    
      ; ----------------------------------------------------------------
      'd3d_1': BEGIN
        user_step   = -0.0001D ; -0.0001D ; -0.0005D
        user_finish =  0.08D ; 0.10D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
;          option.ctr_boundary = [-0.9D*user_finish, user_finish, user_finish, 0.7D*user_finish, -0.15D*user_finish, -0.08D*user_finish]  ; grid_d3d_1a 
          option.ctr_boundary = [-3.0D*user_finish, user_finish, 1.7D*user_finish, 0.9D*user_finish, -0.15D*user_finish, -0.08D*user_finish] 
    
        option.pol_2nd_niter    = 1
        option.ctr_2nd_xpt_slow = 1.0D
    
        option.rad_res_core_n   = 51       ; grid_d3d_1d
        option.rad_res_core_min = 0.0005D 
        option.rad_res_core_max = 0.0200D 
        option.rad_res_core_exp = 5.000D 
        option.rad_res_core_opt = 2     
    
        xrange = [ 0.8, 2.5]  ; full
        yrange = [-1.5, 1.5]
;        xrange = [ 1.15, 1.3]  ; upper PFR
;        yrange = [ 1.1 , 1.4]
;        xrange = [ 1.1, 2.2]  ; upper half
;        yrange = [ 0.2, 1.4]
;        xrange      = [ 1.5, 1.9]  ; lower divertor
;        yrange      = [-1.4,-1.0]
;        xrange      = [ 1.5, 1.6]  ; lower outer target
;        yrange      = [-1.4,-1.3]
    
        xpoint_zone = [ 1.2,1.9,-1.38,1.36]
        wall_path = '/home/ITER/lisgos/divimp/shots/d3d/default/'
        wall_file = 'vessel_wall2.dat'
        change_sign = 1
        ; ERROR grid_GetOrthogonal: Normal vector requested away from vertex on outer wall
        equ_path  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
;        equ_file  = 'g105500.03500.x2.equ' ; LSND, 2nd x-point just inside the vessel
        equ_file  = 'g105500.03500.x8.equ' ; LSND, 2nd x-point just inside the vessel
        END
      ; ----------------------------------------------------------------
      'd3d_2': BEGIN
        user_step   = -0.0005D
        user_finish =  0.10D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.1D*user_finish, user_finish, user_finish, user_finish, -0.3D*user_finish, -0.3D*user_finish]
        xrange = [ 0.8, 2.5]  ; full
        yrange = [-1.5, 1.5]
;       xrange = [ 1.0, 1.6]  ; upper PFR
;       yrange = [ 1.1, 1.4]
        xpoint_zone = [ 1.2,1.9,-1.38,1.36]
        wall_path = '/home/ITER/lisgos/divimp/shots/d3d/default/'
        wall_file = 'vessel_wall2.dat'
        change_sign = 1
        ; works...
        equ_path  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
        equ_file  = 'g131245.02800.x2.equ' ; USND
        END
      ; ----------------------------------------------------------------
      'd3d_3': BEGIN
        user_step   = -0.0005D
        user_finish =  0.10D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.1D*user_finish, user_finish, user_finish, user_finish, -0.3D*user_finish, -0.3D*user_finish]
        xrange = [ 0.8, 2.5]  ; full
        yrange = [-1.5, 1.5]
;        xrange = [ 1.0, 1.4]  ; upper PFR
;        yrange = [ 1.1, 1.4]
        xpoint_zone = [ 1.2,1.9,-1.38,1.36]
        wall_path = '/home/ITER/lisgos/divimp/shots/d3d/default/'
        wall_file = 'vessel_wall2.dat'
        change_sign = 1
        ; works...
        equ_path  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
        equ_file  = 'g134585.03005.x4.equ' ; UDND
        END
      ; ----------------------------------------------------------------
      'd3d_4': BEGIN
        user_step   = -0.0005D
        user_finish =  0.10D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.1D*user_finish, user_finish, user_finish, user_finish, -0.3D*user_finish, -0.3D*user_finish]
        xrange = [ 0.8, 2.5]  ; full
        yrange = [-1.5, 1.5]
;        xrange = [ 1.0, 1.4]  ; upper PFR
;        yrange = [ 1.1, 1.4]
        xpoint_zone = [ 1.2,1.9,-1.38,1.36]
        wall_path = '/home/ITER/lisgos/divimp/shots/d3d/default/'
        wall_file = 'vessel_wall2.dat'
        change_sign = 1
    
        ; Corner at the entrance to the upper plenum causes problems, although the 
        ; code seems to cope, sort of, by acting defensively.  This is related to an 
        ; OUTSTANDING ISSUE where there are multiple valid segments along a contour,
        ; which can happen if the angle of the contours relative to the surface is such
        ; that a tangency point is not associated with the entrance to the secondary
        ; region of the vaccum vessel -- the pump plenums in this case and for C-Mod
        ; Perhaps filter the contour when an unusual number of intersections are found, 
        ; that is, more than two, and try to decide which to keep...
    
        ; I THINK THIS PROBLEM IS PARTLY param.SOLVED ACTUALLY...
        ; From the behaviour of MAST and JET runs... and even this D3D run.  What I'm not 
        ; clear on is why the D3D case stops expanding radially when this issue is encountered.
    
        equ_path  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
        equ_file  = 'g135487.03865.x4.equ' ; LSND
        END
      ; ----------------------------------------------------------------
      'd3d_5': BEGIN
        user_step   = -0.0005D
        user_finish =  0.10D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.1D*user_finish, user_finish, user_finish, user_finish, -0.3D*user_finish, -0.3D*user_finish]
        xrange = [ 0.8, 2.5]  ; full
        yrange = [-1.5, 1.5]
;        xrange = [ 1.0, 1.4]  ; upper PFR
;        yrange = [ 1.1, 1.4]
        xpoint_zone = [ 1.2,1.9,-1.38,1.36]
        wall_path = '/home/ITER/lisgos/divimp/shots/d3d/default/'
        wall_file = 'vessel_wall2.dat'
        change_sign = 1
        ; fails... but for the same reason as the above equilibrium, although this time
        ; for the upper inner plenum entrance...
        equ_path  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
        equ_file  = 'g123417.03500.x2.equ' ; LSND, upper point just outside the vessel
        END
      ; ----------------------------------------------------------------
      'd3d_6': BEGIN
        user_step   = -0.0001D  ; -0.0005D
        user_finish =  0.10D 
        IF (NOT KEYWORD_SET(boundary)) THEN  $
;          option.ctr_boundary = [-0.9D*user_finish, user_finish, user_finish, user_finish, -0.15D*user_finish, -0.08D*user_finish]
;          option.ctr_boundary = [-0.9D*user_finish, user_finish, 1.5D*user_finish, user_finish, -0.175D*user_finish, -0.3D*user_finish] ; original grid for Dave, old wall (no shelf)
;          option.ctr_boundary = [-0.9D*user_finish, user_finish, 1.5D*user_finish, 0.67D*user_finish, -0.16D*user_finish, -0.3D*user_finish] ; grid_d3d_6a
          option.ctr_boundary = [-3.5D*user_finish, user_finish, 1.9D*user_finish, 0.85D*user_finish, -0.16D*user_finish, -0.3D*user_finish] ; grid_d3d_6b,c,d
    
        option.pol_2nd_niter = 1
    
        option.rad_res_core_n   = 51       ; grid_d3d_6d
        option.rad_res_core_min = 0.0005D 
        option.rad_res_core_max = 0.0200D 
        option.rad_res_core_exp = 5.000D 
        option.rad_res_core_opt = 2     
    
    
        xrange = [ 0.8, 2.5]  ; full
        yrange = [-1.5, 1.5]
;        xrange = [ 1.0, 1.4]  ; upper PFR
;        yrange = [ 1.1, 1.4]
        xpoint_zone = [ 1.1,1.9,-1.38,1.36]
        wall_path = '/home/ITER/lisgos/divimp/shots/d3d/default/'
;        wall_file = 'wall_d3d_shelf.dat'
        wall_file = 'vessel_wall4.dat'  ; old wall with new shelf
;        wall_file = 'vessel_wall3.dat'  ; old wall
        change_sign = 1
        equ_path  = '/home/ITER/lisgos/divimp/shots/d3d/146198_3010/'
        equ_file  = 'g146198.03010.x8.equ'  ; screws up on the inside
        END
      ; ----------------------------------------------------------------
      'd3d_7': BEGIN
        user_step   = -0.0005D
        user_finish =  0.10D 
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.9D*user_finish, 0.775D*user_finish, 1.5D*user_finish, 0.67D*user_finish, -0.03D*user_finish, -0.3D*user_finish]
    
        option.pol_2nd_niter = 1
    
        xrange = [ 0.8, 2.5]  ; full
        yrange = [-1.5, 1.5]
    
;       xrange = [1.55,1.65]  
;       yrange = [-1.255,-1.245]
    
    
        xpoint_zone = [ 1.1,1.9,-1.38,1.36]
        wall_path = '/home/ITER/lisgos/divimp/shots/d3d/default/'
        wall_file = 'wall_d3d_shelf_slant.dat'
        change_sign = 1
        equ_path  = '/home/ITER/lisgos/divimp/shots/d3d/140422_3765/'
        equ_file  = 'g140422.03765.x8.equ
        END
      ; ----------------------------------------------------------------
      'd3d_8': BEGIN
        user_step   = -0.0005D
        user_finish =  0.10D 
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.9D*user_finish, user_finish, 1.5D*user_finish, 0.67D*user_finish, -0.175D*user_finish, -0.3D*user_finish]
    
        option.pol_2nd_niter = 1
    
        xrange = [ 0.8, 2.5]  ; full
        yrange = [-1.5, 1.5]
    
        xpoint_zone = [ 1.1,1.9,-1.38,1.36]
        wall_path = '/home/ITER/lisgos/divimp/shots/d3d/default/'
        wall_file = 'wall_d3d_shelf.dat'
        change_sign = 1
        equ_path  = '/home/ITER/lisgos/divimp/shots/d3d/140424_3005/'
        equ_file  = 'g140424.03005.x8.equ'
        END
      ; ----------------------------------------------------------------
      'd3d_9': BEGIN
        user_step   = -0.0001D
        user_finish =  0.10D 
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-1.8D*user_finish, user_finish, 1.8D*user_finish, 0.500D*user_finish, -0.175D*user_finish, -0.3D*user_finish]
    
        option.pol_2nd_niter = 1
    
;        option.ctr_2nd_xpt_slow = 0.005D
        option.ctr_2nd_xpt_frac = 0.01D
    
        option.rad_res_core_n   = 51       ; grid_d3d_9d
        option.rad_res_core_min = 0.0005D 
        option.rad_res_core_max = 0.0200D 
        option.rad_res_core_exp = 5.000D 
        option.rad_res_core_opt = 2     
    
        xrange = [ 0.8, 2.5]  ; full
        yrange = [-1.5, 1.5]
;        xrange = [ 1.3, 1.4]  ; 
;        yrange = [ 1.2, 1.4]
    
        xpoint_zone = [ 1.1,1.9,-1.38,1.36]
        wall_path = '/home/ITER/lisgos/divimp/shots/d3d/default/'
;        wall_file = 'wall_d3d_upper_duct_closed.dat'
;        wall_file = 'vessel_wall4.dat'  ; old wall with new shelf
;        wall_file = 'wall_d3d_shelf_slant2.dat'
        wall_file = 'wall_d3d_shelf_upper_closed.dat'
        change_sign = 1
        equ_path  = '/home/ITER/lisgos/divimp/shots/d3d/146210_2005/'
        equ_file  = 'g146210.02005.x8.equ'
        END
      ; ----------------------------------------------------------------
      'd3d_10': BEGIN
        user_step   = -0.0001D
        user_finish =  0.10D 
        IF (NOT KEYWORD_SET(boundary)) THEN  $
;          option.ctr_boundary = [-0.5D*user_finish, 0.68D*user_finish, 1.4D*user_finish, 0.700D*user_finish, -0.175D*user_finish, -0.3D*user_finish]
          option.ctr_boundary = [-3.3D*user_finish, 0.68D*user_finish, 1.4D*user_finish, 0.700D*user_finish, -0.175D*user_finish, -0.3D*user_finish]
    
        option.pol_2nd_niter = 1
    
;        option.ctr_2nd_xpt_slow = 0.005D
;        option.ctr_2nd_xpt_frac = 0.01D
    
        option.ctr_minimum_length = 0.2D
    
        option.rad_res_core_n   = 51       ; grid_d3d_10d
        option.rad_res_core_min = 0.0005D 
        option.rad_res_core_max = 0.0200D 
        option.rad_res_core_exp = 5.000D 
        option.rad_res_core_opt = 2     
    
        xrange = [ 0.8, 2.5]  ; full  1.60 -0.35  1.70 -0.25
        yrange = [-1.5, 1.5]
;        xrange = [ 1.4, 1.6]  
;        yrange = [-1.2,-0.6]
;        xrange = [ 1.55, 1.65] 
;        yrange = [-0.55,-0.45]
        xrange = [ 1.60, 1.75]
        yrange = [-0.35,-0.15]
    
        xpoint_zone = [ 1.0,1.9,-1.38,1.36]
        wall_path = '/home/ITER/lisgos/divimp/shots/d3d/default/'
;        wall_file = 'wall_d3d_upper_duct_closed.dat'
;        wall_file = 'vessel_wall4.dat'  ; old wall with new shelf
;        wall_file = 'wall_d3d_shelf_slant2.dat'
        wall_file = 'wall_d3d_shelf_upper_closed.dat'
        change_sign = 1
        equ_path  = '/home/ITER/lisgos/divimp/shots/d3d/144981_3990/'
        equ_file  = 'g144981.03990.x8.equ'
        END
      ; ----------------------------------------------------------------
      'd3d_11': BEGIN
        user_step   = -0.0001D
    
;        user_finish =  0.02D   ; for _11a and _11b
;        IF (NOT KEYWORD_SET(boundary)) THEN  $
;          option.ctr_boundary = [-1.8D*user_finish, user_finish, 1.0D*user_finish, 1.000D*user_finish, -0.175D*user_finish, -0.3D*user_finish]
;   
        user_finish =  0.10D    ; for _11c
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.8D*user_finish, user_finish, 1.0D*user_finish, 1.000D*user_finish, -0.045D*user_finish, -0.3D*user_finish]
;          option.ctr_boundary = [-0.8D*user_finish, user_finish, 1.0D*user_finish, 1.000D*user_finish, -0.08D*user_finish, -0.3D*user_finish]  ; epic fail
    
    
        option.pol_2nd_niter = 1
    
;        option.ctr_2nd_xpt_slow = 0.005D
        option.ctr_2nd_xpt_frac = 0.01D
    
;        option.rad_res_core_n   = 51       ; grid_d3d_11a and _11c (good separatrix resolution)
;        option.rad_res_core_min = 0.0005D 
;        option.rad_res_core_max = 0.0200D 
;        option.rad_res_core_exp = 5.000D 
;        option.rad_res_core_opt = 2     
;;;;      option.rad_res_sol_n   = 101   
;        option.rad_res_sol_min = 0.0005D   
;        option.rad_res_sol_max = 0.010D    
;        option.rad_res_sol_exp = 1.000D
;;;;      option.rad_res_sol_opt = 2     
;   
        option.rad_res_core_n   = 51       ; grid_d3d_11b and _11d
        option.rad_res_core_min = 0.010D 
        option.rad_res_core_max = 0.020D 
        option.rad_res_core_exp = 5.000D 
        option.rad_res_core_opt = 2     
;        option.rad_res_sol_n   = 101   
        option.rad_res_sol_min = 0.0005D   
        option.rad_res_sol_max = 0.100D    
        option.rad_res_sol_exp = 1.000D ; 1.000D
;        option.rad_res_sol_opt = 2     
;        option.tar_threshold_size = 0.50D 
    
        option.pol_res_min = 0.001D
        option.pol_res_max = 0.150D ; 0.100D
        option.pol_res_opt = 1     
        option.pol_res_exp = 0.750D
    
        xrange = [ 0.8, 2.5]  ; full
        yrange = [-1.5, 1.5]
    
        xpoint_zone = [ 1.1,1.9,-1.38,1.36]
        wall_path = '/home/ITER/lisgos/divimp/shots/d3d/default/'
;        wall_file = 'wall_d3d_upper_duct_closed.dat'
;        wall_file = 'vessel_wall4.dat'  ; old wall with new shelf
;        wall_file = 'wall_d3d_shelf_slant2.dat'
;        wall_file = 'wall_d3d_shelf_upper_closed.dat'
;        wall_file = 'wall_d3d_shelf_blocky_upper_dome.dat'
        wall_file = 'wall_d3d_shelf_modified.dat'
        change_sign = 1
        equ_path  = '/home/ITER/lisgos/divimp/shots/d3d/134585_3005/'
        equ_file  = 'g134585.03005.x8.equ'
        END
      ; ----------------------------------------------------------------
      'east_1': BEGIN
        user_step   = -0.001D
        user_finish =  0.1000D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.4D*user_finish, user_finish, 0.5D*user_finish, 0.4D*user_finish, -0.15D*user_finish, -0.3D*user_finish]
    
;        option.pol_res_min = 0.010D
;        option.pol_res_max = 0.300D
;        option.pol_res_exp = 1.0D
    
;        option.rad_res_core_n   = 51   
;        option.rad_res_core_min = 0.050D 
;        option.rad_res_core_max = 0.100D 
;        option.rad_res_core_exp = 1.000D 
;        option.rad_res_core_opt = 2      
    
;        option.rad_res_sol_n   = 101   
;        option.rad_res_sol_min = 0.030D
;        option.rad_res_sol_max = 0.100D   ; 0.100D
;        option.rad_res_sol_exp = 2.000D
;        option.rad_res_sol_opt = 2     
    
        xrange      = [ 1.3,2.8]
        yrange      = [-1.3,1.3]
        xpoint_zone = [1.50,2.10,-1.10,1.10]
        wall_path = '/home/ITER/lisgos/fuse_data/east/shots/default/'
        wall_file = 'main_wall.dat'
        change_sign = 0
        equ_path = '/home/ITER/lisgos/fuse_data/east/shots/92101/'
        equ_file = '92102.x4.equ'
        END
      ; ----------------------------------------------------------------
      'east_2': BEGIN
        user_step   = -0.0002D  ; -0.001D
        user_finish =  0.050D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.6D*user_finish, 0.9D*user_finish, 0.5D*user_finish, 0.4D*user_finish, -0.15D*user_finish, -0.10D*user_finish]
    
        option.ctr_minimum_length = 0.5D
    
        option.pol_res_min = 0.001D
        option.pol_res_max = 0.100D
        option.pol_res_exp = 1.0D
    
        xrange      = [ 1.3,2.8]
        yrange      = [-1.3,1.3]
        xpoint_zone = [1.50,2.10,-1.10,1.10]
        wall_path = '/home/ITER/lisgos/fuse_data/east/shots/default/'
        wall_file = 'main_wall.dat'
        change_sign = 1
        equ_path = '/home/ITER/lisgos/fuse_data/east/shots/31729/'
        equ_file = 'east_31729_4500.x4.equ'
        END
      ; ----------------------------------------------------------------
      'east_3': BEGIN
        user_step   = -0.0002D  ; -0.001D
        user_finish =  0.050D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.6D*user_finish, 0.9D*user_finish, 0.5D*user_finish, 0.4D*user_finish, -0.15D*user_finish, -0.10D*user_finish]
    
        option.ctr_minimum_length = 0.5D
    
        xrange      = [ 1.3,2.8]
        yrange      = [-1.3,1.3]
;        xrange      = [ 1.3,1.8] ; lower inner
;        yrange      = [-1.3,0.3]
;        xrange      = [ 1.3,1.5] ; debug
;        yrange      = [ 0.4,0.7]
    
        xpoint_zone = [1.50,2.10,-1.10,1.10]
        wall_path = '/home/ITER/lisgos/fuse_data/east/shots/default/'
        wall_file = 'main_wall.dat'
        change_sign = 1
        equ_path = '/home/ITER/lisgos/fuse_data/east/shots/31729/'
        equ_file = 'east_31729_6000.x4.equ'
        END
      ; ----------------------------------------------------------------
      'east_4': BEGIN
        user_step   = -0.0002D  ; -0.001D
        user_finish =  0.050D
        IF (NOT KEYWORD_SET(boundary)) THEN  $
          option.ctr_boundary = [-0.6D*user_finish, 0.9D*user_finish, 0.5D*user_finish, 0.3D*user_finish, -0.14D*user_finish, -0.10D*user_finish]
    
        option.ctr_minimum_length = 0.5D
    
        option.pol_res_min = 0.001D
        option.pol_res_max = 0.100D
        option.pol_res_exp = 1.0D
    
        xrange      = [ 1.3,2.8]
        yrange      = [-1.3,1.3]
    
        xpoint_zone = [1.50,2.10,-1.10,1.10]
        wall_path = '/home/ITER/lisgos/fuse_data/east/shots/default/'
        wall_file = 'main_wall.dat'
        change_sign = 1
        equ_path = '/home/ITER/lisgos/fuse_data/east/shots/31729/'
        equ_file = 'east_31729_7200.x4.equ'
        END
      ; ----------------------------------------------------------------
       ELSE: BEGIN
         PRINT,'ERROR grid_Input: Unknown machine_shot designation'
         PRINT,'  MACHINE_SHOT = ',machine
         RETURN, -1
         END
    ENDCASE

;   --------------------------------------------------------------------
;   --------------------------------------------------------------------
;   --------------------------------------------------------------------

  ENDELSE

  option.debug  = debug
  option.xrange = xrange
  option.yrange = yrange

  result = grid_Execute( wall_path         , $
                         wall_file         , $
                         equ_path          , $
                         equ_file          , $
                         xpoint_zone       , $
                         psi_zone          , $
                         machine           , $
                         shot              , $
                         step              , $
                         user_step         , $
                         user_finish       , $
                         change_sign       , $
                         save    = save    , $
                         preview = preview , $
                         check   = check   , $
                         all     = all     , $
                         mswin   = mswin   , $
                         limiter = limiter )

END
;
; ======================================================================
;
;
; ======================================================================
;
PRO grid_Batch, step, save=save, debug=debug

  CASE step OF
    1: BEGIN
      a=grid_Main(iter=4 ,step=step,debug=debug,save=save)
      a=grid_Main(iter=10,step=step,debug=debug,save=save)
      a=grid_Main(/mast  ,step=step,debug=debug,save=save)
      a=grid_Main(/cmod  ,step=step,debug=debug,save=save)
      a=grid_Main(/aug   ,step=step,debug=debug,save=save)
      a=grid_Main(d3d=1  ,step=step,debug=debug,save=save)  ; not working...
      a=grid_Main(d3d=2  ,step=step,debug=debug,save=save)  ; USN
      a=grid_Main(d3d=3  ,step=step,debug=debug,save=save)  ; UDN
      a=grid_Main(d3d=4  ,step=step,debug=debug,save=save)  ; LSN, OK, but tangency point on 'inside' a contour near the top is not supported, so a limitation
      a=grid_Main(d3d=5  ,step=step,debug=debug,save=save)  ; LSN (upper x-point just outside the vessel), OK
      a=grid_Main(d3d=6  ,step=step,debug=debug,save=save)  ; LSN, OK, 2nd x-oint just outside, which causes contouring near the inner wall to truncate
      a=grid_Main(/east  ,step=step,debug=debug,save=save)  
      a=grid_Main(/jet   ,step=step,debug=debug,save=save)
      END
    ELSE: BEGIN

    ;;; STEP=2, REFINEMENT PROBLEM   a=grid_Main(iter=4 ,step=step,debug=debug,save=save)
      a=grid_Main(iter=10,step=step,debug=debug,save=save)
      a=grid_Main(/mast  ,step=step,debug=debug,save=save)
      a=grid_Main(/cmod  ,step=step,debug=debug,save=save)
      a=grid_Main(/aug   ,step=step,debug=debug,save=save)
    ;  ;a=grid_Main(d3d=1  ,step=step,debug=debug,save=save) ; not working...
      a=grid_Main(d3d=2  ,step=step,debug=debug,save=save)  ; USN
    ;;;;;step2 fail;;;  a=grid_Main(d3d=3  ,step=step,debug=debug,save=save)  ; UDN
      a=grid_Main(d3d=4  ,step=step,debug=debug,save=save)  ; LSN, OK, but tangency point on 'inside' a contour near the top is not supported, so a limitation
    ;;;;;step2 fail;;;  a=grid_Main(d3d=5  ,step=step,debug=debug,save=save)  ; LSN (upper x-point just outside the vessel), OK
      a=grid_Main(d3d=6  ,step=step,debug=debug,save=save)  ; LSN, OK, 2nd x-oint just outside, which causes contouring near the inner wall to truncate
      a=grid_Main(/east  ,step=step,debug=debug,save=save)  
      a=grid_Main(/jet   ,step=step,debug=debug,save=save)
      END
  ENDCASE

; frozen 21092011
;;  a=grid_Main(iter=4 ,step=step,debug=debug,save=save)
;;  a=grid_Main(iter=10,step=step,debug=debug,save=save)
;;  a=grid_Main(/mast  ,step=step,debug=debug,save=save)
;;  a=grid_Main(/cmod  ,step=step,debug=debug,save=save)
;;  a=grid_Main(/aug   ,step=step,debug=debug,save=save)
;;  ;a=grid_Main(d3d=1  ,step=step,debug=debug,save=save) ; not working...
;;  a=grid_Main(d3d=2  ,step=step,debug=debug,save=save)  ; USN
;;;;;;step2 fail;;;  a=grid_Main(d3d=3  ,step=step,debug=debug,save=save)  ; UDN
;  a=grid_Main(d3d=4  ,step=step,debug=debug,save=save)  ; LSN, OK, but tangency point on 'inside' a contour near the top is not supported, so a limitation
;;;;;;step2 fail;;;  a=grid_Main(d3d=5  ,step=step,debug=debug,save=save)  ; LSN (upper x-point just outside the vessel), OK
;  a=grid_Main(/east  ,step=step,debug=debug,save=save)  
;  a=grid_Main(/jet   ,step=step,debug=debug,save=save)

END
;
; ======================================================================
;
PRO grid_Main, args

  IF (NOT KEYWORD_SET(args)) THEN args = 'none'

  PRINT, N_ELEMENTS(args),args

  machine = args[1]
  shot    = args[2]
  equ     = args[3]


  PRINT,'PREVIEW = ',args[0]
  PRINT,'MACHINE = ',machine
  PRINT,'SHOT    = ',shot
  PRINT,'EQU     = ',equ



  IF (args[0] EQ '1') THEN BEGIN
    a=grid_Run(machine=machine,shot=shot,equ=equ,/debug,/preview)
  ENDIF ELSE BEGIN
    a=grid_Run(machine=machine,shot=shot,equ=equ,/debug,/save,/all,/check)
  ENDELSE

  dummy = ' '
  READ,"Press a key and then ENTER to quit: ",dummy

  EXIT

END
