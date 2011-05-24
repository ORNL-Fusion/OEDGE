;
; ======================================================================
;
; improvements:
;	-the D3D cases that don't work, see below (the C-Mod case as well)
;	-the MAST problem
;	-the TS problems -- see below
;
;
FUNCTION grid_Main, iter=iter, mast=mast, cmod=cmod, west=west, jet=jet, aug=aug, d3d=d3d,  $
                    load=load,  $
                    debug=debug,  $
                    save=save

  IF (KEYWORD_SET(iter)) THEN machine = 'iter'
  IF (KEYWORD_SET(mast)) THEN machine = 'mast'
  IF (KEYWORD_SET(cmod)) THEN machine = 'cmod'
  IF (KEYWORD_SET(west)) THEN machine = 'west'
  IF (KEYWORD_SET(jet )) THEN machine = 'jet' 
  IF (KEYWORD_SET(aug )) THEN machine = 'aug' 
  IF (KEYWORD_SET(d3d )) THEN machine = 'd3d' 

  IF (NOT KEYWORD_SET(machine)) THEN machine = 'iter'


  grid_SetupParameters


; orange - o-point
; purple - primary x-point
; green  - secondary x-point

  
  CASE machine OF
    ; ----------------------------------------------------------------
    'iter': BEGIN
      user_step   = -0.003D; -0.003D
      user_finish = 3.0D

      xrange = [ 4.02,4.06]
      yrange = [ 1.5,4.0]
      xrange = [ 4.02,4.06]
      yrange = [-1.5,1.5]
      xrange = [ 4.0, 6.0]  ; divertor
      yrange = [-4.7,-3.0]
      ;xrange = [ 4.3, 4.6]
      ;yrange = [-4.0,-3.8]

      xrange = [ 6.52,6.65] ; <-- contours are unnecessarily close together here, so add a check that 
      yrange = [ 3.84,3.91]     ; avoids this -- note that a warning already appears on the screen
      ;xrange = [ 5.6,5.9]      ; pretty easy to just flag this and cause the grid_AddContour call to fail
      ;yrange = [ 4.4,4.7]      ; But then again, the rings appear to be well enough resolved, and you will lose
                               ; the well defined gap... think some more.

      xrange = [ 8.2,8.5]  ; notice how the rings don't get into the gaps since they die before then because
      yrange = [-1.0,2.0]  ; they are pretty short -- adjust the shortness before stopping parameter?
                           ; They could also be hitting the outer PSI boundary, but I don't think so (but didn't check)

      xrange = [ 6.0, 7.5] ; some contours close to each other again     
      yrange = [-3.5,-1.5] ; to I need to got to an x8 equilibrium -- makes a difference?

      xrange = [ 4.8,5.1]  ; top
      yrange = [ 4.4,4.9]
      xrange = [ 4.85,4.95]  ; top
      yrange = [ 4.65,4.80]
      ;xrange = [ 4.0,4.95] ; top
      ;yrange = [ 4.3,5.0]

      xrange = [ 4.02,4.06] ; inner wall
      yrange = [ 1.0,2.0]
      xrange = [ 4.035,4.045] ; inner wall lower
      yrange = [-2.05,-1.95]
      xrange = [ 3.5,8.75]     
      yrange = [-5.0,5.0]
      xpoint_zone = [4.0,7.0,-3.8,5.5]
      path_wall = '/home/ITER/lisgos/divimp/shots/iter/1514/'
      file_wall = 'psi_wall_simple.dat'
      change_sign = 0
      ; doesn't appear to be catching the small gap at the top properly anymore
      path_equ = '/home/ITER/lisgos/divimp/shots/iter/equilibria/'
      file_equ = 'feat_001.x4.equ'
      path_equ = '/home/ITER/lisgos/divimp/shots/iter/1514/'
      file_equ = 'Baseline2008-li0.70.x4.equ'
      END
    ; ----------------------------------------------------------------    
    'mast': BEGIN
      ; Gets stuck on inner wall when the equilibrium data peters out and the contour doesn't
      ; make it to the tangency point boundary.
      ; Need to recognize this somehow and switch to using the wall intersection rather than
      ; the tangecy point boundary.
      user_step   = -0.0005D
      user_finish =  0.03D
      xrange      = [0.05 ,2.2]
      yrange      = [-2.5 ,2.5]
      xpoint_zone = [0.4,1.0,-1.5,1.5]
      path_wall = '/home/ITER/lisgos/divimp/shots/mast/default/'
      file_wall = 'main_wall.dat'
      change_sign = 0
      ; trouble with multiple target intersections when far out
      path_equ = '~/fuse_data/mast/shots/24867/'
      file_equ = '24867_335.equ'
      ; Fails on the inner wall because of limited efit data
      ;path_equ = '~/divimp/shots/mast/24860/'
      ;file_equ = 'carre.24860_240.equ'
      ; works, but get into trouble if going to far out(>0.1) as multiple wall interesections appear
      ; can perhaps have a selection rule based on the appearance of an odd number of intersections, i.e.
      ; 3 rather than 4...
      path_equ = '/home/ITER/lisgos/divimp/shots/mast/13018_250/'
      file_equ = 'sonnet_13018_250.equ'
      END
    ; ----------------------------------------------------------------
    'cmod': BEGIN
      user_step   = -0.0005D0
      user_finish =  0.025D
      IF (KEYWORD_SET(load)) THEN BEGIN
        xrange      = [0.439,0.441]
        yrange      = [-0.05, 0.05]
        xrange      = [0.450,0.490]
        yrange      = [-0.5,-0.45]
        ;xrange      = [0.600,0.650]
        ;yrange      = [-0.55,-0.50]
        ;xrange      = [0.600,0.700]
        ;yrange      = [-0.50,-0.40]
        xrange      = [0.565,0.575]
        yrange      = [0.50,0.52]
        xrange      = [0.43,0.50]
        yrange      = [-0.1,0.1]
        xrange      = [0.50,0.60]  ; Upper x-point
        yrange      = [0.35,0.45]
        xrange      = [ 0.50, 0.60]  ; Lower x-point
        yrange      = [-0.45,-0.35]
        xrange      = [0.4,0.7] ; divertor
        yrange      = [-0.6,-0.3]
        xrange      = [0.4,1.1]
        yrange      = [-0.65,0.65]
      ENDIF ELSE BEGIN
        xrange      = [0.3,1.2]
        yrange      = [-0.7,0.7]
      ENDELSE
      xpoint_zone = [0.52,0.75,-0.45,0.45]
      path_wall = '/home/ITER/lisgos/divimp/shots/cmod/1100303017_0138/'
      file_wall = 'vessel_wall4.dat'  ; no pump plenum entrance
;      file_wall = 'vessel_wall.dat'  
      change_sign = 1
      ; works... but need to re-address the problem of "funny business" in the PFZ, see vessel_wall4.dat file
      path_equ  = '/home/ITER/lisgos/divimp/shots/cmod/1100303017_0138/'
      file_equ  = 'cmod.1100303017.01380.x4.equ'
      END
    ; ----------------------------------------------------------------
    'west': BEGIN
      user_step   = -0.01D
      user_finish =  0.80D
      xrange      = [1.4,3.3]
      yrange      = [-1.0,1.0]
      xpoint_zone = [2.15,2.7,-0.73,0.75]
      path_wall = '/home/ITER/lisgos/divimp/shots/ts2/upgrade/'
      file_wall = 'vessel_wall.dat'
      change_sign = 0
      ; dies on corner - might be same as D3D problem
       path_equ = '/home/ITER/lisgos/divimp/shots/ts2/600kA_LN_1cm/'
       file_equ = 'TS2_600kA_LN.x2.cr.equ'
      ; gets pissed off by the second x-point that's just inside the vessel
      ; path_equ = '/home/ITER/lisgos/divimp/shots/ts2/1MA_LN_1cm/'
      ; file_equ = 'TS2_1MA_LN.x2.equ'
      ; pissed off by CDNness... to be expected 
      ;file_equ = '/home/ITER/lisgos/divimp/shots/ts2/600kA_CDN/TS2_600kA.x2.equ'
      ;file_equ = '/home/ITER/lisgos/divimp/shots/ts2/600kA_CDN/TS2_600kA.x2.equ'
      END
    ; ----------------------------------------------------------------
    'jet': BEGIN
      user_step   = -0.001D
      user_finish =  0.80D
      xrange      = [1.6,4.1]
      yrange      = [-2.0,2.3]
      xpoint_zone = [2.4,3.1,-1.8,2.0]
      path_wall = '/home/ITER/lisgos/divimp/shots/jet/default/'
      file_wall = 'vessel_wall.dat'
      change_sign = 1
      ; Gets stuck at the nasty tangency point at the outer midplane.
      ; Need to recognize this and force a step, even though the wall isn't followed
      ; perfectly.
      path_equ = '/home/ITER/lisgos/divimp/shots/jet/68124_49000/'
      file_equ = 'JET_68124_49000.X2.equ' ; LSND
      END
    ; ----------------------------------------------------------------
    'aug': BEGIN
      user_step   = -0.0005D
      user_finish =  0.10D
      xrange      = [ 1.2, 1.7]  ; divertor
      yrange      = [-1.3,-0.8]
      xrange      = [ 1.285, 1.287]  ; super close up
      yrange      = [-0.91,-0.89]
      ;xrange      = [ 1.45,1.55]
      ;yrange      = [ 1.15,1.25]
      xrange      = [ 2.05,2.15] ; top of outer limiter
      yrange      = [ 0.40,0.60]
      xrange      = [ 1.1,1.3]  ; top of inner wall
      yrange      = [ 0.6,1.2]
      xrange      = [ 1.15,1.17]  ; top of inner wall - zoom
      yrange      = [ 0.725,0.775]
      xrange      = [ 0.9,2.4]  ; full
      yrange      = [-1.4,1.4]
      xpoint_zone = [ 1.30,1.75,-1.20,1.35]
      path_wall = '/home/ITER/lisgos/divimp/shots/aug/default/'
      file_wall = 'vessel_wall.dat'
      change_sign = 0
      ; seems to work
      path_equ  = '/home/ITER/lisgos/divimp/shots/aug/22575_3000/'
      file_equ  = '22575_72x18.sonnet.equ'
      END
    ; ----------------------------------------------------------------
    'd3d': BEGIN
      user_step   = -0.0005D
      user_finish =  0.10D
      IF (KEYWORD_SET(load)) THEN BEGIN
        xrange = [ 1.0,1.5]
        yrange = [ 1.0,1.5]   ; *** check more d3d equilibiria!!! ***
        xrange = [ 1.0,1.20]
        yrange = [ 0.8,1.6]
        ;xrange = [ 1.0,1.60]  ; fix the poor curvature on the dome that causes extra rings
        ;yrange = [ 0.8,1.6]   ; fix the incorrect PFZ assignment in the upper plenum entrance (need to 
        ;xrange = [ 1.6,1.8]  ; replace with an OSM-like scheme: o-point, x-point lines + radial connection map information
        ;yrange = [-1.5,-1.0]
        xrange = [ 1.6,1.8]
        yrange = [-1.4,-1.3]
        xrange = [ 1.2,1.5]
        yrange = [ 1.2,1.4]
        xrange = [ 1.0,1.4]  ; Upper inner corner
        yrange = [ 1.1,1.4]
        xrange = [ 0.8,2.5]
        yrange = [-1.5,1.5]
      ENDIF ELSE BEGIN
        xrange = [ 0.8,2.5]
        yrange = [-1.5,1.5]
      ENDELSE
      xpoint_zone = [ 1.2,1.9,-1.38,1.36]
      path_wall = '/home/ITER/lisgos/divimp/shots/d3d/default/'
      file_wall = 'vessel_wall2.dat'
      change_sign = 1
      ; works...
      ;path_equ  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
      ;file_equ  = 'g105500.03500.x2.equ' ; LSND, 2nd x-point just inside the vessel
      ; fails in lower left corner when calculating cross-field contours...
      path_equ  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
      file_equ  = 'g131245.02800.x2.equ' ; USND
      ; works...
      ;path_equ  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
      ;file_equ  = 'g134585.03005.x4.equ' ; UDND
      ; Corner at the entrance to the upper plenum causes problems, although the 
      ; code seems to cope, sort of, by acting defensively.  This is related to an 
      ; OUTSTANDING ISSUE where there are multiple valid segments along a contour,
      ; which can happen if the angle of the contours relative to the surface is such
      ; that a tangency point is not associated with the entrance to the secondary
      ; region of the vaccum vessel -- the pump plenums in this case and for C-Mod
      ; Perhaps filter the contour when an unusual number of intersections are found, 
      ; that is, more than two, and try to decide which to keep...

      ; I THINK THIS PROBLEM IS PARTLY SOLVED ACTUALLY...
      ; From the behaviour of MAST and JET runs... and even this D3D run.  What I'm not 
      ; clear on is why the D3D case stops expanding radially when this issue is encountered.

      ;path_equ  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
      ;file_equ  = 'g135487.03865.x4.equ' ; LSND
      ; fails... but for the same reason as the above equilibrium, although this time
      ; for the upper inner plenum entrance...
      ;path_equ  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
      ;file_equ  = 'g123417.03500.x2.equ' ; LSND, upper point just outside the vessel
      END
    ; ----------------------------------------------------------------
     ELSE: BEGIN
       stop
       END
  ENDCASE




  IF (KEYWORD_SET(load)) THEN BEGIN
;  IF (1 EQ 1 AND process_2nd) THEN BEGIN

    RESTORE,'stored_scan_'+machine+'.sav'
    contour_array.contour1.separatrix = 1

    help,scan_params,/struct

    grid_ZoneWall, wall.pt1, wall.pt2, debug=debug, xrange=xrange, yrange=yrange

    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1

    ; Check if the secondary separatrix needs to be added:
    IF (scan_params.process_2nd NE -1) THEN BEGIN
;      xrange = [4.5,5.5]
;      yrange = [4.1,5.0]
      ; Add the secondary separatrix:
      status = 1
      IF (scan_params.failure_2nd_pfz EQ 1) THEN BEGIN
        psi_shift =  0.000001D 
        ; Check if the secondary x-point is outside the vessel:
        focus_x = b.x[b.null_i[2]]
        focus_y = b.y[b.null_j[2]]
        inside = grid_PointInPolygon(focus_x,focus_y,wall.x,wall.y)
        IF (NOT inside) THEN BEGIN
          ; Get the distance from the 2nd null to the wall:
          dist1 = grid_PerpDistance(focus_x,focus_y,wall.x,wall.y)
          ; Check if the there's a wall tangency point that's very close to the
          ; secondary x-point, and if yes, then don't add contours associated
          ; with the 2nd x-point:    
          tags  = STRUPCASE(TAG_NAMES(contour_array))
          nctr  = N_ELEMENTS(tags)
          FOR ictr = 1, nctr DO BEGIN
            ctr = grid_ExtractStructure(contour_array,tags[ictr-1])      
            IF (ctr.tangent_i EQ -1) THEN CONTINUE
            dist2 = SQRT( (focus_x-ctr.tangent_p1[0])^2 + (focus_y-ctr.tangent_p1[1])^2 )
            IF (dist2 LT 10.0D*dist1) THEN status = 0
            print,'check 2nd',dist1,dist2,status
          ENDFOR
        ENDIF
      ENDIF ELSE  $
        psi_shift = -0.000001D

      IF (status) THEN  $
        status = grid_AddContour(b, wall, scan_params, contour_array, -1,  $
                                 psi_val=b.psi_2nd_xpoint + psi_shift,  $
                                 debug=debug, xrange=xrange, yrange=yrange)
    ENDIF

    grid_TrimContours, b, contour_array, wall,  $
                       debug=debug, xrange=xrange, yrange=yrange

    contour_array = grid_InstallXPoints(b,scan_params,contour_array,  $
                                        debug=debug, xrange=xrange, yrange=yrange)

    print, '---------------------------------------------'
    print, '--- first wall update -----------------------'
    print, '---------------------------------------------'
    grid_UpdateWall, b, contour_array, wall,  $
                     debug=debug, xrange=xrange, yrange=yrange

    print, '--------------------------------------------------'
    print, '--- process wall ---------------------------------'
    print, '--------------------------------------------------'

    grid_ProcessWall, b, wall, scan_params, contour_array, $
                      debug=debug, xrange=xrange, yrange=yrange

    print, '--------------------------------------------------'
    print, '--- subsequent wall update -----------------------'
    print, '--------------------------------------------------'
    grid_UpdateWall, b, contour_array, wall,  $
                     debug=debug, xrange=xrange, yrange=yrange


; add check to amek sure that the contours don intersect the walls except at the ends, other wise there
; which can happen near insufficiently resolved tangency points


;    help,contour_array,/struct
;   xrange = [6.4,7.0]
;   yrange = [3.7,4.0]


    contour_array = grid_BuildRadialMap(contour_array, wall,  $
                                        debug=debug, xrange=xrange, yrange=yrange)

    contour_array = grid_SetPoloidalDomains(contour_array, wall, b, $
                                            debug=debug, xrange=xrange, yrange=yrange)


    contour_array = grid_CreatePoloidalPoints(contour_array, wall,  $
                                            debug=debug, xrange=xrange, yrange=yrange)

stop

    grid_Debug,b,contour_array,wall,xrange,yrange

  ENDIF 









  wall = grid_LoadWall(path_wall+file_wall,/debug,xrange=xrange,yrange=yrange)

  b = grid_ReadEQUFile(path_equ+file_equ,CHANGE_SIGN=change_sign)
      ;shade_surf,b.psi,b.x,b.y; ,ax=0  
  b = grid_FindNullPoints(b,xpoint_zone,1,/debug)

  b = grid_RefineSeparatrices(b)  
;  b = grid_RefineSeparatrices(b, /debug, xrange=xrange, yrange=yrange)  

  IF (0 EQ 1) THEN BEGIN
    plot,wall.x,wall.y,color=Truecolor('Yellow')
    ; help,b,/struct
    OPLOT, [b.x[b.null_i[0]]], [b.y[b.null_j[0]]], PSYM=6, color=TrueColor('Red')      ; o-point
    OPLOT, [b.x[b.null_i[1]]], [b.y[b.null_j[1]]], PSYM=6, color=TrueColor('Green')   ; primary   x-point
    IF (N_ELEMENTS(b.null_i) EQ 3) THEN  $
      OPLOT, [b.x[b.null_i[2]]], [b.y[b.null_j[2]]], PSYM=6, color=TrueColor('Blue' ) ; secondary x-point
    OPLOT, [xpoint_zone[0],xpoint_zone[0],xpoint_zone[1],xpoint_zone[1],xpoint_zone[0]],  $
          [xpoint_zone[2],xpoint_zone[3],xpoint_zone[3],xpoint_zone[2],xpoint_zone[2]],color=TrueColor('Red')
    stop
  ENDIF

;  CONTOUR, b.psi, b.x, b.y, NLEVELS=20, c_labels=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
;stop



  b = grid_AnalyseBoundary(b,wall,user_step,user_finish, machine,  $
                           debug=debug,xrange=xrange,yrange=yrange,save=save)








  result = 'wow'


  RETURN, result
END
;
; ======================================================================
;
PRO grid_Batch

  a=grid_Main(/mast,/debug,/save)
  a=grid_Main(/cmod,/debug,/save)
  a=grid_Main(/iter,/debug,/save)
  a=grid_Main(/aug ,/debug,/save)
  a=grid_Main(/d3d ,/debug,/save)

END
;
; ======================================================================
;

