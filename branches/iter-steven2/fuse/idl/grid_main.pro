;
; ======================================================================
;
; improvements:
;	-the D3D cases that don't work, see below (the C-Mod case as well)
;	-the MAST problem
;	-the TS problems -- see below
;
;
FUNCTION grid_Main, iter=iter, mast=mast, cmod=cmod, west=west, jet=jet, aug=aug, d3d=d3d, east=east,  $
                    debug=debug,  $
                    save=save,    $
                    step=step

  IF (KEYWORD_SET(iter)) THEN machine = 'iter' 
  IF (KEYWORD_SET(mast)) THEN machine = 'mast' 
  IF (KEYWORD_SET(cmod)) THEN machine = 'cmod' 
  IF (KEYWORD_SET(west)) THEN machine = 'west' 
  IF (KEYWORD_SET(jet )) THEN machine = 'jet'  
  IF (KEYWORD_SET(aug )) THEN machine = 'aug'  
  IF (KEYWORD_SET(d3d )) THEN machine = 'd3d'  
  IF (KEYWORD_SET(east)) THEN machine = 'east' 

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
      PRINT,'ERROR grid_Main: MACHINE not set'
      RETURN,-1
      END
  ENDCASE

  IF (NOT KEYWORD_SET(step)) THEN step = 1

  machine = machine + "_" + STRTRIM(STRING(shot),2)

; ----------------------------------------------------------------------
; ----------------------------------------------------------------------
; ----------------------------------------------------------------------

  grid_SetupParameters

; orange - o-point
; purple - primary x-point
; green  - secondary x-point
  
  CASE machine OF
    ; ----------------------------------------------------------------
    'iter_4': BEGIN
      user_step   = -0.003D ; -0.003D
      user_finish =  4.0D ; 3.0D
      xrange      = [ 3.5,8.75]     
      yrange      = [-5.0,5.0]
      xpoint_zone = [4.0,7.0,-3.8,5.5]
      path_wall = '/home/ITER/lisgos/divimp/shots/iter/1514/'
      file_wall = 'psi_wall_simple.dat'
      change_sign =  0
      path_equ = '/home/ITER/lisgos/divimp/shots/iter/equilibria/'
      file_equ = 'feat_001.x4.equ'
      END
    ; ----------------------------------------------------------------
    'iter_10': BEGIN
      user_step   = -0.003D; -0.003D
      user_finish =  3.0D
      xrange      = [ 3.5,8.75]     
      yrange      = [-5.0,5.0]
      xpoint_zone = [4.0,7.0,-3.8,5.5]
      path_wall = '/home/ITER/lisgos/divimp/shots/iter/1514/'
      file_wall = 'psi_wall_simple.dat'
      change_sign =  0
      path_equ = '/home/ITER/lisgos/divimp/shots/iter/1514/'
      file_equ = 'Baseline2008-li0.70.x4.equ'
      END
    ; ----------------------------------------------------------------    
    'mast_1': BEGIN
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
    'cmod_1': BEGIN
      user_step   = -0.0005D0
      user_finish =  0.025D
      xrange      = [0.3,1.2]
      yrange      = [-0.7,0.7]
;      xrange      = [0.5,0.8]   ; upper x-point
;      yrange      = [0.2,0.7]
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
    'west_1': BEGIN
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
    'jet_1': BEGIN
      user_step   = -0.001D
      user_finish =  0.80D
      xrange      = [ 1.6, 4.1]  ; full
      yrange      = [-2.0, 2.3]
      xrange      = [ 2.3, 3.0]  ; divertor
      yrange      = [-1.8,-1.2]
      xpoint_zone = [2.4,3.1,-1.8,2.0]
      path_wall = '/home/ITER/lisgos/divimp/shots/jet/default/'
      file_wall = 'vessel_wall.dat'
      change_sign = 1
      path_equ = '/home/ITER/lisgos/divimp/shots/jet/68124_49000/'
      file_equ = 'JET_68124_49000.X2.equ' ; LSND
      END
    ; ----------------------------------------------------------------
    'aug_1': BEGIN
      user_step   = -0.0005D
      user_finish =  0.10D
      xrange      = [ 0.9, 2.4]  ; full
      yrange      = [-1.4, 1.4]
;      xrange      = [ 1.2, 1.7]  ; lower divertor
;      yrange      = [-1.3,-0.8]
      xpoint_zone = [ 1.30,1.75,-1.20,1.35]
      path_wall = '/home/ITER/lisgos/divimp/shots/aug/default/'
      file_wall = 'vessel_wall.dat'
      change_sign = 0
      path_equ  = '/home/ITER/lisgos/divimp/shots/aug/22575_3000/' 
      file_equ  = '22575_72x18.sonnet.equ'
      END
    ; ----------------------------------------------------------------
    'd3d_1': BEGIN
      user_step   = -0.0001D ; -0.0005D
      user_finish =  0.08D ; 0.10D
      xrange = [ 0.8, 2.5]  ; full
      yrange = [-1.5, 1.5]
;      xrange = [ 1.1, 1.8]  ; upper PFR
;      yrange = [ 0.8, 1.4]
;      xrange = [ 1.1, 2.2]  ; upper half
;      yrange = [ 0.2, 1.4]
;      xrange      = [ 1.5, 1.9]  ; lower divertor
;      yrange      = [-1.4,-1.0]
      xpoint_zone = [ 1.2,1.9,-1.38,1.36]
      path_wall = '/home/ITER/lisgos/divimp/shots/d3d/default/'
      file_wall = 'vessel_wall2.dat'
      change_sign = 1
      ; ERROR grid_GetOrthogonal: Normal vector requested away from vertex on outer wall
      path_equ  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
;      file_equ  = 'g105500.03500.x2.equ' ; LSND, 2nd x-point just inside the vessel
      file_equ  = 'g105500.03500.x8.equ' ; LSND, 2nd x-point just inside the vessel
      END
    ; ----------------------------------------------------------------
    'd3d_2': BEGIN
      user_step   = -0.0005D
      user_finish =  0.10D
      xrange = [ 0.8, 2.5]  ; full
      yrange = [-1.5, 1.5]
;     xrange = [ 1.0, 1.6]  ; upper PFR
;     yrange = [ 1.1, 1.4]
      xpoint_zone = [ 1.2,1.9,-1.38,1.36]
      path_wall = '/home/ITER/lisgos/divimp/shots/d3d/default/'
      file_wall = 'vessel_wall2.dat'
      change_sign = 1
      ; works...
      path_equ  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
      file_equ  = 'g131245.02800.x2.equ' ; USND
      END
    ; ----------------------------------------------------------------
    'd3d_3': BEGIN
      user_step   = -0.0005D
      user_finish =  0.10D
      xrange = [ 0.8, 2.5]  ; full
      yrange = [-1.5, 1.5]
;      xrange = [ 1.0, 1.4]  ; upper PFR
;      yrange = [ 1.1, 1.4]
      xpoint_zone = [ 1.2,1.9,-1.38,1.36]
      path_wall = '/home/ITER/lisgos/divimp/shots/d3d/default/'
      file_wall = 'vessel_wall2.dat'
      change_sign = 1
      ; works...
      path_equ  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
      file_equ  = 'g134585.03005.x4.equ' ; UDND
      END
    ; ----------------------------------------------------------------
    'd3d_4': BEGIN
      user_step   = -0.0005D
      user_finish =  0.10D
      xrange = [ 0.8, 2.5]  ; full
      yrange = [-1.5, 1.5]
;      xrange = [ 1.0, 1.4]  ; upper PFR
;      yrange = [ 1.1, 1.4]
      xpoint_zone = [ 1.2,1.9,-1.38,1.36]
      path_wall = '/home/ITER/lisgos/divimp/shots/d3d/default/'
      file_wall = 'vessel_wall2.dat'
      change_sign = 1

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

      path_equ  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
      file_equ  = 'g135487.03865.x4.equ' ; LSND
      END
    ; ----------------------------------------------------------------
    'd3d_5': BEGIN
      user_step   = -0.0005D
      user_finish =  0.10D
      xrange = [ 0.8, 2.5]  ; full
      yrange = [-1.5, 1.5]
;      xrange = [ 1.0, 1.4]  ; upper PFR
;      yrange = [ 1.1, 1.4]
      xpoint_zone = [ 1.2,1.9,-1.38,1.36]
      path_wall = '/home/ITER/lisgos/divimp/shots/d3d/default/'
      file_wall = 'vessel_wall2.dat'
      change_sign = 1
      ; fails... but for the same reason as the above equilibrium, although this time
      ; for the upper inner plenum entrance...
      path_equ  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
      file_equ  = 'g123417.03500.x2.equ' ; LSND, upper point just outside the vessel
      END
    ; ----------------------------------------------------------------
    'd3d_6': BEGIN
      user_step   = -0.0005D
      user_finish =  0.10D
      xrange = [ 0.8, 2.5]  ; full
      yrange = [-1.5, 1.5]
;      xrange = [ 1.0, 1.4]  ; upper PFR
;      yrange = [ 1.1, 1.4]
      xpoint_zone = [ 1.1,1.9,-1.38,1.36]
      path_wall = '/home/ITER/lisgos/divimp/shots/d3d/default/'
      file_wall = 'vessel_wall3.dat'
      change_sign = 1
      path_equ  = '/home/ITER/lisgos/divimp/shots/d3d/146198_3010/'
      file_equ  = 'g146198.03010.x8.equ'  ; screws up on the inside
      END
    ; ----------------------------------------------------------------
    'east_1': BEGIN
      user_step   = -0.001D
      user_finish =  0.1000D
      xrange      = [ 1.3,2.8]
      yrange      = [-1.3,1.3]
      xpoint_zone = [1.50,2.10,-1.10,1.10]
      path_wall = '/home/ITER/lisgos/fuse_data/east/shots/92101/'
      file_wall = 'main_wall.dat'
      change_sign = 0
      path_equ = '/home/ITER/lisgos/fuse_data/east/shots/92101/'
      file_equ = '92102.equ'
      END
    ; ----------------------------------------------------------------
     ELSE: BEGIN
       PRINT,'ERROR grid_Main: Unknown machine_shot designation'
       PRINT,'  MACHINE_SHOT = ',machine
       RETURN, -1
       END
  ENDCASE

; ----------------------------------------------------------------------
; ----------------------------------------------------------------------
; ----------------------------------------------------------------------

  CASE step OF 
;   --------------------------------------------------------------------
;   Radial processing:
;   --------------------------------------------------------------------
    1: BEGIN
      wall = grid_LoadWall(path_wall+file_wall,/debug,xrange=xrange,yrange=yrange)

      b = grid_ReadEQUFile(path_equ+file_equ,CHANGE_SIGN=change_sign)
          ;shade_surf,b.psi,b.x,b.y; ,ax=0  
      b = grid_FindNullPoints(b,xpoint_zone,1,/debug)
      
      b = grid_RefineSeparatrices(b)  
      ;b = grid_RefineSeparatrices(b, /debug, xrange=xrange, yrange=yrange)  

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

      ;CONTOUR, b.psi, b.x, b.y, NLEVELS=20, c_labels=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

      print,'null_x',b.x[b.null_i]
      print,'null_y',b.y[b.null_j]

      b = grid_AnalyseBoundary(b,wall,user_step,user_finish, machine,  $
                               debug=debug,xrange=xrange,yrange=yrange,save=save)
      END
;   --------------------------------------------------------------------
;   Load in the saved file from the radial processing and do all the 
;   poloidal processing.
;   --------------------------------------------------------------------
    2: BEGIN
   
      RESTORE,'stored_step1_'+machine+'.sav'
      contour_array.contour1.separatrix = 1
   
      help,scan_params,/struct
   
      grid_ZoneWall, wall.pt1, wall.pt2, debug=debug, xrange=xrange, yrange=yrange
   
      PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
   
      ; Check if the secondary separatrix needs to be added:
      IF (scan_params.process_2nd NE -1) THEN BEGIN
;        xrange = [4.5,5.5]
;        yrange = [4.1,5.0]
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
   
;      grid_TrimContours, b, contour_array, wall,  $
;                         debug=debug, xrange=xrange, yrange=yrange
   
      contour_array = grid_InstallXPoints(b,scan_params,contour_array,  $
                                          debug=debug, xrange=xrange, yrange=yrange)
   
      print, '---------------------------------------------'
      print, '--- first wall update -----------------------'
      print, '---------------------------------------------'
      grid_UpdateWall, b, contour_array, wall,  $
                       debug=debug, xrange=xrange, yrange=yrange
;      grid_Debug,b,contour_array,wall,xrange,yrange
      print, '--------------------------------------------------'
      print, '--- process wall ---------------------------------'
      print, '--------------------------------------------------'
      grid_ProcessWall, b, wall, scan_params, contour_array, $
                        debug=debug, xrange=xrange, yrange=yrange
      print, '--------------------------------------------------'
      print, '--- subsequent wall update -----------------------'
      print, '--------------------------------------------------'
      grid_UpdateWall, b, contour_array, wall, $
                       debug=debug, xrange=xrange, yrange=yrange
   
   
      contour_array = grid_BuildRadialMap(contour_array, wall,  $
                                          debug=debug, xrange=xrange, yrange=yrange)



      contour_array = grid_RadialRefinement(contour_array, wall, b,  $
                                            debug=debug, xrange=xrange, yrange=yrange)

;print,'diggin it'
;stop
   
      contour_array = grid_SetPoloidalDomains(contour_array, wall, b, $
                                              debug=debug, xrange=xrange, yrange=yrange)
   
   
      contour_array = grid_CreatePoloidalPoints(contour_array, wall, b,  $
                                              debug=debug, xrange=xrange, yrange=yrange)
   

      IF (KEYWORD_SET(save)) THEN  $
        SAVE,filename='stored_step2_'+machine+'.sav',b,wall,contour_array,scan_params

      print,'taking a rest now'
      RETURN, 1
      stop
      grid_Debug,b,contour_array,wall,xrange,yrange

      END
;   --------------------------------------------------------------------
;   Take the radial and poloidal descretisation and build up the
;   rings/tubes and cells.
;   --------------------------------------------------------------------
    3: BEGIN

      RESTORE,'stored_step2_'+machine+'.sav'
;      help,contour_array,/struct
;      help,contour_array.contour1,/struct
;      print,contour_array.contour1.map_in
;      print,contour_array.contour1.map_out
   
      IF (KEYWORD_SET(debug)) THEN BEGIN
        tags = STRUPCASE(TAG_NAMES(contour_array))
        nctr = N_ELEMENTS(tags)
        PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
        OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
        FOR i = 0, nctr-1 DO BEGIN
          ctr = grid_ExtractStructure(contour_array,tags[i])      
          OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
          OPLOT,ctr.grid_x,ctr.grid_y,color=Truecolor('Pink'),PSYM=6  
          npts = N_ELEMENTS(ctr.x)
          XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('White')
        ENDFOR
      ENDIF
   
      contour_array = grid_AddCoreRings(  $
        contour_array, wall, b, debug=debug, xrange=xrange, yrange=yrange)
   
      tubes = grid_BuildGrid(  $
        contour_array, wall, b, debug=debug, xrange=xrange, yrange=yrange)

      IF (KEYWORD_SET(save)) THEN  $
        SAVE,filename='stored_step3_'+machine+'.sav',b,wall,contour_array,scan_params,tubes

;      grid_Debug,b,contour_array,wall,xrange,yrange

      RETURN, 1
      END
;   --------------------------------------------------------------------
;   Generate the grid file:
;   --------------------------------------------------------------------
    4: BEGIN

      RESTORE,'stored_step3_'+machine+'.sav'
;      help,contour_array,/struct
;      help,contour_array.contour1,/struct
;      print,contour_array.contour1.map_in
;      print,contour_array.contour1.map_out
   
      IF (KEYWORD_SET(debug)) THEN BEGIN
        tags = STRUPCASE(TAG_NAMES(contour_array))
        nctr = N_ELEMENTS(tags)
        PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
        OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
        FOR i = 0, nctr-1 DO BEGIN
          ctr = grid_ExtractStructure(contour_array,tags[i])      
          OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
          OPLOT,ctr.grid_x,ctr.grid_y,color=Truecolor('Pink'),PSYM=3
          npts = N_ELEMENTS(ctr.x)
          XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('White')
        ENDFOR
      ENDIF
   
      help,tubes,/struct

      tags  = STRUPCASE(TAG_NAMES(tubes))
      ntube = N_ELEMENTS(tags)
      FOR itube = 1, ntube DO BEGIN
        tube = grid_ExtractStructure(tubes,tags[itube-1])
        FOR icell = 0, tube.n-1 DO BEGIN
          OPLOT,[tube.x1[icell],tube.x2[icell],tube.x2[icell+1],tube.x1[icell+1],tube.x1[icell]],  $
                [tube.y1[icell],tube.y2[icell],tube.y2[icell+1],tube.y1[icell+1],tube.y1[icell]],  $
                COLOR=Truecolor('Lightblue')
        ENDFOR
      ENDFOR

      fname = 'stored_step4_'+machine+'.grid'

      result = grid_WriteGridFile(fname, tubes, contour_array, wall, b,  $ 
                                  debug=debug, xrange=xrange, yrange=yrange)

      RETURN, 1
      END
;   --------------------------------------------------------------------
;   --------------------------------------------------------------------
    ELSE: BEGIN
      PRINT,'ERROR grid_Main: Unknown step'
      PRINT,' STEP =',step
      RETURN,-1
      END
;   --------------------------------------------------------------------
;   --------------------------------------------------------------------
  ENDCASE

  result = 1

  RETURN, result
END
;
; ======================================================================
;
PRO grid_Batch

  debug = 1
  save  = 1
  step  = 2

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

