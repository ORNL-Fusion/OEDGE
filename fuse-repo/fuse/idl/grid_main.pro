; a=grid_Run(d3d=9,step=4,/debug,/check,/save,/all)
;
; ======================================================================
;
;
FUNCTION grid_Execute, wall_path         , $
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
                       limiter = limiter

  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  debug  = option.debug
  xrange = option.xrange
  yrange = option.yrange

 
  aspect_ratio = .90 * (yrange[1] - yrange[0]) / (xrange[1] - xrange[0])
  ; aspect_ratio = 1.0
  !P.BACKGROUND = TrueColor('White')
  WINDOW,0,XSIZE=900/aspect_ratio,YSIZE=900

  cont = 1


  WHILE (cont) DO BEGIN

    cont = 0

    SWITCH step OF 
;     --------------------------------------------------------------------
;     Radial processing:
;     --------------------------------------------------------------------
      1: BEGIN

        wall = grid_LoadWall(wall_path+wall_file) ; ,/debug,xrange=xrange,yrange=yrange)
      
        b = grid_ReadEQUFile(equ_path+equ_file,CHANGE_SIGN=change_sign,psi_zone=psi_zone)
            ;shade_surf,b.psi,b.x,b.y; ,ax=0  

        b = grid_FindNullPoints(b,xpoint_zone,1,/debug)

        IF (machine EQ 'west_psi') THEN BEGIN
;          b.null_n++
;          b = grid_UpdateStructure(b,'null_i',[(b.null_i[0]+b.null_i[1])/2,b.null_i])
;          b = grid_UpdateStructure(b,'null_j',[N_ELEMENTS(b.y)/2,b.null_j])

;stop

;          wall.x = wall.x / 100.0
;          wall.y = wall.y / 100.0
print,b.null_i
        ENDIF 

;        print,xpoint_zone
;        print,equ_path+equ_file
;        help,b,/struct
        
        b = grid_RefineSeparatrices(b)  
;        b = grid_RefineSeparatrices(b, /debug, xrange=xrange, yrange=yrange)  
      
        IF (0 EQ 1) THEN BEGIN
          plot,wall.x,wall.y,color=Truecolor('Black')
          ; help,b,/struct
          OPLOT, [b.x[b.null_i[0]]], [b.y[b.null_j[0]]], PSYM=6, color=TrueColor('Red')      ; o-point
          OPLOT, [b.x[b.null_i[1]]], [b.y[b.null_j[1]]], PSYM=6, color=TrueColor('Green')   ; primary   x-point
          IF (N_ELEMENTS(b.null_i) EQ 3) THEN  $
            OPLOT, [b.x[b.null_i[2]]], [b.y[b.null_j[2]]], PSYM=6, color=TrueColor('Blue' ) ; secondary x-point
          OPLOT, [xpoint_zone[0],xpoint_zone[0],xpoint_zone[1],xpoint_zone[1],xpoint_zone[0]],  $
                 [xpoint_zone[2],xpoint_zone[3],xpoint_zone[3],xpoint_zone[2],xpoint_zone[2]],color=TrueColor('Red')
          stop
        ENDIF
      
;        print,'null_x',b.x[b.null_i]
;        print,'null_y',b.y[b.null_j]

        IF (KEYWORD_SET(preview)) THEN BEGIN

          first_plot = 1

          print,'option.ctr_boundary',option.ctr_boundary
      
          color = ['Black','Orange','Blue','Purple','Cyan','Red','Green']
          PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
;
;         --------------------------------------------------------------
;
          IF (b.null_n LE 2 OR -option.ctr_boundary[1] GT b.psi_2nd_xpoint) THEN single_null = 1 ELSE single_null = 0

          FOR i = 0, param.PFZ_SECONDARY DO BEGIN

            IF ((i EQ param.SOL_LFS OR i EQ param.SOL_HFS OR i EQ param.PFZ_SECONDARY) AND  $
                b.null_n LE 2) THEN CONTINUE

            IF (b.null_n GE 2) THEN BEGIN
              IF (i EQ param.PFZ_SECONDARY) THEN psi_start = b.psi_2nd_xpoint  $
                                            ELSE psi_start = b.psi_1st_xpoint
            ENDIF ELSE psi_start = 0.0

            IF (i EQ 0) THEN BEGIN
              IF (b.null_n GE 2) THEN ctr = grid_ExtractContour(b.psi, b.x, b.y, b.psi_1st_xpoint) $
                                 ELSE CONTINUE
            ENDIF ELSE  $
              ctr = grid_ExtractContour(b.psi, b.x, b.y, psi_start-option.ctr_boundary[i-1])

            ibrk = WHERE(ctr.dist GT option.ctr_break_distance) 
            nseg = N_ELEMENTS(ibrk) 
            ibrk = [-1,ibrk,ctr.n-1]

            FOR iseg = 0, nseg DO BEGIN
              print,'----------------------------------------'
              print,'contour seg:',iseg,ibrk[iseg]+1,ibrk[iseg+1]

              ; Not a real segment so skip it, may need to strengthen this check:
              IF (ibrk[iseg]+2 EQ ibrk[iseg+1] OR  $
                  ibrk[iseg]+1 EQ ibrk[iseg+1] OR  $
                  ibrk[iseg]   EQ ibrk[iseg+1]) THEN CONTINUE  
            
              x = ctr.x[ibrk[iseg]+1:ibrk[iseg+1]]
              y = ctr.y[ibrk[iseg]+1:ibrk[iseg+1]]

              ; Avoid double-null contours if a single-null grid is specified:
              IF ( ( (i EQ param.SOL_LFS OR i EQ param.SOL_HFS OR i EQ param.PFZ_SECONDARY) AND  $
                          single_null )                                                     OR   $ 
                   ( (i EQ param.SOL                                                      ) AND  $
                      NOT single_null AND b.null_n GE 2 ) ) THEN CONTINUE

              ; Check if the contour crosses the midplane:
              IF ( MIN(y) LT b.y[b.null_j[0]] AND MAX(y) GT b.y[b.null_j[0]] ) THEN cross = 1 ELSE cross = 0
              IF ( (i LE 4 AND NOT cross) OR  $
                   (i GE 5 AND     cross) ) THEN CONTINUE

              meanx = MEAN(x)
              IF ( (i EQ 3 AND meanx LT b.x[b.null_i[0]]) OR  $
                   (i EQ 4 AND meanx GT b.x[b.null_i[0]]) ) THEN CONTINUE       

              IF (b.null_n GE 2) THEN BEGIN
                meany = MEAN(y)
                IF ( (i EQ 5 AND meany GT b.y[b.null_j[0]] AND b.y[b.null_j[1]] LT b.y[b.null_j[0]] ) OR  $
                     (i EQ 5 AND meany LT b.y[b.null_j[0]] AND b.y[b.null_j[1]] GT b.y[b.null_j[0]] ) OR  $
                     (i EQ 6 AND meany LT b.y[b.null_j[0]] AND b.y[b.null_j[1]] LT b.y[b.null_j[0]] ) OR  $
                     (i EQ 6 AND meany GT b.y[b.null_j[0]] AND b.y[b.null_j[1]] GT b.y[b.null_j[0]] ) ) THEN CONTINUE       
              ENDIF 

              IF (first_plot) THEN BEGIN
                first_plot = 0
                PLOT, x, y, COLOR=Truecolor(color[i]), LINESTYLE=2,  $
                      XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1, /NOERASE
              ENDIF ELSE  $
                PLOT, x, y, COLOR=Truecolor(color[i]), PSYM=3,  $
                      XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1, /NOERASE

            ENDFOR

          ENDFOR

          PLOT,wall.x,wall.y,color=Truecolor('Black'),  $
               XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1, /NOERASE
          OPLOT, [b.x[b.null_i[0]]], [b.y[b.null_j[0]]], PSYM=6, color=TrueColor('Red')      ; o-point
          IF (b.null_n GE 2) THEN  $
            OPLOT, [b.x[b.null_i[1]]], [b.y[b.null_j[1]]], PSYM=6, color=TrueColor('Green')   ; primary   x-point
          IF (b.null_n GE 3) THEN  $
            OPLOT, [b.x[b.null_i[2]]], [b.y[b.null_j[2]]], PSYM=6, color=TrueColor('Blue' ) ; secondary x-point
      
          print,'psi at 1st separatrix',b.psi_1st_xpoint
          IF (b.null_n GE 3) THEN  $
            print,'       2nd separatrix',b.psi_2nd_xpoint
          print,'psi at core      ',b.psi[b.null_i[0],b.null_j[0]]
          print,'PSI_SOL=',option.ctr_boundary[1]
          XYOUTS,0.65,0.20,'CORE'          ,color=Truecolor(color[1]),/NORMAL,CHARSIZE=1.2
          XYOUTS,0.65,0.18,'SOL'           ,color=Truecolor(color[2]),/NORMAL,CHARSIZE=1.2 
          XYOUTS,0.65,0.16,'SOL, LFS'      ,color=Truecolor(color[3]),/NORMAL,CHARSIZE=1.2
          XYOUTS,0.65,0.14,'SOL, HFS'      ,color=Truecolor(color[4]),/NORMAL,CHARSIZE=1.2
          XYOUTS,0.65,0.12,'PFZ         '  ,color=Truecolor(color[5]),/NORMAL,CHARSIZE=1.2
          XYOUTS,0.65,0.10,'PFZ, SECONDARY',color=Truecolor(color[6]),/NORMAL,CHARSIZE=1.2

          RETURN, -1
        ENDIF 
      
        IF (debug) THEN PRINT,'Calling grid_AnalyseBoundary',mswin

        b = grid_AnalyseBoundary(b,wall,user_step,user_finish, machine,  $
                                 save=save,mswin=mswin,limiter=limiter)

        IF (NOT KEYWORD_SET(all)) THEN BREAK
        END
;     --------------------------------------------------------------------
;     Distribute radial contours and identify poloidal domains:
;     --------------------------------------------------------------------
      2: BEGIN
      
        RESTORE,'grid_data/stored_step1_'+machine+'.sav'
        IF (KEYWORD_SET(hold_state)) THEN state = hold_state

        contour_array.contour1.separatrix = 1 ; why is this here, and not somewhere else?

        IF (state.version NE 1.0) THEN BEGIN
          PRINT,'ERROR grid_Execute: STATE variable version not correct'
          STOP
        ENDIF
      
        grid_ZoneWall, wall.pt1, wall.pt2, debug=debug, xrange=xrange, yrange=yrange
      
        PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
      
        ; Check if the secondary separatrix needs to be added:
        IF (scan_params.process_2nd NE -1) THEN BEGIN

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
            psi_shift = -0.0000001D ; -0.00000001D

          IF (status) THEN  $
            status = grid_AddContour(b, wall, scan_params, contour_array, -1,  $
                                     psi_val=b.psi_2nd_xpoint + psi_shift,  $
                                     debug=debug, xrange=xrange, yrange=yrange)

        ENDIF
      
        IF (state.geometry NE param.LIMITED) THEN  $
          contour_array = grid_InstallXPoints(b,scan_params,contour_array,  $
                                              debug=debug, xrange=xrange, yrange=yrange)

        grid_UpdateWall, b, contour_array, wall,  $
                         debug=debug, xrange=xrange, yrange=yrange
;        grid_Debug,b,contour_array,wall,xrange,yrange

        grid_ProcessWall, b, wall, scan_params, contour_array, $
                          debug=debug, xrange=xrange, yrange=yrange

        grid_UpdateWall, b, contour_array, wall, $
                         debug=debug, xrange=xrange, yrange=yrange

        contour_array = grid_BuildRadialMap(contour_array, wall,  $
                                            debug=debug, xrange=xrange, yrange=yrange)
      
        contour_array = grid_RadialRefinement(contour_array, wall, b,  $
                                              debug=debug, xrange=xrange, yrange=yrange)

        contour_array = grid_SetPoloidalDomains(contour_array, wall, b, $
                                                debug=debug, xrange=xrange, yrange=yrange)
      
        IF (KEYWORD_SET(save)) THEN  $
          SAVE,filename='grid_data/stored_step2_'+machine+'.sav',b,wall,contour_array,scan_params,state
      
        print,'taking a rest now'
        ;grid_Debug,b,contour_array,wall,xrange,yrange

        IF (NOT KEYWORD_SET(all)) THEN BREAK      
        END
;     --------------------------------------------------------------------
;     Set poloidal descretisation (cell building):
;     --------------------------------------------------------------------
      3: BEGIN
      
        RESTORE,'grid_data/stored_step2_'+machine+'.sav'
        IF (KEYWORD_SET(hold_state)) THEN state = hold_state

        IF (KEYWORD_SET(debug)) THEN BEGIN
          tags = STRUPCASE(TAG_NAMES(contour_array))
          nctr = N_ELEMENTS(tags)
          PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
          OPLOT,wall.x,wall.y,color=Truecolor('Black')
          FOR i = 0, nctr-1 DO BEGIN
            ctr = grid_ExtractStructure(contour_array,tags[i])      
            OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
            npts = N_ELEMENTS(ctr.x)
            XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('Darkgrey')
          ENDFOR
        ENDIF
      
        contour_array = grid_CreatePoloidalPoints(contour_array, wall, b,  $
                                                debug=debug, xrange=xrange, yrange=yrange)
      
        IF (state.geometry NE param.LIMITED) THEN  $
          contour_array = grid_EnforceMidplaneAlignment(contour_array, wall, b)
      
        IF (KEYWORD_SET(save)) THEN  $
          SAVE,filename='grid_data/stored_step3_'+machine+'.sav',b,wall,contour_array,scan_params,state

        IF (debug) THEN PRINT,'STEP 3 COMPLETE'

        IF (NOT KEYWORD_SET(all)) THEN BREAK      
        END
;     --------------------------------------------------------------------
;     Add core rings/cells and generate grid:
;     --------------------------------------------------------------------
      4: BEGIN
      
        RESTORE,'grid_data/stored_step3_'+machine+'.sav'
        IF (KEYWORD_SET(hold_state)) THEN state = hold_state
      
        IF (KEYWORD_SET(debug)) THEN BEGIN
          tags = STRUPCASE(TAG_NAMES(contour_array))
          nctr = N_ELEMENTS(tags)
          PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Darkgray')
          OPLOT,wall.x,wall.y,color=Truecolor('Black')
          FOR i = 0, nctr-1 DO BEGIN
            ctr = grid_ExtractStructure(contour_array,tags[i])      
            OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
            OPLOT,ctr.grid_x,ctr.grid_y,color=Truecolor('Pink'),PSYM=3
            npts = N_ELEMENTS(ctr.x)
            XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('Darkgrey')
          ENDFOR
        ENDIF
      
        contour_array = grid_AddCoreRings(  $
          contour_array, wall, b, debug=debug, xrange=xrange, yrange=yrange)
      
        tubes = grid_BuildGrid(  $
          contour_array, wall, b, debug=debug, xrange=xrange, yrange=yrange)
      
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
      
        IF (KEYWORD_SET(save)) THEN  $
          SAVE,filename='grid_data/stored_step4_'+machine+'.sav',b,wall,contour_array,scan_params,tubes,state

        IF (debug) THEN PRINT,'STEP 4 COMPLETE'
      
        IF (NOT KEYWORD_SET(all)) THEN BREAK      
        END
;     --------------------------------------------------------------------
;     Check the grid for problems and implement corrections or generate
;     the grid file:
;     --------------------------------------------------------------------
      5: BEGIN
      
        RESTORE,'grid_data/stored_step4_'+machine+'.sav'
        IF (KEYWORD_SET(hold_state)) THEN state = hold_state

        IF (KEYWORD_SET(debug)) THEN BEGIN
          tags = STRUPCASE(TAG_NAMES(contour_array))
          nctr = N_ELEMENTS(tags)
          PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Darkgray')
          OPLOT,wall.x,wall.y,color=Truecolor('Black')
          FOR i = 0, nctr-1 DO BEGIN
            ctr = grid_ExtractStructure(contour_array,tags[i])      
            OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
            OPLOT,ctr.grid_x,ctr.grid_y,color=Truecolor('Pink'),PSYM=3
            npts = N_ELEMENTS(ctr.x)
            XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('Darkgrey')
          ENDFOR
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
        ENDIF
 
        desperation:

        IF (KEYWORD_SET(check)) THEN BEGIN
;      
;         ----------------------------------------------------------------
;         CHECK IF THERE IS A PROBLEM WITH THE GRID
;      
          ; Store current parameters related to automated grid refinement:
          state_refinement_rad = state.refinement_rad
          state_refinement_pol = state.refinement_pol
          IF (state.refinement_rad EQ 0) THEN  $
            state_ref_rad_n = 0 ELSE  $
            state_ref_rad_n = N_ELEMENTS(state.ref_rad_psi1)
          IF (state.refinement_pol EQ 0) THEN  $
            state_ref_pol_n = 0 ELSE  $
            state_ref_pol_n = N_ELEMENTS(state.ref_pol_pos )

          result = grid_CheckIntegrity(contour_array, wall, b, tubes)

          IF ((state.refinement_rad GT state_refinement_rad) OR  $
              (state.refinement_pol GT state_refinement_pol)) THEN BEGIN

            cont = 1
            step = 0
       
            print,'n1',state_ref_rad_n,N_ELEMENTS(state.ref_rad_psi1)
            print,'n2',state_ref_pol_n,N_ELEMENTS(state.ref_pol_pos)

            print,'n3',state.refinement_rad,state_refinement_rad
            print,'n4',state.refinement_pol,state_refinement_pol

            IF (state.refinement_rad GT state_refinement_rad) THEN BEGIN
              step = 2 
              IF (state_ref_rad_n EQ N_ELEMENTS(state.ref_rad_psi1)) THEN cont = 0
            ENDIF
       
            IF (state.refinement_pol GT state_refinement_pol) THEN BEGIN
              IF (step EQ 0) THEN step = 3
              IF (state_ref_pol_n EQ N_ELEMENTS(state.ref_pol_pos)) THEN cont = 0
            ENDIF
             
            IF (cont EQ 0) THEN BEGIN
              PRINT,'ERROR grid_Execute: Integrity check failed but no enhanced refinement specified'
              STOP
            ENDIF
       
            hold_state = state

          ENDIF
        ENDIF 

        IF (KEYWORD_SET(check) AND NOT cont) THEN BEGIN
;
;         --------------------------------------------------------------
;         ADJUST GRID CELLS
;
          count = 0
          again = 1
          WHILE (again) DO BEGIN
            count = count + 1
            IF (count EQ 4) THEN BEGIN
              print, 'whoa! 3 passes not enough...'
;              GOTO, desperation
              stop
            ENDIF
            result = grid_MassageCells(contour_array, wall, b, tubes, again)
          ENDWHILE

;          IF (count GT 1) THEN GOTO, desperation

        ENDIF

        IF (NOT cont) THEN BEGIN
;
;         --------------------------------------------------------------
;         WRITE THE GRID FILE
;
          cont = 0
          fname = 'grid_' + machine 
          result = grid_WriteGridFile(fname, tubes, contour_array, wall, b,  $ 
                                      debug=debug, xrange=xrange, yrange=yrange)

          fname = 'grid_output'
          result = grid_WriteGridFile(fname, tubes, contour_array, wall, b,  $ 
                                      debug=debug, xrange=xrange, yrange=yrange)

          PRINT,'----------------------------'
          PRINT,'Grid generated successufully'
          PRINT,'----------------------------'
        ENDIF

        BREAK
        END
;     --------------------------------------------------------------------
;     --------------------------------------------------------------------
      ELSE: BEGIN
        PRINT,'ERROR grid_Execute: Unknown step'
        PRINT,' STEP =',step
        RETURN,-1
        END
;     --------------------------------------------------------------------
;     --------------------------------------------------------------------
    ENDSWITCH

  ENDWHILE

  result = 1

  RETURN, result
END
;
; ======================================================================
;

