;
; ======================================================================
;
FUNCTION grid_MassageCells, contours, wall, b, tubes, again
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  debug  = option.debug
  xrange = option.xrange
  yrange = option.yrange

  again = 0

  ctags = STRUPCASE(TAG_NAMES(contours))
  ttags = STRUPCASE(TAG_NAMES(tubes))

  nctr  = N_ELEMENTS(ctags)
  ntube = N_ELEMENTS(ttags)

  print,'--------------------------------------'
  print,'Let the massage begin...'
  print,' '

  FOR itube = 1, ntube DO BEGIN

    tube = grid_ExtractStructure(tubes,ttags[itube-1])      

    ictr1 = tube.ictr1
    ictr2 = tube.ictr2

    ctr1 = grid_ExtractStructure(contours,ctags[ictr1-1])
    ctr2 = grid_ExtractStructure(contours,ctags[ictr2-1])

    ; Assume that the core is well behaved:
    IF (ctr1.region EQ param.CORE OR ctr2.region EQ param.CORE) THEN CONTINUE

;    print,'--------------------------------------',itube,ntube
;    print,'itube',itube,FORMAT='(A,I10)'
;    print,'ictr1',ictr1,ctr1.map_in,ctr1.map_out,FORMAT='(A,5I10)'
;    print,'ictr2',ictr2,ctr2.map_in,ctr2.map_out,FORMAT='(A,5I10)'
;
;   -------------------------------------------------------------------
;   SCAN OVER THE CELLS IN EACH TUBE
;
    FOR ipass = 1, 2 DO BEGIN

       IF (ipass EQ 1) THEN BEGIN
         istart =  tube.n / 2
         iend   =  0
         istep  = -1
       ENDIF ELSE BEGIN
         istart =  tube.n / 2
         iend   =  tube.n - 1
         istep  = +1
       ENDELSE

      FOR icell = istart, iend, istep DO BEGIN ; 0, tube.n-1 DO BEGIN

       log_message = 1
;      FOR icell = tube.n-1, 0, -1 DO BEGIN ; 0, tube.n-1 DO BEGIN
      
        cell_x = [tube.x1[icell],tube.x2[icell],tube.x2[icell+1],tube.x1[icell+1],tube.x1[icell]]
        cell_y = [tube.y1[icell],tube.y2[icell],tube.y2[icell+1],tube.y1[icell+1],tube.y1[icell]]
      
      
        oplot,cell_x,cell_y,COLOR=Truecolor('Darkgray')
;     
;       -----------------------------------------------------------------
;       CHECK IF ANY INTERIOR CELL ANGLES ARE LARGER THAN 180 DEGREES
;     
        FOR j2 = 0, 3 DO BEGIN
      
          j1 = j2 - 1
          j3 = j2 + 1
          j4 = j2 + 2
          IF (j1 EQ -1) THEN j1 = 3
          IF (j3 EQ  4) THEN j3 = 0
          IF (j4 GE  4) THEN j4 = j4 - 4
      
          hold_x = -999.0D
      
          cont = 1
          WHILE (cont) DO BEGIN
            cont = 0
      
            v1 = [cell_x[j1],cell_y[j1]]
            v2 = [cell_x[j2],cell_y[j2]]
            v3 = [cell_x[j3],cell_y[j3]]
            v4 = [cell_x[j4],cell_y[j4]]
            
            as = SQRT( (v1[0]-v2[0])^2 + (v1[1]-v2[1])^2 )
            bs = SQRT( (v4[0]-v2[0])^2 + (v4[1]-v2[1])^2 )
            cs = SQRT( (v4[0]-v1[0])^2 + (v4[1]-v1[1])^2 )
            angle1 = ACOS( (as^2 + bs^2 - cs^2) / (2.0D * as * bs) ) * 180.0D / 3.1415D

;            as1 = as
;            bs1 = bs
;            cs1 = cs
            
            as = SQRT( (v3[0]-v2[0])^2 + (v3[1]-v2[1])^2 )
            bs = SQRT( (v4[0]-v2[0])^2 + (v4[1]-v2[1])^2 )
            cs = SQRT( (v4[0]-v3[0])^2 + (v4[1]-v3[1])^2 )
            angle2 = ACOS( (as^2 + bs^2 - cs^2) / (2.0D * as * bs) ) * 180.0D / 3.1415D
            
;            as2 = as
;            bs2 = bs
;            cs2 = cs

            angle = angle1 + angle2
      
; if (itube eq 58) then begin
;   print,'checking',icell,j2,angle1,angle2,angle
; endif
      
            IF (angle GT 180.0D) THEN BEGIN

              IF (log_message) THEN BEGIN              
                print,'********** LOGGED BIG ANGLE! **********'
                print,icell,j2,angle1,angle2,angle	          
;                print,'v',v1[0],v2[0],v3[0],v4[0]
;                print,' ',v1[1],v2[1],v3[1],v4[1]
;                print,'a,b,cs1',as1,bs1,cs1
;                print,'a,b,cs2',as2,bs2,cs2
                log_message = 0
              ENDIF

              cont  = 1
              again = 1
              
              deltax = MAX(cell_x) - MIN(cell_x)
              deltay = MAX(cell_y) - MIN(cell_y)
              cen_x = MEAN(cell_x[0:N_ELEMENTS(cell_x)-2])
              cen_y = MEAN(cell_y[0:N_ELEMENTS(cell_y)-2])
;              xrange = [cen_x-1.0D*deltax,cen_x+1.0D*deltax]
;              yrange = [cen_y-1.0D*deltay,cen_y+1.0D*deltay]
              PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
              OPLOT, wall.x, wall.y, COLOR=Truecolor('Black')
              OPLOT,cell_x,cell_y,COLOR=Truecolor('Orange')
              OPLOT,cell_x,cell_y,COLOR=Truecolor('Orange'),PSYM=7

              XYOUTS, [v1[0]], [v1[1]], CHARSIZE=1.1, ALIGNMENT=0.0, 'j1'
              XYOUTS, [v2[0]], [v2[1]], CHARSIZE=1.1, ALIGNMENT=0.0, 'j2'
              XYOUTS, [v3[0]], [v3[1]], CHARSIZE=1.1, ALIGNMENT=0.0, 'j3'
              XYOUTS, [v4[0]], [v4[1]], CHARSIZE=1.1, ALIGNMENT=0.0, 'j4'

              XYOUTS, [cell_x[0]],[cell_y[0]], CHARSIZE=1.1, ALIGNMENT=1.0, 'v1'
              XYOUTS, [cell_x[1]],[cell_y[1]], CHARSIZE=1.1, ALIGNMENT=1.0, 'v2'
              XYOUTS, [cell_x[2]],[cell_y[2]], CHARSIZE=1.1, ALIGNMENT=1.0, 'v3'
              XYOUTS, [cell_x[3]],[cell_y[3]], CHARSIZE=1.1, ALIGNMENT=1.0, 'v4'

              CASE j2 OF
                0: BEGIN
                  iside1 = 1
                  iside2 = 2
                  END
                1: BEGIN
                  iside1 = 0
                  iside2 = 3
                  END
                2: BEGIN
                  iside1 = 3
                  iside2 = 0
;                  stop
                  END
                3: BEGIN  
;                  IF (tube.region EQ param.PFZ OR tube.region EQ param.PFZ_SECONDARY) THEN BEGIN
;                    iside1 = 3 ; 2 ; 3
;                    iside2 = 0 ; 1 ; 0
;                  ENDIF ELSE BEGIN
                    iside1 = 2 ; 3
                    iside2 = 1 ; 0
;                  ENDELSE
                  END
              ENDCASE
      
;              print,j1,j2,j3,j4,iside1
            
              IF (hold_x EQ -999.0D) THEN BEGIN
                hold_x = cell_x[iside1]
                hold_y = cell_y[iside1]
              ENDIF
              
              frac = 0.05D
              new_x = cell_x[iside1] + frac * (cell_x[iside2] - cell_x[iside1])
              new_y = cell_y[iside1] + frac * (cell_y[iside2] - cell_y[iside1])
      
              cell_x[iside1] = new_x 
              cell_y[iside1] = new_y 
      
            ENDIF
          ENDWHILE
      
          IF (hold_x NE -999.0D) THEN BEGIN
      
            print,icell,j2,angle1,angle2,angle
            OPLOT,cell_x,cell_y,COLOR=Truecolor('Lightgreen'),PSYM=7
;     
;           --------------------------------------------------------------
;           UPDATE CELL POINT ON CURRENT FOCUS TUBE
;     
            i1 = WHERE(tube.x1 EQ hold_x AND tube.y1 EQ hold_y,count1)
            i2 = WHERE(tube.x2 EQ hold_x AND tube.y2 EQ hold_y,count2)
            IF (count1 GT 0) THEN BEGIN
              tube.x1[i1] = new_x
              tube.y1[i1] = new_y
            ENDIF
            IF (count2 GT 0) THEN BEGIN
              tube.x2[i2] = new_x
              tube.y2[i2] = new_y
            ENDIF
            tubes = grid_UpdateStructure(tubes,ttags[itube-1],tube)
;     
;           --------------------------------------------------------------
;           FIND CORRESPONDING POINT ON ADJACENT TUBE (NO CONNECTION MAP)
;     
            FOR itube2 = 1, ntube DO BEGIN
              IF (itube EQ itube2) THEN CONTINUE
      
              tube2 = grid_ExtractStructure(tubes,ttags[itube2-1])      
              i1 = WHERE(tube2.x1 EQ hold_x AND tube2.y1 EQ hold_y,count1)
              i2 = WHERE(tube2.x2 EQ hold_x AND tube2.y2 EQ hold_y,count2)
      
              IF (count1 NE 0) THEN BEGIN
                print,'found 1!',itube,itube2
                tube2.x1[i1] = new_x
                tube2.y1[i1] = new_y
                tubes = grid_UpdateStructure(tubes,ttags[itube2-1],tube2)
                BREAK
              ENDIF
              IF (count2 NE 0) THEN BEGIN
                 print,'found 2!',itube,itube2
                tube2.x2[i2] = new_x
                tube2.y2[i2] = new_y
                tubes = grid_UpdateStructure(tubes,ttags[itube2-1],tube2)
                BREAK
              ENDIF 
            ENDFOR          
      
;            FOR ictr = 1, nctr DO BEGIN
;              ctr = grid_ExtractStructure(contours,ctags[ictr-1])
;              i = WHERE(ctr.grid_x EQ hold_x AND ctr.grid_y EQ hold_y, count)              
;              IF (count NE 0) THEN BEGIN
;                print,'found 3!',ictr
;                ctr.grid_x[i] = new_x
;                ctr.grid_y[i] = new_y
;                ctr = grid_UpdateStructure(contours,ctags[ictr-1],ctr)
;              ENDIF
;            ENDFOR

            IF (itube2 EQ ntube+1) THEN BEGIN
              print,'*** could be a problem *** ' ; but maybe not if the tube is 'hanging' somehow
              print,'ctr1',ctr1.psi
              print,' map in ',ctr1.map_in
              print,' map out',ctr1.map_out
              print,'ctr2',ctr2.psi
              print,' map in ',ctr2.map_in
              print,' map out',ctr2.map_out
;              stop
            ENDIF
      
          ENDIF
      
        ENDFOR  ; j2
      ENDFOR  ; icell
    ENDFOR  ; ipasgs
  ENDFOR  ; itube


;stop

END

;
; ======================================================================
;
FUNCTION grid_CheckIntegrity, contours, wall, b, tubes,  $
                              problem_radial   = problem__radial,  $
                              problem_poloidal = problem__poloidal
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  debug  = option.debug
  xrange = option.xrange
  yrange = option.yrange

  ctags = STRUPCASE(TAG_NAMES(contours))
  ttags = STRUPCASE(TAG_NAMES(tubes))

  nctr  = N_ELEMENTS(ctags)
  ntube = N_ELEMENTS(ttags)

  IF (state.refinement_rad EQ 0) THEN BEGIN
    ref_rad_x1   = [-999.0]
    ref_rad_x2   = [-999.0]
    ref_rad_y1   = [-999.0]
    ref_rad_y2   = [-999.0]
    ref_rad_psi1 = [-999.0]
    ref_rad_psi2 = [-999.0]
  ENDIF ELSE BEGIN
    ref_rad_x1   = state.ref_rad_x1  
    ref_rad_y1   = state.ref_rad_y1  
    ref_rad_x2   = state.ref_rad_x2  
    ref_rad_y2   = state.ref_rad_y2  
    ref_rad_psi1 = state.ref_rad_psi1
    ref_rad_psi2 = state.ref_rad_psi2
  ENDELSE

  IF (state.refinement_pol EQ 0) THEN BEGIN
    ref_pol_i1  = [-999  ]
    ref_pol_i2  = [-999  ]
    ref_pol_pos = [-999.0]
  ENDIF ELSE BEGIN
    ref_pol_i1  = state.ref_pol_i1
    ref_pol_i2  = state.ref_pol_i2
    ref_pol_pos = state.ref_pol_pos
  ENDELSE

  bubble_count = 0
;
; ----------------------------------------------------------------------
; LOOP OVER ALL CELL AND SEE IF VERTICES FROM NEIGHBOURING TUBES
; ARE INSIDE THE CURRENT FOCUS CELL, i.e. IF THE GRID IS MALFORMED
;
  problem__radial   = 0
  problem__poloidal = 0

  FOR itube = 1, ntube DO BEGIN

    tube = grid_ExtractStructure(tubes,ttags[itube-1])      

    ictr1 = tube.ictr1
    ictr2 = tube.ictr2

    ctr1 = grid_ExtractStructure(contours,ctags[ictr1-1])
    ctr2 = grid_ExtractStructure(contours,ctags[ictr2-1])

    ; Assume that the core is well behaved:
    IF (ctr1.region EQ param.CORE OR ctr2.region EQ param.CORE) THEN CONTINUE

;oplot,ctr1.grid_x,ctr1.grid_y,COLOR=Truecolor('Lightgreen'), PSYM=6
;oplot,ctr2.grid_x,ctr2.grid_y,COLOR=Truecolor('Lightgreen'), PSYM=6

    print,'--------------------------------------',itube,ntube
    print,'itube',itube,FORMAT='(A,I10)'
    print,'ictr1',ictr1,ctr1.map_in,ctr1.map_out,FORMAT='(A,5I10)'
    print,'ictr2',ictr2,ctr2.map_in,ctr2.map_out,FORMAT='(A,5I10)'

    map_i = [REFORM(ctr1.map_in),REFORM(ctr1.map_out),  $
             REFORM(ctr2.map_in),REFORM(ctr2.map_out)]
;    j = WHERE(map_i NE -1 AND map_i NE ictr1 AND map_i NE ictr2, count)
    j = WHERE(map_i NE -1, count)

    IF (count EQ 0) THEN BEGIN
      print, 'unexpected this is'
      stop
    ENDIF

    map_i = map_i[j]

;    print, 'map_i',map_i
    
    FOR map_j = 0, N_ELEMENTS(map_i)-1 DO BEGIN
 
      ictr = (map_i[map_j])[0]
;      print,'going map_i',map_i
;      print,'ictr   ',ictr
      ctr  = grid_ExtractStructure(contours,ctags[ictr-1])      
;
;     -------------------------------------------------------------------
;     SCAN OVER THE CELLS IN EACH TUBE
;
      FOR icell = 0, tube.n-1 DO BEGIN

        cell_x = [tube.x1[icell],tube.x2[icell],tube.x2[icell+1],tube.x1[icell+1],tube.x1[icell]]
        cell_y = [tube.y1[icell],tube.y2[icell],tube.y2[icell+1],tube.y1[icell+1],tube.y1[icell]]

        cen_x = MEAN(cell_x[0:N_ELEMENTS(cell_x)-2])
        cen_y = MEAN(cell_y[0:N_ELEMENTS(cell_y)-2])

        oplot,cell_x,cell_y,COLOR=Truecolor('Pink')
;
;       -----------------------------------------------------------------
;       CHECK IF THE CELL ITSELF IS MALFORMED: OPPOSITE SIDES CROSS
;
        IF (map_j EQ 0) THEN BEGIN
          result = grid_Intersection([cell_x[0],cell_y[0]],[cell_x[1],cell_y[1]],  $
                                     [cell_x[2],cell_y[2]],[cell_x[3],cell_y[3]], 0, status=status)
          IF (status) THEN BEGIN

            bubble_count++

            deltax = MAX(cell_x) - MIN(cell_x)
            deltay = MAX(cell_y) - MIN(cell_y)
            xrange = [cen_x-0.6D*deltax,cen_x+0.6D*deltax]
            yrange = [cen_y-0.6D*deltay,cen_y+0.6D*deltay]

            PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
            OPLOT, wall.x, wall.y, COLOR=Truecolor('Black')
            FOR itube2 = 1, ntube DO BEGIN
              tube2 = grid_ExtractStructure(tubes,ttags[itube2-1])
              FOR icell2 = 0, tube2.n-1 DO BEGIN
                OPLOT,[tube2.x1[icell2],tube2.x2[icell2],tube2.x2[icell2+1],tube2.x1[icell2+1],tube2.x1[icell2]],  $
                      [tube2.y1[icell2],tube2.y2[icell2],tube2.y2[icell2+1],tube2.y1[icell2+1],tube2.y1[icell2]],  $
                      COLOR=Truecolor('Lightblue')
              ENDFOR
            ENDFOR
            OPLOT,ctr.x,ctr.y,COLOR=Truecolor('Red')
            OPLOT,cell_x,cell_y

            print, '=================================================>>>> bubble',bubble_count,icell,map_j

            print,'********** ENDS CROSSING! **********',i,icell,map_j,itube

            j = WHERE(ref_rad_psi1 EQ MIN([ctr1.psi,ctr2.psi]) AND  $
                      ref_rad_psi2 EQ MAX([ctr1.psi,ctr2.psi]), count)
            IF (count EQ 0) THEN BEGIN

              print,'********** LOGGED! **********'
    
              problem__radial = 1
    
              ref_rad_x1   = [ref_rad_x1  ,tube.x1[0]]
              ref_rad_y1   = [ref_rad_y1  ,tube.y1[0]]
              ref_rad_x2   = [ref_rad_x2  ,tube.x2[0]]
              ref_rad_y2   = [ref_rad_y2  ,tube.y2[0]]
              ref_rad_psi1 = [ref_rad_psi1,MIN([ctr1.psi,ctr2.psi])]
              ref_rad_psi2 = [ref_rad_psi2,MAX([ctr1.psi,ctr2.psi])]

            ENDIF

            CONTINUE

;            if (bubble_count EQ 2) THEN STOP

          ENDIF

          result = grid_Intersection([cell_x[0],cell_y[0]],[cell_x[1],cell_y[1]],  $
                                     [cell_x[2],cell_y[2]],[cell_x[3],cell_y[3]], 0, status=status)
          IF (status) THEN BEGIN
            print, 'no shit!'
            stop
          ENDIF

        ENDIF
;
;       -----------------------------------------------------------------
;       IDENTIFY VERTICES THAT ARE NEAR THE CURRENT FOCUS CELL
;
        x = ctr.grid_x
        y = ctr.grid_y
;        oplot,x,y,COLOR=Truecolor('Black')

        proximity = SQRT( (x - cen_x)^2 + (y- cen_y)^2)
        range     = MAX([SQRT( (cell_x[2] - cell_x[0])^2 + (cell_y[2] - cell_y[0])^2) ,  $
                         SQRT( (cell_x[3] - cell_x[1])^2 + (cell_y[3] - cell_y[1])^2) ]) 
        i = WHERE(proximity LT range, count)  
        IF (count EQ 0) THEN CONTINUE
        x = x[i]
        y = y[i]
;
;       -----------------------------------------------------------------
;       CHECK EACH VERTEX TO SEE IF IT'S INSIDE THE FOCUS CELL
;
        FOR i = 0, N_ELEMENTS(x)-1 DO BEGIN

          inside = grid_PointInPolygon((x[i])[0],(y[i])[0],cell_x,cell_y)

          IF (inside EQ 1) THEN BEGIN
;
;           -------------------------------------------------------------
;           YES, WE HAVE A PROBLEM, STORE LOCATION 
;
;            xspan = 0.1D
;            yspan = 5.0D
;            xspan = 0.02D
;            yspan = 1.0D
;            xspan = 0.725D
;            yspan = 0.725D
;            xrange = [x[i]-range*xspan,x[i]+range*xspan] 
;            yrange = [y[i]-range*yspan,y[i]+range*yspan]

            deltax = MAX(cell_x) - MIN(cell_x)
            deltay = MAX(cell_y) - MIN(cell_y)
            scale = 15.0D
            xrange = [cen_x-scale*deltax,cen_x+scale*deltax]
            yrange = [cen_y-scale*deltay,cen_y+scale*deltay]

            xrange = option.xrange
            yrange = option.yrange

            delta = 0.5D / 2.0D0
            xrange = [cen_x-delta,cen_x+delta]
            yrange = [cen_y-delta,cen_y+delta]


            PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
            OPLOT, wall.x, wall.y, COLOR=Truecolor('Black')
            FOR itube2 = 1, ntube DO BEGIN
              tube2 = grid_ExtractStructure(tubes,ttags[itube2-1])
              FOR icell2 = 0, tube2.n-1 DO BEGIN
                OPLOT,[tube2.x1[icell2],tube2.x2[icell2],tube2.x2[icell2+1],tube2.x1[icell2+1],tube2.x1[icell2]],  $
                      [tube2.y1[icell2],tube2.y2[icell2],tube2.y2[icell2+1],tube2.y1[icell2+1],tube2.y1[icell2]],  $
                      COLOR=Truecolor('Lightblue')
              ENDFOR
            ENDFOR
            OPLOT,ctr.x,ctr.y,COLOR=Truecolor('Red')
            OPLOT,cell_x,cell_y
            OPLOT,x,y,PSYM=6
            OPLOT,[x[i]],[y[i]],PSYM=6,COLOR=Truecolor('Red')
            XYOUTS, [x[i]], [y[i]], CHARSIZE=1.1, ALIGNMENT=0.5, $
                    STRTRIM(STRING(ictr ),2) +  $
                    STRTRIM(STRING(ictr1),2) +  $
                    STRTRIM(STRING(ictr2),2) ,  $
                    COLOR=TrueColor('Darkgray')
;print,'map_i',map_i
;print,'ictr1,2,none',ictr1,ictr2,ictr

            IF (ictr EQ ictr1 OR ictr EQ ictr2) THEN BEGIN
;
;             ----------------------------------------------------------
;             A PROBLEM IS FOUND ON THE TUBE ITSELF, SO NEED HIGHER
;             RADIAL RESOLUTION
;
              print,'********** RADIAL PROBLEM! **********',i,icell,map_j,itube

              j = WHERE(ref_rad_psi1 EQ MIN([ctr1.psi,ctr2.psi]) AND  $
                        ref_rad_psi2 EQ MAX([ctr1.psi,ctr2.psi]), count)
              IF (count GT 0) THEN CONTINUE

              print,'********** LOGGED! **********'

              problem__radial = 1

              ref_rad_x1   = [ref_rad_x1  ,tube.x1[0]]
              ref_rad_y1   = [ref_rad_y1  ,tube.y1[0]]
              ref_rad_x2   = [ref_rad_x2  ,tube.x2[0]]
              ref_rad_y2   = [ref_rad_y2  ,tube.y2[0]]
              ref_rad_psi1 = [ref_rad_psi1,MIN([ctr1.psi,ctr2.psi])]
              ref_rad_psi2 = [ref_rad_psi2,MAX([ctr1.psi,ctr2.psi])]

            ENDIF ELSE BEGIN
;
;             ----------------------------------------------------------
;             A PROBLEM IS FOUND WITH NEIGHBOURING CONTOURS, SO NEED
;             BETTER POLOIDAL RESOLUTION
;
              print,'********** POLOIDAL PROBLEM! **********',i,icell,map_j,itube

              problem__poloidal = 1

              j = WHERE(ctr.grid_x EQ x[i] AND ctr.grid_y EQ y[i], count)
              IF (count NE 1) THEN BEGIN
                PRINT,'ERROR grid_CheckIntegrity: Grid point not found on contour'
                STOP
              ENDIF
              
              j = (j[0])[0]
              FOR j1 = j-1,                        0, -1 DO IF (ctr.grid_s[j1] NE 0) THEN BREAK  
              FOR j2 = j  , N_ELEMENTS(ctr.grid_s)-1     DO IF (ctr.grid_s[j2] NE 0) THEN BREAK  
              IF (j1 EQ -1 OR j2 EQ N_ELEMENTS(ctr.grid_s)) THEN BEGIN
                PRINT,'ERROR grid_CheckIntegrity: Marker points not identified'
                PRINT,' ICTR   = ',ictr
                PRINT,' GRID_S = ',ctr.grid_s 
                PRINT,' J,N    = ',j,N_ELEMENTS(ctr.grid_s)
                PRINT,' J1,J2  = ',j1,j2
                STOP
              ENDIF            
	      
              print,'ictr',ictr
              print,'j1,j,j2',j1,j,j2
              print,'grid_s',ctr.grid_s[j1],ctr.grid_s[j2]
	      
              print,'section1',state.section_i1
              print,'section2',state.section_i2
	      
              ; Just a quick check to see if things will be easy, i.e. not skipping a marker
              ; point somehow:
              span_xpoint = 0
              FOR k = 0, N_ELEMENTS(state.section_i1)-1 DO   $ 
                IF (state.section_i1[k] EQ ctr.grid_s[j1] AND  $
                    state.section_i2[k] EQ ctr.grid_s[j2]) THEN BREAK
              IF (k EQ N_ELEMENTS(state.section_i1)) THEN BEGIN
                ; Check if the issue is the primary x-point: 
                span_separatrix = 0
                ctr_s = grid_ExtractStructure(contours,'contour1') ; separatrix contour
                l1 = WHERE( ctr_s.grid_s EQ ctr.grid_s[j1] , count1)
                l2 = WHERE( ctr_s.grid_s EQ ctr.grid_s[j2] , count2)
                IF (count1 EQ 0 OR count2 EQ 0) THEN BEGIN
                  PRINT,'ERROR grid_CheckIntegrity: This is odd'
                  STOP
                ENDIF
                PRINT,'l1,2=',l1,l2
                l = WHERE( ctr_s.grid_s[l1+1:l2-1] NE 0 , count)
                IF (count NE 1) THEN BEGIN
                  PRINT,'ERROR grid_CheckIntegrity: What does "baffled" mean?'
                  STOP
                ENDIF                
                PRINT,'l=',l,ctr_s.grid_s[l]
                PRINT,'  ',ctr_s.grid_s[l1:l2]

                ; Super lame check to see if this is the separatix point, i.e. see if there
                ; are two of them:
                dummy = WHERE (ctr_s.grid_x EQ (ctr_s.grid_x[l+l1+1])[0] AND  $
                               ctr_s.grid_y EQ (ctr_s.grid_y[l+l1+1])[0], count)
print,'count',count
                IF (count EQ 2) THEN BEGIN
                  span_xpoint = 1
                  ; The region to be refinded spans the primary x-point in GRID_S space,
                  ; so need to find out where it is and adjust frac depending on which
                  ; side of the x-point the region is located:
                  frac_s = grid_GetFrac(ctr_s.grid_x[l1:l2],ctr_s.grid_y[l1:l2],dist=dist)
print,'frac_s',frac_s
                  frac_s = frac_s[l+1]
                ENDIF ELSE BEGIN
                  PRINT,'ERROR grid_CheckIntegrity: Clean segment span not identified, need development'
                  STOP
                ENDELSE
              ENDIF
              ; Another convenience check to make sure that the span is on the primary separatrix
              ; which is the most straighforward situation to address:
              IF (NOT span_xpoint) THEN BEGIN
                l = WHERE(state.section_i2 EQ -2, count) 
                IF (count EQ 0) THEN BEGIN
                  PRINT,'ERROR grid_CheckIntegrity: Separatrix malformed'
                  STOP
                ENDIF
                IF (l LT k) THEN BEGIN
                  PRINT,'ERROR grid_CheckIntegrity: Separatrix segment boundary not found'
                  PRINT,'  L = ',l
                  PRINT,'  K = ',k
                  STOP
                ENDIF
              ENDIF	
      
              print,'l,k',l,k
	      
              frac = grid_GetFrac(ctr.grid_x[j1:j2],ctr.grid_y[j1:j2],dist=dist)
	      
              IF (span_xpoint) THEN BEGIN

print,'frac_s',frac_s,frac[j-j1]
print, 'made it here'
                ; Decide which side of the x-point the focus point is on and adjust
                ; the indices accordingly:
                IF (frac[j-j1] LE frac_s) THEN BEGIN
                  j2 = j1 + (l+1 - l1)
                  grid_s1 = ctr_s.grid_s[j1    ]
                  grid_s2 = ctr_s.grid_s[l+1+l1] ; add x-point _S value
                ENDIF ELSE BEGIN
                  j1 = j2 - (l2 - l-1)
                  grid_s1 = ctr_s.grid_s[l+1+l1]
                  grid_s2 = ctr_s.grid_s[j2    ]
                ENDELSE

print,'updated range',j1,j,j2
print,'ctr.grid_s   ',ctr.grid_s
print,'             ',grid_s1,grid_s2

                frac = grid_GetFrac(ctr.grid_x[j1:j2],ctr.grid_y[j1:j2],dist=dist)

              ENDIF ELSE BEGIN
                grid_s1 = ctr.grid_s[j1]
                grid_s2 = ctr.grid_s[j2]
              ENDELSE

              ref_pol_i1  = [ref_pol_i1 ,grid_s1   ]
              ref_pol_i2  = [ref_pol_i2 ,grid_s2   ]
              ref_pol_pos = [ref_pol_pos,frac[j-j1]]
	      


            ENDELSE

          ENDIF

        ENDFOR  ; I (IDENTIFIED VERTICES)

      ENDFOR  ; ICELL

    ENDFOR  ; MAP_J

  ENDFOR  ; ITUBE



  print,'ref_rad_psi1',ref_rad_psi1
  print,'ref_rad_psi2',ref_rad_psi2

  print,'ref_pol_i1 ',ref_pol_i1
  print,'ref_pol_i2 ',ref_pol_i2
  print,'ref_pol_pos',ref_pol_pos

  IF (problem__radial) THEN BEGIN
    state.refinement_rad++
    n = N_ELEMENTS(ref_rad_x1)-1
    IF (state.refinement_rad EQ 1) THEN i = 1 ELSE i = 0
    state = grid_UpdateStructure(state,'ref_rad_x1'  ,ref_rad_x1  [i:n])
    state = grid_UpdateStructure(state,'ref_rad_y1'  ,ref_rad_y1  [i:n])
    state = grid_UpdateStructure(state,'ref_rad_x2'  ,ref_rad_x2  [i:n])
    state = grid_UpdateStructure(state,'ref_rad_y2'  ,ref_rad_y2  [i:n])
    state = grid_UpdateStructure(state,'ref_rad_psi1',ref_rad_psi1[i:n])
    state = grid_UpdateStructure(state,'ref_rad_psi2',ref_rad_psi2[i:n])
  ENDIF

  IF (problem__poloidal) THEN BEGIN
    state.refinement_pol++
    n = N_ELEMENTS(ref_pol_i1)-1
    IF (state.refinement_pol EQ 1) THEN i = 1 ELSE i = 0
    state = grid_UpdateStructure(state,'ref_pol_i1' ,ref_pol_i1 [i:n])
    state = grid_UpdateStructure(state,'ref_pol_i2' ,ref_pol_i2 [i:n])
    state = grid_UpdateStructure(state,'ref_pol_pos',ref_pol_pos[i:n])
  ENDIF

help,state,/struct
print, 'being integrious',state.refinement_pol
;print, 'being integrious',state.refinement_pol++
;stop

;if (state.refinement EQ 3) THEN stop



END
;
; ======================================================================
;
; ======================================================================
;
FUNCTION grid_BuildGrid, c_array, wall, b, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  IF (debug) THEN BEGIN
    PRINT, ' '
    PRINT, '----------------------------------------------------------------------'
    PRINT, ''
    PRINT, 'BUILDING THE GRID'
    PRINT, ''
    PRINT, '----------------------------------------------------------------------'
    PRINT, ' '
  ENDIF

  tags  = STRUPCASE(TAG_NAMES(c_array))
  nctr  = N_ELEMENTS(tags)

  IF (KEYWORD_SET(debug)) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
    OPLOT,wall.x,wall.y,color=Truecolor('Black')
    FOR i = 0, nctr-1 DO BEGIN
      ctr = grid_ExtractStructure(c_array,tags[i])      
      OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
      OPLOT,ctr.grid_x,ctr.grid_y,color=Truecolor('Pink'),PSYM=6
      npts = N_ELEMENTS(ctr.x)
      XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('Darkgray')
    ENDFOR
  ENDIF


;
; ----------------------------------------------------------------------
; COLLECTION SOME INFORMATION ABOUT SEPARATRICES
;
  ; isep[1] = primary separatrix
  ; isep[2] = secondary separatrix, with x-point located on the vessel wall
  ; isep[3] = secondary separatrix, x-point inside the vessel, first leg
  ; isep[4] = secondary separatrix, x-point inside the vessel, second leg

  isep = MAKE_ARRAY(4,  /INT,VALUE=-1)
  ixpt = MAKE_ARRAY(4,2,/INT,VALUE=-1)

  isearch = [1, 2, 2, 2]
  iresult = [2, 1, 1, 1]

  FOR ictr = 1, nctr DO BEGIN

    ctr = grid_ExtractStructure(c_array,tags[ictr-1])      

    IF (ctr.separatrix EQ 0) THEN CONTINUE

    i = ctr.separatrix - 1

    isep[i] = ictr

    IF (state.geometry EQ param.LIMITED) THEN BEGIN
      ixpt[i,*] = [0,N_ELEMENTS(ctr.grid_x)-1]
    ENDIF ELSE BEGIN
      ixpt[i,*] = WHERE(ctr.grid_x EQ b.x[b.null_i[isearch[i]]] AND  $
                        ctr.grid_y EQ b.y[b.null_j[isearch[i]]], count)
      IF (count NE iresult[i]) THEN BEGIN
        PRINT,'ERROR grid_BuildGrid: Mangled primary x-point search'
        PRINT,'  ICTR      =',ictr
        PRINT,'  TAGS      =',tags[ictr-1]
        PRINT,'  SEPATATRIX=',i+1
        PRINT,'  IXPT      =',ixpt[i]
        STOP 
      ENDIF
    ENDELSE

  ENDFOR

  print,'separatrix data..............'
  print,isep
  print,ixpt[*,0]
  print,ixpt[*,1]

;
; ----------------------------------------------------------------------
; IDENTIFY THE CORE CONTOURS
;
  ; Identify the core contours:
  ilist = -1
  FOR ictr = nctr, 1, -1 DO BEGIN
    ctr = grid_ExtractStructure(c_array,tags[ictr-1])          
    IF (ctr.region EQ param.CORE) THEN ilist = [ilist,ictr]
  ENDFOR
  ilist = ilist[1:N_ELEMENTS(ilist)-1] ; Get rid of first entry

  print,'core contours', ilist

;
; ----------------------------------------------------------------------
; FIND WHERE TO START THE GRID WHEN FOLLOWING THE WALL
;
  ; Start with the contour that's just outside the low-index primary
  ; strike-point (the inner target strike-point for lower single null):
  ctr  = grid_ExtractStructure(c_array,'contour1')  ; separatrix
  ictr = (ctr.map_out[0])[0]  
  iwall1 = WHERE(wall.ptc EQ ictr AND wall.ptt EQ 1, count)
  IF (count NE 1) THEN BEGIN
    PRINT,'ERROR grid_BuildGrid: Something is wrong with the wall'
    PRINT,'  IWALL1=',iwall1
    STOP 
  ENDIF
  istart = iwall1

  grid_n = 0

  core_active = 1

  i1 = 0
  i2 = 0

  icore = 0

  tube_n  = 0
;
; ----------------------------------------------------------------------
; FOLLOW THE WALL AND BUILD THE GRID RINGS
;
  WHILE (1) DO BEGIN  

    IF (core_active) THEN BEGIN
;  
;     --------------------------------------------------------------------
;     PROCESS THE CORE RINGS FIRST
;  
      IF (icore EQ N_ELEMENTS(ilist)-1) THEN BEGIN
        ictr1 = ilist[icore]
        ictr2 = 1
        i2 = ixpt[0,0]  ; Set the outer starting point to the x-point
      ENDIF ELSE BEGIN
        ictr1 = ilist[icore  ]
        ictr2 = ilist[icore+1]
      ENDELSE
      ctr1 = grid_ExtractStructure(c_array,tags[ictr1-1])      
      ctr2 = grid_ExtractStructure(c_array,tags[ictr2-1])      

      print,'new contour ---core---- ',ictr1,'-----------',FORMAT='(A,I6,A)'
      print,'            ----------- ',ictr2,'-----------',FORMAT='(A,I6,A)'

    ENDIF ELSE BEGIN
;  
;     --------------------------------------------------------------------
;     SCAN BACKWARD ALONG THE WALL TO FIND THE NEIGHBOURING CONTOUR
;  
      print,'new contour ----------- ',iwall1,ictr1,'-----------',FORMAT='(A,2I6,A)'

      iwall2 = iwall1
      WHILE (1) DO BEGIN
        iwall2--
        IF (iwall2 LT 0) THEN iwall2 = N_ELEMENTS(wall.ptc) - 1
        IF (wall.ptc[iwall2] NE -1 AND (wall.ptt[iwall2] EQ 1 OR  $
                                        wall.ptt[iwall2] EQ 3)) THEN BREAK
      ENDWHILE
      ictr2 = (wall.ptc[iwall2])[0]
      ctr2  = grid_ExtractStructure(c_array,tags[ictr2-1])      
   
      print,'            ----------- ',iwall2,ictr2,'-----------',FORMAT='(A,2I6,A)'

      ; Special case when the neighbouring starting location is a tangency point:
      IF (wall.ptt[iwall2] EQ 3) THEN BEGIN
        print,'   scramble!'
        dummy = MIN( SQRT( ((wall.pt1[0,iwall2])[0]-ctr2.grid_x)^2 +  $
                           ((wall.pt1[1,iwall2])[0]-ctr2.grid_y)^2), i2)
        print,wall.pt1[0:1,iwall2]
        print,ctr2.grid_x[i2],ctr2.grid_y[i2],ctr2.grid_s[i2]
      ENDIF

    ENDELSE  ; CORE_ACTIVE block
;  
;   ----------------------------------------------------------------------
;   ASSEMBLE RING
;  
    tube_x1 = -1.0D
    tube_y1 = -1.0D
    tube_x2 = -1.0D
    tube_y2 = -1.0D

    WHILE (1) DO BEGIN
 
      IF (ctr1.separatrix NE 0) THEN BEGIN
        i = ctr1.separatrix - 1

;        print, 'data check',i,ixpt[i,*]

        IF (i EQ 0 AND i1 EQ ixpt[i,0]) THEN BEGIN
          print,'----1 blasterizing a new contour dude----'
          i1 = ixpt[i,1]
        ENDIF

        IF (i EQ 2 AND i1 EQ ixpt[i,0]) THEN BEGIN
          print,'----2 blasterizing a new contour dude----'
          ctr1 = grid_ExtractStructure(c_array,tags[isep[i+1]-1])                
          i1 = ixpt[i+1,0]
        ENDIF

        IF (i EQ 3 AND i1 EQ ixpt[i,0]) THEN BEGIN
          print,'----3 blasterizing a new contour dude----'
          ctr1 = grid_ExtractStructure(c_array,tags[isep[i-1]-1])                
          i1 = ixpt[i-1,0]
        ENDIF

      ENDIF   

;      IF (grid_n EQ 0) THEN BEGIN
        grid_n++
;        print,'i1,i2',i1,i2,N_ELEMENTS(ctr1.grid_x),N_ELEMENTS(ctr2.grid_x)
        IF (core_active) THEN BEGIN
          cell_x = [ctr1.grid_x[i1  ], ctr2.grid_x[i2  ], ctr2.grid_x[i2+1], ctr1.grid_x[i1+1]]
          cell_y = [ctr1.grid_y[i1  ], ctr2.grid_y[i2  ], ctr2.grid_y[i2+1], ctr1.grid_y[i1+1]]
        ENDIF ELSE BEGIN
          cell_x = [ctr2.grid_x[i2  ], ctr1.grid_x[i1  ], ctr1.grid_x[i1+1], ctr2.grid_x[i2+1]]
          cell_y = [ctr2.grid_y[i2  ], ctr1.grid_y[i1  ], ctr1.grid_y[i1+1], ctr2.grid_y[i2+1]]
        ENDELSE
;      ENDIF ELSE BEGIN
;      ENDELSE

      IF (KEYWORD_SET(debug)) THEN  $
        OPLOT,[cell_x,cell_x[0]],  $
              [cell_y,cell_y[0]],color=Truecolor('Lightblue')

      ; Collect the cell data:      
      tube_x1 = [tube_x1,cell_x[0]]
      tube_y1 = [tube_y1,cell_y[0]]
      tube_x2 = [tube_x2,cell_x[1]]
      tube_y2 = [tube_y2,cell_y[1]]

      i1++
      i2++

      IF ((i1 EQ N_ELEMENTS(ctr1.grid_x)-1                       ) OR  $
          (i2 EQ N_ELEMENTS(ctr2.grid_x)-1 AND ctr2.region GE param.PFZ)) THEN BEGIN  ; Special exit for broken PFR contours:

        tube_x1 = [tube_x1[1:N_ELEMENTS(tube_x1)-1],cell_x[3]]
        tube_y1 = [tube_y1[1:N_ELEMENTS(tube_y1)-1],cell_y[3]]
        tube_x2 = [tube_x2[1:N_ELEMENTS(tube_x2)-1],cell_x[2]]
        tube_y2 = [tube_y2[1:N_ELEMENTS(tube_y2)-1],cell_y[2]]

        IF (core_active) THEN BEGIN
          ictr0 = ictr1
          ctr0  = ctr1
          ictr1 = ictr2
          ctr1  = ctr2
          ictr2 = ictr0
          ctr2  = ctr0
        ENDIF

        tube = {  $
          region      : ctr2.region           ,  $ 
          ictr1       : ictr2                 ,  $
          ictr2       : ictr1                 ,  $
          psi1        : ctr2.psi              ,  $
          psi2        : ctr1.psi              ,  $
          separatrix1 : ctr2.separatrix       ,  $
          separatrix2 : ctr1.separatrix       ,  $
          map_in      : ctr2.map_in           ,  $
          map_out     : ctr1.map_out          ,  $
          n           : N_ELEMENTS(tube_x1)-1 ,  $
          x1          : tube_x1               ,  $
          y1          : tube_y1               ,  $
          x2          : tube_x2               ,  $
          y2          : tube_y2    }
        tube_n++
        tag = 'tube' + STRING(tube_n,FORMAT='(I0)')
        IF (KEYWORD_SET(tubes)) THEN tubes = CREATE_STRUCT(tubes,tag,tube) ELSE  $
                                     tubes = CREATE_STRUCT(      tag,tube)

        BREAK
      ENDIF

    ENDWHILE
;
;   --------------------------------------------------------------------
;   ADVANCE TO THE NEXT CONTOUR
;
    i1 = 0
    i2 = 0

    IF (core_active) THEN BEGIN
      icore++

;      print, 'checking core completiong',icore,N_ELEMENTS(ilist)

      ; Switch over to following the wall if the core is complete:
      IF (icore EQ N_ELEMENTS(ilist)) THEN BEGIN
        iwall1 = istart
        ictr1  = (wall.ptc[iwall1])[0]
        ctr1   = grid_ExtractStructure(c_array,tags[ictr1-1])      

        core_active = 0
      ENDIF

    ENDIF ELSE BEGIN

      WHILE (1) DO BEGIN
        status = 0
        iwall1++
        IF (iwall1 GT N_ELEMENTS(wall.ptc)-1) THEN iwall1 = 0

;        print,'>>> ------- checking',iwall1
;        print,'>>> pti',wall.pti[iwall1]
;        print,'>>> ptc',wall.ptc[iwall1]
;        print,'>>> ptt',wall.ptt[iwall1]
;help,wall,/struct
;stop

        IF (wall.ptc[iwall1] NE -1 AND wall.ptt[iwall1] EQ 1) THEN BEGIN
          ictr1 = (wall.ptc[iwall1])[0]
          ctr1  = grid_ExtractStructure(c_array,tags[ictr1-1])              
;          print,'ctr region   ',ictr1,ctr1.region
;          print,'ctr in ',ctr1.map_in ,FORMAT='(A,2I6)'
;          print,'ctr out',ctr1.map_out,FORMAT='(A,2I6)'
          SWITCH ctr1.region OF
            param.SOL    : 
            param.SOL_LFS: 
            param.SOL_HFS: BEGIN
              status = 1
              BREAK
              END
            param.PFZ          : 
            param.PFZ_SECONDARY: BEGIN
              IF (ctr1.map_in[0] NE -1) THEN status = 1 
              BREAK
              END
          ENDSWITCH
        ENDIF 

;       Special case for tangency points in private flux regions:
        IF (wall.ptc[iwall1] NE -1 AND wall.ptt[iwall1] EQ 3) THEN BEGIN
          ictr1 = (wall.ptc[iwall1])[0]
          ctr1  = grid_ExtractStructure(c_array,tags[ictr1-1])              
;          print,'ctr region ah!',ictr1,ctr1.region
;          print,'ctr in ',ctr1.map_in ,FORMAT='(A,2I6)'
;          print,'ctr out',ctr1.map_out,FORMAT='(A,2I6)'
          SWITCH ctr1.region OF
            param.SOL: 
            param.SOL_LFS: 
            param.SOL_HFS: 
            param.PFZ          :
            param.PFZ_SECONDARY: BEGIN
              IF (ctr1.map_in[1] NE -1) THEN BEGIN
                ; Set I1 so that the the ring starts at the tangency point that's
                ; part way down the contour:
                dummy = MIN( SQRT( ((wall.pt1[0,iwall1])[0]-ctr1.grid_x)^2 +  $
                                   ((wall.pt1[1,iwall1])[0]-ctr1.grid_y)^2), i1)
                status = 1
              ENDIF ELSE BEGIN
                IF (ctr1.map_in[1] EQ -1) THEN BEGIN
                  ; Do nothing, because there's no ring outside the section of
                  ; the contour past the tangency point:
                ENDIF ELSE BEGIN
                  print, 'case not handled yet'
                  stop
                ENDELSE
              ENDELSE
              BREAK
              END
          ENDSWITCH
        ENDIF 
      
        IF (status) THEN BREAK
      ENDWHILE

      ; Exit the loop if back around to the start of the primary separatrix:
      IF (iwall1 EQ istart) THEN BREAK

    ENDELSE

  ENDWHILE

;help,tubes,/struct

  result = tubes



  RETURN, result

END
;
; ======================================================================
;
; ======================================================================
;


