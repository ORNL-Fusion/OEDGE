;
;
; ======================================================================
;
FUNCTION grid_FindNullPoints, b, xpoint_zone, axis_x, axis_y, mode, debug=debug

  b_x   = b.x
  b_y   = b.y
  b_pol = b.b_pol

  colors = ['White','Red','Green','Blue','Orange','Purple','Silver']

;  shade_surf,b_pol,b_x,b_y,ax=-45

  dim = SIZE(b_pol,/DIMENSIONS)
  
  min_n = 10 ; 20
  min_b_pol= MAKE_ARRAY(min_n,/FLOAT,VALUE=1.0E+10)
  min_i    = MAKE_ARRAY(min_n,/LONG,VALUE=0)
  min_j    = MAKE_ARRAY(min_n,/LONG,VALUE=0)
  min_dist = (MAX(b_y) - MIN(b_y)) / FLOAT(min_n)

  FOR j = 0, dim[1]-1 DO BEGIN
    FOR i = 0, dim[0]-1 DO BEGIN
      IF (b_x[i] LE xpoint_zone[0] OR b_x[i] GE xpoint_zone[1] OR  $
          b_y[j] LE xpoint_zone[2] OR b_y[j] GE xpoint_zone[3]) THEN CONTINUE

;     Check if the local poloidal field component strength is less than those
;     in the current list of minimum values:
      FOR k = 0, min_n-1 DO BEGIN

        IF (min_b_pol[k] EQ 1.0E+10) THEN  $
          dist = -1.0                ELSE  $
          dist = ABS(b_y[min_j[k]] - b_y[j])

;        print,k,i,j,dist,min_dist,b_y[min_j[k]],min_b_pol[k]

        IF (dist LT min_dist) THEN BEGIN
          IF (b_pol[i,j] LT min_b_pol[k]) THEN BEGIN
            min_b_pol[k] = b_pol[i,j]
            min_i    [k] = i
            min_j    [k] = j
          ENDIF
          BREAK
        ENDIF
      ENDFOR

    ENDFOR
  ENDFOR

  IF (KEYWORD_SET(debug)) THEN BEGIN
    CONTOUR, b_pol, b_x, b_y, NLEVELS=20, /c_labels
    OPLOT, [xpoint_zone[0],xpoint_zone[0],xpoint_zone[1],xpoint_zone[1],xpoint_zone[0]],  $
           [xpoint_zone[2],xpoint_zone[3],xpoint_zone[3],xpoint_zone[2],xpoint_zone[2]],color=TrueColor('Red')
    FOR k = 0, min_n-1 DO BEGIN
      IF ( min_b_pol[k] EQ 1.0E+10) THEN CONTINUE
      OPLOT, [b_x[min_i[k]]], [b_y[min_j[k]]], PSYM=6, color=TrueColor('Blue')
      print,b_x[min_i[k]], b_y[min_j[k]], min_b_pol[k]
    ENDFOR
  ENDIF
 
; Take only those points that were assigned:
  i = WHERE(min_b_pol LT 1.0E+10, n)  
  min_b_pol = min_b_pol[i]
  min_i     = min_i    [i]
  min_j     = min_j    [i]

  print,' '
  print,min_b_pol
  print,min_i
  print,min_j

; Check if there's a local minimum (or an approximate one anyway) by 
; evaluating the derivatives at the null points and looking for a 
; local extrema in both the x and y directions:

;  space_range = 0.05 ; take all points on the x,y grid within 10 cm
  space_range = MIN([ABS(b_x[1] - b_x[0]), ABS(b_y[1] - b_y[0])]) * 5.0

  min_check = MAKE_ARRAY(N_ELEMENTS(min_b_pol),/LONG,VALUE=0)

  FOR k = 0, N_ELEMENTS(min_b_pol)-1 DO BEGIN
    i = WHERE( ABS(b_x - b_x[min_i[k]]) LT space_range )
    j = WHERE( ABS(b_y - b_y[min_j[k]]) LT space_range )
    x =        b_x  [ i[0] : i[N_ELEMENTS(i)-1] ]
    y = REFORM(b_pol[ i[0] : i[N_ELEMENTS(i)-1] , min_j[k] ],N_ELEMENTS(i))
    y = DERIV(x,y)
    dum = WHERE(y LE 0.0, xcount_neg)
    dum = WHERE(y GT 0.0, xcount_pos)
    IF (KEYWORD_SET(debug)) THEN BEGIN
      print,'----------------------------------------'
      print,'derivative x',k
      print,' i:  ',i
      print,' k:  ',j
      print,' y:  ',y
    ENDIF
    x =        b_y  [            j[0] : j[N_ELEMENTS(j)-1] ]
    y = REFORM(b_pol[ min_i[k] , j[0] : j[N_ELEMENTS(j)-1] ],N_ELEMENTS(j))
    y = DERIV(x,y)
    dum = WHERE(y LE 0.0, ycount_neg)
    dum = WHERE(y GT 0.0, ycount_pos)
    IF (xcount_neg GT 0 AND xcount_pos GT 0 AND   $
        ycount_neg GT 0 AND ycount_pos GT 0) THEN min_check[k] = 1
    IF (KEYWORD_SET(debug)) THEN BEGIN
      print,'derivative y'
      print,' y:  ',y
      print,' '
      print,' counts   :',xcount_neg,xcount_pos,ycount_neg,ycount_pos
      print,' min_check:',min_check[k]
      print,' '
    ENDIF
  ENDFOR

  i = WHERE(min_check,n)
  min_b_pol = min_b_pol[i]
  min_i     = min_i    [i]
  min_j     = min_j    [i]

  print,min_b_pol
  print,min_i
  print,min_j

  IF (KEYWORD_SET(debug)) THEN BEGIN
    CONTOUR, b_pol, b_x, b_y, NLEVELS=20, /c_labels
    OPLOT, [xpoint_zone[0],xpoint_zone[0],xpoint_zone[1],xpoint_zone[1],xpoint_zone[0]],  $
           [xpoint_zone[2],xpoint_zone[3],xpoint_zone[3],xpoint_zone[2],xpoint_zone[2]],color=TrueColor('Red')
    FOR k = 0, N_ELEMENTS(min_b_pol)-1 DO BEGIN
      IF ( min_b_pol[k] EQ 1.0E+10) THEN CONTINUE
      OPLOT, [b_x[min_i[k]]], [b_y[min_j[k]]], PSYM=6, color=TrueColor('Blue')
      print,b_x[min_i[k]], b_y[min_j[k]], min_b_pol[k]
    ENDFOR
  ENDIF

;  return,-1

; Calculate the difference in PSI between the o-point in the core (peak value) and the 
; null points, and sort from closest to furthest:

  max_psi = -1.0E+10
  max_i = -1
  max_j = -1
  FOR j = 0, dim[1]-1 DO BEGIN
    FOR i = 0, dim[0]-1 DO BEGIN
;    i = imax MOD dim[0]
;    j = imax / dim[0]
      IF (b_x[i] LE xpoint_zone[0] OR b_x[i] GE xpoint_zone[1] OR  $
          b_y[j] LE xpoint_zone[2] OR b_y[j] GE xpoint_zone[3]) THEN CONTINUE
      IF (b.psi[i,j] GT max_psi) THEN BEGIN
;      IF (b.psi_raw[i,j] GT max_psi) THEN BEGIN
        max_i = i
        max_j = j
        max_psi = b.psi[i,j]
;        max_psi = b.psi_raw[i,j]
      ENDIF
    ENDFOR
  ENDFOR

  psi_diff = max_psi - b.psi[min_i,min_j]
;  psi_diff = max_psi - b.psi_raw[min_i,min_j]
;  psi_diff = psi_max - b.psi_raw[min_i,min_j]
  i = SORT(psi_diff)
  min_b_pol = min_b_pol[i]
  min_i     = min_i    [i]
  min_j     = min_j    [i]

  print,max_psi
  print,psi_diff

  print,min_b_pol
  print,min_i
  print,min_j



  result = CREATE_STRUCT( b,'null_n',n    ,  $
                            'null_i',min_i,  $
                            'null_j',min_j)
;                            'psi_1st_xpoint',b.psi[null_i[1],null_j[1],
;                            'psi_2nd_xpoint',b.psi[null_i[2],null_j[2])

  IF (KEYWORD_SET(debug)) THEN BEGIN
    CONTOUR, b_pol, b_x, b_y, NLEVELS=10, c_labels=[1,1,1,1,1,1,1,1,1,1]
    OPLOT, [xpoint_zone[0],xpoint_zone[0],xpoint_zone[1],xpoint_zone[1],xpoint_zone[0]],  $
           [xpoint_zone[2],xpoint_zone[3],xpoint_zone[3],xpoint_zone[2],xpoint_zone[2]],color=TrueColor('Red')
    OPLOT, [b_x[min_i[0]]], [b_y[min_j[0]]], PSYM=6, color=TrueColor('Red')      ; o-point
    OPLOT, [b_x[min_i[1]]], [b_y[min_j[1]]], PSYM=6, color=TrueColor('Green')   ; primary   x-point
    IF (N_ELEMENTS(min_i) EQ 3) THEN  $
      OPLOT, [b_x[min_i[2]]], [b_y[min_j[2]]], PSYM=6, color=TrueColor('Blue' ) ; secondary x-point
  ENDIF

  RETURN, result

END
;
; ======================================================================
;
FUNCTION grid_RefineSeparatrices, b, debug=debug, xrange=xrange, yrange=yrange

  result = b

  result = CREATE_STRUCT(result,'psi_1st_xpoint',1.0E+10)
  result = CREATE_STRUCT(result,'psi_2nd_xpoint',1.0E+10)

help,b,/struct

  IF (NOT KEYWORD_SET(user_xrange)) THEN user_xrange = [ 3.0,9.0]  ; lame
  IF (NOT KEYWORD_SET(user_yrange)) THEN user_yrange = [-6.0,6.0]

  psi   = b.psi 
  psi_x = b.x
  psi_y = b.y

  FOR j = 1, b.null_n-1 DO BEGIN

    psi_xpt =  b.psi[b.null_i[j],b.null_j[j]]
    psi_step  = MAX(psi) * 0.0001D

    FOR i = 0, 3 DO BEGIN
      psi_start = psi_xpt + 20.0D * psi_step
      psi_end   = psi_xpt - 20.0D * psi_step
      psi_save  = -1.0D0
      length_xpt = 0.0D0
      length_last = 0.0D0
      delta_max  = 0.0D0
      delta_last = 0.0D0
      FOR psi_xpt = psi_start, psi_end, -psi_step DO BEGIN
    
        ctr = grid_ExtractContour(psi, psi_x, psi_y, psi_xpt)

        IF (KEYWORD_SET(debug)) THEN BEGIN
          PLOT,ctr.x,ctr.y,color=Truecolor('Green'),  $
               XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1
        ENDIF

        ibrk = WHERE(ctr.dist GT 0.10)
        nseg = N_ELEMENTS(ibrk) 
        ibrk = [-1,ibrk,ctr.n-1]
    
        length_max = 0.0D
        FOR iseg = 0, nseg DO BEGIN
          IF ((ibrk[iseg+1]-ibrk[iseg]) LT 5) THEN CONTINUE
          x = ctr.x[ibrk[iseg]+1:ibrk[iseg+1]]
          y = ctr.y[ibrk[iseg]+1:ibrk[iseg+1]]
          length = grid_CalcLength(x,y)
          IF (length GT length_max) THEN length_max = length
        ENDFOR
     
        IF (length_last GT 0.0D) THEN BEGIN 
          delta = length_max - length_last 
          IF (delta GT delta_max) THEN BEGIN
            delta_max = delta
            psi_save = psi_xpt
          ENDIF
        ENDIF
    
        length_last = length_max
    
        IF (KEYWORD_SET(debug)) THEN print,'xpt',j,psi_xpt,length_max,psi_save
      ENDFOR
      psi_xpt = psi_save
      psi_step = psi_step * 0.1D
    ENDFOR
    
    CASE j OF
      1: result.psi_1st_xpoint = psi_xpt
      2: result.psi_2nd_xpoint = psi_xpt-psi_step
;      1: result = CREATE_STRUCT(result,'psi_1st_xpoint',psi_xpt)
;      2: result = CREATE_STRUCT(result,'psi_2nd_xpoint',psi_xpt-psi_step)
      ELSE: stop
    ENDCASE


  ENDFOR 

  RETURN,result

END
;
; ======================================================================
;
FUNCTION grid_AddContour, b, wall, scan_params, contour_array, mode,  $
                          psi_val=psi_val, focus_x=focus_x, focus_y=focus_y,  $
                          debug=debug, xrange=xrange, yrange=yrange, rage=rage
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------

  IF (NOT KEYWORD_SET(xrange)) THEN xrange = [ 3.0,9.0]  ; lame
  IF (NOT KEYWORD_SET(yrange)) THEN yrange = [-6.0,6.0]

  IF (NOT KEYWORD_SET(mode)) THEN mode = 0

  result = -1

  geometry    = scan_params.geometry
  process_2nd = scan_params.process_2nd

  tags      = STRUPCASE(TAG_NAMES(contour_array))
  contour_n = N_ELEMENTS(tags)

  psi   = b.psi 
  psi_x = b.x
  psi_y = b.y

  SWITCH mode OF 
    -1: BEGIN  ; Special case: secondary separatrix...
      focus_x = b.x[b.null_i[2]]
      focus_y = b.y[b.null_j[2]]
      BREAK
      END
     1: 
     2: BEGIN
      IF (NOT KEYWORD_SET(focus_x) OR KEYWORD_SET(psi_val) OR  $
          NOT KEYWORD_SET(focus_y)) THEN BEGIN
        PRINT, 'ERROR grid_AddContour: Invalid call for MODE=1,2'      
        STOP
      ENDIF
 
      n = N_ELEMENTS(psi_x)
      i = WHERE(focus_x GE psi_x[0:n-2] AND focus_x LT psi_x[1:n-1], count)
      IF (count EQ 0) THEN stop    
      n = N_ELEMENTS(psi_y)
      j = WHERE(focus_y GE psi_y[0:n-2] AND focus_y LT psi_y[1:n-1], count)
      IF (count EQ 0) THEN stop

      print,focus_y
      print,psi_y[j],psi_y[j+1]

      psi_1 = INTERPOL(psi[*,j  ],psi_x,focus_x)
      psi_2 = INTERPOL(psi[*,j+1],psi_x,focus_x)
      frac = (focus_y - psi_y[j]) / (psi_y[j+1] - psi_y[j])
      psi_val = psi_1 + frac * (psi_2 - psi_1)
         
      BREAK
      END
    ELSE: BEGIN
      PRINT, 'ERROR grid_AddContour: MODE not found'      
      PRINT, 'MODE =',mode
      STOP
      END
  ENDSWITCH

  ctr = grid_ExtractContour(psi, psi_x, psi_y, psi_val)
  ibrk = WHERE(ctr.dist GT 0.10)
  nseg = N_ELEMENTS(ibrk) 
  ibrk = [-1,ibrk,ctr.n-1]

  IF (KEYWORD_SET(debug)) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, /NOERASE
  ENDIF

  FOR iseg = 0, nseg DO BEGIN
    print, 'iseg',iseg,nseg

    IF ((ibrk[iseg+1]-ibrk[iseg]) LT 5) THEN CONTINUE
  
    x = ctr.x[ibrk[iseg]+1:ibrk[iseg+1]]
    y = ctr.y[ibrk[iseg]+1:ibrk[iseg+1]]

    ; Special check for the case where the secondary PFZ wasn't properly defined,
    ; almost certainly when the 2nd x-point is very close to the vessel wall:
    IF (mode EQ -1 AND scan_params.failure_2nd_pfz EQ 1) THEN BEGIN
      mean_y = MEAN(y)
      IF ((geometry EQ LOWER_NULL AND mean_y GT b.y[b.null_j[2]]) OR  $
          (geometry EQ UPPER_NULL AND mean_y LT b.y[b.null_j[2]])) THEN CONTINUE
    ENDIF

    ; Check if contour runs close to the focus point:
    proximity = SQRT ( (x-focus_x)^2 + (y-focus_y)^2 )

    IF (KEYWORD_SET(rage)) THEN BEGIN
      PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
      OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
      OPLOT,x,y,color=Truecolor('Purple')
      OPLOT,[x[0]],[y[0]],color=Truecolor('White'),PSYM=6
      print, 'min proximity', min(proximity)
    ENDIF

    ; Contour must approach within 5 cm of the focus point:
    IF (MIN(proximity,min_i) GT 0.05D) THEN CONTINUE

    ; Refine the grid near this point, for cases where sharp corners are a problem:
    IF (mode EQ 1 OR mode EQ 2) THEN BEGIN
      count = 0
      WHILE (count LT 2 OR (count LT 4 AND MIN(proximity) GT 0.001D)) DO BEGIN
        grid_RefineContour, x, y, min_i,  $
                            debug=debug, xrange=xrange, yrange=yrange
        proximity = SQRT ( (x-focus_x)^2 + (y-focus_y)^2 )
        dummy = MIN(proximity,min_i)
        count++
      ENDWHILE
    ENDIF

    IF (KEYWORD_SET(rage)) THEN BEGIN
      print, 'min proximity ---', min(proximity)
      print, 'mode             ', mode
    ENDIF

    ; Make sure that the closest contour vertex to the focus point is
    ; inside the grid:
 ;   CASE mode OF
 ;     -1: 
 ;      1: IF (NOT grid_PointInPolygon(x[min_i],y[min_i],wall.x,wall.y)) THEN min_i = min_i + 1
 ;      2: IF (NOT grid_PointInPolygon(x[min_i],y[min_i],wall.x,wall.y)) THEN min_i = min_i - 1
 ;   ENDCASE

;   print, 'processing 2'


    IF (NOT KEYWORD_SET(rage)) THEN $
      PLOT,x,y,color=Truecolor('Lightgreen'),  $
           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

    IF (mode EQ 2) THEN BEGIN
      x = x[0:min_i]
      y = y[0:min_i]
    ENDIF ELSE BEGIN
      ; Find the first point that's inside the wall:
      FOR j = min_i, N_ELEMENTS(x)-2 DO BEGIN
        inside = grid_PointInPolygon(x[j],y[j],wall.x,wall.y)      
        IF (inside) THEN BREAK
      ENDFOR
      IF (j EQ N_ELEMENTS(x)-1) THEN BEGIN
        PRINT, 'ERROR grid_AddContour: No points inside wall (1)'
        PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
        OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
        OPLOT,x,y,color=Truecolor('Magenta')
        OPLOT,x,y,color=Truecolor('White'), PSYM=3
        STOP
      ENDIF
      ; Search for the wall intersection:
      FOR i = j, N_ELEMENTS(x)-2 DO BEGIN
        inter = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                   wall.pt1, wall.pt2, 1, status=status)
        IF (status) THEN BREAK
      ENDFOR
      IF (NOT status) THEN BEGIN
        PRINT, 'NO INTERSECTION WITH WALL 2 FOUND!'
;        CONTINUE
      ENDIF ELSE BEGIN
        IF (N_ELEMENTS(inter.i) EQ 1) THEN BEGIN
          x = x[0:i+1]
          y = y[0:i+1]
        ENDIF ELSE BEGIN
          inside = grid_PointInPolygon(x[min_i],y[min_i],wall.x,wall.y)
          IF (inside EQ 0 OR N_ELEMENTS(inter.i) EQ 0) THEN BEGIN
            xrange = [inter.x[0] - 0.3D, inter.x[0] + 0.3D]
            yrange = [inter.y[0] - 0.3D, inter.y[0] + 0.3D]
            PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
            OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
            OPLOT,[x[min_i]],[y[min_i]],color=Truecolor('Magenta'), PSYM=6
            OPLOT,[x[i]],[y[i]],color=Truecolor('Orange'), PSYM=6
            OPLOT,x,y,color=Truecolor('Magenta')
            OPLOT,x,y,color=Truecolor('White'), PSYM=3
            print,'problem here as well, again (1)',N_ELEMENTS(inter.i)
            stop
          ENDIF
      
          PRINT,'WARNING grid_AddContour: Working hard 1'
      
          ; Take the intersection point that's closest:
          dist = SQRT( (x[i]-inter.x)^2 + (y[i]-inter.y)^2 )
          dummy = MIN(dist, j)
          ; Extend the length of the line segment so that it's just (and I mean
          ; just) beyond the wall:
          length = 0.1D * (MAX(dist) - MIN(dist))
          x = [x[0:i], inter.x[j] + length / dist[j] * (inter.x[j] - x[i])]
          y = [y[0:i], inter.y[j] + length / dist[j] * (inter.y[j] - y[i])]
        ENDELSE
      ENDELSE

    ENDELSE

    IF (mode EQ 1) THEN BEGIN
      x = x[min_i:N_ELEMENTS(x)-1]
      y = y[min_i:N_ELEMENTS(y)-1]
    ENDIF ELSE BEGIN
      ; Find the second point that's inside the wall (since the first point
      ; will register a wall intersection according to the search method
      ; used below):
      count = 0
      FOR j = min_i, 0, -1 DO BEGIN
        inside = grid_PointInPolygon(x[j],y[j],wall.x,wall.y)      
        IF (inside AND count EQ 1) THEN BREAK
        IF (inside AND count EQ 0) THEN count = 1
      ENDFOR
      IF (KEYWORD_SET(rage)) THEN BEGIN
        print,'j',j
        print,N_ELEMENTS(x),min_i
        OPLOT,[x[min_i]],[y[min_i]],color=Truecolor('Blue'),PSYM=6
      ENDIF
      IF (j EQ -1) THEN BEGIN
        PRINT, 'ERROR grid_AddContour: No points inside wall (2)'
        STOP
      ENDIF
      IF (KEYWORD_SET(rage)) THEN BEGIN
        print,'j',j
        print,N_ELEMENTS(x),min_i
        OPLOT,[x[j]],[y[j]],color=Truecolor('Blue'),PSYM=6
      ENDIF
      FOR i = j, 0, -1 DO BEGIN
        inter = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                   wall.pt1, wall.pt2, 1, status=status)
        IF (status) THEN BREAK
      ENDFOR
      IF (NOT status) THEN BEGIN
        PRINT, 'NO INTERSECTION WITH WALL 1 FOUND!',iseg
        CONTINUE
      ENDIF ELSE BEGIN
        IF (N_ELEMENTS(inter.i) EQ 1) THEN BEGIN
          x = x[i:N_ELEMENTS(x)-1]
          y = y[i:N_ELEMENTS(y)-1]
          min_i = min_i - i
        ENDIF ELSE BEGIN
          inside = grid_PointInPolygon(x[min_i],y[min_i],wall.x,wall.y)

          IF (inside EQ 0 OR N_ELEMENTS(inter.i) EQ 0) THEN BEGIN
            xrange = [inter.x[0] - 0.3D, inter.x[0] + 0.3D]
            yrange = [inter.y[0] - 0.3D, inter.y[0] + 0.3D]
            PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
            OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
            OPLOT,[x[min_i]],[y[min_i]],color=Truecolor('Magenta'), PSYM=6
            OPLOT,[x[i]],[y[i]],color=Truecolor('Orange'), PSYM=6
            OPLOT,x,y,color=Truecolor('Magenta')
            OPLOT,x,y,color=Truecolor('White'), PSYM=6
            print,'problem here as well, again (2)',N_ELEMENTS(inter.i)
            stop
          ENDIF
      
          PRINT,'WARNING grid_AddContour: Working hard 2'
      
          ; Take the intersection point that's closest:
          dist = SQRT( (x[i+1]-inter.x)^2 + (y[i+1]-inter.y)^2 )
          dummy = MIN(dist, j)
          ; Extend the length of the line segment so that it's just (and I mean
          ; just) beyond the wall:
          length = 0.1D * (MAX(dist) - MIN(dist))
          x = [inter.x[j] + length / dist[j] * (inter.x[j] - x[i+1]), x[i+1:N_ELEMENTS(x)-1]]
          y = [inter.y[j] + length / dist[j] * (inter.y[j] - y[i+1]), y[i+1:N_ELEMENTS(y)-1]]
          min_i = min_i - i
        ENDELSE
      ENDELSE

    ENDELSE

    IF (KEYWORD_SET(debug)) THEN BEGIN
;      OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
      IF (NOT KEYWORD_SET(rage)) THEN OPLOT,x,y,color=Truecolor('White')
      OPLOT,[x[0]],[y[0]],color=Truecolor('Lightgreen'), PSYM=6
      OPLOT,[x[N_ELEMENTS(x)-1]],[y[N_ELEMENTS(y)-1]],color=Truecolor('Green'), PSYM=6
      OPLOT,[focus_x],[focus_y],color=Truecolor('Orange'), PSYM=6
    ENDIF


    ; Decide the region:
    mean_x = MEAN(x)
    mean_y = MEAN(y)
    region = SOL
    IF (psi_val GT b.psi_1st_xpoint) THEN BEGIN
      IF ((geometry EQ LOWER_NULL AND mean_y GT b.y[b.null_j[1]]) OR  $
          (geometry EQ UPPER_NULL AND mean_y LT b.y[b.null_j[1]])) THEN  $
        region = CORE ELSE region = PFZ
    ENDIF
    ; Secondary x-point is inside the vacuum vessel:
    IF (process_2nd GE 0) THEN BEGIN               ; *** THIS CHECK IS NOT PERFECT, SINCE SMALL NEAR-WALL RINGS IN
      IF (psi_val GT b.psi_2nd_xpoint) THEN BEGIN  ;  THE SOL COULD QUALIFY FOR PFZ ***
        IF ((geometry EQ LOWER_NULL AND mean_y GT b.y[b.null_j[2]]) OR  $
            (geometry EQ UPPER_NULL AND mean_y LT b.y[b.null_j[2]])) THEN  $
          region = PFZ
      ENDIF        
    ENDIF
    print,'region',region

    tangent_i = -1
    p1 = [0.0D,0.0D]
    p2 = p1

    ; Define radial boundary associated with the symmetry point, for
    ; the case where the x-point is inside the vessel but no secondary
    ; PFZ was generated:
    IF (mode EQ -1) THEN BEGIN
      IF (scan_params.failure_2nd_pfz EQ 1) THEN BEGIN
        tangent_i = min_i
        p1 = [x[tangent_i],y[tangent_i]]
        p2 = grid_GetOrthogonal(x,y,p1,2,1.0D,span=10)
print, 'data', contour_n+1
print,p1
print,p2
        PLOT,[p1[0],p2[0]],[p1[1],p2[1]],color=Truecolor('Lightgreen'),  $
             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
        separatrix = 2
      ENDIF ELSE BEGIN
        mean_x = MEAN(x)
        IF ((geometry EQ LOWER_NULL AND mean_x LT b.x[b.null_i[1]]) OR  $
            (geometry EQ UPPER_NULL AND mean_x GT b.x[b.null_i[1]])) THEN BEGIN
          separatrix = 3 
          save_x3 = x          
          save_y3 = y
        ENDIF ELSE BEGIN 
          separatrix = 4
          save_x4 = x          
          save_y4 = y
        ENDELSE
        ; Bit of a messy job here of calculating a normal vector representing
        ; the trajectory between the secondary x-point and the core, to be used
        ; later when setting up the poloidal sections:
        IF (KEYWORD_SET(save_x3) AND KEYWORD_SET(save_x4)) THEN BEGIN
          focus_x = b.x[b.null_i[2]]
          focus_y = b.y[b.null_j[2]]
          proximity = SQRT ( (save_x3-focus_x)^2 + (save_y3-focus_y)^2 )
          dummy = MIN(proximity,min_i3)
          proximity = SQRT ( (save_x4-focus_x)^2 + (save_y4-focus_y)^2 )
          dummy = MIN(proximity,min_i4)
          ; Patch together the PFZ regions for the two secondary separatrix contours:
          x2 = [save_x4[0:min_i4],save_x3[min_i3:N_ELEMENTS(save_x3)-1]]
          y2 = [save_y4[0:min_i4],save_y3[min_i3:N_ELEMENTS(save_y3)-1]]
          OPLOT,x2,y2,color=Truecolor('Yellow')
          ; Find an approximation to the normal vector:
          p1 = [x2[min_i4],y2[min_i4]]
          p2 = grid_GetOrthogonal(x2,y2,p1,1,0.2D,span=10)
          OPLOT,[p1[0],p2[0]],[p1[1],p2[1]],color=Truecolor('Lightgreen')
        ENDIF
      ENDELSE
    ENDIF ELSE separatrix = 0

    contour_data = {  $
      state        : 0                 ,  $
      origin       : 2                 ,  $
      separatrix   : separatrix        ,  $
      region       : region            ,  $
      psi          : psi_val           ,  $
      tangent_i    : tangent_i         ,  $
      tangent_p1   : p1                ,  $
      tangent_p2   : p2                ,  $
      focus_x      : focus_x           ,  $
      focus_y      : focus_y           ,  $
      x            : x                 ,  $
      y            : y                 }

    contour_n++      
    name = 'contour' + STRING(contour_n,FORMAT='(I0)')
    IF (contour_n EQ 1) THEN  $
      contour_array = CREATE_STRUCT(              name,contour_data) ELSE  $
      contour_array = CREATE_STRUCT(contour_array,name,contour_data)

    result = contour_n

  ENDFOR ; Scan over contour segments


  RETURN, result

END
;
; ======================================================================
;
FUNCTION grid_InstallXPoints, b, scan_params, contour_array,  $
                              debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------

  tags      = STRUPCASE(TAG_NAMES(contour_array))
  contour_n = N_ELEMENTS(tags)      

  IF (KEYWORD_SET(debug)) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
    FOR i = 1, contour_n DO BEGIN
      ctr = grid_ExtractStructure(contour_array,tags[i-1])      
      OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
      OPLOT,ctr.x,ctr.y,color=Truecolor('White'),PSYM=3
    ENDFOR
  ENDIF



  count_2nd = 0

  FOR ictr = 1, contour_n DO BEGIN    

    tag = 'contour' + STRING(ictr,FORMAT='(I0)')
    ctr = grid_ExtractStructure(contour_array,tag)  

    status = 0
    IF (ctr.separatrix EQ 1) THEN BEGIN
;    IF (ABS(ctr.psi-b.psi_1st_xpoint) LT 1.0D-6) THEN BEGIN
      status = 1
      focus_x = b.x[b.null_i[1]]
      focus_y = b.y[b.null_j[1]]
      clean_range = 0.01  ; parameter
    ENDIF
    IF (ctr.separatrix GE 2) THEN BEGIN
;     Secondary x-point on the vessel wall:
      status = 2
      focus_x = b.x[b.null_i[2]]
      focus_y = b.y[b.null_j[2]]
      clean_range = 0.01  ; parameter
    ENDIF
    IF (ctr.separatrix GE 3) THEN BEGIN
;     Secondary x-point was inside the vessel
      status = 2
      focus_x = b.x[b.null_i[2]]
      focus_y = b.y[b.null_j[2]]
      clean_range = 0.02  ; parameter
      count_2nd++
    ENDIF
    IF (status EQ 0) THEN CONTINUE

    print,tag

    print,'before',ctr.separatrix,ctr.tangent_i
    x = ctr.x
    y = ctr.y      

    print,'status',status

    FOR i = N_ELEMENTS(y)-2, 0, -1 DO BEGIN 
      j = 0
      IF (status EQ 1) THEN BEGIN
        IF ( (y[i] LE focus_y AND y[i+1] GT focus_y)  OR    $
             (y[i] GE focus_y AND y[i+1] LT focus_y) ) THEN BEGIN
          print,'found'
          j = i
        ENDIF
      ENDIF ELSE BEGIN
        IF ( (scan_params.failure_2nd_pfz EQ 0 AND           $
              ((y[i] LE focus_y AND y[i+1] GT focus_y) OR    $
               (y[i] GE focus_y AND y[i+1] LT focus_y))) OR  $
             (scan_params.failure_2nd_pfz EQ 1 AND           $ 
              ((x[i] LE focus_x AND x[i+1] GT focus_x) OR    $
               (x[i] GE focus_x AND x[i+1] LT focus_x))) ) THEN BEGIN
;          if (scan_params.failure_2nd_pfz EQ 1) then begin
;            print, '*** outdated code I think, remove ***'
;            stop
;          endif
          print, 'found'
          j = i
        ENDIF
      ENDELSE

      IF (j NE 0) THEN BEGIN
        IF (x[i] EQ focus_x AND y[i] EQ focus_y) THEN BEGIN
        ENDIF ELSE BEGIN
          x = [x[0:i],focus_x,x[i+1:N_ELEMENTS(x)-1]]
          y = [y[0:i],focus_y,y[i+1:N_ELEMENTS(y)-1]]
          j++
        ENDELSE
        save_j = j
      ENDIF

    ENDFOR

    ; Delete points that are really close to the x-point since the contouring
    ; is usually messed up there:
    proximity = SQRT ( (x-focus_x)^2 + (y-focus_y)^2 )
    distance  = clean_range * grid_Length(x,y)
    IF (status EQ 1) THEN BEGIN
      i = WHERE(proximity LT distance AND x LT focus_x) 
      IF (KEYWORD_SET(debug)) THEN BEGIN
        OPLOT,x[i],y[i],color=Truecolor('White'),PSYM=6
      ENDIF
      frac = (DINDGEN(5) + 1.0D) / 6.0D
      pad1_x = x[i[0]-1] + frac * (focus_x - x[i[0]-1])
      pad1_y = y[i[0]-1] + frac * (focus_y - y[i[0]-1])
      pad2_x = focus_x + frac * (x[i[N_ELEMENTS(i)-1]+1] - focus_x)
      pad2_y = focus_y + frac * (y[i[N_ELEMENTS(i)-1]+1] - focus_y)
      x = [ x[ 0 : i[0]-1 ], pad1_x, focus_x, pad2_x, x[ i[N_ELEMENTS(i)-1]+1 : N_ELEMENTS(x)-1 ] ]
      y = [ y[ 0 : i[0]-1 ], pad1_y, focus_y, pad2_y, y[ i[N_ELEMENTS(i)-1]+1 : N_ELEMENTS(y)-1 ] ]
      proximity = SQRT ( (x-focus_x)^2 + (y-focus_y)^2 )
      i = WHERE(proximity LT distance AND x GT focus_x) 
      IF (KEYWORD_SET(debug)) THEN BEGIN
        OPLOT,x[i],y[i],color=Truecolor('White'),PSYM=6
      ENDIF
      pad1_x = x[i[0]-1] + frac * (focus_x - x[i[0]-1])
      pad1_y = y[i[0]-1] + frac * (focus_y - y[i[0]-1])
      pad2_x = focus_x + frac * (x[i[N_ELEMENTS(i)-1]+1] - focus_x)
      pad2_y = focus_y + frac * (y[i[N_ELEMENTS(i)-1]+1] - focus_y)
      x = [ x[ 0 : i[0]-1 ], pad1_x, focus_x, pad2_x, x[ i[N_ELEMENTS(i)-1]+1 : N_ELEMENTS(x)-1 ] ]
      y = [ y[ 0 : i[0]-1 ], pad1_y, focus_y, pad2_y, y[ i[N_ELEMENTS(i)-1]+1 : N_ELEMENTS(y)-1 ] ]
    ENDIF ELSE BEGIN
      i = WHERE(proximity LT distance)
      IF (KEYWORD_SET(debug)) THEN BEGIN
        OPLOT,x[i],y[i],color=Truecolor('White'),PSYM=6
      ENDIF
      n = N_ELEMENTS(i)
      x = [ x[ 0 : MAX([0,i[0]-1]) ], focus_x, x[ MIN([i[n-1]+1,N_ELEMENTS(x)-1]) : N_ELEMENTS(x)-1 ] ]
      y = [ y[ 0 : MAX([0,i[0]-1]) ], focus_y, y[ MIN([i[n-1]+1,N_ELEMENTS(y)-1]) : N_ELEMENTS(y)-1 ] ]
      IF (ctr.separatrix EQ 2) THEN save_j = i[0]
    ENDELSE

    IF (KEYWORD_SET(debug)) THEN BEGIN
      OPLOT,x,y,color=Truecolor('Lightgreen'),PSYM=7
    ENDIF  

    ; Setup an outward facing vector for the secondary x-point so that the 
    ; inner SOL can be divided into inner and outer regions:
    IF (count_2nd EQ 1) THEN BEGIN
      save_x = x
      save_y = y
    ENDIF
    IF (count_2nd EQ 2) THEN BEGIN
      i = WHERE(x      EQ focus_x AND y      EQ focus_y)
      j = WHERE(save_x EQ focus_x AND save_y EQ focus_y) 
      IF (ctr.separatrix EQ 3) THEN BEGIN
        x2 = [save_x[0:j],x     [i+1:N_ELEMENTS(x     )-1]]
        y2 = [save_y[0:j],y     [i+1:N_ELEMENTS(y     )-1]]
      ENDIF ELSE BEGIN 
        x2 = [x     [0:i],save_x[j+1:N_ELEMENTS(save_x)-1]]
        y2 = [y     [0:i],save_y[j+1:N_ELEMENTS(save_y)-1]]
      ENDELSE
      ; Find an approximation to the 'normal' vector:
      p1 = [focus_x,focus_y]
      p2 = grid_GetOrthogonal(x2,y2,p1,1,0.3D,span=10)
      IF (debug) THEN BEGIN
        OPLOT,[p1[0],p2[0]],[p1[1],p2[1]],color=Truecolor('Lightgreen')
        OPLOT,x2,y2                      ,color=Truecolor('Red'),PSYM=6
      ENDIF
      ctr.tangent_p1 = p1
      ctr.tangent_p2 = p2
    ENDIF

    ctr = grid_UpdateStructure(ctr,'x',x)    
    ctr = grid_UpdateStructure(ctr,'y',y)    

    ctr = CREATE_STRUCT(ctr,'null_x',focus_x,'null_y',focus_y)

    ; *** should remove this I think ***
    IF (status EQ 2 AND scan_params.failure_2nd_pfz EQ 1) THEN ctr.tangent_i = save_j

    IF (ctr.separatrix EQ 2) THEN BEGIN
     ctr.tangent_i  = save_j
     ctr.tangent_p1 = [focus_x,focus_y] 
    ENDIF

    print,'after',ctr.separatrix,ctr.tangent_i

    contour_array = grid_UpdateStructure(contour_array,tag,ctr)

  ENDFOR

;stop

  RETURN, contour_array

END
;
; ======================================================================
;
PRO grid_AddWallPoint, wall_pti, wall_ptc, wall_ptt,  $
                       wall_pt1, wall_pt2, i, icontour, itarget, x, y

;help,wall_pti
;help,wall_ptc
;help,wall_ptt
;help,wall_pt1
;help,wall_pt2
;help,i
;help,icontour
;help,itarget
;help,x
;help,y

  n = N_ELEMENTS(wall_pti)

  IF (i LT 0 OR i GT n-1) THEN BEGIN
    PRINT, 'ERROR grid_AddWallPoint: Array bounds violation'
    PRINT, '  I=',i
    PRINT, '  N=',n
    STOP
  ENDIF

  IF (i LT n-1) THEN BEGIN
    wall_pti = [wall_pti[0:i],wall_pti[0],wall_pti[i+1:n-1]]
    wall_ptc = [wall_ptc[0:i],icontour   ,wall_ptc[i+1:n-1]]
    wall_ptt = [wall_ptt[0:i],itarget    ,wall_ptt[i+1:n-1]]
    wall_pt1 = TRANSPOSE([[REFORM(wall_pt1[0,0:i]),x,REFORM(wall_pt1[0,i+1:n-1])],  $
                          [REFORM(wall_pt1[1,0:i]),y,REFORM(wall_pt1[1,i+1:n-1])]])
  ENDIF ELSE BEGIN
    wall_pti = [wall_pti[0:i],wall_pti[0]]
    wall_ptc = [wall_ptc[0:i],icontour   ]
    wall_ptt = [wall_ptt[0:i],itarget    ]
    wall_pt1 = TRANSPOSE([[REFORM(wall_pt1[0,0:i]),x],  $
                          [REFORM(wall_pt1[1,0:i]),y]])
  ENDELSE

  IF (i GT 0) THEN  $
    wall_pt2 = TRANSPOSE([[REFORM(wall_pt2[0,0:i-1]),x,REFORM(wall_pt2[0,i:n-1])],  $
                          [REFORM(wall_pt2[1,0:i-1]),y,REFORM(wall_pt2[1,i:n-1])]]) $
  ELSE  $       
    wall_pt2 = TRANSPOSE([[x,REFORM(wall_pt2[0,0:n-1])],  $
                          [y,REFORM(wall_pt2[1,0:n-1])]])

  grid_ZoneWall, wall_pt1, wall_pt2, nx=2, ny=2

END
;
; ======================================================================
;
FUNCTION grid_FindWallNeighbour, wall_pt1, wall_pt2, wall_pti, index, direction

  result = -1
 
  status = 0

  FOR i = 0, N_ELEMENTS(wall_pti)-1 DO BEGIN
    CASE direction OF
      1: IF (wall_pt2[0,i] EQ wall_pt1[0,index] AND  $
             wall_pt2[1,i] EQ wall_pt1[1,index]) THEN status = 1
      2: IF (wall_pt1[0,i] EQ wall_pt2[0,index] AND  $
             wall_pt1[1,i] EQ wall_pt2[1,index]) THEN status = 1
      ELSE: BEGIN
        PRINT, 'ERROR grid_FindWallNeighbour: Invalid direction specified'
        STOP
        END
    ENDCASE
    IF (status) THEN BREAK
  ENDFOR

  IF (NOT status) THEN BEGIN
    ; This may be because the segments don't close on themselves... make provision...
    PRINT, 'ERROR grid_FindWallNeighbour: No neighbour found'
    STOP
  ENDIF

  IF (ABS(wall_pti[i]) NE ABS(wall_pti[index])) THEN BEGIN
    ; This could happen one day depending on how the wall is loaded in (connected segments
    ; from different files), but I doubt it...
    PRINT, 'ERROR grid_FindWallNeighbour: Neighbour from different segment group'
    STOP
  ENDIF

  result = i

  RETURN, result
END
;
; ======================================================================
;
PRO grid_TrimContours, b, contour_array, wall, kill=kill,  $
                       debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------

  IF (NOT KEYWORD_SET(xrange)) THEN xrange = [ 3.0,9.0]  ; lame
  IF (NOT KEYWORD_SET(yrange)) THEN yrange = [-6.0,6.0]

  result = contour_array

  tags = STRUPCASE(TAG_NAMES(contour_array))
  nctr = N_ELEMENTS(tags)

  IF (KEYWORD_SET(debug)) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
    FOR i = 1, nctr DO BEGIN
      ctr = grid_ExtractStructure(contour_array,tags[i-1])      
      OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
      OPLOT,ctr.x,ctr.y,color=Truecolor('White'),PSYM=3
      npts = N_ELEMENTS(ctr.x)
      XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i),2),color=Truecolor('White')
    ENDFOR
  ENDIF



;***  perhaps need to re-specify wall_pt from wall.x/y here, for a fresh start each time this
;     routine is called

  wall_pt1 = wall.pt1
  wall_pt2 = wall.pt2

  FOR ictr = 0, nctr-1 DO BEGIN

    print, '--------------------'
    print, 'ictr',ictr+1

    tag = tags[ictr]
    ctr = grid_ExtractStructure(contour_array,tag)      

    x = ctr.x
    y = ctr.y
 
;    x[0] = x[1] + 2.0D * (x[0] - x[1])
;    y[0] = y[1] + 2.0D * (y[0] - y[1])

    tangent_i = ctr.tangent_i

    OPLOT,x,y,color=Truecolor('Red'),psym=3
;
;   -------------------------------------------------------------------- 
;   
;
    FOR i = 0, N_ELEMENTS(x)-3 DO BEGIN
      result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                 wall_pt1, wall_pt2, 1, status=status)
      IF (status) THEN BREAK
    ENDFOR
    IF (NOT status) THEN BEGIN
      PRINT,'ERROR grid_TrimContours: No wall intersection found 1'
      STOP
    ENDIF ELSE BEGIN
      print,'shit 1',i,N_ELEMENTS(x)-1
      x = x[i:N_ELEMENTS(x)-1] ; This clean up appears to be required, but I'm not sure why...
      y = y[i:N_ELEMENTS(y)-1]

;      dist = grid_PerpDistance(x[0],y[0],[REFORM(wall_pt1[0,*]),wall_pt1[0,0]],  $
;                                         [REFORM(wall_pt1[1,*]),wall_pt1[1,0]])
;      IF (dist LT 0.001D) THEN BEGIN  ; parameter, threshold distance between the last point on the contour (outside the wall) and the wall
;        length = grid_Length(x[0:1],y[0:1])
;        print, 'point on wall, inner!',x[0],y[0],dist,length
;        frac = (length + 0.001D) / length  ; paremeter, 1 mm extension (see below also)
;;        x[0] = x[1] + frac * (x[0] - x[1])
;;        y[0] = y[1] + frac * (y[0] - y[1])
;;        xspan = 0.01D ; 0.50D ; 0.01D
;;        yspan = 0.05D ; 0.50D ; 0.05D
;;        xrange = [x[0] - xspan, x[0] + xspan]
;;        yrange = [y[0] - yspan, y[0] + yspan]
;;        PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
;;        OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange')
;;        OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange'), PSYM=7
;;        OPLOT,[x[0]],[y[0]],color=Truecolor('Magenta'),PSYM=6
;;        stop
;      ENDIF

      IF (tangent_i NE -1) THEN tangent_i = tangent_i - i
    ENDELSE
;
;   -------------------------------------------------------------------- 
;   
;
;   print,'HUNTING!'
    FOR i = N_ELEMENTS(x)-2, 1, -1 DO BEGIN
      result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                 wall_pt1, wall_pt2, 1, status=status)
      IF (status) THEN BREAK
    ENDFOR
    IF (NOT status) THEN BEGIN
      PRINT,'ERROR grid_TrimContours: No wall intersection found 2'

      xspan = 0.02D ; 0.50D ; 0.01D
      yspan = 0.5D ; 0.50D ; 0.05D
      i = N_ELEMENTS(x)-1
      xrange = [x[i] - xspan, x[i] + xspan]
      yrange = [y[i] - yspan, y[i] + yspan]
      PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
      OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange')
      OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange'), PSYM=7
      OPLOT,x,y,color=Truecolor('Magenta'),PSYM=6

      STOP
    ENDIF ELSE BEGIN
      print,'shit 2',i,N_ELEMENTS(x)-1
      x = x[0:i+1]
      y = y[0:i+1]
;      ; Check that the end point isn't too close to the wall:
;      n = N_ELEMENTS(x)
;      dist = grid_PerpDistance(x[n-1],y[n-1],[REFORM(wall_pt1[0,*]),wall_pt1[0,0]],  $
;                                             [REFORM(wall_pt1[1,*]),wall_pt1[1,0]])
;      IF (dist LT 0.001D) THEN BEGIN  ; parameter, threshold distance between the last point on the contour (outside the wall) and the wall
;        length = grid_Length(x[n-2:n-1],y[n-2:n-1])
;        print, 'point on wall, outer!',x[n-1],y[n-1],dist,length
;        frac = (length + 0.001D) / length
;        x[n-1] = x[n-2] + frac * (x[n-1] - x[n-2])
;        y[n-1] = y[n-2] + frac * (y[n-1] - y[n-2])
;;        stop
;      ENDIF

;      ctr.wall_p2 = [result.x[0],result.y[0]]
    ENDELSE

    OPLOT,x,y,color=Truecolor('Pink'),PSYM=6

    ctr = grid_UpdateStructure(ctr,'x',x)    
    ctr = grid_UpdateStructure(ctr,'y',y)    
    ctr.tangent_i = tangent_i
    contour_array = grid_UpdateStructure(contour_array,tag,ctr)

  ENDFOR

  IF (KEYWORD_SET(kill)) THEN stop

END
;
; ======================================================================
;
PRO grid_UpdateWall, b, contour_array, wall, kill=kill, $
                          debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------


  grid_TrimContours, b, contour_array, wall, kill=kill,  $
                     debug=debug, xrange=xrange, yrange=yrange

;  IF (KEYWORD_SET(kill)) THEN  $
;    grid_Debug,b,contour_array,wall,xrange,yrange



  IF (NOT KEYWORD_SET(xrange)) THEN xrange = [ 3.0,9.0]  ; lame
  IF (NOT KEYWORD_SET(yrange)) THEN yrange = [-6.0,6.0]

  result = contour_array

  tags = STRUPCASE(TAG_NAMES(contour_array))
  nctr = N_ELEMENTS(tags)


;***  perhaps need to re-specify wall_pt from wall.x/y here, for a fresh start each time this
;     routine is called  -- DO NOT DO THIS! since it clashes with some assumptions is PROCESSWALL


  wall_pti = wall.pti
  wall_ptc = wall.ptc
  wall_ptt = wall.ptt
  wall_pt1 = wall.pt1
  wall_pt2 = wall.pt2


  FOR ictr = 1, nctr DO BEGIN

    print, 'ictr targets ',ictr,' ',tags[ictr-1]

    tag = tags[ictr-1]
    ctr = grid_ExtractStructure(contour_array,tag)      

    x = ctr.x
    y = ctr.y
 
    tangent_i = ctr.tangent_i

    OPLOT,x,y,color=Truecolor('Red'),psym=3

    i = WHERE(wall_ptc EQ ictr AND wall_ptt EQ 1, count)
    IF (count EQ 0) THEN BEGIN
      result = grid_Intersection([x[0],y[0]], [x[1],y[1]],  $
                                 wall_pt1, wall_pt2, 1, status=status)
      IF (NOT status) THEN BEGIN
        PRINT,'ERROR grid_UpdateWall: No wall intersection found (1)'

        xrange = [x[0] - 0.1D, x[0] + 0.1D]
        yrange = [y[0] - 0.1D, y[0] + 0.1D]
        PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
        OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
        OPLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Orange')
        OPLOT,x,y,color=Truecolor('Magenta')
        OPLOT,x,y,color=Truecolor('White'), PSYM=6
        OPLOT,[x[0]],[y[0]],color=Truecolor('Red'), PSYM=6

        STOP
      ENDIF ELSE BEGIN
        ;print,'damn 1' ; ,result.i,N_ELEMENTS(wall_pti)-1
        IF (N_ELEMENTS(result.i) GT 1) THEN BEGIN
          PRINT,'ERROR grid_UpdateWall: Multiple intersections detected (1)'

          help,result,/struct
          print,result.i
          print,result.s
          print,result.t

          xrange = [result.x[0] - 0.03D, result.x[0] + 0.03D]
          yrange = [result.y[0] - 0.03D, result.y[0] + 0.03D]
          PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
          OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
          OPLOT,wall.pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange')
          OPLOT,wall.pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange'),PSYM=7
          FOR i = 0, N_ELEMENTS(wall_pt1[0,*])-1 DO  $
            XYOUTS,0.5*(wall_pt1[0,i]+wall_pt2[0,i]),  $
                   0.5*(wall_pt1[1,i]+wall_pt2[1,i]),  $
                   STRTRIM(STRING(i),2),color=Truecolor('White')
          OPLOT,x,y,color=Truecolor('Magenta')
          OPLOT,x,y,color=Truecolor('White'), PSYM=6

          STOP
        ENDIF ELSE BEGIN
          i = result.i
          t = result.t        
          status = 1
          IF (t LT 0.5D) THEN j = i ELSE  $
                              j = grid_FindWallNeighbour(wall_pt1,wall_pt2,wall_pti,i,2)
          dist = SQRT((result.x-wall_pt1[0,j])^2 + (result.y-wall_pt1[1,j])^2)
          IF (dist LT 0.0001D) THEN BEGIN
            ;print, 'close 1'
            ;OPLOT,[wall_pt1[0,j]], [wall_pt1[1,j]],color=Truecolor('Lightblue'), PSYM=4
            ;print,wall_pti[j],wall_ptc[j],wall_ptt[j]
            ;print,i,j,t,dist
            IF (wall_ptc[j] NE -1 AND dist GT 1.0D-6) THEN BEGIN
              IF (wall_ptc[j] EQ ictr) THEN BEGIN
                print,'duplicate 1',wall_ptc[j]
                status = 0 
              ENDIF ELSE BEGIN  ; Assume that this point was assigned to this ring on the previous pass
                PRINT, 'WARNING grid_UpdateWall: Contour wall points very close together (1)'
                PRINT, '  DIST = ',dist
                ;xrange = [result.x[0] - 0.1D, result.x[0] + 0.1D]
                ;yrange = [result.y[0] - 0.1D, result.y[0] + 0.1D]
                ;PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
                ;OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
                ;OPLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Orange')
                ;OPLOT,x,y,color=Truecolor('Magenta')
                ;OPLOT,x,y,color=Truecolor('White'), PSYM=6
                ;stop
              ENDELSE
            ENDIF ELSE BEGIN
              print, 'taking over 1'
              print,ictr
              PRINT, '  DIST = ',dist
              wall_ptc[j] = ictr
              wall_ptt[j] = 1
              k = grid_FindWallNeighbour(wall_pt1,wall_pt2,wall_pti,j,1)
              save_pt1 = wall_pt1[*,k]
              save_pt2 = wall_pt2[*,k]
              wall_pt1[0,j] = result.x[0]
              wall_pt1[1,j] = result.y[0]
              wall_pt2[0,k] = result.x[0]
              wall_pt2[1,k] = result.y[0]
              ; Add new point to preserve the wall shape as much as possible:
              new_x = 0.1D * save_pt1[0] + 0.9D * save_pt2[0]
              new_y = 0.1D * save_pt1[1] + 0.9D * save_pt2[1]
;             grid_AddWallPoint, wall_pti, wall_ptc, wall_ptt, wall_pt1, wall_pt2,  $
;                                k, -1 ,-1, new_x, new_y
              status = 0
            ENDELSE
          ENDIF 
          IF (status) THEN  $
            grid_AddWallPoint, wall_pti, wall_ptc, wall_ptt, wall_pt1, wall_pt2,  $
                               i, ictr, 1, result.x[0], result.y[0]
        ENDELSE
      ENDELSE
    ENDIF

    i = WHERE(wall_ptc EQ ictr AND wall_ptt EQ 2, count)
    IF (count EQ 0) THEN BEGIN
      n = N_ELEMENTS(x)
      result = grid_Intersection([x[n-2],y[n-2]], [x[n-1],y[n-1]],  $
                                 wall_pt1, wall_pt2, 1, status=status)
      IF (NOT status) THEN BEGIN
        PRINT,'ERROR grid_UpdateWall: No wall intersection found (2)'
        xspan = 0.11D ; 0.50D ; 0.01D
        yspan = 0.15D ; 0.50D ; 0.05D
        xrange = [x[n-2] - xspan, x[n-2] + xspan]
        yrange = [y[n-2] - yspan, y[n-2] + yspan]
        PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
        OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
        OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange')
        OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange'), PSYM=7
        OPLOT,x,y,color=Truecolor('Magenta')
        OPLOT,x,y,color=Truecolor('White'), PSYM=6
        STOP
      ENDIF ELSE BEGIN
        ;print,'damn 2' ; ,result.i,N_ELEMENTS(wall_pti)-1
        IF (N_ELEMENTS(result.i) GT 1) THEN BEGIN
          PRINT,'ERROR grid_UpdateWall: Multiple intersections detected 2'
          STOP
        ENDIF ELSE BEGIN
          i = result.i
          t = result.t        
          status = 1
          IF (t LT 0.5D) THEN j = i ELSE  $
                              j = grid_FindWallNeighbour(wall_pt1,wall_pt2,wall_pti,i,2)
          dist = SQRT((result.x-wall_pt1[0,j])^2 + (result.y-wall_pt1[1,j])^2)
          IF (dist LT 0.0001D) THEN BEGIN
            ;print, 'close 2'
            ;OPLOT,[wall_pt1[0,j]], [wall_pt1[1,j]],color=Truecolor('Lightblue'), PSYM=4
            ;print,wall_pti[j],wall_ptc[j],wall_ptt[j]
            ;print,i,j,t,dist
            IF (wall_ptc[j] NE -1 AND dist GT 1.0D-6) THEN BEGIN
              IF (wall_ptc[j] EQ ictr) THEN BEGIN
                print,'duplicate 2',wall_ptc[j]
                status = 0 
              ENDIF ELSE BEGIN  ; Assume that this point was assigned to this ring on the previous pass
                PRINT, 'WARNING grid_UpdateWall: Contour wall points very close together (2)'
                PRINT, '  DIST = ',dist

                ;xrange = [result.x[0] - 0.1D, result.x[0] + 0.1D]
                ;yrange = [result.y[0] - 0.1D, result.y[0] + 0.1D]
                ;PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
                ;OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
                ;OPLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Orange')
                ;OPLOT,x,y,color=Truecolor('Magenta')
                ;OPLOT,x,y,color=Truecolor('White'), PSYM=6
                ;stop

              ENDELSE
            ENDIF ELSE BEGIN
              print, 'taking over 2'
              print,ictr
              wall_ptc[j] = ictr 
              wall_ptt[j] = 2
              k = grid_FindWallNeighbour(wall_pt1,wall_pt2,wall_pti,j,1)
              wall_pt1[0,j] = result.x[0]
              wall_pt1[1,j] = result.y[0]
              wall_pt2[0,k] = result.x[0]
              wall_pt2[1,k] = result.y[0]
              status = 0
            ENDELSE
          ENDIF 
          IF (status) THEN  $
            grid_AddWallPoint, wall_pti, wall_ptc, wall_ptt, wall_pt1, wall_pt2,  $
                               i, ictr, 2, result.x[0], result.y[0]
;          IF (0 EQ 1 AND (t LT 0.01D OR t GT 0.99D)) THEN BEGIN
;            print, 'close 2, stopping'
;            stop
;          ENDIF ELSE BEGIN
;           grid_AddWallPoint, wall_pti, wall_ptc, wall_ptt, wall_pt1, wall_pt2,  $
;                              i, ictr, 2, result.x[0], result.y[0]
;          ENDELSE
        ENDELSE
      ENDELSE
    ENDIF

    OPLOT,wall_pt1[0,*], wall_pt1[1,*],color=Truecolor('Orange'), PSYM=7
    OPLOT,wall.pt1[0,*], wall.pt1[1,*],color=Truecolor('Green') , PSYM=7
    OPLOT,x,y,color=Truecolor('Pink')    


;    contour_array = grid_UpdateStructure(contour_array,tag,ctr)


; *** need to clear out double points, removing the old ones, but keeping _pti -ve ones?
; *** make is so the wall zoning can be sent 

  ENDFOR
;
; ----------------------------------------------------------------------
; LOOK AFTER TANGENCY POINTS
; 
  FOR ictr = 1, nctr DO BEGIN

    print, 'ictr tangency ',ictr,' ',tags[ictr-1]

    tag = tags[ictr-1]
    ctr = grid_ExtractStructure(contour_array,tag)      

    IF (ctr.tangent_i EQ -1) THEN CONTINUE

    ; Check if tangency point has already been assigned:
    i = WHERE(wall_ptc EQ ictr AND wall_ptt EQ 3, count)
    IF (count GT 0) THEN CONTINUE

    x         = ctr.x
    y         = ctr.y
    tangent_i = ctr.tangent_i
    p1        = ctr.tangent_p1
    p2        = ctr.tangent_p2
 
    length = SQRT( (p1[0]-p2[0])^2 + (p1[1]-p2[1])^2 )   ; *** This distance -- looking forward and backward -- should	      
    p3 = p1 + 0.01D / length * (p2 - p1)                 ; be the same as the threshold distance for saying how far	       
    p4 = p1 - 0.01D / length * (p2 - p1)                 ; outside the wall the 2nd x-point can be before being forced inside 
                                                          
    ;OPLOT,[p3[0],p4[0]], [p3[1],p4[1]], color=Truecolor('Green') 
    ;OPLOT,x,y,color=Truecolor('Red')    
    ;OPLOT,x,y,color=Truecolor('White'),PSYM=3

    result = grid_Intersection([p3[0],p3[1]], [p4[0],p4[1]],  $
                               wall_pt1, wall_pt2, 1, status=status)
    IF (NOT status) THEN BEGIN
      PRINT,'ERROR grid_UpdateWall: No tangency point wall intersection found'
      xrange = [p3[0] - 0.20D, p3[0] + 0.20D]
      yrange = [p3[1] - 0.20D, p3[1] + 0.20D]
      PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
      OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
      OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange')
      OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange'), PSYM=7
      OPLOT,x,y,color=Truecolor('Magenta')
      OPLOT,x,y,color=Truecolor('White'), PSYM=6
      OPLOT,[p3[0],p4[0]], [p3[1],p4[1]], color=Truecolor('Green')
      STOP
    ENDIF ELSE BEGIN
      ;print,'hot damn 1' ; ,result.i,N_ELEMENTS(wall_pti)-1
      IF (N_ELEMENTS(result.i) GT 1) THEN BEGIN
        PRINT,'ERROR grid_UpdateWall: Multiple tangency intersections detected'
        PRINT,'  I=',result.i
        PRINT,'  S=',result.s
        PRINT,'  T=',result.t
        xrange = [result.x[0] - 0.05D, result.x[0] + 0.05D]
        yrange = [result.y[0] - 0.05D, result.y[0] + 0.05D]
        PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
        OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
        OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange')
        OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange'), PSYM=7
        OPLOT,x,y,color=Truecolor('Magenta')
        OPLOT,x,y,color=Truecolor('White'), PSYM=6
        OPLOT,[p3[0],p4[0]], [p3[1],p4[1]], color=Truecolor('Green')
        STOP
      ENDIF ELSE BEGIN
        i = result.i
        t = result.t        
        status = 1
        IF (t LT 0.5D) THEN j = i ELSE  $
                            j = grid_FindWallNeighbour(wall_pt1,wall_pt2,wall_pti,i,2)
        dist = SQRT((result.x-wall_pt1[0,j])^2 + (result.y-wall_pt1[1,j])^2)
        IF (dist LT 0.0001D) THEN BEGIN
          print, 'close tangency'
          OPLOT,[wall_pt1[0,j]], [wall_pt1[1,j]],color=Truecolor('Lightblue'), PSYM=4
          print,wall_pti[j],wall_ptc[j],wall_ptt[j]
          print,i,j,t,dist
          IF (wall_ptc[j] NE -1) THEN BEGIN
            IF (wall_ptc[j] EQ ictr) THEN BEGIN
              print,'duplicate tangency',wall_ptc[j]
              status = 0 
            ENDIF ELSE  $  ; Assume that this point was assigned to this ring on the previous pass
              PRINT, 'WARNING grid_UpdateWall: Contour wall points very close together'
          ENDIF ELSE BEGIN
            print, 'taking over wall point with tangency'
            print,x[tangent_i],y[tangent_i]
            wall_ptc[j] = ictr
            wall_ptt[j] = 3
            k = grid_FindWallNeighbour(wall_pt1,wall_pt2,wall_pti,j,1)
            wall_pt1[0,j] = x[tangent_i]
            wall_pt1[1,j] = y[tangent_i]
            wall_pt2[0,k] = x[tangent_i]
            wall_pt2[1,k] = y[tangent_i]
            status = 0
          ENDELSE
        ENDIF 
        IF (status) THEN  $
          grid_AddWallPoint, wall_pti, wall_ptc, wall_ptt, wall_pt1, wall_pt2,  $
                             i, ictr, 3, x[tangent_i], y[tangent_i]
      ENDELSE
    ENDELSE

    IF (status) THEN print, 'point added to wall'

    OPLOT,wall_pt1[0,*], wall_pt1[1,*],color=Truecolor('Orange'), PSYM=7




; *** need to clear out double points, removing the old ones, but keeping _pti -ve ones?
; *** make is so the wall zoning can be sent 

  ENDFOR



  wall = grid_UpdateStructure(wall,'pti',wall_pti)
  wall = grid_UpdateStructure(wall,'ptc',wall_ptc)
  wall = grid_UpdateStructure(wall,'ptt',wall_ptt)
  wall = grid_UpdateStructure(wall,'pt1',wall_pt1)
  wall = grid_UpdateStructure(wall,'pt2',wall_pt2)



  grid_ZoneWall, wall.pt1, wall.pt2, debug=debug, xrange=xrange, yrange=yrange


;  n=N_ELEMENTS(wall_pti)
;  FOR j = 0, n-1 DO BEGIN
;    print,j,wall.pti[j],wall.pt1[0,j],wall.pt1[1,j],wall.pt2[0,j],wall.pt2[1,j]
;  ENDFOR

  RETURN
END
;
; ======================================================================
;
PRO grid_ProcessWall, b, wall, scan_params, contour_array, $
                      debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------

  LO_INDEX = 1
  HI_INDEX = 2
  TANGENCY = 3
  OUTSIDE  = -1

  result = -1

  n = N_ELEMENTS(wall.x)


  wall_pti = wall.pti
  wall_ptc = wall.ptc
  wall_ptt = wall.ptt
  wall_pt1 = wall.pt1
  wall_pt2 = wall.pt2

  ; Find out how many groups of wall segments there are:
  FOR ngrp = 1, 100 DO BEGIN
    j = WHERE(wall_pti EQ ngrp, count)
    IF (count EQ 0) THEN BREAK
  ENDFOR
  ngrp = ngrp - 1

  print, 'number of wall groups', ngrp

  FOR igrp = 1, ngrp DO BEGIN

    IF (igrp EQ 1) THEN BEGIN
      ; Start at the separatrix segment:
      iwall = WHERE(wall_ptc EQ 1 AND wall_ptt EQ 1)
    ENDIF ELSE BEGIN
      print, 'not ready for multiple wall groups'
      stop
    ENDELSE

    location       = LO_INDEX
    tangent_active = 0

    i = WHERE(wall_pti EQ igrp, nseg)
    cnt = 0
    print,'nseg',nseg,N_ELEMENTS(wall_pti)

    WHILE (cnt LT nseg) DO BEGIN
      ; print, 'wall',iwall

      iwall = grid_FindWallNeighbour(wall_pt1,wall_pt2,wall_pti,iwall,2)

      CASE wall_ptt[iwall] OF 
        LO_INDEX: BEGIN
          i = WHERE(wall_ptt EQ -2, count)
          IF (location EQ LO_INDEX OR tangent_active) THEN BEGIN
            IF (count GT 0) THEN wall_ptt[i] = LO_INDEX
          ENDIF ELSE BEGIN
            IF (count GT 0) THEN wall_ptt[i] = OUTSIDE        
          ENDELSE
          location = LO_INDEX
          tangent_active = 0
          END
        HI_INDEX: BEGIN
          i = WHERE(wall_ptt EQ -2, count)
          IF (location EQ HI_INDEX OR tangent_active) THEN BEGIN
            IF (count GT 0) THEN wall_ptt[i] = HI_INDEX
          ENDIF ELSE BEGIN
            IF (count GT 0) THEN wall_ptt[i] = OUTSIDE   
          ENDELSE
          location = HI_INDEX
          tangent_active = 0
          END
        TANGENCY: BEGIN
          i = WHERE(wall_ptt EQ -2, count)
          IF (tangent_active) THEN  $
            IF (location EQ LO_INDEX) THEN location = HI_INDEX ELSE  $
                                           location = LO_INDEX
          IF (count GT 0) THEN wall_ptt[i] = location
          tangent_active = 1
          END
        OUTSIDE: wall_ptt[iwall] = -2
      ENDCASE

      cnt++
    ENDWHILE

;   OPLOT,wall_pt1[0,*], wall_pt1[1,*],color=Truecolor('Orange'), PSYM=7
;   OPLOT,wall.pt1[0,*], wall.pt1[1,*],color=Truecolor('Green') , PSYM=7

  ENDFOR

;  n=N_ELEMENTS(wall_pti)
;  FOR j = 0, n-1 DO BEGIN
;    print,j,wall_ptc[j],wall_ptt[j]
;  ENDFOR




  PLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Yellow'),  $
       XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
  PLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Orange'), PSYM=7,  $
       XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

;  i = WHERE(wall_ptt EQ 2, count) 
;  IF (count GT 0) THEN  $
;    PLOT,wall_pt1[0,i],wall_pt1[1,i],color=Truecolor('Pink'), PSYM=7,  $
;         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

  i = WHERE(wall_ptt EQ 1, count) 
  IF (count GT 0) THEN  $
    PLOT,wall_pt1[0,i],wall_pt1[1,i],color=Truecolor('Pink'), PSYM=7,  $
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

  i = WHERE(wall_ptt EQ 2, count) 
  IF (count GT 0) THEN  $
    PLOT,wall_pt1[0,i],wall_pt1[1,i],color=Truecolor('Magenta'), PSYM=7,  $
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

  i = WHERE(wall_ptt EQ 3, count) 
  IF (count GT 0) THEN  $
    PLOT,wall_pt1[0,i],wall_pt1[1,i],color=Truecolor('Lightgreen'),PSYM=7,  $
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

;
; Calculate the angles associated with wall vertices:
;
  print,'filling in corners'

  n = N_ELEMENTS(wall_pti)
 
  wall_angle = MAKE_ARRAY(n,/DOUBLE,VALUE=0.0D)

  FOR i = 0, n-1 DO BEGIN
    j = grid_FindWallNeighbour(wall_pt1,wall_pt2,wall_pti,i,1)
    k = grid_FindWallNeighbour(wall_pt1,wall_pt2,wall_pti,i,2)

    va = wall_pt1[*,j]
    vb = wall_pt1[*,i]
    vc = wall_pt1[*,k]

;    print,va
;    print,vb
;    print,vc

;    PLOT,[vb[0]],[vb[1]],color=Truecolor('Blue'),PSYM=6,  $
;         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

    sa = SQRT( (vb[0]-vc[0])^2 + (vb[1]-vc[1])^2 )
    sb = SQRT( (va[0]-vc[0])^2 + (va[1]-vc[1])^2 )
    sc = SQRT( (va[0]-vb[0])^2 + (va[1]-vb[1])^2 )

    IF (sa GT 0.01D OR sc GT 0.01D) THEN BEGIN
      cos_theta = (sa^2 + sc^2 - sb^2) / (2.0D * sa * sc)
 
      IF ((1.0D0 - ABS(cos_theta)) GT 1.0E-06) THEN  $
        theta = ABS(ACOS(cos_theta) * 180.0D / DOUBLE(!PI)) ELSE  $
        theta = 180.0D    
    ENDIF ELSE theta = 180.0D

    wall_angle[i] = theta

    ; print,i,'yellow',wall_ptt[i],wall_angle[i],cos_theta
  ENDFOR




; i = WHERE(wall_angle LT 150.0D, count)
; IF (count GT 0) THEN  $
;   PLOT,wall.pt1[0,i],wall.pt1[1,i],color=Truecolor('Blue'), PSYM=7,  $
;        XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE


;
; Add a contour for any sharp angles that are inside the grid:
;
  n = N_ELEMENTS(wall_pti)
  count = 0
  FOR i = 0, n-1 DO BEGIN

    IF (wall_angle[i] GT 150.0D) THEN CONTINUE
    IF (wall_ptc  [i] NE -1    ) THEN CONTINUE
    IF (wall_ptt  [i] EQ -1    ) THEN CONTINUE

    count = count + 1

    PLOT,[wall.pt1[0,i]],[wall.pt1[1,i]],color=Truecolor('Blue'),PSYM=1,  $
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
    XYOUTS, wall.pt1[0,i], wall.pt1[1,i], STRTRIM(STRING(count),2), color=Truecolor('White')
    print, '----- working the wall', count,' -----'

    focus_x = wall.pt1[0,i]
    focus_y = wall.pt1[1,i]


;    if (count eq 8) then rage = 1 else rage = 0


    status = grid_AddContour(b, wall, scan_params, contour_array, wall_ptt[i], $
                             focus_x=focus_x, focus_y=focus_y,  $
                             debug=debug, xrange=xrange, yrange=yrange, rage=rage)

;    if (count eq 8) then stop

    IF (status GE 1) THEN wall_ptc[i] = status

    ; Need to adjust the wall point to match the new contour:
    tags = STRUPCASE(TAG_NAMES(contour_array))
    nctr = N_ELEMENTS(tags)
    ctr  = grid_ExtractStructure(contour_array,tags[nctr-1])      
    j    = grid_FindWallNeighbour(wall_pt1,wall_pt2,wall_pti,i,1)
    CASE wall_ptt[i] OF
      1: BEGIN
        wall_pt1[0  ,i] = 0.9D * ctr.x[0] + 0.1D * ctr.x[1]
        wall_pt1[1  ,i] = 0.9D * ctr.y[0] + 0.1D * ctr.y[1]
        wall_pt2[0:1,j] = wall_pt1[0:1,i]
        END
      2: BEGIN
        n = N_ELEMENTS(ctr.x)
        wall_pt1[0  ,i] = 0.9D * ctr.x[n-1] + 0.1D * ctr.x[n-2]
        wall_pt1[1  ,i] = 0.9D * ctr.y[n-1] + 0.1D * ctr.y[n-2]
        wall_pt2[0:1,j] = wall_pt1[0:1,i]
        END
    ENDCASE

;    if (count eq 4) then stop

  ENDFOR



  wall = grid_UpdateStructure(wall,'ptc',wall_ptc)
  wall = grid_UpdateStructure(wall,'ptt',wall_ptt)
  wall = grid_UpdateStructure(wall,'pt1',wall_pt1)
  wall = grid_UpdateStructure(wall,'pt2',wall_pt2)



;  n=N_ELEMENTS(wall_pti)
;  FOR j = 0, n-1 DO BEGIN
;    print,j,wall.pti[j],wall.pt1[0,j],wall.pt1[1,j],wall.pt2[0,j],wall.pt2[1,j],wall_ptc[j]
;  ENDFOR
;stop


  grid_ZoneWall, wall.pt1, wall.pt2, debug=debug, xrange=xrange, yrange=yrange

  RETURN
END
;
; ======================================================================
;
FUNCTION grid_RadialRefinement,c_array, wall, b,  $
                          debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------

  tags  = STRUPCASE(TAG_NAMES(c_array))
  nctr  = N_ELEMENTS(tags)
  nwall = N_ELEMENTS(wall.ptc)

  IF (KEYWORD_SET(debug)) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
    OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
    OPLOT,[wall.pt1[0,*]],[wall.pt1[1,*]],color=Truecolor('Hotpink'),PSYM=6
    FOR i = 0, nctr-1 DO BEGIN
      ctr = grid_ExtractStructure(c_array,tags[i])      
      OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
      npts = N_ELEMENTS(ctr.x)
      XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('White')
    ENDFOR
  ENDIF

help,c_array,/struct

;
; ----------------------------------------------------------------------
; SCAN OVER THE CONTOURS AND INCREASE THE RADIAL RESOLUTION
;  

;
; ----------------------------------------------------------------------
; FIND WHERE TO START THE GRID WHEN FOLLOWING THE WALL
;
  ; Start with the contour that's just outside the low-index primary
  ; strike-point (the inner target strike-point for lower single null):
  ctr1  = grid_ExtractStructure(c_array,'contour1')  ; separatrix
  ictr1 = (ctr1.map_out[0])[0]                       ; step out radially by one contour
  iwall1 = WHERE(wall.ptc EQ ictr1 AND wall.ptt EQ 1, count)
  IF (count NE 1) THEN BEGIN
    PRINT,'ERROR grid_RadialRefinement: Something is wrong with the wall'
    PRINT,'  IWALL1=',iwall1
    STOP 
  ENDIF
  ctr1   = grid_ExtractStructure(c_array,tags[ictr1-1])  
  istart = iwall1

  ictr_new = 0
;
; ----------------------------------------------------------------------
; FOLLOW THE WALL AND BUILD THE GRID RINGS
;
  WHILE (1) DO BEGIN  
;  
;   --------------------------------------------------------------------
;   SCAN BACKWARD ALONG THE WALL TO FIND THE NEIGHBOURING CONTOUR
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

    IF (KEYWORD_SET(debug)) THEN BEGIN
      OPLOT,[wall.pt1[0,iwall1],wall.pt1[0,iwall2]],  $
            [wall.pt1[1,iwall1],wall.pt1[1,iwall2]],color=Truecolor('Lightblue'),PSYM=1
    ENDIF
;  
;   --------------------------------------------------------------------
;   DECIDE IF A NEW CONTOUR IS REQUIRED  
    print,ctr1.psi,ctr2.psi

    psi_val = 0.5D * (ctr1.psi + ctr2.psi)

    IF (1) THEN BEGIN

      print, 'psi_val = ',psi_val

      ctr = grid_ExtractContour(b.psi, b.x, b.y, psi_val)

      ibrk = WHERE(ctr.dist GT 0.10)  ; PARAMETER
      nseg = N_ELEMENTS(ibrk) 
      ibrk = [-1,ibrk,ctr.n-1]
;
;     ------------------------------------------------------------------
;     LOOP OVER INDIVIDUAL CONTOUR SEGMENTS AND EXAMINE
;  
      new_contour = 0

      FOR iseg = 0, nseg DO BEGIN

        print,'----------------------------------------'
        print,'contour seg:',iseg,ibrk[iseg]+1,ibrk[iseg+1]

        ; Not a real segment so skip it, may need to strengthen this check:
        IF (ibrk[iseg]+2 EQ ibrk[iseg+1] OR  $
            ibrk[iseg]+1 EQ ibrk[iseg+1] OR  $
            ibrk[iseg]   EQ ibrk[iseg+1]) THEN CONTINUE  

        x = ctr.x[ibrk[iseg]+1:ibrk[iseg+1]]
        y = ctr.y[ibrk[iseg]+1:ibrk[iseg+1]]

        IF (KEYWORD_SET(debug)) THEN  $
          PLOT,[x,x[0]],[y,y[0]],color=Truecolor('Silver'), PSYM=3,  $
               XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
;
;       ------------------------------------------------------------------
;       TEST IF THE CONTOUR IS IN THE RIGHT PLACE
;  
        FOR i = 0, N_ELEMENTS(x)-2 DO BEGIN  
;         Find intersections with the wall segments:
          result = grid_Intersection([x[i],y[i]],[x[i+1],y[i+1]],  $
                                     wall.pt1[*,iwall2:iwall1],    $
                                     wall.pt2[*,iwall2:iwall1], 0, status=status)
          IF (KEYWORD_SET(debug)) THEN  $
            PLOT,[x[i]],[y[i]],color=Truecolor('Orange'), PSYM=3,  $
                 XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
          IF (status) THEN BREAK
        ENDFOR
;
;       --------------------------------------------------------------------
;       THE CONTOUR INTERSECTS THE WALL IN THE RIGHT PLACE, SO TALLY-HO!
;  
        IF (status) THEN BEGIN
;
;         ------------------------------------------------------------------
;         TRIM THE OTHER END
;  
          x = x[i:N_ELEMENTS(x)-1]
          y = y[i:N_ELEMENTS(y)-1]
          FOR i = 1, N_ELEMENTS(x)-2 DO BEGIN  
            result = grid_Intersection([x[i],y[i]],[x[i+1],y[i+1]],  $
                                       wall.pt1,wall.pt2,1,status=status)
            IF (status) THEN BREAK
          ENDFOR
          IF (NOT status) THEN BEGIN
            PRINT,'ERROR grid_RadialRefinement: Second wall intersection not found'
            STOP
          ENDIF 
          x = x[0:i+1]
          y = y[0:i+1]

          IF (KEYWORD_SET(debug)) THEN  $
            PLOT,x,y,color=Truecolor('Cyan'), PSYM=3,  $
                 XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
;
;         ------------------------------------------------------------------
;         ADD THE CONTOUR TO THE LIST AND EXIT THE LOOP
;  
          ctr = {  $
            state      : 0           ,  $
            origin     : 0           ,  $
            separatrix : 0           ,  $
            region     : ctr2.region ,  $  ; is this OK? probably not
            psi        : psi_val     ,  $
            tangent_i  : -1          ,  $
            tangent_p1 : [0.0D,0.0D] ,  $
            tangent_p2 : [0.0D,0.0D] ,  $
            x          : x           ,   $
            y          : y          }

          ictr_new++
          tag = 'contour' + STRING(nctr+ictr_new,FORMAT='(I0)')

          c_array = CREATE_STRUCT(c_array,tag,ctr)

          new_contour = 1
          BREAK
        ENDIF

      ENDFOR

      IF (NOT new_contour) THEN BEGIN
        PRINT,'grid_RadialRefinement: New contour not found'
        STOP
      ENDIF

    ENDIF



    IF (ictr_new EQ 2) THEN BREAK  ; ***** EARLY EXIT ****




;
;   --------------------------------------------------------------------
;   ADVANCE TO THE NEXT CONTOUR
;
    WHILE (1) DO BEGIN
      status = 0
      iwall1++
      IF (iwall1 GT N_ELEMENTS(wall.ptc)-1) THEN iwall1 = 0

      IF (wall.ptc[iwall1] NE -1 AND wall.ptt[iwall1] EQ 1) THEN BEGIN
        ictr1 = (wall.ptc[iwall1])[0]
        ctr1  = grid_ExtractStructure(c_array,tags[ictr1-1])              
        CASE ctr1.region OF
          SOL: status = 1
          PFZ: IF (ctr1.map_in[0] NE -1) THEN status = 1 
        ENDCASE
      ENDIF 

;     Special case for tangency points in private flux regions:
      IF (wall.ptc[iwall1] NE -1 AND wall.ptt[iwall1] EQ 3) THEN BEGIN
        ictr1 = (wall.ptc[iwall1])[0]
        ctr1  = grid_ExtractStructure(c_array,tags[ictr1-1])              
        CASE ctr1.region OF
          SOL: 
          PFZ: BEGIN
            IF (ctr1.map_in[1] NE -1) THEN BEGIN
              status = 1
            ENDIF ELSE BEGIN
              IF (ctr1.map_in[1] EQ -1) THEN BEGIN
                ; Do nothing, because there's no ring outside the section of
                ; the contour past the tangency point:
              ENDIF ELSE BEGIN
                print, 'radial refinemet: case not handeled yet'
                stop
              ENDELSE
            ENDELSE
           END
        ENDCASE
      ENDIF 
    
      IF (status) THEN BREAK
    ENDWHILE

    ; Exit the loop if back around to the start of the primary separatrix:
    IF (iwall1 EQ istart) THEN BREAK

  ENDWHILE

  help,c_array,/struct

  grid_UpdateWall, b, c_array, wall, $
                   debug=debug, xrange=xrange, yrange=yrange

  c_array = grid_BuildRadialMap(c_array, wall,  $
                                debug=debug, xrange=xrange, yrange=yrange)

  RETURN, c_array

END
