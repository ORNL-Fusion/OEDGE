;
;
;
;
; ======================================================================
;
FUNCTION grid_AnalyseBoundary, b, xpoint_zone, mode, debug=debug


  psi   = b.psi_raw
  psi_x = b.x
  psi_y = b.y

  psib = b.psi[b.null_i[0],b.null_j[0]]

  print,b.x[b.null_i[0]],b.y[b.null_j[0]]
  print,b.b_pol[b.null_i[0],b.null_j[0]]

  result = CREATE_STRUCT( b,'psi_b',psib)

;stop
  IF (KEYWORD_SET(debug)) THEN BEGIN
    HELP,result,/STRUCT
    CONTOUR, psi, psi_x, psi_y, LEVELS=[psib-0.001,psib,psib+0.001]
  ENDIF

  RETURN, result

END
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
  
  min_n = 20
  min_b_pol= MAKE_ARRAY(min_n,/FLOAT,VALUE=1.0E+10)
  min_i    = MAKE_ARRAY(min_n,/LONG,VALUE=0)
  min_j    = MAKE_ARRAY(min_n,/LONG,VALUE=0)
  min_dist = (MAX(b_y) - MIN(b_y)) / FLOAT(min_n)

;  FOR j = 0, 150 DO BEGIN
  FOR j = 0, dim[1]-1 DO BEGIN
    FOR i = 0, dim[0]-1 DO BEGIN
      IF (b_x[i] LE xpoint_zone[0] OR b_x[i] GE xpoint_zone[1] OR  $
          b_y[j] LE xpoint_zone[2] OR b_y[j] GE xpoint_zone[3]) THEN CONTINUE

;     Check if the local poloidal field component strength is less than those
;     in the current list of minium values:
      FOR k = 0, min_n-1 DO BEGIN

        IF (min_b_pol[k] EQ 1.0E+10) THEN  $
          dist = -1.0 ELSE dist = SQRT( (b_y[min_j[k]]-b_y[j])^2 )

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
    CONTOUR, b_pol, b_x, b_y, NLEVELS=20
    OPLOT, [xpoint_zone[0],xpoint_zone[0],xpoint_zone[1],xpoint_zone[1],xpoint_zone[0]],  $
           [xpoint_zone[2],xpoint_zone[3],xpoint_zone[3],xpoint_zone[2],xpoint_zone[2]],color=TrueColor('Red')
    k = 0
    OPLOT, [b_x[min_i[k]]], [b_y[min_j[k]]], PSYM=6, color=TrueColor('Blue')
    print,b_x[min_i[k]], b_y[min_j[k]]
  ENDIF
 



; Sort the selected B_pol values into ascending order:
  print,min_b_pol
  print,min_i
  print,min_j

;  i = SORT(min_b_pol)
;  min_b_pol = min_b_pol[i]
;  min_i     = min_i    [i]
;  min_j     = min_j    [i]



; Take only those points that were assigned:
  i = WHERE(min_b_pol LT 1.0E+10, n)  
  min_b_pol = min_b_pol[i]
  min_i     = min_i    [i]
  min_j     = min_j    [i]

  print,min_b_pol
  print,min_i
  print,min_j

; Check if there's a local minimum (or an approximate one anyway) by 
; evaluating the derivatives at the null points and looking for a 
; local extrema in both the x and y directions:

  space_range = 0.05 ; take all points on the x,y grid within 10 cm

  min_check = MAKE_ARRAY(N_ELEMENTS(min_b_pol),/LONG,VALUE=0)

  FOR k = 0, N_ELEMENTS(min_b_pol)-1 DO BEGIN
    i = WHERE( ABS(b_x - b_x[min_i[k]]) LT space_range )
    j = WHERE( ABS(b_y - b_y[min_j[k]]) LT space_range )
    x =        b_x  [ i[0] : i[N_ELEMENTS(i)-1] ]
    y = REFORM(b_pol[ i[0] : i[N_ELEMENTS(i)-1] , min_j[k] ],N_ELEMENTS(i))
;    print,SIZE(x,/DIMENSIONS),SIZE(y,/DIMENSIONS)
    y = DERIV(x,y)
    dum = WHERE(y LE 0.0, xcount_neg)
    dum = WHERE(y GT 0.0, xcount_pos)
    IF (KEYWORD_SET(debug)) THEN BEGIN
      print,i,j
      print,'GO derivative x'
      print,y
    ENDIF
    x =        b_y  [            j[0] : j[N_ELEMENTS(j)-1] ]
    y = REFORM(b_pol[ min_i[k] , j[0] : j[N_ELEMENTS(j)-1] ],N_ELEMENTS(j))
;    print,SIZE(x,/DIMENSIONS)
;    print,SIZE(y,/DIMENSIONS)
    y = DERIV(x,y)
    dum = WHERE(y LE 0.0, ycount_neg)
    dum = WHERE(y GT 0.0, ycount_pos)
    IF (KEYWORD_SET(debug)) THEN BEGIN
      print,'   derivative y'
      print,y
      print, 'counts', xcount_neg,xcount_pos, ycount_neg,ycount_pos
    ENDIF
    IF (xcount_neg GT 0 AND xcount_pos GT 0 AND   $
        ycount_neg GT 0 AND ycount_pos GT 0) THEN min_check[k] = 1
    IF (KEYWORD_SET(debug)) THEN BEGIN
      print,'   derivative y'
      print,y
      print, 'counts', xcount_neg,xcount_pos, ycount_neg,ycount_pos
      print, 'min_check',min_check[k]
      print, ''
    ENDIF
  ENDFOR
;  print,min_check

  ilist = WHERE(min_check,n)
  min_b_pol = min_b_pol[ilist]
  min_i     = min_i    [ilist]
  min_j     = min_j    [ilist]

  print,min_b_pol
  print,min_i
  print,min_j

; Calculate the difference in PSI between the o-point in the core (peak value) and the 
; null points, and sort from closest to furthest:

  psi_max = MAX(b.psi_raw,imax)
  i = imax MOD dim[0]
  j = imax / dim[0]
  psi_diff = psi_max - b.psi_raw[min_i,min_j]
  ilist = SORT(psi_diff)
  min_b_pol = min_b_pol[ilist]
  min_i     = min_i    [ilist]
  min_j     = min_j    [ilist]

  result = CREATE_STRUCT( b,'null_n',n,'null_i',min_i,'null_j',min_j )

  IF (KEYWORD_SET(debug)) THEN BEGIN
    OPLOT, [b_x[i]], [b_y[j]], PSYM=6, color=TrueColor('Orange') ; o-point
    OPLOT, [b_x[min_i[0]]], [b_y[min_j[0]]], PSYM=6, color=TrueColor('Purple') ; x-points
    OPLOT, [b_x[min_i[1]]], [b_y[min_j[1]]], PSYM=6, color=TrueColor('Green') ; x-points
;    print,min_b_pol[i],min_i[i],min_j[i]
;    print,b_x[min_i[i]], b_y[min_j[i]]
;    HELP,result,/STRUCT
  ENDIF

  RETURN, result

END

;
; ======================================================================
;
FUNCTION grid_ExtractContour, psi, psi_x, psi_y, psib;, xpoint_zone=xpoint_zone

; Generate the contour data for the separatrix:
; --------------------------------------------------------------------
  CONTOUR, psi, psi_x, psi_y, LEVELS=[psib], /PATH_DATA_COORDS,  $
           PATH_XY=path_xy, PATH_INFO=path_info

  x = path_xy[0,*]
  y = path_xy[1,*]

  IF (KEYWORD_SET(xpoint_zone)) THEN  $
    i = WHERE( x GE xpoint_zone[0] AND x LE xpoint_zone[1]  AND  $
               y GE xpoint_zone[2] AND y LE xpoint_zone[3], count_pts) ELSE  $
    i = WHERE( x NE -999.0 AND y NE -999.0, count_pts)
  
  IF (count_pts EQ 0) THEN BEGIN
    PRINT, 'no contours in x-point zone'
  ENDIF

  n = N_ELEMENTS(x)

; Find the distance between consectutive points:
  dist = MAKE_ARRAY(n,/FLOAT,VALUE=0.0)
  FOR i = 0, n-2 DO dist[i] = SQRT((x[i]-x[i+1])^2 + (y[i]-y[i+1])^2)

  result = { dist : dist, x : x, y : y, n : n, count_pts : count_pts }

  RETURN, result
END
;
;
; ======================================================================
FUNCTION grid_AnalyseBoundary_OLD, psi_x, psi_y, psi, xpoint_zone, mode, plots=plots

  colors = ['White','Red','Green','Blue','Orange','Purple','Silver']

  psib_1 = -MAX(psi) * 0.9D0
  psib_2 = -psib_1
  psib_step = (psib_2 - psib_1) / 100.0

  psib = psib_1

;  psib = -0.3

  print,min(psi),max(psi)

  count_last = 0
  count_max  = 0

  psib_mark1 = -999.0
  psib_mark2 = -999.0

  WHILE ( psib LT psib_2 ) DO BEGIN

;
;   Generate the contour data for the separatrix:
;   --------------------------------------------------------------------
;    CONTOUR, psi, psi_x, psi_y, LEVELS=[psib], /PATH_DATA_COORDS,  $
;             PATH_XY=path_xy, PATH_INFO=path_info
;
;    x = path_xy[0,*]
;    y = path_xy[1,*]
;
;    i = WHERE( x GE xpoint_zone[0] AND x LE xpoint_zone[1]  AND  $
;               y GE xpoint_zone[2] AND y LE xpoint_zone[3], count_pts)
;    
;    IF (count_pts EQ 0) THEN BEGIN
;      PRINT, 'no contours in x-point zone'
;    ENDIF
;
;    n = N_ELEMENTS(x)
;
;   Find the distance between consectutive points:
;    dist = MAKE_ARRAY(n,/FLOAT,VALUE=0.0)
;    FOR i = 0, n-2 DO dist[i] = SQRT((x[i]-x[i+1])^2 + (y[i]-y[i+1])^2)

    result = grid_ExtractContour(psi, psi_x, psi_y, psib, xpoint_zone=xpoint_zone)

    count_pts = result.count_pts
    dist      = result.dist
    n         = result.n   
    x         = result.x
    y         = result.y

;   Use discontinuities in the contours to isolate flux surfaces:
    i = WHERE(dist GT 0.1, count)
    i = [i,n-1]

    print,count,' --- ',i
    print,psib_mark1,psib_mark2,count_max,count,count_last

    IF (count GT 0) THEN BEGIN

      IF (psib EQ psib_1) THEN  $
        PLOT, [xpoint_zone[0],xpoint_zone[0],xpoint_zone[1],xpoint_zone[1],xpoint_zone[0]],  $
              [xpoint_zone[2],xpoint_zone[3],xpoint_zone[3],xpoint_zone[2],xpoint_zone[2]],  $
              XRANGE=[MIN(x),MAX(x)], YRANGE=[MIN(y),MAX(y)]   
          
      FOR j = 0, count DO BEGIN
        IF (j EQ 0) THEN OPLOT, x[0       :i[j]], y[0       :i[j]], COLOR=TrueColor(colors[j]) ELSE  $
                         OPLOT, x[i[j-1]+1:i[j]], y[i[j-1]+1:i[j]], COLOR=TrueColor(colors[j])
;          PLOT, x[0:i[j]], y[0:i[j]],  $
;               XRANGE=[MIN(x),MAX(x)],  $
;               YRANGE=[MIN(y),MAX(y)]   $
      ENDFOR
    ENDIF

    IF (count_pts GT 0) THEN BEGIN
      IF (count GT count_last) THEN BEGIN
        psib_mark1 = psib
        count_max = count
      ENDIF
      IF (count LT count_last AND psib_mark2 EQ -999.0 and psib_mark1 NE -999.0) THEN psib_mark2 = psib
;      IF (count LT count_last AND count EQ count_max-1) THEN psib_mark2 = psib
      count_last = count
    ENDIF  


;    contour,psi,levels=[psib]



    psib = psib + psib_step





    CONTINUE



    RETURN, -1

    x_core = x[0     :i[0]]
    y_core = y[0     :i[0]]
    n_core = N_ELEMENTS(x_core)
    x_ldiv = x[i[0]+1:i[1]]
    y_ldiv = y[i[0]+1:i[1]]
    n_ldiv = N_ELEMENTS(x_ldiv)
    x_udiv = x[i[1]+1:n-1 ]
    y_udiv = y[i[1]+1:n-1 ]
    n_udiv = N_ELEMENTS(x_udiv)

;   Find the spatial extent of the core:
    y_core_cen = MEAN(y_core)
    y_core_min = MIN(y_core)
    y_core_max = MAX(y_core)
;
;   Lower x-point:
;   --------------------------------------------------------------------
    i = WHERE(y_core LT 0.5 * (y_core_cen + y_core_min))  ; Select the lower core points only, to speed things up
  
;   Find the minimum distance between the lower divertor points and the 
;   core points:  
    min_dist = 1.0E+20
    min_j = -1
    min_k = -1
    FOR j = 0, N_ELEMENTS(i)-1 DO BEGIN
      FOR k = 0, n_ldiv-1 DO BEGIN
        dist = SQRT( (x_ldiv[k] - x_core[i[j]])^2 +  $
                     (y_ldiv[k] - y_core[i[j]])^2 )
        IF (dist LT min_dist) THEN BEGIN
          min_dist = dist
          min_j = j
          min_k = k
        ENDIF
      ENDFOR
    ENDFOR
    xpt_x_ldiv = x_core[i[min_j]]
    xpt_y_ldiv = y_core[i[min_j]]
;
;   Upper x-point:
;   --------------------------------------------------------------------

;
;   Make a pretty plot:
;   --------------------------------------------------------------------
    IF (KEYWORD_SET(plots)) THEN BEGIN
;      PLOT, x_core, y_core,  $
;           XRANGE=[MIN(x_ldiv),MAX(x_core)],  $
;           YRANGE=[MIN(y_ldiv),MAX(y_udiv)]
;      OPLOT, x_ldiv, y_ldiv
;      OPLOT, x_udiv, y_udiv
;      OPLOT, [xpt_x_ldiv], [xpt_y_ldiv], PSYM=6
;      OPLOT, [xpt_x_udiv], [xpt_y_udiv], PSYM=6
    ENDIF

;
;
;
    psib = psib + psib_step

  ENDWHILE 


; *** LEF TOFF *** SOMETHIGNS NOT RIGHT!

    CONTOUR, psi, psi_x, psi_y, LEVELS=[psib_mark1], /PATH_DATA_COORDS,  $
             PATH_XY=path_xy, PATH_INFO=path_info
    x = path_xy[0,*]
    y = path_xy[1,*]
    OPLOT, x, y, COLOR=TrueColor('Purple') 

    CONTOUR, psi, psi_x, psi_y, LEVELS=[psib_mark2], /PATH_DATA_COORDS,  $
             PATH_XY=path_xy, PATH_INFO=path_info
    x = path_xy[0,*]
    y = path_xy[1,*]
    OPLOT, x, y, COLOR=TrueColor('Orange') 


    print,psib_mark1,psib_mark2


return,-1
;*** LEFT OFF *** 
;NOT working for MAST case, since psib_mark2 appears to be inside the core still
;May need more sophisticated x-point detection as well, although maybe not...


  psib_1 = psib_mark2
  psib_2 = psib_mark1

  psib = psib_1

  count_last   = 999
  count_adjust = 0

  WHILE ( psib GT psib_2 ) DO BEGIN

    adjust = 0

    result = grid_ExtractContour(psi, psi_x, psi_y, psib, xpoint_zone=xpoint_zone)

    dist = result.dist
    n    = result.n   

;   Use discontinuities in the contours to isolate flux surfaces:
    i = WHERE(dist GT 0.1, count_break)
    i = [0,i,n-1]

    print,'------------------------'

    print,i
    print,count_break

    x = result.x
    y = result.y

    print,x[i[0]  ],x[i[1]],'  ',dist[i[1]-1:i[1]+1]
    print,x[i[1]+1],x[i[2]],'  ',dist[i[2]-1:i[2]+1]
    print,x[i[2]+1]

    print,x[n-1],y[n-1]

    print,count_break,' --groovy-- ',psib

    count = 0

    FOR j = 0, count_break DO BEGIN
 
      IF (j EQ 0) THEN i1 = i[j] ELSE i1 = i[j] + 1
      i2 = i[j+1]

      k = WHERE( x[i1:i2] GE xpoint_zone[0] AND x[i1:i2] LE xpoint_zone[1]  AND  $
                 y[i1:i2] GE xpoint_zone[2] AND y[i1:i2] LE xpoint_zone[3], count_pts)

      IF (count_pts GT 0) THEN count = count + 1

      print,'count_pts',i1,i2,count_pts


    ENDFOR

    print,count,' --      -- ',psib,count_last

    IF (count_adjust MOD 2 EQ 0 AND count LT count_last AND count_last NE 999) THEN adjust = 1
    IF (count_adjust MOD 2 EQ 1 AND count GT count_last AND count_last NE 999) THEN adjust = 1

    print,'adjust',adjust

    count_last = count

    IF (adjust) THEN BEGIN
      count_adjust = count_adjust + 1
      IF (count_adjust EQ 5) THEN BREAK
      psib_step = -0.3 * psib_step
    ENDIF

    psib = psib - psib_step

  ENDWHILE

  CONTOUR, psi, psi_x, psi_y, LEVELS=[psib], /PATH_DATA_COORDS,  $
           PATH_XY=path_xy, PATH_INFO=path_info
  x = path_xy[0,*]
  y = path_xy[1,*]
  OPLOT, x, y, COLOR=TrueColor('Lightgreen') 


return, -1

;
; Collect the results:
; ----------------------------------------------------------------------
  result = {  $
    upper_x : xpt_x_udiv ,  $
    upper_y : xpt_y_udiv ,  $
    lower_x : xpt_x_ldiv ,  $
    lower_y : xpt_y_ldiv }

  RETURN, result

END
;
; ======================================================================
;
FUNCTIOn grid_AnalyseSeparatrix, psi_x, psi_y, psi, psi_sep, mode, plots=plots

;
; Generate the contour data for the separatrix:
; ----------------------------------------------------------------------
  CONTOUR, psi, psi_x, psi_y, levels=[psi_sep], /path_data_coords,  $
           path_xy=path_xy, path_info=path_info

  x = path_xy[0,*]
  y = path_xy[1,*]
  n = N_ELEMENTS(x)

; Find the distance between consectutive points:
  dist = MAKE_ARRAY(n,/FLOAT,VALUE=0.0)
  FOR i = 0, n-2 DO dist[i] = SQRT((x[i]-x[i+1])^2 + (y[i]-y[i+1])^2)

; Divide the contour data into regions - core, lower divertor and upper 
; divertor - using the discontinuities in the distances between points:
  i = WHERE(dist GT 0.1)

  x_core = x[0     :i[0]]
  y_core = y[0     :i[0]]
  n_core = N_ELEMENTS(x_core)
  x_ldiv = x[i[0]+1:i[1]]
  y_ldiv = y[i[0]+1:i[1]]
  n_ldiv = N_ELEMENTS(x_ldiv)
  x_udiv = x[i[1]+1:n-1 ]
  y_udiv = y[i[1]+1:n-1 ]
  n_udiv = N_ELEMENTS(x_udiv)

; Find the spatial extent of the core:
  y_core_cen = MEAN(y_core)
  y_core_min = MIN(y_core)
  y_core_max = MAX(y_core)
;
; Lower x-point:
; ----------------------------------------------------------------------
  i = WHERE(y_core LT 0.5 * (y_core_cen + y_core_min))  ; Select the lower core points only, to speed things up
  
; Find the minimum distance between the lower divertor points and the 
; core points:  
  min_dist = 1.0E+20
  min_j = -1
  min_k = -1
  FOR j = 0, N_ELEMENTS(i)-1 DO BEGIN
    FOR k = 0, n_ldiv-1 DO BEGIN
      dist = SQRT( (x_ldiv[k] - x_core[i[j]])^2 +  $
                   (y_ldiv[k] - y_core[i[j]])^2 )
      IF (dist LT min_dist) THEN BEGIN
        min_dist = dist
        min_j = j
        min_k = k
      ENDIF
    ENDFOR
  ENDFOR
  xpt_x_ldiv = x_core[i[min_j]]
  xpt_y_ldiv = y_core[i[min_j]]
;
; Upper x-point:
; ----------------------------------------------------------------------
  i = WHERE(y_core GT 0.5 * (y_core_cen + y_core_max))

; Find the minimum distance between the lower divertor points and the 
; core points:  
  min_dist = 1.0E+20
  min_j = -1
  min_k = -1
  FOR j = 0, N_ELEMENTS(i)-1 DO BEGIN
    FOR k = 0, n_udiv-1 DO BEGIN
      dist = SQRT( (x_udiv[k] - x_core[i[j]])^2 +  $
                   (y_udiv[k] - y_core[i[j]])^2 )
      IF (dist LT min_dist) THEN BEGIN
        min_dist = dist
        min_j = j
        min_k = k
      ENDIF
    ENDFOR
  ENDFOR
  xpt_x_udiv = x_core[i[min_j]]
  xpt_y_udiv = y_core[i[min_j]]
;
; Make a pretty plot:
; ----------------------------------------------------------------------
  IF (KEYWORD_SET(plots)) THEN BEGIN
    PLOT, x_core, y_core,  $
         XRANGE=[MIN(x_ldiv),MAX(x_core)],  $
         YRANGE=[MIN(y_ldiv),MAX(y_udiv)]
    OPLOT, x_ldiv, y_ldiv
    OPLOT, x_udiv, y_udiv
    OPLOT, [xpt_x_ldiv], [xpt_y_ldiv], PSYM=6
    OPLOT, [xpt_x_udiv], [xpt_y_udiv], PSYM=6
  ENDIF
;
; Collect the results:
; ----------------------------------------------------------------------
  result = {  $
    upper_x : xpt_x_udiv ,  $
    upper_y : xpt_y_udiv ,  $
    lower_x : xpt_x_ldiv ,  $
    lower_y : xpt_y_ldiv }

  RETURN, result

END
;
; ======================================================================
;
PRO grid_Main

  
  xpoint_zone = [0.4,0.8,-1.5,1.5]
  b = grid_ReadEQUFile('~/fuse_data/mast/shots/24867/24867_335.equ')
;  b = grid_ReadEQUFile('~/divimp/shots/mast/24860/carre.24860_240.equ')

;  shade_surf,b.psi_raw,b.x,b.y,ax=0
;stop

  b = grid_FindNullPoints (b,xpoint_zone,1,/debug)
;  b = grid_AnalyseBoundary(b,xpoint_zone,1,/debug)



;  xpoint_zone = [0.4,0.8,-1.5,1.5]
;  b = grid_ReadEQUFile('/home/ITER/lisgos/divimp/shots/mast/13018_250/sonnet_13018_250.equ')
;  b = grid_FindNullPoints(b,xpoint_zone,1,/debug)
;  b = grid_AnalyseBoundary(b,xpoint_zone,1,/plots)

;  xpoint_zone = [2.0,2.4,-0.7,0.7]
;  a = grid_ReadEQUFile('/home/ITER/lisgos/divimp/shots/ts2/600kA_LN_1cm/TS2_600kA_LN.x2.cr.equ')
;  b = grid_AnalyseBoundary(a.x,a.y,a.psi_raw,xpoint_zone,1,/plots)

;  xpoint_zone = [2.0,2.4,-0.7,0.7]
;  a = grid_ReadEQUFile('/home/ITER/lisgos/divimp/shots/ts2/600kA_CDN/TS2_600kA.x2.cr.equ')
;  b = grid_AnalyseBoundary(a.x,a.y,a.psi_raw,xpoint_zone,1,/plots)

;  xpoint_zone = [2.0,2.4,-0.7,0.7]
;  a = grid_ReadEQUFile('/home/ITER/lisgos/divimp/shots/ts/upgrade_1MA_1cm/TS2_1MA_LN.x2.cr.equ')
;  b = grid_AnalyseBoundary(a.x,a.y,a.psi_raw,xpoint_zone,1,/plots)

;  xpoint_zone = [4.0,6.0,-3.8,5.5]
;  a = grid_ReadEQUFile('/home/ITER/lisgos/divimp/shots/iter/i1514/Baseline2008-li0.70.x4.equ')
;  b = grid_AnalyseBoundary(a.x,a.y,a.psi_raw,xpoint_zone,1,/plots)


;  HELP,a,/struct
END