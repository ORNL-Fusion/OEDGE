;
;
;
; ======================================================================
;
FUNCTIOn grid_AnalyseBoundary, psi_x, psi_y, psi, xpoint_zone, mode, plots=plots


  colors = ['White','Red','Green','Blue','Orange','Purple','Silver']

  psib_1 = -MAX(psi) * 0.9D0
  psib_2 = -psib_1
  psib_step = (psib_2 - psib_1) / 100.0

  psib = psib_1

;  psib = -0.3

  print,min(psi),max(psi)

  WHILE ( psib LT psib_2 ) DO BEGIN
;
;   Generate the contour data for the separatrix:
;   --------------------------------------------------------------------
    CONTOUR, psi, psi_x, psi_y, LEVELS=[psib], /PATH_DATA_COORDS,  $
             PATH_XY=path_xy, PATH_INFO=path_info

    x = path_xy[0,*]
    y = path_xy[1,*]

    i = WHERE( x GE xpoint_zone[0] AND x LE xpoint_zone[1]  AND  $
               y GE xpoint_zone[2] AND y LE xpoint_zone[3], count)
    IF (count EQ 0) THEN BEGIN
      PRINT, 'no contours in x-point zone - x scan'
;      STOP
    ENDIF
;    x = x[i]
;    y = y[i]

;    i = WHERE( y GE xpoint_zone[2] AND y LE xpoint_zone[3], count)
;    IF (count EQ 0) THEN BEGIN
;      PRINT, 'no contours in x-point zone - y scan'
;      STOP
;    ENDIF
;    x = x[i]
;    y = y[i]

    n = N_ELEMENTS(x)

;   Find the distance between consectutive points:
    dist = MAKE_ARRAY(n,/FLOAT,VALUE=0.0)
    FOR i = 0, n-2 DO dist[i] = SQRT((x[i]-x[i+1])^2 + (y[i]-y[i+1])^2)

;   Use discontinuities in the contours to isolate flux surfaces:
    i = WHERE(dist GT 0.1, count)
    i = [i,n-1]

    print,count,' --- ',i




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

  
;  xpoint_zone = [0.4,0.8,-1.5,1.5]
;  a = grid_ReadEQUFile('/home/ITER/lisgos/divimp/shots/mast/13018_250/sonnet_13018_250.equ')
;  b = grid_AnalyseBoundary(a.x,a.y,a.psi_raw,xpoint_zone,1,/plots)

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


  HELP,a,/struct
END