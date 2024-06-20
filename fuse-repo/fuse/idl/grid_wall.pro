;
; ======================================================================
;
; ======================================================================
;
FUNCTION grid_RadialDistribution, c_array, wall, b, region, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

;
; ----------------------------------------------------------------------
; GET THE MAPPING BETWEEN x (R) AND psi AT THE OUTER MIDPLANE
;
  i = b.null_i[0]  ; O-point in the core, i.e. the plasma centre
  j = b.null_j[0]
  nx = N_ELEMENTS(b.x)
  psi   = b.psi[i:nx-1,j]
  psi_x = b.x  [i:nx-1  ]  
;
; ----------------------------------------------------------------------
; PROCESS THE DISTRIBUTION OF CONTOURS FOR EACH REGION
;
  SWITCH region OF
;   --------------------------------------------------------------------
    param.CORE: BEGIN

;      res_n   = 51
;      res_min = 0.002D
;      res_max = 0.10D
;;     res_min = 0.001D ; 0.001D
;;     res_max = 0.050D ; 0.050D
;      res_exp = 2.0D

      res_n   = option.rad_res_core_n   ;101
      res_min = option.rad_res_core_min ;0.001D ; 0.001D
      res_max = option.rad_res_core_max ;0.020D ; 0.050D
      res_exp = option.rad_res_core_exp ;1.0D

      res = DINDGEN(res_n) / DOUBLE(res_n-1)

print,'origin',res     

      ; Apply the distribution function of choice:
;     res = [ 0.0D, res[1:res_n-1]^0.7  ]
;      res = [ 0.0D, EXP(res[1:res_n-1])-1.0D ]

      CASE option.rad_res_core_opt OF
        param.TANH: BEGIN
          x = res[1:N_ELEMENTS(res)-1] * 10.0D 
          shift = -5.0D
          res = [res[0],( EXP(0.5D*res_exp*(x+shift)) - 1.0D ) /  $
                        ( EXP(0.5D*res_exp*(x+shift)) + 1.0D ) ]
          res = [res[0], (res[1:res_n-1] - res[1]) * (res_max - res_min) / MAX(res) + res_min]
          END
        ELSE: BEGIN
          PRINT,'ERROR grid_RadialDistribution: Unrecognised distribution'
          PRINT,'  OPT = ', option.rad_res_sol_opt
          STOP
          END
      ENDCASE



;      plot,res,psym=6

      ; Rescale the point nearest the separatrix (RES = 0.0) so that it
      ; is equal to RES_MIN:
print,'before',res     
      scale = (res_min / (res[1])[0])
      res = res * scale
print,'after ',res     

      ; Limit the maxium spatial resolution:
      FOR i = 0, N_ELEMENTS(res)-1 DO res[i] = MIN([res[i],res_max])
print,'after ',res     
      ; Generate distribution of points in space:
      dist = res
      FOR i = 1, N_ELEMENTS(dist)-1 DO dist[i] = dist[i-1] + res[i]
print,'dist  ',dist
      ; Find out where the separatrix is for the X vs. PSI data:
      delta = ABS(psi - b.psi_1st_xpoint)
      dummy = MIN(delta,i)
     
      FOR i = 0, N_ELEMENTS(psi)-2 DO $
        IF (psi[i] GE b.psi_1st_xpoint AND psi[i+1] LT b.psi_1st_xpoint) THEN BREAK
      IF (i EQ N_ELEMENTS(psi)-1) THEN BEGIN
        PRINT,'ERROR grid_RadialDistribution: Separatrix not found'
        STOP
      ENDIF

      frac = (b.psi_1st_xpoint - psi[i]) / (psi[i+1] - psi[i])
      psi_xsep = psi_x[i] + frac * (psi_x[i+1] - psi_x[i])

print,'raw',psi[i:i+1],b.psi_1st_xpoint,i
print,'x  ',psi_x[i:i+1],psi_xsep
;stop
print,'psi',psi
print,'psi_x-psi_x[i]',psi_x - psi_x[i]
      ; Select the PSI_X and PSI data in the core, rescale PSI_X so that
      ; the origin is at the separatrix, and resort:
      psi_x = psi_xsep - REVERSE([psi_x[0:i],psi_xsep        ])
      psi   =            REVERSE([psi  [0:i],b.psi_1st_xpoint])
print,'psi_x',psi_x
print,'psi',psi
;stop
;plot,psi_x,psi
      ; Find the boundary in PSI and PSI_X:
      psi_boundary = b.psi_1st_xpoint - option.ctr_boundary[region-1]
print,psi_boundary
      FOR i = 0, N_ELEMENTS(psi)-2 DO  $
        IF (psi[i] LE psi_boundary AND psi[i+1] GT psi_boundary) THEN BREAK
      IF (i EQ N_ELEMENTS(psi)-1) THEN BEGIN
        PRINT,'ERROR grid_RadialDistribution: Boundary not found'
        PRINT,'  PSI_BOUNDARY = ',psi_boundary
        PRINT,'  PSI_MIN,_MAX = ',psi[0],psi[N_ELEMENTS(psi)-1]
        STOP
      ENDIF
      frac = (psi_boundary - psi[i]) / (psi[i+1] - psi[i])
      psi_x_boundary = psi_x[i] + frac * (psi_x[i+1] - psi_x[i])
;print,i,N_ELEMENTS(psi)-1
;print,psi_x_boundary
      ; Trim the distribution of contours so that it does not extend beyond the
      ; inner boundary:
      i = WHERE(dist LT psi_x_boundary, count)
      IF (count EQ 0) THEN BEGIN
        PRINT,'ERROR grid_RadialDistribution: Inner boundary too close to separatrix'
        PRINT,'  PSI_X_BOUNDARY = ',psi_x_boundary
        PRINT,'  PSI_X          = ',psi_x
        STOP
      ENDIF

print,'i ',i
      IF (N_ELEMENTS(i) LT 2) THEN  $
        psi_dist = [0.0D                      , psi_x_boundary] ELSE  $  ; no contours inside specified boundary
        psi_dist = [dist[i[0:N_ELEMENTS(i)-2]], psi_x_boundary]
      psi_ctr  = SPLINE( psi_x, psi, psi_dist ,/DOUBLE)              

      n = N_ELEMENTS(psi_dist)
      psi_dist = psi_dist[1:n-1]
      psi_ctr  = psi_ctr [1:n-1]
      izero = -1

print,'done',psi_dist
print,'    ',psi_ctr     

      BREAK
      END
;   --------------------------------------------------------------------
    param.SOL          :
    param.SOL_LFS      : 
    param.SOL_HFS      :
    param.PFZ          :
    param.PFZ_SECONDARY: BEGIN

;      res_n   = 101
;      res_min = 0.001D ; 0.001D
;      res_max = 0.020D ; 0.050D
;      res_exp = 1.0D

      res_n   = option.rad_res_sol_n   ;101
      res_min = option.rad_res_sol_min ;0.001D ; 0.001D
      res_max = option.rad_res_sol_max ;0.020D ; 0.050D
      res_exp = option.rad_res_sol_exp ;1.0D

      res = DINDGEN(res_n) / DOUBLE(res_n-1)
print,res

      ; Apply the distribution function of choice:

print,'res1',res


      CASE option.rad_res_sol_opt OF
        param.TANH: BEGIN
          x = res[1:N_ELEMENTS(res)-1] * 10.0D 
          shift = -5.0D
          res = [res[0],( EXP(0.5D*res_exp*(x+shift)) - 1.0D ) /  $
                        ( EXP(0.5D*res_exp*(x+shift)) + 1.0D ) ]
          res = [res[0], (res[1:res_n-1] - res[1]) * (res_max - res_min) / MAX(res) + res_min]
          END
        ELSE: BEGIN
          PRINT,'ERROR grid_RadialDistribution: Unrecognised distribution'
          PRINT,'  OPT = ', option.rad_res_sol_opt
          STOP
          END
      ENDCASE
print,'x',x
print,'res2',res

      ; Rescale the point nearest the separatrix (RES = 0.0) so that it
      ; is equal to RES_MIN:
print,'before',res     
      scale = res_min / res[1]
      res = res * scale
print,'after ',res     
      ; Limit the maxium spatial resolution:
      FOR i = 0, N_ELEMENTS(res)-1 DO res[i] = MIN([res[i],res_max])

;plot,res,psym=6
;stop

      ; Mirror:
      izero = N_ELEMENTS(res)-1
      res = [REVERSE(-res[1:izero]),res]

      ; Generate distribution of points in space:
      dist = res
      FOR i = izero-1, 0                 , -1 DO dist[i] = dist[i+1] + res[i]
      FOR i = izero+1, N_ELEMENTS(dist)-1     DO dist[i] = dist[i-1] + res[i]

      ; Find out where the separatrix is for the X vs. PSI data:
      delta = ABS(psi - b.psi_1st_xpoint)
      dummy = MIN(delta,i)

      ; Rescale PSI_X so that the origin is at the separatrix, and resort:
      psi_x = psi_x - psi_x[i]

      i = WHERE(dist GT psi_x[0] AND dist LT psi_x[N_ELEMENTS(psi_x)-1], count)
      IF (count EQ 0) THEN BEGIN
        PRINT,'ERROR grid_RadialDistribution: No boundary contours found'
        STOP
      ENDIF

      psi_dist = dist[i]
      psi_ctr  = SPLINE( psi_x, psi, psi_dist ,/DOUBLE)              

      izero = WHERE(psi_dist EQ 0.0D)

      BREAK
      END
;   --------------------------------------------------------------------
  ENDSWITCH



  result = { dist : psi_dist, psi : psi_ctr, izero : izero } 

  RETURN, result

END
;
;
; ======================================================================
;
FUNCTION grid_FindNullPoints, b, xpoint_zone, axis_x, axis_y, mode, debug=debug

print,'input',xpoint_zone, axis_x

  b_x   = b.x
  b_y   = b.y
  b_pol = b.b_pol

  colors = ['Darkgrey','Red','Green','Blue','Orange','Purple','Silver']

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
    CONTOUR, b_pol, b_x, b_y, NLEVELS=20, /c_labels, color=Truecolor('Black')
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
    CONTOUR, b_pol, b_x, b_y, NLEVELS=20, /c_labels, color=Truecolor('Black')
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

  print,'max_psi   ',max_psi
  print,'psi_diff  ',psi_diff
  print,'min_b_pol ',min_b_pol
  print,'b.psi     ',b.psi[min_i,min_j]
  print,min_i
  print,min_j


  result = CREATE_STRUCT( b,'null_n',n    ,  $
                            'null_i',min_i,  $
                            'null_j',min_j)
;                            'psi_1st_xpoint',b.psi[null_i[1],null_j[1],
;                            'psi_2nd_xpoint',b.psi[null_i[2],null_j[2])

  IF (KEYWORD_SET(debug)) THEN BEGIN
    CONTOUR, b_pol, b_x, b_y, NLEVELS=10, c_labels=[1,1,1,1,1,1,1,1,1,1], color=Truecolor('Black')
    OPLOT, [xpoint_zone[0],xpoint_zone[0],xpoint_zone[1],xpoint_zone[1],xpoint_zone[0]],  $
           [xpoint_zone[2],xpoint_zone[3],xpoint_zone[3],xpoint_zone[2],xpoint_zone[2]],color=TrueColor('Red')
    OPLOT, [b_x[min_i[0]]], [b_y[min_j[0]]], PSYM=6, color=TrueColor('Red')      ; o-point
    IF (n GE 2) THEN  $
      OPLOT, [b_x[min_i[1]]], [b_y[min_j[1]]], PSYM=6, color=TrueColor('Green') ; primary   x-point
    IF (n GE 3) THEN  $
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

  IF (NOT KEYWORD_SET(user_xrange)) THEN user_xrange = [ 3.0,9.0]  ; lame
  IF (NOT KEYWORD_SET(user_yrange)) THEN user_yrange = [-6.0,6.0]

  psi   = b.psi 
  psi_x = b.x
  psi_y = b.y

  FOR j = 1, b.null_n-1 DO BEGIN

    psi_xpt =  b.psi[b.null_i[j],b.null_j[j]]
    psi_step  = MAX(psi) * 0.0001D

    IF (j EQ 2) THEN BEGIN
      psi_1st_xpt = result.psi_1st_xpoint
    ENDIF

    FOR i = 0, 3 DO BEGIN

      psi_start   = psi_xpt + 20.0D * psi_step
      psi_end     = psi_xpt - 20.0D * psi_step
      psi_save    = -1.0D0
      length_xpt  =  0.0D0
      length_last =  0.0D0
      delta_max   =  0.0D0
      delta_last  =  0.0D0

      FOR psi_xpt = psi_start, psi_end, -psi_step DO BEGIN

        IF (j EQ 2) THEN  $   ; Don't get tangled up with the primary x-point if a close-to-connected double-null equilibrium
          IF (psi_xpt GT psi_1st_xpt - psi_step) THEN CONTINUE

        ctr = grid_ExtractContour(psi, psi_x, psi_y, psi_xpt)

        ibrk = WHERE(ctr.dist GT 0.10)
        nseg = N_ELEMENTS(ibrk) 
        ibrk = [-1,ibrk,ctr.n-1]
    
        length_max = 0.0D
        FOR iseg = 0, nseg DO BEGIN
          IF ((ibrk[iseg+1]-ibrk[iseg]) LT 5) THEN CONTINUE
          x = ctr.x[ibrk[iseg]+1:ibrk[iseg+1]]
          y = ctr.y[ibrk[iseg]+1:ibrk[iseg+1]]

          IF (KEYWORD_SET(debug)) THEN BEGIN
            PLOT,x,y,color=Truecolor('Green'),  $
                 XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,PSYM=3
          ENDIF

          length = grid_CalcLength(x,y)
          IF (length GT length_max) THEN length_max = length
        ENDFOR
     
        IF (length_last GT 0.0D) THEN BEGIN 
          delta = length_max - length_last 
          IF (delta GT delta_max) THEN BEGIN
            delta_max = delta
            psi_save  = psi_xpt
          ENDIF
        ENDIF
    
        length_last = length_max
    
;        IF (KEYWORD_SET(debug)) THEN print,'xpt',j,psi_xpt,length_max,psi_save
      ENDFOR
      psi_xpt  = psi_save
      psi_step = psi_step * 0.1D
    ENDFOR
    
    CASE j OF
      1: result.psi_1st_xpoint = psi_xpt
      2: result.psi_2nd_xpoint = psi_xpt-psi_step
;      2: result.psi_2nd_xpoint = psi_xpt-psi_step

;      1: result = CREATE_STRUCT(result,'psi_1st_xpoint',psi_xpt)
;      2: result = CREATE_STRUCT(result,'psi_2nd_xpoint',psi_xpt-psi_step)
      ELSE: stop
    ENDCASE

  ENDFOR 

;        ctr = grid_ExtractContour(psi, psi_x, psi_y, psi_xpt)
;            PLOT,ctr.x,ctr.y,color=Truecolor('Green'),  $
;                 XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,PSYM=3

;stop

  RETURN,result

END
;
; ======================================================================
;
FUNCTION grid_AddContour, b, wall, scan_params, contour_array, mode,  $
                          psi_val=psi_val, focus_x=focus_x, focus_y=focus_y, no_region=no_region,  $
                          debug=debug, xrange=xrange, yrange=yrange, rage=rage
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  IF (NOT KEYWORD_SET(xrange)) THEN xrange = [ 3.0,9.0]  ; lame
  IF (NOT KEYWORD_SET(yrange)) THEN yrange = [-6.0,6.0]

  IF (NOT KEYWORD_SET(mode)) THEN mode = 0

  result = -1

;  geometry    = scan_params.geometry
  geometry    = state.geometry  ; scan_params.geometry
  process_2nd = scan_params.process_2nd

  tags      = STRUPCASE(TAG_NAMES(contour_array))
  contour_n = N_ELEMENTS(tags)

  IF (state.geometry EQ param.LIMITED) THEN BEGIN
    ctr = grid_ExtractStructure(contour_array,tags[0])      
    contact_x = ctr.tangent_p1[0]     
    contact_y = ctr.tangent_p1[1]
  ENDIF

  psi   = b.psi 
  psi_x = b.x
  psi_y = b.y

; need to fix things for limiter grids, i.e. put the origin of the contour near the LCFS tangency point

  SWITCH mode OF 
    -1: BEGIN  ; Special case: secondary separatrix...
      focus_x = b.x[b.null_i[2]]
      focus_y = b.y[b.null_j[2]]
      BREAK
      END
     1: 
     2:
     3:  BEGIN
      IF (N_ELEMENTS(focus_x) EQ 0 OR N_ELEMENTS(psi_val) NE 0 OR  $
          N_ELEMENTS(focus_y) EQ 0) THEN BEGIN
        PRINT, 'ERROR grid_AddContour: Invalid call for MODE=1,2,3'      
        STOP
      ENDIF
 
      n = N_ELEMENTS(psi_x)
      i = WHERE(focus_x GE psi_x[0:n-2] AND focus_x LT psi_x[1:n-1], count)
      IF (count EQ 0) THEN stop    
      n = N_ELEMENTS(psi_y)
      j = WHERE(focus_y GE psi_y[0:n-2] AND focus_y LT psi_y[1:n-1], count)
      IF (count EQ 0) THEN stop

      IF (KEYWORD_SET(debug)) THEN BEGIN
        print,focus_y
        print,psi_y[j],psi_y[j+1]
      ENDIF

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
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, /NOERASE, color=Truecolor('Black')
  ENDIF

  FOR iseg = 0, nseg DO BEGIN

    IF (KEYWORD_SET(debug)) THEN  $
      print, 'iseg',iseg,nseg

    IF ((ibrk[iseg+1]-ibrk[iseg]) LT 5) THEN CONTINUE
  
    x = ctr.x[ibrk[iseg]+1:ibrk[iseg+1]]
    y = ctr.y[ibrk[iseg]+1:ibrk[iseg+1]]

    ; Special check for the case where the secondary PFZ wasn't properly defined,
    ; almost certainly when the 2nd x-point is very close to the vessel wall:
    IF (mode EQ -1 AND scan_params.failure_2nd_pfz EQ 1) THEN BEGIN
      mean_y = MEAN(y)
      IF (( (geometry EQ param.LSN OR geometry EQ param.LDN) AND mean_y GT b.y[b.null_j[2]]) OR  $
          ( (geometry EQ param.USN OR geometry EQ param.UDN) AND mean_y LT b.y[b.null_j[2]])) THEN CONTINUE
;      IF ((geometry EQ param.LOWER_NULL AND mean_y GT b.y[b.null_j[2]]) OR  $
;          (geometry EQ param.UPPER_NULL AND mean_y LT b.y[b.null_j[2]])) THEN CONTINUE
    ENDIF

    IF (state.geometry EQ param.LIMITED) THEN BEGIN
      inside = grid_PointInPolygon(x,y,wall.x,wall.y)
      i = WHERE(inside EQ 0, count)
      proximity = SQRT ( (x[i]-contact_x)^2 + (y[i]-contact_y)^2 )
      dummy = MIN(proximity,imin)
      j = i[imin]
      hold_x = x
      hold_y = y
      n = N_ELEMENTS(x)
      IF (j GT 0) THEN BEGIN
        x = [ hold_x[j:n-1] , hold_x[0:j-1], hold_x[j] ]
        y = [ hold_y[j:n-1] , hold_y[0:j-1], hold_y[j] ]
      ENDIF
    ENDIF

    ; Check if contour runs close to the focus point:
    proximity = SQRT ( (x-focus_x)^2 + (y-focus_y)^2 )

    IF (KEYWORD_SET(rage)) THEN BEGIN
      PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
      OPLOT,wall.x,wall.y,color=Truecolor('Black')
      OPLOT,x,y,color=Truecolor('Purple')
      OPLOT,[x[0]],[y[0]],color=Truecolor('Darkgrey'),PSYM=6
      print, 'min proximity', min(proximity)
    ENDIF

    ; Contour must approach within 5 cm of the focus point:
    IF (MIN(proximity,min_i) GT 0.05D) THEN CONTINUE  ; parameter

    ; Refine the grid near this point, for cases where sharp corners are a problem:
    IF (mode EQ 1 OR mode EQ 2 OR mode EQ 3) THEN BEGIN
      count = 0
      WHILE (count LT 2 OR (count LT 4 AND MIN(proximity) GT 0.001D)) DO BEGIN
        grid_RefineContour, x, y, min_i ; ,  $
;                            debug=debug, xrange=xrange, yrange=yrange
        proximity = SQRT ( (x-focus_x)^2 + (y-focus_y)^2 )
        dummy = MIN(proximity,min_i)
        count++
      ENDWHILE
    ENDIF

    IF (KEYWORD_SET(rage)) THEN BEGIN
      print, 'min proximity ---', min(proximity)
      print, 'mode             ', mode
    ENDIF

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
        PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
        OPLOT,wall.x,wall.y,color=Truecolor('Black')
        OPLOT,x,y,color=Truecolor('Magenta')
        OPLOT,x,y,color=Truecolor('Darkgrey'), PSYM=3
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
            PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
            OPLOT,wall.x,wall.y,color=Truecolor('Black')
            OPLOT,[x[min_i]],[y[min_i]],color=Truecolor('Magenta'), PSYM=6
            OPLOT,[x[i]],[y[i]],color=Truecolor('Orange'), PSYM=6
            OPLOT,x,y,color=Truecolor('Magenta')
            OPLOT,x,y,color=Truecolor('Darkgrey'), PSYM=3
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
            PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
            OPLOT,wall.x,wall.y,color=Truecolor('Black')
            OPLOT,[x[min_i]],[y[min_i]],color=Truecolor('Magenta'), PSYM=6
            OPLOT,[x[i]],[y[i]],color=Truecolor('Orange'), PSYM=6
            OPLOT,x,y,color=Truecolor('Magenta')
            OPLOT,x,y,color=Truecolor('Darkgrey'), PSYM=6
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
;      OPLOT,wall.x,wall.y,color=Truecolor('Black')
      IF (NOT KEYWORD_SET(rage)) THEN OPLOT,x,y,color=Truecolor('Darkgrey')
      OPLOT,[x[0]],[y[0]],color=Truecolor('Lightgreen'), PSYM=6
      OPLOT,[x[N_ELEMENTS(x)-1]],[y[N_ELEMENTS(y)-1]],color=Truecolor('Green'), PSYM=6
      OPLOT,[focus_x],[focus_y],color=Truecolor('Orange'), PSYM=6
    ENDIF

    ; Decide the region:
    IF (N_ELEMENTS(no_region) EQ 0) THEN BEGIN
      mean_x = MEAN(x)
      mean_y = MEAN(y)
      region = param.SOL
      IF (N_ELEMENTS(b.null_i) GT 2 AND psi_val LE b.psi_2nd_xpoint) THEN  $
;        IF (mean_x GT b.x[b.null_j[1]]) THEN  $  ; bug, SL 18/10/2011
        IF (mean_x GT b.x[b.null_i[2]]) THEN  $
          region = param.SOL_LFS ELSE region = param.SOL_HFS
      IF (psi_val GT b.psi_1st_xpoint) THEN BEGIN
        IF (( (geometry EQ param.LSN OR geometry EQ param.LDN) AND mean_y GT b.y[b.null_j[1]]) OR  $
            ( (geometry EQ param.USN OR geometry EQ param.UDN) AND mean_y LT b.y[b.null_j[1]])) THEN  $
;        IF ((geometry EQ param.LOWER_NULL AND mean_y GT b.y[b.null_j[1]]) OR  $
;            (geometry EQ param.UPPER_NULL AND mean_y LT b.y[b.null_j[1]])) THEN  $
          region = param.CORE ELSE region = param.PFZ
      ENDIF
      ; Secondary x-point is inside the vacuum vessel:
      IF (process_2nd GE 0) THEN BEGIN               ; *** THIS CHECK IS NOT PERFECT, SINCE SMALL NEAR-WALL RINGS IN
        IF (psi_val GT b.psi_2nd_xpoint) THEN BEGIN  ;  THE param.SOL COULD QUALIFY FOR param.PFZ ***
          IF (( (geometry EQ param.LSN OR geometry EQ param.LDN) AND mean_y GT b.y[b.null_j[2]]) OR  $
              ( (geometry EQ param.USN OR geometry EQ param.UDN) AND mean_y LT b.y[b.null_j[2]])) THEN  $
;          IF ((geometry EQ param.LOWER_NULL AND mean_y GT b.y[b.null_j[2]]) OR  $
;              (geometry EQ param.UPPER_NULL AND mean_y LT b.y[b.null_j[2]])) THEN  $
            region = param.PFZ_SECONDARY
        ENDIF        
      ENDIF
    ENDIF ELSE region = -1
    IF (KEYWORD_SET(debug)) THEN print,'region',region

    tangent_i = -1
    p1 = [0.0D,0.0D]
    p2 = p1

    ; Define radial boundary associated with the symmetry point, for
    ; the case where the x-point is inside the vessel but no secondary
    ; param.PFZ was generated:
    IF (mode EQ -1) THEN BEGIN
      IF (scan_params.failure_2nd_pfz EQ 1) THEN BEGIN
        tangent_i = min_i
        p1 = [x[tangent_i],y[tangent_i]]
;print, 'mark 1'
        p2 = grid_GetOrthogonal(x,y,p1,2,1.0D,span=10)
;print, 'data', contour_n+1
;print,p1
;print,p2
        PLOT,[p1[0],p2[0]],[p1[1],p2[1]],color=Truecolor('Lightgreen'),  $
             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
        separatrix = 2
      ENDIF ELSE BEGIN
        mean_x = MEAN(x)
        IF (( (geometry EQ param.LSN OR geometry EQ param.LDN) AND mean_x LT b.x[b.null_i[1]]) OR  $
            ( (geometry EQ param.USN OR geometry EQ param.UDN) AND mean_x GT b.x[b.null_i[1]])) THEN BEGIN
;        IF ((geometry EQ param.LOWER_NULL AND mean_x LT b.x[b.null_i[1]]) OR  $
;            (geometry EQ param.UPPER_NULL AND mean_x GT b.x[b.null_i[1]])) THEN BEGIN
;          IF (mean_x GT b.x[b.null_j[1]]) = param.SOL_LFS ELSE region = param.SOL_HFS
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
;          OPLOT,x2,y2,color=Truecolor('Black')
          ; Find an approximation to the normal vector:
          p1 = [x2[min_i4],y2[min_i4]]
;print, 'mark 2'
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

;  tags      = STRUPCASE(TAG_NAMES(contour_array))
;  contour_n = N_ELEMENTS(tags)
;  IF (contour_n EQ 17) THEN begin
;   oplot,x,y,psym=6
;   stop
;  ENDIF

  RETURN, result

END
;
; ======================================================================
;
FUNCTION grid_InstallXPoints, b, scan_params, contour_array,  $
                              debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  tags      = STRUPCASE(TAG_NAMES(contour_array))
  contour_n = N_ELEMENTS(tags)      

  IF (KEYWORD_SET(debug)) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
    FOR i = 1, contour_n DO BEGIN
      ctr = grid_ExtractStructure(contour_array,tags[i-1])      
      OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
      OPLOT,ctr.x,ctr.y,color=Truecolor('Darkgrey'),PSYM=3
    ENDFOR
  ENDIF


  save_j = 0

  count_2nd = 0

  FOR ictr = 1, contour_n DO BEGIN    

    tag = 'contour' + STRING(ictr,FORMAT='(I0)')
    ctr = grid_ExtractStructure(contour_array,tag)  

    status = 0
    IF (ctr.separatrix EQ 1) THEN BEGIN
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

    print,'-----------------------------------------'
    print,'before',ctr.separatrix,ctr.tangent_i
    x = ctr.x
    y = ctr.y      

    print,'status',status,scan_params.failure_2nd_pfz

    FOR i = N_ELEMENTS(y)-2, 0, -1 DO BEGIN 
      j = 0
      IF (status EQ 1) THEN BEGIN
        IF ( (y[i] LE focus_y AND y[i+1] GT focus_y)  OR    $
             (y[i] GE focus_y AND y[i+1] LT focus_y) ) THEN BEGIN
          j = i
          print,'found 1',status,j,N_ELEMENTS(x),focus_x,focus_y
          print,'       ',x[0              ],y[0              ]
          print,'       ',x[N_ELEMENTS(x)-1],y[N_ELEMENTS(y)-1]
        ENDIF
      ENDIF ELSE BEGIN
        proximity = SQRT ( (x-focus_x)^2 + (y-focus_y)^2 )
        IF ( (scan_params.failure_2nd_pfz EQ 0 AND           $
              ((y[i] LE focus_y AND y[i+1] GT focus_y) OR    $
               (y[i] GE focus_y AND y[i+1] LT focus_y))) OR  $
             (scan_params.failure_2nd_pfz EQ 1 AND           $ 
              proximity[i] LT 0.05D            AND           $   ; parameter
              ((x[i] LE focus_x AND x[i+1] GT focus_x) OR    $
               (x[i] GE focus_x AND x[i+1] LT focus_x))) ) THEN BEGIN
          print, 'found 2',status
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

;        if (status eq 2) then begin
;          oplot,[x[save_j]],[y[save_j]],psym=6,color=truecolor('Black')
;          stop
;        endif

      ENDIF

    ENDFOR


    IF (KEYWORD_SET(debug)) THEN BEGIN
      OPLOT,x,y,color=Truecolor('Lightgreen'),PSYM=5+status
    ENDIF  

    ; Setup an outward facing vector for the secondary x-point so that the 
    ; inner param.SOL can be divided into inner and outer regions:
    IF (count_2nd EQ 1) THEN BEGIN
      save_x = x
      save_y = y
    ENDIF
    IF (count_2nd EQ 2) THEN BEGIN
      i = WHERE(x      EQ focus_x AND y      EQ focus_y, count1)
      j = WHERE(save_x EQ focus_x AND save_y EQ focus_y, count2) 
      IF (count1 EQ 0 OR count2 EQ 0) THEN BEGIN
        PRINT,'ERROR grid_InstallXPoints: Secondary x-point not found'
        PRINT,'  COUNT1,2 = ',count1,count2
        IF (count1 EQ 0) THEN BEGIN
          OPLOT,x,y,color=Truecolor('Orange'),PSYM=6
        ENDIF
        STOP
      ENDIF
      IF (ctr.separatrix EQ 3) THEN BEGIN

print, 'j',j ,N_ELEMENTS(x)
        x2 = [save_x[0:j],x     [i+1:N_ELEMENTS(x     )-1]]
        y2 = [save_y[0:j],y     [i+1:N_ELEMENTS(y     )-1]]
      ENDIF ELSE BEGIN 
       help,x
       help,y
       help,i
       help,j
       help,save_x
       help,save_y

        x2 = [x     [0:i],save_x[j+1:N_ELEMENTS(save_x)-1]]
        y2 = [y     [0:i],save_y[j+1:N_ELEMENTS(save_y)-1]]
      ENDELSE
      ; Find an approximation to the 'normal' vector:
      p1 = [focus_x,focus_y]
print, 'mark 3'
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

oplot,ctr.x,ctr.y




    ctr = CREATE_STRUCT(ctr,'null_x',focus_x,'null_y',focus_y)

    ; *** should remove this I think ***
    IF (status EQ 2 AND scan_params.failure_2nd_pfz EQ 1) THEN ctr.tangent_i = save_j

    IF (ctr.separatrix EQ 2) THEN BEGIN
     ctr.tangent_i  = save_j
     ctr.tangent_p1 = [focus_x,focus_y] 
    ENDIF

;    print,'after',ctr.separatrix,ctr.tangent_i

    contour_array = grid_UpdateStructure(contour_array,tag,ctr)

  ENDFOR

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
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  IF (NOT KEYWORD_SET(xrange)) THEN xrange = [ 3.0,9.0]  ; lame
  IF (NOT KEYWORD_SET(yrange)) THEN yrange = [-6.0,6.0]

  result = contour_array

  tags = STRUPCASE(TAG_NAMES(contour_array))
  nctr = N_ELEMENTS(tags)

  IF (KEYWORD_SET(debug)) THEN BEGIN

    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
    FOR i = 1, nctr DO BEGIN
      ctr = grid_ExtractStructure(contour_array,tags[i-1])      
      OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
      OPLOT,ctr.x,ctr.y,color=Truecolor('Darkgrey'),PSYM=3
      npts = N_ELEMENTS(ctr.x)
      XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i),2),color=Truecolor('Darkgrey')
    ENDFOR

    PRINT, ''
    PRINT, '----------------------------------------------------------------------'
    PRINT, ''
    PRINT, 'TRIMMING CONTOURS'
    PRINT, ''
    PRINT, '----------------------------------------------------------------------'
    PRINT, ''

  ENDIF


;***  perhaps need to re-specify wall_pt from wall.x/y here, for a fresh start each time this
;     routine is called

  wall_pt1 = wall.pt1
  wall_pt2 = wall.pt2

  FOR ictr = 0, nctr-1 DO BEGIN

    ; Don't trip the LCFS for a limiter grid, since it is already tidy: 
    IF (state.geometry EQ param.LIMITED AND ictr EQ 0) THEN CONTINUE

    refine_i = -1
    refine_n =  0

    cont = 1
    WHILE (cont) DO BEGIN
      cont = 0

      print, '--------------------'
      print, 'ictr',ictr+1
      
      tag = tags[ictr]
      ctr = grid_ExtractStructure(contour_array,tag)      

      IF (refine_i EQ -1) THEN BEGIN      
        x = ctr.x
        y = ctr.y
      ENDIF ELSE BEGIN
        ; Refine the contour if this is the second pass through:
        grid_RefineContour, x, y, refine_i
        IF (ctr.tangent_i GT -1) THEN BEGIN
          ctr.tangent_i = WHERE( x EQ ctr.tangent_p1[0] AND y EQ ctr.tangent_p1[1], count) 
          IF (count EQ 0) THEN BEGIN
            PRINT,'ERROR grid_TrimContours: Oh for goodness sake'
            PRINT,'  X = ',x
            PRINT,'  Y = ',y
            PRINT,' P1 = ',ctr.tangent_p1
            STOP
          ENDIF
        ENDIF
      ENDELSE

      tangent_i = ctr.tangent_i
      
      OPLOT,x,y,color=Truecolor('Red'),psym=3
;     
;     -------------------------------------------------------------------- 
;     
;     
      FOR j = 0, 1 DO BEGIN
        FOR i = 0, N_ELEMENTS(x)-3 DO BEGIN
          result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                     wall_pt1, wall_pt2, 1, status=status)
          IF (status) THEN BREAK
        ENDFOR
        IF (NOT status AND j EQ 0) THEN BEGIN
          print, 'TRYING TO SAVE THE DAY! 1'
          n = N_ELEMENTS(x) - 1
          x[0] = 2.0D * (x[0] - x[1]) + x[1]
          y[0] = 2.0D * (y[0] - y[1]) + y[1]
        ENDIF ELSE BREAK
      ENDFOR
      IF (NOT status) THEN BEGIN
        PRINT,'ERROR grid_TrimContours: No wall intersection found 1'
        STOP
      ENDIF ELSE BEGIN
        print,'shit 1',i,N_ELEMENTS(x)-1,ctr.tangent_p1
;        print,'      ',x[0:i+1]
;        print,'      ',y[0:i+1]
        x = x[i:N_ELEMENTS(x)-1] ; This clean up appears to be required, but I'm not sure why...
        y = y[i:N_ELEMENTS(y)-1]
      
        IF (tangent_i NE -1) THEN tangent_i = tangent_i - i
      ENDELSE
;
;     ------------------------------------------------------------------ 
;   
;
      FOR j = 0, 1 DO BEGIN      

        FOR i = N_ELEMENTS(x)-2, 1, -1 DO BEGIN
          result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                     wall_pt1, wall_pt2, 1, status=status)
          IF (status) THEN BREAK
        ENDFOR

        IF (NOT status AND j EQ 0) THEN BEGIN
          print, 'TRYING TO SAVE THE DAY! 2'
          n = N_ELEMENTS(x) - 1
          x[n] = 2.0D * (x[n] - x[n-1]) + x[n-1]
          y[n] = 2.0D * (y[n] - y[n-1]) + y[n-1]
        ENDIF ELSE BREAK

      ENDFOR

      IF (NOT status) THEN BEGIN
        PRINT,'ERROR grid_TrimContours: No wall intersection found 2'
      
        xspan = 0.02D ; 0.50D ; 0.01D
        yspan = 0.5D ; 0.50D ; 0.05D
        i = N_ELEMENTS(x)-1
        xrange = [x[i] - xspan, x[i] + xspan]
        yrange = [y[i] - yspan, y[i] + yspan]
        PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
        OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange')
        OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange'), PSYM=7
        OPLOT,x,y,color=Truecolor('Magenta'),PSYM=6
      
        STOP
      ENDIF ELSE BEGIN
        print,'shit 2',i,N_ELEMENTS(x)-1,ctr.tangent_p1
;        print,'      ',x[i-1:N_ELEMENTS(x)-1]
;        print,'      ',y[i-1:N_ELEMENTS(y)-1]
        x = x[0:i+1]
        y = y[0:i+1]
      ENDELSE
;
;     ------------------------------------------------------------------    
;     CHECK THAT THERE ARE ENOUGH POINTS ON THE CONTOUR (>5)
;
      IF (N_ELEMENTS(x) LE 5) THEN BEGIN

        IF (refine_n EQ 10) THEN BEGIN
          print,'here in the shit, again',N_ELEMENTS(x)
          stop
        ENDIF

        print,'ictr',ictr
        print,N_ELEMENTS(x)
        print, 'nforcing local refinement because not enough points',N_ELEMENTS(x)

        n = N_ELEMENTS(x)

        refine_i = WHERE( ctr.x EQ x[n/2] AND ctr.y EQ y[n/2] , count)
        IF (count EQ 0) THEN BEGIN
          refine_i = n/2
          print,'refine_i guess',refine_i
;          PRINT,'ERROR grid_TrimContours: Another great mystery presents itself'
;          STOP
        ENDIF

        refine_n++

        cont = 1
      ENDIF

    ENDWHILE

    OPLOT,x,y,color=Truecolor('Lightgreen'),PSYM=3

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
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  grid_TrimContours, b, contour_array, wall, kill=kill,  $
                     debug=debug, xrange=xrange, yrange=yrange

;  IF (KEYWORD_SET(kill)) THEN  $
;    grid_Debug,b,contour_array,wall,xrange,yrange

  IF (debug) THEN BEGIN
    PRINT, ' '
    PRINT, '----------------------------------------------------------------------'
    PRINT, ''
    PRINT, 'UPDATE WALL'
    PRINT, ''
    PRINT, '----------------------------------------------------------------------'
    PRINT, ' '
  ENDIF

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

    print, '----------------------------------------------------------------------'
    print, 'TARGETS ',ictr,' ',tags[ictr-1]
    print, '----------------------------------------------------------------------'

    tag = tags[ictr-1]
    ctr = grid_ExtractStructure(contour_array,tag)      

    print,'ctr.psi=',ctr.psi

    IF (state.geometry EQ param.LIMITED AND ctr.separatrix EQ 1) THEN CONTINUE

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
        PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
        OPLOT,wall.x,wall.y,color=Truecolor('Black')
        OPLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Orange')
        OPLOT,x,y,color=Truecolor('Magenta')
        OPLOT,x,y,color=Truecolor('Darkgrey'), PSYM=6
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
          PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
          OPLOT,wall.x,wall.y,color=Truecolor('Black')
          OPLOT,wall.pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange')
          OPLOT,wall.pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange'),PSYM=7
          FOR i = 0, N_ELEMENTS(wall_pt1[0,*])-1 DO  $
            XYOUTS,0.5*(wall_pt1[0,i]+wall_pt2[0,i]),  $
                   0.5*(wall_pt1[1,i]+wall_pt2[1,i]),  $
                   STRTRIM(STRING(i),2),color=Truecolor('Darkgrey')
          OPLOT,x,y,color=Truecolor('Magenta')
          OPLOT,x,y,color=Truecolor('Darkgrey'), PSYM=6

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
                ;OPLOT,wall.x,wall.y,color=Truecolor('Black')
                ;OPLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Orange')
                ;OPLOT,x,y,color=Truecolor('Magenta')
                ;OPLOT,x,y,color=Truecolor('Darkgrey'), PSYM=6
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
        PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
        OPLOT,wall.x,wall.y,color=Truecolor('Black')
        OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange')
        OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange'), PSYM=7
        OPLOT,x,y,color=Truecolor('Magenta')
        OPLOT,x,y,color=Truecolor('Darkgrey'), PSYM=6
        STOP
      ENDIF ELSE BEGIN
        ;print,'damn 2' ; ,result.i,N_ELEMENTS(wall_pti)-1
        IF (N_ELEMENTS(result.i) GT 1) THEN BEGIN
          PRINT,'ERROR grid_UpdateWall: Multiple intersections detected (2)'
          xspan = 0.50D ; 0.50D ; 0.01D
          yspan = 0.50D ; 0.50D ; 0.05D
          xrange = [x[n-2] - xspan, x[n-2] + xspan]
          yrange = [y[n-2] - yspan, y[n-2] + yspan]
          PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
          OPLOT,wall.x,wall.y,color=Truecolor('Black')
          OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange')
          OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange'), PSYM=7
          OPLOT,x,y,color=Truecolor('Magenta')
          OPLOT,x,y,color=Truecolor('Darkgrey'), PSYM=6
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
                ;OPLOT,wall.x,wall.y,color=Truecolor('Black')
                ;OPLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Orange')
                ;OPLOT,x,y,color=Truecolor('Magenta')
                ;OPLOT,x,y,color=Truecolor('Darkgrey'), PSYM=6
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

;  FOR i = 0, N_ELEMENTS(wall_pt1[0,*])-1 DO BEGIN
;     PRINT,'wall i',i,wall_pt1[0:1,i],wall_pt2[0:1,i],wall_ptc[i],FORMAT='(A,I6,2(2X,2F15.9),I6)'
;  ENDFOR
;
; ----------------------------------------------------------------------
; LOOK AFTER TANGENCY POINTS
; 
  FOR ictr = 1, nctr DO BEGIN

    print, '----------------------------------------------------------------------'
    print, 'TANGENCY ',ictr,' ',tags[ictr-1]
    print, '----------------------------------------------------------------------'

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
                                                         ; parameter

    result = grid_Intersection([p3[0],p3[1]], [p4[0],p4[1]],  $
                               wall_pt1, wall_pt2, 1, status=status)
    IF (NOT status) THEN BEGIN
      PRINT,'ERROR grid_UpdateWall: No tangency point wall intersection found'
;      xrange = [p3[0] - 0.20D, p3[0] + 0.20D]
;      yrange = [p3[1] - 0.20D, p3[1] + 0.20D]
      xrange = [MIN(x),MAX(x)]
      yrange = [MIN(y),MAX(y)]
      print,'xrange=',xrange
      print,'yrange=',yrange
      print,'p1=',p1
      print,'p2=',p2
      print,'p3=',p3
      PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
      OPLOT,wall.x,wall.y,color=Truecolor('Black')
      OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange')
      OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange'), PSYM=7
      OPLOT,x,y,color=Truecolor('Magenta')
      OPLOT,x,y,color=Truecolor('Darkgrey'), PSYM=6
      OPLOT,[p3[0],p4[0]], [p3[1],p4[1]], color=Truecolor('Green')
      STOP
    ENDIF ELSE BEGIN
      ;print,'hot damn 1' ; ,result.i,N_ELEMENTS(wall_pti)-1
      IF (N_ELEMENTS(result.i) GT 1) THEN BEGIN
        PRINT,'ERROR grid_UpdateWall: Multiple tangency intersections detected'
        PRINT,'  I=',result.i
        PRINT,'  S=',result.s
        PRINT,'  T=',result.t
 ;       PRINT,'  1x=',wall_pt1(0,result.i[0]:result.i[1]),FORMAT='(A,3F15.9)'
 ;       PRINT,'  1y=',wall_pt1(1,result.i[0]:result.i[1]),FORMAT='(A,3F15.9)'
 ;       PRINT,'  2x=',wall_pt2(0,result.i[0]:result.i[1]),FORMAT='(A,3F15.9)'
 ;       PRINT,'  2y=',wall_pt2(1,result.i[0]:result.i[1]),FORMAT='(A,3F15.9)'
        xrange = [result.x[0] - 0.005D, result.x[0] + 0.005D]
        yrange = [result.y[0] - 0.005D, result.y[0] + 0.005D]
        PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
        OPLOT,wall.x,wall.y,color=Truecolor('Black')
        OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange')
        OPLOT,wall_pt1[0,*],wall_pt1[1,*],color=Truecolor('Orange'), PSYM=7
        OPLOT,x,y,color=Truecolor('Magenta')
        OPLOT,x,y,color=Truecolor('Darkgrey'), PSYM=6
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

print,'tangent_i',tangent_i
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


  n=N_ELEMENTS(wall_pti)
  FOR j = 0, n-1 DO BEGIN
    print,j,wall.pti[j],wall.ptt[j],wall.ptc[j],wall.pt1[0,j],wall.pt1[1,j],wall.pt2[0,j],wall.pt2[1,j]
  ENDFOR

  RETURN
END
;
; ======================================================================
;
PRO grid_ProcessWall, b, wall, scan_params, contour_array, $
                      debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  IF (debug) THEN BEGIN
    PRINT, ' '
    PRINT, '----------------------------------------------------------------------'
    PRINT, ''
    PRINT, 'PROCESS WALL'
    PRINT, ''
    PRINT, '----------------------------------------------------------------------'
    PRINT, ' '
  ENDIF

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
      IF (state.geometry EQ param.LIMITED) THEN  $
        iwall = WHERE(wall_ptc EQ 1 AND wall_ptt EQ 3)  ELSE  $  ; Start at the LCFS tangency point for limiter grids
        iwall = WHERE(wall_ptc EQ 1 AND wall_ptt EQ 1)           ; Start at the primary separatrix segment
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




  PLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Black'),  $
       XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
  PLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Orange'), PSYM=7,  $
       XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

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
; ----------------------------------------------------------------------
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

    IF (sa GT 0.01D AND sc GT 0.01D) THEN BEGIN            ; parameter
;    IF (sa GT 0.01D OR sc GT 0.01D) THEN BEGIN            ; parameter
      cos_theta = (sa^2 + sc^2 - sb^2) / (2.0D * sa * sc)
 
      IF ((1.0D0 - ABS(cos_theta)) GT 1.0E-06) THEN  $
        theta = ABS(ACOS(cos_theta) * 180.0D / DOUBLE(!PI)) ELSE  $
        theta = 180.0D    
    ENDIF ELSE theta = 180.0D

    wall_angle[i] = theta

    ; print,i,'yellow',wall_ptt[i],wall_angle[i],cos_theta
  ENDFOR

;
; Add a contour for any sharp angles that are inside the grid:
;
  n = N_ELEMENTS(wall_pti)
  count = 0
  FOR i = 0, n-1 DO BEGIN

    IF (wall_angle[i] GT 150.0D) THEN CONTINUE ; 150.0D) THEN CONTINUE
    IF (wall_ptc  [i] NE -1    ) THEN CONTINUE
    IF (wall_ptt  [i] EQ -1    ) THEN CONTINUE

    count = count + 1

    PLOT,[wall.pt1[0,i]],[wall.pt1[1,i]],color=Truecolor('Blue'),PSYM=1,  $
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
    XYOUTS, wall.pt1[0,i], wall.pt1[1,i], STRTRIM(STRING(count),2), color=Truecolor('Darkgrey')
    print, '----- working the wall', count,' -----'

    focus_x = wall.pt1[0,i]
    focus_y = wall.pt1[1,i]


    status = grid_AddContour(b, wall, scan_params, contour_array, wall_ptt[i], $
                             focus_x=focus_x, focus_y=focus_y,  $
                             debug=debug, xrange=xrange, yrange=yrange, rage=rage)

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
FUNCTION grid_BuildRadialMap, c_array, wall, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  tags  = STRUPCASE(TAG_NAMES(c_array))
  nctr  = N_ELEMENTS(tags)
  nwall = N_ELEMENTS(wall.ptc)

  IF (state.geometry EQ param.LIMITED) THEN  $
    istart = WHERE(wall.ptc EQ 1 AND wall.ptt EQ 3, count)  ELSE  $  ; Start at the LCFS tangency point for limiter grids
    istart = WHERE(wall.ptc EQ 1 AND wall.ptt EQ 1, count)           ; Start at the primary separatrix segment
;  istart = WHERE(wall.ptc EQ 1 AND wall.ptt EQ 1, count)
  IF (count EQ 0) THEN BEGIN
    PRINT, 'ERROR grid_BuildRadialMap: Starting contour not found'
    STOP
  ENDIF ELSE IF (KEYWORD_SET(debug)) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
    OPLOT,wall.x,wall.y,color=Truecolor('Black')
    OPLOT,[wall.pt1[0,istart]],[wall.pt1[1,istart]],color=Truecolor('Hotpink'),PSYM=6
    FOR i = 0, nctr-1 DO BEGIN
      ctr = grid_ExtractStructure(c_array,tags[i])      
      OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
      npts = N_ELEMENTS(ctr.x)
      XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('Darkgrey'),CHARSIZE=2.0
    ENDFOR
  ENDIF

  ; There should be up to 2 contours mapped to any other contour, which
  ; occurs when there was a tagnecy point:
  map_in  = MAKE_ARRAY(nctr,3,/LONG,VALUE=-1L)
  map_out = map_in

  iwall = istart
  cwall = 1
  ictr_last = 1
  itar_last = 1

  PRINT,'---------- START RADIAL MAP ------------'

  WHILE (cwall LE nwall) DO BEGIN
    cwall++
    iwall++
    IF (iwall GT nwall-1) THEN iwall = 0

    print,nwall,iwall,cwall,FORMAT='(3I8)'    

    ictr = wall.ptc[iwall]  ; contour index of wall point
    itar = wall.ptt[iwall]  ; target 

    ; Cycle to the next wall point if this one isn't associated with a contour:
    IF (ictr EQ -1) THEN CONTINUE

    IF (KEYWORD_SET(debug)) THEN  $
      OPLOT,[wall.pt1[0,iwall]],[wall.pt1[1,iwall]],color=Truecolor('Hotpink'),PSYM=1

    print,'OK---',ictr,itar,ictr_last,itar_last,FORMAT='(A8,4I8)'    

    IF ((itar EQ 1 OR itar EQ 3) AND (itar_last EQ 1 OR itar_last EQ 3)) THEN BEGIN  ; NEW

      tag = 'CONTOUR'+STRTRIM(STRING(ictr),2)
      region = grid_GetRegion(c_array,tag)

      IF (itar NE 3 AND itar_last NE 3) THEN BEGIN
        IF ((region GE param.SOL AND region LE param.SOL_HFS) OR itar_last EQ 1) THEN i = 0 ELSE i = 1
;        IF (region EQ param.SOL OR itar_last EQ 1) THEN i = 0 ELSE i = 1
        IF (map_in[ictr-1,i] NE -1) THEN BEGIN
          PRINT, 'ERROR grid_BuildRadialMap: Inward overload (1)'
          FOR i = 0, nctr-1 DO print,i+1,map_in[i,0:1],format='(3I6)' 
          STOP
        ENDIF
        map_in[ictr-1,i] = ictr_last 
        map_in[ictr-1,2] = region
      ENDIF

      IF ((region GE param.SOL AND  $
           region LE param.SOL_HFS AND itar EQ 1 AND itar_last EQ 3) OR  $
          (region GE param.PFZ     AND itar EQ 3 AND itar_last EQ 1)) THEN BEGIN
        IF (region GE param.PFZ) THEN i = 1 ELSE i = 0
;      IF ((region EQ param.SOL AND itar EQ 1 AND itar_last EQ 3) OR  $
;          (region EQ param.PFZ AND itar EQ 3 AND itar_last EQ 1)) THEN BEGIN
;        IF (region EQ param.PFZ) THEN i = 1 ELSE i = 0
        IF (map_in[ictr-1,i] NE -1) THEN BEGIN
          PRINT, 'ERROR grid_BuildRadialMap: Inward overload (2)'
          FOR i = 0, nctr-1 DO print,i+1,map_in[i,0:1],format='(3I6)' 
          STOP
        ENDIF
        map_in[ictr-1,i] = ictr_last 
        map_in[ictr-1,2] = region
      ENDIF

      tag = 'CONTOUR'+STRTRIM(STRING(ictr_last),2)
      region = grid_GetRegion(c_array,tag)

      IF (itar NE 3 AND itar_last NE 3) THEN BEGIN
        IF (region GE param.PFZ OR itar_last EQ 1) THEN i = 0 ELSE i = 1
;        IF (region EQ param.PFZ OR itar_last EQ 1) THEN i = 0 ELSE i = 1
        IF (map_out[ictr_last-1,i] NE -1) THEN BEGIN
          PRINT, 'ERROR grid_BuildRadialMap: Outward overload'
          FOR i = 0, nctr-1 DO print,i+1,map_out[i,0:1],format='(3I6)' 
          STOP
        ENDIF
        map_out[ictr_last-1,i] = ictr 
        map_out[ictr_last-1,2] = region
      ENDIF

      IF ((region GE param.SOL AND  $
           region LE param.SOL_HFS AND itar EQ 1 AND itar_last EQ 3) OR  $
          (region GE param.PFZ     AND itar EQ 3 AND itar_last EQ 1)) THEN BEGIN
        IF (region GE param.PFZ) THEN i = 0 ELSE i = 1
;      IF ((region EQ param.SOL AND itar EQ 1 AND itar_last EQ 3) OR  $
;          (region EQ param.PFZ AND itar EQ 3 AND itar_last EQ 1)) THEN BEGIN
;        IF (region EQ param.PFZ) THEN i = 0 ELSE i = 1
        IF (map_out[ictr_last-1,i] NE -1) THEN BEGIN
          PRINT, 'ERROR grid_BuildRadialMap: Outward overload (2)'
          FOR i = 0, nctr-1 DO print,i+1,map_out[i,0:1],format='(3I6)' 
          STOP
        ENDIF
        map_out[ictr_last-1,i] = ictr 
        map_out[ictr_last-1,2] = region
      ENDIF

    ENDIF

    print,'  MAP',map_in [ictr-1,0:1],FORMAT='(A8,3I8)'    
    print,'     ',map_out[ictr-1,0:1],FORMAT='(A8,3I8)'    

    ictr_last = ictr
    itar_last = itar

  ENDWHILE


  ; Small update to map for secondary separatrix:
  FOR ictr = 1, nctr DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[ictr-1])      
    IF (ctr.separatrix LT 3) THEN CONTINUE
    ; Scan the wall to find the next point that references a ring and
    ; assign that ring to the connection map:
    iwall = WHERE(wall.ptc EQ ictr AND wall.ptt EQ 2, count)
    IF (count LE 0) THEN BEGIN
      PRINT,'ERROR grid_BuildRadialMap: Unable to find separatrix segment'
      PRINT,'  ICTR=',ictr
      PRINT,'  SEP =',ctr.separatrix 
      STOP      
    ENDIF
    cnt = 0
    WHILE (1 EQ 1) DO BEGIN
      cnt++
      IF (cnt GT nwall) THEN BEGIN
        PRINT,'ERROR grid_BuildRadialMap: Unable to find wall reference'
        STOP      
      ENDIF
      iwall = grid_FindWallNeighbour(wall.pt1,wall.pt2,wall.pti,iwall,2)
      IF (wall.ptc[iwall] GT 0 AND wall.ptt[iwall] EQ 2) THEN BEGIN
        map_in[ictr-1,1] = wall.ptc[iwall]
        IF (ctr.separatrix EQ 4) THEN  map_out[wall.ptc[iwall]-1,1] = ictr
        BREAK
      ENDIF
    ENDWHILE
  ENDFOR


  ; Update the connection map information for each contour:
  FOR i = 0, nctr-1 DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[i])      

    dummy = WHERE(STRUPCASE(TAG_NAMES(ctr)) EQ 'MAP_IN',count)
    IF (count EQ 0) THEN BEGIN
      ctr = CREATE_STRUCT(ctr,'map_in',map_in[i,0:1],'map_out',map_out[i,0:1]) 
    ENDIF ELSE BEGIN
      ctr = grid_UpdateStructure(ctr,'map_in' ,map_in [i,0:1])
      ctr = grid_UpdateStructure(ctr,'map_out',map_out[i,0:1])
    ENDELSE
    c_array = grid_UpdateStructure(c_array,tags[i],ctr)
  ENDFOR


  print,'in'  
  FOR i = 0, nctr-1 DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[i])      
    print,i+1,map_in[i,0:2],'  ',ctr.map_in[0:1],format='(4I6,A,2I6)'
  ENDFOR

  print,'out'  
  FOR i = 0, nctr-1 DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[i])      
    print,i+1,map_out[i,0:2],'  ',ctr.map_out[0:1],format='(4I6,A,2I6)'
  ENDFOR

  result = c_array

  RETURN, result

END
;
; ======================================================================
;
FUNCTION grid_RadialRefinement,c_array, wall, b, $
                          debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  IF (debug) THEN BEGIN
    PRINT, ' '
    PRINT, '----------------------------------------------------------------------'
    PRINT, ''
    PRINT, 'RADIAL REFINEMENT'
    PRINT, ''
    PRINT, '----------------------------------------------------------------------'
    PRINT, ' '
  ENDIF

  FOR ipass = 1, 5 DO BEGIN

    print, '============================================================'
    print, '  ipass ',ipass
    print, '============================================================'

    tags  = STRUPCASE(TAG_NAMES(c_array))
    nctr  = N_ELEMENTS(tags)
    nwall = N_ELEMENTS(wall.ptc)

    add_contour = 0

    IF (state.geometry EQ param.LIMITED) THEN BEGIN
      ctr = grid_ExtractStructure(c_array,tags[0])      
      contact_x = ctr.tangent_p1[0]     
      contact_y = ctr.tangent_p1[1]
    ENDIF

    IF (KEYWORD_SET(debug)) THEN BEGIN
      PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black')
      OPLOT,wall.x,wall.y,color=Truecolor('Black')
      OPLOT,[wall.pt1[0,*]],[wall.pt1[1,*]],color=Truecolor('Hotpink'),PSYM=6
      FOR i = 0, nctr-1 DO BEGIN
        ctr = grid_ExtractStructure(c_array,tags[i])      
        OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
        npts = N_ELEMENTS(ctr.x)
        XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('Darkgrey')
      ENDFOR
    ENDIF

    CASE ipass OF 
      1: BEGIN
;       ----------------------------------------------------------------
;       COLLECT THE DESIRED DISTRIBUTION OF RINGS IN THE SOL AND PFRs
;  
        dist = grid_RadialDistribution(c_array, wall, b, param.SOL,  $
                                       debug=debug, xrange=xrange, yrange=yrange)

        psi_min = ABS(dist.psi[dist.izero] - dist.psi[dist.izero-1])
        print,'min a',dist.dist[dist.izero-1:dist.izero+1]
        print,'min a',dist.psi [dist.izero-1:dist.izero+1],psi_min
        END
;       ----------------------------------------------------------------
      ELSE:
    ENDCASE
;
;   --------------------------------------------------------------------
;   FIND WHERE TO START THE GRID WHEN FOLLOWING THE WALL
;
    ; Start with the contour that's just outside the low-index primary
    ; strike-point (the inner target strike-point for lower single null):

    ctr1  = grid_ExtractStructure(c_array,'contour1')  ; separatrix
    IF (ipass EQ 3) THEN itarget = 2 ELSE itarget = 1
    ictr1 = (ctr1.map_out[0])[0]                       ; step out radially by one contour
    iwall1 = WHERE(wall.ptc EQ ictr1 AND wall.ptt EQ itarget, count)
    IF (count NE 1) THEN BEGIN
      PRINT,'ERROR grid_RadialRefinement: Something is wrong with the wall'
      PRINT,'  IWALL1=',iwall1
      STOP 
    ENDIF
    ctr1   = grid_ExtractStructure(c_array,tags[ictr1-1])  
    istart = iwall1
;
;   --------------------------------------------------------------------
;   FOLLOW THE WALL AND BUILD THE GRID RINGS
;
    ictr_new    = 0
    nskip_total = 0

    WHILE (1) DO BEGIN  
;  
;     ------------------------------------------------------------------
;     SCAN BACKWARD ALONG THE WALL TO FIND THE NEIGHBOURING CONTOUR
;
      print,'new contour ----------- ',iwall1,ictr1,'-----------',FORMAT='(A,2I6,A)'
      print,'      istart,ipass=',istart,ipass,FORMAT='(A,2I6)'

      iwall2 = iwall1
      print,'      ptc,ptt=',wall.ptc[iwall1],wall.ptt[iwall1],FORMAT='(A,2I6)'
      WHILE (1) DO BEGIN
        iwall2--
        IF (iwall2 LT 0) THEN iwall2 = N_ELEMENTS(wall.ptc) - 1
        print,'      iwall2',iwall2,wall.ptc[iwall2],wall.ptt[iwall2],FORMAT='(A,3I6)'
        IF (wall.ptc[iwall2] NE -1 AND (wall.ptt[iwall2] EQ itarget OR  $
                                        wall.ptt[iwall2] EQ 3          )) THEN BREAK
      ENDWHILE

      ictr2 = (wall.ptc[iwall2])[0]
      ctr2  = grid_ExtractStructure(c_array,tags[ictr2-1])      
 
      print,'            ----------- ',iwall2,ictr2,'-----------',FORMAT='(A,2I6,6X,A)'
      print,'      ptc,ptt=',wall.ptc[iwall2],wall.ptt[iwall2],FORMAT='(A,2I6)'

      IF (KEYWORD_SET(debug)) THEN BEGIN
        OPLOT,[wall.pt1[0,iwall1],wall.pt1[0,iwall2]],  $
              [wall.pt1[1,iwall1],wall.pt1[1,iwall2]],color=Truecolor('Lightblue'),PSYM=1
      ENDIF

      count  = 0
      status = 1

      SWITCH ipass OF
;  
;       ----------------------------------------------------------------
;       DECIDE IF A NEW CONTOUR IS REQUIRED  
;
        1: BEGIN
          print,'range',ctr1.psi,ctr2.psi,FORMAT='(A,2F15.7)'

          iselect = WHERE(dist.psi GT (ctr1.psi)[0] AND dist.psi LT (ctr2.psi)[0], count )
          iselect2= WHERE(dist.psi GT (ctr2.psi)[0] AND dist.psi LT (ctr1.psi)[0], count2)

          IF (count EQ 0 AND count2 NE 0) THEN BEGIN
            print,'count',count,count2
            PRINT,'DAMN THE COUNT'
            STOP
          ENDIF
          BREAK

          END
;  
;       ----------------------------------------------------------------
;       CHECK IF TARGET SEGMENTS ARE TOO LONG
;
        2:
        3: BEGIN
          oplot,[wall.pt1[0,iwall1]],[wall.pt1[1,iwall1]],psym=6,color=Truecolor('Lightgreen')

          length = SQRT( (wall.pt1[0,iwall1] - wall.pt1[0,iwall2])^2 +  $
                         (wall.pt1[1,iwall1] - wall.pt1[1,iwall2])^2 )

          v1 = wall.pt1[0:1,iwall1]
          v2 = wall.pt2[0:1,iwall1]
          IF (ipass EQ 2) THEN v3 = [ctr1.x[                   3],ctr1.y[                   3]] ELSE  $
                               v3 = [ctr1.x[N_ELEMENTS(ctr1.x)-3],ctr1.y[N_ELEMENTS(ctr1.x)-3]]

          as = SQRT( (v2[0]-v1[0])^2 + (v2[1]-v1[1])^2 )
          bs = SQRT( (v3[0]-v1[0])^2 + (v3[1]-v1[1])^2 )
          cs = SQRT( (v3[0]-v2[0])^2 + (v3[1]-v2[1])^2 )

          angle = ACOS( (as^2 + bs^2 - cs^2) / (2.0D * as * bs) ) * 180.0D / 3.1415D
          IF (angle GT 90.0D) THEN angle = 180.0D - angle

          IF (angle GT option.tar_threshold_angle) THEN  $ 
            threshold = option.tar_threshold_size ELSE  $ 
            threshold = option.tar_threshold_size *  $
                        MIN([5.0D,option.tar_threshold_angle/angle])  

          IF (length GT threshold) THEN BEGIN  

            oplot,[wall.pt1[0,iwall1]],[wall.pt1[1,iwall1]],psym=6,color=Truecolor('Black')
            print,'   >>>> doh!',ipass,length,threshold
            print,'            ',iwall1,iwall2

            count = MAX( [1, LONG(length / threshold)] )
            dist = DINDGEN(count+2) / DOUBLE(count+1)
            count = N_ELEMENTS(dist)

            ref_psi = ctr1.psi + dist[1] * (ctr2.psi - ctr1.psi)  ; short-hand method didn't work for COUNT=4... IDL bug...
            FOR i = 2, count-2 DO  $
              ref_psi = [ref_psi, ctr1.psi + dist[i] * (ctr2.psi - ctr1.psi)]

            count = N_ELEMENTS(ref_psi)

           IF (length GT 0.4D) THEN BEGIN
                    print,'    ******************* HUGE! ********************'
                    oplot,[wall.pt1[0,iwall1],wall.pt1[0,iwall2]],  $
                          [wall.pt1[1,iwall1],wall.pt1[1,iwall2]],thick=2.0,color=Truecolor('Cyan')
;stop
           ENDIF


; *** LOOK HERE MATE! ***
; if (ipass EQ 2) then count = 0
;  if (ipass EQ 3) then count = 0  

          ENDIF
          BREAK
          END
;  
;       ----------------------------------------------------------------
;       CHECK TO SEE IF THE RADIAL WIDTH RINGS IS UNUSUALLY LARGE, WHICH
;       CAN HAPPEN NEAR AN X-POINT
;
        4: BEGIN

          print,'checking ring width',ictr1,ictr2
         
          count = 0

          IF ( (ctr1.region GE param.PFZ AND ctr2.region LT param.PFZ) OR  $
               (ctr1.region LT param.PFZ AND ctr2.region GE param.PFZ) OR  $
               ctr1.region GE param.PFZ OR ctr2.region GT param.PFZ ) THEN BEGIN
            print,'skipping, in different regions'
            BREAK
          ENDIF

          IF ( (ctr1.separatrix NE 0 AND ctr2.separatrix NE 0) ) THEN BEGIN
            print,'skipping, in different regions'
            BREAK
          ENDIF

          x = ctr2.x
          y = ctr2.y

          FOR i = 0, N_ELEMENTS(ctr1.x)-1 DO BEGIN

            focus_x = ctr1.x[i]
            focus_y = ctr1.y[i]

            proximity = SQRT ( (x - focus_x)^2 + (y - focus_y)^2 )

            IF (MIN(proximity) GT 0.1D0) THEN BEGIN ; 1.0D0) THEN BEGIN ; 0.5D0 100.0D0) THEN BEGIN ; 0.10D) THEN BEGIN  ; parameter  ; dev

              oplot,ctr1.x,ctr1.y,color=Truecolor('Lightgreen')
              oplot,ctr2.x,ctr2.y,color=Truecolor('Lightgreen')
              oplot,[focus_x],[focus_y],color=Truecolor('Lightgreen'),PSYM=6

              psi_new = 0.5D * (ctr1.psi + ctr2.psi)

              print,ctr1.psi
              print,ctr2.psi
              print,psi_new

              IF (count EQ 0) THEN BEGIN
                count = 1
                ref_psi = [psi_new]
              ENDIF ELSE BEGIN

                psi_span = ABS(ctr1.psi - ctr2.psi)

                FOR j = 0, N_ELEMENTS(ref_psi)-1 DO  $
                  IF (ABS(psi_new - ref_psi[j]) LT 0.05D * psi_span) THEN BREAK ; parameter

                IF (j EQ N_ELEMENTS(ref_psi)) THEN BEGIN
                  count++
                  ref_psi = ref_psi[ref_psi,psi_new]
                ENDIF ELSE BEGIN

 print,'caught!'

                ENDELSE

 print, 'count',count

              ENDELSE
              
              BREAK
 
            ENDIF 

          ENDFOR


          BREAK
          END 
;  
;       ----------------------------------------------------------------
;       ALSO CHECK THE REFINEMENT LIST DEFINED WHEN VALIDATING THE GRID
;
        5: BEGIN

          IF (state.refinement_rad GT 0) THEN BEGIN

            ref_psi = 0.5D * (state.ref_rad_psi1 + state.ref_rad_psi2)

            i = WHERE(ref_psi GT MIN([ctr1.psi,ctr2.psi]) AND  $
                      ref_psi LT MAX([ctr1.psi,ctr2.psi]), count)

            print, 'checking.....',MIN([ctr1.psi,ctr2.psi]),MAX([ctr1.psi,ctr2.psi])
            print, 'psi',ref_psi
            print, 'count',count
            IF (count GT 0) THEN ref_psi = ref_psi[i]  
          ENDIF
          BREAK
          END 
;
;     ------------------------------------------------------------------
;
      ENDSWITCH
;  
;     ------------------------------------------------------------------
;     ADD THE NEW CONTOURS, IF ANY
;
      IF (count GT 0) THEN BEGIN

        new_contour = 0
        nskip       = 0

        FOR icount = 0, count-1 DO BEGIN

          IF (debug) THEN BEGIN
            PRINT, ' '
            PRINT, ' ----------------------------------------------------------------------'
            PRINT, ' TRYING TO ADD A CONTOUR',ipass,icount,count
            PRINT, ' ----------------------------------------------------------------------'
            PRINT, ' '
          ENDIF

          status = 1

          CASE ipass OF

            1: BEGIN

              print, 'iselect',iselect
              print, 'psi    ',dist.psi[iselect],dist.izero

              psi_val = dist.psi[iselect[icount]] 
              psi_span = ABS(dist.psi[iselect[icount]] - dist.psi[iselect[icount]-1]) * 0.5D  ; parameter
              IF ((ABS(ctr1.psi - psi_val) LT psi_span) OR  $        
                  (ABS(ctr2.psi - psi_val) LT psi_span)) THEN BEGIN
                nskip++
                print,'************Dropping***************',nctr+ictr_new
                status = 0
              ENDIF
              END
            2: psi_val = ref_psi[icount] 
            3: psi_val = ref_psi[icount] 
            4: psi_val = ref_psi[icount] 
            5: psi_val = ref_psi[icount] 
          ENDCASE

          IF (NOT status) THEN CONTINUE

          print, 'psi_val = ',icount,psi_val
;       
;         ----------------------------------------------------------------
;         COLLECT THE CONTOUR AND DIVIDE INTO SEGMENTS
;       
          ctr = grid_ExtractContour(b.psi, b.x, b.y, psi_val)
        
          ibrk = WHERE(ctr.dist GT 0.10D, nbreak)  ; PARAMETER
          IF (nbreak EQ 0) THEN BEGIN
            nseg = 0
            ibrk = [-1,ctr.n-1]
          ENDIF ELSE BEGIN
            nseg = N_ELEMENTS(ibrk) 
            ibrk = [-1,ibrk,ctr.n-1]
          ENDELSE
;       
;         ----------------------------------------------------------------
;         LOOP OVER INDIVIDUAL CONTOUR SEGMENTS AND CHECK LOCATION
;       
          FOR iseg = 0, nseg DO BEGIN
        
            print,'----------------------------------------'
            print,'contour seg:',iseg,ibrk[iseg]+1,ibrk[iseg+1]
        
            ; Not a real segment so skip it, may need to strengthen this check:
            IF (ibrk[iseg]+2 EQ ibrk[iseg+1] OR  $
                ibrk[iseg]+1 EQ ibrk[iseg+1] OR  $
                ibrk[iseg]   EQ ibrk[iseg+1]) THEN CONTINUE  
        
            x = ctr.x[ibrk[iseg]+1:ibrk[iseg+1]]
            y = ctr.y[ibrk[iseg]+1:ibrk[iseg+1]]

 
            IF (state.geometry EQ param.LIMITED) THEN BEGIN

              inside = grid_PointInPolygon(x,y,wall.x,wall.y)
              i = WHERE(inside EQ 0, count)
              proximity = SQRT ( (x[i]-contact_x)^2 + (y[i]-contact_y)^2 )
              dummy = MIN(proximity,imin)
              j = i[imin]
              hold_x = x
              hold_y = y
              n = N_ELEMENTS(x)
              IF (j GT 0) THEN BEGIN
                x = [ hold_x[j:n-1] , hold_x[0:j-1], hold_x[j] ]
                y = [ hold_y[j:n-1] , hold_y[0:j-1], hold_y[j] ]
              ENDIF

            ENDIF

            IF (KEYWORD_SET(debug)) THEN  $
              PLOT,[x,x[0]],[y,y[0]],color=Truecolor('Silver'), PSYM=3,  $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
;       
;           ----------------------------------------------------------------
;           TEST IF THE CONTOUR IS IN THE RIGHT PLACE
;       
;print, 'n',n_elements(x)
print, 'wall2,1',iwall2,iwall1

            IF (iwall2 LE iwall1) THEN BEGIN

              IF (iwall2 EQ iwall1) THEN BEGIN
                wall_pt1 = wall.pt1[*,iwall1]
                wall_pt2 = wall.pt2[*,iwall1]
              ENDIF ELSE BEGIN
                wall_pt1 = wall.pt1[*,iwall2  :iwall1-1]
                wall_pt2 = wall.pt1[*,iwall2+1:iwall1  ]
              ENDELSE

;              wall_pt1 = wall.pt1[*,iwall2:iwall1]
;              wall_pt2 = wall.pt2[*,iwall2:iwall1]
            ENDIF ELSE BEGIN
 print, 'fix please'

              n = N_ELEMENTS(wall.ptc)
              wall_pt1 = [[REFORM(wall.pt1[0,iwall2:n-1]),REFORM(wall.pt1[0,0:iwall1])],  $
                          [REFORM(wall.pt1[1,iwall2:n-1]),REFORM(wall.pt1[1,0:iwall1])]]
              wall_pt2 = [[REFORM(wall.pt2[0,iwall2:n-1]),REFORM(wall.pt2[0,0:iwall1])],  $
                          [REFORM(wall.pt2[1,iwall2:n-1]),REFORM(wall.pt2[1,0:iwall1])]]

              wall_pt1 = TRANSPOSE(wall_pt1)
              wall_pt2 = TRANSPOSE(wall_pt2)
print,wall_pt1
print,wall_pt2

;stop
;              wall_pt1 = [wall.pt1[*,iwall2:n-1],wall.pt1[*,0:iwall1]]
;              wall_pt2 = [wall.pt2[*,iwall2:n-1],wall.pt2[*,0:iwall1]]
            ENDELSE


              PLOT,[wall_pt1[0,*]], $
                   [wall_pt1[1,*]], $
                   color=Truecolor('Lightblue'),  $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE


            FOR i = 0, N_ELEMENTS(x)-2 DO BEGIN  
;             Find intersections with the wall segments:
              result = grid_Intersection([x[i],y[i]],[x[i+1],y[i+1]],  $
                                         wall_pt1,wall_pt2,0,status=status)
              IF (status) THEN BREAK
            ENDFOR
            IF ( ((ipass NE 3 AND i NE N_ELEMENTS(x)-1) OR  $
                  (ipass EQ 3 AND i NE 0              )) AND KEYWORD_SET(debug)) THEN  $
              PLOT,[x[0:i]],[y[0:i]],color=Truecolor('Orange'), PSYM=3,  $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
;       
;           ------------------------------------------------------------------
;           THE CONTOUR INTERSECTS THE WALL IN THE RIGHT PLACE, SO TALLY-HO!
;       
            IF (status) THEN BEGIN
;       
;             ----------------------------------------------------------------
;             TRIM THE OTHER END
;       
              strip_point = 0

              IF (ipass EQ 3) THEN BEGIN
                ; Searching backwards:
                x = x[0:i+1]
                y = y[0:i+1]
                FOR i = N_ELEMENTS(x)-2, 1, -1 DO BEGIN  
                  result = grid_Intersection([x[i],y[i]],[x[i-1],y[i-1]],  $
                                             wall.pt1,wall.pt2,1,status=status)
                  IF (status AND i EQ N_ELEMENTS(x)-2) THEN BEGIN
                    status      = 0
                    strip_point = 1
                  ENDIF
                  IF (status) THEN BREAK
                ENDFOR
                IF (NOT status) THEN BEGIN
                  PRINT,'ERROR grid_RadialRefinement: Second wall intersection not found (1)'
                  STOP
                ENDIF 
                IF (strip_point) THEN BEGIN
                  x = [x[i-1:N_ELEMENTS(x)-3],x[N_ELEMENTS(x)-1]]
                  y = [y[i-1:N_ELEMENTS(y)-3],y[N_ELEMENTS(y)-1]]
                ENDIF ELSE BEGIN
                  x = x[i-1:N_ELEMENTS(x)-1]
                  y = y[i-1:N_ELEMENTS(y)-1]
                ENDELSE
              ENDIF ELSE BEGIN
                ; Searching forwards:
                x = x[i:N_ELEMENTS(x)-1]
                y = y[i:N_ELEMENTS(y)-1]
                FOR i = 1, N_ELEMENTS(x)-2 DO BEGIN  
                  result = grid_Intersection([x[i],y[i]],[x[i+1],y[i+1]],  $
                                             wall.pt1,wall.pt2,1,status=status)
                  ; If there's an intersection on the first point then it's because it's right on
                  ; the line -- this is a rare but encountered problem:
                  IF (status AND i EQ 1) THEN BEGIN
                    status      = 0
                    strip_point = 1
                  ENDIF
                  IF (status) THEN BREAK
                ENDFOR
                IF (NOT status) THEN BEGIN
                  PRINT,'ERROR grid_RadialRefinement: Second wall intersection not found (2)'
                  STOP
                ENDIF 
                IF (strip_point) THEN BEGIN
                  x = [x[0],x[2:i+1]]
                  y = [y[0],y[2:i+1]]
                ENDIF ELSE BEGIN
                  x = x[0:i+1]
                  y = y[0:i+1]
                ENDELSE
              ENDELSE

              IF (KEYWORD_SET(debug)) THEN  $
                PLOT,x,y,color=Truecolor('Cyan'), PSYM=3,  $
                     XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

IF ( nctr+ictr_new EQ 65) then begin

print,'wall_pt12',iwall1,iwall2
print,'wall_ptc',wall.ptc(iwall1),wall.ptt(iwall1)
print,'1:',wall_pt1
print,'2:',wall_pt2


print,x
print,y
;stop
endif 
;       
;             ----------------------------------------------------------------
;             ADD THE CONTOUR TO THE LIST AND EXIT THE LOOP
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
              add_contour = 1
              BREAK
            ENDIF
        
          ENDFOR

        ENDFOR  ; ICOUNT loop

        nskip_total = nskip_total + nskip

        IF (NOT new_contour AND nskip NE icount) THEN BEGIN
          span = 5.00
          xrange = [wall.pt1[0,iwall1]-span,wall.pt1[0,iwall1]+span]
          yrange = [wall.pt1[1,iwall1]-span,wall.pt1[1,iwall1]+span]
          PLOT,wall.x,wall.y,color=Truecolor('Black'),XRANGE=xrange,YRANGE=yrange
          OPLOT,ctr.x,ctr.y,PSYM=3,COLOR=Truecolor('Cyan')
          OPLOT,[wall_pt1[0,*]],[wall_pt1[1,*]],PSYM=6
          OPLOT,[wall_pt2[0,*]],[wall_pt2[1,*]],PSYM=6,COLOR=Truecolor('Lightblue')
          PRINT,'grid_RadialRefinement: New contour not found'
          STOP
        ENDIF

      ENDIF  ; COUNT.GT.0 check
;
;     ------------------------------------------------------------------
;     ADVANCE TO THE NEXT CONTOUR
;
      WHILE (1) DO BEGIN
        status = 0
        iwall1++
        IF (iwall1 GT N_ELEMENTS(wall.ptc)-1) THEN iwall1 = 0

        print,' ++++ iwall1',iwall1,wall.ptc[iwall1],wall.ptt[iwall1],FORMAT='(A,3I6)'
        oplot,[wall.pt1[0,iwall1]],[wall.pt1[1,iwall1]],psym=7,color=Truecolor('BLUE')

        IF (wall.ptc[iwall1] NE -1 AND wall.ptt[iwall1] EQ itarget) THEN BEGIN
          ictr1 = (wall.ptc[iwall1])[0]
          ctr1  = grid_ExtractStructure(c_array,tags[ictr1-1])              
          SWITCH ctr1.region OF
            param.SOL    :     
            param.SOL_LFS: 
            param.SOL_HFS: BEGIN
              status = 1
              IF (ipass EQ 3 AND ctr1.map_out[0] EQ -1) THEN status = 0
              BREAK
              END
            param.PFZ          :            
            param.PFZ_SECONDARY: BEGIN
              status = 1
              IF (ipass NE 3 AND ctr1.map_in [0] EQ -1) THEN status = 0
              BREAK
              END
          ENDSWITCH
        ENDIF 

;       Special case for tangency points:
        IF (wall.ptc[iwall1] NE -1 AND wall.ptt[iwall1] EQ 3) THEN BEGIN
          ictr1 = (wall.ptc[iwall1])[0]
          ctr1  = grid_ExtractStructure(c_array,tags[ictr1-1])              
          SWITCH ctr1.region OF
            param.SOL    : 
            param.SOL_LFS: 
            param.SOL_HFS: BEGIN
              IF (ipass EQ 3) THEN BEGIN
                IF (ctr1.map_out[0] NE -1) THEN BEGIN
                  status = 1
                ENDIF ELSE BEGIN
                  IF (ctr1.map_out[0] EQ -1) THEN BEGIN
                    ; Do nothing, because there's no ring outside the section of
                    ; the contour past the tangency point:
                  ENDIF ELSE BEGIN
                    print, 'radial refinemet: case not handeled yet'
                    stop
                  ENDELSE
                ENDELSE
              ENDIF
              BREAK
              END
            param.PFZ          :
            param.PFZ_SECONDARY: BEGIN
              IF (ipass NE 3) THEN BEGIN
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
              ENDIF
              BREAK
              END
          ENDSWITCH
        ENDIF 
    
        IF (status) THEN BREAK
      ENDWHILE

      ; Exit the loop if back around to the start of the primary separatrix:
      IF (iwall1 EQ istart) THEN BREAK

    ENDWHILE


;if (ipass eq 4) then stop ; dev

  print,'nskip',nskip_total
;stop

;  help,c_array,/struct


;   if (ipass eq 5) then stop

    IF (add_contour) THEN BEGIN
      grid_UpdateWall, b, c_array, wall, $
                       debug=debug, xrange=xrange, yrange=yrange

      c_array = grid_BuildRadialMap(c_array, wall,  $
                                    debug=debug, xrange=xrange, yrange=yrange)
    ENDIF
;if (ipass eq 5) then stop ; dev


 
  ENDFOR ; IPASS loop

  RETURN, c_array

END
