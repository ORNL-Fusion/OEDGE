;
; ======================================================================
;
FUNCTION grid_LoadWall, file, debug=debug, xrange=xrange, yrange=yrange

;  wall = cortex_LoadAnnotationData(1,'/home/ITER/lisgos/divimp/shots/iter/1514/psi_wall.dat') 
;  wall = cortex_LoadAnnotationData(1,'/home/ITER/lisgos/divimp/shots/mast/default/main_wall.dat') ; MAST

  wall = cortex_LoadAnnotationData(1,file)

  n = N_ELEMENTS(wall.x)

help,wall,/struct
  IF (wall.x[0] NE wall.x[n-1] OR  $
      wall.y[0] NE wall.y[n-1]) THEN BEGIN
    wall = { x : [wall.x,wall.x[0]], y : [wall.y,wall.y[0]] }
    n++
  ENDIF
help,wall,/struct

; Rearrange wall data so that it doesn't need to be continuous: 
  pt1 = MAKE_ARRAY(2,n-1,/DOUBLE,VALUE=0.0D)
  pt2 = pt1
  pt1[0,*] = DOUBLE(wall.x[0:n-2])
  pt1[1,*] = DOUBLE(wall.y[0:n-2])
  pt2[0,*] = DOUBLE(wall.x[1:n-1])
  pt2[1,*] = DOUBLE(wall.y[1:n-1])

;   FOR j = 0, n-2 DO BEGIN
;     print,pt1[0,j],pt1[1,j],pt2[0,j],pt2[1,j]
;   ENDFOR

  IF (KEYWORD_SET(debug)) THEN BEGIN
    IF (NOT KEYWORD_SET(xrange)) THEN xrange = [ 3.0,9.0]
    IF (NOT KEYWORD_SET(yrange)) THEN yrange = [-6.0,6.0]
    PLOT,wall.x,wall.y,color=Truecolor('Yellow'),XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1
  ENDIF


  grid_ZoneWall, wall, pt1, pt2, debug=debug, xrange=xrange, yrange=yrange



   
  result = {  $
    x   : wall.x ,  $
    y   : wall.y ,  $
    pt1 : pt1    ,  $
    pt2 : pt2    }

  RETURN, result
END

;
; ======================================================================
;
FUNCTION grid_Length, x, y

  result = 0.0D

  FOR i = 0, N_ELEMENTS(x)-2 DO BEGIN
    result = result + SQRT( (x[i] - x[i+1])^2 + (y[i] - y[i+1])^2 )
  ENDFOR

  RETURN, result
END
;
; ======================================================================
;
FUNCTION grid_ExtractStructure, struct, tag

  index = WHERE(TAG_NAMES(struct) EQ STRUPCASE(tag), count)

  IF (count NE 0) THEN BEGIN
    result = struct.(index[0]) 
  ENDIF ELSE BEGIN
    PRINT, 'ERROR grid_ExtractStructure: TAG not found'
    PRINT, '  TAG = ',tag
    result = -1
  ENDELSE

  RETURN, result
END
;
; ======================================================================
;
FUNCTION grid_ReplaceStructure, struct_array,struct,tag

  struct_tags = STRUPCASE(TAG_NAMES(struct_array))

  FOR i = 0, N_ELEMENTS(struct_tags)-1 DO BEGIN
    IF (struct_tags[i] EQ STRUPCASE(tag)) THEN BEGIN
      val = struct 
      print, 'replacing! ',tag
      print,val.state
    ENDIF ELSE  $
      val = grid_ExtractStructure(struct_array,struct_tags[i])  

    IF (i EQ 0) THEN  $
      result = CREATE_STRUCT(       struct_tags[i],val) ELSE  $
      result = CREATE_STRUCT(result,struct_tags[i],val)
  ENDFOR

  RETURN, result
END
;
;  ======================================================================
;
PRO grid_Debug, contours, wall, xrange, yrange


  PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
  PLOT,wall.x,wall.y,color=Truecolor('Yellow'),  $
       XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE


;  help,contours,/struct

  tags = STRUPCASE(TAG_NAMES(contours))

  FOR i = 1, N_ELEMENTS(tags) DO BEGIN

    name = 'contour' + STRING(i,FORMAT='(I0)')

    print,name,i, N_ELEMENTS(tags)

    cnt = grid_ExtractStructure(contours,name)  

  
    PLOT,cnt.x,cnt.y,color=Truecolor('Red'),  $
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,PSYM=3,/NOERASE

    IF (i EQ N_ELEMENTS(tags)) THEN BEGIN
print,cnt.boundary1_p1[0],cnt.boundary1_p2[0]
print,cnt.boundary1_p1[1],cnt.boundary1_p2[1]

      OPLOT,[cnt.boundary1_p1[0],cnt.boundary1_p2[0]],  $
            [cnt.boundary1_p1[1],cnt.boundary1_p2[1]],  $
            color=Truecolor('Lightgreen')

      OPLOT,[cnt.boundary2_p1[0],cnt.boundary2_p2[0]],  $
            [cnt.boundary2_p1[1],cnt.boundary2_p2[1]],  $
            color=Truecolor('Lightgreen')

    ENDIF

    XYOUTS, [cnt.x[0]], [cnt.y[0]], CHARSIZE=1.0, ALIGNMENT=0.5, $
            STRTRIM(STRING(i),2), COLOR=TrueColor('White')

  ENDFOR

;  SAVE,filename='contours2.sav',contours

  stop

END
;
; ======================================================================
;
; need to add some efficiency, i.e. domain decomposition (as per usual)
;
FUNCTION grid_Intersection, p1, p2, seg1, seg2, mode, status=status
  ; --------------------------------------------------------------------
  COMMON zone, zone_list, zone_index, zone_nx  , zone_ny  ,  $
               zone_minx, zone_maxx , zone_miny, zone_maxy
  ; --------------------------------------------------------------------

  a1 = p1
  a2 = p2

  CASE mode OF
;   ----------------------------------------------------------------------
    0: BEGIN
      seg_list = LINDGEN(N_ELEMENTS(seg1[0,*]))
      END
;   ----------------------------------------------------------------------
    1: BEGIN

      deltax = (zone_maxx - zone_minx) / DOUBLE(zone_nx)
      deltay = (zone_maxy - zone_miny) / DOUBLE(zone_ny)
     
      i1 = LONG((a1[0] - zone_minx) / deltax)
      j1 = LONG((a1[1] - zone_miny) / deltay)
      i2 = LONG((a2[0] - zone_minx) / deltax)
      j2 = LONG((a2[1] - zone_miny) / deltay)

;     Basic check to see if points are inside the region of interest.  Note that
;     this check can give a false result if the line segments span the region but
;     the end points themselves are outside, but this is (very) unlikely to happen:     
      IF ((i1 LT 0 OR i1 GT zone_nx-1) OR  $
          (j1 LT 0 OR j1 GT zone_ny-1) OR  $
          (i2 LT 0 OR i2 GT zone_nx-1) OR  $
          (j1 LT 0 OR j1 GT zone_ny-1)) THEN BEGIN
        status =  0
        result = -1
        RETURN, result
      ENDIF

;      print,i1,j1,i2,j2
;      print,zone_nx,zone_ny
;      print, '-->',zone_index[0,i1,j1],zone_index[1,i1,j1]

      IF (i1 EQ i2 AND j1 EQ j2) THEN BEGIN
        IF (zone_index[0,i1,j1] EQ -1) THEN BEGIN
;         No wall segments associated with this zone:
          status =  0
          result = -1
          RETURN, result
        ENDIF
        seg_list = zone_list[zone_index[0,i1,j1]:zone_index[1,i1,j1]]
      ENDIF ELSE BEGIN 
;       More than one zone involved or one end point outside zoned region, 
;       so search everywhere:
        seg_list = LINDGEN(N_ELEMENTS(seg1[0,*])) 
      ENDELSE                             
      END
;   --------------------------------------------------------------------
    ELSE: BEGIN
      PRINT,'ERROR grid_Intersect: Unknown MODE'
      PRINT,'  MODE=',mode
      END
  ENDCASE


  FOR j = 0, N_ELEMENTS(seg_list)-1 DO BEGIN

    i = seg_list[j]        

    b1 = seg1[*,i]
    b2 = seg2[*,i]

    IF ( ABS(b1[0]-b2[0]) LT 1.0D-6 ) THEN b1[0] = b1[0] + 0.00001D
  
    Lint, a1, a2, b1, b2, c, d, flag=flag
  
    IF (flag NE 1) THEN CONTINUE
  
    s = (c[0] - a1[0]) / (a2[0] - a1[0])    
    t = (c[0] - b1[0]) / (b2[0] - b1[0])    

    IF ( s GE 0.0D AND s LE 1.0D AND  $
         t GE 0.0D AND t LT 1.0D) THEN BEGIN
  
      IF (NOT KEYWORD_SET(int_i)) THEN BEGIN
        int_i    = [i]
        int_s    = [s]
        int_t    = [t]      
        int_x    = [c[0]]
        int_y    = [c[1]]
      ENDIF ELSE BEGIN
        int_i    = [int_i,i]
        int_s    = [int_s,s]
        int_t    = [int_t,t]
        int_x    = [int_x,c[0]]
        int_y    = [int_y,c[1]]
      ENDELSE
  
    ENDIF
  ENDFOR
  
  IF (NOT KEYWORD_SET(int_t)) THEN BEGIN
;   No intersection found between the line and segment list:
    status =  0 ; Can only be 0 or 1 since used as a logical variable below
    result = -1
  ENDIF ELSE BEGIN
    status = 1  ; Can only be 0 or 1 since used as a logical variable below
    result = {          $
      i    : int_i   ,  $
      s    : int_s   ,  $
      t    : int_t   ,  $ 
      x    : int_x   ,  $
      y    : int_y   }          
  ENDELSE

  RETURN,result

END
;
;
; ======================================================================
;
; From http://www.dfanning.com/tips/point_in_polygon.html.
;
FUNCTION grid_PointInPolygon, x, y, px, py

;  n = N_ELEMENTS(px)
;  IF (px[0] NE px[n-1] OR py[0] NE py[n-1]) THEN BEGIN
;    px = [px,px[0]]
;  ENDIF

  object = Obj_New('IDLanROI', px, py)
  result = object->ContainsPoints(x, y)
  Obj_Destroy, object

  RETURN, result
END
;
;
; ======================================================================
;
FUNCTION grid_ExtractContour, psi, psi_x, psi_y, psi_val  ;, xpoint_zone=xpoint_zone

    CONTOUR, psi, psi_x, psi_y, NLEVELS=20, c_labels=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
stop

; x_range = !X.RANGE  ; Not sure why !X. and !Y. are zero...
; y_range = !Y.RANGE

; Generate the contour data for PSI_VAL:
  CONTOUR, psi, psi_x, psi_y, LEVELS=[psi_val], /PATH_DATA_COORDS,  $
           PATH_XY=path_xy, PATH_INFO=path_info, /PATH_DOUBLE

; !X.RANGE = x_range
; !Y.RANGE = y_range

  x = REFORM(path_xy[0,*])
  y = REFORM(path_xy[1,*])

;  IF (KEYWORD_SET(xpoint_zone)) THEN  $
;    i = WHERE( x GE xpoint_zone[0] AND x LE xpoint_zone[1]  AND  $
;               y GE xpoint_zone[2] AND y LE xpoint_zone[3], count_pts) ELSE  $
;    i = WHERE( x NE -999.0 AND y NE -999.0, count_pts)
;  
;  IF (count_pts EQ 0) THEN BEGIN
;    PRINT, 'no contours in x-point zone'
;  ENDIF

  n = N_ELEMENTS(x)

; Find the distance between consectutive points:
  dist = MAKE_ARRAY(n,/FLOAT,VALUE=0.0)
  FOR i = 0, n-2 DO dist[i] = SQRT((x[i]-x[i+1])^2 + (y[i]-y[i+1])^2)

  result = { dist : dist, x : x, y : y, n : n } ; , count_pts : count_pts }

  RETURN, result
END
;
;
; ======================================================================
;
PRO grid_SetActiveContour, contour_array, int_nlast, region, active_contour, store_forced,  $
                           psi_start, psi_end, psi_val, psi_adjust, psi_step, sht_count, last_good_psi_val,  $
                           boundary1_p1, boundary1_p2, boundary2_p1, boundary2_p2,  $
                           the_search_continues, contour_n, direction, process_2nd,  $
                           psi_1st_separatrix, psi_2nd_separatrix, user_step
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, FORWARD, BACKWARD
  ; ------------------------------------------------------------------

  ; Loop through contours and find those with a tangecy point,
  ; but where the radial progression of the grid on either side of the 
  ; tangecy point (in the poloidal direction) has not been processed yet: 

  n = N_ELEMENTS(TAG_NAMES(contour_array))
  status = 0

  FOR state = 0, 1 DO BEGIN
    FOR i = 1, n DO BEGIN  
      tag = 'contour' + STRING(i,FORMAT='(I0)')

      cnt = grid_ExtractStructure(contour_array,tag)

      IF (cnt.tangent_i EQ -1) THEN CONTINUE

      print,'TAG:-----------> ',state,' ',tag,cnt.state


;     help,cnt,/struct


      IF (state EQ cnt.state) THEN BEGIN
        CASE cnt.state OF
          0: BEGIN
            status = 1
            cnt.state = 1
            psi_val = cnt.psi
            psi_start = psi_val
            direction = cnt.direction
            int_nlast[region] = 2
            boundary1_p1 = cnt.boundary1_p1
            boundary1_p2 = cnt.boundary1_p2
            boundary2_p1 = cnt.tangent_p1
            boundary2_p2 = cnt.tangent_p2
            IF (region EQ SOL AND direction EQ backward) THEN BEGIN
              IF (process_2nd EQ 2) THEN psi_end = psi_2nd_separatrix ELSE $
                                         psi_end = psi_1st_separatrix 
           ;  *** NEED TO KEEP TRACK OF WHICH END POINT RAN A-FOUL AND 
           ;      USE THAT INFORMATION HERE TO DECIDE WHAT TO DO ***
             ENDIF
            END
          1: BEGIN
            status = 1
            cnt.state = 2
            psi_val = cnt.psi
            psi_start = psi_val
            direction = cnt.direction
            int_nlast[region] = 2
            boundary1_p1 = cnt.tangent_p1
            boundary1_p2 = cnt.tangent_p2
            boundary2_p1 = cnt.boundary2_p1
            boundary2_p2 = cnt.boundary2_p2
            IF (region EQ SOL AND direction EQ backward) THEN BEGIN
              direction = FORWARD
              psi_end = psi_1st_separatrix - 3.0D
              boundary1_p1 = [-999.0D,-999.0D]
              boundary1_p2 = [-999.0D,-999.0D]
            ENDIF
            END
          ELSE: BEGIN
            PRINT,'ERROR grid_SetActiveContour: Unknown contour state'
            PRINT,' TAG   = ',tag
            PRINT,' STATE = ',cnt.state
            STOP
            END
        ENDCASE
        IF (status) THEN BEGIN
          contour_array = grid_ReplaceStructure(contour_array,cnt,tag)
          BREAK
        ENDIF
      ENDIF
    ENDFOR
    IF (status) THEN BREAK
  ENDFOR

  ; Final exit condition for PSI_VAL loop:
  IF (NOT status) THEN the_search_continues = 0 

  active_contour = i
  store_forced   = 0

  IF (region EQ SOL AND direction EQ BACKWARD) THEN   $
    psi_step = user_step * direction ELSE  $
    psi_step = user_step * direction

  psi_adjust = psi_step 

  sht_count = 0
  last_good_psi_val = 0.0D

      print,'  boundary1:',FLOAT(boundary1_p1),FLOAT(boundary1_p2),FORMAT='(A,4F10.4)'
      print,'  boundary2:',FLOAT(boundary2_p1),FLOAT(boundary2_p2),FORMAT='(A,4F10.4)'

END
;
;
; ======================================================================
;
PRO grid_RefineContour, x, y, i, xrange, yrange

  n_pts = 10

  i1 = MAX([0            ,i-1 ]  )
  j1 = MAX([0            ,i-10]  )
  i2 = MIN([i+2 ,N_ELEMENTS(x)-1])
  j2 = MIN([i+11,N_ELEMENTS(x)-1])
  ref_x = x[j1:j2]
  ref_y = y[j1:j2]
  frac = (DINDGEN(n_pts-1) + 1.0D) / DOUBLE(n_pts)
  new_x = x[i1] + frac * (x[i2] - x[i1])
  new_y = y[i1] + frac * (y[i2] - y[i1])
  IF (ABS(x[i]-x[i+1]) LT  $
      ABS(y[i]-y[i+1])) THEN BEGIN
    IF (ref_y[0] GT ref_y[1]) THEN BEGIN
      new_y = REVERSE(new_y)
      ref_x = REVERSE(ref_x)
      ref_y = REVERSE(ref_y)
      status = 1
;      print,'reversing y!'
    ENDIF ELSE status = 0
    new_x = SPLINE( ref_y, ref_x, new_y ,/DOUBLE)              
  ENDIF ELSE BEGIN
    IF (ref_x[0] GT ref_x[1]) THEN BEGIN
      new_x = REVERSE(new_x)
      ref_x = REVERSE(ref_x)
      ref_y = REVERSE(ref_y)
      status = 1
;      print,'reversing x!'
    ENDIF ELSE status = 0
    new_y = SPLINE( ref_x, ref_y, new_x ,/DOUBLE)              

;print,i1,i2
;print,x[i1],x[i2]
;print,ref_x
;print,ref_y
;print,new_x
;print,new_y
;PLOT,[new_x],[new_y],color=Truecolor('Green'), PSYM=6,  $
;     XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
;PLOT,x,y,color=Truecolor('Orange'),  $
;     XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,PSYM=3,/NOERASE
;stop
  ENDELSE

  IF (status) THEN BEGIN
    new_x = REVERSE(new_x)
    new_y = REVERSE(new_y)
  ENDIF

  PLOT,[new_x],[new_y],color=Truecolor('Red'), PSYM=6,  $
       XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

  x = [x[0:i1],new_x,x[i2:N_ELEMENTS(x)-1]]
  y = [y[0:i1],new_y,y[i2:N_ELEMENTS(y)-1]]

END

;
; ======================================================================
;
PRO grid_ZoneWall, wall, wall_pt1, wall_pt2, debug=debug, xrange=xrange, yrange=yrange

  COMMON zone, zone_list, zone_index, zone_nx  , zone_ny  ,  $
               zone_minx, zone_maxx , zone_miny, zone_maxy

; Setup zones:
  zone_nx = 10
  zone_ny = 10
  zone_list  = MAKE_ARRAY(2*N_ELEMENTS(wall_pt1[0,*]),/LONG,VALUE=-1L)
  zone_index = MAKE_ARRAY(2,zone_nx,zone_ny,/LONG,VALUE=-1L)
  zone_minx = MIN(wall.x) - 0.1D
  zone_maxx = MAX(wall.x) + 0.1D
  zone_miny = MIN(wall.y) - 0.1D
  zone_maxy = MAX(wall.y) + 0.1D
  deltax = (zone_maxx - zone_minx) / DOUBLE(zone_nx)
  deltay = (zone_maxy - zone_miny) / DOUBLE(zone_ny)

  box_pt1 = MAKE_ARRAY(2,5,/DOUBLE,VALUE=0.0D)
  box_pt2 = box_pt1

  x = zone_minx
  count = -1
  FOR i = 0, zone_nx-1 DO BEGIN
    y = zone_miny
    FOR j = 0, zone_ny-1 DO BEGIN
      box_pt1[0,*] = [x       , x         , x + deltax, x + deltax, x       ]
      box_pt1[1,*] = [y       , y + deltay, y + deltay, y         , y       ] 
      box_pt2[0,*] = [x       , x + deltax, x + deltax, x         , x       ]
      box_pt2[1,*] = [y+deltay, y + deltay, y         , y         , y+deltay] 

      FOR k = 0, N_ELEMENTS(wall_pt1[0,*])-1 DO BEGIN
        status = 0
        ; Check if end points are inside the zone:
        IF (wall_pt1[0,k] GE x AND wall_pt1[0,k] LE x+deltax AND $
            wall_pt1[1,k] GE y AND wall_pt1[1,k] LE y+deltay) THEN status = 1
        IF (status EQ 0) THEN BEGIN
          IF (wall_pt2[0,k] GE x AND wall_pt2[0,k] LE x+deltax AND $
              wall_pt2[1,k] GE y AND wall_pt2[1,k] LE y+deltay) THEN status = 1
          ; If end points not inside the zone, check for intersections of the 
          ; line segment with the zone boundary:
          IF (status EQ 0) THEN  $
            result = grid_Intersection([wall_pt1[0,k],wall_pt1[1,k]],  $
                                       [wall_pt2[0,k],wall_pt2[1,k]],  $
                                       box_pt1, box_pt2, 0, status=status)
        ENDIF
        ; Register the wall segment if it is associated wit the zone, noting
        ; that a given segment can interact with more than one zone:
        IF (status GT 0) THEN BEGIN
          count = count + 1
          zone_list[count] = k
          IF (zone_index[0,i,j] EQ -1) THEN zone_index[0,i,j] = count
          zone_index[1,i,j] = count
        ENDIF
      ENDFOR  

;      print, 'zone_index',i,j,zone_index[*,i,j]    
      IF (KEYWORD_SET(debug)) THEN BEGIN
        PLOT,box_pt1[0,*],box_pt1[1,*],color=Truecolor('Purple'),  $
             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
      ENDIF

      y = y + deltay
    ENDFOR
    x = x + deltax
  ENDFOR
END
;
; ======================================================================
;
PRO grid_ClipToWall, x, y, wall_pt1, wall_pt2, cycle, xrange, yrange

  cycle = 0

  cont  = 1
  WHILE (cont) DO BEGIN
    cont = 0

    ; Find intersections with wall starting at beginning of the wall:
    FOR i = 0, N_ELEMENTS(x)-2 DO BEGIN
      result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                 wall_pt1, wall_pt2, 1, status=status)
      IF (status) THEN BREAK
    ENDFOR
    IF (NOT status) THEN BEGIN
      PRINT, 'NO WALL INTERSECTION FOUND 1!'
      cycle = 1
      BREAK
    ENDIF
    IF (N_ELEMENTS(result.i) EQ 1) THEN BEGIN
      x = x[i:N_ELEMENTS(x)-1]
      y = y[i:N_ELEMENTS(y)-1]
    ENDIF ELSE BEGIN
      grid_RefineContour, x, y, i, xrange, yrange
      cont = 1
    ENDELSE

    ; Find intersections with wall starting at the end of the contour:
    FOR i = N_ELEMENTS(x)-2, 0, -1 DO BEGIN
      result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                 wall_pt1, wall_pt2, 1, status=status)
      IF (status) THEN BREAK
    ENDFOR
    IF (NOT status) THEN BEGIN
      PRINT, 'NO WALL INTERSECTION FOUND 2!'
      cycle = 1
      BREAK
    ENDIF
    IF (N_ELEMENTS(result.i) EQ 1) THEN BEGIN
      x = x[0:i+1]
      y = y[0:i+1]
    ENDIF ELSE BEGIN
      grid_RefineContour, x, y, i, xrange, yrange
      cont = 1
;       oplot,result.x,result.y,psym=6,color=Truecolor('Red')
;       help,result,/struct
;       stop

    ENDELSE
  ENDWHILE

      oplot,x,y,psym=3,color=Truecolor('Purple')


END
;
; ======================================================================
;
; How do I avoid the distorted target cells that plagued the JET grid, i.e. the 
; cells overlapped?  Just keep the resolution high? 

; Need to bring in the other wall segments now I think, and look for intersections, rather 
; than checking if a point is inside or outside (although good to have some appropriate
; spatial refinement before this step? to make sure the entrances and exits are all accounted for?)

;      inside = grid_PointInPolygon(x,y,wall.x,wall.y)
;      i = WHERE(inside EQ 1, count)
;      IF (count EQ 0) THEN CONTINUE  
;      n = N_ELEMENTS(i) - 1
;      m = N_ELEMENTS(x) - 1
;      IF (i[0] GT 0) THEN j = i[0] - 1 ELSE j = 0  ; Should check if contour entirely inside wall, and
;      IF (i[n] LT m) THEN k = i[n] + 1 ELSE k = m  ; if this matches expectations...
;      x      = x[j:k]
;      y      = y[j:k]
;      inside = inside[j:k]
;      inside[0] = 1 & inside[N_ELEMENTS(inside)-1] = 1

            ;  xrange = [ 3.0,9.0]
            ;  yrange = [-6.0,6.0]
            ;  PLOT,[boundary1_p1[0],boundary1_p2[0]],  $
            ;       [boundary1_p1[1],boundary1_p2[1]],  $
            ;       color=Truecolor('Lightgreen'),  $
            ;       XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1
            ;  PLOT,[boundary2_p1[0],boundary2_p2[0]],  $
            ;       [boundary2_p1[1],boundary2_p2[1]],  $
            ;       color=Truecolor('Green'),  $
            ;       XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
            ;  PLOT,x,y,  $
            ;       color=Truecolor('Purple'), PSYM=3, $
            ;       XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
            ;  PLOT,wall.x,wall.y,color=Truecolor('Yellow'),  $
            ;       XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

FUNCTION grid_AnalyseBoundary, b, wall, user_step, finish, xrange=user_xrange, yrange=user_yrange, debug=debug
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, FORWARD, BACKWARD

  FORWARD       =  1.0D
  BACKWARD      = -1.0D
  CORE          = 0
  SOL           = 1
  SOL_LFS       = 2  
  SOL_HFS       = 3  
  PFZ           = 4  
  PFZ_SECONDARY = 5  
  ; ------------------------------------------------------------------

  IF (NOT KEYWORD_SET(user_xrange)) THEN user_xrange = [ 3.0,9.0]
  IF (NOT KEYWORD_SET(user_yrange)) THEN user_yrange = [-6.0,6.0]

  psi   = b.psi_raw
  psi_x = b.x
  psi_y = b.y

  psi_1st_xpoint = b.psi[b.null_i[1],b.null_j[1]]
  psi_2nd_xpoint = b.psi[b.null_i[2],b.null_j[2]]

  region = SOL

  direction = FORWARD

  psi_val    = DOUBLE(psi_1st_xpoint) + 0.001D 

;  psi_val    = DOUBLE(psi_1st_xpoint) - 0.00001D ; MAST
  psi_start  = psi_val
  psi_end    = psi_1st_xpoint - finish
  psi_step   = user_step * direction
  psi_adjust = psi_step
  tol_adjust = 1.0D-5


  ; Setup local wall segment arrays:
  wall_pt1 = wall.pt1
  wall_pt2 = wall.pt2


  ; Initialisation:

  int_nlast = MAKE_ARRAY(6,/LONG,VALUE=-1)

  contour_n = 0

  boundary1_p1 = [-999.0D,-999.0D]
  boundary1_p2 = [-999.0D,-999.0D]
  boundary2_p1 = [-999.0D,-999.0D]
  boundary2_p2 = [-999.0D,-999.0D]

  save_boundary1_p1 = [-999.0D,-999.0D]
  save_boundary1_p2 = [-999.0D,-999.0D]
  save_boundary2_p1 = [-999.0D,-999.0D]
  save_boundary2_p2 = [-999.0D,-999.0D]

  the_search_continues = 1
  active_contour       = 0
  store_forced         = 0

  IF (0 EQ 1) THEN BEGIN
    xrange = [ 3.0,9.0]
    yrange = [-6.0,6.0]
    xrange = [ 4.0,5.5]
    yrange = [ 4.0,5.0]
    PLOT,wall.x,wall.y,color=Truecolor('Yellow'),  $
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1
    ctr = grid_ExtractContour(psi, psi_x, psi_y, psi_1st_xpoint)
    PLOT,ctr.x,ctr.y,color=Truecolor('Green'),  $
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
    ctr = grid_ExtractContour(psi, psi_x, psi_y, psi_2nd_xpoint)
    PLOT,ctr.x,ctr.y,color=Truecolor('Orange'),  $
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
    stop
  ENDIF


  last_psi_val      = -1.0D
  last_good_psi_val =  0.0D

  seg_active = MAKE_ARRAY(100,/LONG,VALUE=0L)

  iteration_count = 0

  ; Decide if the secondary x-point is inside the main wall:
  inside = grid_PointInPolygon(b.x[b.null_i[2]],b.y[b.null_j[2]],wall.x,wall.y)
  i = WHERE(inside EQ 1, count)
  IF (count EQ 0) THEN process_2nd = -1 ELSE  $
                       process_2nd =  0

;
; ----------------------------------------------------------------------
; MAIN LOOP
;  
  WHILE (the_search_continues) DO BEGIN

    iteration_count = iteration_count + 1


    IF (0 EQ 1 AND contour_n EQ 0) THEN BEGIN
      RESTORE,'contours.sav'
      region = SOL
      contour_array = contours
      active_contour = 20
      contour_n = 20
      grid_SetActiveContour, contour_array, int_nlast, region, active_contour, store_forced,  $
                             psi_start, psi_end, psi_val, psi_adjust, psi_step, sht_count, last_good_psi_val,   $
                             boundary1_p1, boundary1_p2, boundary2_p1, boundary2_p2,  $
                             the_search_continues, contour_n, direction, process_2nd,  $
                             psi_1st_xpoint, psi_2nd_xpoint, user_step
      the_search_continues = 1
;      psi_start = psi_2nd_xpoint - 0.1D 
;      psi_end   = psi_1st_xpoint 
;      psi_val   = psi_start
;      direction = BACKWARD
;      psi_step  = -0.005D * direction

      psi_start = psi_2nd_xpoint - 0.0001D
      psi_val = psi_start

      process_2nd = 2
      int_nlast[SOL] = 2
    ENDIF

    print,'===========',contour_n,iteration_count,'============='
    print,'psi_values: ',psi_start,psi_val,psi_end,process_2nd, direction

    IF (process_2nd EQ 0 AND  $
        psi_adjust NE psi_step AND psi_val LT psi_2nd_xpoint) THEN  $
      process_2nd = 1


    ctr = grid_ExtractContour(psi, psi_x, psi_y, psi_val)

    ibrk = WHERE(ctr.dist GT 0.10)
    nseg = N_ELEMENTS(ibrk) 
    ibrk = [-1,ibrk,ctr.n-1]

    IF (1 EQ 1 AND iteration_count EQ 2)  THEN BEGIN
      xrange = user_xrange
      yrange = user_yrange
;      xrange = [ 4.5,5.5]
;      yrange = [ 4.0,5.5]
;      xrange = [ 3.8,4.2]
;      yrange = [-2.0,-1.0]
      PLOT,ctr.x,ctr.y,color=Truecolor('White'), PSYM=3, $ 
           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1
      PLOT,[ctr.x[0]],[ctr.y[0]],color=Truecolor('Red'), PSYM=6, $
      XRANGE=xrange,YRANGE=yrange,/NOERASE,XSTYLE=1,YSTYLE=1
         PLOT,[boundary1_p1[0],boundary1_p2[0]],  $
              [boundary1_p1[1],boundary1_p2[1]],  $
              color=Truecolor('Lightgreen'),  $
              XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
         PLOT,[boundary2_p1[0],boundary2_p2[0]],  $
              [boundary2_p1[1],boundary2_p2[1]],  $
              color=Truecolor('Green'),  $
              XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
         PLOT,wall.x,wall.y,color=Truecolor('Yellow'),  $
              XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

print,'... nseg',nseg
      stop
    ENDIF

    ; Flag to register when at least one segment of the contour is in the region of interest:
    seg_valid = 0  
;
;   --------------------------------------------------------------------
;   LOOP OVER INDIVIDUAL CONTOUR SEGMENTS
;  
    FOR iseg = 0, nseg DO BEGIN

      ; Skip the rest of the segments in the list if a properly processed segment has
      ; already been found:
      IF (seg_valid EQ 1) THEN CONTINUE

      ; Special case to deal with contour trajectory near the 
      ; secondary x-point:
      IF (save_boundary1_p1[0] NE -999.0D0) THEN BEGIN
        IF (seg_active[iseg] EQ 1 AND psi_val LT psi_2nd_xpoint) THEN BEGIN
          boundary1_p1 = new_boundary1_p1
          boundary1_p2 = new_boundary1_p2
          boundary2_p1 = new_boundary2_p1
          boundary2_p2 = new_boundary2_p2
        ENDIF ELSE BEGIN
          boundary1_p1 = save_boundary1_p1
          boundary1_p2 = save_boundary1_p2
          boundary2_p1 = save_boundary2_p1
          boundary2_p2 = save_boundary2_p2
        ENDELSE
      ENDIF

      print,'----------------------------------------'
      print,'contour seg:',iseg,ibrk[iseg]+1,ibrk[iseg+1],seg_valid

      ; Not a real segment, may need to strengthen this check:
      IF (ibrk[iseg]+1 EQ ibrk[iseg+1] OR  $
          ibrk[iseg]   EQ ibrk[iseg+1]) THEN CONTINUE  

      x = ctr.x[ibrk[iseg]+1:ibrk[iseg+1]]
      y = ctr.y[ibrk[iseg]+1:ibrk[iseg+1]]

      ;        xrange = user_xrange
      ;        yrange = user_yrange
      ;        PLOT,[boundary1_p1[0],boundary1_p2[0]],  $
      ;             [boundary1_p1[1],boundary1_p2[1]],  $
      ;             color=Truecolor('Lightgreen'),  $
      ;             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1
      ;        PLOT,[boundary2_p1[0],boundary2_p2[0]],  $
      ;             [boundary2_p1[1],boundary2_p2[1]],  $
      ;             color=Truecolor('Green'),  $
      ;             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
      ;        PLOT,x,y,  $
      ;             color=Truecolor('Red'), PSYM=3, $
      ;             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
      ;        PLOT,wall.x,wall.y,color=Truecolor('Yellow'),  $
      ;             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
      ;       if (contour_n EQ 2 AND iseg EQ 2) then stop


      ; Isolate region of interest along the contour by clipping at the poloidal 
      ; boundaries:

      store_x = x
      store_y = y
      cont    = 1
      WHILE (cont EQ 1 OR cont EQ 3) DO BEGIN
        cont = 0

        print,'  boundary1:',FLOAT(boundary1_p1),FLOAT(boundary1_p2),FORMAT='(A,4F10.4)'
        print,'  boundary2:',FLOAT(boundary2_p1),FLOAT(boundary2_p2),FORMAT='(A,4F10.4)'
              print,'psi...',psi_val,psi_2nd_xpoint,psi_val GT psi_2nd_xpoint
        print,N_ELEMENTS(x)

        failure = 0 

        IF (boundary1_p1[0] NE -999.0D) THEN BEGIN
          ; Find intersections with first boundary:
          FOR i = N_ELEMENTS(x)-2, 1, -1 DO BEGIN
            result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                       boundary1_p1, boundary1_p2, 0, status=status)
            IF (status) THEN BREAK
          ENDFOR
          IF (NOT status) THEN BEGIN
            PRINT, 'NO INTERSECTION WITH BOUNDARY 1 FOUND!'
            failure = failure + 1
          ENDIF ELSE BEGIN
            IF (N_ELEMENTS(result.i) EQ 1) THEN BEGIN
              x = x[i:N_ELEMENTS(x)-1]
              y = y[i:N_ELEMENTS(y)-1]
              x[0] = result.x[0]
              y[0] = result.y[0]
            ENDIF ELSE BEGIN
               print,'problem here'
               stop
            ENDELSE
          ENDELSE
        ENDIF
        
        IF (boundary2_p1[0] NE -999.0D) THEN BEGIN
          ; Find intersections with second boundary:
          FOR i = 0, N_ELEMENTS(x)-2 DO BEGIN
            result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                       boundary2_p1, boundary2_p2, 0, status=status)
            IF (status) THEN BREAK
          ENDFOR
          IF (NOT status) THEN BEGIN
            PRINT, 'NO INTERSECTION WITH BOUNDARY 2 FOUND!'
            failure = failure + 2
          ENDIF ELSE BEGIN
            IF (N_ELEMENTS(result.i) EQ 1) THEN BEGIN
              x = x[0:i]
              y = y[0:i]
              x[i] = result.x[0]
              y[i] = result.y[0]
            ENDIF ELSE BEGIN
               print,'problem here as well'
               stop
            ENDELSE
          ENDELSE
        ENDIF
        
        print, 'failure', failure

        IF ((failure EQ 1 AND boundary2_p1[0] EQ -999.0D) OR  $
            (failure EQ 2 AND boundary1_p1[0] EQ -999.0D) OR  $
             failure EQ 3) THEN  $
          cont = 2  $
        ELSE BEGIN
          IF (failure EQ 1 OR failure EQ 2) THEN BEGIN
        
            FOR i = 0, N_ELEMENTS(x)-2 DO BEGIN
              result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                         wall_pt1, wall_pt2, 1, status=status)
              IF (status) THEN BREAK
            ENDFOR
            IF (NOT status) THEN cont = 2  $
            ELSE BEGIN

              print,'checking...'
              print,'psi...',psi_val,psi_2nd_xpoint

              xrange = user_xrange
              yrange = user_yrange
              PLOT,[boundary1_p1[0],boundary1_p2[0]],  $
                   [boundary1_p1[1],boundary1_p2[1]],  $
                   color=Truecolor('Lightgreen'),  $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1
              PLOT,[boundary2_p1[0],boundary2_p2[0]],  $
                   [boundary2_p1[1],boundary2_p2[1]],  $
                   color=Truecolor('Green'),  $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
              PLOT,x,y,  $
                   color=Truecolor('Red'), PSYM=3, $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
              PLOT,wall.x,wall.y,color=Truecolor('Yellow'),  $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
;              if (contour_n EQ 21) then stop


              IF (cont EQ 3) THEN stop

              IF (     psi_val LT psi_2nd_xpoint AND  $
                  last_psi_val GT psi_2nd_xpoint) THEN BEGIN

                seg_active[iseg] = 1

                save_boundary1_p1 = boundary1_p1
                save_boundary1_p2 = boundary1_p2
                save_boundary2_p1 = boundary2_p1
                save_boundary2_p2 = boundary2_p2
                print,MEAN(x,/DOUBLE),b.x[b.null_i[2]]
              
                IF (MEAN(x,/DOUBLE) LT b.x[b.null_i[2]]) THEN BEGIN
                  print, '>>>>>>>>> turning off second boundary'
                  new_boundary1_p1 = boundary1_p1
                  new_boundary1_p2 = boundary1_p2
                  new_boundary2_p1 = [-999.0D0,-999.0D0]
                  new_boundary2_p2 = [-999.0D0,-999.0D0]
                ENDIF ELSE BEGIN
                  print, '>>>>>>>>> turning off first boundary'
                  new_boundary1_p1 = [-999.0D0,-999.0D0]
                  new_boundary1_p2 = [-999.0D0,-999.0D0]
                  new_boundary2_p1 = boundary2_p1
                  new_boundary2_p2 = boundary2_p2
                ENDELSE

                boundary1_p1 = new_boundary1_p1
                boundary1_p2 = new_boundary1_p2
                boundary2_p1 = new_boundary2_p1
                boundary2_p2 = new_boundary2_p2

                ;print, 'wow!'
                ;stop    

                print,N_ELEMENTS(x)
                x = store_x
                y = store_y
   
                cont = 3
              ENDIF

            ENDELSE

          ENDIF

        ENDELSE

      ENDWHILE

      ; Cycle to the next segment if the current segment is not suitable:
      IF (cont EQ 2) THEN CONTINUE


      CASE contour_n OF
        217: BEGIN 
          xrange = [4.5,5.5] ; top
          yrange = [4.0,5.0] 
          END
        299: BEGIN 
          xrange = [5.0,7.0] ; top
          yrange = [3.6,4.8] 
          END
        ELSE: BEGIN
          xrange = user_xrange
          yrange = user_yrange
;          xrange = [ 3.8,4.2]  ; First intesection on CC
;          yrange = [ 1.0,2.0]
;          xrange = [4.5,5.5] ; top
;          yrange = [4.0,5.0] 
;          xrange = [ 4.5,5.5]
;          yrange = [ 4.0,5.5]
;          xrange = [ 3.8,4.2]
;          yrange = [-2.0,-1.0]
;          xrange = [ 7.0,8.0 ]
;          yrange = [ 3.0,4.0 ]
          END
      ENDCASE
      PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
      PLOT,wall.x,wall.y,color=Truecolor('Yellow'),  $
           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
      PLOT,x,y,color=Truecolor('White'),  $
           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,PSYM=3,/NOERASE

      ; Only take the parts of the contour segment that are inside the main wall, by dropping
      ; all points before and after the first and last wall intersections, respectively:       
      grid_ClipToWall, x, y, wall_pt1, wall_pt2, cycle, xrange, yrange

      print,'CYCLE:',cycle

    IF (1 EQ 1 AND iteration_count EQ 46 AND iseg EQ 0) THEN BEGIN
stop
      xrange = user_xrange
      yrange = user_yrange
;      xrange = [ 4.5,5.5]
;      yrange = [ 4.0,5.5]
;      xrange = [ 3.8,4.2]
;      yrange = [-2.0,-1.0]
         PLOT,x,y,color=Truecolor('White'), PSYM=3, $ 
              XRANGE=xrange,YRANGE=yrange,/NOERASE,XSTYLE=1,YSTYLE=1
         PLOT,[boundary1_p1[0],boundary1_p2[0]],  $
              [boundary1_p1[1],boundary1_p2[1]],  $
              color=Truecolor('Lightgreen'),  $
              XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
         PLOT,[boundary2_p1[0],boundary2_p2[0]],  $
              [boundary2_p1[1],boundary2_p2[1]],  $
              color=Truecolor('Green'),  $
              XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
print,'... nseg',nseg
      stop
    ENDIF






      ; Cycle if outside the core and no wall intersections were found:
      IF (cycle AND region NE CORE) THEN CONTINUE

      ; Check if the contour is radially outward of the secondary x-point,
      ; if that's the current region of interest (PROCESS_2ND = 2):
      print,'...',MEAN(x),b.x[b.null_i[2]],process_2nd
      IF (process_2nd EQ 2 AND MEAN(x) LT b.x[b.null_i[2]]) THEN BEGIN
        print,'...looking...........'
        CONTINUE 
      ENDIF

      seg_valid = 1
;
;     ------------------------------------------------------------------
;     LOOP TO FIND ??? ("TANGENCY") POINTS WITH THE WALL
;  
      cont      = 1

      WHILE (cont) DO BEGIN
        cont  = 0
        int_i = -1
        first_intersection = 1

        FOR i = 0, N_ELEMENTS(x)-2 DO BEGIN  
;         Find intersections with the wall segments:
          result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                     wall_pt1, wall_pt2, 1, status=status)

          IF (status) THEN BEGIN

            IF (N_ELEMENTS(result.i) GT 1) THEN BEGIN
              PLOT,[result.x],[result.y],color=Truecolor('White'), PSYM=6,  $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
              PRINT,'******** MULTIPLE SEGMENT INTERSECTIONS *********'
            ENDIF ELSE  $
              PLOT,[result.x],[result.y],color=Truecolor('Orange'), PSYM=6,  $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

            IF (first_intersection) THEN BEGIN
              first_intersection = 0  ; Can I replace this for some check on the nature of int_i?
              int_i    = i
              int_x    = result.x[0]
              int_y    = result.y[0]
            ENDIF ELSE BEGIN
              int_i    = [int_i,i]
              int_x    = [int_x,result.x[0]]
              int_y    = [int_y,result.y[0]]
            ENDELSE

          ENDIF

        ENDFOR 

        print,'      int_i:',int_i      
  
        int_n = N_ELEMENTS(int_i)

        ; Check if there are multiple wall intersections very close together, and 
        ; spatially refine the contour locally.  This creates overhead as the
        ; contour is passed through more than once, but I can't see a more general
        ; way of doing this:
        FOR j = 0, int_n-2 DO BEGIN
          IF (int_i[j+1]-int_i[j] LE 2) THEN BEGIN
            grid_RefineContour, x, y, int_i[j], xrange, yrange
            cont = 1  ; Loop through the contour again to refine intersection data
            print,'cycling...'
          ENDIF
        ENDFOR
      ENDWHILE ; Wall intersection loop
;
;     ------------------------------------------------------------------
;     PROCESS ??? ('tangecy points')
;  
      store_contour = 0

      CASE region OF
;       ----------------------------------------------------------------
        SOL: BEGIN
          IF (int_nlast[region] EQ -1 ) THEN BEGIN
            int_nlast[region] = int_n
            store_contour = 1
            sht_count = 0
          ENDIF ELSE BEGIN
            PRINT,'     psi_adjust 1:',psi_adjust,int_n,int_nlast[region]

            if (contour_n EQ 99) then begin
              print,'  boundary1:',boundary1_p1,boundary1_p2
              print,'  boundary2:',boundary2_p1,boundary2_p2
              print,'range:',xrange,yrange
              PLOT,[boundary1_p1[0],boundary1_p2[0]],  $
                   [boundary1_p1[1],boundary1_p2[1]],  $
                   color=Truecolor('Lightreen'),  $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
              PLOT,[boundary2_p1[0],boundary2_p2[0]],  $
                   [boundary2_p1[1],boundary2_p2[1]],  $
                   color=Truecolor('Lightgreen'),  $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

              PLOT,xrange,yrange,/NODATA,/NOERASE
              stop
            ENDIF
;            PRINT,'       check 1: ',int_nlast[region] LT int_n AND psi_adjust LT 0.0D
;            PRINT,'       check 2: ',int_nlast[region] GE int_n AND psi_adjust GT 0.0D
;            PRINT,'       check 3: ',psi_val LE psi_end AND int_n EQ 2

            IF (int_nlast[region] LT int_n AND psi_adjust*direction LT 0.0D) THEN BEGIN
              psi_adjust = -0.3D * psi_adjust
              sht_count = 0
            ENDIF ELSE BEGIN          
              IF (int_nlast[region] GE int_n AND psi_adjust*direction GT 0.0D) THEN BEGIN
                psi_adjust = -0.3D * psi_adjust
                sht_count = 0
              ENDIF ELSE BEGIN
                ; Exit conditions:
                IF (((direction EQ forward  AND psi_val LE psi_end)  OR  $
                     (direction EQ backward AND psi_val GE psi_end)) AND $
                    int_n EQ 2) THEN store_contour = 1
                IF (store_forced GT 0) THEN store_contour = 1
                length = grid_Length(x,y) 
;                print,'length:',length,sht_count
                IF (length LT 0.1D) THEN BEGIN
                  sht_count = sht_count + 1
                  IF (sht_count GE 3) THEN BEGIN
                    store_contour = 1
                  ENDIF
                ENDIF ELSE  sht_count = 0
                ; Store valid PSI from last standard progression step:
                IF (psi_adjust*direction LT 0.0D) THEN last_good_psi_val = psi_val  
              ENDELSE
            ENDELSE
          ENDELSE

          PRINT,'     psi_adjust 2:',psi_adjust,int_n,int_nlast[region]

          ; Search for tangency point complete, store the contour:
          broke_contour = 0
          IF (psi_adjust*direction GT 0.0D      AND  $
              psi_adjust*direction LT tol_adjust) THEN BEGIN
            store_contour  = 1
            broke_contour  = 1
          ENDIF

          IF (store_contour) THEN BEGIN
            n = N_ELEMENTS(int_i) - 1        
            dist1 = SQRT( (int_x[0  ] - int_x[1])^2 + (int_y[0  ] - int_y[1])^2 )
            dist2 = SQRT( (int_x[n-1] - int_x[n])^2 + (int_y[n-1] - int_y[n])^2 )
            print,'distances: ',dist1,dist2
            IF (dist1 LT 1.0D-5 OR dist2 LT 1.0D-5) THEN BEGIN
              print, 'funny business!'
 
              psi_start = psi_2nd_xpoint - 0.2D 
              psi_end   = psi_1st_xpoint 
              psi_val   = psi_val - 0.1D
              direction = BACKWARD
              psi_step  = -0.005D * direction
              psi_adjust = psi_step

              store_contour = 0
            ENDIF
          ENDIF

          IF (store_contour) THEN BEGIN

            CASE int_n OF 
              2: tangent_i = -1
              3: BEGIN
                ; Selecting the second point in the intersection list isn't perfect
                ; since you can get two intersections on the same line segment when
                ; resolving corners (I think), which results in two points being 
                ; counted, but they are almost on top of each other so it doesn't 
                ; really matter if you take int_i[1] or int_i[2]:
                tangent_i = int_i[1] + 1
                x = [x[0:int_i[1]],int_x[1],x[int_i[1]+1:N_ELEMENTS(x)-1]]
                y = [y[0:int_i[1]],int_y[1],y[int_i[1]+1:N_ELEMENTS(y)-1]]
                PLOT,[x[tangent_i]],[y[tangent_i]],color=Truecolor('Yellow'), PSYM=6,  $
                     XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
                END
              4: BEGIN
                 ; Take an estimate of the tangency point since a single intersection point (implying
	         ; that this is a 'real' tangecy point between the contour and the wall) was not found:
                tangent_i = int_i[1] + 1 + (int_i[2] - int_i[1]) / 2  ; OK?
                END
              ELSE: BEGIN
	   
                PLOT,[boundary1_p1[0],boundary1_p2[0]],  $
                     [boundary1_p1[1],boundary1_p2[1]],  $
                     color=Truecolor('Lightgreen'),  $
                     XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1
                PLOT,[boundary2_p1[0],boundary2_p2[0]],  $
                     [boundary2_p1[1],boundary2_p2[1]],  $
                     color=Truecolor('Green'),  $
                     XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
                PLOT,x,y,  $
                     color=Truecolor('Red'), PSYM=3, $
                     XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
                PLOT,wall.x,wall.y,color=Truecolor('Yellow'),  $
                     XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
	   
                PRINT,'ERROR grid_AnalyseBoundary: Tangency point ill-defined'
                STOP
                END
            ENDCASE

            ; Define radial boundary associated with the symmetry point:
            p1 = [0.0D,0.0D]
            p2 = p1
            IF (broke_contour) THEN BEGIN
              vx = -(y[tangent_i] - y[tangent_i-1])
              vy =   x[tangent_i] - x[tangent_i-1]
              IF (direction EQ backward) THEN BEGIN
                vx = -vx
                vy = -vy
              ENDIF 
              length = SQRT(vx^2 + vy^2)
              p1 = [x[tangent_i],y[tangent_i]]
              p2 = [p1[0] + 1.0D / length * vx,  $
                    p1[1] + 1.0D / length * vy]
              PLOT,[p1[0],p2[0]],[p1[1],p2[1]],color=Truecolor('Lightgreen'),  $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
            ENDIF

            print, '         <><><><>><><><>'
            print, 'storing contour...',sht_count,store_forced
            contour_n = contour_n + 1
            contour_data = {  $
              state        : 0                 ,  $
              region       : region            ,  $
              direction    : direction         ,  $
              psi          : psi_val           ,  $
              adjust       : psi_adjust        ,  $ 
              tangent_i    : tangent_i         ,  $
              tangent_p1   : p1                ,  $
              tangent_p2   : p2                ,  $
              boundary1_p1 : boundary1_p1      ,  $
              boundary1_p2 : boundary1_p2      ,  $
              boundary2_p1 : boundary2_p1      ,  $
              boundary2_p2 : boundary2_p2      ,  $
              int_n        : int_n             ,  $
              int_nlast    : int_nlast[region] ,  $
              int_i        : int_i             ,  $
              int_x        : int_x             ,  $
              int_y        : int_y             ,  $
              x            : x                 ,  $
              y            : y                 }

            name = 'contour' + STRING(contour_n,FORMAT='(I0)')
            IF (contour_n EQ 1) THEN  $
              contour_array = CREATE_STRUCT(              name,contour_data) ELSE  $
              contour_array = CREATE_STRUCT(contour_array,name,contour_data)

;            IF ((psi_adjust GE 0.0D AND psi_adjust LT tol_adjust) OR  $
            IF (broke_contour OR  $
                store_forced  OR  $
                (direction EQ forward  AND psi_val LE psi_end OR   $
                 direction EQ backward AND psi_val GE psi_end) OR  $
                sht_count GE 3) THEN BEGIN

              grid_SetActiveContour, contour_array, int_nlast, region, active_contour, store_forced,  $
                                     psi_start, psi_end, psi_val, psi_adjust, psi_step, sht_count, last_good_psi_val,   $
                                     boundary1_p1, boundary1_p2, boundary2_p1, boundary2_p2,  $
                                     the_search_continues, contour_n, direction, process_2nd,  $
                                     psi_1st_xpoint, psi_2nd_xpoint, user_step
              last_psi_val = -1.0D
              seg_active = MAKE_ARRAY(100,/LONG,VALUE=0L)
              save_boundary1_p1 = [-999.0D,-999.0D]
              save_boundary1_p2 = [-999.0D,-999.0D]
              save_boundary2_p1 = [-999.0D,-999.0D]
              save_boundary2_p2 = [-999.0D,-999.0D]

              print, 'active_contour', active_contour
              print,int_x
              print,int_y

;              if (0 EQ 1 AND direction EQ backward) then begin
              if (0 EQ 1 AND contour_n EQ 35) then begin
                xrange = user_xrange
                yrange = user_yrange
;                xrange = [ 3.8,4.10]  ; upper CC
;                yrange = [ 2.0, 3.0]
;                xrange = [ 3.8,4.10]  ; lower inner CC
;                yrange = [-2.6,-1.0]
;                xrange = [ 4.5,5.5]  ; top
;                yrange = [ 4.0,5.0]
;                xrange = [5.0,7.0] ; top outer
;                yrange = [3.6,4.8] 
                grid_Debug,contour_array,wall,xrange,yrange
              endif

            ENDIF
          ENDIF
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN
          PRINT,'ERROR grid_AnalyseBoundary: Unrecognized grid region'
          STOP
          END
      ENDCASE

;      PLOT,x,y,color=Truecolor('Orange'),  $
;           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,PSYM=3,/NOERASE
;      PLOT,[x[0]],[y[0]],color=Truecolor('Orange'),  $
;           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,PSYM=6,/NOERASE
;print,seg_valid
;         if (seg_active[iseg] EQ 1) THEN STOP

    ENDFOR  ; Contour segment loop

    PRINT,'seg_valid', seg_valid

    IF (NOT seg_valid) THEN BEGIN
      print,' '
      print,'  PROBLEM !', last_good_psi_val
      print,'           ', active_contour
      print,' '

;              xrange = user_xrange
;              yrange = user_yrange
;              PLOT,[boundary1_p1[0],boundary1_p2[0]],  $
;                   [boundary1_p1[1],boundary1_p2[1]],  $
;                   color=Truecolor('Lightgreen'),  $
;                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1
;              PLOT,[boundary2_p1[0],boundary2_p2[0]],  $
;                   [boundary2_p1[1],boundary2_p2[1]],  $
;                   color=Truecolor('Green'),  $
;                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
;              PLOT,x,y,  $
;                   color=Truecolor('Red'), PSYM=3, $
;                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
;       stop










      IF (active_contour EQ 0) THEN BEGIN
        PRINT,'ERROR'
        STOP
      ENDIF
      IF (last_good_psi_val EQ 0.0D) THEN BEGIN
        PRINT, 'ABANDONING ACTIVE CONTOUR'
        grid_SetActiveContour, contour_array, int_nlast, region, active_contour, store_forced,  $
                               psi_start, psi_end, psi_val, psi_adjust, psi_step, sht_count, last_good_psi_val,   $
                               boundary1_p1, boundary1_p2, boundary2_p1, boundary2_p2,  $
                               the_search_continues, contour_n, direction, process_2nd,  $
                               psi_1st_xpoint, psi_2nd_xpoint, user_step

        seg_active = MAKE_ARRAY(100,/LONG,VALUE=0L)
        save_boundary1_p1 = [-999.0D,-999.0D]
        save_boundary1_p2 = [-999.0D,-999.0D]
        save_boundary2_p1 = [-999.0D,-999.0D]
        save_boundary2_p2 = [-999.0D,-999.0D]
      ENDIF ELSE BEGIN
        PRINT, 'RECOVERING CONTOUR' 
        psi_val      = last_good_psi_val
        psi_adjust   = 0.0D 
        store_forced = 1

        xrange = user_xrange
        yrange = user_yrange
        PLOT,x,y,color=Truecolor('Cyan'), PSYM=3, $
             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
      ENDELSE

    ENDIF

    last_psi_val = psi_val
    cont = 1
    WHILE (cont) DO BEGIN
      cont = 0
      psi_val = psi_val + psi_adjust
      IF ((direction EQ forward  AND psi_val GT psi_start) OR  $
          (direction EQ backward AND psi_val LT psi_start)) THEN BEGIN
        PRINT,'--- PROTECTING PSI_VAL! ---------------'
        print,'psi_values: ',psi_start,psi_val,psi_end

        psi_val = psi_val - psi_adjust 
        IF (psi_val*direction GT psi_start) THEN BEGIN
          PRINT,'ERROR grid_AnalyseBoundary: psi value is out of whack'
          STOP
        ENDIF
        psi_adjust = psi_adjust * 0.3D
        cont = 1
      ENDIF
      print,'psi_values: ',psi_start,psi_val,psi_end
    ENDWHILE

    IF (process_2nd EQ 1 AND the_search_continues EQ 0) THEN BEGIN
      print, ' >>>>>>>>>>>>> 2nd x-point <<<<<<<<<<<< contour_n',contour_n

;      grid_Debug,contour_array,wall,xrange,yrange

      psi_start         = psi_2nd_xpoint - 0.0001D
      psi_val           = psi_start
      psi_adjust        = psi_step
      last_psi_val      = -1.0D
      last_good_psi_val =  0.0D

      the_search_continues = 1
      active_contour       = 0
      store_forced         = 0

      boundary1_p1 = [-999.0D,-999.0D]
      boundary1_p2 = [-999.0D,-999.0D]
      boundary2_p1 = [-999.0D,-999.0D]
      boundary2_p2 = [-999.0D,-999.0D]

      process_2nd = 2
    ENDIF

  ENDWHILE  ; PSIN loop


  print, '......all done......'
  xrange = user_xrange
  yrange = user_yrange
;  xrange = [ 3.8,4.10]  ; upper CC
;  yrange = [ 2.0, 3.0]
;  xrange = [ 3.8,4.10]  ; lower inner CC
;  yrange = [-2.6,-1.0]
;  xrange = [ 4.5,5.5]  ; top
;  yrange = [ 4.0,5.0]
;  xrange = [5.0,7.0] ; top outer
;  yrange = [3.6,4.8] 
;  xrange = [7.0,8.5] ; outer wide angle
;  yrange = [-2.0,4.0]
;  xrange = [3.8,4.4] ; inner
;  yrange = [-3.0,5.0]
;  xrange = [ 4.5,5.5]
;  yrange = [ 4.0,5.5]
;  xrange = [ 7.0,8.0 ]
;  yrange = [ 3.0,4.0 ]
  grid_Debug,contour_array,wall,xrange,yrange


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
      IF (b.psi_raw[i,j] GT max_psi) THEN BEGIN
        max_i = i
        max_j = j
        max_psi = b.psi_raw[i,j]
      ENDIF
    ENDFOR
  ENDFOR

  psi_diff = max_psi - b.psi_raw[min_i,min_j]
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
    OPLOT, [b_x[min_i[0]]], [b_y[min_j[0]]], PSYM=6, color=TrueColor('Orange')   ; o-point
    OPLOT, [b_x[min_i[1]]], [b_y[min_j[1]]], PSYM=6, color=TrueColor('Purple')   ; primary   x-point
    IF (N_ELEMENTS(min_i) EQ 3) THEN  $
      OPLOT, [b_x[min_i[2]]], [b_y[min_j[2]]], PSYM=6, color=TrueColor('Green' ) ; secondary x-point
  ENDIF

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
FUNCTION grid_Main, iter=iter, mast=mast, cmod=cmod

  IF (KEYWORD_SET(iter)) THEN machine = 'iter'
  IF (KEYWORD_SET(mast)) THEN machine = 'mast'
  IF (KEYWORD_SET(cmod)) THEN machine = 'cmod'

  IF (NOT KEYWORD_SET(machine)) THEN machine = 'iter'


; orange - o-point
; purple - primary x-point
; green  - secondary x-point

  
  CASE machine OF
    'iter': BEGIN
      xpoint_zone = [4.0,7.0,-3.8,5.5]
;      b = grid_ReadEQUFile('/home/ITER/lisgos/divimp/shots/iter/equilibria/feat_001.x4.equ');,/CHANGE_SIGN)
      b = grid_ReadEQUFile('/home/ITER/lisgos/divimp/shots/iter/1514/Baseline2008-li0.70.x4.equ');,/CHANGE_SIGN)
      b = grid_FindNullPoints (b,xpoint_zone,1,/debug)
      wall = grid_LoadWall('/home/ITER/lisgos/divimp/shots/iter/1514/psi_wall.dat',/debug)
      result = grid_AnalyseBoundary(b,wall,-0.005D,3.0D,/debug)
      ;shade_surf,b.psi_raw,b.x,b.y; ,ax=0  
      END
    
    'mast': BEGIN
       xpoint_zone = [0.4,1.0,-1.5,1.5]
      ; b = grid_ReadEQUFile('~/fuse_data/mast/shots/24867/24867_335.equ')
       b = grid_ReadEQUFile('~/divimp/shots/mast/24860/carre.24860_240.equ')
      ;shade_surf,b.psi_raw,b.x,b.y,ax=0
      ;stop
      b = grid_FindNullPoints (b,xpoint_zone,1,/debug)
      wall = grid_LoadWall('/home/ITER/lisgos/divimp/shots/mast/default/main_wall.dat', $
                           /debug,xrange=[0.05,2.2],yrange=[-2.5,2.5])
      b = grid_AnalyseBoundary(b,wall,-0.0005D,0.03D,xrange=[0.05,2.2],yrange=[-2.5,2.5],/debug)
      END

    'cmod': BEGIN

      wall = grid_LoadWall('/home/ITER/lisgos/divimp/shots/cmod/1100303017_0138/vessel_wall.dat',  $
                           /debug,xrange=[0.3,1.2],yrange=[-0.8,0.8])
       xpoint_zone = [0.52,0.75,-0.45,0.45]
       b = grid_ReadEQUFile('/home/ITER/lisgos/divimp/shots/'+  $
                            'cmod/1100303017_0138/cmod.1100303017.01380.x2.equ',/CHANGE_SIGN)
      b = grid_FindNullPoints (b,xpoint_zone,1,/debug)
    CONTOUR, b.psi, b.x, b.y, NLEVELS=20, c_labels=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
      b = grid_AnalyseBoundary(b,wall,-0.0005D,0.03D,  $
                           /debug,xrange=[0.3,1.2],yrange=[-0.8,0.8])

      END
     ELSE: BEGIN
       stop
       END
  ENDCASE



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

  RETURN, result
END