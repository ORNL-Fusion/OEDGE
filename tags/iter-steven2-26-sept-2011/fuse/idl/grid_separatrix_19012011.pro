;
; ======================================================================
;
PRO grid_SetupParameters
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
                 
  FORWARD       =  1.0D
  BACKWARD      = -1.0D
  CORE          = 1
  SOL           = 2
  SOL_LFS       = 3  
  SOL_HFS       = 4  
  PFZ           = 5  
  PFZ_SECONDARY = 6  

  HFS = 1
  LFS = 2

  LOWER_NULL    = 1
  UPPER_NULL    = 2
  ; ------------------------------------------------------------------
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
    PRINT, '  TAG = >'+STRTRIM(tag,2)+'<'
    PRINT, '  NAMES = ',TAG_NAMES(struct)
;    HELP, struct, /struct
    HELP, struct.contour12,/struct
    result = -1
  ENDELSE

  RETURN, result
END
;
; ======================================================================
;
FUNCTION grid_UpdateStructure, struct_array,tag,struct

  status = 0

  tags = STRUPCASE(TAG_NAMES(struct_array))

  FOR i = 0, N_ELEMENTS(tags)-1 DO BEGIN
    IF (tags[i] EQ STRUPCASE(tag)) THEN BEGIN
      val = struct
      status = 1
;      print, 'replacing! ',tag
;      print,val.state
    ENDIF ELSE  $
      val = grid_ExtractStructure(struct_array,tags[i])  

    IF (i EQ 0) THEN  $
      result = CREATE_STRUCT(       tags[i],val) ELSE  $
      result = CREATE_STRUCT(result,tags[i],val)
  ENDFOR

  IF (status EQ 0) THEN BEGIN
    PRINT, 'ERROR grid_UpdateStructure: Specified TAG not found'
    PRINT,'  TAG           =  ',STRUPCASE(tag)
    PRINT,'  STRUCTURE TAGS= ',tags    
    STOP
  ENDIF

  RETURN, result
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

  IF (ABS(a1[0]-a2[0]) LT 1.0D-6) THEN a2[0] = a2[0] + 0.0000001D
;  IF (ABS(a1[1]-a2[1]) LT 1.0D-6) THEN a2[1] = a2[1] + 0.00001D

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

;print, 'go'

;  min_t = 1000.0D

  FOR j = 0, N_ELEMENTS(seg_list)-1 DO BEGIN

    i = seg_list[j]        

    b1 = seg1[*,i]
    b2 = seg2[*,i]

    IF ( ABS(b1[0]-b2[0]) LT 1.0D-6 ) THEN b1[0] = b1[0] + 0.0000001D
;    IF ( ABS(b1[1]-b2[1]) LT 1.0D-6 ) THEN b1[1] = b1[1] + 0.00001D
  
    Lint, a1, a2, b1, b2, c, d, flag=flag
  
    IF (flag NE 1) THEN CONTINUE
  
    s = (c[0] - a1[0]) / (a2[0] - a1[0])    
    t = (c[0] - b1[0]) / (b2[0] - b1[0])    

;print,s,t,FORMAT='(2E20.10,I6)',KEYWORD_SET(int_i)
;    IF (ABS(t) LT min_t) THEN BEGIN
;      min_t = t
;      min_s = s
;    ENDIF

    IF ( s GE  0.0D    AND s LE 1.0D AND  $
         t GE -1.0D-10 AND t LT 1.0D) THEN BEGIN
;      print, 'ok'
      IF (KEYWORD_SET(int_x)) THEN BEGIN
        int_i    = [int_i,i]
        int_s    = [int_s,s]
        int_t    = [int_t,t]
        int_x    = [int_x,c[0]]
        int_y    = [int_y,c[1]]
      ENDIF ELSE BEGIN
        int_i    = [i]
        int_s    = [s]
        int_t    = [t]      
        int_x    = [c[0]]
        int_y    = [c[1]]
      ENDELSE
  
    ENDIF
  ENDFOR

; print,'min_s,t',min_s,min_t
  
;print,'set',KEYWORD_SET(int_x)

  IF (KEYWORD_SET(int_x)) THEN BEGIN
    status = 1  ; Can only be 0 or 1 since used as a logical variable below
    result = {          $
      i    : int_i   ,  $
      s    : int_s   ,  $
      t    : int_t   ,  $ 
      x    : int_x   ,  $
      y    : int_y   }          
  ENDIF ELSE BEGIN
;   No intersection found between the line and segment list:
;print,'not',NOT KEYWORD_SET(int_t)
    status =  0 ; Can only be 0 or 1 since used as a logical variable below
    result = -1
  ENDELSE

;help,result
;print,status

  RETURN,result

END
;
; ======================================================================
;
PRO grid_ZoneWall, wall_pt1, wall_pt2, nx=nx, ny=ny,  $
                   debug=debug, xrange=xrange, yrange=yrange

  COMMON zone, zone_list, zone_index, zone_nx  , zone_ny  ,  $
               zone_minx, zone_maxx , zone_miny, zone_maxy

  IF (NOT KEYWORD_SET(nx)) THEN nx = 10
  IF (NOT KEYWORD_SET(ny)) THEN ny = 10


;  wall_pt1 = wall.pt1
;  wall_pt2 = wall.pt2

  wall_n = N_ELEMENTS(wall_pt1[0,*])

; Setup zones:
  zone_nx = nx
  zone_ny = ny
  zone_list  = MAKE_ARRAY(3*wall_n,/LONG,VALUE=-1L)
  zone_index = MAKE_ARRAY(2,zone_nx,zone_ny,/LONG,VALUE=-1L)
  zone_minx = MIN([wall_pt1[0,*],wall_pt2[0,*]]) - 0.1D
  zone_maxx = MAX([wall_pt1[0,*],wall_pt2[0,*]]) + 0.1D
  zone_miny = MIN([wall_pt1[1,*],wall_pt2[1,*]]) - 0.1D
  zone_maxy = MAX([wall_pt1[1,*],wall_pt2[1,*]]) + 0.1D
  deltax = (zone_maxx - zone_minx) / DOUBLE(zone_nx)
  deltay = (zone_maxy - zone_miny) / DOUBLE(zone_ny)

  box_pt1 = MAKE_ARRAY(2,5,/DOUBLE,VALUE=0.0D)
  box_pt2 = box_pt1

  wall_check = MAKE_ARRAY(wall_n,/LONG,VALUE=0L)

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
;          print,count,i,j,k
          zone_list[count] = k
          IF (zone_index[0,i,j] EQ -1) THEN zone_index[0,i,j] = count
          zone_index[1,i,j] = count
          wall_check[k] = 1
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

  i = WHERE(wall_check EQ 0,count)
  IF (count GT 0) THEN BEGIN
    PRINT,'ERROR grid_ZoneWall: Wall segment not assigned to zone'
    STOP
  ENDIF



  IF (0 EQ 1) THEN BEGIN
    i = 0
    j = 1
    box_pt1[0,*] = [zone_minx+(i  )*deltax,zone_minx+(i  )*deltax,zone_minx+(i+1)*deltax,zone_minx+(i+1)*deltax,zone_minx+(i  )*deltax]
    box_pt1[1,*] = [zone_miny+(j  )*deltay,zone_miny+(j+1)*deltay,zone_miny+(j+1)*deltay,zone_miny+(j  )*deltay,zone_miny+(j  )*deltay]
    box_pt2[0,*] = [zone_minx+(i  )*deltax,zone_minx+(i+1)*deltax,zone_minx+(i+1)*deltax,zone_minx+(i  )*deltax,zone_minx+(i  )*deltax]
    box_pt2[1,*] = [zone_miny+(j+1)*deltay,zone_minx+(j+1)*deltax,zone_miny+(j  )*deltay,zone_miny+(j  )*deltay,zone_miny+(j+1)*deltay]

    OPLOT,box_pt1[0,*],box_pt1[1,*],color=Truecolor('Lightgreen')


    OPLOT,wall.x,wall.y,color=Truecolor('Yellow'), PSYM=6
    print,zone_index[0:1,i,j]
    stop
  ENDIF

END
;
; ======================================================================
;
FUNCTION grid_LoadWall, file, debug=debug, xrange=xrange, yrange=yrange

;  wall = cortex_LoadAnnotationData(1,'/home/ITER/lisgos/divimp/shots/iter/1514/psi_wall.dat') 
;  wall = cortex_LoadAnnotationData(1,'/home/ITER/lisgos/divimp/shots/mast/default/main_wall.dat') ; MAST

  wall = cortex_LoadAnnotationData(1,file)


; Check for redundant points:
  x = wall.x
  y = wall.y
  wall_n = N_ELEMENTS(x)
  FOR i = wall_n-2, 0, -1 DO BEGIN
    dist = SQRT( (x[i]-x[i+1])^2 + (y[i]-y[i+1])^2 )    
    IF (dist LT 1.0D-06) THEN BEGIN
      print, 'deleting wall segment',i
      n = N_ELEMENTS(x)
      IF (i EQ wall_n-2) THEN BEGIN
        x = x[0:i]
        y = y[0:i]
      ENDIF ELSE BEGIN
        x = [x[0:i],x[i+2:n-1]]
        y = [y[0:i],y[i+2:n-1]]
      ENDELSE
    ENDIF
  ENDFOR
  wall = grid_UpdateStructure(wall,'x',x)
  wall = grid_UpdateStructure(wall,'y',y)

  n = N_ELEMENTS(wall.x)

  IF (wall.x[0] NE wall.x[n-1] OR  $
      wall.y[0] NE wall.y[n-1]) THEN BEGIN
    wall = { x : [wall.x,wall.x[0]], y : [wall.y,wall.y[0]] }
    n++
  ENDIF

  nwall = 1  ; In future, more than one group of wall segments may be loaded, so 
             ; making some prevision here.  For inside-the-vessel checks will need to isolate
             ; the main wall segments...
  pti = MAKE_ARRAY(N_ELEMENTS(wall.x)-1,/LONG,VALUE= 1)

  ptc = MAKE_ARRAY(N_ELEMENTS(wall.x)-1,/LONG,VALUE=-1)
  ptt = MAKE_ARRAY(N_ELEMENTS(wall.x)-1,/LONG,VALUE=-1)

  wall = CREATE_STRUCT(wall, 'n',nwall, 'i',pti)

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


  wall = CREATE_STRUCT(wall,'pt1',pt1,'pt2',pt2,'pti',pti,'ptc',ptc,'ptt',ptt)


  grid_ZoneWall, pt1, pt2, debug=debug, xrange=xrange, yrange=yrange


;  n=N_ELEMENTS(wall.pti)
;  FOR j = 0, n-1 DO BEGIN
;    print,j,wall.pti[j],wall.pt1[0,j],wall.pt1[1,j],wall.pt2[0,j],wall.pt2[1,j]
;  ENDFOR
;  stop

;  result = {  $
;    x   : wall.x ,  $
;    y   : wall.y ,  $
;    pt1 : pt1    ,  $
;    pt2 : pt2    }

  result = wall

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
;  ======================================================================
;
PRO grid_Debug, b, contour_array, wall, xrange, yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------

  PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1

  OPLOT,wall.x,wall.y,color=Truecolor('Yellow')

  OPLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Orange')
  OPLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Orange'), PSYM=7

;  i = WHERE(wall.ptt EQ 2, count) 
;  IF (count GT 0) THEN  $
;    PLOT,wall.pt1[0,i],wall.pt1[1,i],color=Truecolor('Pink'), PSYM=7,  $
;         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

  i = WHERE(wall.ptc NE -1, count) 
  IF (count GT 0) THEN  $
    OPLOT,wall.pt1[0,i],wall.pt1[1,i],color=Truecolor('White'),PSYM=6

  i = WHERE(wall.ptt EQ 1, count) 
  IF (count GT 0) THEN  $
    OPLOT,wall.pt1[0,i],wall.pt1[1,i],color=Truecolor('Pink'), PSYM=7

  i = WHERE(wall.ptt EQ 2, count) 
  IF (count GT 0) THEN  $
    OPLOT,wall.pt1[0,i],wall.pt1[1,i],color=Truecolor('Magenta'), PSYM=7

  i = WHERE(wall.ptt EQ 3, count) 
  IF (count GT 0) THEN  $
    OPLOT,wall.pt1[0,i],wall.pt1[1,i],color=Truecolor('Lightgreen'),PSYM=7




;  help,contour_array,/struct

  tags = STRUPCASE(TAG_NAMES(contour_array))

  FOR i = 1, N_ELEMENTS(tags) DO BEGIN

    name = 'contour' + STRING(i,FORMAT='(I0)')

;    print,name,i, N_ELEMENTS(tags)

    cnt = grid_ExtractStructure(contour_array,name)  

    XYOUTS, [cnt.x[0]], [cnt.y[0]], CHARSIZE=1.0, ALIGNMENT=0.5, $
            STRTRIM(STRING(i),2), COLOR=TrueColor('White')

    CASE cnt.region OF
      SOL: color = 'Red'
      PFZ: color = 'Orange'
      ELSE: stop
    ENDCASE
  
    OPLOT,cnt.x,cnt.y,color=Truecolor(color)

    OPLOT,cnt.x,cnt.y,color=Truecolor('White'), PSYM=3

;    print,i,cnt.tangent_i

    IF (cnt.tangent_i NE -1) THEN  $
      OPLOT,[cnt.x[cnt.tangent_i]],[cnt.y[cnt.tangent_i]],  $
            color=Truecolor('Lightgreen'), PSYM=6

    IF (cnt.origin EQ 1 AND i EQ N_ELEMENTS(tags)) THEN BEGIN

      OPLOT,[cnt.boundary1_p1[0],cnt.boundary1_p2[0]],  $
            [cnt.boundary1_p1[1],cnt.boundary1_p2[1]],  $
            color=Truecolor('Lightgreen')

      OPLOT,[cnt.boundary2_p1[0],cnt.boundary2_p2[0]],  $
            [cnt.boundary2_p1[1],cnt.boundary2_p2[1]],  $
            color=Truecolor('Lightgreen')

    ENDIF


  ENDFOR



  OPLOT, [b.x[b.null_i[0]]], [b.y[b.null_j[0]]], PSYM=7, color=TrueColor('Red')      ; o-point
  OPLOT, [b.x[b.null_i[1]]], [b.y[b.null_j[1]]], PSYM=7, color=TrueColor('Green')   ; primary   x-point
  IF (N_ELEMENTS(b.null_i) EQ 3) THEN  $
    OPLOT, [b.x[b.null_i[2]]], [b.y[b.null_j[2]]], PSYM=7, color=TrueColor('Blue' ) ; secondary x-point

;  print,'saving contour_array'
;  SAVE,filename='contour_array.sav',contour_array


  print,'done in debug',n_elements(tags)

  stop

END
;
;
; ======================================================================
;
; From http://www.dfanning.com/tips/point_in_polygon.html.
;
FUNCTION grid_PointInPolygon, x, y, px, py, id=id

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

;    CONTOUR, psi, psi_x, psi_y, NLEVELS=20, c_labels=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
;stop

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
                           psi_1st_xpoint, psi_2nd_xpoint, user_step, user_finish
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
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
            IF (region EQ SOL AND direction EQ BACKWARD) THEN BEGIN
              IF (process_2nd EQ 2) THEN psi_end = psi_2nd_xpoint ELSE $
                                         psi_end = psi_1st_xpoint 
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
            IF (region EQ SOL AND direction EQ BACKWARD) THEN BEGIN
              direction = FORWARD
              psi_end = psi_1st_xpoint - user_finish
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
          contour_array = grid_UpdateStructure(contour_array,tag,cnt)
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

;  IF (region EQ SOL AND direction EQ BACKWARD) THEN   $
;    psi_step = user_step * direction ELSE  $
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
PRO grid_RefineContour, x, y, i, n_pts=n_pts,  $
                        debug=debug, xrange=xrange, yrange=yrange

  IF (NOT KEYWORD_SET(n_pts)) THEN n_pts = 10

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
      grid_RefineContour, x, y, i,  $
                          debug=debug, xrange=xrange, yrange=yrange
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
      grid_RefineContour, x, y, i,  $
                          debug=debug, xrange=xrange, yrange=yrange
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
FUNCTION grid_CalcLength, x, y, istart=istart, iend=iend

  result = 0.0D

  IF (NOT KEYWORD_SET(istart)) THEN istart = 0
  IF (NOT KEYWORD_SET(iend  )) THEN iend   = N_ELEMENTS(x)-2

  FOR i = istart, iend DO  $
    result = result + SQRT ( (x[i]-x[i+1])^2 + (y[i] - y[i+1])^2 )

  RETURN, result
END
;
; ======================================================================
;
FUNCTION grid_CheckSecondary, b, wall, geometry, psi_val, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------ 

  print, 'redundant'
  STOP

  result = 0

help,b,/struct

  IF (NOT KEYWORD_SET(user_xrange)) THEN user_xrange = [ 3.0,9.0]
  IF (NOT KEYWORD_SET(user_yrange)) THEN user_yrange = [-6.0,6.0]

  psi   = b.psi 
  psi_x = b.x
  psi_y = b.y

  ctr = grid_ExtractContour(psi, psi_x, psi_y, psi_val)

  ibrk = WHERE(ctr.dist GT 0.10)
  nseg = N_ELEMENTS(ibrk) 
  ibrk = [-1,ibrk,ctr.n-1]
  
  ; Identify the contour assoicaited with the secondary PFZ:
  length_max = 0.0D
  FOR iseg = 0, nseg DO BEGIN
    IF ((ibrk[iseg+1]-ibrk[iseg]) LT 5) THEN CONTINUE
    x = ctr.x[ibrk[iseg]+1:ibrk[iseg+1]]
    y = ctr.y[ibrk[iseg]+1:ibrk[iseg+1]]
    mean_y = MEAN(y)
    IF ((geometry EQ LOWER_NULL AND mean_y LT b.y[b.null_j[2]]) OR  $
        (geometry EQ UPPER_NULL AND mean_y GT b.y[b.null_j[2]])) THEN CONTINUE
    length = grid_CalcLength(x,y)
    IF (length GT length_max) THEN BEGIN
      length_max = length
      save_x = x
      save_y = y
    ENDIF
  ENDFOR

  x = save_x
  y = save_y
  inside = grid_PointInPolygon(x,y,wall.x,wall.y)
  i = WHERE(inside EQ 1, count)
  IF (count GT 0) THEN result = 1 


  IF (KEYWORD_SET(debug)) THEN BEGIN
    print,inside
    print,count
    print,result

    PLOT,wall.x,wall.y,color=Truecolor('Yellow'),  $
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
    PLOT,x,y,color=Truecolor('Green'),  $
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1
  ENDIF

  RETURN,result
END
;
; ======================================================================
;
; How do I avoid the distorted target cells that plagued the JET grid, i.e. the 
; cells overlapped?  Just keep the resolution high? 
;
; Need to bring in the other wall segments now I think, and look for intersections, rather 
; than checking if a point is inside or outside (although good to have some appropriate
; spatial refinement before this step? to make sure the entrances and exits are all accounted for?)
;
;
FUNCTION grid_AnalyseBoundary, b, wall, user_step, user_finish, machine, $
                               xrange=user_xrange, yrange=user_yrange, debug=debug, save=save
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------

  IF (NOT KEYWORD_SET(user_xrange)) THEN user_xrange = [ 3.0,9.0]
  IF (NOT KEYWORD_SET(user_yrange)) THEN user_yrange = [-6.0,6.0]


  region = SOL

  direction = FORWARD


  psi   = b.psi ; psi_raw
  psi_x = b.x
  psi_y = b.y

  ; Set lower or upper null:
  IF (b.y[b.null_j[1]] LT b.y[b.null_j[0]]) THEN  $
    geometry = LOWER_NULL ELSE geometry = UPPER_NULL

  psi_1st_xpoint = b.psi_1st_xpoint


  psi_val    = DOUBLE(psi_1st_xpoint) ; - 0.00001D ; MAST
  psi_start  = psi_val
  psi_end    = psi_1st_xpoint - user_finish
  psi_step   = user_step * direction
  psi_adjust = psi_step
  tol_adjust = 1.0D-7  ; 7 ; -6

  ; Decide what to do about the secondary x-point:
  process_2nd     = -1
  process_2nd_pfz =  0
  psi_span       = 0.0D0
  psi_2nd_xpoint = -999.0D
  IF (N_ELEMENTS(b.null_i) GE 3) THEN BEGIN
    ; Check whether or not the non-primary x-points are inside the vessel wall:
    psi_2nd_xpoint = b.psi_2nd_xpoint
    psi_span       = ABS(psi_1st_xpoint - psi_2nd_xpoint)
    IF (psi_2nd_xpoint GT psi_end) THEN BEGIN
      ; If the secondary x-point is inside the main wall, flag that processing
      ; is required:
      inside = grid_PointInPolygon(b.x[b.null_i[2]],b.y[b.null_j[2]],wall.x,wall.y)
      i = WHERE(inside EQ 1, count)
      IF (count GT 0) THEN BEGIN
        process_2nd = 0
        ; Decide if there's a well defined secondary PFR, i.e. see if the 2nd x-point is
        ; right on the wall:
        process_2nd_pfz = 1
;        process_2nd_pfz = grid_CheckSecondary(b,wall,geometry,psi_2nd_xpoint-user_step*0.01D,  $
;                                              debug=debug,xrange=xrange,yrange=yrange)
      ENDIF
    ENDIF 
  ENDIF



  ; Setup local wall segment arrays:
  wall_pt1 = wall.pt1
  wall_pt2 = wall.pt2

;  wall = CREATE_STRUCT(wall,'wall_pt1',wall_pt1,'wall_pt2',wall_pt2)



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

  active_2nd = 0
  process_pfz = 0
  failure_2nd_pfz = 0
;
; ----------------------------------------------------------------------
; MAIN LOOP
;  
  WHILE (the_search_continues) DO BEGIN

    iteration_count = iteration_count + 1

    IF (0 EQ 1 AND contour_n EQ 0) THEN BEGIN
      RESTORE,'contour_array.sav'
      region = SOL
      active_contour = 20
      contour_n = 20
      grid_SetActiveContour, contour_array, int_nlast, region, active_contour, store_forced,  $
                             psi_start, psi_end, psi_val, psi_adjust, psi_step, sht_count, last_good_psi_val,   $
                             boundary1_p1, boundary1_p2, boundary2_p1, boundary2_p2,  $
                             the_search_continues, contour_n, direction, process_2nd,  $
                             psi_1st_xpoint, psi_2nd_xpoint, user_step, user_finish
      the_search_continues = 1
;      psi_start = psi_2nd_xpoint - 0.1D 
;      psi_end   = psi_1st_xpoint 
;      psi_val   = psi_start
;      direction = BACKWARD
;      psi_step  = -0.005D * direction

      psi_start = psi_2nd_xpoint - 0.001D
      psi_val = psi_start
      active_2nd = HFS
      process_2nd = 2
      int_nlast[SOL] = 2


;      ; Start secondary PFZ processing:
;      process_2nd = 3
;      psi_start         = psi_2nd_xpoint + 0.000001D
;      psi_end           = psi_start + user_finish * 0.3D  ; replace with specific parameter
;      psi_val           = psi_start
;      psi_adjust        = psi_step
;      last_psi_val      = -1.0D
;      last_good_psi_val =  0.0D
;      the_search_continues = 1
;      active_contour       = 0
;      store_forced         = 0
;      direction = BACKWARD
;      boundary1_p1 = [-999.0D,-999.0D]
;      boundary1_p2 = [-999.0D,-999.0D]
;      boundary2_p1 = [-999.0D,-999.0D]
;      boundary2_p2 = [-999.0D,-999.0D]
;      process_pfz = 2
;      region = PFZ
;      int_nlast[region] = 2  ; Rubbish, remove [region] in the code...


;      ; primary pfr processing
;      process_2nd = 3
;      psi_start         = psi_1st_xpoint + user_finish * 0.29D
;      psi_end           = psi_start + user_finish * 0.01D ; replace 0.3 with specific pfz parameter
;      psi_val           = psi_start
;      psi_adjust        = psi_step
;      last_psi_val      = -1.0D
;      last_good_psi_val =  0.0D
;      the_search_continues = 1
;      active_contour       = 0
;      store_forced         = 0
;      direction = BACKWARD
;      boundary1_p1 = [-999.0D,-999.0D]
;      boundary1_p2 = [-999.0D,-999.0D]
;      boundary2_p1 = [-999.0D,-999.0D]
;      boundary2_p2 = [-999.0D,-999.0D]
;      process_pfz = 1
;      region = PFZ
;      int_nlast[region] = 2  ; Rubbish, remove [region] in the code...

    ENDIF

    print,'===========',contour_n,iteration_count,'============='
    print,'psi_values: ',psi_start,psi_val,psi_end,process_2nd, direction

; help,b,/struct

    IF (process_2nd EQ 0 AND psi_val LT psi_2nd_xpoint) THEN process_2nd = 1
      

    ctr = grid_ExtractContour(psi, psi_x, psi_y, psi_val)

    ibrk = WHERE(ctr.dist GT 0.10)
    nseg = N_ELEMENTS(ibrk) 
    ibrk = [-1,ibrk,ctr.n-1]

    IF (0 EQ 1 AND process_2nd EQ 1) THEN BEGIN

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
        IF (seg_active[iseg] EQ 1 AND  $
            psi_val LT psi_2nd_xpoint) THEN BEGIN
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

      ; Not a real segment so skip it, may need to strengthen this check:
      IF (ibrk[iseg]+2 EQ ibrk[iseg+1] OR  $
          ibrk[iseg]+1 EQ ibrk[iseg+1] OR  $
          ibrk[iseg]   EQ ibrk[iseg+1]) THEN CONTINUE  

      x = ctr.x[ibrk[iseg]+1:ibrk[iseg+1]]
      y = ctr.y[ibrk[iseg]+1:ibrk[iseg+1]]


      IF (0 EQ 1 AND iteration_count EQ 1 AND iseg EQ 2) THEN BEGIN
      
        print,ibrk[iseg],ibrk[iseg+1]
        print,ibrk[iseg]+1 EQ ibrk[iseg+1]
        print,ibrk[iseg]   EQ ibrk[iseg+1]

        xrange = user_xrange
        yrange = user_yrange
;        xrange = [ 4.5,5.5]
;        yrange = [ 4.0,5.5]
;        xrange = [ 3.8,4.2]
;        yrange = [-2.0,-1.0]
        PLOT,x,y,color=Truecolor('White'), PSYM=3, $ 
             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1
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

        print,'psi:',psi_val
        print,'1st',psi_1st_xpoint
        print,'2nd',psi_2nd_xpoint

        ;print,x
        ;print,y
        stop
      ENDIF

      ; Isolate region of interest along the contour by clipping at the poloidal 
      ; boundaries:

      caution = 0
      store_x = x
      store_y = y
      cont    = 1
      first_pass = 1
      WHILE (cont EQ 1 OR cont EQ 3) DO BEGIN
        cont = 0

        print,'  boundary1:',FLOAT(boundary1_p1),FLOAT(boundary1_p2),FORMAT='(A,4F10.4)'
        print,'  boundary2:',FLOAT(boundary2_p1),FLOAT(boundary2_p2),FORMAT='(A,4F10.4)'
              print,'psi...',psi_val,psi_2nd_xpoint,psi_val GT psi_2nd_xpoint
        print,N_ELEMENTS(x)

        failure = 0 

        IF (boundary1_p1[0] NE -999.0D) THEN BEGIN
          ; Find intersections with first boundary:
          FOR i = N_ELEMENTS(x)-2, 0, -1 DO BEGIN
;          FOR i = N_ELEMENTS(x)-2, 1, -1 DO BEGIN  -changed 14/12/2010 -- not sure why this was here...
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
            IF (caution EQ 0) THEN  $  ; *** CAN REMOVE THIS CAUTION REFERENCE? ***
              failure = failure + 2  $
            ELSE BEGIN
              boundary1_p1 = save_boundary1_p1
              boundary1_p2 = save_boundary1_p2
              save_boundary1_p1[0] = -999.0D0 ; Turn off boundary updating
              failure = 3
            ENDELSE                  
          ENDIF ELSE BEGIN
            IF (N_ELEMENTS(result.i) EQ 1) THEN BEGIN
              x = x[0:i+1]
              y = y[0:i+1]
              x[i+1] = result.x[0]
              y[i+1] = result.y[0]
            ENDIF ELSE BEGIN
               print,'problem here as well'
               stop
            ENDELSE
          ENDELSE
        ENDIF
        
        print, 'failure', failure


        psi_check = 0
        IF (first_pass AND (process_2nd EQ -1 AND psi_span GT 0.0D0)) THEN BEGIN
print, 'cont check:', cont
          psi_delta = ABS(psi_val - psi_2nd_xpoint)
          psi_check = (psi_delta LT 0.05D * psi_span) AND (psi_val LT psi_2nd_xpoint)
        ENDIF

        IF ((failure EQ 1 AND boundary2_p1[0] EQ -999.0D AND NOT psi_check) OR  $
            (failure EQ 2 AND boundary1_p1[0] EQ -999.0D AND NOT psi_check) OR  $
             failure EQ 3) THEN  $
          cont = 2  $
        ELSE BEGIN
          IF (failure EQ 1 OR failure EQ 2) THEN BEGIN

            FOR i = 0, N_ELEMENTS(x)-2 DO BEGIN
              result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                         wall_pt1, wall_pt2, 1, status=status)
              IF (status) THEN BREAK
            ENDFOR

            print,'checking, inside grid...?',status

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


              IF (cont EQ 3) THEN stop

              IF (process_2nd EQ -1 OR  $
                  (process_2nd GE 0 AND psi_val      LT psi_2nd_xpoint AND  $
                                        last_psi_val GT psi_2nd_xpoint)) THEN BEGIN

                seg_active[iseg] = 1

                save_boundary1_p1 = boundary1_p1
                save_boundary1_p2 = boundary1_p2
                save_boundary2_p1 = boundary2_p1
                save_boundary2_p2 = boundary2_p2

                IF (failure EQ 2) THEN BEGIN
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

                x = store_x
                y = store_y
   
                cont = 3
              ENDIF

            ENDELSE

          ENDIF

        ENDELSE

        first_pass = 0

      ENDWHILE

      ; Cycle to the next segment if the current segment is not suitable:
      IF (cont EQ 2) THEN BEGIN
        IF (seg_active[iseg] EQ 1) THEN seg_active[iseg] = 0
        CONTINUE
      ENDIF 


      CASE contour_n OF
        221: BEGIN 
          xrange = [4.0,5.0] ; top, in a bit
          yrange = [4.0,5.0] 
          END
        220: BEGIN 
          xrange = [4.5,5.5] ; top
          yrange = [4.0,5.0] 
          END
        220: BEGIN 
          xrange = [4.5,5.5] ; top
          yrange = [4.0,5.0] 
          END
        221: BEGIN 
          xrange = [4.5,5.5] ; top
          yrange = [4.0,5.0] 
          END
        222: BEGIN 
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
      OPLOT, [b.x[b.null_i[0]]], [b.y[b.null_j[0]]], PSYM=6, color=TrueColor('Red')      ; o-point
      OPLOT, [b.x[b.null_i[1]]], [b.y[b.null_j[1]]], PSYM=6, color=TrueColor('Green')   ; primary   x-point
      IF (N_ELEMENTS(b.null_i) EQ 3) THEN  $
        OPLOT, [b.x[b.null_i[2]]], [b.y[b.null_j[2]]], PSYM=6, color=TrueColor('Blue' ) ; secondary x-point
      PLOT,[boundary1_p1[0],boundary1_p2[0]],  $
           [boundary1_p1[1],boundary1_p2[1]],  $
           color=Truecolor('Lightgreen'),  $
           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
      PLOT,[boundary2_p1[0],boundary2_p2[0]],  $
           [boundary2_p1[1],boundary2_p2[1]],  $
           color=Truecolor('Green'),  $
           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE


      ; Only take the parts of the contour segment that are inside the main wall, by dropping
      ; all points before and after the first and last wall intersections, respectively:       
      grid_ClipToWall, x, y, wall_pt1, wall_pt2, cycle, xrange, yrange

      print,'CYCLE:',cycle

    IF (0 EQ 1 AND iteration_count EQ 27 AND iseg EQ 1) THEN BEGIN
      
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

      ; Cycle if no wall interesections were found:
      IF (cycle) THEN CONTINUE

      ; Check if the contour is radially outward of the secondary x-point,
      ; if that's the current region of interest (PROCESS_2ND = 2):
      IF ( process_2nd EQ 2 ) THEN BEGIN
        print,'...',MEAN(x),b.x[b.null_i[2]],process_2nd
        mean_x = MEAN(x)
        IF ((active_2nd EQ HFS AND mean_x LT b.x[b.null_i[2]]) OR  $
            (active_2nd EQ LFS AND mean_x GT b.x[b.null_i[2]])) THEN BEGIN
          print,'...looking...........'
          CONTINUE 
        ENDIF
      ENDIF

      IF ( process_pfz EQ 1 ) THEN BEGIN
        mean_y = MEAN(y)
        IF ((geometry EQ LOWER_NULL AND mean_y GT b.y[b.null_j[1]]) OR  $
            (geometry EQ UPPER_NULL AND mean_y LT b.y[b.null_j[1]])) THEN BEGIN
          print,'...looking 1st pfr...........'
          CONTINUE 
        ENDIF
      ENDIF

      IF ( process_pfz EQ 2 ) THEN BEGIN
        mean_y = MEAN(y)
        IF ((geometry EQ LOWER_NULL AND mean_y LT b.y[b.null_j[2]]) OR  $
            (geometry EQ UPPER_NULL AND mean_y GT b.y[b.null_j[2]])) THEN BEGIN
          print,'...looking 2nd pfr...........'
          CONTINUE 
        ENDIF
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
        ; contour is passed through more than once, but I can't see a more efficent
        ; way of making this general check:
        FOR j = 0, int_n-2 DO BEGIN
          IF (int_i[j+1]-int_i[j] LE 2) THEN BEGIN
            grid_RefineContour, x, y, int_i[j],  $
                                debug=debug, xrange=xrange, yrange=yrange
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
      cycle = 0

;      CASE region OF
;       ----------------------------------------------------------------
;        SOL: BEGIN



          IF (int_nlast[region] EQ -1 ) THEN BEGIN
            int_nlast[region] = int_n
            store_contour = 1
            sht_count = 0
          ENDIF ELSE BEGIN
            PRINT,'     psi_adjust 1:',psi_adjust,int_n,int_nlast[region]

            if (contour_n EQ 999) then begin
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
                IF (((direction EQ FORWARD  AND psi_val LE psi_end)  OR  $
                     (direction EQ BACKWARD AND psi_val GE psi_end)) AND $
                    int_n EQ 2) THEN store_contour = 1
                IF (store_forced GT 0) THEN store_contour = 1
                length = grid_Length(x,y) 
;                print,'length:',length,sht_count
                IF (length LT 0.1D) THEN BEGIN
                  sht_count = sht_count + 1
                  IF (sht_count GE 3) THEN store_contour = 1
                ENDIF ELSE  sht_count = 0

                IF (int_n EQ 1 AND int_nlast[region] EQ 2) THEN BEGIN
                  IF (last_good_psi_val NE 0.0D) THEN BEGIN
                    PRINT, 'STRANGE EVENT, PERHAPS OUT OF SPACE, TAKING LAST GOOD CONTOUR'
                    psi_val      = last_good_psi_val
                    psi_adjust   = 0.0D 
                    store_forced = 1

                    xrange = user_xrange
                    yrange = user_yrange
                    PLOT,x,y,color=Truecolor('Cyan'), PSYM=3, $
                         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
                      print,x,y
                    CONTINUE
                  ENDIF ELSE BEGIN
                    PRINT, 'NOT SURE WHAT TO DO'
                    PLOT,[boundary1_p1[0],boundary1_p2[0]],  $
                         [boundary1_p1[1],boundary1_p2[1]],  $
                         color=Truecolor('Lightgreen'), $
                         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
                    PLOT,[boundary2_p1[0],boundary2_p2[0]],  $
                         [boundary2_p1[1],boundary2_p2[1]],  $
                         color=Truecolor('Lightgreen'), $
                         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
                    STOP
                  ENDELSE
                ENDIF

                ; Store valid PSI from last standard progression step:
                IF (psi_adjust*direction LT 0.0D) THEN last_good_psi_val = psi_val  

              ENDELSE
            ENDELSE
          ENDELSE

;          IF (cycle) THEN CONTINUE


          PRINT,'     psi_adjust 2:',psi_adjust,int_n,int_nlast[region],store_contour

          ; Search for tangency point complete, store the contour:
          broke_contour = 0
          IF (psi_adjust*direction GT 0.0D      AND  $
              psi_adjust*direction LT tol_adjust) THEN BEGIN
            store_contour  = 1
            broke_contour  = 1
          ENDIF

          ; Some additional checks for strange conditions:
          reverse_direction = 0
 
          IF (store_contour) THEN BEGIN
            n = N_ELEMENTS(int_i) - 1        
            dist1 = SQRT( (int_x[0  ]-int_x[1])^2 + (int_y[0  ]-int_y[1])^2 )
            dist2 = SQRT( (int_x[n-1]-int_x[n])^2 + (int_y[n-1]-int_y[n])^2 )
            print,'distances: ',dist1,dist2
;            stop
            IF (dist1 LT 1.0D-3 OR dist2 LT 1.0D-3) THEN BEGIN
              print, ' '
              print, ' *********** funny business! ************** ' 
              print, ' '
              reverse_direction = 1
            ENDIF
          ENDIF

          IF (store_contour AND (NOT reverse_direction) AND int_n EQ 3) THEN BEGIN
            inside = grid_PointInPolygon(x,y,wall.x,wall.y)
            dummy = WHERE(inside[int_i[0]+1:int_i[2]-1] EQ 0, count)
            print,'counting:',count
            PLOT,x[int_i[0]+1:int_i[1]-1],y[int_i[0]+1:int_i[1]-1],  $
                 color=Truecolor('Pink'), PSYM=7, $
                 XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
            PLOT,x[int_i[1]+1:int_i[2]-1],y[int_i[1]+1:int_i[2]-1],  $
                 color=Truecolor('Pink'), PSYM=7, $
                 XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
            IF (count GT 0) THEN BEGIN
              print, ' '
              print, ' *********** funny business - counting! ************** ' 
              print, ' '
              reverse_direction = 1
            ENDIF
          ENDIF

          IF (reverse_direction) THEN BEGIN
            IF (region EQ SOL) THEN BEGIN
              psi_start = psi_2nd_xpoint - 0.2D 
              psi_end   = psi_1st_xpoint 
              psi_val   = psi_val - 0.1D
              direction = BACKWARD
              psi_step  = -0.005D * direction
              psi_adjust = psi_step
            ENDIF ELSE BEGIN
              print, 'pfz not done yet'
              STOP
            ENDELSE  
            store_contour = 0
          ENDIF



          IF (store_contour) THEN BEGIN

            CASE int_n OF 
              2: tangent_i = -1
              3: BEGIN
                ; Selecting the second point in the intersection list isn't perfect
                ; since you can get two intersections on the same line segment when
                ; resolving corners (I think), which results in two points being 
                ; counted, but they are almost on top of each other so it doesn't 
                ; really matter if you take int_i[1] or int_i[2] (I hope):
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
              IF (direction EQ BACKWARD) THEN BEGIN
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

            print, '         <><><><><><><><>'
            print, 'storing contour...',sht_count,store_forced
            contour_n = contour_n + 1
            contour_data = {  $
              state        : 0                 ,  $
              origin       : 1                 ,  $
              separatrix   : 0                 ,  $
              region       : region            ,  $
              direction    : direction         ,  $
              active_2nd   : active_2nd        ,  $
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
            print, 'broke_contour', broke_contour
            print, 'store_forced ', store_forced
            print, 'end1         ', psi_val LE psi_end
            print, 'end2         ', psi_val GE psi_end
            print, 'sht_count    ', sht_count

            IF (broke_contour OR  $
                store_forced  OR  $
                (direction EQ FORWARD  AND psi_val LE psi_end) OR  $
                (direction EQ BACKWARD AND psi_val GE psi_end) OR  $
                sht_count GE 3) THEN BEGIN

              PRINT, 'setting active contour'

              grid_SetActiveContour, contour_array, int_nlast, region, active_contour, store_forced,  $
                                     psi_start, psi_end, psi_val, psi_adjust, psi_step, sht_count, last_good_psi_val,   $
                                     boundary1_p1, boundary1_p2, boundary2_p1, boundary2_p2,  $
                                     the_search_continues, contour_n, direction, process_2nd,  $
                                     psi_1st_xpoint, psi_2nd_xpoint, user_step, user_finish
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
              if (0 EQ 1 AND contour_n EQ 20) then begin
                xrange = user_xrange
                yrange = user_yrange
;                xrange = [ 3.8,4.10]  ; upper CC
;                yrange = [ 2.0, 3.0]
;                xrange = [ 3.8,4.10]  ; lower inner CC
;                yrange = [-2.6,-1.0]
                xrange = [ 4.5,5.5]  ; top
                yrange = [ 4.0,5.0]
;                xrange = [5.0,7.0] ; top outer
;                yrange = [3.6,4.8] 
help,contour_array, /struct
                grid_Debug,b,contour_array,wall,xrange,yrange
              endif

            ENDIF
          ENDIF

;          END
;       ----------------------------------------------------------------
;        ELSE: BEGIN
;          PRINT,'ERROR grid_AnalyseBoundary: Unrecognized grid region'
;          STOP
;          END
;      ENDCASE

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

;     IF (active_contour EQ 0) THEN BEGIN
;       PRINT,'ERROR'
;       STOP
;     ENDIF
      IF (last_good_psi_val EQ 0.0D) THEN BEGIN
        PRINT, 'ABANDONING ACTIVE CONTOUR'
        IF (process_2nd EQ 3) THEN failure_2nd_pfz = 1
        ; Choose a new contour and continue:
        grid_SetActiveContour, contour_array, int_nlast, region, active_contour, store_forced,  $
                               psi_start, psi_end, psi_val, psi_adjust, psi_step, sht_count, last_good_psi_val,   $
                               boundary1_p1, boundary1_p2, boundary2_p1, boundary2_p2,  $
                               the_search_continues, contour_n, direction, process_2nd,  $
                               psi_1st_xpoint, psi_2nd_xpoint, user_step, user_finish
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
      IF ((direction EQ FORWARD  AND psi_val GT psi_start) OR  $
          (direction EQ BACKWARD AND psi_val LT psi_start)) THEN BEGIN
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
;      print,'psi_values: ',psi_start,psi_val,psi_end
    ENDWHILE



    IF (process_2nd EQ 1 AND the_search_continues EQ 0) THEN BEGIN
      print, ' '
      print, ' >>>>>>>>>>>>> 2nd x-point <<<<<<<<<<<< contour_n',contour_n
      print, ' '

      ; grid_Debug,b,contour_array,wall,xrange,yrange

      psi_start         = psi_2nd_xpoint + user_step
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

    IF (1 EQ 1 AND process_pfz EQ 0 AND the_search_continues EQ 0) THEN BEGIN
      print, ' '
      print, ' >>>>>>>>>>>>> primary pfz <<<<<<<<<<<< contour_n',contour_n
      print, ' '

      IF (process_2nd EQ 2) THEN process_2nd = 3

;      grid_Debug,b,contour_array,wall,xrange,yrange


      psi_start         = psi_1st_xpoint - user_step
      psi_end           = psi_start + user_finish * 0.3D ; replace 0.3 with specific pfz parameter
      psi_val           = psi_start
      psi_adjust        = psi_step
      last_psi_val      = -1.0D
      last_good_psi_val =  0.0D

      print,user_step
      print,psi_1st_xpoint
      print,psi_start
      print,psi_val
      print,psi_end

      the_search_continues = 1
      active_contour       = 0
      store_forced         = 0

      direction = BACKWARD

      boundary1_p1 = [-999.0D,-999.0D]
      boundary1_p2 = [-999.0D,-999.0D]
      boundary2_p1 = [-999.0D,-999.0D]
      boundary2_p2 = [-999.0D,-999.0D]

      process_pfz = 1

      region = PFZ
      int_nlast[region] = 2  ; Rubbish, remove [region] in the code...
    ENDIF

    IF (1 EQ 1 AND process_pfz     EQ 1 AND  $ ; process_2nd GE 0 AND  $
                   process_2nd_pfz EQ 1 AND the_search_continues EQ 0) THEN BEGIN
      print, ' '
      print, ' >>>>>>>>>>>>> secondary pfz <<<<<<<<<<<< contour_n',contour_n
      print, ' '

      IF (process_2nd EQ 2) THEN process_2nd = 3

      ;xrange = [ 4.5,5.5]  ; top
      ;yrange = [ 4.0,5.0]    
      ;grid_Debug,b,contour_array,wall,xrange,yrange


      psi_start         = psi_2nd_xpoint - user_step * 0.01D
      psi_end           = psi_start + user_finish * 0.3D  ; replace with specific parameter
      psi_val           = psi_start
      psi_adjust        = psi_step
      last_psi_val      = -1.0D
      last_good_psi_val =  0.0D

      the_search_continues = 1
      active_contour       = 0
      store_forced         = 0

      direction = BACKWARD

      boundary1_p1 = [-999.0D,-999.0D]
      boundary1_p2 = [-999.0D,-999.0D]
      boundary2_p1 = [-999.0D,-999.0D]
      boundary2_p2 = [-999.0D,-999.0D]

      process_pfz = 2

      region = PFZ
      int_nlast[region] = 2  ; Rubbish, remove [region] in the code...
    ENDIF

    ; One time check to register which region is active (HFS or LFS) after 
    ; passing the secondary x-point when processing the SOL:
    IF (process_2nd EQ 1 AND active_2nd EQ 0) THEN BEGIN
      mean_x = MEAN(x,/DOUBLE)
      IF (mean_x LT b.x[b.null_i[2]]) THEN active_2nd = HFS ELSE active_2nd = LFS
      print, 'setting active_2nd'
      print,mean_x
      print, b.x[b.null_i[2]]
      print,active_2nd
    ENDIF


  ENDWHILE  ; PSIN loop


  print, '......all done......'


  scan_params = {  $
    geometry        : geometry        ,  $
    process_2nd     : process_2nd     ,  $
    process_2nd_pfz : process_2nd_pfz ,  $
    failure_2nd_pfz : failure_2nd_pfz }
                  
  
  IF (KEYWORD_SET(save)) THEN  $
    SAVE,filename='stored_scan_'+machine+'.sav',b,wall,contour_array,scan_params

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
;  xrange = [ 1.5,2.0 ]
;  yrange = [-1.4,-1.1]
;  xrange = [ 0.9, 1.1]
;  yrange = [-1.2,-0.6]
;  xrange = [ 0.9, 1.6]
;  yrange = [ 1.0, 1.5]
;         xrange = [3.5,6.5] ; top
;         yrange = [2.0,5.0] 
  grid_Debug,b,contour_array,wall,xrange,yrange


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
      1: result = CREATE_STRUCT(result,'psi_1st_xpoint',psi_xpt)
      2: result = CREATE_STRUCT(result,'psi_2nd_xpoint',psi_xpt-psi_step)
      ELSE: stop
    ENDCASE


  ENDFOR 

  RETURN,result

END
;
; ======================================================================
;
FUNCTION grid_GetOrthogonal, x_val,y_val,p1_val,direction,length,ictr=ictr,ctrs=ctrs,span=span
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------


  x = x_val
  y = y_val

  status_2nd = 0

  IF (KEYWORD_SET(ictr) AND KEYWORD_SET(ctrs)) THEN BEGIN
    tags = STRUPCASE(TAG_NAMES(ctrs))
    nctr = N_ELEMENTS(tags)   
    ctr  = grid_ExtractStructure(ctrs,tags[ictr-1])      
    IF (ctr.separatrix GT 1) THEN BEGIN
      PRINT, 'WHOA! SECONDARY SEPARATRIX DETECTED!'
      PRINT, ctr.separatrix

; left off

; the problem, now, naturally, is with the ITER grid, where there's a tangency point that faces
; away from the separatrix... can't exactly see what to do with it just yet...

   ; need to add a sector point corresponding to the secondary null I think, using the
   ; new cross-field vector -- not sure but I think it's required so that there's a proper
   ; division between the LFS and HFS when partioning the poloidal divisions

   ; then, for radial contour crossing the 2nd sep., need to convolve them with the null normal
   ; to make sure they conform, although there's no guarantee, I don't think, that there can't
   ; be cross-over problems for some equil+wall combinations...

   ; need to check both the ctr.sep=3,4 contours to find the one with the vector in it...

   ; and perhaps need to recheck the proximity of P1 to the 2nd null here, since a better
   ; point will have been added since the vector was generated...


      ; Find out which contour has the cross-field vector that's associated 
      ; with the secondary null: 

      FOR i = 1, nctr DO BEGIN
        ctr = grid_ExtractStructure(ctrs,tags[i-1])      

        print,ctr.separatrix,ctr.tangent_p1[0]

        IF ((ctr.separatrix EQ 2) OR  $
            ((ctr.separatrix EQ 3 OR ctr.separatrix EQ 4) AND  $
              ctr.tangent_p1[0] NE 0.0D0)) THEN BREAK
      ENDFOR
      IF (i EQ nctr+1) THEN BEGIN
        PRINT, 'ERROR grid_GetOrthogonal: Secondary null cross-field vector not found'
        STOP
      END

      p1 = ctr.tangent_p1
      p2 = 2.0D0* p1 - ctr.tangent_p2
      vx = p2[0] - p1[0]
      vy = p2[1] - p1[1]
      length_now = SQRT(vx^2 + vy^2)
      p2 = [p1[0] + length / length_now * vx,  $
            p1[1] + length / length_now * vy]
      p1_2nd = p1
      p2_2nd = p2
      OPLOT,[p1_2nd[0],p2_2nd[0]],[p1_2nd[1],p2_2nd[1]],color=Truecolor('Red')

      status_2nd = 1

    ENDIF
  ENDIF

  p1 = p1_val

  ; Find where the point is along the contour:
  n = N_ELEMENTS(x)
  FOR i = 0, n-2 DO BEGIN
    IF (ABS(x[i+1]-x[i]) GT 1.0D-8) THEN s = (p1[0] - x[i]) / (x[i+1] - x[i]) ELSE s = -999.0D
    IF (ABS(y[i+1]-y[i]) GT 1.0D-8) THEN t = (p1[1] - y[i]) / (y[i+1] - y[i]) ELSE t = -999.0D
    ;print,i,s,t
    IF ((ABS(s-t) LT 1.0D-7 AND s LT 0.9999999D) OR  $
        (s EQ -999.0D AND (t GE 0.0D AND t LT 1.0D)) OR  $
        (t EQ -999.0D AND (s GE 0.0D AND s LT 1.0D))) THEN BREAK 
  ENDFOR
  IF (i EQ n-1) THEN BEGIN
    PRINT, 'ERROR grid_GetOrthogonal: Location of point on line not found'
    STOP
  ENDIF
  IF (ABS(s) GT 1.0E-7) THEN BEGIN
    ; This is a non-standard situation at present, so some of the logic below
    ; needs to be modified to make this case valid.
    PRINT, 'ERROR grid_GetOrthogonal: Normal vector requested away from vertex'
    STOP
  ENDIF

;  print,'test',i

  ; Scan a few points in either direction to try and account for wiggle
  ; in the contour, especially for points that are very close together:
  IF (NOT KEYWORD_SET(span)) THEN span = 4
  IF (i LT span OR N_ELEMENTS(x)-1-i LT span) THEN BEGIN
    ; Too close to the end of the ring to do anything other than
    ; use the immediate line segment:
    x1 = x[i-1]
    y1 = y[i-1]
    x2 = x[i+1]
    y2 = y[i+1]
  ENDIF ELSE BEGIN
    i1 = i - span
    i2 = i + span 
    length1 = grid_Length(x[i1:i],y[i1:i])
    length2 = grid_Length(x[i:i2],y[i:i2])
    IF (length2 GT length1) THEN BEGIN
      ;print, 'fixing 2'
      x1 = x[i1]
      y1 = y[i1]
      dist = 0.0D
      FOR i2 = i+1, N_ELEMENTS(x)-1 DO BEGIN
        dist2 = grid_Length(x[i2-1:i2],y[i2-1:i2])
        IF (dist+dist2 GT length1) THEN BEGIN
          frac = (length1 - dist) / dist2
          x2 = x[i2-1] + frac * (x[i2] - x[i2-1])
          y2 = y[i2-1] + frac * (y[i2] - y[i2-1])
          dist = dist + SQRT( (x2-x[i2-1])^2 + (y2-y[i2-1])^2 )
          BREAK
        ENDIF
        dist = dist + dist2
      ENDFOR
    ENDIF ELSE BEGIN
      ;print, 'fixing 1'
      x2 = x[i2]
      y2 = y[i2]
      dist = 0.0D
      FOR i1 = i-1, 0, -1 DO BEGIN
        dist1 = grid_Length(x[i1:i1+1],y[i1:i1+1])
        IF (dist+dist1 GT length2) THEN BEGIN
          frac = (length2 - dist) / dist1
          x1 = x[i1+1] + frac * (x[i1] - x[i1+1])
          y1 = y[i1+1] + frac * (y[i1] - y[i1+1])
          dist = dist + SQRT( (x1-x[i1+1])^2 + (y1-y[i1+1])^2 )
          BREAK
        ENDIF
        dist = dist + dist1
      ENDFOR
    ENDELSE
    print, 'lengths',length1,length2,dist,frac
  ENDELSE

  vx = -(y1 - y2)
  vy =   x1 - x2
  IF (direction EQ 2) THEN BEGIN
    vx = -vx
    vy = -vy
  ENDIF
  length_now = SQRT(vx^2 + vy^2)
  p2 = [p1[0] + length / length_now * vx,  $
        p1[1] + length / length_now * vy]


  IF (status_2nd) THEN BEGIN
    ctr = grid_ExtractStructure(ctrs,tags[ictr-1])      
    separatrix = ctr.separatrix

    focus_x = p1_2nd[0]
    focus_y = p1_2nd[1]

    proximity = SQRT ( (x-focus_x)^2 + (y-focus_y)^2 )
    dummy = MIN(proximity,i_2nd)

    IF ((i GT i_2nd AND separatrix EQ 3) OR  $
        (i LT i_2nd AND separatrix EQ 4)) THEN BEGIN
      PRINT, 'NOT READY, SORRY MAN'
      stop
    ENDIF

    IF (i LE i_2nd) THEN dist1 = grid_Length(x[i:i_2nd],y[i:i_2nd]) ELSE  $
                         dist1 = grid_Length(x[i_2nd:i],y[i_2nd:i]) 

    ; Need to take the length of the primary separatrix as some kind of 
    ; distance scale since the secondary separatrix can be very short if
    ; it hits the wall:
    ctr = grid_ExtractStructure(ctrs,tags[0])      
    IF (ctr.separatrix NE 1) THEN BEGIN
      PRINT, 'ERROR grid_GetOrthogonal: Not the separatrix'
      STOP
    ENDIF
    n = N_ELEMENTS(ctr.x)
    dist2 = grid_Length(ctr.x[0:n-1],ctr.y[0:n-1])

    frac = dist1 / dist2

    print,i,i_2nd,frac,dist1,dist2
    threshold = 0.10D
    IF (frac LT threshold) THEN BEGIN
      frac = frac / threshold
      print,frac
      v1 = p2     - p1 
      v2 = p2_2nd - p1_2nd
      v3 = frac * v1 + (1.0D - frac) * v2
      p2 = p1 + v3
    ENDIF

  ENDIF



  result = p2

  RETURN, result
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

;  IF (mode EQ -1) THEN BEGIN
;    ; Special case: secondary separatrix...
;    focus_x = b.x[b.null_i[2]]
;    focus_y = b.y[b.null_j[2]]
;  ENDIF ELSE BEGIN
;    IF (KEYWORD_SET(focus_x) AND KEYWORD_SET(focus_y) AND  $
;        NOT KEYWORD_SET(psi_val)) THEN BEGIN
;
;      ;PLOT,[focus_x],[focus_y],color=Truecolor('Blue'),PSYM=6,  $
;      ;    XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
;
;      n = N_ELEMENTS(psi_x)
;      i = WHERE(focus_x GE psi_x[0:n-2] AND focus_x LT psi_x[1:n-1], count)
;      IF (count EQ 0) THEN stop      
;      n = N_ELEMENTS(psi_y)
;      j = WHERE(focus_y GE psi_y[0:n-2] AND focus_y LT psi_y[1:n-1], count)
;      IF (count EQ 0) THEN stop
;
;      print,focus_y
;      print,psi_y[j],psi_y[j+1]
;
;      psi_1 = INTERPOL(psi[*,j  ],psi_x,focus_x)
;      psi_2 = INTERPOL(psi[*,j+1],psi_x,focus_x)
;
;      frac = (focus_y - psi_y[j]) / (psi_y[j+1] - psi_y[j])
;
;      psi_val = psi_1 + frac * (psi_2 - psi_1)
;
;;      ctr = grid_ExtractContour(psi, psi_x, psi_y, psi_val)
;;      x = ctr.x
;;      y = ctr.y
;;      PLOT,x,y,color=Truecolor('Purple'),  $
;;           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
;;stop
;
;;test =1 
;
;    ENDIF
;  ENDELSE


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
            print,'problem here as well, again 1',N_ELEMENTS(inter.i)
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
            xrange = [inter.x[0] - 0.1D, inter.x[0] + 0.1D]
            yrange = [inter.y[0] - 0.1D, inter.y[0] + 0.1D]
            PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
            OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
            OPLOT,[x[min_i]],[y[min_i]],color=Truecolor('Magenta'), PSYM=6
            OPLOT,[x[i]],[y[i]],color=Truecolor('Orange'), PSYM=6
            OPLOT,x,y,color=Truecolor('Magenta')
            OPLOT,x,y,color=Truecolor('White'), PSYM=3
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
      IF (psi_val LT b.psi_2nd_xpoint) THEN BEGIN  ;  THE SOL COULD QUALIFY FOR PFZ ***
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
          p2 = grid_GetOrthogonal(x2,y2,p1,1,0.2D,span=3)
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
FUNCTION grid_InstallXPoints, b, scan_params, contour_array
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------

  tags      = STRUPCASE(TAG_NAMES(contour_array))
  contour_n = N_ELEMENTS(tags)      


  FOR icnt = 1, contour_n DO BEGIN    

    tag = 'contour' + STRING(icnt,FORMAT='(I0)')
    cnt = grid_ExtractStructure(contour_array,tag)  

    status = 0
    IF (ABS(cnt.psi-b.psi_1st_xpoint) LT 1.0D-6) THEN BEGIN
      status = 1
      focus_x = b.x[b.null_i[1]]
      focus_y = b.y[b.null_j[1]]
    ENDIF
    IF (cnt.origin EQ 2) THEN BEGIN
      status = 2
      focus_x = b.x[b.null_i[2]]
      focus_y = b.y[b.null_j[2]]
    ENDIF
    IF (status EQ 0) THEN CONTINUE

    print,tag

    x = cnt.x
    y = cnt.y      

    print,'status',status

    FOR i = N_ELEMENTS(y)-2, 0, -1 DO BEGIN 

;      print,icnt,i,x[i],focus_x

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
          print, 'found'
          j = i
        ENDIF
      ENDELSE

      IF (j NE 0) THEN BEGIN
;        print,focus_x,x[i],x[i+1]
;        print,focus_y,y[i],y[i+1]
;        print, 'OK!'
        IF (x[i  ] EQ focus_x AND y[i  ] EQ focus_y) THEN BEGIN
          
        ENDIF ELSE BEGIN
          x = [x[0:i],focus_x,x[i+1:N_ELEMENTS(x)-1]]
          y = [y[0:i],focus_y,y[i+1:N_ELEMENTS(y)-1]]
          j++
        ENDELSE
        save_j = j
      ENDIF

    ENDFOR

    cnt = grid_UpdateStructure(cnt,'x',x)    
    cnt = grid_UpdateStructure(cnt,'y',y)    
;    cnt.separatrix = status

    IF (status EQ 2 AND scan_params.failure_2nd_pfz EQ 1) THEN cnt.tangent_i = save_j

;   NEED TO CLEAN UP AROUND THE X-POINTS, i.e. delete some points near by and then interpolate
;   new ones
;   need to send some kind of warning if the resolution of the equilibirum mesh is poor

    contour_array = grid_UpdateStructure(contour_array,tag,cnt)

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
PRO grid_TrimContours, b, contour_array, wall,  $
                       debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------

  IF (NOT KEYWORD_SET(xrange)) THEN xrange = [ 3.0,9.0]  ; lame
  IF (NOT KEYWORD_SET(yrange)) THEN yrange = [-6.0,6.0]

  result = contour_array

  tags = STRUPCASE(TAG_NAMES(contour_array))
  nctr = N_ELEMENTS(tags)

;***  perhaps need to re-specify wall_pt from wall.x/y here, for a fresh start each time this
;     routine is called

  wall_pt1 = wall.pt1
  wall_pt2 = wall.pt2

  FOR ictr = 0, nctr-1 DO BEGIN

print, 'ictr',ictr

    tag = tags[ictr]
    ctr = grid_ExtractStructure(contour_array,tag)      

    x = ctr.x
    y = ctr.y
 
;    x[0] = x[1] + 2.0D * (x[0] - x[1])
;    y[0] = y[1] + 2.0D * (y[0] - y[1])

    tangent_i = ctr.tangent_i

    OPLOT,x,y,color=Truecolor('Red'),psym=3

    FOR i = 0, N_ELEMENTS(x)-2 DO BEGIN
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
      IF (tangent_i NE -1) THEN tangent_i = tangent_i - i
    ENDELSE

    FOR i = N_ELEMENTS(x)-2, 1, -1 DO BEGIN
      result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                 wall_pt1, wall_pt2, 1, status=status)
      IF (status) THEN BREAK
    ENDFOR
    IF (NOT status) THEN BEGIN
      PRINT,'ERROR grid_TrimContours: No wall intersection found 2'
      STOP
    ENDIF ELSE BEGIN
      x = x[0:i+1]
      y = y[0:i+1]
;      ctr.wall_p2 = [result.x[0],result.y[0]]
    ENDELSE

    OPLOT,x,y,color=Truecolor('Pink')    

    ctr = grid_UpdateStructure(ctr,'x',x)    
    ctr = grid_UpdateStructure(ctr,'y',y)    
    ctr.tangent_i = tangent_i
    contour_array = grid_UpdateStructure(contour_array,tag,ctr)

  ENDFOR

END
;
; ======================================================================
;
PRO grid_UpdateWall, b, contour_array, wall,  $
                          debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------



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


  FOR ictr = 0, nctr-1 DO BEGIN

print, 'ictr targets ',ictr,' ',tags[ictr]

    tag = tags[ictr]
    ctr = grid_ExtractStructure(contour_array,tag)      

    x = ctr.x
    y = ctr.y
 
    tangent_i = ctr.tangent_i

    OPLOT,x,y,color=Truecolor('Red'),psym=3

    i = WHERE(wall_ptc EQ ictr+1 AND wall_ptt EQ 1, count)
    IF (count EQ 0) THEN BEGIN
      result = grid_Intersection([x[0],y[0]], [x[1],y[1]],  $
                                 wall_pt1, wall_pt2, 1, status=status)
      IF (NOT status) THEN BEGIN
        PRINT,'ERROR grid_UpdateWall: No wall intersection found 1'

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
        print,'damn 1' ; ,result.i,N_ELEMENTS(wall_pti)-1
        IF (N_ELEMENTS(result.i) GT 1) THEN BEGIN
          PRINT,'ERROR grid_UpdateWall: Multiple intersections detected 1'

          xrange = [result.x[0] - 0.003D, result.x[0] + 0.003D]
          yrange = [result.y[0] - 0.003D, result.y[0] + 0.003D]
          PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
          OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
          OPLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Orange')
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
          IF (dist LT 0.001D) THEN BEGIN
            ;print, 'close 1'
            ;OPLOT,[wall_pt1[0,j]], [wall_pt1[1,j]],color=Truecolor('Lightblue'), PSYM=4
            ;print,wall_pti[j],wall_ptc[j],wall_ptt[j]
            ;print,i,j,t,dist
            IF (wall_ptc[j] NE -1) THEN BEGIN
              IF (wall_ptc[j] EQ ictr + 1) THEN BEGIN
                print,'duplicate 1',wall_ptc[j]
                status = 0 
              ENDIF ELSE BEGIN  ; Assume that this point was assigned to this ring on the previous pass
                PRINT, 'WARNING grid_UpdateWall: Contour wall points very close together (1)'

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
              print,ictr+1
              wall_ptc[j] = ictr + 1
              wall_ptt[j] = 1
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
                               i, ictr+1, 1, result.x[0], result.y[0]
        ENDELSE
      ENDELSE
    ENDIF

    i = WHERE(wall_ptc EQ ictr+1 AND wall_ptt EQ 2, count)
    IF (count EQ 0) THEN BEGIN
      n = N_ELEMENTS(x)
      result = grid_Intersection([x[n-2],y[n-2]], [x[n-1],y[n-1]],  $
                                 wall_pt1, wall_pt2, 1, status=status)
      IF (NOT status) THEN BEGIN
        PRINT,'ERROR grid_UpdateWall: No wall intersection found 2'
        STOP
      ENDIF ELSE BEGIN
        print,'damn 2' ; ,result.i,N_ELEMENTS(wall_pti)-1
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
          IF (dist LT 0.001D) THEN BEGIN
            ;print, 'close 2'
            ;OPLOT,[wall_pt1[0,j]], [wall_pt1[1,j]],color=Truecolor('Lightblue'), PSYM=4
            ;print,wall_pti[j],wall_ptc[j],wall_ptt[j]
            ;print,i,j,t,dist
            IF (wall_ptc[j] NE -1) THEN BEGIN
              IF (wall_ptc[j] EQ ictr + 1) THEN BEGIN
                print,'duplicate 2',wall_ptc[j]
                status = 0 
              ENDIF ELSE BEGIN  ; Assume that this point was assigned to this ring on the previous pass
                PRINT, 'WARNING grid_UpdateWall: Contour wall points very close together (2)'

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
              print,ictr+1
              wall_ptc[j] = ictr + 1
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
                               i, ictr+1, 2, result.x[0], result.y[0]
;          IF (0 EQ 1 AND (t LT 0.01D OR t GT 0.99D)) THEN BEGIN
;            print, 'close 2, stopping'
;            stop
;          ENDIF ELSE BEGIN
;           grid_AddWallPoint, wall_pti, wall_ptc, wall_ptt, wall_pt1, wall_pt2,  $
;                              i, ictr+1, 2, result.x[0], result.y[0]
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
; Look after tangency points:
;


  FOR ictr = 0, nctr-1 DO BEGIN

    print, 'ictr tangency ',ictr,' ',tags[ictr]

    tag = tags[ictr]
    ctr = grid_ExtractStructure(contour_array,tag)      

    IF (ctr.tangent_i EQ -1) THEN CONTINUE

    ; Check if tangency point has already been assigned:
    i = WHERE(wall_ptc EQ ictr+1 AND wall_ptt EQ 3, count)
    IF (count GT 0) THEN CONTINUE

    x         = ctr.x
    y         = ctr.y
    tangent_i = ctr.tangent_i
    p1        = ctr.tangent_p1
    p2        = ctr.tangent_p2
 
    length = SQRT( (p1[0]-p2[0])^2 + (p1[1]-p2[1])^2 )
    p3 = p1 + 0.001D / length * (p2 - p1)
    p4 = p1 - 0.001D / length * (p2 - p1)

    OPLOT,[p3[0],p4[0]], [p3[1],p4[1]], color=Truecolor('Green') 
    OPLOT,x,y,color=Truecolor('Red')    
    OPLOT,x,y,color=Truecolor('White'),PSYM=3

    result = grid_Intersection([p3[0],p3[1]], [p4[0],p4[1]],  $
                               wall_pt1, wall_pt2, 1, status=status)
    IF (NOT status) THEN BEGIN
      PRINT,'ERROR grid_UpdateWall: No tangency point wall intersection found'
      STOP
    ENDIF ELSE BEGIN
      print,'hot damn 1' ; ,result.i,N_ELEMENTS(wall_pti)-1
      IF (N_ELEMENTS(result.i) GT 1) THEN BEGIN
        PRINT,'ERROR grid_UpdateWall: Multiple tangency intersections detected'
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
            IF (wall_ptc[j] EQ ictr + 1) THEN BEGIN
              print,'duplicate tangency',wall_ptc[j]
              status = 0 
            ENDIF ELSE  $  ; Assume that this point was assigned to this ring on the previous pass
              PRINT, 'WARNING grid_UpdateWall: Contour wall points very close together'
          ENDIF ELSE BEGIN
            print, 'taking over wall point with tangency'
            print,x[tangent_i],y[tangent_i]
            wall_ptc[j] = ictr + 1
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
                             i, ictr+1, 3, x[tangent_i], y[tangent_i]
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

    if (count eq 8) then rage = 1 else rage = 0


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
;    print,j,wall.pti[j],wall.pt1[0,j],wall.pt1[1,j],wall.pt2[0,j],wall.pt2[1,j]
;  ENDFOR


;  n=N_ELEMENTS(wall_pti)
;  FOR j = 0, n-1 DO BEGIN
;    print,j,wall.pti[j],wall.pt1[0,j],wall.pt1[1,j],wall.pt2[0,j],wall.pt2[1,j]
;  ENDFOR

  grid_ZoneWall, wall.pt1, wall.pt2, debug=debug, xrange=xrange, yrange=yrange

  RETURN
END
;
; ======================================================================
;
FUNCTION grid_GetRegion, c_array, tag

  tags  = STRUPCASE(TAG_NAMES(c_array))

;    ctr = grid_ExtractStructure(c_array,'CONTOUR'+STRTRIM(STRING(ictr),2))      
    print,'tag:'+tag
;    help,c_array,/struct
    j = -1 
    FOR i = 0, N_ELEMENTS(tags)-1 DO BEGIN
;      print,tags[i],tag eq tags[i]
      IF (tag EQ tags[i]) THEN BEGIN
        j = i
        BREAK
      ENDIF
    ENDFOR
    IF (j EQ -1) THEN STOP
;    i = WHERE(TAG_NAMES(c_array) EQ STRUPCASE(tag), count)
;    IF (count EQ 0) THEN STOP
;    ctr = c_array.(i[0]) 
    ctr = c_array.(j) 
;help,ctr,/struct
;stop
  result = ctr.region

  RETURN, result
END
;
; ======================================================================
;
FUNCTION grid_BuildRadialMap, c_array, wall, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------

  tags  = STRUPCASE(TAG_NAMES(c_array))
  nctr  = N_ELEMENTS(tags)
  nwall = N_ELEMENTS(wall.ptc)

  istart = WHERE(wall.ptc EQ 1 AND wall.ptt EQ 1, count)
  IF (count EQ 0) THEN BEGIN
    PRINT, 'ERROR grid_BuildRadialMap: Starting contour not found'
    STOP
  ENDIF ELSE IF (KEYWORD_SET(debug)) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
    OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
    OPLOT,[wall.pt1[0,istart]],[wall.pt1[1,istart]],color=Truecolor('Hotpink'),PSYM=6
    FOR i = 0, nctr-1 DO BEGIN
      ctr = grid_ExtractStructure(c_array,tags[i])      
      OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
      npts = N_ELEMENTS(ctr.x)
      XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('White')
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

;    ctr = grid_ExtractStructure(c_array,'CONTOUR'+STRTRIM(STRING(ictr),2))      

    IF (KEYWORD_SET(debug)) THEN  $
      OPLOT,[wall.pt1[0,iwall]],[wall.pt1[1,iwall]],color=Truecolor('Hotpink'),PSYM=1

    print,'OK---',ictr,itar,ictr_last,itar_last,FORMAT='(A8,4I8)'    

    IF (itar EQ 1 AND (itar_last EQ 1 OR itar_last EQ 3)) THEN BEGIN
;      i = WHERE(map_in [ictr-1,*] EQ -1,count)
;      IF (count GT 0) THEN map_in [ictr-1,i[0]] = ictr_last ELSE BEGIN
;        PRINT, 'ERROR grid_BuildRadialMap: Inward overload'
;        FOR i = 0, nctr-1 DO print,i+1,map_in[i,0:1],format='(3I6)'
;        STOP
;      ENDELSE

      tag = 'CONTOUR'+STRTRIM(STRING(ictr),2)
      region = grid_GetRegion(c_array,tag)

      IF (region EQ SOL OR itar_last EQ 1) THEN i = 0 ELSE i = 1
      IF (map_in[ictr-1,i] NE -1) THEN BEGIN
        PRINT, 'ERROR grid_BuildRadialMap: Inward overload'
        FOR i = 0, nctr-1 DO print,i+1,map_in[i,0:1],format='(3I6)' 
        STOP
      ENDIF
      map_in[ictr-1,i] = ictr_last 
      map_in[ictr-1,2] = region

      ; Only assign outward map to previous ring if the previous wall
      ; point was at the end of a contour, and not a tangency point: 
;      i = WHERE(map_out[ictr_last-1,*] EQ -1,count)

;print,ictr_last
;print,count
;print,'start out',map_out[ictr_last-1,*],'end'
;print,i
;help,map_in
      tag = 'CONTOUR'+STRTRIM(STRING(ictr_last),2)
      region = grid_GetRegion(c_array,tag)

      IF (region EQ PFZ OR itar_last EQ 1) THEN i = 0 ELSE i = 1
      IF (map_out[ictr_last-1,i] NE -1) THEN BEGIN
        PRINT, 'ERROR grid_BuildRadialMap: Outward overload'
        FOR i = 0, nctr-1 DO print,i+1,map_out[i,0:1],format='(3I6)' 
        STOP
      ENDIF
      map_out[ictr_last-1,i] = ictr 
      map_out[ictr_last-1,2] = region

;      IF (count GT 0) THEN map_out[ictr_last-1,i[0]] = ictr ELSE BEGIN
;        PRINT, 'ERROR grid_BuildRadialMap: Outward overload'
;        FOR i = 0, nctr-1 DO print,i+1,map_out[i,0:1],format='(3I6)' 
;        STOP
;      ENDELSE


    ENDIF

    print,'  MAP',map_in [ictr-1,0:1],FORMAT='(A8,3I8)'    
    print,'     ',map_out[ictr-1,0:1],FORMAT='(A8,3I8)'    

    ictr_last = ictr
    itar_last = itar

  ENDWHILE



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
FUNCTION grid_SetPoloidalDomains, c_array, wall, debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON params, CORE, SOL, PFZ, FORWARD, BACKWARD, LOWER_NULL, UPPER_NULL, HFS, LFS
  ; ------------------------------------------------------------------

  tags  = STRUPCASE(TAG_NAMES(c_array))
  nctr  = N_ELEMENTS(tags)


  IF (KEYWORD_SET(debug)) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
    OPLOT,wall.x,wall.y,color=Truecolor('Yellow')
    FOR i = 0, nctr-1 DO BEGIN
      ctr = grid_ExtractStructure(c_array,tags[i])      
      OPLOT,ctr.x,ctr.y,color=Truecolor('Red')    
      npts = N_ELEMENTS(ctr.x)
      XYOUTS,ctr.x[npts/2],ctr.y[npts/2],STRTRIM(STRING(i+1),2),color=Truecolor('White')
    ENDFOR
  ENDIF



  ; There should be up to 2 contours mapped to any other contour, which
  ; occurs when there was a tagnecy point:
  map_in  = MAKE_ARRAY(nctr,3,/LONG,VALUE=-1L)
  map_out = map_in




  PRINT,'---------- START DOMAIN IDENTIFICATION ------------'


  ; Process tangency points:

  section = 0

  FOR ictr = 1, nctr DO BEGIN 
    ctr = grid_ExtractStructure(c_array,tags[ictr-1])      

    IF (ctr.tangent_i EQ -1) THEN CONTINUE

    print, 'ictr ----- ',ictr,' ',tags[ictr-1]

    OPLOT,ctr.x,ctr.y,color=Truecolor('Pink')    

    OPLOT,[ctr.tangent_p1[0]],[ctr.tangent_p1[1]],color=Truecolor('Lightgreen'),PSYM=6    

    ; Follow the trajectory into the separatrix:

    p1 =        ctr.tangent_p1
    p2 = 2.0D * ctr.tangent_p1 - ctr.tangent_p2

    OPLOT,[p1[0],p2[0]],[p1[1],p2[1]],color=Truecolor('Lightgreen')    

    section++
    dummy = WHERE(STRUPCASE(TAG_NAMES(ctr)) EQ 'SECTION',count)
    IF (count EQ 0) THEN BEGIN
      ctr = CREATE_STRUCT(ctr,'section',section,'section_x',p1[0],  $
                                                'section_y',p1[1])
    ENDIF ELSE BEGIN
      ctr = grid_UpdateStructure(ctr,'section'  ,[ctr2.section  ,section])        
      ctr = grid_UpdateStructure(ctr,'section_x',[ctr2.section_x,p1[0]  ])
      ctr = grid_UpdateStructure(ctr,'section_y',[ctr2.section_y,p1[1]  ])      
    ENDELSE


;    IF (ctr.region EQ SOL) THEN ictr2 = ctr.map_in [0] ELSE  $
;                                ictr2 = ctr.map_out[0]

    cnt=0
    ictr3 = ictr
    WHILE (1 EQ 1) DO BEGIN

      
      print, 'round =========>',cnt,ictr3

      ; Not trying to use the connection map here ecau
      min_ictr2 = -1
      min_ctr2  = -1
      min_inter = -1
      min_s     = 1.0D+06
      FOR ictr2 = 1, nctr DO BEGIN
        IF (ictr2 EQ ictr3) THEN CONTINUE
        ctr2 = grid_ExtractStructure(c_array,tags[ictr2-1])      
        n = N_ELEMENTS(ctr2.x)
        seg1 = MAKE_ARRAY(2,n-1,/DOUBLE,VALUE=0.0D)
        seg2 = seg1
        seg1[0,*] = ctr2.x[0:n-2]
        seg1[1,*] = ctr2.y[0:n-2]
        seg2[0,*] = ctr2.x[1:n-1]
        seg2[1,*] = ctr2.y[1:n-1]
        inter = grid_Intersection(p1, p2, seg1, seg2, 0, status=status)
        IF (status) THEN BEGIN
          print, 'intersection',ictr2,N_ELEMENTS(inter.s)
          IF (MIN(inter.s) LT min_s) THEN BEGIN
            min_inter = inter
            min_s     = MIN(inter.s)
            min_ictr2 = ictr2
            min_ctr2  = ctr2
          ENDIF
        ENDIF
      ENDFOR
      print,'taken -->',min_ictr2
      ictr2 = min_ictr2
      ctr2  = min_ctr2
      inter = min_inter
      dummy = MIN(inter.s,min_i)
      inter_i = inter.i[min_i]
      inter_x = inter.x[min_i]
      inter_y = inter.y[min_i]

      print, 'separatrix ====>',ictr,ctr2.separatrix

      OPLOT,ctr2.x,ctr2.y,color=Truecolor('White')    

      OPLOT,[inter_x],[inter_y],color=Truecolor('White'),PSYM=6
      OPLOT,[ctr2.x[inter_i  ]],[ctr2.y[inter_i  ]],color=Truecolor('Red'),PSYM=6
      OPLOT,[ctr2.x[inter_i+1]],[ctr2.y[inter_i+1]],color=Truecolor('Red'),PSYM=6


;      vx = -(ctr2.y[inter_i] - ctr2.y[inter_i+1])
;      vy =   ctr2.x[inter_i] - ctr2.x[inter_i+1]
;      IF (ctr.region EQ PFZ) THEN BEGIN
;        vx = -vx
;        vy = -vy
;      ENDIF 
;      length = SQRT(vx^2 + vy^2)
;      p1_new = [inter_x,inter_y]
;      p2_new = [p1_new[0] + 1.50D / length * vx,  $
;                p1_new[1] + 1.50D / length * vy]
                  
      ; Add point to the contour and register the designation of the 
      ; tangency point:
      ; Note that this intersection point is the results of a linear
      ; interpolation, so the curvature of the contour will be a bit
      ; off.  Might be worth refining the contour after the intersection 
      ; is identified and then finding a new intersection, in order
      ; to minimize this effect:
      n = N_ELEMENTS(ctr2.x)
      x = [ctr2.x[0:inter_i],inter_x,ctr2.x[inter_i+1:n-1]]
      y = [ctr2.y[0:inter_i],inter_y,ctr2.y[inter_i+1:n-1]]
      ctr2 = grid_UpdateStructure(ctr2,'x',x)
      ctr2 = grid_UpdateStructure(ctr2,'y',y)
      dummy = WHERE(STRUPCASE(TAG_NAMES(ctr2)) EQ 'SECTION',count)
      IF (count EQ 0) THEN BEGIN
        ctr2 = CREATE_STRUCT(ctr2,'section',section,'section_x',inter_x,  $
                                                    'section_y',inter_y)
      ENDIF ELSE BEGIN
        ctr2 = grid_UpdateStructure(ctr2,'section'  ,[ctr2.section  ,section   ])        
        ctr2 = grid_UpdateStructure(ctr2,'section_x',[ctr2.section_x,inter_x])
        ctr2 = grid_UpdateStructure(ctr2,'section_y',[ctr2.section_y,inter_y])      
      ENDELSE
      c_array = grid_UpdateStructure(c_array,tags[ictr2-1],ctr2)


      ; Quit when the separatrix has been reached:
      IF (ictr2 EQ 1) THEN BREAK


      ; Calculate the new cross-field trajectory:
      p1_new = [inter_x,inter_y]
      IF (ctr.region EQ PFZ) THEN direction = 2 ELSE direction = 1
      p2_new = grid_GetOrthogonal(ctr2.x,ctr2.y,p1_new,direction,1.5D,ictr=ictr2,ctrs=c_array)

      OPLOT,[p1_new[0],p2_new[0]],[p1_new[1],p2_new[1]],color=Truecolor('Blue')

      cnt++
      ;IF (cnt EQ 2) THEN STOP

;      print,'region',ctr2.region, ctr2.map_in [0]
;      IF (ctr2.region EQ SOL) THEN ictr2 = ctr2.map_in [0] ELSE  $
;                                   ictr2 = ctr2.map_out[0]
      p1 = p1_new
      p2 = p2_new              
 
      ictr3 = ictr2

    ENDWHILE

  ENDFOR




  stop


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
      user_step   = -0.003D
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
      xrange = [ 3.0,9.0]     
      yrange = [-6.0,6.0]
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
      ;path_equ = '/home/ITER/lisgos/divimp/shots/mast/13018_250/'
      ;file_equ = 'sonnet_13018_250.equ'
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
        ;xrange      = [0.4,0.7] ; divertor
        ;yrange      = [-0.6,-0.3]
        xrange      = [0.565,0.575]
        yrange      = [0.50,0.52]
        xrange      = [0.43,0.50]
        yrange      = [-0.1,0.1]
        xrange      = [0.3,1.2]
        yrange      = [-0.8,0.8]
      ENDIF ELSE BEGIN
        xrange      = [0.3,1.2]
        yrange      = [-0.8,0.8]
      ENDELSE
      xpoint_zone = [0.52,0.75,-0.45,0.45]
      path_wall = '/home/ITER/lisgos/divimp/shots/cmod/1100303017_0138/'
      file_wall = 'vessel_wall4.dat'  ; no pump plenum entrance
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
      ; dies on corner
      ; path_equ = '/home/ITER/lisgos/divimp/shots/ts2/600kA_LN_1cm/'
      ; file_equ = 'TS2_600kA_LN.x2.cr.equ'
      ; gets pissed off by the second x-point that's just inside the vessel
       path_equ = '/home/ITER/lisgos/divimp/shots/ts2/1MA_LN_1cm/'
       file_equ = 'TS2_1MA_LN.x2.equ'
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
        xrange = [ 1.0,1.2]
        yrange = [ 0.8,1.2]
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
      path_equ  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
      file_equ  = 'g105500.03500.x2.equ' ; LSND, 2nd x-point just inside the vessel
      ; works...
      ;path_equ  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
      ;file_equ  = 'g131245.02800.x2.equ' ; USND
      ; works...
      ;path_equ  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
      ;file_equ  = 'g134585.03005.x4.equ' ; UDND
      ; Corner at the entrance to the upper plenum causes problems:
      ;path_equ  = '/home/ITER/lisgos/divimp/shots/d3d/equilibria/'
      ;file_equ  = 'g135487.03865.x4.equ' ; LSND
      ; works (but had to get rid of the tiny step on the floor at the entrance to the lower pump plenum)...
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
      IF (scan_params.failure_2nd_pfz EQ 1) THEN  $
        psi_shift = 0.000001D ELSE psi_shift = -0.000001D

      status = grid_AddContour(b, wall, scan_params, contour_array, -1,  $
                               psi_val=b.psi_2nd_xpoint + psi_shift,  $
                               debug=debug, xrange=xrange, yrange=yrange)
    ENDIF

    grid_TrimContours, b, contour_array, wall,  $
                       debug=debug, xrange=xrange, yrange=yrange

    contour_array = grid_InstallXPoints(b,scan_params,contour_array)

    print, 'first wall update'
    grid_UpdateWall, b, contour_array, wall,  $
                     debug=debug, xrange=xrange, yrange=yrange

    grid_ProcessWall, b, wall, scan_params, contour_array, $
                      debug=debug, xrange=xrange, yrange=yrange


    print, 'subsequent wall update'
    grid_UpdateWall, b, contour_array, wall,  $
                     debug=debug, xrange=xrange, yrange=yrange


; add check to amek sure that the contours don intersect the walls except at the ends, other wise there
; which can happen near insufficiently resolved tangency points


;    help,contour_array,/struct
;   xrange = [6.4,7.0]
;   yrange = [3.7,4.0]


    contour_array = grid_BuildRadialMap(contour_array, wall,  $
                                        debug=debug, xrange=xrange, yrange=yrange)



    contour_array = grid_SetPoloidalDomains(contour_array, wall,  $
                                            debug=debug, xrange=xrange, yrange=yrange)



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











  RETURN, result
END