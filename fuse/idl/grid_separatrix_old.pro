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
; ======================================================================
;
; need to add some efficiency, i.e. domain decomposition (as per usual)
;
FUNCTION grid_Intersection, p1, p2, seg1, seg2, mode, status=status

  COMMON zone, zone_list, zone_index, zone_nx  , zone_ny  ,  $
               zone_minx, zone_maxx , zone_miny, zone_maxy

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
      IF ((i1 LT 0 OR i1 GT zone_nx) AND  $
          (j1 LT 0 OR j1 GT zone_ny) AND  $
          (i2 LT 0 OR i2 GT zone_nx) AND  $
          (j1 LT 0 OR j1 GT zone_ny)) THEN BEGIN
        status =  0
        result = -1
        RETURN, result
      ENDIF

;      print,i1,j1,i2,j2
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
  
    IF ( ABS(b1[0]-b2[0]) LT 1.0E-7 ) THEN b1[0] = b1[0] + 0.000001D
  
    Lint, a1, a2, b1, b2, c, d, flag=flag
  
    IF (flag NE 1) THEN CONTINUE
  
    s = (c[0] - a1[0]) / (a2[0] - a1[0])    
    t = (c[0] - b1[0]) / (b2[0] - b1[0])    
  
    IF ( s GE 0.0D AND s LE 1.0D AND  $
         t GE 0.0D AND t LT 1.0D) THEN BEGIN
  
      IF (NOT KEYWORD_SET(int_i)) THEN BEGIN
        int_i    = i
        int_s    = s
        int_t    = t      
        int_x    = c[0]
        int_y    = c[1]
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
FUNCTION grid_AnalyseBoundary, b, xpoint_zone, mode, debug=debug

  COMMON zone, zone_list, zone_index, zone_nx  , zone_ny  ,  $
               zone_minx, zone_maxx , zone_miny, zone_maxy


  psi   = b.psi_raw
  psi_x = b.x
  psi_y = b.y

  psib = b.psi[b.null_i[1],b.null_j[1]]

;  print,b.x[b.null_i[1]],b.y[b.null_j[1]]
;  print,b.b_pol[b.null_i[1],b.null_j[0]]

  result = CREATE_STRUCT( b,'psi_b',psib)

;stop
  IF (KEYWORD_SET(debug)) THEN BEGIN
    HELP,result,/STRUCT
    CONTOUR, psi, psi_x, psi_y, LEVELS=[psib-0.001,psib,psib+0.001]
  ENDIF


  a = grid_ExtractContour(psi, psi_x, psi_y, psib)
  help,a,/struct
  i = WHERE(a.dist GT 0.10)
  i = i[0]
;  PLOT,[MIN(a.x[0:i]),MAX(a.x[0:i])], [MIN(a.y[0:i]),MAX(a.y[0:i])] ,/NODATA
;  OPLOT,a.x[0:i],a.y[0:i]  NOT SURE WHAT THE FUCK IS GOING ON HERE!

;  ctr = grid_ExtractContour(psi, psi_x, psi_y, psib)
;  i = WHERE(ctr.dist GT 0.10)
;  x = ctr.x[0:i[0]]
;  y = ctr.y[0:i[0]]
;  OPLOT,x,y,color=Truecolor('White')


  xrange = [ 3.0,9.0]
  yrange = [-6.0,6.0]

  PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1

;  PLOT,[4.0,4.0,8.0,8.0,4.0],[-3.0,3.0,3.0,-3.0,-3.0],color=Truecolor('White'),XRANGE=xrange,YRANGE=yrange,/NOERASE,XSTYLE=1,YSTYLE=1

  wall = cortex_LoadAnnotationData(1,'/home/ITER/lisgos/divimp/shots/iter/1514/psi_wall.dat') 
;  wall = cortex_LoadAnnotationData(1,'/home/ITER/lisgos/divimp/shots/mast/default/main_wall.dat') ; MAST
  PLOT,wall.x,wall.y,color=Truecolor('Yellow'),XRANGE=xrange,YRANGE=yrange,/NOERASE,XSTYLE=1,YSTYLE=1





;  PLOT,[4.0,4.0,8.0,8.0,4.0],[-3.0,3.0,3.0,-3.0,-3.0],color=Truecolor('Red'),XRANGE=xrange,YRANGE=yrange,/NOERASE,XSTYLE=1,YSTYLE=1


;  print,'range:',!x.range,!y.range


;  OPLOT,[x[150]],[y[150]],color=Truecolor('Red'),PSYM=6

;  Rearrange wall data so that it doesn't need to be continuous: 
   n = N_ELEMENTS(wall.x)
   wall_pt1 = MAKE_ARRAY(2,n-1,/DOUBLE,VALUE=0.0D)
   wall_pt2 = wall_pt1
   wall_pt1[0,*] = DOUBLE(wall.x[0:n-2])
   wall_pt1[1,*] = DOUBLE(wall.y[0:n-2])
   wall_pt2[0,*] = DOUBLE(wall.x[1:n-1])
   wall_pt2[1,*] = DOUBLE(wall.y[1:n-1])
;   FOR j = 0, n-2 DO BEGIN
;     print,wall_pt1[1,j],wall_pt2[1,j]
;   ENDFOR



; Setup zones:
  zone_nx = 3
  zone_ny = 15
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

help,zone_list

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

      print, 'zone_index',i,j,zone_index[*,i,j]    
      PLOT,box_pt1[0,*],box_pt1[1,*],color=Truecolor('Purple'),  $
           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

      y = y + deltay
    ENDFOR
    x = x + deltax
  ENDFOR


help,zone_list

;  FOR psi_val = psib, psib, -.5 DO BEGIN
;  FOR psi_val = psib, psib-3.0, -.5 DO BEGIN
;  FOR psi_val = psib+3.0, psib-3.0, -.5 DO BEGIN
;  FOR psi_val = psib+0.015, psib-0.025, -0.001 DO BEGIN  ; MAST


; Need to have a pre-pass that looks to see where the tangency points are between the 
; flux surfaces and the wall (need to iterate to find the exact points)
; Throw away really short segments? 

; How do I avoid the distorted target cells that plagued the JET grid, i.e. the 
; cells overlapped?  Just keep the resolution high? 

; Need to bring in the other wall segments now I think, and look for intersections, rather 
; than checking if a point is inside or outside (although good to have some appropriate
; spatial refinement before this step? to make sure the entrances and exits are all accounted for?)

; Whatś going to be complicated is keeping good track (i.e. intelligently) of how the wall intersections are 
; evolving, so that I can properly collect the tangency contours, after which the assemply of the
; rings should proceed methodically... (famour last words...)

; Need to proceed outward from the separatrix, recognizing a change in the patter / structure of the
; accepted segments, where the segments are determined based on the point-in-wall and wall
; interesection filters, i.e. when the number of valid segments changes then you know to look for a 
; tangency point... although it may get complicated when a surface is exiting and entering the vessel 
; at various points

  int_nlast = MAKE_ARRAY(6,/LONG,VALUE=-1)
  CORE          = 0
  SOL           = 1
  SOL_LFS       = 2  
  SOL_HFS       = 3  
  PFZ           = 4  
  PFZ_SECONDARY = 5  

  psi_val    = DOUBLE(psib)
  psi_start  = psi_val
  psi_end    = psib - 3.0D
  psi_step   = -0.5D
  psi_adjust = psi_step
  tol_adjust = 1.0D-3

  contour_n = 0

  boundary1_p1 = [-999.0D,-999.0D]
  boundary1_p2 = [-999.0D,-999.0D]
  boundary2_p1 = [-999.0D,-999.0D]
  boundary2_p2 = [-999.0D,-999.0D]

  the_search_continues = 1

  WHILE (the_search_continues) DO BEGIN
;  WHILE (psi_val GT psi_end) DO BEGIN

    print,'psi_values: ',psi_val,psi_start,psi_end

    ctr = grid_ExtractContour(psi, psi_x, psi_y, psi_val)

    ibrk = WHERE(ctr.dist GT 0.10)
    nseg = N_ELEMENTS(ibrk) 
    ibrk = [-1,ibrk,ctr.n-1]

    FOR iseg = 0, nseg DO BEGIN

      print,'contour seg:',iseg,ibrk[iseg]+1,ibrk[iseg+1]
;      print,'  boundary1:',boundary1_p1,boundary1_p2
;      print,'  boundary2:',boundary2_p1,boundary2_p2

      ; Not a real segment, may need to strengthen this check:
      IF (ibrk[iseg]+1 EQ ibrk[iseg+1]) THEN CONTINUE  

      x = ctr.x[ibrk[iseg]+1:ibrk[iseg+1]]
      y = ctr.y[ibrk[iseg]+1:ibrk[iseg+1]]

      PLOT,x,y,color=Truecolor('White'),XRANGE=xrange,YRANGE=yrange,/NOERASE,XSTYLE=1,YSTYLE=1

      ; Only take the parts of the contour segment that are (mostly) inside the main wall,
      ; by dropping all points before and after the first and last intersections, respectively:       
      inside = grid_PointInPolygon(x,y,wall.x,wall.y)
      i = WHERE(inside EQ 1, count)
      IF (count EQ 0) THEN CONTINUE  
      n = N_ELEMENTS(i) - 1
      m = N_ELEMENTS(x) - 1
      IF (i[0] GT 0) THEN j = i[0] - 1 ELSE j = 0  ; Should check if contour entirely inside wall, and
      IF (i[n] LT m) THEN k = i[n] + 1 ELSE k = m  ; if this matches expectations...
      x      = x[j:k]
      y      = y[j:k]
      inside = inside[j:k]
      inside[0] = 1 & inside[N_ELEMENTS(inside)-1] = 1







      IF (boundary1_p1[0] NE -999.0D OR boundary2_p1[0] NE -999.0D) THEN BEGIN

;        xrange = [ 4.0,4.2]  ; First intersection on the inner wall
;        yrange = [ 0.4,0.6]
;        xrange = [ 3.8,5.2]  ; lower inner CC
;        yrange = [-3.0,2.0]
        xrange = [ 3.0,9.0]
        yrange = [-6.0,6.0]

        PLOT,wall.x,wall.y,color=Truecolor('Yellow'),  $
             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1
        PLOT,x,y,color=Truecolor('Orange'),  $
             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,PSYM=3,/NOERASE

        IF (boundary1_p1[0] NE -999.0D) THEN BEGIN
stop
          PLOT,[boundary2_p1[0],boundary2_p2[0]],  $
               [boundary2_p1[1],boundary2_p2[1]],color=Truecolor('Lightgreen'),  $
               XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

          FOR i = 0, N_ELEMENTS(x)-2 DO BEGIN
            ; Find intersections with the wall segments:
            result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                       boundary2_p1, boundary2_p2, 0, status=status)
            IF (status) THEN BREAK
          ENDFOR
          IF (inside[i] EQ 0) THEN BEGIN
            FOR j = i, 0, -1 DO IF (inside[j] EQ 1) THEN BREAK
            j = j + 1
            inside[j] = 1
          ENDIF ELSE j = i
          x      = x[0:j]
          y      = y[0:j]
          inside = inside[0:j]

          print,'i,j:',i,j
          PLOT,x,y,color=Truecolor('Lightgreen'),  $
               XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,PSYM=3,/NOERASE
        ENDIF


        IF (boundary2_p1[0] NE -999.0D) THEN BEGIN
          PLOT,[boundary2_p1[0],boundary2_p2[0]],  $
               [boundary2_p1[1],boundary2_p2[1]],color=Truecolor('Lightgreen'),  $
               XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

          FOR i = 0, N_ELEMENTS(x)-2 DO BEGIN
            ; Find intersections with the wall segments:
            result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                       boundary2_p1, boundary2_p2, 0, status=status)
            IF (status) THEN BREAK
          ENDFOR
          IF (NOT status) THEN BEGIN
            PRINT, 'NO INTERSECTION WITH BOUNDARY FOUND!'
pause
; f !d.name eq 'X' and !p.multi(0) eq 0 then begin
;  cursor,xw,yw,/nowait,/device
; 
;  if xw eq -1 then begin
;    if keyword_set(on_screen) then continue
;    if n_elements(xcrs) eq 0 then xcrs=.9*!d.x_vsize 
;    if n_elements(ycrs) eq 0 then ycrs=.9*!d.y_vsize 
;    tvcrs,xcrs,ycrs
;  endif
;  device,cursor_standard=16
;  cursor,xdum,ydum,/wait,/device
;  if keyword_set(up) then cursor,xdum,ydum,/up,/device
;  device,cursor_standard=30
;endif
            CONTINUE  ; ***
          ENDIF
          IF (inside[i] EQ 0) THEN BEGIN
            FOR j = i, 0, -1 DO IF (inside[j] EQ 1) THEN BREAK
            j = j + 1
            inside[j] = 1
          ENDIF ELSE j = i
          x      = x[0:j]
          y      = y[0:j]
          inside = inside[0:j]

          print,'i,j:',i,j
          PLOT,x,y,color=Truecolor('Lightgreen'),  $
               XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,PSYM=3,/NOERASE
        ENDIF

;help,result,/struct

      ENDIF ELSE BEGIN
        xrange = [ 3.8,4.2]
        yrange = [ 1.0,2.0]
;        xrange = [ 3.0,9.0]
;        yrange = [-6.0,6.0]
       
        PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
        PLOT,wall.x,wall.y,color=Truecolor('Yellow'),  $
             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
        PLOT,x,y,color=Truecolor('Orange'),  $
             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,PSYM=3,/NOERASE
     ENDELSE




; Need to :
;   - add the dome perhaps so that the 'floating structure' feature is included
;   - have a quick try with the MAST-U equilibirum

;  19/11/2010
;   - required intelligence
;     -region switching, i.e. from SOL to SOL_INNER, etc.
;     -a decrease in the number of intersections when moving outward, as more of the contour is outside the wall
;     -local refinement along the contour of the tangency point region as the _adjust parameter get small


      cont = 1
      WHILE (cont) DO BEGIN
        cont = 0

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
            ENDIF ELSE  $
              PLOT,[result.x],[result.y],color=Truecolor('Orange'), PSYM=6,  $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

            IF (first_intersection) THEN BEGIN
              first_intersection = 0  ; Can I replace this for some check on the nature of int_i?
              int_i    = i
              int_x    = result.x
              int_y    = result.y
            ENDIF ELSE BEGIN
              int_i    = [int_i,i]
              int_x    = [int_x,result.x]
              int_y    = [int_y,result.y]
            ENDELSE

          ENDIF

        ENDFOR ; Loop over contour looking for wall intersections

        print,'      int_i:',int_i      
  
        int_n = N_ELEMENTS(int_i)

        ; Check if there are multiple wall intersections very close together, and 
        ; spatially refine the contour locally.  This creates overhead as the
        ; contour is passed through (at least) twice, but can't see a more general
        ; way of doing this:
        FOR j = 0, int_n-2 DO BEGIN
          IF (int_i[j+1]-int_i[j] LE 2) THEN BEGIN
            i1 = MAX([0            ,int_i[j]-1 ]    )
            j1 = MAX([0            ,int_i[j]-10]    )
            i2 = MIN([int_i[j+1]+2 ,N_ELEMENTS(x)-1])
            j2 = MIN([int_i[j+1]+10,N_ELEMENTS(x)-1])
            frac = (DINDGEN(9.0) + 1.0D) / 10.0D
            IF (ABS(x[int_i[j]]-x[int_i[j+1]]) LT  $
                ABS(y[int_i[j]]-y[int_i[j+1]])) THEN BEGIN
              new_y = y[i1] + frac * (y[i2] - y[i1])
              new_x = SPLINE( y[j1:j2], x[j1:j2], new_y ,/DOUBLE)              
            ENDIF ELSE BEGIN
 print,x
 print,y
            PLOT,x,y,color=Truecolor('Orange'),  $
                 XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,PSYM=3,/NOERASE
              print,'not ready'
              STOP
            ENDELSE

            PLOT,[new_x],[new_y],color=Truecolor('Red'), PSYM=6,  $
                 XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

            x = [x[0:i1],new_x,x[i2:N_ELEMENTS(x)-1]]
            y = [y[0:i1],new_y,y[i2:N_ELEMENTS(y)-1]]

            cont = 1  ; Loop through the contour again to refine intersection data
            print,'cycling...'
          ENDIF
        ENDFOR

      ENDWHILE ; Wall intersection loop


      region = SOL

      CASE region OF
        SOL: BEGIN
          
          store_contour = 0

          IF (int_nlast[region] EQ -1 ) THEN BEGIN
            int_nlast[region] = int_n
            store_contour = 1
          ENDIF ELSE BEGIN
            PRINT,'     psi_adjust 1:',psi_adjust,int_n,int_nlast[region]
            PRINT,'       check 1: ',int_nlast[region] LT int_n AND psi_adjust LT 0.0D
            PRINT,'       check 2: ',int_nlast[region] GE int_n AND psi_adjust GT 0.0D
            PRINT,'       check 3: ',psi_val LE psi_end AND int_n EQ 2
            IF (int_nlast[region] LT int_n AND psi_adjust LT 0.0D) THEN BEGIN
              psi_adjust = -0.3D * psi_adjust
            ENDIF ELSE BEGIN          
              IF (int_nlast[region] GE int_n AND psi_adjust GT 0.0D) THEN BEGIN
                psi_adjust = -0.3D * psi_adjust
              ENDIF ELSE BEGIN
                ; Exit condition:
                IF (psi_val LE psi_end AND int_n EQ 2) THEN store_contour = 1
              ENDELSE
            ENDELSE
          ENDELSE

          PRINT,'     psi_adjust 2:',psi_adjust,int_n,int_nlast[region]

;         Search for tangency point complete, store the contour:
          IF (psi_adjust GT 0.0D AND psi_adjust LT tol_adjust) THEN store_contour = 1
 
          print,'store_contour: ',store_contour

          IF (store_contour) THEN BEGIN

            CASE int_n OF 
;            CASE int_n-int_nlast[region] OF 
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
                PRINT,'ERROR grid_AnalyseBoundary: Tangency point ill-defined'
                STOP
                END
            ENDCASE

            ; Define radial boundary associated with the symmetry point:
            p1 = [0.0D,0.0D]
            p2 = p1
            IF (psi_adjust GT 0.0D AND psi_adjust LT tol_adjust) THEN BEGIN
               vx = -(y[tangent_i] - y[tangent_i-1])
              vy =   x[tangent_i] - x[tangent_i-1]
              length = SQRT(vx^2 + vy^2)
              p1 = [x[tangent_i],y[tangent_i]]
              p2 = [p1[0] + 0.5D / length * vx,  $
                    p1[1] + 0.5D / length * vy]

              PLOT,[p1[0],p2[0]],[p1[1],p2[1]],color=Truecolor('Lightgreen'),  $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
            ENDIF

            print, 'storing contour...'
            contour_n = contour_n + 1
            contour_data = {                  $
; need to record the contour index that gave birth to this one...
              state        : 0                 ,  $
              region       : region            ,  $
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
              y            : y         }

            name = 'contour' + STRING(contour_n,FORMAT='(I0)')
            IF (contour_n EQ 1) THEN  $
              contour_array = CREATE_STRUCT(              name,contour_data) ELSE  $
              contour_array = CREATE_STRUCT(contour_array,name,contour_data)

            IF ((psi_adjust GT 0.0D AND psi_adjust LT tol_adjust) OR  $
                 psi_val LE psi_end) THEN BEGIN

              ; Loop through contours and find those with a tangecy point,
              ; but where the radial progression of the grid on either side of the 
              ; tangecy point (in the poloidal direction) has not been processed yet: 

              n = N_ELEMENTS(TAG_NAMES(contour_array))
              status = 0

              FOR state = 0, 1 DO BEGIN
                FOR i = 1, n DO BEGIN  
                  tag = 'contour' + STRING(i,FORMAT='(I0)')
print,tag
                  cnt = grid_ExtractStructure(contour_array,tag)

                  IF (cnt.tangent_i EQ -1) THEN CONTINUE

                  help,cnt,/struct

                  IF (state EQ cnt.state) THEN BEGIN
                    CASE cnt.state OF
                      0: BEGIN
                        status = 1
                        cnt.state = 1
                        psi_val = cnt.psi
                        psi_start = psi_val
                        int_nlast[region] = 2
                        boundary1_p1 = cnt.boundary1_p1
                        boundary1_p2 = cnt.boundary1_p2
                        boundary2_p1 = cnt.tangent_p1
                        boundary2_p2 = cnt.tangent_p2
                        END
                      1: BEGIN
                        status = 1
                        cnt.state = 2
                        psi_val = cnt.psi
                        psi_start = psi_val
                        int_nlast[region] = 2
                        boundary1_p1 = cnt.tangent_p1
                        boundary1_p2 = cnt.tangent_p2
                        boundary2_p1 = cnt.boundary2_p1
                        boundary2_p2 = cnt.boundary2_p2
                        END
                      ELSE: BEGIN
                        PRINT,'ERROR grid_AnalyseBoundary: Unknown contour state'
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

              psi_adjust = psi_step * 0.1D
              print,int_x
              print,int_y

              if (contour_n EQ 7) then stop
            ENDIF

          ENDIF
          END
        ELSE: BEGIN
          PRINT,'ERROR grid_AnalyseBoundary: Unrecognized grid region'
          STOP
          END
      ENDCASE

    ENDFOR  ; Contour segment loop

    psi_val = psi_val + psi_adjust

    IF (psi_val GT psi_start) THEN BEGIN
      PRINT,'PROTCECTING PSI_VAL! ---------------'
      psi_val = psi_val + 100.0D * tol_adjust
 stop
    ENDIF

  ENDWHILE  ; PSIN loop
;  ENDFOR  






;  a = grid_ExtractContour(psi, psi_x, psi_y, psib+1.0)
;  i = sort(a.dist)
;;  print,a.dist[i]
;  i = WHERE(a.dist GT 0.10)
;  i = i[0]
;  OPLOT,a.x[0:i],a.y[0:i],color=Truecolor('Green')

;  a = grid_ExtractContour(psi, psi_x, psi_y, psib+2.0)
;  i = sort(a.dist)
;;  print,a.dist[i]
;  i = WHERE(a.dist GT 0.10)
;  i = i[0]
;  OPLOT,a.x[0:i],a.y[0:i],color=Truecolor('Green')
;  OPLOT,a.x     ,a.y     ,color=Truecolor('Green')

;  a = grid_ExtractContour(psi, psi_x, psi_y, psib-0.5)
;  i = sort(a.dist)
;;  print,a.dist[i]
;  i = WHERE(a.dist GT 0.10)
;  i = i[0]
;  OPLOT,a.x[0:i],a.y[0:i],color=Truecolor('Red')

;  a = grid_ExtractContour(psi, psi_x, psi_y, psib-1.0)
;  i = sort(a.dist)
;;  print,a.dist[i]
;  i = WHERE(a.dist GT 0.10)
;  i = i[0]
;  OPLOT,a.x[0:i],a.y[0:i],color=Truecolor('Red')


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
    CONTOUR, b_pol, b_x, b_y, NLEVELS=20
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
    CONTOUR, b_pol, b_x, b_y, NLEVELS=20
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

  result = CREATE_STRUCT( b,'null_n',n,'null_i',min_i,'null_j',min_j )

  IF (KEYWORD_SET(debug)) THEN BEGIN
    CONTOUR, b_pol, b_x, b_y, NLEVELS=20
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

; orange - o-point
; purple - primary x-point
; green  - secondary x-point

  xpoint_zone = [4.0,7.0,-3.8,5.5]
  b = grid_ReadEQUFile('/home/ITER/lisgos/divimp/shots/iter/1514/Baseline2008-li0.70.x4.equ');,/CHANGE_SIGN)
  b = grid_FindNullPoints (b,xpoint_zone,1);,/debug)
  b = grid_AnalyseBoundary(b,xpoint_zone,1,/debug)
;shade_surf,b.psi_raw,b.x,b.y; ,ax=0  

  


; xpoint_zone = [0.4,1.0,-1.5,1.5]
;; b = grid_ReadEQUFile('~/fuse_data/mast/shots/24867/24867_335.equ')
; b = grid_ReadEQUFile('~/divimp/shots/mast/24860/carre.24860_240.equ')
;;shade_surf,b.psi_raw,b.x,b.y,ax=0
;;stop
;  b = grid_FindNullPoints (b,xpoint_zone,1,/debug)
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