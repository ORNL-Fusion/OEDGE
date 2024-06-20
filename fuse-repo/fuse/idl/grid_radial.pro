;
; ======================================================================
;
FUNCTION grid_LoadAnnotations, code, file


;print,code
;print,file


  CASE code OF
;   --------------------------------------------------------------------
    1: BEGIN

      fp = 3
      FREE_LUN, fp
      OPENR, fp, file, error=err
      IF (err NE 0) THEN BEGIN
        PRINT,'ERROR grid_LoadAnnotationData: Unable to open data file
        PRINT,' FILE=',file
        RETURN, -1
      ENDIF

      buffer = ' '
      first_line = 1

      WHILE NOT EOF(fp) DO BEGIN 
        READF,fp,buffer      
        buffer = STRTRIM(buffer,2)
        ; print,'buffer:',buffer
        IF (STRPOS(buffer,'#') EQ 0 OR STRPOS(buffer,'*') EQ 0) THEN CONTINUE
        str = STRSPLIT(buffer,' '+STRING(9B),/EXTRACT)  ; STRING(9B) is the ASCII TAB character
        IF (N_ELEMENTS(str) EQ 1) THEN CONTINUE         ; weak check to avoid blank lines
        IF (first_line) THEN BEGIN
          x = FLOAT(str[0])
          y = FLOAT(str[1])
        ENDIF ELSE BEGIN
          x = [x,FLOAT(str[0])]
          y = [y,FLOAT(str[1])]
        ENDELSE
        first_line = 0
      ENDWHILE

      FREE_LUN, fp

      result = {x : x, y : y}

      END
;   --------------------------------------------------------------------
    2: BEGIN
      file = grid_UpdateFile(file)
      status = inOpenInterface(file)
      IF (status LT 0) THEN BEGIN
        result = CREATE_STRUCT('version',0.0,'file','none')
        RETURN, result
      ENDIF      
      x1 = inGetData('X_1')
      y1 = inGetData('Y_1')
      x2 = inGetData('X_2')
      y2 = inGetData('Y_2')
      inCloseInterface  
      result = {x1 : x1, y1 : y1, x2 : x2, y2 : y2}
      END
;   --------------------------------------------------------------------
    3: BEGIN  ; EIRENE triangles
      file = grid_UpdateFile(file)
      status = inOpenInterface(file)
      IF (status LT 0) THEN BEGIN
        result = CREATE_STRUCT('version',0.0,'file','none')
        RETURN, result
      ENDIF      
      v1 = inGetData('VERTEX_1')
      v2 = inGetData('VERTEX_2')
      v3 = inGetData('VERTEX_3')
;      s1 = inGetData('SURFACE_1')
;      s2 = inGetData('SURFACE_2')
;      s3 = inGetData('SURFACE_3')
      x  = inGetData('VERTEX_X')
      y  = inGetData('VERTEX_Y')
      inCloseInterface  
      v = [[v1],[v2],[v3]]
      s = -1
;      s = [[s1],[s2],[s3]]
      result = { version : 1.0 , file : file, v : v, surface : s, x : x , y : y }
      END
;   --------------------------------------------------------------------
    ELSE:
  ENDCASE
  
  RETURN, result
END

;
; ======================================================================
;
FUNCTION grid_LoadAnnotationData, plot,file

;  annotate = { version : 1.0 } 
;  status_annotate = -1
  annotate_data = -1
  FOR i = 0, plot.annotate_n-1 DO BEGIN
    CASE plot.annotate_code[i] OF
      1: annotate_data = grid_LoadAnnotations(plot.annotate_code[i],       plot.annotate_file[i]) 
      2: annotate_data = grid_LoadAnnotations(plot.annotate_code[i],file + plot.annotate_file[i]) 
      3: annotate_data = grid_LoadAnnotations(plot.annotate_code[i],file + plot.annotate_file[i]) 
      ELSE: BEGIN
        PRINT, 'ERROR grid_LoadAnnotationData: Unknown annotation code'
        PRINT, '  PLOT TAG = ',plot.tag
        PRINT, '  CODE     = ',plot.annotate_code[i]
        annotate_data = -1
        END
    ENDCASE
    name = 'data' + STRING(i+1,FORMAT='(I0)')
    IF (i EQ 0) THEN annotate = CREATE_STRUCT(         name,annotate_data)  ELSE  $
                     annotate = CREATE_STRUCT(annotate,name,annotate_data)
;    status_annotate = 0
  ENDFOR 

  IF (plot.annotate_n EQ 1) THEN annotate = CREATE_STRUCT(annotate,'dummy',0)  ; Can remove in future...

  result = annotate

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
FUNCTION grid_Intersection, p1, p2, seg1, seg2, mode, status=status, debug=debug, show=show
  ; --------------------------------------------------------------------
  COMMON zone, zone_list, zone_index, zone_nx  , zone_ny  ,  $
               zone_minx, zone_maxx , zone_miny, zone_maxy
  ; --------------------------------------------------------------------

  a1 = p1
  a2 = p2

  IF (ABS(a1[0]-a2[0]) LT 1.0D-10) THEN a2[0] = a2[0] + 1.0D-10  ; changed 03/10/11 SL
;  IF (ABS(a1[0]-a2[0]) LT 1.0D-6) THEN a2[0] = a2[0] + 0.0000001D
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

  IF (KEYWORD_SET(show)) THEN print, 'going on now...'

  min_s = MAKE_ARRAY(10,/DOUBLE,VALUE=-1.0D)
  min_t = min_s
  min_i = MAKE_ARRAY(10,/LONG  ,VALUE=-1L  )
  min_c = -1

  FOR j = 0, N_ELEMENTS(seg_list)-1 DO BEGIN

    i = seg_list[j]        

    b1 = seg1[*,i]
    b2 = seg2[*,i]

    IF ( ABS(b1[0]-b2[0]) LT 1.0D-10 ) THEN b1[0] = b1[0] + 1.0D-10
;    IF ( ABS(b1[0]-b2[0]) LT 1.0D-6 ) THEN b1[0] = b1[0] + 0.0000001D
;    IF ( ABS(b1[1]-b2[1]) LT 1.0D-6 ) THEN b1[1] = b1[1] + 0.00001D
  
    Lint, a1, a2, b1, b2, c, d, flag=flag
  
    IF (flag NE 1) THEN CONTINUE
  
    s = (c[0] - a1[0]) / (a2[0] - a1[0])    
    t = (c[0] - b1[0]) / (b2[0] - b1[0])    

    IF (KEYWORD_SET(show)) THEN  $
      print,s,t,FORMAT='(2F20.10,I6,2X,4F15.7)',KEYWORD_SET(int_i),a1,a2,b1,b2

;    IF (ABS(t) LT min_t) THEN BEGIN
;      min_t = t
;      min_s = s
;    ENDIF

;    IF ( s GE  0.0D     AND s LE 1.0D   AND  $
;         t GE -1.0D-10  AND t LT 1.0D ) THEN BEGIN

    IF ( s GE 0.0D AND s LE 1.0D) THEN BEGIN
      IF (ABS(t) LT 1.0D-5 OR ABS(1.0D0-t) LT 1.0D-5) THEN BEGIN
        min_c++
        min_i[min_c] = i
        min_s[min_c] = s
        min_t[min_c] = t
      ENDIF

;      IF ( t GE 0.0D AND t LE 1.0D ) THEN BEGIN
      IF ( t GE 0.0D AND t LT 1.0D ) THEN BEGIN  ; -changed 01/09/2011, SL

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
  
    ENDIF
  ENDFOR

  IF (min_c NE -1) THEN BEGIN
    help,min_s
    print,'min_s',min_s[WHERE(min_s NE -1.0D)],FORMAT='(A,10D16.10)'
    print,'min_t',min_t[WHERE(min_t NE -1.0D)],FORMAT='(A,10D16.10)'
  ENDIF
  
;print,'set',KEYWORD_SET(int_x)

  IF (KEYWORD_SET(int_x)) THEN BEGIN
    ; Check that the multiple intersections are legit:
    IF (N_ELEMENTS(int_x) EQ 2) THEN BEGIN
      IF (int_i[0]+1 EQ int_i[1]) THEN BEGIN
print,'special check for multiple intersections...'
print,int_i
print,int_s
print,int_t
;        plot,seg1[0,*],seg1[1,*],color=Truecolor('Black')
;        oplot,seg1[0,*],seg1[1,*],color=Truecolor('Black'),PSYM=7
;        oplot,[p1[0],p2[0]],[p1[1],p2[1]],color=Truecolor('Lightgreen')

        i = int_i[0]
        sa = SQRT( (seg1[0,i  ]-seg2[0,i  ])^2 + (seg1[1,i  ]-seg2[1,i  ])^2 )
        sb = SQRT( (seg1[0,i  ]-seg2[0,i+1])^2 + (seg1[1,i  ]-seg2[1,i+1])^2 )
        sc = SQRT( (seg1[0,i+1]-seg2[0,i+1])^2 + (seg1[1,i+1]-seg2[1,i+1])^2 )

        cos_b = (sa^2 + sc^2 - sb^2) / (2.0D * sa * sc)
 
        IF ((1.0D0 - ABS(cos_b)) GT 1.0E-06) THEN  $
          theta = ABS(ACOS(cos_b) * 180.0D / DOUBLE(!PI)) ELSE  $
          theta = 180.0D    
   
        ; Check that the intersections are both very close to the common vertex:
        dist = MAX([1.0-int_t[0],int_t[1]])

        print,'theta,dist',theta,dist

        IF (theta GT 179.0 OR dist LT 1.0D-5) THEN BEGIN
          IF (ABS(0.5-int_s[0]) LT ABS(0.5-int_s[1])) THEN i = 0 ELSE i = 1
          int_i = [int_i[i]]
          int_s = [int_s[i]]
          int_t = [int_t[i]]
          int_x = [int_x[i]]
          int_y = [int_y[i]]
        ENDIF
      ENDIF
    ENDIF 
    status = 1  ; Can only be 0 or 1 since used as a logical variable below
    result = {          $
      i    : int_i   ,  $
      s    : int_s   ,  $
      t    : int_t   ,  $ 
      x    : int_x   ,  $
      y    : int_y   }          
  ENDIF ELSE BEGIN
    IF (min_c GT -1) THEN BEGIN
      ; Need to check for intersections that are missed because the line passes
      ; too close to one of the wall segment end points:
      print, 'here',min_c
      print,min_i[0:min_c]
      print,min_s[0:min_c]
      print,min_t[0:min_c]
      CASE min_c OF
;       ----------------------------------------------------------------
        0: BEGIN
          status = 1
          IF (min_t[0] LT 0.5D0) THEN b1 = seg1[*,min_i[0]] ELSE b1 = seg2[*,min_i[0]]
          result = {             $
            i    : min_i[0]   ,  $
            s    : min_s[0]   ,  $
            t    : min_t[0]   ,  $ 
            x    : (b1[0])[0] ,  $
            y    : (b1[1])[0] }
          END
;       ----------------------------------------------------------------
        1: BEGIN
          status = 1
          IF (min_t[0]-1.0D LT -min_t[1]) THEN i = 0 ELSE i = 1
          result = {                   $
            i    : min_i[i]         ,  $
            s    : min_s[i]         ,  $
            t    : min_t[i]         ,  $ 
            x    : seg2[0,min_i[0]] ,  $
            y    : seg2[1,min_i[0]] }
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN
;         No valid intersection found between the line and segment list:
          status =  0 
          result = -1
          END
;       ----------------------------------------------------------------
      ENDCASE
    ENDIF ELSE BEGIN
;     No valid intersection found between the line and segment list:
      status =  0 ; Can only be 0 or 1 since used as a logical variable below
      result = -1
    ENDELSE
  ENDELSE

  ;help,result

  IF (KEYWORD_SET(show)) THEN print,'status=',status

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
      IF (KEYWORD_SET(debug) AND i EQ 0 AND j EQ 0) THEN BEGIN
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


    OPLOT,wall.x,wall.y,color=Truecolor('Black'), PSYM=6
    print,zone_index[0:1,i,j]
    stop
  ENDIF

END
;
; ======================================================================
;
FUNCTION grid_LoadWall, file ; , debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  debug  = option.debug
  xrange = option.xrange
  yrange = option.yrange

  wall = grid_LoadAnnotations(1,file)

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
    PLOT,wall.x,wall.y,color=Truecolor('Black'),XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1
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
PRO grid_Debug, b, contour_array, wall, xrange, yrange, no_stop=no_stop
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  pretty = 0

  IF (pretty) THEN BEGIN
    !P.BACKGROUND=Truecolor('Darkgrey')
    !P.THICK=2.0
    !X.THICK=2.0
    !Y.THICK=2.0
    !X.CHARSIZE=2.0
    !Y.CHARSIZE=2.0
  ENDIF

  PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Darkgrey')

  IF (NOT pretty) THEN BEGIN

    OPLOT,wall.x,wall.y,color=Truecolor('Black')

    OPLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Orange')
    OPLOT,wall.pt1[0,*],wall.pt1[1,*],color=Truecolor('Orange'), PSYM=7

    i = WHERE(wall.ptc NE -1, count) 
    IF (count GT 0) THEN  $
      OPLOT,wall.pt1[0,i],wall.pt1[1,i],color=Truecolor('Darkgrey'),PSYM=6
    
    i = WHERE(wall.ptt EQ 1, count) 
    IF (count GT 0) THEN  $
      OPLOT,wall.pt1[0,i],wall.pt1[1,i],color=Truecolor('Pink'), PSYM=7
    
    i = WHERE(wall.ptt EQ 2, count) 
    IF (count GT 0) THEN  $
      OPLOT,wall.pt1[0,i],wall.pt1[1,i],color=Truecolor('Magenta'), PSYM=7
    
    i = WHERE(wall.ptt EQ 3, count) 
    IF (count GT 0) THEN  $
      OPLOT,wall.pt1[0,i],wall.pt1[1,i],color=Truecolor('Lightgreen'),PSYM=7
  ENDIF ELSE BEGIN
    OPLOT,wall.x,wall.y,color=Truecolor('Blue')
  ENDELSE

  tags = STRUPCASE(TAG_NAMES(contour_array))

  FOR i = 1, N_ELEMENTS(tags) DO BEGIN

    name = 'contour' + STRING(i,FORMAT='(I0)')

;    print,name,i, N_ELEMENTS(tags)

    cnt = grid_ExtractStructure(contour_array,name)  

    IF (NOT pretty) THEN BEGIN
      print,i,cnt.x[0],cnt.y[0]
      XYOUTS, [cnt.x[0]], [cnt.y[0]], CHARSIZE=1.1, ALIGNMENT=0.5, $
              STRTRIM(STRING(i),2), COLOR=TrueColor('Darkgrey')
    ENDIF

    CASE cnt.region OF
      param.CORE         : color = 'Blue'
      param.SOL          : color = 'Red'
      param.SOL_LFS      : color = 'Pink'
      param.SOL_HFS      : color = 'Pink'
      param.PFZ          : color = 'Orange'
      param.PFZ_SECONDARY: color = 'Orange'
      ELSE: BEGIN
        PRINT,'ERROR grid_Debug: Unrecognized REGION'
        PRINT,' REGION = ',ctr.region
        END
    ENDCASE
  
    OPLOT,cnt.x,cnt.y,color=Truecolor(color)

    IF (NOT pretty) THEN OPLOT,cnt.x,cnt.y,color=Truecolor('Darkgrey'), PSYM=3

    IF (NOT pretty AND cnt.tangent_i NE -1) THEN  $
      OPLOT,[cnt.x[cnt.tangent_i]],[cnt.y[cnt.tangent_i]],  $
            color=Truecolor('Lightgreen'), PSYM=6

    IF (NOT pretty AND cnt.origin EQ 1 AND i EQ N_ELEMENTS(tags)) THEN BEGIN
      OPLOT,[cnt.boundary1_p1[0],cnt.boundary1_p2[0]],  $
            [cnt.boundary1_p1[1],cnt.boundary1_p2[1]],  $
            color=Truecolor('Lightgreen')
      OPLOT,[cnt.boundary2_p1[0],cnt.boundary2_p2[0]],  $
            [cnt.boundary2_p1[1],cnt.boundary2_p2[1]],  $
            color=Truecolor('Lightgreen')
    ENDIF

  ENDFOR


  IF (NOT pretty) THEN BEGIN
    OPLOT, [b.x[b.null_i[0]]], [b.y[b.null_j[0]]], PSYM=6, color=TrueColor('Red')      ; o-point
    IF (N_ELEMENTS(b.null_i) GT 1) THEN  $
      OPLOT, [b.x[b.null_i[1]]], [b.y[b.null_j[1]]], PSYM=6, color=TrueColor('Green')  ; primary   x-point
    IF (N_ELEMENTS(b.null_i) GT 2) THEN  $
      OPLOT, [b.x[b.null_i[2]]], [b.y[b.null_j[2]]], PSYM=6, color=TrueColor('Blue' )  ; secondary x-point
  ENDIF

;  print,'saving contour_array'
;  SAVE,filename='contour_array.sav',contour_array


  print,'done in debug',n_elements(tags)

  IF (NOT KEYWORD_SET(no_stop)) THEN stop

END
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
                           psi_1st_xpoint, psi_2nd_xpoint, user_step, user_finish, geometry
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  ; Loop through contours and find those with a tangecy point,
  ; but where the radial progression of the grid on either side of the 
  ; tangecy point (in the poloidal direction) has not been processed yet: 

  n = N_ELEMENTS(TAG_NAMES(contour_array))
  status = 0

  psi_end = psi_1st_xpoint - option.ctr_boundary[param.SOL-1]  ;  user_finish  parameter (need to set region)

  FOR state2 = 0, 1 DO BEGIN
    FOR i = 1, n DO BEGIN  
      tag = 'contour' + STRING(i,FORMAT='(I0)')

      cnt = grid_ExtractStructure(contour_array,tag)

      IF (cnt.tangent_i EQ -1) THEN CONTINUE

      print,'TAG:-----------> ',state2,' ',tag,cnt.state


;     help,cnt,/struct

      IF (state2 EQ cnt.state) THEN BEGIN
        CASE cnt.state OF
          0: BEGIN
            status = 1
            cnt.state = 1
            psi_val = cnt.psi
            psi_start = psi_val
            direction = cnt.direction
            region    = cnt.region
            int_nlast = 2 ; [region] = 2

            ; Special case for the first contour on a limiter grid:
            IF (n EQ 1 AND geometry EQ param.LIMITED) THEN BEGIN
              cnt.state = 2
            ENDIF ELSE BEGIN
              boundary1_p1 = cnt.boundary1_p1
              boundary1_p2 = cnt.boundary1_p2
              boundary2_p1 = cnt.tangent_p2 + 1.01D * (cnt.tangent_p1 - cnt.tangent_p2) 
              boundary2_p2 = cnt.tangent_p2
            ENDELSE

            IF (region EQ param.SOL AND direction EQ param.BACKWARD) THEN BEGIN
              print,'need to check the region' 
              stop
              IF (process_2nd EQ 3) THEN psi_end = psi_2nd_xpoint ELSE $
                                         psi_end = psi_1st_xpoint 
           ;  *** NEED TO KEEP TRACK OF WHICH END POINT RAN A-FOUL AND 
           ;      USE THAT INFORMATION HERE TO DECIDE WHAT TO DO ***
            ENDIF
            IF (region EQ param.PFZ_SECONDARY) THEN psi_begin = psi_2nd_xpoint ELSE  $
                                                    psi_begin = psi_1st_xpoint
            psi_end = psi_begin - option.ctr_boundary[region-1] 
            END
          1: BEGIN
            status = 1
            cnt.state = 2
            psi_val = cnt.psi
            psi_start = psi_val
            direction = cnt.direction
            region    = cnt.region
            int_nlast = 2 ; [region] = 2
            boundary1_p1 = cnt.tangent_p2 + 1.01D * (cnt.tangent_p1 - cnt.tangent_p2) 
;            boundary1_p1 = cnt.tangent_p1 ; temp
            boundary1_p2 = cnt.tangent_p2
            boundary2_p1 = cnt.boundary2_p1
            boundary2_p2 = cnt.boundary2_p2
            IF (region EQ param.SOL AND direction EQ param.BACKWARD) THEN BEGIN
              print, 'need to check region, i.e. _LFS or _HFS'
              stop
              direction = param.FORWARD
              psi_end = psi_1st_xpoint - option.ctr_boundary[param.SOL-1]  ;  user_finish  parameter (need to set region)
              boundary1_p1 = [-999.0D,-999.0D]
              boundary1_p2 = [-999.0D,-999.0D]
            ENDIF
            IF (region EQ param.PFZ_SECONDARY) THEN psi_begin = psi_2nd_xpoint ELSE  $
                                                    psi_begin = psi_1st_xpoint
            psi_end = psi_begin - option.ctr_boundary[region-1] 
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
PRO grid_RefineContour, x, y, i, n_pts=n_pts  ; ,  $
$                        debug=debug, xrange=xrange, yrange=yrange
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  debug  = option.debug
  xrange = option.xrange
  yrange = option.yrange

  IF (N_ELEMENTS(x) LT 3) THEN BEGIN
    IF (N_ELEMENTS(x) EQ 2) THEN BEGIN
      x = [x[0],0.5D*(x[0]+x[1]),x[1]]
      y = [y[0],0.5D*(y[0]+y[1]),y[1]]
      i = 1
      IF (N_ELEMENTS(n_pts) GT 0) THEN n_pts = 3
    ENDIF ELSE BEGIN
      PRINT,'ERROR grid_RefineContour: Well, this routine is being asked to refine a contour'
      PRINT,'                          that only has one point, which is very silly'
      STOP
    ENDELSE
  ENDIF

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

  a = MAX(ref_x,imax,SUBSCRIPT_MIN=imin)
  a = MAX(ref_y,jmax,SUBSCRIPT_MIN=jmin)

;  print,'imax',imax,N_ELEMENTS(ref_x)-1
;  print,'jmax',jmax,N_ELEMENTS(ref_y)-1

  debug_local = 0

  IF (debug_local) THEN BEGIN
    print, '-----REFINING-----'
    print,'new_x',new_x
    print,'new_y',new_y
    print,'deltax,y',ABS(x[i]-x[i+1]),ABS(y[i]-y[i+1])
  ENDIF

  force = -1
  IF ( ((imax NE 0 AND imax NE N_ELEMENTS(ref_x)-1) OR  $
        (imin NE 0 AND imin NE N_ELEMENTS(ref_x)-1)) AND  $
       ( jmax EQ 0 OR  jmax EQ N_ELEMENTS(ref_y)-1)  AND  $
       ( jmin EQ 0 OR  jmin EQ N_ELEMENTS(ref_y)-1) ) THEN force = 1

  IF (debug_local) THEN BEGIN
    print,'imax,min  = ',imax,imin,N_ELEMENTS(ref_x)-1
    print,'jmax,min  = ',jmax,jmin,N_ELEMENTS(ref_y)-1
    print,'force = ',force
  ENDIF

  indep_x = 0
  indep_y = 0
  ecount  = 0

  store_ref_x = ref_x
  store_ref_y = ref_y
  store_new_x = new_x
  store_new_y = new_y

  FOR ipass = 1, 2 DO BEGIN

    IF ((force EQ -1 AND ABS(x[i]-x[i+1]) LT ABS(y[i]-y[i+1])) OR  $
         force EQ  1) THEN BEGIN
      indep_y = 1
      IF (debug_local) THEN print,'indep = y'
      IF (ref_y[0] GT ref_y[1]) THEN BEGIN
        new_y = REVERSE(new_y)
        ref_x = REVERSE(ref_x)
        ref_y = REVERSE(ref_y)
        status = 1
        IF (debug_local) THEN print,'reversing y!'
      ENDIF ELSE status = 0
      j = SORT(ref_y)
      ref_x = ref_x[j]
      ref_y = ref_y[j]
      new_x = SPLINE( ref_y, ref_x, new_y ,/DOUBLE)              
    ENDIF ELSE BEGIN
      indep_x = 1
      IF (debug_local) THEN print,'indep = x'
      IF (ref_x[0] GT ref_x[1]) THEN BEGIN
        new_x = REVERSE(new_x)
        ref_x = REVERSE(ref_x)
        ref_y = REVERSE(ref_y)
        status = 1
        IF (debug_local) THEN print,'reversing x!'
      ENDIF ELSE status = 0
      j = SORT(ref_x)
      ref_x = ref_x[j]
      ref_y = ref_y[j]
      new_y = SPLINE( ref_x, ref_y, new_x ,/DOUBLE)              
    ENDELSE

    IF (status) THEN BEGIN
      new_x = REVERSE(new_x)
      new_y = REVERSE(new_y)

      ref_x = REVERSE(ref_x)
      ref_y = REVERSE(ref_y)
    ENDIF

    IF (debug_local) THEN BEGIN
      print, '------------------',ipass
      print,'ref_x',ref_x
      print,'ref_y',ref_y
      print, '-'
      print,'new_x',new_x
      print,'new_y',new_y
    ENDIF
;
;   ----------------------------------------------------------------------
;   CHECK THAT THE INTERPOLATED POINTS ARE OK
;
    ; Calcualte the perpendicular distance between the interpolated points
    ; and a line segment between the bounding range of the independent
    ; variable:
    dist = 0.0D
    FOR j = 0, N_ELEMENTS(new_x)-1 DO  $
      dist = [dist,grid_PerpDistance(new_x[j],new_y[j],  $
                     [x[i1],x[i2]], [y[i1],y[i2]],/nostop)]  
    ; Check that these distances are not creater than 20% of these end
    ; points, which indicates that the interpolation failed:
    limit = SQRT( (x[i1]-x[i2])^2 + (y[i1]-y[i2])^2 )
    IF (MAX(dist) GT 0.2D * limit) THEN BEGIN
      ecount++
      IF (debug_local) THEN BEGIN
        PRINT,'limit ',limit
        PRINT,'dist ',dist 
        IF (ecount EQ 1) THEN BEGIN
          xspan = 2.0D0 * limit
          yspan = xspan
    
          PLOT, x,y,COLOR=Truecolor('Red'),  $
                XRANGE=[0.5D*(x[i1]+x[i2])-xspan,0.5D*(x[i1]+x[i2])+xspan],  $
                YRANGE=[0.5D*(y[i1]+y[i2])-xspan,0.5D*(y[i1]+y[i2])+xspan],  $
                XSTYLE=1,YSTYLE=1
          OPLOT, ref_x,ref_y,COLOR=Truecolor('Darkgrey'),PSYM=3
          OPLOT, [x[i1],x[i2]],[y[i1],y[i2]],COLOR=Truecolor('Black'),PSYM=7
          OPLOT, new_x,new_y,COLOR=Truecolor('Green'),PSYM=6
        ENDIF
      ENDIF
      IF (ecount EQ 1) THEN BEGIN
;       Restore the interpolation reference data and try again:
        ref_x = store_ref_x
        ref_y = store_ref_y
        new_x = store_new_x
        new_y = store_new_y
        IF (indep_x) THEN force = 1 ELSE force = 2
        ;PRINT, 'PROBLEM DETECTED'
;        STOP
      ENDIF
    ENDIF ELSE BREAK

  ENDFOR ; IPASS loop

  IF (ecount EQ 2) THEN BEGIN
    IF (debug_local) THEN  $
      OPLOT, new_x,new_y,COLOR=Truecolor('Pink'),PSYM=6

    ; Out of desperation, remove any curvature between the points, which
    ; seems to be the best one can do when the interpolation fails for
    ; both passes (x-independent and y-independent):
    new_x = store_new_x
    new_y = store_new_y

    IF (debug_local) THEN BEGIN
      OPLOT, new_x,new_y,COLOR=Truecolor('Blue'),PSYM=6
      STOP
    END
  ENDIF

  IF (ecount EQ 1) THEN BEGIN
    IF (debug_local) THEN  $
      OPLOT, new_x,new_y,COLOR=Truecolor('Silver'),PSYM=6
    ;STOP
  ENDIF

  OPLOT,[new_x],[new_y],color=Truecolor('Red'), PSYM=3

  x = [x[0:i1],new_x,x[i2:N_ELEMENTS(x)-1]]
  y = [y[0:i1],new_y,y[i2:N_ELEMENTS(y)-1]]

END
;
; ======================================================================
;
PRO grid_ClipToWall, x, y, wall_pt1, wall_pt2, cycle


  cycle = 0

  cont  = 1
  WHILE (cont) DO BEGIN
    cont = 0

    ; Find intersections with wall starting at beginning of the wall:

    FOR i = 0, N_ELEMENTS(x)-2 DO BEGIN
;print, 'start of contour',i,N_ELEMENTS(x)
      result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                 wall_pt1, wall_pt2, 1, status=status)
;      IF (status) THEN print, 'result.i',result.i,status
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
      grid_RefineContour, x, y, i ; ,  $
;                          debug=debug, xrange=xrange, yrange=yrange
      cont = 1
    ENDELSE

    ; Find intersections with wall starting at the end of the contour:
    FOR i = N_ELEMENTS(x)-2, 0, -1 DO BEGIN
;print, 'end of contour',N_ELEMENTS(x)
      result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                 wall_pt1, wall_pt2, 1, status=status)
;      IF (status) THEN print, 'result.i',result.i,status
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
      grid_RefineContour, x, y, i ; ,  $
;                          debug=debug, xrange=xrange, yrange=yrange
      cont = 1
;       oplot,result.x,result.y,psym=6,color=Truecolor('Red')
;       help,result,/struct
;       stop

    ENDELSE
  ENDWHILE

  ; Extend the end points out a little further, to avoid problems when the end point is almost
  ; right on the wall, and then the wall is shifted a bit for one reason or another:  -SL, 13/06/2013
;  x[0] = 10.0D * (x[0] - x[1]) + x[1]
;  y[0] = 10.0D * (y[0] - y[1]) + y[1]
;  n = N_ELEMENTS(x) - 1
;  x[n] = 10.0D * (x[n] - x[n-1]) + x[n-1]
;  y[n] = 10.0D * (y[n] - y[n-1]) + y[n-1]


      oplot,x,y,psym=3,color=Truecolor('Purple')

END
;
; ======================================================================
;
FUNCTION grid_LocatePoint,x,y,p1,s=s,t=t

;  print,'point location search',N_ELEMENTS(x)

  ; Find where the point is along the contour:
  n = N_ELEMENTS(x)
  FOR i = 0, n-2 DO BEGIN
    IF (ABS(x[i+1]-x[i]) GT 1.0D-8) THEN s = (p1[0] - x[i]) / (x[i+1] - x[i]) ELSE s = -999.0D
    IF (ABS(y[i+1]-y[i]) GT 1.0D-8) THEN t = (p1[1] - y[i]) / (y[i+1] - y[i]) ELSE t = -999.0D
;    print,'locate',i,n,s,t,x[i+1]-x[i],y[i+1]-y[i],p1[0],x[i],p1[1],y[i],FORMAT='(A,2I6,2F18.12,2F18.12,2X,4F18.12)'
    IF (s EQ -999.0D AND ABS(p1[0]-x[i]) LT 1.0D-8) THEN s = t
    IF (t EQ -999.0D AND ABS(p1[1]-y[i]) LT 1.0D-8) THEN t = s
;    print,'      ',i,n,s,t,x[i+1]-x[i],y[i+1]-y[i],FORMAT='(A,2I6,2F18.12,2F18.12)'
    IF ((ABS(s-t) LT 1.0D-7 AND s GE 0.0D AND s LE 1.0D)) THEN BREAK
;        (s EQ -999.0D AND t EQ -999.0D0                )) THEN BREAK 

;    IF ((ABS(s-t) LT 1.0D-7 AND s LT 0.9999999D) OR  $
;        (s EQ -999.0D AND (t GE -1.0D-10 AND t LT 1.0D)) OR  $
;        (t EQ -999.0D AND (s GE -1.0D-10 AND s LT 1.0D)) OR  $
;        (s EQ -999.0D AND t EQ -999.0D0                )) THEN BREAK 
;;        (s EQ -999.0D AND (t GE 0.0D AND t LT 1.0D)) OR  $   
;;        (t EQ -999.0D AND (s GE 0.0D AND s LT 1.0D))) THEN BREAK 
  ENDFOR
  IF (i EQ n-1) THEN BEGIN
    PRINT, 'ERROR grid_LocatePoint: Location of point on line not found'
    STOP
  ENDIF

  result = i

  RETURN, result
END
;
; ======================================================================
;
FUNCTION grid_GetOrthogonal, x_val,y_val,p1_val,direction,length,ictr=ictr,ctrs=ctrs,span=span
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
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

; the problem, now, naturally, is with the ITER grid, where there's a tangency point that faces
; away from the separatrix... can't exactly see what to do with it just yet...

      ; Find out which contour has the cross-field vector that's associated 
      ; with the secondary null: 

      FOR i = 1, nctr DO BEGIN
        ctr = grid_ExtractStructure(ctrs,tags[i-1])      

;        print,ctr.separatrix,ctr.tangent_p1[0]

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
  i = grid_LocatePoint(x,y,p1,s=s,t=t)

  print, 'i  ',i
  print, 'p1 ',p1

  val1 = ABS(s     ) ; had to move these outside the IF parentheses for some reason...
  val2 = ABS(s-1.0D)

  IF ( val1 GT 1.0D-7 AND val2 GT 1.0D-7) THEN BEGIN
    ; This is a non-standard situation at present, so some of the logic below
    ; needs to be modified to make this case valid.
    val1 = ABS(t     ) ; had to move these outside the IF parentheses for some reason...
    val2 = ABS(t-1.0D)
    IF (val1 LT 1.0D-7 OR val2 LT 1.0D-7) THEN s = t ELSE BEGIN
      PRINT, 'ERROR grid_GetOrthogonal: Normal vector requested away from vertex'
      PRINT, ' S,T = ',s,t
      STOP
    ENDELSE
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
        
      print,'i i_2nd',i,i_2nd
      print,'sep    ',separatrix
      print,'focus_x',focus_x
      print,'focus_y',focus_y
      print,'x,y[i    ]',x[i    ],y[i    ]
      print,'x,y[i_2nd]',x[i_2nd],y[i_2nd]

      OPLOT,x,y,COLOR=Truecolor('Lightgreen')
      OPLOT,[focus_x],[focus_y],PSYM=6,COLOR=Truecolor('Black')

      OPLOT,[x[i    ]],[y[i    ]],PSYM=6,COLOR=Truecolor('Purple')
      OPLOT,[x[i_2nd]],[y[i_2nd]],PSYM=6,COLOR=Truecolor('Orange')

      PRINT, 'NOT READY, SORRY MAN'
;stop

      RETURN, [-999.0D,-999.0D]
;      stop
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
    IF (1 EQ 1 AND frac LT threshold) THEN BEGIN
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
; How do I avoid the distorted target cells that plagued the JET grid, i.e. the 
; cells overlapped?  Just keep the resolution high? 
;
; Need to bring in the other wall segments now I think, and look for intersections, rather 
; than checking if a point is inside or outside (although good to have some appropriate
; spatial refinement before this step? to make sure the entrances and exits are all accounted for?)
;
;
FUNCTION grid_AnalyseBoundary, b, wall, user_step, user_finish, machine, $
                               save=save, mswin=mswin, limiter=limiter
  ; ------------------------------------------------------------------
  COMMON commons, param, option, state
  ; ------------------------------------------------------------------

  debug  = option.debug
  xrange = option.xrange
  yrange = option.yrange

  IF (debug) THEN BEGIN
    PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
    OPLOT,wall.x,wall.y,color=Truecolor('Black')
    OPLOT, [b.x[b.null_i[0]]], [b.y[b.null_j[0]]], PSYM=6, color=TrueColor('Red')     ; o-point
    IF (b.null_n GE 2) THEN  $
      OPLOT, [b.x[b.null_i[1]]], [b.y[b.null_j[1]]], PSYM=6, color=TrueColor('Green')   ; primary   x-point
    IF (b.null_n GE 3) THEN  $
      OPLOT, [b.x[b.null_i[2]]], [b.y[b.null_j[2]]], PSYM=6, color=TrueColor('Blue' ) ; secondary x-point
  ENDIF

  region = param.SOL

  direction = param.FORWARD

  psi   = b.psi ; psi_raw
  psi_x = b.x
  psi_y = b.y

  ; Set lower or upper null:
  IF (b.null_n GE 2) THEN BEGIN
    IF (b.y[b.null_j[1]] LT b.y[b.null_j[0]]) THEN geometry = param.LSN  $  ; Just label as LSN OR USN for now, independent of whether
                                              ELSE geometry = param.USN     ; or not the secondary null is incorporated into the grid.
  ENDIF ELSE geometry = param.LIMITED

  psi_1st_xpoint = b.psi_1st_xpoint

  IF (geometry EQ param.LIMITED) THEN BEGIN
    psi_val   = -option.ctr_boundary[param.CORE-1] * 0.1D ; hard-coded !!! TEMP !!!  
    psi_start = psi_val
    psi_end   = -option.ctr_boundary[param.SOL -1]   
  ENDIF ELSE BEGIN
    ; psi_val   = DOUBLE(psi_1st_xpoint) * 1.0001D ;  changed for JET, not sure what happens to the other grids... sigh... - 13/06/2013  0.9D ; MAST * 0.9999999D ; MAST
    psi_val = DOUBLE(psi_1st_xpoint) ; changed for EAST grids... a better system I hope... -SL, 19/11/2013
   
    psi_start = psi_val
    psi_end   = psi_1st_xpoint - option.ctr_boundary[param.SOL-1]  ;  user_finish  parameter (need to set region)
  ENDELSE
  psi_step   = user_step * direction
  psi_adjust = psi_step
  tol_adjust = 1.0D-7 ; 1.0D-7

  ; Decide what to do about the secondary x-point:
  process_2nd     = -1
  process_2nd_pfz =  0
  psi_span        = 0.0D0
  psi_2nd_xpoint  = -999.0D

  IF (b.null_n GE 3) THEN BEGIN
    xpt2_x = b.x[b.null_i[2]] 
    xpt2_y = b.y[b.null_j[2]] 
    ; Check whether or not the non-primary x-points are inside the vessel wall:
    psi_2nd_xpoint = b.psi_2nd_xpoint
    psi_span       = ABS(psi_1st_xpoint - psi_2nd_xpoint)
    IF (psi_2nd_xpoint GT psi_end) THEN BEGIN
      ; If the secondary x-point is inside the main wall, flag that processing
      ; is required:
      focus_x = b.x[b.null_i[2]]
      focus_y = b.y[b.null_j[2]]
      inside = grid_PointInPolygon(focus_x,focus_y,wall.x,wall.y)
      i = WHERE(inside EQ 1, count)
      IF (count GT 0) THEN BEGIN
        process_2nd     = 0  ; Tell code to process the secondary separatrix
        process_2nd_pfz = 1
        dist = grid_PerpDistance(focus_x,focus_y,wall.x,wall.y,/nostop)
      ENDIF ELSE BEGIN
        ; Double check that the null isn't ever-so-just outside the vessel:
        dist = grid_PerpDistance(focus_x,focus_y,wall.x,wall.y,/nostop)
        IF (dist EQ 1.0D+6) THEN BEGIN
          ; The calculation of the perpendicular distance between the secondary null and the wall
          ; has failed, perhaps due to the point not having a valid projection onto a wall 
          ; segment because it's outside the boundary rectangle that contains all the wall points,
          ; so make a simple estimate:
          dist = MIN(SQRT( (wall.x-focus_x)^2 + (wall.y-focus_y)^2 ))
          PRINT, 'WARNING grid_AnalyseBoundary: Perpendicular distance not valid, so calculating'
          PRINT, '                              the minumum distance from the wall points to the'
          PRINT, '                              secondary null'
          PRINT, '  DISTANCE (m) = ',dist
        ENDIF
        IF (dist LT 0.01D) THEN BEGIN
          PRINT, 'WARNING grid_AnalyseBoundary: Secondary null is outside the vessel, but is very'
          PRINT, '                              close to the wall, so is being brought inside'
          PRINT, '  DISTANCE (m) = ',dist
          process_2nd     = 0  ; Note, the secondary separatrix won't necessarily be added later,
          process_2nd_pfz = 1  ; depending on if a tangency point is found nearby on the wall
        ENDIF
      ENDELSE
    ENDIF 
  ENDIF

  ; Setup local wall segment arrays:
  wall_pt1 = wall.pt1
  wall_pt2 = wall.pt2

;
; ----------------------------------------------------------------------
; INITIALIZATION
;

  int_nlast = -1 ; MAKE_ARRAY(7,/LONG,VALUE=-1)

  contour_n = 0

  boundary1_p1 = [-999.0D,-999.0D]
  boundary1_p2 = [-999.0D,-999.0D]
  boundary2_p1 = [-999.0D,-999.0D]
  boundary2_p2 = [-999.0D,-999.0D]

  save_boundary1_p1 = [-999.0D,-999.0D]
  save_boundary1_p2 = [-999.0D,-999.0D]
  save_boundary2_p1 = [-999.0D,-999.0D]
  save_boundary2_p2 = [-999.0D,-999.0D]

  save_2nd_boundary2_p1 = [-999.0D,-999.0D]
  save_2nd_boundary2_p2 = [-999.0D,-999.0D]

  the_search_continues   = 1
  active_contour         = 0
  store_forced           = 0
  just_passed_2nd_xpoint = 0
  processing_core        = 0

  IF (0 EQ 1) THEN BEGIN
    xrange = [ 3.0,9.0]
    yrange = [-6.0,6.0]
    xrange = [ 4.0,5.5]
    yrange = [ 4.0,5.0]
    PLOT,wall.x,wall.y,color=Truecolor('Black'),  $
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
  IF (geometry EQ param.LIMITED) THEN process_pfz = -1 ELSE process_pfz = 0
  failure_2nd_pfz = 0
  careful_2nd = 0
  careful_relax = 0
  careful_min = 1.0E+6
  immediate_trouble = 0
  local_count = 0
;  last_iseg = -1  ; keeping track of the last valid contour segment


;  silly_count = 0

  second_side = 0

;  be_careful = 0
;
; ----------------------------------------------------------------------
; ----------------------------------------------------------------------
;
; MAIN LOOP
;  
; ----------------------------------------------------------------------
; ----------------------------------------------------------------------
;
  WHILE (the_search_continues) DO BEGIN

    iteration_count = iteration_count + 1
    local_count     = local_count     + 1


    IF (debug) THEN BEGIN
      print,' '
      print,'===========',contour_n,iteration_count,local_count,'============='
      print,'psi_values: ',psi_start,psi_val,psi_end,process_2nd,direction,psi_step,psi_adjust,careful_2nd, careful_relax,region,  $
                           FORMAT='(A,3F17.9,2I6,2F12.7,3I6)'
      print,'psi_end,processing_core=',psi_end,processing_core
    ENDIF


;   IF (iteration_count EQ 5) THEN STOP
; print,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',psi_val,psi_2nd_xpoint,psi_1st_xpoint,psi_val GT psi_2nd_xpoint


    ; Check if the search contour has crossed the secondary separatrix:
    IF (psi_val LT psi_2nd_xpoint AND last_psi_val GT psi_2nd_xpoint) THEN BEGIN
      IF (process_2nd EQ 0) THEN process_2nd = 1
      focus_x = b.x[b.null_i[2]]
      focus_y = b.y[b.null_j[2]]
      proximity = MIN(SQRT( (last_x-focus_x)^2 + (last_y-focus_y)^2 ))
      IF (proximity LT option.ctr_2nd_xpt_proximity AND  $
          save_2nd_boundary2_p1[0] EQ -999.0D) THEN BEGIN
        just_passed_2nd_xpoint = 1
        save_2nd_boundary2_p1  = boundary2_p1
        save_2nd_boundary2_p2  = boundary2_p2
      ENDIF
      IF (geometry EQ param.LSN) THEN region = param.SOL_HFS ELSE  $
                                      region = param.SOL_LFS
;      IF (geometry EQ param.LOWER_NULL) THEN region = param.SOL_HFS ELSE  $
;                                             region = param.SOL_LFS
      psi_end = psi_1st_xpoint - option.ctr_boundary[region-1]

      PRINT,'>>>>>>>>>>>>> PSI BOUNDARY PAST 2nd X-POINT <<<<<<<<<<<<<<<<',region,psi_end

    ENDIF


    ; Generate the contour for the current value of PSI:
    ctr = grid_ExtractContour(psi, psi_x, psi_y, psi_val)


    ibrk = WHERE(ctr.dist GT option.ctr_break_distance) 
    nseg = N_ELEMENTS(ibrk) 
    ibrk = [-1,ibrk,ctr.n-1]

    ; Initialize the flag that registers when at least one segment of 
    ; the contour is inside the current region of interest:
    seg_valid = 0  
;
;   --------------------------------------------------------------------
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

      IF (geometry EQ param.LIMITED AND contour_n GE 1) THEN BEGIN

        PRINT,'before clip'
        ;PRINT,'x:',x
        ;PRINT,'y:',y

        IF (0 EQ 1 AND iseg EQ 1) THEN BEGIN
        
          print,ibrk[iseg],ibrk[iseg+1]
          print,ibrk[iseg]+1 EQ ibrk[iseg+1]
          print,ibrk[iseg]   EQ ibrk[iseg+1]
    
          xrange = [ 3.3,4.0]  ; option.xrange
          yrange = [-0.5,1.1]  ; [0.0,0.6]  ; option.yrange

          PLOT,x,y,color=Truecolor('Darkgrey'), PSYM=6, $ 
               XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1
          PLOT,[boundary1_p1[0],boundary1_p2[0]],  $
               [boundary1_p1[1],boundary1_p2[1]],  $
               color=Truecolor('Lightgreen'),  $
               XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
          PLOT,[boundary2_p1[0],boundary2_p2[0]],  $
               [boundary2_p1[1],boundary2_p2[1]],  $
               color=Truecolor('Green'),  $
               XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
          PLOT,wall.x,wall.y,color=Truecolor('Black'),  $
               XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
          print,'... nseg',nseg
    
          print,'psi:',psi_val
          print,'1st',psi_1st_xpoint
          print,'2nd',psi_2nd_xpoint
          print,'contact_x,y',contact_x,contact_y


        ENDIF

          inside = grid_PointInPolygon(x,y,wall.x,wall.y)


          i = WHERE(inside EQ 0, count)

          proximity = SQRT ( (x[i]-contact_x)^2 + (y[i]-contact_y)^2 )
          dummy = MIN(proximity,imin)

;          OPLOT,[x[i[imin]]],[y[i[imin]]],color=Truecolor('Blue'),PSYM=6

          j = i[imin]
      
        
          i = WHERE(inside EQ 1, count)
;          OPLOT,[x[i]],[y[i]],color=Truecolor('Pink'), PSYM=6


          hold_x = x
          hold_y = y

          n = N_ELEMENTS(x)

          IF (j GT 0) THEN BEGIN
            x = [ hold_x[j:n-1] , hold_x[0:j-1], hold_x[j] ]
            y = [ hold_y[j:n-1] , hold_y[0:j-1], hold_y[j] ]
          ENDIF
    
;          print, 'n',n
;          print, N_ELEMENTS(x)

;          print,x
;          print,y

;          print,contact_x
;          print,contact_y

;          stop
;        ENDIF

      ENDIF

      ; Isolate region of interest along the contour by clipping at the poloidal 
      ; boundaries:

      caution = 0
      store_x = x
      store_y = y
      cont    = 1
      first_pass = 1
;
;     ------------------------------------------------------------------
;     DECIDE IF THIS SEGMENT IS THE ONE INSIDE THE REGION OF INTEREST
;  
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

              xrange = option.xrange
              yrange = option.yrange
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
              PLOT,wall.x,wall.y,color=Truecolor('Black'),  $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

              IF (cont EQ 3) THEN stop

              mean_x = MEAN(x)

              IF (failure EQ 1) THEN test_x = boundary1_p1[0] ELSE  $
                                     test_x = boundary2_p1[0]
;
;             ----------------------------------------------------------
;             TRY TO CORRECT FOR A DOUBLE-NULL GRID
;
              IF (N_ELEMENTS(b.null_i) GE 3) THEN BEGIN

                secondary_x = b.x[b.null_i[2]]
                print, 'testing',secondary_x,mean_x,test_x,failure,process_2nd
                ; This check is not fool proof unfortunately, since it's feasible that
                ; the tangency point for the bounday near the second null is inside
                ; the null and the correct contour is outside...
                IF ((process_2nd EQ -1 AND  $
                     (mean_x LT secondary_x AND 0.95D*test_x LT secondary_x OR     $
                     (mean_x GT secondary_x AND 1.05D*test_x GT secondary_x))) OR  $
                    (process_2nd GE  0 AND psi_val      LT psi_2nd_xpoint AND      $
                                           last_psi_val GT psi_2nd_xpoint)) THEN   $
	         cont = 3 ELSE cont = 2 ; Failure
;
;             ----------------------------------------------------------
;             TRY TO CORRECT FOR A SINGLE-NULL GRID
;
              ENDIF ELSE BEGIN
                print,'failure',failure
                print,'iseg   ',iseg ; ,last_iseg
                cont = 3
              ENDELSE 
;
;             ----------------------------------------------------------
;             TURN OFF ONE OF THE RADIAL CONTOUR BOUNDARIES DEFINED BY
;             THE TANGENCY POINTS AND TRY AGAIN
;
              IF (cont EQ 3) THEN BEGIN

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
;
;     ------------------------------------------------------------------
;     Cycle to the next segment if the current segment is not suitable:
;
      print,'cont',cont
      IF (cont EQ 2) THEN BEGIN
        IF (seg_active[iseg] EQ 1) THEN seg_active[iseg] = 0
        CONTINUE
      ENDIF 

      ; Debugging:
      CASE contour_n OF
        221: BEGIN 
          xrange = [4.0,5.0] ; top, in a bit
          yrange = [4.0,5.0] 
          END
        ELSE: BEGIN
          xrange = option.xrange
          yrange = option.yrange
          END
      ENDCASE
      PLOT, xrange, yrange ,/NODATA, XSTYLE=1, YSTYLE=1
      PLOT,wall.x,wall.y,color=Truecolor('Black'),  $
           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
      PLOT,x,y,color=Truecolor('Darkgrey'),  $
           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,PSYM=3,/NOERASE
      OPLOT, [b.x[b.null_i[0]]], [b.y[b.null_j[0]]], PSYM=6, color=TrueColor('Red')      ; o-point
      IF (b.null_n EQ 2) THEN  $
        OPLOT, [b.x[b.null_i[1]]], [b.y[b.null_j[1]]], PSYM=6, color=TrueColor('Green')  ; primary   x-point
      IF (b.null_n EQ 3) THEN  $
        OPLOT, [b.x[b.null_i[2]]], [b.y[b.null_j[2]]], PSYM=6, color=TrueColor('Blue' )  ; secondary x-point
      PLOT,[boundary1_p1[0],boundary1_p2[0]],  $
           [boundary1_p1[1],boundary1_p2[1]],  $
           color=Truecolor('Lightgreen'),  $
           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
      PLOT,[boundary2_p1[0],boundary2_p2[0]],  $
           [boundary2_p1[1],boundary2_p2[1]],  $
           color=Truecolor('Green'),  $
           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

      PLOT,x,y,color=Truecolor('Black'),  $
           XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,PSYM=3,/NOERASE

      ; Only take the parts of the contour segment that are inside the main wall, by dropping
      ; all points before and after the first and last wall intersections, respectively:       
      hold_x = x
      hold_y = y

;print,'x:',x
;print,'y:',y

      ; Needed to move this check here, instead of below the call to _ClipToWall because
      ; the segments-of-non-interest were confusing the clipping routine when tangency
      ; points were ignored.  -SL, 05/05/14 

      IF ( process_2nd EQ 1 ) THEN BEGIN
        mean_x = MEAN(x)
        ; print,'...',MEAN(x),b.x[b.null_i[2]],process_2nd
        ; print,'...looking 1A...........',mean_x,geometry EQ param.LSN
        IF ((geometry EQ param.LSN AND mean_x GT b.x[b.null_i[2]]) OR  $
            (geometry EQ param.USN AND mean_x LT b.x[b.null_i[2]])) THEN BEGIN
          print,'...looking 1 and skipping...'
          CONTINUE 
        ENDIF
      ENDIF

      grid_ClipToWall, x, y, wall_pt1, wall_pt2, cycle   ; , xrange, yrange

      last_x = x
      last_y = y

      IF (processing_core) THEN BEGIN
        distance = SQRT ( ( hold_x[0] - hold_x[N_ELEMENTS(x)-1] )^2 +  $
                          ( hold_y[0] - hold_y[N_ELEMENTS(y)-1] )^2 )
        length = grid_Length(hold_x,hold_y) 
        IF ( distance / length GT 0.10 ) THEN cycle = 1 ELSE BEGIN  $
          x = hold_x
          y = hold_y
        ENDELSE
        print, '  distance ---------------------------------',distance,length,distance/length*100.0D
;        print,x
;        print,y
      ENDIF 

      print,'CYCLE:',cycle
;
;     ------------------------------------------------------------------
;     CYCLE THE MAIN LOOP AND ADVANCE psi IF NO SEGMENTS ARE INSIDE WALL
;
      IF (cycle) THEN CONTINUE

      mean_x = MEAN(x)
      mean_y = MEAN(y)

      ; Moved above...
      ;IF ( 1 EQ 1 AND  process_2nd EQ 1 ) THEN BEGIN
      ;  print,'...',MEAN(x),b.x[b.null_i[2]],process_2nd
      ;  IF ((geometry EQ param.LSN AND mean_x GT b.x[b.null_i[2]]) OR  $
      ;      (geometry EQ param.USN AND mean_x LT b.x[b.null_i[2]])) THEN BEGIN
      ;    print,'...looking 1...........'
      ;    CONTINUE 
      ;  ENDIF
      ;ENDIF

      ; Check if the contour is radially outward of the secondary x-point, if that's
      ; the current region of interest (PROCESS_2ND = 2):
      IF ( process_2nd EQ 3 ) THEN BEGIN
;       print, 'hey, the code actually gets here'
;       stop 
;      IF ( process_2nd EQ 2 ) THEN BEGIN
        print,'...',MEAN(x),b.x[b.null_i[2]],process_2nd
;        mean_x = MEAN(x)
        IF ((active_2nd EQ param.HFS AND mean_x LT b.x[b.null_i[2]]) OR  $
            (active_2nd EQ param.LFS AND mean_x GT b.x[b.null_i[2]])) THEN BEGIN
          print,'...looking 2...........'
          CONTINUE 
        ENDIF
      ENDIF

      IF ( process_pfz EQ 1 ) THEN BEGIN
;        mean_y = MEAN(y)
        IF ((geometry EQ param.LSN AND mean_y GT b.y[b.null_j[1]]) OR  $
            (geometry EQ param.USN AND mean_y LT b.y[b.null_j[1]])) THEN BEGIN
;        IF ((geometry EQ param.LOWER_NULL AND mean_y GT b.y[b.null_j[1]]) OR  $
;            (geometry EQ param.UPPER_NULL AND mean_y LT b.y[b.null_j[1]])) THEN BEGIN
          print,'...looking 1st pfr...........'
          CONTINUE 
        ENDIF
      ENDIF

      IF ( process_pfz EQ 2 ) THEN BEGIN
;        mean_y = MEAN(y)
        IF ((geometry EQ param.LSN AND mean_y LT b.y[b.null_j[2]]) OR  $
            (geometry EQ param.USN AND mean_y GT b.y[b.null_j[2]])) THEN BEGIN
          print,'...looking 2nd pfr...........'
          CONTINUE 
        ENDIF

      ENDIF
;
;     ------------------------------------------------------------------
;     THE SEGMENT HAS "SURVIVED" AND WILL NOW BE PROCESSED
;  
      print,process_2nd,psi_val, psi_2nd_xpoint,psi_2nd_xpoint+MAX(b.psi)

      IF (geometry NE param.LIMITED) THEN BEGIN

        IF (process_2nd EQ -1 AND careful_2nd EQ 0 AND  $
            psi_val LT (        option.ctr_2nd_xpt_frac  * psi_1st_xpoint +      $
                        (1.0D - option.ctr_2nd_xpt_frac) * psi_2nd_xpoint)  AND  $
            psi_val GT psi_2nd_xpoint) THEN BEGIN
      
          dist = MIN( SQRT( (xpt2_x - x)^2 + (xpt2_y - y)^2 ) )
      
          IF (dist LT option.ctr_2nd_xpt_proximity) THEN BEGIN  ; 0.10D) THEN BEGIN  ; parameter
            careful_2nd   = 1
            careful_relax = 1
            psi_step   = psi_step   * option.ctr_2nd_xpt_slow  ; 0.05D
            psi_adjust = psi_adjust * option.ctr_2nd_xpt_slow  ; 0.05D
          ENDIF
        ENDIF

        IF (process_2nd NE -1 AND careful_2nd NE 1 AND  $
                                  careful_2nd NE 3 AND  $
            psi_val LT (        option.ctr_2nd_xpt_frac  * psi_1st_xpoint +      $
                        (1.0D - option.ctr_2nd_xpt_frac) * psi_2nd_xpoint)  AND  $
            psi_val GT psi_2nd_xpoint) THEN BEGIN

          careful_2nd++
          psi_step   = psi_step   * option.ctr_2nd_xpt_slow  
          psi_adjust = psi_adjust * option.ctr_2nd_xpt_slow
        ENDIF

      ENDIF

      IF (careful_2nd NE 0 AND psi_val LT psi_2nd_xpoint) THEN BEGIN
        careful_2nd   = 0
        careful_relax = 0
        psi_step   = psi_step   / option.ctr_2nd_xpt_slow
        psi_adjust = psi_adjust / option.ctr_2nd_xpt_slow
      ENDIF
;
;     ------------------------------------------------------------------
;     SET THE CORRECT REGION 
;  
      IF (process_pfz EQ 0 AND psi_val LT psi_2nd_xpoint) THEN BEGIN
        IF (mean_x LT b.x[b.null_i[2]]) THEN region = param.SOL_HFS ELSE  $
                                             region = param.SOL_LFS
        psi_end = psi_1st_xpoint - option.ctr_boundary[region-1]
      ENDIF

      seg_valid = 1
;
;     ------------------------------------------------------------------
;     ------------------------------------------------------------------
;     LOOP TO FIND ??? ("TANGENCY") POINTS WITH THE WALL
;  
      cont      = 1
      only_once = 1
      WHILE (cont) DO BEGIN
        cont  = 0
        int_i = -1
        first_intersection = 1

        FOR i = 0, N_ELEMENTS(x)-2 DO BEGIN  
;         Find intersections with the wall segments:
          result = grid_Intersection([x[i],y[i]], [x[i+1],y[i+1]],  $
                                     wall_pt1, wall_pt2, 1, status=status)

          IF (status) THEN BEGIN
;            print,'x,x1',[x[i],y[i]], [x[i+1],y[i+1]]

            IF (N_ELEMENTS(result.i) GT 1) THEN BEGIN

              PLOT,[result.x],[result.y],color=Truecolor('Darkgrey'), PSYM=6,  $
                   XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
              PRINT,'******** MULTIPLE SEGMENT INTERSECTIONS *********',N_ELEMENTS(result.i)

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

        IF (NOT careful_relax) THEN BEGIN ; careful_2nd NE 1 AND careful_2nd NE 3) THEN BEGIN
          ; Check if there are multiple wall intersections very close together, and 
          ; spatially refine the contour locally.  This creates overhead as the
          ; contour is passed through more than once, but I can't see a more efficent
          ; way of making this general check:
          FOR j = 0, int_n-2 DO BEGIN
            IF ((int_i[j+1]-int_i[j] LE 2) OR  $
                (int_n EQ 3 AND j EQ 1 AND only_once)) THEN BEGIN
              grid_RefineContour, x, y, int_i[j] ; ,  $
;                                  debug=debug, xrange=xrange, yrange=yrange
              cont = 1  ; Loop through the contour again to refine intersection data
              print,'refining and cycling...'
              IF (int_n EQ 3 AND j EQ 1) THEN BEGIN
                ; this can happen when concecutive contour points are vertical and the wall segment is
                ; vertical, or more generally when the contour is poorly resolved and tangential, and
                ; also when you're right at the tolerance for contour position refinement and so are
                ; ready to store the contour... this seems to sort things out, but perhaps good to add
                ; some additional checking that it's the right thing to do
                only_once = 0
                print,'handling N_INT=3 special... (new)'
              ENDIF
            ENDIF
          ENDFOR
        ENDIF ELSE BEGIN
          print,'NOT refining and cycling...'
        ENDELSE
      ENDWHILE ; Wall intersection loop
;
;     ------------------------------------------------------------------
;     ------------------------------------------------------------------
;
;     Record the fact that there are more than two intersections on the first pass,
;     which isn't good since this is a small step from the last intersection (although
;     doesn't necessarily mean trouble: 
      IF (local_count EQ 1 AND int_n NE 2) THEN BEGIN
        print, 'REGISTERING THAT THINGS ARE LOOKING SHAKY'
        immediate_trouble = immediate_trouble + 1
      ENDIF ELSE BEGIN
        IF (immediate_trouble GT 0 AND int_n EQ 2) THEN BEGIN
          print, 'EVERYTHING IS OK NOW!'
          immediate_trouble = 0
        ENDIF
      ENDELSE
;
;     ------------------------------------------------------------------
;     ANALYSE ??? ('tangecy points'), AND DECIDE IF CONTOUR IS TO BE
;     BROKEN
;  
      store_contour = 0
      cycle = 0

      IF ( int_nlast EQ -1 AND NOT processing_core ) THEN BEGIN
        PRINT,' *** INT_NLAST = -1, STORING CONTOUR ***'
        int_nlast = int_n 
        store_contour = 1
        sht_count = 0
      ENDIF ELSE BEGIN
        PRINT,'     psi_adjust 1:',psi_adjust,int_n,int_nlast ; [region]

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
;        PRINT,'       check 1: ',int_nlast[region] LT int_n AND psi_adjust LT 0.0D
;        PRINT,'       check 2: ',int_nlast[region] GE int_n AND psi_adjust GT 0.0D
;        PRINT,'       check 3: ',psi_val LE psi_end AND int_n EQ 2

        IF (int_nlast LT int_n AND psi_adjust*direction LT 0.0D) THEN BEGIN
          psi_adjust = -0.3D * psi_adjust
          sht_count = 0
        ENDIF ELSE BEGIN          
          IF (int_nlast GE int_n AND psi_adjust*direction GT 0.0D) THEN BEGIN
            psi_adjust = -0.3D * psi_adjust
            sht_count = 0
          ENDIF ELSE BEGIN
            ; Exit conditions:
            IF (((direction EQ param.FORWARD  AND psi_val LE psi_end)  OR  $
                 (direction EQ param.BACKWARD AND psi_val GE psi_end)) AND $
                int_n EQ 2) THEN BEGIN
              PRINT,' *** PSI BOUNDARY HIT, STORING CONTOUR ***'
              store_contour = 1
            ENDIF
            IF (store_forced GT 0) THEN BEGIN
              PRINT,' *** STORAGE OF THE CONTOUR IS FORCES SOMEWHERE! ***'
              store_contour = 1
            ENDIF
            length = grid_Length(x,y) 
            region_length = param.SOL
            IF (process_pfz EQ 1) THEN region_length = param.PFZ
            IF (process_pfz EQ 2) THEN region_length = param.PFZ_SECONDARY
            IF (psi_adjust EQ psi_step                                  AND  $
                length     LT option.ctr_minimum_length[region_length]) THEN BEGIN
              sht_count = sht_count + 1
              IF (sht_count GE option.ctr_minimum_count) THEN BEGIN
                PRINT,' *** STORING CONTOUR BECAUSE IT IS TOO SHORT ***'
                store_contour = 1  ; 3) THEN store_contour = 1
              ENDIF
            ENDIF ELSE  sht_count = 0

            IF (int_n EQ 1 AND int_nlast EQ 2) THEN BEGIN
              IF (last_good_psi_val NE 0.0D) THEN BEGIN
                PRINT, '*** STRANGE EVENT, PERHAPS OUT OF SPACE, TAKING LAST GOOD CONTOUR ***'
                psi_val      = last_good_psi_val
                psi_adjust   = 0.0D 
                store_forced = 1

                xrange = option.xrange
                yrange = option.yrange
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

      PRINT,'     psi_adjust 2:',psi_adjust,int_n,int_nlast,store_contour

      ; Search for tangency point complete, store the contour:
      broke_contour = 0
      IF (psi_adjust*direction GT 0.0D      AND  $
          psi_adjust*direction LT tol_adjust) THEN BEGIN
        PRINT,' *** STORING CONTOUR BECAUSE TANGENCY POINT IDENTIFIED ***'
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
        print,'need to check region specification, i.e. should be SOL_LFS or SOL_HFS'
        stop
        IF (region EQ param.SOL) THEN BEGIN
          psi_start = psi_2nd_xpoint - 0.2D 
          psi_end   = psi_1st_xpoint 
          psi_val   = psi_val - 0.1D
          direction = param.BACKWARD
          psi_step  = -0.005D * direction
          psi_adjust = psi_step
        ENDIF ELSE BEGIN
          print, 'pfz not done yet'
          STOP
        ENDELSE  
        store_contour = 0
      ENDIF
;
;     ------------------------------------------------------------------
;     ------------------------------------------------------------------
;     STORE THE CONTOUR AND MOVE ON
;
      IF (store_contour) THEN BEGIN

        IF (debug) THEN BEGIN
          print, ' '
          print, ' ------------------------------------------------------'
          print, ' |                                                    |'
          print, ' |                 STORING CONTOUR                    |'
          print, ' |                                                    |'
          print, ' ------------------------------------------------------'
          print, ' '
          print, ' int_n,broke_contour,contour_n:',int_n,broke_contour,contour_n
          print, ' '
        ENDIF
;         
;       ----------------------------------------------------------------
;       CHOSING THE TANGENCY POINT (tangent_i)
;         
        CASE int_n OF 
          2: BEGIN
             ; Do not treat the first limited contour as a tangency point;
             tangent_i = -1

             IF (geometry EQ param.LIMITED AND contour_n EQ 0) THEN  $
               tangent_i = int_i[0] + 1 + (int_i[1] - int_i[0]) / 2 

             END
          3: BEGIN
            ; Selecting the second point in the intersection list isn't perfect
            ; since you can get two intersections on the same line segment when
            ; resolving corners (I think), which results in two points being 
            ; counted, but they are almost on top of each other so it doesn't 
            ; really matter if you take int_i[1] or int_i[2] (I hope):
            tangent_i = int_i[1] + 1
            x = [x[0:int_i[1]],int_x[1],x[int_i[1]+1:N_ELEMENTS(x)-1]]
            y = [y[0:int_i[1]],int_y[1],y[int_i[1]+1:N_ELEMENTS(y)-1]]
            PLOT,[x[tangent_i]],[y[tangent_i]],color=Truecolor('Black'), PSYM=6,  $
                 XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
            END
          4: BEGIN
            ; Take an estimate of the tangency point since a single intersection point (implying
            ; that this is a 'real' tangency point between the contour and the wall) was not found:
            tangent_i = int_i[1] + 1 + (int_i[2] - int_i[1]) / 2  ; OK?
            END
          6: BEGIN
            ; Well, the only case of this happening so far is when the wall and the decrete contour 
            ; really look conformal, so make an approximation:
            print,'LOADS OF TANGENCY POINTS! TAKING A CHANCE!'
            tangent_i = LONG(MEAN(int_i[1:4]))
            print,int_i
	    print,tangent_i
;            stop
            END
          ELSE: BEGIN
       
            print,'int_n',int_n
            print,'int_i',int_i
            xspan = 0.1D
            yspan = 0.1D
            xrange = [int_x[1]-xspan,int_x[1]+xspan]
            yrange = [int_y[1]-xspan,int_y[1]+xspan]
            PLOT,[boundary1_p1[0],boundary1_p2[0]],  $
                 [boundary1_p1[1],boundary1_p2[1]],  $
                 color=Truecolor('Lightgreen'),  $
                 XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1
            PLOT,[boundary2_p1[0],boundary2_p2[0]],  $
                 [boundary2_p1[1],boundary2_p2[1]],  $
                 color=Truecolor('Green'),  $
                 XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
            PLOT,wall.x,wall.y,color=Truecolor('Black'),  $
                 XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
            PLOT,x,y,  $
                 color=Truecolor('Red'), PSYM=3, $
                 XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

            print,int_x
            print,int_y

            PLOT,int_x,int_y,  $
                 color=Truecolor('Lightblue'), PSYM=6, $
                 XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
       
            PRINT,'ERROR grid_AnalyseBoundary: Tangency point ill-defined'
            STOP
            END
        ENDCASE
;         
;       ----------------------------------------------------------------
;       DEFINING THE PERPENDICULAR ASSOICATED WITH THE TANGENCY POINT
;         
        p1 = [0.0D,0.0D]
        p2 = p1
        IF (broke_contour) THEN BEGIN
          side = 2
          IF (direction EQ param.BACKWARD) THEN BEGIN
            side = 1
          ENDIF 

          p1 = [x[tangent_i],y[tangent_i]]

          p2 = grid_GetOrthogonal(x,y,p1,side,1.0D,span=10)

          PLOT,[p1[0],p2[0]],[p1[1],p2[1]],color=Truecolor('Lightgreen'),  $
               XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

          IF (careful_relax) THEN BEGIN ; careful_2nd EQ 1 OR careful_2nd EQ 3) THEN BEGIN
            careful_2nd++
            careful_relax = 0

            dist = MIN( SQRT( (xpt2_x - x)^2 + (xpt2_y - y)^2 ) )

            ; If the tangency point is very close to the secondary null then
            ; use a non-standard contour boundary check ("parallel" to the wall
            ; at the tangency point rather than perpendicular):           

            IF (dist LT option.ctr_2nd_xpt_parallel) THEN BEGIN  ; 0.05D) THEN BEGIN  ; parameter
              careful_n  = contour_n + 1
              careful_p1 = p1
              careful_p2 = p2
              FOR frac = 0.05D, 1.0D, 0.05D DO BEGIN         
                p3 = p1 + frac * (p2 - p1)
                p4 = [p3[0],p3[1]] + [-(p2[1]-p1[1]), (p2[0]-p1[0])]
                p5 = [p3[0],p3[1]] + [ (p2[1]-p1[1]),-(p2[0]-p1[0])]

                PLOT,[p4[0],p5[0]],[p4[1],p5[1]],color=Truecolor('Lightgreen'),  $
                     XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE

                dummy = grid_Intersection([p4[0],p4[1]], [p5[0],p5[1]],  $
                                          wall_pt1, wall_pt2, 0, status=status)
                print,frac,status
                IF (NOT status) THEN BREAK
              ENDFOR
              p1 = p4
              p2 = p5
            ENDIF
          ENDIF
        ENDIF
        IF (debug) THEN BEGIN
          print, 'tangent_i =',tangent_i
          print, 'tangent_p1=',p1
          print, 'tangent_p2=',p2
        ENDIF
        IF (tangent_i NE -1 AND p1[0] EQ 0.0D and p1[1] EQ 0.0D) THEN STOP
;         
;       ----------------------------------------------------------------
;       REGISTERING THE CONTOUR DATA IN THE DATA STRUCTURE
;         
        IF (contour_n EQ 0 AND geometry EQ param.LIMITED) THEN BEGIN

            psi_1st_xpoint  = psi_val
          b.psi_1st_xpoint  = psi_val

          contact_x = p1[0]
          contact_y = p1[1]

          processing_core = 0       

          ; Reorder the points so that the contour starts and ends on the tangency point:
          i = WHERE(x EQ contact_x AND y EQ contact_y, count)
          IF (count EQ 0) THEN BEGIN
            PRINT, 'ERROR grid_AnalyseBoundary: LCFS tangency point not found in list'
            STOP
          ENDIF 
          n = N_ELEMENTS(x)
          x = [ x[i:n-1] , x[0:i-1], x[i] ]
          y = [ y[i:n-1] , y[0:i-1], y[i] ]
          tangent_i = 0            

;     print,'tangent_i',tangent_i
;     print,'n',n,N_ELEMENTS(x)
;     stop

        ENDIF

        contour_n = contour_n + 1
        contour_data = {  $
          state        : 0             ,  $
          origin       : 1             ,  $
          separatrix   : 0             ,  $
          region       : region        ,  $
          direction    : direction     ,  $
          active_2nd   : active_2nd    ,  $
          psi          : psi_val       ,  $
          adjust       : psi_adjust    ,  $ 
          tangent_i    : tangent_i     ,  $
          tangent_p1   : p1            ,  $
          tangent_p2   : p2            ,  $
          boundary1_p1 : boundary1_p1  ,  $
          boundary1_p2 : boundary1_p2  ,  $
          boundary2_p1 : boundary2_p1  ,  $
          boundary2_p2 : boundary2_p2  ,  $
          int_n        : int_n         ,  $
          int_nlast    : int_nlast     ,  $
          int_i        : int_i         ,  $
          int_x        : int_x         ,  $
          int_y        : int_y         ,  $
          x            : x             ,  $
          y            : y             }

        name = 'contour' + STRING(contour_n,FORMAT='(I0)')

;if (contour_n EQ 2) then stop

        IF (contour_n EQ 1) THEN BEGIN
          contour_array = CREATE_STRUCT(name,contour_data)
        ENDIF ELSE BEGIN
          ; Check that the contour is not redundant -- an ugly 'fix' for this problem, but
          ; I don't know what else to do a the moment, see above under 'redundant':
          check_tags  = STRUPCASE(TAG_NAMES(contour_array))
          check_nctr  = N_ELEMENTS(check_tags)          
          FOR check_ictr = 0, check_nctr-1 DO BEGIN
            check_ctr = grid_ExtractStructure(contour_array,check_tags[check_ictr])      
            length = grid_Length(x,y) 
            mean_x       = MEAN(x)
            mean_y       = MEAN(y)
            mean_check_x = MEAN(check_ctr.x)
            mean_check_y = MEAN(check_ctr.y)
            IF ( ABS( psi_val - check_ctr.psi )          LT 1.0D-6  AND  $
                 ABS( mean_x  - mean_check_x  ) / length LT 0.01    AND  $
                 ABS( mean_y  - mean_check_y  ) / length LT 0.01   ) THEN BEGIN
              PRINT, '*** REDUNDANT CONTOUR DETECTED, DISCARDING ***'
              stop

              contour_n = contour_n - 1
              BREAK
            ENDIF
          ENDFOR
          ; Store the contour:
          IF (check_ictr EQ check_nctr) THEN  $   
            contour_array = CREATE_STRUCT(contour_array,name,contour_data)
        ENDELSE

        IF (process_2nd EQ 1) THEN process_2nd = 2

        print, 'broke_contour', broke_contour
        print, 'store_forced ', store_forced
        print, 'end1         ', psi_val LE psi_end
        print, 'end2         ', psi_val GE psi_end
        print, 'sht_count    ', sht_count

        IF (broke_contour OR  $
            store_forced  OR  $
            (direction EQ param.FORWARD  AND psi_val LE psi_end) OR  $
            (direction EQ param.BACKWARD AND psi_val GE psi_end) OR  $
            sht_count GE 3) THEN BEGIN

          PRINT, 'setting active contour'

          grid_SetActiveContour, contour_array, int_nlast, region, active_contour, store_forced,  $
                                 psi_start, psi_end, psi_val, psi_adjust, psi_step, sht_count, last_good_psi_val,   $
                                 boundary1_p1, boundary1_p2, boundary2_p1, boundary2_p2,  $
                                 the_search_continues, contour_n, direction, process_2nd,  $
                                 psi_1st_xpoint, psi_2nd_xpoint, user_step, user_finish, geometry
          last_psi_val = -1.0D
          seg_active = MAKE_ARRAY(100,/LONG,VALUE=0L)
          save_boundary1_p1 = [-999.0D,-999.0D]
          save_boundary1_p2 = [-999.0D,-999.0D]
          save_boundary2_p1 = [-999.0D,-999.0D]
          save_boundary2_p2 = [-999.0D,-999.0D]

          print, 'active_contour       =', active_contour
          print, 'the search continues =', the_search_continues
          print,int_x
          print,int_y
        ENDIF

        immediate_trouble = 0
        local_count = 0
;        last_iseg = -1

        careful_2nd = 0
        careful_relax = 0
        just_passed_2nd_xpoint = 0

        ; ...
        IF (geometry EQ param.LIMITED AND contour_n EQ 1) THEN the_search_continues = 0

      ENDIF  ; STORE_CONTOUR = TRUE block
;
;     ------------------------------------------------------------------
;     ------------------------------------------------------------------
;
    ENDFOR  ; ISEG loop (list of contour segments)
;
;   --------------------------------------------------------------------
;   --------------------------------------------------------------------
;

    PRINT,'seg_valid final', seg_valid

    IF ( NOT second_side AND process_2nd EQ 3 AND NOT seg_valid ) THEN BEGIN
      ; This is here because sometimes the scan of the HFS is gets messed up
      ; depending on the relation between the centre-column tangency point
      ; and the secondary separatrix.  However, in some cases it can cause a
      ; redundant scan on the LFS, and I haven't figured out how to elimiate
      ; this yet -- as a 'temporary' measure a check for redundancy is made
      ; just before the contour is stored (ugly):
      IF ( boundary1_p1[0] NE -999.0D OR boundary2_p1[0] NE -999.0D ) THEN BEGIN
        second_side = 1
        boundary1_p1 = [-999.0D,-999.0D]
        boundary1_p2 = [-999.0D,-999.0D]
        boundary2_p1 = [-999.0D,-999.0D]
        boundary2_p2 = [-999.0D,-999.0D]
        CONTINUE
      ENDIF ELSE STOP
    ENDIF

    IF (NOT seg_valid) THEN BEGIN
      print,' '
      print,'  PROBLEM !', last_good_psi_val
      print,'           ', active_contour
      print,'           ', contour_n
      print,' '

      IF (iteration_count EQ 1 AND geometry EQ param.LIMITED) THEN processing_core = 1

      IF (processing_core) THEN BEGIN
        IF (psi_adjust EQ psi_step) THEN BEGIN
          psi_val = psi_val + psi_adjust
          CONTINUE
        ENDIF ELSE BEGIN
          IF (psi_adjust GT 0.0) THEN BEGIN
            psi_adjust = -psi_adjust * 0.3D
            psi_val = psi_val + psi_adjust
            CONTINUE
          ENDIF
        ENDELSE

      ENDIF ELSE BEGIN IF (just_passed_2nd_xpoint) THEN BEGIN
        just_passed_2nd_xpoint = 0
        boundary2_p1 = [-999.0D0,-999.0D0]  
        boundary2_p2 = [-999.0D0,-999.0D0]
        print, '------------------------here we go-------------------------'
;        stop
        CONTINUE  ; restart the main loop to avoid updating PSI_VAL

      ENDIF ELSE BEGIN IF (last_good_psi_val EQ 0.0D) THEN BEGIN
        PRINT, 'ABANDONING ACTIVE CONTOUR'
        IF (process_2nd EQ 4) THEN failure_2nd_pfz = 1
;        IF (process_2nd EQ 3) THEN failure_2nd_pfz = 1
        ; Choose a new contour and continue:
        grid_SetActiveContour, contour_array, int_nlast, region, active_contour, store_forced,  $
                               psi_start, psi_end, psi_val, psi_adjust, psi_step, sht_count, last_good_psi_val,   $
                               boundary1_p1, boundary1_p2, boundary2_p1, boundary2_p2,  $
                               the_search_continues, contour_n, direction, process_2nd,  $
                               psi_1st_xpoint, psi_2nd_xpoint, user_step, user_finish, geometry
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

        xrange = option.xrange
        yrange = option.yrange
        PLOT,x,y,color=Truecolor('Cyan'), PSYM=3, $
             XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
      ENDELSE
      ENDELSE
      ENDELSE

    ENDIF

    IF (geometry EQ param.LIMITED AND the_search_continues EQ 0 AND contour_n EQ 1) THEN BEGIN

      psi_val           = psi_val ;  + 0.01D * psi_step
      psi_adjust        = psi_step
      last_psi_val      = -1.0D
      last_good_psi_val =  0.0D

;      print,user_step
;      print,psi_1st_xpoint
;      print,psi_start
;      print,psi_val
;      print,psi_end

      the_search_continues = 1
      active_contour       = 0
      store_forced         = 0

      boundary1_p1 = [-999.0D,-999.0D]
      boundary1_p2 = [-999.0D,-999.0D]
      boundary2_p1 = [-999.0D,-999.0D]
      boundary2_p2 = [-999.0D,-999.0D]

      CONTINUE

    ENDIF




    last_psi_val = psi_val
    cont = 1
    WHILE (cont) DO BEGIN
      cont = 0
      psi_val = psi_val + psi_adjust

      IF (ABS(psi_adjust) LE tol_adjust) THEN BEGIN

        PRINT,'--- PSI SUPER KICK! ---------------'
        print,'psi_values kick: ',psi_start,psi_val,psi_end,psi_adjust,immediate_trouble

        IF (immediate_trouble) THEN BEGIN
          print,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
          print,'giving it a kick'
          print,'last_good_psi_val',last_good_psi_val          

          psi_val    = psi_start + DOUBLE(immediate_trouble) * 5.0D * psi_step 
          psi_adjust = psi_step          

          local_count = 0

          BREAK
        ENDIF ELSE BEGIN

        ENDELSE
      ENDIF

      IF ((direction EQ param.FORWARD  AND psi_val GT psi_start) OR  $
          (direction EQ param.BACKWARD AND psi_val LT psi_start)) THEN BEGIN

        PRINT,'--- PROTECTING PSI_VAL! ---------------'
        print,'psi_values protect: ',psi_start,psi_val,psi_end,psi_adjust,immediate_trouble

;        IF (active_contour EQ 16) THEN BEGIN
;          xrange = [3.7,3.9]  ; temp
;          yrange = [0.0,0.5]
;          xrange_user = xrange
;          yrange_user = yrange
;        ENDIF

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





    IF ((process_2nd EQ 1 OR process_2nd EQ 2) AND the_search_continues EQ 0) THEN BEGIN
      print, ' '
      print, ' '
      print, ' '
      print, ' '
      print, ' '
      print, ' '
      print, ' '
      print, ' >>>>>>>>>>>>> other side of 2nd x-point <<<<<<<<<<<< contour_n',contour_n
      print, ' '
      print, ' '
      print, ' '
      print, ' '
      print, ' '
      print, ' '
      print, ' '


      ; grid_Debug,b,contour_array,wall,xrange,yrange

      IF (geometry EQ param.LSN) THEN region = param.SOL_LFS ELSE  $
                                      region = param.SOL_HFS 

      psi_start         = psi_2nd_xpoint + user_step
      psi_end           = psi_1st_xpoint - option.ctr_boundary[region-1]       
      psi_val           = psi_start
      psi_adjust        = psi_step
      last_psi_val      = -1.0D
      last_good_psi_val =  0.0D

      the_search_continues = 1
      active_contour       = 0
      store_forced         = 0

      boundary1_p1 = [-999.0D,-999.0D]
      boundary1_p2 = [-999.0D,-999.0D]
      boundary2_p1 = save_2nd_boundary2_p1 
      boundary2_p2 = save_2nd_boundary2_p2 

      process_2nd = 3
    ENDIF


    IF (process_pfz EQ 0 AND the_search_continues EQ 0) THEN BEGIN
      print, ' '
      print, ' >>>>>>>>>>>>> primary pfz <<<<<<<<<<<< contour_n',contour_n
      print, ' '

      IF (process_2nd EQ 3) THEN process_2nd = 4

      psi_start         = psi_1st_xpoint - user_step
      psi_end           = psi_start - option.ctr_boundary[param.PFZ-1] ; user_finish * 0.3D ; replace 0.3 with specific pfz parameter
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

      direction = param.BACKWARD

      boundary1_p1 = [-999.0D,-999.0D]
      boundary1_p2 = [-999.0D,-999.0D]
      boundary2_p1 = [-999.0D,-999.0D]
      boundary2_p2 = [-999.0D,-999.0D]

      process_pfz = 1

      region = param.PFZ
      int_nlast = 2 ; [region] = 2  ; Rubbish, remove [region] in the code...
    ENDIF

    IF (process_pfz     EQ 1 AND  $
        process_2nd_pfz EQ 1 AND the_search_continues EQ 0) THEN BEGIN
      print, ' '
      print, ' >>>>>>>>>>>>> secondary pfz <<<<<<<<<<<< contour_n',contour_n
      print, ' '

      IF (process_2nd EQ 3) THEN process_2nd = 4

      psi_start         = psi_2nd_xpoint - user_step * 0.01D
      psi_end           = psi_start - option.ctr_boundary[param.PFZ_SECONDARY-1]  ; user_finish * 0.3D  ; replace with specific parameter
      psi_val           = psi_start
      psi_adjust        = psi_step
      last_psi_val      = -1.0D
      last_good_psi_val =  0.0D

      the_search_continues = 1
      active_contour       = 0
      store_forced         = 0

      direction = param.BACKWARD

      boundary1_p1 = [-999.0D,-999.0D]
      boundary1_p2 = [-999.0D,-999.0D]
      boundary2_p1 = [-999.0D,-999.0D]
      boundary2_p2 = [-999.0D,-999.0D]

      process_pfz = 2

      region = param.PFZ_SECONDARY
      int_nlast = 2 ; [region] = 2  ; Rubbish, remove [region] in the code...
    ENDIF

    ; One time check to register which region is active (HFS or LFS) after 
    ; passing the secondary x-point when processing the SOL:
    IF (process_2nd EQ 1 AND active_2nd EQ 0) THEN BEGIN
      mean_x = MEAN(DOUBLE(x))
      IF (mean_x LT b.x[b.null_i[2]]) THEN active_2nd = param.HFS ELSE active_2nd = param.LFS
      print, 'setting active_2nd'
      print,mean_x
      print, b.x[b.null_i[2]]
      print,active_2nd
    ENDIF


  ENDWHILE  ; PSIN loop
; 
; ---------------------------------------------------------------------- 
; ---------------------------------------------------------------------- 
;


; 
; ---------------------------------------------------------------------- 
; CORRECT THE TANGENCY SEGMENT ASSOCIATED WITH THE 'CAREFUL' CONTOUR
; THAT WAS JUST INSIDE THE SECONDARY X-POINT
;
  IF (KEYWORD_SET(careful_p1)) THEN BEGIN

    PRINT,'REPLACING THE TANGENT VECTOR NEAR 2nd X-POINT'

    tag = 'contour'+STRTRIM(STRING(careful_n),2)
    print,'updating careful contour >'+tag+'<'
    ctr = grid_ExtractStructure(contour_array,tag)      
    ctr.tangent_p1 = careful_p1
    ctr.tangent_p2 = careful_p2
    contour_array = grid_UpdateStructure(contour_array,tag,ctr)
    PLOT,[careful_p1[0],careful_p2[0]],[careful_p1[1],careful_p2[1]],color=Truecolor('Lightblue'),  $
         XRANGE=xrange,YRANGE=yrange,XSTYLE=1,YSTYLE=1,/NOERASE
  ENDIF

  print, '......all done......'

; 
; ---------------------------------------------------------------------- 
; SET THE GEOMETRY PROPERLY, i.e. WITH REFERENCE TO THE 2nd X-POINT
;


  IF (process_2nd_pfz GT 0) THEN  $
    IF (geometry EQ param.LSN) THEN geometry = param.LDN ELSE  $
                                    geometry = param.UDN

  scan_params = {  $
    geometry        : geometry        ,  $
    process_2nd     : process_2nd     ,  $
    process_2nd_pfz : process_2nd_pfz ,  $
    failure_2nd_pfz : failure_2nd_pfz }
                  
  state.geometry = geometry

print,mswin
print,mswin + '\idl\grid_data\stored_step1_'+machine+'.sav'

  IF (KEYWORD_SET(save)) THEN BEGIN
    IF (KEYWORD_SET(mswin)) THEN  $
      SAVE,filename=mswin + '\idl\grid_data\stored_step1_'+machine+'.sav',b,wall,contour_array,scan_params,state ELSE  $
      SAVE,filename='grid_data/stored_step1_'+machine+'.sav',b,wall,contour_array,scan_params,state
  ENDIF

  xrange = option.xrange
  yrange = option.yrange
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
  grid_Debug,b,contour_array,wall,xrange,yrange, /no_stop


  RETURN, result

END


;
;

