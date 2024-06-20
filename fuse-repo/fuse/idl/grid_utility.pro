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
;    HELP, struct.contour12,/struct
    result = -1
  ENDELSE

  RETURN, result
END
;
; ======================================================================
;
FUNCTION grid_GetIndex, ctrs, tag, value

  tags = STRUPCASE(TAG_NAMES(ctrs))
  nctr = N_ELEMENTS(tags)   

  result = -1

  FOR ictr = 1, nctr DO BEGIN
    ctr = grid_ExtractStructure(ctrs,tags[ictr-1])
    tags2 = STRUPCASE(TAG_NAMES(ctr))
    i = WHERE(tags2 EQ STRUPCASE(tag), count)
    IF (count NE 1) THEN CONTINUE

    IF (ctr.(i) EQ value) THEN result = ictr 
  ENDFOR

  RETURN, result
END
;
; ======================================================================
;
FUNCTION grid_GetFrac, x, y, dist=dist

  n = N_ELEMENTS(x)
  dist = MAKE_ARRAY(n,/DOUBLE,VALUE=0.0D0)
  FOR i = 1, n-1 DO  $
    dist[i] = dist[i-1] + SQRT( (x[i] - x[i-1])^2 + (y[i] - y[i-1])^2 )

;  FOR i = 1, n-1 DO  $
;    print,'shit',x[i],y[i], SQRT( (x[i] - x[i-1])^2 + (y[i] - y[i-1])^2 ), dist[i]

  frac = dist / dist[n-1]

  result = frac

  RETURN, result
END

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
; ======================================================================
;
;
;
FUNCTION grid_PerpDistance, x, y, px, py, s=min_s, index=index, nostop=nostop

  min_s = -1.0D
  index = -1

  n = N_ELEMENTS(px)  
  min_dist = 1.0D+6
  FOR i = 0, n-2 DO BEGIN
    dist = PNT_LINE( [x,y], [px[i],py[i]], [px[i+1],py[i+1]], p1)

    IF (ABS(px[i]-px[i+1]) GT 1.0D-6) THEN  $
      s = (p1[0] - px[i]) / (px[i+1] - px[i]) ELSE  $
      s = (p1[1] - py[i]) / (py[i+1] - py[i])
         
    IF (s GE 0.0D AND s LT 1.0D AND dist LT min_dist) THEN BEGIN
      min_s    = s
      min_dist = dist
      index    = i
    ENDIF
;    print, 'perp',s,dist,min_dist
;    print, '    ',px[i],py[i],'    ',px[i+1],py[i+1]
;    print, '    ',x,y
  ENDFOR

  IF (min_dist EQ 1.0D+6) THEN BEGIN
    IF (KEYWORD_SET(nostop)) THEN BEGIN
    ENDIF ELSE BEGIN
      PRINT, 'ERROR grid_PerpDistance: Failure'
      STOP
    ENDELSE
  ENDIF

  result = min_dist

  RETURN, result
END
;
; ======================================================================
;
