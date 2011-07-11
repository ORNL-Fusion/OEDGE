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

  frac = dist / dist[n-1]

  result = frac

  RETURN, result
END
