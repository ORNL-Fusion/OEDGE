;; PURPOSE
;; mapping a regular gridded surface to a new regular grid using
;; interpolation
;;
;; INPUT
;;
;;   data 2D (nxm))  array of data values
;;   oldx n elements vector with old xaxis values
;;   oldy m elements vector with old y values
;;   newx n* element vector with new x values
;;   newy m* element vector with new y values
;;
;; KEYWORDS 
;;  cubic parameter for cubic interpolation see IDL manuel (default
;;  -0.5)
;;  noextrap do not extrabolate if values are out of range
;;  log interpolate on logarithmic scale
;;
;; OUTPUT
;;   2D array (n*xm*)with new values 
;;
;; HISTORY 
;; 16/05/2000 creation
;;
function map2d,data,oldx,oldy,newx,newy,cubic=cubic,noextrap=noextrap,log=log
  
                                ;ensure double precision
  data = double(data)
  oldx=double(oldx)
  oldy=double(oldy)
  newx=double(newx)
  newy=double(newy)
                                ; first check dimensions
  s = size(data) 
  n = n_elements(oldx)
  m = n_elements(oldy)
  IF NOT ARG_PRESENT(cubic) THEN cubic=-0.5
  
  IF (s[0] NE 2) THEN BEGIN
    print,'ERROR: no valid input arrays'
    return,0
  ENDIF 
  
  IF (s[1] NE n) OR (s[2] NE m) THEN BEGIN
    print,'ERROR: axis do not match'
    print,s,', n: ',n,', m: ',m
    return,0
  ENDIF    
  
                                ; next check the x range
  IF (MAX(oldx) LT MAX(newx)) OR (MIN(oldx) GT MIN(newx)) THEN BEGIN
    IF KEYWORD_SET(noextrap) THEN BEGIN
      print,'ERROR: new x range values out of range -> exiting'
      return,0
    ENDIF ELSE BEGIN
      print,'WARNING: new x values out of range -> setting to nearest value!'
      print,MIN(oldx),MAX(oldx),format = '("Old range: [",F12.6,", ",F12.6,"]")'
      print,MIN(newx),MAX(newx),format = '("New range: [",F12.6,", ",F12.6,"]")'
    ENDELSE
  ENDIF
                                ; next check the y range
  IF (MAX(oldy) LT MAX(newy)) OR (MIN(oldy) GT MIN(newy)) THEN BEGIN
    IF KEYWORD_SET(noextrap) THEN BEGIN
      print,'ERROR: new y range values out of range -> exiting'
      return,0
    ENDIF ELSE BEGIN
      print,'WARNING: new y values out of range -> setting to nearest value!'
      print,MIN(oldy),MAX(oldy),format = '("Old range: [",F12.6,", ",F12.6,"]")'
      print,MIN(newy),MAX(newy),format = '("New range: [",F12.6,", ",F12.6,"]")'
    ENDELSE
  ENDIF
  
                                ; finding the indizes using linear interpolation
  IF n_elements(oldx) EQ 1 THEN BEGIN
      IF oldx NE newx THEN BEGIN
          print,'ERROR: cannot interpolate from one x element'
          return,0
      ENDIF ELSE x = oldx
  ENDIF ELSE x = INTERPOL(dindgen(n),oldx,newx)
  IF n_elements(oldy) EQ 1 THEN BEGIN
      IF oldy NE newy THEN BEGIN
          print,'ERROR: cannot interpolate from one x element'
          return,0
      ENDIF ELSE y = oldy
  ENDIF ELSE y = INTERPOL(dindgen(m),oldy,newy)
  
                                ; now map the data using cubic
                                ; interpolation
  IF KEYWORD_SET(log) THEN BEGIN
    return,  exp(INTERPOLATE(alog(data),x,y,/grid,cubic=cubic))
  ENDIF ELSE BEGIN
    return, INTERPOLATE(data,x,y,/grid,cubic=cubic)
  ENDELSE
  
END
  
  
  
