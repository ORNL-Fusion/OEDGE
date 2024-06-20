;; Helper function for fluxcoordinates it calculates grad f in
;; cylindrical coordinates from a given 2D flux function perpendicular 
;; on theta
;;
;; Expecting field of the form f(r,z,t) t
;; Note it calculates in (r,theta,z) coordinates
function gradcv,r,z,field
  
  sf = size(field)
  sr = size(r)
  sz = size(z)
  IF sf[0] EQ 3 THEN BEGIN
    IF (sr[0] EQ 1) AND (sr[1] EQ sf[1]) AND ( sz[0] EQ 1) AND (sz[1] EQ sf[2]) THEN BEGIN
      res = dblarr(sf[1],sf[2],sf[3],3)
      FOR j=0,sf[3]-1 DO BEGIN
; slmod begin - bug 25/01/2010
        FOR i=0,sf[2]-1 DO res[*,i,j,0] = DERIV(r,field[*,i,j])
        FOR i=0,sf[1]-1 DO res[i,*,j,2] = DERIV(z,field[i,*,j])
;
;        FOR i=0,sf[1]-1 DO res[*,i,j,0] = DERIV(r,field[*,i,j])
;        FOR i=0,sf[2]-1 DO res[i,*,j,2] = DERIV(z,field[i,*,j])
; slmod end
      ENDFOR
  ENDIF ELSE BEGIN
      print,'ERROR: wrong field dimensions'
    ENDELSE
  ENDIF ELSE BEGIN
    print,'ERROR: No valid field specified'
  ENDELSE
  
  return,res
end
