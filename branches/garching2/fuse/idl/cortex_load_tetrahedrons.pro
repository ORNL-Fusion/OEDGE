;
;
; ======================================================================
;
FUNCTION cortex_LoadTetrahedronData, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  x      = inGetData('X')  
  y      = inGetData('Y')  
  z      = inGetData('Z')  
  dens   = inGetData('NE')  
  te     = inGetData('TE')
  n_D    = inGetData('N_D')
  n_D2   = inGetData('N_D2')
  dalpha = inGetData('D_ALPHA')

  inCloseInterface

  i = STRPOS(file,'/',/REVERSE_SEARCH)

  result = {          $
    desc    : 'tetrahedron grid centroid data' ,  $
    file    : STRMID(file,i+1) ,                  $
    version : 1.0  ,  $
    x       : x    ,  $
    y       : y    ,  $
    z       : z    ,  $
    dens    : dens ,  $
    te      : te   ,  $
    n_D     : n_D  ,  $
    n_D2    : n_D2 ,  $
    dalpha  : dalpha }
  RETURN, result
END
;
; ======================================================================
;
