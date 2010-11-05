;
; ======================================================================
;
FUNCTION cortex_LoadProfile, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  path   = inGetData('PATH')
  signal = inGetData('SIGNAL')
  inCloseInterface  

  result = {  $
    version : 1.0  ,  $
    file    : file ,  $
    path    : path ,  $
    signal  : signal }

  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadIntegrals, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  i      = inGetData('i')
  j      = inGetData('j')
  signal = inGetData('data')
  inCloseInterface  

  result = {  $
    version : 1.0  ,  $
    file    : file ,  $
    xindex  : i      ,  $
    yindex  : j      ,  $
    signal  : signal }

  RETURN,result
END
;
; ======================================================================
;