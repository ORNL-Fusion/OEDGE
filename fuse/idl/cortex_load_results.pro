;
; ======================================================================
;
FUNCTION cortex_LoadIntegrals, file

  inOpenInterface, file
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