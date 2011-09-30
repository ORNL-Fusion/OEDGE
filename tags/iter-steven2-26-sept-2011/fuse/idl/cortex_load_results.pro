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
  signal = inGetData('SIGNAL_01')
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

  n        = inGetData('N_SIGNAL')
  z        = inGetData('ATOMIC_NUMBER')
  a        = inGetData('ATOMIC_MASS')
  charge   = inGetData('CHARGE')
  database = inGetData('DATABASE')
  wlngth   = inGetData('WAVELENGTH')
  i        = inGetData('I')
  j        = inGetData('J')
  x1       = inGetData('X1')
  y1       = inGetData('Y1')
  z1       = inGetData('Z1')
  x2       = inGetData('X2')
  y2       = inGetData('Y2')
  z2       = inGetData('Z2')

  signal = MAKE_ARRAY(N_ELEMENTS(i),n,/DOUBLE,VALUE=0.0D)
  FOR k = 1, n DO BEGIN
    tag = 'SIGNAL_'+STRING(k,FORMAT='(I02)')
    signal[*,k-1] = inGetData(tag)
  ENDFOR

  inCloseInterface  

  result = {  $
    version       : 1.0     ,  $
    file          : file    ,  $
    n_signal      : n       ,  $
    xindex        : i       ,  $
    yindex        : j       ,  $
    atomic_number : z       ,  $
    atomic_mass   : a       ,  $
    charge        : charge  ,  $ 
    wavelength    : wlngth  ,  $ 
    database      : database,  $
    signal        : signal  ,  $
    x1            : x1      ,  $
    y1            : y1      ,  $
    z1            : z1      ,  $
    x2            : x2      ,  $
    y2            : y2      ,  $
    z2            : z2      }

  RETURN,result
END
;
; ======================================================================
;