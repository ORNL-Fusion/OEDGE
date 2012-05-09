;
; ======================================================================
;
FUNCTION cortex_LoadSpectraInfo, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  index         = inGetData('INDEX')
  cell          = inGetData('CELL') 
  path          = inGetData('PATH') 
  delta         = inGetData('DELTA') 

  result = {  $
    verison       : 1.0           ,  $
    file          : file          ,  $

    index         : index         ,  $
    cell          : cell          ,  $
    path          : path          ,  $
    delta         : delta         ,  $

    dummy         : -1            }
  RETURN,result
END;
; ======================================================================
;
FUNCTION cortex_LoadEireneSpectrum, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  bin           = inGetData('BIN')
  flux          = inGetData('FLUX')
  stdev         = inGetData('STDE') 
  integral      = inGetData('INTEGRAL') 
  species_type  = inGetData('SPECIES_TYPE')
  spectrum_type = inGetData('SPECTRUM_TYPE')
  index         = inGetData('INDEX')
  min_value     = inGetData('MIN_VALUE')
  max_value     = inGetData('MAX_VALUE')
  stratum       = inGetData('STRATUM')
  inCloseInterface  

  result = {  $
    verison       : 1.0           ,  $
    file          : file          ,  $
    bin           : bin           ,  $
    flux          : flux          ,  $
    stdev         : stdev         ,  $
    integral      : integral      ,  $
    species_type  : species_type  ,  $
    spectrum_type : spectrum_type ,  $
    index         : index         ,  $
    min_value     : min_value     ,  $
    max_value     : max_value     ,  $
    stratum       : stratum       }
  RETURN,result
END
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

  z        = inGetData('Z')
  a        = inGetData('A')
  charge   = inGetData('CHARGE')
  wlngth   = inGetData('WLNGTH')
  integral = inGetData('INT')
  v1       = inGetData('V1')
  v2       = inGetData('V2')

  path   = inGetData('PATH')
  delta  = inGetData('DELTA')
  dens   = inGetData('NE')
  te     = inGetData('TE')
  ti     = inGetData('TI')
  n_d    = inGetData('N_D')
  n_d2   = inGetData('N_D2')
  signal = inGetData('SIGNAL_01')

  inCloseInterface  

  result = {  $
    version  : 1.0      ,  $
    file     : file     ,  $
    z        : z        ,  $
    a        : a        ,  $
    charge   : charge   ,  $
    wlngth   : wlngth   ,  $
    integral : integral ,  $
    v1       : v1       ,  $
    v2       : v2       ,  $
    path     : path     ,  $
    delta    : delta    ,  $
    dens     : dens     ,  $
    te       : te       ,  $
    ti       : ti       ,  $
    n_d      : n_d      ,  $
    n_d2     : n_d2     ,  $
    signal   : signal   }

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
