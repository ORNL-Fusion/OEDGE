;
; ======================================================================
;
FUNCTION cortex_LoadEireneHistoryData, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  strings   = inGetData('STRINGS')
  time_end  = inGetData('TIME')
  gauge     = inGetData('GAUGE_I1')
  number    = inGetData('GAUGE_I2')
  phi       = inGetData('GAUGE_PHI') * (-180.0 / 3.141592) 
  strata    = inGetData('STRATA')
  iteration = inGetData('FLUID_ITERATION') 
  history   = inGetData('HISTORY') 
  index     = inGetData('GAUGE') 
  x         = inGetData('GAUGE_X')
  y         = inGetData('GAUGE_Y') 
  volume    = inGetData('VOLUME') 
  dens_atm  = inGetData('DENS_ATM') 
  dens_mol  = inGetData('DENS_MOL') 
  Eavg_atm  = inGetData('EAVG_ATM') 
  Eavg_mol  = inGetData('EAVG_MOL') 
  p2_atm    = inGetData('P2_ATM') 
  p2_mol    = inGetData('P2_MOL') 

  result = {  $
    version   : 1.0       ,  $
    file      : file      ,  $
    strings   : strings   ,  $
    time_end  : time_end  ,  $
    gauge     : gauge     ,  $
    number    : number    ,  $
    phi       : phi       ,  $
    strata    : strata    ,  $
    iteration : iteration ,  $
    history   : history   ,  $
    index     : index     ,  $
    dens_atm  : dens_atm  ,  $
    dens_mol  : dens_mol  ,  $
    Eavg_atm  : Eavg_atm  ,  $
    Eavg_mol  : Eavg_mol  ,  $
    p2_atm    : p2_atm    ,  $
    p2_mol    : p2_mol    ,  $
    dummy     : 1.0       }

  RETURN,result

END
;
; ======================================================================
;
FUNCTION cortex_LoadSurfaceTemperatures, file

;  file = cortex_UpdateFile(file)

;  status = inOpenInterface(file)
;  IF (status LT 0) THEN BEGIN
;    result = CREATE_STRUCT('version',0.0,'file','none')
;    RETURN, result
;  ENDIF

  lun = 3
  FREE_LUN, lun
  OPENR, lun, file, error=err
  IF (err NE 0) THEN BEGIN
    PRINT,'ERROR cortex_LoadSurfaceTemperatures: Unable to open file'
    PRINT,'  FILE = ',file 
    RETURN, -1
  ENDIF
  first_pass = 1
  buffer = ' '
  WHILE (NOT EOF(lun)) DO BEGIN
    READF,lun,buffer
    i = STRPOS(buffer,'*')
    IF (i EQ 0) THEN CONTINUE
    str = STRSPLIT(STRTRIM(buffer,2),' ',/EXTRACT)
    IF (first_pass) THEN BEGIN
      dist = FLOAT(str[0])
      temp = FLOAT(str[1])
      first_pass = 0
    ENDIF ELSE BEGIN
      dist = [dist,FLOAT(str[0])]
      temp = [temp,FLOAT(str[1])]
    ENDELSE
  ENDWHILE
  CLOSE,lun

  result = {  $
    version : 1.0  ,  $
    file    : file ,  $
    dist    : dist ,  $
    temp    : temp ,  $
    dummy   : -1   }

  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadSpectraInfo, file

  file = cortex_UpdateFile(file)

print,'spectra info loader',file

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  index = inGetData('INDEX')
  cell  = inGetData('CELL') 
  path  = inGetData('PATH') 
  delta = inGetData('DELTA') 
  v_ind = inGetData('V_IND') 
  v_1x  = inGetData('V_1X') 
  v_1y  = inGetData('V_1Y') 
  v_2x  = inGetData('V_2X') 
  v_2y  = inGetData('V_2Y') 

  result = {  $
    version : 1.0   ,  $
    file    : file  ,  $
    index   : index ,  $
    cell    : cell  ,  $
    path    : path  ,  $
    delta   : delta ,  $
    v_ind   : v_ind ,  $
    v_1x    : v_1x  ,  $
    v_1y    : v_1y  ,  $
    v_2x    : v_2x  ,  $
    v_2y    : v_2y  ,  $
    dummy   : -1    }

  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadEireneStrata, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  type   = inGetData('TYPE')
  target = inGetData('TARGET')
  range1 = inGetData('RANGE1')
  range2 = inGetData('RANGE2') 

  result = {  $
    verison : 1.0    ,  $
    file    : file   ,  $
    type    : type   ,  $
    target  : target ,  $
    range1  : range1 ,  $
    range2  : range2 }

  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadEireneSpectrum, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file,/no_check)
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

print,'profile file loader',file

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
FUNCTION cortex_LoadIntegrals, file, direct=direct

  IF (N_ELEMENTS(direct) EQ 0) THEN file = cortex_UpdateFile(file)

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
FUNCTION cortex_LoadImpurityEmission_DIVIMP, file, direct=direct

  IF (N_ELEMENTS(direct) EQ 0) THEN file = cortex_UpdateFile(file)

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
  pos      = inGetData('POS')
  tube     = inGetData('TUBE')

  signal = MAKE_ARRAY(N_ELEMENTS(pos),n,/DOUBLE,VALUE=0.0D)
  FOR k = 1, n DO BEGIN
    tag = 'SIGNAL_'+STRING(k,FORMAT='(I02)')
    signal[*,k-1] = inGetData(tag)
  ENDFOR

  inCloseInterface  

  result = {  $
    version       : 1.0      ,  $
    file          : file     ,  $
    n_signal      : n        ,  $
    pos           : pos      ,  $
    tube          : tube     ,  $
    atomic_number : z        ,  $
    atomic_mass   : a        ,  $
    charge        : charge   ,  $ 
    wavelength    : wlngth   ,  $ 
    database      : database ,  $
    signal        : signal   }

  RETURN,result
END
;
; ======================================================================
;
