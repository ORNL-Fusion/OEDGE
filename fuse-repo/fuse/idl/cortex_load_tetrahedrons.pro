;
; ======================================================================
;
FUNCTION cortex_LoadTetrahedronSlice, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  r1  = inGetData('R1')
  z1  = inGetData('Z1')
  r2  = inGetData('R2')
  z2  = inGetData('Z2')
  r3  = inGetData('R3')
  z3  = inGetData('Z3')
  r4  = inGetData('R4')
  z4  = inGetData('Z4')
  nv  = inGetData('NV')
  data = inGetData('VAL')
  opt = inGetData('OPT')

  inCloseInterface

  x = [ [r1], [r2], [r3], [r4] ]
  y = [ [z1], [z2], [z3], [z4] ]

  result = {  $
    file    : file ,  $
    version : 1.0  ,  $
    opt     : opt  ,  $
    x       : x    ,  $
    y       : y    ,  $
    nv      : nv   ,  $
    data    : data   }

  RETURN, result
END
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

  iobj     = inGetData('IOBJ')
  region   = inGetData('REGION')
  x        = inGetData('X')  
  y        = inGetData('Y')  
  z        = inGetData('Z')  
  dens     = inGetData('NE')  
;  te       = inGetData('TE')
;  dens_D   = inGetData('D_DENS')
;  dens_D2  = inGetData('D2_DENS')
;  avgeng_D = inGetData('D_AVGENG')
;  dalpha   = inGetData('D_ALPHA')
  s_ion    = inGetData('S_ION')

  inCloseInterface

  i = STRPOS(file,'/',/REVERSE_SEARCH)
  n = N_ELEMENTS(x)

  result = {  $
    desc      : 'tetrahedron data' ,  $
    file      : STRMID(file,i+1)   ,  $
    version   : 1.0      ,  $
    n         : n        ,  $
    iobj      : iobj     ,  $
    region    : region   ,  $
    x         : x        ,  $
    y         : y        ,  $
    z         : z        ,  $
    dens      : dens     ,  $
;    te        : te       ,  $
;    dens_D    : dens_D   ,  $
;    dens_D2   : dens_D2  ,  $
;    avgeng_D  : avgeng_D ,  $
;    dalpha    : dalpha   ,  $
    s_ion     : s_ion    }

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_LoadTriangleGeometry, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  x  = inGetData('X')  
  y  = inGetData('Y')  
  z  = inGetData('Z')
  v0 = inGetData('V1')
  v1 = inGetData('V2')
  v2 = inGetData('V3')

  i = inGetData('ISRF')

  n = N_ELEMENTS(v0)
  v = MAKE_ARRAY(3,n,/LONG,VALUE=0)
  v[0,*] = v0 - 1
  v[1,*] = v1 - 1
  v[2,*] = v2 - 1

  inCloseInterface

  result = {  $
    file : file   ,  $   
    type : 'wall' ,  $
    n    : n      ,  $
    i    : i      ,  $
    v    : v      ,  $
    x    : x      ,  $
    y    : y      ,  $
    z    : z      }

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_LoadTriangleData, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  x  = inGetData('X')  
  y  = inGetData('Y')  
  z  = inGetData('Z')
  v0 = inGetData('V1')
  v1 = inGetData('V2')
  v2 = inGetData('V3')

  i              = inGetData('ISRF')
  in_std         = inGetData('IN_STD')
  in_add         = inGetData('IN_ADD')
  iliin          = inGetData('ILIIN')
  area           = inGetData('AREA')
  nt_par_atm_1_0 = inGetData('NT_PAR_ATM_1_0') * 1.602E-19
  in_par_atm_1_0 = inGetData('IN_PAR_ATM_1_0') * 1.602E-19
  in_ene_atm_1_0 = inGetData('IN_ENE_ATM_1_0') * 1.602E-19
  nt_par_mol_1_0 = inGetData('NT_PAR_MOL_1_0') * 1.602E-19
  in_par_mol_1_0 = inGetData('IN_PAR_MOL_1_0') * 1.602E-19
  in_ene_mol_1_0 = inGetData('IN_ENE_MOL_1_0') * 1.602E-19
;  nt_par_imp_1_0 = inGetData('NT_PAR_IMP_1_0') * 1.602E-19

  n = N_ELEMENTS(in_std)
  v = MAKE_ARRAY(3,n,/LONG,VALUE=0)
  v[0,*] = v0 - 1
  v[1,*] = v1 - 1
  v[2,*] = v2 - 1

  inCloseInterface

  result = {  $
    file           : file           ,  $   
    type           : 'flux'         ,  $
    n              : n              ,  $
    i              : i              ,  $
    in_std         : in_std         ,  $
    in_add         : in_add         ,  $
    iliin          : iliin          ,  $
    nt_par_atm_1_0 : nt_par_atm_1_0 ,  $
    in_par_atm_1_0 : in_par_atm_1_0 ,  $
    in_ene_atm_1_0 : in_ene_atm_1_0 ,  $
    nt_par_mol_1_0 : nt_par_mol_1_0 ,  $
    in_par_mol_1_0 : in_par_mol_1_0 ,  $
    in_ene_mol_1_0 : in_ene_mol_1_0 ,  $
;    nt_par_imp_1_0 : nt_par_imp_1_0 ,  $
    area           : area           ,  $
    v              : v              ,  $
    x              : x              ,  $
    y              : y              ,  $
    z              : z              }

  RETURN, result
END
;
; ======================================================================
;
