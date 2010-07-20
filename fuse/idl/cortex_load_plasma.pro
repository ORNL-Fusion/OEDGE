;
; ======================================================================
;
FUNCTION cortex_LoadWallProfiles, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  index           = inGetData('INDEX')
  index_pin       = inGetData('INDEX_PIN')
  r_cen           = inGetData('R_CEN')
  z_cen           = inGetData('Z_CEN')
  length          = inGetData('LENGTH')
  atom_par_flux   = inGetData('ATOM_PAR_FLUX')
  atom_avg_energy = inGetData('ATOM_AVG_ENERGY')
  inCloseInterface  

  atom_energy_flux = atom_par_flux * atom_avg_energy * 1.603E-19 * 1.0E-6
  toroidal_area    = 2.0 * 3.141592 * r_cen * length

  result = {  $
    version : 1.0  ,  $
    file    : file ,  $
    index            : index            ,  $
    index_pin        : index_pin        ,  $
    r_cen            : r_cen            ,  $
    z_cen            : z_cen            ,  $
    length           : length           ,  $
    toroidal_area    : toroidal_area    ,  $
    atom_par_flux    : atom_par_flux    ,  $
    atom_avg_energy  : atom_avg_energy  ,  $
    atom_energy_flux : atom_energy_flux }

  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadDIVIMPSummary, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  ions_requested        = inGetData('IONS_REQUESTED') 
  neutrals_launched     = inGetData('NEUTRALS_LAUNCHED') 
  neutrals_failed       = inGetData('NEUTRALS_FAILED') 
  ions_created          = inGetData('IONS_CREATED') 
  ions_reaching_core    = inGetData('IONS_REACHING_CORE') 
  ions_lost_wall        = inGetData('IONS_LOST_WALL')    
  ions_lost_targets     = inGetData('IONS_LOST_TARGET')
  ions_lost_state_limit = inGetData('IONS_LOST_STATE_LIMIT')
  ions_lost_recombined  = inGetData('IONS_LOST_RECOMBINED')
  ion_atomic_number     = inGetData('ION_ATOMIC_NUMBER')    
  max_charge_state      = inGetData('MAX_CHARGE_STATE')
  inCloseInterface  

  result = {  $
    verison               : 1.0   ,  $
    file                  : file  ,  $
    ions_requested        : ions_requested        ,  $
    neutrals_launched     : neutrals_launched     ,  $
    neutrals_failed       : neutrals_failed       ,  $
    ions_created          : ions_created          ,  $
    ions_reaching_core    : ions_reaching_core    ,  $
    ions_lost_wall        : ions_lost_wall        ,  $
    ions_lost_targets     : ions_lost_targets     ,  $
    ions_lost_state_limit : ions_lost_state_limit ,  $
    ions_lost_recombined  : ions_lost_recombined  ,  $
    ion_atomic_number     : ion_atomic_number     ,  $
    max_charge_state      : max_charge_state      }

  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadNodeData, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  tube  = inGetData('TUBE') 
  s_opt = inGetData('S_OPT') 
  i_sym = inGetData('M_NODE') 
  i_end = inGetData('N_NODE') 

  n = N_ELEMENTS(tube)
  m = MAX(i_end)

  s    = MAKE_ARRAY(m,n,/FLOAT,VALUE=0.0)
  jsat = MAKE_ARRAY(m,n,/FLOAT,VALUE=0.0)
  dens = MAKE_ARRAY(m,n,/FLOAT,VALUE=0.0)
  pe   = MAKE_ARRAY(m,n,/FLOAT,VALUE=0.0)
  te   = MAKE_ARRAY(m,n,/FLOAT,VALUE=0.0)
  ti   = MAKE_ARRAY(m,n,/FLOAT,VALUE=0.0)
  FOR i = 1, m DO BEGIN
    tag = 'NODE_'+STRING(i,FORMAT='(I1)')+'_'  ; *** BREAKS IF MORE THAN 9 NODES! LOOK UP ZERO PADDING FORMAT SPECIFIER! ***
    s   [i-1,*] = inGetData(tag+'S')
    jsat[i-1,*] = inGetData(tag+'JSAT')
    dens[i-1,*] = inGetData(tag+'DENS')
    pe  [i-1,*] = inGetData(tag+'PE')
    te  [i-1,*] = inGetData(tag+'TE')
    ti  [i-1,*] = inGetData(tag+'TI')
  ENDFOR

  type        = inGetData('TYPE')
  tube_range1 = inGetData('TUBE_RANGE1')
  tube_range2 = inGetData('TUBE_RANGE2')
  x           = inGetData('RAD_X')
  y           = inGetData('RAD_Y')

  inCloseInterface  

  result = {  $
    verison : 1.0   ,  $
    file    : file  ,  $
    tube    : tube  ,  $
    s_opt   : s_opt ,  $
    i_sym   : i_sym ,  $
    i_end   : i_end ,  $
    s       : s     ,  $
    jsat    : jsat  ,  $
    dens    : dens  ,  $
    pe      : pe    ,  $
    te      : te    ,  $
    ti      : ti    ,  $
    type    : type  ,  $
    tube_range1 : tube_range1 ,  $
    tube_range2 : tube_range2 ,  $
    x       : x     ,  $
    y       : y     } 

  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadEireneData, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  index    = inGetData('INDEX') 
  pos      = inGetData('POS')   
  tube     = inGetData('TUBE') 
  s        = inGetData('S') 
  ion_net  = inGetData('ION_NET')
  rec_net  = inGetData('REC_NET')
  mom_net  = inGetData('MOM_NET')
  qe_net   = inGetData('QE_NET')
  qi_net   = inGetData('QI_NET')
  atm_dens = inGetData('ATM_DENS')
  mol_dens = inGetData('MOL_DENS')
  balmer_alpha = inGetData('BALMER_ALPHA')
  balmer_gamma = inGetData('BALMER_GAMMA')
  inCloseInterface  

  result = {  $
    verison      : 1.0          ,  $
    file         : file         ,  $
;    index        : index        ,  $
    pos          : pos          ,  $  
    tube         : tube         ,  $
    s            : s            ,  $
    ion_net      : ion_net      ,  $
    rec_net      : rec_net      ,  $
    mom_net      : mom_net      ,  $
    qe_net       : qe_net       ,  $
    qi_net       : qi_net       ,  $
    atm_dens     : atm_dens     ,  $
    mol_dens     : mol_dens     ,  $
    balmer_alpha : balmer_alpha ,  $
    balmer_gamma : balmer_gamma }
  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadEnergySpectrum, file

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
  surface_index = inGetData('SURFACE_INDEX')
  min_energy    = inGetData('MIN_ENERGY')
  max_energy    = inGetData('MAX_ENERGY')
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
    surface_index : surface_index ,  $
    min_energy    : min_energy    ,  $
    max_energy    : max_energy    ,  $
    stratum       : stratum       }
  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadSourceData, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  index   = inGetData('INDEX') 
  tube    = inGetData('TUBE' ) 
  s       = inGetData('S'    ) 
  par_net = inGetData('PAR_NET')
  par_ion = inGetData('PAR_ION')
  par_rec = inGetData('PAR_REC')
  par_usr = inGetData('PAR_USR')
  par_ano = inGetData('PAR_ANO')
  mom_net = inGetData('MOM_NET')
  mom_vol = inGetData('MOM_VOL')
  mom_usr = inGetData('MOM_USR')
  mom_ano = inGetData('MOM_ANO')
  ene_net = inGetData('ENE_NET')
  ene_ion = inGetData('ENE_ION')
  ene_fit = inGetData('ENE_FIT')
  eni_net = inGetData('ENI_NET')
  inCloseInterface  

  result = {             $
    verison : 1.0     ,  $
    file    : file    ,  $
    index   : index   ,  $
    tube    : tube    ,  $
    s       : s       ,  $
    par_net : par_net ,  $
    par_ion : par_ion ,  $ 
    par_rec : par_rec ,  $
    par_usr : par_usr ,  $
    par_ano : par_ano ,  $
    mom_net : mom_net ,  $
    mom_vol : mom_vol ,  $
    mom_usr : mom_usr ,  $
    mom_ano : mom_ano ,  $
    ene_net : ene_net ,  $
    ene_ion : ene_ion ,  $
    ene_fit : ene_fit ,  $
    eni_net : eni_net }
  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadPlasmaData, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  index = inGetData('INDEX') 
  tube  = inGetData('TUBE') 
  s     = inGetData('S') 
  dens  = inGetData('NE')
  vi    = inGetData('VI')
  cs    = inGetData('CS')
  M     = inGetData('MACHNO')
  pe    = inGetData('PE')
  pi    = inGetData('PI')
  te    = inGetData('TE')
  ti    = inGetData('TI')
  inCloseInterface  

  result = {           $
    verison : 1.0   ,  $
    file    : file  ,  $
    index   : index ,  $
    tube    : tube  ,  $
    s       : s     ,  $
    dens    : dens  ,  $
    vi      : vi    ,  $
    cs      : cs    ,  $
    M       : M     ,  $
    pe      : pe    ,  $
    pi      : pi    ,  $
    te      : te    ,  $
    ti      : ti    }
  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadPedestalModel, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  a        = inGetData('PED_A'       )
  volume   = inGetData('CORE_VOLUME' )
  cross_ne = inGetData('PED_CROSS_NE') 
  cross_te = inGetData('PED_CROSS_TE')
  cross_ti = inGetData('PED_CROSS_TI')
  r        = inGetData('PED_R'       )
  volfr    = inGetData('PED_VOLFR'   )
  dens     = inGetData('PED_NE'      )
  te       = inGetData('PED_TE'      )
  ti       = inGetData('PED_TI'      )
  inCloseInterface  

  new_r = (FINDGEN(100.0) + 0.5) / 100.0 * a
  new_deltar = r[1] - r[0]

;print,new_r / a

  new_dens = INTERPOL(dens,r,new_r)
  new_te   = INTERPOL(te  ,r,new_r)
  new_ti   = INTERPOL(ti  ,r,new_r)

  new_volume = new_r * new_deltar
  new_volume = new_volume / TOTAL(new_volume) * volume

  PRINT,new_volume
  print,TOTAL(new_volume)

  volume = volfr * volume

  result = {                $
    verison    : 1.0        ,  $
    file       : file       ,  $
    a          : a          ,  $
    cross_ne   : cross_ne   ,  $
    cross_te   : cross_te   ,  $
    cross_ti   : cross_ti   ,  $
    new_r      : new_r      ,  $
    new_volume : new_volume ,  $
    new_dens   : new_dens   ,  $
    new_te     : new_te     ,  $
    new_ti     : new_ti     ,  $
    r          : r          ,  $
    volume     : volume     ,  $
    dens       : dens       ,  $
    te         : te         ,  $
    ti         : ti    }
  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadTargetData, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  tube        = inGetData('TAR_TUBE')
  ring        = inGetData('TAR_RING')
  rho         = inGetData('TAR_RHO')  
  psin        = inGetData('TAR_PSIN')  
  lo_location = inGetData('TAR_LO_LOCATION')
  lo_s    = inGetData('TAR_LO_S')
  lo_p    = inGetData('TAR_LO_P')
  lo_jsat = inGetData('TAR_LO_JSAT')
  lo_dens = inGetData('TAR_LO_NE')
  lo_M    = inGetData('TAR_LO_MACHNO')
  lo_pe   = inGetData('TAR_LO_PE')
  lo_pi   = inGetData('TAR_LO_PI')
  lo_te   = inGetData('TAR_LO_TE')
  lo_ti   = inGetData('TAR_LO_TI')
  hi_location = inGetData('TAR_HI_LOCATION')
  hi_s    = inGetData('TAR_HI_S')
  hi_p    = inGetData('TAR_HI_P')
  hi_jsat = inGetData('TAR_HI_JSAT')
  hi_dens = inGetData('TAR_HI_NE')
  hi_M    = inGetData('TAR_HI_MACHNO')
  hi_pe   = inGetData('TAR_HI_PE')
  hi_pi   = inGetData('TAR_HI_PI')
  hi_te   = inGetData('TAR_HI_TE')
  hi_ti   = inGetData('TAR_HI_TI')
  inCloseInterface

  location = MAKE_ARRAY(N_ELEMENTS(lo_location),2,/FLOAT,VALUE=0.0)
  s    = MAKE_ARRAY(N_ELEMENTS(lo_s   ),2,/FLOAT,VALUE=0.0)
  p    = MAKE_ARRAY(N_ELEMENTS(lo_p   ),2,/FLOAT,VALUE=0.0)
  jsat = MAKE_ARRAY(N_ELEMENTS(lo_jsat),2,/FLOAT,VALUE=0.0)
  dens = MAKE_ARRAY(N_ELEMENTS(lo_dens),2,/FLOAT,VALUE=0.0)
  M    = MAKE_ARRAY(N_ELEMENTS(lo_M   ),2,/FLOAT,VALUE=0.0)
  pe   = MAKE_ARRAY(N_ELEMENTS(lo_pe  ),2,/FLOAT,VALUE=0.0)
  pi   = MAKE_ARRAY(N_ELEMENTS(lo_pi  ),2,/FLOAT,VALUE=0.0)
  te   = MAKE_ARRAY(N_ELEMENTS(lo_te  ),2,/FLOAT,VALUE=0.0)
  ti   = MAKE_ARRAY(N_ELEMENTS(lo_ti  ),2,/FLOAT,VALUE=0.0)

  location[*,0] = lo_location[*]
  location[*,1] = hi_location[*]

  s   [*,0] = lo_s   [*]
  s   [*,1] = hi_s   [*]
  p   [*,0] = lo_p   [*]
  p   [*,1] = hi_p   [*]
  jsat[*,0] = ABS(lo_jsat[*])  ; *** Shorter way to do this remapping..? ***
  jsat[*,1] = ABS(hi_jsat[*])
  dens[*,0] = lo_dens[*]
  dens[*,1] = hi_dens[*]
  M   [*,0] = lo_M   [*]
  M   [*,1] = hi_M   [*]
  pe  [*,0] = lo_pe  [*]
  pe  [*,1] = hi_pe  [*]
  pi  [*,0] = lo_pi  [*]
  pi  [*,1] = hi_pi  [*]
  te  [*,0] = lo_te  [*]
  te  [*,1] = hi_te  [*]
  ti  [*,0] = lo_ti  [*]
  ti  [*,1] = hi_ti  [*]

  result = {          $
    verison : 1.0  ,  $
    file    : file ,  $
    tube    : tube ,  $
    ring    : ring ,  $
    rho     : rho  ,  $
    psin    : psin ,  $
    location : location,  $
    s       : s    ,  $
    p       : p    ,  $
    jsat    : jsat ,  $
    dens    : dens ,  $
    M       : M    ,  $
    pe      : pe   ,  $
    pi      : pi   ,  $
    te      : te   ,  $
    ti      : ti   }
  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadMidplaneProfiles, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  ring   = inGetData('MID_RING')  
  r      = inGetData('MID_R')  
  rho    = inGetData('MID_RHO')  
  psin   = inGetData('MID_PSIN')  
  l      = inGetData('MID_L')
  dens   = inGetData('MID_NE')
  vb     = inGetData('MID_VB')
  machno = inGetData('MID_M' )
  p      = inGetData('MID_P' )
  te     = inGetData('MID_TE')
  ti     = inGetData('MID_TI')
  inCloseInterface

  result = {  $
    version : 1.0      ,  $
    file    : file     ,  $
    ring    : ring     ,  $
    r       : r        ,  $
    rho     : rho      ,  $
    psin    : psin     ,  $
    l       : l        ,  $
    dens    : dens     ,  $
    vb      : vb       ,  $
    machno  : machno   ,  $
    p       : p        ,  $
    te      : te       ,  $
    ti      : ti  }

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_LoadEIRENEImpurityData, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  grid_isep = inGetData('GRID_ISEP')
  grid_ipfz = inGetData('GRID_IPFZ')
  pos  = inGetData('POS')  
  tube = inGetData('TUBE')  
  s    = inGetData('S')
  imp_dens  = MAKE_ARRAY(N_ELEMENTS(tube),1,/FLOAT,VALUE=0.0)
  imp_ioniz = imp_dens
  imp_dens [*,0] = inGetData('IMP_DENS_00') 
  imp_ioniz[*,0] = inGetData('IMP_IONIZ_00') 
  inCloseInterface

  result = {  $
    version   : 1.0       ,  $
    file      : file      ,  $
    grid_isep : grid_isep ,  $
    grid_ipfz : grid_ipfz ,  $
    pos       : pos       ,  $
    tube      : tube      ,  $
    s         : s         ,  $
    imp_dens  : imp_dens  ,  $
    imp_ioniz : imp_ioniz }

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_LoadDIVIMPImpurityData_Density, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  div_influx  = inGetData('DIV_IMPURITY_INFLUX')  
  eir_influx  = inGetData('EIR_IMPURITY_INFLUX')  
  div_init_iz = inGetData('IMP_INITIAL_IZ')  
  div_max_iz  = inGetData('IMP_MAX_IZ')
  z           = inGetData('IMP_Z')
  a           = inGetData('IMP_A')

  grid_isep = inGetData('GRID_ISEP')
  grid_ipfz = inGetData('GRID_IPFZ')

  pos  = inGetData('POS')  
  tube = inGetData('TUBE')  
  s    = inGetData('S')

  imp_data = MAKE_ARRAY(N_ELEMENTS(tube),z+2,/FLOAT,VALUE=0.0)
  FOR i = 0, z DO imp_data[*,i] = inGetData('IMP_DENS_'+STRING(i,FORMAT='(I02)'))
  FOR j = 0, N_ELEMENTS(ring)-1 DO imp_data[j,z+1] = TOTAL(imp_data[j,1:z])

  inCloseInterface

  result = {  $
    version     : 1.0         ,  $
    file        : file        ,  $
    grid_isep   : grid_isep   ,  $
    grid_ipfz   : grid_ipfz   ,  $
    div_influx  : div_influx  ,  $
    eir_influx  : eir_influx  ,  $
    div_init_iz : div_init_iz ,  $
    div_max_iz  : div_max_iz  ,  $
    z           : z           ,  $
    a           : a           ,  $
    pos         : pos         ,  $
    tube        : tube        ,  $
    s           : s           ,  $
    imp_dens    : imp_data    }

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_LoadDIVIMPImpurityData_Ionisation, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  div_influx  = inGetData('DIV_IMPURITY_INFLUX')  
  eir_influx  = inGetData('EIR_IMPURITY_INFLUX')  
  div_init_iz = inGetData('IMP_INITIAL_IZ')  
  div_max_iz  = inGetData('IMP_MAX_IZ')
  z           = inGetData('IMP_Z')
  a           = inGetData('IMP_A')

  grid_isep = inGetData('GRID_ISEP')
  grid_ipfz = inGetData('GRID_IPFZ')

  pos  = inGetData('POS')  
  tube = inGetData('TUBE')  
  s    = inGetData('S')

  imp_data = MAKE_ARRAY(N_ELEMENTS(tube),z+2,/FLOAT,VALUE=0.0)
  FOR i = 0, z DO imp_data[*,i] = inGetData('IMP_IONIZ_'+STRING(i,FORMAT='(I02)'))
  FOR j = 0, N_ELEMENTS(ring)-1 DO imp_data[j,z+1] = TOTAL(imp_data[j,1:z])

  inCloseInterface

  result = {  $
    version     : 1.0         ,  $
    file        : file        ,  $
    grid_isep   : grid_isep   ,  $
    grid_ipfz   : grid_ipfz   ,  $
    div_influx  : div_influx  ,  $
    eir_influx  : eir_influx  ,  $
    div_init_iz : div_init_iz ,  $
    div_max_iz  : div_max_iz  ,  $
    z           : z           ,  $
    a           : a           ,  $
    pos         : pos         ,  $
    tube        : tube        ,  $
    s           : s           ,  $
    imp_ioniz   : imp_data    }

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_LoadCoreProfiles, file    

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  div_influx  = inGetData('DIV_IMPURITY_INFLUX')  
  eir_influx  = inGetData('EIR_IMPURITY_INFLUX')  
  div_init_iz = inGetData('IMP_INITIAL_IZ')  
  div_max_iz  = inGetData('IMP_MAX_IZ')
  z           = inGetData('IMP_Z')
  a           = inGetData('IMP_A')

  ring     = inGetData('RING')  
  r        = inGetData('MID_R')  
  r_a      = inGetData('MID_R/A')  
  rho      = inGetData('RHO')  
  psin     = inGetData('PSIN')  
  l        = inGetData('L')
  dens     = inGetData('AVG_NE')
  p        = dens
  te       = inGetData('AVG_TE')
  ti       = inGetData('AVG_TI')
  volume   = inGetData('IMP_VOL')
  i_frac   = inGetData('IMP_FRAC_%')
  i_frac_e = inGetData('IMP_FRAC_E_%')
  zeff     = inGetData('IMP_ZEFF')

  ni = MAKE_ARRAY(z+2,N_ELEMENTS(ring),/FLOAT,VALUE=0.0)
  FOR i = 1, z DO ni[i,*] = inGetData('IMP_AVG_DENS_'+STRING(i,FORMAT='(I02)'))
  FOR i = 0, N_ELEMENTS(ring)-1 DO ni[z+1,i] = TOTAL(ni[1:z,i])

  inCloseInterface

  result = {  $
    version     : 1.0         ,  $
    file        : file        ,  $
    div_influx  : div_influx  ,  $
    eir_influx  : eir_influx  ,  $
    div_init_iz : div_init_iz ,  $
    div_max_iz  : div_max_iz  ,  $
    z           : z           ,  $
    a           : a           ,  $
    ring        : ring        ,  $
    r           : r           ,  $
    r_a         : r_a         ,  $
    rho         : rho         ,  $
    psin        : psin        ,  $
    l           : l           ,  $
    dens        : dens        ,  $
    p           : p           ,  $
    te          : te          ,  $
    ti          : ti          ,  $
    volume      : volume      ,  $
    i_frac      : i_frac      ,  $
    i_frac_e    : i_frac_e    ,  $
    zeff        : zeff        ,  $
    ni          : ni          }

  RETURN, result
END
;
; ======================================================================
;