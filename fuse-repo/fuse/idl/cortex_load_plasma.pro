;
; ======================================================================
;
FUNCTION cortex_SolverDebugData, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

;  index = inGetData('INDEX') 
;  tube  = inGetData('TUBE') 
  s     = inGetData('ZMESH') 
;  r     = inGetData('R') 
;  z     = inGetData('Z') 
  dens  = inGetData('NE')
  vi    = inGetData('VB')
;  cs    = inGetData('CS')
;  M     = inGetData('MACHNO')
;  pe    = inGetData('PE')
;  pi    = inGetData('PI')
  te    = inGetData('TE')
  ti    = inGetData('TI')
  inCloseInterface  

  result = {           $
    verison : 1.0   ,  $
    file    : file  ,  $
;    index   : index ,  $
;    tube    : tube  ,  $
    s       : s     ,  $
;    r       : r     ,  $
;    z       : z     ,  $
    dens    : dens  ,  $
    vi      : vi    ,  $
;    cs      : cs    ,  $
;    M       : M     ,  $
;    pe      : pe    ,  $
;    pi      : pi    ,  $
    te      : te    ,  $
    ti      : ti    }

  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadSputterData_DIVIMP, file

  print, 'Loading DIVIMP sputter data...'

  file = cortex_UpdateFile(file)



  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  absfac_total  = inGetData('ABSFAC_TOTAL')   
  absfac        = inGetData('ABSFAC')   
  absfac_in     = inGetData('ABSFAC_IN')     ; scaling factor for the incident particle flux
  bomb_z        = inGetData('BOMB_Z')        ; bombarding species atomic number
  bomb_a        = inGetData('BOMB_A')        ; bombarding species atomic mass
  target_z      = inGetData('TARGET_Z')      ; target material atomic number   

  modifier      = inGetData('FACT')          ; wall segment surface sputtering modifier in DIVIMP
  launch_prob   = inGetData('PROB')          ; wall segment impurity launch probability

  index         = inGetData('ID')            ; wall segment index in DIVIMP
  length        = inGetData('LENGTH')        ; wall segment length (m) in poloidal plane
  pos_r         = inGetData('POS_R')         ; wall segment R position in machine coordinates (poloidal plane)
  pos_z         = inGetData('POS_Z')         ; wall segment Z position in machine coordinates
  modifier      = inGetData('FACT')          ; yield modifier set in DIVIMP

  n = N_ELEMENTS(index)

  distance = MAKE_ARRAY(n,/FLOAT,VALUE=0.0)
  distance[0] = 0.5 * length[0]
  FOR i = 1, n-1 DO BEGIN
    distance[i] = distance[i-1] + 0.5 * length[i-1]
    distance[i] = distance[i  ] + 0.5 * length[i  ]
  ENDFOR

  result = {  $
    version   : 1.0      ,  $
    file      : file     ,  $
    absfac    : absfac   ,  $
    absfac_in : absfac_in,  $
    bomb_z    : bomb_z   ,  $
    bomb_a    : bomb_a   ,  $
    target_z  : target_z ,  $
    index     : index    ,  $
    pos_r     : pos_r    ,  $
    pos_z     : pos_z    ,  $
    length    : length   ,  $
    modifier  : modifier ,  $
    distance  : distance }
 
  FOR i = 1, N_ELEMENTS(absfac) DO BEGIN

    t1 = '_' + STRING(i,FORMAT='(I02)')

;    print,'t1=',t1
    max_charge = inGetData('MAX_CHARGE'+t1)
;    print,'max charge',max_charge

    data_set = { max_charge : max_charge }

    flux  = MAKE_ARRAY(n,max_charge,/FLOAT,VALUE=0.0)
    e0    = flux
    yield = flux

    FOR j = 1, max_charge DO BEGIN
      t2 = t1 + '_' + STRING(j,FORMAT='(I02)')
      flux [*,j-1] = inGetData('FLUX' +t2)
      e0   [*,j-1] = inGetData('E0'   +t2)
      yield[*,j-1] = inGetData('YIELD'+t2)
    ENDFOR
    data_set = CREATE_STRUCT(data_set,'flux',flux,'e0',e0,'yield',yield)  

    name = 'data' + STRTRIM(STRING(i),2)

;    print,'name ',name

    result = CREATE_STRUCT(result,name,data_set)  

  ENDFOR

  inCloseInterface  

;help,result,/struct
;help,result.data1,/struct
;print,result.data5.flux[*,1]
;psclose
;stop

  print, 'done.'

  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadWallErosion_DIVIMP, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  absfac        = inGetData('DIV_IMPURITY_INFLUX')   
  index         = inGetData('INDEX_ID')
  atomic_number = inGetData('IMP_Z')   
  atomic_mass   = inGetData('IMP_A')   
  r0            = inGetData('R0')    
  z0            = inGetData('Z0')    
  r             = inGetData('R')     
  z             = inGetData('Z')    
  dist_1        = inGetData('DIST1')    
  dist_2        = inGetData('DIST2')    
  tot_ero       = inGetData('TOT_ERO2') * absfac
  tot_dep       = inGetData('TOT_DEP2') * absfac
  tot_net       = inGetData('TOT_NET2') * absfac
  atm_ero       = inGetData('ATM_ERO')    
  atm_dep       = inGetData('ATM_DEP') 
  atm_net       = inGetData('ATM_NET')
  ion_flux      = inGetData('SURFACE_ION_FLUX')  ; perpendicular flux density
  atm_flux      = inGetData('SURFACE_ATM_FLUX')

  inCloseInterface  

  length = dist_2 - dist_1
  dist   = 0.5 * (dist_1 + dist_2)
  area   = 2.0 * !PI * r * length 

  ion_ero = tot_ero - atm_ero
  ion_dep = tot_dep - atm_dep
  ion_net = tot_net - atm_net

  result = {  $
    version       : 1.0           ,  $
    file          : file          ,  $
    index         : index         ,  $
    absfac        : absfac        ,  $
    atomic_number : atomic_number ,  $
    atomix_mass   : atomic_mass   ,  $
    r0            : r0            ,  $
    z0            : z0            ,  $
    r             : r             ,  $
    z             : z             ,  $
    dist_1        : dist_1        ,  $
    dist_2        : dist_2        ,  $
    dist          : dist          ,  $
    length        : length        ,  $
    area          : area          ,  $
    tot_ero       : tot_ero       ,  $
    tot_dep       : tot_dep       ,  $
    tot_net       : tot_net       ,  $
    ion_ero       : ion_ero       ,  $
    ion_dep       : ion_dep       ,  $
    ion_net       : ion_net       ,  $
    atm_ero       : atm_ero       ,  $
    atm_dep       : atm_dep       ,  $
    atm_net       : atm_net       ,  $
    ion_flux      : ion_flux      ,  $
    atm_flux      : atm_flux      } 

  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_LoadWallProfiles_EIRENE, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  tot_em_par_atm_2_1  = inGetData('TOT_EM_IMP_1')  ; SPECIES=2 (impurity usually), SOURCE=1 (bulk ions)
  tot_em_par_atm_2_2  = inGetData('TOT_EM_IMP_2')  ; SPECIES=2 (impurity usually), SOURCE=1 (test atoms)
  length              = inGetData('LENGTH')
  area                = inGetData('AREA')

  in_par_atm_1    = inGetData('IN_PAR_ATM_1_0')  ; SPECIES=1 (deuterium usually), SOURCE=0 (all, i.e. bulk ions, atoms, molecules, and test ions)
  in_ene_atm_1    = inGetData('IN_ENE_ATM_1_0')  

  in_par_atm_2    = inGetData('IN_PAR_ATM_2_0')  ; SPECIES=2 (impurity usually), ...see above...
  in_ene_atm_2    = inGetData('IN_ENE_ATM_2_0')  

  em_par_atm_2_1  = inGetData('EM_PAR_ATM_2_1')  ; SPECIES=1 (deuterium usually), SOURCE=1 (bulk ions)
  em_ene_atm_2_1  = inGetData('EM_ENE_ATM_2_1')  
  em_par_atm_2_2  = inGetData('EM_PAR_ATM_2_2')  ; SPECIES=2 (impurity usually), SOURCE=1 (test atoms)
  em_ene_atm_2_2  = inGetData('EM_ENE_ATM_2_2')  

  inCloseInterface  

  result = {  $
    version        : 1.0            ,  $
    file           : file           ,  $
    tot_em_par_atm_2_1 : tot_em_par_atm_2_1,  $
    tot_em_par_atm_2_2 : tot_em_par_atm_2_2,  $
    length         : length         ,  $
    area           : area           ,  $
    in_par_atm_1   : in_par_atm_1   ,  $
    in_ene_atm_1   : in_ene_atm_1   ,  $
    in_par_atm_2   : in_par_atm_2   ,  $
    in_ene_atm_2   : in_ene_atm_2   ,  $
    em_par_atm_2_1 : em_par_atm_2_1 ,  $
    em_ene_atm_2_1 : em_ene_atm_2_1 ,  $
    em_par_atm_2_2 : em_par_atm_2_2 ,  $
    em_ene_atm_2_2 : em_ene_atm_2_2 }

  RETURN,result
END
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
  index_target    = inGetData('INDEX_TARGET')
  index_cell      = inGetData('INDEX_CELL')
  index_ring      = inGetData('INDEX_RING')
  index_tube      = inGetData('INDEX_TUBE')
  r_cen           = inGetData('R_CEN')
  z_cen           = inGetData('Z_CEN')
  r_vertex1       = inGetData('R_VERTEX1')
  z_vertex1       = inGetData('Z_VERTEX1')
  r_vertex2       = inGetData('R_VERTEX2')
  z_vertex2       = inGetData('Z_VERTEX2')
  length          = inGetData('LENGTH')
  atom_par_flux   = inGetData('ATOM_PAR_FLUX')
  atom_avg_energy = inGetData('ATOM_AVG_ENERGY')
  mol_par_flux    = inGetData('MOL_PAR_FLUX')
  mol_avg_energy  = inGetData('MOL_AVG_ENERGY')
  psin            = inGetData('PSIN')
  rho             = inGetData('RHO')
;  l               = inGetData('L')
  jsat            = inGetData('JSAT')
  bratio          = inGetData('BRATIO')
  costet          = inGetData('COSTET')
  dens            = inGetData('NE')
  vb              = inGetData('VB')
  te              = inGetData('TE')
  ti              = inGetData('TI')

  total_rad_volume  = inGetData('TOT_RAD_VOL')
  total_rad_surface = inGetData('TOT_RAD_SUR')
  rad_impurity      = inGetData('RAD_IMP')
  rad_hydrogen      = inGetData('RAD_H')
  g_perp            = inGetData('G_PERP')

  inCloseInterface  

  atom_energy_flux = atom_par_flux * atom_avg_energy * 1.603E-19 * 1.0E-6
  mol_energy_flux  = mol_par_flux  * mol_avg_energy  * 1.603E-19 * 1.0E-6
  toroidal_area    = 2.0 * 3.141592 * r_cen * length

  dist = length
  dist[0] = 0.0
  FOR i = 1, N_ELEMENTS(length)-1 DO  $
    dist[i] = dist[i-1] + 0.5 * (length[i-1] + length[i])

  result = {  $
    version          : 1.0              ,  $
    file             : file             ,  $
    index            : index            ,  $
    index_target     : index_target     ,  $
    index_pin        : index_pin        ,  $
    index_cell       : index_cell       ,  $
    index_ring       : index_ring       ,  $
    index_tube       : index_tube       ,  $
    r_cen            : r_cen            ,  $
    z_cen            : z_cen            ,  $
    r_vertex1        : r_vertex1        ,  $
    z_vertex1        : z_vertex1        ,  $
    r_vertex2        : r_vertex2        ,  $
    z_vertex2        : z_vertex2        ,  $
    length           : length           ,  $
    dist             : dist             ,  $
    toroidal_area    : toroidal_area    ,  $
    atom_par_flux    : atom_par_flux    ,  $
    atom_avg_energy  : atom_avg_energy  ,  $
    atom_energy_flux : atom_energy_flux ,  $
    mol_par_flux     : mol_par_flux     ,  $  ; keeping units of D m-2 s-1 for now, rather than D2 m-2 s-1
    mol_avg_energy   : mol_avg_energy   ,  $  ;  (corrected depending on the plot)
    mol_energy_flux  : mol_energy_flux  ,  $  ; 
    psin             : psin             ,  $  ; 
    rho              : rho              ,  $  ; 
;    l                : l                ,  $  ; 
    jsat             : jsat             ,  $  ; 
    bratio           : bratio           ,  $  ; 
    costet           : costet           ,  $  ; 
    dens             : dens             ,  $  ; 
    vb               : vb               ,  $  ; 
    te               : te               ,  $  ; 
    ti               : ti               ,  $  ;
    total_rad_volume  : total_rad_volume  ,  $  ;
    total_rad_surface : total_rad_surface ,  $  ;
    rad_impurity      : rad_impurity      ,  $  ;
    rad_hydrogen      : rad_hydrogen      ,  $  ;
    g_perp            : g_perp            }

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

print,'file ',file

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
print,tube
print,n,m

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
  r        = inGetData('R') 
  z        = inGetData('Z') 
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
    r            : r            ,  $
    z            : z            ,  $
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
FUNCTION cortex_LoadSourceData, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  index   = inGetData('INDEX') 
  tube    = inGetData('TUBE' ) 
  s       = inGetData('S') 
  r       = inGetData('R') 
  z       = inGetData('Z') 
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
  eni_ion = inGetData('ENI_ION')
  eni_fit = inGetData('ENI_FIT')
  inCloseInterface  

  result = {             $
    verison : 1.0     ,  $
    file    : file    ,  $
    index   : index   ,  $
    tube    : tube    ,  $
    s       : s       ,  $
    r       : r       ,  $
    z       : z       ,  $
    par_net : par_net ,  $
    par_ion : par_ion ,  $ 
    par_rec : par_rec ,  $
    par_usr : par_usr ,  $
    par_ano : par_ano ,  $
    mom_net : mom_net ,  $
    mom_vol : mom_vol ,  $
    mom_usr : mom_usr ,  $ ; hello
    mom_ano : mom_ano ,  $
    ene_net : ene_net ,  $
    ene_ion : ene_ion ,  $
    ene_fit : ene_fit ,  $
    eni_net : eni_net ,  $
    eni_ion : eni_ion ,  $
    eni_fit : eni_fit }

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
  r     = inGetData('R') 
  z     = inGetData('Z') 
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
    r       : r     ,  $
    z       : z     ,  $
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
    ring    : ring     ,  $  ; index of the ring on the grid
    r       : r        ,  $  ; major radius of the centre of the midplane cell on the grid (approximate)
    rho     : rho      ,  $  ; radial distance from the LFS separatrix
    psin    : psin     ,  $  ; normalized magnetic flux (psin=1.0 at the separatrix)
    l       : l        ,  $  ; connection length (target-to-target distance) for the ring
    dens    : dens     ,  $  ; electron density
    vb      : vb       ,  $  ; parallel plasma flow velocity (m s-1)
    machno  : machno   ,  $  ; parallel Mach number
    p       : p        ,  $  ; total plasma pressure
    te      : te       ,  $  ; electron temperature
    ti      : ti  }          ; D+ ion temperature 

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
FUNCTION cortex_LoadEmissionData, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  wlngth = inGetData('WLNGTH')  
  cell   = inGetData('CELL')
  tube   = inGetData('TUBE')
  s      = inGetData('KSS')
  p      = inGetData('KPS')

  lines = MAKE_ARRAY(N_ELEMENTS(cell),N_ELEMENTS(wlngth),/FLOAT,VALUE=0.0)
  FOR i = 1, N_ELEMENTS(wlngth)  $
    DO lines[*,i-1] = inGetData('LINE_'+STRING(i,FORMAT='(I02)'))

  inCloseInterface

  result = {  $
    version : 1.0    ,  $
    file    : file   ,  $
    wlngth  : wlngth ,  $
    cell    : cell   ,  $
    tube    : tube   ,  $
    s       : s      ,  $
    p       : p      ,  $
    lines   : lines  }

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
  cell = inGetData('TUBE')  
  s    = inGetData('S')

  imp_data = MAKE_ARRAY(N_ELEMENTS(cell),div_max_iz+2,/FLOAT,VALUE=0.0)
  FOR i = 0, div_max_iz DO imp_data[*,i] = inGetData('IMP_DENS_'+STRING(i,FORMAT='(I02)'))
  FOR j = 0, N_ELEMENTS(cell)-1 DO imp_data[j,div_max_iz+1] = TOTAL(imp_data[j,1:div_max_iz])

  PRINT, 'WARNING cortex_LoadDIVIMPImpurityData_Density: Not summing neutral state, due to artifact with time-dependent runs with zero sputtering in the current time-step, where the dummy launch in DIVIMP gets counted, and it should not be'


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
    cell        : cell        ,  $
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

  imp_data = MAKE_ARRAY(N_ELEMENTS(tube),div_max_iz+2,/FLOAT,VALUE=0.0)
  FOR i = 0, z DO imp_data[*,i] = inGetData('IMP_IONIZ_'+STRING(i,FORMAT='(I02)'))
  FOR j = 0, N_ELEMENTS(tube)-1 DO imp_data[j,z+1] = TOTAL(imp_data[j,1:z])

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
