; b=grid_readequfile('/home/ITER/lisgos/fuse_data/mast/shots/13018_250/13018_250_mod.equ',/debug)
; b=grid_readequfile('/home/ITER/lisgos/divimp/shots/iter/i1514/Baseline2008-li0.70.x4.equ',/debug)
; b=grid_readequfile('/home/ITER/lisgos/divimp/shots/mast/24860/carre.24860_240.equ',/debug)
;
; ======================================================================
;
FUNCTION grid_GetComponents, b, debug=debug

; A lot of this is overdone in a vain attempt to get something useful out of
; fluxcoordinates -- which I've abandoned for now, but am leaving all the
; infrastructure in for some future attempt to re-visit what's going on...

  nr = N_ELEMENTS(b.x)
  nz = N_ELEMENTS(b.y)
  nt = 2

  r2D = MAKE_ARRAY(nr,nz,/DOUBLE,VALUE=0.0D)  
  z2D = MAKE_ARRAY(nr,nz,/DOUBLE,VALUE=0.0D)  
  FOR i = 0, nz-1 DO r2D[*,i] = b.x
  FOR i = 0, nr-1 DO z2D[i,*] = b.y

  r       = MAKE_ARRAY(nr      ,/DOUBLE,VALUE=0.0D)
  z       = MAKE_ARRAY(nz      ,/DOUBLE,VALUE=0.0D)
  t       = MAKE_ARRAY(      nt,/DOUBLE,VALUE=0.0D)
  flux    = MAKE_ARRAY(nr,nz,nt,/DOUBLE,VALUE=0.0D)
  lcfs    = MAKE_ARRAY(      nt,/DOUBLE,VALUE=0.0D)
  axisval = MAKE_ARRAY(      nt,/DOUBLE,VALUE=0.0D)
  axispos = MAKE_ARRAY(2,    nt,/DOUBLE,VALUE=0.0D)

  psi_max = MAX(b.psi,imax)

  r = b.x  
  z = b.y  

  FOR it = 0, nt-1 DO BEGIN
    t      [    it] = 0.0
    flux   [*,*,it] = b.psi[*,*]
    lcfs   [    it] = b.psi_boundary
    axisval[    it] = b.psi[imax]
    axispos[0,  it] = r2D  [imax]
    axispos[1,  it] = z2D  [imax]
  ENDFOR

  B_phi = b.B_phi_ref * b.R_phi_ref / r2D

; For some reason, this doesn't work when calculating Bratio.  I checked the calculated
; B_pol and B_phi against b=read_flux(13018) at Culham and they seem OK, but I can't seem
; to get from there to B_total -- note, comparing the _fc quanities with the R,Z ones
; calculated below seems to be OK on the midplane, but funny things happen elsewhere --
; I dunno, but I "wasted" a lot of time on this... - SL, 27/01/2010
;  a = fluxcoordinates(r,z,t,flux,lcfs,axisval,axispos,/unit)
;  btheta = a.norm.cov[*,*,*,0] * a.norm.cov[*,*,*,2] * a.jacobian  
;  for i=0,n_elements(btheta[0,0,*])-1 do btheta[*,*,i]=btheta[*,*,i]*a.norm.psi[i]
;  B_pol_fc    = a.metriccov[*,*,*,1,1] * btheta
;  B_psi_fc    = a.metriccov[*,*,*,0,1] * btheta
;  B_total_fc  = SQRT(B_pol_fc[*,*,0]^2 + B_psi_fc[*,*,0]^2 + B_phi^2)
;  B_ratio_fc  = ABS(B_pol_fc) / B_total_fc


  grad_flux = gradcv(r,z,flux)  ; gradcv: 0=r, 1=psi, 2=z

  B_r = MAKE_ARRAY(nr,nz,nt,/DOUBLE,VALUE=0.0D)
  B_z = MAKE_ARRAY(nr,nz,nt,/DOUBLE,VALUE=0.0D)
  FOR k = 0, nr-1 DO BEGIN 
    B_r[k,*,*] = -grad_flux[k,*,*,2] / r[k]  ; -1/R(dpsi/dZ)
    B_z[k,*,*] =  grad_flux[k,*,*,0] / r[k]  ;  1/R(dpsi/dR)
  ENDFOR

  B_pol = SQRT(B_r^2 + B_z^2)
   
  B_total = SQRT(B_r[*,*,0]^2 + B_z[*,*,0]^2 + B_phi^2)

  B_ratio = ABS(B_pol[*,*,0]) / B_total

  result = CREATE_STRUCT(b ,  $
;    'B_pol_fc'   , B_pol_fc  [*,*,0] ,  $
;    'B_psi_fc'   , B_psi_fc  [*,*,0] ,  $
;    'B_total_fc' , B_total_fc        ,  $
    'B_phi'      , B_phi             ,  $
    'B_r'        , B_r       [*,*,0] ,  $
    'B_z'        , B_z       [*,*,0] ,  $
    'B_pol'      , B_pol     [*,*,0] ,  $
    'B_total'    , B_total           ,  $
    'B_ratio'    , B_ratio           )

  RETURN, result

END
;
; ======================================================================
;
FUNCTION grid_ReadEquFile, filename, debug=debug, change_sign=change_sign

  IF (KEYWORD_SET(debug)) THEN debug = 1 ELSE debug = 0
;
;   --------------------------------------------------------------------
;
  fp = 3
  FREE_LUN, fp
  OPENR, fp, filename, error=err
  IF (err NE 0) THEN BEGIN
    PRINT,'ERROR: Unable to open EQU file -' + filename
    RETURN, -1
  ENDIF

  buffer = ' '
  
  WHILE (1) DO BEGIN

    READF,fp,buffer

    IF (STRMATCH(buffer,'*jm    =*') EQ 1) THEN BEGIN
      str = STRSPLIT(buffer,' ',/EXTRACT)
      jm = LONG(str[2])
      IF (debug) THEN PRINT,'found jm',jm
    ENDIF

    IF (STRMATCH(buffer,'*km    =*') EQ 1) THEN BEGIN
      str = STRSPLIT(buffer,' ',/EXTRACT)
      km = LONG(str[2])
      IF (debug) THEN PRINT,'found km',km
    ENDIF

    IF (STRMATCH(buffer,'*psib  =*') EQ 1) THEN BEGIN
      str = STRSPLIT(buffer,' ',/EXTRACT)
      psib = DOUBLE(str[2])
      IF (debug) THEN PRINT,'found psib',psib
    ENDIF

    IF (STRMATCH(buffer,'*btf   =*') EQ 1) THEN BEGIN
      str = STRSPLIT(buffer,' ',/EXTRACT)
      btf = DOUBLE(str[2])
      IF (debug) THEN PRINT,'found btf',btf
    ENDIF

    IF (STRMATCH(buffer,'*rtf   =*') EQ 1) THEN BEGIN
      str = STRSPLIT(buffer,' ',/EXTRACT)
      rtf = DOUBLE(str[2])
      IF (debug) THEN PRINT,'found rtf',rtf
    ENDIF

    IF (STRMATCH(buffer,'*r(1:jm)*') EQ 1) THEN BEGIN
      x = MAKE_ARRAY(jm,/DOUBLE,VALUE=0.0)
      FOR i = 0, jm-1, 5 DO BEGIN
        READF,fp,buffer
        str = STRSPLIT(buffer,' ',/EXTRACT)                   
        x[i:MIN([i+4,jm-1])] = DOUBLE(str)
      ENDFOR
      IF (debug) THEN PRINT,'found R',x[0],x[jm-1]
    ENDIF

    IF (STRMATCH(buffer,'*z(1:km)*') EQ 1) THEN BEGIN
      y = MAKE_ARRAY(km,/DOUBLE,VALUE=0.0)
      FOR i = 0, km-1, 5 DO BEGIN
        READF,fp,buffer
        str = STRSPLIT(buffer,' ',/EXTRACT)                   
        y[i:MIN([i+4,km-1])] = DOUBLE(str)
      ENDFOR
      IF (debug) THEN PRINT,'found Z',y[0],y[km-1]
    ENDIF

    IF (STRMATCH(buffer,'*psi(j,k)*') EQ 1) THEN BEGIN
      psi = MAKE_ARRAY(jm,km,/DOUBLE,VALUE=0.0)
      j = -1
      k =  0
      WHILE (1) DO BEGIN                
        READF,fp,buffer
        str = STRSPLIT(buffer,' ',/EXTRACT)                   
        FOR i = 0, N_ELEMENTS(str)-1 DO BEGIN
          j = j + 1 
          IF (j EQ jm) THEN BEGIN
            j = 0
            k = k + 1
          ENDIF
          psi[j,k] = DOUBLE(str[i])
        ENDFOR
        IF (j EQ jm -1 AND k EQ km-1) THEN BREAK
      ENDWHILE
      IF (debug) THEN PRINT,'found psi',psi[0,0],psi[jm-1,km-1]
      BREAK
    ENDIF

  ENDWHILE

  FREE_LUN, fp
;
;   --------------------------------------------------------------------
;

  psi_raw = psi

  IF (KEYWORD_SET(change_sign)) THEN psi = -psi

;  IF (psib EQ 0.0D) THEN BEGIN
;    psib = 1.0D
;    psi = psi + psib
;  ENDIF
;  psi_min = MIN(psi) + psib
;  psin = 1.0D / ((psi+psib-psi_min) / (psib-psi_min))  ; DIII-D definition (I think)

; psin: normalized psi sqrt((psi-psi_0)/(psi_s-psi_0))  Culham definition in fluxcoordinates.pro from H. Meyer
  psi_0 = MAX(psi) 
  psi_s = psib
  psin = ((psi - psi_0) / (psi_s - psi_0))

;  print,min(psi+psib-psi_min),max(psi+psib-psi_min)

  b = {                     $       
    psi_boundary : psib       ,  $
;    psi_min      : psi_min    ,  $
    x            : x          ,  $
    y            : y          ,  $
    B_phi_ref    : btf        ,  $
    R_phi_ref    : rtf        ,  $
    psi_raw      : psi_raw    ,  $
    psi          : psi        ,  $
    psin         : psin          }

  result = grid_GetComponents(b,debug=debug)

  RETURN, result 
END
