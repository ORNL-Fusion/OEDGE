;
; getb,'m-cry-0009a.celldata',machine='MAST',equ='~/divimp/dg/culham/longpx_2d_mastu_ctf4_0_snm_2.equ'
; getb,'j-bfg-0001h.celldata',machine='MAST',equ='~/divimp/dg/jet/JET_68124_49000.X4.equ'
; getb,'i-iwl-0002a.dat.cell',machine='EQU',equ='Baseline2008-li0.70.x4.equGRID_ITER_21JAN2001.equ'
; getb,'t-new-0001a.dat.cell',machine='EQU',equ='~/divimp/shots/ts2/600kA_LN_1cm/TS2_600kA_LN.x2.cr.equ'
; getb,'c-new-0000a.dat.cell',machine='EQU',equ='~/divimp/shots/cmod/1100303017_0138/cmod.1100303017.01380.x8.equ',/CHANGE_SIGN
; getb,'j-bfg-0006d.dat.cell',machine='EQU',equ='~/divimp/shots/jet/68124_49000/JET_68124_49000.equ',/CHANGE_SIGN
; getb,'j-bfg-0007d.dat.cell',machine='EQU',equ='~/divimp/shots/jet/68124_49000/JET_68124_49000.equ',/CHANGE_SIGN
; getb,'j-bfg-0006e.dat.cell',machine='EQU',equ='~/divimp/shots/jet/68124_49000/JET_68124_49000.equ',/CHANGE_SIGN
; getb,'m-det-0003a.dat.cell',machine='EQU',equ='~/divimp/shots/mast/24860/carre.24860_240.equ'
; getb,'m-ch4-0003a.objects.centre',machine='EQU',equ='~/divimp/shots/mast/13948/efm0139.48.280.equ'
;
; ======================================================================
;
FUNCTION ReadEquFile_OLD, filename, debug=debug

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
;
;   --------------------------------------------------------------------
;
  IF (psib EQ 0.0) THEN BEGIN
    PRINT,'Fising...'
    psib = 1.0
    psi = psi + psib
  ENDIF

  psi_min = MIN(psi) + psib

  psin = 1.0 / ((psi+psib-psi_min) / (psib-psi_min))  ; DIII-D definition (I think)

 print,min(psi+psib-psi_min),max(psi+psib-psi_min)

  result = {                     $       
    psi_boundary : psib       ,  $
    psi_min      : psi_min    ,  $
    x            : x          ,  $
    y            : y          ,  $
    psi_raw      : psi        ,  $
    psi          : psi + psib ,  $
    psin         : psin          }

  RETURN, result 
END
;
; ======================================================================
;
FUNCTION getb_LoadGridData, fpin, file_name, debug=debug

  result = -1

;
; ----------------------------------------------------------------------
;
  buffer = ' '
  OPENR,fpin,file_name, error=err
  WHILE NOT EOF(fpin) DO BEGIN 
    READF,fpin,buffer
    IF (STRPOS(buffer,'*') EQ -1) THEN BEGIN
      file_format = LONG(buffer)
      BREAK        
    ENDIF
  ENDWHILE
  n = 0
  WHILE NOT EOF(fpin) DO BEGIN 
    READF,fpin,buffer
    IF (STRPOS(buffer,'*' ) EQ -1 AND STRLEN(STRTRIM(buffer)) NE 0) THEN n = n + 1
    IF (STRPOS(buffer,'**') NE -1) THEN BEGIN
      PRINT,'ERROR getb_LoadGridData: Number format error'
      STOP
    ENDIF
  ENDWHILE 
  CLOSE,fpin
  IF (KEYWORD_SET(debug)) THEN PRINT,'FILE FORMAT, N =',file_format,n
;
; ----------------------------------------------------------------------
;
  CASE (file_format) OF
;   --------------------------------------------------------------------
    3: BEGIN
      ik   = MAKE_ARRAY(n,/LONG  ,VALUE=0)
      ir   = MAKE_ARRAY(n,/LONG  ,VALUE=0)
      rcen = MAKE_ARRAY(n,/DOUBLE,VALUE=0.0)
      zcen = MAKE_ARRAY(n,/DOUBLE,VALUE=0.0)
      rmid = MAKE_ARRAY(n,/DOUBLE,VALUE=0.0)
      zmid = MAKE_ARRAY(n,/DOUBLE,VALUE=0.0)
      r1   = MAKE_ARRAY(n,/DOUBLE,VALUE=0.0)
      z1   = MAKE_ARRAY(n,/DOUBLE,VALUE=0.0)
      OPENR,fpin,file_name, error=err
;     Read header:
      WHILE NOT EOF(fpin) DO BEGIN 
        READF,fpin,buffer
        IF (STRPOS(buffer,'*') EQ -1) THEN BREAK
      ENDWHILE
      i = -1
      WHILE NOT EOF(fpin) DO BEGIN 
        READF,fpin,buffer
        IF (STRPOS(buffer,'*') EQ -1 AND STRLEN(STRTRIM(buffer)) NE 0) THEN BEGIN
          buffer_array = STRSPLIT(buffer,/EXTRACT)
          i = i + 1
          ik[i] = (buffer_array[0])
          ir[i] = (buffer_array[1])
          rcen[i] = DOUBLE(buffer_array[2])
          zcen[i] = DOUBLE(buffer_array[3])
          rmid[i] = DOUBLE(buffer_array[4])
          zmid[i] = DOUBLE(buffer_array[5])
          r1  [i] = DOUBLE(buffer_array[6])
          z1  [i] = DOUBLE(buffer_array[7])
        ENDIF
      ENDWHILE
      CLOSE,fpin
      result = CREATE_STRUCT('n',n,'ik',ik,'ir',ir,      $
                             'x_cen',rcen,'y_cen',zcen,  $
                             'x_mid',rmid,'y_mid',zmid,  $
                             'x_1'  ,r1  ,'y_1'  ,z1     )
      END
;   --------------------------------------------------------------------
    ELSE: BEGIN
      PRINT,'ERROR GetB: Unknown file format'
      PRINT,' FILE_FORMAT= ',file_format
      STOP
      END
  ENDCASE

  RETURN, result
END
;
; ======================================================================
;
;  getb,'../../divimp/results/i-iwl-0003b.dat.cell',equ='../../divimp/shots/iter/111/GRID_ITER_21JAN2001.equ',machine='EQU'
;  getb,'../../divimp/results/a-new-0000b.dat.cell',equ='../../divimp/shots/asdex/22575_3000/22575_72x18.sonnet.equ',machine='EQU'
;
PRO getb,input_file_name,machine=machine,path=path,fast=fast,equ=equ,change_sign=change_sign

IF (NOT KEYWORD_SET(machine)) THEN machine = 'MAST'

; IDL routine to find Btotal and PSI_n for a list of given R,Z.
;
; from DIII-D code written by Adam McClean (U. Toronto)
;
; Declare variables.
fpin = 1
fpout = 2

; Make sure units are free.
free_lun,fpout
free_lun,fpin

; Open input file name.
  OPENR, fpin, input_file_name, error=err
  IF (err NE 0) THEN BEGIN
    print,"Error opening input file: ",input_file_name,err
    printf,-2,!error_state.msg
    free_lun,fpout
    free_lun,fpin
    stop
  ENDIF
  CLOSE,fpin
;
; ----------------------------------------------------------------------
;
  grid = getb_LoadGridData(fpin,input_file_name)
;
; ----------------------------------------------------------------------
;
CASE (machine) OF
; ----------------------------------------------------------------------
  'EQU': b = grid_ReadEquFile(equ,CHANGE_SIGN=change_sign)
; ----------------------------------------------------------------------
  'MAST': BEGIN
    print,shot
    print,time

    ndat = 0
    WHILE not eof(fpin) do begin
      ndat = ndat + 1
      READF,fpin,buffer
    ENDWHILE
    close, fpin
    print,ndat

;    ndat = 2

    ik = MAKE_ARRAY(ndat,/LONG ,VALUE=0)
    ir = MAKE_ARRAY(ndat,/LONG ,VALUE=0)
    r  = MAKE_ARRAY(ndat,/FLOAT,VALUE=0.0)
    z  = MAKE_ARRAY(ndat,/FLOAT,VALUE=0.0)
;    t  = time/1000.0
    t  = MAKE_ARRAY(ndat,/FLOAT,VALUE=time/1000.0)

    OPENR, fpin, input_file_name, error=err
    READF,fpin,buffer
    READF,fpin,buffer

    FOR i = 0, ndat-1 DO BEGIN
      READF,fpin,buffer
      buffer_array = strsplit(buffer,/extract)
      ik[i] = (buffer_array[0])
      ir[i] = (buffer_array[1])
      r [i] = float(buffer_array[2])
      z [i] = float(buffer_array[3])
; I think this was supposed to trigger the code to output BRATIO, etc. as well, see 0.99999 below...
;      IF (r[i] LT  0.06 OR r[i] GT 2.00 OR  $
;          z[i] LT -2.00 OR z[i] GT 2.00) THEN BEGIN  
;        r[i] = 0.99999
;        z[i] = 0.0
;      ENDIF
    ENDFOR

;    psi_n = get_psin(shoti,r,z,t,/mast)                  ; MAST
;    path = '/home/slisgo/divimp/dg/culham/'  ; For special improved TRANSP EFITS


;    ENDIF ELSE BEGIN
;      res = flux_coord_sl(shoti,r,z,t,/mast,path=path,fast=fast)
;      help,res,/struct
;      psi_n   = res.normflux
;      IF (NOT KEYWORD_SET(fast)) THEN BEGIN
;        btotal  = SQRT(res.br^2 + res.bphi^2 + res.bz^2)  
;        bpol_cl = res.btheta
;        bpol    = SQRT(res.br^2 + res.bz^2)
;        bratio  = bpol / btotal  
;        bfield  = [0.0,0.0]
;      ENDIF
;    ENDELSE

    FOR i = 0, ndat-1 DO BEGIN       
      IF (r[i] EQ 0.99999 AND z[i] EQ 0.0) THEN BEGIN
        psi_n [i] = 0.0
        bratio[i] = 0.0
        btotal[i] = 0.0
        res.br  [i] = 0.0
        res.bphi[i] = 0.0
        res.bz  [i] = 0.0
      ENDIF

      IF (KEYWORD_SET(fast) OR KEYWORD_SET(equ)) THEN BEGIN
        printf,fpout,r[i],z[i],0.0,psi_n[i],ik[i],ir[i],  $
                       format='(4F10.6,2I6)'
      ENDIF ELSE BEGIN
        printf,fpout,r[i],z[i],bratio[i],psi_n[i],ik[i],ir[i],  $
                       bpol_cl[i],bpol[i],btotal[i],              $
                       res.br[i],res.bphi[i],res.bz[i],           $
                       format='(4F10.6,2I6,2X,3F12.8,2X,3F12.7)'
      ENDELSE
    ENDFOR

    END
; ----------------------------------------------------------------------
  'D3D' : BEGIN
; Loop until EOF.
    while not eof(fpin) do begin

	; Begin reading cell/ring number and central R,Z.
	readf,fpin,buffer
	buffer_array = strsplit(buffer,/extract)
	ik = (buffer_array[0])
	ir = (buffer_array[1])
	r_m = float(buffer_array[2])
	z_m = float(buffer_array[3])
	
	; Determine Bpol and Btor at gizen r,z.
        bfield = [0.0,0.0]
;	bfield = calculate_bprz(r_m,z_m,a,g,grid=grid)   ! DIII-D ???
	
	; Calculate Bratio = Bpol / Btot
        bratio = 0.0
;	bratio = abs(bfield[0]) / abs(bfield[1])

	; Calculate normalized PSI at r,z.
        rin = [ 0.613]
        zin = [-0.973]
        tin = [ 0.220] 
        rin = [r_m]
        zin = [z_m]
        tin = [time/1000.0]
;        psi_n = efit_rz2rho(rin,zin,tin,shot=980116027,/PSINORM)   ; C-Mod
;        psi_n = get_psin(shoti,rin,zin,tin,/mast)                  ; MAST
;    	 psi_n = calculate_psi(r_m,z_m,a,g,/norm)                   ; DIII-D

        print,"psin",psi_n

	; Write output line to file.
	printf,fpout,r_m,z_m,bratio,psi_n,ik,ir,abs(bfield[0]),abs(bfield[1]),format='(4F10.6,2I6,2X,2F10.6)'
    endwhile

          END
ENDCASE


;help,grid,/struct
;help,b,/struct
; contour,b.psin,b.x,b.y,levels=[0.9,1.0,1.1],c_annotation=''
;stop
;
; ----------------------------------------------------------------------
;
  psin = MAKE_ARRAY(grid.n,/DOUBLE,VALUE=0.0)
  FOR i = 0, grid.n-1 DO BEGIN
    psin[i] = Map2D(b.psin,b.x,b.y,grid.x_mid[i],grid.y_mid[i],cubic=-0.5)
  ENDFOR
  grid = CREATE_STRUCT(grid,'psin',psin)

  b_ratio = MAKE_ARRAY(grid.n,/DOUBLE,VALUE=0.0)
  FOR i = 0, grid.n-1 DO BEGIN
    b_ratio[i] = Map2D(b.b_ratio,b.x,b.y,grid.x_cen[i],grid.y_cen[i],cubic=-0.5)
  ENDFOR
  grid = CREATE_STRUCT(grid,'b_ratio',b_ratio)

  psi = MAKE_ARRAY(grid.n,/DOUBLE,VALUE=0.0)
  FOR i = 0, grid.n-1 DO BEGIN
    psi[i] = Map2D(b.psi,b.x,b.y,grid.x_mid[i],grid.y_mid[i],cubic=-0.5)
  ENDFOR
  grid = CREATE_STRUCT(grid,'psi',psi)

;help,grid
;
; ----------------------------------------------------------------------
;


;help,grid,/struct



  outfile = input_file_name+'.out'
  OPENW, fpout, outfile, ERROR=err
  IF (err NE 0) THEN BEGIN
    PRINT,'ERROR GetB: Problem opening output file'
    PRINT,'FILE_NAME= ',outfile
    PRINTF,-2,!error_state.msg
    FREE_LUN,fpin
    FREE_LIN,fpout
    STOP
  ENDIF
; Write version number:
  PRINTF,fpout,'1.00'
; Write new table descriptions in fpout.
  IF (KEYWORD_SET(fast)) THEN BEGIN
    printf,fpout,"R(m)","Z(m)","PSI_n",  $
                   format='(3A10)'
  ENDIF ELSE BEGIN
;    printf,fpout,"R(m)","Z(m)","Bratio","PSI_n","IK","IR",  $
;                   "Bpol_cl","Bpol","Btot","Br","Bphi","Bz",   $
;                   format='(/4A10,2A6,2X,3A12,2X,3A12)'
    PRINTF,fpout,'* ','IK','IR','PSIn','B_ratio',  $
                   format='(A2,2A6,2A14)'


    FOR i = 0, grid.n-1 DO BEGIN
      PRINTF,fpout,grid.ik[i],grid.ir[i],grid.psin[i], grid.b_ratio[i], $
                   grid.psi[i],  $  ; *** TEMP ***
                   FORMAT='(2X,2I6,3F14.7)'
    ENDFOR

  ENDELSE






; Close output file unit.
FREE_LUN, fpin
FREE_LUN, fpout
;exit
END
