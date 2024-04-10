;
; ======================================================================
;
PRO adas_FractionalAbundance_Weighted

  color = [ 'Red', 'Blue' ]
  path  = '/work/projects/adas/'

  uid  = 'adas'
  elem = [ 'he' , 'c'  , 'fe' , 'w' ]
  iz0  = [  2   ,  6   ,  26  ,  74 ]
  year = [  96  ,  96  ,  96  ,  88 ]

  temp = adas_vector(low=1, high=20, num=100)
  dens = fltarr(100) + 1e12

  atom  = 'c'
  state = 2

  i = (WHERE(elem EQ atom,count))[0]
  IF (count EQ 0) THEN BEGIN
    PRINT,'ERROR adas_FractionalAbundance_Weighted: Unknown element'
    PRINT,'  ATOM = ',atom
    STOP
  ENDIF

  CASE elem[i] OF
    ; ------------------------------------------------------------------
    'he': BEGIN
      CASE state OF 
        0: BEGIN
          file = path + 'adas/adf15/pec96#he/pec96#he_pju#he0.dat'
          line  = 668
          CASE line OF  ; pju96   pju08   pjr96? ; pju96
            447: block = 11    ; 6     ; 26	  11    
            588: block = 9     ; 1     ; 24	  9     
            668: block = 13    ; 4     ; 28	  13    
            707: block = 8     ; 2     ; 23	  8     
            728: block = 12    ; 5     ; 27	  12    
            ELSE: BEGIN
              PRINT,'ERROR adas_LoadPEC: Unrecognised line'
              PRINT,'  LINE = ',line
              STOP
              END
          ENDCASE
          END
        1: BEGIN
          file = path + 'adas/adf15/pec96#he/pec96#he_pju#he1.dat'
          line = 468
          CASE line OF  ; pju96   pju08   pjr96
            468: block = 8      ; 6     ; 26
            ELSE: BEGIN
              PRINT,'ERROR adas_LoadPEC: Unrecognised line'
              PRINT,'  LINE = ',line
              STOP
              END
          ENDCASE
          END
        ELSE: BEGIN
          PRINT,'ERROR adas_LoadPEC: Invalid charge state'
          PRINT,'  STATE = ',state
          STOP
          END
      ENDCASE
      END
    ; ------------------------------------------------------------------
    'c': BEGIN
      CASE state OF 
        1: BEGIN
          file = path + 'adas/adf15/pec93#c/pec93#c_pju#c1.dat'
          line  = 514
          CASE line OF  ; pju93
            514: block = 16
            ELSE: BEGIN
              PRINT,'ERROR adas_LoadPEC: Unrecognised line'
              PRINT,'  LINE = ',line
              STOP
              END
          ENDCASE
          END
        2: BEGIN
          file = path + 'adas/adf15/pec93#c/pec93#c_pju#c2.dat'
          line = 465
          CASE line OF  ; pju93 
            465: block = 15  
            ELSE: BEGIN
              PRINT,'ERROR adas_LoadPEC: Unrecognised line'
              PRINT,'  LINE = ',line
              STOP
              END
          ENDCASE
          END
        ELSE: BEGIN
          PRINT,'ERROR adas_LoadPEC: Invalid charge state'
          PRINT,'  STATE = ',state
          STOP
          END
      ENDCASE
      END
    ; ------------------------------------------------------------------
    ELSE: BEGIN
      PRINT,'ERROR adas_LoadPEC: Unrecognised element'
      PRINT,'  ELEM = ',elem[i]
      STOP
      END
  ENDCASE


;  read_adf11,uid='adas',year=96,iz0=1,iz1=1,class='scd',te=temp,dens=dens*1.E-06,data=adas_data

;  file = '/home/adas/adas/adf15/pec96#h/pec96#h_pju#h0.dat'
;  file = './pec08_pju#he0.dat'
;  file = '/home/adas/adas/adf15/pec96#he/pec96#he_pju#he0.dat'
;  file = '/home/adas/adas/adf15/pec96#he/pec96#he_pjr#he0.dat'


  read_adf15,file=file,block=block,te=temp,dens=dens,data=data,wlngth=wlngth

  print,'ADAS:',file,block,wlngth

  PLOT, temp, data / MAX(data), XSTYLE=1, YSTYLE=1

  IF (atom EQ 'w') THEN BEGIN
    run_adas405, uid=uid, year=year[i], elem=elem[i], te=temp, dens=dens, $
                 frac=frac, files={acd:'/home/ITER/omullam/adas/adf11/acd88/acd88_w.dat',  $
                                   scd:'/home/ITER/omullam/adas/adf11/scd88/scd88_w.dat',  $
                                   plt:'/home/ITER/omullam/adas/adf11/plt88/plt88_w.dat',  $
                                   prb:'/home/ITER/omullam/adas/adf11/prb88/prb88_w.dat'}
  ENDIF ELSE  $
    run_adas405, uid=uid, year=year[i], elem=elem[i], te=temp, dens=dens, frac=frac

;  PLOT, [1,20000], [0.01,1.1], /nodata, XSTYLE=1, ystyle=1, CHARSIZE=2.0, /XLOG

  FOR j = 0, iz0[i] DO BEGIN

    OPLOT, temp, frac.ion[*,j], color=Truecolor(color[j MOD 2])

    IF ((state EQ 1 AND j EQ 1 AND (atom EQ 'he' OR atom EQ 'c')) OR  $
        (state EQ 2 AND j EQ 2 AND (                atom EQ 'c'))) THEN  $
      OPLOT, temp, (frac.ion[*,j] * data) / (MAX(frac.ion[*,j] * data)), color=Truecolor('Yellow')

  ENDFOR

  XYOUTS, 0.5, 0.9, elem[i], /NORMAL, CHARSIZE=2.0, color=Truecolor('White')

stop



END