;
; ======================================================================
;
PRO cortex_PageTitle, plot, ps, file, dev_xsize, dev_ysize, title, notes, charsize

  ypos = 0.940
  IF (notes NE 'default') THEN ypos = 0.965

  IF (file NE 'none') THEN BEGIN
    str = file                       ; extract the case name from the data file
    str = STRSPLIT(str,'/',/EXTRACT)
    str = str[N_ELEMENTS(str)-1]     ; take the last sub-string
    str = STRSPLIT(str,'.',/EXTRACT)
    str = str[0]                     ; take the first one
  ENDIF

  n = STRLEN(title) 
  nmax = LONG(60.0 * 2.0 / (1.5 * charsize))
  IF (n GT nmax) THEN BEGIN
    i = STRPOS(STRMID(title,0,nmax),' ',/REVERSE_SEARCH)
    str1 = STRMID(title,0  ,i   )
    str2 = STRMID(title,i+1,nmax)
    XYOUTS, 0.02 * dev_xsize, (ypos + 0.045 * charsize / 1.7) * dev_ysize, CHARSIZE=1.5*charsize, str1, /DEVICE
    XYOUTS, 0.02 * dev_xsize, (ypos                         ) * dev_ysize, CHARSIZE=1.5*charsize, str2, /DEVICE
  ENDIF ELSE BEGIN
    XYOUTS, 0.02 * dev_xsize, ypos * dev_ysize, CHARSIZE=1.5*charsize, title, /DEVICE
  ENDELSE

  IF (notes NE 'unknown') THEN BEGIN
    n = STRLEN(notes) 
    nmax = 105
    IF (n GT nmax) THEN BEGIN
      i = STRPOS(STRMID(notes,0,nmax),' ',/REVERSE_SEARCH)
      str1 = STRMID(notes,0,nmax)
      str2 = STRMID(notes,nmax+1,n)
      XYOUTS, 0.02 * dev_xsize, (ypos - 0.025) * dev_ysize, CHARSIZE=1.0, str1, /DEVICE
      XYOUTS, 0.02 * dev_xsize, (ypos - 0.050) * dev_ysize, CHARSIZE=1.0, str2, /DEVICE
      ypos = ypos - 0.020
    ENDIF ELSE BEGIN
      XYOUTS, 0.02 * dev_xsize, (ypos - 0.025) * dev_ysize, CHARSIZE=1.0, notes, /DEVICE
    ENDELSE              
  ENDIF

  RETURN
END
;
; ======================================================================
;
FUNCTION cortex_ExtractStructure, struct, index

  result = {error : -1} 

; Improve this routine (vastly) with this IDL feature...
;  struct = cortex_ExtractStructure(plot_array,idata)                 
;  struct = struct.node
;  index = WHERE(TAG_NAMES(struct) EQ STRUPCASE(tag), count)
;  IF (count NE 0) THEN BEGIN
;    val = struct.(index[0]) 
;  ENDIF ELSE BEGIN
;    PRINT, 'ERROR cortex_PlotNodes: TAG not found'
;    PRINT, '  TAG = ',tag
;    RETURN
;  ENDELSE

  IF (index LT 1 OR index GT 40) THEN BEGIN
    PRINT,'ERROR cortrex_ExtractStructure: Index out of bounds'
    PRINT,' INDEX  = ',index
    PRINT,' STRUCT = '
    HELP,struct,/STRUCT
    STOP
  ENDIF

  IF (index EQ 1 ) THEN result = struct.data1
  IF (index EQ 2 ) THEN result = struct.data2
  IF (index EQ 3 ) THEN result = struct.data3
  IF (index EQ 4 ) THEN result = struct.data4
  IF (index EQ 5 ) THEN result = struct.data5
  IF (index EQ 6 ) THEN result = struct.data6
  IF (index EQ 7 ) THEN result = struct.data7
  IF (index EQ 8 ) THEN result = struct.data8
  IF (index EQ 9 ) THEN result = struct.data9
  IF (index EQ 10) THEN result = struct.data10
  IF (index EQ 11) THEN result = struct.data11
  IF (index EQ 12) THEN result = struct.data12
  IF (index EQ 13) THEN result = struct.data13
  IF (index EQ 14) THEN result = struct.data14
  IF (index EQ 15) THEN result = struct.data15
  IF (index EQ 16) THEN result = struct.data16
  IF (index EQ 17) THEN result = struct.data17
  IF (index EQ 18) THEN result = struct.data18
  IF (index EQ 19) THEN result = struct.data19
  IF (index EQ 20) THEN result = struct.data20
  IF (index EQ 21) THEN result = struct.data21
  IF (index EQ 22) THEN result = struct.data22
  IF (index EQ 23) THEN result = struct.data23
  IF (index EQ 24) THEN result = struct.data24
  IF (index EQ 25) THEN result = struct.data25
  IF (index EQ 26) THEN result = struct.data26
  IF (index EQ 27) THEN result = struct.data27
  IF (index EQ 28) THEN result = struct.data28
  IF (index EQ 29) THEN result = struct.data29
  IF (index EQ 30) THEN result = struct.data30
  IF (index EQ 31) THEN result = struct.data31
  IF (index EQ 32) THEN result = struct.data32
  IF (index EQ 33) THEN result = struct.data33
  IF (index EQ 34) THEN result = struct.data34
  IF (index EQ 35) THEN result = struct.data35
  IF (index EQ 36) THEN result = struct.data36
  IF (index EQ 37) THEN result = struct.data37
  IF (index EQ 38) THEN result = struct.data38
  IF (index EQ 39) THEN result = struct.data39
  IF (index EQ 40) THEN result = struct.data40

  IF (N_ELEMENTS(TAG_NAMES(result)) EQ 1) THEN BEGIN
    PRINT, 'ERROR Cortex ExtractStructure: Structure element not found'
    PRINT,'  INDEX  = ',index
    PRINT,'  STRUCT = '
    HELP,struct,/struct     
    STOP
  ENDIF

  RETURN, result
END
