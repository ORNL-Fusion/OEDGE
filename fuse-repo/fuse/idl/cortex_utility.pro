;
; ======================================================================
; 
;  from here: http://www.idlcoyote.com/tips/variable_undefine.html
;
   PRO cortex_Undefine, varname  
   tempvar = SIZE(TEMPORARY(varname))
   END
;
; ======================================================================
; 
FUNCTION cortex_GetValues, str, values

  IF (str EQ 'none') THEN RETURN, 0

  result = 0

  str_comma = STRSPLIT(str,',',/EXTRACT)
  
  values = [-1]

  FOR i = 0, N_ELEMENTS(str_comma)-1 DO BEGIN
    j = STRPOS(str_comma[i],'-')
    IF (j EQ -1) THEN BEGIN
      values = [values,FLOAT(str_comma[i])]
    ENDIF ELSE BEGIN
      str_dash = STRSPLIT(str_comma[i],'-',/EXTRACT)
      print, 'not ready'
      print,str_comma[i]
      print,str_dash
      stop
;       IF (LONG(val) GE LONG(str_dash[0]) AND  $
;          LONG(val) LE LONG(str_dash[1])) THEN result = 1
    ENDELSE
  ENDFOR

  IF (N_ELEMENTS(values) GT 1) THEN BEGIN
    values = values[1:N_ELEMENTS(values)-1]
    RETURN, 1
  ENDIF ELSE RETURN, 0

END
;
; ======================================================================
; 
FUNCTION cortex_CheckIndex, val, str

  IF (str EQ 'none') THEN RETURN, 0
  IF (str EQ 'all' ) THEN RETURN, 1

  result = 0

  str_comma = STRSPLIT(str,',',/EXTRACT)
  
  FOR i = 0, N_ELEMENTS(str_comma)-1 DO BEGIN
    j = STRPOS(str_comma[i],'-')
    IF (j EQ -1) THEN BEGIN
      IF (LONG(val) EQ LONG(str_comma[i])) THEN result = 1
    ENDIF ELSE BEGIN
      str_dash = STRSPLIT(str_comma[i],'-',/EXTRACT)
      IF (LONG(val) GE LONG(str_dash[0]) AND  $
          LONG(val) LE LONG(str_dash[1])) THEN result = 1
    ENDELSE
    IF (result EQ 1) THEN BREAK
  ENDFOR

  RETURN, result
END

;
; ======================================================================
;
FUNCTION cortex_UpdateFile, file

 COMMON options, dir_structure


 IF (dir_structure EQ 0) THEN RETURN, file


 i = STRPOS(file,'/',/REVERSE_SEARCH)

 str1 = STRMID(file,0,i+1)
 str2 = STRMID(file,i+1)

 str = STRSPLIT(str2,'-',/EXTRACT)


 family = str[0] + '-' + str[1] + '/'
; family = STRMID(str2,0,5) + '/'
 child  = STRMID(str2,0,7) + '/'

; PRINT,str1
; PRINT,family
; PRINT,child
; PRINT,str2

 CASE dir_structure OF
   1: file = str1 + family + str2
   2: file = str1 + family + child + str2
 ENDCASE

; PRINT,file

; PRINT, 'HERE in UIPDATE FILE'
; STOP
 RETURN,file

END
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
  nmax = LONG(75.0 * 2.0 / (1.5 * charsize))
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
    nmax = 125
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
FUNCTION cortex_ExtractStructure, struct, tag, name=name

  IF (KEYWORD_SET(name)) THEN BEGIN
    tag_local = tag
  ENDIF ELSE BEGIN
    tag_local = 'data' + STRTRIM(STRING(tag),2)
  ENDELSE  

  index = WHERE(TAG_NAMES(struct) EQ STRUPCASE(tag_local), count)

  IF (count NE 0) THEN BEGIN
    result = struct.(index[0]) 
  ENDIF ELSE BEGIN
    PRINT, 'ERROR grid_ExtractStructure: TAG not found'
    PRINT, '  TAG = >'+STRTRIM(tag_local,2)+'<'
    PRINT, '  NAMES = ',TAG_NAMES(struct)
;    HELP, struct, /struct
;    HELP, struct.contour12,/struct
    result = {error:-1}
  ENDELSE

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_ExtractStructure_old, struct, index

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

  IF (index LT 1 OR index GT 109) THEN BEGIN                        ; *** FIX THIS MINDLESS IMPLEMENTATION! ***
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
  IF (index EQ 41) THEN result = struct.data41
  IF (index EQ 42) THEN result = struct.data42
  IF (index EQ 43) THEN result = struct.data43
  IF (index EQ 44) THEN result = struct.data44
  IF (index EQ 45) THEN result = struct.data45
  IF (index EQ 46) THEN result = struct.data46
  IF (index EQ 47) THEN result = struct.data47
  IF (index EQ 48) THEN result = struct.data48
  IF (index EQ 49) THEN result = struct.data49
  IF (index EQ 50) THEN result = struct.data50
  IF (index EQ 51) THEN result = struct.data51
  IF (index EQ 52) THEN result = struct.data52
  IF (index EQ 53) THEN result = struct.data53
  IF (index EQ 54) THEN result = struct.data54
  IF (index EQ 55) THEN result = struct.data55
  IF (index EQ 56) THEN result = struct.data56
  IF (index EQ 57) THEN result = struct.data57
  IF (index EQ 58) THEN result = struct.data58
  IF (index EQ 59) THEN result = struct.data59
  IF (index EQ 60) THEN result = struct.data60
  IF (index EQ 61) THEN result = struct.data61
  IF (index EQ 62) THEN result = struct.data62
  IF (index EQ 63) THEN result = struct.data63
  IF (index EQ 64) THEN result = struct.data64
  IF (index EQ 65) THEN result = struct.data65
  IF (index EQ 66) THEN result = struct.data66
  IF (index EQ 67) THEN result = struct.data67
  IF (index EQ 68) THEN result = struct.data68
  IF (index EQ 69) THEN result = struct.data69
  IF (index EQ 70) THEN result = struct.data70
  IF (index EQ 71) THEN result = struct.data71
  IF (index EQ 72) THEN result = struct.data72
  IF (index EQ 73) THEN result = struct.data73
  IF (index EQ 74) THEN result = struct.data74
  IF (index EQ 75) THEN result = struct.data75
  IF (index EQ 76) THEN result = struct.data76
  IF (index EQ 77) THEN result = struct.data77
  IF (index EQ 78) THEN result = struct.data78
  IF (index EQ 79) THEN result = struct.data79
  IF (index EQ 80) THEN result = struct.data80
  IF (index EQ 81) THEN result = struct.data81
  IF (index EQ 82) THEN result = struct.data82
  IF (index EQ 83) THEN result = struct.data83
  IF (index EQ 84) THEN result = struct.data84
  IF (index EQ 85) THEN result = struct.data85
  IF (index EQ 86) THEN result = struct.data86
  IF (index EQ 87) THEN result = struct.data87
  IF (index EQ 88) THEN result = struct.data88
  IF (index EQ 89) THEN result = struct.data89
  IF (index EQ 90) THEN result = struct.data90
  IF (index EQ 91) THEN result = struct.data91
  IF (index EQ 92) THEN result = struct.data92
  IF (index EQ 93) THEN result = struct.data93
  IF (index EQ 94) THEN result = struct.data94
  IF (index EQ 95) THEN result = struct.data95
  IF (index EQ 96) THEN result = struct.data96
  IF (index EQ 97) THEN result = struct.data97
  IF (index EQ 98) THEN result = struct.data98
  IF (index EQ 99) THEN result = struct.data99
  IF (index EQ 100) THEN result = struct.data100
  IF (index EQ 101) THEN result = struct.data101
  IF (index EQ 102) THEN result = struct.data102
  IF (index EQ 103) THEN result = struct.data103
  IF (index EQ 104) THEN result = struct.data104
  IF (index EQ 105) THEN result = struct.data105
  IF (index EQ 106) THEN result = struct.data106
  IF (index EQ 107) THEN result = struct.data107
  IF (index EQ 108) THEN result = struct.data108
  IF (index EQ 109) THEN result = struct.data109

  IF (N_ELEMENTS(TAG_NAMES(result)) EQ 1) THEN BEGIN
    PRINT, 'ERROR Cortex ExtractStructure: Structure element not found'
    PRINT,'  INDEX  = ',index
    PRINT,'  STRUCT = '
    HELP,struct,/struct     
    STOP
  ENDIF

  RETURN, result
END
