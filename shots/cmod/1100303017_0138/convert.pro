;
; ======================================================================
;
PRO ConvertThomson

  filename = 'thomson'

  RESTORE, filename+'.sav'

  fp = 3
  FREE_LUN, fp
  OPENW, fp, filename+'.dat', error=err
  IF (err NE 0) THEN BEGIN
    PRINT,'ERROR: Unable to open output file
  ENDIF

  PRINTF,fp,'* IDL Conversion from .sav to ASCII file'
  PRINTF,fp,'*'
  FOR i = 0, N_ELEMENTS(comment)-1 DO PRINTF,fp,'* '+comment[i] ; 
  PRINTF,fp,'*'
  PRINTF,fp,'*        R_e','       psi_e','         n_e','         T_e',  $
         FORMAT='(4A12)'
  FOR i = 1, N_ELEMENTS(R_e)-1 DO  $
     PRINTF,fp, R_e[i], psi_e[i], n_e[i], T_e[i],  $
            FORMAT='(2F12.6,E12.4,F12.4)'

  FREE_LUN, fp

END
