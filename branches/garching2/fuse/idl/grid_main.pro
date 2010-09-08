
;
; ======================================================================
;
PRO grid_Main, args

  IF (NOT KEYWORD_SET(args)) THEN args = 'none'

  PRINT, N_ELEMENTS(args),args

  CASE (args[0]) OF
    'suppliment': BEGIN
      getb,args[1],machine='EQU',equ=args[2]
      END
    ELSE: BEGIN
      PRINT,'grid_Main: Unknown command'
      PRINT,'  COMMAND = ',args[0]
      END
  ENDCASE
  EXIT
END


