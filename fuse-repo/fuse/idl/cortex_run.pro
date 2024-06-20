;
; ======================================================================
;
args = command_line_args()
;@cortex_make
RESTORE,FILENAME='cortex.sav'
cortex_main, args
exit
;
; ======================================================================
;
