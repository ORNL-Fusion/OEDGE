;---------------------------------------------------------------------------
; Function: IS_STRING
; Concept Author: Julich-People
; Date: 15.10.03
;---------------------------------------------------------------------------
; IS_STRING
;
; Function returns TRUE if parameter is a string
;
; Format:
;   flag=IS_STRING(var)
;
; flag - BYTE*1, TRUE(1) if var is a string 
;                FALSE(0) otherwise
;
;---------------------------------------------------------------------------
;

function is_string, var

  return, size(var, /type) eq 7

end
