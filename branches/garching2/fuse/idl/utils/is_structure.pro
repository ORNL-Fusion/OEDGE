;--------------------------------------------------------------------------
; Function: IS_STRUCTURE
; Date: 15.10.03
; Author: R.Martin
;--------------------------------------------------------------------------
; IS_STRUCTURE
; 
; Function returns true if a variable is not of type STRUCTURE
;
; Calling Sequence
;
; result=IS_STRUCTURE(var)
;
; var      : any IDL-variable type
; result   : TRUE if var is a STRUCTURE, FALSE(1B) otherwise
;
;--------------------------------------------------------------------------
;

function is_structure, var

  return, size(var, /type) eq 8

end
  
