;--------------------------------------------------------------------------
; Function: NOT_STRUCTURE
; Date: 15.10.03
; Author: R.Martin
;--------------------------------------------------------------------------
; NOT_STRUCTURE
; 
; Function returns true if a variable is not of type STRUCTURE
;
; Calling Sequence
;
; result=NOT_STRUCTURE(var)
;
; var      : any IDL-variable type
; result   : FALSE(0B) if var is a STRUCTURE, TRUE(1B) otherwise
;
;--------------------------------------------------------------------------
;

function not_structure, var

  return, size(var, /type) ne 8

end
  
