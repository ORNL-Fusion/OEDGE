;--------------------------------------------------------------------------
; Function: NOT_STRING
; Date: 15.10.03
; Author: R.Martin
;--------------------------------------------------------------------------
; NOT_STRING
; 
; Function returns true if a variable is not of type STRING
;
; Calling Sequence
;
; result=NOT_STRING(var)
;
; var      : any IDL-variable type
; result   : FALSE(0B) if var is a STRING, TRUE(1B) otherwise
;
;--------------------------------------------------------------------------
;

function not_string, var

  return, size(var, /type) ne 7

end
  
