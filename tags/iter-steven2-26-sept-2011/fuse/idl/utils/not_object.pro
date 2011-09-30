;--------------------------------------------------------------------------
; Function: NOT_OBJECT
; Date: 15.10.03
; Author: R.Martin
;--------------------------------------------------------------------------
; NOT_OBJECT
; 
; Function returns true if a variable is not and IDL-Object
;
; Calling Sequence
;
; result=NOT_OBJECT(var)
;
; var      : any IDL-variable type
; result   : FALSE(0B) if var is an IDL-Object, TRUE(1B) otherwise
;
;--------------------------------------------------------------------------
;

function not_object, var

  return, size(var, /type) ne 11

end
  
