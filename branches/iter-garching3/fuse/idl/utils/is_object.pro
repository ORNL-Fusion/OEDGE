;--------------------------------------------------------------------------
; Function: IS_OBJECT
; Date: 02.09.05
; Author: R.Martin
;--------------------------------------------------------------------------
; IS_OBJECT
; 
; Function returns true if a variable is of type OBJECT
;
; Calling Sequence
;
; result=IS_OBJECT(var)
;
; var      : any IDL-variable type
; result   : True(1B) if var is an IDL-Object, False(0B) otherwise
;
;--------------------------------------------------------------------------
;

function is_object, var

  return, size(var, /type) eq 11

end
  
