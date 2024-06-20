;--------------------------------------------------------------------------
; Function: IS_POINTER
; Date: 15.10.03
; Author: R.Martin
;--------------------------------------------------------------------------
; NOT_POINTER
; 
; Function returns ture if a variable is not a pointer.
;
; Calling Sequence
;
; result=IS_POINTER(var)
;
; var      : any IDL-variable type
; result   : TRUE if var is a valid IDL pointer,
;            FALSE(1B) otherwise
;
;--------------------------------------------------------------------------
;

function is_pointer, var

  return, size(var, /type) eq 10

end
  
