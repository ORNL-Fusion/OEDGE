;--------------------------------------------------------------------------
; Function: NOT_POINTER
; Date: 15.10.03
; Author: R.Martin
;--------------------------------------------------------------------------
; NOT_POINTER
; 
; Function returns ture if a variable is not a pointer.
;
; Calling Sequence
;
; result=NOT_POINTER(var)
;
; var      : any IDL-variable type
; result   : FALSE if var is an IDL pointer,
;            TRUE(1B) otherwise
;
;--------------------------------------------------------------------------
;

function not_pointer, var

  return, size(var, /type) ne 10

end
  
