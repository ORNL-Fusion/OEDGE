;--------------------------------------------------------------------------
; Function: UNDEFINED
; Date: 15.10.03
; Author: R.Martin
;--------------------------------------------------------------------------
; UNDEFINED
; 
; Returns TRUE(1B) if a variable is undefined

; Calling Sequence
;
; result=UNDEFINED(var)
;
; var      : any IDL-variable type
; result   : Returns TRUE(1B) if var exists
;            FALSE(1) otherwise
;
;--------------------------------------------------------------------------
;

function undefined, var

  return, (n_elements(var) eq 0)

end
  
