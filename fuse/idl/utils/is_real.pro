;--------------------------------------------------------------------------
; Function: IS_REAL
; Date: 15.10.03
; Author: R.Martin
;--------------------------------------------------------------------------
; IS_REAL
; 
; Function returns ture if a variable is any non-Complex number type
;
; Calling Sequence
;
; result=IS_REAL(var)
;
; var      : any IDL-variable type
; result   : TRUE(1B) if var is any IDL-number other than COMPLEX
;            FALSE(0B) otherwise
;
;--------------------------------------------------------------------------
;

function is_real, var

  mask=[0,1,1,1,1,1,0,0, 0,0,0,0,1,1,1,1]
  return, byte(mask(size(var, /type)))

end
  
