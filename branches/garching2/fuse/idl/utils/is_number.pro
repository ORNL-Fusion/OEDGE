;--------------------------------------------------------------------------
; Function: IS_NUMBER
; Date: 02.02.04
; Author: R.Martin
;--------------------------------------------------------------------------
; IS_NUMBER
; 
; Function returns ture if a variable is any numerical type
;
; Calling Sequence
;
; result=IS_NUMBER(var)
;
; var      : any IDL-variable type
; result   : TRUE(1B) if var is any IDL-number other than COMPLEX
;            FALSE(0B) otherwise
;
;--------------------------------------------------------------------------
;

function is_number, var

  mask=[0,1,1,1,1,1,1,0, 0,1,0,0,1,1,1,1]
  return, byte(mask(size(var, /type)))

end
  
