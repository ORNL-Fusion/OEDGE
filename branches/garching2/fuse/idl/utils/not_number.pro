;--------------------------------------------------------------------------
; Function: NOT_NUMBER
; Date: 02.02.04
; Author: R.Martin
;--------------------------------------------------------------------------
; NOT_NUMBER
; 
; Function returns ture if a variable is not a numerical data-type
;
; Calling Sequence
;
; result=NOT_NUMBER(var)
;
; var      : any IDL-variable type
; result   : FLASE(0B) if var is any IDL-numerical data-type
;            TRUE(1B) otherwise
;
;--------------------------------------------------------------------------
;

function not_number, var

  mask=[1,0,0,0,0,0,0,1, 1,0,1,1,0,0,0,0]
  return, byte(mask(size(var, /type)))

end
  
