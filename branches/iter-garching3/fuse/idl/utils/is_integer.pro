;--------------------------------------------------------------------------
; Function: IS_INTEGER
; Date: 15.10.03
; Author: R.Martin
;--------------------------------------------------------------------------
; IS_INTEGER
; 
; Function returns ture if a variable can be used as an interger.
; In this case an integer is any of the IDL types BYTE, INT, UINT, LONG
; ULONG
;
; Calling Sequence
;
; result=IS_INTEGER(var)
;
; var      : any IDL-variable type
; result   : Returns TRUE(1B) if var is any INTEGER-type, 
;            FALSE(1) otherwise
;
;--------------------------------------------------------------------------
;

function is_integer, var

  mask=[0,1,1,1,0,0,0,0,   0,0,0,0,1,1,1,1]
  return, byte(mask(size(var, /type)))

end
  
