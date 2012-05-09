;--------------------------------------------------------------------------
; Function: NOT_INTEGER
; Date: 15.10.03
; Author: R.Martin
;--------------------------------------------------------------------------
; NOT_INTEGER
; 
; Function returns ture if a variable is not of type of an integer type.
; In this case an integer is any of the IDL types BYTE, INT, UINT, LONG
; ULONG
;
; Calling Sequence
;
; result=NOT_INTEGER(var)
;
; var      : any IDL-variable type
; result   : Returns FALSE(0) if var is any INTEGER-type, 
;            TRUE(1) otherwise
;
;--------------------------------------------------------------------------
;

function not_integer, var

  mask=[1,0,0,0,1,1,1,1,   1,1,1,1,0,0,0,0]
  return, byte(mask(size(var, /type)))

end
  
