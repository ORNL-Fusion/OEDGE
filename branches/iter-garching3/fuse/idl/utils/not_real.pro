;--------------------------------------------------------------------------
; Function: NOT_REAL
; Date: 15.10.03
; Author: R.Martin
;--------------------------------------------------------------------------
; NOT_REAL
; 
; Function returns ture if a variable cannot be used as a real-number.
;
; Calling Sequence
;
; result=NOT_REAL(var)
;
; var      : any IDL-variable type
; result   : FLASE(0B) if var is any IDL-number other than COMPLEX
;            TRUE(1B) otherwise
;
; /ptr     : Optional keyword, will dereference a pointer and report on
;            the contents 
;
;--------------------------------------------------------------------------
;

function not_real, var, ptr=ptrkey

  mask=[1,0,0,0,0,0,1,1, 1,1,1,1,0,0,0,0]

  if keyword_set(ptrkey) then begin
    if (ptr_valid(var) eq 0) then return, 1B
    return, byte(mask(size(*var(0), /type)))
  endif else $
    return, byte(mask(size(var, /type)))

end
  
