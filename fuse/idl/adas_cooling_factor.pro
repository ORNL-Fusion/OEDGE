;
; from M. O'Mullane, via email - SL, 22/04/2010
;

PRO cooling_curve

; Carbon in a fixed density and varying temperature plasma

iz0   = 6; 14 ; 26 ; 74 ; 6
uid   = 'adas'
elem  = 'c' ; 'ar'; 'fe'; 'w'; 'c'
year  = 89 ; 89 ; 89 ; 97 ; 96
;f_plt = '/work/projects/adas/adas/adf11/plt89/plt89_ar.dat'
;f_prb = '/work/projects/adas/adas/adf11/prb89/prb89_ar.dat'
;f_plt = '/work/projects/adas/adas/adf11/plt89/plt89_fe.dat'
;f_prb = '/work/projects/adas/adas/adf11/prb89/prb89_fe.dat'
;f_plt = '/home/ITER/lisgos/divimp/adas/adf11/plt97/plt97_w.dat'
;f_prb = '/home/ITER/lisgos/divimp/adas/adf11/prb97/prb97_w.dat'
f_plt = '/work/projects/adas/adas/adf11/plt89/plt89_c.dat'
f_prb = '/work/projects/adas/adas/adf11/prb89/prb89_c.dat'
;f_plt = '/home/adas/adas/adf11/plt96/plt96_c.dat'
;f_prb = '/home/adas/adas/adf11/prb96/prb96_c.dat'

itval = 200
te    = adas_vector(low=1, high=1e5, num=itval)
dens  = fltarr(itval) + 1.0e13

; Or build it up from plt and prb components - we need the equilibrium
; fractional abundance - frac.ion - from above.



; The easy way

run_adas405, uid=uid, year=year, elem=elem, te=te, dens=dens, $
             frac=frac, power=power

plot_oo, te, power.total*1.0E-6, charsize=2

;help,power,/struct
;print,power.stage

;stop


pow = fltarr(itval)

for j = 0, iz0-1 do begin    ; neutral to H-like
  
  iz1 = j + 1
  
  read_adf11, file=f_plt, class='plt', te=te, dens=dens, iz0=iz0, iz1=iz1, data=plt
  read_adf11, file=f_prb, class='prb', te=te, dens=dens, iz0=iz0, iz1=iz1, data=prb     
  
  pow = pow + frac.ion[*,j] * plt + frac.ion[*,j+1] * prb  ; nb: j+1 for prb
    
endfor

oplot, te, pow*1E-6, psym=2

stop




END
