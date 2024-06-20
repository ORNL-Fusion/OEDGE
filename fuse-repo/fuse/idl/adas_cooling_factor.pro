;
; from M. O'Mullane, via email - SL, 22/04/2010
;

PRO adas_CoolingFactor

  color = [ 'White', 'Red', 'Green', 'Blue', 'Magenta', 'Cyan', 'Lightgreen' ]

  uid   = 'adas'
  iz0   = [2   , 4   , 6  , 7  ,  14  , 26  , 74 ]
  elem  = ['he', 'be', 'c', 'n',  'ar', 'fe', 'w']
  year  = [89  , 89  , 89 , 89 ,  89  , 89  , 97 ]

;  f_plt = [ '/work/projects/adas/adas/adf11/plt89/plt89_c.dat'      ,  $
;            '/work/projects/adas/adas/adf11/plt89/plt89_ar.dat'     ,  $
;            '/work/projects/adas/adas/adf11/plt89/plt89_fe.dat'     ,  $
;            '/home/ITER/lisgos/divimp/adas/adf11/plt97/plt97_w.dat' ]
;  f_prb = [ '/work/projects/adas/adas/adf11/prb89/prb89_c.dat'      ,  $
;            '/work/projects/adas/adas/adf11/prb89/prb89_ar.dat'     ,  $
;            '/work/projects/adas/adas/adf11/prb89/prb89_fe.dat'     ,  $
;            '/home/ITER/lisgos/divimp/adas/adf11/prb97/prb97_w.dat' ]
  ;f_plt = '/home/adas/adas/adf11/plt96/plt96_c.dat'
  ;f_prb = '/home/adas/adas/adf11/prb96/prb96_c.dat'


  itval = 200
  te    = adas_vector(low=1, high=1e5, num=itval)
  dens  = fltarr(itval) + 1.0e13


  PLOT, [8.0, 40000.0], [5.0E-37, 1.0E-30], /NODATA, /XLOG, /YLOG, XSTYLE=1, YSTYLE=1, CHARSIZE=1.5

  pow = fltarr(itval)

  FOR i = 0, N_ELEMENTS(iz0)-1 DO BEGIN

;   The easy way

    IF (elem[i] EQ 'w') THEN BEGIN
      ; w_path = '/home/ITER/omullam/adas/adf11/'
      ; w_year = '88'
       w_path = '/home/ITER/lisgos/divimp/adas/adf11/'
       w_year = '08'
      run_adas405, uid=uid, year=year[i], elem=elem[i], te=te, dens=dens, $
                   frac=frac, power=power, files={acd:w_path+'acd'+w_year+'/acd'+w_year+'_w.dat',  $
                                                  scd:w_path+'scd'+w_year+'/scd'+w_year+'_w.dat',  $
                                                  plt:w_path+'plt'+w_year+'/plt'+w_year+'_w.dat',  $
                                                  prb:w_path+'prb'+w_year+'/prb'+w_year+'_w.dat'}
    ENDIF ELSE  $
      run_adas405, uid=uid, year=year[i], elem=elem[i], te=te, dens=dens, $
                   frac=frac, power=power

print, i,elem[i],' ',color[i]
    OPLOT, te, power.total*1.0E-6, color=Truecolor(color[i])

    XYOUTS, 0.25 + FLOAT(i) * 0.1, 0.9, elem[i], /NORMAL, CHARSIZE=2.0, color=Truecolor(color[i])

;    Or build it up from plt and prb components - we need the equilibrium
;    fractional abundance - frac.ion - from above.
;    FOR j = 0, iz0[i]-1 DO BEGIN    ; neutral to H-like
;  
;      iz1 = j + 1
;      read_adf11, file=f_plt[i], class='plt', te=te, dens=dens, iz0=iz0[i], iz1=iz1, data=plt
;      read_adf11, file=f_prb[i], class='prb', te=te, dens=dens, iz0=iz0[i], iz1=iz1, data=prb     
;      pow = pow + frac.ion[*,j] * plt + frac.ion[*,j+1] * prb  ; nb: j+1 for prb
;    ENDFOR
;    print,min(pow)*1E-6,max(pow)*1E-6
;    OPLOT, te, pow*1E-6, color=Truecolor(color[i])
   
  ENDFOR

END
