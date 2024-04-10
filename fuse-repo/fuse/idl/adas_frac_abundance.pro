PRO adas_FractionalAbundance

  color = [ 'Red', 'Blue' ]

  uid  = 'adas'
  elem = [ 'he', 'fe' , 'w' ]
  iz0  = [  2  ,  26  ,  74 ]
  year = [  96 ,  96  ,  88 ]

  te=adas_vector(low=1, high=20000, num=100)
  dens=fltarr(100)+1e12

  i = 0

print,elem[i]
print,year[i]
print,uid

  IF (elem[i] EQ 'w') THEN BEGIN
    run_adas405, uid=uid, year=year[i], elem=elem[i], te=te, dens=dens, $
                 frac=frac, files={acd:'/home/ITER/omullam/adas/adf11/acd88/acd88_w.dat',  $
                                   scd:'/home/ITER/omullam/adas/adf11/scd88/scd88_w.dat',  $
                                   plt:'/home/ITER/omullam/adas/adf11/plt88/plt88_w.dat',  $
                                   prb:'/home/ITER/omullam/adas/adf11/prb88/prb88_w.dat'}
  ENDIF ELSE  $
    run_adas405, uid=uid, year=year[i], elem=elem[i], te=te, dens=dens, frac=frac

  PLOT, [1,20000], [0.01,1.1], /nodata, XSTYLE=1, ystyle=1, CHARSIZE=2.0, /XLOG

  FOR j = 0, iz0[i] DO OPLOT, te, frac.ion[*,j], color=Truecolor(color[j MOD 2])

  XYOUTS, 0.5, 0.9, elem[i], /NORMAL, CHARSIZE=2.0, color=Truecolor('White')

END