PRO adas_go

elem = 'fe'
iz0  = 26

te=adas_vector(low=1, high=20000, num=100)
dens=fltarr(100)+1e12
run_adas405, uid='adas', elem=elem, year=96, te=te, dens=dens, frac=frac

plot, [1,20000], [0.01,1.1], /nodata, XSTYLE=1, ystyle=1, CHARSIZE=2.0, /XLOG

for j = 0, iz0 do oplot, te, frac.ion[*,j]

END