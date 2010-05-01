PRO adas_go

te=adas_vector(low=1, high=20000, num=100)
dens=fltarr(100)+1e12
run_adas405, uid='adas', elem='fe', year=96, te=te, dens=dens, frac=frac
plot_oo, [1,20000], [0.01,1.1], /nodata, ystyle=1
for j = 0, 26 do oplot, te, frac.ion[*,j]

END