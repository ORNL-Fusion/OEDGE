;
; ======================================================================
;
FUNCTION cortex_RadiationLoss, dens, te, A, impurity_fraction

  CASE A OF
     1: elem = 'h'
     2: elem = 'he'
     4: elem = 'be'
     6: elem = 'c'
    18: elem = 'ar'
    26: elem = 'fe'
    ELSE: BEGIN
      PRINT,'ERROR cortex_AverageChargeState: Unrecognised element'
      PRINT,'  A =',A
      STOP
      END
  ENDCASE

  colors = ['Black','Red','Green','Blue','Orange','Purple', 'Hotpink', 'Darkseagreen', 'Silver']

; Mnemonic Class Data prefix
; ACD Coll.-diel. recom. coefft. R & U (or none)
; SCD Coll.-diel. ionis. coefft. R & U (or none)
; CCD Coll.-rad. charge exch. coefft. R & U (or none)
; PRB Coll.-diel. recom./brems. power coefft. R & U (or none)          ***
; PRC Coll.-rad. charge exch. recom. power coefft. R & U (or none)
; QCD Coll.-rad. metastable cross coupling coefft. R
; XCD Coll.-diel. parent meta. cross-coupling coefft. R
; PLT Coll.-rad. excit. line power coefft. R & U (or none)             ***
; PLS Coll.-rad. specific line excit. power coefft. R & U (or none)

PRINT,' ADAS CALL COMMENTED OUT! '
stop
; te=adas_vector(low=10, high=1000, num=100) 
 dens=MAKE_ARRAY(100,/FLOAT,VALUE=1.0E+18)

  run_adas405, uid='adas', elem=elem, year=96, te=te, dens=dens*1.0E-06, frac=frac

  FOR j = 1, A DO BEGIN
print,'j=',j
    read_adf11,uid='adas',year=93,iz0=A,iz1=j,class='plt',te=te,dens=dens*1.0E-06,data=data

    IF (j EQ 1) THEN PLOT, te, data*1.0E-6 , /XLOG, /YLOG , yrange=[1.0E-37, 1.0E-30] , title = 'beryllium adf11 plt  +1-black, +2-red, +3-green, +4-blue' , ytitle = 'W m3', xtitle = 'te (eV)'  $
                ELSE OPLOT, te, data*1.0E-6, color=Truecolor(colors[j-1])
  ENDFOR



return,0

END
;
; ======================================================================
;
FUNCTION cortex_AverageChargeState, dens, te, A

  CASE A OF
     1: elem = 'h'
     2: elem = 'he'
     4: elem = 'be'
     6: elem = 'c'
    18: elem = 'ar'
    26: elem = 'fe'
    ELSE: BEGIN
      PRINT,'ERROR cortex_AverageChargeState: Unrecognised element'
      PRINT,'  A =',A
      STOP
      END
  ENDCASE

;  te=adas_vector(low=1, high=1000, num=100) 
;  dens=fltarr(100)+1e12
;help,dens 
  run_adas405, uid='adas', elem=elem, year=96, te=te, dens=dens*1.0E-06, frac=frac

  colors = ['Black','Red','Green','Blue','Orange','Purple', 'Hotpink', 'Darkseagreen', 'Silver']

;  plot_oo, [1,20000], [0.01,1.1], /nodata, ystyle=1
;  for j = 0, A do oplot, te, frac.ion[*,j] ; , color=TRUECOLOR(colors[j])

  avg_state = MAKE_ARRAY(N_ELEMENTS(te),/FLOAT,VALUE=0.0)
  FOR j = 0, A DO BEGIN
    avg_state = avg_state + frac.ion[*,j] * FLOAT(j) 
  ENDFOR

  RETURN, avg_state

END
;
; ======================================================================
;
FUNCTION cortex_FusionRate, temperature

;  T = temperature 
  T = DOUBLE(temperature / 1.0E+3)



  C = [6.6610D+0, 643.41D-16, 15.13600D+3, 75.1890D+3 ,  $
       4.6064D+3, 13.500D+3 , -0.10675D+3, 0.01366D+3]

  zeta = 1.0D0 - (        C[2]*T + C[4]*T^2 + C[6]*T^3) /  $
                 (1.0D0 + C[3]*T + C[5]*T^2 + C[7]*T^3)

  xi = C[0] / (T^(1.0D0/3.0D0))

;print,T
;print,C

  result = C[1] * zeta^(-5.0D0/6.0D0) * xi^2.0D0 *  $
           EXP(-3.0D0 * zeta^(1.0D0/3.0D0) * xi)

  result = result * 1.0D-6  ; convert from cm^3 to m^3

  i = WHERE(T LT 0.2D0 OR T GT 100.0D0, n)
  IF (n GT 0) THEN result[i] = 0.0D0

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_PlotPedestalModel, plot, data_array, ps=ps

  PRINT
  PRINT,'----------------------- NEW PLOT -----------------------'
  PRINT

  MAXNYDATA = 7
  MAXNXDATA = 1000

  focus = plot.focus

  !P.BACKGROUND = TrueColor('White')


  dev_xsize = !D.X_SIZE
  dev_ysize = !D.Y_SIZE

  notes    = plot.notes
  charsize = plot.charsize
  charsize_labels = charsize

  xy_label = [0.96,0.04,0.12,0.88]

  IF (focus) THEN BEGIN
    xy_label = [0.93,0.07,0.13,0.87]
    charsize_labels = charsize_labels * 1.2
  ENDIF

  option = plot.option
  IF (option EQ 100 OR option EQ 101) THEN option = 999

  plot_peak = plot.peak

  CASE option OF
;   --------------------------------------------------------------------
    1: BEGIN
       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot = 6
       plot_xn = 2
       plot_yn = 3
       title = plot.title 
       subtitle = ['1','2','3','4','5','6']
       xtitle   = 'r (m)'
       ytitle   = ['1','2','3','4','5','6']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace = [1,1,1,1,2,1]
 
       impurity_species  = [2    ,4   ,18    ,26    ]
       impurity_fraction = [0.01 ,0.01,0.01  ,0.01  ]

       ntrace[1] = N_ELEMENTS(impurity_species)
       ntrace[2] = N_ELEMENTS(impurity_species) + 2
       ntrace[5] = N_ELEMENTS(impurity_species) 

       END
;   --------------------------------------------------------------------
    ELSE: BEGIN  
      PRINT, 'ERROR cortex_PlotPedestalModel: Unrecognised plot option'
      PRINT, '  OPTION = ',option,' (',plot.option,')'
      RETURN, -1
      END
  ENDCASE

  colors = ['Black','Red','Green','Blue','Orange','Purple', 'Hotpink', 'Darkseagreen', 'Silver']
;
; Setup plot:
; ----------------------------------------------------------------------
;
;
; Setup plot area:
; ----------------------------------------------------------------------

  size = TOTAL(plot.frame_bnds[1] - plot.frame_bnds[0])

  IF (focus) THEN BEGIN
    plot_xn = 1
    plot_yn = 1
  ENDIF
  xsize = (size - 2.0 * plot_xboarder - FLOAT(plot_xn-1) * plot_xspacing) / FLOAT(plot_xn)
  ysize = (1.0  - 2.0 * plot_yboarder - FLOAT(plot_yn-1) * plot_yspacing) / FLOAT(plot_yn) 
  xi = 1
  yi = 0

  FOR iplot = 1, nplot DO BEGIN
    IF (focus NE 0 AND focus NE iplot) THEN CONTINUE

    ndata = N_ELEMENTS(TAG_NAMES(data_array))
    IF (ndata LE 0) THEN BEGIN
      PRINT, 'ERROR cortex_PlotPedestalModel: No data found'
      RETURN, -1
    ENDIF

    yi = yi + 1
    IF (yi EQ plot_yn+1) THEN BEGIN
      xi = xi + 1
      yi = 1
    ENDIF
    
    xcen  =       xsize * (0.5 + FLOAT(xi-1)) + 1.5 * plot_xboarder + plot_xspacing * FLOAT(xi-1) + plot.frame_bnds[0]
    ycen  = 1.0 - ysize * (0.5 + FLOAT(yi-1)) - 1.0 * plot_yboarder - plot_yspacing * FLOAT(yi-1)
    xpos  = [xcen - 0.5 * xsize, xcen + 0.5 * xsize]
    ypos  = [ycen - 0.5 * ysize, ycen + 0.5 * ysize]

    xmin =  1.0E+35
    xmax = -1.0E+35
    ymin =  1.0E+35
    ymax = -1.0E+35
    FOR idata = 1, ndata DO BEGIN
      val = cortex_ExtractStructure(data_array,idata)
      CASE option OF
;       ----------------------------------------------------------------
        1: BEGIN
          file = val.pedestal.file
          integral = ' '
          str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
          str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
          labels[0] = labels[0] + STRING(idata-1) + '/' + str[0] + integral + ' :'

          xdata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)
          ydata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)

          i = WHERE(val.pedestal.r LE val.pedestal.a,n)

          xdata[0:n-1,0] = val.pedestal.r[i] / val.pedestal.a

          stored_energy = 1.602E-19 * 1.5 * TOTAL(val.pedestal.volume * val.pedestal.dens *  $
                                                  (val.pedestal.te + val.pedestal.ti))
     

          PRINT,'STORED_ENERGY (MJ)= ',stored_energy * 1.0E-6
          PRINT,'VOLUME=             ',TOTAL(val.pedestal.volume)
          PRINT,'VOLUME=             ',TOTAL(val.pedestal.new_volume)





          CASE iplot OF
            1: BEGIN
              xdata[0:n-1,0] = val.pedestal.new_r / val.pedestal.a
              ydata[0:n-1,0] = val.pedestal.new_volume
              END
            2: BEGIN
              n = N_ELEMENTS(val.pedestal.new_r)
              impurity_charge = MAKE_ARRAY(n,N_ELEMENTS(impurity_species),/FLOAT,VALUE=0.0)
              FOR j = 0, N_ELEMENTS(impurity_species)-1 DO BEGIN
                xdata[0:n-1,j] = val.pedestal.new_r / val.pedestal.a
                ydata[0:n-1,j] = cortex_AverageChargeState(val.pedestal.new_dens,val.pedestal.new_te,impurity_species[j])
                i = WHERE(ydata[0:n-1,j] LT 0.0, i_count)
                IF (i_count GT 0) THEN ydata[i,j] = 0.0
                impurity_charge[*,j] = ydata[0:n-1,j]
              ENDFOR
              END
            3: BEGIN
              xdata[0:n-1,0] = val.pedestal.new_r / val.pedestal.a
              ydata[0:n-1,0] = val.pedestal.new_dens
              impurity_dilution = MAKE_ARRAY(n,/FLOAT,VALUE=0.0)
              FOR j = 0, N_ELEMENTS(impurity_species)-1 DO BEGIN
                ydata[0:n-1,j+1] = val.pedestal.new_dens * impurity_charge[*,j] * impurity_fraction[j]
                impurity_dilution = impurity_dilution + impurity_charge[*,j] * impurity_fraction[j]
              ENDFOR
              fuel_density = val.pedestal.new_dens * (1.0 - impurity_dilution)
              ydata[0:n-1,ntrace[iplot-1]-1] = fuel_density
              END
            4: BEGIN
              xdata[0:n-1,0] = val.pedestal.new_r / val.pedestal.a
              ydata[0:n-1,0] = val.pedestal.new_te
              END
            5: BEGIN
              n = N_ELEMENTS(val.pedestal.new_r)
              xdata[0:n-1,0] = val.pedestal.new_r / val.pedestal.a
              ydata[0:n-1,0] = ((0.5*fuel_density)^2) * cortex_FusionRate(val.pedestal.new_te) 
              ydata[0:n-1,1] = ydata[0:n-1,0] * val.pedestal.new_volume
              PRINT,'FUSION POWER= ',TOTAL(ydata[0:n-1,1]) / 1.0E+6 * 17.6E+6 * 1.602E-19
              END
            6: BEGIN
              xdata[0:n-1,0] = val.pedestal.new_r / val.pedestal.a
              FOR j = 0, N_ELEMENTS(impurity_species)-1 DO BEGIN                                
                a= cortex_RadiationLoss(val.pedestal.new_dens,val.pedestal.new_te,4,0.01)
return,0
                ydata[0:n-1,j] = val.pedestal.new_dens * impurity_fraction[j]
              ENDFOR
              END
          ENDCASE
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotPedestalModel: Unrecognised plot option'
          PRINT, '  OPTION = ',plot.option
          RETURN, -1
          END
      ENDCASE

;     Package up the data for plotting:
      name = 'data' + STRING(idata,FORMAT='(I0)')
      data = { n : ntrace[iplot-1], x : xdata, y : ydata, file : file } 
      IF (idata EQ 1) THEN data_store = CREATE_STRUCT(           name,data) ELSE  $
                           data_store = CREATE_STRUCT(data_store,name,data)
      IF (N_ELEMENTS(WHERE(plot.xrange EQ 0.0)) NE 2) THEN BEGIN
        i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1])        
        IF (N_ELEMENTS(i) EQ 1) THEN BEGIN
          PRINT,'ERROR cortex_PlotPedestalModel: No data within XRANGE'
          PRINT,'  IPLOT          = ',iplot
          PRINT,'  IDATA          = ',idata
          PRINT,'  XRANGE         = ',plot.xrange
          PRINT,'  XDATA MIN,MAX  = ',MIN(xdata),MAX(xdata)
          RETURN, -1
        ENDIF
      ENDIF ELSE i = WHERE(xdata NE -999.0)
      xmin = MIN([xmin,xdata[i]])
      xmax = MAX([xmax,xdata[i]])
;     Reform YDATA for finding the maximum value:    *** ALMOST CERTAINLY A BETTER WAY TO DO THIS ***
      i = WHERE(ydata[*,0] NE -999.0)
      ydata1 = REFORM(ydata[i,0])
      FOR j = 1, ntrace[iplot-1]-1 DO BEGIN
        i = WHERE(ydata[*,j] NE -999.0)
        ydata1 = [ydata1,REFORM(ydata[i,j])] 
      ENDFOR
      ymin = MIN([ymin,ydata1])
      ymax = MAX([ymax,ydata1])
    ENDFOR
    IF (ymin GT 0.0 AND ymax GT 0.0) THEN ymin = 0.0  ; Makes things a bit clearer on the plots I think...
    IF (ymin LT 0.0 AND ymax LT 0.0) THEN ymax = 0.0
    deltay = ymax - ymin
    ymin = ymin - 0.05 * deltay
    ymax = ymax + 0.05 * deltay

    PRINT,'YMIN,YMAX=',ymin,ymax,ntrace[iplot-1]

;   Axes:
    xrange = [xmin,xmax]
    yrange = [ymin,ymax]
    position = [xpos[0],ypos[0],xpos[1],ypos[1]]

    plot_type = default_plot_type                                           
    IF (focus NE 0 OR yi EQ plot_yn OR iplot EQ nplot) THEN plot_type = 2   ; Show x-axis label

    CASE (plot_type) OF
      1: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTICKFORMAT='(A1)',/NOERASE                                  
      2: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTITLE=xtitle,/NOERASE                                 
    ENDCASE

;   Write sub-title for each plot:
    XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
            (xy_label[2] * ypos[0] + xy_label[3] * ypos[1]) * dev_ysize,  $
            subtitle[iplot-1], CHARSIZE=charsize_labels, /DEVICE

;   Data:
    FOR idata = 1, ndata DO BEGIN
      val = cortex_ExtractStructure(data_store,idata)
      CASE option OF
;       ----------------------------------------------------------------
        1: BEGIN
          cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                          dev_xsize, dev_ysize, charsize_labels, colors

          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 

          i = WHERE(val.x[*,0] NE -999.0)

;              print,'data: ',val.x[i,0], val.y[i,0]
          OPLOT, val.x[i,0], val.y[i,0], COLOR=TrueColor(colors[idata-1]); , psym=6

          CASE iplot OF
            1: 
            2: BEGIN
              FOR j = 1, N_ELEMENTS(impurity_species)-1 DO BEGIN
;               print,'data: ',val.x[i,0], val.y[i,j]
                OPLOT, val.x[i,0], val.y[i,j], COLOR=TrueColor(colors[j])
              ENDFOR
              END
            3: BEGIN
              FOR j = 1, N_ELEMENTS(impurity_species) DO BEGIN
                OPLOT, val.x[i,0], val.y[i,j], COLOR=TrueColor(colors[j])
              ENDFOR
              OPLOT, val.x[i,0], val.y[i,ntrace[iplot-1]-1], COLOR=TrueColor('Black'), LINESTYLE=2
              END
            4:
            5: BEGIN
              i = WHERE(val.x[*,0] NE -999.0)
              OPLOT, val.x[i,0], val.y[i,1], COLOR=TrueColor(colors[idata-1])
              END
            6: BEGIN
              FOR j = 1, N_ELEMENTS(impurity_species)-1 DO  $
                OPLOT, val.x[i,0], val.y[i,j], COLOR=TrueColor(colors[j])
              END
          ENDCASE
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotPedestalModel: Unrecognised plot option'
          PRINT, '  OPTION = ',plot.option,' (',option,')'
          RETURN, -1
          END
      ENDCASE

    ENDFOR  ; idata loop

  ENDFOR  ; iplots loop
;
; Put main title on the page
; ----------------------------------------------------------------------
  cortex_PageTitle, plot, ps, 'none', dev_xsize, dev_ysize, title, notes, charsize

  RETURN, 0

END
;
; ======================================================================
;




;
; ======================================================================
;
FUNCTION cortex_PlotPedestalModel_OLD, plot, mid, ps=ps

  tags = TAG_NAMES(mid)

  ndata = N_ELEMENTS(TAG_NAMES(mid))
  IF (ndata LE 0) THEN BEGIN
    PRINT, 'ERROR Cortex PlotMidplaneProfiles: No data found'
    RETURN, -1
  ENDIF
  PRINT,'NDATA= ',ndata

  focus = plot.focus

;  return

;  window_id = 0
;  window_xsize = 700
;  window_ysize = 900

  !P.BACKGROUND = TrueColor('White')

  plot_title  = 'PEDESTAL PLOT'
  plot_xtitle = 'psin'
  plot_ytitle = '??? (???)'

  colors = ['Black','Red','Green','Blue']

  plot_xboarder = 0.1
  plot_yboarder = 0.1
;
; Setup plot:
; ----------------------------------------------------------------------
;
  IF (focus NE 0) THEN BEGIN
    plot_xi = [1,1]
    plot_yi = [1,1]
    plot_xn = 1
    plot_yn = 1
  ENDIF ELSE BEGIN
    plot_xi = [1,2]
    plot_yi = [1,1]
    plot_xn = 2
    plot_yn = 1
  ENDELSE

  xsize = (1.0 - 2.0 * plot_xboarder) / FLOAT(plot_xn)
  ysize = (1.0 - 2.0 * plot_yboarder) / FLOAT(plot_yn)

  FOR iplot = 1, 2 DO BEGIN

    IF (focus NE 0 AND focus NE iplot) THEN CONTINUE

    xi = plot_xi[iplot-1]
    yi = plot_yi[iplot-1]

    PRINT,'XI,YI=',xi,yi

    xcen  =       xsize * (0.5 + FLOAT(xi - 1)) + plot_xboarder
    ycen  = 1.0 - ysize * (0.5 + FLOAT(yi - 1)) - plot_yboarder
    xpos  = [xcen - 0.5  * xsize, xcen + 0.35 * xsize]
    ypos  = [ycen - 0.45 * ysize, ycen + 0.5  * ysize]

    xmin =  1.0E+35
    xmax = -1.0E+35
    ymin =  1.0E+35
    ymax = -1.0E+35
    FOR idata = 1, ndata DO BEGIN
      val = cortex_ExtractStructure(mid,idata)

      xdata = val.r
      IF (iplot EQ 1) THEN ydata = val.dens
      IF (iplot EQ 2) THEN ydata = [val.te,val.ti]

      xmin = MIN([xmin,xdata])
      xmax = MAX([xmax,xdata])
      ymin = MIN([ymin,ydata])
      ymax = MAX([ymax,ydata])
    ENDFOR
    PRINT, 'XMIN,MAX=',xmin,xmax
    PRINT, 'YMIN,MAX=',ymin,ymax

;    IF (iplot EQ 2) THEN ymax = 500.0

;   Axes:
    xrange = [xmin,xmax]
    yrange = [ymin,ymax]
    position = [xpos[0],ypos[0],xpos[1],ypos[1]]

    CASE (iplot) OF
      focus: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position, YTITLE=plot_ytitle
      1    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position, YTITLE=plot_ytitle
      2    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position, YTITLE=plot_ytitle,/NOERASE
    ENDCASE

    XYOUTS, 0.95*xmin+0.05*xmax, 0.1*ymin+0.9*ymax, plot_title
;   Data:
    FOR idata = 1, ndata DO BEGIN
      val = cortex_ExtractStructure(mid,idata)
      xdata = val.r
      IF (iplot EQ 1) THEN ydata = val.dens
      IF (iplot EQ 2) THEN ydata = val.te  

      OPLOT, [val.a,val.a], [ymin,ymax], LINESTYLE=3 
      IF (iplot EQ 1) THEN OPLOT, [val.cross_ne,val.cross_ne], [ymin,ymax], LINESTYLE=3, COLOR=TrueColor('red')
      IF (iplot EQ 2) THEN OPLOT, [val.cross_te,val.cross_te], [ymin,ymax], LINESTYLE=3, COLOR=TrueColor('red')

      OPLOT, xdata, ydata, COLOR=TrueColor(colors[idata-1])
 
      IF (iplot EQ 2) THEN BEGIN
        ydata = val.ti
        OPLOT, xdata, ydata, COLOR=TrueColor(colors[idata-1]), LINESTYLE=2 
      ENDIF

      step = 0.05
      frac = step * FLOAT(idata) + 0.1
      PRINT,frac
      IF (iplot EQ 1) THEN  $
        XYOUTS, 0.95*xmin+0.05*xmax,frac*ymin+(1.0-frac)*ymax,  $
                val.file, COLOR=TrueColor(colors[idata-1])
    ENDFOR

  ENDFOR

  RETURN, 0

END
;
; ======================================================================
;



