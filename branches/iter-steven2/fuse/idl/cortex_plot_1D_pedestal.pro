;
; ======================================================================
;
FUNCTION cortex_PlotPedestalModel, plot, mid, ps=ps

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



