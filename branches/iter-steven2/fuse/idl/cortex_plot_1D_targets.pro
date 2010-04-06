;
; ======================================================================
;
FUNCTION cortex_PlotTargetProfiles, target, ps=ps

;.resolve_all
;@script
;help,/source_files

  tags = TAG_NAMES(target)
  print,tags

  help,target,/struct
;  help,target.data1,/struct
;  help,target.data2,/struct

  ndata = N_ELEMENTS(TAG_NAMES(target))
  IF (ndata LE 0) THEN BEGIN
    PRINT, 'ERROR Cortex PlotTargetProfiles: No data found'
    RETURN, -1
  ENDIF
  PRINT,'NDATA= ',ndata


  window_id = 0
  window_xsize = 700
  window_ysize = 700

  !P.BACKGROUND = TrueColor('White')

  plot_xn = 2
  plot_yn = 3
  plot_xboarder = 0.1
  plot_yboarder = 0.1


  plot_title  = 'TARGET PLOT'
  plot_xtitle = 'psin (m)'
  plot_ytitle = '??? (???)'

  colors = ['Black','Red','Green','Blue']

;
; Setup plot:
; ----------------------------------------------------------------------
;
  plot_xi = [1 , 1, 1, 2, 2, 2]
  plot_yi = [1 , 2, 3, 1, 2, 3]

  xsize = (1.0 - 2.0 * plot_xboarder) / FLOAT(plot_xn)
  ysize = (1.0 - 2.0 * plot_yboarder) / FLOAT(plot_yn)

  FOR iplot = 1, 6 DO BEGIN

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
      val = cortex_ExtractStructure(target,idata)
      xdata = val.psin
      IF (iplot EQ 1) THEN ydata = val.jsat[*,0]
      IF (iplot EQ 2) THEN ydata = val.te  [*,0]
      IF (iplot EQ 3) THEN ydata = val.ti  [*,0]
      IF (iplot EQ 4) THEN ydata = val.jsat[*,1]
      IF (iplot EQ 5) THEN ydata = val.te  [*,1]
      IF (iplot EQ 6) THEN ydata = val.ti  [*,1]

      IF (iplot EQ 2) THEN ydata = [ydata,val.ti[*,0]]
      IF (iplot EQ 5) THEN ydata = [ydata,val.ti[*,1]]

      xmin = MIN([xmin,xdata])
      xmax = MAX([xmax,xdata])
      ymin = MIN([ymin,ydata])
      ymax = MAX([ymax,ydata])
    ENDFOR
;    xmin =  4.1
;    xmax =  4.2
;    ymin = -4.0
;    ymax = 0.5E+19
    PRINT, 'XMIN,MAX=',xmin,xmax
    PRINT, 'YMIN,MAX=',ymin,ymax

;   Axes:
    xrange = [xmin,xmax]
    yrange = [ymin,ymax]
    position = [xpos[0],ypos[0],xpos[1],ypos[1]]

    CASE (iplot) OF
      1: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=plot_ytitle,XTICKFORMAT='(A1)'                                    
      2: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=plot_ytitle,XTICKFORMAT='(A1)',/NOERASE
      3: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=plot_ytitle,XTITLE=plot_xtitle,/NOERASE
      4: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=plot_ytitle,XTICKFORMAT='(A1)',/NOERASE
      5: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=plot_ytitle,XTICKFORMAT='(A1)',/NOERASE
      6: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=plot_ytitle,XTITLE=plot_xtitle,/NOERASE                                 
    ENDCASE

    XYOUTS, 0.95*xmin+0.05*xmax, 0.1*ymin+0.9*ymax, plot_title
;   Data:
    FOR idata = 1, ndata DO BEGIN
      val = cortex_ExtractStructure(target,idata)
      xdata = val.psin
      IF (iplot EQ 1) THEN ydata = val.jsat[*,0]
      IF (iplot EQ 2) THEN ydata = val.te  [*,0]
      IF (iplot EQ 3) THEN ydata = val.ti  [*,0]
      IF (iplot EQ 4) THEN ydata = val.jsat[*,1]
      IF (iplot EQ 5) THEN ydata = val.te  [*,1]
      IF (iplot EQ 6) THEN ydata = val.ti  [*,1]

      OPLOT, xdata, ydata, COLOR=TrueColor(colors[idata-1])
 
      IF (iplot EQ 2) THEN BEGIN
        ydata = val.ti[*,0]
        OPLOT, xdata, ydata, COLOR=TrueColor(colors[idata-1]),LINESTYLE=2 
      ENDIF
      IF (iplot EQ 5) THEN BEGIN
        ydata = val.ti[*,1]
        OPLOT, xdata, ydata, COLOR=TrueColor(colors[idata-1]),LINESTYLE=2 
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


