;
; ======================================================================
;
FUNCTION cortex_PlotCoreProfiles, plot, core, ps=ps

;.resolve_all
;@script
;help,/source_files

  vdata = core

  tags = TAG_NAMES(vdata)
;  tags = TAG_NAMES(mid)

  PRINT
  PRINT,'----------------------- NEW PLOT -----------------------'
  PRINT

  ndata = N_ELEMENTS(TAG_NAMES(vdata))
;  ndata = N_ELEMENTS(TAG_NAMES(mid))
  IF (ndata LE 0) THEN BEGIN
    PRINT, 'ERROR Cortex PlotCoreProfiles: No data found'
    RETURN, -1
  ENDIF
  PRINT,'NDATA= ',ndata

  focus = plot.focus

  !P.BACKGROUND = TrueColor('White')

  plot_title  = 'Core average radial plot'
  plot_subtitle = ['ELECTRON DENSITY','PLASMA PRESSURE'       ,  $
                   'PLASMA TEMPERATURES (Te-solid, Ti-dashed)',  $
                   'Zeff','IMPURITY CONCENTRATION','IMPURITY DENSITY'] 
  plot_xtitle = 'r/a' ; 'rho (m)'
  plot_ytitle = ['AVG ne (m-3)','AVG p (Pa)','AVG Te,i (eV)','Zeff','AVG C','AVG n_I+ (m-3)']

  colors = ['Black','Red','Green','Blue','Orange','Purple', 'Hotpink', 'Lime']

  plot_xboarder = 0.05
  plot_yboarder = 0.1

  plot_xspacing = 0.125
  plot_yspacing = 0.025
;
; Setup plot:
; ----------------------------------------------------------------------
;
  IF (focus NE 0) THEN BEGIN
    plot_xi = [1,1,1,1,1,1]
    plot_yi = [1,1,1,1,1,1]
    plot_xn = 1
    plot_yn = 1
  ENDIF ELSE BEGIN
    plot_xi = [1,1,1,2,2,2]
    plot_yi = [1,2,3,1,2,3]
    plot_xn = 2
    plot_yn = 3
  ENDELSE

  xsize = (1.0 - 2.0 * plot_xboarder - 1.0 * plot_xspacing) / FLOAT(plot_xn) 
  ysize = (1.0 - 2.0 * plot_yboarder - 2.0 * plot_yspacing) / FLOAT(plot_yn)

  FOR iplot = 1, 6 DO BEGIN
    IF (focus NE 0 AND focus NE iplot) THEN CONTINUE

    xi = plot_xi[iplot-1]
    yi = plot_yi[iplot-1]

    xcen  =       xsize * (0.5 + FLOAT(xi-1)) + 1.5 * plot_xboarder + plot_xspacing * FLOAT(xi-1)
    ycen  = 1.0 - ysize * (0.5 + FLOAT(yi-1)) - 1.0 * plot_yboarder - plot_yspacing * FLOAT(yi-1)
    xpos  = [xcen - 0.5 * xsize, xcen + 0.5 * xsize]
    ypos  = [ycen - 0.5 * ysize, ycen + 0.5 * ysize]
;    xpos  = [xcen - 0.5  * xsize, xcen + 0.35 * xsize]
;    ypos  = [ycen - 0.45 * ysize, ycen + 0.5  * ysize]

    xmin =  1.0E+35
    xmax = -1.0E+35
    ymin =  1.0E+35
    ymax = -1.0E+35
    FOR idata = 1, ndata DO BEGIN
      val = cortex_ExtractStructure(vdata,idata)
;      val = cortex_ExtractStructure(mid,idata)
      xdata = val.r_a

      IF (iplot EQ 1) THEN ydata = val.dens
      IF (iplot EQ 2) THEN ydata = val.p
      IF (iplot EQ 3) THEN ydata = [val.te,val.ti]
      IF (iplot EQ 4) THEN ydata = val.zeff
      IF (iplot EQ 5) THEN ydata = val.i_frac
      IF (iplot EQ 6) THEN BEGIN
        ydata = REFORM(val.ni[0,*])
        FOR i = 1, val.z+1 DO ydata = [ydata,REFORM(val.ni[i,*])]
      ENDIF

      IF (N_ELEMENTS(WHERE(plot.xrange EQ 0.0)) NE 2) THEN BEGIN
        i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1])        
        IF (N_ELEMENTS(i) EQ 1) THEN BEGIN
          PRINT,'ERROR cortex_PlotCoreProfiles: No data within XRANGE'
          PRINT,'  IPLOT  = ',iplot
          PRINT,'  IDATA  = ',idata
          PRINT,'  XRANGE = ',plot.xrange
          PRINT,'  XDATA  = '
          PRINT,xdata
          RETURN, -1
        ENDIF
        xmin = MIN([xmin,xdata[i]])
        xmax = MAX([xmax,xdata[i]])
        ymin = MIN([ymin,ydata[i]])
        ymax = MAX([ymax,ydata[i]])
      ENDIF ELSE BEGIN
        xmin = MIN([xmin,xdata])
        xmax = MAX([xmax,xdata])
        ymin = MIN([ymin,ydata])
        ymax = MAX([ymax,ydata])
      ENDELSE
    ENDFOR
    deltay = ymax - ymin
    ymin = ymin - 0.05 * deltay
    ymax = ymax + 0.05 * deltay

    PRINT, 'XMIN,MAX=',xmin,xmax
    PRINT, 'YMIN,MAX=',ymin,ymax

;   Axes:
    xrange = [xmin,xmax]
    yrange = [ymin,ymax]
    position = [xpos[0],ypos[0],xpos[1],ypos[1]]

    CASE (iplot) OF
      focus: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=plot_ytitle[iplot-1],XTITLE=plot_xtitle
      1    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=plot_ytitle[iplot-1],XTICKFORMAT='(A1)'                                    
      2    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=plot_ytitle[iplot-1],XTICKFORMAT='(A1)',/NOERASE
      3    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=plot_ytitle[iplot-1],XTITLE=plot_xtitle,/NOERASE
      4    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=plot_ytitle[iplot-1],XTICKFORMAT='(A1)',/NOERASE
      5    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=plot_ytitle[iplot-1],XTICKFORMAT='(A1)',/NOERASE
      6    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=plot_ytitle[iplot-1],XTITLE=plot_xtitle,/NOERASE                                 
    ENDCASE

    XYOUTS, 0.95*xmin+0.05*xmax, 0.1*ymin+0.9*ymax, plot_subtitle[iplot-1]
;   Data:
    FOR idata = 1, ndata DO BEGIN
      val = cortex_ExtractStructure(vdata,idata)
;      val = cortex_ExtractStructure(mid,idata)
      xdata = val.r_a
      IF (iplot EQ 1) THEN ydata = val.dens
      IF (iplot EQ 2) THEN ydata = val.p
      IF (iplot EQ 3) THEN ydata = val.te
      IF (iplot EQ 4) THEN ydata = val.zeff    
      IF (iplot EQ 5) THEN ydata = val.i_frac
      IF (iplot EQ 6) THEN ydata = val.ni[val.z+1,*]
;      IF (iplot EQ 1) THEN ydata = val.l   
;      IF (iplot EQ 2) THEN ydata = val.dens
;      IF (iplot EQ 3) THEN ydata = val.vb
;      IF (iplot EQ 4) THEN ydata = val.machno  
;      IF (iplot EQ 5) THEN ydata = val.p
;      IF (iplot EQ 6) THEN ydata = val.te  

;     Print labels first so they are below the lines:
      IF (iplot EQ 1 or iplot EQ focus) THEN BEGIN 
        IF (focus NE 0) THEN BEGIN
          step = 0.03
        ENDIF ELSE BEGIN
          step = 0.07
        ENDELSE
        frac = step * FLOAT(idata) + 0.1
        str = val.file                   ; extract the case name from the data file
        str = STRSPLIT(str,'/',/EXTRACT)
        str = str[N_ELEMENTS(str)-1]     ; take the last sub-string
        str = STRSPLIT(str,'.',/EXTRACT)
        str = str[0]                     ; take the first one
        XYOUTS, 0.95*xmin+0.05*xmax,frac*ymin+(1.0-frac)*ymax,  $
                str, COLOR=TrueColor(colors[idata-1])
      ENDIF

      OPLOT, [0.0,0.0], [ymin,ymax], LINESTYLE=3 

      IF (iplot EQ 6 AND ndata EQ 1) THEN BEGIN
        LOADCT, 3
        FOR i = 0, val.z DO BEGIN
          frac = FLOAT(i-1) / FLOAT(val.z-1) * 0.9
          color = LONG(frac * 255.0)
          OPLOT, xdata, val.ni[i,*], COLOR=color 
;          XYouts, xdata[10], val.ni[i,10], 'a', COLOR=color, ALIGN=0.5, CHARSIZE=0.75
        ENDFOR
      ENDIF

      OPLOT, xdata, ydata, COLOR=TrueColor(colors[idata-1])
 
      IF (iplot EQ 3) THEN BEGIN
        ydata = val.ti
        OPLOT, xdata, ydata, COLOR=TrueColor(colors[idata-1]),LINESTYLE=2 
      ENDIF

    ENDFOR

  ENDFOR

  RETURN, 0

END
;
; ======================================================================
;


