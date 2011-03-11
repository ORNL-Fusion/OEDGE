;
; ======================================================================
;
FUNCTION cortex_PlotSummary, plot, plot_array, ps=ps



  PRINT
  PRINT,'----------------------- NEW PLOT -----------------------'
  PRINT


  ndata = N_ELEMENTS(TAG_NAMES(plot_array))
  IF (ndata LE 0) THEN BEGIN
    PRINT, 'ERROR cortex_PlotSummary: No data found'
    RETURN, -1
  ENDIF

  USERSYM, [-2, -2, 2, 2], [-2, 2, 2, -2], /FILL

  dev_xsize = !D.X_SIZE
  dev_ysize = !D.Y_SIZE

  xy_label = [0.95,0.05,0.10,0.90]

  xlabels = plot.xlabels

  psym = 6

  IF (plot.ylog EQ -1) THEN ylog = 1 ELSE ylog = plot.ylog

  nset = MAX(plot.case_set)

  focus = 2  ; plot.focus

  !P.BACKGROUND = TrueColor('White')
; !P.FONT = 1
; DEVICE, SET_FONT='Courier', /TT_FONT

  title    = plot.title
  notes    = plot.notes
  charsize = plot.charsize

  subtitle = ['Zeff','CONCENTRATION (%)','3','4','5','6']
  xtitle = 'bin'
  ytitle = [' ','impurity concentration (%)','3','4','5','6']

  colors = ['Black','Red','Green','Blue','Orange','Purple', 'Hotpink', 'Lime']

  plot_xboarder = 0.05
  plot_yboarder = 0.10

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
;...
    xi = plot_xi[iplot-1]
    yi = plot_yi[iplot-1]
    xcen  =       xsize * (0.5 + FLOAT(xi-1)) + 1.5 * plot_xboarder + plot_xspacing * FLOAT(xi-1)
    ycen  = 1.0 - ysize * (0.5 + FLOAT(yi-1)) - 1.0 * plot_yboarder - plot_yspacing * FLOAT(yi-1)
    xpos  = [xcen - 0.5 * xsize, xcen + 0.5 * xsize]
    ypos  = [ycen - 0.5 * ysize, ycen + 0.5 * ysize]
;...
    imax = 0
    FOR iset = 1, nset DO imax = MAX([imax,N_ELEMENTS(WHERE(plot.case_set EQ iset))])
    xmin = -0.5
    xmax = FLOAT(imax-1)+0.5

    IF (N_ELEMENTS(WHERE(plot.yrange EQ 0.0)) NE 2) THEN BEGIN
      ymin = plot.yrange[0]
      ymax = plot.yrange[1]
    ENDIF ELSE BEGIN
      ymin =  1.0E+35
      ymax = -1.0E+35
      FOR idata = 1, ndata DO BEGIN
        val = cortex_ExtractStructure(plot_array,idata)
        IF (val.core.file EQ 'none') THEN CONTINUE
        IF (iplot EQ 1) THEN ydata = MEAN(val.core.zeff)
        IF (iplot EQ 2) THEN ydata = MEAN(val.core.i_frac)
        ymin = MIN([ymin,ydata])
        ymax = MAX([ymax,ydata])
      ENDFOR
; 
      IF (N_ELEMENTS(WHERE(plot.ylimit EQ 0.0)) NE 4) THEN ymax = MAX([ymax,plot.ylimit])
      deltay = ymax - ymin
      ymin = ymin - 0.05 * deltay
      ymax = ymax + 0.05 * deltay
    ENDELSE
;    PRINT, 'XMIN,MAX=',xmin,xmax
;    PRINT, 'YMIN,MAX=',ymin,ymax

;   Axes:
    xrange = [xmin,xmax]
    yrange = [ymin,ymax]
    position = [xpos[0],ypos[0],xpos[1],ypos[1]]

    CASE (iplot) OF
      focus: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1, YLOG=ylog, XMINOR=1,  $
                   POSITION=position,YTITLE=ytitle[iplot-1],XTICKFORMAT='(A1)',CHARSIZE=charsize
      1    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=ytitle[iplot-1],XTICKFORMAT='(A1)'                                    
      2    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=ytitle[iplot-1],XTICKFORMAT='(A1)',/NOERASE
      3    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=ytitle[iplot-1],XTITLE=xtitle,/NOERASE
      4    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=ytitle[iplot-1],XTICKFORMAT='(A1)',/NOERASE
      5    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=ytitle[iplot-1],XTICKFORMAT='(A1)',/NOERASE
      6    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=ytitle[iplot-1],XTITLE=xtitle,/NOERASE                                 
    ENDCASE

;...
    FOR i = 0, N_ELEMENTS(plot.ylimit)-1 DO BEGIN
      IF (plot.ylimit[i] NE 0.0) THEN  $
        OPLOT, [xmin,xmax], [plot.ylimit[i],plot.ylimit[i]], LINESTYLE=3 
    ENDFOR
;
    XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
            (xy_label[2] * ypos[0] + xy_label[3] * ypos[1]) * dev_ysize,  $
            subtitle[iplot-1], /DEVICE
;   Print labels:
    FOR iset = 1, nset DO BEGIN
      IF (iplot EQ 1 OR iplot EQ focus) THEN BEGIN 
        IF (focus NE 0) THEN BEGIN
          step = 0.03
        ENDIF ELSE BEGIN
          step = 0.07
        ENDELSE
        frac = step * FLOAT(iset) + xy_label[2]
        str = plot.case_set_name[iset]
        XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
                (       frac * ypos[0] +  (1.0-frac) * ypos[1]) * dev_ysize,  $
                str, COLOR=TrueColor(colors[iset-1]), /DEVICE
      ENDIF
    ENDFOR
;   Data:
    xdata = FINDGEN(imax+1) - 1.0
    FOR iset = 1, nset DO BEGIN
      ydata = 0.0
      edata = 0
      FOR idata = 1, ndata DO BEGIN
        IF (plot.case_set[idata-1] NE iset) THEN CONTINUE
        val = cortex_ExtractStructure(plot_array,idata)
        IF (val.core.file EQ 'none') THEN BEGIN
          ydata = [ydata,-999.0] 
          edata = [edata,0     ] 
          CONTINUE      
        ENDIF
        IF (iplot EQ 1) THEN ydata = [ydata,MEAN(val.core.zeff  )]
        IF (iplot EQ 2) THEN ydata = [ydata,MEAN(val.core.i_frac)]
;        IF (iplot EQ 3) THEN ydata = val.core.te
;        IF (iplot EQ 4) THEN ydata = val.core.zeff    
;        IF (iplot EQ 5) THEN ydata = val.core.i_frac
;        IF (iplot EQ 6) THEN ydata = val.core.ni[val.core.z+1,*]

;       Check for problems with the simulation: 
        ecode = 0
        IF (plot.warnings) THEN BEGIN
          IF (val.summary.ions_reaching_core LT 1000.0) THEN ecode = TrueColor('LightYellow')
          IF (val.summary.ions_reaching_core LT  100.0) THEN ecode = TrueColor('Yellow')
          IF (val.summary.ions_reaching_core LT   10.0) THEN ecode = TrueColor('Greenyellow')
          IF (FLOAT(val.summary.ions_created  ) /  $
              FLOAT(val.summary.ions_requested) LT 0.99) THEN ecode = TrueColor('Lightblue')
        ENDIF
        edata = [edata,ecode]

      ENDFOR

      ydata = ydata * plot.scale_factor

      n = N_ELEMENTS(ydata) - 1
      i = 1
      j = 0
      WHILE (i LT n) DO BEGIN
        FOR j = i, n DO IF (ydata[j] EQ -999.0) THEN BREAK
        IF (i LT j) THEN BEGIN
          FOR k = i, j-1 DO IF (edata[k] NE 0) THEN OPLOT, [xdata[k]], [ydata[k]], PSYM=8, COLOR=edata[k]
          OPLOT, xdata[i:j-1], ydata[i:j-1], COLOR=TrueColor(colors[iset-1])
          OPLOT, xdata[i:j-1], ydata[i:j-1], COLOR=TrueColor(colors[iset-1]), PSYM=psym
        ENDIF
        i = j + 1
      ENDWHILE

    ENDFOR  ; iset

;   x-axis label:
    IF (focus NE 0 AND xlabels[0] NE 'unknown') THEN BEGIN
      FOR j = 0, imax-1 DO BEGIN
        xlimit = ((xpos[1] - xpos[0]) / FLOAT(imax)) * 0.95
        i =  0          
        l =  0
        k =  0.0
        WHILE (1) DO BEGIN
          k = k + 0.7 * charsize
          str = STRMID(xlabels[j],l,STRLEN(xlabels[j]))
          xwidth = FLOAT(STRLEN(STRTRIM(str,2))) * 0.0085 * charsize
          cont = 0
          WHILE (xwidth GT xlimit) DO BEGIN
            cont = 1
            i = STRPOS(str,' ',/REVERSE_SEARCH)
            IF (i EQ -1) THEN i = LONG(14.0 / charsize)
            str = STRMID(str,0,i)
            xwidth = FLOAT(STRLEN(STRTRIM(str,2))) * 0.0085 * charsize
          END
          xshift = 0.5 * xwidth / (xpos[1] - xpos[0]) * (xmax - xmin)
          xfrac  = (xdata[j+1] - xshift + 0.5) / (xdata[imax] + 1.0) 
          XYOUTS, ((1.0 - xfrac) * xpos[0] + xfrac * xpos[1]) * dev_xsize,  $
                  (ypos[0] - 0.01 - 0.03 * k                ) * dev_ysize,  $
                  STRTRIM(str,2), CHARSIZE=charsize, /DEVICE
          IF (NOT cont) THEN BREAK
          l = l + i + 1
        END

      ENDFOR
    ENDIF

  ENDFOR  ; iplot
;
;
; ----------------------------------------------------------------------  (PUT IN A FUNCTION)
  ypos = 0.940
  IF (notes NE 'default') THEN ypos = 0.965

  n = STRLEN(title) 
  nmax =  LONG(60.0 * 2.0 / (1.5 * charsize)) ; LONG(53.0 * 2.0 / (1.5 * charsize)) ;  LONG(60.0 * 2.0 / (1.5 * charsize))
  IF (n GT nmax) THEN BEGIN
    i = STRPOS(STRMID(title,0,nmax),' ',/REVERSE_SEARCH)
    str1 = STRMID(title,0,i)
    str2 = STRMID(title,i+1)
    XYOUTS, 0.02 * dev_xsize, (ypos + 0.045 * charsize / 1.7) * dev_ysize, CHARSIZE=1.5*charsize, str1, /DEVICE
    XYOUTS, 0.02 * dev_xsize, (ypos                         ) * dev_ysize, CHARSIZE=1.5*charsize, str2, /DEVICE
  ENDIF ELSE BEGIN
    XYOUTS, 0.02 * dev_xsize, ypos * dev_ysize, CHARSIZE=1.5*charsize, title, /DEVICE
  ENDELSE

  IF (notes NE 'unknown') THEN BEGIN
    n = STRLEN(notes) 
    nmax = 125; 110 ; 125;
    IF (n GT nmax) THEN BEGIN
      i = STRPOS(STRMID(notes,0,nmax),' ',/REVERSE_SEARCH)
      str1 = STRMID(notes,0,i)
      str2 = STRMID(notes,i+1)
      XYOUTS, 0.02 * dev_xsize, (ypos - 0.025) * dev_ysize, CHARSIZE=1.0, str1, /DEVICE
      XYOUTS, 0.02 * dev_xsize, (ypos - 0.050) * dev_ysize, CHARSIZE=1.0, str2, /DEVICE
      ypos = ypos - 0.020
    ENDIF ELSE BEGIN
      XYOUTS, 0.02 * dev_xsize, (ypos - 0.025) * dev_ysize, CHARSIZE=1.0, notes, /DEVICE
    ENDELSE              
  ENDIF

  RETURN, 0

END
;
; ======================================================================
;


