; 
;
; ======================================================================
;
FUNCTION cortex_TargetProfiles, plot, data_array, ps=ps

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
       nplot = 3
       plot_xn = 1
       plot_yn = 3
       title = plot.title 
       subtitle = ['1','2','3']
       xtitle   = 'psi_n'
       ytitle   = ['1','2','3']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace = [1]
       END
;   --------------------------------------------------------------------
    ELSE: BEGIN  
      PRINT, 'ERROR cortex_PlotTargetProfile: Unrecognised plot option'
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
      PRINT, 'ERROR cortex_PlotTargetProfile: No data found'
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
          file = val.target.file
          integral = ' '
          str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
          str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
          labels[0] = labels[0] + STRING(idata-1) + '/' + str[0] + integral + ' :'
          ntrace = [1, 1, 1]  ; Number of data lines on each plot

          i = WHERE(val.target.location EQ 2, count_i)

          ndata = MAKE_ARRAY(MAXNYDATA,/LONG,VALUE=0)
          xdata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=0.0)      
          ydata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=0.0)      

          CASE iplot OF
            1: BEGIN
              ndata[  0] = count_i
              xdata[*,0] = val.target.psin[i]
              ydata[*,0] = val.target.psin[i]
              END
            2: BEGIN
              ndata[  0] = count_i
              xdata[*,0] = val.target.psin[i]
              ydata[*,0] = val.target.psin[i]
              END
            3: BEGIN
              ndata[  0] = count_i
              xdata[*,0] = val.target.psin[i]
              ydata[*,0] = val.target.psin[i]
              END
          ENDCASE
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotTargetProfile: Unrecognised plot option'
          PRINT, '  OPTION = ',plot.option
          RETURN, -1
          END
      ENDCASE

;     Package up the data for plotting:
      name = 'data' + STRING(idata,FORMAT='(I0)')
      data = { n : ndata, x : xdata, y : ydata, file : file } 
      IF (idata EQ 1) THEN data_store = CREATE_STRUCT(           name,data) ELSE  $
                           data_store = CREATE_STRUCT(data_store,name,data)
      IF (N_ELEMENTS(WHERE(plot.xrange EQ 0.0)) NE 2) THEN BEGIN
        i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1])        
        IF (N_ELEMENTS(i) EQ 1) THEN BEGIN
          PRINT,'ERROR cortex_PlotTargetProfile: No data within XRANGE'
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
      ydata1 = REFORM(ydata[i,0])
      FOR j = 1, ntrace[iplot-1]-1 DO ydata1 = [ydata1,REFORM(ydata[i,j])] 
      ymin = MIN([ymin,ydata1])
      ymax = MAX([ymax,ydata1])
    ENDFOR
    IF (ymin GT 0.0 AND ymax GT 0.0) THEN ymin = 0.0  ; Makes things a bit clearer on the plots I think...
    IF (ymin LT 0.0 AND ymax LT 0.0) THEN ymax = 0.0
    deltay = ymax - ymin
    ymin = ymin - 0.05 * deltay
    ymax = ymax + 0.05 * deltay

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
 
          OPLOT, val.x[*,0], val.y[*,0], COLOR=TrueColor(colors[idata-1])
          CASE iplot OF
            1: 
            2: 
            3: 
            ELSE:
          ENDCASE
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotTargetProfile: Unrecognised plot option'
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
FUNCTION cortex_PlotTargetProfiles_OLD, target, ps=ps

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


