;
; *** ADD {SHOW GRID} OPTION ***
; 
;
; ======================================================================
;
FUNCTION cortex_PlotImage, plot, data_array, ps=ps

  PRINT
  PRINT,'----------------------- NEW PLOT -----------------------'
  PRINT

  MAXNYDATA = 7

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

;  LOADCT,0
  nlevels = 20
  c_colors = LONG(FINDGEN(nlevels) * (255.0 / FLOAT(nlevels-1)) ) 

;  IF (plot.zlog NE 0) THEN zlog = 1

  save_ypos = [-999.0,-999.0]

  CASE option OF
;   --------------------------------------------------------------------
    1: BEGIN
       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot    = 6
       plot_xn  = 2
       plot_yn  = 3
       title    = plot.title 
       subtitle = ['DATA / unknown','DATA / unknown','DATA / unknown',  $
                   'DATA / unknown','DATA / unknown','DATA / unknown']
       xtitle   = ['x index','x index','x index','x index','x index','x index']
       ytitle   = ['y index','y index','y index','y index','y index','y index']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace   = [ 1,1,1,2,1,1]
       type     = [-1,4,1,1,0,0]
       END
;   --------------------------------------------------------------------
    ELSE: BEGIN  
      PRINT, 'ERROR cortex_PlotImage: Unrecognised plot option'
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
    IF ((focus NE 0 AND focus NE iplot) OR (type[iplot-1] EQ 0)) THEN CONTINUE

    thick = !P.THICK

    ndata = N_ELEMENTS(TAG_NAMES(data_array))
    IF (ndata LE 0) THEN BEGIN
      PRINT, 'ERROR cortex_PlotImage: No data found'
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

    ; Fancy trick for breaking the regularity of the plots, if you don't mind
    ; me saying, so that you can have a tall plot next to some short plots:
    IF (type[iplot-1] EQ -1) THEN BEGIN
      IF (save_ypos[0] EQ -999.0) THEN save_ypos = ypos
      CONTINUE
    ENDIF ELSE IF (save_ypos[0] NE -999.0) THEN BEGIN
      ypos[1] = save_ypos[1]
      save_ypos = [-999.0,-999.0]
    ENDIF

    xmin =  1.0E+35
    xmax = -1.0E+35
    ymin =  1.0E+35
    ymax = -1.0E+35

    FOR idata = 1, ndata DO BEGIN
      val = cortex_ExtractStructure(data_array,idata)
      CASE option OF
;       ----------------------------------------------------------------
        1: BEGIN
          integral = ' '
          nintegral = N_ELEMENTS(TAG_NAMES(val.integral))

          val_data = cortex_ExtractStructure(val.integral,1)  ; *** the 1 is temporary, or should be 1 plot? ***

          file = val_data.file

          IF (nintegral EQ 1) THEN BEGIN
            str = STRSPLIT(file,'/',/EXTRACT)                   
            str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
            labels[0] = labels[0] + STRING(idata-1) + '\' + str[0] + integral + ' :'
          ENDIF ELSE BEGIN
            IF (ndata GT 1) THEN BEGIN
              PRINT,'cortex_PlotImage: Sorry, can only plot one case at the moment '
              PRINT,'                  if plotting multiple LOS integrals'
              STOP
            ENDIF
            IF (iplot EQ 1) THEN BEGIN
              str = STRSPLIT(file,'/',/EXTRACT)                 
              str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
              title = title + ' : CASE ' + str[0] 
            ENDIF
          ENDELSE

          xdata = val_data.xindex 
          ydata = val_data.yindex 
          zdata = 0.0
          nx = MAX(xdata)
          ny = MAX(ydata)
          val_data = cortex_ExtractStructure(val.integral,1)
          CASE iplot OF
            2: BEGIN
              zdata = ROTATE(TRANSPOSE(REFORM(val_data.y2,nx,ny)),3)  ; Outch!
              END
            3: BEGIN
              x   = val_data.x2
              y   = val_data.y2
              z   = val_data.z2
              r   = SQRT(x^2 + z^2)
              phi = ATAN(z / x) * 180.0 / !PI
              phi = ATAN(z,x) * 180.0 / !PI
;print,phi
;              i = WHERE(phi GT -85.0 AND phi LT -80.0)
              i = WHERE(r GT 5.55 AND r LT 5.75)
              zdata = MAKE_ARRAY(N_ELEMENTS(x),/FLOAT,VALUE=0.0)
              zdata[i] = phi[i]
;              zdata[i] = 1.0
              
;              zdata = y


;              x   = ROTATE(TRANSPOSE(REFORM(val_data.x2,nx,ny)),3)  ; Outch!
;              y   = ROTATE(TRANSPOSE(REFORM(val_data.y2,nx,ny)),3)  ; Outch!
;              z   = ROTATE(TRANSPOSE(REFORM(val_data.z2,nx,ny)),3)  ; Outch!
;              r   = SQRT(x^2 + z^2)
;              phi = ATAN(z / x) * 180.0 / !PI

              
;              zdata = SQRT(x^2 + z^2)
;;              zdata = ROTATE(TRANSPOSE(REFORM(val_data.y2,nx,ny)),3)  ; Outch!
              END
            4: BEGIN
              ; Assume vertical stripes only:
              x   = val_data.x2
              y   = val_data.y2
              z   = val_data.z2
              r   = SQRT(x^2 + z^2)
              phi = ATAN(z,x) * 180.0 / !PI
              dist = SQRT( (val_data.x1[0]-x)^2 + (val_data.y1[0]-y)^2 + (val_data.z1[0]-z)^2 )

	      ;5.55432 -4.55619  bottom of target, from i-ref-0001b
              ;5.55423 -4.07974  top

              i = WHERE(r GT  5.550 AND r LT  5.560 AND  $  ; bottom of the vertical target
                        y GT -4.556 AND y LT -4.536)
              j = WHERE(r GT  5.550 AND r LT  5.560 AND  $  ; top
                        y GT -4.100 AND y LT -4.080)
              zdata = MAKE_ARRAY(N_ELEMENTS(x),2,/FLOAT,VALUE=0.0)

              threshold = 0.5 ; degree
              FOR k = 0, N_ELEMENTS(i)-1 DO BEGIN
                FOR l = 0, N_ELEMENTS(j)-1 DO BEGIN
                  IF (ABS(phi[i[k]]-phi[j[l]]) LT threshold) THEN BEGIN
                    zdata[i[k],0] = 1.0
                    zdata[j[l],1] = 1.0
                    print,k,l,phi[i[k]],phi[j[l]],ABS(phi[i[k]]-phi[j[l]]),  $
                          xdata[i[k]],ydata[i[k]],  $
                          xdata[j[l]],ydata[j[l]],  $
                          format='(2I6,3F10.3,2X,2I6,2X,2I6)'
                  ENDIF
                ENDFOR
              ENDFOR

              save_xdata = xdata
              save_ydata = ydata
              save_zdata = zdata
              i = WHERE(save_zdata[*,0] NE 0.0, count)
              j = WHERE(save_zdata[*,1] NE 0.0)
              xdata = MAKE_ARRAY(count,2,/FLOAT,VALUE=0.0)
              ydata = xdata
              zdata = xdata
              max_y = MAX(save_ydata)
              xdata[*,0] = save_xdata[i]
              ydata[*,0] = max_y - save_ydata[i]
              zdata[*,0] = save_zdata[i,0]
              xdata[*,1] = save_xdata[j]
              ydata[*,1] = max_y - save_ydata[j]
              zdata[*,1] = save_zdata[j,1]
              dist_pixel  = SQRT( (xdata[*,0]-xdata[*,1])^2 + (ydata[*,0]-ydata[*,1])^2 )
              dist_target = SQRT( (r    [i  ]-r    [j  ])^2 + (y    [i  ]-y    [j  ])^2 )
              FOR k = 0, count-1 DO BEGIN
                print,k,  $
                      xdata[k,0],ydata[k,0],  $
                      xdata[k,1],ydata[k,1],  $
                      phi[i[k]],phi[j[k]],    $
                      dist_pixel [k],  $
                      dist_target[k],  $
                      dist_pixel[k]/dist_target[k], $ 
                      (dist_pixel[k]/dist_target[k])/(MAX(dist_pixel/dist_target)),  $
                      format='(I6,4F10.3,2X,2F10.3,2X,4F10.3)'
              ENDFOR
              PRINT,'mean distance',MEAN(dist)

              END
;           5: BEGIN
              ; Don't assume the strip is vertical, just take something close by
              ; in PHI but closer in terms of distance on the chip:
          ENDCASE
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotImage: Unrecognised plot option'
          PRINT, '  OPTION = ',plot.option
          RETURN, -1
          END
      ENDCASE

;     Package up the data for plotting:
      name = 'data' + STRING(idata,FORMAT='(I0)')
      data = { x : xdata, y : ydata, z : zdata, file : file } 
      IF (idata EQ 1) THEN data_store = CREATE_STRUCT(           name,data) ELSE  $
                           data_store = CREATE_STRUCT(data_store,name,data)
      IF (N_ELEMENTS(WHERE(plot.xrange EQ 0.0)) NE 2) THEN BEGIN
        i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1], count)        
        IF (count EQ 0) THEN BEGIN
          PRINT,'ERROR cortex_PlotImage: No data within XRANGE'
          PRINT,'  IPLOT          = ',iplot
          PRINT,'  IDATA          = ',idata
          PRINT,'  XRANGE         = ',plot.xrange
          PRINT,'  XDATA MIN,MAX  = ',MIN(xdata),MAX(xdata)
          RETURN, -1
        ENDIF
      ENDIF ELSE i = WHERE(xdata NE -999.0)
      xmin = MIN([xmin,xdata[i]])
      xmax = MAX([xmax,xdata[i]])
      IF (ntrace[iplot-1] GT 1) THEN BEGIN
        ydata1 = ydata[*,0]
        FOR j = 1, ntrace[iplot-1]-1 DO ydata1 = [ydata1,REFORM(ydata[i,j])] 
      ENDIF ELSE ydata1=ydata
      ymin = MIN([ymin,ydata1])
      ymax = MAX([ymax,ydata1])
    ENDFOR

;   Axes:
    xrange = [xmin,xmax]
    yrange = [ymin,ymax] 
    position = [xpos[0],ypos[0],xpos[1],ypos[1]]

    plot_type = default_plot_type                                           
    IF (focus NE 0 OR yi EQ plot_yn OR iplot EQ nplot) THEN plot_type = 2   ; Show x-axis label

    IF (type[iplot-1] EQ 3) THEN plot_type = 3
    IF (type[iplot-1] EQ 4) THEN plot_type = 4

    CASE (plot_type) OF
      1: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTICKFORMAT='(A1)',/NOERASE, YLOG=ylog                                  
      2: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTITLE=xtitle[iplot-1],/NOERASE, YLOG=ylog                                 
      3: CONTOUR, zdata, NLEVELS=nlevels,  $ ; MIN_VALUE=0.1, /FILL, C_COLORS=c_colors,  $ ; /FILL,  $
                  C_LABELS=[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0], $
                  POSITION=position,YTITLE=ytitle[iplot-1],XTITLE=xtitle[iplot-1],/NOERASE, YLOG=ylog                                 
      4: SURFACE, zdata, AZ=-10, CHARSIZE=2.0,  $ 
         POSITION=position,YTITLE=ytitle[iplot-1],XTITLE=xtitle[iplot-1],/NOERASE, YLOG=ylog                                 
    ENDCASE

;   Write sub-title for each plot:
    XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
            (xy_label[2] * ypos[0] + xy_label[3] * ypos[1]) * dev_ysize,  $
            subtitle[iplot-1], CHARSIZE=charsize_labels, /DEVICE

;   Data:
    FOR idata = 1, ndata DO BEGIN
      val = cortex_ExtractStructure(data_store,idata)

      cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                      dev_xsize, dev_ysize, charsize_labels, colors

      CASE option OF
;       ----------------------------------------------------------------
       1: BEGIN
         ;OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
         ;OPLOT, val.x, val.y[*,0], COLOR=TrueColor(colors[idata-1]), THICK=thick
         CASE iplot OF
;          -------------------------------------------------------------
           3: BEGIN
             i = WHERE(val.z NE 0.0 AND val.z GT -75.0 AND val.z LT -73.0)
             i = WHERE(val.z NE 0.0 AND val.z GT -45.0 AND val.z LT -43.0)
             color = MAKE_ARRAY(N_ELEMENTS(val.x[i]),/LONG,VALUE=0)
             frac = (val.z[i] - MIN(val.z[i])) / (MAX(val.z[i]) - MIN(val.z[i]))
             frac_color = LONG(frac * 255.0)
;print,MIN(val.z[i]),MAX(val.z[i])
;help,frac
;help,frac_color
;help,val.x
;help,val.z
;print,val.z[i]
             LOADCT,3
             n = N_ELEMENTS(val.x[i])
             FOR j = 0, n-1 DO  $
;              OPLOT, [val.x[i[j]]], [MAX(val.y)-val.y[i[j]]], COLOR=Truecolor('Black'), PSYM=7
;             OPLOT, [val.x[i[j]]], [val.y[i[j]]], COLOR=Truecolor('Black'), PSYM=7
               OPLOT, [val.x[i[j]]], [MAX(val.y)-val.y[i[j]]], COLOR=frac_color[j], PSYM=7
;             OPLOT, val.x[i], val.y[i], COLOR=TrueColor(colors[idata-1]), PSYM=7
             END
;          -------------------------------------------------------------
           4: BEGIN
             OPLOT, [val.x], [val.y], COLOR=Truecolor('Black'), PSYM=3
             OPLOT, [val.x], [val.y], COLOR=Truecolor('Red')  , PSYM=3
;             i = WHERE(val.z[*,0] NE 0.0)
;             OPLOT, [val.x[i]], [MAX(val.y)-val.y[i]], COLOR=Truecolor('Black'), PSYM=3
;             i = WHERE(val.z[*,1] NE 0.0)
;             OPLOT, [val.x[i]], [MAX(val.y)-val.y[i]], COLOR=Truecolor('Red')  , PSYM=3
             END
;            -----------------------------------------------------------
           ELSE:
         ENDCASE
         END
       111: BEGIN
          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
          OPLOT, val.x, val.y[*,0], COLOR=TrueColor(colors[idata-1]), THICK=thick
          CASE iplot OF
            1: 
            ELSE:
          ENDCASE
          END
;       ----------------------------------------------------------------
         2: BEGIN
;          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
;          OPLOT, val.x, val.y, COLOR=TrueColor(colors[idata-1]), THICK=thick
          i = WHERE(val.z EQ 1.0)
print,i
          OPLOT, val.x[i], val.y[i], COLOR=TrueColor(colors[idata-1]), PSYM=6
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotImage: Unrecognised plot option'
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


