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

  MAXNYDATA = 20

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

  IF (plot.ylog EQ 1) THEN ylog = 1

  save_xpos = [-999.0,-999.0]
  save_ypos = [-999.0,-999.0]

  CASE option OF
;   --------------------------------------------------------------------
    1: BEGIN
       plot_type = [1]
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot = 1
       plot_xn = 1
       plot_yn = 1
       title = plot.title 
       subtitle = ['EMISSION LINE INTEGRAL / photons ster-1 m-2 s-1']
       xtitle   = ['view index']
       ytitle   = ['integral (photons ster-1 m-2 s-1)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace   = [1]
       focus    =  1
       END
    2: BEGIN
       plot_type = [3,-1,1]
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot = 3
       plot_xn = 3
       plot_yn = 1
       title = plot.title 
       subtitle = ['none'   ,'none','EMISSION LINE INTEGRAL / photons ster-1 m-2 s-1']
       xtitle   = ['x index','none','view index']
       ytitle   = ['y index','none','integral (photons m-2 s-1 ster-1)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace   = [1,1,1]
       END
;   --------------------------------------------------------------------
    ELSE: BEGIN  
      PRINT, 'ERROR cortex_PlotImage: Unrecognised plot option'
      PRINT, '  OPTION = ',option,' (',plot.option,')'
      RETURN, -1
      END
  ENDCASE

  colors = ['Black','Red','Green','Blue','Orange','Purple', 'Hotpink', 'Darkseagreen', 'Silver']
;  colors = ['Black','Red','Blue','Orange','Purple', 'Hotpink', 'Green']
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
    IF (plot_type[iplot-1] EQ -1) THEN BEGIN
      IF (save_xpos[0] EQ -999.0) THEN BEGIN
        save_xpos = xpos
        save_ypos = ypos
      ENDIF
      CONTINUE
    ENDIF ELSE IF (save_xpos[0] NE -999.0) THEN BEGIN
      xpos[0] = save_xpos[0]
      ypos[1] = save_ypos[1]
      save_xpos = [-999.0,-999.0]
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
          print,'die'
          stop
          END
;       ----------------------------------------------------------------
        2: BEGIN
          integral = ' '
          nintegral = N_ELEMENTS(TAG_NAMES(val.integral))

          val_data = cortex_ExtractStructure(val.integral,1)      ; *** the 1 is temporary, or should be 1 plot? ***

          file = val_data.file

          j = plot.signal

          IF (ndata GT 1) THEN BEGIN
            PRINT,'cortex_PlotImage: Sorry, can only plot one case at the moment '
            PRINT,'                  if plotting multiple LOS integrals'
            STOP
          ENDIF
          IF (iplot EQ 3) THEN BEGIN
            str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
            str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
            case_name = str
            labels[iplot-1] = ''
            FOR i = 0, nintegral-1 DO BEGIN
              file =  plot.data_file[2*i]
              str = STRSPLIT(file,'.',/EXTRACT)
              str = str[N_ELEMENTS(str)-1]
              labels[iplot-1] = labels[iplot-1] + case_name +  ', ' + str + ' ' + plot.data_file[2*i+1] + ' :'
            ENDFOR
;            stop
          ENDIF

          ; Set scaling of results, based on whether a specific concentration is being forced:
          IF (plot.concentration NE 0.0) THEN BEGIN
            help,val,/struct
            help,val.core,/struct
            n = N_ELEMENTS(val.core.i_frac)
            scale = plot.concentration / val.core.i_frac[n-1] 
          ENDIF ELSE scale = 1.0
          scale = scale / (4.0 * !PI)          

          xdata = [0.0,1.0] ; val_data.xindex 
;          FOR i = 0, N_ELEMENTS(xdata)-1 DO xdata[i] = MAX([xdata[i],val_data.yindex[i]])  ; *** TEMP *** 
;          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      

          CASE iplot OF
;           ------------------------------------------------------------
            1: BEGIN
              i = 1
              val_data = cortex_ExtractStructure(val.integral,i)
              xdata = val_data.xindex 
              ydata = val_data.yindex 
              nx = MAX(xdata)
              ny = MAX(ydata)
              xdata = INDGEN(nx) + 1
              ydata = INDGEN(ny) + 1
              zdata = ROTATE(TRANSPOSE(REFORM(val_data.signal[*,j-1],nx,ny)),3)  ; Outch!
              END
;           ------------------------------------------------------------
            2: ydata = [0.0,1.0]
;           ------------------------------------------------------------
            3: BEGIN
               j = plot.signal
               ntrace[iplot-1] = nintegral
               FOR i = 1, nintegral DO BEGIN
                 val_data = cortex_ExtractStructure(val.integral,i)

                 k = LONG(STRMID(val_data.stripe,1))
                 CASE STRMID(val_data.stripe,0,1) OF
                   'r': k = WHERE(val_data.yindex EQ MAX(val_data.yindex) - k + 1, n)
                   'c': k = WHERE(val_data.xindex EQ k                           , n)  ; *** MIGHT NEED SOMETHING LIKE THE LINE ABOVE!? ***
                   ELSE: BEGIN
                     PRINT,'cortex_PlotImage: Trouble understanding specification'
                     PRINT,' SELECTOR=',STRMID(val_data.stripe,0,1)
                     END
                 ENDCASE

                 IF (i EQ 1) THEN BEGIN
                   xdata = DINDGEN(n) + 1.0
                   ydata = MAKE_ARRAY(n,MAXNYDATA,/FLOAT,VALUE=0.0)      
                 ENDIF ELSE BEGIN
                   IF (N_ELEMENTS(xdata) NE n) THEN BEGIN
;                    print,'shit!',n,(n / N_ELEMENTS(xdata))
;                    print,'k',k
                     n = N_ELEMENTS(xdata)
                     l = FINDGEN(n) / FLOAT(n-1)
;                    print,'l',l
;                    print,'l',N_ELEMENTS(k)
                     l = LONG(l * FLOAT(N_ELEMENTS(k)-1))
;                    print,'l',l
                     k = k[l]
;                    print,'k',k
                   ENDIF
                 ENDELSE

                 ydata[*,i-1] = val_data.signal[k,j-1] * scale

                 CASE STRMID(val_data.stripe,0,1) OF
                   'r': 
                   'c': ydata[*,i-1] = REVERSE(ydata[*,i-1])
                 ENDCASE

                 IF (i EQ 1) THEN BEGIN
                   IF (idata EQ 1) THEN BEGIN
                     atomic_number = val_data.atomic_number[j-1]
                     CASE atomic_number OF 
                       1: element_name = 'DEUTERIUM'
                       2: element_name = 'HELIUM'
                       4: element_name = 'BERYLLIUM'
                       ELSE: BEGIN
                         PRINT,'ERROR cortex_PlotImage: Unknown element'
                         STOP
                         END
                     ENDCASE
                     wavelength = val_data.wavelength[j-1]
                     charge     = val_data.charge    [j-1]
                     title = title + ', ' + element_name + ' +' + STRING(charge,FORMAT='(I0)') +  $
                             ', WAVELENGTH=' + STRING(wavelength,FORMAT='(F5.1)') + ' nm'

                     CASE FIX(wavelength) OF
                       656: title = title + ' (D_alpha)'
                       ELSE:
                     ENDCASE

                     IF (plot.concentration NE 0.0) THEN  $
                       title = title + ' (ION CONC.=' + STRING(plot.concentration,FORMAT='(F6.2)') + ' %)'

                   ENDIF ELSE BEGIN

                   ENDELSE
                 ENDIF ELSE BEGIN
;           stop
;                   ydata = [0.0,1.0]
                 ENDELSE



               ENDFOR
               END
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
      data = { x : xdata, y : ydata, file : file } 
      IF (idata EQ 1) THEN data_store = CREATE_STRUCT(           name,data) ELSE  $
                           data_store = CREATE_STRUCT(data_store,name,data)

      IF (N_ELEMENTS(WHERE(plot.xrange EQ 0.0)) NE 2) THEN BEGIN
        i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1])        
        IF (N_ELEMENTS(i) EQ 1) THEN BEGIN
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
;     Reform YDATA for finding the maximum value:    *** ALMOST CERTAINLY A BETTER WAY TO DO THIS ***
      ydata1 = REFORM(ydata[i,0])
      FOR j = 1, ntrace[iplot-1]-1 DO ydata1 = [ydata1,REFORM(ydata[i,j])] 
      ymin = MIN([ymin,ydata1])
      ymax = MAX([ymax,ydata1])
    ENDFOR
    IF (NOT KEYWORD_SET(ylog)) THEN BEGIN
      IF (ymin GT 0.0 AND ymax GT 0.0) THEN ymin = 0.0  ; Makes things a bit clearer on the plots I think...
      IF (ymin LT 0.0 AND ymax LT 0.0) THEN ymax = 0.0
      deltay = ymax - ymin
      ymin = ymin - 0.05 * deltay
      ymax = ymax + 0.05 * deltay
    ENDIF

;   Axes:
    xrange = [xmin,xmax]
    IF (plot.yrange[1] NE 0.0) THEN yrange = plot.yrange ELSE yrange = [ymin,ymax] 
    position = [xpos[0],ypos[0],xpos[1],ypos[1]]

    type = plot_type[iplot-1]

    IF (type NE 3 AND (focus NE 0 OR yi EQ plot_yn OR iplot EQ nplot)) THEN type = 2   ; Show x-axis label

    CASE (type) OF
      1: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTICKFORMAT='(A1)',/NOERASE, YLOG=ylog                                  
      2: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTITLE=xtitle[iplot-1],/NOERASE, YLOG=ylog                                 
      3: BEGIN
        nlevels = 100
        LOADCT,3
        c_colors = LONG(FINDGEN(nlevels) * (255.0 / FLOAT(nlevels-1)) )

        aspect_ratio = FLOAT(N_ELEMENTS(xdata)) / FLOAT(N_ELEMENTS(ydata))
        xwidth = xpos[1] - xpos[0]
        ywidth = ypos[1] - ypos[0]

       display_ratio = 0.69 ; 3.25 / 3.35  ; Small correction based on measurements from the GhostView display of the plot

        ypos[0] = ypos[1] - xwidth / (aspect_ratio * display_ratio)

    position = [xpos[0],ypos[0],xpos[1],ypos[1]]

; print,ypos
;stop

        CONTOUR, ALOG10(zdata), xdata, ydata, NLEVELS=nlevels, /FILL, C_COLORS=c_colors,  $ ; /FILL,  $
                 XRANGE=[1,MAX(xdata)], XSTYLE=1, $ 
                 YRANGE=[1,MAX(ydata)], YSTYLE=1, $
;                 C_LABELS=[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0], $
                 POSITION=position,YTITLE=ytitle[iplot-1],XTITLE=xtitle[iplot-1],/NOERASE ; , YLOG=ylog
        END
    ENDCASE

    IF (type EQ 3) THEN CONTINUE

;   Write sub-title for each plot:
    XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
            (xy_label[2] * ypos[0] + xy_label[3] * ypos[1]) * dev_ysize,  $
            subtitle[iplot-1], CHARSIZE=charsize_labels, /DEVICE

    IF (plot_peak NE 0) THEN BEGIN
      FOR idata = 1, ndata DO BEGIN
        val = cortex_ExtractStructure(data_store,idata)
        IF (idata EQ 1) THEN BEGIN
          xdata = val.x
          ydata = val.y
        ENDIF ELSE BEGIN
          FOR i = 0, N_ELEMENTS(xdata)-1 DO BEGIN
            FOR j = 0, N_ELEMENTS(ydata[0,*])-1 DO BEGIN
              IF (xdata[i] NE val.x[i]) THEN BEGIN
                print, 'trouble'
                STOP
              ENDIF
              ydata[i,j] = MAX([ydata[i,j],val.y[i,j]])
            ENDFOR
          ENDFOR
        ENDELSE
      ENDFOR
      IF (plot_peak EQ 2) THEN BEGIN
        labels[0] = 'peak values:'
        file  = 'PEAK'
        ndata = 1  
        name = 'data' + STRING(ndata,FORMAT='(I0)')
        data = { x : xdata, y : ydata, file : file } 
        data_store = CREATE_STRUCT(name,data)
      ENDIF ELSE BEGIN
        labels[0] = labels[0] + 'peak values :'
        file  = 'peak values'
        ndata = ndata + 1  
        name = 'data' + STRING(ndata,FORMAT='(I0)')
        data = { x : xdata, y : ydata, file : file } 
        data_store = CREATE_STRUCT(data_store,name,data)
      ENDELSE
    ENDIF

;   Data:
    FOR idata = 1, ndata DO BEGIN
      val = cortex_ExtractStructure(data_store,idata)

      cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                      dev_xsize, dev_ysize, charsize_labels, colors, step = 0.04

      CASE option OF
;       ----------------------------------------------------------------
        1: BEGIN
          END
;       ----------------------------------------------------------------
        2: BEGIN
          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
          IF (plot_peak EQ 1 AND idata EQ ndata) THEN thick = 2.0
          OPLOT, val.x, val.y[*,0], COLOR=TrueColor(colors[idata-1]), THICK=thick
          CASE iplot OF
            3: BEGIN
               FOR i = 1, ntrace[iplot-1]-1 DO BEGIN
                 OPLOT, val.x, val.y[*,i], COLOR=TrueColor(colors[i]), THICK=thick; , LINESTYLE=i
               ENDFOR
               END
            ELSE:
          ENDCASE
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


