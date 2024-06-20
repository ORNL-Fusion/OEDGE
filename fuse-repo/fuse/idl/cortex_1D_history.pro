;
; ======================================================================
;
FUNCTION cortex_PlotHistory, plot, data_array, grid=grid, wall=wall, annotate=annotate, ps=ps

  PRINT
  PRINT,'----------------------- NEW PLOT -----------------------'
  PRINT,'PLOT OPTION=',plot.option

  MAXNYDATA = 20
  MAXNXDATA = 1000

  focus = plot.focus

  !P.BACKGROUND = TrueColor('White')

  dev_xsize = !D.X_SIZE
  dev_ysize = !D.Y_SIZE

  notes    = plot.notes
  charsize = plot.charsize
  charsize_labels = charsize

  !P.CHARSIZE  = charsize
  !P.CHARTHICK = plot.thick
  !P.THICK     = plot.thick
  !X.THICK     = plot.thick
  !Y.THICK     = plot.thick
  !Z.THICK     = plot.thick

  xy_label = [0.96,0.04,0.12,0.88]

  IF (focus) THEN BEGIN
    xy_label = [0.93,0.07,0.13,0.87]
    charsize_labels = charsize_labels * 1.2
  ENDIF

  option = plot.option

  plot_peak = plot.peak

  CASE option OF
;   --------------------------------------------------------------------
    1: BEGIN
       val = data_array.data1

       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot = 4
       plot_xn = 2
       plot_yn = [2,2]
       plot_type = [1,1,1,1]
       ntrace = [1,1,1,1] ; Number of data lines on each plot
;       ntrace = [4,4,4,4] ; Number of data lines on each plot
       title = plot.title 
       plot_subtitle = val.strings
       subtitle = plot_subtitle  ; over-written below when data is being assigned
;       subtitle = ['gas puff location / Pa','top of main chamber / Pa',  $
;                   'main / Pa','divertor / Pa']
       xtitle   = 'time (s)'
       ytitle   = ['... (unit)','... (unit)','... (unit)','... (unit)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       END
;   --------------------------------------------------------------------
    2: BEGIN

       val = data_array.data1
       str = STRSPLIT(val.file,'/',/EXTRACT)                   ; Extract case name to STR
       str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot = MAX(val.gauge)
       plot_xn = 1
       plot_yn = nplot
       plot_type = REPLICATE(1,nplot)

       ntrace = REPLICATE(MAX(val.history),nplot) ; Number of data lines on each plot

       title = plot.title + ': ' + str[0] + ': ' + STRING(val.time_end*1.0E+3,FORMAT='(F0.1)') + ' ms'
       subtitle = val.strings
       xtitle   = 'toroidal position (degrees)'
       ytitle   = REPLICATE('p_D2 (Pa)',nplot)
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)

       IF (plot.show_grid) THEN BEGIN
         plot_xn     = 2
         plot_yn     = [1,nplot]
         subtitle    = ['blank',subtitle]
         ytitle      = ['blank',ytitle]
         plot_type   = [3,REPLICATE(1,nplot)]
         ntrace      = [1,ntrace]
         iplot_start = 3
         nplot       = nplot + 1
       ENDIF

       END
;   --------------------------------------------------------------------
    3: BEGIN

       val = data_array.data1

       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot = MAX(val.gauge)
       plot_xn = 1
       plot_yn = nplot
       plot_type = REPLICATE(1,nplot)
       ntrace = REPLICATE(2,nplot) ; Number of data lines on each plot
;       ntrace = REPLICATE(MAX(val.history),nplot) ; Number of data lines on each plot
       title = plot.title 
       subtitle = val.strings
       xtitle   = 'iteration number'
       ytitle   = REPLICATE('p_D2 (Pa)',nplot)
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)

       IF (plot.show_grid) THEN BEGIN
         plot_xn     = 2
         plot_yn     = [1,nplot]
         subtitle    = ['blank',subtitle]
         ytitle      = ['blank',ytitle]
         plot_type   = [3,REPLICATE(1,nplot)]
         ntrace      = [1,ntrace]
         iplot_start = 3
         nplot       = nplot + 1
       ENDIF

       END
;   --------------------------------------------------------------------
    4: BEGIN

       val = data_array.data1

       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot = 6
       plot_xn = 2
       plot_yn = [3,3]
       plot_type = [1,1,1,1,1,1]
       ntrace = [1,1,1,1,1,1] ; Number of data lines on each plot
       title = plot.title 
       subtitle = ['ATOM DENSITY / m-3','ATOM AVERAGE ENERGY / eV','ATOM PRESSURE / Pa',  $
                   'MOLECULE DENSITY / m-3','MOLECULE AVERAGE ENERGY / eV','MOLECULE PRESSURE / Pa']
       xtitle   = 'iteration'
       ytitle   = ['n_atm (m-3)','E_avg_atm (eV)','p_atm (Pa)',  $
                   'n_mol (m-3)','E_avg_mol (eV)','p_mol (Pa)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)

       IF (plot.show_grid) THEN BEGIN
         plot_xn     = 2
         plot_yn     = [1,nplot]
         subtitle    = ['blank',subtitle]
         ytitle      = ['blank',ytitle]
         plot_type   = [3,REPLICATE(1,nplot)]
         ntrace      = [1,ntrace]
         iplot_start = 3
         nplot       = nplot + 1
       ENDIF

       END
;   --------------------------------------------------------------------
    ELSE: BEGIN  
      PRINT, 'ERROR cortex_PlotHistory: Unrecognised plot option'
      PRINT, '  OPTION = ',option,' (',plot.option,')'
      RETURN, -1
      END
  ENDCASE

  colors = ['Black','Red','Green','Blue','Darkseagreen', 'Hotpink', 'Orange', 'Silver']
;
; Setup plot:
; ----------------------------------------------------------------------
;
;
; Setup plot area:
; ----------------------------------------------------------------------

  size = TOTAL(plot.frame_bnds[1] - plot.frame_bnds[0])
;  size = TOTAL(0.5 - plot.frame_bnds[0])

  IF (focus) THEN BEGIN
    plot_xn = 1
    plot_yn = [1];    plot_yn = 1
  ENDIF

  xsize = (size - 2.0 * plot_xboarder - FLOAT(plot_xn-1) * plot_xspacing) / FLOAT(plot_xn)
;  ysize = (1.0  - 2.0 * plot_yboarder - FLOAT(plot_yn-1) * plot_yspacing) / FLOAT(plot_yn) 
  xi = 1
  yi = 0

  first_plot = 1

  FOR iplot = 1, nplot DO BEGIN
    IF (focus NE 0 AND focus NE iplot) THEN CONTINUE

    ndata = N_ELEMENTS(TAG_NAMES(data_array))

    IF (ndata LE 0) THEN BEGIN
      PRINT, 'ERROR cortex_PlotRadialProfile: No data found'
      RETURN, -1
    ENDIF

    yi = yi + 1
;    IF (yi EQ plot_yn+1) THEN BEGIN
    IF (yi EQ plot_yn[xi-1]+1) THEN BEGIN
      xi = xi + 1
      yi = 1
    ENDIF
    
    ysize = (1.0  - 2.0 * plot_yboarder - FLOAT(plot_yn[xi-1]-1) * plot_yspacing) / FLOAT(plot_yn[xi-1]) 
 
print, 'xsize,ysize',xsize,ysize,ndata

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

      strata_max = MAX(val.strata)

      CASE option OF
;       ----------------------------------------------------------------
        1: BEGIN  ; Time series

          val = data_array.data1

          file = val.file
          integral = ' '
          str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
          str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
          labels[0] = labels[0] + STRING(idata-1) + '/' + str[0] + integral + ' :'

          xdata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      
          ydata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      

;          FOR i = 1, 5 DO BEGIN
;            gauge_i = i + (iplot-1)*8
;;            if (iplot EQ 3) then gauge_i = gauge_i + 8
;            strata_i = WHERE(val.gauge EQ gauge_i AND val.strata EQ 6, count)
;;            strata_i = WHERE(val.gauge EQ gauge_i AND val.strata EQ 3, count)

          gauge = iplot ; *** FULL GEOMETRY ***

          IF (gauge EQ 4) THEN gauge = 7  ; *** FULL GEOMETRY ***

          subtitle[iplot-1] = plot_subtitle[gauge-1]

          FOR i = 1, 1 DO BEGIN


;           number = 2*i ; -1
            number = 3


            strata_i = WHERE(val.gauge EQ gauge AND val.number EQ number AND val.strata EQ strata_max, count)

print,'iplot,gauge_i',gauge,number,count


            xdata[0,i-1] = 0.0
            ydata[0,i-1] = 0.0

            xdata[1:count,i-1] = val.history[strata_i] * 0.02 ;  * 0.1
            ydata[1:count,i-1] = val.p2_mol [strata_i]

            CASE iplot OF
              1: BEGIN
                END
              2: BEGIN
                END
              3: BEGIN
                END
              4: BEGIN
                END
            ENDCASE
 
          ENDFOR

          END
;       ----------------------------------------------------------------
        2: BEGIN  ; Toroidal distribution

          IF (first_plot EQ 1) THEN BEGIN
            first_plot = 0
            file = val.file
            integral = ' '
            str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
            str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
            labels[0] = ' '
;            labels[0] = labels[0] + STRING(idata-1) + '/' + str[0] + integral + ' :'
          ENDIF

          xdata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      
          ydata = xdata

          IF (plot.show_grid) THEN igauge = iplot - 1 ELSE igauge = iplot

          IF (igauge EQ 0) THEN BREAK

          FOR itrace = 1, ntrace[iplot-1] DO BEGIN

            i = WHERE(val.history EQ itrace AND val.gauge EQ igauge AND val.strata EQ strata_max, count)

            print, 'count',itrace,igauge,count
            print, 'i',i
            print, 'number',val.number[i]

            xdata[1:count,itrace-1] = val.phi[i] ; val.number[i]
            ydata[1:count,itrace-1] = val.p2_mol[i]
 
          ENDFOR

          END
;       ----------------------------------------------------------------
        3: BEGIN  ; Iteration series - one value only

          val = data_array.data1

          file = val.file
          integral = ' '
          str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
          str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
          labels[0] = labels[0] + STRING(idata-1) + '/' + str[0] + integral + ' :'

          xdata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      
          ydata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      

;          FOR i = 1, 5 DO BEGIN
;            gauge_i = i + (iplot-1)*8
;;            if (iplot EQ 3) then gauge_i = gauge_i + 8
;            strata_i = WHERE(val.gauge EQ gauge_i AND val.strata EQ 6, count)
;;            strata_i = WHERE(val.gauge EQ gauge_i AND val.strata EQ 3, count)

          gauge = iplot ; *** FULL GEOMETRY ***

          number = 1

          strata_i = WHERE(val.gauge EQ gauge AND val.number EQ number AND val.strata EQ strata_max, count)

          print,'iplot,gauge_i',gauge,number,count

          xdata[0:count-1,0] = FINDGEN(count) + 1.0
          ydata[0:count-1,0] = val.p2_mol[strata_i]

;          CASE iplot OF
;            1: BEGIN
;              END
;            2: BEGIN
;              END
;            3: BEGIN
;              END
;            4: BEGIN
;              END
;          ENDCASE
 
          END
;       ----------------------------------------------------------------
        4: BEGIN  ; Iteration series - multi-plot    

          file = val.file
          integral = ' '
          str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
          str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
          labels[0] = labels[0] + str[0] + integral + ' :'
;          labels[0] = labels[0] + STRING(idata-1) + '/' + str[0] + integral + ' :'

          xdata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      
          ydata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      

          gauge = 4

          number = 1

          strata_i = WHERE(val.gauge EQ gauge AND val.number EQ number AND val.strata EQ strata_max, count)

          print,'iplot,gauge_i',gauge,number,count,iplot

          xdata[0:count-1,0] = FINDGEN(count) + 1.0

          CASE iplot OF
            1: ydata[0:count-1,0] = val.dens_atm[strata_i]
            2: ydata[0:count-1,0] = val.Eavg_atm[strata_i]
            3: ydata[0:count-1,0] = val.p2_atm  [strata_i]
            4: ydata[0:count-1,0] = val.dens_mol[strata_i]
            5: ydata[0:count-1,0] = val.Eavg_mol[strata_i]
            6: ydata[0:count-1,0] = val.p2_mol  [strata_i]
          ENDCASE
 
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotRadialProfile: Unrecognised plot option'
          PRINT, '  OPTION = ',plot.option
          RETURN, -1
          END
      ENDCASE

;     Package up the data for plotting:

;print, 'polot',iplot,nplot
;help,val,/struct
      name = 'data' + STRING(idata,FORMAT='(I0)')
      data = { n : ndata, x : xdata, y : ydata, file : file }  ;  *** NEED A WAY TO ONLY STORE THE NON -999.0 DATA ***

      IF (idata EQ 1) THEN data_store = CREATE_STRUCT(           name,data) ELSE  $
                           data_store = CREATE_STRUCT(data_store,name,data)

      IF (plot_type[iplot-1] EQ 1 AND N_ELEMENTS(WHERE(plot.xrange EQ 0.0)) NE 2) THEN BEGIN
        i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1])        
        IF (N_ELEMENTS(i) EQ 1) THEN BEGIN
          PRINT,'ERROR cortex_PlotWallProfile: No data within XRANGE'
          PRINT,'  IPLOT          = ',iplot
          PRINT,'  IDATA          = ',idata
          PRINT,'  XRANGE         = ',plot.xrange
          PRINT,'  XDATA MIN,MAX  = ',MIN(xdata),MAX(xdata)
          RETURN, -1
        ENDIF
      ENDIF ELSE i = WHERE(xdata NE -999.0, count)
      IF (count NE 0) THEN BEGIN
        xmin = MIN([xmin,xdata[i]])
        xmax = MAX([xmax,xdata[i]])
;       Reform YDATA for finding the maximum value:    *** ALMOST CERTAINLY A BETTER WAY TO DO THIS ***
        ydata1 = REFORM(ydata[i,0])
        IF (ndata EQ 1) THEN  $
          FOR j = 1, ntrace[iplot-1]-1 DO ydata1 = [ydata1,REFORM(ydata[i,j])]  $
        ELSE  $
          FOR j = 1, ndata-1 DO ydata1 = [ydata1,REFORM(ydata[i,j])] 
        j = WHERE(ydata1 NE -999.0,count)
        IF (count NE 0) THEN BEGIN
          ymin = MIN([ymin,ydata1[j]])
          ymax = MAX([ymax,ydata1[j]])
        ENDIF  
      ENDIF

;      FOR itrace = 1, ntrace[iplot-1] DO BEGIN
;        IF (N_ELEMENTS(WHERE(plot.xrange EQ 0.0)) NE 2) THEN BEGIN
;          i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1],count)        
;          IF (count EQ 0) THEN BEGIN
;            PRINT,'ERROR cortex_PlotRadialProfile: No data within XRANGE'
;            PRINT,'  IPLOT          = ',iplot
;            PRINT,'  IDATA          = ',idata
;            PRINT,'  XRANGE         = ',plot.xrange
;            PRINT,'  XDATA MIN,MAX  = ',MIN(xdata),MAX(xdata)
;            RETURN, -1
;          ENDIF
;        ENDIF ELSE i = WHERE(xdata[*,itrace-1] NE -999.0)
;        xmin = MIN([xmin,xdata[i,itrace-1]])
;        xmax = MAX([xmax,xdata[i,itrace-1]])
;        ymin = MIN([ymin,ydata[i,itrace-1]])
;        ymax = MAX([ymax,ydata[i,itrace-1]])
;      ENDFOR
    ENDFOR
    IF (ymin GT 0.0 AND ymax GT 0.0) THEN ymin = 0.0  ; Makes things a bit clearer on the plots I think...
    IF (ymin LT 0.0 AND ymax LT 0.0) THEN ymax = 0.0
    deltay = ymax - ymin
    ymin = ymin - 0.05 * deltay
    ymax = ymax + 0.05 * deltay

;   Axes:


    xrange = [xmin,xmax]
    IF (plot.yrange[0] NE 0.0 OR plot.yrange[1] NE 0.0) THEN  $
      yrange = plot.yrange ELSE yrange = [ymin,ymax]

    position = [xpos[0],ypos[0],xpos[1],ypos[1]]

    type = plot_type[iplot-1]
;    IF (focus NE 0 OR yi EQ plot_yn OR iplot EQ nplot) THEN plot_type = 2   ; Show x-axis label
    IF (focus NE 0 OR yi EQ plot_yn[xi-1] OR iplot EQ nplot) THEN type = 2   ; Show x-axis label

    IF (plot_type[iplot-1] EQ 3) THEN type = 3

print, 'type',iplot,nplot,type
    CASE (type) OF
      1: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black'),  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTICKFORMAT='(A1)',/NOERASE                                  
      2: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1, color=Truecolor('Black'),  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTITLE=xtitle,/NOERASE                                 
      3: BEGIN
        plot_grid        = plot
        plot_grid.size   = 0.7
        plot_grid.center = [0.6*xpos[0]+0.4*xpos[1],0.5*(ypos[0]+ypos[1])]
        status = cortex_PlotFluidGrid(plot_grid, grid, wall, 0, annotate, 'subordinate', 'outline', ps='on')
        END
    ENDCASE
    IF (type EQ 3) THEN CONTINUE

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
          OPLOT, val.x[i,0], val.y[i,0], COLOR=TrueColor(colors[idata-1])

              OPLOT, val.x[i,1], val.y[i,1], COLOR=TrueColor(colors[1]) ;, LINESTYLE=1
              OPLOT, val.x[i,2], val.y[i,2], COLOR=TrueColor(colors[2]) ;, LINESTYLE=1
              OPLOT, val.x[i,3], val.y[i,3], COLOR=TrueColor(colors[3]) ;, LINESTYLE=1


          CASE iplot OF
            1: BEGIN

              END
            2:
            3: 
            ELSE:
          ENDCASE
          END
;       ----------------------------------------------------------------
        2: BEGIN
          cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                          dev_xsize, dev_ysize, charsize_labels, colors

          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 



          LOADCT, 33


          FOR itrace = 1, ntrace[iplot-1] DO BEGIN



            frac = FLOAT(itrace-1) / FLOAT(ntrace[iplot-1]-1)

            frac_color = LONG(frac * 255.0)





            i = WHERE(val.x[*,itrace-1] NE -999.0)
            OPLOT, val.x[i,itrace-1], val.y[i,itrace-1], COLOR=frac_color
          ENDFOR

          END
;       ----------------------------------------------------------------
        3: BEGIN
          cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                          dev_xsize, dev_ysize, charsize_labels, colors

          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 

          i = WHERE(val.x[*,0] NE -999.0)
          OPLOT, val.x[i,0], val.y[i,0], COLOR=TrueColor(colors[idata-1])

              OPLOT, val.x[i,1], val.y[i,1], COLOR=TrueColor(colors[1]) ;, LINESTYLE=1
              OPLOT, val.x[i,2], val.y[i,2], COLOR=TrueColor(colors[2]) ;, LINESTYLE=1
              OPLOT, val.x[i,3], val.y[i,3], COLOR=TrueColor(colors[3]) ;, LINESTYLE=1


          CASE iplot OF
            1: BEGIN

              END
            2:
            3: 
            ELSE:
          ENDCASE
          END
;       ----------------------------------------------------------------
        4: BEGIN
          cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                          dev_xsize, dev_ysize, charsize_labels, colors

          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 

          i = WHERE(val.x[*,0] NE -999.0)
          OPLOT, val.x[i,0], val.y[i,0], COLOR=TrueColor(colors[idata-1])

              OPLOT, val.x[i,1], val.y[i,1], COLOR=TrueColor(colors[1]) ;, LINESTYLE=1
              OPLOT, val.x[i,2], val.y[i,2], COLOR=TrueColor(colors[2]) ;, LINESTYLE=1
              OPLOT, val.x[i,3], val.y[i,3], COLOR=TrueColor(colors[3]) ;, LINESTYLE=1


          CASE iplot OF
            1: BEGIN

              END
            2:
            3: 
            ELSE:
          ENDCASE
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotHistory: Unrecognised plot option'
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


