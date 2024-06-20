;
; *** ADD {SHOW GRID} OPTION ***
; 
;
; ======================================================================
;
FUNCTION cortex_PlotIntegrals, plot, data_array, ps=ps

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

  CASE option OF
;   --------------------------------------------------------------------
    1: BEGIN
       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot = 1
       plot_xn = 1
       plot_yn = 1
       title = plot.title 
       subtitle = ['EMISSION LINE INTEGRAL / photons m-2 s-1 ster-1']
       xtitle   = 'VIEW INDEX (TOROIDAL DISPLACEMENT)'
       ytitle   = ['integral (photons m-2 s-1 ster-1)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace = [1]
       focus = 1
       END
    2: BEGIN
       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot = 1
       plot_xn = 1
       plot_yn = 1
       title = plot.title 
       subtitle = ['PROFILE ALONG THE DIAGNOSTIC LINE-OF-SIGHT / ' + plot.id]
       xtitle   = 'distance from last lens (m)' 
       ytitle   = ['signal (arb)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace = [1]
       END
;   --------------------------------------------------------------------
    ELSE: BEGIN  
      PRINT, 'ERROR cortex_PlotIntegrals: Unrecognised plot option'
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


  plot.frame_bnds[0] = 0.0  ; *** HACK ***
  plot.frame_bnds[1] = 1.0


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
      PRINT, 'ERROR cortex_PlotIntegrals: No data found'
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
;          integral = '   PARTICLE FLUX INTEGRAL= ' + STRING(val.integral,FORMAT='(E12.4)')
          integral = ' '

          nintegral = N_ELEMENTS(TAG_NAMES(val.integral))

          val_data = cortex_ExtractStructure(val.integral,1)      ; *** the 1 is temporary, or should be 1 plot? ***

          file = val_data.file
;          print,'file= ',file

          IF (nintegral EQ 1) THEN BEGIN
            str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
            str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
            labels[0] = labels[0] + STRING(idata-1) + '\' + str[0] + integral + ' :'
          ENDIF ELSE BEGIN
            IF (ndata GT 1) THEN BEGIN
              PRINT,'cortex_1D_integrals: sorry, can only plot one case at the moment '
              PRINT,'                     if plotting multiple LOS integrals'
              STOP
            ENDIF
            IF (iplot EQ 1) THEN BEGIN
              str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
              str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
              case_name = str
;              labels[0] = labels[0] + STRING(idata-1) + '\' + str[0] + integral + ' :'
;              title = title + ': CASE ' + str[0] 

              labels[0] = ''
              FOR i = 0, nintegral-1 DO BEGIN
                file =  plot.data_file[i]
                str = STRSPLIT(file,'.',/EXTRACT)
                str = str[N_ELEMENTS(str)-1]
                labels[0] = labels[0] + case_name +  ', ' + str + ' :'
              ENDFOR

            ENDIF
          ENDELSE

          ntrace = [1, 1, 1]  ; Number of data lines on each plot
;          ntrace = [2, 1, 1]  ; C-Mod
          xdata = val_data.xindex 
          FOR i = 0, N_ELEMENTS(xdata)-1 DO xdata[i] = MAX([xdata[i],val_data.yindex[i]])  ; *** TEMP *** 
          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      

          ; Set scaling of results, based on whether a specific concentration is being forced:
          IF (plot.concentration NE 0.0) THEN BEGIN
;            help,val,/struct
;            help,val.core,/struct
            n = N_ELEMENTS(val.core.i_frac)
            scale = plot.concentration / val.core.i_frac[n-1] 
          ENDIF ELSE scale = 1.0
          scale = scale / (4.0 * !PI)          

          CASE iplot OF
            1: BEGIN
;help,val_data,/struct
;help,ydata
;print,'scale= ',scale
               j = plot.signal
               ntrace[iplot-1] = nintegral
               FOR i = 1, nintegral DO BEGIN
                 val_data = cortex_ExtractStructure(val.integral,i)
                 ydata[*,i-1] = val_data.signal[*,j-1] * scale
 
                 IF (i EQ 1) THEN BEGIN
;                 IF (nintegral EQ 1) THEN BEGIN
;print,val_data.atomic_number
;print,val_data.atomic_number[j-1]
                   IF (idata EQ 1) THEN BEGIN
                     atomic_number = val_data.atomic_number[j-1]
                     CASE atomic_number OF 
                       1: element_name = 'DEUTERIUM'
                       2: element_name = 'HELIUM'
                       4: element_name = 'BERYLLIUM'
                       6: element_name = 'CARBON'
                       10: element_name = 'NEON'
                       ELSE: BEGIN
                         PRINT,'ERROR cortex_PlotIntegrals: Unknown element'
                         STOP
                         END
                     ENDCASE
                     wavelength = val_data.wavelength[j-1]
                     charge     = val_data.charge    [j-1]
                     title = title + ', ' + element_name + ' +' + STRING(charge,FORMAT='(I0)') +  $
                             ', WAVELENGTH=' + STRING(wavelength,FORMAT='(F5.1)') + ' nm'

                     CASE FIX(wavelength) OF
                       486: title = title + ' (D_beta)'
                       656: title = title + ' (D_alpha)'
                       ELSE:
                     ENDCASE

                     IF (plot.concentration NE 0.0) THEN  $
                       title = title + ' (ION CONC.=' + STRING(plot.concentration,FORMAT='(F6.2)') + ' %)'

                   ENDIF ELSE BEGIN

                   ENDELSE
                 ENDIF ELSE BEGIN
;                   print, 'not sure what to do here with multiple lines'
;                   stop
                 ENDELSE


;                 labels[0] = labels[0] + STRING(i-1) + '\' + plot.data_file[i-1] + ':' 
               ENDFOR
;               ydata[*,1] = [2.5,4.5,1.65,1.6,0.1,0.028,0.02,0.02,0.02] * 1E+21 * 4.0 * 3.1415  ; C-Mod

;   print,xdata
;   print,ydata[*,0]

               END

          ENDCASE

    
          IF (1 EQ 1) THEN BEGIN           

            fp1 = 3
            FREE_LUN,fp1        
            str = STRSPLIT(plot.data_file[0],'.',/EXTRACT)  
            file_name = 'cortex_data/'+plot.case_name[0] + '.los_' + str[1]
            print,file_name

            OPENW,fp1,file_name, ERROR=err
            IF (err NE 0) THEN BEGIN
              PRINT,'ERROR cortex_Plot1DIntegrals: Problem opening data stream (to file)'
              PRINT,'FILE_NAME= ',file_name
              PRINTF,-2,!err.msg
              FREE_LUN,fp1
              STOP
            ENDIF
     	  
            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
            PRINTF,fp1,'* LINE-OF-SITE INTEGRAL DATA FROM RAY (THIS DATA FILE GENERATED IN CORTEX)'
            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
            PRINTF,fp1,'* CASE             ',plot.case_name[0],FORMAT='(2A)'
            PRINTF,fp1,'* TITLE            ',FORMAT='(A)'
            PRINTF,fp1,'* DATE AND TIME    ',FORMAT='(A)
            PRINTF,fp1,'* VIEW IDENTIFIER  ',str(1),FORMAT='(2A)'
            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
            PRINTF,fp1,'{DATA FILE VERSION}',FORMAT='(A)'
            PRINTF,fp1,'     1.0',FORMAT='(A)'
            PRINTF,fp1,'*',FORMAT='(A)
            PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
            PRINTF,fp1,'{EMISSION SIGNALS}',FORMAT='(A)'
            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* Z          - atomic number of particle',FORMAT='(A)'
            PRINTF,fp1,'* A          - atomic weight',FORMAT='(A)'
            PRINTF,fp1,'* CHARGE     - electric charge',FORMAT='(A)'
            PRINTF,fp1,'* WAVELENGTH - of the signal/line [nm]',FORMAT='(A)'
            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* ','SIGNAL','Z'   ,'A'   ,'CHARGE','WAVELENGTH',FORMAT='(A2,A6,3A8,A12)'
            PRINTF,fp1,'* ','     ','    ','     ','      ','[nm]'      ,FORMAT='(A2,A6,3A8,A12)'
            FOR i = 0, N_ELEMENTS(val_data.atomic_number)-1 DO  $
              PRINTF,fp1,i+1,  $
                         val_data.atomic_number[i],  $
                         val_data.atomic_mass  [i],  $
                         val_data.charge       [i],  $
                         val_data.wavelength   [i],  $
                         FORMAT='(2X,I6,3I8,F12.2)'

            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
            PRINTF,fp1,'{LINES-OF-SIGHT}',FORMAT='(A)'
            PRINTF,fp1,N_ELEMENTS(val_data.signal),FORMAT='(I12)'
            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* max,view_25/max = ',MAX(val_data.signal),val_data.signal[25]/MAX(val_data.signal),FORMAT='(A,E10.2,F10.4)'  
            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* INDEX    - index number of the viewing chord',FORMAT='(A)'
            PRINTF,fp1,'* (X,Y,Z)1 - location of the pupil in machine coordinates (X=R,Y=Z) [m]',FORMAT='(A)'
            PRINTF,fp1,'* (X,Y,Z)2 - end point of the chord [m]',FORMAT='(A)'
            PRINTF,fp1,'* SIGNAL   - line-of-sight emission integral in units of [photons ster-1 m-2 s-1]',FORMAT='(A)'
            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* ','INDEX','X1' ,'Y1' ,'Z1' ,'X2' ,'Y2' ,'Z2' ,'SIGNAL_1',FORMAT='(A2,A6,2(3A10,2X),A17)'
            PRINTF,fp1,'* ',' '    ,'[m]','[m]','[m]','[m]','[m]','[m]','[ph st-1 m-2 s-1]',FORMAT='(A2,A6,2(3A10,2X),A17)'
	  
            FOR i2 = 0, N_ELEMENTS(val_data.signal)-1 DO BEGIN
              PRINTF,fp1,  $
                i2,  $
                val_data.x1[i2],  $
                val_data.y1[i2],  $
                val_data.z1[i2],  $
                val_data.x2[i2],  $
                val_data.y2[i2],  $
                val_data.z2[i2],  $
                val_data.signal[i2],  $
                FORMAT='(2X,I6,2(3F10.4,2X),E17.2)'
            ENDFOR
	  
            CLOSE,fp1
            FREE_LUN,fp1

          ENDIF


          END
;       ----------------------------------------------------------------
        2: BEGIN
          integral = ' '

;          help,val_array.integral,/struct
;stop
;          ntrace = N_ELEMENTS(TAG_NAMES(data_array))
;          IF (ntrace LE 0) THEN BEGIN
;            PRINT, 'ERROR cortex_PlotIntegrals: No data found'
;            RETURN, -1
;          ENDIF

          val = cortex_ExtractStructure(val.profile,1)      ; *** the 1 is temporary, or should be 1 plot? ***

;    help,val,/struct

          file = val.file
          str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
          str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
          labels[0] = labels[0] + STRING(idata-1) + '/' + str[0] + integral + ' :'

          ntrace = [1, 1, 1]  ; Number of data lines on each plot

          xdata = val.path
          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
          CASE iplot OF
            1: ydata[*,0] = val.signal
          ENDCASE


          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotIntegrals: Unrecognised plot option'
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
          PRINT,'ERROR cortex_PlotIntegrals: No data within XRANGE'
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

    plot_type = default_plot_type                                           
    IF (focus NE 0 OR yi EQ plot_yn OR iplot EQ nplot) THEN plot_type = 2   ; Show x-axis label

    CASE (plot_type) OF
      1: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTICKFORMAT='(A1)',/NOERASE, YLOG=ylog                                  
      2: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTITLE=xtitle,/NOERASE, YLOG=ylog                                 
    ENDCASE

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
                      dev_xsize, dev_ysize, charsize_labels, colors

      CASE option OF
;       ----------------------------------------------------------------
        1: BEGIN

          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
 
          IF (plot_peak EQ 1 AND idata EQ ndata) THEN thick = 2.0

          OPLOT, val.x, val.y[*,0], COLOR=TrueColor(colors[idata-1]), THICK=thick
          CASE iplot OF
            1: BEGIN
;               print,'ntrace:',ntrace[iplot-1]
               FOR i = 1, ntrace[iplot-1]-1 DO BEGIN
                 OPLOT, val.x, val.y[*,i], COLOR=TrueColor(colors[i]), THICK=thick; , LINESTYLE=i
; print,i,val.y[*,i]
               ENDFOR
               END
            2: 
            3: 
            ELSE:
          ENDCASE
          END
;       ----------------------------------------------------------------
        2: BEGIN

          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
 
          OPLOT, val.x, val.y[*,0], COLOR=TrueColor(colors[idata-1]), THICK=thick
          CASE iplot OF
            1: 
            2: 
            3: 
            ELSE:
          ENDCASE
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotIntegrals: Unrecognised plot option'
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


