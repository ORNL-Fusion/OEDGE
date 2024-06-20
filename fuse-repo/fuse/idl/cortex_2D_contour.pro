;
; ======================================================================
;
FUNCTION cortex_SetColour, range, frac, lastct

  IF (0 EQ 1 AND range[0] LT 0.0 AND range[1] GT 0.0) THEN BEGIN
    frac = 2.0 * (frac - 0.5)
    IF (frac LT 0.0 AND lastct NE 1) THEN BEGIN
      LOADCT, 1, /SILENT  
      lastct = 1
    ENDIF 
    IF (frac GE 0.0 AND lastct NE 3) THEN BEGIN
      LOADCT, 3, /SILENT 
      lastct = 3
    ENDIF
  ENDIF

  result = LONG(ABS(frac) * 255.0)

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_GetSideVertices, grid, iobj, iside

  isrf = grid.obj_iside[iside,iobj]
  IF (isrf EQ 0) THEN BEGIN
    PRINT, 'ERROR cortex_GetSideVertices: Surface index is zero'
    STOP
  ENDIF

  i = WHERE(grid.srf_index EQ ABS(isrf))                              ; Should replace with a mapping array?
  IF (N_ELEMENTS(i) NE 1 OR i EQ -1) THEN BEGIN
    PRINT, 'ERROR cortex_GetSideVertices: Unable to map surface index'
    STOP
  ENDIF

  IF (isrf GT 0) THEN BEGIN
    j = grid.vtx_map[grid.srf_ivtx[0                 ,i]]
  ENDIF ELSE BEGIN
    j = grid.vtx_map[grid.srf_ivtx[grid.srf_nvtx[i]-1,i]]
  ENDELSE

  result = grid.vtx[*,j]

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_Plot2DContour, plot, plot_data, grid, wall, annotate, ps=ps, time_slice=time_slice

  flip = plot.flip  ; for ribbon grid plots

  window_id = 0
  window_xsize = 700
  window_ysize = 700

  dev_xsize = !D.X_SIZE
  dev_ysize = !D.Y_SIZE
  !P.BACKGROUND = TrueColor('White')

  title    = plot.title
  notes    = plot.notes
  charsize = plot.charsize

  !P.CHARSIZE  = charsize
  !P.CHARTHICK = plot.thick
  !P.THICK     = plot.thick
  !X.THICK     = plot.thick
  !Y.THICK     = plot.thick
  !Z.THICK     = plot.thick

  xy_label = [0.08,0.92,0.07,0.93]

  xtitle = 'R (m)'
  ytitle = 'Z (m)'
  IF (flip) THEN BEGIN
    xtitle = 'distance parallel to field, s (m)'
    ytitle = 'radial distance (m)'
  ENDIF

;
; Setup plot area:
; ----------------------------------------------------------------------
  size = 0.8
  IF (plot.size NE 0.0) THEN size = plot.size

  IF (N_ELEMENTS(WHERE(plot.center EQ 0.0)) EQ 2) THEN BEGIN 
    xcen = 0.5 * TOTAL(plot.frame_bnds) 
    ycen = 0.5
  ENDIF ELSE BEGIN
    IF (plot.center[0] NE 0.0) THEN xcen = plot.center[0]  $
                               ELSE xcen = 0.5 * TOTAL(plot.frame_bnds)
    IF (plot.center[1] NE 0.0) THEN ycen = plot.center[1]  $
                               ELSE ycen = 0.5
  ENDELSE  

;
; Setup plot zoom:
; ----------------------------------------------------------------------
  IF (N_ELEMENTS(WHERE(plot.zoom EQ 0.0)) EQ 4) THEN BEGIN 
    IF (plot.option EQ 900) THEN BEGIN
      xmin = MIN(grid.x)
      xmax = MAX(grid.x)
      ymin = MIN(grid.x)
      ymax = MAX(grid.x)
    ENDIF ELSE IF (plot.option EQ 901) THEN BEGIN
      i = WHERE(plot_data.x NE 0.0)
      xmin = MIN(plot_data.x[i]) - 0.1
      xmax = MAX(plot_data.x[i]) + 0.1
      i = WHERE(plot_data.x NE 0.0)
      ymin = MIN(plot_data.y[i]) - 0.1
      ymax = MAX(plot_data.y[i]) + 0.1
    ENDIF ELSE BEGIN
      xmin = MIN([[grid.vtx[0,*]],[wall.v1[0,*]],[wall.v2[0,*]]]) - 0.1
      xmax = MAX([[grid.vtx[0,*]],[wall.v1[0,*]],[wall.v2[0,*]]]) + 0.1
      ymin = MIN([[grid.vtx[1,*]],[wall.v1[1,*]],[wall.v2[1,*]]]) - 0.1
      ymax = MAX([[grid.vtx[1,*]],[wall.v1[1,*]],[wall.v2[1,*]]]) + 0.1
    ENDELSE
  ENDIF ELSE BEGIN
    xmin = plot.zoom[0]
    xmax = plot.zoom[2]
    ymin = plot.zoom[1]
    ymax = plot.zoom[3]
  ENDELSE
  xdelta = xmax - xmin
  ydelta = ymax - ymin
;
;
; ----------------------------------------------------------------------
  IF (KEYWORD_SET(ps)) THEN BEGIN
    xsize = 297.0  ; Landscape A4
    ysize = 210.0  
    display_ratio = 3.25 / 3.35  ; Small correction based on measurements from the GV display
    aspect_ratio = (xdelta / xsize) / (ydelta / ysize) * display_ratio
  ENDIF ELSE BEGIN
    PRINT, 'NEED TO FIX WINDOW ASPECT RATIO'
    RETURN, -1
    display_ratio = 4.4
    aspect_ratio = xdelta / ydelta / display_ratio *            $ 
                   (FLOAT(window_ysize) / FLOAT(window_xsize))
    WINDOW, window_id, XSIZE=window_xsize, YSIZE=window_ysize
  ENDELSE
  IF (plot.aspect_ratio GT 0.0) THEN aspect_ratio = plot.aspect_ratio
;
;
  IF ((xdelta / xsize) LT (ydelta / ysize)) THEN BEGIN
    xpos = [xcen - 0.5 * size * aspect_ratio,  $
            xcen + 0.5 * size * aspect_ratio]
    ypos = [ycen - 0.5 * size, ycen + 0.5 * size]
  ENDIF ELSE BEGIN
    xpos = [xcen - 0.5 * size, xcen + 0.5 * size]
    ypos = [ycen - 0.5 * size / aspect_ratio,  $
            ycen + 0.5 * size / aspect_ratio]
  ENDELSE
;
; Check if the plot fits in the allocated portion of the page, and scale the plot
; as necessary:
  IF (NOT plot.rigid AND (xpos[0] LT plot.frame_bnds[0] OR xpos[1] GT plot.frame_bnds[1])) THEN BEGIN
    ratio = 0.8 * (plot.frame_bnds[1] - plot.frame_bnds[0]) / (xpos[1] - xpos[0]) 
    xpos = [xcen - ratio * 0.5 * (xpos[1] - xpos[0]),  $
            xcen + ratio * 0.5 * (xpos[1] - xpos[0])]
    ypos = [ycen - ratio * 0.5 * (ypos[1] - ypos[0]),  $
            ycen + ratio * 0.5 * (ypos[1] - ypos[0])]
  ENDIF
;
; Store the plot information
; ----------------------------------------------------------------------
  plot.zoom = [xmin,ymin,xmax,ymax]
  plot_xrange    = [xmin,xmax]
  plot_yrange    = [ymin,ymax]  
;  plot.xrange   = [xmin,xmax]
;  plot.yrange   = [ymin,ymax]
  plot.position = [xpos[0],ypos[0],xpos[1],ypos[1]]
;
; Draw the plot:
; ----------------------------------------------------------------------
;
  color_table = plot.color_table

  state = plot.state
  IF (state EQ -1) THEN state = plot_data.div_max_iz + 1

  plot_labels = MAKE_ARRAY(1000,/STRING)

  plot_labels[1  ] = 'ELECTRON DENSITY          : n_e (m-3)'
  plot_labels[2  ] = 'PARALLEL PLASMA FLOW      : v_parallel (m s-1)'
  plot_labels[3  ] = 'PARALLEL FLOW MACH NUMBER : Mach no.'
  plot_labels[4  ] = 'ELECTRON PRESSURE         : p_e (Pa)'
  plot_labels[5  ] = 'ION PRESSURE              : p_i (Pa)'
  plot_labels[6  ] = 'TOTAL PLASMA PRESSURE     : p (Pa)'
  plot_labels[7  ] = 'ELECTRON TEMPERATURE      : T_e (eV)'
  plot_labels[8  ] = 'ION TEMPERATURE           : T_i (eV)'

  plot_labels[200] = 'EIRENE ATOM DENSITY : n_D (m-3)' 
  plot_labels[201] = 'EIRENE MOLECULAR DENSITY : n_D2 (m-3)' 
  plot_labels[220] = 'EIRENE BALMER ALPHA : D_alpha (ph m-3 s-1)' 
  plot_labels[222] = 'EIRENE BALMER GAMMA : D_gamma (ph m-3 s-1)' 

  plot_labels[300] = 'EIRENE IMPURITY ATOM DENSITY: n_I+0 (m-3)' 

  IF (plot.option EQ 400) THEN BEGIN
    IF (state EQ plot_data.div_max_iz+1) THEN  $
      plot_labels[400] = 'DIVIMP IMPURITY DENSITY: n_i all states (m-3)'  $
    ELSE  $
      plot_labels[400] = 'DIVIMP IMPURITY DENSITY: n_i +'+STRTRIM(STRING(state),2)+' (m-3)'

    IF (state GT plot_data.div_max_iz AND state NE plot_data.div_max_iz+1) THEN BEGIN
      print,'!!! HACK in 2D contour: NO DATA SO SETTING STATE=0 !!!', state, plot_data.div_max_iz
      state = 0
    ENDIF
  ENDIF

  plot_labels[900] = 'WALLDYN REDISTRIBUTION MATRIX ' + plot.case_name[0] + ' : destination probability'
  plot_labels[901] = 'GRID SLICE, PORT 18, ' + plot.case_name[0] 

  str = STRSPLIT(plot_data.file,'/',/EXTRACT)                   ; Extract case name to STR
  case_name = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)

  CASE plot.option OF
;   ------------------------------------------------------------------
    1: data = plot_data.dens
    2: data = plot_data.vi
    3: BEGIN
       data = plot_data.vi / plot_data.cs
       IF (flip) THEN data = -data
       END
    4: data = plot_data.pe
    5: data = plot_data.pi
    6: data = plot_data.pe + plot_data.pi
    7: data = plot_data.te
    8: data = plot_data.ti
;   ------------------------------------------------------------------
    200: data = plot_data.atm_dens
    201: data = plot_data.mol_dens
    220: BEGIN
        data = plot_data.balmer_alpha
        i = WHERE(plot_data.tube LE 12, count)
;        i = WHERE(data GT 1.E+25, count)
        IF (count GT 0) THEN BEGIN
          data[i] = 0.0
          print,'*** zeroing out super-high Halpha values ***'
        ENDIF
        END
    222: data = plot_data.balmer_gamma
;   ------------------------------------------------------------------
    300: data = plot_data.imp_dens[*,0]  ; EIRENE
;   ------------------------------------------------------------------
    400: data = plot_data.imp_dens[*,state] * plot_data.div_influx   ; DIVIMP
;   ------------------------------------------------------------------
    500: BEGIN
      plot_labels[500] = 'DIVIMP IMPURITY EMISSION : ' + STRTRIM(STRING(plot_data.wavelength[plot.signal-1]),2)+'nm (ph m-3 s-1)'
      data = plot_data.signal[*,plot.signal-1] 
      END
;   ------------------------------------------------------------------
    900: BEGIN
      xtitle = 'SOURCE WALL INDEX'
      ytitle = 'DESTINATION WALL INDEX'
      ticklen = -0.01
      data_array = plot_data.redep_matrix
      ndata = N_ELEMENTS(TAG_NAMES(data_array))
      FOR idata = 1, ndata DO BEGIN
        val = cortex_ExtractStructure(data_array,idata)
        print,val.file
        IF (STRPOS(val.file, plot.case_name[0]) NE -1) THEN BEGIN
          data = val.data
          BREAK
        ENDIF        
      ENDFOR
      IF (idata EQ ndata+1) THEN BEGIN
        PRINT, 'cortex_Plot2DContour: WALLDYN redistribution matrix data not found'
        PRINT, ' CASE = ' + plot.case_name[0]
      ENDIF
      END
;   ------------------------------------------------------------------
    901: BEGIN
      xtitle = 'R (m)'
      ytitle = 'Z (m)' 

      plot_labels[901] = plot_labels[901] + ': D2 pressure (Pa)'

      ndata = 1

      val = plot_data
      print,val.file

      data = val.data

      print,'plot option',val.opt

      grid = { obj_n : N_ELEMENTS(data) } 

      END
;   ------------------------------------------------------------------
    ELSE: BEGIN
      PRINT, 'ERROR cortex_PlotContour: Unrecognised option'
      PRINT, '  PLOT OPTION = ',plot.option
      RETURN, -1
      END
  ENDCASE

;  *** FIND A BETTER WAY TO SORT THROUGH THESE ***
  plot_index = plot.option
;  IF (plot.option GE 100 AND plot.option LE 199) THEN plot_index = plot_index -  99
;  IF (plot.option GE 200 AND plot.option LE 299) THEN plot_index = plot_index - 199
;  IF (plot.option GE 300 AND plot.option LE 399) THEN plot_index = plot_index - 299
;  IF (plot.option GE 400 AND plot.option LE 499) THEN plot_index = plot_index - 399

  labels = STRSPLIT(plot_labels[plot_index],':',/EXTRACT)
  IF (plot.plot_title NE 'unknown') THEN plot_title = plot.plot_title ELSE   $
                                         plot_title = STRTRIM(labels[0],2)
  IF (plot.no_title) THEN plot_title = ' '

  scale_label = STRTRIM(labels[1],2)

  IF (plot.log) THEN FOR i = 0L, N_ELEMENTS(data)-1 DO  $                    ; A nicer way to do this..?
    IF (data[i] GT 0.0) THEN data[i] = ALOG10(data[i]) ELSE data[i] = 0.0

  xrange = [MIN(data),MAX(data)]
  IF (N_ELEMENTS(WHERE(plot.xrange EQ 0.0)) NE 2) THEN BEGIN
    IF (plot.xrange[0] NE -99.0) THEN xrange[0] = plot.xrange[0]
    IF (plot.xrange[1] NE -99.0) THEN xrange[1] = plot.xrange[1]
    IF (plot.xrange[0] EQ -88.0) THEN xrange[0] = xrange[1] - 4.0
  ENDIF

  xstyle = plot.xstyle
  ystyle = plot.ystyle

  IF (flip) THEN BEGIN
    plot_xrange_save = plot_xrange
    plot_xrange = plot_yrange
    plot_yrange = plot_xrange_save
  ENDIF

  IF (plot.frame_bnds[0] EQ 0.0) THEN BEGIN
    ytickformat = ''
  ENDIF ELSE BEGIN
    ytickformat = '(A1)'
    ytitle = ''
  ENDELSE

  PLOT, plot_xrange, plot_yrange, /NODATA, XSTYLE=xstyle, YSTYLE=ystyle, /NOERASE,  $
        TICKLEN =ticklen      ,                                           $
        POSITION=plot.position,                                           $
        TITLE=plot_title, XTITLE=xtitle, YTITLE=ytitle,                   $
        COLOR=TrueColor('Black'), CHARSIZE=1.0, YTICKFORMAT=ytickformat

;   Write sub-title for each plot:

  IF (plot.frame_bnds[0] EQ 0.0 AND time_slice NE -1) THEN BEGIN
    ; aps setting
    ; time_str = STRTRIM(STRING(time_slice*10),2) + ' us'
    time_str = STRTRIM(STRING(time_slice*2),2) + ' ms'
    XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
            (xy_label[2] * ypos[0] + xy_label[3] * ypos[1]) * dev_ysize,  $
            time_str, CHARSIZE=charsize_labels, /DEVICE, ALIGNMENT=1.0
  ENDIF ELSE BEGIN
    xy_label = [0.05,0.95,0.97,0.03]
    XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
            (xy_label[2] * ypos[0] + xy_label[3] * ypos[1]) * dev_ysize,  $
            plot.case_name[0], CHARSIZE=0.95, /DEVICE, ALIGNMENT=1.0
  ENDELSE



  LOADCT, color_table

  IF (N_ELEMENTS(data) NE grid.obj_n) THEN BEGIN
    PRINT
    PRINT, '------------------------------------------'
    PRINT, '  WARNING cortex_Plot2D: NDATA NE NOBJ'
    PRINT, '------------------------------------------'
    PRINT
  ENDIF

  lastct = -1

  FOR iobj = 0L, grid.obj_n-1 DO BEGIN
;
;   Filter(s):
;    IF (grid.obj_tube[iobj] LT 13 OR grid.obj_tube[iobj] GT 13) THEN CONTINUE
;    IF (grid.obj_tube[iobj] LT 7 OR grid.obj_tube[iobj] GT 7) THEN CONTINUE
;    IF (grid.obj_tube[iobj] NE 8) THEN CONTINUE
;
    IF (plot.option EQ 900) THEN BEGIN
      v = MAKE_ARRAY(3,4,/FLOAT,VALUE=0.0)
      n = SQRT(N_ELEMENTS(data))
      i = FIX(iobj / n)
      j = FIX(iobj MOD n)
      data_val = data[i,j]
      v = [ [grid.x[i+1],grid.x[j  ]], [grid.x[i   ],grid.x[j  ]],  $
            [grid.x[i  ],grid.x[j+1]], [grid.x[i+1 ],grid.x[j+1]] ]
    ENDIF ELSE IF (plot.option EQ 901) THEN BEGIN

      n = val.nv[iobj]
      data_val = data[iobj]
      v = [ val.x[iobj,0:n-1] , val.y[iobj,0:n-1] ]

    ENDIF ELSE BEGIN
      v = MAKE_ARRAY(3,grid.obj_nside[iobj],/FLOAT,VALUE=0.0)
      FOR iside = 0, grid.obj_nside[iobj]-1 DO  $
        v[*,iside] = cortex_GetSideVertices(grid,iobj,iside)
      data_val = data[iobj]
    ENDELSE

    frac = (data_val - xrange[0]) / (xrange[1] - xrange[0])
    frac = MAX([frac,0.0])
    frac = MIN([frac,1.0])

    IF (plot.log AND data_val EQ 0.0) THEN frac = 0.0
    
    color = cortex_SetColour(xrange,frac,lastct)

    IF (flip) THEN cortex_FlipArray, v
    POLYFILL, v[0,*], v[1,*], /DATA, COLOR=color, NOCLIP=0 
  ENDFOR
;
; Overlay the grid and wall:
; ----------------------------------------------------------------------
  IF (plot.option NE 900 AND plot.option NE 901) THEN  $
    status = cortex_PlotFluidGrid(plot, grid, wall, 0, annotate, 'overlay', 'outline', ps=ps)  $
  ELSE  $
    OPLOT, [0.9999*xmax,xmin,xmin,0.9999*xmax,0.9999*xmax], [ymin,ymin,0.9999*ymax,0.9999*ymax,ymin], COLOR=Truecolor('Black')
;
; Plot the scale: 
; ----------------------------------------------------------------------
  nsteps = 100.0

  IF (1 EQ 1 OR aspect_ratio LE 1.0) THEN BEGIN
    x0 = (0.95 * xpos[0] + 0.05 * xpos[1]) * dev_xsize
    x1 = (0.05 * xpos[0] + 0.95 * xpos[1]) * dev_xsize
    y0 = (ypos[0] - 0.11 * (ypos[1] - ypos[0])) * dev_ysize
    y1 = (ypos[0] - 0.09 * (ypos[1] - ypos[0])) * dev_ysize
;    y0 = (ypos[0] - 0.15 * (ypos[1] - ypos[0])) * dev_ysize
;    y1 = (ypos[0] - 0.13 * (ypos[1] - ypos[0])) * dev_ysize
;    y0 = (0.10 * ypos[0]                 ) * dev_ysize
;    y1 = (0.30 * ypos[0]                 ) * dev_ysize
    LOADCT, color_table
    lastct = -1
    FOR i = 0, nsteps-1 DO BEGIN
      frac1 = FLOAT(i  ) / nsteps
      frac2 = FLOAT(i+1) / nsteps
      x = [x0 + frac1 * (x1 - x0), x0 + frac2 * (x1 - x0),  $
           x0 + frac2 * (x1 - x0), x0 + frac1 * (x1 - x0)]
      y = [y0, y0, y1 ,y1]
      color = cortex_SetColour(xrange,0.5 * (frac1 + frac2),lastct)
      POLYFILL, x, y, /DEVICE, COLOR=color
    ENDFOR
;   This is nasty...
    log_prefix = ''
    IF (plot.log NE 0) THEN log_prefix = 'LOG10 '
    PLOT, xrange, [0.0,0.0], POSITION=[x0,y0,x1,y1], YTICKFORMAT='(A1)',  $
          TICKLEN=0, /NODATA, /DEVICE, /NOERASE,  $
          XTITLE=log_prefix+scale_label,COLOR=TrueColor('Black'), XSTYLE=1
  ENDIF ELSE BEGIN
    x0 = (xpos[1] + 0.09 * (xpos[1] - xpos[0])) * dev_xsize
    x1 = (xpos[1] + 0.11 * (xpos[1] - xpos[0])) * dev_xsize
    y0 = (0.95 * ypos[0] + 0.05 * ypos[1]) * dev_ysize
    y1 = (0.05 * ypos[0] + 0.95 * ypos[1]) * dev_ysize
    LOADCT, color_table
    lastct = -1
    FOR i = 0, nsteps-1 DO BEGIN
      frac1 = FLOAT(i  ) / nsteps
      frac2 = FLOAT(i+1) / nsteps
      x = [x0, x0, x1 ,x1]
      y = [y0 + frac1 * (y1 - y0), y0 + frac2 * (y1 - y0),  $
           y0 + frac2 * (y1 - y0), y0 + frac1 * (y1 - y0)]
      color = cortex_SetColour(xrange,0.5 * (frac1 + frac2),lastct)
      POLYFILL, x, y, /DEVICE, COLOR=color
    ENDFOR
    log_prefix = ''
    IF (plot.log NE 0) THEN log_prefix = 'LOG10 '
    PLOT, [0.0,0.0], xrange, POSITION=[x0,y0,x1,y1], XTICKFORMAT='(A1)',  $
          TICKLEN=0, /NODATA, /DEVICE, /NOERASE,  $
          YTITLE=log_prefix+scale_label,COLOR=TrueColor('Black'), YSTYLE=1
;    OPLOT, [x0,x0,x1,x1], [y0,y1,y1,y0], COLOR=TrueColor('Black')
  ENDELSE
;
;
; ----------------------------------------------------------------------  (PUT IN A FUNCTION)
  IF (plot.frame EQ 1) THEN BEGIN
    ypos = 0.940
    IF (notes NE 'default') THEN ypos = 0.965

    n = STRLEN(title) 
    nmax = LONG(60.0 * 2.0 / (1.5 * charsize))
    IF (n GT nmax) THEN BEGIN
      i = STRPOS(STRMID(title,0,nmax),' ',/REVERSE_SEARCH)
      str1 = STRMID(title,0  ,i   )
      str2 = STRMID(title,i+1,nmax)
      XYOUTS, 0.02 * dev_xsize, (ypos + 0.045 * charsize / 1.7) * dev_ysize, CHARSIZE=1.5*charsize, str1, /DEVICE
      XYOUTS, 0.02 * dev_xsize, (ypos                         ) * dev_ysize, CHARSIZE=1.5*charsize, str2, /DEVICE
    ENDIF ELSE BEGIN
      XYOUTS, 0.02 * dev_xsize, ypos * dev_ysize, CHARSIZE=1.5*charsize, title, /DEVICE
    ENDELSE

    IF (notes NE 'unknown') THEN BEGIN
      n = STRLEN(notes) 
      nmax = 105
      IF (n GT nmax) THEN BEGIN
        i = STRPOS(STRMID(notes,0,nmax),' ',/REVERSE_SEARCH)
        str1 = STRMID(notes,0,nmax)
        str2 = STRMID(notes,nmax+1,n)
        XYOUTS, 0.02 * dev_xsize, (ypos - 0.025) * dev_ysize, CHARSIZE=1.0, str1, /DEVICE
        XYOUTS, 0.02 * dev_xsize, (ypos - 0.050) * dev_ysize, CHARSIZE=1.0, str2, /DEVICE
        ypos = ypos - 0.020
      ENDIF ELSE  $
        XYOUTS, 0.02 * dev_xsize, (ypos - 0.025) * dev_ysize, CHARSIZE=1.0, notes, /DEVICE
    ENDIF
  ENDIF

  IF (NOT KEYWORD_SET(ps)) THEN BEGIN
    PRINT, 'STOPPING SO YOU CAN SEE THE PLOT'
    RETURN, -1
  ENDIF

  RETURN, 0
END
;
; ======================================================================
;


