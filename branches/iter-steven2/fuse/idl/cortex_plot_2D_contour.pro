;
; ======================================================================
;
FUNCTION cortex_GetVertex, grid, iobj, iside

  isrf = grid.obj_iside[iside,iobj]
  IF (isrf EQ 0) THEN BEGIN
    PRINT, 'ERROR cortex_GetVertex: Surface index is zero'
    STOP
  ENDIF

  i = WHERE(grid.srf_index EQ ABS(isrf))                              ; Should replace with a mapping array?
  IF (N_ELEMENTS(i) NE 1 OR i EQ -1) THEN BEGIN
    PRINT, 'ERROR cortex_GetVertex: Unable to map surface index'
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
FUNCTION cortex_Plot2DContour, plot, plot_data, grid, wall, ps=ps

  window_id = 0
  window_xsize = 700
  window_ysize = 700

  dev_xsize = !D.X_SIZE
  dev_ysize = !D.Y_SIZE
  !P.BACKGROUND = TrueColor('White')

  title    = plot.title
  notes    = plot.notes
  charsize = plot.charsize

  xtitle = 'R (m)'
  ytitle = 'Z (m)'
;
; Setup plot area:
; ----------------------------------------------------------------------
  IF (N_ELEMENTS(WHERE(plot.center EQ 0.0)) EQ 2) THEN BEGIN 
    size = 0.8
    xcen = 0.5 * TOTAL(plot.frame_bnds)
    ycen = 0.5
  ENDIF ELSE BEGIN
    size = 0.8
    IF (plot.size NE 0.0) THEN size = plot.size
    xcen = plot.center[0]
    ycen = plot.center[1]
  ENDELSE  
;
; Setup plot zoom:
; ----------------------------------------------------------------------
  IF (N_ELEMENTS(WHERE(plot.zoom EQ 0.0)) EQ 4) THEN BEGIN 
    xmin = MIN([[grid.vtx[0,*]],[wall.v1[0,*]],[wall.v2[0,*]]]) - 0.2
    xmax = MAX([[grid.vtx[0,*]],[wall.v1[0,*]],[wall.v2[0,*]]]) + 0.2
    ymin = MIN([[grid.vtx[1,*]],[wall.v1[1,*]],[wall.v2[1,*]]]) - 0.2
    ymax = MAX([[grid.vtx[1,*]],[wall.v1[1,*]],[wall.v2[1,*]]]) + 0.2
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
  IF (xpos[0] LT plot.frame_bnds[0] OR xpos[1] GT plot.frame_bnds[1]) THEN BEGIN
    ratio = 0.8 * (plot.frame_bnds[1] - plot.frame_bnds[0]) / (xpos[1] - xpos[0]) 
    print,'ratio',ratio
    xpos = [xcen - ratio * 0.5 * (xpos[1] - xpos[0]),  $
            xcen + ratio * 0.5 * (xpos[1] - xpos[0])]
    ypos = [ycen - ratio * 0.5 * (ypos[1] - ypos[0]),  $
            ycen + ratio * 0.5 * (ypos[1] - ypos[0])]
  ENDIF
;
; Store the plot information
; ----------------------------------------------------------------------
  plot.xrange   = [xmin,xmax]
  plot.yrange   = [ymin,ymax]
  plot.position = [xpos[0],ypos[0],xpos[1],ypos[1]]
;
; Draw the plot:
; ----------------------------------------------------------------------
;
  color_table = plot.color_table

  IF (plot.option GE   1 AND plot.option LE  99) THEN  $
    plot_labels = ['ELECTRON DENSITY          : n_e (m-3)         ' ,  $
                   'PARALLEL PLASMA FLOW      : v_parallel (m s-1)' ,  $
                   'PARALLEL FLOW MACH NUMBER : Mach no.          ' ,  $
                   'ELECTRON PRESSURE         : p_e (?)           ' ,  $
                   'ION PRESSURE              : p_i (?)           ' ,  $
                   'TOTAL PLASMA PRESSURE     : p (?)             ' ,  $
                   'ELECTRON TEMPERATURE      : T_e (eV)          ' ,  $
                   'ION TEMPERATURE           : T_i (eV)          '] 
  IF (plot.option GE 100 AND plot.option LE  199) THEN  $
    plot_labels = ['?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?']
  IF (plot.option GE 200 AND plot.option LE  299) THEN  $
    plot_labels = ['EIRENE ATOM DENSITY : n_D (m-3)' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?']
  IF (plot.option GE 300 AND plot.option LE  499) THEN  $
    plot_labels = ['EIRENE IMPURITY ATOM DENSITY: n_I+0 (m-3)' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?']
  IF (plot.option GE 400 AND plot.option LE  499) THEN  BEGIN
    state = 0
    plot_labels = ['DIVIMP IMPURITY DENSITY: n_I +'+STRTRIM(STRING(state),2)+' (m-3)',  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?' ,  $
                   '?:?']
  ENDIF

  CASE plot.option OF
;   ------------------------------------------------------------------
    1: data = plot_data.dens
    2: data = plot_data.vi
    3: data = plot_data.vi / plot_data.cs
    4: data = plot_data.pe
    5: data = plot_data.pi
    6: data = plot_data.pe + plot_data.pi
    7: data = plot_data.te
    8: data = plot_data.ti
;   ------------------------------------------------------------------
    200: data = plot_data.atm_dens
;   ------------------------------------------------------------------
    300: data = plot_data.imp_dens[*,0]  ; EIRENE
;   ------------------------------------------------------------------
    400: data = plot_data.imp_dens[*,0] * plot_data.div_influx   ; DIVIMP
;   ------------------------------------------------------------------
    ELSE: BEGIN
      PRINT, 'ERROR cortex_Plot2D: Unrecognised option'
      PRINT, '  PLOT OPTION = ',plot.option
      RETURN, -1
      END
  ENDCASE

;  *** FIND A BETTER WAY TO SORT THROUGH THESE ***
  plot_index = plot.option
  IF (plot.option GE 100 AND plot.option LE 199) THEN plot_index = plot_index -  99
  IF (plot.option GE 200 AND plot.option LE 299) THEN plot_index = plot_index - 199
  IF (plot.option GE 300 AND plot.option LE 399) THEN plot_index = plot_index - 299
  IF (plot.option GE 400 AND plot.option LE 499) THEN plot_index = plot_index - 399

  labels = STRSPLIT(plot_labels[plot_index-1],':',/EXTRACT)
  IF (plot.plot_title NE 'unknown') THEN plot_title = plot.plot_title ELSE   $
                                         plot_title = STRTRIM(labels[0],2)
  scale_label = STRTRIM(labels[1],2)


  IF (plot.log) THEN FOR i = 0, N_ELEMENTS(data)-1 DO  $                    ; A nicer way to do this..?
    IF (data[i] GT 0.0) THEN data[i] = ALOG10(data[i]) ELSE data[i] = 0.0

  xrange = [MIN(data),MAX(data)]

  PLOT, plot.xrange, plot.yrange, /NODATA, XSTYLE=1, YSTYLE=1, /NOERASE,  $
        POSITION=plot.position,                                           $
        TITLE=plot_title, XTITLE=xtitle, YTITLE=ytitle,                   $
        COLOR=TrueColor('Black')

  LOADCT, color_table

  IF (N_ELEMENTS(data) NE grid.obj_n) THEN BEGIN
    PRINT
    PRINT, '------------------------------------------'
    PRINT, '  WARNING cortex_Plot2D: NDATA NE NOBJ'
    PRINT, '------------------------------------------'
    PRINT
  ENDIF

print,'xrange=',xrange[0:1]

  FOR iobj = 0, grid.obj_n-1 DO BEGIN
;
;   Filter(s):
;    IF (grid.obj_tube[iobj] LT 13 OR grid.obj_tube[iobj] GT 13) THEN CONTINUE
;    IF (grid.obj_tube[iobj] LT 7 OR grid.obj_tube[iobj] GT 7) THEN CONTINUE
;    IF (grid.obj_tube[iobj] NE 8) THEN CONTINUE
;
    v = MAKE_ARRAY(3,grid.obj_nside[iobj],/FLOAT,VALUE=0.0)

    FOR iside = 0, grid.obj_nside[iobj]-1 DO  $
      v[*,iside] = cortex_GetVertex(grid,iobj,iside)

    frac = (data[iobj] - xrange[0]) / (xrange[1] - xrange[0])
    color = LONG(frac * 255.0)
    POLYFILL, v[0,*], v[1,*], /DATA, COLOR=color, NOCLIP=0 
  ENDFOR
print,plot.xrange[0:1],plot.yrange[0:1]
;
; Overlay the grid and wall:
; ----------------------------------------------------------------------
  status = cortex_PlotFluidGrid(plot, grid, wall, 0, 0, 'overlay', 'outline', ps=ps)
;
; Plot the scale: 
; ----------------------------------------------------------------------
  nsteps = 100.0


print,'ypos',ypos

  IF (1) THEN BEGIN
    x0 = (0.95 * xpos[0] + 0.05 * xpos[1]) * dev_xsize
    x1 = (0.05 * xpos[0] + 0.95 * xpos[1]) * dev_xsize
    y0 = (ypos[0] - 0.12 * (ypos[1] - ypos[0])) * dev_ysize
    y1 = (ypos[0] - 0.10 * (ypos[1] - ypos[0])) * dev_ysize
;    y0 = (ypos[0] - 0.15 * (ypos[1] - ypos[0])) * dev_ysize
;    y1 = (ypos[0] - 0.13 * (ypos[1] - ypos[0])) * dev_ysize

;    y0 = (0.10 * ypos[0]                 ) * dev_ysize
;    y1 = (0.30 * ypos[0]                 ) * dev_ysize
    LOADCT, color_table
    FOR i = 0, nsteps-1 DO BEGIN
      frac1 = FLOAT(i  ) / nsteps
      frac2 = FLOAT(i+1) / nsteps
      x = [x0 + frac1 * (x1 - x0), x0 + frac2 * (x1 - x0),  $
           x0 + frac2 * (x1 - x0), x0 + frac1 * (x1 - x0)]
      y = [y0, y0, y1 ,y1]
      color = LONG( 0.5 * (frac1 + frac2) * 255.0) 
      POLYFILL, x, y, /DEVICE, COLOR=color
    ENDFOR
;   This is nasty...
    log_prefix = ''
    IF (plot.log NE 0) THEN log_prefix = 'LOG '
    PLOT, xrange, [0.0,0.0], POSITION=[x0,y0,x1,y1], YTICKFORMAT='(A1)',  $
          TICKLEN=0, /NODATA, /DEVICE, /NOERASE,  $
          XTITLE=log_prefix+scale_label,COLOR=TrueColor('Black')
  ENDIF ELSE BEGIN
    x0 = 0.10 * xmin + 0.90 * xmax
    x1 = 0.05 * xmin + 0.95 * xmax
    y0 = 0.70 * ymin + 0.30 * ymax
    y1 = 0.05 * ymin + 0.95 * ymax
    LOADCT, color_table
    FOR i = 0, nsteps-1 DO BEGIN
      frac1 = FLOAT(i  ) / nsteps
      frac2 = FLOAT(i+1) / nsteps
      x = [x0, x0, x1 ,x1]
      y = [y0 + frac1 * (y1 - y0), y0 + frac2 * (y1 - y0),  $
           y0 + frac2 * (y1 - y0), y0 + frac1 * (y1 - y0)]
      color = LONG( 0.5 * (frac1 + frac2) * 255.0) 
print,color
      POLYFILL, x, y, /DATA, COLOR=color
    ENDFOR
    OPLOT, [x0,x0,x1,x1], [y0,y1,y1,y0], COLOR=TrueColor('Black')
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


