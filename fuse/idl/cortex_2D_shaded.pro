;
; ======================================================================
;
FUNCTION cortex_Plot2D, plot, grid, wall, plasma, source, eirene, ps=ps

  window_id = 0
  window_xsize = 700
  window_ysize = 700

  !P.BACKGROUND = TrueColor('White')

  xtitle = 'R (m)'
  ytitle = 'Z (m)'
;
; Setup plot area:
; ----------------------------------------------------------------------
  IF (N_ELEMENTS(WHERE(plot.center EQ 0.0)) EQ 2) THEN BEGIN 
    size = 0.8
    xcen = 0.5
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
; Store the plot information
; ----------------------------------------------------------------------
  plot.xrange   = [xmin,xmax]
  plot.yrange   = [ymin,ymax]
  plot.position = [xpos[0],ypos[0],xpos[1],ypos[1]]
;
; Draw the plot:
; ----------------------------------------------------------------------
;
  color_table = 3

  plot_title = ['ELECTRON DENSITY'         ,  $
                'PARALLEL PLASMA FLOW'     ,  $
                'PARALLEL FLOW MACH NUMBER',  $
                'ELECTRON PRESSURE'        ,  $
                'ION PRESSURE'             ,  $
                'TOTAL PLASMA PRESSURE'    ,  $
                'ELECTRON TEMPERATURE'     ,  $
                'ION TEMPERATURE'] 

  CASE plot.option OF
    1: data = plasma.dens
    2: data = plasma.vi
    3: data = plasma.vi / plasma.cs
    4: data = plasma.pe
    5: data = plasma.pi
    6: data = plasma.pe + plasma.pi
    7: data = plasma.te
    8: data = plasma.ti
;   ------------------------------------------------------------------
    ELSE: BEGIN
      PRINT, 'ERROR cortex_Plot2D: Unrecognised option'
      PRINT, '  PLOT OPTION = ',plot.option
      RETURN, -1
      END
  ENDCASE

  IF (plot.xlog NE 0) THEN data = ALOG10(data>0.0)
;  print,data

  xrange = [MIN(data),MAX(data)]
  PRINT,'XRANGE=',xrange

  PLOT, plot.xrange, plot.yrange, /NODATA, XSTYLE=1, YSTYLE=1,     $
        POSITION=plot.position,                                    $
        TITLE=plot_title[plot.option-1], XTITLE=xtitle, YTITLE=ytitle,  $
        COLOR=TrueColor('Black')

  LOADCT, color_table

  IF (N_ELEMENTS(data) NE grid.obj_n) THEN BEGIN
    PRINT
    PRINT, '------------------------------------------'
    PRINT, '  WARNING cortex_Plot2D: NDATA NE NOBJ'
    PRINT, '------------------------------------------'
    PRINT
  ENDIF

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
    POLYFILL, v[0,*], v[1,*], /DATA, COLOR=color
  ENDFOR
;
; Overlay the grid and wall:
; ----------------------------------------------------------------------
  cortex_PlotFluidGrid, plot, grid, wall, 'overlay', 'outline', ps=ps
;
; Plot the scale: 
; ----------------------------------------------------------------------
  nsteps = 100.0

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
    POLYFILL, x, y, /DATA, COLOR=color
  ENDFOR
  OPLOT, [x0,x0,x1,x1], [y0,y1,y1,y0], COLOR=TrueColor('Black')


  IF (NOT KEYWORD_SET(ps)) THEN BEGIN
    PRINT, 'STOPPING SO YOU CAN SEE THE PLOT'
    RETURN, -1
  ENDIF

  RETURN, 0

END
;
; ======================================================================
;


