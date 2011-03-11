;
; ======================================================================
;
PRO cortex_FlipPoints, v1, v2

  v1_save = v1
  v2_save = v2
  v1[0] = v1_save[1]
  v1[1] = v1_save[0]
  v2[0] = v2_save[1]
  v2[1] = v2_save[0]

END
;
; ======================================================================
;
PRO cortex_FlipArray, p

  p_save = p

  p[0,*] = p_save[1,*]
  p[1,*] = p_save[0,*]

END
;
; ======================================================================
;
FUNCTION cortex_PlotFluidGrid, plot, grid, wall, node, annotation, mode, type, ps=ps

; A better system than the one used here is to collect all the data at the beginning of the 
; routine, and then find XMIN, etc. since then the data is only being processed once.  It also helps
; with storing the data later for redrawing, a bit... this scheme should be implemented for
; most of the plots... for another day...

  flip = 0 ; 1

  IF (NOT KEYWORD_SET(mode)) THEN mode = 'main'
  IF (NOT KEYWORD_SET(type)) THEN type = 'full'

  IF (plot.outline         ) THEN type = 'outline'
  IF (plot.flux_surfaces   ) THEN type = 'flux_surfaces'
  IF (plot.no_grid         ) THEN type = 'no_grid'
  IF (plot.equ NE 'default') THEN type = 'equ'
  IF (plot.xticks NE -1    ) THEN xticks = plot.xticks
  IF (plot.yticks NE -1    ) THEN yticks = plot.yticks
  IF (plot.xminor NE -1    ) THEN xminor = plot.xminor
  IF (plot.yminor NE -1    ) THEN yminor = plot.yminor

  window_id = 0
  window_xsize = 700
  window_ysize = 700

  !P.CHARSIZE = plot.charsize
  !P.CHARTHICK = plot.thick
  !P.THICK    = plot.thick
  !X.THICK    = plot.thick
  !Y.THICK    = plot.thick
  !Z.THICK    = plot.thick

  !P.BACKGROUND = TrueColor('White')

  IF (mode EQ 'subordinate' OR mode EQ 'overlay') THEN BEGIN
    plot_title  = ' '
  ENDIF ELSE BEGIN
    plot_title  = ' '
    IF (type EQ 'equ') THEN plot_title  = 'EQUILIBRIUM PLOT (DG/CARRE FORMAT)'
  ENDELSE
  plot_xtitle = 'R (m)'
  plot_ytitle = 'Z (m)'
  IF (flip) THEN BEGIN
    plot_xtitle = 'distance parallel to field, s (m)'
    plot_ytitle = 'radial distance (m)'
  ENDIF

  colors = ['Red','Green','Blue']
;
; Setup plot zoom:
; ----------------------------------------------------------------------
  IF (N_ELEMENTS(WHERE(plot.zoom EQ 0.0)) EQ 4) THEN BEGIN 
    xmin =  1.0E+20
    xmax = -1.0E+20
    ymin =  1.0E+20
    ymax = -1.0E+20
    CASE type OF
;     ------------------------------------------------------------------
      'no_grid': 
;     ------------------------------------------------------------------
      'equ': BEGIN
        xmin = MIN([xmin,MIN(grid.x[*])])
        xmax = MAX([xmax,MAX(grid.x[*])])
        ymin = MIN([ymin,MIN(grid.y[*])])
        ymax = MAX([ymax,MAX(grid.y[*])])
        END
;     ------------------------------------------------------------------
      ELSE: BEGIN      
        xmin = MIN([xmin,MIN(grid.vtx[0,*])])
        xmax = MAX([xmax,MAX(grid.vtx[0,*])])
        ymin = MIN([ymin,MIN(grid.vtx[1,*])])
        ymax = MAX([ymax,MAX(grid.vtx[1,*])])
        END
    ENDCASE
    IF (NOT plot.no_wall) THEN BEGIN
      xmin = MIN([xmin,MIN(wall.v1[0,*]),MIN(wall.v2[0,*])])
      xmax = MAX([xmax,MAX(wall.v1[0,*]),MAX(wall.v2[0,*])])
      ymin = MIN([ymin,MIN(wall.v1[1,*]),MIN(wall.v2[1,*])])
      ymax = MAX([ymax,MAX(wall.v1[1,*]),MAX(wall.v2[1,*])])
    ENDIF
    IF (plot.annotate_n GT 0) THEN BEGIN
      FOR i = 0, plot.annotate_n-1 DO BEGIN
        val = cortex_ExtractStructure(annotation,i+1)   
        CASE plot.annotate_code[i] OF
;         ----------------------------------------------------------------
          1: BEGIN
            xmin = MIN([xmin,MIN(val.x[*])])
            xmax = MAX([xmax,MAX(val.x[*])])
            ymin = MIN([ymin,MIN(val.y[*])])
            ymax = MAX([ymax,MAX(val.y[*])])
            END
;         ----------------------------------------------------------------
          2: 
;         ----------------------------------------------------------------
          ELSE: BEGIN
            PRINT, 'ERROR cortex_PlotFluidGrid: Unrecognised annotation code'
            PRINT, '  TAG  = ',plot.tag
            PRINT, '  INDEX= ',i
            PRINT, '  TYPE = ',plot.annotate_code[i]
            RETURN, -1
            END           
        ENDCASE          
      ENDFOR
    ENDIF
    xdelta = xmax - xmin
    ydelta = ymax - ymin
    boarder = MIN([xdelta,ydelta]) * 0.05
    xmin = xmin - boarder
    xmax = xmax + boarder
    ymin = ymin - boarder
    ymax = ymax + boarder
  ENDIF ELSE BEGIN
    xmin = plot.zoom[0]
    xmax = plot.zoom[2]
    ymin = plot.zoom[1]
    ymax = plot.zoom[3]
  ENDELSE
  PRINT, 'XMIN,MAX=',xmin,xmax
  PRINT, 'YMIN,MAX=',ymin,ymax

  xdelta = xmax - xmin
  ydelta = ymax - ymin
;
; Setup plot area:
; ----------------------------------------------------------------------
;
  IF (KEYWORD_SET(ps)) THEN BEGIN
    xsize = 297.0  ; Landscape A4
    ysize = 210.0  
    display_ratio = 3.25 / 3.35  ; Small correction based on measurements from the GhostView display of the plot
    aspect_ratio = (xdelta / xsize) / (ydelta / ysize) * display_ratio
    plot_xboarder = 0.05
  ENDIF ELSE BEGIN
    PRINT, 'NEED TO FIX WINDOW ASPECT RATIO'
    RETURN, -1
    display_ratio = 4.4
    aspect_ratio = xdelta / ydelta / display_ratio *            $ 
                   (FLOAT(window_ysize) / FLOAT(window_xsize))
    plot_xboarder = 0.05
    WINDOW, window_id, XSIZE=window_xsize, YSIZE=window_ysize
  ENDELSE
  IF (plot.aspect_ratio GT 0.0) THEN aspect_ratio = plot.aspect_ratio

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

  IF ((xdelta / xsize) LT (ydelta / ysize)) THEN BEGIN
    IF (N_ELEMENTS(WHERE(plot.center EQ 0.0)) EQ 2 AND  $
        (plot.frame_bnds[0] NE 0.0 OR plot.frame_bnds[1] NE 1.0)) THEN BEGIN 
      plot.frame_bnds[1] = plot.frame_bnds[0] + size * aspect_ratio + plot_xboarder
      xcen = 0.5 * TOTAL(plot.frame_bnds)
    ENDIF 
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
; ----------------------------------------------------------------------
  IF (xpos[0] LT plot.frame_bnds[0] OR xpos[1] GT plot.frame_bnds[1]) THEN BEGIN
    ratio = size * (plot.frame_bnds[1] - plot.frame_bnds[0]) / (xpos[1] - xpos[0]) 
    print,'ratio',ratio
    xpos = [xcen - ratio * 0.5 * (xpos[1] - xpos[0]),  $
            xcen + ratio * 0.5 * (xpos[1] - xpos[0])]
    ypos = [ycen - ratio * 0.5 * (ypos[1] - ypos[0]),  $
            ycen + ratio * 0.5 * (ypos[1] - ypos[0])]
  ENDIF

  IF (mode EQ 'overlay') THEN BEGIN
;   Do nothing since the axes have been set somewhere else:
  ENDIF ELSE BEGIN
;
;   Store the plot information
;   ----------------------------------------------------------------------
;    plot.xrange   = [xmin,xmax]
;    plot.yrange   = [ymin,ymax]
    plot.zoom = [xmin,ymin,xmax,ymax]
    xrange = [xmin,xmax]
    yrange = [ymin,ymax]
    plot.position = [xpos[0],ypos[0],xpos[1],ypos[1]]
;
;   Create the axes:
;   ----------------------------------------------------------------------
    IF (flip) THEN BEGIN
      xrange_save = xrange
      xrange = yrange
      yrange = xrange_save
    ENDIF

    IF (type NE 'equ') THEN  $
      PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1, /NOERASE,    $
            POSITION=plot.position,                                    $
            TITLE=plot_title, XTITLE=plot_xtitle, YTITLE=plot_ytitle,  $
            COLOR=TrueColor('Black'),  $
            XTICKS=xticks,YTICKS=yticks,XMINOR=xminor,YMINOR=yminor
  ENDELSE
;
; Plot the grid:
; ----------------------------------------------------------------------
  CASE (type) OF
;   ----------------------------------------------------------------
    'no_grid': 
;   ----------------------------------------------------------------
    'equ': BEGIN
;      levels = 0.10 * (FINDGEN(25) / 25.0 - 0.5) + 1.1
      IF (plot.equ_levels[0] NE -999.0) THEN  $
        levels = plot.equ_levels[WHERE(plot.equ_levels NE -999.0)]  $
      ELSE  $
        levels = plot.equ_params[1] *   $
                 (FINDGEN(plot.equ_params[0]) / plot.equ_params[0] - 0.5) +  $
                 plot.equ_params[2]

;      levels = (plot.equ_params[2] - plot.equ_params[1]) *   $
;               (FINDGEN(plot.equ_params[0]) / plot.equ_params[0] - 0.5) +  $
;               plot.equ_params[1]
print,levels

      CONTOUR, grid.psin, grid.x, grid.y, LEVELS=levels, XSTYLE=1, YSTYLE=1, /NOERASE, $
               POSITION=plot.position,  XRANGE=xrange,YRANGE=yrange,  $
               TITLE=plot_title, XTITLE=plot_xtitle, YTITLE=plot_ytitle,  $
               COLOR=TrueColor('Black'),  $
               XTICKS=xticks,YTICKS=yticks,XMINOR=xminor,YMINOR=yminor
      END
;   ----------------------------------------------------------------
    ELSE: BEGIN
      FOR iobj = 0, grid.obj_n-1 DO BEGIN
;     
;       Filter(s):
;        IF (grid.obj_tube[iobj] LT 13 OR grid.obj_tube[iobj] GT 13) THEN CONTINUE
;        IF (grid.obj_tube[iobj] LT 7 OR grid.obj_tube[iobj] GT 7) THEN CONTINUE
;        IF (grid.obj_tube[iobj] NE 8) THEN CONTINUE
;     
        FOR iside = 0, grid.obj_nside[iobj]-1 DO BEGIN
          color = 'Black'
;          IF (plot.show_n GT 0) THEN  $
;            FOR i = 0, plot.show_n-1 DO IF (plot.show_tube[i] EQ grid.obj_tube[iobj]) THEN  $
;              color = plot.show_colour[i]
          CASE (type) OF
;           --------------------------------------------------------
            'full': BEGIN
              isrf = grid.obj_iside[iside,iobj]
              IF (isrf GT 0) THEN BEGIN
                i = WHERE(grid.srf_index EQ isrf)                              ; Should replace with a mapping array?
                IF (N_ELEMENTS(i) NE 1 OR i EQ -1) THEN BEGIN
                  PRINT, 'ERROR PlotFluidGrid: Surface index not defined'
                  RETURN, -1
                ENDIF
                j  = grid.vtx_map[grid.srf_ivtx[0,i]]
                k  = grid.vtx_map[grid.srf_ivtx[1,i]]
                v1 = grid.vtx[*,j]
                v2 = grid.vtx[*,k]
                IF (flip) THEN cortex_FlipPoints, v1, v2
                OPLOT, [v1[0],v2[0]], [v1[1],v2[1]], COLOR=TrueColor(color)
              ENDIF
              END
;           --------------------------------------------------------
            'flux_surfaces': BEGIN
              isrf = grid.obj_iside[iside,iobj]
              IF (iside EQ 1 OR (iside EQ 3 AND grid.obj_omap[3,iobj] EQ -1)) THEN BEGIN
                i = WHERE(grid.srf_index EQ ABS(isrf))                              ; Should replace with a mapping array?
                IF (N_ELEMENTS(i) NE 1 OR i EQ -1) THEN BEGIN
                  PRINT, 'ERROR PlotFluidGrid: Surface index not defined'
                  RETURN, -1
                ENDIF
                j  = grid.vtx_map[grid.srf_ivtx[0,i]]
                k  = grid.vtx_map[grid.srf_ivtx[1,i]]
                v1 = grid.vtx[*,j]
                v2 = grid.vtx[*,k]
                IF (flip) THEN cortex_FlipPoints, v1, v2
                OPLOT, [v1[0],v2[0]], [v1[1],v2[1]], COLOR=TrueColor(color)
              ENDIF
              END
;           --------------------------------------------------------
            'outline': BEGIN
              IF (grid.obj_omap[iside,iobj] LE -1) THEN BEGIN
                isrf = ABS(grid.obj_iside[iside,iobj])
                i  = WHERE(grid.srf_index EQ isrf)                              ; Should replace with a mapping array?
                j  = grid.vtx_map[grid.srf_ivtx[0,i]]
                k  = grid.vtx_map[grid.srf_ivtx[1,i]]
                v1 = grid.vtx[*,j]
                v2 = grid.vtx[*,k]
                IF (flip) THEN cortex_FlipPoints, v1, v2
                OPLOT, [v1[0],v2[0]], [v1[1],v2[1]], COLOR=TrueColor('Black')
              ENDIF
              END
;           --------------------------------------------------------
            ELSE: BEGIN
              PRINT, 'ERROR cortex_PlotFluidGrid: Unrecognised TYPE'
              PRINT, '  TYPE = ',type
              RETURN, -1
              END         
          ENDCASE
        ENDFOR
      ENDFOR
    END
  ENDCASE
;
; Overlay the separatrix:
; ----------------------------------------------------------------------
  IF (type NE 'equ' and NOT plot.no_separatrix) THEN BEGIN

    FOR i1 = -1, plot.show_n-1 DO BEGIN

      IF (i1 EQ -1) THEN BEGIN
        tube  = grid.isep
        color = 'Red'
        side = 3
        iv0 = 0
        iv1 = 1
      ENDIF ELSE BEGIN
        tube  = plot.show_tube  [i1]
        color = plot.show_colour[i1]
        side = 1
        iv0 = 1
        iv1 = 0
      ENDELSE
    
      i = WHERE(grid.obj_tube EQ tube, n)
      IF (n EQ 0) THEN CONTINUE
      v = MAKE_ARRAY(2,n+1,/FLOAT,VALUE=0.0)
      FOR j = 0, n-1 DO BEGIN
        isrf = grid.obj_iside[side,i[j]]
        k = WHERE(grid.srf_index EQ ABS(isrf))
        IF (N_ELEMENTS(k) NE 1 OR k EQ -1) THEN BEGIN
          PRINT, 'ERROR PlotFluidGrid: Surface index not defined for sepatartrix'
          RETURN, -1
        ENDIF
        IF (isrf GT 0) THEN m = grid.vtx_map[grid.srf_ivtx[iv1,k]] ELSE  $
                            m = grid.vtx_map[grid.srf_ivtx[iv0,k]]
        v[0:1,j] = grid.vtx[0:1,m]
        IF (j EQ n-1) THEN v[0:1,j+1] = grid.vtx[0:1,m] 
      ENDFOR
      IF (flip) THEN cortex_FlipArray, v
      OPLOT, v[0,*], v[1,*], COLOR=TrueColor(color)
    ENDFOR
  ENDIF
;
; Show interpolation nodes:
; ----------------------------------------------------------------------
  IF (plot.nodes) THEN BEGIN
    FOR i = 0, N_ELEMENTS(node.type)-2 DO BEGIN
      IF (node.type[i] EQ 0.0 OR node.type[i] NE node.type[i+1]) THEN CONTINUE
      v1 = [node.x[i  ],node.y[i  ]]
      v2 = [node.x[i+1],node.y[i+1]]
      IF (flip) THEN cortex_FlipPoints, v1, v2      
      OPLOT, [v1[0],v2[0]], [v1[1],v2[1]], THICK=4.0, COLOR=TrueColor('Green')
    ENDFOR
  ENDIF
;
; Plot the wall:
; ----------------------------------------------------------------------
  IF (NOT plot.no_wall) THEN BEGIN
    FOR iwall = 0, wall.n-1 DO BEGIN
      IF (wall.class[iwall] NE 1) THEN CONTINUE
;      IF (wall.class[iwall] NE 1 OR wall.target[iwall] NE 0) THEN CONTINUE
      v1 = wall.v1[*,iwall]
      v2 = wall.v2[*,iwall]
      IF (flip) THEN cortex_FlipPoints, v1, v2
      OPLOT, [v1[0],v2[0]], [v1[1],v2[1]], COLOR=TrueColor('Blue')
;
;     Label segments:
;      IF (N_ELEMENTS(WHERE(plot.label_index EQ 0)) NE 2 AND  $
;          iwall GE plot.label_index[0] AND  $                    ; Shouldn't use indexing based on the array index
;          iwall LE plot.label_index[1]) THEN BEGIN               ; Replace with a wall index array in the structure
    ENDFOR
  ENDIF

  IF (plot.label_index NE 'none') THEN BEGIN
    FOR iwall = 0, wall.n-1 DO BEGIN
      IF (wall.class[iwall] NE 1) THEN CONTINUE
      IF (cortex_CheckIndex(iwall+1,plot.label_index)) THEN BEGIN   ; This IWALL+1 business is poor
        v1 = wall.v1[*,iwall]
        v2 = wall.v2[*,iwall]
        vlabel = 0.5 * (v1 + v2)
        IF (flip) THEN cortex_FlipPoints, vlabel, v2      
        IF (plot.label_index NE 'all') THEN BEGIN
          OPLOT, [vlabel[0]], [vlabel[1]], COLOR=TrueColor('Black'),  $
                 PSYM=6, SYMSIZE=0.4
          charsize = 1.0
print, 'vlable!',vlabel[0],xmin,xmax,xmin+0.80*(xmax-xmin)
          IF (vlabel[0] GT xmin+0.80*(xmax-xmin)) THEN  $
            alignment = 1.0 ELSE  $
            alignment = 0.0            
          XYOUTS, vlabel[0], vlabel[1], CHARSIZE=charsize, ALIGNMENT=alignment, $
                  ' '+STRTRIM(STRING(iwall+1),2)+' ', COLOR=TrueColor('Black')   ; This IWALL+1 business is poor
        ENDIF ELSE BEGIN
          charsize = 0.5
          XYOUTS, vlabel[0], vlabel[1], CHARSIZE=charsize, ALIGNMENT=0.5, $
                  STRTRIM(STRING(iwall+1),2), COLOR=TrueColor('Black')        ; This IWALL+1 business is poor
        ENDELSE
      ENDIF           
    ENDFOR
  ENDIF
;
; Annotate the grid (add pretty stuff):
; ----------------------------------------------------------------------
  IF (plot.annotate_n GT 0) THEN BEGIN
;help,plot.annotate,/struct
;print,shit
;stop
    FOR i = 0, plot.annotate_n-1 DO BEGIN
; check if the data is there... somehow...
      val = cortex_ExtractStructure(annotation,i+1)   
      CASE plot.annotate_code[i] OF
;       ----------------------------------------------------------------
        1: BEGIN
          OPLOT, val.x, val.y, COLOR=TrueColor(plot.annotate_colour[i])
          FOR j = 0, N_ELEMENTS(val.x)-2 DO BEGIN
            x = [val.x[j  ],val.y[j  ]]
            y = [val.x[j+1],val.y[j+1]]
          ENDFOR
          END

;       ----------------------------------------------------------------
        2: BEGIN
          FOR j = 0, N_ELEMENTS(val.x1)-1 DO BEGIN
            IF ((j+1) MOD plot.annotate_step[i] NE 0 AND  $
                j NE 0 AND j NE N_ELEMENTS(val.x1)-1) THEN CONTINUE
            x = [val.x1[j],val.x2[j]]
            y = [val.y1[j],val.y2[j]]
            OPLOT, x, y, COLOR=TrueColor(plot.annotate_colour[i])
            XYOUTS, x[1], y[1], CHARSIZE=charsize, ALIGNMENT=1.0, $
                    ' '+STRTRIM(STRING(j+1),2), COLOR=TrueColor(plot.annotate_colour[i])
          ENDFOR
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN
          PRINT, 'ERROR cortex_PlotFluidGrid: Unrecognised annotation code'
          PRINT, '  TAG  = ',plot.tag
          PRINT, '  INDEX= ',i
          PRINT, '  TYPE = ',plot.annotate_code[i]
          RETURN, -1
          END           
      ENDCASE
    ENDFOR
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


