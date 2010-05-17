FUNCTION cortex_SliceGrid, line, grid, status
;
; Slide the grid with the nasty line segment:
; ----------------------------------------------------------------------
;
  a1 = line[0:1]
  a2 = line[2:3]

 print,a1
 print,a2
;
; Find intersection with the separatrix, if any:
;
 

;
; Find intersections between the line segment and the flux-tubes:
;
  FOR iobj = 0, grid.obj_n-1 DO BEGIN

;   if (grid.obj_tube[iobj] ne 1) then RETURN, 1
;
;   Setup the line segment subtending the current grid cell:
;
    isrf = ABS(grid.obj_iside[0,iobj])
    i = WHERE(grid.srf_index EQ isrf)                              ; Should replace with a mapping array?
    IF (N_ELEMENTS(i) NE 1 OR i EQ -1) THEN BEGIN
      PRINT, 'ERROR cortex_SliceGrid: Surface index not defined'
      STOP
    ENDIF
    j  = grid.vtx_map[grid.srf_ivtx[0,i]]
    k  = grid.vtx_map[grid.srf_ivtx[1,i]]
    b1 = 0.5 * (grid.vtx[0:1,j] + grid.vtx[0:1,k])

    isrf = ABS(grid.obj_iside[2,iobj])
    i = WHERE(grid.srf_index EQ isrf)                              ; Should replace with a mapping array?
    IF (N_ELEMENTS(i) NE 1 OR i EQ -1) THEN BEGIN
      PRINT, 'ERROR cortex_SliceGrid: Surface index not defined'
      STOP
    ENDIF
    j  = grid.vtx_map[grid.srf_ivtx[0,i]]
    k  = grid.vtx_map[grid.srf_ivtx[1,i]]
    b2 = 0.5 * (grid.vtx[0:1,j] + grid.vtx[0:1,k])

    Lint, a1, a2, b1, b2, c, d, flag=flag

    IF (flag NE 1) THEN CONTINUE

    s = (c[0] - a1[0]) / (a2[0] - a1[0])    
    t = (c[0] - b1[0]) / (b2[0] - b1[0])    

    IF ( s GE 0.0 AND s LE 1.0 AND  $
         t GE 0.0 AND t LT 1.0) THEN BEGIN

      i = grid.obj_tube[iobj] - 1

      IF (NOT KEYWORD_SET(int_cell)) THEN BEGIN
        int_cell = grid.obj_cell[iobj]
        int_dist = s
        int_psin = grid.tube_psin[i]
        int_rho  = grid.tube_rho [i]
        int_l    = grid.tube_l   [i]
        int_t    = t      
      ENDIF ELSE BEGIN
        int_cell = [int_cell,grid.obj_cell[iobj]]
        int_dist = [int_dist,s]
        int_psin = [int_psin,grid.tube_psin[i]]
        int_rho  = [int_rho ,grid.tube_rho [i]]
        int_l    = [int_l   ,grid.tube_l   [i]]
        int_t    = [int_t   ,t]      
      ENDELSE

;       print,iobj
;       print,b1
;       print,b2

;       print,s,t
;       print,c,d
    ENDIF

  ENDFOR

  IF (NOT KEYWORD_SET(int_cell)) THEN BEGIN
;   No intersection found between the line segment and the grid:
    status = -1
    result = -1
  ENDIF ELSE BEGIN
    status = 0
    result = {          $
      cell : int_cell,  $
      dist : int_dist,  $
      psin : int_psin,  $
      rho  : int_rho ,  $
      l    : int_l   ,  $
      t    : int_t   }          
  ENDELSE

  RETURN,result
END
;
; ======================================================================
;
FUNCTION cortex_PlotRadialProfile, plot, plot_array, ps=ps

  ndata = N_ELEMENTS(TAG_NAMES(plot_array))
  IF (ndata LE 0) THEN BEGIN
    PRINT, 'ERROR Cortex PlotMidplaneProfiles: No data found'
    RETURN, -1
  ENDIF
  PRINT,'NDATA= ',ndata

  focus = plot.focus

  charsize = plot.charsize

  !P.CHARSIZE = charsize
  !P.CHARTHICK = plot.thick
  !P.THICK    = plot.thick
  !X.THICK    = plot.thick
  !Y.THICK    = plot.thick
  !Z.THICK    = plot.thick

  !P.BACKGROUND = TrueColor('White')

  plot_title  = plot.title
  plot_xtitle = 'psi_n'
  plot_ytitle = ['1','n_e (m-3)','3','v_parallel / M','5','T_e (eV)']

  colors = ['Black',   'Red','Green','Blue','Orange','Purple', 'Hotpink', 'Darkseagreen', 'Silver',  $
            'Darkred', 'Greenyellow']
;  colors = ['Red','Green','Blue']
;
; Plot 2D grid with the line segment for the radial profile:
; ----------------------------------------------------------------------
;
  plot_grid = plot
  IF (focus NE 0) THEN BEGIN
    plot_grid.size   = 0.5
    plot_grid.center = [0.88,0.62]
  ENDIF ELSE BEGIN
    plot_grid.size   = 0.5
    plot_grid.center = [0.11,0.50]
  ENDELSE
  plot_data = cortex_ExtractStructure(plot_array,1)
  grid = plot_data.grid
  wall = plot_data.wall
  status = cortex_PlotFluidGrid(plot_grid, grid, wall, 0, 0, 'subordinate', 'outline', ps='on')
;
; Line segment annotation on the grid plot:
; ----------------------------------------------------------------------
;
  FOR i = 1, plot_grid.line_seg_n DO BEGIN
    v1 = plot.line_seg[0+(i-1)*4:1+(i-1)*4]
    v2 = plot.line_seg[2+(i-1)*4:3+(i-1)*4]
    xrange = [plot_grid.zoom[0],plot_grid.zoom[2]]
    yrange = [plot_grid.zoom[1],plot_grid.zoom[3]]
    PLOT, [v1[0],v2[0]], [v1[1],v2[1]], XSTYLE=5, YSTYLE=5, /NOERASE,  $
          XRANGE=xrange, YRANGE=yrange,                                $
          POSITION=plot_grid.position, COLOR=TrueColor(colors[i-1])
  ENDFOR
;
; Setup plot:
; ----------------------------------------------------------------------
;
  IF (focus NE 0) THEN BEGIN
    plot_xi = [1,1,1,1,1,1]
    plot_yi = [1,1,1,1,1,1]
    plot_xn = 1
    plot_yn = 1
    plot_xboarder = 0.0
    plot_yboarder = 0.0
    plot_xspacing = 0.0
    plot_yspacing = 0.0
    xspan  = 0.9
    yspan  = 0.8
    xshift = 0.1
    yshift = -0.1
  ENDIF ELSE BEGIN
    plot_xi = [1,1,1,2,2,2]
    plot_yi = [1,2,3,1,2,3]
    plot_xn = 2
    plot_yn = 3
    plot_xboarder = 0.05
    plot_yboarder = 0.1
    plot_xspacing = 0.125
    plot_yspacing = 0.025
    xspan  = 0.8
    yspan  = 1.0
    xshift = 0.25
    yshift = 0.0
  ENDELSE

  xsize = (xspan - 2.0 * plot_xboarder - 1.0 * plot_xspacing) / FLOAT(plot_xn) 
  ysize = (yspan - 2.0 * plot_yboarder - 2.0 * plot_yspacing) / FLOAT(plot_yn)

  FOR iplot = 1, 6 DO BEGIN
    IF (focus NE 0 AND focus NE iplot) THEN CONTINUE

    xi = plot_xi[iplot-1]
    yi = plot_yi[iplot-1]

    xcen  =       xsize * (0.5 + FLOAT(xi-1)) + 1.5 * plot_xboarder + plot_xspacing * FLOAT(xi-1) + xshift
    ycen  = 1.0 - ysize * (0.5 + FLOAT(yi-1)) - 1.0 * plot_yboarder - plot_yspacing * FLOAT(yi-1) + yshift
    xpos  = [xcen - 0.5 * xsize, xcen + 0.5 * xsize]
    ypos  = [ycen - 0.5 * ysize, ycen + 0.5 * ysize]

    xmin =  1.0E+35
    xmax = -1.0E+35
    FOR idata = 1, ndata DO BEGIN
      plot_data = cortex_ExtractStructure(plot_array,idata)
      plasma = plot_data.plasma 
      FOR iseg = 1, plot.line_seg_n DO BEGIN
        IF (plot_grid.line_seg_s[iseg] EQ -1) THEN CONTINUE  ; No intersection with the grid was found
        inter  = cortex_ExtractStructure(plot_data.inter,iseg)
        xdata = inter.psin ; inter.dist
        xmin = MIN([xmin,xdata])
        xmax = MAX([xmax,xdata])
      ENDFOR
    ENDFOR 
    IF (N_ELEMENTS(WHERE(plot.xrange EQ 0.0)) NE 2) THEN BEGIN
      xmin = MAX([plot.xrange[0],xmin])
      xmax = MIN([plot.xrange[1],xmax])
    ENDIF

    PRINT, plot.xrange
    PRINT, '1 XMIN,MAX=',iplot,xmin,xmax

    ymin =  1.0E+35
    ymax = -1.0E+35
    FOR idata = 1, ndata DO BEGIN
      plot_data = cortex_ExtractStructure(plot_array,idata)
      plasma = plot_data.plasma 
      FOR iseg = 1, plot.line_seg_n DO BEGIN
        IF (plot_grid.line_seg_s[iseg] EQ -1) THEN CONTINUE  ; No intersection with the grid was found
        inter  = cortex_ExtractStructure(plot_data.inter,iseg)
        xdata = inter.psin ; inter.dist
        k = WHERE(xdata GE xmin AND xdata LE xmax)
        i = inter.cell[k]
;       This is a check to make sure that the 1:1 assumption of OBJ to CELL/PLASMA in OSM
;       is valid for this particular data set.  The index information is passed along with
;       the simulation results, so if this assumption is violated there will have to be
;       some additional index searching/mapping here to get the correct data.
        FOR j = 0, N_ELEMENTS(i)-1 DO BEGIN
          IF (plasma.index[i[j]] NE i[j]+1) THEN BEGIN
            PRINT,'ERROR PlotRadialProfile: 1:1 mapping between the GRID and PLASMA arrays has failed'
            PRINT,'  I            = ',i[j]
            PRINT,'  PLASMA.INDEX = ',plasma.index[i[j]]
            RETURN, -1
          ENDIF
        ENDFOR
        IF (iplot EQ 1) THEN ydata = inter.l[k]
        IF (iplot EQ 2) THEN ydata = plasma.dens[i]
        IF (iplot EQ 3) THEN ydata = plasma.vi[i]
        IF (iplot EQ 4) THEN ydata = plasma.vi[i] / plasma.cs[i]
        IF (iplot EQ 5) THEN ydata = plasma.pe[i] + plasma.pi[i]
        IF (iplot EQ 6) THEN ydata = plasma.te[i]
;        IF (iplot EQ 6) THEN ydata = [plasma.te[i],plasma.ti[i]]
        ymin = MIN([ymin,ydata])
        ymax = MAX([ymax,ydata])
      ENDFOR ; iseg loop
    ENDFOR ; idata loop

    ydelta = ymax - ymin
    ymin = ymin - 0.05 * ydelta
    ymax = ymax + 0.05 * ydelta

    PRINT, '2 XMIN,MAX=',iplot,xmin,xmax
    PRINT, '2 YMIN,MAX=',iplot,ymin,ymax

;   Axes:
    xrange = [xmin,xmax]
    yrange = [ymin,ymax]
    position = [xpos[0],ypos[0],xpos[1],ypos[1]]

;    IF (iplot EQ 6) THEN yrange = [0.0,500.0]

    ytitle = STRTRIM(plot_ytitle[iplot-1],2)

    CASE (iplot) OF
      focus: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=ytitle, /NOERASE
      1    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=ytitle,XTICKFORMAT='(A1)',/NOERASE                                    
      2    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=ytitle,XTICKFORMAT='(A1)',/NOERASE
      3    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=ytitle,XTITLE=plot_xtitle,/NOERASE
      4    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=ytitle,XTICKFORMAT='(A1)',/NOERASE
      5    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=ytitle,XTICKFORMAT='(A1)',/NOERASE
      6    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
                   POSITION=position,YTITLE=ytitle,XTITLE=plot_xtitle,/NOERASE                                 
    ENDCASE

;    XYOUTS, 0.95*xmin+0.05*xmax, 0.1*ymin+0.9*ymax, plot_title
;   Data:
    FOR idata = 1, ndata DO BEGIN
      plot_data = cortex_ExtractStructure(plot_array,idata)
      plasma = plot_data.plasma 
      FOR iseg = 1, plot.line_seg_n DO BEGIN
        IF (plot_grid.line_seg_s[iseg] EQ -1) THEN CONTINUE 
        inter  = cortex_ExtractStructure(plot_data.inter,iseg)
        xdata = inter.psin ; inter.dist
        i = inter.cell
        IF (iplot EQ 1) THEN ydata = inter.l
        IF (iplot EQ 2) THEN ydata = plasma.dens[i]
        IF (iplot EQ 3) THEN ydata = plasma.vi[i]
        IF (iplot EQ 4) THEN ydata = plasma.vi[i] / plasma.cs[i]
        IF (iplot EQ 5) THEN ydata = plasma.pe[i] + plasma.pi[i]
        IF (iplot EQ 6) THEN ydata = plasma.te[i]

        OPLOT, [1.0,1.0], [ymin,ymax], LINESTYLE=3 

        OPLOT, xdata, ydata, COLOR=TrueColor(colors[idata+iseg-2])

        IF (iplot EQ 5) THEN BEGIN
          OPLOT, xdata, plasma.pe[i], COLOR=TrueColor(colors[idata+iseg-2]), LINESTYLE=2 
          OPLOT, xdata, plasma.pi[i], COLOR=TrueColor(colors[idata+iseg-2]), LINESTYLE=1
        ENDIF
;        IF (iplot EQ 6) THEN BEGIN
;          OPLOT, xdata, plasma.ti[i], COLOR=TrueColor(colors[idata+iseg-2]), LINESTYLE=2 
;        ENDIF

      ENDFOR  ; iseg loop
    ENDFOR  ; idata loop
  ENDFOR  ; iplot loop

  RETURN, 0

END
;
; ======================================================================
;



