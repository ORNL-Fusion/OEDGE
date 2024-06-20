FUNCTION cortex_SliceGrid, line, grid, status
;
; Slide the grid with the nasty line segment:
; ----------------------------------------------------------------------
;
  a1 = line[0:1]
  a2 = line[2:3]

  if ( ABS(a1[0]-a2[0]) LT 1.0E-7 ) THEN a1[0] = a1[0] + 0.000001

  length = SQRT( (a1[0] - a2[0])^2 + (a1[1] - a2[1])^2 )

 print, '  ------slicing-------'

; print,grid.isep
; print,grid.tube_psin
; print,a1
; print,a2,length
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

    IF ( ABS(b1[0]-b2[0]) LT 1.0E-7 ) THEN b1[0] = b1[0] + 0.000001D

    Lint, a1, a2, b1, b2, c, d, flag=flag

    IF (flag NE 1) THEN CONTINUE

    s = (c[0] - a1[0]) / (a2[0] - a1[0])    
    t = (c[0] - b1[0]) / (b2[0] - b1[0])    

;   IF (grid.obj_tube[iobj] EQ 1) THEN BEGIN
;     print, ' -->',b1,b2,iobj,s,t   
;   ENDIF

    IF ( s GE 0.0 AND s LE 1.0 AND  $
         t GE 0.0 AND t LT 1.0) THEN BEGIN

;      print,'carve!',s,t,iobj

      i = grid.obj_tube[iobj] - 1

      IF (NOT KEYWORD_SET(int_cell)) THEN BEGIN
        int_cell = grid.obj_cell[iobj]
        int_pos  = grid.obj_pos [iobj]
        int_tube = grid.obj_tube[iobj]
        int_frac = s
        int_dist = s * length
        int_psin = grid.tube_psin[i]
        int_rho  = grid.tube_rho [i]
        int_l    = grid.tube_l   [i]
        int_t    = t      
      ENDIF ELSE BEGIN
        int_cell = [int_cell,grid.obj_cell[iobj]]
        int_pos  = [int_pos ,grid.obj_pos [iobj]]
        int_tube = [int_tube,grid.obj_tube[iobj]]
        int_frac = [int_frac,s]
        int_dist = [int_dist,s * length]
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
    result = {           $
      cell : int_cell ,  $
      pos  : int_pos  ,  $
      tube : int_tube ,  $
      frac : int_frac ,  $
      dist : int_dist ,  $
      psin : int_psin ,  $
      rho  : int_rho  ,  $
      l    : int_l    ,  $
      t    : int_t    }          
  ENDELSE

  RETURN,result
END



;
; ======================================================================
;
FUNCTION cortex_PlotRadialProfile_NEW, plot, data_array, ps=ps

  PRINT
  PRINT,'----------------------- NEW PLOT -----------------------'
  PRINT,'PLOT OPTION=',plot.option

;help,data_array.inter,/struct

  MAXNYDATA = 7
  MAXNXDATA = 1000

  focus = plot.focus

  !P.BACKGROUND = TrueColor('White')

  dev_xsize = !D.X_SIZE
  dev_ysize = !D.Y_SIZE

  notes    = plot.notes
  charsize = plot.charsize
  charsize_labels = charsize

  !P.CHARSIZE = charsize
  !P.CHARTHICK = plot.thick
  !P.THICK    = plot.thick
  !X.THICK    = plot.thick
  !Y.THICK    = plot.thick
  !Z.THICK    = plot.thick

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
       PRINT, 'PlotRadialProfile: Invalid option, OPTION = 1'
       STOP
       END
    2: BEGIN
       PRINT, 'PlotRadialProfile: Invalid option, OPTION = 2'
       STOP
       END
    2: BEGIN
       PRINT, 'PlotRadialProfile: Invalid option, OPTION = 3'
       STOP
       END
    4: BEGIN
       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot = 4
       plot_xn = 2
       plot_yn = 2
       title = plot.title 
       subtitle = ['1 / unit','2 / unit',  $
                   '3 / unit','4 / unit']
       xtitle   = 'distance (m)'
       ytitle   = ['1 (unit)','2 (unit)','3 (unit)','4 (unit)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       END
    5: BEGIN
       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot = 4
       plot_xn = 2
       plot_yn = 2
       title = plot.title 
       subtitle = ['1 / unit','2 / unit',  $
                   '3 / unit','4 / unit']
       xtitle   = 'flux-tube No.'
       ytitle   = ['1 (unit)','2 (unit)','3 (unit)','4 (unit)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       END
;   --------------------------------------------------------------------
    ELSE: BEGIN  
      PRINT, 'ERROR cortex_PlotRadialProfile: Unrecognised plot option'
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
      PRINT, 'ERROR cortex_PlotRadialProfile: No data found'
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
          END
;       ----------------------------------------------------------------
        2: BEGIN
          END
;       ----------------------------------------------------------------
        3: BEGIN
          END
;       ----------------------------------------------------------------
        4: BEGIN  ; Plot for C-Mod inner wall midplane gas puff

;          help,val,/struct
          nseg = plot.line_seg_n
          iseg = nseg
          inter = cortex_ExtractStructure(val.inter,iseg)
;          help,inter,/struct
;          i = inter.cell
;          print,inter.psin
;          print,i
;          stop

          file = val.plasma.file
          integral = ' '
          str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
          str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
          labels[0] = labels[0] + STRING(idata-1) + '/' + str[0] + integral + ' :'

          ntrace = [1,2,1,1] ; Number of data lines on each plot

          xdata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      
          ydata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      

          i = inter.cell
          count_i = N_ELEMENTS(i)

          CASE iplot OF
            1: BEGIN
              xdata[0:count_i-1,0] = inter.dist
              ydata[0:count_i-1,0] = val.plasma.dens[i]
              END
            2: BEGIN
              xdata[0:count_i-1,0] = inter.dist
              ydata[0:count_i-1,0] = val.plasma.te[i]
              xdata[0:count_i-1,1] = xdata[0:count_i-1,0]
              ydata[0:count_i-1,1] = 2.0*val.plasma.te[i]
              END
            3: BEGIN
              xdata[0:count_i-1,0] = inter.dist
              ydata[0:count_i-1,0] = val.plasma.te[i]
              END
            4: BEGIN
              xdata[0:count_i-1,0] = inter.dist
              ydata[0:count_i-1,0] = val.plasma.dens[i]
              END
          ENDCASE
          END
;       ----------------------------------------------------------------
        5: BEGIN  ; Average quantities on a flux tube (debugging plot)

          file = val.plasma.file
          integral = ' '
          str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
          str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
          labels[0] = labels[0] + STRING(idata-1) + '/' + str[0] + integral + ' :'

          ntrace = [1,1,1,1] ; Number of data lines on each plot

          ntube = MAX(val.plasma.tube)
          xdata = FINDGEN(ntube) + 1.0
          ydata = xdata
          CASE iplot OF
            1: FOR i = 1, ntube DO ydata[i-1] = ALOG10(MEAN(val.plasma.dens[WHERE(val.plasma.tube EQ i)]))
            2: FOR i = 1, ntube DO ydata[i-1] =        MEAN(val.plasma.vi  [WHERE(val.plasma.tube EQ i)])
            3: FOR i = 1, ntube DO ydata[i-1] = ALOG10(MEAN(val.plasma.te  [WHERE(val.plasma.tube EQ i)]))
            4: FOR i = 1, ntube DO ydata[i-1] = ALOG10(MEAN(val.plasma.ti  [WHERE(val.plasma.tube EQ i)]))
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
      name = 'data' + STRING(idata,FORMAT='(I0)')
      data = { n : ndata, x : xdata, y : ydata, file : file }  ;  *** NEED A WAY TO ONLY STORE THE NON -999.0 DATA ***

      IF (idata EQ 1) THEN data_store = CREATE_STRUCT(           name,data) ELSE  $
                           data_store = CREATE_STRUCT(data_store,name,data)

      FOR itrace = 1, ntrace[iplot-1] DO BEGIN
        IF (N_ELEMENTS(WHERE(plot.xrange EQ 0.0)) NE 2) THEN BEGIN
          i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1],count)        
          IF (count EQ 0) THEN BEGIN
            PRINT,'ERROR cortex_PlotRadialProfile: No data within XRANGE'
            PRINT,'  IPLOT          = ',iplot
            PRINT,'  IDATA          = ',idata
            PRINT,'  XRANGE         = ',plot.xrange
            PRINT,'  XDATA MIN,MAX  = ',MIN(xdata),MAX(xdata)
            RETURN, -1
          ENDIF
        ENDIF ELSE i = WHERE(xdata[*,itrace-1] NE -999.0)
        xmin = MIN([xmin,xdata[i,itrace-1]])
        xmax = MAX([xmax,xdata[i,itrace-1]])
        ymin = MIN([ymin,ydata[i,itrace-1]])
        ymax = MAX([ymax,ydata[i,itrace-1]])
      ENDFOR
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
        4: BEGIN
          cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                          dev_xsize, dev_ysize, charsize_labels, colors

          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 

          i = WHERE(val.x[*,0] NE -999.0)
          OPLOT, val.x[i,0], val.y[i,0], COLOR=TrueColor(colors[idata-1])

          CASE iplot OF
            1: 
            2: BEGIN
              OPLOT, val.x[i,1], val.y[i,1], COLOR=TrueColor(colors[idata-1]), LINESTYLE=2
              END
            3: 
            ELSE:
          ENDCASE
          END
;       ----------------------------------------------------------------
        5: BEGIN
          cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                          dev_xsize, dev_ysize, charsize_labels, colors
          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
          OPLOT, val.x, val.y, COLOR=TrueColor(colors[idata-1])
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotRadialProfile: Unrecognised plot option'
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
FUNCTION cortex_PlotRadialProfile, plot, plot_array, ps=ps

  ndata = N_ELEMENTS(TAG_NAMES(plot_array))
  IF (ndata LE 0) THEN BEGIN
    PRINT, 'ERROR Cortex PlotRadialProfiles: No data found'
    RETURN, -1
  ENDIF

  IF (plot.option GE 4) THEN BEGIN
    PRINT, 'Redirecting to new radial plot code!'
    result = cortex_PlotRadialProfile_NEW(plot, plot_array, ps=ps)
    RETURN, result
  ENDIF

  focus = plot.focus

  notes    = plot.notes
  charsize = plot.charsize

  title  = plot.title

  !P.CHARSIZE = charsize
  !P.CHARTHICK = plot.thick
  !P.THICK    = plot.thick
  !X.THICK    = plot.thick
  !Y.THICK    = plot.thick
  !Z.THICK    = plot.thick

  !P.BACKGROUND = TrueColor('White')

  dev_xsize = !D.X_SIZE
  dev_ysize = !D.Y_SIZE

  CASE plot.xdata OF
    'psi_n'   : plot_xtitle = 'psi_n'
    'distance': plot_xtitle = 'distance along line (m)'
    'rho'     : plot_xtitle = 'rho'
  ENDCASE


  plot_subtitle = ['CONNECTION LENGTH'        ,'ELECTRON DENSITY'     ,'PARALLEL FLOW VELOCITY',  $
                   'PARALLEL FLOW MACH NUMBER','TOTAL PLASMA PRESSURE','PLASMA TEMPERATURES','(Te-solid, Ti-dashed)'] 
  plot_xtitle = 'rho (m)'
  plot_ytitle = ['L (m)','ne (m-3)','v|| (m s-1)','Mach no.','p (Pa)','Te,i (eV)']
  colors = ['Black','Red','Green','Blue','Darkseagreen', 'Hotpink', 'Orange', 'Silver']

;  plot_ytitle = ['1','n_e (m-3)','3','v_parallel / M','5','T_e (eV)']
;  colors = ['Black',   'Red','Green','Blue','Orange','Purple', 'Hotpink', 'Darkseagreen', 'Silver',  $
;            'Darkred', 'Greenyellow']
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
        CASE plot.xdata OF
          'psi_n'   : xdata = inter.psin 
          'distance': xdata = inter.dist 
          'rho'     : xdata = inter.rho
        ENDCASE
        xmin = MIN([xmin,xdata])
        xmax = MAX([xmax,xdata])
      ENDFOR
    ENDFOR 
    IF (N_ELEMENTS(WHERE(plot.xrange EQ 0.0)) NE 2) THEN BEGIN
      xmin = MAX([plot.xrange[0],xmin])
      xmax = MIN([plot.xrange[1],xmax])
    ENDIF

;    PRINT, plot.xrange
;    PRINT, '1 XMIN,MAX=',iplot,xmin,xmax

    ymin =  1.0E+35
    ymax = -1.0E+35
    FOR idata = 1, ndata DO BEGIN
      plot_data = cortex_ExtractStructure(plot_array,idata)
      plasma = plot_data.plasma 
      FOR iseg = 1, plot.line_seg_n DO BEGIN
        IF (plot_grid.line_seg_s[iseg] EQ -1) THEN CONTINUE  ; No intersection with the grid was found
        inter  = cortex_ExtractStructure(plot_data.inter,iseg)
        CASE plot.xdata OF
          'psi_n'   : xdata = inter.psin 
          'distance': xdata = inter.dist 
          'rho'     : xdata = inter.rho
        ENDCASE
        k = WHERE(xdata GE xmin AND xdata LE xmax, count)
        IF (count EQ 0) THEN BEGIN
          PRINT,'ERROR PlotRadialProfile: No data within x-range'
          STOP
        ENDIF
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
;        IF (iplot EQ 6) THEN ydata = plasma.te[i]
        IF (iplot EQ 6) THEN ydata = [plasma.te[i],plasma.ti[i]]
        ymin = MIN([ymin,ydata])
        ymax = MAX([ymax,ydata])
      ENDFOR ; iseg loop
    ENDFOR ; idata loop

    ydelta = ymax - ymin
    ymin = ymin - 0.05 * ydelta
    ymax = ymax + 0.05 * ydelta

;    PRINT, '2 XMIN,MAX=',iplot,xmin,xmax
;    PRINT, '2 YMIN,MAX=',iplot,ymin,ymax

;   Axes:
;    xrange = [xmin,xmax]
;    yrange = [ymin,ymax]
    xrange = [xmin,xmax]
    IF (plot.yrange[0] NE 0.0 OR plot.yrange[1] NE 0.0) THEN  $
      yrange = plot.yrange ELSE yrange = [ymin,ymax]

    position = [xpos[0],ypos[0],xpos[1],ypos[1]]

;    IF (iplot EQ 6) THEN yrange = [5.0,3000.0]
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
      6    : PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $ ; /YLOG, $
                   POSITION=position,YTITLE=ytitle,XTITLE=plot_xtitle,/NOERASE                                 
    ENDCASE

    XYOUTS, 0.95*xmin+0.05*xmax, 0.1*ymin+0.9*ymax, plot_subtitle[iplot-1]
    IF (iplot EQ 6) THEN XYOUTS, 0.95*xmin+0.05*xmax, 0.2*ymin+0.8*ymax, plot_subtitle[iplot]
;   Data:
    FOR idata = 1, ndata DO BEGIN
      plot_data = cortex_ExtractStructure(plot_array,idata)
      plasma = plot_data.plasma 
      FOR iseg = 1, plot.line_seg_n DO BEGIN
        IF (plot_grid.line_seg_s[iseg] EQ -1) THEN CONTINUE 
        inter  = cortex_ExtractStructure(plot_data.inter,iseg)
        CASE plot.xdata OF
          'psi_n'   : xdata = inter.psin 
          'distance': xdata = inter.dist 
          'rho'     : xdata = inter.rho
        ENDCASE
        i = inter.cell
        IF (iplot EQ 1) THEN ydata = inter.l
        IF (iplot EQ 2) THEN ydata = plasma.dens[i]
        IF (iplot EQ 3) THEN ydata = plasma.vi[i]
        IF (iplot EQ 4) THEN ydata = plasma.vi[i] / plasma.cs[i]
        IF (iplot EQ 5) THEN ydata = plasma.pe[i] + plasma.pi[i]
        IF (iplot EQ 6) THEN ydata = plasma.te[i]

;       Print labels first so they are below the lines:
        IF (iplot EQ 1 or iplot EQ focus) THEN BEGIN 
          IF (focus NE 0) THEN BEGIN
            step = 0.03
          ENDIF ELSE BEGIN
            step = 0.10 ; 0.07
          ENDELSE
          frac = step * FLOAT(idata) + 0.1
          val = plasma
          str = val.file                   ; extract the case name from the data file
          str = STRSPLIT(str,'/',/EXTRACT)
          str = str[N_ELEMENTS(str)-1]     ; take the last sub-string
          str = STRSPLIT(str,'.',/EXTRACT)
          str = str[0]                     ; take the first one
          XYOUTS, 0.95*xmin+0.05*xmax,frac*ymin+(1.0-frac)*ymax,  $
                  str, COLOR=TrueColor(colors[idata-1])
        ENDIF

        OPLOT, [1.0,1.0], [ymin,ymax], LINESTYLE=3 

        OPLOT, xdata, ydata, COLOR=TrueColor(colors[idata+iseg-2])

        IF (iplot EQ 5) THEN BEGIN
          OPLOT, xdata, plasma.pe[i], COLOR=TrueColor(colors[idata+iseg-2]), LINESTYLE=2 
          OPLOT, xdata, plasma.pi[i], COLOR=TrueColor(colors[idata+iseg-2]), LINESTYLE=1
        ENDIF
        IF (iplot EQ 6) THEN BEGIN
          OPLOT, xdata, plasma.ti[i], COLOR=TrueColor(colors[idata+iseg-2]), LINESTYLE=2 
        ENDIF

      ENDFOR  ; iseg loop
    ENDFOR  ; idata loop


    ; Discharge data:
    IF (iplot EQ 6 AND plot.id NE 'unknown') THEN BEGIN

      fp1 = 3
      FREE_LUN,fp1        
      file_name = 'cortex_data/' + plot.id + '.dat'
      OPENW,fp1,file_name, ERROR=err
      IF (err NE 0) THEN BEGIN
        PRINT,'ERROR cortex_PlotRadialProfile: Problem opening data stream (to file)'
        PRINT,'  FILE_NAME= ',file_name
        PRINTF,-2,!err.msg
        FREE_LUN,fp1
        STOP
      ENDIF
      
;      PRINTF,fp1,'{NUMBER OF COLUMNS}  8',FORMAT='(A)'
      PRINTF,fp1,'*',FORMAT='(A)'
      PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
      PRINTF,fp1,'* RADIAL PLASMA PROFILE FROM CORTEX - UPSTREAM'
      PRINTF,fp1,'*',FORMAT='(A)'
      PRINTF,fp1,'* CASE ',plot.case_name[0],FORMAT='(2A)'
      PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
      PRINTF,fp1,'*',FORMAT='(A)'
      PRINTF,fp1,'*      psi','rho','ne'   ,'vi'     ,'M','pe'  ,'Te'  ,'Ti'  ,FORMAT='(8A10)'
      PRINTF,fp1,'*         ',' '  ,'(m-3)','(m s-1)',' ','(Pa)','(eV)','(eV)',FORMAT='(8A10)'

      FOR j = 0, N_ELEMENTS(i)-1 DO BEGIN
        PRINTF,fp1,  $
          inter.psin   [j] ,                    $
          inter.rho    [j] ,                    $
          plasma.dens[i[j]],                    $
          plasma.vi  [i[j]],                    $
          plasma.vi  [i[j]] / plasma.cs[i[j]],  $
          plasma.pe  [i[j]],                    $
          plasma.te  [i[j]],                    $
          plasma.ti  [i[j]],                    $
          inter.pos    [j] ,                    $
          inter.tube   [j] ,                    $
          FORMAT='(2F10.5,2E10.2,F10.2,E10.2,2F10.2,I10,I6)'        
      ENDFOR

      CLOSE,fp1
      FREE_LUN,fp1

    ENDIF

  ENDFOR  ; iplot loop


  print,'title ',title

;
; Put main title on the page
; ----------------------------------------------------------------------
  cortex_PageTitle, plot, ps, 'none', dev_xsize, dev_ysize, title, notes, charsize

  RETURN, 0

END
;
; ======================================================================
;



