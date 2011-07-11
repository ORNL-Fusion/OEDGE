;
; ======================================================================
;
FUNCTION cortex_PlotWallProfiles, plot, data_array, grid=grid, wall=wall, annotate=annotate, ps=ps

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
  plot_sum  = plot.sum

  max_value   = REPLICATE(-1.0E+10,10)
  max_erosion = max_value

  default_type = 1
  plot_xboarder = 0.05
  plot_yboarder = 0.1
  plot_xspacing = 0.100
  plot_yspacing = 0.025

  save_ypos = [-999.0,-999.0]

  first_plot = 1

  CASE option OF
;   --------------------------------------------------------------------
    1: BEGIN
       nplot = 3
       plot_xn = 1
       plot_yn = 3
       title = plot.title 
       subtitle = ['ATOMIC PARTICLE FLUX DENSITY / D m-2 s-1'       ,  $
                   'AVERAGE ATOM ENERGY / eV'                       ,  $
                   'ATOMIC ENERGY FLUX DENSITY ON THE WALL / MW m-2']
       xtitle   = 'WALL SEGMENT INDEX'
       ytitle   = ['flux_p (D m-2 s-1)','E_avg (eV)','flux_e (MW m-2)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace    = [1,1,1] 
       plot_type = [1,1,1]  
       END
;   --------------------------------------------------------------------
    2: BEGIN
       nplot = 3
       plot_xn = 1
       plot_yn = 3
       title = plot.title 
       subtitle = ['MOLECULAR PARTICLE FLUX DENSITY / D2 m-2 s-1'      ,  $
                   'AVERAGE MOLECULE ENERGY / eV (per D2)'             ,  $
                   'MOLECULAR ENERGY FLUX DENSITY ON THE WALL / MW m-2']
       xtitle   = 'WALL SEGMENT INDEX'
       ytitle   = ['flux_p (D2 m-2 s-1)','E_avg (eV)','flux_e (MW m-2)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace = [1, 1, 1] 
       plot_type = [1,1,1]  
       END
;   --------------------------------------------------------------------
    3: BEGIN
       nplot = 6

       IF (plot.title EQ 'unknown') THEN title = 'WALL EROSION RATE SUMMARY' ELSE  $
                                         title = plot.title 
       subtitle = ['none','none',  $
                   'ATOM FLUX / m-2 s-1',        $
                   'AVERAGE ATOM ENERGY / eV',   $
                   'IMPURITY INFLUX / m-2 s-1',  $
;                   'EROSION / mm s-1']
                   'EROSION / mm per 1000 shots']   ; For the 1000 shot plot
       xtitle = 'WALL INDEX'
       ytitle = ['none','none',  $
                 'atom flux density (D m-2 s-1)',        $
                 'average indicent atom energy (eV)',    $
                 'impurity influx (particles m-2 s-1)',  $
;                 'erosion rate (mm s-1)']
                 'erosion rate (mm / 1000 shots)']  ; For the 1000 shot plot
       labels = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace = [1,1,1,1,1,1]  ; Number of data lines on each plot
       IF (plot.show_grid) THEN BEGIN
         plot_xn = 3
         plot_yn = 2
         plot_type = [-1,3,1,1,1,1]  ; Type of plot
       ENDIF ELSE BEGIN
         plot_xn = 2
         plot_yn = 2
         plot_type = [ 0,0,1,1,1,1]  
       ENDELSE       
       END
;   --------------------------------------------------------------------
    4: BEGIN
       nplot = 3
       plot_xn = 1
       plot_yn = 3
       title = plot.title 
       subtitle = ['ATOMIC PARTICLE FLUX DENSITY / D m-2 s-1 (converted to D2 for Pa calculation)',  $
                   'MOLECULAR PARTICLE FLUX DENSITY / D2 m-2 s-1',  $
                   'PRESSURE / Pa'                               ]
       xtitle   = 'WALL SEGMENT INDEX'
       ytitle   = ['flux_atm (D m-2 s-1)','flux_mol (D2 m-2 s-1)','p_D2 (Pa)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace    = [1,1,1] 
       plot_type = [1,1,1]  
       END
;   --------------------------------------------------------------------
    5: BEGIN
       IF (plot.title EQ 'unknown') THEN  $
         title = 'DIVIMP WALL EROSION' ELSE  $
         title = plot.title 
       subtitle = ['none','none',  $
                   'PARTICLE FLUX TO SURFACE / m-2 s-1'  ,  $
                   'EROSION / s-1 m-2'                   ,  $
                   'DEPOSITION (as ions only) / s-1 m-2' ,  $
                   'NET EROSION / mm per 1000 pulses, 400 s each']
       xtitle = 'WALL INDEX'
       ytitle = ['none','none',  $
                 'particle flux (m-2 s-1)',  $
                 'erosion (m-2 s-1)'      ,  $
                 'deposition (m-2 s-1)'   ,  $
                 'net erosion (m-2 s-1)']
       labels = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace = [1,1,1,1,1,1]  ; Number of data lines on each plot
       IF (plot.show_grid) THEN BEGIN
         plot_xn = 3
         plot_yn = 2
         plot_type = [-1,3,1,1,1,1]  
       ENDIF ELSE BEGIN
         plot_xn = 2
         plot_yn = 2
         plot_type = [ 0,0,1,1,1,1]  
       ENDELSE       
       nplot = 6
       ndata = N_ELEMENTS(TAG_NAMES(data_array))
       IF (ndata EQ 1) THEN BEGIN
         labels[0] = 'ions:atoms'
         labels[1] = 'total:atoms'
         ntrace[2] = 2
         ntrace[3] = 2
       ENDIF
       END
;   --------------------------------------------------------------------
    ELSE: BEGIN  
      PRINT, 'ERROR cortex_PlotWallProfile: Unrecognised plot option'
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

  title_mod = 1

  FOR iplot = 1, nplot DO BEGIN
    IF ((focus NE 0 AND focus NE iplot) OR (plot_type[iplot-1] EQ 0)) THEN CONTINUE

    thick = !P.THICK

    ndata = N_ELEMENTS(TAG_NAMES(data_array))
    IF (ndata LE 0) THEN BEGIN
      PRINT, 'ERROR cortex_PlotWallProfile: No data found'
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

      file = val.wall.file
      integral = ' '
      str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
      str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
      IF (ndata GT 1 OR option NE 5) THEN BEGIN
        IF (first_plot) THEN  $
          labels[0] = labels[0] + STRING(idata-1) + '\' + str[0] + integral + ' :'
      ENDIF ELSE BEGIN
        IF (title_mod) THEN title = title + ', CASE= ' + str[0]
        title_mod = 0
      ENDELSE

      CASE option OF
;       ----------------------------------------------------------------
        1: BEGIN
          xdata = FIX(val.wall.index)
          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
          CASE iplot OF
            1: BEGIN
              ydata[*,0] = val.wall.atom_par_flux

              j = [107, 116, 167, 290]
              FOR i = 0, 3 DO BEGIN
                k = WHERE(val.wall.index EQ j[i])
                flux = ( val.wall.atom_par_flux[k] + val.wall.mol_par_flux[k]) / 2
                
                temp = 300 
 
                fact = 2.0 * SQRT(!PI / 2.0)

                mass = 2.0 * 1.67E-27

                pressure = fact * SQRT( mass ) * flux * SQRT( 1.38E-23 * temp )

                print, j[i], val.wall.atom_par_flux[k],  $
                             val.wall.mol_par_flux [k],  $
                             flux, temp, pressure, pressure * 7.5,  $
                       FORMAT='(I6,2E10.2,2X,E10.2,2X,I6,2X,2E10.2)'                     
              ENDFOR

              END
            2: BEGIN
              ydata[*,0] = val.wall.atom_avg_energy
              ;mid_atom_flux = 0.5 * (MIN(val.wall.atom_par_flux) + MAX(val.wall.atom_par_flux))
              ;i = WHERE(val.wall.atom_par_flux LT 0.01 * mid_atom_flux, count)
              ;IF (count GE 1) THEN ydata[i,0] = 0.0
              END
            3: ydata[*,0] = val.wall.atom_energy_flux  ; W m-2
          ENDCASE
          END
;       ----------------------------------------------------------------
        2: BEGIN
          xdata = FIX(val.wall.index)
          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
          CASE iplot OF
            1: ydata[*,0] = val.wall.mol_par_flux 
            2: ydata[*,0] = val.wall.mol_avg_energy * 2.0 
            3: ydata[*,0] = val.wall.mol_energy_flux  ; W m-2
          ENDCASE
          END
;       ----------------------------------------------------------------
        3: BEGIN
          xdata = INDGEN(N_ELEMENTS(val.wall.length)) + 1
          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
          CASE iplot OF
            1:
            2:
            3: ydata[*,0] = val.wall.in_par_atm_1
            4: ydata[*,0] = val.wall.in_ene_atm_1
            5: BEGIN
              ydata[*,0] = val.wall.em_par_atm_2_2
              ; Calculate total sputtered source:
              total_source = TOTAL(val.wall.em_par_atm_2_2 * val.wall.length)
              print,'total_source',total_source,val.wall.tot_em_par_atm_2_2[0],val.core.div_influx
              END
            6: BEGIN
              ; Calculate thickness eroded per second:
              CASE FIX(val.summary.ion_atomic_number) OF
                ;                    Density                          Mass per atom          Assume simple 
                ;          (in m)    (g/cm^3)   (convert to kg/m^3)   (kg / atom)            cubic structure  
                 4: atom_diameter = ( 1.85    * 1.0E+3              / (  9.012 * 1.67E-27) )^(-1.0/3.0)
                26: atom_diameter = ( 7.87    * 1.0E+3              / ( 55.845 * 1.67E-27) )^(-1.0/3.0)
                74: atom_diameter = (19.35    * 1.0E+3              / (183.840 * 1.67E-27) )^(-1.0/3.0)
                ELSE: BEGIN
                  PRINT,'ERROR cortex_PlotWallProfiles: Unrecognised element'
                  PRINT,'  A =',FIX(val.summary.ion_atomic_number)
                  END
              ENDCASE
              ; Atomic diameter from http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
              ;   Be=2.1E-10, Fe=2.5E-10, W=2.7E-10
              ; which is about the same as what I get from the above estimate.
              print,'atom_diameter',atom_diameter, FIX(val.summary.ion_atomic_number)

              particles_per_m2 = 1.0 / (atom_diameter^2)

              ;              impurity influx           surface atom density     layer thickness
              ;              (atoms / s / m^2)         (atoms / m^2 / layer)    (m / layer)       convert to mm  
              ydata[*,0] = ( val.wall.em_par_atm_2_2 / particles_per_m2     ) * atom_diameter   * 1.0E+3         

              ;                          pulse time   number of pulses 
              ydata[*,0] =  ydata[*,0] * 400.0      * 1000.0             ; For the 1000 shot plot

              ; Integrate up to get the total erosion:
              CASE FIX(val.summary.ion_atomic_number) OF
                 4: mass = 9.0122
                26: mass = 55.845
                74: mass = 183.84
                ELSE: BEGIN
                  PRINT,'ERROR cortex_PlotWallProfiles: Unrecognised element'
                  PRINT,'  A =',FIX(val.summary.ion_atomic_number)
                  END
              ENDCASE
              i = WHERE(ydata[*,0] GT 0.00001, count) 
              tot_area    = TOTAL(val.wall.area[i]) 
              tot_erosion = TOTAL(val.wall.em_par_atm_2_2 * val.wall.area) *  $
                            mass * 1.67E-27 * 1000.0 *  $
                            14.0 * 3600.0
              print,'debug',count,237-5,tot_area,tot_erosion

              END
          ENDCASE
          END
;       ----------------------------------------------------------------
        4: BEGIN
          xdata = FIX(val.wall.index)
          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
          CASE iplot OF
            1: ydata[*,0] = val.wall.atom_par_flux
            2: ydata[*,0] = val.wall.mol_par_flux 
            3: BEGIN
              atm_flux = val.wall.atom_par_flux / 2.0 ; (assuming a chamber where all D is convereted to D2)
              mol_flux = val.wall.mol_par_flux
             
              flux = (atm_flux + mol_flux) 
              temp = 300.0 
              fact = 2.0 * SQRT(!PI / 2.0)
              mass = 2.0 * 1.67E-27
              pressure = fact * SQRT( mass ) * flux * SQRT( 1.38E-23 * temp )

              ;j = [107, 116, 167, 290]
              j = [30]
              PRINT,'Pressure calculation'
              FOR i = 0, N_ELEMENTS(j)-1 DO BEGIN
                k = WHERE(val.wall.index EQ j[i])
                print, j[i], atm_flux[k],  $
                             mol_flux[k],  $
                             flux[k], temp, pressure[k], pressure[k] * 7.5,  $
                       FORMAT='(I6,2E10.2,2X,E10.2,2X,I6,2X,2E10.2)'                     
              ENDFOR

              ydata[*,0] = pressure

              END
          ENDCASE
          END
;       ----------------------------------------------------------------
        5: BEGIN
          xdata = INDGEN(N_ELEMENTS(val.wall.index)) + 1
          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
          CASE iplot OF
            1:
            2:
            3: BEGIN
              IF (ndata EQ 1) THEN BEGIN
                ydata[*,0] = val.wall.ion_flux
                ydata[*,1] = val.wall.atm_flux
              ENDIF ELSE  $
                ydata[*,0] = val.wall.ion_flux + val.wall.atm_flux
              END
            4: BEGIN
              ydata[*,0] = val.wall.tot_ero
              ydata[*,1] = val.wall.atm_ero
              END
            5: ydata[*,0] = val.wall.tot_dep
            6: BEGIN
              ; Calculate thickness eroded per second:
              CASE FIX(val.wall.atomic_number) OF
                ;                    Density                          Mass per atom          Assume simple 
                ;          (in m)    (g/cm^3)   (convert to kg/m^3)   (kg / atom)            cubic structure  
                 4: atom_diameter = ( 1.85    * 1.0E+3              / (  9.012 * 1.67E-27) )^(-1.0/3.0)
                26: atom_diameter = ( 7.87    * 1.0E+3              / ( 55.845 * 1.67E-27) )^(-1.0/3.0)
                74: atom_diameter = (19.35    * 1.0E+3              / (183.840 * 1.67E-27) )^(-1.0/3.0)
                ELSE: BEGIN
                  PRINT,'ERROR cortex_PlotWallProfiles: Unrecognised element'
                  PRINT,'  A =',FIX(val.summary.ion_atomic_number)
                  END
              ENDCASE
              ; Atomic diameter from http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
              ;   Be=2.1E-10, Fe=2.5E-10, W=2.7E-10
              ; which is about the same as what I get from the above estimate.
              print,'atom_diameter',atom_diameter, FIX(val.summary.ion_atomic_number)

              particles_per_m2 = 1.0 / (atom_diameter^2)

              ;              impurity influx     surface atom density     layer thickness
              ;              (atoms / s / m^2)   (atoms / m^2 / layer)    (m / layer)       convert to mm  
              ydata[*,0] = ( val.wall.tot_net  / particles_per_m2     ) * atom_diameter   * 1.0E+3         

              ;                          pulse time   number of pulses 
              ydata[*,0] =  ydata[*,0] * 400.0      * 1000.0             ; For the 1000 shot plot

              ; Integrate up to get the total erosion:
              CASE FIX(val.summary.ion_atomic_number) OF
                 4: mass = 9.0122
                26: mass = 55.845
                74: mass = 183.84
                ELSE: BEGIN
                  PRINT,'ERROR cortex_PlotWallProfiles: Unrecognised element'
                  PRINT,'  A =',FIX(val.summary.ion_atomic_number)
                  END
              ENDCASE

              tot_area = TOTAL(val.wall.area) 

              length   = val.wall.dist_2 - val.wall.dist_1


              circ = 2.0 * !PI * val.wall.r0
              print,'r0',val.wall.r0

              i = WHERE(plot.xrange EQ 0.0, count)
              IF (count EQ 2) THEN i = INDGEN(N_ELEMENTS(val.wall.index)) ELSE  $
                                   i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1], count)        
              IF (count EQ 0) THEN BEGIN
                PRINT,'cortex_PlotWallProfiles','No data within XRANGE'
                STOP
              ENDIF

              tot_erosion = TOTAL(val.wall.tot_ero[i] * length[i] * circ) *  $
                            mass * 1.67E-27 * 1000.0 *  $
                            14.0 * 3600.0
              tot_deposit = TOTAL(val.wall.tot_dep[i] * length[i] * circ) *  $
                            mass * 1.67E-27 * 1000.0 *  $
                            14.0 * 3600.0
              net_erosion = TOTAL(val.wall.tot_net[i] * length[i] * circ) *  $
                            mass * 1.67E-27 * 1000.0 *  $
                            14.0 * 3600.0

              print,'debug',count,tot_area,tot_erosion,tot_deposit,net_erosion
              print,'debug',ndata

              length = val.wall.dist_2 - val.wall.dist_1
              absfac1 = TOTAL(val.wall.ion_ero * length)
              absfac2 = TOTAL(val.wall.atm_ero * length)
              print,'debug: absfact=',absfac1,absfac2,absfac1+absfac2,val.wall.absfac

              END
          ENDCASE
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotWallProfile: Unrecognised plot option'
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
          PRINT,'ERROR cortex_PlotWallProfile: No data within XRANGE'
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
    IF (ymin GT 0.0 AND ymax GT 0.0) THEN ymin = 0.0  ; Makes things a bit clearer on the plots I think...
    IF (ymin LT 0.0 AND ymax LT 0.0) THEN ymax = 0.0
    deltay = ymax - ymin
    ymin = ymin - 0.12 * deltay
    ymax = ymax + 0.05 * deltay

;   Axes:
    xrange = [xmin,xmax]
    yrange = [ymin,ymax]
    position = [xpos[0],ypos[0],xpos[1],ypos[1]]

    type = default_type                                           
    IF (focus NE 0 OR yi EQ plot_yn OR iplot EQ nplot) THEN type = 2   ; Show x-axis label

    IF (plot_type[iplot-1] EQ 3) THEN type = 3

    first_plot = 0

    CASE (type) OF
      1: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTICKFORMAT='(A1)',/NOERASE                                  
      2: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTITLE=xtitle,/NOERASE                                 
      3: BEGIN  
        plot_grid        = plot
        plot_grid.size   = 0.6
        plot_grid.center = [0.6*xpos[0]+0.4*xpos[1],0.5*(ypos[0]+ypos[1])]
        status = cortex_PlotFluidGrid(plot_grid, grid, wall, 0, annotate, 'subordinate', 'outline', ps='on')
        END
    ENDCASE

    IF (type EQ 3) THEN CONTINUE

;   Write sub-title for each plot:
    XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
            (xy_label[2] * ypos[0] + xy_label[3] * ypos[1]) * dev_ysize,  $
            subtitle[iplot-1], CHARSIZE=charsize_labels, /DEVICE

;   Add a trace that gives the peak value for each wall segment:
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

;   Add a trace that gives the sum overall traces:
    IF (plot_sum NE 0) THEN BEGIN
      FOR idata = 1, ndata DO BEGIN
        val = cortex_ExtractStructure(data_store,idata)
        IF (idata EQ 1) THEN BEGIN
          CASE option OF
;           --------------------------------------------------------
            1: BEGIN
;               CASE iplot OF              ; *** LEFT OFF *** comenting out for now so can @cortex_make
;                 1: BEGIN
;                    xdata = val.x
;                    ydata = val.y
;                    save_ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
;                    END
;                 2: xdata = val.x  &  ydata = vay.y * save_ydata[*,0] 
;                 3: xdata = val.x  &  ydata = val.y
;                ENDCASE
               END
;           --------------------------------------------------------
          ENDCASE
        ENDIF ELSE BEGIN
          FOR i = 0, N_ELEMENTS(xdata)-1 DO BEGIN
            ydata[i,*] = 0.0
            FOR j = 0, N_ELEMENTS(ydata[0,*])-1 DO BEGIN
              IF (xdata[i] NE val.x[i]) THEN BEGIN
                print, 'trouble'
                STOP
              ENDIF
              CASE option OF
;               --------------------------------------------------------
                1: BEGIN
                   CASE iplot OF              
                     1: ydata[i,j] = ydata[i,j] + val.y[i,j]
                     2: ydata[i,j] = ydata[i,j] + val.y[i,j] * save_ydata[j,idata-1]
                     3: ydata[i,j] = ydata[i,j] + val.y[i,j]
                    ENDCASE
                   END
;               --------------------------------------------------------
              ENDCASE
            ENDFOR
          ENDFOR
        ENDELSE
        IF (iplot EQ 1) THEN save_ydata[*,idata-1] = val.y[*,0]
      ENDFOR
      IF (iplot EQ 1) THEN save_ydata[*,ndata] = ydata[*,0]
      IF (iplot EQ 2) THEN ydata[*,0] = ydata[*,0] / save_ydata[*,ndata]
      IF (plot_sum EQ 2) THEN BEGIN
        labels[0] = 'sum values:'
        file  = 'SUM'
        ndata = 1  
        name = 'data' + STRING(ndata,FORMAT='(I0)')
        data = { x : xdata, y : ydata, file : file } 
        data_store = CREATE_STRUCT(name,data)
      ENDIF ELSE BEGIN
        labels[0] = labels[0] + 'sum values :'
        file  = 'sum values'
        ndata = ndata + 1  
        name = 'data' + STRING(ndata,FORMAT='(I0)')
        data = { x : xdata, y : ydata, file : file } 
        data_store = CREATE_STRUCT(data_store,name,data)
      ENDELSE
;     Need to store the particle flux data so that it can be used
;     to average the average energy flux on the next pass of IPLOT:
    ENDIF

    ; Mark vertical lines on plot:
    IF (cortex_GetValues(plot.xmark,values)) THEN BEGIN
      FOR i = 0, N_ELEMENTS(values)-1 DO  $
         OPLOT, [values[i],values[i]], [ymin,ymax], LINESTYLE=1, COLOR=Truecolor('Grey')
    ENDIF

;   Data:
    FOR idata = 1, ndata DO BEGIN
      val = cortex_ExtractStructure(data_store,idata)
      SWITCH option OF
;       ----------------------------------------------------------------
        1: BEGIN
          IF (idata EQ 1) THEN  $
            cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                            dev_xsize, dev_ysize, charsize_labels, colors

          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
 
          IF (plot_peak EQ 1 AND idata EQ ndata) THEN thick = 2.0
          IF (plot_sum  EQ 1 AND idata EQ ndata) THEN thick = 2.0

          OPLOT, val.x, val.y[*,0], COLOR=TrueColor(colors[idata-1]), THICK=thick
          CASE iplot OF
            1: 
            2: 
            3: 
            ELSE:
          ENDCASE
          BREAK
          END
;       ----------------------------------------------------------------
        2:
        3: 
        5: BEGIN
          IF (idata EQ 1) THEN  $
            cortex_DrawKey, iplot-2, focus, labels, xy_label, xpos, ypos,  $
                            dev_xsize, dev_ysize, charsize_labels, colors
          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
          IF (plot_peak EQ 1 AND idata EQ ndata) THEN thick = 2.0
          IF (plot_sum  EQ 1 AND idata EQ ndata) THEN thick = 2.0
          OPLOT, val.x, val.y[*,0], COLOR=TrueColor(colors[idata-1]), THICK=thick

          i = WHERE(val.x GE xmin AND val.x LE xmax)
          max_value[iplot] = MAX([max_value[iplot],val.y[i,0]])
          IF (idata EQ ndata) THEN BEGIN
            ;XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
            ;        (xy_label[3] * ypos[0] + xy_label[2] * ypos[1]) * dev_ysize,  $
            ;XYOUTS, (0.38        * xpos[0] + 0.62        * xpos[1]) * dev_xsize,  $
            ;        (xy_label[2] * ypos[0] + xy_label[3] * ypos[1]) * dev_ysize,  $
            XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
                    (0.97        * ypos[0] + 0.03        * ypos[1]) * dev_ysize,  $
                    'MAX = '+STRTRIM(STRING(max_value[iplot]),2), CHARSIZE=charsize_labels, /DEVICE
            IF (option EQ 5 AND iplot EQ 6) THEN  $
              XYOUTS, (0.50 * xpos[0] + 0.50 * xpos[1]) * dev_xsize,  $
                      (0.97 * ypos[0] + 0.03 * ypos[1]) * dev_ysize,  $
                      'NET = '+STRTRIM(STRING(net_erosion),2) + ' g in 14 h',  $
                      CHARSIZE=charsize_labels, /DEVICE
          ENDIF

          ; Integrate up to get the total erosion:
          IF (iplot EQ 4) THEN BEGIN
;            tot_erosion = 0.0
          ENDIF        

          IF (option EQ 5 AND ndata EQ 1) THEN BEGIN
             print,'trying..',iplot
            SWITCH iplot OF
              3: 
              4: BEGIN
                print,'herererere!',max(val.y[*,1])
                OPLOT, val.x, val.y[*,1], COLOR=TrueColor(colors[1]), THICK=thick
                max_value[iplot] = MAX([max_value[iplot],val.y[i,1]])
                BREAK
                END
              ELSE:
            ENDSWITCH
          ENDIF
          BREAK
          END
;       ----------------------------------------------------------------
        4: BEGIN
          IF (idata EQ 1) THEN  $
            cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                            dev_xsize, dev_ysize, charsize_labels, colors

          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
          IF (plot_peak EQ 1 AND idata EQ ndata) THEN thick = 2.0
          IF (plot_sum  EQ 1 AND idata EQ ndata) THEN thick = 2.0
          OPLOT, val.x, val.y[*,0], COLOR=TrueColor(colors[idata-1]), THICK=thick
          BREAK
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotWallProfile: Unrecognised plot option'
          PRINT, '  OPTION = ',plot.option,' (',option,')'
          RETURN, -1
          END
      ENDSWITCH

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


