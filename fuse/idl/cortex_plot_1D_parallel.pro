;
; *** ADD {SHOW GRID} OPTION ***
; 

;
; ======================================================================
;
PRO cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos, dev_xsize, dev_ysize,  $
                    charsize_labels, colors

  step = 0.1
  IF (focus) THEN step = 0.05
  str = STRSPLIT(labels[iplot-1],':',/EXTRACT)
  FOR i = 0, N_ELEMENTS(str)-1 DO BEGIN       
    icolor = i
    str_label = str[i]
    str_color = STRSPLIT(str_label,'/',/EXTRACT)
    IF (N_ELEMENTS(str_color) GT 1) THEN BEGIN
      icolor = LONG(str_color[0])
      str_label = str_color[1]
    ENDIF
    frac = step * FLOAT(i) + 0.1
    XYOUTS, ( xy_label[0]         * xpos[0] +  xy_label[1]         * xpos[1]) * dev_xsize,  $
            ((xy_label[2] + frac) * ypos[0] + (xy_label[3] - frac) * ypos[1]) * dev_ysize,  $
           str_label, CHARSIZE=charsize_labels,COLOR=TrueColor(colors[icolor]), /DEVICE
  ENDFOR

END
;
; ======================================================================
;
PRO cortex_PlotNodes, plot_array, idata, tag, tube, color

  struct = cortex_ExtractStructure(plot_array,idata)                 
  struct = struct.node

  index = WHERE(TAG_NAMES(struct) EQ STRUPCASE(tag), count)

  IF (count NE 0) THEN BEGIN
    val = struct.(index[0]) 
  ENDIF ELSE BEGIN
    PRINT, 'ERROR cortex_PlotNodes: TAG not found'
    PRINT, '  TAG = ',tag
    RETURN
  ENDELSE

  i = WHERE(struct.tube EQ tube)
  s   = struct.s [*,i]
  val = val      [*,i]

; Some temporary tweaks...
  IF (tag EQ 'jsat') THEN BEGIN 
    val = val / 1.602E-19
    val[0] = -val[0]
  ENDIF
  IF (tag EQ 'pe') THEN val = val * 1.602E-19

  FOR j = 1, (struct.i_end[i])[0] DO BEGIN
    IF (val[j-1] EQ 0.0) THEN CONTINUE
    psym = 4
    IF (j EQ struct.i_sym[i]) THEN psym = 6
    OPLOT, [s[j-1]], [val[j-1]], PSYM=psym, COLOR=TrueColor(color)
  ENDFOR
END
;
; ======================================================================
;
FUNCTION cortex_PlotParallelProfiles, plot, tube, plot_array, ps=ps

  PRINT
  PRINT,'----------------------- NEW PLOT -----------------------'
  PRINT

  ndata = N_ELEMENTS(TAG_NAMES(plot_array))
  IF (ndata LE 0) THEN BEGIN
    PRINT, 'ERROR cortex_PlotParallelProfiles: No data found'
    RETURN, -1
  ENDIF

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

  CASE option OF
;   --------------------------------------------------------------------
    1: BEGIN
       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot = 11
       plot_xn = 3
       plot_yn = 4
       title = plot.title + ': TUBE = ' + STRTRIM(STRING(tube),2)
       subtitle = ['PLASMA DENSITY','PLASMA FLOW','PLASMA PRESSURE','PLASMA TEMPERATURE',  $
                   'PARTICLE SOURCES','PARTICLE FLUX','MOMENTUM SOURCES','PRESSURE','9','10','11','12']
       xtitle   = 's (m)'
       ytitle   = ['n (m-3)','Mach no.','p (?)','Te,i','particles (m-3 s-1)',  $
                  'parallel flux','momentum (?)','pressure (?)',  $
                  'D_a (ph m-3 s-1)','D_g (ph m-3 s-1)','D_a / D_g','12']
       labels   = ['n_e:n_D:n_D2',  $
                   'M',             $
                   'p:p_e:p_i',     $ 
                   'T_e:T_i',       $
                  'solver_net:solver_ion:solver_rec:solver_usr:solver_fit:1/eirene_ion (dashed):2/eirene_rec (dashed)',  $
                  'solver',         $
                  'solver_net:solver_vol:solver_usr:solver_fit',    $
                  'solver',         $
                  'eirene',         $
                  'eirene',         $
                  'eirene']
       END
;   --------------------------------------------------------------------
    999: BEGIN
       default_plot_type = 2
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.075
       nplot = N_ELEMENTS(WHERE(tube NE 0))
       IF (nplot GT 32) THEN BEGIN
         PRINT,'ERROR cortex_PlotParallelProfiles: Sorry, too many TUBES specified (MAX=32)'
         RETURN, -1
       ENDIF
       IF (nplot LE 32) THEN set_nx_ny = [4,8]  ; Move to a function, or resolve using an expression here...
       IF (nplot LE 16) THEN set_nx_ny = [4,4]
       IF (nplot LE 12) THEN set_nx_ny = [3,4]
       IF (nplot LE 9 ) THEN set_nx_ny = [3,3]
       IF (nplot LE 4 ) THEN set_nx_ny = [2,2]
       IF (nplot EQ 2 ) THEN set_nx_ny = [2,1]
       IF (nplot EQ 1 ) THEN set_nx_ny = [1,1]
       plot_xn = set_nx_ny[0]
       plot_yn = set_nx_ny[1]
       IF (plot.state EQ -1) THEN BEGIN
         PRINT,'ERROR cortex_PlotParallelProfiles: Charge STATE not specified'
         RETURN, -1
       ENDIF
       title = plot.title + ': COMPARING EIRENE AND DIVIMP IMPURITY PROFILES'
       state = plot.state
       IF (plot.option EQ 100) THEN subtitle = REPLICATE(['DENSITY +'+STRTRIM(STRING(state),2)],nplot)
       IF (plot.option EQ 101) THEN subtitle = REPLICATE(['IONISATION +'+STRTRIM(STRING(state),2)],nplot)
       xtitle   = 's (m)'
       ytitle   = REPLICATE(['n_I +'+STRTRIM(STRING(state),2)+' (m-3)'],nplot)
       labels   = REPLICATE(['eirene:divimp'],nplot)
       END
;   --------------------------------------------------------------------
    ELSE: BEGIN  
      PRINT, 'ERROR cortex_PlotParallelProfile: Unrecognised plot option'
      PRINT, '  OPTION = ',option,' (',plot.option,')'
      RETURN, -1
      END
  ENDCASE

  colors   = ['Black','Red','Blue','Orange','Purple', 'Hotpink', 'Green']

;
; Setup plot:
; ----------------------------------------------------------------------
;
  IF (focus) THEN BEGIN
    plot_xn = 1
    plot_yn = 1
  ENDIF

  xsize = (1.0 - 2.0 * plot_xboarder - FLOAT(plot_xn-1) * plot_xspacing) / FLOAT(plot_xn) 
  ysize = (1.0 - 2.0 * plot_yboarder - FLOAT(plot_yn-1) * plot_yspacing) / FLOAT(plot_yn)

  xi = 1
  yi = 0

  FOR iplot = 1, nplot DO BEGIN
    IF (focus NE 0 AND focus NE iplot) THEN CONTINUE

    yi = yi + 1
    IF (yi EQ plot_yn+1) THEN BEGIN
      xi = xi + 1
      yi = 1
    ENDIF
    
    xcen  =       xsize * (0.5 + FLOAT(xi-1)) + 1.5 * plot_xboarder + plot_xspacing * FLOAT(xi-1)
    ycen  = 1.0 - ysize * (0.5 + FLOAT(yi-1)) - 1.0 * plot_yboarder - plot_yspacing * FLOAT(yi-1)
    xpos  = [xcen - 0.5 * xsize, xcen + 0.5 * xsize]
    ypos  = [ycen - 0.5 * ysize, ycen + 0.5 * ysize]

    xmin =  1.0E+35
    xmax = -1.0E+35
    ymin =  1.0E+35
    ymax = -1.0E+35
    FOR idata = 1, ndata DO BEGIN
      val = cortex_ExtractStructure(plot_array,idata)

      CASE option OF
;       ----------------------------------------------------------------
        1: BEGIN
          file = val.target.file
          i = WHERE(val.target.tube EQ tube)
          j = WHERE(val.plasma.tube EQ tube)
          k = WHERE(val.source.tube EQ tube)
          l = WHERE(val.eirene.tube EQ tube)
          ntrace = [3, 1, 3, 2, 7, 1, 4, 1, 1, 1, 1]  ; Number of data lines on each plot
          CASE iplot OF
            1 : xdata = [val.target.s[i,0],val.plasma.s[j],val.target.s[i,1]]
            2 : xdata = [val.target.s[i,0],val.plasma.s[j],val.target.s[i,1]]
            3 : xdata = [val.target.s[i,0],val.plasma.s[j],val.target.s[i,1]]
            4 : xdata = [val.target.s[i,0],val.plasma.s[j],val.target.s[i,1]]
            5 : xdata = val.source.s[j]
            6 : xdata = val.source.s[j]
            7 : xdata = val.source.s[j]
            8 : xdata = val.source.s[j]
            9 : xdata = val.eirene.s[j]
            10: xdata = val.eirene.s[j]
            11: xdata = val.eirene.s[j]
          ENDCASE
          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
          CASE iplot OF
            1 : BEGIN
                ydata[*,0] = [val.target.dens[i,0],val.plasma.dens    [j],val.target.dens[i,1]]
                ydata[*,1] = [0.0                 ,val.eirene.atm_dens[l],0.0                 ]
                ydata[*,2] = [0.0                 ,val.eirene.mol_dens[l],0.0                 ]
                END
            2 : BEGIN
                ydata[*,0] = [val.target.M   [i,0],val.plasma.M   [j],val.target.M   [i,1]]
                END
            3 : BEGIN
                ydata[*,1] = [val.target.pe[i,0],val.plasma.pe[j],val.target.pe[i,1]]
                ydata[*,2] = [val.target.pi[i,0],val.plasma.pi[j],val.target.pi[i,1]]
                ydata[*,0] = ydata[*,1] + ydata[*,2]
                END
            4 : BEGIN
                ydata[*,0] = [val.target.te[i,0],val.plasma.te[j],val.target.te[i,1]]
                ydata[*,1] = [val.target.ti[i,0],val.plasma.ti[j],val.target.ti[i,1]]
                END
            5 : BEGIN
                ydata[*,0] = val.source.par_net[k]
                ydata[*,1] = val.source.par_ion[k]
                ydata[*,2] = val.source.par_rec[k]
                ydata[*,3] = val.source.par_usr[k]
                ydata[*,4] = val.source.par_ano[k]
                ydata[*,5] = val.eirene.ion_net[l]
                ydata[*,6] = val.eirene.rec_net[l]
                END
            6 : ydata[*,0] = val.plasma.dens[j] * val.plasma.vi[j]  
            7 : BEGIN
                ydata[*,0] = val.source.mom_net[k]
                ydata[*,1] = val.source.mom_vol[k]
                ydata[*,2] = val.source.mom_usr[k]
                ydata[*,3] = val.source.mom_ano[k]
                END
            8 : ydata[*,0] = val.plasma.pe  [j] + val.plasma.pi[j]
            9 : ydata[*,0] = val.eirene.balmer_alpha[k]
            10: ydata[*,0] = val.eirene.balmer_gamma[k]
            11: BEGIN
               dalpha = val.eirene.balmer_alpha[k]
               dgamma = val.eirene.balmer_gamma[k]
               m = WHERE(dgamma EQ 0.0,count)
               IF (count GT 0            ) THEN dgamma[m] = dalpha[m] ELSE  $
               IF (count EQ N_ELEMENTS(k)) THEN dgamma[m] = 1.0
               ydata[*,0] = dalpha / dgamma
               END
          ENDCASE
          END
;       ----------------------------------------------------------------
        999: BEGIN
          subtitle[iplot-1] = subtitle[iplot-1] + ' TUBE='+STRTRIM(STRING(tube[iplot-1]),2)
          file = val.eirene.file
          i = WHERE(val.eirene.tube EQ tube[iplot-1],count_i)
          j = WHERE(val.divimp.tube EQ tube[iplot-1],count_j)
          IF (count_i EQ 0 OR count_j EQ 0) THEN BEGIN
            PRINT, 'ERROR cortex_PlotParallelProfile: Data not found'
            PRINT, '  OPTION   =',plot.option
            PRINT, '  COUNT_I,J=',count_i,count_j
            PRINT, '  IPLOT    =',iplot
            PRINT, '  TUBE     =',tube[iplot-1]
            RETURN, -1
          ENDIF
          ntrace = MAKE_ARRAY(99,/LONG,VALUE=2)
          xdata = val.eirene.s[i]
          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
          CASE (plot.option) OF
            100: BEGIN
              ydata[*,0] = val.eirene.imp_dens[i,state]
              ydata[*,1] = val.divimp.imp_dens[j,state] * val.divimp.div_influx
              END
            101 : BEGIN
              ydata[*,0] = val.eirene.imp_ioniz[i,state]
              ydata[*,1] = val.divimp.imp_ioniz[j,state] * val.divimp.div_influx
              END
          ENDCASE
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotParallelProfile: Unrecognised plot option'
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
        xmax1 = MAX(xdata)
        i = WHERE(xdata/xmax1 GE plot.xrange[0] AND xdata/xmax1 LE plot.xrange[1])        
        IF (N_ELEMENTS(i) EQ 1) THEN BEGIN
          PRINT,'ERROR cortex_PlotCoreProfiles: No data within XRANGE'
          PRINT,'  IPLOT  = ',iplot
          PRINT,'  IDATA  = ',idata
          PRINT,'  XRANGE = ',plot.xrange
          PRINT,'  XDATA  = '
          PRINT,xdata
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
    ymin = ymin - 0.05 * deltay
    ymax = ymax + 0.05 * deltay

;   Axes:
    xrange = [xmin,xmax]
    yrange = [ymin,ymax]
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
        1: BEGIN
          cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                          dev_xsize, dev_ysize, charsize_labels, colors

          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 

          OPLOT, val.x, val.y[*,0], COLOR=TrueColor(colors[0]) 
          CASE iplot OF
            1: BEGIN
               OPLOT, val.x, val.y[*,1], COLOR=TrueColor(colors[1])
               OPLOT, val.x, val.y[*,2], COLOR=TrueColor(colors[2])
               IF (plot.nodes) THEN cortex_PlotNodes, plot_array, idata, 'dens', tube, colors[0]
               END
            3: BEGIN
               OPLOT, val.x, val.y[*,1], COLOR=TrueColor(colors[1])
               OPLOT, val.x, val.y[*,2], COLOR=TrueColor(colors[2])
               IF (plot.nodes) THEN cortex_PlotNodes, plot_array, idata, 'pe', tube, colors[1]
               END
            4: BEGIN
               OPLOT, val.x, val.y[*,1], COLOR=TrueColor(colors[1])
               IF (plot.nodes) THEN BEGIN
                 cortex_PlotNodes, plot_array, idata, 'te', tube, colors[0]
                 cortex_PlotNodes, plot_array, idata, 'ti', tube, colors[1]
               ENDIF
               END
            5: BEGIN
               OPLOT, val.x, val.y[*,1], COLOR=TrueColor(colors[1])
               OPLOT, val.x, val.y[*,2], COLOR=TrueColor(colors[2])
               OPLOT, val.x, val.y[*,3], COLOR=TrueColor(colors[3])
               OPLOT, val.x, val.y[*,4], COLOR=TrueColor(colors[4])
               OPLOT, val.x, val.y[*,5], COLOR=TrueColor(colors[1]), LINESTYLE=1
               OPLOT, val.x, val.y[*,6], COLOR=TrueColor(colors[2]), LINESTYLE=1
               END
            6: BEGIN
               IF (plot.nodes) THEN cortex_PlotNodes, plot_array, idata, 'jsat', tube, 'Black'
               END
            7: BEGIN
               OPLOT, val.x, val.y[*,1], COLOR=TrueColor(colors[1])
               OPLOT, val.x, val.y[*,2], COLOR=TrueColor(colors[2])
               OPLOT, val.x, val.y[*,3], COLOR=TrueColor(colors[3])
               END
            ELSE:
          ENDCASE
          END
;       ----------------------------------------------------------------
        999: BEGIN
          IF (focus OR iplot EQ 1) THEN  $
            cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                            dev_xsize, dev_ysize, charsize_labels, colors
          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
          OPLOT, val.x, val.y[*,0], COLOR=TrueColor(colors[0]) 
          OPLOT, val.x, val.y[*,1], COLOR=TrueColor(colors[1])
          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotParallelProfile: Unrecognised plot option'
          PRINT, '  OPTION = ',plot.option,' (',option,')'
          RETURN, -1
          END
      ENDCASE

    ENDFOR  ; idata loop

  ENDFOR  ; iplots loop

;
;
; ----------------------------------------------------------------------  (PUT IN A FUNCTION)
  ypos = 0.940
  IF (notes NE 'default') THEN ypos = 0.965

  str = val.file                   ; extract the case name from the data file
  str = STRSPLIT(str,'/',/EXTRACT)
  str = str[N_ELEMENTS(str)-1]     ; take the last sub-string
  str = STRSPLIT(str,'.',/EXTRACT)
  str = str[0]                     ; take the first one

  title = title + ': CASE ' + STRTRIM(str,2)

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
    ENDIF ELSE BEGIN
      XYOUTS, 0.02 * dev_xsize, (ypos - 0.025) * dev_ysize, CHARSIZE=1.0, notes, /DEVICE
    ENDELSE              
  ENDIF

  RETURN, 0

END
;
; ======================================================================
;


