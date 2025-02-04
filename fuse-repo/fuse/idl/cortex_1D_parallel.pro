;
; *** ADD {SHOW GRID} OPTION ***
; 

;
; ======================================================================
;
PRO cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos, dev_xsize, dev_ysize,  $
                    charsize_labels, colors, step=step, offset=offset, starting_icolour=starting_icolour

  IF (NOT KEYWORD_SET(step)) THEN BEGIN
    step = 0.105 ; 0.1
    IF (focus) THEN step = 0.05
  ENDIF
  IF (NOT KEYWORD_SET(offset)) THEN BEGIN
    offset = 0.11 ; 0.10
  ENDIF
  str = STRSPLIT(labels[iplot-1],':',/EXTRACT)
  FOR i = 0, N_ELEMENTS(str)-1 DO BEGIN       
    IF (ARG_PRESENT(starting_icolour)) THEN BEGIN
      IF (i EQ 0) THEN icolor = starting_icolour + 1 ELSE icolor++
;print,'icolor',icolor
    ENDIF ELSE BEGIN
      icolor = i
    ENDELSE
    str_label = str[i]
    str_color = STRSPLIT(str_label,'\',/EXTRACT)
    IF (N_ELEMENTS(str_color) GT 1) THEN BEGIN
      icolor = LONG(str_color[0])
      str_label = str_color[1]
    ENDIF
    frac = step * FLOAT(i) + offset
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
FUNCTION cortex_PlotParallelProfiles, plot, tube, plot_array, ps=ps, time_slice=time_slice

  PRINT
  PRINT,'----------------------- NEW PLOT -----------------------'

  ndata = N_ELEMENTS(TAG_NAMES(plot_array))
  IF (ndata LE 0) THEN BEGIN
    PRINT, 'ERROR cortex_PlotParallelProfiles: No data found'
    RETURN, -1
  ENDIF


;help,plot_array,/struct
;help,plot_array.data1.plasma,/struct
;psclose
;stop

  MAXNYDATA = 7

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

  option = plot.option
  IF (option EQ 100 OR option EQ 101 OR option EQ 102) THEN option = 999

  CASE option OF
;   --------------------------------------------------------------------
    1: BEGIN
       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot   = 8  
       plot_xn = 2
       plot_yn = 4
       title = plot.title + ': TUBE ' + STRTRIM(STRING(tube),2)
       ntrace = [3, 1, 3, 2, 7, 1, 4, 1, 1, 1, 1, 2]  ; Number of data lines on each plot
       subtitle = ['PLASMA DENSITY',          $
                   'PLASMA FLOW',             $
                   'PLASMA PRESSURE',         $
                   'PLASMA TEMPERATURE',      $
                   'PARTICLE SOURCES',        $
                   'PARALLEL PARTICLE FLUX',  $
                   'MOMENTUM SOURCES',        $
                   'PRESSURE']
       xtitle   = 's (m)'
       IF (plot.xdata EQ 'tar_dist') THEN xtitle = 'dist from target (m)'
       ytitle   = ['n (m-3)',              $ 
                   'Mach no.',             $ 
                   'p (Pa)',               $
                   'Te,i (eV)',            $
                   'particles (m-3 s-1)',  $
                   'par. flux (m-2 s-1)',  $
                   'momentum (Pa m-1)',    $
                   'pressure (Pa)']
       labels   = ['n_e:n_D:n_D2',  $
                   'M',             $
                   'p:p_e:p_i',     $ 
                   'T_e:T_i',       $
                   'solver_net:solver_ion:solver_rec:solver_usr:solver_fit:1\eirene_ion (dashed):2\eirene_rec (dashed)',  $
                   'solver',   $
                   'solver_net:solver_vol:solver_usr:solver_fit',    $
                   'solver']

       IF (0 EQ 1) THEN BEGIN
         nplot   = 12
         plot_xn =  3
         plot_yn =  4
         ntrace  = [ntrace, 1, 1, 1, 2]  ; Number of data lines on each plot

         subtitle = [subtitle,  $
                     'EIRENE',  $
                     'EIRENE',  $
                     'EIRENE',  $
                     'EIRENE']

         ytitle   = [ytitle        ,  $
                     '(ph m-3 s-1)',  $
                     '(ph m-3 s-1)',  $
                     ' ',             $
                     '(ph m-3 s-1)']

         labels   = [labels,     $
                     'D_alpha',  $
                     'D_gamma',  $
                     'D_alpha / D_gamma',  $
                     'D_alpha:D_gamma (scaled)']
       ENDIF

       IF (1 EQ 1) THEN BEGIN
         nplot   = 12
         plot_xn =  3
         plot_yn =  4
         ntrace  = [ntrace, 2, 2, 1, 2]  ; Number of data lines on each plot

         subtitle = [subtitle,  $
                     'ELECTRON ENERGY SOURCES',  $
                     'ION ENERGY SOURCES',       $
                     'EIRENE',  $
                     'EIRENE']

         ytitle = [ytitle ,  $
                   '(???)',  $
                   '(???)',  $
                   '(ph m-3 s-1)',  $
                   ' ']

         labels = [labels   ,  $
                   'eirene:solver_net',  $
                   'eirene:solver_net',  $
                   'D_alpha',  $
                   'D_alpha / D_gamma']

       ENDIF

       END
;   --------------------------------------------------------------------
    2: BEGIN
       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot   = 8
       plot_xn = 2
       plot_yn = 4
       title = 'PARALLEL PLASMA PROFILES' 
       subtitle = ['PLASMA DENSITY / m-3'    ,'PLASMA FLOW / Mach no.'    ,'PARALLEL FLOW / m s-1'    ,'ELECTRON STATIC PRESSURE / Pa'  , $
                   'ION STATIC PRESSURE / Pa','TOTAL STATIC PRESSURE / Pa','ELECTRON TEMPERATURE / eV','ION TEMPERATURE / eV']
       xtitle   = 's (m)'
       ytitle   = ['n (m-3)','Mach no.','v_par (m s-1)','p_e (Pa)','p_i (Pa)','p_tot (Pa)','T_e (eV)','T_i (eV)']
       labels   = ['','','','','','','','','','','','']
print,title
;       nplot   = 7
;      plot_xn = 3
;      plot_yn = 4
;      title = plot.title + ': TUBE = ' + STRTRIM(STRING(tube),2)
;      subtitle = ['PLASMA DENSITY','PLASMA FLOW','ELECTRON PRESSURE','ION PRESSURE', 'TOTAL PRESSURE', $
;                  'ELECTRON TEMPERATURE','ION TEMPERATURE','8','9','10','11','12']
;      xtitle   = 's (m)'
;      ytitle   = ['n (m-3)','Mach no.','p_e (?)','p_i (?)','p (?)','T_e (eV)','T_i (eV)','8','9','10','11','12']
;      labels   = ['','','','','','','','','','','','']
       END
;   --------------------------------------------------------------------
    3: BEGIN
       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot   = 1
       plot_xn = 1
       plot_yn = 1
       title = plot.title + ': TUBE = ' + STRTRIM(STRING(tube),2)
       subtitle = ['shit']
       xtitle   = 'p (m)'
       ytitle   = ['sit']
       labels   = ['none']
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
       CASE plot.option OF
         100 : title = plot.title + ': IMPURITY DENSITY'
         ELSE: title = plot.title + ': LINE EMISSION'
       ENDCASE
       
;       title = plot.title + ': COMPARING EIRENE AND DIVIMP IMPURITY PROFILES'
       state = plot.state
       IF (plot.option EQ 100) THEN subtitle = REPLICATE(['CHARGE STATE +'   +STRTRIM(STRING(state),2)],nplot)
       IF (plot.option EQ 101) THEN subtitle = REPLICATE(['IONISATION +'+STRTRIM(STRING(state),2)],nplot)
       IF (plot.option EQ 102) THEN BEGIN
         wlngth = plot_array.data1.emission.wlngth[state-1]
         CASE FIX(wlngth) OF
           514:  charge_state = 'CII'
           465:  charge_state = 'CIII'
           ELSE: charge_state = 'unknown'
         ENDCASE
         subtitle = REPLICATE(['LINE ' + STRING(STRTRIM(wlngth,2)) + ' nm ' + charge_state],nplot)  ; 'EMISSION +'  +STRTRIM(STRING(state),2)],nplot)
       ENDIF
;       IF (plot.option EQ 102) THEN xtitle   = 'p (m)' ELSE $
                                    xtitle   = 's (m)'
       
       CASE plot.option OF
         100 : ytitle = REPLICATE(['n_I +'+STRTRIM(STRING(state),2)+' (arb.)'],nplot)
         102 : ytitle = REPLICATE(['emission (arb.)'],nplot)
         ELSE: ytitle = REPLICATE([''],nplot)
       ENDCASE

       CASE plot.option OF
         100 : labels = REPLICATE([''             ],nplot)
         102 : labels = REPLICATE([''             ],nplot)
         ELSE: labels = REPLICATE(['eirene:divimp'],nplot)
       ENDCASE
       END
;   --------------------------------------------------------------------
    ELSE: BEGIN  
      PRINT, 'ERROR cortex_PlotParallelProfile: Unrecognised plot option'
      PRINT, '  OPTION = ',option,' (',plot.option,')'
      RETURN, -1
;   --------------------------------------------------------------------
      END
  ENDCASE

  colors   = ['Black','Red','Darkgreen','Blue','Orange', 'Purple', 'Silver', 'Hotpink']

  xy_label = [0.96,0.04,0.12,0.88]

  IF (focus OR nplot EQ 1) THEN BEGIN
    xy_label = [0.93,0.07,0.13,0.87]
    charsize_labels = charsize_labels * 1.2
  ENDIF

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
            12: xdata = val.eirene.s[j]
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
                ydata[*,0] =  val.source.par_net[k]
                ydata[*,1] =  val.source.par_ion[k]
                ydata[*,2] = -val.source.par_rec[k]
                ydata[*,3] =  val.source.par_usr[k]
                ydata[*,4] =  val.source.par_ano[k]
                ydata[*,5] =  val.eirene.ion_net[l]
                ydata[*,6] = -val.eirene.rec_net[l]
                END
            6 : ydata[*,0] = val.plasma.dens[j] * val.plasma.vi[j]  
            7 : BEGIN
                ydata[*,0] = val.source.mom_net[k]
                ydata[*,1] = val.source.mom_vol[k]
                ydata[*,2] = val.source.mom_usr[k]
                ydata[*,3] = val.source.mom_ano[k]
                END
            8 : ydata[*,0] = val.plasma.pe[j] + val.plasma.pi[j]
            9 : BEGIN
                IF (0 EQ 1) THEN BEGIN
                  ydata[*,0] = val.eirene.balmer_alpha[k]
                ENDIF

                IF (1 EQ 1) THEN BEGIN
                  ydata[*,0] = val.eirene.qe_net [k]
                  ydata[*,1] = val.source.ene_net[k]
                ENDIF

                END
            10: BEGIN
                IF (0 EQ 1) THEN BEGIN
                  ydata[*,0] = val.eirene.balmer_gamma[k]
                ENDIF

                IF (1 EQ 1) THEN BEGIN
                  ydata[*,0] = val.eirene.qi_net [k]
                  ydata[*,1] = val.source.eni_net[k]
                ENDIF

                END
            11: BEGIN
               IF (0 EQ 1) THEN BEGIN
                 dalpha = val.eirene.balmer_alpha[k]
                 dgamma = val.eirene.balmer_gamma[k]
                 m = WHERE(dgamma EQ 0.0,count)
                 IF (count GT 0            ) THEN dgamma[m] = dalpha[m] ELSE  $
                 IF (count EQ N_ELEMENTS(k)) THEN dgamma[m] = 1.0
                 ydata[*,0] = dalpha / dgamma
               ENDIF

               IF (1 EQ 1) THEN BEGIN
                 ydata[*,0] = val.eirene.balmer_alpha[k]
               ENDIF

               END
            12: BEGIN

               IF (0 EQ 1) THEN BEGIN
                 dalpha = val.eirene.balmer_alpha[k]
                 dgamma = val.eirene.balmer_gamma[k]
                 m = N_ELEMENTS(dalpha) / 2
                 ydata[*,0] = dalpha 
                 ydata[*,1] = dgamma * (dalpha[m] / dgamma[m])
               ENDIF

               IF (1 EQ 1) THEN BEGIN
                 dalpha = val.eirene.balmer_alpha[k]
                 dgamma = val.eirene.balmer_gamma[k]
                 m = WHERE(dgamma EQ 0.0,count)
                 IF (count GT 0            ) THEN dgamma[m] = dalpha[m] ELSE  $
                 IF (count EQ N_ELEMENTS(k)) THEN dgamma[m] = 1.0
                 ydata[*,0] = dalpha / dgamma
               ENDIF

               END
          ENDCASE
          IF (plot.xdata EQ 'tar_dist') THEN BEGIN
            CASE iplot OF
              1 : xdata2 = [0.28,val.plasma.r[j],2.0]
              2 : xdata2 = [0.28,val.plasma.r[j],2.0]
              3 : xdata2 = [0.28,val.plasma.r[j],2.0]
              4 : xdata2 = [0.28,val.plasma.r[j],2.0]
              5 : xdata2 = val.source.r[j]
              6 : xdata2 = val.source.r[j]
              7 : xdata2 = val.source.r[j]
              8 : xdata2 = val.source.r[j]
              9 : xdata2 = val.eirene.r[j]
              10: xdata2 = val.eirene.r[j]
              11: xdata2 = val.eirene.r[j]
              12: xdata2 = val.eirene.r[j]
            ENDCASE             
            ydata2 = ydata
            i = WHERE(xdata LT 4.5)
            xdata = xdata2[i] - 0.28
            ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
            FOR j = 0, MAXNYDATA-1 DO ydata[*,j] = ydata2[i,j]
          ENDIF
          END
;       ----------------------------------------------------------------
        2: BEGIN
          file = val.target.file
          IF (time_slice NE -1) THEN BEGIN          
            labels = MAKE_ARRAY(100,VALUE='',/STRING)
            labels[0] = '                                ' + STRTRIM(STRING(time_slice),2) + ' us'
          ENDIF ELSE BEGIN
            str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
            str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
            labels[iplot-1] = labels[iplot-1] + str[0] + ', ' + STRTRIM(STRING(tube[idata-1]),2) +  ' :'
          ENDELSE
          i = WHERE(val.target.tube EQ tube[idata-1])
          j = WHERE(val.plasma.tube EQ tube[idata-1])
          ntrace = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  ; Number of data lines on each plot
print,'>>>>>>',tube,tube[0]
print,val.target.tube
print,i
;print,  val.target.s[i,0]
;print,  val.plasma.s[j]
;print,  val.target.s[i,1]
help,val.target.s
help,val.plasma.s
          xdata_target    = [val.target.s[i,0],val.plasma.s[j],val.target.s[i,1]]
          xdata_no_target = val.plasma.s[j]
          CASE iplot OF
            1 : xdata = xdata_target
            2 : xdata = xdata_target
            3 : xdata = xdata_no_target
            4 : xdata = xdata_target
            5 : xdata = xdata_target
            6 : xdata = xdata_target
            7 : xdata = xdata_target
            8 : xdata = xdata_target
          ENDCASE
          IF (plot.xdata EQ 'norm') THEN xdata = xdata / MAX(xdata)
          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
          CASE iplot OF
            1 : ydata[*,0] = [val.target.dens[i,0],val.plasma.dens[j],val.target.dens[i,1]]
            2 : ydata[*,0] = [val.target.M   [i,0],val.plasma.M   [j],val.target.M   [i,1]]
            3 : ydata[*,0] =                       val.plasma.vi  [j]                      
            4 : ydata[*,0] = [val.target.pe  [i,0],val.plasma.pe  [j],val.target.pe  [i,1]] 
            5 : ydata[*,0] = [val.target.pi  [i,0],val.plasma.pi  [j],val.target.pi  [i,1]] 
            6 : ydata[*,0] = [val.target.pe  [i,0],val.plasma.pe  [j],val.target.pe  [i,1]] +  $
                             [val.target.pi  [i,0],val.plasma.pi  [j],val.target.pi  [i,1]] 
            7 : ydata[*,0] = [val.target.te  [i,0],val.plasma.te  [j],val.target.te  [i,1]]
            8 : ydata[*,0] = [val.target.ti  [i,0],val.plasma.ti  [j],val.target.ti  [i,1]]

          ENDCASE
;          ntrace = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  ; Number of data lines on each plot
;          xdata_target    = [val.target.s[i,0],val.plasma.s[j],val.target.s[i,1]]
;          xdata_no_target = val.plasma.s[j]
;          CASE iplot OF
;            1 : xdata = xdata_target
;            2 : xdata = xdata_target
;            3 : xdata = xdata_target
;            4 : xdata = xdata_target
;            5 : xdata = xdata_target
;            6 : xdata = xdata_target
;            7 : xdata = xdata_target
;            8 : xdata = xdata_no_target
;            9 : xdata = xdata_no_target
;            10: xdata = xdata_no_target
;            11: xdata = xdata_no_target
;            12: xdata = xdata_no_target
;          ENDCASE
;          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
;          CASE iplot OF
;            1 : ydata[*,0] = [val.target.dens[i,0],val.plasma.dens[j],val.target.dens[i,1]]
;            2 : ydata[*,0] = [val.target.M   [i,0],val.plasma.M   [j],val.target.M   [i,1]]
;            3 : ydata[*,0] = [val.target.pe  [i,0],val.plasma.pe  [j],val.target.pe  [i,1]]
;            4 : ydata[*,0] = [val.target.pi  [i,0],val.plasma.pi  [j],val.target.pi  [i,1]]
;            5 : ydata[*,0] = [val.target.pe  [i,0],val.plasma.pe  [j],val.target.pe  [i,1]] +  $
;                             [val.target.pi  [i,0],val.plasma.pi  [j],val.target.pi  [i,1]] 
;            6 : ydata[*,0] = [val.target.te  [i,0],val.plasma.te  [j],val.target.te  [i,1]]
;            7 : ydata[*,0] = [val.target.ti  [i,0],val.plasma.ti  [j],val.target.ti  [i,1]]
;            8 :
;            9 : 
;            10:
;            11:
;            12:
;          ENDCASE
          END
;       ----------------------------------------------------------------
        3: BEGIN
          file = val.target.file
          str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
          str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
          labels[iplot-1] = labels[iplot-1] + str[0] + ' :'
          i = WHERE(val.ir EQ tube)
          ntrace = [1]
          xdata = [val.p[i]]
          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
          CASE iplot OF
            1 : ydata[*,0] = val.lines[i,state]
          ENDCASE
          END
;       ----------------------------------------------------------------
        999: BEGIN
          IF (idata EQ 1) THEN  $
            subtitle[iplot-1] = subtitle[iplot-1] + ', TUBE '+STRTRIM(STRING(tube[iplot-1]),2)

          CASE (plot.option) OF
;           ------------------------------------------------------------
            100: BEGIN
              ntrace = MAKE_ARRAY(99,/LONG,VALUE=1)
;              ntrace = MAKE_ARRAY(99,/LONG,VALUE=2)
              file = val.divimp.file

              str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
              str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
              labels[iplot-1] = labels[iplot-1] + str[0] + ':'

              i = WHERE(val.divimp.tube EQ tube[iplot-1],count_i)
;              j = WHERE(val.eirene.tube EQ tube[iplot-1],count_j)
              count_j = -1
              IF (count_i EQ 0 OR count_j EQ 0) THEN BEGIN
                PRINT, 'ERROR cortex_PlotParallelProfile: Data not found (1)'
                PRINT, '  OPTION   =',plot.option
                PRINT, '  COUNT_I,J=',count_i,count_j
                PRINT, '  IPLOT    =',iplot
                PRINT, '  TUBE     =',tube[iplot-1]
                RETURN, -1
              ENDIF
              xdata = val.divimp.s[i]
              ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
;              ydata[*,0] = val.divimp.imp_dens[i,state]
              ydata[*,0] = val.divimp.imp_dens[i,state] * val.divimp.div_influx
;              ydata[*,1] = val.eirene.imp_dens[j,state]
              END
;           ------------------------------------------------------------
            101 : BEGIN
              ntrace = MAKE_ARRAY(99,/LONG,VALUE=2)
              file = val.eirene.file
              i = WHERE(val.eirene.tube EQ tube[iplot-1],count_i)
              j = WHERE(val.divimp.tube EQ tube[iplot-1],count_j)
              IF (count_i EQ 0 OR count_j EQ 0) THEN BEGIN
                PRINT, 'ERROR cortex_PlotParallelProfile: Data not found (2)'
                PRINT, '  OPTION   =',plot.option
                PRINT, '  COUNT_I,J=',count_i,count_j
                PRINT, '  IPLOT    =',iplot
                PRINT, '  TUBE     =',tube[iplot-1]
                RETURN, -1
              ENDIF
              xdata = val.eirene.s[i]
              ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
              ydata[*,0] = val.eirene.imp_ioniz[i,state]
              ydata[*,1] = val.divimp.imp_ioniz[j,state] * val.divimp.div_influx
              END
;           ------------------------------------------------------------
            102 : BEGIN

              file = val.emission.file
              str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
              str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
              labels[iplot-1] = labels[iplot-1] + str[0] + ':'

              ntrace = MAKE_ARRAY(99,/LONG,VALUE=1)
              i = WHERE(val.emission.tube EQ tube[iplot-1],count_i)
              IF (count_i EQ 0 ) THEN BEGIN
                PRINT, 'ERROR cortex_PlotParallelProfile: Data not found (3)'
                PRINT, '  OPTION   =',plot.option
                PRINT, '  COUNT_I,J=',count_i,count_j
                PRINT, '  IPLOT    =',iplot
                PRINT, '  TUBE     =',tube[iplot-1]
                RETURN, -1
              ENDIF              
              xdata = val.emission.s[i]
;              xdata = val.emission.p[i]
              ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
              ydata[*,0] = val.emission.lines[i,state-1]
;              ydata[*,0] = ydata[*,0] / MAX(ydata[*,0])
              END
;           ------------------------------------------------------------
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
      IF (plot.flip) THEN BEGIN
        xdata =  MAX(xdata) - xdata
        i = SORT(xdata)
        xdata = xdata[i]
        FOR j = 0, ntrace[iplot-1]-1 DO ydata[*,j] = ydata[i,j]
      ENDIF
      data = { x : xdata, y : ydata, file : file } 
      IF (idata EQ 1) THEN data_store = CREATE_STRUCT(           name,data) ELSE  $
                           data_store = CREATE_STRUCT(data_store,name,data)
      IF (N_ELEMENTS(WHERE(plot.xrange EQ 0.0)) NE 2) THEN BEGIN
        IF (plot.xabsolute) THEN xmax1 = 1.0        ELSE  $
                                 xmax1 = MAX(xdata)
        i = WHERE(xdata/xmax1 GE plot.xrange[0] AND xdata/xmax1 LE plot.xrange[1])        
        IF (N_ELEMENTS(i) EQ 1) THEN BEGIN
          PRINT,'ERROR cortex_PlotParallelProfiles: No data within XRANGE'
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
    IF (plot.xabsolute AND (plot.xrange[0] NE 0.0 OR plot.xrange[1] NE 0.0)) THEN xrange = plot.xrange
    yrange = [ymin,ymax]
    IF (                   (plot.yrange[0] NE 0.0 OR plot.yrange[1] NE 0.0)) THEN yrange = plot.yrange
    position = [xpos[0],ypos[0],xpos[1],ypos[1]]

    plot_type = default_plot_type                                           
    IF (focus NE 0 OR yi EQ plot_yn OR iplot EQ nplot) THEN plot_type = 2   ; Show x-axis label

    ; *** HACK ***
;    IF (iplot EQ 1) THEN yrange = [0.0, 1.0E+21]
;    IF (iplot EQ 2) THEN yrange = [-0.5, 0.5]

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

            9: BEGIN
               IF (1 EQ 1) THEN BEGIN
                 OPLOT, val.x, val.y[*,1], COLOR=TrueColor(colors[1])
               ENDIF    
               END
           10: BEGIN
               IF (1 EQ 1) THEN BEGIN
                 OPLOT, val.x, val.y[*,1], COLOR=TrueColor(colors[1])
               ENDIF    
               END

           12: BEGIN
               OPLOT, val.x, val.y[*,1], COLOR=TrueColor(colors[1])
               END
            ELSE:
          ENDCASE
          END
;       ----------------------------------------------------------------
        2: BEGIN
          IF (focus OR iplot EQ 1) THEN  $
            cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                            dev_xsize, dev_ysize, charsize_labels, colors
          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
          OPLOT, val.x, val.y[*,0], COLOR=TrueColor(colors[idata-1]) 
;          CASE iplot OF
;            4: BEGIN
;               OPLOT, val.x, val.y[*,1], COLOR=TrueColor(colors[1])
;              END
;            ELSE:
;          ENDCASE
          END
;       ----------------------------------------------------------------
        999: BEGIN
          IF (focus OR iplot EQ 1) THEN  $
            cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                            dev_xsize, dev_ysize, charsize_labels, colors
          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
          OPLOT, val.x, val.y[*,0], COLOR=TrueColor(colors[idata-1]) 
;          IF (ntrace[iplot-1] GT 1) THEN  $
;            OPLOT, val.x, val.y[*,1], COLOR=TrueColor(colors[1])
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

  IF (ndata EQ 1) THEN BEGIN
    str = val.file                   ; extract the case name from the data file
    str = STRSPLIT(str,'/',/EXTRACT)
    str = str[N_ELEMENTS(str)-1]     ; take the last sub-string
    str = STRSPLIT(str,'.',/EXTRACT)
    str = str[0]                     ; take the first one

    title = title + ': CASE ' + STRTRIM(str,2) ; + ', TUBE ' + STRTRIM(STRING(tube),2)
  ENDIF

;    title = ''  ; *** HACK ***


  n = STRLEN(title) 
  nmax = LONG(70.0 * 2.0 / (1.5 * charsize))
help,title
help,n
help,nmax
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
    nmax = 125
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


