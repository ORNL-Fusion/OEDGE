;
; ======================================================================
;
FUNCTION cortex_LocatePoint,x,y,p1,s=s,t=t,no_first=no_first,no_last=no_last

;  print,'point location search',N_ELEMENTS(x)

  IF (KEYWORD_SET(no_first)) THEN first_point = 1.0D-06     ELSE first_point = -1.0E-6
  IF (KEYWORD_SET(no_last )) THEN last_point  = 9.99999D-01 ELSE last_point  = 1.000001D0

;print,first_point,last_point,ARG_PRESENT(no_first),ARG_PRESENT(no_last ),KEYWORD_SET(no_first),KEYWORD_SET(no_last )

  ; Find where the point is along the contour:
  n = N_ELEMENTS(x)
  FOR i = 0, n-2 DO BEGIN

    IF (ABS(x[i+1]-x[i]) GT 1.0D-8) THEN s = (p1[0] - x[i]) / (x[i+1] - x[i]) ELSE s = -999.0D
    IF (ABS(y[i+1]-y[i]) GT 1.0D-8) THEN t = (p1[1] - y[i]) / (y[i+1] - y[i]) ELSE t = -999.0D

;   print,'locate',i,n,s,t,x[i+1]-x[i],y[i+1]-y[i],p1[0],x[i],p1[1],y[i],FORMAT='(A,2I6,2F18.12,2F18.12,2X,4F18.12)'

    IF (s EQ -999.0D AND ABS(p1[0]-x[i]) LT 1.0D-8) THEN s = t
    IF (t EQ -999.0D AND ABS(p1[1]-y[i]) LT 1.0D-8) THEN t = s

;    print,'      ',i,n,s,t,x[i+1]-x[i],y[i+1]-y[i],FORMAT='(A,2I6,2F18.12,2F18.12)'


    IF ((ABS(s-t) LT 1.0D-7 AND s GE first_point AND s LE last_point)) THEN BREAK




;        (s EQ -999.0D AND t EQ -999.0D0                )) THEN BREAK 

;    IF ((ABS(s-t) LT 1.0D-7 AND s LT 0.9999999D) OR  $
;        (s EQ -999.0D AND (t GE -1.0D-10 AND t LT 1.0D)) OR  $
;        (t EQ -999.0D AND (s GE -1.0D-10 AND s LT 1.0D)) OR  $
;        (s EQ -999.0D AND t EQ -999.0D0                )) THEN BREAK 
;;        (s EQ -999.0D AND (t GE 0.0D AND t LT 1.0D)) OR  $   
;;        (t EQ -999.0D AND (s GE 0.0D AND s LT 1.0D))) THEN BREAK 
  ENDFOR
  IF (i EQ n-1) THEN BEGIN
    PRINT, 'ERROR grid_LocatePoint: Location of point on line not found'
    STOP
  ENDIF

  result = i

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_SpectrumGeo, spectra, wall, strata, flux, spectrum_str


;help,wall,/struct

  n = N_ELEMENTS(spectrum_str)


  first_spectrum = 1


  FOR i = 0, n-1 DO BEGIN
    FOR j = 1, 100 DO BEGIN

      IF (cortex_CheckIndex(j,spectrum_str[i]) EQ 0) THEN CONTINUE



      k = WHERE(spectra.v_ind EQ j, count)

      IF (count EQ 0) THEN BEGIN
        PRINT,'ERROR cortex_SpectrumGeo: Spectrum wall segment data not found' 
        PRINT,'  SPECTRUM NUMBER=',j
      ENDIF


      v1x = spectra.v_1x[k]
      v1y = spectra.v_1y[k]
      v2x = spectra.v_2x[k] 
      v2y = spectra.v_2y[k] 



      ; find out which segments the points are on

;print,v1x
;print,v1y
;print,v2x
;print,v2y
    
       x = [REFORM(wall.v1[0,*]),wall.v2[0,wall.n-1]]
       y = [REFORM(wall.v1[1,*]),wall.v2[0,wall.n-1]]
;help,x
;help,y

       index_hold = -1

       FOR k = 0, count-1 DO BEGIN

         index = cortex_LocatePoint(x,y,[v2x[k],v2y[k]],/no_last)  ; use the second point since the wall is reversed when calling EIRENE



;print,'--------->',index

         IF (index EQ index_hold) THEN CONTINUE

         IF (k EQ 0) THEN BEGIN
           index_store = index + 1
           vr1 = x[index  ]
           vz1 = y[index  ]
           vr2 = x[index+1]
           vz2 = y[index+1]
         ENDIF ELSE BEGIN
           vr1 = [vr1,x[index  ]]
           vz1 = [vz1,y[index  ]]
           vr2 = [vr2,x[index+1]]
           vz2 = [vz2,y[index+1]]
         ENDELSE

         index_hold = index

       ENDFOR

       index_store = [index_store,index + 1]

       k = SORT(index_store)
       index_store = index_store[k]

;help,vr1


print,'===',i+1,index_store
;print,vr1
;print,vz1
;print,vr2
;print,vz2



         
       data_struct = { spectrum : j, index : index_store, vr1 : vr1, vz1 : vz1, vr2 : vr2, vz2 : vz2 }
       name = 'data' + STRING(j,FORMAT='(I0)')
       IF (first_spectrum EQ 1) THEN  $
         spectrum_geo = CREATE_STRUCT(             name,data_struct) ELSE  $
         spectrum_geo = CREATE_STRUCT(spectrum_geo,name,data_struct)

       first_spectrum = 0


    ENDFOR
  ENDFOR


;help,spectrum_geo,/struct
;stop

  RETURN, spectrum_geo








  FOR i = 0, n-1 DO BEGIN
    FOR j = 1, 100 DO BEGIN
      IF (cortex_CheckIndex(j,spectrum_str[i]) EQ 1) THEN BEGIN

        spectrum = cortex_ExtractStructure(spectra,j)

print, 'local ', i,j, '  ',spectrum_str[i]
;help,spectrum,/struct

         print,'spectra',spectra.cell[j-1]

         istratum = -spectra.cell[j-1] - 1

         print,'strata',strata.target[istratum],strata.range1[istratum],strata.range2[istratum]

         target = strata.target[istratum]
         range = [strata.range1[istratum],strata.range2[istratum]]
         
         IF (target EQ 1) THEN k = WHERE( flux.index_cell EQ 1 AND flux.index_ring EQ range[0], count) ELSE  $
                               k = WHERE( flux.index_cell NE 1 AND flux.index_ring EQ range[0], count)

         print, 'start', target, range[0], k, flux.index_pin[k]

         i1 = flux.index_pin[k]

         IF (target EQ 1) THEN k = WHERE( flux.index_cell EQ 1 AND flux.index_ring EQ range[1], count) ELSE  $
                               k = WHERE( flux.index_cell NE 1 AND flux.index_ring EQ range[1], count)

         print, 'start', target, range[1], k, flux.index_pin[k]

         i2 = flux.index_pin[k]

         index = [i1,i2]
         k = SORT(index)
         index = index[k]

print,'===',index
print,'===',i1,i2

         vr1 = flux.r_vertex1[index[0]-1:index[1]-1]
         vz1 = flux.z_vertex1[index[0]-1:index[1]-1]
         vr2 = flux.r_vertex2[index[0]-1:index[1]-1]
         vz2 = flux.z_vertex2[index[0]-1:index[1]-1]
         
         data_struct = { spectrum : j, index : index, vr1 : vr1, vz1 : vz1, vr2 : vr2, vz2 : vz2 }
         name = 'data' + STRING(j,FORMAT='(I0)')
         IF (first_spectrum EQ 1) THEN  $
           spectrum_geo = CREATE_STRUCT(             name,data_struct) ELSE  $
           spectrum_geo = CREATE_STRUCT(spectrum_geo,name,data_struct)

         first_spectrum = 0

      ENDIF
    ENDFOR
  ENDFOR









  END
;
; ======================================================================
;
FUNCTION cortex_PlotWallProfiles, plot, data_array, grid=grid, wall=wall, annotate=annotate, ps=ps, wdyn=wdyn

  PRINT
  PRINT,'----------------------- NEW PLOT -----------------------'
  PRINT

  MAXNXDATA = 1000
  MAXNYDATA = 20

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

  colors = ['Black','Red','Green','Blue','Orange','Purple','Hotpink', 'Darkseagreen','Silver','Gold','Indigo']
;  colors = ['Black','Red','Blue','Orange','Purple', 'Hotpink', 'Green']

  xy_label = [0.96,0.04,0.12,0.88]

  IF (focus NE 0) THEN BEGIN
    xy_label = [0.93,0.07,0.13,0.87]
    charsize_labels = charsize_labels * 1.2
    !P.CHARSIZE  = charsize * 1.5
  ENDIF

  option = plot.option
  IF (option EQ 100 OR option EQ 101) THEN option = 999

  plot_peak = plot.peak
  plot_sum  = plot.sum

  max_value   = REPLICATE(-1.0E+10,50)
  max_erosion = -1.0E+10
  max_deposit = -1.0E+10
  tot_retain  = max_value

  default_type = 1
  plot_xboarder = 0.05
  plot_yboarder = 0.1
  plot_xspacing = 0.100
  plot_yspacing = 0.025

  save_ypos = [-999.0,-999.0]

  first_plot = 1
  iplot_start = 1

  ; pulse_time = 400.0  ; one ITER pulse  CHANGE
  pulse_time = 60.0 * 60.0 * 24.0 * 365.0  ; one year continuous

  CASE option OF
;   --------------------------------------------------------------------
    1: BEGIN
       nplot = 3
       plot_xn = 1
       plot_yn = [3]
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
       plot_yn = [3]
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
                   'EROSION / mm s-1']
;                   'EROSION / mm per 1000 shots']   ; For the 1000 shot plot
       xtitle = 'wall index'
       ytitle = ['none','none',  $
                 'atom flux density (D m-2 s-1)',        $
                 'average indicent atom energy (eV)',    $
                 'impurity influx (particles m-2 s-1)',  $
                 'erosion rate (mm s-1)']
;                 'erosion rate (mm / 1000 shots)']  ; For the 1000 shot plot
       labels = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace = [1,1,1,1,1,1]  ; Number of data lines on each plot
       IF (plot.show_grid) THEN BEGIN
         plot_xn = 3
         plot_yn = [2]
         plot_type = [-1,3,1,1,1,1]  ; Type of plot
       ENDIF ELSE BEGIN
         plot_xn = 2
         plot_yn = [2]
         plot_type = [ 0,0,1,1,1,1]  
       ENDELSE       
       END
;   --------------------------------------------------------------------
    4: BEGIN
       nplot = 3
       plot_xn = 1
       plot_yn = [3]
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
                   'ION FLUX TO SURF. / m-2 s-1'       ,  $
                   'TARGET ELECTRON TEMP. / eV'    ,  $
                   'EROSION / m-2 s-1'                   ,  $
                   'DEPOSITION (ions) / m-2 s-1' ,  $
;                   'NET DEPOSITION / A-2 s-1'            ,  $  CHANGE
                   'GROSS / mm per year'            ,  $
;                   'NET DEP. / mm (see y-axis)']
;                   'NET DEP. / microns per shot'] ; CHANGE
                   'NET DEP. / mm per year']
       xtitle = 'wall index'
       ytitle = ['none','none',  $
                 'perp. particle flux (m-2 s-1)',  $
                 'Te_target (eV)'         ,  $
                 'erosion (m-2 s-1)'      ,  $
                 'deposition (m-2 s-1)'   ,  $
                 ; 'net deposition (A-2 s-1)'  ,  $  CHANGE
                 'gross (mm per year)'  ,  $
;                 'net deposition (mm s-1)']
;                 'net dep. (mm per 1000 400 s pulses)']
;                 'net dep. (microns per 400 s pulse)']  CHANGE
                 'net dep. (mm per year)']
       labels = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace = [1,1,1,1,1,1,1,1]  ; Number of data lines on each plot
       IF (plot.show_grid) THEN BEGIN
         plot_xn = 4
         plot_yn = [2]
         plot_type = [-1,3,1,1,1,1,1,1]  
       ENDIF ELSE BEGIN
         plot_xn = 3
         plot_yn = [2]
         plot_type = [ 0,0,1,1,1,1,1,1]  
       ENDELSE       
       nplot = 8
       ndata = N_ELEMENTS(TAG_NAMES(data_array))
       IF (ndata EQ 1) THEN BEGIN
         labels[0] = 'ions:atoms'
         labels[2] = 'total:atoms'
         labels[4] = 'erosion:deposition' 
         ntrace[2] = 2
         ntrace[4] = 2
         ntrace[6] = 2  ; CHANGE
       ENDIF
       END
;   --------------------------------------------------------------------
    6: BEGIN

       IF (plot.data_set EQ 0) THEN BEGIN
         PRINT, 'ERROR cortex_PlotWallProfile: Sputtering DATA_SET not specified'
         PRINT, '  OPTION = ',option
         RETURN, -1
       ENDIF

       IF (plot.title EQ 'unknown') THEN  $
         title = 'DIVIMP COMPOSITE WALL EROSION DATA' ELSE  $
         title = plot.title 
       subtitle = ['none','none',  $
                   'BOMBARDING PARTICLE FLUX / m-2 s-1',  $
                   'AVERAGE INDICENT ENERGY / eV'       ,  $
                   'SPUTTERED IMPURITY INFLUX / m-2 s-1',  $
                   'YIELD'                              ]
       xtitle = 'wall segment index'
       ytitle = ['none','none', $
                 'indicent flux (m-2 s-1)',  $
                 'incident energy (eV)'    ,  $
                 'sputtered flux (m-2 s-1)',  $
                 'yield'                   ]
       labels = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace = [1,1,1,1,1,1]  
       IF (plot.show_grid) THEN BEGIN
         plot_xn = 3
         plot_yn = [2]
         plot_type = [-1,3,1,1,1,1]  
       ENDIF ELSE BEGIN
         plot_xn = 2
         plot_yn = [2]
         plot_type = [ 0,0,1,1,1,1]  
       ENDELSE       
       nplot = 6
       ndata = N_ELEMENTS(TAG_NAMES(data_array))
       IF (ndata EQ 1) THEN BEGIN
         labels[0] = 'ions:atoms'
         labels[1] = 'total:atoms'
         labels[2] = 'total:atoms'
         ntrace[3] = 1
         ntrace[4] = 1
         ntrace[5] = 1
       ENDIF
       END
;   --------------------------------------------------------------------
    7: BEGIN  ; background plasma fluxes
       IF (plot.title EQ 'unknown') THEN  $
         title = 'DIVIMP WALL EROSION' ELSE  $
         title = plot.title 
       subtitle = ['none','none',  $
                   'ION FLUX TO SURF. / m-2 s-1'       ,  $
                   'TARGET ELECTRON TEMP. / eV'    ,  $
                   '???     / m-2 s-1'                   ,  $
                   '???        (ions) / m-2 s-1' ,  $
                   '???            / m-2 s-1'               ,  $
                   '???    . / mm (see y-axis)']
       xtitle = 'wall index'
       ytitle = ['none','none',  $
                 'perp. particle flux (m-2 s-1)',  $
                 'Te_target (eV)'         ,  $
                 'erosion (m-2 s-1)'      ,  $
                 'deposition (m-2 s-1)'   ,  $
                 'net deposition (m-2 s-1)'  ,  $
                 'net dep. (mm per 1000 400 s pulses)']
       labels = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace = [1,1,1,1,1,1,1,1]  ; Number of data lines on each plot
       IF (plot.show_grid) THEN BEGIN
         plot_xn = 4
         plot_yn = [2]
         plot_type = [-1,3,1,1,1,1,1,1]  
       ENDIF ELSE BEGIN
         plot_xn = 3
         plot_yn = [2]
         plot_type = [ 0,0,1,1,1,1,1,1]  
       ENDELSE       
       nplot = 8
       ndata = N_ELEMENTS(TAG_NAMES(data_array))
       IF (ndata EQ 1) THEN BEGIN
;         labels[0] = 'ions:atoms'
;         labels[2] = 'total:atoms'
;         ntrace[2] = 2
;         ntrace[4] = 2
       ENDIF
       END
;   --------------------------------------------------------------------
    8: BEGIN

      IF (plot.title EQ 'unknown') THEN title = 'WALLDYN DATA' ELSE title = plot.title 
      subtitle = ['none','none',  $
                  'SURFACE TEMPERATURE / K'  ,  $
                  'NET Be DEP. / atoms m-2 s-1'  ,  $
                  'CO DEP. (50% DT) / T m-1 s-1',  $
                  'NET Be DEP. / microns per shot']
      xtitle = 'wall index'
      ytitle = ['none','none',  $
                'temperature (K)'      ,  $
                'net dep. (Be m-2 s-1)',  $
                'retention (T m-1 s-1)',  $
                'net dep. (microns per 400 s shot)']
      labels = MAKE_ARRAY(100,VALUE=' ',/STRING)
      ntrace = [1,1,1,1,1,1]  ; Number of data lines on each plot
      IF (plot.show_grid) THEN BEGIN
        plot_xn = 3
        plot_yn = [2]
        plot_type = [-1,3,1,1,1,1]
      ENDIF ELSE BEGIN
        plot_xn = 2
        plot_yn = [2]
        plot_type = [ 0, 0,1,1,1,1]  
      ENDELSE       
      nplot = 6
      ndata = N_ELEMENTS(TAG_NAMES(data_array))
      END
;   --------------------------------------------------------------------
    9: BEGIN ; energy resolved atomic fluxes to wall elements

       ; set the number of plots to make
       spectrum_str = STRSPLIT(plot.spectrum_index,' ',/EXTRACT)
       nplot = N_ELEMENTS(spectrum_str)
 
       ; set the number of traces to appear on each plot
       ntrace = REPLICATE(0,nplot)
       FOR i = 0, nplot-1 DO BEGIN
         FOR j = 1, 100 DO IF (cortex_CheckIndex(j,spectrum_str[i]) EQ 1) THEN ntrace[i]++
       ENDFOR

       IF (plot.title EQ 'unknown') THEN  $
         title = 'ENERGY DISTRIBUTION OF ATOMIC FLUX TO THE WALL' ELSE  $
         title = plot.title 

       subtitle = ['INCIDENT ATOM PARTICLE FLUX / m-2 s-1',  $
                   'AVERAGE ATOM INDICENT ENERGY / eV'      ,  $
                   REPLICATE('',nplot)]

       xtitle = ['wall segment index','energy at the centre of each bin (eV)']

       IF (plot.normalize EQ 0) THEN ytitle_dist = 'atoms (s-1)' ELSE ytitle_dist = 'atoms (arb.)'

       ytitle = ['atom flux (m-2 s-1)',  $
                 'avg. atom energy (eV)'   ,  $
                 REPLICATE(ytitle_dist,nplot)]

       labels = MAKE_ARRAY(100,VALUE='',/STRING)

       IF (plot.show_grid) THEN BEGIN
         plot_xn     = 2
         plot_yn     = [4,nplot]
         subtitle    = [REPLICATE('blank',2),subtitle]
         ytitle      = [REPLICATE('blank',2),ytitle]
         plot_type   = [-1,3,1,1,REPLICATE(1,nplot)]
         ntrace      = [ 1,1,1,1,ntrace]
         iplot_start = 3
         nplot       = nplot + 4
       ENDIF ELSE BEGIN
         plot_xn = 1
         plot_yn = nplot
         plot_type = [0,0,1,1,1,1]  
       ENDELSE       



;       ndata = 1 ; N_ELEMENTS(TAG_NAMES(data_array))
;       IF (ndata EQ 1) THEN BEGIN
;         labels[0] = 'ions:atoms'
;         labels[1] = 'total:atoms'
;         labels[2] = 'total:atoms'
;         ntrace[3] = 1
;         ntrace[4] = 1
;         ntrace[5] = 1
;       ENDIF


        ndata = N_ELEMENTS(TAG_NAMES(data_array))

        IF (ndata EQ 1) THEN BEGIN
          str = STRSPLIT(data_array.data1.flux.file,'/',/EXTRACT)                   ; Extract case name to STR
          str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
          title = title + ': ' + str[0] 
        ENDIF


        ; loop over the data set selection and identify the wall segments that are affected

        val_local = cortex_ExtractStructure(data_array,1)
        flux     = val_local.flux
        strata   = val_local.strata
        spectra  = val_local.spectra

;help,data_array,/struct
;help,val_local,/struct
;help,wall,/struct
;help,flux,/struct
;help,strata,/struct
;help,spectra,/struct



        spectrum_geo = cortex_SpectrumGeo(spectra, wall, strata, flux, spectrum_str)






;help,spectrum_geo,/struct
;help,spectrum_geo.data5,/struct
;stop
       END
;   --------------------------------------------------------------------
    10: BEGIN
       nplot = 4
       plot_xn = 1
       plot_yn = [4]
       title = plot.title 
       subtitle = ['PLASMA / MW m-2',  $
                   'CHARGE-EXCANGE ATOMS / MW m-2',  $
                   'LINE EMISSION / MW m-2',  $ 
                   'TOTAL / MW m-2' ]
       xtitle   = 'WALL SEGMENT INDEX'
       ytitle   = ['g_perp (MW m-2)','CX (MW m-2)','rad. (MW m-2)','total (MW m-2)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       labels[2] = '0\hydrogen (solid):0\impurity (dashed)'
       ntrace    = [1,1,2,1] 
       plot_type = [1,1,1,1]  

       END
;   --------------------------------------------------------------------
    ELSE: BEGIN  
      PRINT, 'ERROR cortex_PlotWallProfile: Unrecognised plot option'
      PRINT, '  OPTION = ',option,' (',plot.option,')'
      RETURN, -1
      END
  ENDCASE

;
; Setup plot:
; ----------------------------------------------------------------------
;
;
; Setup plot area:
; ----------------------------------------------------------------------

  frame_bnds = plot.frame_bnds

  IF (focus NE 0) THEN BEGIN
    frame_bnds = [0.0,1.0]
    plot_xn = 1
    plot_yn = [1]
  ENDIF

  size = TOTAL(frame_bnds[1] - frame_bnds[0])

  xsize = (size - 2.0 * plot_xboarder - FLOAT(plot_xn-1) * plot_xspacing) / FLOAT(plot_xn)
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
    IF (yi EQ plot_yn[xi-1]+1) THEN BEGIN
      xi = xi + 1
      yi = 1
    ENDIF

    ysize = (1.0  - 2.0 * plot_yboarder - FLOAT(plot_yn[xi-1]-1) * plot_yspacing) / FLOAT(plot_yn[xi-1]) 
    
    xcen =       xsize * (0.5 + FLOAT(xi-1)) + 1.5 * plot_xboarder + plot_xspacing * FLOAT(xi-1) + frame_bnds[0]
    ycen = 1.0 - ysize * (0.5 + FLOAT(yi-1)) - 1.0 * plot_yboarder - plot_yspacing * FLOAT(yi-1)
    xpos = [xcen - 0.5 * xsize, xcen + 0.5 * xsize]
    ypos = [ycen - 0.5 * ysize, ycen + 0.5 * ysize]

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

      IF (KEYWORD_SET(wall) AND plot.option NE 8) THEN file = wall.file ELSE file = val.wall.file

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

             ; j = [107, 116, 167, 290]
             ; FOR i = 0, 3 DO BEGIN
             ;   k = WHERE(val.wall.index EQ j[i])
             ;   flux = ( val.wall.atom_par_flux[k] + val.wall.mol_par_flux[k]) / 2
             ;   temp = 300 
             ;   fact = 2.0 * SQRT(!PI / 2.0)
             ;   mass = 2.0 * 1.67E-27
             ;   pressure = fact * SQRT( mass ) * flux * SQRT( 1.38E-23 * temp )
             ;   print, j[i], val.wall.atom_par_flux[k],  $
             ;                val.wall.mol_par_flux [k],  $
             ;                flux, temp, pressure, pressure * 7.5,  $
             ;          FORMAT='(I6,2E10.2,2X,E10.2,2X,I6,2X,2E10.2)'                     
             ; ENDFOR
              END
            2: BEGIN
              ydata[*,0] = val.wall.atom_avg_energy
              ;mid_atom_flux = 0.5 * (MIN(val.wall.atom_par_flux) + MAX(val.wall.atom_par_flux))
              ;i = WHERE(val.wall.atom_par_flux LT 0.01 * mid_atom_flux, count)
              ;IF (count GE 1) THEN ydata[i,0] = 0.0
              END
            3: BEGIN
              ydata[*,0] = val.wall.atom_energy_flux  ; W m-2

;              help,val.wall,/struct

              i = WHERE(plot.xrange EQ 0.0, count)
              IF (count EQ 2) THEN i = INDGEN(N_ELEMENTS(val.wall.index)) ELSE  $
                                   i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1], count)        
              IF (count EQ 0) THEN BEGIN
                PRINT,'cortex_PlotWallProfiles','No data within XRANGE'
                STOP
              ENDIF

              print, 'total power:',TOTAL(val.wall.atom_energy_flux[i] *  $
                                          val.wall.toroidal_area   [i]),  $
                                    TOTAL(val.wall.toroidal_area   [i])
 
              END
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
                 6: atom_diameter = ( 2.267   * 1.0E+3              / ( 12.011 * 1.67E-27) )^(-1.0/3.0)  ; graphite
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

              ; print,'atom_diameter',atom_diameter, FIX(val.summary.ion_atomic_number)

              particles_per_m2 = 1.0 / (atom_diameter^2)

              ;              impurity influx           surface atom density     layer thickness
              ;              (atoms / s / m^2)         (atoms / m^2 / layer)    (m / layer)       convert to mm  
              ydata[*,0] = ( val.wall.em_par_atm_2_2 / particles_per_m2     ) * atom_diameter   * 1.0E+3         

              ;                          pulse time   number of pulses 
;              ydata[*,0] =  ydata[*,0] * pulse_time  * 1000.0             ; For the 1000 shot plot

              ; Integrate up to get the total erosion:
              CASE FIX(val.summary.ion_atomic_number) OF
                 4: mass = 9.0122
                 6: mass = 12.011  ; graphite
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

          IF (idata EQ 1) THEN BEGIN
            IF (plot.xorigin EQ 0) THEN BEGIN
              xdata = val.wall.index
            ENDIF ELSE BEGIN
              i = ABS(plot.xorigin) - 1
              xdata = -1.0 * (val.wall.dist - val.wall.dist[i])
              IF (plot.xorigin LT 0) THEN xdata = -1.0 * xdata
              xtitle = 'distance along wall (m)'
            ENDELSE

            ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
          ENDIF

          ; Calculate the distance along the target from the strike-points:
          IF (iplot EQ 3) THEN BEGIN
          ENDIF

          itrace = idata - 1

          CASE iplot OF
            1:
            2:
            3: BEGIN
              ydata[*,itrace] = val.wall.ion_flux
              IF (ndata EQ 1) THEN ydata[*,1] = val.wall.atm_flux  ; Only plotted if data from only one case has been requested
              END
            4: BEGIN
              FOR i = 0, N_ELEMENTS(xdata)-1 DO  $
                IF (val.flux.dens[i] NE -999.0) THEN ydata[i,itrace] = val.flux.te[i] ELSE ydata[i,itrace] = -999.0
              END
            5: BEGIN
              ydata[*,itrace] = val.wall.tot_ero
              IF (ndata EQ 1) THEN ydata[*,1] = val.wall.atm_ero
              END
            6: ydata[*,itrace] = val.wall.tot_dep
            7: BEGIN
               ; ydata[*,itrace] = -val.wall.tot_net * 1.0E-20  ; TOT_NET is the net erosion, but showing net deposition on the plot  CHANGE

              ; Calculate thickness eroded per second:
              CASE FIX(val.wall.atomic_number) OF    
                ;                    Density                          Mass per atom          Assume simple 
                ;          (in m)    (g/cm^3)   (convert to kg/m^3)   (kg / atom)            cubic structure  
                 4: atom_diameter = ( 1.85    * 1.0E+3              / (  9.012 * 1.67E-27) )^(-1.0/3.0)
                 6: atom_diameter = ( 2.267   * 1.0E+3              / ( 12.011 * 1.67E-27) )^(-1.0/3.0)  ; graphite
                26: atom_diameter = ( 7.87    * 1.0E+3              / ( 55.845 * 1.67E-27) )^(-1.0/3.0)
                74: atom_diameter = (19.35    * 1.0E+3              / (183.840 * 1.67E-27) )^(-1.0/3.0)
                ELSE: BEGIN
                  PRINT,'ERROR cortex_PlotWallProfiles: Unrecognised element'
                  PRINT,'  A =',FIX(val.summary.ion_atomic_number)
                  END
              ENDCASE
              particles_per_m2 = 1.0 / (atom_diameter^2)
              ;               impurity influx     surface atom density     layer thickness
              ;               (atoms / s / m^2)   (atoms / m^2 / layer)    (m / layer)       convert to mm  
              ydata[*,0] = ( -val.wall.tot_ero  / particles_per_m2     ) * atom_diameter   * 1.0E+3         
              ydata[*,0] =  ydata[*,0] * pulse_time        
              ydata[*,1] = ( -val.wall.tot_dep  / particles_per_m2     ) * atom_diameter   * 1.0E+3         
              ydata[*,1] =  ydata[*,1] * pulse_time        
              END
            8: BEGIN
              ; Calculate thickness eroded per second:
              CASE FIX(val.wall.atomic_number) OF    
                ;                    Density                          Mass per atom          Assume simple 
                ;          (in m)    (g/cm^3)   (convert to kg/m^3)   (kg / atom)            cubic structure  
                 4: atom_diameter = ( 1.85    * 1.0E+3              / (  9.012 * 1.67E-27) )^(-1.0/3.0)
                 6: atom_diameter = ( 2.267   * 1.0E+3              / ( 12.011 * 1.67E-27) )^(-1.0/3.0)  ; graphite
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

              ; print,'atom_diameter',atom_diameter, FIX(val.summary.ion_atomic_number)

              particles_per_m2 = 1.0 / (atom_diameter^2)

              ; print,'particles_per_m2',particles_per_m2

              ;                   impurity influx     surface atom density     layer thickness
              ;                   (atoms / s / m^2)   (atoms / m^2 / layer)    (m / layer)       convert to mm  
              ydata[*,itrace] = ( -val.wall.tot_net  / particles_per_m2     ) * atom_diameter   * 1.0E+3         

              ;                                    pulse time        convert to microns
              ; ydata[*,itrace] =  ydata[*,itrace] * pulse_time      * 1.0E+3  ; CHANGE
              ydata[*,itrace] =  ydata[*,itrace] * pulse_time                             
              
              ;                                     pulse time   number of pulses 
              ; ydata[*,itrace] =  ydata[*,itrace] * pulse_time * 1000.0             ; For the 1000 shot plot

              ; Integrate up to get the total erosion:
              CASE FIX(val.summary.ion_atomic_number) OF
                 4: mass = 9.0122
                 6: mass = 12.011  ; graphite
                26: mass = 55.845
                74: mass = 183.84
                ELSE: BEGIN
                  PRINT,'ERROR cortex_PlotWallProfiles: Unrecognised element'
                  PRINT,'  A =',FIX(val.summary.ion_atomic_number)
                  END
              ENDCASE

              tot_area = TOTAL(val.wall.area) 

              length   = val.wall.dist_2 - val.wall.dist_1

;              circ = 2.0 * !PI * val.wall.r
              circ = 2.0 * !PI * val.wall.r0                           ; make sure the CIRC array is referenced below if an array is being used!
;              print,'r0',val.wall.r0,MEAN(val.wall.r)
;              print,'r',val.wall.r

              IF (plot.integral[0] NE -999.0) THEN BEGIN
                i = WHERE(xdata GE plot.integral[0] AND xdata LE plot.integral[1], count)
              ENDIF ELSE BEGIN
                i = WHERE(plot.xrange EQ 0.0, count)
                IF (count EQ 2) THEN i = INDGEN(N_ELEMENTS(xdata)) ELSE BEGIN
                  IF ( plot.deposition_only ) THEN  $
                    i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1] AND -val.wall.tot_net GT 0.0, count)  ELSE  $
                    i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1]                             , count)        
                ENDELSE
              ENDELSE
              IF (count EQ 0) THEN BEGIN
                PRINT,'cortex_PlotWallProfiles','No data within XRANGE or INTEGRAL range (5)'
                STOP
              ENDIF

              tot_erosion = TOTAL( val.wall.tot_ero[i] * length[i] * circ) *  $
                            mass * 1.67E-27 * 1000.0                             ; the "1000" is to convert from kg to g
              tot_deposit = TOTAL( val.wall.tot_dep[i] * length[i] * circ) *  $
                            mass * 1.67E-27 * 1000.0
              net_erosion = TOTAL(-val.wall.tot_net[i] * length[i] * circ) *  $  ; -ve sign is because we want net deposition and TOT_NET is erosion
                            mass * 1.67E-27 * 1000.0

              tot_erosion =  tot_erosion * pulse_time  ; * 1000.0  ; For the 1000 shot plot
              tot_deposit =  tot_deposit * pulse_time  ; * 1000.0
              net_erosion =  net_erosion * pulse_time  ; * 1000.0

              max_erosion = MAX([max_erosion,net_erosion])

              print, 'net_erosion = ',idata,net_erosion,max_erosion

              length = val.wall.dist_2 - val.wall.dist_1
              absfac1 = TOTAL(val.wall.ion_ero * length)
              absfac2 = TOTAL(val.wall.atm_ero * length)
              END
          ENDCASE

          IF (ndata EQ 1) THEN BEGIN
            ; Store the data for a big dump
            i = WHERE(plot.xrange EQ 0.0, count)
            IF (count EQ 2) THEN i = INDGEN(N_ELEMENTS(xdata)) ELSE  $
                                 i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1], count)        
            CASE iplot OF
              3: BEGIN
                 dump_data = CREATE_STRUCT(           'xdata' , xdata          [i]) 
                 dump_data = CREATE_STRUCT(dump_data, 'index' , val.wall.index [i]) 
                 dump_data = CREATE_STRUCT(dump_data, 'dens'  , val.flux.dens  [i]) 
                 dump_data = CREATE_STRUCT(dump_data, 'length', val.wall.length[i]) 
                 dump_data = CREATE_STRUCT(dump_data, 'r'     , val.wall.r     [i]) 
                 dump_data = CREATE_STRUCT(dump_data, 'psin'  , val.flux.psin  [i]) 
                 dump_data = CREATE_STRUCT(dump_data, 'rho'   , val.flux.rho   [i]) 
                 dump_data = CREATE_STRUCT(dump_data, 'ti'    , val.flux.ti    [i]) 
                 dump_data = CREATE_STRUCT(dump_data, 'jsat'  , val.flux.jsat  [i])
                 END
              4: BEGIN
                 dump_data = CREATE_STRUCT(dump_data, 'te'    , val.flux.te[i])
                 END
              7: BEGIN 
                 dump_data = CREATE_STRUCT(dump_data, 'net_dep', -val.wall.tot_net[i] )  
                 fp1 = 3
                 FREE_LUN,fp1        
                 file_name = plot.case_name[0] + '.deposition_' + plot.id 
                 OPENW,fp1,'cortex_deposition/'+file_name, ERROR=err
                 IF (err NE 0) THEN BEGIN
                   PRINT,'ERROR cortex_PlotWallProfiles: Problem opening data stream (to file)'
                   PRINT,'  FILE_NAME= ',file_name
                   PRINTF,-2,!err.msg
                   FREE_LUN,fp1
                   STOP
                 ENDIF
                 PRINTF,fp1,'*',FORMAT='(A)'
                 PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
                 PRINTF,fp1,'* WALL DEPOSITION DATA FROM DIVIMP'
                 PRINTF,fp1,'*',FORMAT='(A)'
                 PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
                 PRINTF,fp1,'* CASE            ',plot.case_name[0],FORMAT='(2A)'
                 PRINTF,fp1,'* TITLE           ',FORMAT='(A)'
                 PRINTF,fp1,'* DATE AND TIME   ',FORMAT='(A)'
                 PRINTF,fp1,'*',FORMAT='(A)'
                 PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
                 PRINTF,fp1,'{DATA FILE VERSION}',FORMAT='(A)'
                 PRINTF,fp1,'     1.0',FORMAT='(A)'
                 PRINTF,fp1,'*',FORMAT='(A)'
                 PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
                 PRINTF,fp1,'{ATOMIC NUMBER OF IMPURITY}',FORMAT='(A)'
                 PRINTF,fp1,FIX(val.summary.ion_atomic_number)
                 PRINTF,fp1,'*',FORMAT='(A)'
                 PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
                 PRINTF,fp1,'{WALL SEGMENT DATA}',FORMAT='(A)'
                 PRINTF,fp1,N_ELEMENTS(dump_data.xdata),FORMAT='(I12)'
                 PRINTF,fp1,'*',FORMAT='(A)'
                 PRINTF,fp1,'* dist    - distance from strike-point (usually -- check rho value for dist=0.0 to confirm)',FORMAT='(A)'
                 PRINTF,fp1,'* length  - length of wall segment',FORMAT='(A)'
                 PRINTF,fp1,'* R       - radius of the centre of the wall segment',FORMAT='(A)'
                 PRINTF,fp1,'* jsat    - parallel ion saturation current to the wall/target segment',FORMAT='(A)'
                 PRINTF,fp1,'* n_e     - electron density at the entrance to the sheath',FORMAT='(A)'
                 PRINTF,fp1,'* T_e     - electron temperature at the entrance to the sheath',FORMAT='(A)'
                 PRINTF,fp1,'* T_i     - background hydrogenic ion temperature at the entrance to the sheath',FORMAT='(A)'
                 PRINTF,fp1,'* net_dep - net impurity deposition flux density',FORMAT='(A)'
                 PRINTF,fp1,'* rho     - distance of the wall/target segment from the separatrix, mapped to the at the outer midplane (approx.)',FORMAT='(A)'
                 PRINTF,fp1,'* psi_n   - normalized magnetic flux coordinate',FORMAT='(A)'
                 PRINTF,fp1,'*',FORMAT='(A)'
                 PRINTF,fp1,'* ','index','dist','length','R','jsat' ,'n_e','T_e','T_i','net_dep','rho','psi_n',FORMAT='(A2,A6,2X,5A10,2A6,A13,2A10)'
                 PRINTF,fp1,'* ',' '    ,'m'   ,'m'     ,'m','A m-2','m-3','eV' ,'eV' ,'par m-2 s-1','m'  ,' ',FORMAT='(A2,A6,2X,5A10,2A6,A13,2A10)'
                 j = SORT(dump_data.xdata)
                 FOR i = 0, N_ELEMENTS(dump_data.xdata)-1 DO BEGIN
                   IF (dump_data.jsat[j[i]] EQ -999.0) THEN CONTINUE
                   PRINTF,fp1,  $
                     dump_data.index [j[i]],  $
                     dump_data.xdata [j[i]],  $
                     dump_data.length [j[i]],  $
                     dump_data.r      [j[i]],  $
                     dump_data.jsat   [j[i]],  $
                     dump_data.dens   [j[i]],  $
                     dump_data.te     [j[i]],  $
                     dump_data.ti     [j[i]],  $
                     dump_data.net_dep[j[i]],  $
                     dump_data.rho    [j[i]],  $
                     dump_data.psin   [j[i]],  $
                     FORMAT= '(2X,I6,2X,3F10.6,2E10.2,2F6.1,E13.2,2F10.5)'
                 ENDFOR
                 CLOSE,fp1
                 FREE_LUN,fp1
                 END
              ELSE:
            ENDCASE
            ENDIF

          END
;       ----------------------------------------------------------------
        6: BEGIN

          nstates = 0
          sputter = val.sputter
          data_set = cortex_ExtractStructure(sputter,plot.data_set)
          FOR istate = 1, data_set.max_charge DO BEGIN
            IF (cortex_CheckIndex(istate,plot.states)) THEN nstates++
          ENDFOR
          IF (nstates EQ 0) THEN BEGIN
            PRINT, 'ERROR cortex_PlotWallProfile: No valid charge states selected'
            PRINT, '  OPTION = ',plot.option
            PRINT, '  IDATA  = ',idata
            RETURN, -1
          ENDIF

print,'nstates',nstates,data_set.max_charge

          labels[iplot-1] = ''

          ntrace = [1,1,nstates,nstates,nstates,nstates]  

          xdata = FINDGEN(N_ELEMENTS(sputter.index)) + 1.0  ; sputter.index
          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),nstates,/FLOAT,VALUE=0.0)      
          zdata = ydata

          IF (iplot EQ 3) THEN net_flux = 0.0

          itrace = -1
          FOR istate = 1, data_set.max_charge DO BEGIN
            IF (NOT cortex_CheckIndex(istate,plot.states)) THEN CONTINUE
            itrace++
            labels[iplot-1] = labels[iplot-1] + '+' + STRTRIM(STRING(istate),2) + ':'
            CASE iplot OF
              1: 
              2:
              3: BEGIN
                zdata[*,itrace] = data_set.flux[*,itrace] * sputter.absfac_in[plot.data_set-1]
                ydata[*,itrace] = zdata[*,itrace] * sputter.modifier
                END
              4: BEGIN
                zdata[*,itrace] = data_set.e0  [*,itrace]
                ydata[*,itrace] = zdata[*,itrace] * sputter.modifier 
                END
              5: BEGIN
                ydata[*,itrace] = sputter.modifier * data_set.flux [*,itrace] * sputter.absfac_in[plot.data_set-1] *  $
                                                     data_set.yield[*,itrace]
                net_flux = net_flux + TOTAL( ydata[*,itrace] * sputter.length ) ; * 2.0 * !PI * sputter.pos_r )
                END
              6: BEGIN
                ydata[*,itrace] = sputter.modifier * data_set.yield[*,itrace]
                END
            ENDCASE
          ENDFOR  

          ; Store the data for a big dump
          i = WHERE(plot.xrange EQ 0.0, count)
          IF (count EQ 2) THEN i = INDGEN(N_ELEMENTS(xdata)) ELSE  $
                               i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1], count)        
          CASE iplot OF
            3: BEGIN
              dump_data = CREATE_STRUCT(          'xdata'        ,xdata[i]            ) 
              dump_data = CREATE_STRUCT(dump_data,'index'        ,val.erosion.index[i]) 
              dump_data = CREATE_STRUCT(dump_data,'incident_flux',zdata[i,0:itrace]   )
              

;               dump_data = CREATE_STRUCT(dump_data, 'costet', val.erosion.costet[i]       ) 
;               dump_data = CREATE_STRUCT(dump_data, 'bratio', val.erosion.bratio[i]       ) 
;               dump_data = CREATE_STRUCT(dump_data, 'dens'  , val.flux.dens  [i]       ) 
;               dump_data = CREATE_STRUCT(dump_data, 'length', val.erosion.length[i]       ) 
;               dump_data = CREATE_STRUCT(dump_data, 'r'     , val.erosion.r     [i]       ) 
;               dump_data = CREATE_STRUCT(dump_data, 'r'     , val.flux.r_cen [i]       ) 
;               dump_data = CREATE_STRUCT(dump_data, 'z'     , val.flux.z_cen [i]       ) 
;               dump_data = CREATE_STRUCT(dump_data, 'psin'  , val.flux.psin  [i]       ) 
;               dump_data = CREATE_STRUCT(dump_data, 'rho'   , val.flux.rho   [i]       ) 
;               dump_data = CREATE_STRUCT(dump_data, 'ti'    , val.flux.ti    [i]       ) 
;               dump_data = CREATE_STRUCT(dump_data, 'jsat'  , val.flux.jsat    [i]  )
               END
            4: BEGIN
              dump_data = CREATE_STRUCT(dump_data, 'incident_energy' , zdata[i,0:itrace]   )
               END
            6: BEGIN 

               fp1 = 3
               FREE_LUN,fp1        
               file_name = plot.case_name[0] + '.walldyn'
               IF (plot.id NE 'unknown') THEN file_name = file_name + '_' + plot.id
               OPENW,fp1,'cortex_deposition/'+file_name, ERROR=err
               IF (err NE 0) THEN BEGIN
                 PRINT,'ERROR cortex_PlotWallProfiles: Problem opening data stream (to file)'
                 PRINT,'  FILE_NAME= ',file_name
                 PRINTF,-2,!err.msg
                 FREE_LUN,fp1
                 STOP
               ENDIF
      	     
               n = N_ELEMENTS(dump_data.xdata)

               zeros = MAKE_ARRAY(itrace+1,/FLOAT,VALUE=0.0)      

               CASE sputter.bomb_z[plot.data_set-1] OF
                 1 : BEGIN
                   CASE FIX(sputter.bomb_a[plot.data_set-1]) OF
                     1: element = 'H'
                     2: element = 'D'
                     3: element = 'T'
                     ELSE: BEGIN
                       PRINT,'ERROR cortex_PlotWallProfiles: Hydrogen isotope not recognized'
                       PRINT,'  A = ',INT(sputter.bomb_a[plot.data_set-1])
                       STOP
                       END
                     ENDCASE
                   END
                 2 : element = 'He'
                 10: element = 'Ne'
                 ELSE: BEGIN
                   PRINT,'ERROR cortex_PlotWallProfiles: Atomic number not recognized'
                   PRINT,'  Z = ',sputter.bomb_z[plot.data_set-1]
                   STOP
                   END
               ENDCASE

               ; Map surface temperature profile data onto the wall:

               min_ring = MIN(ABS(val.flux.index_ring))
               j_inner = WHERE(val.flux.index_ring EQ min_ring AND val.flux.index_cell EQ 1)
               j_outer = WHERE(val.flux.index_ring EQ min_ring AND val.flux.index_cell NE 1)

               dist_inner = val.flux.dist - (val.flux.dist[j_inner])[0]  ; No idea why this ()[0] is necessary, but it is...
               dist_outer = val.flux.dist - (val.flux.dist[j_outer])[0]

               flux_temp = MAKE_ARRAY(N_ELEMENTS(val.flux.index),/FLOAT,VALUE=400.0)

               print,'temp_inner',plot.temp_inner

               IF (plot.temp_inner NE 'none') THEN BEGIN
                 help,val.temp_inner,/struct

                 ; Set distances so that the target is near the centre of the full range of distances:
                 delta =             dist_inner[N_ELEMENTS(dist_inner)-1] - dist_inner     [0] +  $
                         0.5 * (val.flux.length[N_ELEMENTS(dist_inner)-1] + val.flux.length[0])
                 j = WHERE(dist_inner LT -0.5 * delta)  
                 dist_inner[j] = dist_inner[j] + delta  

                 j = WHERE(dist_inner GE MIN(val.temp_inner.dist) AND  $
                           dist_inner LE MAX(val.temp_inner.dist))

                 FOR k = 0, N_ELEMENTS(j)-1 DO  $
                   flux_temp[j[k]] = INTERPOL(val.temp_inner.temp,val.temp_inner.dist,dist_inner[j[k]]) + 273.0

                 ;print,val.temp_inner.dist
                 ;print,val.temp_inner.temp
                 ;FOR k = 0, N_ELEMENTS(j)-1 DO  $
                 ;  print,dist_inner[j[k]],flux_temp[j[k]]
               ENDIF

               IF (plot.temp_outer NE 'none') THEN BEGIN
                 delta =             dist_outer[N_ELEMENTS(dist_outer)-1] - dist_outer     [0] +  $
                         0.5 * (val.flux.length[N_ELEMENTS(dist_outer)-1] + val.flux.length[0])
                 j = WHERE(dist_outer LT -0.5 * delta)
                 dist_outer[j] = dist_outer[j] + delta

                 ; Change the direction of the along-the-wall distance for the outer target, since
                 ; the temperature data (at the moment) always has points in the SOL as +ve and
                 ; PFZ points as -ve, which is fine for the inner target, but not for the outer:
                 dist_outer = -dist_outer

                 j = WHERE(dist_outer GE MIN(val.temp_outer.dist) AND  $
                           dist_outer LE MAX(val.temp_outer.dist))
                 FOR k = 0, N_ELEMENTS(j)-1 DO  $
                   flux_temp[j[k]] = INTERPOL(val.temp_outer.temp,val.temp_outer.dist,dist_outer[j[k]]) + 273.0

                 ;print,val.temp_outer.dist
                 ;print,val.temp_outer.temp
                 ;FOR k = 0, N_ELEMENTS(j)-1 DO  $
                 ;  print,dist_outer[j[k]],flux_temp[j[k]]
               ENDIF

               PRINTF,fp1,'*',FORMAT='(A)'
               PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
               PRINTF,fp1,'* DIVIMP-WALLDYN data file (via CORTEX)'
               PRINTF,fp1,'*',FORMAT='(A)'
               PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
               PRINTF,fp1,'* CASE            ',plot.case_name[0],FORMAT='(2A)'
               PRINTF,fp1,'* TITLE           ',FORMAT='(A)'
               PRINTF,fp1,'* DATE AND TIME   ',FORMAT='(A)'
               PRINTF,fp1,'*',FORMAT='(A)'
               PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
               PRINTF,fp1,'{DATA FILE VERSION}',FORMAT='(A)'
               PRINTF,fp1,'     1.1',FORMAT='(A)'
               PRINTF,fp1,'*',FORMAT='(A)'
               PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
               PRINTF,fp1,'{DATA}'
               PRINTF,fp1,n,FORMAT='(I12)'
               PRINTF,fp1,'*',FORMAT='(A)'
               PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
               PRINTF,fp1,'*',FORMAT='(A)'					       
               PRINTF,fp1,'* index   - wall segment index in DIVIMP'       
               PRINTF,fp1,'* r1,z1   - starting point of segment in the R,Z plane'
               PRINTF,fp1,'* r2,z2   - end point'			       
               PRINTF,fp1,'* T_surf  - surface temperature'			       
               PRINTF,fp1,'* flux_D+ - perpendicular (not parallel!) background ion flux density on wall (multiply by 2*PI*R*delta'
               PRINTF,fp1,'*           to get the total flux to the vessel, where R is midpoint of the wall segment and ''delta'' '
               PRINTF,fp1,'*           is the length)'
               PRINTF,fp1,'* T_e     - electron temperature'		
               PRINTF,fp1,'* T_i     - ion temperature'			
               PRINTF,fp1,'* flux_D  - atom flux density from CX and the dissociation of D2 (as calculated by EIRENE)'
               PRINTF,fp1,'* T_D     - average energy of atom flux'
               PRINTF,fp1,'* E_dist  - atom energy distribution index (corresponding to the distribution data included farther' 
               PRINTF,fp1,'*           down in the file), with -1 indicating that the distribution data for that segment is not'
               PRINTF,fp1,'*           available'
               PRINTF,fp1,'*',FORMAT='(A)'
               PRINTF,fp1,'* index','r1' ,'z1' ,'r2' ,'z2' ,'T_surf','flux_D+'  ,'T_e' ,'T_i' ,'flux_D'   ,'T_D' ,'E_dist',  $
                          FORMAT='(A7,1X,2(2A9,1X),A8,2X,A10,2A8,A12,A10,2X,A6)'
               PRINTF,fp1,'*      ','(m)','(m)','(m)','(m)',   '(K)','(m-2 s-1)','(eV)','(eV)','(m-2 s-1)','(eV)',  $
                          FORMAT='(A7,1X,2(2A9,1X),A8,2X,A10,2A8,A12,A10,2X,6X)'

               ; Check through the energy distribution spectra for D, if available:

               edist = MAKE_ARRAY(N_ELEMENTS(val.flux.index),/INT,VALUE=0)
               IF (plot.spectra EQ 1) THEN BEGIN
                 elist = [-1]
                 FOR j = 0, N_ELEMENTS(val.spectra.spectra.cell)-1 DO BEGIN
                   k = -val.spectra.spectra.cell[j]
                   k1 = val.spectra.strata.range1[k-1]
                   k2 = val.spectra.strata.range2[k-1]
                   IF (val.spectra.strata.target[k-1] EQ 1) THEN  $
                     l = WHERE(val.flux.index_ring GE k1 AND val.flux.index_ring LE k2 AND  $
                               val.flux.index_cell EQ 1, count)  $
                   ELSE  $
                     l = WHERE(val.flux.index_ring GE k1 AND val.flux.index_ring LE k2 AND  $
                               val.flux.index_cell GT 1, count)
                   IF (count EQ 0) THEN BEGIN
                     PRINT,'ERROR cortex_PlotWallProfiles: Wall elements not identified'
                     PRINT,'  k    = ,k
                     PRINT,'  k1,k2= ,k1,k2
                   ENDIF
                   k = WHERE(edist[l] NE 0, count)
                   IF (count NE 0) THEN BEGIN
                     PRINT,'ERROR cortex_PlotWallProfiles: Double assingment of a wall segment to a CX energy distribution'
                     PRINT,'  l        = ,l
                     PRINT,'  edist[l] = ,edist[l]
                   ENDIF
                   edist[l] = j + 1
                   elist = [elist,l]
                 ENDFOR
               ENDIF
        
               FOR j = 0, n-1 DO BEGIN
                 k = i[j] 
                 IF (val.flux.index_target[k] GT 0) THEN  $
                   PRINTF,fp1,  $
                     dump_data.index    [j] ,  $
                     val.flux.r_vertex1 [k] ,  $
                     val.flux.z_vertex1 [k] ,  $
                     val.flux.r_vertex2 [k] ,  $
                     val.flux.z_vertex2 [k] ,  $
                     flux_temp          [k] ,  $
                     val.flux.dens[k]*ABS(val.flux.vb[k])*val.flux.costet[k]*val.flux.bratio[k],  $
                     val.flux.te [k] ,  $
                     val.flux.ti [k] ,  $
                     val.flux.atom_par_flux  [k] ,  $
                     val.flux.atom_avg_energy[k] ,  $
                     edist[k], $ ; k, $
                   FORMAT= '(I7,1X,2(2F9.5,1X),F8.0,2X,E10.2,2F8.2,E12.2,F10.2,2X,I6)'  $
                 ELSE  $
                   PRINTF,fp1,  $
                     dump_data.index [j] , $ 
                     val.flux.r_vertex1 [k] ,  $
                     val.flux.z_vertex1 [k] ,  $
                     val.flux.r_vertex2 [k] ,  $
                     val.flux.z_vertex2 [k] ,  $
                     flux_temp          [k] ,  $
                     0.0,  $
                     0.0,  $
                     0.0,  $
                     val.flux.atom_par_flux  [k] ,  $
                     val.flux.atom_avg_energy[k] ,  $
                     edist[k],  $
                   FORMAT= '(I7,1X,2(2F9.5,1X),F8.0,2X,E10.2,2F8.2,E12.2,F10.2,2X,I6)'
               ENDFOR
               IF (plot.spectra EQ 1) THEN BEGIN
                 n = N_ELEMENTS(val.spectra.spectra.cell)
                 PRINTF,fp1,'*',FORMAT='(A)'
                 PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
                 PRINTF,fp1,'{NUMBER OF CX ENERGY SPECTRA}',FORMAT='(A)'
                 PRINTF,fp1,n,FORMAT='(I6)'                 
                 PRINTF,fp1,'*',FORMAT='(A)'					       
                 PRINTF,fp1,'* energy    - the energy at the centre of the bin'
                 PRINTF,fp1,'* flux      - the particle flux to the wall segment for each energy bin'
                 PRINTF,fp1,'* flux_norm - just the flux for each bin divided by the integral over all bins'
                 PRINTF,fp1,'* stdev     - not available at the moment'
                 FOR k = 1, n DO BEGIN
                   spectrum = cortex_ExtractStructure(val.spectra,k)
                   PRINTF,fp1,'*',FORMAT='(A)'
                   PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
                   PRINTF,fp1,'* spectrum no., number of energy bins, min. range, max. range =',FORMAT='(A)'
                   PRINTF,fp1,k,N_ELEMENTS(spectrum.bin),spectrum.min_value,spectrum.max_value,FORMAT='(I6,I8,2F12.2)'
                   PRINTF,fp1,'*',FORMAT='(A)'
                   PRINTF,fp1,'*','bin','energy','  flux' ,'flux_norm','stdev',FORMAT='(A1,A6,A12,3A14)'
                   PRINTF,fp1,'*','   ','  (eV)','(D s-1)','         ','     ',FORMAT='(A1,A6,A12,3A14)'
                   total_flux = TOTAL(spectrum.flux)

                   print,total_flux/1.602E-19,elist[k],val.flux.atom_par_flux[elist[k]] , val.flux.atom_par_flux[elist[k]]/(total_flux/1.602E-19)

                   FOR j = 0, N_ELEMENTS(spectrum.bin)-1 DO BEGIN
                     PRINTF,fp1,j+1,spectrum.bin[j],spectrum.flux[j]/1.602E-19,spectrum.flux[j]/total_flux,spectrum.stdev[j],FORMAT='(1X,I6,F12.2,3E14.4)'
                   ENDFOR
                 ENDFOR
               ENDIF
               PRINTF,fp1,'*',FORMAT='(A)'
               PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
               PRINTF,fp1,'{NUMBER OF IMPURITY SPECIES}'
               PRINTF,fp1,1,FORMAT='(I12)'
               PRINTF,fp1,'*',FORMAT='(A)'
               PRINTF,fp1,'* For each impurity, the flux to the wall segment, the incident energy (after sheath acceleration), and'
               PRINTF,fp1,'* the velocity of the impurity at the entrance to the sheath are provided.  (Although at the moment, '
               PRINTF,fp1,'* not the veloicty.)'
               n = N_ELEMENTS(dump_data.xdata)
               imp_tag = element + '+0'
               FOR j = 1, itrace+1 DO imp_tag = [imp_tag,'+'+STRTRIM(STRING(j),2)]
               spacer = '          '
               FOR j = 1, 74 DO spacer = [spacer,'          ']
               FOR k = 1, 1 DO BEGIN
                 PRINTF,fp1,'*',FORMAT='(A)'
                 PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
                 PRINTF,fp1,'* species no., atomic number, atomic mass, max. charge state ='
                 PRINTF,fp1,k,sputter.bomb_z[plot.data_set-1],sputter.bomb_a[plot.data_set-1],sputter.bomb_z[plot.data_set-1],FORMAT='(I6,I6,F8.1,I6)'

                 PRINTF,fp1,'*',FORMAT='(A)'
                 PRINTF,fp1,'* index',imp_tag[0:itrace+1],imp_tag[0:itrace+1],imp_tag[0:itrace+1],  $
                            FORMAT='(A7,2X,222A10)'
                 PRINTF,fp1,'*      ','(m-2 s-1)',spacer[0:itrace],'(eV)',spacer[0:itrace],'(m s-1)',spacer[0:itrace],  $
                            FORMAT='(A7,2X,222A10)'

                 FOR j = 0, n-1 DO BEGIN
                   k = i[j] 
                   IF (val.flux.index_target[k] GT 0) THEN  $
                     PRINTF,fp1,  $
                       dump_data.index    [j] ,  $
                       0.0,dump_data.incident_flux  [j,0:itrace],  $
                       0.0,dump_data.incident_energy[j,0:itrace],  $
                       0.0,zeros                    [  0:itrace],  $
                     FORMAT= '(I7,2X,222E10.2)'  $
                   ELSE  $
                     PRINTF,fp1,  $
                       dump_data.index [j] , $ 
                       0.0,zeros[0:itrace],  $
                       0.0,zeros[0:itrace],  $
                       0.0,zeros[0:itrace],  $
                     FORMAT= '(I7,2X,222E10.2)'
                 ENDFOR

               ENDFOR




               CLOSE,fp1
               FREE_LUN,fp1
               END
            ELSE:
          ENDCASE


          END
;       ----------------------------------------------------------------
        7: BEGIN   

          IF (plot.xorigin EQ 0) THEN BEGIN
            xdata = val.wall.index
          ENDIF ELSE BEGIN
            i = ABS(plot.xorigin) - 1
            xdata = -1.0 * (val.wall.dist - val.wall.dist[i])
            IF (plot.xorigin LT 0) THEN xdata_ref = -1.0 * xdata_ref  
            xtitle = 'distance along wall (m)'
          ENDELSE

          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
          CASE iplot OF
            1:
            2:
            3: BEGIN
              ydata[*,0] = val.wall.ion_flux
              END
            4: BEGIN
              FOR i = 0, N_ELEMENTS(xdata)-1 DO  $
                IF (val.flux.dens[i] NE -999.0) THEN ydata[i,0] = val.flux.te[i] ELSE ydata[i,0] = -999.0
              END
            5: BEGIN
              ydata[*,0] = val.wall.ion_flux
              END
            6: BEGIN
              ydata[*,0] = val.wall.ion_flux
              END
            7: BEGIN
              ydata[*,0] = val.wall.ion_flux
              END
            8: BEGIN
              ydata[*,0] = val.wall.ion_flux
              END
          ENDCASE
          ; Store the data for a big dump
          i = WHERE(plot.xrange EQ 0.0, count)
          IF (count EQ 2) THEN i = INDGEN(N_ELEMENTS(xdata)) ELSE  $
                               i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1], count)        
          END
;       ----------------------------------------------------------------
        8: BEGIN   
;          print,val.wall.file
;          help,dwyn,/struct
;          help,dwyn.netbedep,/struct
;          print,dwyn.netbedep.run

          FOR i = 0, wdyn.netbedep.n-1 DO BEGIN
            IF (STRPOS(wdyn.netbedep.run[i],val.wall.file) NE -1) THEN BREAK 
          ENDFOR
          IF (i EQ wdyn.netbedep.n) THEN BEGIN
            PRINT, 'ERROR cortex_PlotWallProfile: WALLDYN case not found'
            PRINT, 'CASE = ',val.wall.file
          ENDIF

          netbedep = REFORM(wdyn.netbedep.data[i,*])
          codep    = REFORM(wdyn.codep   .data[i,*])
          walltemp = REFORM(wdyn.walltemp.data[i,*])

          IF (idata EQ 1) THEN BEGIN
            xdata  = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      
            xdata1 = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      
            ydata  = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      
          ENDIF

          ; Choose wall:
          wdyn_wall = wdyn.wall_0903

          n = wdyn_wall.n

          itrace = idata - 1

          IF (plot.xorigin EQ 0) THEN BEGIN
            FOR j = 0, n-1 DO BEGIN 
              xdata [3*j  ,itrace] = FLOAT(j) - 0.5
              xdata [3*j+1,itrace] = FLOAT(j)
              xdata [3*j+2,itrace] = FLOAT(j) + 0.5
              xdata1[  j  ,itrace] = xdata[3*j+1,itrace]
            ENDFOR
          ENDIF ELSE BEGIN
            i = ABS(plot.xorigin) - 1
            FOR j = 0, n-1 DO BEGIN 
              xdata [3*j  ,itrace] = -1.0 * (wdyn_wall.dist[j] - 0.5 * wdyn_wall.length[j] - wdyn_wall.dist[i])
              xdata [3*j+1,itrace] = -1.0 * (wdyn_wall.dist[j]                             - wdyn_wall.dist[i])
              xdata [3*j+2,itrace] = -1.0 * (wdyn_wall.dist[j] + 0.5 * wdyn_wall.length[j] - wdyn_wall.dist[i])
              xdata1[  j  ,itrace] = xdata[3*j+1,itrace]
            ENDFOR
            i = WHERE(xdata [*,itrace] NE -999.0,count)
            IF (plot.xorigin LT 0 AND count NE 0) THEN xdata [i,itrace] = -1.0 * xdata [i,itrace]  
            i = WHERE(xdata1[*,itrace] NE -999.0,count)
            IF (plot.xorigin LT 0 AND count NE 0) THEN xdata1[i,itrace] = -1.0 * xdata1[i,itrace]  
            xtitle = 'distance along wall (m)'
          ENDELSE

          CASE iplot OF
            1:
            2:
            3: BEGIN
              FOR j = 0, n-1 DO BEGIN 
                ydata[3*j  ,itrace] = walltemp[j]
                ydata[3*j+1,itrace] = walltemp[j]
                ydata[3*j+2,itrace] = walltemp[j]
              ENDFOR
              END
            4: BEGIN
              FOR j = 0, n-1 DO BEGIN 
                ydata[3*j  ,itrace] = netbedep[j]
                ydata[3*j+1,itrace] = netbedep[j]
                ydata[3*j+2,itrace] = netbedep[j]
              ENDFOR

              atomic_number = 4

              ; Calculate thickness eroded per second:
              CASE atomic_number OF    
                ;                    Density                          Mass per atom          Assume simple 
                ;          (in m)    (g/cm^3)   (convert to kg/m^3)   (kg / atom)            cubic structure  
                 4: atom_diameter = ( 1.85    * 1.0E+3              / (  9.012 * 1.67E-27) )^(-1.0/3.0)
                 6: atom_diameter = ( 2.267   * 1.0E+3              / ( 12.011 * 1.67E-27) )^(-1.0/3.0)  ; graphite
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
              
              particles_per_m2 = 1.0 / (atom_diameter^2)

              ;           impurity influx     surface atom density     layer thickness
              ;           (atoms / s / m^2)   (atoms / m^2 / layer)    (m / layer)       convert to mm  
              tot_bedep = netbedep          / particles_per_m2       * atom_diameter   * 1.0E+3         

              ;                       pulse time   convert to microns
              tot_bedep = tot_bedep * pulse_time * 1.0E+3
              
              ; Integrate up to get the total erosion:
              CASE atomic_number OF
                 4: mass = 9.0122
                 6: mass = 12.011  ; graphite
                26: mass = 55.845
                74: mass = 183.84
                ELSE: BEGIN
                  PRINT,'ERROR cortex_PlotWallProfiles: Unrecognised element'
                  PRINT,'  A =',FIX(val.summary.ion_atomic_number)
                  END
              ENDCASE

              length = wdyn_wall.length

              ; circ = 2.0 * !PI * 0.5 * (wdyn_wall.R_start + wdyn_wall.R_end)   
              circ = 2.0 * !PI * 6.0                                           ; make sure the CIRC array is referenced below if an array is being used!

              IF (plot.integral[0] NE -999.0) THEN BEGIN
                i = WHERE(xdata1[*,itrace] GE plot.integral[0] AND xdata1[*,itrace] LE plot.integral[1], count)
              ENDIF ELSE BEGIN
                i = WHERE(plot.xrange EQ 0.0, count)
                IF (count EQ 2) THEN i = WHERE(xdata1[*,itrace] NE -999.0        ,count) ELSE  $
                                     i = WHERE(xdata1[*,itrace] GE plot.xrange[0] AND  $
                                               xdata1[*,itrace] LE plot.xrange[1],count)        
              ENDELSE
              IF (count EQ 0) THEN BEGIN
                PRINT,'cortex_PlotWallProfiles: No data within XRANGE or INTEGRAL range (8/5)'
                PRINT,'  XRANGE = ',plot.xrange
                PRINT,'  XDATA1 = ',xdata1[*,itrace]
                STOP
              ENDIF
              
              tot_deposit = TOTAL( netbedep[i] * length[i] * circ ) *  $
                            mass * 1.67E-27 * 1000.0                      ; the 1000 converts from kg to g

              tot_migrate = TOTAL( ABS(netbedep[i]) * length[i] * circ ) *  $
                            mass * 1.67E-27 * 1000.0

              tot_deposit =  tot_deposit * pulse_time  ; * 1000.0
              tot_migrate =  tot_migrate * pulse_time  ; * 1000.0

              max_deposit = MAX([max_deposit,tot_deposit])

              print, 'net_/max_deposit,tot_migrate = ',tot_deposit,max_deposit,tot_migrate

              END
            5: BEGIN
              FOR j = 0, n-1 DO BEGIN 
                ydata[3*j  ,itrace] = codep[j]
                ydata[3*j+1,itrace] = codep[j]
                ydata[3*j+2,itrace] = codep[j]
              ENDFOR

              circ = 2.0 * !PI * 6.0                                           ; make sure the CIRC array is referenced below if an array is being used!

              IF (plot.integral[0] NE -999.0) THEN BEGIN
                i = WHERE(xdata1[*,itrace] GE plot.integral[0] AND xdata1[*,itrace] LE plot.integral[1], count)
              ENDIF ELSE BEGIN
                i = WHERE(plot.xrange EQ 0.0, count)
                IF (count EQ 2) THEN i = WHERE(xdata1[*,itrace] NE -999.0        ,count) ELSE  $
                                     i = WHERE(xdata1[*,itrace] GE plot.xrange[0] AND  $
                                               xdata1[*,itrace] LE plot.xrange[1],count)        
              ENDELSE
              IF (count EQ 0) THEN BEGIN
                PRINT,'cortex_PlotWallProfiles: No data within XRANGE or INTEGRAL range (8/6)'
                PRINT,'  XRANGE = ',plot.xrange
                PRINT,'  XDATA1 = ',xdata1[*,itrace]
                STOP
              ENDIF
              
              mass = 3.0

              tot_retain[itrace] = TOTAL( codep[i] * circ ) * mass * 1.67E-27 * 1000.0  ; the 1000 converts from kg to g
                            
              tot_retain[itrace] = tot_retain[itrace] * pulse_time * 0.5  ; the 0.5 is for 50% DT mix

              print, 'tot_retain = ',tot_retain[itrace]

              END
            6: BEGIN
              FOR j = 0, n-1 DO BEGIN 
                ydata[3*j  ,itrace] = netbedep[j]
                ydata[3*j+1,itrace] = netbedep[j]
                ydata[3*j+2,itrace] = netbedep[j]
              ENDFOR

              atomic_number = 4

              ; Calculate thickness eroded per second:
              CASE atomic_number OF    
                ;                    Density                          Mass per atom          Assume simple 
                ;          (in m)    (g/cm^3)   (convert to kg/m^3)   (kg / atom)            cubic structure  
                 4: atom_diameter = ( 1.85    * 1.0E+3              / (  9.012 * 1.67E-27) )^(-1.0/3.0)
                 6: atom_diameter = ( 2.267   * 1.0E+3              / ( 12.011 * 1.67E-27) )^(-1.0/3.0)  ; graphite
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
              
              particles_per_m2 = 1.0 / (atom_diameter^2)

              i = WHERE(ydata[*,itrace] NE -999.0,count)
              IF (count EQ 0) THEN BEGIN
                PRINT, 'ERROR cortex_PlotWallProfile: YDATA not found'
                RETURN, -1
              ENDIF

              ;                 impurity influx    surface atom density     layer thickness
              ;                 (atoms / s / m^2)  (atoms / m^2 / layer)    (m / layer)       convert to mm  
              ydata[i,itrace] = ydata[i,itrace]  / particles_per_m2       * atom_diameter   * 1.0E+3         

              ;                                   pulse time   convert to microns
              ydata[i,itrace] = ydata[i,itrace] * pulse_time * 1.0E+3

              END
          ENDCASE

          END
;       ----------------------------------------------------------------
        9: BEGIN

;help,val,/struct
;help,val.flux,/struct
;help,val.spectra,/struct
;help,val.spectra.strata,/struct
;help,val.spectra.spectra,/struct

print,iplot,ntrace[iplot-1]

          xdata = MAKE_ARRAY(10000,10,/FLOAT,VALUE=-999.0)      
          ydata = xdata

;          labels[iplot-1] = labels[iplot-1] + '+' + STRTRIM('shit',2) + ':'

          flux = val.flux

          itrace = 0

          n = N_ELEMENTS(flux.index)

          CASE iplot OF
            iplot_start+0: BEGIN
              xdata[0:n-1,itrace] = flux.index
              ydata[0:n-1,itrace] = flux.atom_par_flux
              END
            iplot_start+1: BEGIN
              xdata[0:n-1,itrace] = flux.index
              ydata[0:n-1,itrace] = flux.atom_avg_energy
              i = WHERE(flux.atom_par_flux LT 1.E-08*MAX(flux.atom_par_flux), count)
;              i = WHERE(flux.atom_par_flux LT 1.E-04*MAX(flux.atom_par_flux), count)
              IF (count GE 1) THEN ydata[i,itrace] = 0.0
              END
            ELSE: BEGIN
              IF (iplot GT iplot_start) THEN BEGIN
                ; find the next spectrum to plot
                j = iplot-iplot_start-2

                print, val.flux.file
                print, 'idata ',idata, j, '  ', spectrum_str[j]

                IF (iplot EQ iplot_start+2) THEN subtitle[iplot-1] = 'ENERGY DISTRIBUTION' ELSE subtitle[iplot-1] = ''

                IF (ndata NE 1 AND iplot EQ iplot_start+2) THEN subtitle[iplot-1] = subtitle[iplot-1] + ', '

                IF (ndata NE 1                           ) THEN subtitle[iplot-1] = subtitle[iplot-1] + 'DATA SET '+ spectrum_str[j] 

                IF (ndata GT 1 AND iplot EQ iplot_start+2) THEN BEGIN
                  str = STRSPLIT(val.flux.file,'/',/EXTRACT)                   ; Extract case name to STR
                  str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
                  labels[iplot-1] = labels[iplot-1] + str[0] + ':'
                ENDIF

                itrace = 0
                FOR i = 1, 100 DO BEGIN
                  IF (cortex_CheckIndex(i,spectrum_str[j]) EQ 1) THEN BEGIN
                    IF (ndata EQ 1) THEN BEGIN
                      labels[iplot-1] = labels[iplot-1] + STRING(i,FORMAT='(I0)') + ':'
                    ENDIF
                    spectrum = cortex_ExtractStructure(val.spectra,i)
                    n = N_ELEMENTS(spectrum.bin)
                    xdata[0:n-1,itrace] = spectrum.bin
                    ydata[0:n-1,itrace] = spectrum.flux
                    IF (plot.normalize EQ 1) THEN ydata[0:n-1,itrace] = ydata[0:n-1,itrace] / MAX(ydata[0:n-1,itrace]) ELSE  $
                                                  ydata[0:n-1,itrace] = ydata[0:n-1,itrace] / 1.602E-19
                    itrace++
                  ENDIF
                ENDFOR

                print,itrace,ntrace[iplot-1]

              ENDIF ELSE BEGIN
                xdata[0,itrace] = -1.0
              ENDELSE
            END
          ENDCASE

          IF (ndata EQ 1) THEN BEGIN

            fp1 = 3
            FREE_LUN,fp1        
            file_name = plot.case_name[0] + '.wall_energy_spectra'
            IF (plot.id NE 'unknown') THEN file_name = file_name + '_' + plot.id
            OPENW,fp1,'cortex_data/'+file_name, ERROR=err
            IF (err NE 0) THEN BEGIN
              PRINT,'ERROR cortex_PlotWallProfiles: Problem opening data stream (to file)'
              PRINT,'  FILE_NAME= ',file_name
              PRINTF,-2,!err.msg
              FREE_LUN,fp1
              STOP
            ENDIF
      	    
            n = N_ELEMENTS(flux.index)

            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
            PRINTF,fp1,'* EIRENE energy spectra for atomic fluxes to surfaces'
            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
            PRINTF,fp1,'* CASE            ',plot.case_name[0],FORMAT='(2A)'
            PRINTF,fp1,'* TITLE           ',FORMAT='(A)'
            PRINTF,fp1,'* DATE AND TIME   ',FORMAT='(A)'
            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
            PRINTF,fp1,'{DATA FILE VERSION}',FORMAT='(A)'
            PRINTF,fp1,'     1.0',FORMAT='(A)'
            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
            PRINTF,fp1,'{WALL DATA}'
            PRINTF,fp1,n,FORMAT='(I12)'
            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
            PRINTF,fp1,'*',FORMAT='(A)'					       
            PRINTF,fp1,'* index   - wall segment index in EIRENE'       
            PRINTF,fp1,'* r1,z1   - starting point of segment in the R,Z plane'
            PRINTF,fp1,'* r2,z2   - end point'			       
            PRINTF,fp1,'* T_surf  - surface temperature'			       
            PRINTF,fp1,'* flux_D+ - perpendicular (not parallel!) background ion flux density on wall (multiply by 2*PI*R*delta'
            PRINTF,fp1,'*           to get the total flux to the vessel, where R is midpoint of the wall segment and ''delta'' '
            PRINTF,fp1,'*           is the length)'
            PRINTF,fp1,'* T_e     - electron temperature'		
            PRINTF,fp1,'* T_i     - ion temperature'			
            PRINTF,fp1,'* flux_D  - atom flux density from CX and the dissociation of D2 (as calculated by EIRENE)'
            PRINTF,fp1,'* T_D     - average energy of atom flux'
            PRINTF,fp1,'* E_dist  - atom energy distribution index (corresponding to the distribution data included farther' 
            PRINTF,fp1,'*           down in the file), with -1 indicating that the distribution data for that segment is not'
            PRINTF,fp1,'*           available'
            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* index','r1' ,'z1' ,'r2' ,'z2' ,'T_surf','flux_D+'  ,'T_e' ,'T_i' ,'flux_D'   ,'T_D' ,'E_dist',  $
                       FORMAT='(A7,1X,2(2A9,1X),A8,2X,A10,2A8,A12,A10,2X,A6)'
            PRINTF,fp1,'*      ','(m)','(m)','(m)','(m)',   '(K)','(m-2 s-1)','(eV)','(eV)','(m-2 s-1)','(eV)',  $
                       FORMAT='(A7,1X,2(2A9,1X),A8,2X,A10,2A8,A12,A10,2X,6X)'

            ; Check through the energy distribution spectra for D, if available:

            flux = val.flux

            ; check if surface assigned to a spectrum 

            str = '1-'+STRTRIM(STRING(N_ELEMENTS(spectra.index)),2)
print, 'str ',str
;help,spectra,/struct
;help,strata,/struct
            spec_geo = spectrum_geo
;            spec_geo = cortex_SpectrumGeo(spectra, wall, strata, flux, str)
            spec_str = TAG_NAMES(spec_geo)
            spec_ind = MAKE_ARRAY(n,VALUE=-1,/LONG)
            FOR j = 0, N_ELEMENTS(spec_str)-1 DO BEGIN
              ind = WHERE(spec_str EQ STRUPCASE(spec_str[j]))
              geo = spec_geo.(ind[0]) 
              spec_ind[geo.index[0]-1:geo.index[1]-1] = geo.spectrum
            ENDFOR

            FOR j = 0, n-1 DO BEGIN
              IF (flux.index_target[j] GT 0) THEN  $
                PRINTF,fp1,  $
                  j+1                ,  $
                  flux.r_vertex1 [j] ,  $
                  flux.z_vertex1 [j] ,  $
                  flux.r_vertex2 [j] ,  $
                  flux.z_vertex2 [j] ,  $
                  -1.0               ,  $
                  flux.dens[j]*ABS(flux.vb[j])*flux.costet[j]*flux.bratio[j],  $
                  flux.te [j]             ,  $
                  flux.ti [j]             ,  $
                  flux.atom_par_flux  [j] ,  $
                  flux.atom_avg_energy[j] ,  $
                  spec_ind[j]       , $
                FORMAT= '(I7,1X,2(2F9.5,1X),F8.0,2X,E10.2,2F8.2,E12.2,F10.2,2X,I6)'  $
              ELSE  $
                PRINTF,fp1,  $
                  j+1                    ,  $ 
                  val.flux.r_vertex1 [j] ,  $
                  val.flux.z_vertex1 [j] ,  $
                  val.flux.r_vertex2 [j] ,  $
                  val.flux.z_vertex2 [j] ,  $
                  -1.0                   ,  $
                  0.0,  $
                  0.0,  $
                  0.0,  $
                  val.flux.atom_par_flux  [j] ,  $
                  val.flux.atom_avg_energy[j] ,  $
                  spec_ind[j]           ,  $
                FORMAT= '(I7,1X,2(2F9.5,1X),F8.0,2X,E10.2,2F8.2,E12.2,F10.2,2X,I6)'
            ENDFOR


            ; output energy distribution data
            spectra  = val.spectra
            spectrum = spectra.data1

            n = N_ELEMENTS(spectra.index)

            spectrum_data = REPLICATE(spectrum,n)

            FOR i = 0, n-1 DO BEGIN
               spectrum = cortex_ExtractStructure(spectra,i+1)
               spectrum_data[i] = spectrum
            ENDFOR

            PRINTF,fp1,'*',FORMAT='(A)'					       
            PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
            PRINTF,fp1,'{CHARGE EXCHANGE ATOM WALL ENERGY DISTRIBUTION DATA}',FORMAT='(A)'
            PRINTF,fp1,'* number of spectra, number of energy bins, min. range (eV), max. range (eV) =',FORMAT='(A)'
            PRINTF,fp1,n,N_ELEMENTS(spectrum.bin),spectrum.min_value,spectrum.max_value,FORMAT='(I6,I8,2F12.2)'
            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'* energy    - the energy at the centre of the bin'
            PRINTF,fp1,'* flux      - the particle flux to the wall segment for each energy bin'
;            PRINTF,fp1,'* flux_norm - just the flux for each bin divided by the integral over all bins'
;            PRINTF,fp1,'* stdev     - not available at the moment'
            PRINTF,fp1,'*',FORMAT='(A)'
            PRINTF,fp1,'*','bin','energy','  flux' ,FORMAT='(A1,A6,A10,2X,A10)'
            PRINTF,fp1,'*','   ','  (eV)','(D s-1)',FORMAT='(A1,A6,A10,2X,A10)'
	    
;              total_flux = TOTAL(spectrum.flux)
	    
;              print,total_flux/1.602E-19,elist[k],val.flux.atom_par_flux[elist[k]] , val.flux.atom_par_flux[elist[k]]/(total_flux/1.602E-19)
	    
            FOR j = 0, N_ELEMENTS(spectrum.bin)-1 DO BEGIN
	    
              PRINTF,fp1,j+1,  $
                spectrum_data[0    ].bin[j]            ,  $
                spectrum_data[0:n-1].flux[j]/1.602E-19 ,  $
                FORMAT='(1X,I6,F10.1,2X,100E10.2)'
	    
;              PRINTF,fp1,j+1,spectrum.bin[j],spectrum.flux[j]/1.602E-19,spectrum.flux[j]/total_flux,spectrum.stdev[j],FORMAT='(1X,I6,F12.2,3E14.4)'
	    
	    
            ENDFOR
	    
            flux_total = MAKE_ARRAY(n,VALUE=0.0,/FLOAT)
            ene_avg    = flux_total
            FOR j = 0, n-1 DO BEGIN
              flux_total[j] = TOTAL(spectrum_data[j].flux) / 1.602E-19
              ene_avg   [j] = TOTAL(spectrum_data[j].flux * spectrum_data[j].bin) / (flux_total[j] * 1.602E-19)
            ENDFOR
            ; totals	    
            PRINTF,fp1,  $
              flux_total[0:n-1] ,  $
              FORMAT='(1X,6X,10X,2X,100E10.2)'
            PRINTF,fp1,  $
              ene_avg[0:n-1] ,  $
              FORMAT='(1X,6X,10X,2X,100F10.2)'
	     
	     
            CLOSE,fp1
            FREE_LUN,fp1

          ENDIF

          END
;       ----------------------------------------------------------------
        10: BEGIN
          xdata = FIX(val.wall.index)
          ydata = MAKE_ARRAY(N_ELEMENTS(xdata),MAXNYDATA,/FLOAT,VALUE=0.0)      
          CASE iplot OF
            1: ydata[*,0] = val.wall.g_perp / 1.0E+6  ; MW m-2 
            2: ydata[*,0] = val.wall.atom_energy_flux ; MW m-2
            3: BEGIN
              ydata[*,0] = val.wall.rad_hydrogen / 1.0E+6 ; MW m-2
              ydata[*,1] = val.wall.rad_impurity / 1.0E+6 ; MW m-2
              END
            4: BEGIN
              ydata[*,0] = (val.wall.g_perp       / 1.0E+6) + (val.wall.atom_energy_flux         ) +  $
                           (val.wall.rad_hydrogen / 1.0E+6) + (val.wall.rad_impurity     / 1.0E+6)
              
              i = WHERE(plot.xrange EQ 0.0, count)
              IF (count EQ 2) THEN i = INDGEN(N_ELEMENTS(val.wall.index)) ELSE  $
                                   i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1], count)        
              IF (count EQ 0) THEN BEGIN
                PRINT,'cortex_PlotWallProfiles','No data within XRANGE'
                STOP
              ENDIF

              print, 'total power:',TOTAL(val.wall.atom_energy_flux[i] *  $
                                          val.wall.toroidal_area   [i]),  $
                                    TOTAL(val.wall.toroidal_area   [i])
 
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

      IF ((plot_type[iplot-1] EQ 1 AND N_ELEMENTS(WHERE(plot.xrange EQ 0.0)) NE 2) AND  $
          NOT (option EQ 9 AND iplot GT iplot_start + 1) ) THEN BEGIN
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
      IF (ndata EQ 1) THEN  $
        FOR j = 1, ntrace[iplot-1]-1 DO ydata1 = [ydata1,REFORM(ydata[i,j])]  $
      ELSE  $
        FOR j = 1, ndata-1 DO ydata1 = [ydata1,REFORM(ydata[i,j])] 
      j = WHERE(ydata1 NE -999.0,count)
      IF (count NE 0) THEN BEGIN
        ymin = MIN([ymin,ydata1[j]])
        ymax = MAX([ymax,ydata1[j]])
      ENDIF
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
    IF (focus NE 0 OR yi EQ plot_yn[xi-1] OR iplot EQ nplot) THEN type = 2   ; Show x-axis label

    IF (plot_type[iplot-1] EQ 3) THEN type = 3

    first_plot = 0

print,xrange
print,yrange
print,ytitle[iplot-1]

    CASE (type) OF
      1: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTICKFORMAT='(A1)',/NOERASE,YLOG=0
      2: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTITLE=xtitle[xi-1],/NOERASE,YLOG=0                         
      3: BEGIN  
        plot_grid        = plot
        plot_grid.size   = 0.3 ; 0.6
        plot_grid.center = [0.6*xpos[0]+0.4*xpos[1],0.5*(ypos[0]+ypos[1])]
        status = cortex_PlotFluidGrid(plot_grid, grid, wall, 0, annotate, 'subordinate', 'outline', ps='on',  $
                                      spectrum_geo=spectrum_geo, spectrum_colors=colors, spectrum_ndata=ndata)
        END
    ENDCASE

    IF (type EQ 3) THEN CONTINUE

;   Write sub-title for each plot:
    IF (option NE 9) THEN  $
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

          ; Show the peak average quantity over a specificed plot range:
          IF ((iplot EQ 1 OR iplot EQ 3) AND plot.show_avg AND cortex_GetValues(plot.xmark,values)) THEN BEGIN
            IF (idata EQ 1) THEN BEGIN
              val_wall = cortex_ExtractStructure(data_array,idata)
              max_average = MAKE_ARRAY(N_ELEMENTS(values),/FLOAT,VALUE=0.0)
            ENDIF
            print, '       >>>>>---------<<<<<'
            FOR i = 0, N_ELEMENTS(values)-1, 2 DO BEGIN
              j = WHERE(xdata GE values[i] AND xdata LE values[i+1])
              print,'trying',values[i:i+1]
              print,'check ',xdata[j]
              print,'check ',val_wall.wall.length[j]
              print,'ydata ',val.y[j,0]
              weight = val_wall.wall.length[j] / TOTAL(val_wall.wall.length[j])
              print,'weight',weight
              print,'weight',TOTAL(weight)
              print,'-->   ',TOTAL(val.y[j,0]*weight)
              max_average[i] = MAX([max_average[i],TOTAL(val.y[j,0]*weight)])
              print,'max_avg=',i,idata,max_average[i]
            ENDFOR
            IF (idata EQ ndata) THEN BEGIN
              FOR i = 0, N_ELEMENTS(values)-1, 2 DO  $
                XYOUTS, MIN([0.75*xmax,values[i]]) + 0.005 * (xmax - xmin),  $
                        0.0                        - 0.075 * (ymax - ymin),  $
                        'MAX_AVG = '+STRTRIM(STRING(MAX(max_average[i]),FORMAT='(E10.2)'),2),  $
                        CHARSIZE=charsize_labels
            ENDIF
          ENDIF

          BREAK
          END
;       ----------------------------------------------------------------
        2: BEGIN
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

          ; Show the peak average quantity over a specificed plot range:
          IF (iplot EQ 1 AND plot.show_avg AND cortex_GetValues(plot.xmark,values)) THEN BEGIN
            IF (idata EQ 1) THEN BEGIN
              val_wall = cortex_ExtractStructure(data_array,idata)
              max_average = MAKE_ARRAY(N_ELEMENTS(values),/FLOAT,VALUE=0.0)
            ENDIF
            FOR i = 0, N_ELEMENTS(values)-1, 2 DO BEGIN
              j = WHERE(xdata GE values[i] AND xdata LE values[i+1])
              weight = val_wall.wall.length[j] / TOTAL(val_wall.wall.length[j])
              max_average[i] = MAX([max_average[i],TOTAL(val.y[j,0]*weight)])
            ENDFOR
            IF (idata EQ ndata) THEN BEGIN
              FOR i = 0, N_ELEMENTS(values)-1, 2 DO  $
                XYOUTS, values[i] + 0.005 * (xmax - xmin),  $
                        0.0       - 0.075 * (ymax - ymin),  $
                        'MAX_AVG = '+STRTRIM(STRING(MAX(max_average[i]),FORMAT='(E10.2)'),2),  $
                        CHARSIZE=charsize_labels
            ENDIF
          ENDIF

          BREAK
          END
;       ----------------------------------------------------------------
        3: 
;       ----------------------------------------------------------------
        5: BEGIN
          IF (idata EQ 1) THEN  $
            cortex_DrawKey, iplot-2, focus, labels, xy_label, xpos, ypos,  $
                            dev_xsize, dev_ysize, charsize_labels, colors, step = 0.065

          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
          IF (plot_peak EQ 1 AND idata EQ ndata) THEN thick = 2.0
          IF (plot_sum  EQ 1 AND idata EQ ndata) THEN thick = 2.0

          IF (plot.integral[0] NE -999.0) THEN  $
             OPLOT, plot.integral, [0.0,0.0], COLOR=TrueColor('Black'), PSYM=6, SYMSIZE=0.75 

          i = WHERE(val.y[*,0] NE -999.0,count_i)
          IF (count_i GT 0) THEN  $
            OPLOT, val.x[i], val.y[i,idata-1], COLOR=TrueColor(colors[idata-1]), THICK=thick

          i = WHERE(val.x GE xmin AND val.x LE xmax)
          max_value[iplot] = MAX([max_value[iplot],val.y[i,idata-1]])

          IF (idata EQ ndata) THEN BEGIN
            ;XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
            ;        (xy_label[3] * ypos[0] + xy_label[2] * ypos[1]) * dev_ysize,  $
            ;XYOUTS, (0.38        * xpos[0] + 0.62        * xpos[1]) * dev_xsize,  $
            ;        (xy_label[2] * ypos[0] + xy_label[3] * ypos[1]) * dev_ysize,  $
            XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
                    (0.97        * ypos[0] + 0.03        * ypos[1]) * dev_ysize,  $
                    'MAX = '+STRTRIM(STRING(max_value[iplot]),2), CHARSIZE=charsize_labels, /DEVICE
            IF (option EQ 5 AND iplot EQ 8 AND (ndata EQ 1 OR plot.integral[0] NE -999.0)) THEN BEGIN
              IF (plot.deposition_only) THEN tag = 'DEP = ' ELSE tag = 'NET = '
              XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
                      (0.92        * ypos[0] + 0.08        * ypos[1]) * dev_ysize,  $
                      ; tag+STRTRIM(STRING(max_erosion),2) + ' g in ' + STRTRIM(STRING(FIX(pulse_time)),2) + ' s',  $  CHANGE
                      tag + STRTRIM(STRING(max_erosion / 1000.0),2) + ' kg ',  $
                      CHARSIZE=charsize_labels, /DEVICE
            ENDIF
          ENDIF

          IF (option EQ 5 AND ndata EQ 1) THEN BEGIN
;            print,'trying..',iplot
            SWITCH iplot OF
              3: 
              5: 
              7: BEGIN
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
        6: BEGIN
          IF (iplot EQ 3 AND idata EQ 1) THEN BEGIN
;print, 'labels ',labels[iplot-1]
;print, ypos,dev_ysize
            step = 0.1
            IF (ntrace[iplot-1] GT 5) THEN step = 0.07
            cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                            dev_xsize, dev_ysize, charsize_labels, colors, step=step
            OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
          ENDIF 

;          print, 'ntrace',ntrace[iplot-1],iplot

          FOR itrace = 0, ntrace[iplot-1]-1 DO  $
            OPLOT, val.x, val.y[*,itrace], COLOR=TrueColor(colors[itrace]), THICK=thick

          CASE iplot OF
            1: 
            2: 
            3: 
            ELSE:
          ENDCASE
          BREAK
          END
;       ----------------------------------------------------------------
        7: BEGIN
         IF (idata EQ 1) THEN  $
           cortex_DrawKey, iplot-2, focus, labels, xy_label, xpos, ypos,  $
                           dev_xsize, dev_ysize, charsize_labels, colors

          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 

          i = WHERE(val.y[*,0] NE -999.0,count_i)
          IF (count_i GT 0) THEN  $
            OPLOT, val.x[i], val.y[i,0], COLOR=TrueColor(colors[idata-1]), THICK=thick

          BREAK
          END
;       ----------------------------------------------------------------
        8: BEGIN
          IF (idata EQ 1) THEN  $
            cortex_DrawKey, iplot-2, focus, labels, xy_label, xpos, ypos,  $
                            dev_xsize, dev_ysize, charsize_labels, colors,  step = 0.070

          itrace = idata - 1

          i = WHERE(val.y[*,0] NE -999.0,count)
          IF (count GT 0) THEN  $
            OPLOT, val.x[i,itrace], val.y[i,itrace], COLOR=TrueColor(colors[itrace]), THICK=thick

          i = WHERE(val.x[i,itrace] GE xmin AND val.x[i,itrace] LE xmax)
          max_value[iplot] = MAX([max_value[iplot],val.y[i,itrace]])

          IF (idata EQ ndata) THEN BEGIN
            XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
                    (0.97        * ypos[0] + 0.03        * ypos[1]) * dev_ysize,  $
                    'MAX_VAL = '+STRTRIM(STRING(max_value[iplot]),2), CHARSIZE=charsize_labels, /DEVICE

            IF (iplot EQ 4) THEN  $
              XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
                      (0.88        * ypos[0] + 0.12        * ypos[1]) * dev_ysize,  $
                      'MAX_DEP = ' + STRTRIM(STRING(max_deposit),2) + ' g in ' + STRTRIM(STRING(FIX(pulse_time)),2) + ' s',  $
                      CHARSIZE=charsize_labels, /DEVICE

          ENDIF

          IF (iplot EQ 5) THEN BEGIN

            IF (itrace EQ 0) THEN suffix = ' in ' + STRTRIM(STRING(FIX(pulse_time)),2) + ' s' ELSE suffix = ' '

            frac = FLOAT(itrace)

            print,'frac',frac,itrace,idata

            XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
                    ((0.22 +        frac  * 0.07      ) * ypos[0] +   $
                     (0.78 -        frac  * 0.07      ) * ypos[1]) * dev_ysize,  $
                    STRTRIM(STRING(tot_retain[itrace]),2) + ' g' + suffix,  $
                    CHARSIZE=charsize_labels, /DEVICE, COLOR=Truecolor(colors[itrace])

          ENDIF

          BREAK
          END
;       ----------------------------------------------------------------
        9: BEGIN

          IF (iplot LE iplot_start+1 AND idata EQ 1) THEN BEGIN
            ; loop over the data set selection and identify the wall segments that are affected
            icolour = -1

            color_str = 'LightGrey' 

            FOR i = 0, N_ELEMENTS(spectrum_str)-1 DO BEGIN

print,'spectrum=================== ',spectrum_str[i]

              IF (ndata GT 1) THEN IF (color_str EQ 'Gray') THEN color_str = 'LightGrey' ELSE color_str = 'Gray'

              FOR j = 1, 100 DO BEGIN
                IF (cortex_CheckIndex(j,spectrum_str[i]) EQ 1) THEN BEGIN

;help,spectrum_geo,/struct

                  geo = cortex_ExtractStructure(spectrum_geo,j)

                  y_avg = 0.5 * (ymin + ymax)
                  y_del = 0.5 * (ymax - ymin) * 0.98
   
                  v = [ [FLOAT(geo.index[0]), y_avg-y_del], [FLOAT(geo.index[0]), y_avg+y_del],  $
                        [FLOAT(geo.index[1]), y_avg+y_del], [FLOAT(geo.index[1]), y_avg-y_del] ] 

                  IF (ndata EQ 1) THEN BEGIN
                    icolour++
                    color_str = colors[icolour] 
                  ENDIF

print,'icolour====================',icolour,v

                  IF (geo.index[0] EQ geo.index[1]) THEN  $
                    OPLOT   ,v[0,0:1],v[1,0:1], COLOR=TrueColor(color_str), THICK=3.4  ELSE  $
                    POLYFILL,v[0,*  ],v[1,*  ], /DATA, COLOR=TrueColor(color_str), NOCLIP=0, THICK=3.4

                ENDIF

              ENDFOR
            ENDFOR
          ENDIF

          IF ((iplot LE iplot_start+2 OR ndata GT 1) AND idata EQ 1) THEN icolour = -1
          
          IF (idata EQ 1) THEN BEGIN
print,'icolour--------------------',icolour
            cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                            dev_xsize, dev_ysize, charsize_labels, colors, starting_icolour=icolour
            ;   Write sub-title for each plot:
            XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
              (xy_label[2] * ypos[0] + xy_label[3] * ypos[1]) * dev_ysize,  $
              subtitle[iplot-1], CHARSIZE=charsize_labels, /DEVICE
          ENDIF

          IF (idata EQ 1) THEN BEGIN
            FOR itrace = 0, ntrace[iplot-1]-1 DO BEGIN
              icolour++
              i = WHERE(val.y[*,itrace] NE -999.0,count)
              IF (count GT 0) THEN  $
                OPLOT, val.x[i], val.y[i,itrace], COLOR=TrueColor(colors[icolour]), THICK=thick
            ENDFOR
          ENDIF ELSE BEGIN
            itrace = 0
            icolour++
            i = WHERE(val.y[*,itrace] NE -999.0,count)
            IF (count GT 0) THEN  $
              OPLOT, val.x[i], val.y[i,itrace], COLOR=TrueColor(colors[icolour]), THICK=thick
          ENDELSE

          BREAK
          END
;       ----------------------------------------------------------------
        10: BEGIN
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
            3: OPLOT, val.x, val.y[*,1], COLOR=TrueColor(colors[idata-1]), THICK=thick, LINESTYLE=1
            4: 
            ELSE:
          ENDCASE

          ; Show the peak average quantity over a specificed plot range:
          IF ((iplot EQ 1 OR iplot EQ 3) AND plot.show_avg AND cortex_GetValues(plot.xmark,values)) THEN BEGIN
            IF (idata EQ 1) THEN BEGIN
              val_wall = cortex_ExtractStructure(data_array,idata)
              max_average = MAKE_ARRAY(N_ELEMENTS(values),/FLOAT,VALUE=0.0)
            ENDIF
            print, '       >>>>>---------<<<<<'
            FOR i = 0, N_ELEMENTS(values)-1, 2 DO BEGIN
              j = WHERE(xdata GE values[i] AND xdata LE values[i+1])
              print,'trying',values[i:i+1]
              print,'check ',xdata[j]
              print,'check ',val_wall.wall.length[j]
              print,'ydata ',val.y[j,0]
              weight = val_wall.wall.length[j] / TOTAL(val_wall.wall.length[j])
              print,'weight',weight
              print,'weight',TOTAL(weight)
              print,'-->   ',TOTAL(val.y[j,0]*weight)
              max_average[i] = MAX([max_average[i],TOTAL(val.y[j,0]*weight)])
              print,'max_avg=',i,idata,max_average[i]
            ENDFOR
            IF (idata EQ ndata) THEN BEGIN
              FOR i = 0, N_ELEMENTS(values)-1, 2 DO  $
                XYOUTS, MIN([0.75*xmax,values[i]]) + 0.005 * (xmax - xmin),  $
                        0.0                        - 0.075 * (ymax - ymin),  $
                        'MAX_AVG = '+STRTRIM(STRING(MAX(max_average[i]),FORMAT='(E10.2)'),2),  $
                        CHARSIZE=charsize_labels
            ENDIF
          ENDIF

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


