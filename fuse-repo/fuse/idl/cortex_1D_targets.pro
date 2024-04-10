; 
;
; ======================================================================
;
FUNCTION cortex_PlotTargetProfiles, plot, data_array, ps=ps, slice=slice, spike=spike, time_slice=time_slice

  PRINT
  PRINT,'----------------------- NEW PLOT -----------------------'
  PRINT,'PLOT OPTION=',plot.option

  MAXNYDATA = 7
  MAXNXDATA = 1000

; print,'spike',spike,plot.spike

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
  IF (option EQ 100 OR option EQ 101) THEN option = 999

  plot_peak = plot.peak
;
; Setup plot:
; ----------------------------------------------------------------------
;
  CASE option OF
;   --------------------------------------------------------------------
    1: BEGIN
       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot    = 6
       plot_xn  = 2
       plot_yn  = 3
       title = plot.title 
       subtitle = ['LOW INDEX TARGET jsat / A m-2' ,'LOW INDEX TARGET Mach No.' ,'LOW INDEX TARGET Te,i / eV' ,  $
                   'HIGH INDEX TARGET jsat / A m-2','HIGH INDEX TARGET Mach No.','HIGH INDEX TARGET Te,i / eV']
       xtitle   = 'rho (m)'
;       xtitle   = 'psi_n'
       ytitle   = ['jsat (A m-2)','M','Te,i (eV)','jsat (A m-2)','M','Te,i (eV)']
       labels   = MAKE_ARRAY(100,VALUE='',/STRING)
       ntrace   = [1,1,2,1,1,2]
       END
;   --------------------------------------------------------------------
    2: BEGIN
       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot = 4
       plot_xn = 2
       plot_yn = 2
       subtitle = ['INNER TARGET ION SATURATION CURRENT / A m-2','INNER TARGET T_e / eV',  $
                   'OUTER TARGET ION SATURATION CURRENT / A m-2','OUTER TARGET T_e / eV']
       xtitle   = 'psi_n'
       ytitle   = ['j_sat (A m-2)','T_e (eV)','j_sat (A m-2)','T_e (eV)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace = [1]
       END
;   --------------------------------------------------------------------
    3: BEGIN
       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot = 9
       plot_xn = 3
       plot_yn = 3
       title = plot.title 
       subtitle = ['CONNECTION LENGTH / m','RING INDEX','PSIn',  $
                   'ELECTRON DENSITY / m-3','ELECTRON PRESSURE / eV m-3','TARGET T / eV',  $
                   'PARA. CURRENT / A m-2','PARA. HEAT FLUX / MW m-2','PERP. HEAT FLUX / MW m-2' ]
       xtitle   = 'wall index'
       ytitle   = ['L (m)'   ,'ring index','psi_n',  $
                   'n_e (m-3)','p (eV m-3)','T (eV)',  $
                   'jsat_para (A m-2)','q_para (MW m-2)','q_perp (MW m-2)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       labels[4] = 'target:upstream'
       labels[5] = 'Te:Ti'
       ntrace = [1]
       END
;   --------------------------------------------------------------------
    4: BEGIN
       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot    = 6
       plot_xn  = 3
       plot_yn  = 2
       title    = plot.title 
       subtitle = ['OUTER MIDPLANE ne / m-3' ,  $
                   'OUTER TARGET ne / m-3'   ,  $
                   'OUTER MIDPLANE Te / eV'  ,  $
                   'OUTER TARGET Te / eV'    ,  $
                   'OUTER MIDPLANE Ti / eV'  ,  $ 
                   'OUTER TARGET Ti / eV'    ]
       xtitle   = 'psi_n'
       ytitle   = ['ne (m-3)' ,  $
                   'ne (m-3)' ,  $
                   'Te (eV)'  ,  $
                   'Te (eV)'  ,  $
                   'Ti (eV)'  ,  $
                   'Ti (eV)'  ]
       labels   = MAKE_ARRAY(100,VALUE='',/STRING)
       ntrace   = [1,1,2,1,1,2]
       END
;   --------------------------------------------------------------------
    ELSE: BEGIN  
      PRINT, 'ERROR cortex_PlotTargetProfile: Unrecognised plot option'
      PRINT, '  OPTION = ',option,' (',plot.option,')'
      RETURN, -1
      END
  ENDCASE

  colors   = ['Black','Red','Darkgreen','Blue','Orange', 'Purple', 'Silver', 'Hotpink']
;  colors = ['Black','Red','Green','Blue','Orange','Purple', 'Hotpink', 'Darkseagreen', 'Silver']
;
; Setup plot area:
; ----------------------------------------------------------------------
;
  size = TOTAL(plot.frame_bnds[1] - plot.frame_bnds[0])

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
      PRINT, 'ERROR cortex_PlotTargetProfile: No data found'
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
;
;   Collect and organize the data:
;   --------------------------------------------------------------------
;
    xmin =  1.0E+35
    xmax = -1.0E+35
    ymin =  1.0E+35
    ymax = -1.0E+35

    data_store_set = 0



    FOR idata = 1, ndata DO BEGIN

      IF ((slice GT 1) AND (idata NE 1 AND idata NE slice AND idata NE spike)) THEN CONTINUE
 
;      print, 'selecting', idata, ndata

      val = cortex_ExtractStructure(data_array,idata)
      CASE option OF
;       ----------------------------------------------------------------
        1: BEGIN
          file = val.target.file
          integral = ' '
          str = STRSPLIT(file,'/',/EXTRACT)                    ; Extract case name to STR
          str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
          labels[0] = labels[0] + str[0] + integral + ' :'

          xdata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      
          ydata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      

          n = N_ELEMENTS(val.target.psin)

          xdata[0:n-1,0] = val.target.rho[0:n-1]
;          xdata[0:n-1,0] = val.target.psin[0:n-1]

          CASE iplot OF
            1: BEGIN
              ydata[0:n-1,0] = val.target.jsat[0:n-1,0]
              END
            2: BEGIN
              ydata[0:n-1,0] = val.target.M   [0:n-1,0]
              END
            3: BEGIN
              ydata[0:n-1,0] = val.target.te  [0:n-1,0]
              ydata[0:n-1,1] = val.target.ti  [0:n-1,0]
              END
            4: BEGIN
              ydata[0:n-1,0] = val.target.jsat[0:n-1,1]
              END
            5: BEGIN
              ydata[0:n-1,0] = val.target.M   [0:n-1,1]
              END
            6: BEGIN
              ydata[0:n-1,0] = val.target.te  [0:n-1,1]
              ydata[0:n-1,1] = val.target.ti  [0:n-1,1]
              END
          ENDCASE

          ; Write the target data to a file:
          IF (iplot EQ 6 AND plot.id NE 'unknown') THEN BEGIN
          
            FOR i = 0, 1 DO BEGIN

              IF (i EQ 0) THEN BEGIN
                target_tag = 'low' 
                target_loc = 'INNER'
              ENDIF ELSE BEGIN
                target_tag = 'high'
                target_loc = 'OUTER'
              ENDELSE

              fp1 = 3
              FREE_LUN,fp1        
              file_name = 'cortex_data/' + plot.id + '_' + target_tag + '.dat'
              OPENW,fp1,file_name, ERROR=err
              IF (err NE 0) THEN BEGIN
                PRINT,'ERROR cortex_PlotRadialProfile: Problem opening data stream (to file)'
                PRINT,'  FILE_NAME= ',file_name
                PRINTF,-2,!err.msg
                FREE_LUN,fp1
                STOP
              ENDIF
              
;              PRINTF,fp1,'{NUMBER OF COLUMNS}  8',FORMAT='(A)'
              PRINTF,fp1,'*',FORMAT='(A)'
              PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
              PRINTF,fp1,'* RADIAL PLASMA PROFILE FROM CORTEX - ' + STRUPCASE(target_tag) + ' INDEX TARGET (' + target_loc +  ' FOR LSN)'
              PRINTF,fp1,'*',FORMAT='(A)'
              PRINTF,fp1,'* CASE ',plot.case_name[0],FORMAT='(2A)'
              PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
              PRINTF,fp1,'*',FORMAT='(A)'
              PRINTF,fp1,'*      psi','rho','jsat','ne','M','pe','Te','Ti',FORMAT='(8A10)'
              PRINTF,fp1,'*         ','(m)','(A m-2)','(m-3)',' ','(Pa)','(eV)','(eV)',FORMAT='(8A10)'
            
              FOR j = 0, N_ELEMENTS(val.target.psin)-1 DO BEGIN
                PRINTF,fp1,  $
                  val.target.psin[j  ],  $
                  val.target.rho [j  ],  $
                  val.target.jsat[j,i],  $
                  val.target.dens[j,i],  $
                  val.target.M   [j,i],  $
                  val.target.pe  [j,i],  $
                  val.target.te  [j,i],  $
                  val.target.ti  [j,i],  $
                  FORMAT='(2F10.5,2E10.2,F10.2,E10.2,2F10.2)'        
              ENDFOR
            
              CLOSE,fp1
              FREE_LUN,fp1
 
            ENDFOR
          
          ENDIF

          END
;       ----------------------------------------------------------------
        2: BEGIN
          file = val.target.file
          integral = ' '
          str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
          str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
          labels[0] = labels[0] + STRING(idata-1) + '/' + str[0] + integral + ' :'
          ntrace = [1,1,1,1] ; Number of data lines on each plot

          xdata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      
          ydata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      

          CASE iplot OF
            1: BEGIN
              i = WHERE(val.target.location[*,0] EQ 2, count_i)
              xdata[0:count_i-1,0] = val.target.psin[i]
              ydata[0:count_i-1,0] = val.target.jsat[i,0]
              END
            2: BEGIN
              i = WHERE(val.target.location[*,0] EQ 2, count_i)
              xdata[0:count_i-1,0] = val.target.psin[i]
              ydata[0:count_i-1,0] = val.target.te[i,0]
              END
            3: BEGIN
              i = WHERE(val.target.location[*,1] EQ 4, count_i)
              xdata[0:count_i-1,0] = val.target.psin[i]
              ydata[0:count_i-1,0] = val.target.jsat[i,1]
              END
            4: BEGIN
              i = WHERE(val.target.location[*,1] EQ 4, count_i)
              xdata[0:count_i-1,0] = val.target.psin[i]
              ydata[0:count_i-1,0] = val.target.te[i,1]
              END
          ENDCASE
          END
;       ----------------------------------------------------------------
        3: BEGIN  ; Target profile along the wall

          file = val.wall.file   
          integral = ' '
          str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
          str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
;          labels[0] = labels[0] + STRING(idata-1) + '/' + str[0] + integral + ' :'
          ntrace = [1,1,1, 1,2,2 ,1,1,4]  ; Number of data lines on each plot

          xdata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      
          ydata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      

          IF (plot.xorigin EQ 0) THEN BEGIN
            xdata_ref = val.wall.index
          ENDIF ELSE BEGIN
            i = ABS(plot.xorigin) - 1
            xdata_ref = val.wall.dist - val.wall.dist[i]
            i = WHERE(xdata_ref LT 0.0, count)
            IF (count GE 1) THEN  $
              xdata_ref[i] = -1.0 * (MIN(xdata_ref[i]) - xdata_ref[i]) + MAX(xdata_ref) + 0.5 * val.wall.length[i[0]]
            IF (plot.xorigin LT 0) THEN xdata_ref = -1.0 * xdata_ref  
            xtitle = 'distance along wall (m)'
          ENDELSE

          i = WHERE(plot.xrange EQ 0.0, count)
          IF (count EQ 2) THEN i = INDGEN(N_ELEMENTS(val.wall.index)) ELSE  $
                               i = WHERE(xdata_ref GE plot.xrange[0] AND  $
                                         xdata_ref LE plot.xrange[1], count)        

          n = N_ELEMENTS(i)
          xdata[0:n-1,0] = xdata_ref[i      ]
          xdata[0:n-1,1] = xdata    [0:n-1,0]

          CASE iplot OF
            1: BEGIN
              FOR j = 0, n-1 DO  $
                IF (val.wall.jsat[i[j]] NE -999.0) THEN  $
                  ydata[j,0] = val.wall.l[i[j]]
              ;FOR j = 0, n-1 DO BEGIN
              ;  ring = val.wall.index_ring[i[j]]
              ;  k = WHERE(val.midplane.ring EQ ring, count)
              ;  IF (count EQ 1) THEN ydata[j,0] = val.midplane.l[k]
              ;ENDFOR
              END
            2: BEGIN
              FOR j = 0, n-1 DO  $
                IF (val.wall.jsat[i[j]] NE -999.0) THEN  $
                  ydata[j,0] = val.wall.index_ring[i[j]]
              END
            3: BEGIN
              FOR j = 0, n-1 DO  $
                IF (val.wall.jsat[i[j]] NE -999.0) THEN  $
                  ydata[j,0] = val.wall.psin[i[j]]
              END
            4: BEGIN
              FOR j = 0, n-1 DO  $
                IF (val.wall.jsat[i[j]] NE -999.0) THEN  $
                  ydata[j,0] = val.wall.dens[i[j]]
              END
            5: BEGIN
              FOR j = 0, n-1 DO  $
                IF (val.wall.dens[i[j]] NE -999.0) THEN  $
                  ydata[j,0] = val.wall.dens[i[j]] *   $
                                 ( ( val.wall.te[i[j]] + val.wall.te[i[j]] )  + $ 
                                   ( (2.0 * 1.67E-27 / 1.602E-19) * val.wall.vb[i[j]]^2 ) )  ; *** mass hardcoded ***
              FOR j = 0, n-1 DO BEGIN
                ring = val.wall.index_ring[i[j]]
                k = WHERE(val.midplane.ring EQ ring, count)
                IF (count EQ 1) THEN ydata[j,1] = val.midplane.dens[k] * (val.midplane.te[k] + val.midplane.ti[k])
              ENDFOR
              END
            6: BEGIN
              FOR j = 0, n-1 DO BEGIN
                IF (val.wall.te[i[j]] NE -999.0) THEN BEGIN
                  ydata[j,0] = val.wall.te[i[j]]
                  ydata[j,1] = val.wall.ti[i[j]]
                ENDIF
              ENDFOR
              END


            7: BEGIN
              FOR j = 0, n-1 DO  $
                IF (val.wall.jsat[i[j]] NE -999.0) THEN  $
                  ydata[j,0] = val.wall.jsat[i[j]] 
              END
            8: BEGIN
              gamma = 9.0
              FOR j = 0, n-1 DO  $
                IF (val.wall.jsat[i[j]] NE -999.0) THEN  $
                  ydata[j,0] = val.wall.jsat[i[j]] *  $
                               val.wall.te  [i[j]] * gamma / 1.0E+6
              END
            9: BEGIN
              FOR j = 0, n-1 DO  $
                IF (val.wall.jsat[i[j]] NE -999.0) THEN  $
                  ydata[j,0] = val.wall.jsat[i[j]] * val.wall.bratio[i[j]] * val.wall.costet[i[j]] *  $
                               val.wall.te  [i[j]] * gamma / 1.0E+6
              ; CX heat flux               
              ydata[0:n-1,1] = val.wall.atom_par_flux[i] * val.wall.atom_avg_energy[i] * 1.602E-19 / 1.0E+6
              ; CX heat flux               
              ydata[0:n-1,2] = val.wall.mol_par_flux [i] * val.wall.mol_avg_energy [i] * 1.602E-19 / 1.0E+6

              ydata[0:n-1,3] = ydata[0:n-1,1] + ydata[0:n-1,2]

              END

          ENDCASE

          ; Store the data for a big dump
          CASE iplot OF
            1: BEGIN
               dump_data = CREATE_STRUCT(           'xdata' , xdata          [0:n-1,0] ) 
               dump_data = CREATE_STRUCT(dump_data, 'index' , val.wall.index [i]       ) 
               dump_data = CREATE_STRUCT(dump_data, 'costet', val.wall.costet[i]       ) 
               dump_data = CREATE_STRUCT(dump_data, 'bratio', val.wall.bratio[i]       ) 
               dump_data = CREATE_STRUCT(dump_data, 'dens'  , val.wall.dens  [i]       ) 
               dump_data = CREATE_STRUCT(dump_data, 'length', val.wall.length[i]       ) 
               dump_data = CREATE_STRUCT(dump_data, 'r'     , val.wall.r_cen [i]       ) 
               dump_data = CREATE_STRUCT(dump_data, 'z'     , val.wall.z_cen [i]       ) 
               dump_data = CREATE_STRUCT(dump_data, 'psin'  , val.wall.psin  [i]       ) 
               END
            3: dump_data = CREATE_STRUCT(dump_data, 'rho'   , ydata          [0:n-1,0] )  
            6: BEGIN
               dump_data = CREATE_STRUCT(dump_data, 'te'    , ydata          [0:n-1,0] )  
               dump_data = CREATE_STRUCT(dump_data, 'ti'    , ydata          [0:n-1,1] )  
               END
            7: dump_data = CREATE_STRUCT(dump_data, 'jsat'  , ydata          [0:n-1,0] )  
            8: dump_data = CREATE_STRUCT(dump_data, 'q_para', ydata          [0:n-1,0] )  
            9: BEGIN
               dump_data = CREATE_STRUCT(dump_data, 'q_perp', ydata          [0:n-1,0] )  

               fp1 = 3
               FREE_LUN,fp1        
               file_name = plot.case_name[0] + '.heat_flux_' + plot.id 
               OPENW,fp1,'heat_flux/'+file_name, ERROR=err
               IF (err NE 0) THEN BEGIN
                 PRINT,'ERROR cortex_PlotTargetProfiles: Problem opening data stream (to file)'
                 PRINT,'  FILE_NAME= ',file_name
                 PRINTF,-2,!err.msg
                 FREE_LUN,fp1
                 STOP
               ENDIF
      	 
               PRINTF,fp1,'*',FORMAT='(A)'
               PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
               PRINTF,fp1,'* WALL HEAT FLUX DATA FROM OSM'
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
               PRINTF,fp1,'{GAMMA}',FORMAT='(A)'
               PRINTF,fp1,gamma
               PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
               PRINTF,fp1,'{WALL SEGMENT DATA}',FORMAT='(A)'
               PRINTF,fp1,N_ELEMENTS(dump_data.xdata),FORMAT='(I12)'
               PRINTF,fp1,'*',FORMAT='(A)'
               PRINTF,fp1,'* dist    - distance from strike-point (usually -- check rho value for dist=0.0 to confirm)',FORMAT='(A)'
               PRINTF,fp1,'* length  - length of wall segment',FORMAT='(A)'
               PRINTF,fp1,'* jsat    - parallel ion saturation current to the wall/target segment',FORMAT='(A)'
               PRINTF,fp1,'* n_e     - electron density at the entrance to the sheath',FORMAT='(A)'
               PRINTF,fp1,'* T_e     - electron temperature at the entrance to the sheath',FORMAT='(A)'
               PRINTF,fp1,'* T_i     - background hydrogenic ion temperature at the entrance to the sheath',FORMAT='(A)'
               PRINTF,fp1,'* q_para  - parallel heat flux density at the wall/target segment',FORMAT='(A)'
               PRINTF,fp1,'* q_perp  - perpendicular hear flux density to the wall/target segment',FORMAT='(A)'
               PRINTF,fp1,'* cosine  - consine of the angle between the incident field line and the surface in the poloidal plane',FORMAT='(A)'
               PRINTF,fp1,'* B_ratio - ratio of the poloidal field to the toroidal field',FORMAT='(A)'
               PRINTF,fp1,'* rho     - distance of the wall/target segment from the separatrix, mapped to the at the outer midplane (approx.)',FORMAT='(A)'
               PRINTF,fp1,'* psi_n   - normalized magnetic flux coordinate',FORMAT='(A)'
               PRINTF,fp1,'*',FORMAT='(A)'
               PRINTF,fp1,'* The text string in the right-most column gives the name of the file that contains the velocity '
               PRINTF,fp1,'* distribution for each cell, as calculated along the line-of-sight in EIRENE.'
               PRINTF,fp1,'*',FORMAT='(A)'
               PRINTF,fp1,'* ','index','dist','length','jsat' ,'n_e','T_e','T_i','q_para','q_perp','cosine','B_ratio','rho','psi_n',FORMAT='(A2,A6,2X,4A10,2A6,6A10)'
               PRINTF,fp1,'* ',' '    ,'m'   ,'m'     ,'A m-2','m-3','eV' ,'eV' ,'MW m-2','MW m-2',' '     ,' '      ,'m'  ,' '    ,FORMAT='(A2,A6,2X,4A10,2A6,6A10)'

               j = SORT(dump_data.xdata)

               FOR i = 0, N_ELEMENTS(dump_data.xdata)-1 DO BEGIN

                 IF (dump_data.jsat[j[i]] EQ -999.0) THEN CONTINUE

                 PRINTF,fp1,  $
                   dump_data.index [j[i]],  $
                   dump_data.xdata [j[i]],  $
		   dump_data.length[j[i]],  $
		   dump_data.jsat  [j[i]],  $
		   dump_data.dens  [j[i]],  $
		   dump_data.te    [j[i]],  $
		   dump_data.ti    [j[i]],  $
		   dump_data.q_para[j[i]],  $
		   dump_data.q_perp[j[i]],  $
		   dump_data.costet[j[i]],  $
		   dump_data.bratio[j[i]],  $
		   dump_data.rho   [j[i]],  $
		   dump_data.psin  [j[i]],  $
                   FORMAT= '(2X,I6,2X,2F10.6,2E10.2,2F6.1,2F10.4,2E10.2,2F10.5)'
               ENDFOR

               CLOSE,fp1
               FREE_LUN,fp1

               END
            ELSE:
          ENDCASE

          END
;       ----------------------------------------------------------------
        4: BEGIN
          file = val.target.file
          integral = ' '

          IF (time_slice NE -1) THEN BEGIN          
            labels = MAKE_ARRAY(100,VALUE='',/STRING)
            labels[0] = '               ' + STRTRIM(STRING(time_slice),2) + ' us'
          ENDIF ELSE BEGIN
            str = STRSPLIT(file,'/',/EXTRACT)                    ; Extract case name to STR
            str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
            labels[0] = labels[0] + str[0] + integral + ' :'
          ENDELSE

          xdata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      
          ydata = MAKE_ARRAY(MAXNXDATA,MAXNYDATA,/FLOAT,VALUE=-999.0)      

          n = N_ELEMENTS(val.target.psin  )
          o = N_ELEMENTS(val.midplane.psin)

          CASE iplot OF
            1: BEGIN
              xdata[0:o-1,0] = val.midplane.psin[0:o-1]
              ydata[0:o-1,0] = val.midplane.dens[0:o-1]
              END
            2: BEGIN
              xdata[0:n-1,0] = val.target.psin[0:n-1  ]
              ydata[0:n-1,0] = val.target.dens[0:n-1,1]
              END
            3: BEGIN
              xdata[0:o-1,0] = val.midplane.psin[0:o-1]
              ydata[0:o-1,0] = val.midplane.te  [0:o-1]
              END
            4: BEGIN
              xdata[0:n-1,0] = val.target.psin[0:n-1  ]
              ydata[0:n-1,0] = val.target.te  [0:n-1,1]
              END
            5: BEGIN
              xdata[0:o-1,0] = val.midplane.psin[0:o-1]
              ydata[0:o-1,0] = val.midplane.ti  [0:o-1]
              END
            6: BEGIN
              xdata[0:n-1,0] = val.target.psin[0:n-1  ]
              ydata[0:n-1,0] = val.target.ti  [0:n-1,1]
; print,'ti',val.target.ti  [0:n-1,1]
              END
          ENDCASE

          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotTargetProfile: Unrecognised plot option'
          PRINT, '  OPTION = ',plot.option
          RETURN, -1
          END
      ENDCASE

;     Package data for plotting:
      name = 'data' + STRING(idata,FORMAT='(I0)')
      data = { n : ndata, x : xdata, y : ydata, file : file } 
      IF (NOT data_store_set) THEN data_store = CREATE_STRUCT(           name,data) ELSE  $
                                   data_store = CREATE_STRUCT(data_store,name,data)
      data_store_set = 1
      IF (N_ELEMENTS(WHERE(plot.xrange EQ 0.0)) NE 2) THEN BEGIN
        i = WHERE(xdata GE plot.xrange[0] AND xdata LE plot.xrange[1])        
        IF (N_ELEMENTS(i) EQ 1) THEN BEGIN
          PRINT,'ERROR cortex_PlotTargetProfile: No data within XRANGE'
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
      k = WHERE(ydata1 NE -999.0, count)
      IF (count EQ 0) THEN BEGIN
       PRINT,'ERROR cortex_PlotTargetProfile: No valid Y data found'
        PRINT,'  IPLOT          = ',iplot
        PRINT,'  IDATA          = ',idata
        RETURN, -1
      ENDIF
      ymin = MIN([ymin,ydata1[k]])
      ymax = MAX([ymax,ydata1[k]])
    ENDFOR
;    IF (ymin GT 0.0 AND ymax GT 0.0) THEN ymin = 0.0  ; Makes things a bit clearer on the plots I think...
;    IF (ymin LT 0.0 AND ymax LT 0.0) THEN ymax = 0.0
    deltay = ymax - ymin
    ymin = ymin - 0.05 * deltay
    ymax = ymax + 0.05 * deltay

;
;   Plot away:
;   --------------------------------------------------------------------
;
;   Axes:
    IF (plot.flip EQ 1) THEN xrange = [xmax,xmin] ELSE xrange = [xmin,xmax]
    yrange = [ymin,ymax]
    position = [xpos[0],ypos[0],xpos[1],ypos[1]]

    plot_type = default_plot_type                                           
    IF (focus NE 0 OR yi EQ plot_yn OR iplot EQ nplot) THEN plot_type = 2   ; Show x-axis label

    CASE (plot_type) OF
      1: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTICKFORMAT='(A1)',/NOERASE,  $
               XRANGE=xrange,YRANGE=yrange
      2: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[iplot-1],XTITLE=xtitle,/NOERASE,  $
               XRANGE=xrange,YRANGE=yrange
    ENDCASE

;   Write sub-title for each plot:
    XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
            (xy_label[2] * ypos[0] + xy_label[3] * ypos[1]) * dev_ysize,  $
            subtitle[iplot-1], CHARSIZE=charsize_labels, /DEVICE

;   Data:
    FOR idata = 1, ndata DO BEGIN

      IF ((slice GT 1) AND (idata NE 1 AND idata NE slice AND idata NE spike)) THEN CONTINUE

      val = cortex_ExtractStructure(data_store,idata)
      CASE option OF
;       ----------------------------------------------------------------
        1: BEGIN
          cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                          dev_xsize, dev_ysize, charsize_labels, colors

          OPLOT, [xmin,xmax], [0.0,0.0    ], LINESTYLE=1, COLOR=TrueColor('Black') 
          OPLOT, [0.0 ,0.0 ], [-1.0E+20,1.0E+20], LINESTYLE=1, COLOR=TrueColor('Black') 
;          OPLOT, [1.0 ,1.0 ], [0.0,1.0E+20], LINESTYLE=1, COLOR=TrueColor('Black') 
          i = WHERE(val.y[*,0] NE -999.0,count_i)

          IF (count_i GT 0) THEN BEGIN
            j = SORT(val.x[i,0])
            OPLOT, val.x[i[j],0], val.y[i[j],0], COLOR=TrueColor(colors[idata-1])
            CASE iplot OF
              1: 
              2:
              3: OPLOT, val.x[i[j],0], val.y[i[j],1], COLOR=TrueColor(colors[idata-1]), LINESTYLE=1
              4: 
              5: 
              6: OPLOT, val.x[i[j],0], val.y[i[j],1], COLOR=TrueColor(colors[idata-1]), LINESTYLE=1
              ELSE:
            ENDCASE
          ENDIF

;          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
;          OPLOT, val.x[*,0], val.y[*,0], COLOR=TrueColor(colors[idata-1])
;          CASE iplot OF
;            1: 
;            2: 
;            3: 
;            ELSE:
;          ENDCASE

          END
;       ----------------------------------------------------------------
        2: BEGIN
          cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                          dev_xsize, dev_ysize, charsize_labels, colors

          OPLOT, [xmin,xmax], [0.0,0.0    ], LINESTYLE=1, COLOR=TrueColor('Black') 
          OPLOT, [1.0 ,1.0 ], [0.0,1.0E+20], LINESTYLE=1, COLOR=TrueColor('Black') 
          i = WHERE(val.y[*,0] NE -999.0,count_i)
          IF (count_i GT 0) THEN BEGIN
            x = val.x[i,0]
            y = val.y[i,0]
            j = SORT(x)
            OPLOT, x[j], y[j], COLOR=TrueColor(colors[idata-1])
            CASE iplot OF
              1: 
              2: 
              3: 
              ELSE:
            ENDCASE
          ENDIF
          END
;       ----------------------------------------------------------------
        3: BEGIN
          cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                          dev_xsize, dev_ysize, charsize_labels, colors

          OPLOT, [xmin,xmax], [0.0,0.0    ], LINESTYLE=1, COLOR=TrueColor('Black') 
          OPLOT, [1.0 ,1.0 ], [0.0,1.0E+20], LINESTYLE=1, COLOR=TrueColor('Black') 
     
          FOR itrace = 0, ntrace[iplot-1]-1 DO BEGIN

            i = 0
            n = N_ELEMENTS(val.y[*,itrace])
            l = SORT(val.x[*,itrace])
            x = val.x[l,itrace]
            y = val.y[l,itrace]
            WHILE (i LT n) DO BEGIN
              FOR j = i  , n-1 DO IF (y[j] NE -999.0) THEN BREAK
              FOR k = j+1, n-1 DO IF (y[k] EQ -999.0) THEN BREAK            
              k = MIN([k,n])
              IF (j EQ n) THEN BREAK
              IF (k GT j+1) THEN OPLOT,x[j:k-1],y[j:k-1],COLOR=TrueColor(colors[itrace])  $
                            ELSE OPLOT,x[j:k-1],y[j:k-1],COLOR=TrueColor(colors[itrace]),PSYM=3  ; single point
              i = k + 1
            ENDWHILE              
;            i = WHERE(val.y[*,itrace] NE -999.0,count_i)
;            IF (count_i GT 0) THEN BEGIN
;              x = val.x[i,itrace]
;              y = val.y[i,itrace]
;              j = SORT(x)
;              OPLOT, x[j], y[j], COLOR=TrueColor(colors[itrace+1])
;            ENDIF
;            return, 0
          ENDFOR
          END
;       ----------------------------------------------------------------
        4: BEGIN
          cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                          dev_xsize, dev_ysize, charsize_labels, colors

          OPLOT, [xmin,xmax], [0.0,0.0    ], LINESTYLE=1, COLOR=TrueColor('Black') 
          OPLOT, [1.0 ,1.0 ], [0.0,1.0E+20], LINESTYLE=1, COLOR=TrueColor('Black') 
          i = WHERE(val.y[*,0] NE -999.0,count_i)

          time = (plot.time_span[0] + plot.time_span[2]*(idata-1))*10-1000 
;          print,'out of time', time

;          IF (idata LT ndata) THEN BEGIN
;            IF (time LE 400) THEN colour = TrueColor('Red')  $
;                             ELSE colour = TrueColor('Blue')
;          ENDIF ELSE BEGIN

           colour = TrueColor('Black')
           IF (KEYWORD_SET(slice)) THEN BEGIN
             IF (idata EQ 1     AND idata NE slice) THEN colour = TrueColor('Blue')
             IF (idata EQ spike AND idata NE slice) THEN colour = TrueColor('Red')
           ENDIF

;          colour = TrueColor(colors[idata-1])

          IF (count_i GT 0) THEN BEGIN
            j = SORT(val.x[i,0])
            OPLOT, val.x[i[j],0], val.y[i[j],0], COLOR=colour
          ENDIF

;          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 
;          OPLOT, val.x[*,0], val.y[*,0], COLOR=TrueColor(colors[idata-1])
;          CASE iplot OF
;            1: 
;            2: 
;            3: 
;            ELSE:
;          ENDCASE

          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotTargetProfile: Unrecognised plot option'
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
FUNCTION cortex_PlotTargetProfiles_OLD, target, ps=ps

;.resolve_all
;@script
;help,/source_files

  tags = TAG_NAMES(target)
  print,tags

  help,target,/struct
;  help,target.data1,/struct
;  help,target.data2,/struct

  ndata = N_ELEMENTS(TAG_NAMES(target))
  IF (ndata LE 0) THEN BEGIN
    PRINT, 'ERROR Cortex PlotTargetProfiles: No data found'
    RETURN, -1
  ENDIF
  PRINT,'NDATA= ',ndata


  window_id = 0
  window_xsize = 700
  window_ysize = 700

  !P.BACKGROUND = TrueColor('White')

  plot_xn = 2
  plot_yn = 3
  plot_xboarder = 0.1
  plot_yboarder = 0.1

  plot_title  = 'TARGET PLASMA PROFILES'
  plot_xtitle = 'psin (m)'
  plot_ytitle = '??? (???)'

  colors = ['Black','Red','Green','Blue']

;
; Setup plot:
; ----------------------------------------------------------------------
;
  plot_xi = [1 , 1, 1, 2, 2, 2]
  plot_yi = [1 , 2, 3, 1, 2, 3]

  xsize = (1.0 - 2.0 * plot_xboarder) / FLOAT(plot_xn)
  ysize = (1.0 - 2.0 * plot_yboarder) / FLOAT(plot_yn)

  FOR iplot = 1, 6 DO BEGIN

    xi = plot_xi[iplot-1]
    yi = plot_yi[iplot-1]

    PRINT,'XI,YI=',xi,yi

    xcen  =       xsize * (0.5 + FLOAT(xi - 1)) + plot_xboarder
    ycen  = 1.0 - ysize * (0.5 + FLOAT(yi - 1)) - plot_yboarder
    xpos  = [xcen - 0.5  * xsize, xcen + 0.35 * xsize]
    ypos  = [ycen - 0.45 * ysize, ycen + 0.5  * ysize]

    xmin =  1.0E+35
    xmax = -1.0E+35
    ymin =  1.0E+35
    ymax = -1.0E+35
    FOR idata = 1, ndata DO BEGIN
      val = cortex_ExtractStructure(target,idata)
      xdata = val.psin
      IF (iplot EQ 1) THEN ydata = val.jsat[*,0]
      IF (iplot EQ 2) THEN ydata = val.te  [*,0]
      IF (iplot EQ 3) THEN ydata = val.ti  [*,0]
      IF (iplot EQ 4) THEN ydata = val.jsat[*,1]
      IF (iplot EQ 5) THEN ydata = val.te  [*,1]
      IF (iplot EQ 6) THEN ydata = val.ti  [*,1]

      IF (iplot EQ 2) THEN ydata = [ydata,val.ti[*,0]]
      IF (iplot EQ 5) THEN ydata = [ydata,val.ti[*,1]]

      xmin = MIN([xmin,xdata])
      xmax = MAX([xmax,xdata])
      ymin = MIN([ymin,ydata])
      ymax = MAX([ymax,ydata])
    ENDFOR
;    xmin =  4.1
;    xmax =  4.2
;    ymin = -4.0
;    ymax = 0.5E+19
;    PRINT, 'XMIN,MAX=',xmin,xmax
;    PRINT, 'YMIN,MAX=',ymin,ymax

;   Axes:
    xrange = [xmin,xmax]
    yrange = [ymin,ymax]
    position = [xpos[0],ypos[0],xpos[1],ypos[1]]

    CASE (iplot) OF
      1: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=plot_ytitle,XTICKFORMAT='(A1)'                                    
      2: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=plot_ytitle,XTICKFORMAT='(A1)',/NOERASE
      3: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=plot_ytitle,XTITLE=plot_xtitle,/NOERASE
      4: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=plot_ytitle,XTICKFORMAT='(A1)',/NOERASE
      5: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=plot_ytitle,XTICKFORMAT='(A1)',/NOERASE
      6: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=plot_ytitle,XTITLE=plot_xtitle,/NOERASE                                 
    ENDCASE

    XYOUTS, 0.95*xmin+0.05*xmax, 0.1*ymin+0.9*ymax, plot_title
;   Data:
    FOR idata = 1, ndata DO BEGIN
      val = cortex_ExtractStructure(target,idata)
      xdata = val.psin
      IF (iplot EQ 1) THEN ydata = val.jsat[*,0]
      IF (iplot EQ 2) THEN ydata = val.te  [*,0]
      IF (iplot EQ 3) THEN ydata = val.ti  [*,0]
      IF (iplot EQ 4) THEN ydata = val.jsat[*,1]
      IF (iplot EQ 5) THEN ydata = val.te  [*,1]
      IF (iplot EQ 6) THEN ydata = val.ti  [*,1]

      OPLOT, xdata, ydata, COLOR=TrueColor(colors[idata-1])
 
      IF (iplot EQ 2) THEN BEGIN
        ydata = val.ti[*,0]
        OPLOT, xdata, ydata, COLOR=TrueColor(colors[idata-1]),LINESTYLE=2 
      ENDIF
      IF (iplot EQ 5) THEN BEGIN
        ydata = val.ti[*,1]
        OPLOT, xdata, ydata, COLOR=TrueColor(colors[idata-1]),LINESTYLE=2 
      ENDIF

      step = 0.05
      frac = step * FLOAT(idata) + 0.1
      PRINT,frac
      IF (iplot EQ 1) THEN  $
        XYOUTS, 0.95*xmin+0.05*xmax,frac*ymin+(1.0-frac)*ymax,  $
                val.file, COLOR=TrueColor(colors[idata-1])
    ENDFOR

  ENDFOR

  RETURN, 0

END
;
; ======================================================================
;


