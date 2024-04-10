
; ======================================================================
;
FUNCTION cortex_PlotEireneSpectra, plot, data_array, ps=ps

  PRINT
  PRINT,'----------------------- NEW PLOT -----------------------'
  PRINT

  ndata = N_ELEMENTS(TAG_NAMES(data_array))
  IF (ndata LE 0) THEN BEGIN
    PRINT, 'ERROR cortex_PlotEnergySpectrum: No data found'
    RETURN, -1
  ENDIF

;help,data_array,/struct
;help,data_array.data1,/struct
;help,data_array.data1.SPECTRUM_41_0000281,/struct
;help,data_array.data1.SPECTRUM_02_0000242,/struct
;stop

  MAXNYDATA = 7

  nplot_max = 999

  focus = plot.focus
;  focus = MAX([1,plot.focus])

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
       nplot = 1
       plot_xn = 1
       plot_yn = 1
       title = plot.title 
       subtitle = ['CASE NAME']
       xtitle   = 'BIN Te (eV)'
       ytitle   = ['COUNTS (Amps / BIN_eV)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace = [1]
       END
;   --------------------------------------------------------------------
    2: BEGIN
;      Count the SPECTRUM entries in the data structure:
       FOR idata = 1, ndata DO BEGIN
         val = cortex_ExtractStructure(data_array,idata)
         spectra = val.(1)  ; should always be the case
         tags = TAG_NAMES(val)
         count = 0L

print, 'tags', N_ELEMENTS(tags)

;print,  tags
         FOR itag = 0L, N_ELEMENTS(tags)-1 DO BEGIN
           IF (STRPOS(tags[itag],'SPECTRUM') EQ -1) THEN CONTINUE
           spectrum = val.(itag)
           index = WHERE(tags EQ STRUPCASE(spectrum.profile_name), check)
           IF (check LE 0) THEN BEGIN
             PRINT, 'ERROR cortex_PlotEireneSpectra: Profile not found'
             PRINT, '  profile name = ',spectrum.profile_name
             STOP
           ENDIF
           profile = val.(index[0]) 
;           print,'index',spectrum.spectra_index,spectrum.profile_index
           count++
         ENDFOR



         ; Check that all cases in the list have the same number of spectra:
         IF (idata EQ 1) THEN  $
           count_master = count  $
         ELSE BEGIN
           IF (count NE count_master) THEN BEGIN
             PRINT, 'ERROR cortex_PlotEireneSpectra: Case spectrum counts do not match'
             STOP
           ENDIF
         ENDELSE
       ENDFOR


       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot = count
       nplot_max = MIN([12,count])
       plot_xn = [1, 1, 1, 2, 3, 3, 4, 4, 4, 4, 4, 4] 
       plot_yn = [1, 2, 3, 2, 2, 2, 2, 2, 3, 3, 3, 3]
       plot_xn = (plot_xn[nplot_max-1])[0]
       plot_yn = (plot_yn[nplot_max-1])[0]
       title = plot.title 
       subtitle = ['CASE NAME']
       xtitle   = 'BIN v (m/s)'
       ytitle   = ['COUNTS (Amps / BIN_m/s)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace = [1]

       END

;   --------------------------------------------------------------------
    3: BEGIN
;      Plot the spectra along a line-of-sight on the same plot:

       val = cortex_ExtractStructure(data_array,1)


       tags = TAG_NAMES(val)
        

       spectra = val.(1)  ; should always be the case - but, weak
       profile = val.(2)  ; same


; help,val,/struct
help,spectra,/struct
help,profile,/struct

       spectra_i = WHERE(spectra.index EQ val.index, count)

;print, 'spectra_i', spectra_i 

print,count

spectra_n = N_ELEMENTS(profile.path)

spectrum_n = 30 ; 1 ; 30

overlap = 0.5

       file = spectra.file
       str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
       str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)

       default_plot_type = 1
       plot_xboarder = 0.05
       plot_yboarder = 0.1
       plot_xspacing = 0.100
       plot_yspacing = 0.025
       nplot     = 3
       nplot_max = 3
       plot_xn = [1]
       plot_yn = [3]
       title = 'EIRENE VDF along LOS: case ' + str[0] + ', ' + STRTRIM(plot.id,2) + ' view' 
       subtitle = [' ']
       xtitle   = 'BIN (overlapped)'
       ytitle   = ['COUNTS (Amps / BIN)']
;       ytitle   = ['COUNTS (Amps / BIN)']
       labels   = MAKE_ARRAY(100,VALUE=' ',/STRING)
       ntrace = [spectrum_n]

       END
;   --------------------------------------------------------------------
    ELSE: BEGIN  
      PRINT, 'ERROR cortex_PlotEireneSpectra: Unrecognised plot option'
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
      val = cortex_ExtractStructure(data_array,idata)

      IF (option EQ 2) THEN BEGIN

        ; Unfortunately, have to assemble the list every time through this list:
        tags = TAG_NAMES(val)
        list = [-1]
        FOR itag = 0, N_ELEMENTS(tags)-1 DO BEGIN
          IF (STRPOS(tags[itag],'SPECTRUM') EQ -1) THEN CONTINUE
          IF (list[0] EQ -1) THEN list = [itag] ELSE list = [list,itag]
        ENDFOR

        spectrum = val.(list[iplot-1])
        i = WHERE(tags EQ STRUPCASE(spectrum.profile_name))
        profile = val.(i[0]) 

        ; Prepare output files:
        IF (idata EQ 1 AND iplot EQ 1) THEN BEGIN

          fp1 = 3
          FREE_LUN,fp1        
          file_name = plot.case_name[0] + '.profile_' + plot.id + '_' + STRING(FIX(plot.chord_index),FORMAT='(I3.3)')  ; STRING(spectrum.chord_index,FORMAT='(I3.3)')
          OPENW,fp1,'spectra/'+file_name, ERROR=err
          IF (err NE 0) THEN BEGIN
            PRINT,'ERROR cortex_PlotEireneSpectra: Problem opening data stream (to file)'
            PRINT,'FILE_NAME= ',file_name
            PRINTF,-2,!err.msg
            FREE_LUN,fp1
            STOP
          ENDIF
     
          PRINTF,fp1,'*',FORMAT='(A)'
          PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
          PRINTF,fp1,'* LINE-OF-SITE PROFILE DATA FROM RAY (THIS DATA FILE GENERATED IN CORTEX)'
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
          PRINTF,fp1,'{VIEW IDENTIFIER AND CHORD INDEX}',FORMAT='(A)'
          PRINTF,fp1,plot.id,spectrum.chord_index,FORMAT='(3X,A,I12)'
          PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
          PRINTF,fp1,'{R,Z START AND END POINTS OF THE CHORD (m)}',FORMAT='(A)'
          PRINTF,fp1,profile.v1,FORMAT='(3F10.4)'
          PRINTF,fp1,profile.v2,FORMAT='(3F10.4)'
          PRINTF,fp1,'*',FORMAT='(A)'
          PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
          PRINTF,fp1,'{EMISSION SIGNALS}',FORMAT='(A)'
          PRINTF,fp1,N_ELEMENTS(profile.z),FORMAT='(I12)'
          PRINTF,fp1,'*',FORMAT='(A)'
          PRINTF,fp1,'* Z          - atomic number of particle',FORMAT='(A)'
          PRINTF,fp1,'* A          - atomic weight',FORMAT='(A)'
          PRINTF,fp1,'* CHARGE     - electric charge',FORMAT='(A)'
          PRINTF,fp1,'* WAVELENGTH - of the line',FORMAT='(A)'
          PRINTF,fp1,'* INTEGRAL   - signal strength integrated along the chord, in units of [photons ster-1 m-2 s-1]',FORMAT='(A)'
          PRINTF,fp1,'*',FORMAT='(A)'
          PRINTF,fp1,'* ','SIGNAL','Z'   ,'A'   ,'CHARGE','WAVELENGTH','INTEGRAL',FORMAT='(A2,A6,3A8,2A12)'
          PRINTF,fp1,'* ','     ','    ','     ','      ','(nm)'      ,'        ',FORMAT='(A2,A6,3A8,2A12)'
          FOR i = 0, N_ELEMENTS(profile.z)-1 DO  $
            PRINTF,fp1,i+1                ,  $
                       profile.z       [i],  $
                       profile.a       [i],  $
                       profile.charge  [i],  $
                       profile.wlngth  [i],  $
                       profile.integral[i],  $
                       FORMAT='(2X,I6,3I8,F12.2,E12.2)'
          PRINTF,fp1,'*',FORMAT='(A)'
          PRINTF,fp1,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
          PRINTF,fp1,'{INTEGRATION VOLUMES}',FORMAT='(A)'
          PRINTF,fp1,N_ELEMENTS(list),FORMAT='(I12)'
          PRINTF,fp1,'*',FORMAT='(A)'
          PRINTF,fp1,'* DIST   - distance from the starting point along the chord (line-of-sight)',FORMAT='(A)'
          PRINTF,fp1,'* DELTA  - length of the path through the current integation volume (triangular "cell")',FORMAT='(A)'
          PRINTF,fp1,'* WEIGHT - set to unity for now, but will be used when reflections are included in the analysis',FORMAT='(A)'
          PRINTF,fp1,'* n_e    - electron density',FORMAT='(A)'
          PRINTF,fp1,'* T_e    - electron temperature',FORMAT='(A)'
          PRINTF,fp1,'* T_i    - background hydrogenic ion temperature',FORMAT='(A)'
          PRINTF,fp1,'* n_D    - deuterium atom density',FORMAT='(A)'
          PRINTF,fp1,'* n_D2   - deuterium molecule density',FORMAT='(A)'
          PRINTF,fp1,'* n_i+X  - density of impurity with charge X',FORMAT='(A)'
          PRINTF,fp1,'* SIGNAL - local emission in units of [photons ster-1 m-3 s-1]',FORMAT='(A)'
          PRINTF,fp1,'*',FORMAT='(A)'
          PRINTF,fp1,'* The text string in the right-most column gives the name of the file that contains the velocity '
          PRINTF,fp1,'* distribution for each cell, as calculated along the line-of-sight in EIRENE.'
          PRINTF,fp1,'*',FORMAT='(A)'
          PRINTF,fp1,'* ','INDEX','DIST','DELTA','WEIGHT','n_e'  ,'T_e' ,'T_i','n_D','n_D2','n_i+0','n_i+1','n_i+2','n_i+3','SIGNAL_1',FORMAT='(A2,A6,4A10,2A8,7A10)'
          PRINTF,fp1,'* ',' '    ,'[m]' ,'[m]'  ,' '     ,'[m-3]','[eV]','[eV]','[m-3]','[m-3]','[m-3]','[m-3]','[m-3]','[m-3]',FORMAT='(A2,A6,4A10,2A8,6A10)'

;          fact = 4.0 * PI
;          DO i2 = 1, opt%int_num
;            WRITE(fp,'(2X,I6,3I8,F12.2,1P,E12.2,0P)') i2,
;     .        opt%int_z         (i2)     ,
;     .        opt%int_a         (i2)     ,
;     .        opt%int_charge    (i2)     ,
;     .        opt%int_wlngth    (i2)     ,
;     .        pixel(i1)%integral(i2)/fact
;          ENDDO

;        CLOSE,fp1
;        FREE_LUN,fp1
        ENDIF
      ENDIF

      CASE option OF
;       ----------------------------------------------------------------
        1: BEGIN
          file = val.file
          integral = '   PARTICLE FLUX INTEGRAL= ' + STRING(val.integral,FORMAT='(E12.4)')
          str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
          str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
          labels[0] = labels[0] + STRING(idata-1) + '/' + str[0] + integral + ' :'

          xdata = val.bin
          ydata = val.flux 
          adata = [0]
          IF (MAX(ydata) GT 0.0) THEN ydata = ydata / TOTAL(ydata)

          END
;       ----------------------------------------------------------------
        2: BEGIN

          fp2 = 4
          FREE_LUN,fp2

          file_name = plot.case_name[0] + '.' + STRLOWCASE(tags[list[iplot-1]]) + '_' + plot.id + '_' + STRING(FIX(plot.chord_index),FORMAT='(I3.3)')

          OPENW,fp2,'spectra/'+file_name,ERROR=err
          IF (err NE 0) THEN BEGIN
            PRINT,'ERROR cortex_PlotEireneSpectra: Problem opening spectrum data stream'
             PRINT,'FILE_NAME= ',file_name
            PRINTF,-2,!error_option.msg
            FREE_LUN,fp2
            STOP
          ENDIF

          PRINTF,fp2,'*',FORMAT='(A)'
          PRINTF,fp2,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
          PRINTF,fp2,'* HYDROGEN ATOM VELOCITY DISTRIBUTION FROM EIRENE (DATA FILE GENERATED IN CORTEX)'
          PRINTF,fp2,'*',FORMAT='(A)'
          PRINTF,fp2,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
          PRINTF,fp2,'* CASE            ',plot.case_name[0],FORMAT='(2A)'
          PRINTF,fp2,'* TITLE           ',FORMAT='(A)'
          PRINTF,fp2,'* DATE AND TIME   ',FORMAT='(A)'
          PRINTF,fp2,'*',FORMAT='(A)'
          PRINTF,fp2,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
          PRINTF,fp2,'{DATA FILE VERSION}',FORMAT='(A)'
          PRINTF,fp2,'     1.0',FORMAT='(A)'
          PRINTF,fp2,'*',FORMAT='(A)'
          PRINTF,fp2,'* -----------------------------------------------------------------------------------------------------',FORMAT='(A)'
          PRINTF,fp2,'{VELOCITY DISTRIBUTION}',FORMAT='(A)'
          PRINTF,fp2,N_ELEMENTS(spectrum.bin),FORMAT='(I12)'
          PRINTF,fp2,'*',FORMAT='(A)'
          PRINTF,fp2,'* BIN    - central value of the velicity bin when sampling in EIRENE'
          PRINTF,fp2,'* FLUX   - particle count in each bin [Amps / bin_delta_v], i.e. flux is normalized to the bin width'
          PRINTF,fp2,'*',FORMAT='(A)'
          PRINTF,fp2,'* ','INDEX','BIN','FLUX',FORMAT='(A,A8,2A14)'
          PRINTF,fp2,'* ',' ','(m/s)'         ,FORMAT='(A,A8, A14)'
          FOR i = 0, N_ELEMENTS(spectrum.bin)-1 DO  $
            PRINTF,fp2,i+1,spectrum.bin [i],  $
                           spectrum.flux[i],  $
                   FORMAT='(2X,I8,2E14.5)'

          CLOSE,fp2
          FREE_LUN,fp2

          i = spectrum.profile_index - 1

          fact = 1.0 / (4.0 * !PI) 
          PRINTF,fp1,i+1,  $
                     profile.path  [i]      ,  $
                     profile.delta [i]      ,  $
                     1.0                    ,  $
                     profile.dens  [i]      ,  $
                     profile.te    [i]      ,  $
                     profile.ti    [i]      ,  $
                     profile.n_d   [i]      ,  $
                     profile.n_d2  [i]      ,  $
                     0.0                    ,  $
                     0.0                    ,  $
                     0.0                    ,  $
                     0.0                    ,  $
                     profile.signal[i]*fact ,  $
                     file_name              ,  $  ; plot.case_name[0]+'.'+STRLOWCASE(tags[list[iplot-1]])    ,  $
                  FORMAT='(2X,I6,F10.5,3E10.2,2F8.1,7E10.2,2X,A)'


          IF (iplot LE nplot_max) THEN BEGIN
            file = spectrum.file
;            integral = '   PARTICLE FLUX INTEGRAL= ' + STRING(val.integral,FORMAT='(E12.4)')
            integral = ' '
            str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
            str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
            labels[0] = labels[0] + STRING(idata-1) + '/' + str[0] + integral + ' :'

            IF (iplot EQ 1 AND idata EQ 1) THEN BEGIN
              xdata = MAKE_ARRAY(N_ELEMENTS(spectrum.bin),MAXNYDATA,/FLOAT,VALUE=0.0)  
              ydata = xdata
              adata = [0]
            ENDIF

            xdata[*,idata-1] = spectrum.bin
            ydata[*,idata-1] = spectrum.flux 
;            IF (MAX(ydata) GT 0.0) THEN ydata = ydata / TOTAL(ydata)
          ENDIF

          END
;       ----------------------------------------------------------------
        3: BEGIN


          step = 1.0 / FLOAT(nplot)
          frac = FLOAT(iplot-1) / FLOAT(nplot)

;print,'stepping',step,frac

          i1 = LONG( frac         * spectrum_n)
          i2 = LONG((frac + step) * spectrum_n)

          step_multiplier = LONG( FLOAT(spectra_n) / FLOAT(spectrum_n) + 0.5 ) ; + 1

;          step_multiplier = LONG( FLOAT(spectra_n) / FLOAT(spectrum_n) - 0.5 ) 

;print,'stepping',iplot,step,frac,i1,i2,spectrum_n,step_multiplier

          FOR i = i1, i2-1 DO BEGIN
;          FOR i = 0, spectrum_n-1 DO BEGIN

            IF (1 EQ 1) THEN BEGIN
              j = i * step_multiplier + 0
;              j = spectra_n-1 - i * step_multiplier

              j = MAX([MIN([spectra_n-1,j]),0])

;              print,spectra_i[j]
;              print,spectra.cell[spectra_i[j]]

              spectrum_tag = STRING('SPECTRUM_',spectra_i[j]+1,'_',spectra.cell[spectra_i[j]],FORMAT='(A,I04,A,I07)')

              spectrum = cortex_ExtractStructure(val,spectrum_tag,/name)

            ENDIF ELSE BEGIN

              j = spectra_n-1 - i

              spectrum_tag = STRING('SPECTRUM_',spectra_i[j]+1,'_',spectra.cell[spectra_i[j]],FORMAT='(A,I04,A,I07)')

              spectrum = cortex_ExtractStructure(val,spectrum_tag,/name)

              spectrum.flux = spectrum.flux * 0.025

              FOR k = 1, 39 DO BEGIN

                l = j - k

                spectrum_tag = STRING('SPECTRUM_',spectra_i[l]+1,'_',spectra.cell[spectra_i[l]],FORMAT='(A,I04,A,I07)')
  
; print,'adding ',spectrum_tag

                spectrum1 = cortex_ExtractStructure(val,spectrum_tag,/name)
  
                spectrum.flux = spectrum.flux + spectrum1.flux * 0.025

              ENDFOR

            ENDELSE


print,'----------------',j,spectra_n,MAX(spectrum.flux),' ',spectrum.file

            file = spectrum.file


            n = N_ELEMENTS(spectrum.flux)

            IF (i EQ i1) THEN BEGIN

              xdata = MAKE_ARRAY((1.0-overlap)*n*(i2-i1)+overlap*n,i2-i1+1,/FLOAT,VALUE=0.0)  
;              xdata = MAKE_ARRAY((1.0-overlap)*n*spectrum_n+overlap*n,spectrum_n,/FLOAT,VALUE=0.0)  
;              xdata = MAKE_ARRAY(n*spectrum_n,spectrum_n,/FLOAT,VALUE=0.0)  
              ydata = xdata
              adata = MAKE_ARRAY(3,i2-i1+1,/FLOAT,VALUE=0.0)


            ENDIF

            IF (0 EQ 1 AND i EQ i1 AND iplot EQ 1) THEN BEGIN

               integral = ' '
               str = STRSPLIT(file,'/',/EXTRACT)                   ; Extract case name to STR
               str = STRSPLIT(str[N_ELEMENTS(str)-1],'.',/EXTRACT)
               labels[0] = labels[0] + STRING(idata-1) + '/' + str[0] + integral + ' :'

            ENDIF

            IF (iplot EQ 1 AND idata EQ 1) THEN BEGIN

            ENDIF

            adata[*,i-i1] = [FLOAT(spectra_i[j]+1),FLOAT(spectra.cell[spectra_i[j]]),spectra.path[spectra_i[j]]]

            xdata[*                                                ,i-i1] = FINDGEN((1.0-overlap)*n*(i2-i1)+overlap*n) ; FINDGEN(n*spectrum_n)

; print,'MAX',MAX(xdata[*,i-i1]),n,i2-i1+1

            ydata[(1.0-overlap)*n*(i-i1):(1.0-overlap)*n*(i-i1)+n-1,i-i1] = spectrum.flux
 
; print,'RANGE',(1.0-overlap)*n*(i-i1),(1.0-overlap)*n*(i-i1)+n-1

;            xdata[*                                      ,i] = FINDGEN((1.0-overlap)*n*spectrum_n+overlap*n) ; FINDGEN(n*spectrum_n)
;            ydata[(1.0-overlap)*n*i:(1.0-overlap)*n*i+n-1,i] = spectrum.flux
;            ydata[(1.0-overlap)*n*i:(1.0-overlap)*n*i+n-1,i] = SMOOTH(spectrum.flux,21)

            print,max(spectrum.flux)

          ENDFOR

;          IF (MAX(ydata) GT 0.0) THEN ydata = ydata / TOTAL(ydata)

          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotEnergySpectrum: Unrecognised plot option'
          PRINT, '  OPTION = ',plot.option
          RETURN, -1
          END
      ENDCASE

      ; Just in case there's a request for a very large number of plots:
      IF (iplot GT nplot_max) THEN CONTINUE

      IF (plot.smooth GT 0) THEN ydata = SMOOTH(ydata,plot.smooth)

;     Package up the data for plotting:
      name = 'data' + STRING(idata,FORMAT='(I0)')
      data = { x : xdata, y : ydata, a : adata, file : file } 
      IF (idata EQ 1) THEN data_store = CREATE_STRUCT(           name,data) ELSE  $
                           data_store = CREATE_STRUCT(data_store,name,data)
      IF (N_ELEMENTS(WHERE(plot.xrange EQ 0.0)) NE 2) THEN BEGIN
        xmax1 = MAX(xdata)
        i = WHERE(xdata/xmax1 GE plot.xrange[0] AND xdata/xmax1 LE plot.xrange[1])        
        IF (N_ELEMENTS(i) EQ 1) THEN BEGIN
          PRINT,'ERROR cortex_PlotEnergySpectrum: No data within XRANGE'
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
;      ydata1 = REFORM(ydata[i,0])
;      FOR j = 1, ntrace[iplot-1]-1 DO ydata1 = [ydata1,REFORM(ydata[i,j])] 
;      ymin = MIN([ymin,ydata1])
;      ymax = MAX([ymax,ydata1])
      ymin = MIN([ymin,ydata[i]])
      ymax = MAX([ymax,ydata[i]])
    ENDFOR

    ; Just in case there's a request for a very large number of plots:
    IF (iplot GT nplot_max) THEN CONTINUE

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

print,'position',position,plot_type


    CASE (plot_type) OF
      1: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[0],XTICKFORMAT='(A1)',/NOERASE                                  
      2: PLOT, xrange, yrange, /NODATA, XSTYLE=1, YSTYLE=1,  $
               POSITION=position,YTITLE=ytitle[0],XTITLE=xtitle,/NOERASE                                 
    ENDCASE

;   Write sub-title for each plot:
    XYOUTS, (xy_label[0] * xpos[0] + xy_label[1] * xpos[1]) * dev_xsize,  $
            (xy_label[2] * ypos[0] + xy_label[3] * ypos[1]) * dev_ysize,  $
            subtitle[0], CHARSIZE=charsize_labels, /DEVICE
;
;   Set the plot trace colours when there are lots of traces:
;
    IF (option EQ 3) THEN BEGIN

      IF (spectrum_n GT 1) THEN BEGIN

          step = 1.0 / FLOAT(nplot)
          frac = FLOAT(iplot-1) / FLOAT(nplot)

print,'stepping',step,frac

          i1 = LONG( frac         * spectrum_n)
          i2 = LONG((frac + step) * spectrum_n)

        frac = FINDGEN(i2-i1+1) / FLOAT(i2-i1)
;        frac = FINDGEN(spectrum_n) / FLOAT(spectrum_n-1)
        frac_color = LONG(frac * 255.0)

 print,'frac_color',i1,i2,frac,frac_color

      ENDIF ELSE BEGIN
        frac_color = [0]
      ENDELSE

;      LOADCT, 33 ; 37 ; waves 3 ; red

    ENDIF




;   Data:
    FOR idata = 1, ndata DO BEGIN
      val = cortex_ExtractStructure(data_store,idata)
      CASE option OF
;       ----------------------------------------------------------------
        1: BEGIN
          cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                          dev_xsize, dev_ysize, charsize_labels, colors

          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 

          val_y = val.y[*,0]
          IF (N_ELEMENTS(val.y[*,0]) GT 100) THEN val_y = SMOOTH(val_y,10)

          OPLOT, val.x, val_y, COLOR=TrueColor(colors[idata-1]) 
          CASE iplot OF
            1: 
            ELSE:
          ENDCASE
          END
;       ----------------------------------------------------------------
        2: BEGIN
          cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                          dev_xsize, dev_ysize, charsize_labels, colors



          OPLOT, [xmin,xmax], [0.0,0.0], LINESTYLE=1, COLOR=TrueColor('Black') 

;          val_y = val.y[*,idata-1]
;          IF (N_ELEMENTS(val.y[*,0]) GT 100) THEN val_y = SMOOTH(val_y,10)

          OPLOT, val.x[*,idata-1], val.y[*,idata-1], COLOR=TrueColor(colors[idata-1]) 
          CASE iplot OF
            1: 
            ELSE:
          ENDCASE
          END
;       ----------------------------------------------------------------
        3: BEGIN
          cortex_DrawKey, iplot, focus, labels, xy_label, xpos, ypos,  $
                          dev_xsize, dev_ysize, charsize_labels, colors



          step = 1.0 / FLOAT(nplot)
          frac = FLOAT(iplot-1) / FLOAT(nplot)

;print,'stepping',step,frac

          i1 = LONG( frac         * spectrum_n)
          i2 = LONG((frac + step) * spectrum_n) - 1

          IF (iplot EQ 1) THEN bin1 = val.a[0,0]

          LOADCT, 33

          FOR i = 0, i2-i1 DO BEGIN

            OPLOT, val.x[*,i], val.y[*,i], COLOR=frac_color[i]

            frac = FLOAT(i + 1) / FLOAT (i2 - i1 + 2) 

print,frac, i, i1 , i2
print,val.a[*,i]

            str = STRING(val.a[0,i],'(',val.a[0,i]-bin1+1,'),',val.a[2,i],FORMAT='(I0,A1,I0,A2,F0.2)')

print,str

            XYOUTS, ((1.0 - frac) * xpos[0] + frac        * xpos[1]) * dev_xsize,  $
                    ( xy_label[2] * ypos[0] + xy_label[3] * ypos[1]) * dev_ysize,  $
                    str, CHARSIZE=charsize*0.8, /DEVICE, ALIGNMENT=0.2

          ENDFOR



          OPLOT, [xmin,xmax], [0.0,0.0], COLOR=TrueColor('Black') ; LINESTYLE=1, COLOR=TrueColor('Black') 

          END
;       ----------------------------------------------------------------
        ELSE: BEGIN  
          PRINT, 'ERROR cortex_PlotEnergySpectrum: Unrecognised plot option'
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


  IF (option NE 3) THEN BEGIN
    CLOSE,fp1
    FREE_LUN,fp1
  ENDIF

  RETURN, 0
END
;
; ======================================================================
;


