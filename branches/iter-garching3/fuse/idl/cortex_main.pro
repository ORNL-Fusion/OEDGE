; example : cortex_main,['t-new-0000a','i-osm.ctx']
; 
PRO cortex_batch

  cortex_main,['i-fwp-4200a','i-fwp.ctx']
  cortex_main,['i-fwp-4201a','i-fwp.ctx']
  cortex_main,['i-fwp-4202a','i-fwp.ctx']
  cortex_main,['i-fwp-4203a','i-fwp.ctx']
  cortex_main,['i-fwp-4204a','i-fwp.ctx']
  cortex_main,['i-fwp-4205a','i-fwp.ctx']
  cortex_main,['i-fwp-4206a','i-fwp.ctx']
  cortex_main,['i-fwp-4207a','i-fwp.ctx'] ; hello
  cortex_main,['i-fwp-4208a','i-fwp.ctx']
  cortex_main,['i-fwp-4209a','i-fwp.ctx']

END
;
; ======================================================================
;
PRO cortex_Exit, nargs
  IF (nargs EQ 4) THEN EXIT, status=-1
  STOP
END
;
; ======================================================================
;
FUNCTION cortex_CheckTubes, tube, plot_tag, plot_option

  i = WHERE(tube NE 0, count)   ; *** Move to a function ***  (same below)
  IF (count EQ 0) THEN BEGIN
    PRINT,'ERROR cortex_CheckTubes: TUBE(S) not specified'
    PRINT,'  PLOT TAG=',STRTRIM(plot_tag,2)
    PRINT,'  OPTION  =',plot_option
    RETURN, 1
  ENDIF

  RETURN, 0
END
;
; ======================================================================
;
PRO cortex_GeneratePlots, args

  COMMON options, dir_structure
  
  option = 3

  ps = 'on'

  nargs = N_ELEMENTS(args)

  case_name  = args[0]
  input_file = args[1]

  IF (nargs EQ 4) THEN BEGIN
    dir_structure = 0
    family = ''
    child  = ''
    input_file = args[2] + '/' + input_file
    data_path  = args[3] + '/'
  ENDIF ELSE BEGIN
    dir_structure = 1

    str = STRSPLIT(case_name,'-',/EXTRACT)
    family = str[0] + '-' + str[1] + '/'
;    family = STRMID(case_name,0,5) + '/'
    child  = STRMID(case_name,0,7) + '/'

;    input_file = '/home/slisgo/fuse/input/' + input_file
;    data_path  = '/home/slisgo/divimp/results/'
    input_file = '/$HOME/fuse/input/' + input_file
;    data_path  = '/home/ITER/lisgos/fuse_data/results/'   ; + family + '/' + child + '/'
;    input_file = '/home/ITER/lisgos/divimp/data/' + input_file
    data_path  = '/$HOME/divimp/results/'

  ENDELSE


;  PRINT,'INPUT FILE  ',input_file
;  PRINT,'RESULTS DIR ',path

  status = cortex_LoadPlotData(case_name,input_file,plot_array)
  IF (status EQ -1) THEN cortex_Exit, nargs

  nplots = N_ELEMENTS(TAG_NAMES(plot_array))

; Check to see whether to draw more than one plot on a page:
  frame = MAKE_ARRAY(nplots,/LONG,VALUE=0)
  FOR iplot = 1, nplots DO BEGIN
    plot = cortex_ExtractStructure(plot_array,iplot)
    frame[iplot-1] = plot.frame
  ENDFOR
  FOR iplot = 1, nplots DO BEGIN
    i1 = 1
    i2 = nplots
    FOR i = iplot-1, 1     , -1 DO IF (i1 EQ 1      AND frame[i-1] EQ 1) THEN i1 = i + 1  ; Count backward to identify start of plot group
    FOR i = iplot  , nplots, +1 DO IF (i2 EQ nplots AND frame[i-1] EQ 1) THEN i2 = i      ; Count forward
    IF (i1 NE iplot OR i2 NE iplot) THEN BEGIN
      frac1 = FLOAT(iplot - i1    ) / FLOAT(i2 - i1 + 1)
      frac2 = FLOAT(iplot - i1 + 1) / FLOAT(i2 - i1 + 1)
      tag = 'data' + STRTRIM(STRING(iplot),2)
      index = WHERE(TAG_NAMES(plot_array) EQ STRUPCASE(tag), count)
      IF (count NE 0) THEN BEGIN
        plot_array.(index[0]).frame_bnds = [frac1,frac2]
print,'frac',frac1,frac2
      ENDIF ELSE BEGIN
        PRINT, 'ERROR cortex_GeneratePlots: TAG not found'
        PRINT, '  TAG = ',tag
        RETURN
      ENDELSE
    ENDIF
  ENDFOR
print,'frame',frame


; Setup PostScript printing:
  IF (ps EQ 'on') THEN BEGIN
    CASE dir_structure OF
      0: file = data_path + case_name + '.idl.ps'
      1: file = data_path + family + case_name + '.idl.ps'
      2: file = data_path + family + child + case_name + '.idl.ps'
    ENDCASE
    PRINT, 'POSTSCRIPT FILE=',file
    PSOPEN, filename=file
  ENDIF

;  HELP,plot_array,/struct
;
; ----------------------------------------------------------------------  
  FOR iplot = 1, nplots DO BEGIN

;    HELP,plot_array,/struct

    plot = cortex_ExtractStructure(plot_array,iplot)

    IF (plot.data_path NE 'default') THEN path = plot.data_path ELSE path = data_path

;    PRINT,' '
;    PRINT,' ',iplot
    PRINT,'---- PLOT STRUCTURE:',iplot,' ',plot.tag
;    HELP,plot,/struct    

    IF (plot.option EQ 0) THEN CONTINUE

;    PRINT,'NCASES=',N_ELEMENTS(plot.case_name),N_ELEMENTS(WHERE(plot.case_name NE 'unknown'))

    ncase = N_ELEMENTS(WHERE(plot.case_name NE 'unknown'))

    status = 0
;
;   For more than one plot on a page, check if the spacing has been adjusted by the
;   previous plot:
;
print,'adjusting befor',plot.frame_bnds,plot.frame_bnds[1]-plot.frame_bnds[0]

    IF (iplot GT 1) THEN BEGIN
      plot_last = cortex_ExtractStructure(plot_array,iplot-1)

print,'adjusting',plot_last.frame_bnds[1],plot.frame_bnds[0]

      IF (plot_last.frame EQ 0 AND plot_last.frame_bnds[1] NE plot.frame_bnds[0]) THEN BEGIN

        frame_delta = plot.frame_bnds[1] - plot.frame_bnds[0]

        plot.frame_bnds[0] = plot_last.frame_bnds[1] 

        ; Check if the existing frame width should be preserved because the
        ; next plot is on this page and is of the same type:
        IF (iplot LT nplots) THEN BEGIN
          plot_next = cortex_ExtractStructure(plot_array,iplot+1)        
          IF (plot.frame EQ 0 AND plot.tag EQ plot_next.tag) THEN   $
            plot.frame_bnds[1] = plot.frame_bnds[0] + frame_delta
        ENDIF 
        IF (iplot EQ nplots) THEN BEGIN
          IF (plot.tag EQ plot_last.tag) THEN   $
            plot.frame_bnds[1] = plot.frame_bnds[0] + frame_delta
        ENDIF 

      ENDIF
    ENDIF

print,'adjusting after',plot.frame_bnds,plot.frame_bnds[1]-plot.frame_bnds[0]
;
;
;
    CASE plot.tag OF
;     ------------------------------------------------------------------
      'PLOT 2D LOS INTEGRAL': BEGIN
        CASE plot.option OF
;         --------------------------------------------------------------         
          1: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              dummy = WHERE(plot.data_file NE 'unknown',count)
              FOR j = 0, count-1 DO BEGIN
                integral = cortex_LoadIntegrals(file_path + plot.data_file[j] + '_signal')
                name = 'data' + STRING(j+1,FORMAT='(I0)')
                IF (j EQ 0) THEN integral_array = CREATE_STRUCT(               name,integral) ELSE  $
                                 integral_array = CREATE_STRUCT(integral_array,name,integral)
              ENDFOR
              plot_data = { integral : integral_array, dummy : 1.0 }  ; Dummy is there to calm down the error check in ExtractStructure, which is lame, and needs fixing...
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
;            status = cortex_PlotImage(plot, data_array, ps=ps)
            END
;         --------------------------------------------------------------         
          ELSE: BEGIN  
            PRINT, 'ERROR cortex_GeneratePlots: Unrecognised image plot option'
            PRINT, '  OPTION = ',option
            status = -1
            END
        ENDCASE
        END
;     ------------------------------------------------------------------
      'PLOT 1D LOS INTEGRAL': BEGIN
        CASE plot.option OF
;         --------------------------------------------------------------         
          1: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              ; Load core impurtiy concentration data, if required:
              IF (plot.concentration NE 0.0) THEN  $
                   core = cortex_LoadCoreProfiles(file_path + 'idl.core_impurities')  $
              ELSE core = -1
              j = WHERE(plot.data_file NE 'unknown',count)
              FOR j = 0, count-1 DO BEGIN
                integral = cortex_LoadIntegrals(file_path + plot.data_file[j] + '_signal')
                name = 'data' + STRING(j+1,FORMAT='(I0)')
                IF (j EQ 0) THEN integral_array = CREATE_STRUCT(               name,integral) ELSE  $
                                 integral_array = CREATE_STRUCT(integral_array,name,integral)
              ENDFOR
              plot_data = { integral : integral_array, core : core }
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotIntegrals(plot, data_array, ps=ps)
            END
;         --------------------------------------------------------------         
          2: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              count = 1
              FOR j = 0, count-1 DO BEGIN
                profile = cortex_LoadProfile(file_path + plot.data_file[0] + plot.id)
                name = 'data' + STRING(j+1,FORMAT='(I0)')
                IF (j EQ 0) THEN profile_array = CREATE_STRUCT(              name,profile) ELSE  $
                                 profile_array = CREATE_STRUCT(profile_array,name,profile)
              ENDFOR
              plot_data = { profile : profile_array, dummy : 1.0 }  ; Dummy is there to calm down the error check in ExtractStructure, which is lame, and needs fixing...
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotIntegrals(plot, data_array, ps=ps)
            END
;         --------------------------------------------------------------         
          ELSE: BEGIN  
            PRINT, 'ERROR cortex_GeneratePlots: Unrecognised integral plot option'
            PRINT, '  OPTION = ',option
            status = -1
            END
        ENDCASE
        END
;     ------------------------------------------------------------------
      'PLOT 1D WALL PROFILE': BEGIN
        CASE plot.option OF
;         --------------------------------------------------------------         
          1: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              wall = cortex_LoadWallProfiles(file_path + plot.data_file[0])
              plot_data = { wall : wall, dummy : 1.0 }  ; Dummy is there to calm down the error check in ExtractStructure, which is lame, and needs fixing...
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotWallProfiles(plot, data_array, ps=ps)
            END
;         --------------------------------------------------------------         
          2: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              wall = cortex_LoadWallProfiles(file_path + plot.data_file[0])
              plot_data = { wall : wall, dummy : 1.0 }  ; Dummy is there to calm down the error check in ExtractStructure, which is lame, and needs fixing...
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotWallProfiles(plot, data_array, ps=ps)
            END
;         --------------------------------------------------------------         
          3: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              IF (i EQ 0 AND plot.show_grid EQ 1) THEN BEGIN
                grid_data = cortex_LoadFluidGrid(file_path + plot.data_file[3])
                wall_data = cortex_LoadWall     (file_path + plot.data_file[4])
                IF (plot.annotate_n NE 0) THEN  $
                  annotate = cortex_LoadAnnotationData(plot,file) ELSE annotate = 0 
              ENDIF
              wall    = cortex_LoadWallProfiles_EIRENE(file_path + plot.data_file[0])
              summary = cortex_LoadDIVIMPSummary      (file_path + plot.data_file[1])
              core    = cortex_LoadCoreProfiles       (file_path + plot.data_file[2])
              plot_data = { wall : wall, summary : summary, core : core }
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotWallProfiles(plot, data_array, grid=grid_data, wall=wall_data, annotate=annotate, ps=ps)
            END
;         --------------------------------------------------------------         
          4: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              wall = cortex_LoadWallProfiles(file_path + plot.data_file[0])
              plot_data = { wall : wall, dummy : 1.0 }  ; Dummy is there to calm down the error check in ExtractStructure, which is lame, and needs fixing...
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotWallProfiles(plot, data_array, ps=ps)
            END
;         --------------------------------------------------------------         
          5: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              IF (i EQ 0 AND plot.show_grid EQ 1) THEN BEGIN
                grid_data = cortex_LoadFluidGrid(file_path + plot.data_file[2])
                wall_data = cortex_LoadWall     (file_path + plot.data_file[3])
                IF (plot.annotate_n NE 0) THEN  $
                  annotate = cortex_LoadAnnotationData(plot,file) ELSE annotate = 0 
              ENDIF
              wall    = cortex_LoadWallErosion_DIVIMP(file_path + plot.data_file[0])
              summary = cortex_LoadDIVIMPSummary     (file_path + plot.data_file[1])
              plot_data = { wall : wall, summary : summary }
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotWallProfiles(plot, data_array, grid=grid_data, wall=wall_data, annotate=annotate, ps=ps)
            END
;         --------------------------------------------------------------         
          6: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              IF (i EQ 0 ) THEN BEGIN
                wall_data = cortex_LoadWall(file_path + plot.data_file[4])
                IF (plot.show_grid EQ 1) THEN BEGIN
                  grid_data = cortex_LoadFluidGrid(file_path + plot.data_file[3])
                  IF (plot.annotate_n NE 0) THEN  $
                    annotate = cortex_LoadAnnotationData(plot,file) ELSE annotate = 0 
                ENDIF
              ENDIF
              sputter = cortex_LoadSputterData_DIVIMP(file_path + plot.data_file[0])
              erosion = cortex_LoadWallErosion_DIVIMP(file_path + plot.data_file[1])
              summary = cortex_LoadDIVIMPSummary     (file_path + plot.data_file[2])
              plot_data = { sputter : sputter, erosion : erosion, summary : summary }
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotWallProfiles(plot, data_array, grid=grid_data, wall=wall_data, annotate=annotate, ps=ps)
            END
;         --------------------------------------------------------------         
          ELSE: BEGIN  
            PRINT, 'ERROR cortex_GeneratePlots: Unrecognised wall plot option'
            PRINT, '  OPTION = ',option
            status = -1
            END
        ENDCASE
        END
;     ------------------------------------------------------------------
      'PLOT 1D PARALLEL PROFILE': BEGIN
        CASE plot.option OF
;         --------------------------------------------------------------         
          1: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              plasma = cortex_LoadPlasmaData(file_path + plot.data_file[0])
              source = cortex_LoadSourceData(file_path + plot.data_file[1])
              eirene = cortex_LoadEireneData(file_path + plot.data_file[2])
              target = cortex_LoadTargetData(file_path + plot.data_file[3])
              plot_data = { plasma : plasma, source : source, eirene : eirene, target : target }
              IF (plot.nodes) THEN BEGIN
                node = cortex_LoadNodeData  (path + plot.case_name[i] + '.' + plot.data_file[4])
                plot_data = CREATE_STRUCT(plot_data, 'node', node) 
              ENDIF
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            FOR i = 0, N_ELEMENTS(plot.tubes)-1 DO BEGIN
              IF (plot.tubes[i] EQ 0) THEN CONTINUE
              tube = plot.tubes[i]
              IF (tube LT 0) THEN BEGIN
                print,'automatic tube selection not working yet, need to load grid data'
                stop
              ENDIF
help,plot.tubes
help,tube
print,tube
              status = cortex_PlotParallelProfiles(plot, tube, data_array, ps=ps)
              IF (i LT N_ELEMENTS(plot.tubes)-1) THEN ERASE
            ENDFOR
            END
;         --------------------------------------------------------------         
          2: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              plasma = cortex_LoadPlasmaData(file_path + plot.data_file[0])
              target = cortex_LoadTargetData(file_path + plot.data_file[1])
              plot_data = { plasma : plasma, target : target }
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            FOR i = 0, N_ELEMENTS(plot.tubes)-1 DO BEGIN
              IF (plot.tubes[i] EQ 0) THEN CONTINUE
              tube = plot.tubes[i]
              IF (tube LT 0) THEN BEGIN
                print,'automatic tube selection not working yet, need to load grid data'
                stop
              ENDIF
              status = cortex_PlotParallelProfiles(plot, tube, data_array, ps=ps)
              IF (i LT N_ELEMENTS(plot.tubes)-1) THEN ERASE
            ENDFOR
            END
;         --------------------------------------------------------------         
          100: BEGIN
            tube = plot.tubes
            IF (cortex_CheckTubes(tube,plot.tag,plot.option)) THEN BREAK
            file_path = path + plot.case_name[0] + '.'
            eirene = cortex_LoadEIRENEImpurityData        (file_path + plot.data_file[0])
            divimp = cortex_LoadDIVIMPImpurityData_Density(file_path + plot.data_file[1])
            plot_data = { eirene : eirene, divimp : divimp }
            data_array = CREATE_STRUCT('data1',plot_data)
            status = cortex_PlotParallelProfiles(plot, tube, data_array, ps=ps)
            END
;         --------------------------------------------------------------         
          101: BEGIN
            tube = plot.tubes
            IF (cortex_CheckTubes(tube,plot.tag,plot.option)) THEN BREAK
            file_path = path + plot.case_name[0] + '.'
            eirene = cortex_LoadEIRENEImpurityData           (file_path + plot.data_file[0])
            divimp = cortex_LoadDIVIMPImpurityData_Ionisation(file_path + plot.data_file[1])
            plot_data = { eirene : eirene, divimp : divimp }
            data_array = CREATE_STRUCT('data1',plot_data)
            status = cortex_PlotParallelProfiles(plot, tube, data_array, ps=ps)
            END
;         --------------------------------------------------------------         
          ELSE: BEGIN  
            PRINT, 'ERROR cortex_GeneratePlots: Unrecognised parallel plot option'
            PRINT, '  OPTION = ',option
            status = -1
            END
        ENDCASE
        END
;     ------------------------------------------------------------------
      'PLOT 1D EIRENE SPECTRUM': BEGIN
 
        found_first = 0

        FOR i = 0, ncase-1 DO BEGIN
          file = path + plot.case_name[i] + '.idl.eirene_spectra'
          print, 'spectra file ',file
          spectra = cortex_LoadSpectraInfo(file)          
; help,spectra,/struct
          plot_data = CREATE_STRUCT('case_name',plot.case_name[i],'spectra',spectra) ; if this assignment changes, i.e. 'spectra' is not the second element, then need to change the val.(1) assignment in cortex_1D_spectrum

          ichord = FIX(plot.chord_index)

          ; Load profile data for this chord:
          profile_name = 'profile_' + plot.id + '_' + STRING(ichord,FORMAT='(I3.3)')
          file = path + plot.case_name[i] + '.idl.' + profile_name
          print, '  profile file ',file
          profile = cortex_LoadProfile(file)
          profile = CREATE_STRUCT(profile,'spectra_chord',ichord)
          plot_data = CREATE_STRUCT(plot_data,profile_name,profile)

          found  = 0
          nchord = MAX(spectra.index)
          FOR ichord = 1L, nchord DO BEGIN
;            IF (NOT cortex_CheckIndex(ichord,plot.chord_index)) THEN CONTINUE

            j = WHERE(spectra.index EQ ichord, count)

            ; Need to compare the data in PROFILE with the sub-set of data in SPECTRA
            ; to match the chords specified for EIRENE with the selected chord from RAY:

            ; First test: see if the number of elements are the same:
            IF (count NE N_ELEMENTS(profile.path)) THEN CONTINUE
            ; Second test: see if the path data are the same:
            found = 1
            FOR k = 0, count-1 DO BEGIN
              IF ((ABS(profile.path [k] - spectra.path [j[k]]) GT 5.0E-6) OR  $
                  (ABS(profile.delta[k] - spectra.delta[j[k]]) GT 5.0E-6)) THEN found = 0
            ENDFOR
            IF (NOT found) THEN CONTINUE

            IF (found_first) THEN BEGIN
              PRINT, 'ERROR cortex_GeneratePlots: More than one spectra chord identified'
              status = -1
            ENDIF
            found_first = 1

            FOR k = 0L, count-1 DO BEGIN
              icell = spectra.cell[j[k]]
              IF (NOT cortex_CheckIndex(icell,plot.cell_index)) THEN CONTINUE              
              name = 'spectrum_' + STRING(j[k]+1,FORMAT='(I4.4)') + '_' +  $
                                   STRING(icell ,FORMAT='(I7.7)')
              file = path + plot.case_name[i] + '.idl.eirene_' + name + '_sum'
;              print, '  spectrum file ',file
              spectrum = cortex_LoadEireneSpectrum(file)
              spectrum = CREATE_STRUCT(spectrum,'profile_name' ,profile_name,  $
                                                'profile_index',k+1         ,  $
                                                'chord_index'  ,ichord      ,  $
                                                'spectra_index',j[k]+1      ,  $
                                                'spectra_cell' ,icell)
              plot_data = CREATE_STRUCT(plot_data,name,spectrum)
            ENDFOR

          ENDFOR
;          plot_data = cortex_LoadEnergySpectrum(file)
          IF (found_first) THEN BEGIN
            name = 'data' + STRING(i+1,FORMAT='(I0)')
            IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                             data_array = CREATE_STRUCT(data_array,name,plot_data)
          ENDIF
        ENDFOR
        IF (NOT found_first) THEN BEGIN
          PRINT, 'ERROR cortex_GeneratePlots: Spectra chord not identified'
          status = -1        
        ENDIF ELSE  $
          status = cortex_PlotEireneSpectra(plot,data_array, ps=ps)
        END
;     ------------------------------------------------------------------
      'PLOT 1D TARGET PROFILE': BEGIN
        FOR i = 0, ncase-1 DO BEGIN
          file = path + plot.case_name[i] + '.'
          target   = cortex_LoadTargetData      (file + plot.data_file[0])
          midplane = cortex_LoadMidplaneProfiles(file + plot.data_file[1])
          plot_data = { target : target, midplane : midplane }
          name = 'data' + STRING(i+1,FORMAT='(I0)')
          IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                           data_array = CREATE_STRUCT(data_array,name,plot_data)
        ENDFOR
        status = cortex_PlotTargetProfiles(plot,data_array, ps=ps)
        END
;     ------------------------------------------------------------------
      'PLOT 1D RADIAL PROFILE': BEGIN
        CASE plot.option OF
;         --------------------------------------------------------------         
          1: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file = path + plot.case_name[i] + '.'
              mid = cortex_LoadMidplaneProfiles(file + plot.data_file[0])
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN mid_array = CREATE_STRUCT(          name,mid) ELSE  $
                               mid_array = CREATE_STRUCT(mid_array,name,mid)
            ENDFOR
;            help, mid_array, /struct
            status = cortex_PlotMidplaneProfiles(plot, mid_array, ps=ps)
            END
;         --------------------------------------------------------------         
          2: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              grid   = cortex_LoadFluidGrid (path + plot.case_name[i] + '.' + plot.data_file[0])
              wall   = cortex_LoadWall      (path + plot.case_name[i] + '.' + plot.data_file[1])
              plasma = cortex_LoadPlasmaData(path + plot.case_name[i] + '.' + plot.data_file[2])
;              target = cortex_LoadTargetData(path + plot.case_name[i] + '.' + plot.data_file[3])
              target = -1
              FOR j = 1, plot.line_seg_n DO BEGIN
                inter = cortex_SliceGrid(plot.line_seg[0+(j-1)*4:3+(j-1)*4],grid,status)
                IF (status EQ -1) THEN BEGIN
                  plot.line_seg_s = -1
                ENDIF ELSE BEGIN
                  name = 'data' + STRING(j,FORMAT='(I0)')
                  IF (j EQ 1) THEN inter_array = CREATE_STRUCT(            name,inter) ELSE  $
                                   inter_array = CREATE_STRUCT(inter_array,name,inter)
                ENDELSE
              ENDFOR
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              plot_data = { grid : grid, wall : wall, inter : inter_array, plasma : plasma, target : target }
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotRadialProfile(plot, data_array, ps=ps)
            END
;         --------------------------------------------------------------         
          3: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              core = cortex_LoadCoreProfiles(path + plot.case_name[i] + '.' + plot.data_file[0])
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN core_array = CREATE_STRUCT(           name,core) ELSE  $
                               core_array = CREATE_STRUCT(core_array,name,core)
            ENDFOR
            status = cortex_PlotCoreProfiles(plot, core_array, ps=ps)
            END
;         --------------------------------------------------------------         
          4: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              grid   = cortex_LoadFluidGrid (path + plot.case_name[i] + '.' + plot.data_file[0])
              wall   = cortex_LoadWall      (path + plot.case_name[i] + '.' + plot.data_file[1])
              plasma = cortex_LoadPlasmaData(path + plot.case_name[i] + '.' + plot.data_file[2])
              target = -1
              FOR j = 1, plot.line_seg_n DO BEGIN
                inter = cortex_SliceGrid(plot.line_seg[0+(j-1)*4:3+(j-1)*4],grid,status)
                IF (status EQ -1) THEN BEGIN
                  plot.line_seg_s = -1
                ENDIF ELSE BEGIN
                  name = 'data' + STRING(j,FORMAT='(I0)')
                  IF (j EQ 1) THEN inter_array = CREATE_STRUCT(            name,inter) ELSE  $
                                   inter_array = CREATE_STRUCT(inter_array,name,inter)
                ENDELSE
              ENDFOR

              name = 'data' + STRING(i+1,FORMAT='(I0)')
              plot_data = { grid : grid, wall : wall, inter : inter_array, plasma : plasma, target : target }
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotRadialProfile(plot, data_array, ps=ps)
            END
;         --------------------------------------------------------------         
          5: BEGIN  ; Average plasma quantities for each flux-tube
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              plasma = cortex_LoadPlasmaData(file_path + plot.data_file[1])
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              plot_data = { plasma : plasma, dummy : -1 }
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotRadialProfile(plot, data_array, ps=ps)
            END
;         --------------------------------------------------------------         
          ELSE: BEGIN  
            PRINT, 'ERROR cortex_GeneratePlots: Unrecognised radial plot option'
            PRINT, '  OPTION = ',option
            status = -1
            END
        ENDCASE

        END
;     ------------------------------------------------------------------
      'PLOT 1D PEDESTAL MODEL': BEGIN
        FOR i = 0, ncase-1 DO BEGIN
          pedestal = cortex_LoadPedestalModel(path + plot.case_name[i] + '.' + plot.data_file[0])
          plot_data = { pedestal : pedestal, dummy : 0 }
          name = 'data' + STRING(i+1,FORMAT='(I0)')
          IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                           data_array = CREATE_STRUCT(data_array,name,plot_data)
        ENDFOR
        status = cortex_PlotPedestalModel(plot, data_array, ps=ps)
        END
;     ------------------------------------------------------------------
      'PLOT 1D SUMMARY': BEGIN
        CASE plot.option OF
;         --------------------------------------------------------------         
          1: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              IF (plot.case_name[i] EQ 'none') THEN BEGIN
                core    = CREATE_STRUCT('version',0.0,'file','none')
                summary = CREATE_STRUCT('version',0.0,'file','none')
              ENDIF ELSE BEGIN
                file = path + plot.case_name[i] + '.'
                core = cortex_LoadCoreProfiles (file + plot.data_file[0])
                IF (plot.warnings) THEN summary = cortex_LoadDIVIMPSummary(file + plot.data_file[1]) ELSE summary = 0
;                plot_data = { core:core, summary:summary }
              ENDELSE
              plot_data = { core:core, summary:summary }
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotSummary(plot, data_array, ps=ps)
            END
          ELSE: BEGIN  
            PRINT, 'ERROR cortex_GeneratePlots: Unknown summary plot option'
            PRINT, '  OPTION = ',grid.option
            status = -1
            END
        ENDCASE
        END
;     ------------------------------------------------------------------
      'PLOT 2D FLUID GRID': BEGIN
        IF (plot.load AND KEYWORD_SET(grid_save)) THEN BEGIN
          plot_save = plot
          grid = grid_save.grid     
          wall = grid_save.wall     
          node = grid_save.node     
          plot = grid_save.plot
          plot.frame      = plot_save.frame
          plot.frame_bnds = plot_save.frame_bnds
          plot.save = plot_save.save  
          plot.load = plot_save.load  
          plot.zoom = plot_save.zoom  
          plot.show_n       = plot_save.show_n
          plot.show_tube    = plot_save.show_tube 
          plot.show_colour  = plot_save.show_colour
          plot.aspect_ratio = plot_save.aspect_ratio  
          IF (plot_save.annotate_n NE 0) THEN BEGIN            
            n = plot_save.annotate_n
            i = plot.annotate_n
            j = i + n - 1
            plot.annotate_n           = plot.annotate_n + n
            plot.annotate_code  [i:j] = plot_save.annotate_code  [0:n-1]
	    plot.annotate_colour[i:j] = plot_save.annotate_colour[0:n-1]
	    plot.annotate_file  [i:j] = plot_save.annotate_file  [0:n-1]
	    plot.annotate_step  [i:j] = plot_save.annotate_step  [0:n-1]
	    plot.annotate_label [i:j] = plot_save.annotate_label [0:n-1]
          ENDIF
          IF (plot.annotate_n NE 0) THEN  $
            annotate = cortex_LoadAnnotationData(plot,file) ELSE annotate = 0 
        ENDIF ELSE BEGIN
          file = path + plot.case_name[0] + '.'
          IF (plot.equ NE 'default') THEN grid = grid_ReadEQUFile    (plot.equ                ) ELSE  $
                                          grid = cortex_LoadFluidGrid(file + plot.data_file[0])
          wall = cortex_LoadWall(file + plot.data_file[1])
          IF (plot.nodes) THEN node = cortex_LoadNodeData(file + plot.data_file[2]) ELSE node = 0
          IF (plot.annotate_n NE 0) THEN  $
            annotate = cortex_LoadAnnotationData(plot,file) ELSE annotate = 0 
        ENDELSE
        IF (plot.save) THEN BEGIN
          grid_save = {            $
            grid     : grid     ,  $
            wall     : wall     ,  $
            node     : node     ,  $
            plot     : plot     }
        ENDIF ELSE BEGIN
          status = cortex_PlotFluidGrid(plot, grid, wall, node, annotate, 'main', 'full', ps=ps)

          frame_bnds = plot.frame_bnds          
          IF (plot.load) THEN BEGIN
            plot            = plot_save
            plot.frame_bnds = frame_bnds
          ENDIF
        ENDELSE
        END
;     ------------------------------------------------------------------
      'PLOT 2D FLUID GRID - DEBUG': BEGIN
        grid = cortex_LoadFluidGrid_Debug(path + plot.case_name[0] + '.' + plot.data_file[0])
        wall = cortex_LoadFluidWall_Debug(path + plot.case_name[0] + '.' + plot.data_file[1])
;        wall   = {               $
;          n      : 0          ,  $
;          class  : 0          ,  $
;          group  : 0          ,  $
;          index  : 0          ,  $
;          tube   : 0          ,  $
;          target : 0          ,  $
;          v1     : [0.0,0.0]  ,  $
;          v2     : [0.0,0.0]     }
        status = cortex_PlotFluidGrid(plot, grid, wall, 0, 0, 'main', 'full', ps=ps)
        END
;     ------------------------------------------------------------------
      'PLOT 2D CONTOUR': BEGIN
        file_path = path + plot.case_name[0] + '.'
        grid = cortex_LoadFluidGrid(file_path + plot.data_file[0])
        wall = cortex_LoadWall     (file_path + plot.data_file[1]) 
        IF (plot.annotate_n NE 0) THEN  $
          annotate = cortex_LoadAnnotationData(plot,file_path) ELSE annotate = 0 
        IF (plot.option GE   1 AND plot.option LE  99) THEN plot_data = cortex_LoadPlasmaData(file_path + plot.data_file[2])
        IF (plot.option GE 100 AND plot.option LE 199) THEN plot_data = cortex_LoadSourceData(file_path + plot.data_file[3])
        IF (plot.option GE 200 AND plot.option LE 299) THEN plot_data = cortex_LoadEIRENEData(file_path + plot.data_file[4])        
        IF (plot.option GE 300 AND plot.option LE 399) THEN plot_data = cortex_LoadEIRENEImpurityData(file_path + plot.data_file[5])        
        IF (plot.option GE 400 AND plot.option LE 400) THEN plot_data = cortex_LoadDIVIMPImpurityData_Density(file_path + plot.data_file[6])        
        IF (plot.option LE   0 OR  plot.option GE 401) THEN BEGIN
          PRINT, 'ERROR cortex_GeneratePlots: Unrecognised 2D contour plot option'
          PRINT, '  PLOT TAG = ',plot.tag
          cortex_Exit, nargs
        ENDIF
        status = cortex_Plot2DContour(plot, plot_data, grid, wall, annotate, ps=ps)
        END
;     ------------------------------------------------------------------
      'PLOT 2D IMAGE': BEGIN
        CASE plot.option OF
;         --------------------------------------------------------------         
          1: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              j = WHERE(plot.data_file NE 'unknown',count)
              FOR j = 0, count-1 DO BEGIN
                integral = cortex_LoadIntegrals(file_path + plot.data_file[j] + '_signal')
                name = 'data' + STRING(j+1,FORMAT='(I0)')
                IF (j EQ 0) THEN integral_array = CREATE_STRUCT(               name,integral) ELSE  $
                                 integral_array = CREATE_STRUCT(integral_array,name,integral)
              ENDFOR
              plot_data = { integral : integral_array, dummy :  0 }
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotImage(plot, data_array, ps=ps)
            END
;         --------------------------------------------------------------         
          2: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              j = WHERE(plot.data_file NE 'unknown',count)
              FOR j = 0, count-1, 2 DO BEGIN
                integral = cortex_LoadIntegrals(file_path + plot.data_file[j] + '_signal')
                integral = CREATE_STRUCT(integral,'stripe',plot.data_file[j+1])
                name = 'data' + STRING(j/2+1,FORMAT='(I0)')
                IF (j EQ 0) THEN integral_array = CREATE_STRUCT(               name,integral) ELSE  $
                                 integral_array = CREATE_STRUCT(integral_array,name,integral)
              ENDFOR
              plot_data = { integral : integral_array, dummy :  0 }
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotImage(plot, data_array, ps=ps)
            END
;         --------------------------------------------------------------         
          ELSE: BEGIN  
            PRINT, 'ERROR cortex_GeneratePlots: Unrecognised image plot option'
            PRINT, '  OPTION = ',option
            status = -1
            END
        ENDCASE
        END
;     ------------------------------------------------------------------
      'PLOT 3D TEST': BEGIN
        file_path = path + plot.case_name[0] + '.'
        plot_data = cortex_LoadTetrahedronData(file_path + plot.data_file[0])
        help,plot_data,/struct
        tet = plot_data
        SAVE, filename='tetrahedron_test.sav', tet
;        status = cortex_Plot2DContour(plot, plot_data, ps=ps)
        END
;     ------------------------------------------------------------------
      ELSE: BEGIN
        PRINT, 'ERROR cortex_GeneratePlots: Unrecognised plot'
        PRINT, '  PLOT TAG = ',plot.tag
        cortex_Exit, nargs
        END
    ENDCASE
;
;   Update PLOT structure in the PLOT_ARRAY structure in case data has been
;   modified or added when drawing the plot:
;   --------------------------------------------------------------------
    tag = 'data' + STRTRIM(STRING(iplot),2)
    index = WHERE(TAG_NAMES(plot_array) EQ STRUPCASE(tag), count)
    plot_array.(index[0]) = plot
;
;
;
    IF (plot.frame EQ 1) THEN ERASE

    IF (status NE 0) THEN cortex_Exit, nargs

  ENDFOR

  PSCLOSE

END

;
; ======================================================================
;
PRO cortex_Main, args

  IF (NOT KEYWORD_SET(args)) THEN args = 'none'

  PRINT, N_ELEMENTS(args),args

  cortex_GeneratePlots, args

 ;  IF (N_ELEMENTS(args) EQ 4) THEN EXIT, status=0
END
;
;
;
