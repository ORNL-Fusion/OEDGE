;
; example : cortex_main,['t-new-0000a','i-osm.ctx']
;
;
; 3d gas puff: cortex_main,['i-cxd-0010b','users_guide_dev.ctx']
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

  COMMON in_commons, param, master
  COMMON options, dir_structure
  
  IF (N_ELEMENTS(master) NE 0) THEN UNDEFINE, master

  option = 3

  ps = 'on'

  nargs = N_ELEMENTS(args)

  SWITCH STRTRIM(args[0],2) OF
;   ------------------------------------------------------------         
    'sd_no': 
    'sd_yes': BEGIN
      IF (args[0] EQ 'sd_yes') THEN dir_structure = 1 ELSE dir_structure = 0
      case_name  = args[1]
      input_file = args[3] + '/' + args[2]
      data_path  = args[4] + '/'
      BREAK
      END
;   ------------------------------------------------------------         
    'windows': BEGIN
      dir_structure = 0
      case_name  = args[1]
      input_file = args[2]
      input_file = 'c:\Users\Steven\Desktop\share\fuse\input\' + input_file
      data_path  = 'c:\Users\Steven\Desktop\share\results\'
      BREAK
      END
;   ------------------------------------------------------------         
    ELSE: BEGIN  
      dir_structure = 1 ; 0 ; 1 
      case_name  = args[0]
      input_file = args[1]
      input_file = '/$HOME/fuse/input/' + input_file
      data_path  = '/$HOME/divimp/results/'
      BREAK
      END
;   ------------------------------------------------------------         
  ENDSWITCH

;  IF (args[0] EQ 'sd_no' OR args[0] EQ 'sd_yes' ) THEN BEGIN
;    IF (args[0] EQ 'sd_yes') THEN dir_structure = 1
;    case_name  = args[1]
;    input_file = args[3] + '/' + args[2]
;    data_path  = args[4] + '/'
;;print,case_name
;;print,input_file
;;print,data_path
;;stop
;  ENDIF ELSE BEGIN
;    dir_structure = 1 ; 0 ; 1 
;    case_name  = args[0]
;    input_file = args[1]
;    input_file = '/$HOME/fuse/input/' + input_file
;    data_path  = '/$HOME/divimp/results/'
;;print,case_name
;;print,input_file
;;print,data_path
;;stop
;  ENDELSE

  IF (dir_structure GE 1) THEN BEGIN
    str = STRSPLIT(case_name,'-',/EXTRACT)
    family = str[0] + '-' + str[1] + '/'
;    family = STRMID(case_name,0,5) + '/'
    child  = STRMID(case_name,0,7) + '/'
  ENDIF


;  case_name  = args[0]
;  input_file = args[1]
;  IF (nargs EQ 4) THEN BEGIN
;    dir_structure = 0
;    family = ''
;    child  = ''
;    input_file = args[2] + '/' + input_file
;    data_path  = args[3] + '/'
;  ENDIF ELSE BEGIN
;    dir_structure = 1 ; 1
;    str = STRSPLIT(case_name,'-',/EXTRACT)
;    family = str[0] + '-' + str[1] + '/'
;;    family = STRMID(case_name,0,5) + '/'
;    child  = STRMID(case_name,0,7) + '/'
;;    input_file = '/home/slisgo/fuse/input/' + input_file
;;    data_path  = '/home/slisgo/divimp/results/'
;    input_file = '/$HOME/fuse/input/' + input_file
;;    data_path  = '/home/ITER/lisgos/fuse_data/results/'   ; + family + '/' + child + '/'
;;    input_file = '/home/ITER/lisgos/divimp/data/' + input_file
;    data_path  = '/$HOME/divimp/results/'
;  ENDELSE


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

;        IF (plot.frame EQ 0) THEN BEGIN
;          frac1 = 0.0
;          frac2 = 0.5
;        ENDIF ELSE BEGIN
;          frac1 = 0.0
;          frac2 = 1.0
;        ENDELSE

;        frac1 = 0.24  ; screwing around with 3D plots? 
;        frac2 = 0.43

        plot_array.(index[0]).frame_bnds = [frac1,frac2]

        print,'frac',frac1,frac2,i1,i2,plot.frame

      ENDIF ELSE BEGIN
        PRINT, 'ERROR cortex_GeneratePlots: TAG not found'
        PRINT, '  TAG = ',tag
        RETURN
      ENDELSE
    ENDIF
  ENDFOR
print,'frame',frame


  ps_active       = 0
  ps_active_local = 0

;  HELP,plot_array,/struct
;
; ----------------------------------------------------------------------  
  FOR iplot = 1, nplots DO BEGIN


; Setup PostScript printing:
    IF (ps_active EQ 0 and plot.ps_local EQ 0) THEN BEGIN
      IF (ps EQ 'on') THEN BEGIN
        CASE dir_structure OF
          0: file = data_path + case_name + '.idl.ps'
          1: file = data_path + family + case_name + '.idl.ps'
          2: file = data_path + family + child + case_name + '.idl.ps'
        ENDCASE
        PRINT, 'POSTSCRIPT FILE=',file
        PSOPEN, filename=file
        ps_active = 1
      ENDIF
    ENDIF


;    HELP,plot_array,/struct

    plot = cortex_ExtractStructure(plot_array,iplot)

    IF (plot.data_path NE 'default') THEN path = plot.data_path ELSE path = data_path

    PRINT,'---- PLOT STRUCTURE:',iplot,' ',plot.tag

    IF (plot.option EQ 0 OR plot.skip EQ 1) THEN CONTINUE

    ncase = N_ELEMENTS(WHERE(plot.case_name NE 'unknown'))

    status = 0
;
;   For more than one plot on a page, check if the spacing has been adjusted by the
;   previous plot:
;
    IF (iplot GT 1) THEN BEGIN
      plot_last = cortex_ExtractStructure(plot_array,iplot-1)
      IF (plot_last.frame EQ 0 AND plot_last.frame_bnds[1] NE plot.frame_bnds[0]) THEN BEGIN

        print,'screwing around with the frame, probably due to adjustment in cortex_3D_grid'
        print,plot_last.frame_bnds
        print,plot.frame_bnds
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


;    IF (iplot GT 1) THEN BEGIN
;      plot1 = cortex_ExtractStructure(plot_array,iplot-1)
;      print,'frame_bnds i-1',plot.frame_bnds
;    ENDIF
;    print,'frame_bnds i  ',plot.frame_bnds
;    IF (iplot LT nplots) THEN BEGIN
;      plot1 = cortex_ExtractStructure(plot_array,iplot+1)
;      print,'frame_bnds i+1',plot1.frame_bnds
;    ENDIF
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
            time_s = -1
            time_e = -1
            time_d =  1
            path_local = path
            time_suffix = ''
            IF (plot.time_span[0] NE -1) THEN BEGIN
              time_s = plot.time_span[0]
              time_e = plot.time_span[1]
              time_d = plot.time_span[2]
              dir_structure = 0
            ENDIF
            k = -1
            FOR j = time_s, time_e, time_d DO BEGIN
              FOR i = 0, ncase-1 DO BEGIN
                IF (j NE -1) THEN BEGIN
                   time_suffix = '_' + (STRTRIM(STRING(j,FORMAT='(I8.8)'),2))
                   path_local = '/$HOME/divimp/shots/iter/aps_elm/results/' + plot.case_name[i] + '/'
                ENDIF                  
                file_path = path_local + plot.case_name[i] + time_suffix + '.'
print,'file_path ',file_path
                IF (i EQ 0 AND plot.show_grid EQ 1) THEN BEGIN
                  grid_data = cortex_LoadFluidGrid(file_path + plot.data_file[4])
                  wall_data = cortex_LoadWall     (file_path + plot.data_file[5])
                  IF (plot.annotate_n NE 0) THEN  $
                    annotate = cortex_LoadAnnotationData(plot,file) ELSE annotate = 0 
                ENDIF
                wall     = cortex_LoadWallErosion_DIVIMP(file_path + plot.data_file[0])
                summary  = cortex_LoadDIVIMPSummary     (file_path + plot.data_file[1])
                flux     = cortex_LoadWallProfiles      (file_path + plot.data_file[2])
                midplane = cortex_LoadMidplaneProfiles  (file_path + plot.data_file[3])
                plot_data = { wall : wall, summary : summary, flux : flux, midplane : midplane }
                k++
                name = 'data' + STRING(k+1,FORMAT='(I0)')
                IF (k EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                                 data_array = CREATE_STRUCT(data_array,name,plot_data)
              ENDFOR
            ENDFOR
            status = cortex_PlotWallProfiles(plot, data_array, grid=grid_data, wall=wall_data, annotate=annotate, ps=ps)
            END
;         --------------------------------------------------------------         
          6: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              IF (i EQ 0 ) THEN BEGIN
                wall_data = cortex_LoadWall(file_path + plot.data_file[5])
                IF (plot.show_grid EQ 1) THEN BEGIN
                  grid_data = cortex_LoadFluidGrid(file_path + plot.data_file[4])
                  IF (plot.annotate_n NE 0) THEN  $
                    annotate = cortex_LoadAnnotationData(plot,file) ELSE annotate = 0 
                ENDIF
                IF (plot.spectra EQ 1) THEN BEGIN
                  ; 
                  strata = cortex_LoadEireneStrata(file_path + plot.data_file[6])
                  spectra_data = CREATE_STRUCT('strata',strata)
                  ; 
                  spectra = cortex_LoadSpectraInfo(file_path + plot.data_file[7])
                  spectra_data = CREATE_STRUCT(spectra_data,'spectra',spectra) 
                  ; 
                  FOR k = 0, N_ELEMENTS(spectra.index)-1 DO BEGIN
                    icell = spectra.cell[k]
                    IF (icell LT 0) THEN itag = -icell + 1  $
                                    ELSE itag =  icell
                    name = 'spectrum_' + STRING(k+1 ,FORMAT='(I4.4)') + '_' +  $
                                         STRING(itag,FORMAT='(I7.7)')
                    file = path + plot.case_name[i] + '.idl.eirene_' + name + '_sum'
                    print,'  spectrum file ',file
                    spectrum = cortex_LoadEireneSpectrum(file)
                    spectrum = CREATE_STRUCT(spectrum,'cell',icell,'tag',itag)
                    name = 'data' + STRING(k+1,FORMAT='(I0)') 
                    spectra_data = CREATE_STRUCT(spectra_data,name,spectrum)
                  ENDFOR
                ENDIF ELSE spectra_data = -1
                IF (plot.temp_inner NE 'none') THEN temp_inner = cortex_LoadSurfaceTemperatures(plot.temp_inner)  $
                                               ELSE temp_inner = -1
                IF (plot.temp_outer NE 'none') THEN temp_outer = cortex_LoadSurfaceTemperatures(plot.temp_outer)  $
                                               ELSE temp_outer = -1
              ENDIF
              sputter = cortex_LoadSputterData_DIVIMP(file_path + plot.data_file[0])
              erosion = cortex_LoadWallErosion_DIVIMP(file_path + plot.data_file[1])
              summary = cortex_LoadDIVIMPSummary     (file_path + plot.data_file[2])
              flux    = cortex_LoadWallProfiles      (file_path + plot.data_file[3])
              plot_data = { sputter : sputter, erosion : erosion, summary : summary , flux : flux, spectra : spectra_data ,  $
                            temp_inner : temp_inner , temp_outer : temp_outer }
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotWallProfiles(plot, data_array, grid=grid_data, wall=wall_data,  $
                                             annotate=annotate, ps=ps)
            END
;         --------------------------------------------------------------         
          7: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              IF (i EQ 0 AND plot.show_grid EQ 1) THEN BEGIN
                grid_data = cortex_LoadFluidGrid(file_path + plot.data_file[4])
                wall_data = cortex_LoadWall     (file_path + plot.data_file[5])
                IF (plot.annotate_n NE 0) THEN  $
                  annotate = cortex_LoadAnnotationData(plot,file) ELSE annotate = 0 
              ENDIF
              wall     = cortex_LoadWallErosion_DIVIMP(file_path + plot.data_file[0])
              flux     = cortex_LoadWallProfiles      (file_path + plot.data_file[1])
              target   = cortex_LoadTargetData        (file_path + plot.data_file[2])
              midplane = cortex_LoadMidplaneProfiles  (file_path + plot.data_file[3])
              plot_data = { wall : wall, flux : flux, target : target, midplane : midplane }
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
            status = cortex_PlotWallProfiles(plot, data_array, grid=grid_data, wall=wall_data, annotate=annotate, ps=ps)
            END
;         --------------------------------------------------------------         
          8: BEGIN

            FOR i = 0, ncase-1 DO BEGIN

              file_path = path + plot.case_name[i] + '.'

              IF (i EQ 0 AND plot.show_grid EQ 1) THEN BEGIN
                grid_data = cortex_LoadFluidGrid(file_path + plot.data_file[0])
                wall_data = cortex_LoadWall     (file_path + plot.data_file[1])
                IF (plot.annotate_n NE 0) THEN  $
                  annotate = cortex_LoadAnnotationData(plot,file) ELSE annotate = 0 
              ENDIF

              wall = { file : plot.case_name[i] }

              plot_data = { wall : wall, dummy : -1 }

              name = 'data' + STRING(i+1,FORMAT='(I0)')

              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)

            ENDFOR

            walldyn = cortex_LoadWalldyn(plot)

            status = cortex_PlotWallProfiles(plot,data_array,wdyn=walldyn,grid=grid_data, wall=wall_data, annotate=annotate, ps=ps)

            END
;         --------------------------------------------------------------         
          9: BEGIN
             ; is in PlotWallProfiles, but not here... strange...
            END
;         --------------------------------------------------------------         
          10: BEGIN
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
              status = cortex_PlotParallelProfiles(plot, tube, data_array, ps=ps)
              IF (i LT N_ELEMENTS(plot.tubes)-1) THEN ERASE
            ENDFOR
            END
;         --------------------------------------------------------------         
          2: BEGIN  ; APS 2014 
            time_s = -1
            time_e = -1
            time_d =  1
            path_local = path
            time_suffix = ''
            IF (plot.time_span[0] NE -1) THEN BEGIN
              time_s = plot.time_span[0]
              time_e = plot.time_span[1]
              time_d = plot.time_span[2]
              dir_structure = 0
            ENDIF

            FOR j = time_s, time_e, time_d DO BEGIN

              FOR i = 0, ncase-1 DO BEGIN

                IF (j NE -1) THEN BEGIN
                   time_suffix = '_' + (STRTRIM(STRING(j,FORMAT='(I8.8)'),2))
                   path_local = '/$HOME/divimp/shots/iter/aps_elm/results/' + plot.case_name[i] + '/'
                ENDIF                  
                file_path = path_local + plot.case_name[i] + time_suffix + '.'
;                file_path = path_local + plot.case_name[i] + '.'
                plasma = cortex_LoadPlasmaData(file_path + plot.data_file[0])
                target = cortex_LoadTargetData(file_path + plot.data_file[1])
                plot_data = { plasma : plasma, target : target }
                name = 'data' + STRING(i+1,FORMAT='(I0)')
                IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                                 data_array = CREATE_STRUCT(data_array,name,plot_data)
              ENDFOR

              dummy = WHERE(plot.tubes NE 0,ntube)
              IF (ntube EQ ncase) THEN BEGIN

                print,'going for it...'
                status = cortex_PlotParallelProfiles(plot, plot.tubes, data_array, ps=ps, time_slice=-1)        

              ENDIF ELSE BEGIN

                FOR i = 0, N_ELEMENTS(plot.tubes)-1 DO BEGIN
                  IF (plot.tubes[i] EQ 0) THEN CONTINUE
                  tube = plot.tubes[i]
                  IF (tube LT 0) THEN BEGIN
                    print,'automatic tube selection not working yet, need to load grid data'
                    stop
                  ENDIF
            
                  IF (plot.ps_local EQ 1 and ps_active_local EQ 0) THEN BEGIN
                    IF (ps_active EQ 1) THEN PSCLOSE
                    IF (ps EQ 'on') THEN BEGIN
                      time_suffix = '_' + (STRTRIM(STRING(j,FORMAT='(I8.8)'),2))
                      path_local  = '/$HOME/divimp/shots/iter/aps_elm/results/' + plot.case_name[i] + '/'
                      file_path   = path_local + plot.case_name[i] + '_' + plot.id + time_suffix + '.'
                      file = file_path + 'idl.ps'
                      PRINT, 'POSTSCRIPT FILE=',file
                      PSOPEN, filename=file
                      ps_active       = 0
                      ps_active_local = 1
                    ENDIF
                  ENDIF

                  IF (j NE -1) THEN time_slice = j * 10 ELSE time_slice = -1

                  status = cortex_PlotParallelProfiles(plot, tube, data_array, ps=ps, time_slice=time_slice)

                  IF (i LT N_ELEMENTS(plot.tubes)-1) THEN ERASE

                ENDFOR

              ENDELSE

              IF (plot.ps_local EQ 1 AND plot.frame EQ 1) THEN BEGIN
                ps_active_local = 0
                ERASE
                PSCLOSE
              ENDIF

            ENDFOR

            END
;         --------------------------------------------------------------         
          100: BEGIN
            tube = plot.tubes
            IF (cortex_CheckTubes(tube,plot.tag,plot.option)) THEN BREAK
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              divimp = cortex_LoadDIVIMPImpurityData_Density(file_path + plot.data_file[0])
              ;eirene = cortex_LoadEIRENEImpurityData        (file_path + plot.data_file[1])
              eirene = -1
              plot_data = { eirene : eirene, divimp : divimp }
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
;            plot_data = { eirene : eirene, divimp : divimp }
;            data_array = CREATE_STRUCT('data1',plot_data)
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
          102: BEGIN
            tube = plot.tubes
            IF (cortex_CheckTubes(tube,plot.tag,plot.option)) THEN BREAK
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              emission = cortex_LoadEmissionData(file_path + plot.data_file[0])  
              plot_data = { emission : emission, dummy : -1 }
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
;            plot_data = cortex_LoadEmissionData(file_path + plot.data_file[0])        
;            data_array = CREATE_STRUCT('data1',plot_data)
            status = cortex_PlotParallelProfiles(plot, tube, data_array, ps=ps)
            END
;         --------------------------------------------------------------         
          103: BEGIN  ; debugging PT solver
            j = -1
print,plot.iter
            FOR i = plot.iter[0], plot.iter[1], plot.iter[2] DO BEGIN
              file_path = path + plot.case_name[0] + '_' + STRING(i,FORMAT='(I8.8)') + '.'

print,file_path + plot.data_file[0]

              plasma = cortex_SolverDebugData(file_path + plot.data_file[0])
              plot_data = { plasma : plasma, dummy : -1 }
              j++
              name = 'data' + STRING(j+1,FORMAT='(I0)')
              IF (j EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                               data_array = CREATE_STRUCT(data_array,name,plot_data)
            ENDFOR
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

              PRINT, '  count           = ',count
              PRINT, '  N(profile.path) = ',N_ELEMENTS(profile.path)

            IF (count NE N_ELEMENTS(profile.path)) THEN BEGIN
              PRINT, 'WARNING cortex_GeneratePlots: Spectra count and profile path have different numbers of elements, skipping'
              PRINT, '  count           = ',count
              PRINT, '  N(profile.path) = ',N_ELEMENTS(profile.path)
              CONTINUE
            ENDIF
            ; Second test: see if the path data are the same:
            found = 1
            FOR k = 0, count-1 DO BEGIN
              IF ((ABS(profile.path [k] - spectra.path [j[k]]) GT 1.0E-5) OR  $
                  (ABS(profile.delta[k] - spectra.delta[j[k]]) GT 1.0E-5)) THEN found = 0
;              IF ((ABS(profile.path [k] - spectra.path [j[k]]) GT 5.0E-6) OR  $
;                  (ABS(profile.delta[k] - spectra.delta[j[k]]) GT 5.0E-6)) THEN found = 0
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
        status = 0

        time_s = -1
        time_e = -1
        time_d =  1
        path_local = path
        time_suffix = ''
        IF (plot.time_span[0] NE -1) THEN BEGIN
          time_s = plot.time_span[0]
          time_e = plot.time_span[1]
          time_d = plot.time_span[2]
          dir_structure = 0
        ENDIF
        k = -1
        FOR j = time_s, time_e, time_d DO BEGIN

          FOR i = 0, ncase-1 DO BEGIN
            k++
            IF (j NE -1) THEN BEGIN
               time_suffix = '_' + (STRTRIM(STRING(j,FORMAT='(I8.8)'),2))
               path_local = '/$HOME/divimp/shots/iter/aps_elm/results/' + plot.case_name[i] + '/'
            ENDIF                  
            file_path = path_local + plot.case_name[i] + time_suffix + '.'
print,'file_path ',file_path

            SWITCH plot.option OF
;             ------------------------------------------------------------         
              1: 
              2: 
              4: BEGIN
                target   = cortex_LoadTargetData      (file_path + plot.data_file[0])
                midplane = cortex_LoadMidplaneProfiles(file_path + plot.data_file[1])
                plot_data = { target : target, midplane : midplane }
                name = 'data' + STRING(k+1,FORMAT='(I0)')
                IF (k EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                                 data_array = CREATE_STRUCT(data_array,name,plot_data)
                BREAK
                END
;             ------------------------------------------------------------         
              3: BEGIN
                wall     = cortex_LoadWallProfiles    (file_path + plot.data_file[0])
                midplane = cortex_LoadMidplaneProfiles(file_path + plot.data_file[1])
                plot_data = { wall : wall, midplane : midplane }
                name = 'data' + STRING(k+1,FORMAT='(I0)')
                IF (k EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                                 data_array = CREATE_STRUCT(data_array,name,plot_data)
                BREAK
                END
;             ------------------------------------------------------------         
              ELSE: BEGIN  
                PRINT, 'ERROR cortex_GeneratePlots: Unrecognised target plot option'
                PRINT, '  OPTION = ',option
                status = -1
                BREAK
                END
;             ------------------------------------------------------------         
            ENDSWITCH
          ENDFOR
        ENDFOR

        i = 0
        k = 0
        spike = -1
        FOR j = time_s, time_e, time_d DO BEGIN
          k++
          IF (j EQ plot.spike) THEN spike = k
          IF (plot.ps_local EQ 1 and ps_active_local EQ 0) THEN BEGIN
            IF (ps_active EQ 1) THEN PSCLOSE
            IF (ps EQ 'on') THEN BEGIN
              time_suffix = '_' + (STRTRIM(STRING(j,FORMAT='(I8.8)'),2))
              path_local  = '/$HOME/divimp/shots/iter/aps_elm/results/' + plot.case_name[i] + '/'
              file_path   = path_local + plot.case_name[i] + '_' + plot.id + time_suffix + '.'
              file = file_path + 'idl.ps'
              PRINT, 'POSTSCRIPT FILE=',file
              PSOPEN, filename=file
              ps_active       = 0
              ps_active_local = 1
            ENDIF
          ENDIF

          IF (time_s NE -1) THEN time_slice = j * 10 ELSE time_slice = -1

;          print, 'going for a plot...',k,time_slice
;          help,data_array,/struct

          IF (status eq 0) THEN  $
            status = cortex_PlotTargetProfiles(plot,data_array, ps=ps,slice=k,spike=spike,time_slice=time_slice)

          IF (time_s NE -1) THEN ERASE

          IF (plot.ps_local EQ 1 AND plot.frame EQ 1) THEN BEGIN
            ps_active_local = 0
            ERASE
            PSCLOSE
          ENDIF
      
        ENDFOR

        END
;     ------------------------------------------------------------------
      'PLOT 1D RADIAL PROFILE': BEGIN
        CASE plot.option OF
;         --------------------------------------------------------------         
          1: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file = path + plot.case_name[i] + '.'
              mid = cortex_LoadMidplaneProfiles(file+plot.data_file[0])
              name = 'data' + STRING(i+1,FORMAT='(I0)')
              IF (i EQ 0) THEN mid_array = CREATE_STRUCT(          name,mid) ELSE  $
                               mid_array = CREATE_STRUCT(mid_array,name,mid)
            ENDFOR
;            help, mid_array, /struct
            status = cortex_PlotMidplaneProfiles(plot,mid_array,ps=ps)
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
print, 'status', status
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
;         --------------------------------------------------------------         
;         Time evolution of DIVIMP impurity concentration just inside the separatrix
          2: BEGIN

            time_s = -1
            time_e = -1
            time_d =  1
            path_local = path
            time_suffix = ''
            IF (plot.time_span[0] NE -1) THEN BEGIN
              time_s = plot.time_span[0]
              time_e = plot.time_span[1]
              time_d = plot.time_span[2]
              dir_structure = 0
            ENDIF
            plot_hold = plot
            k = -1
            FOR i = 0, ncase-1 DO BEGIN

              FOR n = 0, N_ELEMENTS(plot.case_set)-1 DO BEGIN
                IF (plot.case_set[n] EQ i+1) THEN BREAK
              ENDFOR
;print,n
              FOR j = time_s, time_e, time_d DO BEGIN
                k++
                IF (time_s NE -1) THEN BEGIN
                  time_suffix = '_' + (STRTRIM(STRING(j,FORMAT='(I8.8)'),2))
                  path_local = '/$HOME/divimp/shots/iter/aps_elm/results/' + plot.case_name[i] + '/'
                  IF (j GT time_s) THEN BEGIN
                    plot.case_set = [plot.case_set[0:n],i+1,plot.case_set[n+1:N_ELEMENTS(plot.case_set)-2]]
;print,plot.case_set[0:10]
                  ENDIF
                ENDIF
                IF (plot.case_name[i] EQ 'none') THEN BEGIN
                  core    = CREATE_STRUCT('version',0.0,'file','none')
                  summary = CREATE_STRUCT('version',0.0,'file','none')
                ENDIF ELSE BEGIN
                  file = path_local + plot.case_name[i] + time_suffix + '.'
                  core = cortex_LoadCoreProfiles(file + plot.data_file[0])
print,'file',file
                  IF (plot.warnings) THEN summary = cortex_LoadDIVIMPSummary(file + plot.data_file[1]) ELSE summary = 0
                ENDELSE
                plot_data = { core:core, summary:summary, time_slice:j }
                name = 'data' + STRING(k+1,FORMAT='(I0)')
                IF (k EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                                 data_array = CREATE_STRUCT(data_array,name,plot_data)
              ENDFOR
            ENDFOR

            status = cortex_PlotSummary(plot, data_array, ps=ps)

            END
;         --------------------------------------------------------------         
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
        IF (plot.no_wall EQ 0) THEN BEGIN 
          wall = cortex_LoadFluidWall_Debug(path + plot.case_name[0] + '.' + plot.data_file[1])
        ENDIF ELSE  $
          wall   = {               $
            n      : 0          ,  $
            class  : 0          ,  $
            group  : 0          ,  $
            index  : 0          ,  $
            tube   : 0          ,  $
            target : 0          ,  $
            v1     : [0.0,0.0]  ,  $
            v2     : [0.0,0.0]     }
        status = cortex_PlotFluidGrid(plot, grid, wall, 0, 0, 'main', 'full', ps=ps)
        END
;     ------------------------------------------------------------------
      'PLOT 2D CONTOUR': BEGIN
        
        time_s = -1
        time_e = -1
        time_d =  1
        file_path = path + plot.case_name[0] + '.'        
        file_suffix = ''
        plot_local = plot
        IF (plot.time_span[0] NE -1) THEN BEGIN
          time_s = plot.time_span[0]
          time_e = plot.time_span[1]
          time_d = plot.time_span[2]
          dir_structure = 0
        ENDIF

        FOR time_i = time_s, time_e, time_d DO BEGIN

          k = -1
          WHILE (1) DO BEGIN
            IF (time_s NE -1) THEN BEGIN
              k++
              IF (k EQ 0) THEN BEGIN
                plot_local = plot
              ENDIF ELSE BEGIN
                plot_local = cortex_ExtractStructure(plot_array,iplot+k)
              ENDELSE
              file_suffix = '_' + STRTRIM(STRING(time_i,FORMAT='(I8.8)'),2) + '.'
              file_path = '/$HOME/divimp/shots/iter/bfm/results/' + plot_local.case_name[0] + '/' + plot_local.case_name[0] + file_suffix
              file_ps   = '/$HOME/divimp/shots/iter/bfm/results/' + plot_local.case_name[0] + '/' + plot_local.case_name[0] + '_' + plot.id + file_suffix
;              file_path = '/$HOME/divimp/shots/iter/aps_elm/results/' + plot_local.case_name[0] + '/' + plot_local.case_name[0] + file_suffix
;              file_ps   = '/$HOME/divimp/shots/iter/aps_elm/results/' + plot_local.case_name[0] + '/' + plot_local.case_name[0] + '_' + plot.id + file_suffix
            ENDIF

;            print, '2d contour: path ',file_path,time_i

            IF (plot_local.option LT 900) THEN BEGIN
              grid = cortex_LoadFluidGrid(file_path + plot_local.data_file[0])
              wall = cortex_LoadWall     (file_path + plot_local.data_file[1]) 
              IF (plot_local.annotate_n NE 0) THEN  $
                annotate = cortex_LoadAnnotationData(plot_local,file_path) ELSE annotate = 0 
              IF (plot_local.option GE   1 AND plot_local.option LE  99) THEN plot_data = cortex_LoadPlasmaData(file_path + plot_local.data_file[2])
              IF (plot_local.option GE 100 AND plot_local.option LE 199) THEN plot_data = cortex_LoadSourceData(file_path + plot_local.data_file[3])
              IF (plot_local.option GE 200 AND plot_local.option LE 299) THEN plot_data = cortex_LoadEIRENEData(file_path + plot_local.data_file[4])        
              IF (plot_local.option GE 300 AND plot_local.option LE 399) THEN plot_data = cortex_LoadEIRENEImpurityData(file_path + plot_local.data_file[5])        
              IF (plot_local.option GE 400 AND plot_local.option LE 400) THEN plot_data = cortex_LoadDIVIMPImpurityData_Density(file_path + plot_local.data_file[6])        
;              IF (plot_local.option GE 500 AND plot_local.option LE 500) THEN plot_data = cortex_LoadEmissionData(file_path + plot_local.data_file[8])        
              IF (plot_local.option GE 500 AND plot_local.option LE 500) THEN plot_data = cortex_LoadImpurityEmission_DIVIMP(file_path + plot_local.data_file[9])        
            ENDIF
            IF (plot_local.option EQ 900) THEN BEGIN
              ; WALLDYN redistribution matix plot - semi-hack:
              plot_data = cortex_LoadWalldyn(plot_local)
              n = plot_data.redep_matrix.data1.n
              x = FINDGEN(n+1) + 1.0  ;  FINDGEN(n+1) / FLOAT(n) 
              grid     = {obj_n : n*n, n : n, x : x }
              wall     = -1.0
              annotate = -1.0
            ENDIF
            IF (plot_local.option LE   0 OR  plot_local.option GE 901) THEN BEGIN
              PRINT, 'ERROR cortex_GeneratePlots: Unrecognised 2D contour plot option'
              PRINT, '  PLOT TAG = ',plot_local.tag
              cortex_Exit, nargs
            ENDIF

            IF (plot.ps_local EQ 1 and ps_active_local EQ 0) THEN BEGIN
              IF (ps_active EQ 1) THEN BEGIN
                ERASE
                PSCLOSE
              ENDIF
              IF (ps EQ 'on') THEN BEGIN
                file = file_ps + 'idl.ps'
                PRINT, 'POSTSCRIPT FILE=',file
                PSOPEN, filename=file
                ps_active       = 0
                ps_active_local = 1
              ENDIF
            ENDIF
          
            status = cortex_Plot2DContour(plot_local, plot_data, grid, wall, annotate, ps=ps,time_slice=time_i)
          
            IF (plot.ps_local EQ 1 AND plot_local.frame EQ 1) THEN BEGIN
              ps_active_local = 0
              ERASE
              PSCLOSE
              BREAK
            ENDIF

            IF (time_s EQ -1) THEN BREAK

          ENDWHILE

        ENDFOR

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
      'PLOT 3D TRIANGLES': BEGIN
        file_path = path + plot.case_name[0] + '.'
        IF (plot.quick EQ 1) THEN BEGIN
          RESTORE, FILENAME='output/'+STRTRIM(plot.case_name[0],2)+'_1.sav'
        ENDIF ELSE BEGIN
          data1 = cortex_LoadTriangleGeometry(file_path + plot.data_file[0])
          data2 = cortex_LoadTriangleData    (file_path + plot.data_file[1])
          data3 = cortex_LoadTetrahedronData (file_path + plot.data_file[2])
          SAVE, FILENAME='output/'+STRTRIM(plot.case_name[0],2)+'_1.sav', data1, data2, data3
        ENDELSE


print,4
print,file_path + plot.data_file[3]
        data4 = cortex_LoadAnnotations(4,file_path + plot.data_file[3])

;        plot_data = { data1 : data1 }
;        plot_data = { data1 : data1, data2 : data2 }
;        plot_data = { data1 : data1, data2 : data2, data3 : data3 }
        plot_data = { data1 : data1, data2 : data2, data3 : data3, data4 : data4 }

        IF (plot.annotate_n NE 0) THEN  $
          annotate = cortex_LoadAnnotationData(plot,file_path) ELSE annotate = 0 

; help,plot_data,/struct

;        PSCLOSE
;        SET_PLOT, 'Z' 
;        DEVICE, Set_Resolution=[800,800],DECOMPOSED=0, SET_PIXEL_DEPTH=24

        status = cortex_Plot3DTriangles(plot, plot_data, annotate=annotate, ps=ps)

;        buffer = TVRD(TRUE=1)
;        SET_PLOT, 'X'
;        WINDOW, 0, XSIZE=1000, YSIZE=1000
;        TV, buffer, TRUE=1
;        PSOPEN, filename=file
;        TV, buffer, TRUE=1

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

    IF (status NE 0) THEN BEGIN
      IF (ps_active EQ 1) THEN PSCLOSE
      cortex_Exit, nargs
    ENDIF

  ENDFOR

  IF (ps_active EQ 1) THEN PSCLOSE

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
