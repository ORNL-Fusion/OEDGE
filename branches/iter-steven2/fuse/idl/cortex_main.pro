; example : cortex_main,['t-new-0000a','i-osm.ctx']
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
  
  option = 3

  ps = 'on'

  nargs = N_ELEMENTS(args)

  case_name  = args[0]
  input_file = args[1]

  IF (nargs EQ 4) THEN BEGIN
    input_file = args[2] + '/' + input_file
    data_path  = args[3] + '/'
  ENDIF ELSE BEGIN
    input_file = '/home/ITER/lisgos/divimp/data/' + input_file
    data_path  = '/home/ITER/lisgos/divimp/results/'
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
    file = data_path + case_name + '.idl.ps'
;      PRINT, 'POSTSCRIPT FILE=',file
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
    PRINT,'---- PLOT STRUCTURE:',iplot,plot.tag
;    HELP,plot,/struct    

    IF (plot.option EQ 0) THEN CONTINUE

;    PRINT,'NCASES=',N_ELEMENTS(plot.case_name),N_ELEMENTS(WHERE(plot.case_name NE 'unknown'))

    ncase = N_ELEMENTS(WHERE(plot.case_name NE 'unknown'))

    status = 0
;
;   For more than one plot on a page, check if the spacing has been adjusted by the
;   previous plot:
;
    IF (iplot GT 1) THEN BEGIN
      plot_last = cortex_ExtractStructure(plot_array,iplot-1)
      IF (plot_last.frame EQ 0 AND plot_last.frame_bnds[1] NE plot.frame_bnds[0]) THEN  $
        plot.frame_bnds[0] = plot_last.frame_bnds[1] 
    ENDIF
;
;
;
    CASE plot.tag OF
;     ------------------------------------------------------------------
      'PLOT 1D LOS INTEGRAL': BEGIN
        CASE plot.option OF
;         --------------------------------------------------------------         
          1: BEGIN
            FOR i = 0, ncase-1 DO BEGIN
              file_path = path + plot.case_name[i] + '.'
              j = WHERE(plot.data_file NE 'unknown',count)
              FOR j = 0, count-1 DO BEGIN
                integral = cortex_LoadIntegrals(file_path + plot.data_file[0])
                name = 'data' + STRING(j+1,FORMAT='(I0)')
                IF (j EQ 0) THEN integral_array = CREATE_STRUCT(               name,integral) ELSE  $
                                 integral_array = CREATE_STRUCT(integral_array,name,integral)
              ENDFOR
              plot_data = { integral : integral_array, dummy : 1.0 }  ; Dummy is there to calm down the error check in ExtractStructure, which is lame, and needs fixing...
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
      'PLOT 1D EIRENE ENERGY SPECTRUM': BEGIN
        FOR i = 0, ncase-1 DO BEGIN
          file = path + plot.case_name[i] + '.idl.eirene_spectrum_001_sum'
          plot_data = cortex_LoadEnergySpectrum(file)
          name = 'data' + STRING(i+1,FORMAT='(I0)')
; IMPROVE THIS BY INITIALISING data_array AT THE TOP OF THE LOOP AND STANDARDIZING THE PACKAGING
; OF PLOT DATA FOR ALL PLOTS, ie THEY ALL USE data_array
          IF (i EQ 0) THEN data_array = CREATE_STRUCT(           name,plot_data) ELSE  $
                           data_array = CREATE_STRUCT(data_array,name,plot_data)
        ENDFOR
        status = cortex_PlotEnergySpectrum(plot,data_array, ps=ps)
        END
;     ------------------------------------------------------------------
      'PLOT 1D TARGET PROFILE': BEGIN
        FOR i = 0, ncase-1 DO BEGIN
          tar = cortex_LoadTargetData(path + plot.case_name[i] + '.' + plot.data_file[0])
          name = 'data' + STRING(i+1,FORMAT='(I0)')
          IF (i EQ 0) THEN tar_array = CREATE_STRUCT(          name,tar) ELSE  $
                           tar_array = CREATE_STRUCT(tar_array,name,tar)
        ENDFOR
        status = cortex_PlotTargetProfiles(tar_array, ps=ps)
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
              target = cortex_LoadTargetData(path + plot.case_name[i] + '.' + plot.data_file[3])

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
          mid = cortex_LoadPedestalModel(path + plot.case_name[i] + '.' + plot.data_file[0])
          name = 'data' + STRING(i+1,FORMAT='(I0)')
          IF (i EQ 0) THEN mid_array = CREATE_STRUCT(          name,mid) ELSE  $
                           mid_array = CREATE_STRUCT(mid_array,name,mid)
        ENDFOR
        status = cortex_PlotPedestalModel(plot, mid_array, ps=ps)
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
                plot_data = { core:core, summary:summary }
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
        file = path + plot.case_name[0] + '.'
        IF (plot.equ NE 'default') THEN grid = grid_ReadEQUFile    (plot.equ                ) ELSE  $
                                        grid = cortex_LoadFluidGrid(file + plot.data_file[0])
        wall = cortex_LoadWall(file + plot.data_file[1])
        IF (plot.nodes) THEN node = cortex_LoadNodeData(file + plot.data_file[2]) ELSE node = 0
        IF (plot.annotate_n NE 0) THEN BEGIN
          annotate = -1
          FOR i = 0, plot.annotate_n-1 DO BEGIN
            CASE plot.annotate_code[i] OF
              1: annotate_data = cortex_LoadAnnotationData(plot.annotate_code[i],plot.annotate_file[i]) 
              ELSE: BEGIN
                PRINT, 'ERROR cortex_GeneratePlots: Unknown annotation code'
                PRINT, '  PLOT TAG = ',plot.tag
                PRINT, '  CODE     = ',plot.annotate_code[i]
                annotate_data = -1
                END
            ENDCASE
            name = 'data' + STRING(i+1,FORMAT='(I0)')
            IF (annotate EQ -1) THEN annotate = CREATE_STRUCT(         name,annotate_data)  ELSE  $
                                     annotate = CREATE_STRUCT(annotate,name,annotate_data)
          ENDFOR 
        ENDIF ELSE annotate = 0 
        status = cortex_PlotFluidGrid(plot, grid, wall, node, annotate, 'main', 'full', ps=ps)
        END
;     ------------------------------------------------------------------
      'PLOT 2D FLUID GRID - DEBUG': BEGIN
        grid = cortex_LoadFluidGrid_Debug(path + plot.case_name[0] + '.' + plot.data_file[0])
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
        file_path = path + plot.case_name[0] + '.'
        grid = cortex_LoadFluidGrid (file_path + plot.data_file[0])
        wall = cortex_LoadWall      (file_path + plot.data_file[1])
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
        status = cortex_Plot2DContour(plot, plot_data, grid, wall, ps=ps)
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