;
; ======================================================================
; 
; Based on inReadLine in interface.pro, so if improved here, improve there as well...
;
FUNCTION cortex_ReadLine, fp, buffer
  status = -1
  WHILE NOT EOF(fp) DO BEGIN
;  WHILE (1) DO BEGIN
    READF,fp,buffer
    i = STRPOS(buffer,'$')
    IF (i EQ  0) THEN CONTINUE 
    IF (i NE -1) THEN BEGIN
      str    = STRSPLIT(buffer,'$',/EXTRACT)
      buffer = str[0]
    ENDIF
    i = STRPOS(buffer,'*')               ; *** Combine with above check somehow... ***
    IF (i EQ  0) THEN CONTINUE 
    IF (i NE -1) THEN BEGIN
      str    = STRSPLIT(buffer,'*',/EXTRACT)
      buffer = str[0]
    ENDIF
    i = STRPOS(buffer,';') 
    IF (i EQ  0) THEN CONTINUE 
    IF (i NE -1) THEN BEGIN
      str    = STRSPLIT(buffer,';',/EXTRACT)
      buffer = str[0]
    ENDIF
;    i = STRPOS(buffer,'!') 
;    IF (i EQ  0) THEN CONTINUE 
;    IF (i NE -1) THEN BEGIN
;      str    = STRSPLIT(buffer,'!',/EXTRACT)
;      buffer = str[0]
;    ENDIF
    IF (STRLEN(STRTRIM(buffer)) EQ 0) THEN CONTINUE
    IF (STRMATCH(buffer,'*{EXIT}*') EQ 1) THEN RETURN, 0
    IF (STRMATCH(buffer, '*{*'    ) EQ 1) THEN RETURN, 1
    RETURN, 2
  ENDWHILE
  RETURN, status
END 
;
; ======================================================================
;
FUNCTION cortex_ProcessPlotStruct,plot_struct,plot_array,default,n

  n = n + 1

  IF (n GT 1) THEN BEGIN
    plot_name = 'data' + STRING(n-1,FORMAT='(I0)')   
    IF (n EQ 2) THEN BEGIN
      plot_array = CREATE_STRUCT(plot_name,plot_struct)
    ENDIF ELSE BEGIN
      plot_array = CREATE_STRUCT(plot_array,plot_name,plot_struct)
    ENDELSE
  ENDIF

  long_array   = MAKE_ARRAY(100,/LONG  ,VALUE=0        )
  string_array = MAKE_ARRAY(100,/STRING,VALUE='unknown')
  float_array  = MAKE_ARRAY(100,/FLOAT ,VALUE='-999.0' )

  line_seg = MAKE_ARRAY(4*5,/FLOAT,VALUE=0.0)

  struct = {                               $
    verison         : 1.0               ,  $
    tag             : 'unknown'         ,  $
    option          : -1                ,  $
    id              : 'unknown'         ,  $
    title           : 'unknown'         ,  $
    notes           : 'unknown'         ,  $
    plot_title      : 'unknown'         ,  $
    xtitle          : 'default'         ,  $
    ytitle          : 'default'         ,  $
    xlabels         : string_array      ,  $
    charsize        : 1.0               ,  $
    thick           : 1.0               ,  $
    xstyle          : 1                 ,  $
    ystyle          : 1                 ,  $
    color_table     : 5                 ,  $
    default         : 'default'         ,  $
    case_name       : string_array      ,  $
    case_set        : [0,long_array]    ,  $
    case_set_name   : string_array      ,  $
    data_path       : default.data_path ,  $
    data_file       : string_array      ,  $
    tubes           : long_array        ,  $
    state           : 0                 ,  $
    nodes           : 0                 ,  $
    annotate_n      : 0                 ,  $  ; Number of annotations requested for 2D fluid grid plot
    annotate_code   : long_array        ,  $  ; Type of annotation
    annotate_colour : string_array      ,  $  ; Colour
    annotate_file   : string_array      ,  $  ; Location of data describing the annotation
    annotate_step   : [1,1,1,1,1,1,1,1] ,  $  ; Index step when plotting
    annotate_label  : long_array        ,  $  ; Label annotations
    equ             : 'default'         ,  $
    equ_params      : [20,1.0,1.0]      ,  $
    equ_levels      : float_array       ,  $
    frame           : 1                 ,  $
    frame_bnds      : [0.0,1.0]         ,  $
    outline         : 0                 ,  $
    flux_surfaces   : 0                 ,  $  ; PlotFluidGrid - show tube radial boundaries
    show_grid       : 0                 ,  $
    no_grid         : 0                 ,  $
    no_wall         : 0                 ,  $
    no_separatrix   : 0                 ,  $
    area            : [0.0,0.0,0.0,0.0] ,  $
    center          : [0.0,0.0]         ,  $
    size            : 0.0               ,  $
    zoom            : [0.0,0.0,0.0,0.0] ,  $
    xticks          : -1                ,  $
    yticks          : -1                ,  $
    xminor          : -1                ,  $
    yminor          : -1                ,  $
    xrange          : [0.0,0.0]         ,  $
    yrange          : [0.0,0.0]         ,  $
    ylimit          : [0.0,0.0,0.0,0.0] ,  $
    xmark           : 'none'            ,  $
    peak            : 0                 ,  $
    sum             : 0                 ,  $
    scale_factor    : 1.0               ,  $
    log             : 0                 ,  $
    xdata           : 'psi_n'           ,  $
    xlog            : 0                 ,  $
    ylog            : -1                ,  $
    smooth          : 0                 ,  $
    position        : [0.0,0.0,0.0,0.0] ,  $
    focus           :  0                ,  $
    warnings        :  0                ,  $
    aspect_ratio    : -1.0              ,  $
    line_seg_n      : 0                 ,  $
    line_seg_s      : [0,0,0,0,0]       ,  $
    line_seg        : line_seg          ,  $
    label_index     : 'none'            ,  $
    show_n          : 0                 ,  $
    show_tube       : long_array        ,  $
    show_colour     : string_array      ,  $
    dummy           : -1         }

;  help,struct,/struct

  RETURN, struct
END

;
; ======================================================================
;
FUNCTION cortex_ApplySubstitutions, buffer

  debug = 0

  str = STRSPLIT(STRTRIM(buffer,2),'{',/EXTRACT)

  ilist = 0
  FOR i = 0, N_ELEMENTS(str)-1 DO BEGIN
    IF (STRUPCASE(STRMID(str[i],0,6)) EQ 'STRING') THEN BEGIN
      i1 = STRPOS(str[i],'}') + 1
      i2 = STRLEN(str[i])
      str1 = STRTRIM(STRMID(str[i],i1,i2),2)
      i1 = STRPOS(str1,' ')
      str2 = STRTRIM(STRMID(str1,0   ,i1          ),2)
      str3 = STRTRIM(STRMID(str1,i1+1,STRLEN(str1)),2)
      IF (ilist EQ 0) THEN BEGIN
        list_tag = str2
        list_str = str3
      ENDIF ELSE BEGIN
        list_tag = [list_tag,str2]
        list_str = [list_str,str3]
      ENDELSE
      ilist = ilist + 1
      str[i] = 'empty'
    ENDIF
  ENDFOR

  IF (debug) THEN FOR i = 0, ilist-1 DO print,list_tag[i]
  IF (debug) THEN FOR i = 0, ilist-1 DO print,list_str[i]

  FOR i = 0, ilist-1 DO BEGIN
    FOR j = 0, N_ELEMENTS(str)-1 DO BEGIN 

      status = 1
      WHILE (status) DO BEGIN
        status = 0
        str1 = str[j]
        str2 = list_tag[i]
        str3 = list_str[i]

        IF (debug) THEN BEGIN
          print,str1
          print,str2
          print,str3
        ENDIF

        len = STRLEN(str2)
        FOR k = 0, STRLEN(str1)-len DO BEGIN 
          IF (debug) THEN PRINT,'>'+str2+'<>'+STRMID(str1,k,len)+'<'
          IF STRCMP(str2,STRMID(str1,k,len),/FOLD_CASE) THEN BEGIN
            str_sub = STRSPLIT(STRTRIM(STRMID(str1,k),2),' ',/EXTRACT)
            IF (N_ELEMENTS(str_sub) GE 4) THEN BEGIN
              IF (str_sub[1] EQ 'sub') THEN BEGIN
                IF (debug) THEN print,str3
                str_tmp = str3
                FOR l = 0, 99 DO str_tmp = STR_REPLACE(str_tmp,str_sub[2],str_sub[3])
                str3 = str_tmp
                IF (debug) THEN print,str3
              ENDIF
              IF (debug) THEN print,str1
              str_tmp = STRMID(str1,0,k-1) + ' ' + str_sub[0]
              FOR l = 4, N_ELEMENTS(str_sub)-1 DO str_tmp = str_tmp + ' ' + str_sub[l] 
              str1 = str_tmp
            ENDIF
            IF (debug) THEN print,'working'
            IF (debug) THEN  print,str1
            str1 = STR_REPLACE(str1,str2,str3)
            IF (debug) THEN print,str1
            IF (debug) THEN print, 'done'
            str[j] = str1
            status = 1
            BREAK
          ENDIF
        ENDFOR
      ENDWHILE

    ENDFOR
  ENDFOR

; Rebuild buffer:
  buffer = ' '
  FOR i = 0, N_ELEMENTS(str)-1 DO IF (str[i] NE 'empty') THEN buffer = buffer + ' {' + str[i]

;print,buffer
;stop

  RETURN, STRTRIM(buffer,2)

END
;
; ======================================================================
;
FUNCTION cortex_LoadPlotData,case_name,input_file,result

  file_name = input_file

  default = { data_path : 'default' }

  fp = 1
  FREE_LUN,fp
  OPENR, fp, file_name, ERROR=err
  IF (err NE 0) THEN BEGIN
    PRINT,'ERROR cortex_LoadPlotData: Unable to open input file'
    PRINT,' FILE_NAME= ',file_name
    PRINTF,-2,!error_state.msg
    FREE_LUN,fp  
    RETURN, -1
  ENDIF

  buffer = ' '
  monster_buffer = buffer

; Load in the entire input file:

  WHILE (1) DO BEGIN
    status = cortex_ReadLine(fp,buffer)    
    IF (status EQ -1) THEN BREAK
    monster_buffer = monster_buffer + STRTRIM(buffer,2) + ' '
  ENDWHILE
  buffer = monster_buffer + ' {EXIT} '

;print,monster_buffer

  buffer = cortex_ApplySubstitutions(buffer)

  n = 0
  nchop = 0
  ncase = 0
  plot_struct = -1
  plot_array  = -1

  WHILE (1) DO BEGIN

    IF (nchop EQ 0) THEN BEGIN
;      status = cortex_ReadLine(fp,buffer)
      buffer_array = STRSPLIT(STRTRIM(buffer,2),'{',/EXTRACT)
      nchop = N_ELEMENTS(buffer_array)
      IF (nchop GT 1) THEN BEGIN
        ichop = 0
        nchop = nchop - 1
        buffer = '{' + buffer_array[ichop]
      ENDIF ELSE BEGIN
        nchop = 0
      ENDELSE
    ENDIF ELSE BEGIN
      ichop = ichop + 1
      nchop = nchop - 1
      buffer = '{' + buffer_array[ichop]
    ENDELSE

    IF (status EQ 2) THEN BEGIN
      PRINT,'ERROR cortex_LoadPlotData: Tag expected but not found'
      PRINT,'  FILE_NAME= ',file_name
      PRINT,'  BUFFER   = ',STRTRIM(buffer)
      RETURN,-1
    ENDIF

;   Extract tag and isolate data:
    i1 = STRPOS(buffer,'{')
    i2 = STRPOS(buffer,'}')
    tag  = STRUPCASE(STRMID(buffer,i1+1,i2-i1-1))
    data = STRTRIM(STRMID(buffer,i2+1),2)
    data_array = STRSPLIT(data,/EXTRACT)

;    PRINT,'TAG       :',tag
;    PRINT,'DATA      :',data
;    PRINT,'DATA_ARRAY:',data_array[0],data_array[1]

    CASE tag OF

; LEFT OFF
; need to set new 2D integral plot that show the image... and does some analysis...
;     ------------------------------------------------------------------
      'PLOT 2D LOS INTEGRAL': BEGIN
        ncase = 1
        nset  = 0
        plot_struct = cortex_ProcessPlotStruct(plot_struct,plot_array,default,n)
        plot_struct.tag          = tag
        plot_struct.option       = FIX(data)
        plot_struct.title        = '2D LOS INTEGRALS'
        plot_struct.default      = case_name
        plot_struct.case_name[0] = case_name
        CASE plot_struct.option OF
          0:
          1: 
          ELSE: BEGIN
            PRINT,'ERROR cortex_LoadPlotData: Unknown 2D LOS integral plot option'
            PRINT,'  FILE_NAME= ',file_name
            PRINT,'  TAG      = ',tag
            PRINT,'  OPTION   = ',plot_struct.option
            RETURN,-1
            END
          ENDCASE
        END
;     ------------------------------------------------------------------
      'PLOT 1D LOS INTEGRAL': BEGIN
        ncase = 1
        nset  = 0
        plot_struct = cortex_ProcessPlotStruct(plot_struct,plot_array,default,n)
        plot_struct.tag          = tag
        plot_struct.option       = FIX(data)
        plot_struct.title        = 'LOS INTEGRALS'
        plot_struct.default      = case_name
        plot_struct.case_name[0] = case_name
        CASE plot_struct.option OF
          0:
          1: 
          2: plot_struct.data_file[0] = 'idl.ray_profile_' 
          ELSE: BEGIN
            PRINT,'ERROR cortex_LoadPlotData: Unknown 1D LOS integral plot option'
            PRINT,'  FILE_NAME= ',file_name
            PRINT,'  TAG      = ',tag
            PRINT,'  OPTION   = ',plot_struct.option
            RETURN,-1
            END
          ENDCASE
        END
;     ------------------------------------------------------------------
      'PLOT 1D WALL PROFILE': BEGIN
        ncase = 1
        nset  = 0
        plot_struct = cortex_ProcessPlotStruct(plot_struct,plot_array,default,n)
        plot_struct.tag          = tag
        plot_struct.option       = FIX(data)
        plot_struct.default      = case_name
        plot_struct.case_name[0] = case_name
        CASE plot_struct.option OF
          0:
          1: plot_struct.data_file[0] = 'idl.divimp_flux_wall'
          2: plot_struct.data_file[0] = 'idl.divimp_flux_wall'
          3: BEGIN
            plot_struct.data_file[0] = 'idl.eirene_flux_wall'   ; Wall fluxes from EIRENE
            plot_struct.data_file[1] = 'idl.divimp_summary'
            plot_struct.data_file[2] = 'idl.core_impurities'
            plot_struct.data_file[3] = 'idl.fluid_grid'
            plot_struct.data_file[4] = 'idl.fluid_wall'
            END
          4: plot_struct.data_file[0] = 'idl.divimp_flux_wall'
          ELSE: BEGIN
            PRINT,'ERROR cortex_LoadPlotData: Unknown 1D wall profile plot option'
            PRINT,'  FILE_NAME= ',file_name
            PRINT,'  TAG      = ',tag
            PRINT,'  OPTION   = ',plot_struct.option
            RETURN,-1
            END
          ENDCASE
        END
;     ------------------------------------------------------------------
      'PLOT 1D PARALLEL PROFILE': BEGIN
        ncase = 1
        nset  = 0
        plot_struct = cortex_ProcessPlotStruct(plot_struct,plot_array,default,n)
        plot_struct.tag          = tag
        plot_struct.option       = FIX(data)
        plot_struct.title        = 'PROFILES PARALLEL TO THE MAGNETIC FIELD (FLUX-TUBES)'
        plot_struct.default      = case_name
        plot_struct.case_name[0] = case_name
        CASE plot_struct.option OF
          0:
          1: BEGIN
             plot_struct.data_file[0] = 'idl.fluid_plasma'
             plot_struct.data_file[1] = 'idl.fluid_sources'
             plot_struct.data_file[2] = 'idl.fluid_eirene'
             plot_struct.data_file[3] = 'idl.fluid_targets'
             plot_struct.data_file[4] = 'idl.osm_nodes'
             END
          2: BEGIN
             plot_struct.data_file[0] = 'idl.fluid_plasma'
             plot_struct.data_file[1] = 'idl.fluid_targets'
             plot_struct.data_file[2] = 'idl.osm_nodes'
             END
          100: BEGIN
             plot_struct.data_file[0] = 'idl.eirene_imp'
             plot_struct.data_file[1] = 'idl.divimp_imp_density'
             END
          101: BEGIN
             plot_struct.data_file[0] = 'idl.eirene_imp'
             plot_struct.data_file[1] = 'idl.divimp_imp_ionisation'
             END
          ELSE: BEGIN
            PRINT,'ERROR cortex_LoadPlotData: Unknown 1D parallel profile plot option'
            PRINT,'  FILE_NAME= ',file_name
            PRINT,'  TAG      = ',tag
            PRINT,'  OPTION   = ',plot_struct.option
            RETURN,-1
            END
          ENDCASE
        END
;     ------------------------------------------------------------------
      'PLOT 1D EIRENE ENERGY SPECTRUM': BEGIN
        ncase = 1
        nset  = 0
        plot_struct = cortex_ProcessPlotStruct(plot_struct,plot_array,default,n)
        plot_struct.tag          = tag
        plot_struct.option       = FIX(data)
        plot_struct.title        = 'EIRENE ENERGY SPECTRUM'
        plot_struct.default      = case_name
        plot_struct.case_name[0] = case_name
        END
;     ------------------------------------------------------------------
      'PLOT 1D TARGET PROFILE': BEGIN
        ncase = 1
        nset  = 0
        plot_struct = cortex_ProcessPlotStruct(plot_struct,plot_array,default,n)
        plot_struct.tag          = tag
        plot_struct.option       = FIX(data)
        plot_struct.title        = 'TARGET PLOT'
        plot_struct.default      = case_name
        plot_struct.case_name[0] = case_name
        plot_struct.data_file[0] = 'idl.fluid_targets'
        plot_struct.data_file[1] = 'idl.midplane'
        END
;     ------------------------------------------------------------------
      'PLOT 1D RADIAL PROFILE': BEGIN
        ncase = 1
        nset  = 0
        plot_struct = cortex_ProcessPlotStruct(plot_struct,plot_array,default,n)
        plot_struct.tag     = tag
        plot_struct.option  = FIX(data)
        plot_struct.default = case_name
        CASE plot_struct.option OF
          1: BEGIN
            plot_struct.title        = 'RADIAL MIDPLANE PLOT'
            plot_struct.case_name[0] = case_name
            plot_struct.data_file[0] = 'idl.midplane'
            END
          2: BEGIN
            plot_struct.title        = 'GENERAL RADIAL PLOT'
            plot_struct.case_name[0] = case_name
            plot_struct.data_file[0] = 'idl.fluid_grid'
            plot_struct.data_file[1] = 'idl.fluid_wall'
            plot_struct.data_file[2] = 'idl.fluid_plasma'
            plot_struct.data_file[3] = 'idl.fluid_targets'
            END
          3: BEGIN
            plot_struct.title        = 'RADIAL CORE PLOT'
            plot_struct.case_name[0] = case_name
            plot_struct.data_file[0] = 'idl.core_impurities'
            END
          4: BEGIN
            plot_struct.title        = 'NEW (!) GENERAL RADIAL PLOT'
            plot_struct.case_name[0] = case_name
            plot_struct.data_file[0] = 'idl.fluid_grid'
            plot_struct.data_file[1] = 'idl.fluid_wall'
            plot_struct.data_file[2] = 'idl.fluid_plasma'
            plot_struct.data_file[3] = 'idl.fluid_targets'
            END
          5: BEGIN
            plot_struct.title        = 'NEW (!) GENERAL RADIAL PLOT: Average plasma quantities for debugging'
            plot_struct.case_name[0] = case_name
            plot_struct.data_file[1] = 'idl.fluid_plasma'
            END
          ELSE: BEGIN
            PRINT,'ERROR cortex_LoadPlotData: Unknown 1D radial plot option'
            PRINT,'  FILE_NAME= ',file_name
            PRINT,'  TAG      = ',tag
            PRINT,'  OPTION   = ',plot_struct.option
            RETURN,-1
            END
        ENDCASE
        END
;     ------------------------------------------------------------------
      'PLOT 1D SUMMARY': BEGIN
        ncase = 0
        nset  = 0
        plot_struct = cortex_ProcessPlotStruct(plot_struct,plot_array,default,n)
        plot_struct.tag     = tag
        plot_struct.option  = FIX(data)
        plot_struct.default = case_name
        CASE plot_struct.option OF
          0:
          1: BEGIN
            plot_struct.title        = 'CORE DATA SUMMARY PLOT'
            plot_struct.data_file[0] = 'idl.core_impurities'
            plot_struct.data_file[1] = 'idl.divimp_summary'
            END
          ELSE: BEGIN
            PRINT,'ERROR cortex_LoadPlotData: Unknown summary plot option'
            PRINT,'  FILE_NAME= ',file_name
            PRINT,'  TAG      = ',tag
            PRINT,'  OPTION   = ',plot_struct.option
            RETURN,-1
            END
        ENDCASE
        END
;     ------------------------------------------------------------------
      'PLOT 1D PEDESTAL MODEL': BEGIN
        ncase = 1
        nset  = 0
        plot_struct = cortex_ProcessPlotStruct(plot_struct,plot_array,default,n)
        plot_struct.tag          = tag
        plot_struct.option       = 1
        plot_struct.title        = 'PEDESTAL MODEL'
        plot_struct.default      = case_name
        plot_struct.case_name[0] = case_name
        plot_struct.data_file[0] = 'idl.pedestal'
        END
;     ------------------------------------------------------------------
      'PLOT 2D FLUID GRID': BEGIN
        ncase = 1
        nset  = 0
        plot_struct = cortex_ProcessPlotStruct(plot_struct,plot_array,default,n)
        plot_struct.tag          = tag
        plot_struct.option       = 1
        plot_struct.title        = '2D GRID PLOT'
        plot_struct.default      = case_name
        plot_struct.case_name[0] = case_name
        plot_struct.data_file[0] = 'idl.fluid_grid'
        plot_struct.data_file[1] = 'idl.fluid_wall'
        plot_struct.data_file[2] = 'idl.osm_nodes'
        END
      'PLOT 2D FLUID GRID - DEBUG': BEGIN
        ncase = 1
        nset  = 0
        plot_struct = cortex_ProcessPlotStruct(plot_struct,plot_array,default,n)
        plot_struct.tag          = tag
        plot_struct.option       = 1
        plot_struct.title        = '2D GRID PLOT - DEBUG'
        plot_struct.default      = case_name
        plot_struct.case_name[0] = case_name
        plot_struct.data_file[0] = 'idl.fluid_grid_debug'
        END
;     ------------------------------------------------------------------
      'PLOT 2D CONTOUR': BEGIN
        ncase = 1
        nset  = 0
        plot_struct = cortex_ProcessPlotStruct(plot_struct,plot_array,default,n)
        plot_struct.tag          = tag
        plot_struct.option       = data
        plot_struct.title        = '2D CONTOUR / SHADED PLOT'
        plot_struct.default      = case_name
        plot_struct.case_name[0] = case_name
        plot_struct.data_file[0] = 'idl.fluid_grid'
        plot_struct.data_file[1] = 'idl.fluid_wall'
        plot_struct.data_file[2] = 'idl.fluid_plasma'
        plot_struct.data_file[3] = 'idl.fluid_sources'
        plot_struct.data_file[4] = 'idl.fluid_eirene'
        plot_struct.data_file[5] = 'idl.eirene_imp'
        plot_struct.data_file[6] = 'idl.divimp_imp_density'
        plot_struct.data_file[7] = 'idl.divimp_imp_ionisation'
        END
;     ------------------------------------------------------------------
      'PLOT 3D TEST': BEGIN
        ncase = 1
        nset  = 0
        plot_struct = cortex_ProcessPlotStruct(plot_struct,plot_array,default,n)
        plot_struct.tag          = tag
        plot_struct.option       = data
        plot_struct.title        = '3D'
        plot_struct.default      = case_name
        plot_struct.case_name[0] = case_name
        plot_struct.data_file[0] = 'idl.tet_centroid'
        END
;     ------------------------------------------------------------------
      'ID'           : plot_struct.id           = STRTRIM(data,2)
      'TITLE'        : plot_struct.title        = STRTRIM(data,2)
      'PLOT TITLE'   : plot_struct.plot_title   = STRTRIM(data,2)
      'NOTES'        : plot_struct.notes        = STRTRIM(data,2)
      'DATA PATH'    : default.data_path        = data
      'DATA FILE'    : BEGIN
        FOR j = 0, N_ELEMENTS(data_array)-1 DO BEGIN
          i = WHERE(plot_struct.data_file NE 'unknown',count)  ; Not sure I like this scheme...
          plot_struct.data_file[count] = data_array[j]
        ENDFOR
        END
      'NO FRAME'     : plot_struct.frame        = 0
      'NO ERASE'     : plot_struct.frame        = 0
      'OUTLINE'      : plot_struct.outline      = 1
      'FLUX SURFACES': plot_struct.flux_surfaces= 1
      'SHOW GRID'    : plot_struct.show_grid    = 1
      'NO GRID'      : plot_struct.no_grid      = 1
      'NO WALL'      : plot_struct.no_wall      = 1
      'NO SEPARATRIX': plot_struct.no_separatrix= 1
      'EQU'          : plot_struct.equ          = STRTRIM(data,2)
      'EQU PARAMS'   : plot_struct.equ_params   = FLOAT(data_array)
      'EQU LEVELS'   : plot_struct.equ_levels   = FLOAT(data_array)
      'TUBES'        : plot_struct.tubes        = LONG(data_array)
      'STATE'        : plot_struct.state        = LONG(data)
      'NODES'        : plot_struct.nodes        = 1
      'XLABELS'      : plot_struct.xlabels      = STRSPLIT(data,':',/EXTRACT)
      'CHARSIZE'     : plot_struct.charsize     = FLOAT(data)
      'THICK'        : plot_struct.thick        = FLOAT(data)
      'XSTYLE'       : plot_struct.xstyle       = FIX(data)
      'YSTYLE'       : plot_struct.ystyle       = FIX(data)
      'COLOR TABLE'  : plot_struct.color_table  = FIX(data)
      'FOCUS'        : plot_struct.focus        = LONG(STRTRIM(data,2))
      'WARNINGS'     : plot_struct.warnings     = 1
      'ZOOM'         : plot_struct.zoom         = FLOAT(data_array)
      'XTICKS'       : plot_struct.xticks       = FIX(data)
      'YTICKS'       : plot_struct.yticks       = FIX(data)
      'XMINOR'       : plot_struct.xminor       = FIX(data)
      'YMINOR'       : plot_struct.yminor       = FIX(data)
      'XRANGE'       : plot_struct.xrange       = FLOAT(data_array)
      'YRANGE'       : plot_struct.yrange       = FLOAT(data_array)
      'YLIMIT'       : plot_struct.ylimit       = FLOAT(data_array)
      'XMARK'        : plot_struct.xmark        = data
      'PEAK'         :  $
        IF (data_array[0] EQ 'only') THEN plot_struct.peak = 2 ELSE  $
                                          plot_struct.peak = 1
      'SUM'          :  $
        IF (data_array[0] EQ 'only') THEN plot_struct.peak = 2 ELSE  $
                                          plot_struct.peak = 1
      'SCALE FACTOR' : plot_struct.scale_factor = FLOAT(data)
      'LOG'          : plot_struct.log          = 1
      'XDATA'        : plot_struct.xdata        = STRTRIM(data,2)
      'XLOG'         : plot_struct.xlog         = 1
      'YLOG'         : plot_struct.ylog         = 1
      'NO_YLOG'      : plot_struct.ylog         = 0
      'SMOOTH'       : plot_struct.smooth       = FIX(data)
      'INDICES'      : plot_struct.label_index  = data
      'ASPECT RATIO' : plot_struct.aspect_ratio = FLOAT(data_array)
      'ANNOTATE'     : BEGIN
        i = plot_struct.annotate_n
        j = LONG(data_array[0])
        plot_struct.annotate_code  [i] = j
        plot_struct.annotate_colour[i] = STRTRIM(data_array[1],2)
        plot_struct.annotate_file  [i] = STRTRIM(data_array[2],2)        
        CASE j OF
          1:
          2: BEGIN
            IF (N_ELEMENTS(data_array) GE 4) THEN plot_struct.annotate_step[i] = STRTRIM(data_array[3],2)        
            IF (N_ELEMENTS(data_array) GE 5) THEN BEGIN
              IF (STRTRIM(data_array[4],2) EQ 'label') THEN plot_struct.annotate_label[i] = 1       
            ENDIF
            END
          ELSE:
        ENDCASE
        plot_struct.annotate_n++
        END
      'SET'          : BEGIN
        nset = nset + 1
        plot_struct.case_set_name[nset] = STRTRIM(data,2)
        END
      'CASE CLEAR'   : ncase = 0
      'CASE DEFAULT' : BEGIN
         plot_struct.default      = STRTRIM(data,2)
         plot_struct.case_name[0] = STRTRIM(data,2)  
         END                                         
      'CASE'         : BEGIN
        case_default = plot_struct.default
        FOR i = 0, N_ELEMENTS(data_array)-1 DO BEGIN
          new_case = STRTRIM(data_array[i],2)
          IF (STRMATCH(new_case,'*>*',/FOLD_CASE)) THEN BEGIN
            str = STRSPLIT(data_array[i],'>',/EXTRACT)
            nstr = N_ELEMENTS(str) - 1
            IF (STRLEN(str[nstr]) EQ 0) THEN BEGIN
              new_case = 'none'
            ENDIF ELSE BEGIN
;              new_case = STRMID(plot_struct.case_name[ncase-1],0,  $
;                                STRLEN(case_default)-STRLEN(str[nstr])) + str[nstr]
              new_case = STRMID(case_default,0,STRLEN(case_default)-STRLEN(str[nstr])) + str[nstr]
              case_default = new_case
            ENDELSE
          ENDIF ELSE BEGIN
            case_default = new_case
            IF (i EQ 0 AND nset EQ 0) THEN ncase = 0  ;  NEW RULE! WILL SCREW SOME THINGS UP! -SL, 14/06.2010
          ENDELSE
          plot_struct.case_name[ncase] = new_case
          plot_struct.case_set [ncase] = nset
          ncase = ncase + 1
        ENDFOR
        END
      'LINE SEGMENT' : BEGIN
        plot_struct.line_seg_n++
        IF (plot_struct.line_seg_n GT 5) THEN BEGIN
          PRINT,'ERROR cortex_LoadPlotData: Too many line segments requested (5 max.)'
          RETURN, -1
        ENDIF
        i = plot_struct.line_seg_n
        plot_struct.line_seg[0+(i-1)*4:3+(i-1)*4] = FLOAT(data_array)
        END
      'SHOW'         : BEGIN
         plot_struct.show_tube  [plot_struct.show_n] = LONG(data_array[0])
         plot_struct.show_colour[plot_struct.show_n] =      data_array[1]
         plot_struct.show_n++
         END
;     ------------------------------------------------------------------
      'EXIT': BEGIN
        dummy = cortex_ProcessPlotStruct(plot_struct,plot_array,default,n)
        result = plot_array
        RETURN, 0
        END
;     ------------------------------------------------------------------
      ELSE: BEGIN
        PRINT,'ERROR cortex_LoadPlotData: Unknown plot tag'
        PRINT,'  FILE_NAME= ',file_name
        PRINT,'  TAG      = ',tag
        RETURN,-1
        END
    ENDCASE

  ENDWHILE

END

