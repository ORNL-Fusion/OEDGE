;
; ======================================================================
;
FUNCTION cortex_LoadAnnotations, code, file

;  file = cortex_UpdateFile(file)

;print,code
;print,file

  CASE code OF
;   --------------------------------------------------------------------
    1: BEGIN

      fp = 3
      FREE_LUN, fp
      OPENR, fp, file, error=err
      IF (err NE 0) THEN BEGIN
        PRINT,'ERROR cortex_LoadAnnotationData: Unable to open data file
        PRINT,' FILE=',file
        RETURN, -1
      ENDIF

      buffer = ' '
      first_line = 1

      WHILE NOT EOF(fp) DO BEGIN 
        READF,fp,buffer      
        buffer = STRTRIM(buffer,2)
; print,'buffer:',buffer
        IF (STRPOS(buffer,'#') EQ 0 OR STRPOS(buffer,'*') EQ 0) THEN CONTINUE
        str = STRSPLIT(buffer,' '+STRING(9B),/EXTRACT)  ; STRING(9B) is the ASCII TAB character
        IF (N_ELEMENTS(str) EQ 1) THEN CONTINUE         ; weak check to avoid blank lines
        IF (first_line) THEN BEGIN
          x = FLOAT(str[0])
          y = FLOAT(str[1])
        ENDIF ELSE BEGIN
          x = [x,FLOAT(str[0])]
          y = [y,FLOAT(str[1])]
        ENDELSE
        first_line = 0
      ENDWHILE

      FREE_LUN, fp

      result = {x : x, y : y}

      END
;   --------------------------------------------------------------------
    2: BEGIN
      file = cortex_UpdateFile(file)
      status = inOpenInterface(file)
      IF (status LT 0) THEN BEGIN
        result = CREATE_STRUCT('version',0.0,'file','none')
        RETURN, result
      ENDIF      
      x1 = inGetData('X_1')
      y1 = inGetData('Y_1')
      x2 = inGetData('X_2')
      y2 = inGetData('Y_2')
      inCloseInterface  
      result = {x1 : x1, y1 : y1, x2 : x2, y2 : y2}
      END
;   --------------------------------------------------------------------
    3: BEGIN  ; EIRENE triangles
      file = cortex_UpdateFile(file)
      status = inOpenInterface(file)
      IF (status LT 0) THEN BEGIN
        result = CREATE_STRUCT('version',0.0,'file','none')
        RETURN, result
      ENDIF      
      v1 = inGetData('VERTEX_1')
      v2 = inGetData('VERTEX_2')
      v3 = inGetData('VERTEX_3')
;      s1 = inGetData('SURFACE_1')
;      s2 = inGetData('SURFACE_2')
;      s3 = inGetData('SURFACE_3')
      x  = inGetData('VERTEX_X')
      y  = inGetData('VERTEX_Y')
      inCloseInterface  
      v = [[v1],[v2],[v3]]
      s = -1
;      s = [[s1],[s2],[s3]]
      result = { version : 1.0 , file : file, v : v, surface : s, x : x , y : y }
      END
;   --------------------------------------------------------------------
    4: BEGIN  ; field line
      file = cortex_UpdateFile(file)
      status = inOpenInterface(file)
      IF (status LT 0) THEN BEGIN
        result = CREATE_STRUCT('version',0.0,'file','none')
        RETURN, result
      ENDIF      
      index = inGetData('INDEX')
      x1    = inGetData('X_1')
      y1    = inGetData('Y_1')
      z1    = inGetData('Z_1')
      inCloseInterface  
      v = [[x1],[y1],[z1]]
      result = { version : 1.0, file : file, index : index, v : v }
      END
;   --------------------------------------------------------------------
    ELSE:
  ENDCASE
  
  RETURN, result
END

;
; ======================================================================
;
FUNCTION cortex_LoadAnnotationData, plot,file

;  annotate = { version : 1.0 } 
;  status_annotate = -1
  annotate_data = -1
  FOR i = 0, plot.annotate_n-1 DO BEGIN
    CASE plot.annotate_code[i] OF
      1: annotate_data = cortex_LoadAnnotations(plot.annotate_code[i],       plot.annotate_file[i]) 
      2: annotate_data = cortex_LoadAnnotations(plot.annotate_code[i],file + plot.annotate_file[i]) 
      3: annotate_data = cortex_LoadAnnotations(plot.annotate_code[i],file + plot.annotate_file[i]) 
      ELSE: BEGIN
        PRINT, 'ERROR cortex_LoadAnnotationData: Unknown annotation code'
        PRINT, '  PLOT TAG = ',plot.tag
        PRINT, '  CODE     = ',plot.annotate_code[i]
        annotate_data = -1
        END
    ENDCASE
    name = 'data' + STRING(i+1,FORMAT='(I0)')
    IF (i EQ 0) THEN annotate = CREATE_STRUCT(         name,annotate_data)  ELSE  $
                     annotate = CREATE_STRUCT(annotate,name,annotate_data)
;    status_annotate = 0
  ENDFOR 

  IF (plot.annotate_n EQ 1) THEN annotate = CREATE_STRUCT(annotate,'dummy',0)  ; Can remove in future...

  result = annotate

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_LoadWall, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  class  = inGetData('WALL_CLASS')  
  group  = inGetData('WALL_GROUP')  
  index  = inGetData('WALL_INDEX')  
  tube   = inGetData('WALL_TUBE')  
  target = inGetData('WALL_TARGET')
  n = N_ELEMENTS(class)
  v1 = MAKE_ARRAY(2,n,/FLOAT,VALUE=0.0)
  v2 = v1
  v1[0,*] = inGetData('WALL_V1_X')  
  v1[1,*] = inGetData('WALL_V1_Y')  
  v2[0,*] = inGetData('WALL_V2_X')  
  v2[1,*] = inGetData('WALL_V2_Y')  

  inCloseInterface

  result = {           $
    file   : file   ,  $
    n      : n      ,  $
    class  : class  ,  $
    group  : group  ,  $
    index  : index  ,  $
    tube   : tube   ,  $
    target : target ,  $
    v1     : v1     ,  $
    v2     : v2     }
  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_LoadFluidGrid, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  grd_isep = inGetData('GRD_ISEP')  
  grd_ipfz = inGetData('GRD_IPFZ')  

  grd_r0   = inGetData('GRD_R0')  
  grd_z0   = inGetData('GRD_Z0')  

  grd_rxpt = inGetData('GRD_RXPT')  
  grd_zxpt = inGetData('GRD_ZXPT')  

  tube_n    = inGetData('TUBE_N')  
  tube_psin = inGetData('TUBE_PSIN')  
  tube_rho  = inGetData('TUBE_RHO')  
  tube_l    = inGetData('TUBE_L')  

  obj_index = inGetData('OBJ_INDEX')  
  obj_n = N_ELEMENTS(obj_index)
  obj_cell  = inGetData('OBJ_IND_CELL')  
  obj_pos   = inGetData('OBJ_IND_POS' )  
  obj_tube  = inGetData('OBJ_IND_TUBE')  
  obj_nside = inGetData('OBJ_NSIDE'   )  

  obj_iside = MAKE_ARRAY(4,obj_n,/LONG,VALUE=0)
  obj_iside[0,*] = inGetData('OBJ_ISIDE_1')  
  obj_iside[1,*] = inGetData('OBJ_ISIDE_2')  
  obj_iside[2,*] = inGetData('OBJ_ISIDE_3')  
  obj_iside[3,*] = inGetData('OBJ_ISIDE_4')  

  obj_omap = MAKE_ARRAY(4,obj_n,/LONG,VALUE=0)
  obj_omap[0,*] = inGetData('OBJ_OMAP_1')  
  obj_omap[1,*] = inGetData('OBJ_OMAP_2')  
  obj_omap[2,*] = inGetData('OBJ_OMAP_3')  
  obj_omap[3,*] = inGetData('OBJ_OMAP_4')  

  srf_index = inGetData('SRF_INDEX')
  nsrf = N_ELEMENTS(srf_index)
  srf_type = inGetData('SRF_TYPE')
  srf_nvtx = inGetData('SRF_NVTX')
  srf_ivtx = MAKE_ARRAY(2,nsrf,/LONG,VALUE=0)
  srf_ivtx[0,*] = inGetData('SRF_IVTX_1')
  srf_ivtx[1,*] = inGetData('SRF_IVTX_2')

  vtx_index = inGetData('VTX_INDEX')
  nvtx = N_ELEMENTS(vtx_index)
  vtx = MAKE_ARRAY(3,nvtx,/DOUBLE,VALUE=0.0)
  vtx[0,*] = inGetData('VTX_1')
  vtx[1,*] = inGetData('VTX_2')
  vtx[2,*] = inGetData('VTX_3')

  inCloseInterface

;  help,nsrf

; Setup index mapping for surfaces:   *** NOT CURRENTLY IN USE ***
  maxindex = MAX(srf_index)
  srf_map = MAKE_ARRAY(maxindex+1,/LONG,VALUE=0)
  FOR i = 0L, nsrf-1L DO srf_map[srf_index[i]] = i

; Setup index mapping for vertices:
  maxindex = MAX(vtx_index)
  vtx_map = MAKE_ARRAY(maxindex+1,/LONG,VALUE=0)
  FOR i = 0, nvtx-1 DO vtx_map[vtx_index[i]] = i

  result = {                 $
    version   : 1.0       ,  $
    file      : file      ,  $
    isep      : grd_isep  ,  $
    ipfz      : grd_ipfz  ,  $
    r0        : grd_r0    ,  $
    z0        : grd_z0    ,  $
    rxpt      : grd_rxpt  ,  $
    zxpt      : grd_zxpt  ,  $
    tube_n    : tube_n    ,  $
    tube_psin : tube_psin ,  $
    tube_rho  : tube_rho  ,  $
    tube_l    : tube_l    ,  $
    obj_n     : obj_n     ,  $
    obj_index : obj_index ,  $
    obj_cell  : obj_cell  ,  $
    obj_pos   : obj_pos   ,  $
    obj_tube  : obj_tube  ,  $
    obj_nside : obj_nside ,  $
    obj_iside : obj_iside ,  $
    obj_omap  : obj_omap  ,  $
    nsrf      : nsrf      ,  $
    srf_index : srf_index ,  $
    srf_map   : srf_map   ,  $
    srf_type  : srf_type  ,  $ 
    srf_nvtx  : srf_nvtx  ,  $
    srf_ivtx  : srf_ivtx  ,  $
    nvtx      : nvtx      ,  $
    vtx_index : vtx_index ,  $
    vtx_map   : vtx_map   ,  $
    vtx       : vtx }

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_LoadFluidGrid_Debug, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  ik   = inGetData('IK')  
  ir   = inGetData('IR')  
  nv   = inGetData('NV')
  rv_1 = inGetData('RV_1')
  zv_1 = inGetData('ZV_1')
  rv_2 = inGetData('RV_2')
  zv_2 = inGetData('ZV_2')
  rv_3 = inGetData('RV_3')
  zv_3 = inGetData('ZV_3')
  rv_4 = inGetData('RV_4')
  zv_4 = inGetData('ZV_4')

  inCloseInterface

  obj_n = N_ELEMENTS(ik)
  nside = MAX(nv)  

  rv = MAKE_ARRAY(nside,obj_n,/FLOAT,VALUE=0.0)
  rv[0,*] = rv_1
  rv[1,*] = rv_2
  rv[2,*] = rv_3
  IF (nside GE 4) THEN rv[3,*] = rv_4

  zv = MAKE_ARRAY(nside,obj_n,/FLOAT,VALUE=0.0)
  zv[0,*] = zv_1
  zv[1,*] = zv_2
  zv[2,*] = zv_3
  IF (nside GE 4) THEN zv[3,*] = zv_4

;  help,rv
;  help,zv
;stop

  grd_isep = 1
  grd_ipfz = 1

  tube_n    = MAX(ir)
  tube_psin = MAKE_ARRAY(1,tube_n,/FLOAT,VALUE=1.0)
  tube_rho  = MAKE_ARRAY(1,tube_n,/FLOAT,VALUE=1.0)
  tube_l    = MAKE_ARRAY(1,tube_n,/FLOAT,VALUE=1.0)

  obj_index = MAKE_ARRAY(obj_n,/LONG,VALUE=0)
  obj_cell  = ik
  obj_pos   = ik
  obj_tube  = ir
  obj_nside = nv
  obj_iside = MAKE_ARRAY(nside,obj_n,/LONG,VALUE=0)
  obj_omap  = MAKE_ARRAY(nside,obj_n,/LONG,VALUE=0)

  nsrf = obj_n * nside
  nvtx = nsrf * 2L

;print,tube_n,obj_n,nside,nsrf,nvtx


  srf_index = MAKE_ARRAY(  nsrf,/LONG  ,VALUE=0  )
  srf_type  = MAKE_ARRAY(  nsrf,/LONG  ,VALUE=0  )
  srf_nvtx  = MAKE_ARRAY(  nsrf,/LONG  ,VALUE=0  )
  srf_ivtx  = MAKE_ARRAY(2,nsrf,/LONG  ,VALUE=0  )
  vtx_index = MAKE_ARRAY(  nvtx,/LONG  ,VALUE=0  )
  vtx       = MAKE_ARRAY(3,nvtx,/DOUBLE,VALUE=0.0)

;  print,'nobj,nside:',obj_n,nside

;  print,'nsrf,vtx:',nsrf,nvtx

  nvtx = -1L
  nsrf = -1L

;  help,srf_ivtx

  FOR i = 0L, LONG(obj_n-1) DO BEGIN
    FOR j = 0L, LONG(nside-1) DO BEGIN
      nsrf = nsrf + 1

      obj_iside[j,i] = nsrf + 1L
      obj_omap [j,i] = 0

      srf_index[nsrf] = nsrf + 1L
      srf_type [nsrf] = 0
      srf_nvtx [nsrf] = 2

      i1 = j
      i2 = j + 1L
      IF (i2 GT nside-1) THEN i2 = 0L

;print,'nvtx',nsrf,nvtx,i,j,obj_n,nside,obj_n*nside,obj_n*nside*2
;help,vtx_index
      nvtx = nvtx + 1L
      vtx_index[nvtx] = nvtx
      vtx[0,nvtx] = rv[i1,i]
      vtx[1,nvtx] = zv[i1,i]
      srf_ivtx[0,nsrf] = nvtx

      nvtx = nvtx + 1L
      vtx_index[nvtx] = nvtx
      vtx[0,nvtx] = rv[i2,i]
      vtx[1,nvtx] = zv[i2,i]
      srf_ivtx[1,nsrf] = nvtx
    ENDFOR
  ENDFOR

  nsrf = nsrf + 1L
  nvtx = nvtx + 1L

;  print,'nsrf,vtx:',nsrf,nvtx

; Setup index mapping for surfaces:   *** NOT CURRENTLY IN USE ***
  maxindex = MAX(srf_index)
  srf_map = MAKE_ARRAY(maxindex+1,/LONG,VALUE=0)
  FOR i = 0L, nsrf-1 DO BEGIN
    srf_map[srf_index[i]] = i
  ENDFOR

; Setup index mapping for vertices:
  maxindex = MAX(vtx_index)
  vtx_map = MAKE_ARRAY(maxindex+1,/LONG,VALUE=0)
  FOR i = 0L, LONG(nvtx-1) DO BEGIN
    vtx_map[vtx_index[i]] = i
  ENDFOR

  result = {                 $
    isep      : grd_isep  ,  $
    ipfz      : grd_ipfz  ,  $
    tube_n    : tube_n    ,  $
    tube_psin : tube_psin ,  $
    tube_rho  : tube_rho  ,  $
    tube_l    : tube_l    ,  $
    obj_n     : obj_n     ,  $
    obj_index : obj_index ,  $
    obj_cell  : obj_cell  ,  $
    obj_pos   : obj_pos   ,  $
    obj_tube  : obj_tube  ,  $
    obj_nside : obj_nside ,  $
    obj_iside : obj_iside ,  $
    obj_omap  : obj_omap  ,  $
    nsrf      : nsrf      ,  $
    srf_index : srf_index ,  $
    srf_map   : srf_map   ,  $
    srf_type  : srf_type  ,  $ 
    srf_nvtx  : srf_nvtx  ,  $
    srf_ivtx  : srf_ivtx  ,  $
    nvtx      : nvtx      ,  $
    vtx_index : vtx_index ,  $
    vtx_map   : vtx_map   ,  $
    vtx       : vtx }

  RETURN, result
END
;
; ======================================================================
;
FUNCTION cortex_LoadFluidWall_Debug, file

  file = cortex_UpdateFile(file)

  status = inOpenInterface(file)
  IF (status LT 0) THEN BEGIN
    result = CREATE_STRUCT('version',0.0,'file','none')
    RETURN, result
  ENDIF

  dummy = inGetData('x1')  

  n = N_ELEMENTS(dummy)

  v1 = MAKE_ARRAY(2,n,/FLOAT,VALUE=0.0)
  v2 = v1
  v1[0,*] = inGetData('x1')  
  v1[1,*] = inGetData('y1')  
  v2[0,*] = inGetData('x2')  
  v2[1,*] = inGetData('y2')  

  inCloseInterface

  cl = MAKE_ARRAY(n,/LONG,VALUE=1)  

  result = {        $
     n      : n  ,  $
     class  : cl ,  $
     group  : 0  ,  $
     index  : 0  ,  $
     tube   : 0  ,  $
     target : 0  ,  $
     v1     : v1 ,  $
     v2     : v2    }

  RETURN, result
END
;
; ======================================================================
;
PRO cortex_load_geometry
;  interface
;  path = '/home/ITER/lisgos/divimp/results/'
;  file = path + 'i-new-0006a.idl.fluid_grid'
;  geo = LoadFluidGridGeometry(file)
END
;
; ======================================================================
;

