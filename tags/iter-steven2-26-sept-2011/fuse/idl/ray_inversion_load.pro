;
;
;
; ======================================================================
;
FUNCTION LoadInversion,file,mode,tag_data

  CASE mode OF
;   --------------------------------------------------------------------
    3: BEGIN
;     Load output from RAY/OUT:

      status = inOpenInterface(file)
      index = inGetData('i')
      nvtx  = inGetData('npts')
      ndata = LONG(N_ELEMENTS(index))
      nmax = MAX(nvtx)
      xvtx = MAKE_ARRAY(ndata,nmax,/FLOAT,VALUE=0.0)
      yvtx = MAKE_ARRAY(ndata,nmax,/FLOAT,VALUE=0.0)
      FOR i = 1, nmax DO BEGIN
         tag_x = 'x' + STRTRIM(STRING(i),2)
         tag_y = 'y' + STRTRIM(STRING(i),2)
         xvtx[*,i-1] = inGetData(tag_x)
         yvtx[*,i-1] = inGetData(tag_y)
      ENDFOR
      data = inGetData(tag_data)
      inCloseInterface

      x = MAKE_ARRAY(ndata,/FLOAT,VALUE=0.0)
      y = MAKE_ARRAY(ndata,/FLOAT,VALUE=0.0)      
      FOR i = 0L, ndata-1 DO BEGIN
        x[i] = MEAN(xvtx[i,0:nvtx[i]-1])
        y[i] = MEAN(yvtx[i,0:nvtx[i]-1])
      ENDFOR

      type_inv = {         $
        version : 1.00  ,  $  ;
        file    : file  ,  $  ;
        index   : index ,  $  ;
        n       : ndata ,  $  ;
        x       : x     ,  $  ;
        y       : y     ,  $  ;
        data    : data  ,  $  ;
        nvtx    : nvtx  ,  $  ;
        xvtx    : xvtx  ,  $  ;
        yvtx    : yvtx  }
      END
;   --------------------------------------------------------------------
    1: BEGIN
;     Load output from OSM:

      status = inOpenInterface(file)
      a = inGetData(datatag)
      IF (a[0] EQ -1) THEN BEGIN
        PRINT,'ERROR LoadInversion: Tag not found in data file'
        PRINT,'  TAG =',datatag
        PRINT,'  FILE=',file
      ENDIF
        PRINT,'  TAG =',datatag
        PRINT,'  FILE=',file
      n = N_ELEMENTS(a)
      type_inv = {                                       $
        n     : -1L      ,  $  ;
        x     : fltarr(n),  $  ;
        y     : fltarr(n),  $  ;
        data  : fltarr(n)}     ;
      type_inv.n = n
      type_inv.x   [0:n-1] = inGetData('x')
      type_inv.y   [0:n-1] = inGetData('y')
      type_inv.data[0:n-1] = a
      inCloseInterface
      END
;   --------------------------------------------------------------------
    0: BEGIN
 
      print,'Load method not active'
      stop 

      status = inOpenInterface(file)

      a = inGetData('cell_data')

      n = N_ELEMENTS(a)

      inv.inv.n = n
      inv.inv.data[0:n-1] = a

      xpts = MAKE_ARRAY(4,n,/FLOAT,value=0.0)  ; Number of polygon sides hardcoded...
      ypts = MAKE_ARRAY(4,n,/FLOAT,value=0.0)
      
      xpts[0,*] = inGetData('cell_vtx_x1')  ; Put into a loop...
      xpts[1,*] = inGetData('cell_vtx_x2')
      xpts[2,*] = inGetData('cell_vtx_x3')
      xpts[3,*] = inGetData('cell_vtx_x4')

      ypts[0,*] = inGetData('cell_vtx_y1')
      ypts[1,*] = inGetData('cell_vtx_y2')
      ypts[2,*] = inGetData('cell_vtx_y3')
      ypts[3,*] = inGetData('cell_vtx_y4')

      FOR i = 0, n-1 DO BEGIN
        inv.inv.x[i] = MEAN(xpts[0:3,i])
        inv.inv.y[i] = MEAN(ypts[0:3,i])
      ENDFOR

      inCloseInterface

      END
;   --------------------------------------------------------------------
    2: BEGIN
;     Load output from RAY - OLD DATA FILE FORMAT:

      fp = 2
      FREE_LUN,fp
      OPENR,fp,file,error=error
      IF (error NE 0) THEN BEGIN
        PRINT,'Troubles mate, with ',file
        STOP
      ENDIF

      READF,fp,buffer
      print,buffer

      buffer_array = STRSPLIT(buffer,/extract)
      n = LONG(buffer_array[0])

      type_inv = {                                       $
        n     : -1L      ,  $  ;
        x     : fltarr(n),  $  ;
        y     : fltarr(n),  $  ;
        data  : fltarr(n)}     ;

      xpts = MAKE_ARRAY(4,/FLOAT,value=0.0)
      ypts = MAKE_ARRAY(4,/FLOAT,value=0.0)

      FOR i = 1L, n DO BEGIN
        buffer = ''
        readf,fp,buffer
;        print,buffer
        buffer_array = strsplit(buffer,/extract)
        ind = LONG(buffer_array[0])
        val = MAX([0.0,FLOAT(buffer_array[1])])
        num  = LONG(buffer_array[2])
        FOR i1 = 0, num-1 DO BEGIN
          xpts[i1] = FLOAT(buffer_array[i1*2+3])
          ypts[i1] = FLOAT(buffer_array[i1*2+4])
        ENDFOR
        xcen = MEAN(xpts[0:num-1])
        ycen = MEAN(ypts[0:num-1])
        type_inv.n = type_inv.n + 1
        type_inv.x   [type_inv.n] = xcen
        type_inv.y   [type_inv.n] = ycen
        type_inv.data[type_inv.n] = val
      ENDFOR
      END

    3: BEGIN
;     Load output from RAY:
;     ------------------------------------------------------------------
      status = inOpenInterface(file)
      a = inGetData(datatag)
      IF (a[0] EQ -1) THEN BEGIN
        PRINT,'ERROR LoadInversion: Tag not found in data file'
        PRINT,'  TAG =',datatag
        PRINT,'  FILE=',file
      ENDIF
        PRINT,'  TAG =',datatag
        PRINT,'  FILE=',file
      n = N_ELEMENTS(a)
      type_inv = {                                       $
        n     : -1L      ,  $  ;
        x     : fltarr(n),  $  ;
        y     : fltarr(n),  $  ;
        data  : fltarr(n)}     ;
      type_inv.n = n
      type_inv.x   [0:n-1] = inGetData('x')
      type_inv.y   [0:n-1] = inGetData('y')
      type_inv.data[0:n-1] = a
      inCloseInterface

      END
  ENDCASE

;  result = CREATE_STRUCT(inv, 'inv', type_inv)
  result = type_inv

;  dim1 = 112  ; Just for debugging, image not longer passed back
;  dim2 = 102
;  image = MAKE_ARRAY(dim1,dim2,/FLOAT,value=0.0)
;  window,1,xsize=dim1,ysize=dim2,retain=2
;  tvscl,image,/order  ; Image drawn from the top-down

  RETURN, result

END
;
; ======================================================================
;
FUNCTION InterpolateToGrid, inv, grid

  mode = 1

  x = inv.x   [0:inv.n-1]
  y = inv.y   [0:inv.n-1]
  z = inv.data[0:inv.n-1]

  CASE mode OF
    1: BEGIN
      TRIANGULATE, x, y, tr, b
      xout = grid.x[0:grid.nx-1]
      yout = grid.y[0:grid.ny-1]
      result = TRIGRID(x, y, z, tr, xout=xout, yout=yout)  
      grid = CREATE_STRUCT(grid,'data',result)
      END
    2: BEGIN
      n = 500
      x0 = grid.x[0]
      y0 = grid.y[0]
      dx = (grid.x[grid.nx-1] - x0) / FLOAT(n-1)
      dy = (grid.y[grid.ny-1] - y0) / FLOAT(n-1)
      print,x0,y0
      print,dx,dy
      result = GRID_TPS(x,y,z,ngrid=[n,n],start=[x0,y0],delta=[dx,dy])
      window,0
      loadct,5
      tvscl, result
      stop
      END
  ENDCASE

  result = { inv : inv, grid : grid }
    
  RETURN, result
END
;
; ======================================================================
;
FUNCTION InterpolateInversion,fname,mode,region, datatag
;
; Load inversion data:
; ----------------------------------------------------------------------
;  help,fname
;  help,mode
;  help,inv
;  help,datatag
  inv = LoadInversion(fname,mode,datatag)
;
; Setup inversion grid interpolation mesh
; ----------------------------------------------------------------------
  n = inv.n
  x1 = MIN(inv.x[0:n-1]) - 0.001
  x2 = MAX(inv.x[0:n-1]) + 0.001
  y1 = MIN(inv.y[0:n-1]) - 0.001
  y2 = MAX(inv.y[0:n-1]) + 0.001

  CASE region OF
   -1: BEGIN  ; Custom
      x1 = 0.23
      x2 = 1.28
      y1 = 0.87
      y2 = 1.92
      END
   -2: BEGIN  ; Custom - upper inner divertor
      x1 = 0.25
      x2 = 0.90
      y1 = 0.95
      y2 = 1.60
      END
    1: BEGIN  ; Upper inner
      x1 = MAX([x1, 0.25])  
      x2 = MIN([x2, 0.70])
      y1 = MAX([y1, 1.05])
      y2 = MIN([y2, 1.40])
      END
    2: BEGIN  ; Lower inner
      x1 = MAX([x1, 0.25])  
      x2 = MIN([x2, 0.80])
      y1 = MAX([y1,-1.60])
      y2 = MIN([y2,-0.95])
      END
    3: BEGIN  ; Upper outer
      x1 = 0.40  ; MAX([x1,0.40])
      x2 = 1.30  ; MIN([x2,1.30])
      y1 = 1.05  ; MAX([y1,1.05])
      y2 = 1.90  ; MIN([y2,1.90])
      END
    4: BEGIN  ; Lower outer
      x1 = MAX([x1, 0.40])  
      x2 = MIN([x2, 1.30])
      y1 = MAX([y1,-1.90])
      y2 = MIN([y2,-1.05])
      END
    ELSE: 
  ENDCASE

  print,'INVERSION BOUNDS:',x1,x2,y1,y2

  MAXGRDPTSX = 500 ; 750

  nx = LONG(MAXGRDPTSX)
  ny = LONG(MAXGRDPTSX)
; Account for "aspect ratio":
  IF ((x2 - x1) GE (y2 - y1)) THEN BEGIN
    ny = LONG(FLOAT(nx) * (y2 - y1) / (x2 - x1))
  ENDIF ELSE BEGIN
    nx = LONG(FLOAT(ny) * (x2 - x1) / (y2 - y1))
  ENDELSE

;  type_grid = {                                       $
;    nx    : nx            ,  $  ;
;    ny    : ny            ,  $  ;
;    x     : fltarr(nx)    ,  $  ;
;    y     : fltarr(ny)    ,  $  ;
;    data  : fltarr(nx,ny) }     ;
;  type_grid.x = (x2 - x1) * FINDGEN(nx) / FLOAT(nx-1) + x1
;  type_grid.y = (y2 - y1) * FINDGEN(ny) / FLOAT(ny-1) + y1
;  result = CREATE_STRUCT(inv, 'grid', type_grid)  

  grid = {   $
    nx : nx                                         ,  $  ;
    ny : ny                                         ,  $  ;
    x  : (x2 - x1) * FINDGEN(nx) / FLOAT(nx-1) + x1 ,  $  ;
    y  : (y2 - y1) * FINDGEN(ny) / FLOAT(ny-1) + y1 }

  result = InterpolateToGrid(inv, grid)

  RETURN, result
END
;
; ======================================================================
;
FUNCTION GetInversion,    $
    file                   ,  $  ;
    camera=camera          ,  $  ;
    shot=shot              ,  $  ;
    frame=frame            ,  $  ;
    channel=channel        ,  $  ;
    suffix=suffix          ,  $  ;
    ext=ext                ,  $  ;
    maxscale=maxscale      ,  $  ;
    casename=casename      ,  $  ;
    datatag=datatag        ,  $  ;
    region=region          ,  $  ;  Region of interest
    full=full              ,  $  ;
    path=path              ,  $  ;
    chisq_limit=chisq_limit      ;
;
;
; 
  IF (NOT KEYWORD_SET(ext)) THEN ext = ['cgm']

  IF (NOT KEYWORD_SET(region     )) THEN region  = 0
  IF (NOT KEYWORD_SET(camera     )) THEN camera  = 'DIVCAM'
  IF (NOT KEYWORD_SET(shot       )) THEN shot    = -1
  IF (NOT KEYWORD_SET(frame      )) THEN frame   = -1
  IF (NOT KEYWORD_SET(channel    )) THEN channel = -1
  IF (NOT KEYWORD_SET(machine    )) THEN machine = 'MAST'
  IF (NOT KEYWORD_SET(chisq_limit)) THEN chisq_limit = 0.2
  IF (NOT KEYWORD_SET(path       )) THEN path = './data/'
;
; Define inversion data structure:
; ----------------------------------------------------------------------
  type_inversion = {                                       $
    machine : machine                      ,  $  ;
    camera  : camera                       ,  $  ;
    shot    : 0                            ,  $  ;
    frame   : 0                            ,  $  ;
    channel : 0                            ,  $  ;
    time    : 0.0                          ,  $  ;
    region  : region                       ,  $  ;
    method  : strarr(512)                  }

; Load data:
; ----------------------------------------------------------------------
  
  inv = type_inversion




  blah = 1

  CASE 1 OF
;;   --------------------------------------------------------------------
;    (ext EQ 'ipx'): BEGIN
;      prefix[0] = '/net/fuslsa/data/MAST_IMAGES/rda/rda0'
;      prefix[1] = '/net/fuslsa/data/MAST_IMAGES/rda/rdb0'
;      file_name = prefix[channel]+STRTRIM(STRING(shot),1)+'.ipx'
;      desc = ipx_open(file_name)
;      image = ipx_frame(desc,frame,TIME=t)
;      inv.grid.nx = 1000
;      inv.grid.ny = 751 
;      END
;;   --------------------------------------------------------------------
;    (ext EQ 'sav'): BEGIN
;      print,'SAV not ready'
;      stop
;      RESTORE, filename=fname[i], /verbose
;      dim_x = N_ELEMENTS(image_data[*,0])
;      dim_y = N_ELEMENTS(image_data[0,*])
;      image1 = MAKE_ARRAY([dim_x,dim_y],/FLOAT,value=0.0)
;      FOR j = 0, dim_y-1 DO image1[*,j] = image_data[*,dim_y-1-j]
;      END
;;   --------------------------------------------------------------------
;    (ext EQ 'mxe'): BEGIN
;;      Hack for getting linear Dalpha camera data...
;;      a=loadinversion('/home/slisgo/divimp/shots/mast/images/17469_127_1.cgm',2,inv) 
;;      return,a
;      inv.method = ext
;      IF (STRLEN(ext) EQ 3) THEN ext = '.' + ext
;      fname = path + camera + '_' +                   $
;              STRTRIM(STRING(shot   ),1) + '_' +  $
;              STRTRIM(STRING(frame  ),1) + '_' +  $
;              STRTRIM(STRING(channel),1) + ext
;;      PRINT,fname
;      inv = InterpolateInversion(fname,2,inv)
;;      inv.shot    = shot
;;      inv.frame   = frame
;;      inv.channel = channel
;;      window,0,xsize=inv.grid.nx,ysize=inv.grid.ny
;;      loadct, 6 ; 3
;;      tvscl,inv.grid.data
;;      STOP 
;      END
;   --------------------------------------------------------------------
    (ext EQ 'cgm'): BEGIN

;      PRINT,file

      result = InterpolateInversion(file,3,region,'data')   

;      window,0,xsize=inv.grid.nx,ysize=inv.grid.ny
;      loadct, 6
;      tvscl,inv.grid.data ; <2.0E+21
      END
;;   --------------------------------------------------------------------
;    (ext EQ 'osm'): BEGIN
;      path = '~/fuse_data/mast/images/'
;      IF (STRLEN(ext) EQ 3) THEN ext = '.' + ext
;      IF (NOT KEYWORD_SET(casename)) THEN BEGIN
;        PRINT,'ERROR GetInversion: OSM case name not specified'
;        STOP
;      ENDIF
;;      fname = path+STRTRIM(STRING(shot   ),1)+'_'+  $
;;                   STRTRIM(STRING(frame  ),1)+'_'+  $
;;                   STRTRIM(STRING(channel),1)+ext
;      fname = path + casename + '.osm'
;;      PRINT,fname,' ',datatag
;      inv = InterpolateInversion(fname,1,inv,datatag)   
;;      window,0,xsize=inv.grid.nx,ysize=inv.grid.ny
;;      loadct, 6
;;      tvscl,inv.grid.data ; <2.0E+21
;      END
;   ------------------------------------------------------------------------
    ELSE: BEGIN 
      PRINT,'ERROR ray_GetInversion: Unrecognized file extension'
      STOP
      END
  ENDCASE

;  inv = ProcessImage(inv, chisq_limit)

  IF (region LT 0) THEN BEGIN  ; For a hardcoded / custom region of interest...
    LOADCT,5
    IF (ext EQ 'osm') THEN BEGIN
      fname = './images/MODEL_'+casename+'_'+datatag
    ENDIF ELSE BEGIN
      fname = path+STRTRIM(STRING(shot   ),1)+'_'+  $
                   STRTRIM(STRING(frame  ),1)+'_'+  $
                   STRTRIM(STRING(channel),1)
    ENDELSE
    image_data = inv.grid.data
    maxval= MAX(image_data)
    IF (KEYWORD_SET(maxscale)) THEN BEGIN 
      IF (maxscale LT 0.0) THEN maxscale = ABS(maxscale) * maxval
      i = WHERE(image_data GE maxscale)
      image_data[i] = maxscale
      print,'PNG: maxval:',maxscale
    ENDIF ELSE BEGIN
      print,'PNG: maxval:',maxval
    ENDELSE
    image_byte = BYTSCL(image_data)
    TV,image_byte
    TVLCT, red, green, blue, /GET
    imageRGB = BYTARR(3, inv.grid.nx, inv.grid.ny)
    imageRGB[0, *, *] = red  [image_byte]  
    imageRGB[1, *, *] = green[image_byte]  
    imageRGB[2, *, *] = blue [image_byte] 
    fname=fname+'_'+ext+'.png'
    print,fname
    WRITE_PNG, fname, imageRGB
;    WRITE_JPEG, fname+'_inv.jpg', imageRGB, TRUE = 1, QUALITY = 100

    RETURN, -1.0
  ENDIF

  grid = result.grid
  raw  = result.inv

  result = {  $
    version : 1.00                               ,  $
    machine : machine                            ,  $
    camera  : camera                             ,  $
    shot    : shot                               ,  $
    frame   : frame                              ,  $
    channel : channel                            ,  $
    region  : region                             ,  $
    x       : grid.x   [0:grid.nx-1]             ,  $
    y       : grid.y   [0:grid.ny-1]             ,  $
    data    : grid.data[0:grid.nx-1,0:grid.ny-1] }

  IF (KEYWORD_SET(full)) THEN result = CREATE_STRUCT(result,'raw',raw,'grid',grid)

  RETURN, result
END
;
; ======================================================================
;


