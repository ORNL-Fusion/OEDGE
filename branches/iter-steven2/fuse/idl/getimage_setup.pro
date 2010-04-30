;
; ======================================================================
; 
FUNCTION GetFilterDataStructure

  struct = {                    $
    verison      : 1.0         ,  $
    id           : 0           ,  $
    line         : -1.0        ,  $
    tag          : '         ' ,  $
    supplier     : '         ' ,  $
    sn           : '         ' ,  $  ; Serial no
    cwl          : 0.0         ,  $ 
    fwhm         : 0.0         ,  $ 
    transmission : 0.0         ,  $  ; 
    nd           : 0.0         ,  $  ; Neutral density filter(s)
    cavities     : 0           ,  $  ; 
    index        : 0.0         ,  $
    nshape       : 0           ,  $
    xshape       : FLTARR(10)  ,  $
    yshape       : FLTARR(10)  ,  $
    n            : 20          ,  $
    nr           : FLTARR(20)  ,  $
    tr           : FLTARR(20)}

  RETURN, struct
END
;
; ======================================================================
;
FUNCTION GetImageCalibrationStructure

  n       = 64
  degrees = 8      ; needs a better fitting scheme... always wobbles near the top

  struct ={                             $
    version  : 1.0                   ,  $  ;

    file     : STRING(' ',FORMAT='(A256)'),  $  ;  Name of image file used for calibration

;   Absolute:
    serial_no   : 0                   ,  $  ;
    mode        : -1                  ,  $  ;
    power       : -1.0                ,  $  ;
    power_units : 'Watts st-1 m-2'    ,  $  ;
    counts      : -1.0                ,  $  ;
    gain        : -1.0                ,  $  ;
    shutter     : -1.0                ,  $  ;
    cube        : -1.0                ,  $  ;
    fudge       : -1.0                ,  $  ;
    date        : -1L                 ,  $  ;

;   Window:
    window       : -1.0               ,  $  ;
    window_tag   : STRING(' ',FORMAT='(A256)') ,  $  ;  
    w_n          : -1                 ,  $  ;
    w_shot_range : [-1L,-1L]          ,  $  ;
    w_interpol   : -1                 ,  $  ; Interpolation scheme
    w_wlngth     : FLTARR(100)        ,  $  ;
    w_trans      : FLTARR(100,2)      ,  $  ;
    w_tr         : 1.0                ,  $  ;  *** remove ***

    cube_n        : -1                 ,  $  ;
    cube_wlngth   : FLTARR(100)        ,  $  ;
    cube_straight : FLTARR(100)        ,  $  ;
    cube_orthog   : FLTARR(100)        ,  $  ;
    cube_tr       : 1.0                ,  $  ; *** remove ***

;   Mask:

;   Vignette: (restored from vinetting calibration file)                                                 
    v_date   :  0L                   ,  $  ;
    v_file   : STRING(' ',FORMAT='(A256)') ,  $  ;  
    path     : STRING(' ',FORMAT='(A256)') ,  $  ;  
    n        : n                     ,  $  ;
    angle    : FLTARR(1000)          ,  $  ;
    max_r    : FLTARR(1000)          ,  $  ;
    ndeg     : degrees               ,  $  ;
    A        : FLTARR(1000,12)          $  ;
    }    ; 

  RETURN, struct
END
;
; ======================================================================
;
FUNCTION GetImageDataStructure, xdim, ydim

  print,'getting image data structure',xdim,ydim

  xthumb = MAX([1,xdim/5])
  ythumb = MAX([1,ydim/5])

  calibration = GetImageCalibrationStructure()
  filter      = GetFilterDataStructure()

  struct ={                             $
    version  : 1.0                   ,  $  ;  
    warnings : 0                     ,  $  ;  

    date      : -1L                   ,  $  ;  Data that data was recorded
    date_time : STRING(' ',FORMAT='(A256)') ,  $  ;  Data that data was recorded
    device    : STRING(' ',FORMAT='(A256)') ,  $  ;  Device name
    window    : STRING(' ',FORMAT='(A256)') ,  $  ;  Port / window, ie. "HU12A"
    line      : -1.0                        ,  $  ;  Emission line to the nearest (for file access) 
    shot      : -1L                         ,  $  ;  Plasma discharge
    time      : -1.0                        ,  $  ;  Time stamp for image  (...need to make an interval...)
    frame     : -1L                         ,  $  ;  Frame number in series
                                                   
    camera    : STRING(' ',FORMAT='(A256)') ,  $  ;  Camera identification code      
    channel   : -1                          ,  $  ;  Channel number for multi-camera systems     
    shutter   : -1.0                        ,  $  ;  
    gain      : [-1.0,-1.0]                 ,  $  ;  
    file      : 'none'                      ,  $  ;  Name of image file if file specified directly
			           
    registration : 0                 ,  $  ;  View registration information

    filter   : filter                ,  $  ;  Filter structure
                                                   
    format   : 'ext'                 ,  $  ;  File type of image source
    depth    : -1                    ,  $  ;  Image bit depth          
			           
    scl_pho  : -1.0                  ,  $  ;  Scale factor to give photons m-2 s-1 (default)
    scl_pow  : -1.0                  ,  $  ;  Scale factor to give W m-2
                                                   
    xdim     : xdim                  ,  $  ;        
    ydim     : ydim                  ,  $  ; 

    xbin     : 1                     ,  $  ; 
    ybin     : 1                     ,  $  ; 
                                                   
    xwin     : [1,xdim]              ,  $  ;  Pixel x-range of image within [xdim,ydim] image space
    ywin     : [1,ydim]              ,  $  ; 
                                                   
    xcen     : xdim / 2              ,  $  ;  Pixel x-location of the centre of the view
    ycen     : ydim / 2              ,  $  ; 
    radius   : (xdim + ydim) / 4     ,  $  ;  Radial extent of valid image data about XCEN,YCEN
                                               
    raw       : fltarr(xdim  ,ydim  ),  $  ;  Raw image data
    mask      : fltarr(xdim  ,ydim  ),  $  ;  Mask for raw image data
    map_r     : fltarr(xdim  ,ydim  ),  $  ;        
    map_theta : fltarr(xdim  ,ydim  ),  $  ;        
;    data      : fltarr(xdim  ,ydim  ),  $  ;    *** THIS WAS COMMENT OUT, WHY? 12.5.09 ***
    vignette  : fltarr(xdim  ,ydim  ),  $  ; 
    blueshift : fltarr(xdim  ,ydim  ),  $  ; 

    cal       : calibration          ,  $  ;

    thumb     : fltarr(xthumb,ythumb)   $  ; 
    }    ; 

  struct.mask = 1

  print,'getting image data structure: done'

  RETURN, struct

END
