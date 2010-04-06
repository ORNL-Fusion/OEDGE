;
;
;
;
;
;
;
;
;
;
;
;=======================================================================
;
FUNCTION GetImage,           $
   calibrate  = calibrate,   $
   binary     = binary,      $
   debug      = debug,       $
   background = background,  $ 
   maskleft   = maskleft,    $
   maskright  = maskright,   $
   masktop    = masktop,     $
   maskbottom = maskbottom,  $
   maskradius = maskradius,  $
   maskpoly   = maskpoly,    $
   maskbox    = maskbox,     $
   png        = png,         $
   flip       = flip,        $
   centre     = centre,      $
   colour     = colour,      $
   brighten   = brighten,    $
   clean      = clean,       $
   rotate     = rotate,      $
   shot       = shot,        $
   channel    = channel,     $   ; channel
   frame      = frame,       $
   device     = device,      $
   camera     = camera,      $   ; MAST: 'divcam' (default), 'zebra', 'rgb'
   scale      = scale ,      $
   file       = file,        $
   analyse    = analyse,     $
   shift      = shift,       $
   radius     = radius,      $
   path       = path,        $
   save       = save,        $
   sname      = sname,       $
   png_save   = png_save,    $
   raw        = raw,         $
   date       = date,        $
   line       = line,        $
   type       = type,        $
   filter     = filter,      $
   depth      = depth,       $
   window     = window,      $
   nocal      = nocal,       $
   plots      = plots
;
; ----------------------------------------------------------------------
; Process keywords and assign defaults, as necesary:
; ----------------------------------------------------------------------
;
  IF (N_ELEMENTS(binary   ) EQ 0) THEN binary     = 0  
  IF (N_ELEMENTS(debug    ) EQ 0) THEN debug      = 0  
  IF (N_ELEMENTS(png      ) EQ 0) THEN png        = 0  
  IF (N_ELEMENTS(flip     ) EQ 0) THEN flip       = 0  
  IF (N_ELEMENTS(centre   ) EQ 0) THEN centre     = 0  

  IF (NOT KEYWORD_SET(date      )) THEN date       = -1
  IF (NOT KEYWORD_SET(line      )) THEN line       = -1.0
  IF (NOT KEYWORD_SET(shot      )) THEN shot       = -1
  IF (NOT KEYWORD_SET(frame     )) THEN frame      = -1
  IF (NOT KEYWORD_SET(scale     )) THEN scale      = 1.0
  IF (NOT KEYWORD_SET(device    )) THEN device     = 'MAST'
  IF (NOT KEYWORD_SET(path      )) THEN path       = '$FUSEHOME_DATA/'+STRLOWCASE(device)+'/images'
  IF (NOT KEYWORD_SET(calibrate )) THEN calibrate  = 0
  IF (NOT KEYWORD_SET(clean     )) THEN clean      = 0
  IF (NOT KEYWORD_SET(background)) THEN background = 0

  IF (KEYWORD_SET(file)) THEN BEGIN
    device = 'FILE'
    IF (NOT KEYWORD_SET(window )) THEN window = 'none'
    IF (NOT KEYWORD_SET(channel)) THEN channel = -1
  ENDIF ELSE BEGIN
;  ENDIF
;  IF (NOT KEYWORD_SET(file)) THEN BEGIN
    IF (NOT KEYWORD_SET(shot) AND NOT calibrate) THEN BEGIN
      PRINT,'Error GetImage: Shot number not specified'
      STOP
    ENDIF
    IF (NOT KEYWORD_SET(frame) AND NOT calibrate) THEN BEGIN
      PRINT,'Error GetImage: Frame number not specified'
      STOP
    ENDIF
  ENDELSE
;  ENDIF

  IF (NOT KEYWORD_SET(channel)) THEN BEGIN
    PRINT,'Error GetImage: CHANNEL not specified'
    RETURN,-1
  ENDIF

  IF (NOT KEYWORD_SET(camera)) THEN BEGIN
    CASE device OF
      'FILE': camera = 'unknown'
      'MAST': camera = 'DIVCAM'
      'CMOD': camera = 'CMOD_CAMERA01'
      ELSE: BEGIN
        PRINT,'ERROR GetImage: CAMERA not specified' 
        RETURN,-1
        END
    ENDCASE
  ENDIF

  CASE camera OF
    'DIVCAM': BEGIN
      CASE (channel) OF
       -1: channel_tag = 'unknown'
        1: channel_tag = 'a'
        2: channel_tag = 'b'
        3: channel_tag = 'c'
        4: channel_tag = 'd'
      ENDCASE
      END
    'RGB': BEGIN
      CASE (channel) OF
        1: channel_tag = 'L BREM R'
        2: channel_tag = 'L BREM G'
        3: channel_tag = 'L BREM B'
        4: channel_tag = 'U Dalpha'
        5: channel_tag = 'U C5+'
        6: channel_tag = 'U He+'
      ENDCASE
      END
    ELSE: channel_tag = 'unknown' 
  ENDCASE  
;
; Echo inputs:
;
  IF (KEYWORD_SET(debug)) THEN BEGIN
    PRINT,'shot       = ',shot
    PRINT,'frame      = ',frame
    PRINT,'channel    = ',channel
    PRINT,'device     = ',device
    PRINT,'camera     = ',camera
;    PRINT,'file	      = ',file
    PRINT,'date	      = ',date
    PRINT,'line	      = ',line
    PRINT,'shot	      = ',shot
    PRINT,'channel    = ',channel
    PRINT,'frame      = ',frame
    PRINT,'path       = ',path
    PRINT,'calibrate  = ',calibrate
;    PRINT,'shift      = ',shift
;    PRINT,'radius     = ',radius
;    PRINT,'filter     = ',filter
;    PRINT,'window     = ',window
    PRINT,'clean      = ',clean
    PRINT,'background = ',background
    PRINT,'plots      = ',plots     
  ENDIF		     
;		     
; Initialise colour scheme:
  IF (KEYWORD_SET(plots)) THEN LOADCT,4
;
; ----------------------------------------------------------------------
; Load the image data:
; ----------------------------------------------------------------------
;
  image = LoadImage(device,camera,file,date,line,shot,channel,             $
                    frame,path,calibrate,shift,radius,filter,window,  $
                    clean,background,plots)     ; Make some of these optional...
;
; ----------------------------------------------------------------------
; Manipulate the image:
; ----------------------------------------------------------------------
;
; ----------------------------------------------------------------------
; 
  AdjustImage, image, brighten, flip
;
; ----------------------------------------------------------------------
;
; 
;
;
  IF (centre) THEN BEGIN                       
    IF (KEYWORD_SET(debug)) THEN BEGIN
      PRINT,'radius = ',radius
      PRINT,'shift  = ',shift
      PRINT,'depth  = ',depth
    ENDIF

    FindImageCentre,image,radius,shift,depth,clean

    RETURN, image
  ENDIF
;
; ----------------------------------------------------------------------
;
  IF (calibrate) THEN BEGIN                       

    GetImageMaps, image

;    IF (clean) THEN CleanImage,image_cal

    cal = TraceRadii(image)

    print,cal.angle[WHERE(cal.max_r NE 0)]
    print,cal.max_r[WHERE(cal.max_r NE 0)]

    image.vignette = 1.0

    FOR i1 = 0, cal.n-1 DO BEGIN

      angle1 = cal.angle[i1]
      IF (i1 LT cal.n-1) THEN BEGIN
        i2 = i1 + 1
        angle2 = cal.angle[i2]
      ENDIF ELSE BEGIN
        i2 = 0
        angle2 = cal.angle[i2] + 360.0
      ENDELSE

      j = WHERE(image.map_theta GE angle1 AND  $
                image.map_theta LE angle2 AND  $
                image.map_r LE image.radius)

      frac = (image.map_theta[j] - angle1) / (angle2 - angle1)

      print,angle1,angle2,N_ELEMENTS(j),MIN(frac),MAX(frac)
;      print,i1,cal.A[i1,0:cal.ndeg]
;      print,i2,cal.A[i2,0:cal.ndeg]

      cval1 = POLY(image.map_r[j],cal.A[i1,0:cal.ndeg])
      cval2 = POLY(image.map_r[j],cal.A[i2,0:cal.ndeg])

      sf = cval1 + frac * (cval2 - cval1) 

      i = WHERE(sf LT 1.0E-5)
      IF (i[0] NE -1) THEN sf[i] = 1.0

      image.vignette[j] = 1.0 / sf

      image.cal.mode = -1

    ENDFOR

    RemoveArtificialBlack, image

    image_data = image.raw * image.vignette

;   Make sure that noise near the edge of the image doesn't spoil the
;   plot:
    saturation = 2 ^ image.depth - 2    
    i = WHERE(image_data GT saturation)
    IF (i[0] NE -1) THEN image_data[i] = saturation
  ENDIF 
;
; ----------------------------------------------------------------------
;
  IF (KEYWORD_SET(raw)) THEN BEGIN
    RemoveArtificialBlack, image
  ENDIF
;
; ----------------------------------------------------------------------
;
  image_data = image.raw 
;
; ----------------------------------------------------------------------
;
  IF (KEYWORD_SET(rotate)) THEN BEGIN
    rotate_code = 0
    IF (rotate EQ 90 ) THEN rotate_code = 1
    IF (rotate EQ 180) THEN rotate_code = 2
    IF (rotate EQ 270) THEN rotate_code = 3
    image_data = ROTATE(image_data,rotate_code)
  ENDIF
;

; *** Need to expand the image out to the same size as the calibration image...

;
; ----------------------------------------------------------------------
;
  IF (NOT calibrate AND NOT KEYWORD_SET(raw)) THEN BEGIN

    CASE image.camera OF
;     ------------------------------------------------------------------    
      'unknown': BEGIN
        image_data = image.raw
        GetImageMaps, image
        IF (clean) THEN image_data = CleanImage(image_data)
        END
;     ------------------------------------------------------------------    
      'LINCAM': BEGIN
        GetImageMaps, image
        END
;     ------------------------------------------------------------------    
      'ZEBRA': BEGIN
        GetImageMaps, image
        IF (clean) THEN image_data = CleanImage(image_data)
        END
;     ------------------------------------------------------------------    
      'PHOTRON': BEGIN
        image_data = image.raw
        GetImageMaps, image
        IF (clean) THEN image_data = CleanImage(image_data)
;        image = CREATE_STRUCT(image,'data',image_data)
        dim = SIZE(image_data,/DIMENSIONS)
        image_data[dim[0]/2,dim[1]-1] = 2^image.depth-1
        END
;     ------------------------------------------------------------------    
      'RGB': BEGIN
        GetImageMaps, image
        IF (clean) THEN image_data = CleanImage(image_data)
        END
;     ------------------------------------------------------------------    
      'MWIR': BEGIN
        GetImageMaps, image
        IF (clean) THEN image_data = CleanImage(image_data)
        image_min = MIN([image_data])
        image_max = MAX([image_data])

        max_val = MAX(image_data)
        bottom = 0.135 * max_val
        top    = 0.185 * max_val
;        bottom = 0.50 * max_val  ; for 19374, 68
;        top    = 0.55 * max_val
        print,'BOTTOM:',bottom
        FOR ix = 0, image.xdim-1 DO BEGIN
          FOR iy = 0, image.ydim-1 DO BEGIN
            IF (image_data[ix,iy] LT bottom ) THEN image_data[ix,iy] = bottom
            IF (image_data[ix,iy] GT top    ) THEN image_data[ix,iy] = top
          ENDFOR
        ENDFOR
        END
;     ------------------------------------------------------------------    
      'LWIR': BEGIN
        GetImageMaps, image
        IF (clean) THEN image_data = CleanImage(image_data)
        print,'min,max:',MIN(image_data),MAX(image_data)
        max_val = 2^image.depth - 1 ; MAX(image_data)
        bottom = 0.999 * max_val
        top    = 1.00 * max_val
;        bottom = 0.00 * max_val ; for 19374, 42
;        top    = 1.00 * max_val
        print,'BOTTOM,TOP:',bottom,top
        FOR ix = 0, image.xdim-1 DO BEGIN
          FOR iy = 0, image.ydim-1 DO BEGIN
            IF (image_data[ix,iy] LT bottom) THEN image_data[ix,iy] = bottom
            IF (image_data[ix,iy] GT top   ) THEN image_data[ix,iy] = top
          ENDFOR
        ENDFOR
        print,'min,max:',MIN(image_data),MAX(image_data)
        END
;     ------------------------------------------------------------------    
      'FFC': BEGIN
;        image_data = image.raw
;        IF (KEYWORD_SET(rotate)) THEN BEGIN
;          rotate_code = 0
;          IF (rotate EQ 90 ) THEN rotate_code = 1
;          IF (rotate EQ 180) THEN rotate_code = 2
;          IF (rotate EQ 270) THEN rotate_code = 3
;          image_data = ROTATE(image_data,rotate_code)
;        ENDIF
        GetImageMaps, image
        IF (clean) THEN image_data = CleanImage(image_data)
;        image = CREATE_STRUCT(image,'data',image_data)
;        dim = SIZE(image_data,/DIMENSIONS)
;        image.xdim = dim[0]
;        image.ydim = dim[1]
;        image_data[dim[0]/2,dim[1]-1] = 2^image.depth-1
        END
;     ------------------------------------------------------------------    
      'DIVCAM': BEGIN
        IF (KEYWORD_SET(nocal)) THEN BEGIN
;          image = CREATE_STRUCT(image,'data',image.raw)
;          image_data = image.raw
        ENDIF ELSE BEGIN
          LoadCalibrationData, image
          image_data = image_data * image.vignette
;          image = CREATE_STRUCT(image,'data',image_data)
        ENDELSE

        RemoveArtificialBlack, image

;       Clean:
        IF (clean) THEN image_data = CleanImage(image_data)

        gain_ref = 10.0

        print,'cal mode:',image.cal.mode

        CASE image.cal.mode OF      
          1: BEGIN  ; camera + filter + cube (full system)

            cal_gf = 10.0 ^ ((image.cal.gain - gain_ref) / 20.0)
            cf = image.cal.power * image.cal.shutter * cal_gf / image.cal.counts

            print,'window tr:',image.cal.window
            print,'filter nd:',image.filter.nd
            print,'fudge  tr:',image.cal.fudge

            losses = image.cal.window * image.filter.nd * image.cal.fudge
            gf = 10.0 ^ ((image.gain[0] - gain_ref) / 20.0)
            scale_factor = (4.0 * !PI * cf * 0.001) /  $
                           (losses * image.shutter * gf)

            print,image.cal.power,cal_gf,image.cal.shutter,image.cal.counts
            print,cf
            print,image.cal.window,image.filter.nd,image.cal.cube,image.cal.fudge
            print,image.shutter,gf
            print,scale_factor

            scale_factor = scale_factor * image.line * 1.0E-09 /  $
                           (6.63E-34 * 3.0E+08)

            print,'SCALE FACTOR:',scale_factor
;           Scale the image:
            image_data = image_data * scale_factor
            END

          3: BEGIN  ; camera + filter

            cal_gf = 10.0 ^ ((image.cal.gain - gain_ref) / 20.0)
            cf = image.cal.power * image.cal.shutter * cal_gf / image.cal.counts

            print,'window tr:',image.cal.window
            print,'filter nd:',image.filter.nd
            print,'cube   tr:',image.cal.cube
            print,'fudge  tr:',image.cal.fudge

            losses = image.cal.window * image.filter.nd * image.cal.cube * image.cal.fudge
            gf = 10.0 ^ ((image.gain[0] - gain_ref) / 20.0)
            scale_factor = (4.0 * !PI * cf * 0.001) /  $
                           (losses * image.shutter * gf)

            print,image.cal.power,cal_gf,image.cal.shutter,image.cal.counts
            print,cf
            print,image.cal.window,image.filter.nd,image.cal.cube,image.cal.fudge
            print,image.shutter,gf
            print,scale_factor

            scale_factor = scale_factor * image.line * 1.0E-09 /  $
                           (6.63E-34 * 3.0E+08)

            print,scale_factor

;           Scale the image:
            image_data = image_data * scale_factor
            END

          -1:  ; none

          ELSE: BEGIN
            PRINT,'Error ...: Unrecognized calibration mode'
            PRINT,'  MODE = ',image.cal.mode
            END
        ENDCASE
        END
;     ------------------------------------------------------------------    
      'CMOD_CAMERA01': BEGIN
        IF (NOT KEYWORD_SET(maskradius)) THEN maskradius = 275
        image_data = image.raw 
        GetImageMaps, image

        i = WHERE(image.map_r GT 168 AND image.map_r LE maskradius )
;        a = 0.00012645*image.map_r[i]^3-0.0925*image.map_r[i]^2+20.325*image.map_r[i]-1193.5
        a = 0.0001045*image.map_r[i]^3-0.0742*image.map_r[i]^2+15.727*image.map_r[i]-834.67
        a = a / 210.0
        image_data[i] = image_data[i] / a

        IF (clean) THEN image_data = CleanImage(image_data)
        END
    ENDCASE

    ; Filter blue-shift correction:
    IF (image.camera EQ 'DIVCAM' AND image.filter.tag NE 'empty' AND  $
        NOT KEYWORD_SET(nocal)) THEN BEGIN
      xdat = image.filter.nr * image.radius
;      print,'MAX XDAT:  ',max(xdat)
;      print,'IMAGE.RAD: ',image.radius
;      print,'IMA MAP:' ,N_ELEMENTS(WHERE(image.map_r GE image.radius))
;      print,image.map_r

      ydat = image.filter.tr

      image.blueshift = INTERPOL(ydat,xdat,image.map_r)
      image.blueshift[WHERE(image.map_r GE image.radius)] = 1.0

      image_data = image_data / image.blueshift
    ENDIF

  ENDIF
;
;
;
;
;
;
  image = CREATE_STRUCT(image,'data',image_data)
  dim = SIZE(image_data,/DIMENSIONS)
  image.xdim = dim[0]
  image.ydim = dim[1]
;
; Hide parts of the image by setting the pixel values to zero.  All pixel
; values that were naturally zero are set to one instead:
;
  MaskImage, image, maskradius, maskleft, maskright, masktop,  $
                    maskbottom, maskbox , maskpoly
;
; Some debugging output:
;
  IF (KEYWORD_SET(debug)) THEN BEGIN
    PRINT,'device  = ',image.device
    PRINT,'date    = ',image.date
    PRINT,'window  = ',image.window
    PRINT,'line    = ',image.line
    PRINT,'shot    = ',image.shot
    PRINT,'frame   = ',image.frame
    PRINT,'time    = ',image.time
    PRINT,'channel = ',image.channel
    PRINT,'radius  = ',image.radius
    PRINT,'xcen    = ',image.xcen
    PRINT,'ycen    = ',image.ycen
  ENDIF
;
; ----------------------------------------------------------------------
;
  IF (KEYWORD_SET(plots)) THEN BEGIN
    IF (KEYWORD_SET(colour)) THEN BEGIN
      DEVICE, DECOMPOSED=0
      IF (colour EQ 1) THEN colour = 4
      LOADCT,colour
      IF (calibrate) THEN safe_colors,/first
    ENDIF ELSE BEGIN
      loadct,0
    ENDELSE

    dim = SIZE(image.raw,/DIMENSIONS)
    max_val = MAX(image.raw)
    image_raw = image.raw * scale
    FOR ix = 0, dim[0]-1 DO BEGIN
      FOR iy = 0, dim[1]-1 DO BEGIN
        IF (image_raw[ix,iy] GT max_val) THEN image_raw[ix,iy] = max_val
      ENDFOR
    ENDFOR

    IF (image.ydim EQ 1) THEN BEGIN
      window,0
      safe_colors,/first
      plot,image.data
    ENDIF ELSE BEGIN
      WINDOW,0,xsize=image.xdim,ysize=image.ydim,retain=2      
      TVSCL,image_raw,/ORDER
    ENDELSE

    max_val = MAX(image.data)
    image_data = image.data * scale
;    i = WHERE(image_data GT max_val)
;    IF (i NE -1) THEN image_rdata[i] = max_val
    FOR ix = 0, image.xdim-1 DO BEGIN
      FOR iy = 0, image.ydim-1 DO BEGIN
        IF (image_data[ix,iy] GT max_val) THEN image_data[ix,iy] = max_val
      ENDFOR
    ENDFOR

    IF (image.ydim EQ 1) THEN BEGIN
    ENDIF ELSE BEGIN
      window,1,xsize=image.xdim,ysize=image.ydim,retain=2
      TVSCL,image_data,/ORDER
    ENDELSE

;    safe_colors,/first
;    window,2,retain=2
;    plot,image.data[500,*]
;    oplot,image.raw[500,*],color=2
  ENDIF
;
; ----------------------------------------------------------------------
;
  IF (KEYWORD_SET(analyse)) THEN AnalyseImage,image.data
;
; ----------------------------------------------------------------------
;
  IF (KEYWORD_SET(save) OR KEYWORD_SET(save_png)) THEN   $
    SaveImageData, image, calibrate, path, scale, sname, save_png
;
; ----------------------------------------------------------------------
;
  RETURN, image

END

