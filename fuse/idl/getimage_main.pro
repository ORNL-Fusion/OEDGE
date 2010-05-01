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
   reg        = reg      ,   $  ; registration image shot number
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
   full       = full,        $
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

  IF (NOT KEYWORD_SET(colour)) THEN colour = 5

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
    file = 'none'
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
    'FFC': BEGIN
      CASE (channel) OF
       -1: channel_tag = 'unknown'
        1: channel_tag = 'a'
        2: channel_tag = 'b'
        3: channel_tag = 'c'
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
  IF (KEYWORD_SET(plots) AND KEYWORD_SET(colour)) THEN BEGIN
    DEVICE, DECOMPOSED=0
    IF (colour EQ 1) THEN LOADCT, 3 ELSE LOADCT, colour
  ENDIF
;
; ----------------------------------------------------------------------
; Load the image data:
; ----------------------------------------------------------------------
;
  image = LoadImage(device,camera,file,date,line,shot,channel,             $
                    frame,path,calibrate,reg,shift,radius,filter,window,  $
                    clean,background,plots)     ; Make some of these optional...


  IF (KEYWORD_SET(debug)) THEN BEGIN
    PRINT,'filter     = ',filter
    PRINT,'window     = ',window
  ENDIF		     
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
  IF (centre NE 0) THEN BEGIN                       
    IF (KEYWORD_SET(debug)) THEN BEGIN
      PRINT,'radius = ',radius
      PRINT,'shift  = ',shift
      PRINT,'depth  = ',depth
    ENDIF
    FindImageCentre,image,centre,radius,shift,depth,clean,background,rotate
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
                image.map_r LE image.radius, count)
      frac = (image.map_theta[j] - angle1) / (angle2 - angle1)

      print,angle1,angle2,count,MIN(frac),MAX(frac)
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
;    IF (rotate EQ 90 ) THEN rotate_code = 1
;    IF (rotate EQ 180) THEN rotate_code = 2
;    IF (rotate EQ 270) THEN rotate_code = 3
;    top  = image.ywin[0]
;    left = image.xwin[0]
;    width  = image.xwin[1] - image.xwin[0] + 1
;    height = image.ywin[1] - image.ywin[0] + 1
    xwin = image.xwin
    ywin = image.ywin
    CASE rotate OF
       90: BEGIN
        rotate_code = 1 
        image.xwin[0] = ywin[1] ; bottom becomes left
        image.xwin[1] = ywin[0] ; top becomes right
        image.ywin[0] = xwin[0] ; left becomes top
        image.ywin[1] = xwin[1] ; right become bottom
        END
      180: BEGIN
        rotate_code = 2 
        image.xwin[0] = xwin[1] ; right becomes left  
        image.xwin[1] = xwin[0] ; left become right
        image.ywin[0] = ywin[1] ; bottom to top
        image.ywin[1] = ywin[0] ; top to bottom
       END
      270: BEGIN
        rotate_code = 3 
        image.xwin[0] = ywin[0] ; top becomes left
        image.xwin[1] = ywin[1] ; bottom becomes right
        image.ywin[0] = image.xdim - xwin[1] + 1 ; right becomes top
        image.ywin[1] = image.xdim - xwin[0] + 1 ; left becomes bottom
       END
    ENDCASE

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
;        GetImageMaps, image
;        IF (clean) THEN image_data = CleanImage(image_data)
;        image = CREATE_STRUCT(image,'data',image_data)
;        dim = SIZE(image_data,/DIMENSIONS)
;        image.xdim = dim[0]
;        image.ydim = dim[1]
;        image_data[dim[0]/2,dim[1]-1] = 2^image.depth-1


        IF (KEYWORD_SET(nocal)) THEN BEGIN
;          image = CREATE_STRUCT(image,'data',image.raw)
;          image_data = image.raw
        ENDIF ELSE BEGIN
          LoadCalibrationData, image
          IF (KEYWORD_SET(shift)) THEN BEGIN
            image_vignette  = image.vignette
            image_map_r     = image.map_r
            image_map_theta = image.map_theta
            ShiftImage,shift[0],shift[1],image_vignette
            ShiftImage,shift[0],shift[1],image_map_r
            ShiftImage,shift[0],shift[1],image_map_theta
            i = WHERE(image_vignette  EQ -1.0, count_i)
            j = WHERE(image_map_r     EQ -1.0, count_j)
            k = WHERE(image_map_theta EQ -1.0, count_k)
            IF (count_i GT 0) THEN image_vignette [i] = 0.0
            IF (count_j GT 0) THEN image_map_r    [j] = 0.0
            IF (count_k GT 0) THEN image_map_theta[k] = 0.0
            image.vignette  = image_vignette 
            image.map_r     = image_map_r    
            image.map_theta = image_map_theta
          ENDIF
          image_data = image_data * image.vignette
          image_data[WHERE(image.map_r GT image.radius)] = -1  ; Artifical black outside calibrated region
        ENDELSE

        RemoveArtificialBlack, image

;       Clean:
        IF (clean) THEN image_data = CleanImage(image_data)

        CASE image.cal.mode OF      
         -1:  ; none
          1: BEGIN  ; camera + filter + cube (full system)
            PRINT,'Absolute FFC calibration not applied'
            END
          ELSE: BEGIN
            PRINT,'ERROR GetImage: Unrecognized FFC calibration mode'
            PRINT,'  MODE = ',image.cal.mode
            END
        ENDCASE

        END
;     ------------------------------------------------------------------    
      'DIVCAM': BEGIN
        IF (KEYWORD_SET(nocal)) THEN BEGIN
;          image = CREATE_STRUCT(image,'data',image.raw)
;          image_data = image.raw
        ENDIF ELSE BEGIN
          LoadCalibrationData, image
          IF (KEYWORD_SET(shift)) THEN BEGIN
            image_vignette  = image.vignette
            image_map_r     = image.map_r
            image_map_theta = image.map_theta
            ShiftImage,shift[0],shift[1],image_vignette
            ShiftImage,shift[0],shift[1],image_map_r
            ShiftImage,shift[0],shift[1],image_map_theta
            i = WHERE(image_vignette  EQ -1.0, count_i)
            j = WHERE(image_map_r     EQ -1.0, count_j)
            k = WHERE(image_map_theta EQ -1.0, count_k)
            IF (count_i GT 0) THEN image_vignette [i] = 0.0
            IF (count_j GT 0) THEN image_map_r    [j] = 0.0
            IF (count_k GT 0) THEN image_map_theta[k] = 0.0
            image.vignette  = image_vignette 
            image.map_r     = image_map_r    
            image.map_theta = image_map_theta
          ENDIF
          image_data = image_data * image.vignette
          image_data[WHERE(image.map_r GT image.radius)] = -1  ; Artifical black outside calibrated region
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
    IF ((image.camera EQ 'DIVCAM' OR image.camera EQ 'FFC') AND  $
        image.filter.tag NE 'empty' AND  $
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

      PRINT,'Blueshift filter correciton applied'
    ENDIF

  ENDIF
;
;
;
;
  size_adjust = 1.0

  IF (NOT KEYWORD_SET(full)) THEN BEGIN
    CASE camera OF
      'FFC': BEGIN
        print,image.xwin
        print,image.ywin

        top    = image.ywin[0]
        left   = image.xwin[0]
        width  = image.xwin[1] - image.xwin[0] + 1
        height = image.ywin[1] - image.ywin[0] + 1
        print,top,left,width,height
        
        image_new = MAKE_ARRAY(width,height,/INTEGER,VALUE=0)
        print,width-1
        print,left-1
        print,left-1+width-1
        print,left-1+width-1
        print,top-1
        FOR iy = 0, height-1 DO BEGIN
          image_new[0:width-1,iy] = image_data[left-1:left-1+width-1,iy+top-1]
        ENDFOR
        image_data = image_new
        END
      'DIVCAM': BEGIN
print,image.xbin,image.ybin
        IF (image.xbin EQ 2 AND image.ybin EQ 2) THEN BEGIN
          image_new = MAKE_ARRAY(500,500,/LONG,VALUE=0L)
          FOR i = 0, 499 DO BEGIN
            FOR j = 0, 499 DO BEGIN
              image_new[i,j] = image_data[2*i,2*j]  ; No need to rescale signal strength since the calibration (should have) been
            ENDFOR                                  ; done without binning...
          ENDFOR 
        ENDIF
        image_data = image_new
        size_adjust = 0.5
        END
      ELSE:
    ENDCASE
  ENDIF
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
    PRINT,'xdim    = ',image.xdim
    PRINT,'ydim    = ',image.ydim
  ENDIF
;
; ----------------------------------------------------------------------
;
  IF (KEYWORD_SET(plots)) THEN BEGIN

    extension = STRMID(image.file,2,/REVERSE_OFFSET) 
    CASE extension OF
      'bmp': order = 0
      'ipx': order = 1
      ELSE: BEGIN
        PRINT,'ERROR GetImage: Unknown file type'
        RETURN, -1
        END
    ENDCASE

    IF (KEYWORD_SET(colour)) THEN BEGIN
      DEVICE, DECOMPOSED=0
      IF (colour EQ 1) THEN colour = 3
      LOADCT,colour
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
      limit = (2 ^ image.depth - 2)
      IF (KEYWORD_SET(shift)) THEN BEGIN
        angle = FINDGEN(180)*2.0       ; Dotted circle showing where the image boundary is assumed to be
        reference_radius = image.radius - 30
        ix = image.xcen + image.radius * SIN (angle*!PI/180.0) 
        iy = image.ycen + image.radius * COS (angle*!PI/180.0) 
        image_raw[ix,iy] = max_val
        angle = FINDGEN(3*360) / 3.0   ; Finer circle showing where the vignette mask is being applied
        ix = shift[0] + image.xcen + image.radius * SIN (angle*!PI/180.0) 
        iy = shift[1] + image.ycen + image.radius * COS (angle*!PI/180.0) 
        image_raw[ix,iy] = max_val
      ENDIF
      WINDOW,0,xsize=dim[0],ysize=dim[1],retain=2      
      TVSCL,image_raw,ORDER=order
    ENDELSE

    max_val = MAX(image.data)
    image_data = image.data * scale
;    max_val = (2 ^ image.depth - 2)
    i = WHERE(image_data GT max_val,count)
    IF (count GT 0) THEN image_data[i] = max_val

;    FOR ix = 0, image.xdim-1 DO BEGIN
;      FOR iy = 0, image.ydim-1 DO BEGIN
;        IF (image_data[ix,iy] GT max_val) THEN image_data[ix,iy] = max_val
;      ENDFOR
;    ENDFOR

    IF (image.ydim EQ 1) THEN BEGIN
    ENDIF ELSE BEGIN
      IF (KEYWORD_SET(shift) AND (camera NE 'FFC' OR KEYWORD_SET(full))) THEN BEGIN
        angle = FINDGEN(180)*2.0       ; Dotted circle showing where the image boundary is assumed to be
        reference_radius = image.radius - 30
        ix = image.xcen * size_adjust + image.radius * size_adjust * SIN (angle*!PI/180.0) 
        iy = image.ycen * size_adjust + image.radius * size_adjust * COS (angle*!PI/180.0) 
        image_data[ix,iy] = max_val
        angle = FINDGEN(3*360) / 3.0   ; Finer circle showing where the vignette mask is being applied
        ix = shift[0] * size_adjust + image.xcen * size_adjust + image.radius * size_adjust * SIN (angle*!PI/180.0) 
        iy = shift[1] * size_adjust + image.ycen * size_adjust + image.radius * size_adjust * COS (angle*!PI/180.0) 
        image_data[ix,iy] = max_val
      ENDIF
      window,1,xsize=image.xdim,ysize=image.ydim,retain=2
      TVSCL,image_data,ORDER=order
    ENDELSE

    IF (KEYWORD_SET(calibrate)) THEN BEGIN
      WINDOW,3,XSIZE=image.xdim,YSIZE=image.ydim,RETAIN=2
      TVSCL,image.vignette,ORDER=order
      WINDOW,4,XSIZE=image.xdim,YSIZE=image.ydim,RETAIN=2
      TVSCL,image_data*image.vignette,ORDER=order
    ENDIF

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
    SaveImageData, image, path, scale, sname, save_png, calibrate=calibrate, order=order
;
; ----------------------------------------------------------------------
;
  RETURN, image

END

