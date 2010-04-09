;
; ========================================================================
;
PRO FindImageCentre, image, centre, radius, shift, depth, clean, background, rotate
;PRO FindImageCentre, image, radius=radius, shift=shift, depth=depth
  
  image_centre = image.raw

  extension = STRMID(image.file,2,/REVERSE_OFFSET) 
  PRINT,'FILE EXTENSION = ',extension
  CASE extension OF
    'bmp': order = 0
    'ipx': order = 1
    ELSE: BEGIN
      PRINT,'ERROR getimage_FindImageCentre: Unknown file type'
      RETURN
      END
  ENDCASE

  PRINT,'ORDER=',order

  IF (KEYWORD_SET(rotate)) THEN BEGIN
    rotate_code = 0
    IF (rotate EQ 90 ) THEN rotate_code = 1
    IF (rotate EQ 180) THEN rotate_code = 2
    IF (rotate EQ 270) THEN rotate_code = 3
    image_centre = ROTATE(image_centre,rotate_code)
  ENDIF

; Subtrack non-zero background by sampling around the edge of the image:
  IF (KEYWORD_SET(background)) THEN BEGIN
    IF (background EQ 1) THEN BEGIN
      xwidth = N_ELEMENTS(image.raw[*,0])
      ywidth = N_ELEMENTS(image.raw[0,*])
      idelta = xwidth / 10
      jdelta = ywidth / 10

      idata =         TOTAL(image.raw[0              :idelta ,0       ])
      idata = idata + TOTAL(image.raw[xwidth-idelta-1:xwidth-1,0       ])
      idata = idata + TOTAL(image.raw[0              :idelta ,ywidth-1])
      idata = idata + TOTAL(image.raw[xwidth-idelta-1:xwidth-1,ywidth-1])
      idata = idata / (4.0 * FLOAT(idelta+1))
      
      jdata =         TOTAL(image.raw[0       ,0              :jdelta  ])
      jdata = jdata + TOTAL(image.raw[xwidth-1,0              :jdelta  ])
      jdata = jdata + TOTAL(image.raw[0       ,ywidth-jdelta-1:ywidth-1])
      jdata = jdata + TOTAL(image.raw[xwidth-1,ywidth-jdelta-1:ywidth-1])
      jdata = jdata / (4.0 * FLOAT(jdelta+1))
      
      PRINT,'X,Y BACKGROUND LEVELS = ',idata,jdata
      PRINT,'X,YWIDTH=',xwidth,ywidth
     
      image_centre = image_centre - (idata + jdata) / 2
    ENDIF ELSE BEGIN
      image_centre = image_centre - background
    ENDELSE
    image_centre[WHERE(image_centre LT 0)] = 0     
  ENDIF

  IF (KEYWORD_SET(radius)) THEN image.radius = radius
  IF (KEYWORD_SET(shift) ) THEN ShiftImage,-shift[0],-shift[1],image_centre
  IF (KEYWORD_SET(depth) ) THEN image.depth = depth

  RemoveArtificialBlack, image

  IF (clean) THEN image_centre = CleanImage(image_centre)

  
  IF (FLOAT(centre) NE 1.0) THEN cut_off = centre ELSE cut_off = 0.10

  limit = (2 ^ image.depth - 2) * cut_off

  PRINT,'CUT OFF=',cut_off
  PRINT,'LIMIT=',limit

  i = WHERE(image_centre GT limit, count)
  IF (count GT 0) THEN image_centre[WHERE(image_centre GT limit)] = limit ELSE BEGIN  $
    PRINT, 'ERROR getimage_FindImageCentre: No data above cut off level'
    STOP
  ENDELSE

  angle = FINDGEN(360) 
  reference_radius = image.radius - 30
  ix = image.xdim / 2 + reference_radius * SIN (angle*!PI/180.0) 
  iy = image.ydim / 2 + reference_radius * COS (angle*!PI/180.0) 
  image_centre[ix,iy] = 1.5 * limit

  FOR rad = image.radius, image.radius+30 DO BEGIN
    ix = image.xdim / 2 + rad * SIN (angle*!PI/180.0) 
    iy = image.ydim / 2 + rad * COS (angle*!PI/180.0) 
    FOR i = 0, 359 DO ix[i] = MAX([0.0,MIN([image.xdim,ix[i]])])
    FOR i = 0, 359 DO iy[i] = MAX([0.0,MIN([image.ydim,iy[i]])])
    image_centre[ix,iy] = 1.5 * limit
  ENDFOR

  print,min(image_centre),max(image_centre)

  WINDOW,0,XSIZE=image.xdim,YSIZE=image.ydim,RETAIN=2
  TVSCL,image_centre,ORDER=order

  image = CREATE_STRUCT(image,'data',image_centre)
END
;
; ========================================================================
;
PRO SetCubeTransmission, image

  line = image.line  ; should perhaps be image.filter.line (or whatever...)? 

  IF (image.cal.cube_n LT 1) THEN BEGIN
    PRINT,'Error SetCubeTransmission: Data not found'
    STOP
  ENDIF

  i = WHERE(image.cal.cube_wlngth NE 0.0)

  IF (image.cal.cube_wlngth[i[0              ]] GT line OR  $
      image.cal.cube_wlngth[i[N_ELEMENTS(i)-1]] LT line) THEN BEGIN
    PRINT,'Error SetCubeTransmission: Wavelength out of interpolation range'
    STOP
  ENDIF 

  CASE image.channel OF
    1: tr = image.cal.cube_straight
    2: tr = image.cal.cube_orthog
    ELSE: BEGIN
      PRINT,'Error SetCubeTransmission: Unknown channel'
      PRINT,'  CHANNEL = ',image.channel
      STOP
      END
  ENDCASE

  a = FINDGEN(400) + 400.0

  image.cal.cube = INTERPOL(tr[i],image.cal.cube_wlngth[i],line)

END
;
; ========================================================================
;
PRO SetWindowTransmission, image

  line = image.line  ; should perhaps be image.filter.line (or whatever...)? 

  IF (image.window EQ 'none') THEN RETURN

  IF (image.cal.w_n EQ -1) THEN BEGIN
    PRINT,'Error SetWindowTransmission: Transmission data not found'
    STOP
  ENDIF

  i = WHERE(image.cal.w_wlngth NE 0.0)

  IF (image.cal.w_wlngth[i[0              ]] GT line OR  $
      image.cal.w_wlngth[i[N_ELEMENTS(i)-1]] LT line) THEN BEGIN
    PRINT,'Error SetWindowTransmission: Wavelength out of bounds'
    STOP
  ENDIF 

  CASE image.cal.w_interpol OF

    1: BEGIN

       image.cal.window = INTERPOL(image.cal.w_trans[i,0],image.cal.w_wlngth[i],line)
       END
 
    2: BEGIN

       image.cal.window = INTERPOL(image.cal.w_trans[i,1],image.cal.w_wlngth[i],line)

;      a = FINDGEN(400) + 400.0
;       result = INTERPOL(image.cal.w_trans[i,1],image.cal.w_wlngth[i],a)
;       window,2,retain=2
;       safe_colors,/first
;       plot,a,result
;       oplot,image.cal.w_wlngth[i],image.cal.w_trans[i,1],psym=6
;       stop

      END

    ELSE: BEGIN
      PRINT,'Error SetWindowTransmission: Unrecognized interpolation scheme'
      PRINT,'  W_INTERPOL = ',image.cal.w_interpol
      END
  ENDCASE

END
;
; ========================================================================
;
FUNCTION GetLine, fp, buffer

  WHILE (1) DO BEGIN

    READF,fp,buffer

    buffer = STRTRIM(buffer,2)

    i = STRPOS(buffer,'$')
    IF (i EQ  0) THEN CONTINUE 
    IF (i NE -1) THEN BEGIN
      str    = STRSPLIT(buffer,'$',/extract)
      buffer = str[0]
    ENDIF

;    print,buffer

    IF (STRMATCH(buffer,'*{end}*') EQ 1) THEN return, 0
    IF (STRMATCH(buffer, '*{*'   ) EQ 1) THEN return, 1
    return, 2
  ENDWHILE

END 
;
; ========================================================================
;
PRO LoadCalibrationData, image

  date    = image.date
  window  = image.window
  shot    = image.shot
  channel = image.channel
  line    = image.line

  path = '~/fuse_data/mast/camera_calibration/'

  fp = 2
  FREE_LUN,fp
  fname = path+'divcam_cal.txt'
  OPENR,fp,fname,error=error
  IF (error NE 0) THEN BEGIN
    print,'Troubles mate, unable to find ',fname
    STOP
  ENDIF  

  buffer = ' '

  cal_actual = PTR_NEW(GetImageCalibrationStructure())
  cal_dummy  = PTR_NEW(GetImageCalibrationStructure())

  cont = 1  
  WHILE (cont GE 1) DO BEGIN
    IF (cont EQ 1) THEN status = GetLine(fp, buffer)

    i1 = STRPOS(buffer,'{')
    i2 = STRPOS(buffer,'}')
    tag = STRMID(buffer,i1+1,i2-i1-1)

;    print,'TAG:',tag

    cal = cal_dummy

    CASE tag OF
       
      'vignetting': BEGIN
         WHILE (GetLine(fp,buffer) EQ 2) DO BEGIN
           str = STRSPLIT(buffer,/extract)           
           IF (shot GE LONG(str[0]) AND shot LE LONG(str[1])) THEN cal = cal_actual
         ENDWHILE
         
         WHILE (GetLine(fp,buffer) EQ 2) DO BEGIN
           str = STRSPLIT(buffer,/extract)           
           print,line,channel
           print,str
           IF ((FLOAT(str[0]) EQ line OR FLOAT(str[0]) EQ -1.0) AND  $
               LONG (str[1]) EQ channel ) THEN BEGIN
             (*cal).v_date = LONG(str[2])
             (*cal).v_file = str[3]
           ENDIF
         ENDWHILE
         cont = 2
         END

      'absolute calibration': BEGIN
         WHILE (GetLine(fp,buffer) EQ 2) DO BEGIN
           str = STRSPLIT(buffer,/extract)           
           IF (shot GE LONG(str[0]) AND shot LE LONG(str[1])) THEN cal = cal_actual
         ENDWHILE

         WHILE (GetLine(fp,buffer) EQ 2) DO BEGIN
           str = STRSPLIT(buffer,/extract)           
           IF (LONG (str[0]) EQ channel AND  $
               FLOAT(str[1]) EQ line    ) THEN BEGIN
             (*cal).serial_no = UINT (str[3])
             (*cal).mode      = UINT (str[4])
             (*cal).power     = FLOAT(str[5])  ;  Watts
             (*cal).counts    = FLOAT(str[6])
             (*cal).gain      = FLOAT(str[7])
             (*cal).shutter   = FLOAT(str[8])
             (*cal).fudge     = FLOAT(str[9])
             (*cal).date      = LONG (str[10])
             (*cal).file      = str[11]
           ENDIF
         ENDWHILE
         cont = 2
         END

      'window transmission': BEGIN
         WHILE (GetLine(fp,buffer) EQ 2) DO BEGIN
           str = STRSPLIT(buffer,/extract)           
           IF (window NE 'none' AND shot GE LONG(str[0]) AND  $
                                    shot LE LONG(str[1])) THEN BEGIN
             cal = cal_actual
             (*cal).w_shot_range[0] = LONG(str[0])
             (*cal).w_shot_range[1] = LONG(str[1])
             (*cal).w_interpol      = LONG(str[2])
           ENDIF
         ENDWHILE

         status = GetLine(fp,buffer)
         str = STRSPLIT(buffer,/extract)           
         i = -1
         FOR j = 0, N_ELEMENTS(str)-1 DO  $
           IF (STRMATCH(window,'*'+str[j]+'*') EQ 1) THEN i = 2 * (j - 1) + 1
 
         IF (i EQ -1) THEN BEGIN
           cal = cal_dummy
           i = 1
           print,'turning off absolute calibration ???'
         ENDIF

         (*cal).window_tag = str[(i-1)/2+1]
         WHILE (GetLine(fp,buffer) EQ 2) DO BEGIN
           str = STRSPLIT(buffer,/extract)                      
           (*cal).w_n = (*cal).w_n + 1
           (*cal).w_wlngth[(*cal).w_n]   = FLOAT(str[0])
           (*cal).w_trans [(*cal).w_n,0] = FLOAT(str[i])
           (*cal).w_trans [(*cal).w_n,1] = FLOAT(str[i+1])
         ENDWHILE
         cont = 2           
         END

      'cube transmission': BEGIN
         WHILE (GetLine(fp,buffer) EQ 2) DO BEGIN
           str = STRSPLIT(buffer,/extract)           
           IF (shot GE LONG(str[0]) AND  $
               shot LE LONG(str[1])) THEN cal = cal_actual
         ENDWHILE

         IF (GetLine(fp,buffer) EQ 2) THEN BEGIN
           str = STRSPLIT(buffer,/extract)           
           scale_factor = FLOAT(str)
         ENDIF ELSE BEGIN
           PRINT,'Error LoadCalibrationData: Bad cube transmission format'
           PRINT,'  BUFFER = ',STRTRIM(buffer)
           STOP
         ENDELSE

         WHILE (GetLine(fp,buffer) EQ 2) DO BEGIN
           str = STRSPLIT(buffer,/extract)                      
           (*cal).cube_n = (*cal).cube_n + 1
           (*cal).cube_wlngth  [(*cal).cube_n] = FLOAT(str[0])
           (*cal).cube_straight[(*cal).cube_n] = FLOAT(str[1]) * scale_factor[0]
           (*cal).cube_orthog  [(*cal).cube_n] = FLOAT(str[2]) * scale_factor[1]
         ENDWHILE
         cont = 2           
         END

      ELSE: cont = 1
    ENDCASE

    IF (STRMATCH(buffer,'*{end}*') EQ 1) THEN BREAK

  ENDWHILE


; Not good enough ... need to do selection based on shot# in {absolute calibration}
; code above:

  image.cal = *cal_actual

  LoadVignetteMap, image, date, line

  IF (image.cal.mode NE -1) THEN BEGIN 

    IF (window NE 'none' AND image.cal.w_shot_range[0] EQ -1) THEN BEGIN
      PRINT,'Error LoadCalibrationData: window transmission data not found'
      STOP
    ENDIF

    SetWindowTransmission,image

    IF (image.cal.mode EQ 3 OR image.cal.mode EQ 4) THEN  $
      SetCubeTransmission, image
  ENDIF

  CLOSE, fp

;  help,image.cal,/struct
;  stop

END
;
; ======================================================================
;
PRO LoadVignetteMap, image, date, line

  shot    = image.shot
  channel = image.channel

  status = -1

;print,image.device
;stop

  CASE image.device OF

    'FILE' : BEGIN

        CASE channel OF
          -1: BEGIN
            image_radius = 490
            image_xshift = 0
            image_yshift = 0 
            status = 0
            END
          1: BEGIN
             image_radius = 470 ; from ...? 
;             image_xshift = +24  ; from 20060731 000100
;             image_yshift =  +0  
;             IF (NOT center) THEN RESTORE,filename=path+'/calibration/20060731_1.cal.sav'
             image_xshift = +12  ; from 20060929 000720
             image_yshift =  +5  ;                     
             IF (NOT center) THEN RESTORE,filename=path+'/calibration/20060929_1.cal.sav'
             END
          2: BEGIN
             image_xshift = +18  ; from ...
             image_yshift = +5  
             IF (NOT center) THEN RESTORE,filename=path+'/calibration/20060731_2.cal.sav'
             END
        ENDCASE

      END

    'MAST' : BEGIN

      IF ( 0 EQ 1 AND shot GE 14993 AND shot LE 15021) THEN BEGIN
;       frame = 17
      
        CASE channel OF
          1: BEGIN
             image_radius = 450
             image_xshift = +27  ; from rda015188.ipx  frame 17 (okay, but not great)
             image_yshift = +11
             IF (NOT center) THEN RESTORE,filename=path+'/calibration/original_1.sav'  ; Loads NFIT, YFIT and AFIT (and overwrites IMAGE_RADIUS at the moment...)
             END
          2: BEGIN
             image_radius = 470 
             image_xshift = +23  ; from 14993 frame 1, but difficult, no dedicated registration image available
             image_yshift = +5   
             IF (NOT center) THEN RESTORE,filename=path+'/calibration/original_2.sav'
             END
         ENDCASE
      
      ENDIF
      
      IF ( 0 EQ 1 AND shot GE 15115 AND shot LE 15189) THEN BEGIN
      
        CASE channel OF
          1: BEGIN
             image_radius = 450
             image_xshift = +27  ; rda015188.ipx, frame 17 (okay, but not great)
             image_yshift = +11
             IF (NOT center) THEN RESTORE,filename=path+'/calibration/original_1.sav'
             END
          2: BEGIN
             image_radius  = 470 
             image_xshift = +23  ; rdb015186.ipx, frame 13
             image_yshift = +15   
             IF (NOT center) THEN RESTORE,filename=path+'/calibration/original_2.sav'
             END
         ENDCASE
      
      ENDIF
;      channel_tag = 'b'
;      shot = '/net/fuslsa/data/MAST_IMAGES/rdb/rdb015116.ipx'  ; HU12
;      shot = '/net/fuslsa/data/MAST_IMAGES/rdb/rdb015163.ipx'  ; HM01 lower x-pt: maskview=390,maskleft=175,maskbox=[370,462,858,999]
;      shot = '/net/fuslsa/data/MAST_IMAGES/rdb/rdb015169.ipx'  ; M01 midplane: maskview=500,maskleft=375,maskbox=[1,553,1,209,1,550,798,999]
;      shot = '/net/fuslsa/data/MAST_IMAGES/rdb/rdb015189.ipx'  ; HL01: maskview=500,masktop=74,maskbox=[613,999,0,122]
;      shot = '/net/fuslsa/data/MAST_IMAGES/rdb/rdb015154.ipx'  ; HM01 lower x-pt: 
;      frame = 13
      
      
      IF ( shot GE 16249 AND shot LE 16480) THEN BEGIN
;      IF ( shot GE 16472 AND shot LE 16480) THEN BEGIN
;       frame = 13
      
        CASE channel OF
          1: BEGIN
             image_radius = 470 
             image_xshift = +12  ; from 20060929 000720
             image_yshift =  +5  ;                     
;             image_xshift = +24  ; from 20060731 000100  +12 ; from rda016472.ipx  frame 2
;             image_yshift =  +0  ;                       -14 ;
             IF (NOT center) THEN RESTORE,filename=path+'/calibration/20060731_1.cal.sav'
             END
          2: BEGIN
             image_radius = 480 
;             image_xshift = +23  ; from 20060929 000731  
;             image_yshift =  +3  ;                    
             image_xshift = +18  ; from 20060731 000101  +23 ; from rdb016472.ipx  frame 2
             image_yshift =  +5  ;                        -5 ; 0
             IF (NOT center) THEN RESTORE,filename=path+'/calibration/20060731_2.cal.sav'
             END
         ENDCASE
      ENDIF

      IF ( (shot GE 13943 AND shot LE 13950) OR  $
           (shot GE 15021 AND shot LE 15021) OR  $
           (shot GE 15115 AND shot LE 15186) OR  $
           (shot GE 16746 AND shot LE 23000) ) THEN BEGIN  ; M6B
        path_cal = './calibration/'

        CASE channel OF
          1: BEGIN

             file = 'cal_'+STRING(image.cal.v_date,format='(I8)')+'_'+image.cal.v_file+'.sav'

;             IF (line EQ 434.0) THEN file = '20061215_7281_1.cal.sav'  ; 728 image
;             IF (line EQ 728.1) THEN file = 'cal_20061219_000208.ipx.sav'  ; from 728 image
;             IF (line EQ  -1.0) THEN file = '20061215_7281_1.cal.sav'  ; from 728 image
               
           
;             image_radius = 500
;             image_xshift = +25  ; from 20061215 000211_728
;             image_yshift =  +5  ;                     
;            IF (NOT center) THEN RESTORE,filename=path+'/calibration/20060731_1.cal.sav'
             END
          2: BEGIN

             file = 'cal_'+STRING(image.cal.v_date,format='(I8)')+'_'+image.cal.v_file+'.sav'

;             IF (line EQ 656.3) THEN file = '20060929_2.cal.sav'  ; from ???

;             image_radius = 480 
;             image_xshift = +23  ; from 20060929 000731  
;             image_yshift =  +3  ;                    
;             image_xshift = +18  ; from 20060731 000101  +23 ; from rdb016472.ipx  frame 2
;             image_yshift =  +5  ;                        -5 ; 0
;             IF (NOT center) THEN RESTORE,filename=path+'/calibration/20060731_2.cal.sav'
             END
         ENDCASE
 
         PRINT,'Vignette file'
         PRINT,file
         RESTORE,filename=path_cal+file,/verbose  ; from 728 image

         IF (KEYWORD_SET(cal_xcen)) THEN BEGIN
           status = 1
           image.xcen      = cal_xcen
           image.ycen      = cal_ycen
           image.radius    = cal_radius
           image.vignette  = cal_vignette
           image.map_r     = cal_map_r
           image.map_theta = cal_map_theta
           image.raw[WHERE(image.map_r GT image.radius)] = -1  ; Artifical black outside calibrated region
         ENDIF 

      ENDIF

      IF (shot GE 24800 AND shot LE 24869) THEN BEGIN  ; Detachment

        path_cal = '~/fuse_data/mast/camera_calibration/'

        CASE channel OF
          1: BEGIN
            file = 'cal_'+STRING(image.cal.v_date,format='(I8)')+'_'+image.cal.v_file+'.sav'
            END
          2: BEGIN
            file = 'cal_'+STRING(image.cal.v_date,format='(I8)')+'_'+image.cal.v_file+'.sav'
            END
          3: BEGIN
            file = 'cal_'+STRING(image.cal.v_date,format='(I8)')+'_'+image.cal.v_file+'.sav'
            END
          4: BEGIN
            file = 'cal_'+STRING(image.cal.v_date,format='(I8)')+'_'+image.cal.v_file+'.sav'
            END
         ENDCASE
 
         PRINT,'Vignette file'
         PRINT,file
         RESTORE,filename=path_cal+file,/verbose  

         IF (KEYWORD_SET(cal_xcen)) THEN BEGIN
           status = 1
           image.xcen      = cal_xcen
           image.ycen      = cal_ycen
           image.radius    = cal_radius
           image.vignette  = cal_vignette
           image.map_r     = cal_map_r
           image.map_theta = cal_map_theta
         ENDIF 

      ENDIF

      
;      shot = '/net/fuslsa/data/MAST_IMAGES/rda/rda013796.ipx'
;      frame = 2 ; 12
;      shot = '/net/fuslsa/data/MAST_IMAGES/rda/rda013807.ipx'
;      frame = 12
;      shot = '/net/fuslsa/data/MAST_IMAGES/rda/rda013816.ipx'
;      frame = 12
;      shot = '/net/fuslsa/data/MAST_IMAGES/rda/rda013948.ipx'
;      frame = 14
;      shot = '/net/fuslsa/data/MAST_IMAGES/rda/rda013803.ipx'
;      frame = 11
;      shot = '/net/fuslsa/data/MAST_IMAGES/rda/rda013802.ipx'
;      frame = 12 ; 13796 image processed so far 


      END
  ENDCASE

  IF (status EQ -1) THEN BEGIN
    PRINT,'Error LoadVignetteMap: map not found'
    PRINT,'  shot    : ',image.shot
    PRINT,'  channel : ',image.channel
    PRINT,'  line    : ',image.line
    STOP
  ENDIF

END
;
;----------------------------------------------------------------------
;
FUNCTION TraceRadii,image

  cal = GetImageCalibrationStructure()

  rfit = FINDGEN(image.radius)

  x = 1
  CASE x OF  
    1: BEGIN

      !P.MULTI = [0, 4, 2]
       window,2,retain=2,xsize=1200,ysize=800

       dtheta = 360.0 / FLOAT(cal.n)
       rtheta = MIN([2.0,0.5*dtheta])

       max_val = MEAN(image.raw(WHERE(image.map_r LT 5)))

       xzero = MAKE_ARRAY(100,/FLOAT,VALUE=0.0)
       yzero = MAKE_ARRAY(100,/FLOAT,VALUE=max_val)

       FOR i1 = 0, cal.n-1 DO BEGIN
         theta = FLOAT(i1) * dtheta

         theta1 = theta - rtheta
         theta2 = theta + rtheta

         map_theta = image.map_theta

         IF (theta1 LT 0.0 OR theta2 GT 360.0) THEN BEGIN
           IF (theta2 GT 360.0) THEN BEGIN
             theta1 = theta1 - 360.0
             theta2 = theta2 - 360.0
           ENDIF  
           i = WHERE(map_theta GT 180.0)
           map_theta[i] = map_theta[i] - 360.0
         ENDIF

         i = WHERE(map_theta GE theta1 AND map_theta LE theta2 AND  $
                   image.map_r LT image.radius)

         print, theta1, theta2, N_ELEMENTS(i)

         xdat = [xzero,image.map_r[i]]
         ydat = [yzero,image.raw  [i]] / max_val

         A = POLY_FIT(xdat,ydat,cal.ndeg)            
         vfit = POLY(rfit,A)

         plot,xdat,ydat,psym=6,xrange=[0.0,MAX(xdat)],xstyle=1
         oplot,rfit,vfit,color=2

         cal.angle[i1] = theta
         cal.max_r[i1] = MAX(xdat)
         cal.A    [i1,0:cal.ndeg] = A
         print,MAX(xdat)
       ENDFOR
 
       !P.MULTI = 0

       END

  ENDCASE  

  RETURN, cal

END
;
;----------------------------------------------------------------------
;
PRO ShiftImage,xshift,yshift,image

  xwidth = N_ELEMENTS(image[*,0])
  ywidth = N_ELEMENTS(image[0,*])

print,xwidth,ywidth
print,xshift,yshift

; x-axis shift:
  IF (xshift NE 0) THEN BEGIN
    FOR iy = 0, ywidth-1 DO BEGIN    
      IF ( xshift LT 0 ) THEN BEGIN
        FOR ix = 0, xwidth-1+xshift DO image[ix,iy] = image[ix-xshift,iy]
        image[xwidth+xshift:xwidth-1,iy] = -1  ; Blank the end of the image with artificial black
      ENDIF ELSE BEGIN
        FOR ix = xwidth-1, xshift, -1 DO image[ix,iy] = image[ix-xshift,iy]
        image[0:xshift-1,iy] = -1 
      ENDELSE
    ENDFOR
  ENDIF

; y-axis shift:
  IF (yshift NE 0) THEN BEGIN
    FOR ix = 0, xwidth-1 DO BEGIN
       IF ( yshift LT 0 ) THEN BEGIN
         FOR iy = 0, ywidth-1+yshift DO image[ix,iy] = image[ix,iy-yshift]
         image[ix,ywidth+yshift:ywidth-1] = -1 
       ENDIF ELSE BEGIN
         IF ( yshift NE 0 ) THEN BEGIN
           FOR iy = ywidth-1, yshift, -1 DO image[ix,iy] = image[ix,iy-yshift] 
           image[ix,0:yshift-1] = -1 
         ENDIF
       ENDELSE
    ENDFOR
  ENDIF

END
;
;----------------------------------------------------------------------
;
FUNCTION gfunct, x, a, pder

;  F = A[0] - EXP( A[1] * X )

;  IF N_PARAMS() GT 4 THEN $
;    pder = [[ replicate(1.0,N_ELEMENTS(X))  ], $
;            [ -1.0 * A[1] * EXP( A[1] * X ) ]]

  F = A[0] - A[1] * X ^ A[2] 
;  F = A[1] - A[2] * COS( X * !PI / 180.0 ) ^ 4

  IF N_PARAMS() GT 4 THEN $
    pder = [ [ replicate(1.0,N_ELEMENTS(X))        ], $
             [ -1.0 * X ^ A[2]                     ], $
             [ -1.0 * A[1] * A[2] * X ^ (A[2] - 1) ] ] 


  RETURN, f

END
;
;----------------------------------------------------------------------
;
PRO CURVEFIT_SL, x, y, a, fit

  
  fit=gfunct(x,a)

  nelements = N_ELEMENTS(x)

  b = a

  print,nelements

  minrms = 1.0E+10

  miny = MIN( y )

  FOR p1 = 20.0, 40.0, 1.0 DO BEGIN
    FOR p2 =  4.0,  6.0, 0.01 DO BEGIN
;  FOR p1 = 20.0, 40.0, 1.0 DO BEGIN
;    FOR p2 =  0.0,  45.0, 1.0 DO BEGIN


      a[1] = p1
      a[2] = p2

      fit=gfunct(x,a)

      IF ( MIN(fit) GT miny ) THEN BEGIN

        rms = 0.0
        FOR i1 = 0, nelements-1 DO BEGIN
          rms = rms + (y[i1] - fit[i1])^2
        ENDFOR
        rms = rms / FLOAT( nelements - 1)

;        print,'RMS:',a[1],a[2],rms

        IF ( rms LT minrms ) THEN BEGIN
          minrms = rms
          b = a
        ENDIF

      ENDIF

    ENDFOR
  ENDFOR 

  print,'MINRMS:',minrms,b

  a = b

  fit=gfunct(x,a)

;  RETURN,fit

END
