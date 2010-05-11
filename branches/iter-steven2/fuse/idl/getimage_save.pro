;
; ======================================================================
;
PRO SaveImageData, image, path, scale, sname, save_png, calibrate=calibrate, order=order, colour=colour

  binary = 0

  IF (NOT KEYWORD_SET(colour)) THEN colour = 3

  IF (KEYWORD_SET(calibrate)) THEN BEGIN

    IF (NOT KEYWORD_SET(path)) THEN BEGIN
      PRINT,'ERROR getimage_SaveImageData: Calibration save path not set'
      RETURN
    ENDIF

    str = STRSPLIT(image.cal.file,'/',/extract)

    i = N_ELEMENTS(str)

    fname = path + 'cal_' + str[i-1] + '.sav'  
;    fname = path + 'cal_' + str[i-2] + '_' + str[i-1] + '.sav'  

    cal_xcen	  = image.xcen     
    cal_ycen	  = image.ycen    
    cal_radius	  = image.radius  
    cal_vignette  = image.vignette 
    cal_map_r	  = image.map_r  
    cal_map_theta = image.map_theta 

    SAVE, cal_xcen,cal_ycen,cal_radius,cal_vignette,cal_map_r,cal_map_theta, $
          filename=fname,/verbose,/compress

  ENDIF ELSE BEGIN

    IF ( binary EQ 0 ) THEN BEGIN

      output_path = path + '/'

      IF (NOT KEYWORD_SET(sname)) THEN BEGIN
        fname = output_path+                           $
                STRTRIM(       image.camera  ,1)+'_'+  $
                STRTRIM(STRING(image.shot   ),1)+'_'+  $
                STRTRIM(STRING(image.frame  ),1)+'_'+  $
                STRTRIM(STRING(image.channel),1)
      ENDIF ELSE BEGIN
        fname = sname
      ENDELSE
 
      ; IDL save file
      image_shot    = image.shot
      image_frame   = image.frame
      image_time    = image.time
      image_channel = image.channel
      image_data    = image.data
      IF (NOT KEYWORD_SET(save_png)) THEN  $
        SAVE, image_shot, image_frame, image_time, image_channel, image_data,  $
              filename=fname+'.sav',/verbose,/compress
            
      ; JPEG
      IF ( image.ydim GT 1 ) THEN BEGIN     

;        max_val = LONG(MAX(image_data))
;        bottom = max_val * 20 / 100
;        print,'BOTTOM:',bottom,max_val,scale
;        image_data = image.data * scale
;        FOR ix = 0, image.xdim-1 DO BEGIN
;          FOR iy = 0, image.ydim-1 DO BEGIN
;            IF (image_data[ix,iy] GT max_val) THEN image_data[ix,iy] = max_val
;            IF (image_data[ix,iy] LT bottom ) THEN image_data[ix,iy] = bottom
;          ENDFOR
;        ENDFOR

;        WINDOW,6,xsize=image.xdim,ysize=image.ydim,retain=2
;        DEVICE, DECOMPOSED=0
;        LOADCT, 5
;        TVSCL,image_data,ORDER=order
;        WSET,6
;return
;        image_capture = TVRD(True=1)


; For IR images I think...
;        max_val = MAX(image_data)
;        bottom = 0.9* max_val
;        print,'BOTTOM:',bottom,scale
;        image_data = image.data * scale
;        FOR ix = 0, image.xdim-1 DO BEGIN
;          FOR iy = 0, image.ydim-1 DO BEGIN
;            IF (image_data[ix,iy] GT max_val) THEN image_data[ix,iy] = max_val
;            IF (image_data[ix,iy] LT bottom ) THEN image_data[ix,iy] = bottom
;          ENDFOR
;        ENDFOR

        max_val = MAX(image_data)
        image_data = image_data * scale
        i = WHERE(image_data GT max_val,count)
        IF (count GT 0) THEN image_data[i] = max_val

        image_byte = BYTSCL(image_data)
        WINDOW,6,xsize=image.xdim,ysize=image.ydim,retain=2,/PIXMAP
        DEVICE, DECOMPOSED=0
        WSET,6
        LOADCT, colour
        TV,image_byte,/order
        TVLCT, red, green, blue, /GET
        imageRGB = BYTARR(3, image.xdim, image.ydim)
        imageRGB[0, *, *] = red  [image_byte]  
        imageRGB[1, *, *] = green[image_byte]  
        imageRGB[2, *, *] = blue [image_byte] 
        print,fname+'.jpg'

        IF (KEYWORD_SET(save_png)) THEN BEGIN
          WRITE_PNG, fname+'.png', imageRGB, /ORDER
          WDELETE,6
          RETURN
        ENDIF ELSE BEGIN
          help,red
          WRITE_JPEG, fname+'.jpg', imageRGB, TRUE = 1, QUALITY = 100, /ORDER

          WDELETE,6
        ENDELSE
      ENDIF 


      image_data = image.store2
      max_val = MAX(image_data)
      image_data = image_data * scale
      i = WHERE(image_data GT max_val,count)
      IF (count GT 0) THEN image_data[i] = max_val

      image_byte = BYTSCL(image_data)
      dim = SIZE(image_data,/DIMENSIONS)
      WINDOW,6,xsize=dim[0],ysize=dim[1],retain=2,/PIXMAP
      DEVICE, DECOMPOSED=0
      WSET,6
      LOADCT, colour
      TV,image_byte,/order
      TVLCT, red, green, blue, /GET
      imageRGB = BYTARR(3, dim[0], dim[1])

print,image.xdim,image.ydim
print,dim[0],dim[1]

      imageRGB[0, *, *] = red  [image_byte]  
      imageRGB[1, *, *] = green[image_byte]  
      imageRGB[2, *, *] = blue [image_byte] 
      print,fname+'_raw.jpg'
      WRITE_JPEG, fname+'_raw.jpg', imageRGB, TRUE = 1, QUALITY = 100, /ORDER



      ; RAY file:
      OPENW,unit,fname+'.idl',/GET_LUN,ERROR=err    

      IF (err NE 0) THEN BEGIN
        PRINTF,-2,!ERR_STRING
      ENDIF ELSE BEGIN

        PRINTF,unit,'{version}   : ',1.1                  ,FORMAT='(A,F10.2)'

        PRINTF,unit,'ID          : ','blank'              ,FORMAT='(2A)'
        PRINTF,unit,'SIZE        : ',-1                   ,FORMAT='(A,I12)'
        PRINTF,unit,'CODEC       : ','blank'              ,FORMAT='(2A)'
        PRINTF,unit,'DATA/TIME   : ','blank'              ,FORMAT='(2A)'
        PRINTF,unit,'SHOT        : ',image.shot           ,FORMAT='(A,I12)'
        PRINTF,unit,'FRAME       : ',image.frame          ,FORMAT='(A,I12)'
        PRINTF,unit,'TIME        : ',image.time           ,FORMAT='(A,F12.6)'
        PRINTF,unit,'CHANNEL     : ',image.channel        ,FORMAT='(A,I12)'
        PRINTF,unit,'TRIGGER     : ',-1.0                 ,FORMAT='(A,F12.8)'
        PRINTF,unit,'LENS        : ','blank'              ,FORMAT='(2A)'
        PRINTF,unit,'FILTER      : ',image.filter.tag     ,FORMAT='(2A)'
        PRINTF,unit,'VIEW        : ',image.cal.window_tag ,FORMAT='(2A)'
        PRINTF,unit,'NUMFRAMES   : ',-1                   ,FORMAT='(A,I12)'
        PRINTF,unit,'CAMERA      : ','blank'              ,FORMAT='(2A)'
        PRINTF,unit,'WIDTH       : ',image.xdim           ,FORMAT='(A,I12)'
        PRINTF,unit,'HEIGHT      : ',image.ydim           ,FORMAT='(A,I12)'
        PRINTF,unit,'DEPTH       : ',image.depth          ,FORMAT='(A,I12)'
        PRINTF,unit,'ORIENTATION : ',-1                   ,FORMAT='(A,I12)'
        PRINTF,unit,'TAPS        : ',-1                   ,FORMAT='(A,I12)'
        PRINTF,unit,'COLOUR      : ',-1                   ,FORMAT='(A,I12)'
        PRINTF,unit,'HBIN        : ',-1                   ,FORMAT='(A,I12)'
        PRINTF,unit,'LEFT        : ',-1                   ,FORMAT='(A,I12)'
        PRINTF,unit,'RIGHT       : ',-1                   ,FORMAT='(A,I12)'
        PRINTF,unit,'VBIN        : ',-1                   ,FORMAT='(A,I12)'
        PRINTF,unit,'TOP         : ',-1                   ,FORMAT='(A,I12)'
        PRINTF,unit,'BOTTOM      : ',-1                   ,FORMAT='(A,I12)'
        PRINTF,unit,'OFFSET      : ',-1,-1                ,FORMAT='(A,2I12)'
        PRINTF,unit,'GAIN        : ',-1.0,-1.0            ,FORMAT='(A,2F12.8)'
        PRINTF,unit,'PREEXP      : ',-1                   ,FORMAT='(A,I12)'
        PRINTF,unit,'EXPOSURE    : ',-1                   ,FORMAT='(A,I12)'
        PRINTF,unit,'STROBE      : ',-1                   ,FORMAT='(A,I12)'
        PRINTF,unit,'CCD_TEMP    : ',-1.0                 ,FORMAT='(A,F12.8)'

;        PRINTF,unit,'ID          : ',desc.header.id       ,FORMAT='(2A)'
;        PRINTF,unit,'SIZE        : ',desc.header.size     ,FORMAT='(A,I12)'
;        PRINTF,unit,'CODEC       : ',desc.header.codec    ,FORMAT='(2A)'
;        PRINTF,unit,'DATA/TIME   : ',desc.header.date_time,FORMAT='(2A)'
;        PRINTF,unit,'SHOT        : ',desc.header.shot     ,FORMAT='(A,I12)'
;        PRINTF,unit,'FRAME       : ',frame                ,FORMAT='(A,I12)'
;        PRINTF,unit,'TIME        : ',time                 ,FORMAT='(A,F12.6)'
;        PRINTF,unit,'CHANNEL     : ',channel              ,FORMAT='(A,I12)'
;        PRINTF,unit,'TRIGGER     : ',desc.header.trigger  ,FORMAT='(A,F12.8)'
;        PRINTF,unit,'LENS        : ',desc.header.lens     ,FORMAT='(2A)'
;        PRINTF,unit,'FILTER      : ',desc.header.filter   ,FORMAT='(2A)'
;        PRINTF,unit,'VIEW        : ',desc.header.view     ,FORMAT='(2A)'
;        PRINTF,unit,'NUMFRAMES   : ',desc.header.numframes,FORMAT='(A,I12)'
;        PRINTF,unit,'CAMERA      : ',desc.header.camera   ,FORMAT='(2A)'
;        PRINTF,unit,'WIDTH       : ',desc.header.width    ,FORMAT='(A,I12)'
;        PRINTF,unit,'HEIGHT      : ',desc.header.height   ,FORMAT='(A,I12)'
;        PRINTF,unit,'DEPTH       : ',desc.header.depth    ,FORMAT='(A,I12)'
;        PRINTF,unit,'ORIENTATION : ',desc.header.orient   ,FORMAT='(A,I12)'
;        PRINTF,unit,'TAPS        : ',desc.header.taps     ,FORMAT='(A,I12)'
;        PRINTF,unit,'COLOUR      : ',desc.header.color    ,FORMAT='(A,I12)'
;        PRINTF,unit,'HBIN        : ',desc.header.hBin     ,FORMAT='(A,I12)'
;        PRINTF,unit,'LEFT        : ',desc.header.left     ,FORMAT='(A,I12)'
;        PRINTF,unit,'RIGHT       : ',desc.header.right    ,FORMAT='(A,I12)'
;        PRINTF,unit,'VBIN        : ',desc.header.vBin     ,FORMAT='(A,I12)'
;        PRINTF,unit,'TOP         : ',desc.header.top      ,FORMAT='(A,I12)'
;        PRINTF,unit,'BOTTOM      : ',desc.header.bottom   ,FORMAT='(A,I12)'
;        PRINTF,unit,'OFFSET      : ',desc.header.offset   ,FORMAT='(A,2I12)'
;        PRINTF,unit,'GAIN        : ',desc.header.gain     ,FORMAT='(A,2F12.8)'
;        PRINTF,unit,'PREEXP      : ',desc.header.preExp   ,FORMAT='(A,I12)'
;        PRINTF,unit,'EXPOSURE    : ',desc.header.exposure ,FORMAT='(A,I12)'
;        PRINTF,unit,'STROBE      : ',desc.header.strobe   ,FORMAT='(A,I12)'
;        PRINTF,unit,'CCD_TEMP    : ',desc.header.ccd_temp ,FORMAT='(A,F12.8)'

        FOR iy = 0, image.ydim-1 DO BEGIN
           FOR ix = 0, image.xdim-1, 10 DO BEGIN
              ix1 = ix
              ix2 = ix1 + 9
              IF ( ix2 GT image.xdim-1 ) THEN ix2 = image.xdim-1
;              IF ( ix2 GT image.ydim-1 ) THEN ix2 = image.xdim-1  ...BUG!
              PRINTF,unit,image.data[ix1:ix2,iy],FORMAT='(10E12.4)'        
           ENDFOR
        ENDFOR

        CLOSE,unit

      ENDELSE

    ENDIF ELSE BEGIN
;     Attempting to write the image file as a binary, but this just seems to
;     create loads of problems when reading the file into OUT (FORTRAN):

      fname = path+'/test-data'
      OPENW,unit,fname,/GET_LUN,ERROR=err
      FOR ix = 0, 999 DO BEGIN
        FOR iy = 0, 999 DO BEGIN
          WRITEU,unit,image2(ix,iy)           
        ENDFOR
      ENDFOR
;      WRITEU,unit,image2
      CLOSE,unit      


      fname = path+'/test-header'
      OPENW,unit,fname,/GET_LUN,ERROR=err
      IF (err NE 0) THEN BEGIN
        PRINTF,-2,!ERR_STRING
      ENDIF ELSE BEGIN

        desc.header.id = STRING(desc.header.id,FORMAT='(A8)')
        desc.header.codec = STRING(desc.header.codec,FORMAT='(A8)')
        desc.header.date_time = STRING(desc.header.date_time,FORMAT='(A20)')
        desc.header.lens = STRING(desc.header.lens,FORMAT='(A24)')
        desc.header.filter = STRING(desc.header.filter,FORMAT='(A24)')
        desc.header.view = STRING(desc.header.view,FORMAT='(A64)')
        desc.header.camera = STRING(desc.header.camera,FORMAT='(A64)')

        desc.header.numFrames = 1

        PRINT,'Writing file test.ipx'

        WRITEU,unit,desc.header

        CLOSE,unit

      ENDELSE

    ENDELSE

  ENDELSE

  print,'IMAGE DATA SAVED'


END

