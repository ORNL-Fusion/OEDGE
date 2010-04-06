;
; ======================================================================
;
; subroutine: PositionImageData
;
; It is possible for cameras to sub-window, i.e. use only a sub-region
; of the CCD when recording data.  If this is the case then the image
; data may need to be shifted within the image.raw space so that the data
; maps properly to the full view calibration data.
;
PRO PositionImageData, image

  IF (image.xwin[0] NE 1) THEN BEGIN
    PRINT,'ERROR: IMAGE X CORRECTION NOT READY'
    STOP
  ENDIF

; Check if image data needs to be shifted vertically:     
  IF (image.ywin[0] NE 1) THEN BEGIN
;   Setup:  
    height = image.ywin[1] - image.ywin[0] + 1 
    yshift = image.ywin[0] - 1
;   Move data:
    FOR i = height+yshift-1, yshift, -1 DO $
      image.raw[0:image.xdim-1,i] = image.raw[0:image.xdim-1,i-yshift]
;   Blank region above image data:
    image.raw[0:image.xdim-1,0:yshift-1] = -1  ; The -1 indicates artificial black
  ENDIF                                        ; (don't know how to assign global parameters yet...)

END
;
;
; ======================================================================
;
;
FUNCTION LoadImage,device,camera,file,date,line,shot,channel,frame,path,  $
                   calibrate,shift,radius,filter,window,clean,background, $
                   plots


  IF (KEYWORD_SET(file)) THEN BEGIN

    extension = STRMID(file,2,/REVERSE_OFFSET) 
  
    PRINT,'FILE EXTENSION = ',extension

    CASE extension OF
;     ------------------------------------------------------------------
      'bmp': image_raw = READ_BMP(file)
;     ------------------------------------------------------------------
      'ipx': BEGIN
        desc=ipx_open(file)
        IF (KEYWORD_SET(calibrate)) THEN BEGIN
;         Load all frames in file and average:       
          FOR frame = 0, N_ELEMENTS(desc.frametime) DO BEGIN
            image_raw = ipx_frame(desc,frame,TIME=time)
            IF (frame EQ 0) THEN BEGIN
              image_sum = image_raw
            ENDIF ELSE BEGIN
              image_sum = image_sum + image_raw
            ENDELSE        
          ENDFOR
          image_raw = image_sum / N_ELEMENTS(desc.frametime)
    
          print,'Calibration file:'
          print,file
    
          shot   = -1
          frame  = frame - 1
          time   = -1.0
        ENDIF ELSE BEGIN
          IF (NOT KEYWORD_SET(frame)) THEN frame = 1
          image_raw = ipx_frame(desc,frame,TIME=time)
        ENDELSE
        END
;     ------------------------------------------------------------------
      ELSE: BEGIN
        PRINT,'Error LoadImage: Unrecognized file extension'
        PRINT,'  EXTENSION = ',extension
        STOP
        END
    ENDCASE

    IF (NOT KEYWORD_SET(camera)) THEN camera = 'unknown'
    format = extension

  ENDIF ELSE BEGIN

    CASE device OF
;     ------------------------------------------------------------------
;     ------------------------------------------------------------------
      'CMOD': BEGIN
;       Assign dummy structure:
        file = '/net/fuslsa/data/MAST_IMAGES/rda/rda013950.ipx'
        frame = 1
        desc=ipx_open(file)
        CASE 2 OF
          1: BEGIN
            shot = './calib_dg_image_for_steve.gif'
            frame = 1
            READ_GIF,shot,image_byte
;           For the image header:
            image_xwidth  = 640
            image_ywidth  = 480
            image_xcentre = 320
            image_ycentre = 240   
            image = MAKE_ARRAY(image_xwidth,image_ywidth,/INTEGER,value=-1)
            FOR i = image_xwidth-1, 33, -1 DO BEGIN
              FOR j = 0, image_ywidth-1-10 DO BEGIN
                image[i,image_ywidth-1-j] = image_byte[i-33,j+10]
;                FOR j = image_ywidth-1, 10, -1 DO BEGIN
;                image[i,image_ywidth-1-j] = image_byte[i-33,j-10]
              ENDFOR
            ENDFOR
            pixel_depth = 8
            END
;         --------------------------------------------------------------
          2: BEGIN
;            file = path+'/images/1060725023_32_1.png'
            file = '/home/slisgo/divimp/shots/cmod/1060725023/1060725023_32_1.png'
            shot = 1060725023
            frame = 32
            time = 1.087
            ; for determining the centre of the image
            ; image=getimage(shot=1060725023,frame=32,channel=1,device='CMOD',/plots,/centre,radius=275,depth=11,shift=[65,25]) 
            ;file = '/home/slisgo/divimp/shots/cmod/1060725023/1060725017_lisgo_5.png'
            ;shot = 1060725017
            ;frame = 1
            ;time = 0.1
            image_raw = READ_PNG(file,/ORDER)
;           From the image header:
            image_xwidth  = 480
            image_ywidth  = 640
            image_xcentre = 240
            image_ycentre = 320   
            pixel_depth = 14
            format = 'png'
            IF (NOT KEYWORD_SET(window)) THEN window = 'UNKNOWN'
            END
        ENDCASE
        image_radius  = 200
        image_xshift =   0 
        image_yshift =   0
;        help,image,/STRUCT   
;        window,0,xsize=image_xwidth,ysize=image_ywidth,retain=2
;        tvscl,image
        desc.header.shot   = shot
        desc.header.filter = 'Dgamma 434 nm'
        desc.header.width  = image_xwidth
        desc.header.height = image_ywidth
        desc.header.right  = image_xwidth
        desc.header.bottom = image_ywidth
        END
;     ------------------------------------------------------------------
;     ------------------------------------------------------------------
      'MAST': BEGIN
        CASE camera OF 
;         --------------------------------------------------------------
          'MWIR': BEGIN
            format = 'ir'
            directory = STRING(shot/1000,FORMAT='(I3.3)')
            file = "$MAST_IMAGES/" + directory + "/" +  $
                   STRTRIM(STRING(shot),2) + "/rir0" +  $ 
                   STRTRIM(STRING(shot,FORMAT='(i5.5)'),2) + ".ipx"
            desc = ipx_open(file) 
            image_raw = ipx_frame(desc,frame,time=time)
            END
;         --------------------------------------------------------------
          'LWIR': BEGIN
            format = 'ir'
            directory = STRING(shot/1000,FORMAT='(I3.3)')
            file = "$MAST_IMAGES/" + directory + "/" +  $
                   STRTRIM(STRING(shot),2) + "/rit0" +  $ 
                   STRTRIM(STRING(shot,FORMAT='(i5.5)'),2) + ".ipx"
;            image_path = '/net/fuslsa/data/MAST_IMAGES/rir/rir0'  ; Poor hardcoding...
;            file = image_path+STRTRIM(STRING(shot),1)+'.ipx'
;...        Need to subtract the NUC manually:
            desc = ipx2open(file) 
            image_nuc = ipx2frame(desc,1,/ref)
            image_raw = ipx2frame(desc,frame,time=time) 
            image_raw = image_nuc - image_raw
            image_raw = 2^desc.fileinfo.depth - image_raw
            image_raw = image_raw < 2^desc.fileinfo.depth-1     
;...        Need to modify the DESC structure to conform with the old IPX convention:
            desc.fileinfo.left = desc.fileinfo.left + 1
            desc.fileinfo.top  = desc.fileinfo.top  + 1
            header = desc.fileinfo
            header = CREATE_STRUCT(header,'right' ,desc.fileinfo.left+desc.fileinfo.width -1)            
            header = CREATE_STRUCT(header,'bottom',desc.fileinfo.top +desc.fileinfo.height-1)            
            desc = CREATE_STRUCT(desc,'header',header)
            END
;         --------------------------------------------------------------
          'ZEBRA': BEGIN
            format = 'zeb'
            desc = read_rzz(0,0,shot)
            IF (NOT KEYWORD_SET(window)) THEN window = 'HM10'
            END
;         --------------------------------------------------------------
          'RGB': BEGIN
            PRINT,'RGB NOT SUPPORTED AT THE MOMENT'
            STOP
            format = 'rgb'
;            dum = getrgbdata(s=shot)  ; CAUSING A COMPILE ERROR AT THE MOMENT
            CASE channel OF
              1: desc = REFORM(dum.ls.r[*,*,frame])  ; Brem (TS=toroidally symmetric) + BES (not TS)
              2: desc = REFORM(dum.ls.g[*,*,frame])  ; Brem (TS), noisy
              3: desc = REFORM(dum.ls.b[*,*,frame])  ; Brem (TS), noisy
              4: desc = REFORM(dum.us.r[*,*,frame])  ; Dalpha (TS)
              5: desc = REFORM(dum.us.g[*,*,frame])  ; CVI (TS) + C-CX BES (not TS)
              6: desc = REFORM(dum.us.b[*,*,frame])  ; He+ (same as above?)
              ELSE: BEGIN
                PRINT,'ERROR LoadImage: Unrecognized RGB channel'
                PRINT,'  CHANNEL = ',channel
                END              
            ENDCASE
            IF (channel GE 1 and channel LE 3) THEN desc_t = dum.llt[frame]  ; Store the frame time
            IF (channel GE 4 and channel LE 6) THEN desc_t = dum.uut[frame]  
            desc_filter_id = 'NOT LISTED AT PRESENT'
            END
;         --------------------------------------------------------------
          'PHOTRON': BEGIN
            CASE channel OF
              1: image_path = '/net/fuslsa/data/MAST_IMAGES/rba/rba0'  ; Poor hardcoding...
              2: image_path = '/net/fuslsa/data/MAST_IMAGES/rbb/rbb0'
              3: image_path = '/net/fuslsa/data/MAST_IMAGES/rbc/rbc0'
;              3: image_path = '/net/fuslsa/data/MAST_NEW/rbc0'
            ENDCASE
            file = image_path+STRTRIM(STRING(shot),1)+'.ipx'
            print,file
            desc = ipx_open(file) 
            image_raw = ipx_frame(desc,frame,time=time)
            format = 'ipx_photron'
            END
;         --------------------------------------------------------------
          'LINCAM': BEGIN
            format = 'lnc'
;           IDL> restore, 'write_radius.fnc'
;           ms    - inverted data that is put into 'ADA_Dalpha inverted'
;           tcam  - time series
;           rscam - radial poistion
;           frcam - radial position for raw data (?)
;           fcam  - raw data that is put into 'ADA_Dalpha raw_full' (in Volts)
            show,shot,'m4',ms,tcam,rscam,cam,traw,acam,rcam,fcam,facam,frcam,tint
;            print,acam*180/3.14,N_ELEMENTS(acam)
;            stop
            help,fcam,/struct
            END
;         --------------------------------------------------------------
          'DIVCAM': BEGIN
            image_path = "$MAST_IMAGES/" + STRING(shot/1000,FORMAT='(I3.3)') + $ 
                         "/" +  STRTRIM(STRING(shot,FORMAT='(I5.5)'),2) 
            CASE channel OF
              1: image_path = image_path + "/rda0"  
              2: image_path = image_path + "/rdb0"  
              3: image_path = image_path + "/rdc0"  
              4: image_path = image_path + "/rdd0"
            ENDCASE

            file = image_path + STRTRIM(STRING(shot,FORMAT='(I5.5)'),2) + ".ipx"
; file = './rda000724.ipx'
; file = './rda015912.ipx'
; file = './rdd024808.ipx'
print,file
            desc = ipx_open(file) 
;help,desc.header,/struct
            image_raw = ipx_frame(desc,frame,time=time)
;
            IF (background) THEN BEGIN
              file = image_path+STRTRIM(STRING(background),1)+'.ipx'
              PRINT,file,' (background subtraction)'
              background_desc = ipx_open(file) 
              background_image_raw = ipx_frame(background_desc,frame,time=time)
              FOR i = 0, N_ELEMENTS(image_raw[*,0])-1 DO BEGIN
                FOR j = 0, N_ELEMENTS(image_raw[0,*])-1 DO BEGIN
                  IF (image_raw[i,j] GT background_image_raw[i,j]) THEN BEGIN
                    image_raw[i,j] = image_raw[i,j] - background_image_raw[i,j]
                  ENDIF ELSE BEGIN
                    image_raw[i,j] = 0
                  ENDELSE
                ENDFOR
              ENDFOR
            ENDIF
            format = 'ipx'
            END
;         --------------------------------------------------------------
          'FFC': BEGIN
            image_path = "$MAST_IMAGES/" + STRING(shot/1000,FORMAT='(I3.3)') + $ 
                         "/" +  STRTRIM(STRING(shot),2) 
            CASE channel OF
              1: image_path = image_path + "/rba0"  
              2: image_path = image_path + "/rbb0"  
              3: image_path = image_path + "/rbc0"  
            ENDCASE
            file = image_path + STRTRIM(STRING(shot,FORMAT='(i5.5)'),2) + ".ipx"

;image=getimage(shot=24815,frame=350,channel=3,camera='FFC',/plots)
; file = './rda000724.ipx'
; file = './rba023925.ipx'
            PRINT,file
            desc = ipx_open(file) 
            image_raw = ipx_frame(desc,frame,time=time)
;stop
            format = 'ipx_photron'
            END
        ENDCASE        
        END
;   --------------------------------------------------------------------
;   --------------------------------------------------------------------
    ENDCASE
  ENDELSE
;
;
;
  CASE format OF
;   ------------------------------------------------------------------
    'ipx_photron': BEGIN
      dim = SIZE(image_raw,/DIMENSIONS)
      image = GetImageDataStructure(dim[0],dim[1])                        ; Will need to convert to standard size for vignetting correction, as for DIVCAM
      image.raw[0:dim[0]-1,0:dim[1]-1] = image_raw[0:dim[0]-1,0:dim[1]-1] ; Need to remove exclusion call PositionImageData below
      image.depth   = desc.header.depth
      image.shutter = desc.header.exposure
      image.gain    = desc.header.gain[0]
      image.xwin = [desc.header.left,desc.header.right ]
      image.ywin = [desc.header.top ,desc.header.bottom]
      IF (NOT KEYWORD_SET(window)) THEN window = desc.header.view
      IF (NOT KEYWORD_SET(filter)) THEN filter = 'unknown'
      image.cal.window_tag = window
      END
;   ------------------------------------------------------------------
    'ipx': BEGIN
      image = GetImageDataStructure(1000,1000)
;     Need to make sure the data is mapped properly to the image space that has
;     been setup:
      dim = SIZE(image_raw,/DIMENSIONS)
;      print,'dim:',dim
      IF (dim[0] NE 1000 OR dim[1] NE 1000) THEN BEGIN
        IF (dim[0] EQ 500 AND dim[1] EQ 500) THEN BEGIN
          IF (clean) THEN BEGIN                ; Attempt to clean-up image here rather then after
            image_raw = CleanImage(image_raw)  ; calibration since the current 'cleaning' algorithm
            clean = 0                          ; fails...
          ENDIF
          FOR j = 0, dim[1]-1 DO BEGIN
            FOR i = 0, dim[0]-1 DO BEGIN
              image.raw[2*i  ,2*j  ] = image_raw[i,j]
              image.raw[2*i  ,2*j+1] = image_raw[i,j]
              image.raw[2*i+1,2*j  ] = image_raw[i,j]
              image.raw[2*i+1,2*j+1] = image_raw[i,j]
            ENDFOR
          ENDFOR
          image.raw = image.raw / 4.0  
        ENDIF ELSE BEGIN
          PRINT,'ERROR: UNRECOGNIZED DivCam IMAGE SIZE'
          PRINT,'  DIM:',dim[0],dim[1]
          STOP
        ENDELSE
;       For some reason desc.header.left,right,top and bottom are set for a 
;       1000x1000 image, so no need to reassign these at the moment...
      ENDIF ELSE  $
        image.raw[0:dim[0]-1,0:dim[1]-1] = image_raw[0:dim[0]-1,0:dim[1]-1]
      image.depth   = desc.header.depth
      image.shutter = desc.header.exposure
      image.gain    = desc.header.gain
      image.xwin = [desc.header.left,desc.header.right ]
      image.ywin = [desc.header.top ,desc.header.bottom]
      IF (NOT KEYWORD_SET(window)) THEN window = desc.header.view
      IF (NOT KEYWORD_SET(filter)) THEN filter = desc.header.filter
      END
;   ------------------------------------------------------------------
    'lnc': BEGIN
      image = GetImageDataStructure(1024,1)
      image.raw  = REVERSE(fcam[frame,*],2)
      image.data = REVERSE(fcam[frame,*],2) * 0.9 / 0.82  ; End of campaign window transmission correction from "MEB04 ViewportsList-JD.doc"
      image.depth   = -1                                ; calibration, see divimp/data/m-ohm-notes.txt
      image.shutter = -1.0
      image.gain    = -1.0
      image.xwin = 1
      image.ywin = 1
      IF (NOT KEYWORD_SET(window)) THEN window = 'HM10'
      IF (NOT KEYWORD_SET(filter)) THEN filter = 'lnc_Dalpha'
      time = tcam[frame]
      image.cal.window_tag = window
      END
;   ------------------------------------------------------------------
    'zeb': BEGIN
      image = GetImageDataStructure(128,128)
      dim = SIZE(desc.image_raw,/DIMENSIONS)
      FOR j = 0, dim[1]-1 DO image.raw [j,*] = desc.image_raw [dim[1]-1-j,*,frame]  ; Which way is up?
      FOR j = 0, dim[1]-1 DO image.data[j,*] = desc.image_cal2[dim[1]-1-j,*,frame]
;      image.raw  = desc.image_raw [*,*,frame]
;      image.data = desc.image_cal2[*,*,frame]
      image.depth   = -1
      image.shutter = -1.0
      image.gain    = -1.0
      image.xwin = 1
      image.ywin = 1
      IF (NOT KEYWORD_SET(window)) THEN window = 'HM10'
      IF (NOT KEYWORD_SET(filter)) THEN filter = STRING(desc.filter_id)
      time = desc.t[frame]
      image.cal.window_tag = window
      END
;   ------------------------------------------------------------------
    'rgb': BEGIN
      image = GetImageDataStructure(640,480)
      FOR j = 0, 479 DO image.data[*,j] = desc[*,479-j]
      image.raw    = image.data
      image.depth   = -1
      image.shutter = -1.0
      image.gain    = -1.0
      image.xwin    = 1
      image.ywin    = 1
      IF (NOT KEYWORD_SET(window)) THEN window = 'HM??'
      IF (NOT KEYWORD_SET(filter)) THEN filter = STRING(desc_filter_id)
      time = desc_t
      image.cal.window_tag = window
      END
;   ------------------------------------------------------------------
    'ir': BEGIN
      image = GetImageDataStructure(desc.header.width,desc.header.height)
      image = CREATE_STRUCT(image,'data',image_raw)
;      image.data    = image_raw
      image.raw     = image.data
      image.depth   = desc.header.depth
      image.shutter = desc.header.exposure
      image.gain    = desc.header.gain
      image.xwin    = [desc.header.left,desc.header.right ]
      image.ywin    = [desc.header.top ,desc.header.bottom]
      IF (NOT KEYWORD_SET(window)) THEN window = desc.header.view
      IF (NOT KEYWORD_SET(filter)) THEN filter = desc.header.filter
      time = desc.frametime[frame]
      END
;   ------------------------------------------------------------------
    ELSE: BEGIN
      time = -1.0
      dim = SIZE(image_raw,/DIMENSIONS)
      print,dim
      image = GetImageDataStructure(dim[0],dim[1])
      IF (N_ELEMENTS(dim) EQ 3) THEN BEGIN
        image.raw[*,*] = image_raw[0,*,*] + image_raw[1,*,*] +  $
                         image_raw[2,*,*]
      ENDIF ELSE BEGIN
        image.raw = image_raw
      ENDELSE
;      print,'width,height:',dim
;      print,'            :',image.xdim,image.ydim
      image.depth  = 8              ; *** BAD ***
;      window,0,xsize=image_xwidth,ysize=image_ywidth,retain=2
;      tvscl,image
      END
  ENDCASE

  image.device  = device
  image.camera  = camera
  image.window  = window
  image.format  = format
  image.shot    = shot
  image.date    = date
  image.channel = channel
  image.frame   = frame
  image.time    = time
;
;
;
  IF (KEYWORD_SET(calibrate)) THEN BEGIN
    IF (NOT KEYWORD_SET(shift)) THEN shift = [0,0]
    image.cal.file = file
    image.xcen     = image.xcen + shift[0]
    image.ycen     = image.ycen + shift[1]
    IF (KEYWORD_SET(radius)) THEN image.radius = radius
  ENDIF ELSE BEGIN
    IF (KEYWORD_SET(filter)) THEN AssignFilterData, image, filter, plots
  ENDELSE
;
;
;
  IF (camera NE 'PHOTRON' AND camera NE 'FFC') THEN PositionImageData, image

  RETURN,image
END
;
; ======================================================================
;
