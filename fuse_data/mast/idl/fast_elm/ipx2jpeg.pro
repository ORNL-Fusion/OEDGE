; ELM movie
;
; ipx2jpeg,15622,startframe=382,lastframe=494
; 

PRO ipx2jpeg, shotnr, startframe=startframe, lastframe=lastframe, rotate=rotate

	; SPecify filename
; slmod begin
;   filename = "$MAST_IMAGES/rba/rba0"+$ 
;                   strtrim(STRING(shotnr,FORMAT='(i5.5)'),2)+".ipx"
;
   directory = STRING(shotnr/1000,FORMAT='(I3.3)')

   filename = "$MAST_IMAGES/" + directory + "/" +  $
              STRTRIM(STRING(shotnr),2) + "/rbb0" +  $ 
              STRTRIM(STRING(shotnr,FORMAT='(i5.5)'),2) + ".ipx"
   print,filename

   IF (KEYWORD_SET(rotate)) THEN BEGIN
     rotate_code = 0
     IF (rotate EQ 90 ) THEN rotate_code = 1
     IF (rotate EQ 180) THEN rotate_code = 2
     IF (rotate EQ 270) THEN rotate_code = 3
   ENDIF
; slmod end


;   shotnr  = strtrim(STRING(i,FORMAT='(i5.5)'),2)
; slmod begin
;   cmd = 'ls $MAST_IMAGES/rba/ | grep -ni '+$
;           'rba0'+strtrim(shotnr, 2)+'.ipx'
;
;   cmd = 'ls $MAST_IMAGES/rbb/ | grep -ni '+$
;           'rbb0'+strtrim(shotnr, 2)+'.ipx'
;
   cmd = 'ls $MAST_IMAGES/' + directory + '/' +  $
         STRTRIM(STRING(shotnr),2) +  $ 
         ' | grep -ni '+$
         'rbb0' + STRTRIM(STRING(shotnr,FORMAT='(i5.5)'),2) + ".ipx"

   print,cmd
; slmod end
   spawn, cmd , monster, count
        
  ;check length of file by using STRLEN(file)
                
   IF  (strlen(monster)) EQ  0 THEN BEGIN
          PRINT, 'Error, no camera images for shot'+'   '+strtrim(shotnr,2)

   ENDIF   ELSE     BEGIN

 ; Pick a frame - this is to be the middle frame
       fd = IPX_OPEN(filename, Rate=rate)
       shotnr  = fd.header[0].shot
       nframes = fd.header[0].numframes
       width   = fd.header.width
       height  = fd.header.height

IF NOT keyword_set(startframe) THEN startframe = 0
IF NOT keyword_set(lastframe) THEN lastframe = nframes

 
;       startframe= 0
;       startframe=round(startframe)
;       lastframe=10

       readnumber = lastframe - startframe + 1
       times = FLTARR(readnumber)
;slmod begin
       count = 0
       loadct,3  ; 3
; slmod end
       FOR n = startframe, lastframe-1 DO BEGIN
            ; Catch frame by specifying fd and n
            frame = ipx_frame(fd, n, time=ftime)
            times[n-startframe] = ftime

         
; slmod begin
            IF (KEYWORD_SET(rotate_code)) THEN frame = ROTATE(frame,rotate_code)


            dim = SIZE(frame,/DIMENSIONS)

;            frame = REBIN(frame,dim[0]*2,dim[1]*2)

            dim = SIZE(frame,/DIMENSIONS)

            print,'MAX 1',MAX(frame)

            frame[0,0] = 1023
            print,'MAX 2',MAX(frame)
            image_byte = BYTSCL(frame)
            window,6,xsize=dim[0],ysize=dim[1],retain=2
            image_byte = image_byte*1.5
            tvscl,image_byte,/order

;stop
;            tv,image_byte
            TVLCT, red, green, blue, /GET
            imageRGB = BYTARR(3, dim[0], dim[1])
            imageRGB[0, *, *] = red  [image_byte]
            imageRGB[1, *, *] = green[image_byte]
            imageRGB[2, *, *] = blue [image_byte]
            fname = './'+strtrim(shotnr,2)$
                        +'_frame_'+strtrim(n,2)+'.jpg'
            WRITE_JPEG, fname, imageRGB, TRUE = 1, QUALITY = 100, /ORDER
            wdelete,6
            PRINT,  shotnr, n
;          ENDELSE
;
;       convert data from 10 to 8 bit as write_jpeg only uses 8 bit numbers
;
;            frame= REVERSE(frame, 2)
;            frame = frame*255./1024.
;           WRITE_JPEG,'./'+strtrim(shotnr,2)$
;                      +'_frame_'+strtrim(n,2)+'.jpg', frame, QUALITY=100
; slend
        ENDFOR
    ENDELSE
END
