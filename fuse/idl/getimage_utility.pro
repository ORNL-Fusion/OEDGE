;
; ======================================================================
;
; subroutine: RemoveArtificialBlack
;
PRO RemoveArtificialBlack, image

  i = WHERE(image.raw LE -1)
  j = WHERE(image.raw EQ  0)

  IF (i[0] NE -1) THEN image.raw[i] = 0
  IF (j[0] NE -1) THEN image.raw[j] = 1

END
;
; ======================================================================
;
; subroutine: MaskImage
;
PRO MaskImage, image, maskradius, maskleft, maskright, masktop, maskbottom, maskbox, maskpoly

  image.mask = 1.0

  IF (KEYWORD_SET(maskradius)) THEN BEGIN
    PRINT, 'MASKRADIUS IS BROKEN FOR FFC...'
    STOP 
    i = WHERE(image.map_r GE maskradius)
    IF (i[0] NE -1) THEN image.mask[i] = 0.0
  ENDIF

  xshift = image.xshift
  yshift = image.yshift

  IF (KEYWORD_SET(maskleft  )) THEN image.mask[0                 :maskleft-1  +xshift,*] = 0.0
  IF (KEYWORD_SET(maskright )) THEN image.mask[maskright-1+xshift:image.xdim-1       ,*] = 0.0

  IF (KEYWORD_SET(masktop   )) THEN image.mask[*,0                  :masktop-1+yshift] = 0.0
  IF (KEYWORD_SET(maskbottom)) THEN image.mask[*,maskbottom-1+yshift:image.ydim-1    ] = 0.0

  IF (NOT KEYWORD_SET(maskbox)) THEN maskbox = 0
  IF (N_ELEMENTS(maskbox) GT 1) THEN image.mask[maskbox[0]-1+xshift:maskbox[2 ]-1+xshift,       $
                                                maskbox[1]-1+yshift:maskbox[3 ]-1+yshift] = 0.0
  IF (N_ELEMENTS(maskbox) GT 4) THEN image.mask[maskbox[4]-1+xshift:maskbox[6 ]-1+xshift,       $
                                                maskbox[5]-1+yshift:maskbox[7 ]-1+yshift] = 0.0
  IF (N_ELEMENTS(maskbox) GT 8) THEN image.mask[maskbox[8]-1+xshift:maskbox[10]-1+xshift,       $
                                                maskbox[9]-1+yshift:maskbox[11]-1+yshift] = 2.0

  IF (KEYWORD_SET(maskpoly)) THEN BEGIN
    image_pixels = MAKE_ARRAY(2,LONG(image.xdim)*LONG(image.ydim),/LONG)
    FOR ix = 1, image.xdim  DO BEGIN
      FOR iy = 1, image.ydim DO BEGIN
        image_pixels[0,ix-1+(iy-1)*LONG(image.xdim)] = ix
        image_pixels[1,ix-1+(iy-1)*LONG(image.xdim)] = iy
      ENDFOR
    ENDFOR
    pt1 = 0
    WHILE (pt1 LT N_ELEMENTS(maskpoly)) DO BEGIN
      npts = maskpoly[pt1]
      pts = MAKE_ARRAY(2,npts,/LONG)
      FOR i1 = 1, npts DO BEGIN
        pts[0,i1-1] = maskpoly[pt1+1+2*(i1-1)  ] + xshift
        pts[1,i1-1] = maskpoly[pt1+1+2*(i1-1)+1] + yshift
      ENDFOR
      o=obj_new('IDLanROI',pts)
      p=o->containspoints(image_pixels)
      i = WHERE(p GT 0)
      image.mask[image_pixels[0,i],image_pixels[1,i]] = 0.0
      obj_destroy,o
      pt1 = pt1 + 2*npts + 1
    ENDWHILE
  ENDIF
    
; Add 1.0 to avoid real black being mistaken for artifical black, i.e.
; make sure no no-masked pixels have a value 0.0;
  image.data = (1.0 + image.data) * image.mask[0:image.xdim-1,0:image.ydim-1] 

END
;
; ======================================================================
;
; subroutine: GetImageMaps

PRO GetImageMaps, image

  image.map_r     = -1.0
  image.map_theta = -1.0

  FOR ix = 0, image.xdim-1 DO BEGIN
    FOR iy = 0, image.ydim-1 DO BEGIN

      IF (image.raw[ix,iy] NE -1) THEN BEGIN

        x = FLOAT(ix+1-image.xcen)
        y = FLOAT(iy+1-image.ycen)

        radius = SQRT( x^2 + y^2 ) 

        IF (ABS(radius GT 0.0)) THEN BEGIN  ; floating point error if this check isn't made...?
          theta = ACOS( x / radius) * 180.0 / 3.1415
          IF (y GT 0.0) THEN theta = 360.0 - theta  
        ENDIF ELSE BEGIN
          theta = 0.0
        ENDELSE
        IF (theta EQ 360.0) THEN theta = 0.0

        image.map_r    [ix,iy] = radius
        image.map_theta[ix,iy] = theta

      ENDIF
    ENDFOR
  ENDFOR

END

;
; ======================================================================
;
; subroutine: AdjustImage
;
; It is possible for cameras to sub-window, i.e. use only a sub-region
; of the CCD when recording data.  If this is the case then the image
; data may need to be shifted within the image.raw space so that the data
; maps properly to the full view calibration data.
;
PRO AdjustImage, image, brighten, flip

  IF (KEYWORD_SET(brighten)) THEN BEGIN
    image.raw = UINT(FLOAT(image.raw) * brighten)
  ENDIF

  IF (KEYWORD_SET(flip)) THEN BEGIN
    image_flip = image.raw
    FOR ix = 0, image.xdim-1 DO BEGIN
      FOR iy = 0, image.ydim-1 DO BEGIN
        image_flip[ix,image.ydim-1-iy] = image.raw[ix,iy]
      ENDFOR
    ENDFOR
    image.raw = image_flip
  ENDIF

END
;
;
; ====================================================================================================
;
FUNCTION CleanImage, image


;  cleaned = MEDIAN(image,2)
;  RETURN, cleaned


  dim = SIZE(image,/DIMENSIONS)
  image_xwidth = dim[0]
  image_ywidth = dim[1]

  CASE 1 OF  
    1: BEGIN
      average = MAKE_ARRAY(image_xwidth,image_ywidth,/FLOAT,value=0.0)
      cleaned = MAKE_ARRAY(image_xwidth,image_ywidth,/UINT, value=0)

      FOR ix = 0, image_xwidth-1 DO BEGIN
        FOR iy = 0, image_ywidth-1 DO BEGIN
          IF ( image[ix,iy] GT 0 AND  $
               ix GT 0 AND ix LT image_xwidth-1 AND  $
               iy GT 0 AND iy LT image_ywidth-1) THEN BEGIN

            square = FLOAT([image[ix-1,iy-1],image[ix-1,iy  ],image[ix-1,iy+1],image[ix  ,iy-1],  $
                            image[ix  ,iy+1],image[ix+1,iy-1],image[ix+1,iy  ],image[ix+1,iy+1]])

            s1 = MEAN(square)
 
            average[ix,iy] = s1

            diff = ABS(FLOAT(image[ix,iy] - s1)) / s1  

;            IF (diff GT 2.0) THEN BEGIN
            IF (diff GT 1.05) THEN BEGIN
;              print,ix,iy
              cleaned[ix,iy] = s1
            ENDIF ELSE BEGIN
              cleaned[ix,iy] = image[ix,iy]
            ENDELSE

          ENDIF
        ENDFOR
      ENDFOR

;      window,4,xsize=image_xwidth,ysize=image_ywidth,retain=2
;      tvscl,average,/order

;      window,5,xsize=image_xwidth,ysize=image_ywidth,retain=2
;      tvscl,cleaned,/order

      END

  ENDCASE

  RETURN, cleaned
 
END

;
;----------------------------------------------------------------------
;
PRO MAXAVG, y, span, value

  nelements = N_ELEMENTS(y)

  sumavg = FINDGEN(nelements)
  sumavg[0:nelements-1] = 0.0

  FOR i1 = 0+span, nelements-span-1 DO BEGIN
     sumavg[i1] = TOTAL( y[i1-span:i1+span] ) 
  ENDFOR

  value = MAX(sumavg) / ( 2.0 * FLOAT(span) + 1 )

END
