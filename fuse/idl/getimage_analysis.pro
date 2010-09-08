;
; ======================================================================
;
PRO AnalyseImage, image

  xdim = image.xdim
  ydim = image.ydim

  val_y = MAKE_ARRAY(xdim,/FLOAT,VALUE=0)


  iy1 = 445  ; 895 ; 465
  iy2 = 535  ; 895 ; 535

  FOR ix = 0, xdim-1 DO BEGIN
    val_y[ix] = FLOAT(TOTAL(image.data[ix,iy1:iy2])) / FLOAT(iy2 - iy1 + 1)
  ENDFOR

  max_val = MAX(val_y,imax)

;  print,val_y,max_val

;  print,"MAXIMUM:",FLOAT(ABS(imax-(xdim/2-image_xshift)))/FLOAT(xdim/2),  $
;                   FLOAT(max_val)/5687.0,max_val

  image_data=image.data
  image_data[*,iy1] = max_val
  image_data[*,iy2] = max_val

  val_y = val_y / max_val


  loadct,4
  window,3,xsize=xdim,ysize=ydim,retain=2
  tvscl,image_data,/order

  xcen = 480

  safe_colors,/first
  window,0,retain=2
  plot,REVERSE(val_y[0:xcen])
  oplot,val_y[xcen:xdim-1],color=2

  filter = GetFilterTransmissionProfile(image.filter)

  oplot,filter.nr*500,filter.tr,psym=6,color=3

  STOP

END
;
; ======================================================================
;
; subroutine: GetImageStatistics

PRO GetImageStatistics

; Initialization:
  saturation_count = 0.0
  saturation_total = 0.0
  saturation_value = 2 ^ image.depth - 2
  saturation_fraction = 0.0

  FOR ix = 0, image.xdim-1 DO BEGIN
    FOR iy = 0, image.ydim-1 DO BEGIN

      IF (1) THEN BEGIN

        IF ( pixel_radius LE FLOAT(image.radius) ) THEN BEGIN

;         Saturation fraction:
          saturation_total = saturation_total + 1.0
          IF ( image[ix,iy] EQ saturation_value ) THEN saturation_count = saturation_count + 1.0

        ENDIF ELSE BEGIN

;         Blank pixel:
          IF ( maskradius NE 0) THEN image[ix,iy] = -1   ; need a flag to turn this off when centring the image on the lens

        ENDELSE

;...    ...
        IF ( centre ) THEN BEGIN
;          image[ix,iy] = image[ix,iy] ^ 0.3

          IF ( image[ix,iy] GT 0.05 * saturation_value ) THEN image[ix,iy] = 0.05 * saturation_value
        ENDIF


      ENDIF  ; End of check for artifical black

    ENDFOR
  ENDFOR

END
;
; ======================================================================
;
