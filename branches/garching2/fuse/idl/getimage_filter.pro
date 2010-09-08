;
; ======================================================================
;
FUNCTION FilterBlueShift,wlngth,index,angle

  shift = wlngth * SQRT(1.0 - (1.0 / index^2) * (SIN(angle*3.141592/180.0))^2)  ; Replace with PI global

  RETURN, (shift - wlngth)
END
;
; ======================================================================
;
FUNCTION FilterTransmission,filter,shift

  delta = ABS(filter.cwl + shift - filter.line) 

  delta = delta / (filter.fwhm / 2.0)  ;   ??? THIS FACTOR 2 IS CORRECT ???

   IF ( 0 ) THEN BEGIN
    print,filter.line
    print,filter.cwl
    print,'2:',shift
    print,'3:',ABS(filter.cwl+shift-filter.line)
    print,'d:',delta
  ENDIF

  n = N_ELEMENTS(WHERE(filter.yshape NE 0.0))

  tr = INTERPOL(ALOG10(filter.yshape[0:n-1]),filter.xshape[0:n-1],delta)

  RETURN, (10.0^tr)

END
;
; ======================================================================
;
FUNCTION PixelTransmission,distance,radius,height,filter

; Parameters:
  d = DOUBLE(distance)  ; 310.0       ; distance from field lens to filter (different for each camera?)
  r = DOUBLE(radius)    ; 10.0 ; 16.7 / 2.0  ; radius of image at filter 
  h = DOUBLE(height)    ; 0.0; 30.0        ; height of image at field lens


  PI = 3.1415926536D

; Distribution of points on filter:
  nrad = 20.0D
  nang = 60.0D
  drad = r / nrad
  dang = 360.0D / nang
  num = UINT(nrad*nang)
;  print,num
  x = MAKE_ARRAY(num,/DOUBLE,value=0.0D)
  y = MAKE_ARRAY(num,/DOUBLE,value=0.0D)
  w = MAKE_ARRAY(num,/DOUBLE,value=0.0D)

; Define ray intersections:
  i = -1
  FOR rad = 0.5D*drad, r, drad DO BEGIN
    FOR ang = 0.0D, 359.0D, dang DO BEGIN
      i = i + 1
      x[i] = rad * COS(ang*PI/180.0D)
      y[i] = rad * SIN(ang*PI/180.0D)
      w[i] = drad * (2.0D * PI * rad) / nang  ; not quite right, especially near the centre
    ENDFOR
;    print,rad,r,w[i]
  ENDFOR
  w = w / (PI * r^2)  

;  print,'x:',x
;  print,'y:',y

; Distance from image to lens:
  a = SQRT( x^2 + (y-h)^2 + d^2 )

; Angle of ray incidence (relative to normal):
  angle = ACOS( d / a ) * 180.0D / PI

  shift = FilterBlueShift(filter.cwl,filter.index,angle)

  tr = FilterTransmission(filter,shift)

  average_tr = TOTAL( tr * w ) 



  IF ( 0 ) THEN BEGIN
    print,'a:',a
    print,'d/a:',d/a
    print,'angle:',angle
    print,'shift:',shift
    print,'tr:',tr
    print,'w:',w
    print,'i/num:',i,num
    print,TOTAL(w)
    print,PI*r^2
    print,'TRANS:',d,r,h
    print,'LINE/CWL:',filter.line,filter.cwl
    print,average_tr
  ENDIF
  
  RETURN,average_tr

END
;
; ======================================================================
;
FUNCTION GetFilterTransmissionProfile, filter, camera, plots


  optics ={id           : 0  ,  $
           distance     : 0.0,  $ 
           radius       : 0.0,  $
           height       : 0.0 }


  CASE camera OF
    'DIVCAM': BEGIN
      optics.id       = 0
      optics.distance = 310.0 - 22.5 ; 290.0 - 22.5      
      optics.radius   = 10.0 ; 10.0    
      optics.height   = 30.0 ; 33.0
      END
    'FFC': BEGIN
      optics.id       = 0
      optics.distance = 310.0 - 22.5 + 51.0  ; Extra tube length, but not quite enough... measured by James on 08/04/2010
      optics.radius   = 10.0 ; 10.0          ; am I making an error here, in the geometry assumptions/setup, and 
      optics.height   = 30.0 ; 33.0          ; the fact that I'm not filling the CCD...
      END
    ELSE: BEGIN
      PRINT,'ERROR getimage_GetFilterTransmissionProfile: Camera not recognised'
      PRINT,'  CAMERA= ',camera
      END
  ENDCASE

  print,'LINE:',filter.line

  n = filter.n
  
  nr = MAKE_ARRAY(n,/FLOAT,VALUE=0.0)
  tr = MAKE_ARRAY(n,/FLOAT,VALUE=0.0)

  FOR i = 0, n-1 DO BEGIN
    nr[i] = FLOAT(i) / FLOAT(n-1)
    height = nr[i] * optics.height
    tr[i] = PixelTransmission(optics.distance,optics.radius,height,filter)
  ENDFOR

  filter.nr = nr
  filter.tr = tr

;  filter.tr = 1.0

  IF (KEYWORD_SET(plots)) THEN BEGIN
    window,3,retain=2
    plot,filter.nr*500,filter.tr,psym=6
  ENDIF

  RETURN, filter

END
;
; ======================================================================
; 
PRO GetDivCamFilterData, image, wavelength, plots

  CASE wavelength OF
    411.2: BEGIN
      image.line                = 411.2
      image.filter.line         = 410.174 ; air?
      image.filter.tag          = 'Ddelta'
      image.filter.cwl          = 411.18
      image.filter.fwhm         = 5.35
      image.filter.transmission = 0.4502
      image.filter.cavities     = 2
      image.filter.supplier     = 'Andover surplus'
      image.filter.index        = 1.45
      END

    434.0: BEGIN
      image.line                = 434.0
      image.filter.line         = 434.0462  ; air
      image.filter.tag          = 'Dgamma'
      image.filter.cwl          = 434.01
      image.filter.fwhm         = 1.49
      image.filter.transmission = 0.4641
      image.filter.cavities     = 2
      image.filter.supplier     = 'Andover surplus'
      image.filter.index        = 1.45
      image.filter.xshape = [0.0,0.47,1.0,1.84,3.46,6.3   ,15.0  ]  ; from Andover transmission curve
      image.filter.yshape = [1.0,0.9 ,0.5,0.1 ,0.01,1.0E-3,1.0E-4]  ;  (RHS)
      END

    447.2: BEGIN   
      image.line                = 447.2
      image.filter.line         = 447.15
      image.filter.tag          = 'HeI'
      image.filter.cwl          = 447.22
      image.filter.fwhm         = 2.0
      image.filter.transmission = 0.5745
      image.filter.cavities     = 3   
      image.filter.supplier     = 'Andover custom'
      image.filter.index        = 1.45
      END

    465.0: BEGIN   
      image.line                = 465.0
      image.filter.line         = 465.20  ; This is a triplet -- not sure what I should use...
      image.filter.tag          = 'CIII'
      image.filter.cwl          = 465.2
      image.filter.fwhm         = 10.00 ; 1.42  Useless to use proper FWHM since the line is broad...
      image.filter.transmission = 0.4099
      image.filter.cavities     = 2   
      image.filter.supplier     = 'Andover surplus'
      image.filter.index        = 1.45
      END

    465.0: BEGIN   
      image.line                = 465.0
      image.filter.line         = 465.20  ; This is a triplet -- not sure what I should use...
      image.filter.tag          = 'CIII'
      image.filter.cwl          = 465.29
      image.filter.fwhm         = 10.0 ; 1.42 -- real FWHM, but a problem since this is a triplet and too broad for this filter...
      image.filter.transmission = 0.4099
      image.filter.cavities     = 2   
      image.filter.supplier     = 'Andover surplus'
      image.filter.index        = 1.45
      END

    468.7: BEGIN 
      image.line                = 468.7
      image.filter.line         = 468.56  ; an average of 4 lines from 468.537 to 468.580
      image.filter.tag          = 'HeII'
      image.filter.cwl          = 468.96     
      image.filter.fwhm         = 1.8
      image.filter.transmission = 0.5183
      image.filter.cavities     = 2         ; ???
      image.filter.supplier     = 'Andover custom'
      image.filter.index        = 1.45      ; ???
      END

    486.2: BEGIN 
      image.line                = 468.7
      image.filter.line         = 468.132918  ; average of 3 lines from NIST, weighted by strength
      image.filter.tag          = 'Dbeta'     ;    30 P   4861.2786
      image.filter.cwl          = 468.15      ;    10 P	  4861.2870
      image.filter.fwhm         = 1.64        ;    60 P	  4861.3615
      image.filter.transmission = 0.4500
      image.filter.cavities     = 3
      image.filter.supplier     = 'Andover custom'
      image.filter.index        = 2.05
      END

    514.0: BEGIN 
      image.line                = 514.0
      image.filter.line         = 514.516  ; This is a doublet? triplet?
      image.filter.tag          = 'CII'
      image.filter.cwl          = 514.20
      image.filter.fwhm         = 3.17 
      image.filter.transmission = 0.5658
      image.filter.cavities     = 2     
      image.filter.supplier     = 'Andover surplus'
      image.filter.index        = 2.05
      END

    588.1: BEGIN 
      image.line                = 588.1
      image.filter.line         = 587.56   
      image.filter.tag          = 'HeI'
      image.filter.cwl          = 588.08     
      image.filter.fwhm         = 2.0
      image.filter.transmission = 0.8010
      image.filter.cavities     = 3
      image.filter.supplier     = 'Andover custom'
      image.filter.index        = 2.0
      END

    589.0: BEGIN   ; *** BOGUS ***
      image.line                = 589.0
      image.filter.line         = 589.0   
      image.filter.tag          = 'D-D'
      image.filter.cwl          = 589.0     
      image.filter.fwhm         = 100.0     ; band emission so no blue-shift issues?
      image.filter.transmission = 0.6997
      image.filter.cavities     = 2         ; ???
      image.filter.supplier     = 'Andover surplus'
      image.filter.index        = 1.45      ; ???
      END

    601.2: BEGIN
      image.line                = 601.2
      image.filter.line         = 601.2     ; Fulcher band
      image.filter.tag          = 'D-D'
      image.filter.cwl          = 601.68    ; wrong in John's table?
      image.filter.fwhm         = 100.0     ; band emission so no blue-shift issues?
      image.filter.transmission = 0.6997
      image.filter.cavities     = 2         ; ???
      image.filter.supplier     = 'Andover surplus'
      image.filter.index        = 1.45      ; ???
      END

    656.0: BEGIN
      image.line                = 656.3
      image.filter.line         = 656.280  ; air
      image.filter.tag          = 'Dalpha'
      image.filter.cwl          = 656.70
      image.filter.fwhm         = 3.36
      image.filter.transmission = 0.6983
      image.filter.cavities     = 2 
      image.filter.supplier     = 'Andover stock'
      image.filter.index        = 2.05
      END

    667.8: BEGIN  ; 667 nm  (which one, we have 2...)
      image.line                = 667.8
      image.filter.line         = 667.81517  ; air
      image.filter.tag          = 'HeI'
      image.filter.cwl          = 667.83
      image.filter.fwhm         = 4.67
      image.filter.transmission = 0.7527
      image.filter.cavities     = 4 
      image.filter.supplier     = 'Andover surplus'
      image.filter.index        = 2.05  ; ???
      END

    668.2: BEGIN  
      image.line                = 668.2
      image.filter.line         = 667.81517  ; air
      image.filter.tag          = 'HeI'
      image.filter.cwl          = 668.21
      image.filter.fwhm         = 2.9
      image.filter.transmission = 0.8279
      image.filter.cavities     = 3 
      image.filter.supplier     = 'Andover custom'
      image.filter.index        = 2.05
      END

    706.5: BEGIN
      image.line                = 706.5
      image.filter.line         = 706.52  ; air  ; bug, was 706.68 - SL, 010607
      image.filter.tag          = 'HeI'
      image.filter.cwl          = 706.68
      image.filter.fwhm         = 1.23
      image.filter.transmission = 0.7266
      image.filter.cavities     = 2 
      image.filter.supplier     = 'Andover stock'
      image.filter.index        = 2.05
      image.filter.xshape = [0.0,0.44,1.0,1.86,3.47,6.3   ,15.0  ]  ; from Andover transmission curve
      image.filter.yshape = [1.0,0.9 ,0.5,0.1 ,0.01,1.0E-3,1.0E-4]  ;   (RHS, symmetry assumed)
      END

    728.1: BEGIN
      image.line                = 728.1
      image.filter.line         = 728.1351  ; air
      image.filter.tag          = 'HeI'
      image.filter.cwl          = 728.19
      image.filter.fwhm         =   1.25
      image.filter.transmission = 0.6678
      image.filter.cavities     = 2  
      image.filter.supplier     = 'Andover surplus'
      image.filter.index        = 2.05
      image.filter.xshape = [0.0,0.42,1.0,1.83,3.42,6.3   ,15.0  ]  ; from Andorver transmission curve
      image.filter.yshape = [1.0,0.9 ,0.5,0.1 ,0.01,1.0E-3,1.0E-4]  ;  (RHS)
      END

    734.0: BEGIN
      image.line                = 734.0 ; quasi-continuum measurement for 728 background
      image.filter.line         = 734.0
      image.filter.tag          = 'HeI'
      image.filter.cwl          = 734.22
      image.filter.fwhm         = 1000 ; 5.01 nm in reality but (hopefully) constant illumination
      image.filter.transmission = 0.7797
      image.filter.cavities     = 4
      image.filter.supplier     = 'Andover custom'
      image.filter.index        = 2.05
      END

    ELSE: BEGIN
      PRINT,'Error GetDivCamFilterData: filter not identified'
      PRINT,'  filter: ',wavelength
      STOP
      END   
  ENDCASE

; Apply generic filter transmission profile if measured profile
; not assigned above:
  IF (image.filter.yshape[0] EQ 0.0) THEN BEGIN
    supplier = STRSPLIT(image.filter.supplier,/extract)
    CASE supplier[0] OF 
      'Andover': BEGIN  ; from www.andcorp.com
        CASE image.filter.cavities OF
          2: BEGIN
            image.filter.xshape = [0.0,0.5 ,1.0,2.0 ,3.5 ,6.3   ,15.0  ]
            image.filter.yshape = [1.0,0.9 ,0.5,0.1 ,0.01,1.0E-3,1.0E-4]
            END
          3: BEGIN
            image.filter.xshape = [0.0,0.70,1.0,1.50,2.20,3.20  ,5.40  ,15.0  ]
            image.filter.yshape = [1.0,0.9 ,0.5,0.1 ,0.01,1.0E-3,1.0E-4,1.0E-5]
            END
          4: BEGIN
            image.filter.xshape = [0.0,0.85,1.0,1.25,1.65,2.25  ,4.25  ,12.0  ]
            image.filter.yshape = [1.0,0.9 ,0.5,0.1 ,0.01,1.0E-3,1.0E-4,1.0E-5]
            END
          ELSE: BEGIN
            PRINT,'Error GetDivCamFilterData: Unknown filter cavity type'
            PRINT,'  CAVITIES = ',filter.cavitites
            STOP
            END
        ENDCASE
        END
      ELSE: BEGIN
        PRINT,'Error GetDivCamFilterData: Filter supplier not recognized'
        PRINT,'  SUPPLIER = ',filter.supplier
        STOP
        END
    ENDCASE
  ENDIF

  ; Radial transmission profile which may be important for narrow filters:
  image.filter = GetFilterTransmissionProfile(image.filter,image.camera,plots)

END
;
; ======================================================================
; 
PRO AssignFilterData, image, filter, plots

  check = STRUPCASE(filter)
  IF (check EQ 'EMPTY' OR check EQ 'LNC_DALPHA' OR check EQ 'UNKNOWN') THEN BEGIN
    image.filter.tag = filter
    RETURN
  ENDIF

  PRINT,'FILTER ID=',filter,image.shot
  IF (image.shot GE 24800 AND image.shot LE 25037) THEN BEGIN
    filter = STRTRIM(filter,2)
    CASE filter OF
      '434.0/1.5/25D (AC-M116-': filter = 'Dgamma 434.01/1.49 nm'
      '486/2/25D (AC-0101-01)' : filter = 'Dbeta 486.15/1.64 nm'
      '465/1.5/25D (AC-M116-01': filter = 'CIII 465.29/1.42 nm'
      ELSE:
    ENDCASE
    PRINT,'NEW FILTER ID=',filter    
  ENDIF

  CASE image.device OF

    'MAST': BEGIN

      CASE image.camera OF

        'DIVCAM': BEGIN
          str = STRSPLIT(filter,/extract)                          
          ; Check if a neutral density filter is in series with the the filter:
          nd = [str[0],'1.0']
          IF (STRMATCH(str[0],'*/*') EQ 1) THEN nd = STRSPLIT(str[0],'/',/extract)
          ; ------------------------------------------------
          ; HACK!   *** LOOK HERE ***
          ; ------------------------------------------------
          IF (nd[1] EQ '16') THEN BEGIN
            nd[1] = '19.4'
            PRINT,' ***'
            PRINT,' *** OVER-RIDRING NEUTRAL DENSITY FILTER SCALING! ***'
            PRINT,' ***'
          ENDIF
          image.filter.nd = 1.0 / FLOAT(nd[1])
          ; Get wavelength identifier:
          wlngth = STRSPLIT(str[1],'/',/extract)                          
          PRINT,'WLNGTH[0] =',STRTRIM(STRING(wlngth[0]),2)
          CASE STRTRIM(STRING(wlngth[0]),2) OF
            '411.18': wlngth[0] = 411.2
            '434.01': wlngth[0] = 434.0
            '465.29': wlngth[0] = 465.0
            '486.15': wlngth[0] = 486.2
            '514.20': wlngth[0] = 100.0
            '656.70': wlngth[0] = 656.0
            ELSE:
          ENDCASE
          PRINT,'WLNGTH[0] REMAPPED=',wlngth[0]

          GetDivCamFilterData, image, FLOAT(wlngth[0]), plots
          END

        'FFC': BEGIN
          str = STRSPLIT(filter,/extract)                          
          ; Check if a neutral density filter is in series with the the filter:
          nd = [str[0],'1.0']
          IF (STRMATCH(str[0],'*/*') EQ 1) THEN nd = STRSPLIT(str[0],'/',/extract)
          ; ------------------------------------------------
          ; HACK!   *** LOOK HERE ***
          ; ------------------------------------------------
          IF (nd[1] EQ '16') THEN BEGIN
            nd[1] = '19.4'
            PRINT,' ***'
            PRINT,' *** OVER-RIDRING NEUTRAL DENSITY FILTER SCALING! ***'
            PRINT,' ***'
          ENDIF
          image.filter.nd = 1.0 / FLOAT(nd[1])
          ; Get wavelength identifier:
          wlngth = STRSPLIT(str[1],'/',/extract)                          
          PRINT,'WLNGTH[0] =',STRTRIM(STRING(wlngth[0]),2)
          CASE STRTRIM(STRING(wlngth[0]),2) OF
            '434.01': wlngth[0] = 434.0
            '465.29': wlngth[0] = 465.0
            '514.20': wlngth[0] = 514.0
            '656.70': wlngth[0] = 656.0
            ELSE:
          ENDCASE
          PRINT,'WLNGTH[0] REMAPPED=',wlngth[0]
          GetDivCamFilterData, image, FLOAT(wlngth[0]), plots
          END

        'ZEBRA': image.filter.tag = filter
        'RGB'  : image.filter.tag = filter
        'MWIR' : image.filter.tag = filter
        'LWIR' : image.filter.tag = filter

        ELSE: BEGIN
          PRINT,'Error AssignFilterData: Camera not identified'
          PRINT,'  DEVICE = ',image.camera
          STOP
          END    

      ENDCASE
      END

    'FILE': BEGIN
      IF (image.format EQ 'ipx') THEN BEGIN
        str = STRSPLIT(filter,/extract)                          

        ; Check if a neutral density filter is in series with the the filter:
        nd = [str[0],'1.0']
        IF (STRMATCH(str[0],'*/*') EQ 1) THEN nd = STRSPLIT(str[0],'/',/extract)
        image.filter.nd = 1.0 / FLOAT(nd[1])

        ; Get wavelength identifier:
        wlngth = STRSPLIT(str[1],'/',/extract)                          

        GetDivCamFilterData, image, FLOAT(wlngth[0])
      ENDIF
      END

    ELSE: BEGIN
      PRINT,'Error AssignFilterData: device not identified'
      PRINT,'  DEVICE = ',image.device
      STOP
      END    
  ENDCASE

END

;
;
;
;
;

