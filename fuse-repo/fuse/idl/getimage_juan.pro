

;@getimage_make

; Haven't checked if the images shoudl be shifted by the same amount as in 20100318
; 

PRO getimage_juan,option,plots=plots,scale=scale

  CASE option OF
   1: BEGIN    ; ref                     o/b puff    density ramp
     l_shots  = [22426]
     l_frames = [  150]

     channel = 1
     camera  = 'FFC'
     window  = 'HM07'
     filter  = 'Dalpha 656.0/10.0'
     nocal   = 1

;     rotate  = 270
;     shift   = [21,19]

     masktop  = 100
     maskright = 450
     maskpoly = [17, 0  ,296,50 ,296,106,83 ,234,48 ,268,86 ,291,86 ,291,35 ,450,22 ,529,24 ,529, 71, $
                     559,71 ,570,38 ,820,111,860,296,900,296,900,0  ,0  ,0  ,  $
                 21, 0  ,296,50 ,296,108,544,237,581,267,539,292,542,292,594,465,605,537,601,537,548,  $
                     563,548,576,586,718,558,721,545,810,530,810,446,850,446,860,296,900,296,900,900,0,900]
     END
   2: BEGIN
     l_shots  = [24861,24861,24861,      24862,24869,24866,24866,24866]
     l_frames = [   11,   14,   17,         14,   14,   17,   19,   22]
     l_shots  = [24866]
     l_frames = [   14]

     channel = 2
     rotate  = 270
     shift   = [20,-35]

     maskleft   = 125
     maskbottom = 250
     maskpoly   = [12, 0,55,284,55,294,70,327,70,350,63,377,71,387,92,416,92,417,180,499,302,499,0,0,0]

;     FOR i = 0, N_ELEMENTS(l_shots)-1 DO  $
;       image=getimage(shot=l_shots[i],frame=l_frames[i],channel=2,  $
;                      rotate=270,shift=[20,-35],scale=scale,    $
;                      /clean,/save,plots=plots,colour=5,        $
;                      maskleft=125,maskbottom=250,maskpoly=maskpoly)
     END
   3: BEGIN    ; ref                     o/b puff    density ramp
     l_shots  = [24861,24861,24861,24861,24862,24869,24866,24866,24866]
     l_frames = [  343,  553,  762, 1108,  554,  552,  766,  903, 1109]
;     l_shots  = [24867,24868]  ; attached plasma, inter-ELM H-mode
;     l_frames = [ 1225, 1225] 
     l_shots  = [24866]
     l_frames = [  553]

     channel = 1
     camera  = 'FFC'
     window  = 'HL07'
     shift   = [43,36]

     maskleft   = 100
     maskbottom = 225
     maskpoly = [11, 0,10,173,10,367,44,370,52,396,56,399,139,475,149,475,447,512,447,511,0,0,0]

;     FOR i = 0, N_ELEMENTS(l_shots)-1 DO  $
;       image=getimage(shot=l_shots[i],frame=l_frames[i],channel=1,camera='FFC',window='HL07',  $
;                      shift=[43,36],scale=scale,                                  $
;                      /save,/clean,plots=plots,colour=5,                          $
;                      maskleft=100,maskbottom=225,maskpoly=maskpoly)
     END

  ENDCASE


  FOR i = 0, N_ELEMENTS(l_shots)-1 DO  $
    image=getimage(shot=l_shots[i],frame=l_frames[i],channel=channel,  $
                   camera=camera,window=window,filter=filter,nocal=nocal,  $
                   rotate=rotate,shift=shift,scale=scale,              $
                   /clean,/save,plots=plots,colour=5,                  $
                   maskleft=maskleft,masktop=masktop,maskbottom=maskbottom,  $
                   maskright=maskright,maskpoly=maskpoly,maskbox=maskbox)

END

PRO getimage_juan_batch,option,plots=plots,scale=scale
 
  getimage_juan,1,plots=plots,scale=scale

END

