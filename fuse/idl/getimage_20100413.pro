;@getimage_make

; Haven't checked if the images shoudl be shifted by the same amount as in 20100318
; 

PRO getimage_20100413,option,plots=plots,scale=scale

  CASE option OF
   1: BEGIN
     l_shots  = [25028,25028,25029,25029,25029,25029]
     l_frames = [  603,  673,  392,  603,  810, 1158]
     l_shots  = [25028,25028,25028]
     l_frames = [  392,  810, 1158]

     channel = 3
     camera  = 'FFC'
     window  = 'HL01'
     rotate  = 270
     shift   = [21,19]

     maskleft   = 60
     maskbottom = 260
     maskbox    = [271,125,287,145,349,130,364,153]
     maskpoly   = [14, 0,49,180,47,249,53,254,62,289,65,336,61,351,83,381,88,386,176,430,187,431,218,495,268,495,0,0,0]

;     FOR i = 0, N_ELEMENTS(l_shots)-1 DO  $
;       image=getimage(shot=l_shots[i],frame=l_frames[i],channel=3,camera='FFC',window='HL01',  $
;                      rotate=270,shift=[21,19],scale=scale,                                $
;                      /clean,/save,plots=plots,colour=5,                                   $
;                      maskleft=60,maskbottom=260,maskpoly=maskpoly,maskbox=maskbox)
     END
   2: BEGIN
     l_shots  = [25028,25028,25029,25029,25029]
     l_frames = [   14,   15,   11,   14,   17]
     l_shots  = [25028,25028]
     l_frames = [   11,   17]

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
   3: BEGIN
     maskpoly = [14, 0,49,146,49,263,57,329,66,332,73,362,76,366,85,392,  $
                    92,394,175,424,178,461,218,510,241,511,0,0,0]

     l_shots  = [25028,25028,25029,25029,25029,25029]
     l_frames = [  603,  673,  392,  603,  810, 1158]
     l_shots  = [25028,25028,25028]
     l_frames = [  392,  810, 1158]

             ; 25028, 673 has a strong RHS artifact, and some strangeness near the target
	     ; choose different frame perhaps

     channel = 1
     camera  = 'FFC'
     window  = 'HL07'
     shift   = [43,36]

     maskleft   = 90
     maskbottom = 260

;     FOR i = 0, N_ELEMENTS(l_shots)-1 DO  $
;       image=getimage(shot=l_shots[i],frame=l_frames[i],channel=1,camera='FFC',window='HL07',  $
;                      shift=[43,36],scale=scale,                                  $
;                      /save,/clean,plots=plots,colour=5,                          $
;                      maskleft=90,maskbottom=260,maskpoly=maskpoly)
     END
   4: BEGIN
     l_shots  = [25028,25028,25029,25029,25029]
     l_frames = [   14,   15,   11,   14,   17]
     l_shots  = [25028,25028]
     l_frames = [   11,   17]

     channel = 4
     shift   = [0,-63]

     maskleft   = 115
     maskbottom = 300
     maskpoly   = [12, 0,90,154,83,291,91,384,109,389,119,413,124,415,207,458,212,463,249,499,271,499,0,0,0]
     maskbox    = [294,169,308,193,367,176,380,198]

;     FOR i = 0, N_ELEMENTS(l_shots)-1 DO  $
;       image=getimage(shot=l_shots[i],frame=l_frames[i],channel=4,  $
;                      shift=[0,-63],scale=scale,              $
;                      /clean,/save,plots=plots,colour=5,      $
;                      maskleft=115,maskbottom=300,maskbox=maskbox,maskpoly=maskpoly)     
     END
  ENDCASE

  FOR i = 0, N_ELEMENTS(l_shots)-1 DO  $
    image=getimage(shot=l_shots[i],frame=l_frames[i],channel=channel,  $
                   camera=camera,window=window,                        $
                   rotate=rotate,shift=shift,scale=scale,              $
                   /clean,/save,plots=plots,colour=5,                  $
                   maskleft=maskleft,maskbottom=maskbottom,maskpoly=maskpoly,maskbox=maskbox)
END


PRO getimage_20100413_batch,option,plots=plots,scale=scale
 
  getimage_20100413,1,plots=plots,scale=scale
  getimage_20100413,2,plots=plots,scale=scale
  getimage_20100413,3,plots=plots,scale=scale
  getimage_20100413,4,plots=plots,scale=scale

END

;
; ----------------------------------------------------------------------
; James:
;
; --- 24861 ---
;
; 0.200,0.209, 0.214,0.217, 0.221,0.227, 0.231,0.235, 0.239,0.245, 0.248,0.253,
; 0.259,0.262, 0.267,0.274, 0.278,0.284, 0.288,0.295, 0.312,0.317, 0.321,0.328,
; 0.331,0.338, 0.342,0.348
;
; 201, (215), (229),  243 , (257),  271 ,  285 , (298), (312),  326 , (340)
;
; rbc  67 us  343 0.158673   554 0.200873  762 0.242606  1109 0.311873
; rdb   2 ms   11 0.158725    14 0.200648   17 0.242583
; rba 200 us  343 0.158673   553 0.200806  762 0.242606  1108 0.311806
; rdd 600 us   11 0.158725    14 0.200648   17 0.242583
;
; times=[0.231,0.235, 0.239,0.245, 0.248,0.253]
; lp=getlp([24861],times,[2,8])
; lp = getlp_Plots(lp,/smooth)
;
; --- 25028 ---
;
; 0.200,0.209, 0.211,0.216, 0.221,0.226, 0.229,0.233, 0.235,0.241, 0.244,0.251, 
; 0.253,0.260, 0.263,0.267, 0.270,0.275, 0.279,0.284, 0.286,0.295, 0.298,0.305,
; 0.308,0.314, 0.318,0.324, 0.328,0.334, 0.338,0.344  AIM_DA/TO10
;
; 201,  215 , (229), (243),  257 , (271), (285), (298),  312 , (326), (340)
;
; rbc   50 us  603 0.200656  673 0.214773  
; rdb    3 ms   14 0.200648   15 0.214647
; rba  167 us  603 0.200773  673 0.214656
; rdd  200 us   14 0.200648   15 0.214647	
;
; --- 25029 ---
;
; Repeat of 25028, but with C filters.  Dalpha from AIM_DA/TO10 decreasing from
; 0.220-0.260 -- due to the drop in lower i/b gas puffing?
;
; 0.200,0.210, 0.213,0.218, 0.221,0.226, 0.236,0.242, 0.252,0.257, 0.261,0.265,
; 0.268,0.273, 0.276,0.282, 0.285,0.291, 0.294,0.301, 0.304,0.311, 0.314,0.321, 
; 0.325,0.331, 0.334,0.341, 0.345,0.351 
;
; 201, (215), (229), (243),  257 ,  271 , (285),  298 , {312}, (326),  340
;
; rbc 200 us  392 0.158606   603 0.200806  810 0.242206  1158 0.311806
; rdb   3 ms   11 0.158725    14 0.200648   17 0.242583
; rba 200 us  392 0.158673   603 0.200806  810 0.242206  1158 0.311806
; rdd 200 us   11 0.158725    14 0.200648   17 0.242583
;
; --- 25030 ---
;
; 0.200,0.208, 0.210,0.216, 0.219,0.225, 0.228,0.232, 0.236,0.241, 0.224,0.250,
; 0.253,0.256, 0.260,0.265, 0.268,0.272, 0.276,0.279, 0.282,0.288, 0.294,0.300,
; 0.305,0.311, 0.315,0.321, 0.325,0.331, 0.335,0.341, 0.344,0.351
;
; rdb (3 ms integration), rdd (200 us) 
; 201,  215 , (229), (243), (257),  271 ,  285 ,  298 , {312}, (326),  340
;
;
; --- 24862 and 24869 ---
; 
; Outboard fuelling, H-mode from 250 ms or so, but good match to 24861 before that.
; CII and CIII, Dbeta and Ddelta on 24869.  24869 goes into H-mode a little better than
; 24862, but some good ELMs on both I think.  May be some dithering as early as 210 ms.  
; Sawtooth driven ELMs?  Doesn't seem so, since no real sawtooth activity prior (or 
; after) transition? 
; 
; Interestingly, and in constrast to 24867 and 24868, the plasma seems to return to the semi-detached
; state between ELMs, in particular 24862 where the ELMs are smaller.  Rather intesting movie between 320 and 330 ms, 
; where there's evolution of the inner leg in CIII inter-ELM, but also a brief period of x-point MARFE, and ELMs.
; Need to have a more careful look at it.
;
; 24862 and 24869 are excellent repeats around 201 ms.
;
; 24862
;
; rbc   67 us  554 0.200873 
; rdb    3 ms   14 0.200648 
; rba  200 us  554 0.201006 
; rdd  300 us   14 0.200648 
;
; 24869
;
; 0.200,0.208, 0.211,0.216, 0.220,0.227, 0.229,0.235, 0.237,0.245, 0.249,0.255, (small events except for the last one)
; 0.259,0.275, 0.279,0.284, 0.288,0.292, 0.296,0.300, 0.304,0.307,     H-mode
; 0.312,0.319, 0.324,0.332, 0.337,0.343,                               dithering, H-L transition at 0.309
; 0.346,0.350
;
; 201,  215 , (229),  243 , (257),  271, (285),  298 , (312),  326 ,  340
;
; rdc 200 us  552 0.200606  662 0.214606 
; rbd   3 ms   14 0.200644   15 0.214624 
; rba 200 us  552 0.200606  662 0.214606 (frames near there a bit cleaner)
; rdd 600 us   14 0.200644   15 0.214624 
;
; --- 24866 ---
;
; Density ramp, but not ane_density data, sadly.  Ip flat-top ends at 320 ms.  Da/Dg is lower than 
; reference after 250 ms.  Target profile is identicle to 24861 up to (at least) 210 ms.
;
; 0.200,0.210, 0.230,0.215, 0.220,0.226, 0.230,0.234, 0.237,0.243, 0.248,0.252,
; 0.256,0.261, 0.265,0.268, 0.273,0.278, 0.280,0.286, 0.289,0.295, 0.299,0.304,
; 0.308,0.312, 0.315,0.320
;
; 201,  215 , (229),  243 , (257),  271,   285 , (298),  312 
;
; rbc   50 us  766 0.243256   903 0.270656   1109 0.312056
; rdb    1 ms   17 0.242605    19 0.270562     22 0.312498
; rba  100 us  766 0.243306   903 0.270706   1109 0.311906
; rdd  200 us   17 0.242605    19 0.270562     22 0.312498
;
; --- 24867 and 24868 ---
;
; Used some frames for relative calibration under attached plasma conditions... not liste here, yet.
;
; 24867: rbc (100 us), rba(200 us)
; 24868:
;
; Lower density, H-mode, C filters on 24868, poor match to anything, but perhaps a good reference attached shot
;
; 24868
;
; 0.200,0.219, 0.223,0.236, 0.239,0.250, 0.250,0.257, 0.261,0.267, 0.271,0.278,
; 0.282,0.290, 0.294,0.312, 0.318,0.323, 0.328,0.336, 0.343,0.350
;
; rbc ( ) , rba ( )
; rdb (2 ms integration), rdd (600 us)
; 201,  215 ,  229 ,  243 ,  257 , (271),  285 ,  298 ,  312 , (326), (340)
;



