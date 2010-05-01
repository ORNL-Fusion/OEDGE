;@getimage_make

; Haven't checked if the images shoudl be shifted by the same amount as in 20100318
; 

PRO batch_20100413,option,plots=plots,scale=scale

  CASE option OF
   1: BEGIN
     maskbox  = [271,125,287,145,349,130,364,153]
     maskpoly = [14, 0,49,180,47,249,53,254,62,289,65,336,61,351,83,381,88,386,176,430,187,431,218,495,268,495,0,0,0]
     image=getimage(shot=25028,frame=673,channel=3,camera='FFC',window='HL01',               $
                    plots=plots,rotate=270,colour=5,shift=[21,19],/clean,/save,              $
                    maskleft=60,maskbottom=260,maskpoly=maskpoly,maskbox=maskbox,scale=scale)
     END
   2: BEGIN
     maskpoly=[12, 0,55,284,55,294,70,327,70,350,63,377,71,387,92,416,92,417,180,499,302,499,0,0,0]
     image=getimage(shot=25028,frame=15,channel=2,plots=plots,colour=5,rotate=270,shift=[20,-35],  $
                    /clean,/save,maskleft=125,maskbottom=250,maskpoly=maskpoly,scale=scale)
     END
   3: BEGIN
     maskpoly = [14, 0,49,146,49,263,57,329,66,332,73,362,76,366,85,392,92,394,175,424,178,461,218,510,241,511,0,0,0]
     image=getimage(shot=25028,frame=673,channel=1,camera='FFC',window='HL07',  $
                    plots=plots,shift=[43,36],/save,/clean,scale=scale,         $
                    maskleft=90,maskbottom=260,maskpoly=maskpoly)
     END
   4: BEGIN
     maskpoly = [12, 0,90,154,83,291,91,384,109,389,119,413,124,415,207,458,212,463,249,499,271,499,0,0,0]
     maskbox  = [294,169,308,193,367,176,380,198]
     image=getimage(shot=25028,frame=15,channel=4,plots=plots,colour=5,shift=[0,-63],  $
                    /clean,/save,scale=scale,                                          $
                    maskleft=115,maskbottom=300,maskbox=maskbox,maskpoly=maskpoly)     
     END
  ENDCASE


END


