PRO getimage_psi2010_fast_divertor, plots=plots


  FOR frame = 2300, 2500 DO BEGIN
    image=getimage(shot=22112,frame=frame,channel=1,camera='FFC',plots=plots,  $
                   /nocal,filter='empty',colour=3,scale=2.0,/no_square,/save)
  ENDFOR

END

PRO getimage_psi2010_main_chamber, plots=plots 

  FOR frame = 1, 400 DO BEGIN
    image=getimage(shot=15622,frame=frame,channel=2,camera='FFC',plots=plots,  $
                   /nocal,filter='empty',colour=3,scale=2.0,/no_square,/save,/kill)
  ENDFOR

END