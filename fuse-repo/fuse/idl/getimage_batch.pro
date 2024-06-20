@getimage_make

plots=0


maskbox = [282,128,300,155,363,138,378,159]
maskpoly=[14, 0,58,259,58,271,67,303,70,315,62,354,68,365,87,396,87,400,182,448,182,448,220,479,220,479,0,0,0]
image=getimage(shot=24861,frame=343,channel=3,camera='FFC',window='HL01',plots=plots,rotate=270,colour=5,shift=[21,19],/clean,/save,maskleft=65,maskbottom=280,maskpoly=maskpoly,maskbox=maskbox)
image=getimage(shot=24861,frame=762,channel=3,camera='FFC',window='HL01',plots=plots,rotate=270,colour=5,shift=[21,19],/clean,/save,maskleft=65,maskbottom=280,maskpoly=maskpoly,maskbox=maskbox)
image=getimage(shot=24869,frame=552,channel=3,camera='FFC',window='HL01',plots=plots,rotate=270,colour=5,shift=[21,19],/clean,/save,maskleft=65,maskbottom=280,maskpoly=maskpoly,maskbox=maskbox)
image=getimage(shot=24869,frame=662,channel=3,camera='FFC',window='HL01',plots=plots,rotate=270,colour=5,shift=[21,19],/clean,/save,maskleft=65,maskbottom=280,maskpoly=maskpoly,maskbox=maskbox)


maskpoly=[12, 0,55,284,55,294,70,327,70,350,63,377,71,387,92,416,92,417,180,499,302,499,0,0,0]
image=getimage(shot=24861,frame=11,channel=2,plots=plots,colour=5,rotate=270,shift=[20,-35],/clean,/save,maskleft=125,maskbottom=250,maskpoly=maskpoly)
image=getimage(shot=24861,frame=17,channel=2,plots=plots,colour=5,rotate=270,shift=[20,-35],/clean,/save,maskleft=125,maskbottom=250,maskpoly=maskpoly)
image=getimage(shot=24869,frame=14,channel=2,plots=plots,colour=5,rotate=270,shift=[20,-35],/clean,/save,maskleft=125,maskbottom=250,maskpoly=maskpoly)
image=getimage(shot=24869,frame=15,channel=2,plots=plots,colour=5,rotate=270,shift=[20,-35],/clean,/save,maskleft=125,maskbottom=250,maskpoly=maskpoly)


maskpoly = [11, 0,10,173,10,367,44,370,52,396,56,399,139,475,149,475,447,512,447,511,0,0,0]
image=getimage(shot=24861,frame=343,channel=1,camera='FFC',window='HL07',plots=plots,shift=[43,36],/save,/clean,maskleft=100,maskbottom=225,maskpoly=maskpoly)
image=getimage(shot=24861,frame=762,channel=1,camera='FFC',window='HL07',plots=plots,shift=[43,36],/save,/clean,maskleft=100,maskbottom=225,maskpoly=maskpoly)
image=getimage(shot=24869,frame=552,channel=1,camera='FFC',window='HL07',plots=plots,shift=[43,36],/save,/clean,maskleft=100,maskbottom=225,maskpoly=maskpoly)
image=getimage(shot=24869,frame=662,channel=1,camera='FFC',window='HL07',plots=plots,shift=[43,36],/save,/clean,maskleft=100,maskbottom=225,maskpoly=maskpoly)


maskpoly=[12, 0,90,154,83,291,91,384,109,389,119,413,124,415,207,458,212,463,249,499,271,499,0,0,0]
maskbox = [294,169,308,193,367,176,380,198]
image=getimage(shot=24861,frame=11,channel=4,plots=plots,colour=5,shift=[0,-63],/clean,/save,maskleft=115,maskbottom=310,maskbox=maskbox,maskpoly=maskpoly)
image=getimage(shot=24861,frame=17,channel=4,plots=plots,colour=5,shift=[0,-63],/clean,/save,maskleft=115,maskbottom=310,maskbox=maskbox,maskpoly=maskpoly)
image=getimage(shot=24869,frame=14,channel=4,plots=plots,colour=5,shift=[0,-63],/clean,/save,maskleft=115,maskbottom=310,maskbox=maskbox,maskpoly=maskpoly)
image=getimage(shot=24869,frame=15,channel=4,plots=plots,colour=5,shift=[0,-63],/clean,/save,maskleft=115,maskbottom=310,maskbox=maskbox,maskpoly=maskpoly)
