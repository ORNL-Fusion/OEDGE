;
; a=ray(camera='FFC',shot='24860',frame='760',channel='3',region=2,/sp_peak,fit_cutoff=0.5,fit_sample=100,spt=[0.280,-1.435])
;
; file_a='FFC_24860_760_3_HL01_rbc_Dalpha.cgm'
; file_g='FFC_24860_760_1_HL07_rba_Dgamma.cgm'
;
; a=ray(file=file_a,shot=24860,region=2,fit_cutoff=0.5,fit_sample=50,spt=[0.280,-1.435],xpt=[0.67,-1.13],/plots)
; g=ray(file=file_g,shot=24860,region=2,fit_cutoff=0.5,fit_sample=10,spt=[0.280,-1.435],xpt=[0.67,-1.13],/plots)
; plot,SMOOTH(g.contour.fit_data[10:*],5)*0.05/SMOOTH(a.contour.fit_data[10:*],5),color=Truecolor('Black')
;
; a=ray(file=file_a,shot=24860,region=2,fit_cutoff=0.5,fit_sample=10,spt=[0.280,-1.435],xpt=[0.67,-1.13],/plots)
; g=ray(file=file_g,shot=24860,region=2,fit_cutoff=0.5,fit_sample=10,spt=[0.280,-1.423],xpt=[0.67,-1.13],/plots)
; plot,SMOOTH(g.contour.fit_data,5)*0.05/SMOOTH(a.contour.fit_data,5),color=Truecolor('Black'),yrange=[0.0,0.2]
;
; file_a='FFC_24860_342_3_HL01_rbc_Dalpha.cgm'
; a=ray(file=file_a,shot=24860,region=2,fit_cutoff=0.5,fit_sample=10,spt=[0.280,-1.55],xpt=[0.64,-1.25],/plots)
;
; file_g='FFC_24860_342_1_HL07_rba_Dgamma.cgm'
; g=ray(file=file_g,shot=24860,region=2,fit_cutoff=0.5,fit_sample=10,spt=[0.280,-1.55],xpt=[0.64,-1.25],/plots)
;
; plot,SMOOTH(g.contour.fit_data,5)*0.05/SMOOTH(a.contour.fit_data,5),color=Truecolor('Black')
;
;
; ======================================================================
;
;FUNCTION LoadImage, filename
;
;  path='~slisgo/divimp/shots/uls/images/'
;
;  filename = path + filename + '.sav'
;
;  PRINT,'Restoring ',filename
;
;  RESTORE, filename
;
;  RETURN, image_data
;
;END
;
; ======================================================================
;
;FUNCTION LoadReconstruction, file,ext=ext,path=path
;FUNCTION LoadReconstruction, camera=camera,shot=shot,frame=frame,channel=channel,ext=ext,path=path
FUNCTION LoadReconstruction, shot,frame,channel,camera=camera,ext=ext,path=path


;  result = GetInversion(file,channel=channel,ext=ext,path=path,/full)
                        

result = GetInversion(camera,shot,frame,channel,  $
                      ext=ext,path=path,/full)

;  inv = GetInversion(shot=shot,channel=channel,frame=frame,   $
;                     ext=ext,casename=filename,region=region, $
;                     datatag=datatag,/full)

  RETURN, result

END
;
; ======================================================================
;

FUNCTION ray,  $
  camera=camera,   $ ;
  shot=shot,       $ ;
  frame=frame,     $ ;
  time=time,       $ ;
  channel=channel, $ ;
  suffix=suffix  , $ ;
  ext=ext        , $ ;
  full=full      , $ ;
  region =region , $ ;
  sp_peak=sp_peak, $ ;
  fit_sample = fit_sample, $ ;
  fit_cutoff = fit_cutoff, $ ;
  progress=progress, $ ;
  file   =file   , $ ;
  xpt    =xpt    , $ ;
  spt    =spt    , $ ;
  plots  =plots  , $ ;
  ifile=ifile,     $ ;
  refresh=refresh    ;

  IF (NOT KEYWORD_SET(camera )) THEN camera  = 'DIVCAM'
  IF (NOT KEYWORD_SET(ext    )) THEN ext     = 'cgm'
  IF (NOT KEYWORD_SET(nsample)) THEN nsample = 20
  IF (NOT KEYWORD_SET(path   )) THEN path = '~/fuse_data/mast/images/'

  print,'file= ',file
  print,'path= ',path

  CASE 1 OF
    (KEYWORD_SET(shot) AND KEYWORD_SET(frame) AND KEYWORD_SET(channel)): BEGIN
      path = '~/fuse_data/mast/images/'
      file_name = path + camera + '_' +                   $
                  STRTRIM(STRING(shot   ),2) + '_' +  $
                  STRTRIM(STRING(frame  ),2) + '_' +  $
                  STRTRIM(STRING(channel),2) 
      IF (KEYWORD_SET(suffix)) THEN file_name = file_name + '_' + STRTRIM(suffix,2)
      file_name = file_name + '.' + STRTRIM(ext,2)
      END
    (KEYWORD_SET(file)): file_name = path + file 
    ELSE: BEGIN
      PRINT,'ERROR ray: Unsufficient input data to know what to do'
      RETURN,-1
      END
  ENDCASE
    
  print,'file= ',file_name

  result = GetInversion(file_name,shot=shot,frame=frame,channel=channel,region=region,full=full)
;  result = GetInversion('FFC',24860,760,3,region=2)
  result = ContourImage(result, sp_peak=sp_peak,spt=spt,xpt=xpt,plots=plots)
  result = ProcessImage(result, 100.0, fit_sample=fit_sample, fit_cutoff=fit_cutoff, progress=progress,plots=plots)

  return, result


; Check if image has already been processed, or force the image 
; to be processed if the 'refresh' keyword is set:

  home = '/home/slisgo/idl/ray/'

  machine = 'MAST'
  shot    = 18472
  frame   = 29
  camera  = 'ZEBRA'
  channel = 1

  case_name = 'm-zebra'

;  PRINT,shot,frame,'  '+camera,channel


;  Check if image has already been processed / calibrated, or reprocess
;  if REFRESH_IMAGE is set:
;  .r ../dc/dc_load ../dc/dc_save ../dc/dc_utility ../dc/dc_calibrate ../dc/dc_filter ../dc/dc
;  image=GetImage(shot=shot,frame=frame,channel=channel,camera=camera,  $
;                 maskradius=62,maskleft=5,/clean,/save,path=path)

; Run RAY script that calls the external programme:
; ----------------------------------------------------------------------

;  Check if the geometry matrix has already been calcualted, or whether 
;  REFRESH_MATRIX is set:

  image_name =  $
    camera                     + '_' +  $
    STRTRIM(STRING(shot   ),1) + '_' +  $
    STRTRIM(STRING(frame  ),1) + '_' +  $
    STRTRIM(STRING(channel),1) + '.idl'

  command =  home +  $
    'scripts/run_ray'          + ' ' +  $
    machine                    + ' ' +  $
    STRTRIM(STRING(shot),1)    + ' ' +  $
    STRTRIM(STRING(frame),1)   + ' ' +  $
    camera                     + ' ' +  $
    STRTRIM(STRING(channel),1) + ' ' +  $
    case_name                  + ' ' +  $
    image_name

  PRINT, command

  SPAWN, command

  path = home + 'images/'

  a = LoadReconstruction(camera,shot,frame,channel,ext='cgm',path=path)

return,a
stop




  RETURN, -1

END
;
;
;
;
; ======================================================================
;

FUNCTION ray_new,  $
  camera=camera,   $ ;
  shot=shot,       $ ;
  frame=frame,     $ ;
  time=time,       $ ;
  channel=channel, $ ;
  suffix=suffix  , $ ;
  ext=ext        , $ ;
  region =region , $ ;
  sp_peak=sp_peak, $ ;
  fit_sample = fit_sample, $ ;
  fit_cutoff = fit_cutoff, $ ;
  progress=progress, $ ;
  file   =file   , $ ;
  xpt    =xpt    , $ ;
  spt    =spt    , $ ;
  plots  =plots  , $ ;
  ifile=ifile,     $ ;
  refresh=refresh    ;

  IF (NOT KEYWORD_SET(camera )) THEN camera  = 'DIVCAM'
  IF (NOT KEYWORD_SET(ext    )) THEN ext     = 'cgm'
  IF (NOT KEYWORD_SET(nsample)) THEN nsample = 20
  IF (NOT KEYWORD_SET(path   )) THEN path = '~/fuse_data/mast/images/'

  print,'file= ',file
  print,'path= ',path

  CASE 1 OF
    (KEYWORD_SET(shot) AND KEYWORD_SET(frame) AND KEYWORD_SET(channel)): BEGIN
      path = '~/fuse_data/mast/images/'
      file_name = path + camera + '_' +                   $
                  STRTRIM(STRING(shot   ),2) + '_' +  $
                  STRTRIM(STRING(frame  ),2) + '_' +  $
                  STRTRIM(STRING(channel),2) 
      IF (KEYWORD_SET(suffix)) THEN file_name = file_name + '_' + STRTRIM(suffix,2)
      file_name = file_name + '.' + STRTRIM(ext,2)
      END
    (KEYWORD_SET(file)): file_name = path + file 
    ELSE: BEGIN
      PRINT,'ERROR ray: Unsufficient input data to know what to do'
      RETURN,-1
      END
  ENDCASE
;    
; Run RAY script that calls the external programme:
; ----------------------------------------------------------------------

;  Check if the geometry matrix has already been calcualted, or whether 
;  REFRESH_MATRIX is set:

   command =  home +  $
    'scripts/run_ray'          + ' ' +  $
    machine                    + ' ' +  $
    STRTRIM(STRING(shot),1)    + ' ' +  $
    STRTRIM(STRING(frame),1)   + ' ' +  $
    camera                     + ' ' +  $
    STRTRIM(STRING(channel),1) + ' ' +  $
    case_name                  + ' ' +  $
    image_name

  PRINT, command

  SPAWN, command

  path = home + 'images/'

  a = LoadReconstruction(camera,shot,frame,channel,ext='cgm',path=path)

return,a
stop




  RETURN, -1

END
;
;
;