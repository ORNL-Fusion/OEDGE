; Function for IPX file access
; S. Shibaev Aug 2005
;
; Returns a descriptor for IPX file or 0 in case of any error
; Usage:
;   desc = ipx_open(IPX_file_name,RATE=rate)
; RATE is the maximum frame rate in Hz
; The descriptor is used to get a frame: image = ipx_frame(desc,frame_number)
; All camera parameters can be obtained from the descriptor returned:
;     help,/st,desc
;
function ipx_open, filename, RATE=rate
  if(n_params() lt 1) then begin
    print,'Opens IPX file, usage: desc=ipx_open(filename,RATE=rate)'
    print,'rate is the maximum frame rate in the file'
    return,0
  end
  filename=expand_path(filename)
  openr,unit,filename,/GET_LUN,ERROR=err
  if(err ne 0) then begin
    printf,-2,!ERR_STRING
    return,0
  end
  header={header, $
    id:          string(' ',FORMAT='(A8)'), $
    size:        0UL, $
    codec:       string(' ',FORMAT='(A8)'), $
    date_time:   string(' ',FORMAT='(A20)'), $
    shot:        0L,  $
    trigger:     0E,  $
    lens:        string(' ',FORMAT='(A24)'), $
    filter:      string(' ',FORMAT='(A24)'), $
    view:        string(' ',FORMAT='(A64)'), $
    numFrames:   0UL, $
    camera:      string(' ',FORMAT='(A64)'), $
    width:       0U,  $
    height:      0U,  $
    depth:       0U,  $
    orient:      0UL, $
    taps:        0U,  $
    color:       0U,  $
    hBin:        0U,  $
    left:        0U,  $
    right:       0U,  $
    vBin:        0U,  $
    top:         0U,  $
    bottom:      0U,  $
    offset:      intarr(2), $
    gain:        fltarr(2), $
    preExp:      0UL, $
    exposure:    0UL, $
    strobe:      0UL, $
    board_temp:  0E,  $
    ccd_temp:    0E   $
  }
  readu,unit,header
  if((strpos(header.id,'IPX 01') ne 0)||(header.numFrames lt 1)) then begin
    printf,-2,filename+' is not an IPX file'
    free_lun,unit
    return,0
  end
  if(header.size lt 286) then begin
    ccd_temp=0
  end
  framepos=ulonarr(header.numFrames)
  frametime=dblarr(header.numFrames)
  frame={frame,   $
    fsize:   0UL, $
    ftime:   0D   $
  }
  dt=1e10
  pos=header.size
  for i=0,header.numFrames-1 do begin
    framepos[i]=pos
    point_lun,unit,pos
    readu,unit,frame
    pos=pos+frame.fsize
    frametime[i]=frame.ftime
    if(i gt 0) then begin
      if(frametime[i]-frametime[i-1]) then dt=frametime[i]-frametime[i-1]
    end
  end
  free_lun,unit
  if(arg_present(rate)) then begin
    rate=0L
    if(header.numFrames gt 1) then begin
      if(dt gt 0) then rate=round(1.0/dt)
    end
  end
  desc = {              $   
    filename:filename,  $
    header:header,      $
    framepos:framepos,  $
    frametime:frametime $
  }
  return,desc
end
