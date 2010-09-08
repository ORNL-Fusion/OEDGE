; Function for IPX file access
; S. Shibaev Aug 2005
; Revisions:
;   2 - Reading of compressed IPX files (JP2 and JPC), Nov 2005
;
; Returns a frame pixel array or 0 in case of any error
; Usage: 
;   image = ipx_frame(desc,frame_number,TIME=frame_time)
; or
;   image = ipx_frame(desc,TIME=frame_time,N=frame_number)
; Here desc is the descriptor returned by ipx_open function,
; frame numbers start from 0.
; If the function receives one argument (descriptor) only it
; returns first frame after 'frame_time', frame number is
; returned as 'frame_number'.
;
function ipx_frame, desc, num, TIME=ftime, N=frame
  if(n_params() lt 1) then begin
    print,'Returns a frame of IPX file, usage:'
    print,'  image=ipx_frame(desc,frame_number,TIME=frame_time)'
    print,'or'
    print,'  image=ipx_frame(desc,TIME=frame_time,N=frame_number)'
    return,0
  end
  tags=TAG_NAMES(desc)
  if((N_TAGS(desc) ne 4)||(tags[0] ne 'FILENAME')||(tags[1] ne 'HEADER')|| $
     (tags[2] ne 'FRAMEPOS')||(tags[3] ne 'FRAMETIME')) then begin
    print,'Invalid argument - file descriptor'
    return,0
  end
  fn=0L
  if(n_params() gt 1) then begin
    fn=num
    if(fn lt 0) then fn=0
  end else begin
    if(arg_present(ftime)) then begin
      for fn=0,desc.header.numFrames-1 do $
        if(desc.frametime[fn] ge ftime) then break
    end
  end
  if(fn ge desc.header.numFrames) then fn=desc.header.numFrames-1
  ftm=0D
  if(desc.header.depth gt 8) then $ 
    image=uintarr(desc.header.width,desc.header.height,/NOZERO) $
  else $
    image=bytarr(desc.header.width,desc.header.height,/NOZERO)
  ret = CALL_EXTERNAL('libjp2idl.so','ipx_frame', $
          desc.filename,desc.framepos[fn],image,ftm)
  if(ret ne 0) then begin
    print,'Cannot read image, error ',ret
    return,0
  end
  if(arg_present(ftime)) then ftime=ftm
  if(arg_present(frame)) then frame=fn
  return,image
end
