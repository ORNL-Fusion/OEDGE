; Function for IPX1 and IPX2 file access.
; Sergei.Shibaev@ccfe.ac.uk
;
; Returns a frame pixel array or 0 in case of any error.
; 0 is returned also if there is no requested reference frame.
; Usage, three variants:
; 1) image = ipx2frame(desc, frame_number, TIME=frame_time, EXP=exposure)
;   The frame is accessed by number; function returns frame (image), frame time, and frame exposure.
; 2) image = ipx2frame(desc, TIME=frame_time, NUMBER=frame_number, EXP=exposure) 
;   The frame is accessed by frame time (frame_time on input); 
;   function returns frame (image), actual frame time, frame number, and frame exposure.
; 3) image = ipx2frame(desc, frame_number, /REF)
;   Function returns a reference frame with /REF switch.
;   Version 2 of IPX file can contain 3 reference frames:
;     0 - bad pixels;
;     1 - one point NUC correction;
;     2 - two point NUC correction (gain frame).
;   Ref frame 0 is one byte array. Frames 1 and 2 have the same pixel depth as main frames.
;
; Here 'desc' is the descriptor returned by ipx2open function; frame numbers start from 0.
;
function ipx2frame, desc, num, TIME=ftime, NUMBER=frame, EXP=exposure, REF=ref
  if(n_params() lt 2) then begin
    print,'Returns a frame of IPX file, usage:'
    print,'  image = ipx2frame(desc, frame_number, TIME=frame_time, EXP=exposure)'
    print,'or'
    print,'  image = ipx2frame(desc, TIME=frame_time, NUMBER=frame_number, EXP=exposure)'
    print,'or reference frame'
    print,'  image = ipx2frame(desc, ref_number, /REF)'
    return,0
  end
  tags = TAG_NAMES(desc)
  if((N_TAGS(desc) ne 8) or (tags[0] ne 'FILENAME') or (tags[1] ne 'FILEINFO') or $
     (tags[2] ne 'REFPOS') or (tags[3] ne 'REFSIZE') or (tags[4] ne 'FRAMEPOS') or $
     (tags[5] ne 'FRAMESIZE') or (tags[6] ne 'FRAMETIME') or (tags[7] ne 'FRAMEEXP')) then begin
    print, 'Invalid argument - wrong file descriptor'
    return, 0
  end
  fnum = 0UL
  ftm = 0D
  image = 0
  if(n_params() gt 1) then begin
    fnum = num
  end else begin
    if(not arg_present(ftime)) then begin
      print, 'Invalid arguments, no frame number, no frame time'
      return, 0
    end
    for fnum = 0, desc.fileinfo.numFrames - 1 do begin
      if(desc.frametime[fnum] ge ftime) then break
    end
    if(fnum ge desc.fileinfo.numFrames) then fnum = desc.fileinfo.numFrames - 1
  end
  if(keyword_set(ref)) then begin
    if((fnum lt 0) or (fnum gt 2)) then begin
      print, 'Invalid argument - reference frame number ', fnum
      return, 0
    end
    if((desc.refpos[fnum] lt 1) or (desc.refsize[fnum] lt 16)) then begin
      print, 'Reference frame ', fnum, ' not found'
      return, 0
    end
    imgsize = 0UL
    if(fnum eq 0) then begin
      imgsize = desc.fileinfo.width * desc.fileinfo.height
      image = bytarr(desc.fileinfo.width, desc.fileinfo.height, /NOZERO)
    end else begin
      if(desc.fileinfo.depth le 8) then begin
        imgsize = desc.fileinfo.width * desc.fileinfo.height
        image = bytarr(desc.fileinfo.width, desc.fileinfo.height, /NOZERO)
      end else if(desc.fileinfo.depth le 16) then begin
        imgsize = desc.fileinfo.width * desc.fileinfo.height * 2
        image = uintarr(desc.fileinfo.width, desc.fileinfo.height, /NOZERO)
      end else if(desc.fileinfo.depth le 32) then begin
        imgsize = desc.fileinfo.width * desc.fileinfo.height * 4
        image = ulonarr(desc.fileinfo.width, desc.fileinfo.height, /NOZERO)
      end else begin
        print, 'Cannot handle data with depth ', desc.fileinfo.depth
        return, 0
      end
    end
    if(desc.fileinfo.codec gt 0) then begin
      codec = 1 ; lossless compression for ref frames
      ret = CALL_EXTERNAL('libipx2idl.so', 'ipxframe', $
        desc.filename, codec, desc.refpos[fnum], desc.refsize[fnum], image)
      if(ret ne 0) then begin
        print, 'ipxframe error ', ret
        return, 0
      end
    end else begin
      if(desc.refsize[fnum] ne imgsize) then begin
        print, 'Wrong ref frame ', fnum, ' size ', desc.fileinfo.refsize[fnum]
        return, 0
      end
      openr, unit, desc.filename, /GET_LUN, ERROR = err
      if(err ne 0) then begin
        printf, -2, !ERR_STRING
        return, 0
      end
      ; read frame
      point_lun, unit, desc.refpos[fnum]
      readu, unit, image
      free_lun, unit
    end
  end else begin
    if((fnum lt 0) or (fnum ge desc.fileinfo.numFrames)) then begin
      print, 'There is no frame with number ', fnum
      return, 0
    end
    imgsize = 0UL
    if(desc.fileinfo.depth le 8) then begin
      imgsize = desc.fileinfo.width * desc.fileinfo.height
      image = bytarr(desc.fileinfo.width, desc.fileinfo.height, /NOZERO)
    end else if(desc.fileinfo.depth le 16) then begin
      imgsize = desc.fileinfo.width * desc.fileinfo.height * 2
      image = uintarr(desc.fileinfo.width, desc.fileinfo.height, /NOZERO)
    end else if(desc.fileinfo.depth le 32) then begin
      imgsize = desc.fileinfo.width * desc.fileinfo.height * 4
      image = ulonarr(desc.fileinfo.width, desc.fileinfo.height, /NOZERO)
    end else begin
      print, 'Cannot handle data with depth ', desc.fileinfo.depth
      return, 0
    end
    if(desc.fileinfo.codec gt 0) then begin
; slmod begin
      ret = CALL_EXTERNAL('./libipx2idl.so', 'ipxframe', $
        desc.filename, desc.fileinfo.codec, desc.framepos[fnum], desc.framesize[fnum], image)
;
;      ret = CALL_EXTERNAL('/home/sshibaev/IDL/ipx/libipx2idl.so', 'ipxframe', $
;        desc.filename, desc.fileinfo.codec, desc.framepos[fnum], desc.framesize[fnum], image)
; slmod end
      if(ret ne 0) then begin
        print, 'ipxframe error ', ret
        return, 0
      end
    end else begin
      if(desc.framesize[fnum] ne imgsize) then begin
        print, 'Wrong frame ', fnum, ' size ', desc.framesize[fnum]
        return, 0
      end
      openr, unit, desc.filename, /GET_LUN, ERROR = err
      if(err ne 0) then begin
        printf, -2, !ERR_STRING
        return, 0
      end
      ; read frame
      point_lun, unit, desc.framepos[fnum]
      readu, unit, image
      free_lun, unit
    end
  end
  if(not arg_present(ref)) then begin
    if(arg_present(frame)) then frame = fnum
    if(arg_present(ftime)) then ftime = desc.frametime[fnum]
    if(arg_present(exposure)) then begin
      exposure = desc.fileinfo.exposure
      if(exposure eq 0) then begin
        exposure = desc.frameexp[fnum]
      end
    end
  end
  return, image
end
