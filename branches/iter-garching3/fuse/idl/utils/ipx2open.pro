; Function for IPX file access, both versions 1 and 2
; Sergei.Shibaev@ccfe.ac.uk
;
; Returns a descriptor for IPX file or 0 in case of any error
; Usage:
;   desc = ipx2open(IPX_file_name, RATE=rate)
; RATE is the maximum frame rate in Hz
; The descriptor is used to get a frame:
;   image = ipx2frame(desc, frame_number)
; Full file info can be obtained from the descriptor returned:
;     help, desc.fileinfo, /st
;
function ipx2open, filename, RATE=rate
  if(n_params() lt 1) then begin
    print, 'Opens IPX file, usage: desc = ipx2open(filename, RATE=rate)'
    print, 'rate is the maximum frame rate in the file'
    return, 0
  end
  filename = expand_path(filename)
  openr, unit, filename, /GET_LUN, ERROR = err
  if(err ne 0) then begin
    printf, -2, !ERR_STRING
    return, 0
  end
  ; read file header
  fileid = string(' ', format='(A8)')
  readu, unit, fileid
  if(strpos(fileid, 'IPX ') ne 0) then begin
    free_lun, unit
    printf, -2, filename, ' is not an IPX file'
    return, 0
  end
  ver = 0
  reads, strmid(fileid, 4, 2), ver, format='(I2)'
  if((ver lt 1) or (ver gt 2)) then begin
    free_lun, unit
    printf, -2, filename, ' has unknown version ', ver
    return, 0
  end
  fileinfo = {fileinfo, $
    version:     0,   $
    codec:       0,   $
    datetime:    string(' ', format='(A20)'), $
    shot:        0L,  $
    lens:        string(' ', format='(A24)'), $
    filter:      string(' ', format='(A24)'), $
    view:        string(' ', format='(A80)'), $
    camera:      string(' ', format='(A80)'), $
    numFrames:   0UL, $
    width:       0L,  $
    height:      0L,  $
    depth:       0,   $
    taps:        0,   $
    color:       0,   $
    left:        0L,  $
    top:         0L,  $
    hBin:        0,   $
    vBin:        0,   $
    offset:      fltarr(2), $
    gain:        fltarr(2), $
    preExp:      0E,  $
    exposure:    0E,  $             ; =0 for variable (per frame) exposure
    board_temp:  0E,  $
    ccd_temp:    0E   $
  }
  fpos = 0UL
  if(ver eq 1) then begin
    header = {header, $
      size:        0UL, $
      codec:       string(' ', format='(A8)'), $
      date_time:   string(' ', format='(A20)'), $
      shot:        0L,  $
      trigger:     0E,  $
      lens:        string(' ', format='(A24)'), $
      filter:      string(' ', format='(A24)'), $
      view:        string(' ', format='(A64)'), $
      numFrames:   0UL, $
      camera:      string(' ', format='(A64)'), $
      width:       0U,  $
      height:      0U,  $
      depth:       0U,  $
      orient:      0UL, $             ; not used
      taps:        0U,  $
      color:       0U,  $
      hBin:        0U,  $
      left:        0U,  $
      right:       0U,  $             ; not used
      vBin:        0U,  $
      top:         0U,  $
      bottom:      0U,  $             ; not used
      offset:      intarr(2), $
      gain:        fltarr(2), $
      preExp:      0UL, $
      exposure:    0UL, $
      strobe:      0UL, $             ; not used
      board_temp:  0E,  $
      ccd_temp:    0E   $
    }
    readu, unit, header
    fileinfo.version = 1
    if(strcmp(header.codec, 'JP2', 3)) then fileinfo.codec = 1 $
    else if(strcmp(header.codec, 'JPC', 3)) then fileinfo.codec = 2 $
    else fileinfo.codec = -1
    fileinfo.datetime = strmid(header.date_time, 6, 4) + "-" + strmid(header.date_time, 3, 2) + $
      "-" + strmid(header.date_time, 0, 2) + "T" + strmid(header.date_time, 11)
    fileinfo.shot = header.shot
    fileinfo.lens = header.lens
    fileinfo.filter = header.filter
    fileinfo.view = header.view
    fileinfo.camera = header.camera
    fileinfo.numFrames = header.numFrames
    fileinfo.width = long(header.width)
    fileinfo.height = long(header.height)
    fileinfo.depth = fix(header.depth)
    fileinfo.taps = fix(header.taps)
    fileinfo.color = fix(header.color)
    fileinfo.left = long(header.left)
    fileinfo.top = long(header.top)
    fileinfo.hBin = fix(header.hBin)
    fileinfo.vBin = fix(header.vBin)
    fileinfo.offset[0] = header.offset[0]
    fileinfo.offset[1] = header.offset[1]
    fileinfo.gain[0] = header.gain[0]
    fileinfo.gain[1] = header.gain[1]
    fileinfo.preExp = header.preExp
    fileinfo.exposure = header.exposure
    fileinfo.board_temp = header.board_temp
    if(header.size ne 286) then header.ccd_temp = 0
    fileinfo.ccd_temp = header.ccd_temp
    fpos = header.size
  end else begin
    ; version 2
    fileinfo.version = 2
    sizestr = string(' ', format='(A4)')
    readu, unit, sizestr
    hsize = 0L
    reads, sizestr, hsize, format='(z4)'
    if((hsize le 4) or (hsize gt 512)) then begin
      free_lun, unit
      printf, -2, 'Invalid file header size ', hsize
      return, 0
    end
    bhead = bytarr(hsize - 12)
    readu, unit, bhead
    header = string(bhead)
    fields = strsplit(header, '&', /EXTRACT)
    numf = N_ELEMENTS(fields) - 1
    for i = 0L, numf do begin
      if(strcmp(fields[i], 'codec=', 6)) then begin
        codstr = strlowcase(strmid(fields[i], 6))
        if(strcmp(codstr, 'jp2', 3)) then fileinfo.codec = 1 $
        else if(strcmp(codstr, 'jpc', 3)) then fileinfo.codec = 2 $
        else fileinfo.codec = -1
      end else if(strcmp(fields[i], 'datetime=', 9)) then begin
        fileinfo.datetime = strmid(fields[i], 9, 19)
      end else if(strcmp(fields[i], 'shot=', 5)) then begin
        fileinfo.shot = long(strmid(fields[i], 5))
      end else if(strcmp(fields[i], 'lens=', 5)) then begin
        fileinfo.lens = strmid(fields[i], 5, 23)
      end else if(strcmp(fields[i], 'filter=', 7)) then begin
        fileinfo.filter = strmid(fields[i], 7, 23)
      end else if(strcmp(fields[i], 'view=', 5)) then begin
        fileinfo.view = strmid(fields[i], 5, 63)
      end else if(strcmp(fields[i], 'camera=', 7)) then begin
        fileinfo.camera = strmid(fields[i], 7, 63)
      end else if(strcmp(fields[i], 'frames=', 7)) then begin
        fileinfo.numFrames = ulong(strmid(fields[i], 7))
      end else if(strcmp(fields[i], 'width=', 6)) then begin
        fileinfo.width = long(strmid(fields[i], 6))
      end else if(strcmp(fields[i], 'height=', 7)) then begin
        fileinfo.height = long(strmid(fields[i], 7))
      end else if(strcmp(fields[i], 'depth=', 6)) then begin
        fileinfo.depth = fix(strmid(fields[i], 6))
      end else if(strcmp(fields[i], 'taps=', 5)) then begin
        fileinfo.taps = fix(strmid(fields[i], 5))
      end else if(strcmp(fields[i], 'color=', 6)) then begin
        color = strlowcase(strmid(fields[i], 6))
        if(strcmp(color, 'gbrg/rggb')) then fileinfo.color = 1 $
        else if(strcmp(color, 'gr/bg')) then fileinfo.color = 2
      end else if(strcmp(fields[i], 'left=', 5)) then begin
        fileinfo.left = long(strmid(fields[i], 5))
      end else if(strcmp(fields[i], 'top=', 4)) then begin
        fileinfo.top = long(strmid(fields[i], 4))
      end else if(strcmp(fields[i], 'hbin=', 5)) then begin
        fileinfo.hBin = fix(strmid(fields[i], 5))
      end else if(strcmp(fields[i], 'vbin=', 5)) then begin
        fileinfo.vBin = fix(strmid(fields[i], 5))
      end else if(strcmp(fields[i], 'offset=', 7)) then begin
        ostr = strmid(fields[i], 7)
        if((byte(ostr[0]) eq byte('"')) or (byte(ostr[0]) eq byte("'"))) then begin
          ostr = strmid(ostr, 1, strlen(ostr) - 2)
          f0 = 0E
          f1 = 0E
          reads, ostr, f0, f1
          fileinfo.offset[0] = f0
          fileinfo.offset[1] = f1
        end else begin
          fileinfo.offset[0] = float(ostr)
        end
      end else if(strcmp(fields[i], 'gain=', 5)) then begin
        gstr = strmid(fields[i], 5)
        if((byte(gstr[0]) eq byte('"')) or (byte(gstr[0]) eq byte("'"))) then begin
          gstr = strmid(gstr, 1, strlen(gstr) - 2)
          g0 = 0E
          g1 = 0E
          reads, gstr, g0, g1
          fileinfo.gain[0] = g0
          fileinfo.gain[1] = g1
        end else begin
          fileinfo.gain[0] = float(gstr)
        end
      end else if(strcmp(fields[i], 'preexp=', 7)) then begin
        fileinfo.preExp = ulong(strmid(fields[i], 7))
      end else if(strcmp(fields[i], 'exposure=', 9)) then begin
        fileinfo.exposure = ulong(strmid(fields[i], 9))
      end else if(strcmp(fields[i], 'boardtemp=', 10)) then begin
        fileinfo.board_temp = float(strmid(fields[i], 10))
      end else if(strcmp(fields[i], 'ccdtemp=', 8)) then begin
        fileinfo.ccd_temp = float(strmid(fields[i], 8))
      end
    end
    fpos = hsize
  end
  if(fileinfo.numFrames lt 1) then begin
    free_lun, unit
    printf, -2, 'No frames in ', filename
    return, 0
  end
  if((fileinfo.width lt 1) or (fileinfo.height lt 1) or (fileinfo.depth lt 1) or (fileinfo.codec < 0) or (fileinfo.codec > 2)) then begin
    free_lun, unit
    printf, -2, 'File ', filename, ' is wrong: width=', fileinfo.width, $
      ', height=', fileinfo.height, ', depth=', fileinfo.depth, codec=', fileinfo.codec
    return, 0
  end
  refpos  = ulonarr(3)
  refsize = ulonarr(3)
  framepos  = ulonarr(fileinfo.numFrames)
  framesize = ulonarr(fileinfo.numFrames)
  frametime = dblarr(fileinfo.numFrames)
  frameexp = 0
  if(fileinfo.exposure eq 0) then begin
    frameexp = fltarr(fileinfo.numFrames)
  end
  inter = 1e10
  ; scan whole file
  if(fileinfo.version eq 1) then begin
    nf = 0UL
    while(nf lt fileinfo.numFrames) do begin
      point_lun, unit, fpos
      fhead = {fhead, fsize: 0UL, ftime: 0D}
      readu, unit, fhead
      frametime[nf] = fhead.ftime
      if(nf gt 0) then begin
        dt = frametime[nf] - frametime[nf - 1]
        if(inter gt dt) then inter = dt
      end
      framepos[nf] = fpos + 12
      framesize[nf] = fhead.fsize - 12
      fpos = fpos + fhead.fsize
      nf = nf + 1
    end
  end else begin
    framesz = fileinfo.width * fileinfo.height * ((fileinfo.depth + 7) / 8)
    nf = 0UL
    while(nf lt fileinfo.numFrames) do begin
      point_lun, unit, fpos
      sizestr = string(' ', format='(A2)')
      readu, unit, sizestr
      hsize = 0
      reads, sizestr, hsize, format='(z2)'
      if((hsize le 2) or (hsize ge 128)) then begin
        free_lun, unit
        printf, -2, 'Invalid frame ', nf, ' header: size = ', hsize
        return, 0
      end
      bhead = bytarr(hsize - 2)
      readu, unit, bhead
      shead = string(bhead)
      pos = strpos(shead, 'ftime=')
      if(pos gt 0) then begin            
        ; image frame
        frametime[nf] = double(strmid(shead, pos + 6))
        if(nf gt 0) then begin
          dt = frametime[nf] - frametime[nf - 1]
          if(inter gt dt) then inter = dt
        end
        pos = strpos(shead, 'fsize=')
        if(pos gt 0) then begin
          framesize[nf] = ulong(strmid(shead, pos + 6))
        end else begin
          if(fileinfo.codec gt 0) then begin
            free_lun, unit
            printf, -2, 'Invalid frame ', nf, ' header: no size'
            return, 0
          end
          framesize[nf] = framesz
        end
        if(fileinfo.exposure eq 0) then begin
          pos = strpos(shead, 'fexp=')
          if(pos gt 0) then begin
            frameexp[nf] = float(strmid(shead, pos + 5))
          end
        end
        framepos[nf] = fpos + hsize
        fpos = framepos[nf] + framesize[nf]
        nf = nf + 1
      end else begin
        ; ref frame
        pos = strpos(shead, 'ref=')
        if(pos lt 0) then begin
          free_lun, unit
          printf, -2, 'Invalid frame ', nf, ' header: no time, no ref'
          return, 0
        end
        ref = fix(strmid(shead, pos + 4))
        if((ref lt 0) or (ref gt 2)) then begin
          free_lun, unit
          printf, -2, 'Invalid ref frame number (', ref, ')'
          return, 0
        end
        pos = strpos(shead, 'fsize=')
        if(pos gt 0) then begin
          refsize[ref] = ulong(strmid(shead, pos + 6))
        end else begin
          if(fileinfo.codec gt 0) then begin
            free_lun, unit
            printf, -2, 'Invalid ref ', ref, ' frame header: no size'
            return, 0
          end
          if(ref eq 0) then begin
            ; ref 0 pixel is always 1 byte
            refsize[ref] = fileinfo.width * fileinfo.height
          end else begin
            refsize[ref] = framesz
          end
        end
        refpos[ref] = fpos + hsize
        fpos = refpos[ref] + refsize[ref]
      end
    end
  end
  free_lun,unit
  if(arg_present(rate)) then begin
    rate = 0L
    if(fileinfo.numFrames gt 1) then begin
      if(inter gt 0) then rate = round(1.0 / inter)
    end
  end
  desc = {               $   
    filename :filename,  $
    fileinfo :fileinfo,  $
    refpos   :refpos,    $
    refsize  :refsize,   $
    framepos :framepos,  $
    framesize:framesize, $
    frametime:frametime, $
    frameexp :frameexp   $
  }
  return, desc
end
