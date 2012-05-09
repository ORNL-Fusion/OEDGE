;-------------------------------------------------------------------------
; Procedure: PSCLOSE (Version 1.00)			R.Martin Jan.2008
;
; Procedure to be used with PSOPEN, closes the PS device and
; optionally prints a PS file.
;
; Calling Sequence:
;
;   PSCLOSE {, queue=queue}
;
;   queue: Optional parameter, if defined PSCLOSE closes the PS file
;   and sends the file to the specified print queue.
;
; Example:
;
;   psopen, file='sine.ps'
;   x=findgen(2000)*0.1
;   plot, x, sin(x)
;   psclose, 'd3_209_photocopier'
;
;

pro psclose, queue=queue

  common loc_psoutput, info

  if undefined(info) then info={open:0B, filename:'psoutput.ps', eps:0B}

  if (info.open eq 0) then return

  set_plot, 'ps'
  device, /close
  device, encaps=0
  set_plot, 'x'

  print, 'Plot stored in file: [', info.filename, ']'
  if is_string(queue) and not(info.eps) then begin
    print, 'File sent to ['+queue[0]+']'
    dclc='lpr -P'+queue[0]+' '+info.filename
    spawn, dclc
  endif

  info={open:0B, filename:'psoutput.ps', eps:0B}

end

;-------------------------------------------------------------------------
; 07.01.2008
; - New code based on the HARDCOPY routine
;
