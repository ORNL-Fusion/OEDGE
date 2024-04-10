;-------------------------------------------------------------------------
; Procedure: PSOPEN (Version 1.00)			R.Martin Jan.2008
;
; A standardised routine for creating PS or EPS output (see also
; PSCLOSE)
;
; Opening output file:
;
;   PSOPEN, filename=filename {, /portrait} {, /eps}
;
;   filename: A string containing the name of the PS file to be written,
;	if no filename is given then then the default output filename is 
;	used [sys$scratch:idl.ps] this can be changed using the DEFFILE 
;	option.
;
;   /portarit: Optional parameter, plots by default are landscape, using
;	this option portarte plots can be produced.
;
;   psopen, file='sine.ps'
;   x=findgen(2000)*0.1
;   plot, x, sin(x)
;   psclose, 'd3_209_photocopier'
;
;
;

pro psopen, filename=filename, $
            eps=eps,	       $ 
	    portrait=portrait

  common loc_psoutput, info

  if undefined(info) then info={open:0B, filename:'psoutput.ps', eps:0B}

  if (info.open) then begin
    print, 'PSOPEN: Warning PS output file already open ['+info.filename+']'
    return
  endif

  info.eps=keyword_set(eps)
  if is_string(filename) then info.filename=filename[0] else info.filename=''
  if (info.filename eq '') then info.filename=info.eps ? 'output.eps' : 'output.ps'

  openw, outfile, info.filename, error=error, /get_lun
  if (error ne 0) then begin
    print, 'PSOPEN: Write access denied for file ['+info.filename+']'
    return
  endif
  free_lun, outfile    

  set_plot, 'ps'
  device, filename=expand_path(info.filename), 		$
 	 /color, bits=8, /helvet, encaps=info.eps 

  if keyword_set(eps) then paper=[29.7, 21.0] else paper=[27.7, 19.0]

  if keyword_set(portrait) then begin
    device, /portrait, yoffset=1.0, ysize=paper[0], xoffset=1.0, xsize=paper[1]
  endif else begin
    device, /landscape, yoffset=28.7, ysize=paper[1], xoffset=1.0, xsize=paper[0] 
  endelse

  info.open=1B

end

;-------------------------------------------------------------------------
; Changes
;
; 07.01.08 Version 1.00
; - Split HARDCOPY code into two routines PSOPEN and PSCLOSE
;
