; Added double precision. -15/12/2010, SL
;
; $Id: //depot/idl/IDL_70/idldir/lib/pnt_line.pro#1 $
;
; Copyright (c) 1993-2007, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
Function Pnt_Line_D, p0, l0, l1, pl, INTERVAL = interval
;+
; NAME:
;	PNT_LINE
;
; PURPOSE:
;	This function returns the perpendicular distance between the
;	point P0 and the line between points L0 and L1.
;
; CATEGORY:
;	Geometry.
;
; CALLING SEQUENCE:
;	Result = PNT_LINE(P0, L0, L1 [, Pl])
;
; INPUTS:
;	P0: The location of the point. P0 may have 2 to N elements,
;    	    for N dimensions.
;	L0: One end-point of the line. L0 must have same number of
;	    elements as P0.
;	L1: The other end-point of the line. L1 must have the same
;	    number of elements as LO.
;
; KEYWORD PARAMETERS:
;	INTERVAL: If set, and if the point on the line between L0
;		  and L1 that is closest to PO is not within the
;		  interval [L0, L1], causes the function to return
;		  the distance from P0 to the closer of the two
;		  endpoints L0 and L1.
;
; OUTPUTS:
;	This function returns the distance from point P0 to the line
;	between L0 and L1, unless the closest point on the line is
;	not in the interval [L0, L1] and the keyword INTERVAL is set.
;	In this case, the function returns the distance between P0
;	and the closer of the two end-points.
;
; OPTIONAL OUTPUTS:
;	Pl: The point on the line between L0 and L1 that is closest to P0.
;	    Pl is not necessarily in the interval [L0, L1].
;
; RESTRICTIONS:
;	This function is limited by the machine accuracy of single
;	precision floating point.
;
; PROCEDURE:
;	Solve equations of perpendicular, etc.
;
; EXAMPLE:
;	To print the distance between the point (2,3) and the line
;	from (-3,3) to (5,12), and also the location of the point on
;	the line closest to (2,3), enter the following command:
;
;	  PRINT, PNT_LINE([2,3], [-3,3], [5,12], Pl), Pl
;
; MODIFICATION HISTORY:
; 	DMS, RSI, Jan, 1993.	Written.
;-

; slmod begin
lv = double(l1 - l0)
;
;lv = float(l1 - l0)
; slmod end
l = sqrt(total(lv*lv))
if l eq 0 then begin	;Line is a point
	pl = p0
	return, sqrt(total((l0-float(p0))^2))
	endif
lv = lv / l			;Normal to line
t = - (-total(lv * p0) + total(l0 * lv))/(total(lv * lv))
pl = t * lv + l0		;Closest point on line
out = (t lt 0) or (t gt l)	;Within interval?

if keyword_set(interval) and out then begin  ;Outside interval?
    d1 = sqrt(total((p0-l1)^2))
    d2 = sqrt(total((p0-l0)^2))
    return, d1 < d2
endif else return, sqrt(total((p0 - pl)^2))

end

