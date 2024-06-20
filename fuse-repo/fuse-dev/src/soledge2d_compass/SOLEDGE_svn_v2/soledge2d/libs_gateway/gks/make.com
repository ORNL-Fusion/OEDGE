$ if p1 .nes. ""
$ then
$   opt = f$edit(p1,"UPCASE")
$   if opt .eqs. "INSTALL"
$   then
$     copy/protection=(s:rwed,o:rwed,g:rwed,w:re) libgks.olb, gksfont.dat -
	sys$library
$   else
$     write sys$output "Usage: @make [install]"
$     stop
$   endif
$ endif
$!
$ define/nolog sys sys$library
$!
$ if f$getsyi("arch_name") .eqs. "Alpha" then cc :== cc/standard=vaxc
$!
$ cc compress.c
$ cc gksio.c
$ cc gkscbnd.c
$ fortran gks.for
$ fortran gksinq.for
$ fortran gkserror.for
$ fortran gksroot.for
$ fortran gksmisc.for
$ fortran gksdidd.for
$ fortran gksdps.for
$ fortran gksdtek.for
$ fortran gksdtek2.for
$ cc gksdx11.c
$ cc gksduis.c
$ cc gksdcgm.c
$ cc gksdgksm.c
$ cc gkswiss.c
$ fortran gksdwiss.for
$ fortran gksdhpgl.for
$ fortran gksdvt.for
$ fortran gksdpbm.for
$ fortran gksforio.for
$ fortran gksafm.for
$!
$ library/create libgks.olb -
 compress.obj, gksio.obj, gkscbnd.obj, gks.obj, gksinq.obj, gkserror.obj, -
 gksroot.obj, gksmisc.obj, gksdidd.obj, gksdps.obj, gksdtek.obj, gksdtek2.obj,-
 gksdx11.obj, gksduis.obj, gksdcgm.obj, gkswiss.obj, gksdgksm.obj,-
 gksdwiss.obj, gksdhpgl.obj, gksdvt.obj, gksdpbm.obj, gksforio.obj, gksafm.obj
$!
$ if f$search("demo.for") .nes. ""
$ then
$   fortran demo.for
$ else
$   fortran [.gks]demo.for
$ endif
$!
$ if f$getsyi("arch_name") .eqs. "VAX"
$ then
$   link demo, libgks.olb/library, sys$input/options
	sys$share:decw$dwtlibshr/shareable
	sys$share:decw$xlibshr/shareable
	sys$share:uisshr/shareable
	sys$share:vaxcrtl/shareable
$ else
$   link demo, libgks.olb/library, sys$input/options
	sys$share:decw$xtshr/shareable
	sys$share:decw$xlibshr/shareable
$ endif
$!
$ exit
