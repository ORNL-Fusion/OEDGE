/*
*   Some def's (which may appear in Xlib.h)
*   The defines of the FORTRAN callable code is to
*   deal with those machines which link to lower
*   case c routines only e.g sun.
*
*/
#ifndef Bool
#define Bool int
#endif
#ifndef False
#define False 0
#endif
#ifndef True
#define True 1
#endif

#ifndef SIGIO
#define SIGIO SIGPOLL
#endif

#ifdef apollo
#define SysV
#define link_to_lowercase
#define streams
#endif

#ifdef ardent
#define SysV
#define streams
#endif

#ifdef convex
#define link_to_lowercase
#define clock time
#endif

#ifdef dec
#define link_to_lowercase
#endif

#ifdef hp
#define SysV
#define linklc_nu
#endif

#ifdef ibm
#define linklc_nu
#endif

#ifdef mistral
#define link_to_lowercase
#endif

#ifdef rm
#define SysV
#define link_to_lowercase
#define streams
#endif

#ifdef sg
#define SysV
#define link_to_lowercase
#define streams
#endif

#ifdef sun
#define link_to_lowercase
#define streams
#endif

#ifdef link_to_lowercase
#define XINIT  xinit_
#define XPOLYL xpolyl_
#define XPOLYM xpolym_
#define XTEXT  xtext_
#define XCHARH xcharh_
#define XCHARO xcharo_
#define XCHARB xcharb_
#define XLINET xlinet_
#define XFILLP xfillp_
#define XCOLTB xcoltb_
#define XSETFC xsetfc_
#define XGETIN xgetin_
#define XCLRWN xclrwn_
#define XWNDIM xwndim_
#define XCLOSE xclose_
#define XFLSBF xflsbf_
#endif

#ifdef linklc_nu
#define XINIT  xinit
#define XPOLYL xpolyl
#define XPOLYM xpolym
#define XTEXT  xtext
#define XCHARH xcharh
#define XCHARO xcharo
#define XCHARB xcharb
#define XLINET xlinet
#define XFILLP xfillp
#define XCOLTB xcoltb
#define XSETFC xsetfc
#define XGETIN xgetin
#define XCLRWN xclrwn
#define XWNDIM xwndim
#define XCLOSE xclose
#define XFLSBF xflsbf
#endif


#define EVENTSWANTED ExposureMask | KeyPressMask | ButtonPressMask |StructureNotifyMask

#define MAXPOINTS 1000
#define MAXCOLS 256

#ifndef Max
#define Max(x, y)       (((x) > (y)) ? (x) : (y))
#endif
#ifndef Min
#define Min(x, y)       (((x) < (y)) ? (x) : (y))
#endif
