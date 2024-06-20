/*
 * Copyright @ 1984 - 1996   Josef Heinen
 *
 * Permission to use, copy, and distribute this software and its
 * documentation for any purpose with or without fee is hereby granted,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 *
 * Permission to modify the software is granted, but not the right to
 * distribute the modified code.  Modifications are to be distributed
 * as patches to released version.
 *
 * This software is provided "as is" without express or implied warranty.
 *
 * Send your comments or suggestions to
 *  J.Heinen@kfa-juelich.de.
 *
 *
 * FACILITY:
 *
 *	GLI GKS V4.5
 *
 * ABSTRACT:
 *
 *	This module contains a logical device driver for XWindow
 *	displays.
 *
 * AUTHOR:
 *
 *	Josef Heinen
 *
 * VERSION:
 *
 *	V1.0
 *
 */

#include <stdio.h>

#include "gksdefs.h"

#if !defined(NO_X11)

#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>

#if !defined(VMS) && !defined(MSDOS) && !defined(_WIN32)
#include <unistd.h>
#endif

#ifdef __osf__
int usleep(useconds_t);
#endif

#include <sys/types.h>

#ifdef XSHM
#include <sys/ipc.h>
#include <sys/shm.h>
#endif

#if defined(cray) || defined(__SVR4) || defined(_WIN32)
#include <fcntl.h>
#else
#include <sys/file.h>
#endif

#include <sys/stat.h>

#if !defined(VMS) && !defined(__sun) && !defined(linux) && !defined(_WIN32)
#include <signal.h>
#include <sys/ioctl.h>
#endif

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>
#include <X11/keysym.h>

#include <X11/Intrinsic.h>

#ifdef XSHM
#include <X11/extensions/XShm.h>
#endif

#endif

#ifdef VMS
#include <descrip.h>
#endif

#ifdef cray
#include <fortran.h>
#endif

#ifdef VMS
#define CHARARG(a) struct dsc$descriptor *a
#else
#ifdef cray
#define CHARARG(a) _fcd a
#else
#if defined(_WIN32) && !defined(__GNUC__)
#define CHARARG(a) char *(a), unsigned short a##_len
#else
#define CHARARG(a) char *a
#endif /* _WIN32 */
#endif /* cray */
#endif /* VMS */

#if defined(_WIN32) && !defined(__GNUC__)
#define STDCALL __stdcall
#else
#define STDCALL
#endif

#if !defined(NO_X11)

#include "icon.bm"

#ifndef min
#define min(a,b)	(((a)<(b)) ? (a) : (b))
#endif
#ifndef max
#define max(a,b)	(((a)>(b)) ? (a) : (b))
#endif
#define nint(a)		((int)(a + 0.5))

#define WindowName "GLIgks V4.5.28"

#define DrawBorder	0
#define Undefined	0xffff

#define Pi		3.141592654

#define MAX_TNR         9
#define MAX_PIXMAP      512
#define MAX_SIZE        100

#define PRIVATE_COLORS	8
#define SHARED_COLORS	72
#define CMAP_EXTENT     84
#define MAX_COLORS	256
#define MAX_POINTS	2048
#define PATTERNS	120
#define HATCH_STYLE     108
#define MAX_COLORIND	1024

#define WHITE 255
#define THRESH 127
#define BLACK 0
#define INITERR(X,Y)    (X - (Y ? WHITE : BLACK) + (THRESH - X)/2)

#define LEFT   (1<<0)
#define RIGHT  (1<<1)
#define BOTTOM (1<<2)
#define TOP    (1<<3)

#define CTRL_C 3
#define CTRL_D 4
#define CTRL_Z 26

#define WC_to_NDC(xw, yw, tnr, xn, yn) \
    xn = a[tnr] * (xw) + b[tnr]; \
    yn = c[tnr] * (yw) + d[tnr]

#define WC_to_NDC_rel(xw, yw, tnr, xn, yn) \
    xn = a[tnr] * (xw); \
    yn = c[tnr] * (yw)

#define NDC_to_WC(xn, yn, tnr, xw, yw) \
    xw = ((xn) - b[tnr]) / a[tnr]; \
    yw = ((yn) - d[tnr]) / c[tnr]

#define NDC_to_DC(xn, yn, xd, yd) \
    xd = sint(p->a * (xn) + p->b + 0.5); \
    yd = sint(p->c * (yn) + p->d + 0.5)

#define DC_to_NDC(xd, yd, xn, yn) \
    xn = ((xd) - p->b) / p->a; \
    yn = ((yd) - p->d) / p->c;

#define CharXform(xrel, yrel, x, y) \
    x = cos_f[p->path] * (xrel) - sin_f[p->path] * (yrel); \
    y = sin_f[p->path] * (xrel) + cos_f[p->path] * (yrel);

#define Color8Bit(c) \
    (c) >= 588 ? 80 + (c-588)/56 * 12 + nint((c-588)%56 * 11.0/56.0) : \
    (c) >= 257 ? 8 + nint((c-257)/330.0 * (72-1)) : (c)

static char *fonts[] = {
    "-adobe-times-medium-r-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-times-medium-i-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-times-bold-r-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-times-bold-i-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-helvetica-medium-r-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-helvetica-medium-o-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-helvetica-bold-r-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-helvetica-bold-o-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-courier-medium-r-normal--*-%d0-%d-%d-m-*-*-*",
    "-adobe-courier-medium-o-normal--*-%d0-%d-%d-m-*-*-*",
    "-adobe-courier-bold-r-normal--*-%d0-%d-%d-m-*-*-*",
    "-adobe-courier-bold-o-normal--*-%d0-%d-%d-m-*-*-*",
    "-adobe-symbol-medium-r-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-itc lubalin graph-book-r-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-itc lubalin graph-book-o-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-itc lubalin graph-demi-r-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-itc lubalin graph-demi-o-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-new century schoolbook-medium-r-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-new century schoolbook-medium-i-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-new century schoolbook-bold-r-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-new century schoolbook-bold-i-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-itc avant garde gothic-book-r-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-itc avant garde gothic-book-o-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-itc avant garde gothic-demi-r-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-itc avant garde gothic-demi-o-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-itc souvenir-light-r-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-itc souvenir-light-i-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-itc souvenir-demi-r-normal--*-%d0-%d-%d-p-*-*-*",
    "-adobe-itc souvenir-demi-i-normal--*-%d0-%d-%d-p-*-*-*"
    };
static int n_font = 29;

static int map[32] = {
    22,  9,  5, 14, 18, 26, 13,  1,
    24, 11,  7, 16, 20, 28, 13,  3,
    23, 10,  6, 15, 19, 27, 13,  2,
    25, 12,  8, 17, 21, 29, 13,  4
    };

static float capheights[29] = {
    0.72, 0.72, 0.72, 0.72,
    0.76, 0.76, 0.76, 0.76,
    0.6, 0.6, 0.6, 0.6,
    0.76,
    0.76, 0.76, 0.76, 0.76,
    0.76, 0.76, 0.76, 0.76,
    0.76, 0.76, 0.76, 0.76,
    0.76, 0.76, 0.76, 0.76
    };

static char patterns[PATTERNS][33];
static Bool have_patterns = False;

static float xfac[4] = {0, 0, -0.5, -1};
static float yfac[6] = {0, -1.2, -1, -0.5, 0, 0.2};

static float sin_f[] = {0, 1, 0, -1};
static float cos_f[] = {1, 0, -1, 0};

static int predef_font[] = {1, 1, 1, -2, -3, -4};
static int predef_prec[] = {0, 1, 2, 2, 2, 2};
static int predef_ints[] = {0, 1, 3, 3, 3};
static int predef_styli[] = {1, 1, 1, 2, 3};

static XPoint *points = NULL;
static int max_points = MAX_POINTS;

#if defined(hpux) && !defined(NAGware)
static void (*line_routine_a)();
static void (*fill_routine_a)();
static void (*text_routine_a)();
#endif

typedef enum {
    TypeNone, TypeLocal, TypeCrosshair, TypeCross, TypeRubberband,
    TypeRectangle, TypeDigital, TypeCircle} pe_type;

typedef unsigned char byte;

typedef enum {
    AnyServer, SunServer, SGIServer, DECServer, IBMServer, TEKServer
    } server_type;

typedef struct ws_state_list_struct {
    int wkid;
    int gif, rf;
    FILE *uil;
    Bool packed_ca;
    Bool async_io;
    Widget widget;
    int conid, wstype;
    Display *dpy;
    Bool new_dpy;
    int fd;
    Bool async;
    Screen *screen;
    Bool backing_store;
    unsigned long fg, bg;
    Visual *vis;
    int depth;
    Colormap cmap;
    Window win;
    Bool new_win;
    Pixmap pixmap, icon_pixmap;
    Bool double_buf;
    int shape;
    XImage *shmimage;
#ifdef XSHM
    XShmSegmentInfo shminfo;
#endif
    GC gc, clear, invert;
    long event_mask;
    Cursor cursor, textcursor;
    int swidth, sheight, dpi, x, y, width, height;
    float mwidth, mheight, magnification, window[4], viewport[4];
    int state, mapped;
    Bool empty;
    int path;
    float ratio;
    int cached_size;
    XFontStruct *fstr[29][MAX_SIZE], *cfont;
    int capheight;
    Pixmap tile[MAX_COLORS][PATTERNS];
    Pixmap stipple[MAX_COLORS][PATTERNS];
    Bool ored_patterns;
    XColor color[MAX_COLORS];
    unsigned long pixels[MAX_COLORS], ccolor, ncolors, pcolors, scolors;
    Bool mono_flag, min_colors, dark_bg;
    float red[MAX_COLORS], green[MAX_COLORS], blue[MAX_COLORS];
    float gray[MAX_COLORS];
    int ltype;
    unsigned int lwidth;
    float a, b, c, d;
    pe_type type;
    int px, py;
    char *error;
    server_type server;
    Bool xshm, gr_bc;
    Pixmap *frame;
    int nframes;
    } ws_state_list;

typedef struct { 
    int ch;
    char seq[3];
    char alt_seq[3];
    } compose_keys;

static compose_keys key_bindings[] = { 
    {  34, "\" ",""},       /* quotation mark */
    {  35, "++",""},        /* number sign */
    {  39, "' ",""},        /* apostrophe */    
    {  64, "AA",""},        /* commercial at */
    {  91, "((",""},        /* opening bracket */
    {  92, "//","/<"},      /* backslash */
    {  93, "))",""},        /* closing bracket */
    {  94, "^ ",""},        /* circumflex accent */
    {  96, "` ",""},        /* grave accent */
    { 123, "(-",""},        /* opining brace */
    { 124, "/^",""},        /* vertical line */
    { 125, ")-",""},        /* closing brace */
    { 126, "~ ",""},        /* tilde */     
    { 160, "  ", "" },      /* no break space*/
    { 161, "!!", "" },      /* inverted ! */
    { 162, "C/", "C|" },    /* cent sign */
    { 163, "L-", "L=" },    /* pound sign */
    { 164, "XO", "X0" },    /* currency sign  */
    { 165, "Y-", "Y=" },    /* yen sign  */
    { 166, "||", "!^" },    /* broken vertical bar */
    { 167, "SO", "S!" },    /* section sign */
    { 168, "\"\"", "" },    /* diaeresis */
    { 169, "CO", "C0" },    /* copyright sign */
    { 170, "A_", "" },      /* feminine ordinal */ 
    { 171, "<<", "" },      /* open angle brackets */ 
    { 172, "-,", "" },      /* logical not */
    { 173, "-^", "" },      /* macron */
    { 174, "RO", "" },      /* registered trademark */
    { 175, "--", "" },      /* soft (syllable) hyphen */
    { 176, "0^", "" },      /* degree sign */ 
    { 177, "+-", "" },      /* plus or minus sign */ 
    { 178, "2^", "" },      /* superscript 2 */
    { 179, "3^", "" },      /* superscript 3 */
    { 180, "''", "" },      /* acute accent */
    { 181, "/U", "" },      /* micro sign */
    { 182, "P!", "" },      /* paragraph sign */
    { 183, ".^", "" },      /* middle dot  */
    { 184, ", ", "" },      /* cedilla */
    { 185, "1^", "" },      /* superscript 1 */
    { 186, "O_", "" },      /* masculine ordinal */    
    { 187, ">>", "" },      /* closed angle brackets */
    { 188, "14", "" },      /* fraction one-quarter */  
    { 189, "12", "" },      /* fraction one-half */
    { 190, "34", "" },      /* three quarters */
    { 191, "??", "" },      /* inverted ? */ 
    { 192, "`A", "" },      /* A grave  */ 
    { 193, "'A", "" },      /* A acute  */
    { 194, "^A", "" },      /* A circumflex */
    { 195, "~A", "" },      /* A tilde */
    { 196, "\"A", "" },     /* A umlaut */  
    { 197, "A*", "" },      /* A ring */ 
    { 198, "AE", "" },      /* A E diphthong */
    { 199, "C,", "" },      /* C cedilla */
    { 200, "`E", "" },      /* E grave */
    { 201, "'E", "" },      /* E acute */
    { 202, "^E", "" },      /* E circumflex */
    { 203, "\"E", "" },     /* E umlaut */
    { 204, "`I", "" },      /* I grave */
    { 205, "'I", "" },      /* I acute */
    { 206, "^I", "" },      /* I circumflex */ 
    { 207, "\"I", "" },     /* I umlaut */
    { 208, "-D", "" },      /* capital Icelandic Eth */
    { 209, "~N", "" },      /* N tilde */
    { 210, "`O", "" },      /* O grave */
    { 211, "'O", "" },      /* O acute */
    { 212, "^O", "" },      /* O circumflex */
    { 213, "~O", "" },      /* O tilde */
    { 214, "\"O", "" },     /* O umlaut */
    { 215, "xx", "" },      /* multiplication sign */
    { 216, "o/", "" },      /* O slash */
    { 217, "`U", "" },      /* U grave */
    { 218, "'U", "" },      /* U acute */
    { 219, "^U", "" },      /* U circumflex */
    { 220, "\"U", "" },     /* U umlaut */
    { 221, "'Y", "" },      /* Y acute */
    { 222, "TH", "" },      /* capital Icelandic thorn */
    { 223, "ss", "" },      /* German small sharp s */
    { 224, "`a", "" },      /* a grave */
    { 225, "'a", "" },      /* a acute */
    { 226, "^a", "" },      /* a circumflex */  
    { 227, "~a", "" },      /* a tilde */
    { 228, "\"a", "" },     /* a umlaut */
    { 229, "a*", "" },      /* a ring */
    { 230, "ae", "" },      /* a e diphthong */
    { 231, "c,", "" },      /* c cedilla  */
    { 232, "`e", "" },      /* e grave */
    { 233, "'e", "" },      /* e acute */
    { 234, "^e", "" },      /* e circumflex */
    { 235, "\"e", "" },     /* e umlaut */
    { 236, "`i", "" },      /* i grave */
    { 237, "'i", "" },      /* i acute */ 
    { 238, "^i", "" },      /* i circumflex */
    { 239, "\"i", "" },     /* i umlaut */
    { 240, "-d", "" },      /* small Icelandic Eth */
    { 241, "~n", "" },      /* n tilde  */ 
    { 242, "`o", "" },      /* o grave */
    { 243, "'o", "" },      /* o acute */
    { 244, "^o", "" },      /* o circumflex */ 
    { 245, "~o", "" },      /* o tilde */
    { 246, "\"o", "" },     /* o umlaut */
    { 247, "-:", "" },      /* division sign */
    { 248, "o/", "" },      /* o slash  */
    { 249, "`u", "" },      /* u grave */ 
    { 250, "'u", "" },      /* u acute */
    { 251, "^u", "" },      /* u circumflex */
    { 252, "\"u", "" },     /* u umlaut */
    { 253, "'y", "" },      /* y acute */
    { 254, "th", "" },      /* small Icelandic thorn */
    { 255, "\"y", "" }      /* y umlaut */
    };
static int n_key = sizeof(key_bindings)/sizeof(key_bindings[0]);

static gks_state_list *gksl;
static float a[MAX_TNR], b[MAX_TNR], c[MAX_TNR], d[MAX_TNR];

static ws_state_list *p;
static int error_code, request_code, function_id;


static
int *handler (Display *dpy, XErrorEvent *event)
{
    char str[80], request[40];

    if (event->error_code != error_code || event->request_code != request_code)
	{
	XGetErrorText (dpy, event->error_code, str, sizeof(str));
	gks_fprintf (stderr, "GKS: X Protocol error detected by server: %s\n",
	    str);

	sprintf (request, "XRequest.%d", event->request_code);
	XGetErrorDatabaseText (dpy, "", request, "unknown", str, sizeof(str));
	gks_fprintf (stderr, "Failed request major op code %d (%s)\n",
	    event->request_code, str);

	gks_fprintf (stderr, "Invoked from within GKS function id %d\n",
	    function_id);

	error_code = event->error_code;
	request_code = event->request_code;
	}

    return (NULL);
} 


static
int sint (double a)
{
    if (a > 65535)
	return 65535;
    else if (a < -65535)
	return -65535;
    else
	return (int)(a + 0.5);
}


static
void seg_xform (float *x, float *y)
{
    float xx;

    xx = *x * gksl->mat[0][0] + *y * gksl->mat[0][1] + gksl->mat[2][0];
    *y = *x * gksl->mat[1][0] + *y * gksl->mat[1][1] + gksl->mat[2][1]; 
    *x = xx;
}


static
void seg_xform_rel (float *x, float *y)
{
    float xx;

    xx = *x * gksl->mat[0][0] + *y * gksl->mat[0][1];
    *y = *x * gksl->mat[1][0] + *y * gksl->mat[1][1]; 
    *x = xx;
}


static
void set_clipping (Bool state)
{
    float clrt[4];
    int i, j;
    XRectangle rt;

    if (state && gksl->clip == GCLIP)
        {
        memcpy(clrt, gksl->viewport[gksl->cntnr], 4*sizeof(float));
	if (p->gr_bc == 0)
	    {
	    seg_xform (&clrt[0], &clrt[2]);
	    seg_xform (&clrt[1], &clrt[3]);
	    }
	i = clrt[0] < clrt[1] ? 0 : 1;
	j = clrt[2] < clrt[3] ? 2 : 3;

        rt.x = (int) (p->a * clrt[i]   + p->b);
        rt.y = (int) (p->c * clrt[5-j] + p->d);
        rt.width  = (int) (p->a * (clrt[1-i]-clrt[i]  )) + 2;
        rt.height = (int) (p->c * (clrt[j]  -clrt[5-j])) + 2;

        XSetClipRectangles (p->dpy, p->gc, 0, 0, &rt, 1, Unsorted);
        }
    else
        XSetClipMask (p->dpy, p->gc, None);

    rt.x = 0;
    rt.y = 0;
    rt.width = p->width;
    rt.height = p->height;

    XSetClipRectangles (p->dpy, p->invert, 0, 0, &rt, 1, Unsorted);
}


static
void expose_event (Widget widget, ws_state_list *p, XExposeEvent *event,
    Boolean *continue_to_dispatch)

/*
 *  Handle expose events
 */

{
    if (p->pixmap) {
        set_clipping (False);
	XCopyArea (p->dpy, p->pixmap, p->win, p->gc, event->x, event->y,
	    event->width, event->height, event->x, event->y);
        set_clipping (True);
        }
}


static
Display *open_display (void)

/*
 *  Open display
 */

{
    char *env, *ep;
    char s[80];

    env = (char *) getenv ("GLI_CONID");
    if (!env) env = (char *) getenv ("GKSconid");

    if (p->wstype == 213)
	{
	if (env == NULL)
	    {
	    gks_fprintf (stderr, "GKS: can't obtain widget id\n");
	    return (NULL);
	    }
	else
	    sscanf (env, "%ld", &p->widget); /* may be a 64-Bit pointer */
	}

    if (p->widget == NULL)
	{
	if (p->wstype == 212)
	    {
	    if (env == NULL)
		{
		gks_fprintf (stderr,
		    "GKS: can't obtain pre-existing drawable\n");
		return (NULL);
		}
	    else {
		if (sscanf (env, "%d!%d", &p->dpy, &p->win) != 2)
		    {
		    ep = strchr (env, '!');
		    if (ep != NULL) {
			if (strncmp (++ep, "0x", 2) == 0)
			    sscanf (ep + 2, "%x", &p->win);
			else
			    sscanf (ep, "%d", &p->win);
			}
#ifdef _WIN32
		    if (*env == ':')
			sprintf (s, "localhost%s", env);
		    else
			strcpy (s, env);
#else
		    strcpy (s, env);
#endif
		    strtok (s, "!");
		    p->dpy = XOpenDisplay (s);
		    p->new_dpy = True;
		    }
		}
	    }
	else {
	    if (!env) env = (char *) getenv ("DISPLAY");
	    if (env != NULL)
		{
#ifdef _WIN32
		if (*env == ':')
		    sprintf (s, "localhost%s", env);
		else
		    strcpy (s, env);
#else
		strcpy (s, env);
#endif
		env = s;
		}
	    p->dpy = XOpenDisplay (env);
	    p->new_dpy = True;
	    }

	if (p->dpy == NULL)
            {
	    if (!env) env = "";
	    gks_fprintf (stderr, "GKS: can't open display on \"%s\"\n", env);
	    return (NULL);
	    }

        p->screen = XDefaultScreenOfDisplay(p->dpy);
	}
    else {
	p->dpy = XtDisplay (p->widget);
	p->new_dpy = False;
        p->screen = XtScreenOfObject (p->widget);
        }

    p->fd = ConnectionNumber(p->dpy);
    p->async = False;

    XSetErrorHandler ((XErrorHandler) handler);
    error_code = request_code = 0;

    p->backing_store = (XDoesBackingStore(p->screen) == Always) ||
        ((char *) getenv ("GLI_GKS_BS") != NULL);

    p->mwidth  = XWidthMMOfScreen(p->screen) * 0.001;
    p->mheight = XHeightMMOfScreen(p->screen) * 0.001;
    p->swidth = XWidthOfScreen(p->screen);
    p->sheight = XHeightOfScreen(p->screen);

    p->magnification = 1;

    if ((env = (char *) getenv ("GLI_GKS_DPI")) != NULL)
        p->dpi = atoi(env);
    else
	p->dpi = 75;

    p->ored_patterns = (char *) getenv("GLI_GKS_TRANSPARENT_PATTERNS") != NULL;

    return (p->dpy);
}


static
void set_colors (void)
{
    int i, coli;

    p->dark_bg = (char *) getenv("GLI_GKS_DARK_BG") != NULL;

    for (i=0; i<MAX_COLORS; i++) {
	coli = (i < 2 && p->dark_bg) ? 1 - i : i;
	GQRGB (&coli, &p->red[i], &p->green[i], &p->blue[i]);
	p->gray[i] = 0.3*p->red[i] + 0.59*p->green[i] + 0.11*p->blue[i];
 	}
}


static
void configure_colors (void)
{
    char *env;

    if ((env = (char *) getenv("GLI_GKS_CMAP_EXTENT")) != NULL)
        {
        if (*env)
            p->scolors += atoi(env);
        else
            p->scolors += CMAP_EXTENT;

        if (PRIVATE_COLORS + p->scolors > MAX_COLORS)
            p->scolors = MAX_COLORS;
        }
}


static
Bool allocate_color_cells (int ncolors)
{
    static unsigned long plane_masks[] = {0};
    Bool contig;

    contig = True;
    if (XAllocColorCells (p->dpy, p->cmap, contig, plane_masks, 0, p->pixels,
	ncolors))
	{
	p->pcolors = p->ncolors = ncolors;
	return True;
	}
    else
	return False;
}


static
void allocate_private_colors (void)
{
    int i;
    Bool have_colors = False;

    have_colors = allocate_color_cells (PRIVATE_COLORS + p->scolors);
    if (!have_colors)
	have_colors = allocate_color_cells (PRIVATE_COLORS);
    if (!have_colors) {
	gks_fprintf (stderr, "GKS: unable to allocate private colors\n");
	return;
	}

    for (i=0; i<p->ncolors; i++)
	{
	p->color[i].pixel = p->pixels[i];
	p->color[i].flags = DoRed | DoGreen | DoBlue;

	p->color[i].red   = (unsigned short)(p->red[i]*65535);
	p->color[i].green = (unsigned short)(p->green[i]*65535);
	p->color[i].blue  = (unsigned short)(p->blue[i]*65535);
        }

    XStoreColors (p->dpy, p->cmap, p->color, p->ncolors);
}


static
void pick_closest_colors (XColor *color, int ncolors)
{
    XColor ctab[256];
    int i, j, maxcol, closest, d, mdist;

    maxcol = XCellsOfScreen (p->screen);
    maxcol = (maxcol < 256) ? maxcol : 256;

    p->error = "GKS: failed to allocate at least one pixel; use closest color";

    for (i=0; i<maxcol; i++)
	ctab[i].pixel = i;

    XQueryColors (p->dpy, p->cmap, ctab, maxcol);

    for (i=0; i<ncolors; i++)
	{
	if (color[i].pixel == 0xffff)
	    {
	    mdist = 65536*3;
	    closest = 0;

	    for (j=0; j<maxcol; j++) {
		d = abs(color[i].red - ctab[j].red) +
		    abs(color[i].green - ctab[j].green) +
		    abs(color[i].blue - ctab[j].blue);
		if (d < mdist) {
		    mdist = d;
		    closest = j;
		    }
		}

	    color[i].pixel = ctab[closest].pixel;
	    }
	}
}


static
void allocate_pseudo_colors (void)
{
    int i, k;

    for (i=k=0; i<PRIVATE_COLORS; i++)
	{
	p->color[i].red   = (unsigned short)(p->red[i]*65535);
	p->color[i].green = (unsigned short)(p->green[i]*65535);
	p->color[i].blue  = (unsigned short)(p->blue[i]*65535);

	if (!XAllocColor (p->dpy, p->cmap, &p->color[i]))
	    {
	    p->color[i].pixel = 0xffff;
	    k++;
	    }
	}

    if (k) {
        gks_fprintf (stderr, "GKS: unable to allocate %d of %d colors\n", k,
	    PRIVATE_COLORS);
	pick_closest_colors (p->color, PRIVATE_COLORS);
	}

    p->ncolors = PRIVATE_COLORS;
}


/*******************************************************/
/* 24/32-bit TrueColor display color 'allocation' code */
/*******************************************************/

static
int highbit(unsigned long ul)
{
  /* returns position of highest set bit in 'ul' as an integer (0-31),
   or -1 if none */

  int i;  unsigned long hb;
  hb = 0x8000;  hb = (hb<<16);  /* hb = 0x80000000UL */
  for (i=31; ((ul & hb) == 0) && i>=0;  i--, ul<<=1);
  return i;
}


static
void alloc_color(XColor *color)
{
    unsigned long r, g, b, rmask, gmask, bmask;
    int rshift, gshift, bshift;
    
    /* shift r,g,b so that high bit of 16-bit color specification is 
     * aligned with high bit of r,g,b-mask in visual, 
     * AND each component with its mask,
     * and OR the three components together
     */

    r = color->red;  g = color->green;  b = color->blue;

    rmask = p->vis->red_mask;
    gmask = p->vis->green_mask;
    bmask = p->vis->blue_mask;

    rshift = 15 - highbit(rmask);
    gshift = 15 - highbit(gmask);
    bshift = 15 - highbit(bmask);

    /* shift the bits around */
    if (rshift<0) r = r << (-rshift);
             else r = r >> rshift;

    if (gshift<0) g = g << (-gshift);
             else g = g >> gshift;

    if (bshift<0) b = b << (-bshift);
             else b = b >> bshift;

    r = r & rmask;
    g = g & gmask;
    b = b & bmask;

    color->pixel = r | g | b;

    /* put 'exact' colors into red,green,blue fields */
    /* shift the bits BACK to where they were, now that they've been masked */
    if (rshift<0) r = r >> (-rshift);
             else r = r << rshift;

    if (gshift<0) g = g >> (-gshift);
             else g = g << gshift;

    if (bshift<0) b = b >> (-bshift);
             else b = b << bshift;

    color->red = r;  color->green = g;  color->blue = b;
}


static
void allocate_shared_colors (void)
{
    int i, j, k;

    j = p->ncolors;
    
    for (i=k=0; i<p->scolors; i++)
	{
	p->color[j].red   = (unsigned short)(p->red[j]*65535);
	p->color[j].green = (unsigned short)(p->green[j]*65535);
	p->color[j].blue  = (unsigned short)(p->blue[j]*65535);

#if defined(__cplusplus) || defined(c_plusplus)
	if (p->vis->c_class == TrueColor)
#else
	if (p->vis->class == TrueColor)
#endif
	    {
	    alloc_color(&p->color[j]);
	    }
	else if (!XAllocColor (p->dpy, p->cmap, &p->color[j]))
	    {
	    p->color[j].pixel = 0xffff;
            k++;
	    }

	p->ncolors++;
	j++;
	}

    if (k)
        pick_closest_colors (p->color, p->ncolors);
}


static
void allocate_named_colors (char *rgb)
{
    int i, j, k, coli;
    FILE *file;
    char spec[80], *get;
    Status parse;
    Bool failed = False;

    if ((file = fopen (rgb, "r")) != NULL)
        {
        j = p->ncolors;

        for (i=k=0; i<p->scolors; i++)
	    {
            if ((get = fgets (spec, 80, file)) != NULL);
                {
                strtok (spec, "\n");

                if ((parse = XParseColor (p->dpy, p->cmap, spec, &p->color[j]))
		    != 0)
                    {
	            p->red[j]   = (float)p->color[j].red/65535.0;
	            p->green[j] = (float)p->color[j].green/65535.0;
	            p->blue[j]  = (float)p->color[j].blue/65535.0;
                    p->gray[j]  = 0.3*p->red[j] + 0.59*p->green[j] +
                        0.11*p->blue[j];

		    coli = j;
		    GSRGB (&coli, &p->red[j], &p->green[j], &p->blue[j]);
                    }
                else
                    failed = True;
                }

            if (!get || !parse)
                {
        	p->color[j].red   = (unsigned short)(p->red[j]*65535);
        	p->color[j].green = (unsigned short)(p->green[j]*65535);
        	p->color[j].blue  = (unsigned short)(p->blue[j]*65535);
                }

            if (!XAllocColor (p->dpy, p->cmap, &p->color[j]))
                {
                p->color[j].pixel = 0xffff;
                k++;
                }

            p->ncolors++;
            j++;
            }

        if (k)
            pick_closest_colors (p->color, p->ncolors);

        if (failed)
            gks_fprintf (stderr,
		"GKS: couldn't parse at least one color value\n");

        fclose (file);
        }
    else
        gks_fprintf (stderr, "GKS: can't open RGB file\n");
}


static
void allocate_colors (Bool min_colors)

/*
 *  Allocate colors
 */

{
    char *rgb;
    int i, pix;
    XVisualInfo vinfo;

    if (min_colors)
        {
	if (p->widget == NULL && p->wstype != 212 &&
	    XMatchVisualInfo(p->dpy, DefaultScreen(p->dpy), 24, TrueColor,
	    		     &vinfo) != 0)
	    {
	    p->vis = vinfo.visual;
	    p->depth = vinfo.depth;
	    p->cmap = XCreateColormap (p->dpy, RootWindow(p->dpy, vinfo.screen),
		p->vis, AllocNone);
	    }
	else
	    {
	    p->vis = XDefaultVisualOfScreen(p->screen);
	    p->depth = XDefaultDepthOfScreen(p->screen);
	    p->cmap = XDefaultColormapOfScreen(p->screen);
	    }

	i = p->dark_bg ? 1 : 0;
        p->color[i].pixel = XWhitePixelOfScreen(p->screen);
        p->color[1-i].pixel = XBlackPixelOfScreen(p->screen);
        p->ncolors = 2;

        p->pcolors = 0;
        }

    rgb = (char *) getenv("GLI_GKS_RGB");

    if (rgb == NULL)
        {
#if defined(__cplusplus) || defined(c_plusplus)
        switch (p->vis->c_class)
#else
        switch (p->vis->class)
#endif
            {	
            case PseudoColor:
	    case DirectColor:
	    case StaticGray:
                if (min_colors)
                    {
                    if (p->depth >= 4)
                        allocate_private_colors ();
                    }
                else {
	            if (p->depth >= 8)
			if (p->pcolors < PRIVATE_COLORS + p->scolors)
			    {
			    if (p->pcolors)
				allocate_shared_colors ();
			    else
				allocate_pseudo_colors ();
			    }
		    }
	        break;

	    case StaticColor:
	    case TrueColor:
	    case GrayScale:
                if (min_colors)
	            allocate_pseudo_colors ();
                else
                    allocate_shared_colors ();
	        break;
            }
        }
    else
        {
	if (min_colors)
            allocate_pseudo_colors ();
        else
            allocate_named_colors (rgb);
        }

    if (min_colors)
        {
        p->ccolor = Undefined;

        p->bg = p->color[0].pixel;
        p->fg = p->color[1].pixel;
        }
    else
        p->mono_flag = p->ncolors == 2;

    p->min_colors = min_colors;

    for (i=0; i<p->ncolors; i++)
        {
	pix = p->color[i].pixel;
	GSPIX (&i, &pix);
        }
}


static
void free_colors (void)
{
    unsigned long planes = 0;

    if (p->pcolors)
	{
	XFreeColors (p->dpy, p->cmap, p->pixels, p->pcolors, planes);
        XSync (p->dpy, False);
	}
}


static
void create_window (int win)

/*
 *  Create a window
 */

{
    XSetWindowAttributes xswa;
    XWindowAttributes xwa;
    char icon_name[40];
    char *env, **argv = NULL;
    int argc = 0;
    XSizeHints *hints = NULL;
    unsigned long valuemask;

    /* Set up the window attributes. We want to set the event mask,
     * and the background pixel */

    xswa.background_pixel = p->bg;
    xswa.event_mask = StructureNotifyMask | ExposureMask;

    if (p->backing_store && ((char *) getenv ("GLI_GKS_BS") == NULL))
	xswa.backing_store = Always;
    else
	xswa.backing_store = NotUseful;

    xswa.colormap = p->cmap;
    xswa.border_pixel = p->bg;

    if (p->widget == NULL && p->wstype != 212)
	{
	p->new_win = True;

	p->x =   5 + win*25;
	p->y = 100 + win*25;

	if (!p->uil)
	    {
	    if ((env = (char *) getenv ("GLI_GKS_MAGSTEP")) != NULL)
		p->magnification = pow(1.2, atof(env));

            p->width = p->height = (p->dpi == 100) ?
		(int)(667 * p->magnification) : (int)(500 * p->magnification);
	    }
        else
            p->width = p->height = 16;

	/* Create a window whose parent is the root window */

	p->win = XCreateWindow (p->dpy, XRootWindowOfScreen(p->screen),
	    p->x, p->y, p->width, p->height, 0, p->depth, InputOutput, p->vis,
	    CWBackPixel | CWEventMask | CWBackingStore | CWColormap |
	    CWBorderPixel, &xswa);

        XSelectInput (p->dpy, p->win, xswa.event_mask);

	p->icon_pixmap = XCreatePixmapFromBitmapData (p->dpy,
	    XRootWindowOfScreen(p->screen), (char *) icon_bits,
	    icon_width, icon_height, XBlackPixelOfScreen(p->screen),
	    XWhitePixelOfScreen(p->screen), 1);

	if (p->conid)
	    sprintf (icon_name, "GKSwk %d", p->conid);
	else
	    strcpy (icon_name, "GKSterm");

        XSetStandardProperties (p->dpy, p->win, WindowName, icon_name,
	    p->icon_pixmap, argv, argc, hints);
	
	XStoreName (p->dpy, p->win, WindowName);
	}
    else
	{
	p->new_win = False;

	if (p->wstype != 212)
	    p->win = XtWindow (p->widget);

	XGetWindowAttributes (p->dpy, p->win, &xwa);
	p->x      = xwa.x;
	p->y      = xwa.y;
	p->width  = xwa.width;
	p->height = xwa.height;

	xswa.event_mask |= (xwa.all_event_masks | ButtonPressMask);

	valuemask = CWBackingStore;
	if (p->wstype != 212)
	    valuemask = CWBackPixel | CWEventMask | CWBackingStore | CWColormap;

	XChangeWindowAttributes (p->dpy, p->win, valuemask, &xswa);
	}

    p->event_mask = xswa.event_mask;
}


static 
void set_WM_hints (void)
{
    XSizeHints hints;
    XWMHints wmhints;

    if (p->new_win) {
	hints.flags  = PPosition | PSize;
	hints.x	 = p->x;
        hints.y	 = p->y;
        hints.width  = p->width;
        hints.height = p->height;
    
        XSetNormalHints (p->dpy, p->win, &hints);

        if (p->gif >= 0 || p->rf >= 0) {
            wmhints.initial_state = IconicState;
            wmhints.flags = StateHint;

            XSetWMHints (p->dpy, p->win, &wmhints);
            }
	}
}


static
void create_GC (void)

/*
 *  Create graphics context
 */

{
    XGCValues xgcv;

    xgcv.foreground = p->fg;
    xgcv.background = p->bg;
    p->gc = XCreateGC (p->dpy, p->win, GCForeground | GCBackground, &xgcv);

    p->invert = XCreateGC (p->dpy, p->win, GCForeground | GCBackground, &xgcv);
    XSetFunction (p->dpy, p->invert, GXinvert);
    XSetForeground (p->dpy, p->invert, p->fg ^ p->bg);

    xgcv.foreground = p->bg;
    p->clear = XCreateGC (p->dpy, p->win, GCForeground | GCBackground, &xgcv);
}


static
void free_GC (void)

/*
 *  Free graphics context
 */

{
    XFreeGC (p->dpy, p->clear);
    XFreeGC (p->dpy, p->invert);
    XFreeGC (p->dpy, p->gc);
}


static
void create_pixmap (void)

/*
 *  Create a pixmap
 */

{
    if (!p->backing_store || p->gif >= 0 || p->rf >= 0 || p->uil || p->frame ||
	p->double_buf) {
	p->pixmap = XCreatePixmap (p->dpy, XRootWindowOfScreen(p->screen),
	    p->width, p->height, p->depth);

	XFillRectangle (p->dpy, p->pixmap, p->clear, 0, 0, p->width, p->height);
	}
    else
	p->pixmap = 0;
}


#ifdef XSHM

int XShmQueryExtension (Display *);

static
void create_shared_memory (void)

/*
 *  Create X shared memory
 */

{
    if (p->server != DECServer || !p->xshm)
	{
	p->shmimage = NULL;
	return;
	}

    if (XShmQueryExtension (p->dpy))
	{
	p->shmimage = XShmCreateImage (p->dpy, p->vis, p->depth,
	    p->mono_flag ? XYBitmap : ZPixmap, 0, &p->shminfo,
	    p->width, p->height);

	p->shminfo.shmid = shmget (IPC_PRIVATE, p->shmimage->bytes_per_line *
	    p->shmimage->height, IPC_CREAT | 0777);
        if (p->shminfo.shmid >= 0) {
            p->shminfo.shmaddr = (char *) shmat (p->shminfo.shmid, 0, 0);
	    }

	p->shminfo.readOnly = False;
        XShmAttach (p->dpy, &p->shminfo);
        XSync (p->dpy, False);

        shmctl (p->shminfo.shmid, IPC_RMID, 0);
        p->shmimage->data = p->shminfo.shmaddr;
	}
    else
	p->shmimage = NULL;
}


static
void free_shared_memory (void)

/*
 *  Free X shared memory
 */

{
    if (p->shmimage != NULL)
	{
	XShmDetach (p->dpy, &p->shminfo);
	XDestroyImage (p->shmimage);
	shmdt (p->shminfo.shmaddr);
	}
}

#endif /* XSHM */


static
void initialize_arrays (void)

/*
 *  Initialize the arrays
 */

{
    register int i, j;
    int pat, pa[33];

    if (!have_patterns) {
	for (i = 0; i < PATTERNS; i++) {
	    pat = i;
	    GKQPA (&pat, pa);
	    patterns[i][0] = (char)(*pa);
	    for (j = 1; j <= *pa; j++)
		patterns[i][j] = (char)(~pa[j]);
	    }
	have_patterns = True;
        }

    memset((void *)p->fstr, 0, n_font * MAX_SIZE * sizeof(XFontStruct *));
    memset((void *)p->tile, 0, MAX_COLORS * PATTERNS * sizeof(Pixmap));
    memset((void *)p->stipple, 0, MAX_COLORS * PATTERNS * sizeof(Pixmap));

    p->cached_size = 0;
}


static
void free_tile_patterns (int color)

/*
 *  Free tile patterns
 */

{
    int style;

    for (style=0; style<PATTERNS; style++)
        {
        if (p->tile[color][style] != 0) {
            XFreePixmap (p->dpy, p->tile[color][style]);
            XFreePixmap (p->dpy, p->stipple[color][style]);
	    p->tile[color][style] = p->stipple[color][style] = 0;
            }
        }
}


static
void create_cursor (void)

/*
 *  Create cursor
 */

{
    char *env;
    unsigned int shape = 0;

    if ((env = (char *) getenv ("GLI_XC")) != NULL)
	shape = atoi(env);
    if (!shape)
	shape = XC_draft_small;

    p->cursor = XCreateFontCursor (p->dpy, shape);
    p->textcursor = XCreateFontCursor (p->dpy, XC_xterm);
}


static
void set_color_repr (int i, float r, float g, float b)
{
    unsigned long pixel;
    XSetWindowAttributes xswa;
    int pix;

    if (i < 2 && p->dark_bg)
	return;

    if (i >= 8 && p->min_colors)
        allocate_colors (False);

    p->red[i]   = r;
    p->green[i] = g;
    p->blue[i]  = b;
    p->gray[i]  = 0.3*r + 0.59*g + 0.11*b;

    if (i < 2 && p->ncolors == 2 && p->wstype != 212)
	{
	if (p->gray[i] > 0.5)
	    pixel = XWhitePixelOfScreen(p->screen);
	else
	    pixel = XBlackPixelOfScreen(p->screen);

	if (pixel != p->color[i].pixel) {
	    p->color[i].pixel = pixel;

	    if (i == 0) {
		xswa.background_pixel = pixel;
		XChangeWindowAttributes (p->dpy, p->win, CWBackPixel, &xswa);

		XSetForeground (p->dpy, p->clear, pixel);
		if (p->pixmap)
		    XFillRectangle (p->dpy, p->pixmap, p->clear, 0, 0,
			p->width, p->height);
	    
		XClearWindow (p->dpy, p->win);
		}
	    }
	}

    else if (i <= p->ncolors)
	{
	p->color[i].red   = (unsigned short)(r*65535);
	p->color[i].green = (unsigned short)(g*65535);
	p->color[i].blue  = (unsigned short)(b*65535);

	if (i >= p->pcolors)
	    {
            if (!XAllocColor (p->dpy, p->cmap, &p->color[i]))
                {
                p->color[i].pixel = 0xffff;
                pick_closest_colors (p->color, p->ncolors);
	        }
	    }
	else
	    {
	    XStoreColor (p->dpy, p->cmap, &p->color[i]);
	    XSync (p->dpy, False);
	    }

	if (i == p->ncolors)
	    p->ncolors++;
	}

    if (i < 2)
	{
	p->bg = p->color[0].pixel;
	p->fg = p->color[1].pixel;
	XSetForeground (p->dpy, p->invert, p->fg ^ p->bg);
	}

    if (i < MAX_COLORS)
        {
	pix = p->color[i].pixel;
	GSPIX (&i, &pix);
        }

    p->ccolor = Undefined;
}


static
void set_color (int color)
{
    int i;

    i = Color8Bit(color);

    if (i >= 8 && p->min_colors)
        allocate_colors (False);

    if (i >= p->ncolors)
	{
	if (i >= 8)
	    i = 0;
	else
	    i = 1;
	}
    if (i != p->ccolor)
	{
	XSetForeground (p->dpy, p->gc, p->color[i].pixel);
	p->ccolor = i;
	}
}


static
void set_pattern (int color, int style)
{
    unsigned int w, h;
    char *pattern;

    if (color >= p->ncolors || color >= MAX_COLORS) color = 1;
    if (style >= PATTERNS) style = 1;

    if (style) {
	if (p->tile[color][style] == 0) {
	    pattern = patterns[style];
	    w = h = (*pattern == 32) ? 16 : *pattern;
	    pattern++;
	    p->tile[color][style] = XCreatePixmapFromBitmapData (p->dpy, p->win,
		pattern, w, h, p->color[color].pixel, p->bg, p->depth);
	    p->stipple[color][style] = XCreatePixmapFromBitmapData (p->dpy,
                p->win, pattern, w, h, p->color[color].pixel, p->bg, 1);
            }

        if (p->ored_patterns) {
            XSetFillStyle(p->dpy, p->gc, FillStippled);
            XSetStipple(p->dpy, p->gc, p->stipple[color][style]);
            }
        else {
            XSetFillStyle(p->dpy, p->gc, FillTiled);
            XSetTile(p->dpy, p->gc, p->tile[color][style]);
            }
	}
    else
	XSetFillStyle(p->dpy, p->gc, FillSolid);
}


static
void set_intensity (int color)
{
    int style;
    unsigned int w, h;
    char *pattern;

    if (color) {
	style = (int)(9*(1-p->gray[color-PRIVATE_COLORS])) + 1;
	color = 1;
	set_color (color);

	if (p->tile[color][style] == 0) {
	    pattern = patterns[style];
	    pattern = patterns[style];
	    w = h = (*pattern == 32) ? 16 : *pattern;
	    pattern++;
	    p->tile[color][style] = XCreatePixmapFromBitmapData (p->dpy, p->win,
		pattern, w, h, p->color[color].pixel, p->bg, p->depth);
	    p->stipple[color][style] = XCreatePixmapFromBitmapData (p->dpy,
                p->win, pattern, w, h, p->color[color].pixel, p->bg, 1);
            }

	XSetFillStyle(p->dpy, p->gc, FillTiled);
	XSetTile(p->dpy, p->gc, p->tile[color][style]);
	}
    else
	XSetFillStyle(p->dpy, p->gc, FillSolid);
}


static void configure_event (XConfigureEvent *event);

#if defined (ultrix) || defined(__osf__) || defined(aix) || defined(cray)

static
Bool predicate (Display *display, XEvent *event, char *arg)
{
    return (event->type == ConfigureNotify);
}

static
void process_async_events (void)
{
#ifndef VMS
    XEvent event;

    if (XCheckIfEvent (p->dpy, &event, predicate, NULL))
        configure_event ((XConfigureEvent *) &event);
#endif
}

#endif


static
void async_io (Bool yes)
{
#if defined (ultrix) || defined(__osf__) || defined(aix) || defined(cray)
    struct sigaction    act;
    pid_t               pid = getpid();
    int                 one = 1, zero = 0;

    act.sa_flags        = 0;
    act.sa_handler      = yes ? (void (*)(int)) process_async_events : SIG_IGN;
    (void) sigemptyset (&act.sa_mask);

    sigaction (SIGIO, &act, (struct sigaction *) NULL);

    if (yes) ioctl (p->fd, SIOCSPGRP, (char *) &pid);
    ioctl (p->fd, FIOASYNC, yes ? (char *) &one : (char *) &zero);

    p->async = yes;
#else
    p->async = False;
#endif
}


static
void setup_xform (float *window, float *viewport)
{
    p->a = (p->width-1)/(window[1]-window[0]);
    p->b = -window[0]*p->a;
    p->c = (p->height-1)/(window[2]-window[3]);
    p->d = p->height-1 - window[2]*p->c;
}


static
void configure_event (XConfigureEvent *event)

/*
 *  Handle configure events
 */

{
    float req_aspect_ratio, cur_aspect_ratio;
    int width, height;
    Bool async;

    if (p->widget || p->gif >= 0 || p->rf >= 0 || p->uil || p->frame)
        return;

    p->x = event->x;
    p->y = event->y;
    if (event->width == p->width && event->height == p->height)
        return;

    width  = event->width;
    height = event->height;

    p->viewport[0] = p->x*p->mwidth/p->swidth;
    p->viewport[1] = p->viewport[0] + width*p->mwidth/p->swidth;
    p->viewport[2] = (p->sheight-(p->y+height))*p->mheight/p->sheight;
    p->viewport[3] = p->viewport[2] + height*p->mheight/p->sheight;

    req_aspect_ratio = (p->window[1]-p->window[0])/
        (p->window[3]-p->window[2]);
    cur_aspect_ratio = (p->viewport[1]-p->viewport[0])/
        (p->viewport[3]-p->viewport[2]);

    if (cur_aspect_ratio > req_aspect_ratio) {
        width = (int)(height * req_aspect_ratio);
        p->viewport[1] = p->viewport[0] + (p->viewport[3]-p->viewport[2])*
            req_aspect_ratio;
        }
    else {
        height = (int)(width / req_aspect_ratio);
        p->viewport[3] = p->viewport[2] + (p->viewport[1]-p->viewport[0])/
            req_aspect_ratio;
        }

    if (width != p->width || height != p->height)
        {
        p->width  = width;
        p->height = height;

        if (p->pixmap) {
            XFreePixmap (p->dpy, p->pixmap);
            p->pixmap = XCreatePixmap (p->dpy, XRootWindowOfScreen(p->screen),
                p->width, p->height, p->depth);
            XFillRectangle (p->dpy, p->pixmap, p->clear, 0, 0,
                p->width, p->height);
            }
#ifdef XSHM
	free_shared_memory ();
	create_shared_memory ();
#endif
        setup_xform (p->window, p->viewport);
        set_clipping (True);

        async = p->async;
        GRSGWK (&p->wkid);
        if (async) {
            XSync (p->dpy, False);
            async_io (True);
            }

        return;
        }
    else
        return;
}


static
void async_expose_event (ws_state_list *p)

/*
 *  Handle asynchronous expose events
 */

{
    if (p->pixmap) {
        set_clipping (False);
	XCopyArea (p->dpy, p->pixmap, p->win, p->gc, 0, 0, p->width, p->height,
            0, 0);
        set_clipping (True);
	XSync (p->dpy, False);
	}
}


static
void select_async_input (int flag)
{
    long event_mask;
    
    if (p->wstype != 212) 
	{
	event_mask = p->event_mask;

	if (flag) {
#ifdef VMS
	    if (!p->widget && !p->backing_store)
		XSelectAsyncEvent (p->dpy, p->win, Expose, async_expose_event,
		    p);
#endif
	    }
	else {
#ifdef VMS
	    if (!p->widget && !p->backing_store)
		XSelectAsyncEvent (p->dpy, p->win, Expose, 0, 0);
#endif
	    event_mask |= ButtonPressMask | PointerMotionMask | KeyPressMask |
		KeyReleaseMask;
	    }

	XSelectInput (p->dpy, p->win, event_mask);
	}
}


static
void wait_for_expose (void)
{
    XEvent event;

    if (p->new_win) {
	do
	    XWindowEvent (p->dpy, p->win, StructureNotifyMask, &event);
	while (event.xany.type != MapNotify);
	while (XCheckTypedWindowEvent (p->dpy, p->win, Expose, &event))
	    ;
	}
}


static
void map_window (void)

/*
 *  Map window
 */

{
    /* Windows are not visible until they are mapped - map this window */

    if (!p->mapped) {
	XMapWindow (p->dpy, p->win);
	p->mapped = True;

	if (p->gif < 0 && p->rf < 0)
	    wait_for_expose ();
	
	select_async_input (True);
        if (p->widget && !p->backing_store)
	    XtAddEventHandler (p->widget, ExposureMask, False,
		(XtEventHandler) expose_event, p);
	}
}


static
void unmap_window (void)

/*
 *  Unmap window
 */

{
    if (p->mapped) {
	select_async_input (False);
	if (p->widget && !p->backing_store)
	    XtRemoveEventHandler (p->widget, ExposureMask, False,
		(XtEventHandler) expose_event, p);

	if (!p->widget)
	    XUnmapWindow (p->dpy, p->win);
	p->mapped = False;
	}
}


static
void setup_norm_xform (int tnr, float *wn, float *vp)
{
    a[tnr] = (vp[1]-vp[0])/(wn[1]-wn[0]);
    b[tnr] = vp[0]-wn[0]*a[tnr];
    c[tnr] = (vp[3]-vp[2])/(wn[3]-wn[2]);
    d[tnr] = vp[2]-wn[2]*c[tnr];
}


static
void init_norm_xform (void)
{
    int tnr;

    for (tnr=0; tnr<MAX_TNR; tnr++)
	setup_norm_xform (tnr, gksl->window[tnr], gksl->viewport[tnr]);
}


static
void draw_marker (float xn, float yn, int mtype, float mscale)
{
    int r, d, x, y, i;
    int pc, op;
    XPoint points[13];
    float scale, xr, yr;

    static int marker[26][57] = {
	{5,9,-4,7,4,7,7,4,7,-4,	    /* omark */
	4,-7,-4,-7,-7,-4,-7,4,
	-4,7, 3,9,-4,7,4,7,7,4,
	7,-4,4,-7,-4,-7,-7,-4,
	-7,4,-4,7, 0},
	{5,13,-2,8,2,8,2,2,8,2,	    /* hollow plus */
	 8,-2,2,-2,2,-8,-2,-8,
	 -2,-2,-8,-2,-8,2,-2,2,
	 -2,8, 3,13,-2,8,2,8,
	 2,2,8,2,8,-2,2,-2,2,-8,
	 -2,-8,-2,-2,-8,-2,-8,2,
	 -2,2,-2,8, 0},
	{4,4,-8,0,4,7,4,-7,	    /* solid triangle right */
	 -8,0, 0},
	{4,4,8,0,-4,-7,-4,7,	    /* solid triangle left */
	 8,0, 0},
	{5,4,0,8,7,-4,-7,-4,0,8,    /* triangle up down */
	 5,4,0,-8,-7,4,7,4,0,-8,
	 3,4,0,8,7,-4,-7,-4,0,8,
	 3,4,0,-8,-7,4,7,4,0,-8,
	 0},
	{4,11,0,9,2,2,9,3,3,-1,	    /* solid star */
	 6,-8,0,-3,-6,-8,-3,-1,
	 -9,3,-2,2,0,9, 0},
	{5,11,0,9,2,2,9,3,3,-1,	    /* hollow star */
	 6,-8,0,-3,-6,-8,-3,-1,
	 -9,3,-2,2,0,9,
	 3,11,0,9,2,2,9,3,3,-1,
	 6,-8,0,-3,-6,-8,-3,-1,
	 -9,3,-2,2,0,9, 0},
	{4,5,0,9,9,0,0,-9,-9,0,	    /* solid diamond */
	 0,9, 0},
	{5,5,0,9,9,0,0,-9,-9,0,	    /* hollow diamond */
	 0,9, 3,5,0,9,9,0,0,-9,
	 -9,0,0,9, 0},
	{4,5,9,9,-9,-9,9,-9,-9,9,   /* solid hourglass */
	 9,9, 0},
	{5,5,9,9,-9,-9,9,-9,-9,9,   /* hollow hourglass */
	 9,9, 3,5,9,9,-9,-9,9,-9,
	 -9,9,9,9, 0},
	{4,5,9,9,9,-9,-9,9,-9,-9,   /* solid bowtie */
	 9,9, 0},
	{5,5,9,9,9,-9,-9,9,-9,-9,   /* hollow bowtie */
	 9,9, 3,5,9,9,9,-9,-9,9,
	 -9,-9,9,9, 0},
	{4,5,9,9,9,-9,-9,-9,-9,9,   /* solid square */
	 9,9, 0},
	{5,5,9,9,9,-9,-9,-9,-9,9,   /* hollow square */
	 9,9, 3,5,9,9,9,-9,-9,-9,
         -9,9,9,9, 0},
	{4,4,-9,9,9,9,0,-9,-9,9,    /* solid triangle down */
	 0},
	{5,4,-9,9,9,9,0,-9,-9,9,    /* hollow triangle down */
	 3,4,-9,9,9,9,0,-9,-9,9,
	 0},
	{4,4,0,9,9,-9,-9,-9,0,9,    /* solid triangle up */
	 0},
	{5,4,0,9,9,-9,-9,-9,0,9,    /* hollow triangle up */
	 3,4,0,9,9,-9,-9,-9,0,9, 0},
	{7,0,360, 0},		    /* solid circle */
	{0},			    /* not used */
	{1, 0},			    /* dot */
	{2, 0,0,0,9, 2,0,0,9,0,	    /* plus */
	 2,0,0,0,-9, 2,0,0,-9,0, 
	 0},
	{2,0,0,0,9, 2,0,0,9,3,	    /* asterisk */
	 2,0,0,6,-9, 2,0,0,-6,-9, 
	 2,0,0,-9,3,  0},
	{8,0,360, 6,0,360, 0},	    /* circle */
	{2,0,0,9,9, 2,0,0,9,-9,	    /* diagonal cross */
	 2,0,0,-9,-9, 2,0,0,-9,9,
	 0}};

    r = (int)(3*mscale);
    d = 2*r;
    scale = mscale/3.0;

    NDC_to_DC(xn, yn, x, y);

    pc = 0;    
    mtype = (d > 1) ? mtype + 20 : 21;

    do 
	{
	op = marker[mtype][pc];
	switch (op) {

	    case 1: /* point */
		if (p->pixmap)
		    XDrawPoint (p->dpy, p->pixmap, p->gc, x, y);
		if (!p->double_buf)
		    XDrawPoint (p->dpy, p->win, p->gc, x, y);
		break;

	    case 2: /* line */
		for (i = 0; i < 2; i++)
		    {
                    xr =  scale * marker[mtype][pc+2*i+1];
                    yr = -scale * marker[mtype][pc+2*i+2];
                    seg_xform_rel (&xr, &yr);
                    points[i].x = nint(x-xr);
                    points[i].y = nint(y+yr);
		    }
		if (p->pixmap)
		    XDrawLines (p->dpy, p->pixmap, p->gc, points, 2,
                        CoordModeOrigin);
		if (!p->double_buf)
		    XDrawLines (p->dpy, p->win, p->gc, points, 2,
			CoordModeOrigin);
		pc += 4;
		break;

	    case 3: /* polygon */
		for (i = 0; i < marker[mtype][pc+1]; i++)
		    {
		    xr =  scale * marker[mtype][pc+2+2*i];
		    yr = -scale * marker[mtype][pc+3+2*i];
                    seg_xform_rel (&xr, &yr);
                    points[i].x = nint(x-xr);
                    points[i].y = nint(y+yr);
		    }
		if (p->pixmap)
		    XDrawLines (p->dpy, p->pixmap, p->gc, points,
			marker[mtype][pc+1], CoordModeOrigin);
		if (!p->double_buf)
		    XDrawLines (p->dpy, p->win, p->gc, points,
			marker[mtype][pc+1], CoordModeOrigin);
		pc += 1 + 2 * marker[mtype][pc+1];
		break;

	    case 4: /* filled polygon */
		for (i = 0; i < marker[mtype][pc+1]; i++)
		    {
		    xr =  scale * marker[mtype][pc+2+2*i];
		    yr = -scale * marker[mtype][pc+3+2*i];
                    seg_xform_rel (&xr, &yr);
                    points[i].x = nint(x-xr);
                    points[i].y = nint(y+yr);
		    }
		if (p->pixmap)
		    XFillPolygon (p->dpy, p->pixmap, p->gc, points,
			marker[mtype][pc+1], Complex, CoordModeOrigin);
		if (!p->double_buf)
		    XFillPolygon (p->dpy, p->win, p->gc, points,
			marker[mtype][pc+1], Complex, CoordModeOrigin);
		pc += 1 + 2 * marker[mtype][pc+1];
		break;

	    case 5: /* hollow polygon */
		for (i = 0; i < marker[mtype][pc+1]; i++)
		    {
		    xr =  scale * marker[mtype][pc+2+2*i];
		    yr = -scale * marker[mtype][pc+3+2*i];
                    seg_xform_rel (&xr, &yr);
                    points[i].x = nint(x-xr);
                    points[i].y = nint(y+yr);
		    }
		if (p->pixmap)
		    XFillPolygon (p->dpy, p->pixmap, p->clear, points,
			marker[mtype][pc+1], Complex, CoordModeOrigin);
		if (!p->double_buf)
		    XFillPolygon (p->dpy, p->win, p->clear, points,
			marker[mtype][pc+1], Complex, CoordModeOrigin);
		pc += 1 + 2 * marker[mtype][pc+1];
		break;

	    case 6: /* arc */
		if (p->pixmap)
		    XDrawArc (p->dpy, p->pixmap, p->gc, x-r, y-r, d, d, 
			marker[mtype][pc+1]*64, marker[mtype][pc+2]*64);
		if (!p->double_buf)
		    XDrawArc (p->dpy, p->win, p->gc, x-r, y-r, d, d, 
			marker[mtype][pc+1]*64, marker[mtype][pc+2]*64);
		pc += 2;
		break;

	    case 7: /* filled arc */
		if (p->pixmap)
		    XFillArc (p->dpy, p->pixmap, p->gc, x-r, y-r, d, d, 
			marker[mtype][pc+1]*64, marker[mtype][pc+2]*64);
		if (!p->double_buf)
		    XFillArc (p->dpy, p->win, p->gc, x-r, y-r, d, d, 
			marker[mtype][pc+1]*64, marker[mtype][pc+2]*64);
		pc += 2;
		break;

	    case 8: /* hollow arc */
		if (p->pixmap)
		    XFillArc (p->dpy, p->pixmap, p->clear, x-r, y-r, d, d,
			marker[mtype][pc+1]*64, marker[mtype][pc+2]*64);
		if (!p->double_buf)
		    XFillArc (p->dpy, p->win, p->clear, x-r, y-r, d, d,
			marker[mtype][pc+1]*64, marker[mtype][pc+2]*64);
		pc += 2;
		break;
	    }
	pc++;
	}
    while (op != 0);
}


static
void draw_points (int *n, float *px, float *py, int *tnr)
{
    register int i;
    float xn, yn;

    if (*n > max_points)
	{
	points = (XPoint *) realloc (points, *n * sizeof(XPoint));
	max_points = *n;
	}

    for (i=0; i<*n; i++)
	{
        WC_to_NDC(px[i], py[i], *tnr, xn, yn);
        seg_xform(&xn, &yn);
	NDC_to_DC(xn, yn, points[i].x, points[i].y);
	}

    if (p->pixmap)
	XDrawPoints (p->dpy, p->pixmap, p->gc, points, *n, CoordModeOrigin);
    if (!p->double_buf)
	XDrawPoints (p->dpy, p->win, p->gc, points, *n, CoordModeOrigin);
}


static
void marker_routine (int *n, float *px, float *py, int *tnr, int mtype,
    float mscale)
{
    float clrt[4], x, y;
    register int i;
    register Bool draw;

    if (gksl->clip == GCLIP || mtype != GPOINT)
      {
      if (gksl->clip == GCLIP)
	{
        memcpy(clrt, gksl->viewport[gksl->cntnr], 4*sizeof(float));
	if (p->gr_bc == 0)
	    {
	    seg_xform (&clrt[0], &clrt[2]);
	    seg_xform (&clrt[1], &clrt[3]);
	    }
	}
      set_clipping (False);
      for (i=0; i<*n; i++)
	{
	WC_to_NDC(px[i], py[i], *tnr, x, y);
        seg_xform(&x, &y);

	if (gksl->clip == GCLIP)
	  draw = (x >= clrt[0] && x <= clrt[1] && y >= clrt[2] && y <= clrt[3]);
	else
	  draw = True;

	if (draw)
	  draw_marker (x, y, mtype, mscale);
	}
      set_clipping (True);
      }
    else
      draw_points (n, px, py, tnr);
}


static
void set_line_attr (int linetype, float linewidth)
{
    static char dash_list[35][10] = {
	{8,  4, 2, 4, 2, 4, 2, 4, 6, 0},
	{6,  4, 2, 4, 2, 4, 6, 0, 0, 0},
	{4,  4, 2, 4, 6, 0, 0, 0, 0, 0},
	{8,  3, 2, 3, 2, 3, 2, 3, 6, 0},
	{6,  3, 2, 3, 2, 3, 6, 0, 0, 0},
	{4,  3, 2, 3, 6, 0, 0, 0, 0, 0},
	{8,  3, 2, 3, 2, 3, 2, 3, 4, 0},
	{6,  3, 2, 3, 2, 3, 4, 0, 0, 0},
	{4,  3, 2, 3, 4, 0, 0, 0, 0, 0},
	{2,  1, 1, 0, 0, 0, 0, 0, 0, 0},
	{2,  1, 2, 0, 0, 0, 0, 0, 0, 0},
	{2,  1, 6, 0, 0, 0, 0, 0, 0, 0},
	{2,  1, 8, 0, 0, 0, 0, 0, 0, 0},
	{6,  1, 3, 1, 3, 1, 6, 0, 0, 0},
	{4,  1, 3, 1, 6, 0, 0, 0, 0, 0},
	{8,  6, 2, 1, 2, 1, 2, 1, 2, 0},
	{6,  6, 2, 1, 2, 1, 2, 0, 0, 0},
	{4,  6, 2, 1, 2, 0, 0, 0, 0, 0},
	{4,  9, 3, 5, 3, 0, 0, 0, 0, 0},
	{2,  9, 3, 0, 0, 0, 0, 0, 0, 0},
	{2,  5, 5, 0, 0, 0, 0, 0, 0, 0},
	{2,  5, 3, 0, 0, 0, 0, 0, 0, 0},
	{6,  1, 4, 1, 4, 1, 8, 0, 0, 0},
	{4,  1, 4, 1, 8, 0, 0, 0, 0, 0},
	{2,  1, 1, 0, 0, 0, 0, 0, 0, 0},
	{2,  8, 1, 0, 0, 0, 0, 0, 0, 0},
	{4, 16, 5, 8, 5, 0, 0, 0, 0, 0},
	{2, 16, 5, 0, 0, 0, 0, 0, 0, 0},
	{8,  8, 4, 1, 4, 1, 4, 1, 4, 0},
	{6,  8, 4, 1, 4, 1, 4, 0, 0, 0},
	{0,  0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0,  0, 0, 0, 0, 0, 0, 0, 0, 0},
	{2,  8, 5, 0, 0, 0, 0, 0, 0, 0},
	{2,  1, 2, 0, 0, 0, 0, 0, 0, 0},
	{4,  8, 4, 1, 4, 0, 0, 0, 0, 0}};

    unsigned int width;
    int n;

    if (linewidth > 1)
        width = (unsigned int) (linewidth);
    else
        width = 0;

    if (linetype != p->ltype || width != p->lwidth)
        {
        if (linetype != GLSOLI) {
            XSetLineAttributes (p->dpy, p->gc, width, LineOnOffDash,
                CapNotLast, JoinRound);
	    n = (int) dash_list[linetype+30][0];
            XSetDashes (p->dpy, p->gc, 0, &dash_list[linetype+30][1], n);
	    }
        else
            XSetLineAttributes (p->dpy, p->gc, width, LineSolid, CapNotLast,
	        JoinRound);

        p->ltype = linetype;
        p->lwidth = width;
        }
}


static
int clip_code (int x, int y)
{
    register int code = 0;

    if (x < 0)
        code = LEFT;
    else if (x > p->width)
        code = RIGHT;

    if (y < 0)
        code |= BOTTOM;
    else if (y > p->height)
        code |= TOP;

    return code;
}


static
void clip_line (int *x0, int *x1, int *y0, int *y1, Bool *visible, Bool *clip)
{
    register int c, c0, c1;
    register int x, y;

    c0 = clip_code (*x0, *y0);
    c1 = clip_code (*x1, *y1);

    *clip = c1;
    *visible = False;

    while (c0 | c1) {
        if (c0 & c1) return;
        c = c0 ? c0 : c1;

        if (c & LEFT) {
            x = 0;
            y = (int) (*y0 - (*y1-*y0) * (float) (*x0) / (*x1-*x0));
            }
        else if (c & RIGHT) {
            x = p->width;
            y = (int) (*y0 + (*y1-*y0) * (float) (p->width-*x0) / (*x1-*x0));
            }
        else if (c & BOTTOM) {
            x = (int) (*x0 - (*x1-*x0) * (float) (*y0) / (*y1-*y0));
            y = 0;
            }
        else if (c & TOP) {
            x = (int) (*x0 + (*x1-*x0) * (float) (p->height-*y0) / (*y1-*y0));
            y = p->height;
            }

        if (c == c0) {
            *x0 = x;
            *y0 = y;
            c0 = clip_code (x, y);
            }
        else {
            *x1 = x;
            *y1 = y;
            c1 = clip_code (x, y);
            }
        }
    *visible = True;
}


static
void line_routine (int *n, float *px, float *py, int *linetype, int *tnr)
{
    float x1, y1;
    register int i, j, npoints, m;
    int ix0, iy0, ix1, iy1, x, y;
    Bool visible, clip;

    if (*n > max_points)
	{
	points = (XPoint *) realloc (points, *n * sizeof(XPoint));
	max_points = *n;
	}

    WC_to_NDC(px[0], py[0], *tnr, x1, y1);
    seg_xform(&x1, &y1);
    NDC_to_DC(x1, y1, ix1, iy1);

    npoints = 0;
    m = *linetype ? *n : *n+1;

    for (j = 1;  j < m;  j++) {
        i = j<*n ? j : 0;

	ix0 = ix1;
        iy0 = iy1;

        WC_to_NDC(px[i], py[i], *tnr, x1, y1);
        seg_xform(&x1, &y1);
        NDC_to_DC(x1, y1, ix1, iy1);

        x = ix1;
        y = iy1;
        clip_line (&ix0, &ix1, &iy0, &iy1, &visible, &clip);

        if (visible) {
            if (!npoints) {
                points[0].x = ix0;
                points[0].y = iy0;
                npoints = 1;
                }

            points[npoints].x = ix1;
            points[npoints].y = iy1;
            npoints++;

            if (clip) {
                if (p->pixmap)
                    XDrawLines (p->dpy, p->pixmap, p->gc, points, npoints, 0);
		if (!p->double_buf)
		    XDrawLines (p->dpy, p->win, p->gc, points, npoints, 0);
                npoints = 0;
                }

            if (npoints == MAX_POINTS) {
                if (p->pixmap)
                    XDrawLines (p->dpy, p->pixmap, p->gc, points, npoints, 0);
		if (!p->double_buf)
		    XDrawLines (p->dpy, p->win, p->gc, points, npoints, 0);

		--npoints;
                points[0].x = points[npoints].x;
                points[0].y = points[npoints].y;
                npoints = 1;
                }
            }

        ix1 = x;
        iy1 = y;
        }

    if (npoints > 1) {
	if (p->pixmap)
	    XDrawLines (p->dpy, p->pixmap, p->gc, points, npoints, 0);
	if (!p->double_buf)
	    XDrawLines (p->dpy, p->win, p->gc, points, npoints, 0);
	}
}


static
void polyline (int *n, float *px, float *py)
{
    int ln_type, ln_color;
    float ln_width;

    ln_type  = gksl->asf[0] ? gksl->ltype  : gksl->lindex;
    ln_width = gksl->asf[1] ? gksl->lwidth : 1;
    ln_color = gksl->asf[2] ? gksl->plcoli : 1;

    set_color (ln_color);
    set_line_attr (ln_type, ln_width);

    line_routine (n, px, py, &ln_type, &gksl->cntnr);
}


static
void polymarker (int *n, float *px, float *py)
{
    int mk_type, mk_color;
    float mk_size;

    mk_type  = gksl->asf[3] ? gksl->mtype  : gksl->mindex;
    mk_size  = gksl->asf[4] ? gksl->mszsc  : 1;
    mk_color = gksl->asf[5] ? gksl->pmcoli : 1;

    set_color (mk_color);
    set_line_attr (GLSOLI, 1.0);

    marker_routine (n, px, py, &gksl->cntnr, mk_type, mk_size);
}


static
void fill_routine (int *n, float *px, float *py, int *tnr)
{
    float x, y;
    register int i, npoints;

    if (*n > max_points)
	{
	points = (XPoint *) realloc (points, *n * sizeof(XPoint));
	max_points = *n;
	}

    npoints = *n;
    for (i=0; i<*n; i++)
	{
        WC_to_NDC(px[i], py[i], *tnr, x, y);
        seg_xform(&x, &y);
        NDC_to_DC(x, y, points[i].x, points[i].y);
	}

    if (npoints > 1)
        {
	if (p->pixmap)
	    XFillPolygon (p->dpy, p->pixmap, p->gc, points, npoints,
                p->shape, CoordModeOrigin);
	if (!p->double_buf)
	    XFillPolygon (p->dpy, p->win, p->gc, points, npoints,
		p->shape, CoordModeOrigin);
	}
}


static
void fill_area (int *n, float *px, float *py)
{
    int fl_inter, fl_style, fl_color;
    int ln_type;

    fl_inter = gksl->asf[10] ? gksl->ints   : predef_ints[gksl->findex - 1];
    fl_style = gksl->asf[11] ? gksl->styli  : predef_styli[gksl->findex - 1];
    fl_color = gksl->asf[12] ? gksl->facoli : 1;

    set_color (fl_color);
    set_line_attr (GLSOLI, 1.0);

    if (fl_inter == GSOLID) {
	if (fl_color < p->ncolors || fl_color >= MAX_COLORS)
	    fill_routine (n, px, py, &gksl->cntnr);
	else if (fl_color > 7) {
	    set_intensity (fl_color);
	    fill_routine (n, px, py, &gksl->cntnr);
	    set_intensity (0);
	    }
	else {
	    ln_type = DrawBorder;
	    line_routine (n, px, py, &ln_type, &gksl->cntnr);
	    }
	}
    else if (fl_inter == GPATTR || fl_inter == GHATCH) {
        if (fl_inter == GHATCH)
            set_pattern (fl_color, fl_style + HATCH_STYLE);
        else
            set_pattern (fl_color, fl_style);
	fill_routine (n, px, py, &gksl->cntnr);
	set_pattern (fl_color, 0);
	}
    else {
        ln_type = DrawBorder;
        line_routine (n, px, py, &ln_type, &gksl->cntnr);
        }
}


static
void draw_string (int x, int y, int width, int height, char *chars, int nchars)
{
    Pixmap src, dest;
    register XImage *from, *to;
    register int descent, w, h;
    register int i, j, ii, jj;
    register unsigned long pixel;

    height += 8;
    descent = p->cfont->descent;

    switch (p->path) {
	case 0: break;
	case 1: x = x-height+descent; y -= width; w = height; h = width; break;
	case 2: x -= width; y -= descent; w = width; h = height; break;
	case 3: x -= descent; w = height; h = width; break;
	}

    if (p->path != 0) {
        set_clipping (False);

	src = XCreatePixmap (p->dpy, XRootWindowOfScreen(p->screen),
	    width, height, p->depth);

	XFillRectangle (p->dpy, src, p->clear, 0, 0, width, height);
        XDrawString (p->dpy, src, p->gc, 0, height - descent, chars, nchars);

	dest = XCreatePixmap (p->dpy, XRootWindowOfScreen(p->screen),
	    w, h, p->depth);
	XCopyArea (p->dpy, p->double_buf ? p->pixmap : p->win, dest, p->gc,
	    x, y, w, h, 0, 0);

	from = XGetImage (p->dpy, src, 0, 0, width, height, AllPlanes, ZPixmap);
	to = XGetImage (p->dpy, dest, 0, 0, w, h, AllPlanes, ZPixmap);

        for (i=0; i<width; i++) 
	    for (j=0; j<height; j++) {
		switch (p->path) {
		    case 1: ii = j; jj = h-i-1; break;
		    case 2: ii = w-i-1; jj = h-j-1; break;
		    case 3: ii = w-j-1; jj = i; break;
		    }
		pixel = XGetPixel (from, i, j);
		if (pixel != p->bg)
		    XPutPixel (to, ii, jj, pixel);
		}

        set_clipping (True);

	if (p->pixmap)
	    XPutImage (p->dpy, p->pixmap, p->gc, to, 0, 0, x, y, w, h);
	if (!p->double_buf)
	    XPutImage (p->dpy, p->win, p->gc, to, 0, 0, x, y, w, h);

	XDestroyImage (to);
	XFreePixmap (p->dpy, dest);
	XDestroyImage (from);
	XFreePixmap (p->dpy, src);
	}
    else
	{
	if (p->pixmap)
	    XDrawString (p->dpy, p->pixmap, p->gc, x, y, chars, nchars);
	if (!p->double_buf)
	    XDrawString (p->dpy, p->win, p->gc, x, y, chars, nchars);
	}
}


static
void text_routine (float *x, float *y, int *nchars,
#ifdef VMS
    struct dsc$descriptor *text)
#else
#ifdef cray
    _fcd text)
#else
    char *text)
#endif /* cray */
#endif /* VMS */

{
    char *chars;
    int xorg, yorg, width, height;
    float xrel, yrel, ax, ay;
    int tx_prec;

#ifdef VMS
    chars = text->dsc$a_pointer;
    *nchars = text->dsc$w_length;
#else
#ifdef cray
    chars = _fcdtocp(text);
    *nchars = _fcdlen(text);
#else
    chars = text;
#endif /* cray */
#endif /* VMS */

    NDC_to_DC(*x, *y, xorg, yorg);

    /* Compute text extent */

    width = XTextWidth (p->cfont, chars, *nchars);
    height = p->capheight;

    tx_prec  = gksl->asf[6] ? gksl->txprec : predef_prec[gksl->tindex - 1];

    if (tx_prec == GSTRP) {

	/* Align the text */

	xrel = width * xfac[gksl->txal[0]];
	if (gksl->txal[1] != GABOTT)
	    yrel = height * yfac[gksl->txal[1]];
	else
	    yrel = p->cfont->descent;

	CharXform(xrel, yrel, ax, ay);

	xorg += (int) ax;
	yorg -= (int) ay;
	}

    draw_string (xorg, yorg, width, height, chars, *nchars);
}


static
void set_font (int font)
{
    int family, size, angle;
    char fontname[256];
    float scale, ux, uy, rad;
    float width, height, capheight, points;
    int fontsize;

    font = abs(font);
    if (font >= 101 && font <= 129)
        font -= 100;
    else if (font >= 1 && font <= 32)
        font = map[font-1];
    else
        font = 9;
    family = font-1;

    WC_to_NDC_rel(gksl->chup[0], gksl->chup[1], gksl->cntnr, ux, uy);
    seg_xform_rel(&ux, &uy);

    rad = -atan2(ux,uy);
    angle = (int)(rad*180/Pi + 0.5);
    if (angle < 0) angle += 360;
    p->path = ((angle+45) / 90) % 4;

    scale = sqrt(gksl->chup[0]*gksl->chup[0] + gksl->chup[1]*gksl->chup[1]);
    ux = gksl->chup[0] / scale * gksl->chh;
    uy = gksl->chup[1] / scale * gksl->chh;
    WC_to_NDC_rel(ux, uy, gksl->cntnr, ux, uy);

    width = 0;
    height = sqrt(ux*ux + uy*uy);
    seg_xform_rel(&width, &height);

    height = sqrt(width*width + height*height);
    capheight = height * (fabs(p->c) + 1);

    if (p->server == DECServer || p->server == IBMServer)
	points = 10 * capheight / capheights[family];
    else
	points = 10 * capheight * 100.0/p->dpi;

    p->capheight = nint(capheight);

    if (p->server == DECServer || p->server == IBMServer) {
	size = nint(points/10);
	if (size == 0)
	    size = 1;
	p->ratio = points/(size * 10);
	if (size >= MAX_SIZE)
	    {
	    fontsize = size;
	    size = MAX_SIZE - 1;
	    if (p->cached_size != fontsize)
		{
		p->fstr[family][size] = NULL;
		p->cached_size = fontsize;
		}
	    }
	else
	    fontsize = size;
	}
    else {
    	if (points <= 90) {
	    size = 8; p->ratio = points/80;
    	} else if (points <= 110) {
	    size = 10; p->ratio = points/100;
    	} else if (points <= 130) {
	    size = 12; p->ratio = points/120;
    	} else if (points <= 160) {
	    size = 14; p->ratio = points/140;
    	} else if (points <= 210) {
	    size = 18; p->ratio = points/180;
    	} else {
	    size = 24; p->ratio = points/240;
    	}
	fontsize = size;
    }
    
    if (p->fstr[family][size] == NULL) {
	sprintf (fontname, fonts[family], fontsize, p->dpi, p->dpi);
	p->fstr[family][size] = XLoadQueryFont (p->dpy, fontname);
	
	if (p->fstr[family][size] == NULL) {
	    gks_fprintf (stderr, "GKS: unable to load font %s\n", fontname);
	    p->fstr[family][size] = XLoadQueryFont (p->dpy, "variable");

	    if (p->fstr[family][size] == NULL)
		p->fstr[family][size] = XLoadQueryFont (p->dpy,
		    "fixed");
	    }
	}

    if (p->fstr[family][size] != NULL) {
	p->cfont = p->fstr[family][size];
	XSetFont (p->dpy, p->gc, p->cfont->fid);
	}
}


static
void text (float *px, float *py, int nchars,
#ifdef VMS
    struct dsc$descriptor *chars)
#else
#ifdef cray
    _fcd chars)
#else
    char *chars)
#endif /* cray */
#endif /* VMS */

{
    int tx_font, tx_prec, tx_color;
    Bool avail = True;
    float x, y;

    tx_font  = gksl->asf[6] ? gksl->txfont : predef_font[gksl->tindex - 1];
    tx_prec  = gksl->asf[6] ? gksl->txprec : predef_prec[gksl->tindex - 1];
    tx_color = gksl->asf[9] ? gksl->txcoli : 1;

    if (tx_prec != GSTRKP) {
	set_font (tx_font);
	avail = (p->ratio > 0.8 && p->ratio < 1.3);
	}

    set_color (tx_color);
    set_line_attr (GLSOLI, 1.0);

    if (tx_prec == GSTRP && avail)
	{
        WC_to_NDC(*px, *py, gksl->cntnr, x, y);
        seg_xform(&x, &y);

	text_routine (&x, &y, &nchars, chars);
	}
    else
#if defined(hpux) && !defined(NAGware)
	GTEXTS (px, py, &nchars, chars, &avail,
#if defined(__hp9000s700) || defined(__hp9000s300)
	    line_routine_a, fill_routine_a, text_routine_a);
#else
	    &line_routine_a, &fill_routine_a, &text_routine_a);
#endif /* __hp9000s700 || __hp9000s300 */
#else
	GTEXTS (px, py, &nchars, chars, &avail, line_routine, fill_routine,
            text_routine);
#endif
}


static
void update (void)
{
    if (p->state == GACTIV)
        {
#ifndef VMS
        if (!p->widget && !p->backing_store)
            {
            XEvent event;
            while (XPending (p->dpy)) {
                XNextEvent (p->dpy, &event);
                if (event.type == Expose)
                    expose_event (p->widget, p, (XExposeEvent *) &event, NULL);
                }
            }
        else
#endif
            XSync (p->dpy, False);
	}
}


static
void message (int nchars, char *chars)
{
    XDrawString (p->dpy, p->win, p->invert, 10, 20, chars, nchars);
}


static
void open_gif (int stream)
{
    char *gif_name;
    struct stat buf;
    int fd;

    fd = stream > 100 ? stream - 100 : stream;

    if (isatty(fd) || fstat(fd, &buf) != 0)
        {
        gif_name = (char *) getenv ("GLI_GIF");
        if (!gif_name)
            gif_name = "gli.gif";

        p->gif = open (gif_name, O_CREAT | O_TRUNC | O_WRONLY, 0644);
        if (p->gif < 0)
            gks_fprintf (stderr, "GKS: can't open GIF file\n");
        }
    else
        p->gif = fd;
}


static
void write_gif_word (int w)
{
    byte c;

    c = (w & 0xff);
    write (p->gif, &c, 1);
    c = ((w >> 8) & 0xff);
    write (p->gif, &c, 1);
}


static
void pixmap_to_gif (void)
{
    int size, besize;
    byte c, r, g, b, *pix, *ppix, *beimage;
    register int i, j, k, coli, mcolor;
    XImage *image;
    int BitsPerPixel, ColorMapSize, InitCodeSize;
    unsigned long pixel;

    image = XGetImage (p->dpy, p->pixmap, 0, 0, p->width, p->height, AllPlanes,
        ZPixmap);

    size = p->width * p->height;
    pix = ppix = (byte *) malloc (sizeof(byte) * size);
    beimage = (byte *) malloc (sizeof(byte) * size * 3 / 2);   /* worst case */

    if (pix != NULL && beimage != NULL)
	{
	mcolor = 0;
        for (j=0; j<p->height; j++) {
            for (i=0; i<p->width; i++) {
                pixel = XGetPixel (image, i, j);
                coli = 0;
                for (k=0; k<p->ncolors; k++) {
                    if (pixel == p->color[k].pixel) {
                        coli = k;
                        break;
                        }
                    }
                *ppix++ = coli;
		if (coli > mcolor)
		    mcolor = coli;
                }
            }

	for (BitsPerPixel = 1; BitsPerPixel < 8; BitsPerPixel++)
	    if ((1 << BitsPerPixel) > mcolor) break;

	/* write the GIF header */

	write (p->gif, p->wstype == 218 ? "GIF89a" : "GIF87a", 6);

	write_gif_word (p->width);	/* screen descriptor */
	write_gif_word (p->height);

	c = 0x80;			/* yes, there is a color map */
	c |= (8 - 1) << 4;		/* OR in the color resolution (8) */
	c |= (0) << 3;			/* no, the colors are not sorted */
	c |= (BitsPerPixel - 1);	/* OR in the # of bits per pixel */
	write (p->gif, &c, 1);

	c = 0x0;
	write (p->gif, &c, 1);		/* background color index */
	write (p->gif, &c, 1);		/* pixel aspect ratio */

	/* write colormap */

	ColorMapSize = 1 << BitsPerPixel;

	for (i = 0; i < ColorMapSize; i++) {
	    r = (byte) (255 * p->red[i]);
	    g = (byte) (255 * p->green[i]);
	    b = (byte) (255 * p->blue[i]);
	    write (p->gif, &r, 1);
	    write (p->gif, &g, 1);
	    write (p->gif, &b, 1);
	    }

	/* write extension block */

	if (p->wstype == 218) {		/* transparent GIF? */
	    c = 0x21;
	    write (p->gif, &c, 1);	/* extension block ID (fixed value) */
	    c = 0xf9;
	    write (p->gif, &c, 1);	/* graphic control label (fixed) */
	    c = 0x4;
	    write (p->gif, &c, 1);	/* block size (fixed) */
	    c = 0x0;			/* no disposal method, no user flag */
	    c |= 0x1;			/* transparent color on */
	    write (p->gif, &c, 1);
	    write_gif_word (0);		/* delay time */
	    c = 0x0;
	    write (p->gif, &c, 1);	/* transparent color index = 0 */
	    write (p->gif, &c, 1);	/* block terminator = 0 */
	    }
	
	/* write image block */

	write (p->gif, ",", 1);		/* image separator */

	write_gif_word (0);		/* image left offset */
	write_gif_word (0);		/* image top offset */
	write_gif_word (p->width);	/* image width */
	write_gif_word (p->height);	/* image height */

	c = 0x0;
	write (p->gif, &c, 1);		/* no interlace, sort, color table */

	InitCodeSize = max(BitsPerPixel, 2);
        gkscompress (InitCodeSize + 1, pix, size, beimage, &besize);

	c = InitCodeSize;
	write (p->gif, &c, 1);
	if (write (p->gif, beimage, besize) != besize) {
	    gks_fprintf (stderr, "GKS: can't write GIF file\n");
	    perror ("write");
	    }

        free (beimage);
        free (pix);
        }
    else
        gks_fprintf (stderr, "GKS: can't allocate temporary storage\n");

    write (p->gif, "\0", 1);	/* write out a zero-length packet (EOF) */
    write (p->gif, ";", 1);	/* terminator */

    XDestroyImage (image);
}


static
void open_rf (int stream)
{
    char *rf_name;
    struct stat buf;
    int fd;

    fd = stream > 100 ? stream - 100 : stream;

    if (isatty(fd) || fstat(fd, &buf) != 0)
        {
        rf_name = (char *) getenv ("GLI_RF");
        if (!rf_name)
            rf_name = "gli.rf";

        p->rf = open (rf_name, O_CREAT | O_TRUNC | O_WRONLY, 0644);
        if (p->rf < 0)
            gks_fprintf (stderr, "GKS: can't open rasterfile\n");
        }
    else
        p->rf = fd;
}


static
int compress_rle (byte *pix, int size, byte *beimage)
{
    register byte c, pc;
    register int besize, count, i, j;

    besize = 0;
    count = 0;
    for (i = 0; i < size; ++i)
        {
        c = *pix++;
        if (count > 0)
            {
            if (pc != c)
                {
                if (count == 1 && pc == 128)
                    {
                    beimage[besize++] = 128;
                    beimage[besize++] = 0;
                    count = 0;
                    }
                else if (count > 2 || pc == 128)
                    {
                    beimage[besize++] = 128;
                    beimage[besize++] = count - 1;
                    beimage[besize++] = pc;
                    count = 0;
                    }
                else
                    {
                    for (j = 0; j < count; ++j)
                        beimage[besize++] = pc;
                    count = 0;
                    }
                }
            }
        pc = c;
        ++count;
        if (count == 256)
            {
            beimage[besize++] = 128;
            beimage[besize++] = count - 1;
            beimage[besize++] = c;
            count = 0;
            }
        }
    if (count > 0)
        {
        if (count == 1 && c == 128)
            {
            beimage[besize++] = 128;
            beimage[besize++] = 0;
            }
        if (count > 2 || c == 128)
            {
            beimage[besize++] = 128;
            beimage[besize++] = count - 1;
            beimage[besize++] = c;
            }
        else
            {
            for (j = 0; j < count; ++j)
                beimage[besize++] = c;
            }
        }

    return besize;
}


static
void write_rf_long (long l)
{
    byte c;

    c = ((l >> 24) & 0xff);
    write (p->rf, &c, 1);
    c = ((l >> 16) & 0xff);
    write (p->rf, &c, 1);
    c = ((l >> 8) & 0xff);
    write (p->rf, &c, 1);
    c = (l & 0xff);
    write (p->rf, &c, 1);
}


#define RAS_MAGIC       0x59a66a95
#define RT_BYTE_ENCODED 2               /* Run-length compression of bytes */
#define RMT_EQUAL_RGB   1

static
void pixmap_to_rf (void)
{
    XImage *image;
    int linesize, size, besize, depth = 8;
    byte *pix, *ppix, *beimage;
    register int i, j, k, coli;
    byte rmap[255], gmap[255], bmap[255];
    unsigned long pixel;

    image = XGetImage (p->dpy, p->pixmap, 0, 0, p->width, p->height, AllPlanes,
        ZPixmap);

    linesize = p->width;
    if (linesize % 2) linesize++;

    size = linesize * p->height;

    pix = ppix = (byte *) malloc (sizeof(byte) * size);
    beimage = (byte *) malloc (sizeof(byte) * size * 3 / 2);   /* worst case */

    if (pix != NULL && beimage != NULL)
        {
        for (j = 0; j < p->height; j++) {
            for (i = 0; i < p->width; i++) {
                pixel = XGetPixel (image, i, j);
                coli = 0;
                for (k = 0; k < p->ncolors; k++) {
                    if (pixel == p->color[k].pixel) {
                        coli = k;
                        break;
                        }
                    }
                *ppix++ = coli;
                }
            if (linesize != p->width)
                *ppix++ = 0;
            }

        besize = compress_rle (pix, size, beimage);

        /* write the header */

        write_rf_long (RAS_MAGIC);
        write_rf_long (p->width);
        write_rf_long (p->height);
        write_rf_long (depth);
        write_rf_long (besize);
        write_rf_long (RT_BYTE_ENCODED);
        write_rf_long (RMT_EQUAL_RGB);
        write_rf_long (3 * p->ncolors);

        for (i = 0; i < p->ncolors; i++) {
            rmap[i] = (byte) (255 * p->red[i]);
            gmap[i] = (byte) (255 * p->green[i]);
            bmap[i] = (byte) (255 * p->blue[i]);
            }

        /* write the colormap */

        write (p->rf, rmap, p->ncolors);
        write (p->rf, gmap, p->ncolors);
        write (p->rf, bmap, p->ncolors);

        /* write the image */

        if (write (p->rf, beimage, besize) != besize) {
            gks_fprintf (stderr, "GKS: can't write Sun rle rasterfile\n");
	    perror ("write");
	    }

        free (beimage);
        free (pix);
        }
    else
        gks_fprintf (stderr, "GKS: can't allocate temporary storage\n");

    XDestroyImage (image);
}


#define UIL_HEADER "\
value\n\
    white : color ('white');\n\
    black : color ('black');\n\
    red : color ('red');\n\
    green : color ('green');\n\
    blue : color ('blue');\n\
    cyan : color ('cyan');\n\
    yellow : color ('yellow');\n\
    magenta : color ('magenta');\n\
\n\
value\n\
    color_map : color_table (\n\
        white = '`',\n\
        black = 'd',\n\
        red = 'r',\n\
        green = 'g',\n\
        blue = 'b',\n\
        cyan = 'c',\n\
        yellow = 'y',\n\
        magenta = 'm'\n\
    );\n"


static
void open_uil (void)
{
    char *uil_name;

    uil_name = (char *) getenv ("GLI_UIL");
    if (!uil_name)
        uil_name = "gli.uil";

    p->uil = fopen (uil_name, "w");
    if (p->uil)
        fprintf (p->uil, UIL_HEADER);
    else
        gks_fprintf (stderr, "GKS: can't open UIL file\n");
}


static
void pixmap_to_uil (void)
{
    static char *icon_name;
    XImage *image;
    int i, j, k, n, pix;
    unsigned long pixel;

    static char letter[] = {'`','d','r','g','b','c','y','m'};

    icon_name = (char *) getenv ("GLI_ICON");
    if (!icon_name)
        icon_name = "(unknown)";

    image = XGetImage (p->dpy, p->pixmap, 0, 0, p->width, p->height,
        AllPlanes, ZPixmap);
    n = p->ncolors < PRIVATE_COLORS ? p->ncolors : PRIVATE_COLORS;

    fprintf (p->uil, "\n%s : icon (color_table = color_map", icon_name);
    for (j=0; j<p->height; j++) {
        fprintf (p->uil, ",\n    '");
        for (i=0; i<p->width; i++) {
            pixel = XGetPixel (image, i, j);
            pix = 0;
            for (k=0; k<n; k++) {
                if (pixel == p->color[k].pixel) {
                    pix = k;
                    break;
                    }
                }
            fprintf (p->uil, "%c", letter[pix]);
            }
        fprintf (p->uil, "'");
        }
    if (fprintf (p->uil, "\n    );\n") < 0)
        gks_fprintf (stderr, "GKS: can't write UIL file\n");

    XDestroyImage (image);
}


static
void set_frame_header (int frame)
{
    char header[32];

    sprintf (header, "Frame #%d\n", frame);
    XStoreName (p->dpy, p->win, header);
}


static
void pixmap_loop (void)
{
    int this_frame = 0, inc = 1;
    XEvent event;
    Bool run = True, step = False;

    XSelectInput (p->dpy, p->win, ButtonPressMask);
    XSetClipMask (p->dpy, p->gc, None);
    XSynchronize (p->dpy, True);

    /* Be sure that the window is mapped */

    XMapWindow (p->dpy, p->win);

    for (; p->nframes > 0;) {
	if (run || step)
            {
	    XCopyArea (p->dpy, p->frame[this_frame], p->win, p->gc, 0, 0,
		p->width, p->height, 0, 0);
            this_frame += inc;
            if (this_frame == 0 || this_frame == p->nframes-1)
                inc = -inc;
            step = False;
            set_frame_header (this_frame);
	    }

        while (XPending (p->dpy)) {
            XNextEvent (p->dpy, &event);
            if (event.type == ButtonPress) {
                if (event.xbutton.button == Button1)
		    run = !run;
                else if (event.xbutton.button == Button2)
		    step = True;
		else
		    goto stop;
		}
            }
        }
stop:
    this_frame = p->nframes;
    while (this_frame--)
        XFreePixmap (p->dpy, p->frame[this_frame]);
    free (p->frame);

    p->pixmap = 0;
}



#define COPY_BODY(type) \
    register int i, j, ix, iy, nbytes; \
    register int *ilptr, *ipptr; \
    register byte *blptr, *bpptr, *packed_colia; \
    register type *elptr, *epptr, tmp, *tmpptr; \
    type pixel[MAX_COLORIND]; \
\
    for (i=0;  i<MAX_COLORIND;  i++) { \
	if (!p->mono_flag) { \
	    j = Color8Bit(i); \
	    pixel[i] = (type) p->color[j].pixel; \
	} \
	else \
	    pixel[i] = i; \
	} \
\
    if (p->packed_ca) \
        { \
        packed_colia = (byte *) colia; \
\
        if (dx != dimx || w != dx || h != dy || w != bytes_per_line) \
            { \
            elptr = epptr = ba; \
\
            for (j=0;  j<h;  j++, elptr+=bytes_per_line) { \
                iy = (dy * j) / h; \
                epptr = elptr; \
                blptr = packed_colia + (iy * dimx); \
                for (i=0;  i<w;  i++, epptr++) { \
                    ix = (dx * i) / w; \
                    bpptr = blptr + ix; \
                    *epptr = pixel[*bpptr]; \
                    } \
                } \
            } \
        else \
            { \
            nbytes = w*h; \
	    epptr = ba; \
	    bpptr = packed_colia; \
	    for (i=0;  i<nbytes;  i++) \
                *epptr++ = pixel[*bpptr++]; \
            } \
        } \
    else \
        { \
        if (dx != dimx || w != dx || h != dy || w != bytes_per_line) \
            { \
            elptr = epptr = ba; \
\
            for (j=0;  j<h;  j++, elptr+=bytes_per_line) { \
                iy = (dy * j) / h; \
                epptr = elptr; \
                ilptr = colia + (iy * dimx); \
                for (i=0;  i<w;  i++, epptr++) { \
                    ix = (dx * i) / w; \
                    ipptr = ilptr + ix; \
                    *epptr = pixel[*ipptr & 0x3ff]; \
                    } \
                } \
            } \
        else \
            { \
            nbytes = w*h; \
	    epptr = ba; \
	    ipptr = colia; \
	    for (i=0;  i<nbytes;  i++, ipptr++) \
                *epptr++ = pixel[*ipptr & 0x3ff]; \
            } \
        } \
\
    if (swapx) { \
        ix = 0; \
	w /= 2; \
	for (j=0;  j<h;  j++) { \
            for (i=0;  i<w;  i++) { \
                tmp = ba[i+ix]; \
                ba[i+ix] = ba[ix+w-i]; \
                ba[ix+w-i] = tmp; \
                } \
            } \
	    ix += bytes_per_line; \
        } \
\
    if (swapy) { \
	tmpptr = (type *) malloc (w*sizeof(type)); \
\
        elptr = ba; \
        epptr = ba + h*bytes_per_line; \
        h /=2; \
\
        for (j=0;  j<h;  j++) { \
	    epptr -= bytes_per_line; \
	    memcpy(tmpptr, elptr, w*sizeof(type)); \
	    memcpy(elptr, epptr, w*sizeof(type)); \
	    memcpy(epptr, tmpptr, w*sizeof(type)); \
	    elptr += bytes_per_line; \
	    } \
\
	free (tmpptr); \
	}


static
void copy32 (
    int dx, int dy, int dimx, int *colia,
    int w, int h, int bytes_per_line, int *ba, Bool swapx, Bool swapy)
{
    COPY_BODY(int)
}



static
void copy16 (
    int dx, int dy, int dimx, int *colia,
    int w, int h, int bytes_per_line, short int *ba, Bool swapx, Bool swapy)
{
    COPY_BODY(short int)
}



static
void copy8 (
    int dx, int dy, int dimx, int *colia,
    int w, int h, int bytes_per_line, byte *ba, Bool swapx, Bool swapy)
{
    COPY_BODY(byte)
}



static
void pixmap_to_bitmap (int w, int h, byte *ba)
{
    register byte *pix, *mbuffer, *bbuffer, mvalue, *first;
    register int i, j, k, graylevel, error, bit, row_size;
    int *lerr, *cerr, *terr, *error1, *error2;

    static unsigned char bit_flag[] = {1, 2, 4, 8, 16, 32, 64, 128};

    pix = ba;
    for (j=0; j < h; j++)
	for (i=0; i < w; i++) {
	    *pix = (byte) (p->gray[*pix] * (WHITE - BLACK));
	    pix++;
	    }

    /* Allocate space for error arrays */

    error1 = (int *) calloc (w+2, sizeof(int));
    error2 = (int *) calloc (w+2, sizeof(int));
    bbuffer = mbuffer = (byte *) calloc (w*h, sizeof(byte));
    
    cerr = &error1[1];
    lerr = &error2[1];
    for (error=0, i=0; i < w; ) {        /* The top border */
        mvalue = 0x00;
        for (j=0; (j < 8) && (i < w); j++, i++) {
            graylevel = (int)ba[i] + error;
            bit = graylevel > THRESH ? WHITE : BLACK;
            if (bit)
                mvalue |= (0x01 << j);
            error = graylevel - bit;
            lerr[i] = (THRESH - (int)bit) / 2;
            }
        *mbuffer++ = ~mvalue;
        }

    /*  Process the rest.            1 5 3
     *  Error distribution:          7 x
     */
    for (j=1; j < h; j++) {
        pix = &ba[j*w];
        first = mbuffer;
        for (i=0; i < w; ) {
            mvalue = 0x00;
            for (bit=0; (bit < 8) && (i < w); bit++, i++) {
                graylevel = (int)pix[i] + (int)(lerr[i-1] + (5 * lerr[i]) +
                    (3 * lerr[i+1]) + (7 * cerr[i-1])) / 16;
                if (graylevel > THRESH) {
                    mvalue |= (0x01 << bit);
                    cerr[i] = graylevel - WHITE;
            	    }
		else
                    cerr[i] = graylevel;
                }
            *mbuffer++ = ~mvalue;
            }
        cerr[-1] = INITERR((int)pix[-1], (*first & 0x01));

        /* Swap error buffers */
        terr = error1;
        error1 = error2;
        error2 = terr;
        cerr = &error1[1];
        lerr = &error2[1];
	}

    row_size = (w + 7) / 8;
    for (j=0; j < h; j++)
        for (i=0; i < w; i++) {
	    k = row_size*j + i/8;
	    if (*(bbuffer + k) & bit_flag[i % 8])
		*(ba + k) |= bit_flag[i % 8];
	    else
		*(ba + k) &= ~bit_flag[i % 8];
	    }

    free (bbuffer);
    free (error2);
    free (error1);
}



static
void cell_array (
    float xmin, float xmax, float ymin, float ymax,
    int dx, int dy, int dimx, int *colia)
{
    XImage *image = NULL;
 
    float x1, y1, x2, y2;
    int ix1, ix2, iy1, iy2;
    int	x, y, w, h, bytes_per_line, bitmap_pad;
    byte *ba;
    Bool swapx, swapy;

    WC_to_NDC(xmin, ymax, gksl->cntnr, x1, y1);
    seg_xform(&x1, &y1);
    NDC_to_DC(x1, y1, ix1, iy1);

    WC_to_NDC(xmax, ymin, gksl->cntnr, x2, y2);
    seg_xform(&x2, &y2);
    NDC_to_DC(x2, y2, ix2, iy2);

    w = abs(ix2-ix1) + 1;
    h = abs(iy2-iy1) + 1;
    x = min(ix1, ix2);
    y = min(iy1, iy2);

#ifdef XSHM
    if (p->shmimage != NULL && w == p->width && h == p->height)
	{
	image = p->shmimage;
	ba = (byte *) p->shminfo.shmaddr;
	bytes_per_line = p->shmimage->bytes_per_line;
	}
#endif
    bitmap_pad = (p->depth > 16 ? 32 : (p->depth > 8 ? 16 : 8));
    if (image == NULL)
	{
	ba = (byte *) malloc (w*h * bitmap_pad/8);
	bytes_per_line = w;
	}

    if (ba != NULL)
        {
        swapx = (xmin > xmax) ? True : False;
        swapy = (ymin < ymax) ? True : False;

        if (bitmap_pad == 32)
	    copy32 (dx, dy, dimx, colia, w, h, bytes_per_line, (int *) ba,
		swapx, swapy);
	else if (bitmap_pad == 16)
	    copy16 (dx, dy, dimx, colia, w, h, bytes_per_line, (short int *) ba,
		swapx, swapy);
	else
	    copy8 (dx, dy, dimx, colia, w, h, bytes_per_line, ba,
		swapx, swapy);

        if (p->mono_flag)
	    pixmap_to_bitmap (w, h, ba);

#ifdef XSHM
	if (image != NULL)
	    {
	    XShmPutImage (p->dpy, p->win, p->gc, image, 0, 0, x, y, w, h, True);
            XSync (p->dpy, False);
	    return;
	    }
#endif
	image = XCreateImage (p->dpy, p->vis, p->mono_flag ? 1 : p->depth,
	    p->mono_flag ? XYBitmap : ZPixmap, 0, (char *) ba, w, h,
	    bitmap_pad, 0);
	if (image)
	    {
            if (p->pixmap)
	        XPutImage (p->dpy, p->pixmap, p->gc, image, 0, 0, x, y, w, h);
	    if (!p->double_buf)
		XPutImage (p->dpy, p->win, p->gc, image, 0, 0, x, y, w, h);
            XSync (p->dpy, False);
	    /*
	     * Note: `XDestroyImage' frees both the image structure and the
	     * data pointed to by the image structure (ba)
	     */
#ifdef VMS
            free (ba);
#endif
            XDestroyImage (image);
            }
        else
	    gks_fprintf (stderr, "GKS: unable to create a %dx%d image\n", w, h);
        }
    else
	gks_fprintf (stderr, "GKS: can't allocate %dx%d data array\n", w, h);
}



static
void resize_window (void)
{
    int x, y, width, height, k;
    Arg arg[2];

    if (!p->uil) {
        width = nint((p->viewport[1]-p->viewport[0])/
	    p->mwidth*p->swidth*p->magnification);
        height = nint((p->viewport[3]-p->viewport[2])/
	    p->mheight*p->sheight*p->magnification);
        }
    else {
        width = nint((p->viewport[1]-p->viewport[0])*100);
        height = nint((p->viewport[3]-p->viewport[2])*100);
        }

    x = 4 + nint(p->viewport[0]/p->mwidth*p->swidth);
    y = p->sheight-height-4 - nint(p->viewport[2]/p->mheight*p->sheight);

    if (width != p->width || height != p->height || x != p->x || y != p->y)
        {
	p->x      = x;
	p->y      = y;
        p->width  = width;
	p->height = height;

	if (p->new_win)
            {
	    if (p->mapped)
		unmap_window ();

	    XMoveWindow (p->dpy, p->win, p->x, p->y);
    	    XResizeWindow (p->dpy, p->win, p->width, p->height);
	    }
        else {
            k = 0;
            XtSetArg (arg[k], "width", p->width); k++;
            XtSetArg (arg[k], "height", p->height); k++;
            XtSetValues (p->widget, arg, k);
            }

        if (p->pixmap) {
	    XFreePixmap (p->dpy, p->pixmap);
	    p->pixmap = XCreatePixmap (p->dpy, XRootWindowOfScreen(p->screen),
	        p->width, p->height, p->depth);
	    XFillRectangle (p->dpy, p->pixmap, p->clear, 0, 0,
	        p->width, p->height);
	    }
#ifdef XSHM
	free_shared_memory ();
	create_shared_memory ();
#endif
	}
}


static
void display_cursor (int x, int y)
{
    int xorg, yorg, width, height;
    int dx, dy, r, d;
    char str[16];

    if (x == Undefined && y == Undefined)
        return;

    switch (p->type)
	{
	case TypeNone:
	case TypeCross:
	    break;

	case TypeLocal:
	case TypeCrosshair:
	    XDrawLine (p->dpy, p->win, p->invert, 0, y, p->width, y);
	    XDrawLine (p->dpy, p->win, p->invert, x, 0, x, p->height);
	    break;

	case TypeRubberband:
	    XDrawLine (p->dpy, p->win, p->invert, p->px, p->py, x, y);
	    break;

	case TypeRectangle:
	    xorg = min(p->px, x);
	    yorg = min(p->py, y);
	    width = abs(p->px - x);
	    height = abs(p->py - y);
	    XDrawRectangle (p->dpy, p->win, p->invert, xorg, yorg,
		width, height);
	    break;

	case TypeDigital:
	    sprintf (str, "(%d %d)", x, y);
	    XDrawString (p->dpy, p->win, p->invert, p->px, p->py,
		str, strlen(str));
	    break;

	case TypeCircle:
	    dx = p->px - x;
	    dy = p->py - y;
	    r = (int) (sqrt((double)(dx*dx + dy*dy)) + 0.5);
	    d = 2*r;
	    if (r != 0)
		XDrawArc (p->dpy, p->win, p->invert, p->px-r, p->py-r,
		    d, d, 0, 360*64);
	    break;
	}
}


static
void get_pointer (int *n, float *x, float *y, int *state, int *term)
{
    Window focus, root_win, child_win;
    XEvent event;
    int np, revert, xold, yold;
    KeySym keysym;
    char str[10];
    int inc, xcur, ycur, xwin, ywin;
    static XComposeStatus status = { NULL, 0 };
    unsigned int mask, old_mask;
    
    XGetInputFocus (p->dpy, &focus, &revert);

    XDefineCursor (p->dpy, p->win, p->cursor);
    XRaiseWindow (p->dpy, p->win);

    select_async_input (False);

    if (p->new_win) {
	while (XCheckTypedWindowEvent (p->dpy, p->win, ConfigureNotify, &event))
	    ;
	}

    np = 0;

    xold = p->px;
    yold = p->py;
    old_mask = 0;
    if (xold == Undefined || yold == Undefined)
	{
	XQueryPointer (p->dpy, p->win, &root_win, &child_win, &xcur, &ycur,
	    &xold, &yold, &mask);
	}
    display_cursor (xold, yold);

    *term = 0;

    do
	{
	*state = Undefined;
        inc = 1;

	if (p->wstype == 212)
#ifndef _WIN32
#ifdef __osf__
	    usleep (200000); /* 200 ms */
#else
	    sleep (1);
#endif
#endif

	do
	    {
	    if (p->wstype == 212)
		{
		while (True)
		    {
		    old_mask = mask;
		    XQueryPointer (p->dpy, p->win, &root_win, &child_win,
			&xcur, &ycur, &xwin, &ywin, &mask);
		    if (xwin != xold || ywin != yold || mask != old_mask)
			break;
		    }

		switch (mask)
		    {
		    case Button1Mask:
		    case Button2Mask:
			event.xany.type = ButtonPress;
			event.xbutton.button = (mask == Button1Mask) ?
			    Button1 : Button2;
			event.xbutton.x = xwin;
			event.xbutton.y = ywin;
			break;
		    default:
			event.xany.type = MotionNotify;
			event.xmotion.x = xwin;
			event.xmotion.y = ywin;
		    }
		}
	    else
		 XWindowEvent (p->dpy, p->win, ButtonPressMask |
			PointerMotionMask | KeyPressMask | KeyReleaseMask |
			StructureNotifyMask | ExposureMask, &event);

	    switch (event.xany.type)
	      {
	      case Expose:
		async_expose_event (p);
		xold = yold = Undefined;
                break;

	      case ButtonPress:
		xcur = event.xbutton.x;
		ycur = event.xbutton.y;

		DC_to_NDC(xcur, ycur, *x++, *y++);

		if (event.xbutton.button == Button1) {
		    np++;
		    *state = GOK;
		    }
		else
		    *state = GNONE;
	        break;

	      case MotionNotify:
		display_cursor (xold, yold);

		xcur = event.xmotion.x;
		ycur = event.xmotion.y;
		display_cursor (xcur, ycur);

		xold = xcur;
		yold = ycur;
	        break;

	      case KeyPress:
                xcur = xold;
                ycur = yold;
		display_cursor (xold, yold);

		XLookupString ((XKeyEvent *) &event, str, 9, &keysym, &status);

		switch (keysym)
		    {
		    case XK_Shift_L:
		    case XK_Shift_R:	inc = 10; break;
		    case XK_Left:	xcur -= inc; break;
		    case XK_Right:	xcur += inc; break;
		    case XK_Up:		ycur -= inc; break;
		    case XK_Down:	ycur += inc; break;
		    case XK_Control_L:
		    case XK_Control_R:
                    case XK_Caps_Lock:
                    case XK_Shift_Lock:
                    case XK_Meta_L:
                    case XK_Meta_R:
                    case XK_Alt_L:
                    case XK_Alt_R:
                    case XK_Multi_key:  break;

		    default:
                        if (*str != CTRL_C && *str != CTRL_D && *str != CTRL_Z)
                            {
                            DC_to_NDC(xcur, ycur, *x++, *y++);
                            np++;
                            *state = GOK;
                            }
                        else
                            *state = GNONE;

                        *term = *str;
			break;
		    }

		XWarpPointer (p->dpy, None, p->win, 0, 0, 0, 0, xcur, ycur);

		display_cursor (xcur, ycur);

		xold = xcur;
		yold = ycur;
	        break;

	      case KeyRelease:
		XLookupString ((XKeyEvent *) &event, str, 9, &keysym, &status);

		if (keysym == XK_Shift_L || keysym == XK_Shift_R) inc = 1;
	        break;

	      case ConfigureNotify:
                p->empty = False;
                configure_event ((XConfigureEvent *) &event);
		if (p->empty)
                    xold = yold = Undefined;

                *state = Undefined;
	        break;
	      }
	    }
	while (*state < 0);
	}
    while (np < *n && *state != GNONE);

    display_cursor (xold, yold);

    select_async_input (True);
    
    XSetInputFocus (p->dpy, focus, revert, CurrentTime);
    XRaiseWindow (p->dpy, p->win);

    XUndefineCursor (p->dpy, p->win);
    XSync (p->dpy, False);

    *n = np;
    if (*n > 1)
	*state = GOK;

    p->px = xcur;
    p->py = ycur;
}


static
int lookup_string (char *str)
{
    register int i;
    char s1[3], s2[3];

    s1[0] = str[0]; s1[1] = str[1]; s1[2] = '\0';
    s2[0] = str[1]; s2[1] = str[0]; s2[2] = '\0';

    for (i=0; i<n_key; i++) {
        if (strcmp(s1, key_bindings[i].seq) == 0 ||
            strcmp(s2, key_bindings[i].seq) == 0 ||
            strcmp(s1, key_bindings[i].alt_seq) == 0 ||
            strcmp(s2, key_bindings[i].alt_seq) == 0)
            return (key_bindings[i].ch);
        }

    return (0);
}


static
int dispatch_character (XKeyEvent *event, char *text)
{
    KeySym keysym;
    static char str[10], seq[3];
    static XComposeStatus compose_status = { NULL, 0 };
    static char recall_buffer[256] = "";

    XDrawString (p->dpy, p->win, p->invert, p->px, p->py, text, strlen(text));

    if (event) {
        XLookupString (event, str, 9, &keysym, &compose_status);

#ifndef VMS
        if (keysym == XK_Multi_key)
            compose_status.chars_matched = 1;

        else if (keysym >= XK_space && keysym <= XK_asciitilde)
            {
            switch (compose_status.chars_matched)
                {
                case 1 :
                    seq[0] = (char) keysym;
                    compose_status.chars_matched++;
                    keysym = 0;
                    break;

                case 2 :
                    seq[1] = (char) keysym;
                    seq[2] = '\0';
                    keysym = lookup_string (seq);

                default :
                    compose_status.chars_matched = 0;
                    break;
                }
            }
#endif
        if ((keysym >= XK_space && keysym <= XK_asciitilde) ||
            (keysym >= XK_nobreakspace && keysym <= XK_ydiaeresis)) {
            str[0] = (char) keysym;
            str[1] = '\0';
            strcat (text, str);
            }

        else if (keysym == XK_Delete) {
            if (*text)
                text[strlen(text)-1] = '\0';
            }
        else if (keysym == XK_Up)
            strcpy (text, recall_buffer);

        XDrawString (p->dpy, p->win, p->invert, p->px, p->py, text,
            strlen(text));

        if (keysym == XK_Return)
            strcpy (recall_buffer, text);
        }
    else
        keysym = XK_space;

    return (keysym);
}


static
void draw_text_box (void)

{
    int x, y;
    
    x = p->px - 5;
    y = p->py + 5;
    
    XDrawLine (p->dpy, p->win, p->invert, x, y, x, y-20);
    XDrawLine (p->dpy, p->win, p->invert, x, y, x+100, y);
}


static
void get_string (int *n, char *chars, int *state)
{
    XEvent event;
    Window focus;
    int done, revert;
    char text[256];

    XGetInputFocus (p->dpy, &focus, &revert);
    
    XDefineCursor (p->dpy, p->win, p->textcursor);
    XRaiseWindow (p->dpy, p->win);

    select_async_input (False);

    draw_text_box ();
    
    strcpy (text, "");
    done = FALSE;

    do
	{
	XWindowEvent (p->dpy, p->win, ButtonPressMask | PointerMotionMask |
	  KeyPressMask | ExposureMask, &event);

	switch (event.type)
	  {
	  case Expose:
	    async_expose_event (p);
	    break;

	  case KeyPress :
	    if (dispatch_character ((XKeyEvent *) &event, text) == XK_Return) {
	      strcpy (chars, text);
	      done = TRUE;
	      }
	    break;

	  case ButtonPress :
	    if (event.xbutton.button == Button1) {
	      strcpy (chars, text);
	      done = TRUE;
	      }
	    break;
	  }
	}

    while (!done && event.type != ButtonPress);

    *n = strlen(text);
    dispatch_character (NULL, text);

    draw_text_box ();
    
    select_async_input (True);
    
    XSetInputFocus (p->dpy, focus, revert, CurrentTime);
    XRaiseWindow (p->dpy, p->win);

    XUndefineCursor (p->dpy, p->win);
    XSync (p->dpy, False);

    if (done)
      *state = GOK;
    else
      *state = GNONE;
}

#endif


#if !defined(NO_X11)

void STDCALL GKDXW (
    int *fctid, int *dx, int *dy, int *dimx, int *ia, int *lr1, float *r1,
    int *lr2, float *r2, int *lc, CHARARG(chars), ws_state_list **ptr)
{
    static int win = 0;

    p = *ptr;

    switch (function_id = *fctid) {

      case 2:
/*
 *  Open workstation
 */
	p = (ws_state_list *) malloc (sizeof(ws_state_list));

        p->wkid = ia[0];

        p->packed_ca = getenv("GLI_GKS_PACKED_CELL_ARRAY") ? True : False;
        p->async_io = getenv("GLI_GKS_ASYNC_IO") ? True : False;
        p->double_buf = getenv("GLI_GKS_DOUBLE_BUF") ? True : False;
        p->shape = getenv("GLI_GKS_CONVEX_SHAPE") ? Convex : Complex;
	p->widget = (Widget) NULL;
	p->conid = ia[1];
	p->wstype = ia[2];
        p->gif = -1;
        p->rf = -1;
	p->uil = NULL;
        p->frame = NULL;

	switch (p->wstype) {

	    case 230:
	    case 231:
	    case 232:
	    case 233:
		p->wstype -= 20;
	    case 210:
	    case 211:
	    case 212:
	    case 213:
		if ((unsigned) ia[1] >= 100 + 100) {
		    p->widget = (Widget) *ptr;
                    p->async_io = False;
                    }
		break;

	    case 214:
	    case 215:
	    case 218:
                p->async_io = False;
                break;

	    case 216:
                open_uil ();
                p->async_io = False;
		break;

	    case 217:
                p->frame = (Pixmap *) malloc (MAX_PIXMAP * sizeof(Pixmap));
                p->nframes = 0;
                p->async_io = False;
		break;
	    }

	p->ccolor = Undefined;
	p->ncolors = 0;
        p->scolors = SHARED_COLORS;

	p->error = NULL;

	if (open_display () == NULL) {
            free (p);
            ia[0] = ia[1] = 0;
            return;
            }

	max_points = MAX_POINTS;
	points = (XPoint *) malloc (max_points * sizeof(XPoint));

	set_colors ();
        configure_colors ();
        allocate_colors (True);

        if (p->wstype == 215 || p->wstype == 218)
            open_gif (p->conid);
        else if (p->wstype == 214)
            open_rf (p->conid);

	create_window (win);
	set_WM_hints ();
	create_GC ();
	create_pixmap ();
#ifdef XSHM
	create_shared_memory ();
#endif
	initialize_arrays ();
	create_cursor ();

	p->state = GINACT;
	p->mapped = False;

        p->ltype = GSOLID;
        p->lwidth = 0;

	/* Setup default device transformation */

	p->window[0] = 0; p->window[1] = 1;
	p->window[2] = 0; p->window[3] = 1;

	p->viewport[0] = 0; p->viewport[1] = p->mwidth*p->width/p->swidth;
	p->viewport[2] = 0; p->viewport[3] = p->mheight*p->height/p->sheight;

	setup_xform (p->window, p->viewport);

	p->type = TypeLocal;
	p->px = Undefined;
	p->py = Undefined;

	/* Return state list address, screen width and height */

	if (sizeof(char *) > sizeof(int)) {
	    long *la = (long *)ia;
	    *la = (long)p;
	    }
	else
	    *ia = (long)p;

	r1[0] = p->mwidth; r2[0] = p->mheight;
	r1[1] = p->swidth; r2[1] = p->sheight;

	if (!strncmp(ServerVendor(p->dpy),
	    "X11/NeWS - Sun Microsystems Inc.", 32))
	    p->server = SunServer;
        else if (!strncmp(ServerVendor(p->dpy), "Silicon Graphics", 16))
	    p->server = SGIServer;
        else if (!strncmp(ServerVendor(p->dpy),
	    "DECWINDOWS Digital Equipment Corporation", 40))
	    p->server = DECServer;
        else if (!strncmp(ServerVendor(p->dpy),
	    "International Business Machines", 31))
	    p->server = IBMServer;
        else if (!strncmp(ServerVendor(p->dpy), "TEKtronix, Inc.", 15))
	    p->server = TEKServer;
	else
	    p->server = AnyServer;

	p->xshm = ((char *) getenv ("GLI_GKS_XSHM") != NULL) ? True : False;
	p->gr_bc = ((char *) getenv ("GR_BC") != NULL) ? True : False;
	
#if defined(hpux) && !defined(NAGware)
        line_routine_a = (void (*)()) line_routine;
        fill_routine_a = (void (*)()) fill_routine;
        text_routine_a = (void (*)()) text_routine;
#endif

	if (p->new_win)
	    win++;
	else
	    XClearWindow (p->dpy, p->win);

        /* Get GKS state list address */

        gksl = (gks_state_list *)(ia+4);

        init_norm_xform ();
        set_clipping (True);

	break;

      case 3:
/*
 *  Close workstation
 */
        if (p->async)
            async_io (False);

	if (p->gif >= 0)
            pixmap_to_gif ();
	else if (p->rf >= 0)
            pixmap_to_rf ();

        else if (p->frame)
            if (p->nframes > 1)
		pixmap_loop ();

	if (p->new_win)
	    unmap_window ();
#ifdef XSHM
	free_shared_memory ();
#endif
        if (p->pixmap)
	    XFreePixmap (p->dpy, p->pixmap);

	free_GC ();
	free_colors ();

        if (p->new_win) {
	    XDestroyWindow (p->dpy, p->win);
	    win--;
	    
	    XFreePixmap (p->dpy, p->icon_pixmap);
	    }
	if (p->new_dpy)
	    XCloseDisplay (p->dpy);

        if (p->gif >= 0)
            close (p->gif);
        else if (p->rf >= 0)
            close (p->rf);
        else if (p->uil)
            fclose (p->uil);

	free (points);
	free (p);
	break;

      case 4:
/*
 *  Activate workstation
 */
	p->state = GACTIV;
	break;

      case 5:
/*
 *  Deactivate workstation
 */
	p->state = GINACT;
	break;

      case 6:
/*
 *  Clear workstation
 */
        if (p->async)
            async_io (False);

	if (p->uil)
            pixmap_to_uil ();

        else if (p->frame) {
            if (p->nframes != MAX_PIXMAP) {
		p->frame[p->nframes++] = p->pixmap;
                set_frame_header (p->nframes);
		create_pixmap ();
		}
            }

	if (p->pixmap)
	    XFillRectangle (p->dpy, p->pixmap, p->clear, 0, 0, p->width,
                p->height);
	if (!p->double_buf)
	    XClearWindow (p->dpy, p->win);

	p->empty = True;

        if (!p->frame)
            XSync (p->dpy, False);
	break;

      case 8:
/*
 *  Update workstation
 */
        if (p->async)
            async_io (False);

        if (p->double_buf && ia[1] == GPERFO)
	    async_expose_event (p);

	update ();

	if (p->state == GACTIV) {
	    if (p->error) {
		gks_fprintf (stderr, "%s\n", p->error);
		p->error = NULL;
		}
	    }

        if (p->async_io)
            async_io (True);
        break;

      case 10:
/*
 *  Message
 */
        if (p->async)
            async_io (False);
	if (!p->mapped)
	    map_window ();

	if (p->state == GACTIV)
#ifdef VMS
	    message (chars->dsc$w_length, chars->dsc$a_pointer);
#else
#ifdef cray
	    message (_fcdlen(chars), _fcdtocp(chars));
#else
	    message (strlen(chars), chars);
#endif /* cray */
#endif /* VMS */
	break;

      case 12:
/*
 *  Polyline
 */
        if (p->async)
            async_io (False);
	if (!p->mapped)
	    map_window ();

	if (p->state == GACTIV)
	    polyline (ia, r1, r2);
	break;

      case 13:
/*
 *  Polymarker
 */
        if (p->async)
            async_io (False);
	if (!p->mapped)
	    map_window ();

	if (p->state == GACTIV)
	    polymarker (ia, r1, r2);
	break;

      case 14:
/*
 *  Text
 */
        if (p->async)
            async_io (False);
	if (!p->mapped)
	    map_window ();

	if (p->state == GACTIV)
#ifdef VMS
	    text (r1, r2, chars->dsc$w_length, chars);
#else
#ifdef cray
	    text (r1, r2, _fcdlen(chars), chars);
#else
	    text (r1, r2, strlen(chars), chars);
#endif /* cray */
#endif /* VMS */
	break;

      case 15:
/*
 *  Fill Area
 */
        if (p->async)
            async_io (False);
	if (!p->mapped)
	    map_window ();

	if (p->state == GACTIV)
	    fill_area (ia, r1, r2);
	break;

      case 16:
/*
 *  Cell Array
 */
        if (p->async)
            async_io (False);
	if (!p->mapped)
	    map_window ();

	if (p->state == GACTIV)
            {
	    if (p->min_colors)
                allocate_colors (False);
            cell_array (r1[0], r1[1], r2[0], r2[1], *dx, *dy, *dimx, ia);
            }
	break;

      case 48:
/*
 *  Set color representation
 */
        if (p->async)
            async_io (False);

	if (ia[1] >= 0 && ia[1] < MAX_COLORS) {
	    set_color_repr (ia[1], r1[0], r1[1], r1[2]);
            free_tile_patterns (ia[1]);
            }
	break;

      case 49:
/*
 *  Set window
 */
        if (p->async)
            async_io (False);

	setup_norm_xform (*ia, gksl->window[*ia], gksl->viewport[*ia]);
        set_clipping (True);
	break;

      case 50:
/*
 *  Set viewport
 */
        if (p->async)
            async_io (False);

	setup_norm_xform (*ia, gksl->window[*ia], gksl->viewport[*ia]);
        set_clipping (True);
	break;

      case 52:
/*
 *  Select normalization transformation
 */
        if (p->async)
            async_io (False);

        set_clipping (True);
	break;

      case 53:
/*
 *  Set clipping indicator
 */
        if (p->async)
            async_io (False);

        set_clipping (True);
	break;

      case 54:
/*
 *  Set workstation window
 */
        if (p->async)
            async_io (False);

	p->window[0] = r1[0];
	p->window[1] = r1[1];
	p->window[2] = r2[0];
	p->window[3] = r2[1];

	setup_xform (p->window, p->viewport);
        set_clipping (True);
	break;

      case 55:
/*
 *  Set workstation viewport
 */
	{
	float max_width, max_height;

        if (p->async)
            async_io (False);

	p->viewport[0] = r1[0];
	p->viewport[1] = r1[1];
	p->viewport[2] = r2[0];
	p->viewport[3] = r2[1];

        if (p->gif >= 0 || p->rf >= 0) {
	    max_width  = p->mwidth  / p->swidth  * 1280;
	    max_height = p->mheight / p->sheight * 1024;
        } else if (p->widget == NULL) {
	    max_width  = p->mwidth  * 0.925;
	    max_height = p->mheight * 0.925;
        } else {
	    max_width  = 10.0;
	    max_height = 10.0;
        }
	if (!p->uil)
	    GKFVP (p->viewport, &max_width, &max_height);

	resize_window ();
	set_WM_hints ();

	setup_xform (p->window, p->viewport);
        set_clipping (True);
	break;
	}

      case 69:
/*
 *  Initialize locator
 */
	p->type = (pe_type) ia[3];
        NDC_to_DC(r1[0], r2[0], p->px, p->py);
	break;

      case 81:
/*
 *  Request locator
 */
	{
	int n;

        if (p->async)
            async_io (False);
	if (!p->mapped)
	    map_window ();

	n = 1;
	get_pointer (&n, r1, r2, &ia[0], &ia[3]);
	break;
	}

      case 82:
/*
 *  Request stroke
 */
        if (p->async)
            async_io (False);
	if (!p->mapped)
	    map_window ();

	get_pointer (&ia[2], r1, r2, &ia[0], &ia[3]);
	break;

      case 86:
/*
 *  Request string
 */
        if (p->async)
            async_io (False);
	if (!p->mapped)
	    map_window ();
#ifdef VMS
	get_string (&ia[1], chars->dsc$a_pointer, &ia[0]);
#else
#ifdef cray
	get_string (&ia[1], _fcdtocp(chars), &ia[0]);
#else
	get_string (&ia[1], chars, &ia[0]);
#endif /* cray */
#endif /* VMS */
	break;
      }
}

#else

void STDCALL GKDXW (
    int *fctid, int *dx, int *dy, int *dimx, int *ia, int *lr1, float *r1,
    int *lr2, float *r2, int *lc, CHARARG(chars), void **ptr)
{
    if (*fctid == 2)
        {
        gks_fprintf (stderr,
	    "GKS: X11 not supported for this type of machine\n");
        ia[0] = 0;
        }
}

#endif
