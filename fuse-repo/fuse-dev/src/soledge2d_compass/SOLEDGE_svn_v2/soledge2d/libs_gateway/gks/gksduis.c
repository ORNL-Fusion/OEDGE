/*
 * Copyright @ 1984 - 1993   Josef Heinen
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
 *      GLI GKS V4.5
 *
 * ABSTRACT:
 *
 *      This module contains a logical device driver for the VMS UIS Software
 *
 * AUTHOR:
 *
 *      Josef Heinen
 *
 * VERSION:
 *
 *      V1.0
 *
 */


#include <stdio.h>
#include <math.h>

#ifdef VMS
#include <descrip.h>
#ifdef UIS
#include <uisentry.h>
#include <uisusrdef.h>
#endif /* UIS */
#endif /* VMS */

#ifdef cray
#include <fortran.h>
#endif

#include "gksdefs.h"

#ifdef UIS

#pragma builtins

#define WindowName "GLIgks V4.5"
#define DrawBorder 0
#define Undefined -99

#define BOOL unsigned

#define PRIVATE_COLORS  8
#define SHARED_COLORS   72
#define MAX_COLORS	80
#define MAX_POINTS	512

#define WHITE 255
#define THRESH 127
#define BLACK 0
#define INITERR(X,Y)    (X - (Y ? WHITE : BLACK) + (THRESH - X)/2)

typedef unsigned char byte;


static $DESCRIPTOR(devnam, "SYS$WORKSTATION");
static $DESCRIPTOR(title, WindowName);
static $DESCRIPTOR(font, "UIS$FILL_PATTERNS");

static int state;
static float window[4], viewport[4];
static float width, height, a, b, c, d;

static unsigned int wd_id, vd_id, vcm_id, map_size;
static float retwidth, retheight, retresolx, retresoly;
static unsigned int retpwidth, retpheight;

static int pwidth, pheight;
static unsigned int atb;
static float linewidth;
static float xpoint[MAX_POINTS+1], ypoint[MAX_POINTS+1];
static int npoints;
static float gray[MAX_COLORS];


static set_colors ()
{
    int i, j;
    unsigned int old_atb, new_atb;
    float r, g, b;

    static float rint[] = {1, 0, 1, 0, 0, 0, 1, 1};
    static float gint[] = {1, 0, 0, 1, 0, 1, 1, 0};
    static float bint[] = {1, 0, 0, 0, 1, 1, 0, 1};

    for (i=0; i<PRIVATE_COLORS; i++)
	gray[i] = 0.3*rint[i] + 0.59*gint[i] + 0.11*bint[i];

    old_atb = 0;
    new_atb = (map_size <= 8) ? map_size : 8;
    uis$set_colors (&vd_id, &old_atb, &new_atb, rint, gint, bint);

    for (i=0; i<SHARED_COLORS; i++)
	{
        j = i;
	GQRGB (&j, &r, &g, &b);
	gray[i+PRIVATE_COLORS] = 0.3*r + 0.59*g + 0.11*b;

	if (map_size > PRIVATE_COLORS) {
	    uis$set_color (&vd_id, &new_atb, &r, &g, &b);
	    new_atb++;
	    }
	}
}


static uis_setup ()
{
    unsigned int old_atb, new_atb;
    float x1, y1, x2, y2;
    int i, attributes[3];
 
    uis$get_display_size (&devnam, &retwidth, &retheight,
	&retresolx, &retresoly, &retpwidth, &retpheight);

    uis$get_hw_color_info (&devnam, NULL, &map_size, NULL, NULL, NULL, NULL,
	NULL, NULL, NULL, NULL);

    if (map_size >= 256)
	map_size = PRIVATE_COLORS + SHARED_COLORS;
    else if (map_size >= 16)
	map_size = PRIVATE_COLORS;
    else
	map_size = 2;

    vcm_id = uis$create_color_map (&map_size, NULL, NULL);

    x1 = 0.0;
    y1 = 0.0;
    x2 = pwidth = 485;
    y2 = pheight = 485;
    width = pwidth * retwidth / retpwidth;
    height = pheight * retheight / retpheight;

    vd_id = uis$create_display (&x1, &y1, &x2, &y2, &width, &height, &vcm_id);

    for (i=0; i<map_size; i++) {
	old_atb = 0;
	new_atb = i + 1;
	uis$set_writing_index (&vd_id, &old_atb, &new_atb, &i);
	uis$set_font (&vd_id, &new_atb, &new_atb, &font);
	}

    set_colors ();

    attributes[0] = WDPL$C_PLACEMENT;
    attributes[1] = WDPL$M_TOP | WDPL$M_LEFT;
    attributes[2] = WDPL$C_END_OF_LIST;

    wd_id = uis$create_window (&vd_id, &devnam, &title, &x1, &y1, &x2, &y2,
	&width, &height, attributes);

    uis$push_viewport (&wd_id);

    atb = 0;
    linewidth = 0.0;
}


static set_up_xform (window, viewport)
    float *window, *viewport;
{
    a = (pwidth-1)/(window[1]-window[0]);
    b = window[0]*a + 0.5;
    c = (pheight-1)/(window[3]-window[2]);
    d = window[2]*c + 0.5;
}


static resize_window ()
{
    float new_abs_x, new_abs_y;
    float x1, y1, x2, y2;

    new_abs_x = viewport[0] * 100;
    width = (viewport[1] - viewport[0]) * 100;
    new_abs_y = viewport[2] * 100;
    height = (viewport[3] - viewport[2]) * 100;

    x1 = 0.0;
    x2 = width * retpwidth / retwidth;
    pwidth = (int) x2;
    y1 = 0.0;
    y2 = height * retpheight / retheight;
    pheight = (int) y2;

    uis$resize_window (&vd_id, &wd_id, &new_abs_x, &new_abs_y, &width, &height,
        &x1, &y1, &x2, &y2);

    uis$push_viewport (&wd_id);
}


static change_gc (color, width)
    int color;
    float width;
{
    if (++color != atb || width != linewidth) {
	end_line ();
	atb = (color <= map_size) ? color : map_size;
	linewidth = width;
	uis$set_line_width (&vd_id, &atb, &atb, &linewidth, NULL);
	}
}


static start_line (x, y)
    float *x, *y;
{
    if (npoints > 1)
        uis$plot_array (&vd_id, &atb, &npoints, xpoint, ypoint);

    xpoint[0] = a * *x + b;
    ypoint[0] = c * *y + d;
    npoints = 1;
}


static draw_line (x, y)
    float *x, *y;
{
    xpoint[npoints] = a * *x + b;
    ypoint[npoints] = c * *y + d;

    if (npoints == MAX_POINTS) {
        uis$plot_array (&vd_id, &atb, &npoints, xpoint, ypoint);

        xpoint[0] = xpoint[npoints];
        ypoint[0] = ypoint[npoints];
        npoints = 0;
        }
    npoints++;
}


static end_line ()
{
    if (npoints > 1)
        uis$plot_array (&vd_id, &atb, &npoints, xpoint, ypoint);
    npoints = 0;
}


static uis_move (x, y)
    float *x, *y;
{
    GMOVE (x, y, start_line);
}


static uis_draw (x, y, linetype)
    float *x, *y;
    int *linetype;
{
    GDASH (x, y, start_line, draw_line);
}


static polyline (n, px, py)
    int *n;
    float *px, *py;
{
    int errind, tnr, color;
    float linewidth;
    int linetype;

    /* Inquire current normalization transformation number, linewidth
     * scale factor, line type, and polyline color index */

    GQCNTN (&errind, &tnr);
    GSDT (window, viewport);

    GQLWSC (&errind, &linewidth);
    GQLN (&errind, &linetype);

    GQPLCI (&errind, &color);
    change_gc (color, linewidth);

    GPOLIN (n, px, py, &linetype, &tnr, uis_move, uis_draw);
}


static line_routine (n, px, py, linetype, tnr)
    int *n;
    float *px, *py;
    int *linetype, *tnr;
{
    GPOLIN (n, px, py, linetype, tnr, start_line, draw_line);
}


static polymarker (n, px, py)
    int *n;
    float *px, *py;
{
    int errind, color;

    GSDT (window, viewport);

    GQPMCI (&errind, &color);
    change_gc (color, 1.0);

    GPOLMK (n, px, py, line_routine);
}


static uis_draw_lines (n, px, py, tnr)
    int *n;
    float *px, *py;
    int *tnr;
{
    int i;
    float x, y;

    npoints = (*n < MAX_POINTS) ? *n : MAX_POINTS;

    for (i=0; i<npoints; i++) {
	x = px[i];
	y = py[i];
	GNT (&x, &y, tnr);

	xpoint[i] = a * x + b;
	ypoint[i] = c * y + d;
	}

    uis$plot_array (&vd_id, &atb, &npoints, xpoint, ypoint);
    npoints = 0;
}


static fillarea (n, px, py)
    int *n;
    float *px, *py;
{
    int errind, tnr, style, color;
    float resolution;
    unsigned int pattern;
    int linetype;

    GQCNTN (&errind, &tnr);
    GSDT (window, viewport);

    GQFAIS (&errind, &style);
    GQFACI (&errind, &color);
    change_gc (color, 1.0);

    if (style == GSOLID)
	{
	if (color < map_size) {
	    pattern = PATT$C_FOREGROUND;
	    uis$set_fill_pattern (&vd_id, &atb, &atb, &pattern);

	    uis_draw_lines (n, px, py, &tnr);

	    uis$set_fill_pattern (&vd_id, &atb, &atb, NULL);
	    }
	else if (color >= PRIVATE_COLORS) {
	    pattern = (int)(15*(1-gray[color-PRIVATE_COLORS])) +
		PATT$C_GREY1_16;
	    uis$set_fill_pattern (&vd_id, &atb, &atb, &pattern);

	    uis_draw_lines (n, px, py, &tnr);

	    uis$set_fill_pattern (&vd_id, &atb, &atb, NULL);
 	    }
	else {
	    linetype = DrawBorder;
	    line_routine (n, px, py, &linetype, &tnr);
	    }
	}
    else
	{
	resolution = retpheight;
	GFILLA (n, px, py, line_routine, &resolution);
	}
}


static fill_routine (n, px, py, tnr)
    int *n;
    float *px, *py;
    int *tnr;
{
    unsigned int pattern;

    pattern = PATT$C_FOREGROUND;
    uis$set_fill_pattern (&vd_id, &atb, &atb, &pattern);

    uis_draw_lines (n, px, py, tnr);

    uis$set_fill_pattern (&vd_id, &atb, &atb, NULL);
}


static text (px, py, chars)
    float *px, *py;
    struct dsc$descriptor *chars;
{
    int errind, color, nchars;

    GSDT (window, viewport);

    GQTXCI (&errind, &color);
    change_gc (color, 1.0);

    nchars = strlen(chars);
    GSIMTX (px, py, nchars, chars, line_routine, fill_routine);
}


static pack_data (dx, dy, dimx, colia, w, h, ba, swapx, swapy)
    int dx, dy, dimx;
    int *colia;
    int w, h;
    byte *ba;
    BOOL swapx, swapy;
{
    int i, j, ix, iy, pix;
    int *ilptr, *ipptr;
    byte *elptr, *epptr, tmp;

    elptr = epptr = ba;

    for (j=0;  j<h;  j++, elptr+=w) {
        iy = (dy * j) / h;
        epptr = elptr;
        ilptr = colia + (iy * dimx);
        for (i=0;  i<w;  i++, epptr++) {
            ix = (dx * i) / w;
            ipptr = ilptr + ix;
            if (*ipptr >= PRIVATE_COLORS + SHARED_COLORS)
                *epptr = PRIVATE_COLORS + SHARED_COLORS - 1;
            else if (*ipptr < 0)
                *epptr = 0;
            else
                *epptr = *ipptr;
            }
        }

    if (swapx) {
        for (j=0;  j<h;  j++) {
	    ix = j*w;
            for (i=0;  i<w/2;  i++) {
                tmp = ba[i+ix];
                ba[i+ix] = ba[ix+w-i];
                ba[ix+w-i] = tmp;
                }
            }
        }

    if (swapy) {
        iy = (h-1)*w;
        for (j=0;  j<h/2;  j++) {
	    ix = j*w;
            for (i=0;  i<w;  i++) {
                tmp = ba[i+ix];
                ba[i+ix] = ba[i+iy];
                ba[i+iy] = tmp;
                }
            iy -= w;
            }
        }
}


static pixmap_to_bitmap (w, h, ba)
    int w, h;
    byte *ba;
{
    register byte *pix, *mbuffer, *bbuffer, mvalue, *first;
    register int i, j, graylevel, error, bit, row_size, bdone;
    int *lerr, *cerr, *terr, *error1, *error2;

    static unsigned char bit_flag[] = {1, 2, 4, 8, 16, 32, 64, 128};

    pix = ba;
    for (j=0; j < h; j++)
	for (i=0; i < w; i++)
	    *pix++ = (byte) (gray[*pix] * (WHITE - BLACK));

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

    bdone = 0;
    row_size = (w + 7) / 8;
    for (j=0; j < h; j++)
        for (i=0; i < w; i++) {
	    if (*(bbuffer + row_size*j + i/8) & bit_flag[i % 8])
		_BBSSI (bdone++, ba);
	    else
		_BBCCI (bdone++, ba);
	    }

    free (bbuffer);
    free (error2);
    free (error1);
}


static cell_array (xmin, xmax, ymin, ymax, dx, dy, dimx, colia)
    float xmin, xmax, ymin, ymax;
    int dx, dy, dimx;
    int *colia;
{
    int errind, tnr;
    float x1, y1, x2, y2, tmp;
    byte *rasteraddr;
    int rasterwidth, rasterheight, bitsperpixel;
    BOOL swapx, swapy;

    GQCNTN (&errind, &tnr);
    GSDT (window, viewport);

    x1 = xmin;
    y1 = ymax;
    GNT (&x1, &y1, &tnr);
    GST (&x1, &y1);
    x1 = a * x1 + b;
    y1 = c * y1 + d;

    x2 = xmax;
    y2 = ymin;
    GNT (&x2, &y2, &tnr);
    GST (&x2, &y2);
    x2 = a * x2 + b;
    y2 = c * y2 + d;

    rasterwidth = (int)(fabs(x2 - x1));
    rasterheight = (int)(fabs(y2 - y1));
    rasteraddr = (byte *) malloc (sizeof(byte) * rasterwidth * rasterheight);

    if (swapx = (xmin > xmax)) {
        tmp = x1;
        x1 = x2;
        x2 = tmp;
        }
    if (swapy = (ymin < ymax)) {
        tmp = y1;
        y1 = y2;
        y2 = tmp;
        }

    pack_data (dx, dy, dimx, colia, rasterwidth, rasterheight, rasteraddr,
        swapx, swapy);

    if (map_size <= PRIVATE_COLORS) {
	pixmap_to_bitmap (rasterwidth, rasterheight, rasteraddr);
	bitsperpixel = 1;
	}
    else
	bitsperpixel = 8;

    atb = 0;
    uis$image (&vd_id, &atb, &x1, &y1, &x2, &y2, &rasterwidth, &rasterheight,
	&bitsperpixel, rasteraddr);

    free (rasteraddr);
}



static set_color (color, r, g, b)
    int color;
    float r, g, b;
{
    if (color >= 0 && color < map_size)
	uis$set_color (&vd_id, &color, &r, &g, &b);

    if (color >= PRIVATE_COLORS)
	gray[color - PRIVATE_COLORS] = 0.3*r + 0.59*g + 0.11*b;
}


static button_ast ()
{
    sys$wake (NULL, NULL);
}


static get_pointer (n, x, y, state)
    int *n;
    float *x, *y;
    int *state;
{
    unsigned int kb_id, pattern_count, activex, activey;
    int np, retstate;
    float xpointer, ypointer;

    static unsigned short pattern[] = {
	0x0100, 0x0100, 0x0100, 0x0100, 0x0100, 0x0100, 0x0100, 0xffff,
	0x0100, 0x0100, 0x0100, 0x0100, 0x0100, 0x0100, 0x0100, 0x0100
	};

    kb_id = uis$create_kb (&devnam);
    uis$enable_kb (&kb_id, &wd_id);
    uis$set_button_ast (&vd_id, &wd_id, button_ast, NULL, NULL, NULL, NULL,
	NULL);
    uis$set_pointer_pattern (&vd_id, &wd_id, pattern, &pattern_count,
	&activex, &activey);

    np = 0;
    do
	{
	*state = Undefined;
	do
	    {
	    sys$hiber ();

	    uis$get_buttons (&wd_id, &retstate);

	    if (!(retstate & UIS$M_POINTER_BUTTON_1)) {
		uis$get_pointer_position (&vd_id, &wd_id, &xpointer, &ypointer);

	 	*x++ = (xpointer - b) / a;
		*y++ = (ypointer - d) / c;
		np++;

		*state = GOK;
		}
	    else if (!(retstate & UIS$M_POINTER_BUTTON_2))
		*state = GNONE;
	    }
	while (*state == Undefined);
	}
    while (np < *n && *state != GNONE);

    uis$set_pointer_pattern (&vd_id, &wd_id, NULL, NULL, NULL, NULL);
    uis$set_button_ast (&vd_id, &wd_id, NULL, NULL, NULL, NULL, NULL, NULL);
    uis$disable_kb (&kb_id);
    uis$delete_kb (&kb_id);

    *n = np;
    if (*n > 1)
	*state = GOK;
}

#endif /* VMS */


void GKDGA (int *fctid, int *dx, int *dy, int *dimx, int *ia,
    int *lr1, float *r1, int *lr2, float *r2, int *lc,
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
#ifdef UIS
    switch (*fctid) {

      case 2:
/*
 *  Open workstation
 */
	uis_setup ();

	window[0] = window[2] = 0.0;
	window[1] = window[3] = 1.0;

        viewport[0] = 0; viewport[1] = width / 100;
        viewport[2] = 0; viewport[3] = height / 100;

        set_up_xform (window, viewport);
	break;

      case 3:
/*
 *  Close workstation
 */
	uis$delete_display (&vd_id);
	break;

      case 4:
/*
 *  Activate workstation
 */
	state = GACTIV;
	break;

      case 5:
/*
 *  Deactivate workstation
 */
	state = GINACT;
	break;

      case 6:
/*
 *  Clear workstation
 */
	uis$erase (&vd_id, NULL, NULL, NULL, NULL);
	break;

      case 8:
/*
 *  Update workstation
 */
	end_line ();
	break;

      case 12:
/*
 *  Polyline
 */
        if (state == GACTIV)
            polyline (ia, r1, r2);
	break;

      case 13:
/*
 *  Polymarker
 */
	if (state == GACTIV)
            polymarker (ia, r1, r2);
	break;

      case 14:
/*
 *  Text
 */
        if (state == GACTIV)
            text (r1, r2, chars);
	break;

      case 15:
/*
 *  Fill Area
 */
	if (state == GACTIV)
            fillarea (ia, r1, r2);
	break;

      case 16:
/*
 *  Cell Array
 */
        if (state == GACTIV)
            cell_array (r1[0], r1[1], r2[0], r2[1], *dx, *dy, *dimx, ia);
	break;

      case 48:
/*
 *  Set color representation
 */
	set_color (ia[1], r1[0], r1[1], r1[2]);
	break;

      case 54:
/*
 *  Set workstation window
 */
        window[0] = r1[0];
        window[1] = r1[1];
        window[2] = r2[0];
        window[3] = r2[1];

        set_up_xform (window, viewport);
	break;

      case 55:
/*
 *  Set workstation viewport
 */
        viewport[0] = r1[0];
        viewport[1] = r1[1];
        viewport[2] = r2[0];
        viewport[3] = r2[1];

        resize_window ();

        set_up_xform (window, viewport);
	break;

      case 81:
/*
 *  Request locator
 */
        {
        int n = 1;

        get_pointer (&n, r1, r2, &ia[0]);
	}
	break;

      case 82:
/*
 *  Request stroke
 */
        get_pointer (&ia[2], r1, r2, &ia[0]);
	break;

      case 86:
/*
 *  Request string
 */
	break;
      }
#else
    if (*fctid == 2)
        {
        gks_fprintf (stderr,
	    "GKS: UIS not supported for this type of machine\n");
        ia[0] = 0;
        }
#endif
}

