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
 *      This module contains the C language binding for GLI GKS.
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
#include <string.h>
#include <stdlib.h>

#if !defined(VMS) && !defined(MSDOS) && !defined(_WIN32)
#include <unistd.h>
#endif

#ifndef MSDOS
#include <sys/types.h>
#endif

#if defined (cray) || defined (__SVR4) || defined(MSDOS) || defined(_WIN32)
#include <fcntl.h>
#else
#include <sys/file.h>
#endif

#include "gks.h"
 
#ifdef VMS
#include <descrip.h>
#endif

#ifdef cray
#include <fortran.h>
#endif

#define OK          0
#define MAX_POINTS  2048

#ifdef VAXC
noshare extern int gks_errno;
#else
extern int gks_errno;
#endif

static float *x = NULL, *y = NULL;
static int max_points = 0;


static
void gksrealloc (int n)
{
    if (n > max_points)
    {
        x = (float *) realloc (x, sizeof(float) * n);
        y = (float *) realloc (y, sizeof(float) * n);
        max_points = n;
    }
}


int gopengks (Gfile *errfile, Glong memory)
{
    int errfil, buffer = memory;

    errfil = (errfile != NULL) ? fileno(errfile) : 0;

    GOPKS (&errfil, &buffer);

    if (gks_errno == 0)
    {
        x = (float *) malloc (sizeof(float) * MAX_POINTS);
        y = (float *) malloc (sizeof(float) * MAX_POINTS);
        max_points = MAX_POINTS;
    }

    return gks_errno;
}


int gclosegks (void)
{
    GCLKS ();

    if (gks_errno == 0)
    {
        free (y);
        free (x);
        max_points = 0;
    }

    return gks_errno;
}
 

int gopenws (Gint workstation_id, Gconn *connection, Gwstype *type)
{
    int wkid = workstation_id, conid, wstype;

    conid = (connection != NULL) ?
        open (connection, O_CREAT | O_TRUNC | O_WRONLY, 0644) : 0;
    wstype = (type != NULL) ? *type : 0;

    GOPWK (&wkid, &conid, &wstype);

    return gks_errno;
}
 

int gclosews (Gint workstation_id)
{
    int wkid = workstation_id, errind, conid = 0, wtype;

    GQWKC (&wkid, &errind, &conid, &wtype);
    GCLWK (&wkid);

    if (conid > 0)
        close (conid);

    return gks_errno;
}
 

int gactivatews (Gint workstation_id)
{
    int wkid = workstation_id;

    GACWK (&wkid);

    return gks_errno;
}

 
int gdeactivatews (Gint workstation_id)
{
    int wkid = workstation_id;

    GDAWK (&wkid);

    return gks_errno;
}
 

int gclearws (Gint workstation_id, Gclrflag clearflag)
{
    int wkid = workstation_id, cofl = clearflag;

    GCLRWK (&wkid, &cofl);

    return gks_errno;
}
 

int gupdatews (Gint workstation_id, Gregen regenflag)
{
    int wkid = workstation_id, refl;

    refl = regenflag == GPOSTPONE ? 0 : 1;
    GUWK (&wkid, &refl);

    return gks_errno;
}
 

/*
int gescape (void)
{
} */


int gmessage (Gint workstation_id, Gchar *string)
{
    int wkid = workstation_id, nchars;
    unsigned short chars_len;
    char *chars;
#ifdef VMS
    struct dsc$descriptor_s text;
#endif
#ifdef cray
    _fcd text;
#endif

    chars = string;
    nchars = strlen(chars);
    chars_len = (unsigned short)nchars;
#ifdef VMS
    text.dsc$b_dtype = DSC$K_DTYPE_T;
    text.dsc$b_class = DSC$K_CLASS_S;
    text.dsc$w_length = nchars;
    text.dsc$a_pointer = chars;

    GMSGS (&wkid, &nchars, &text);
#else
#ifdef cray
    text = _cptofcd(chars, nchars);

    GMSGS (&wkid, &nchars, text);
#else
    GMSGS (&wkid, &nchars, chars, chars_len);
#endif /* cray */
#endif /* VMS */

    return gks_errno;
}


int gpolyline (Gint n, Gpoint *points)
{
    int i, npoints = n;
 
    gksrealloc (n);

    for (i = 0; i < n; i++)
    {
        x[i] = points[i].x;
        y[i] = points[i].y;
    }
    GPL (&npoints, x, y);

    return gks_errno;
}
 

int gpolymarker (Gint n, Gpoint *points)
{
    int i, npoints = n;
 
    gksrealloc (n);

    for (i = 0; i < n; i++)
    {
        x[i] = points[i].x;
        y[i] = points[i].y;
    }
    GPM (&npoints, x, y);

    return gks_errno;
}
 

int gtext (Gpoint *position, Gchar *string)
{
    float qx, qy;
    int nchars;
    unsigned short chars_len;
    char *chars;
#ifdef VMS
    struct dsc$descriptor_s text;
#endif
#ifdef cray
    _fcd text;
#endif

    qx = position->x;
    qy = position->y;
    chars = string;
    nchars = strlen(chars);
    chars_len = (unsigned short)nchars;
#ifdef VMS
    text.dsc$b_dtype = DSC$K_DTYPE_T;
    text.dsc$b_class = DSC$K_CLASS_S;
    text.dsc$w_length = nchars;
    text.dsc$a_pointer = chars;

    GTXS (&qx, &qy, &nchars, &text);
#else
#ifdef cray
    text = _cptofcd(chars, nchars);

    GTXS (&qx, &qy, &nchars, text);
#else
    GTXS (&qx, &qy, &nchars, chars, chars_len);
#endif /* cray */
#endif /* VMS */

    return gks_errno;
}
 

int gfillarea (Gint n, Gpoint *points)
{
    int i, npoints = n;
 
    gksrealloc (n);
 
    for (i = 0; i < n; i++)
    {
        x[i] = points[i].x;
        y[i] = points[i].y;
    }
    GFA (&npoints, x, y);

    return gks_errno;
}
 

int gcellarray (Grect *rectangle, Gidim *dimensions, Gint *color)
{
    float qx, qy, rx, ry;
    int dx, dy, scol, srow, ncol, nrow, *colia = color;

    qx = rectangle->ul.x;
    qy = rectangle->ul.y;
    rx = rectangle->lr.x;
    ry = rectangle->lr.y;
    dx = dimensions->x_dim;
    dy = dimensions->y_dim;
    scol = 1;
    srow = 1;
    ncol = dx;
    nrow = dy;

    GCA (&qx, &qy, &rx, &ry, &dx, &dy, &scol, &srow, &ncol, &nrow, colia);

    return gks_errno;
}


int gsetasf (Gasfs *asfs)
{
    int flag[13];

    flag[0]  = asfs->ln_type;
    flag[1]  = asfs->ln_width;
    flag[2]  = asfs->ln_colour;
    flag[3]  = asfs->mk_type;
    flag[4]  = asfs->mk_size;
    flag[5]  = asfs->mk_colour;
    flag[6]  = asfs->tx_fp;
    flag[7]  = asfs->tx_exp;
    flag[8]  = asfs->tx_space;
    flag[9]  = asfs->tx_colour;
    flag[10] = asfs->fl_inter;
    flag[11] = asfs->fl_style;
    flag[12] = asfs->fl_colour;

    GSASF (flag);

    return gks_errno;
}
 

int gsetlineind (Gint index)
{
    int pli = index;

    GSPLI (&pli);

    return gks_errno;
}


int gsetlinetype (Gint type)
{
    int ltype = type;

    GSLN (&ltype);

    return gks_errno;
}
 

int gsetlinewidth (Gfloat width)
{
    float lwidth = width;

    GSLWSC (&lwidth);

    return gks_errno;
}
 

int gsetlinecolourind (Gint colour)
{
    int coli = colour;

    GSPLCI (&coli);

    return gks_errno;
}
 

int gsetmarkerind (Gint index)
{
    int pmi = index;

    GSPMI (&pmi);

    return gks_errno;
}


int gsetmarkertype (Gint type)
{
    int mtype = type;

    GSMK (&mtype);

    return gks_errno;
}
 

int gsetmarkersize (Gfloat size)
{
    float mszsc = size;

    GSMKSC (&mszsc);

    return gks_errno;
}
 

int gsetmarkercolourind (Gint colour)
{
    int coli = colour;

    GSPMCI (&coli);

    return gks_errno;
}
 

int gsettextind (Gint index)
{
    int txi = index;

    GSTXI (&txi);

    return gks_errno;
}


int gsettextfontprec (Gtxfp *txfp)
{
    int font = txfp->font, prec = txfp->prec;

    GSTXFP (&font, &prec);

    return gks_errno;
}
 

int gsetcharexpan (Gfloat exp)
{
    float chxp = exp;

    GSCHXP (&chxp);

    return gks_errno;
}
 

int gsetcharspace (Gfloat spacing)
{
    float chsp = spacing;

    GSCHSP (&chsp);

    return gks_errno;
}
 

int gsettextcolourind (Gint colour)
{
    int coli = colour;

    GSTXCI (&coli);

    return gks_errno;
}
 

int gsetcharheight (Gfloat height)
{
    float chh = height;

    GSCHH (&chh);

    return gks_errno;
}
 

int gsetcharup (Gpoint *charup)
{
    float chux, chuy;

    chux = charup->x;
    chuy = charup->y;

    GSCHUP (&chux, &chuy);

    return gks_errno;
}
 

int gsettextpath (Gtxpath text_path)
{
    int txp = text_path;

    GSTXP (&txp);

    return gks_errno;
}
 

int gsettextalign (Gtxalign *txalign)
{
    int alh, alv;

    alh = txalign->hor;
    alv = txalign->ver;

    GSTXAL (&alh, &alv);

    return gks_errno;
}
 

int gsetfillind (Gint index)
{
    int fai = index;

    GSFAI (&fai);

    return gks_errno;
}


int gsetfillintstyle (Gflinter style)
{
    int ints = style;

    GSFAIS (&ints);

    return gks_errno;
}
 

int gsetfillstyle (Gint index)
{
    int styli = index;

    GSFASI (&styli);

    return gks_errno;
}
 

int gsetfillcolourind (Gint colour)
{
    int coli = colour;

    GSFACI (&coli);

    return gks_errno;
}
 

int gsetcolourrep (Gint workstation_id, Gint index, Gcobundl *rep)
{
    int wkid = workstation_id, coli = index;
    float r, g, b;

    r = rep->red;
    g = rep->green;
    b = rep->blue;

    GSCR (&wkid, &coli, &r, &g, &b);

    return gks_errno;
}
 

int gsetwindow (Gint transform, Glimit *window)
{
    int tnr = transform;
    float xmin, xmax, ymin, ymax;

    xmin = window->xmin;
    xmax = window->xmax;
    ymin = window->ymin;
    ymax = window->ymax;

    GSWN (&tnr, &xmin, &xmax, &ymin, &ymax);

    return gks_errno;
}
 

int gsetviewport (Gint transform, Glimit *viewport)
{
    int tnr = transform;
    float xmin, xmax, ymin, ymax;

    xmin = viewport->xmin;
    xmax = viewport->xmax;
    ymin = viewport->ymin;
    ymax = viewport->ymax;

    GSVP (&tnr, &xmin, &xmax, &ymin, &ymax);

    return gks_errno;
}

 
int gselntran (Gint transform)
{
    int tnr = transform;

    GSELNT (&tnr);

    return gks_errno;
}
 

int gsetclip (Gclip indicator)
{
    int clsw = indicator;

    GSCLIP (&clsw);

    return gks_errno;
}
 

int gsetwswindow (Gint workstation_id, Glimit *window)
{
    int wkid = workstation_id;
    float xmin, xmax, ymin, ymax;

    xmin = window->xmin;
    xmax = window->xmax;
    ymin = window->ymin;
    ymax = window->ymax;

    GSWKWN (&wkid, &xmin, &xmax, &ymin, &ymax);

    return gks_errno;
}

 
int gsetwsviewport (Gint workstation_id, Glimit *viewport)
{
    int wkid = workstation_id;
    float xmin, xmax, ymin, ymax;

    xmin = viewport->xmin;
    xmax = viewport->xmax;
    ymin = viewport->ymin;
    ymax = viewport->ymax;

    GSWKVP (&wkid, &xmin, &xmax, &ymin, &ymax);

    return gks_errno;
}
 

/*
ginitloc (void)
{
} */


int greqloc (Gint workstation_id, Gint device_number, Gqloc *response)
{
    int wkid = workstation_id, lcdnr = device_number;
    int stat, tnr;
    float qx, qy;

    GRQLC (&wkid, &lcdnr, &stat, &tnr, &qx, &qy);

    response->status = (Gistat)stat;
    response->loc.transform = tnr;
    response->loc.position.x = qx;
    response->loc.position.y = qy;

    return gks_errno;
}
 

/*
int greqstroke (void)
{
} */


int greqstring (Gint workstation_id, Gint device_number, Gqstring *response)
{
    int wkid = workstation_id, stdnr = device_number;
    int stat, lostr;

    lostr = strlen(response->string);
#if !defined(VMS) && !defined(cray)
    GRQST (&wkid, &stdnr, &stat, &lostr, response->string,
	(unsigned short)lostr);
#else
    GRQST (&wkid, &stdnr, &stat, &lostr, response->string);
#endif

    response->status = (Gistat)stat;
    response->string[lostr] = '\0';

    return gks_errno;
}
 

int gcreateseg (Gint segment_name)
{
    int segn = segment_name;

    GCRSG (&segn);

    return gks_errno;
}

 
int gcopysegws (Gint workstation_id, Gint segment_name)
{
    int wkid = workstation_id, segn = segment_name;

    GCSGWK (&wkid, &segn);

    return gks_errno;
}
 

int gredrawsegws (Gint workstation_id)
{
    int wkid = workstation_id;

    GRSGWK (&wkid);

    return gks_errno;
}
 

int gcloseseg (void)
{
    GCLSG ();

    return gks_errno;
}
 

int gevaltran (Gpoint *ppoint, Gpoint *pshift, Gfloat angle, Gscale *pscale,
    Gcsw coord, Gfloat result[3][2])
{
    float x0, y0, tx, ty, phi, fx, fy;
    int isw, i, j;    
    float tran[3][2];

    x0  = ppoint->x;
    y0  = ppoint->y;
    tx  = pshift->x;
    ty  = pshift->y;
    phi = angle;
    fx  = pscale->x_scale;
    fy  = pscale->y_scale;
    isw = coord;

    GEVTM (&x0, &y0, &tx, &ty, &phi, &fx, &fy, &isw, tran);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 3; j++)
            result[i][j] = tran[i][j];

    return gks_errno;
}


int gsetsegtran (Gint segment_name, Gfloat segtran[3][2])
{
    int segn = segment_name, i, j;
    float tran[3][2];

    for (i = 0; i < 2; i++)
        for (j = 0; j < 3; j++)
            tran[i][j] = segtran[i][j];

    GSSGT (&segn, tran);

    return gks_errno;
}


int ginqopst (Gint *state)
{
    int opsta;

    GQOPS (&opsta);

    *state = opsta;

    return OK;
}
 

int ginqlevelgks (Gint *level, Gint *error_status)
{
    int errind, lev;

    GQLVKS (&errind, &lev);

    *level = lev;
    *error_status = errind;

    return OK;
}
 

int ginqmaxntrannum (Gint *maxtran, Gint *error_status)
{
    int errind, maxtnr;

    GQMNTN (&errind, &maxtnr);

    *maxtran = maxtnr;
    *error_status = errind;

    return OK;
}
 

/*
int ginqopenws (void)
{
} */


/*
int ginqactivews (void)
{
} */


int ginqcharheight (Gfloat *height, Gint *error_status)
{
    int errind;
    float chh;

    GQCHH (&errind, &chh);

    *height = chh;
    *error_status = errind;

    return OK;
}
 

int ginqcharup (Gpoint *up, Gint *error_status)
{
    float chux, chuy;
    int errind;

    GQCHUP (&errind, &chux, &chuy);

    up->x = chux;
    up->y = chuy;
    *error_status = errind;

    return OK;
}
 

int ginqtextpath (Gtxpath *text_path, Gint *error_status)
{
    int errind, path;

    GQTXP (&errind, &path);

    *text_path = (Gtxpath)path;
    *error_status = errind;

    return OK;
}


int ginqtextalign (Gtxalign *align, Gint *error_status)
{
    int errind, alh, alv;

    GQTXAL (&errind, &alh, &alv);

    align->hor = (Gtxhor)alh;
    align->ver = (Gtxver)alv;

    return OK;
}

 
int ginqasf (Gasfs *asflist, Gint *error_status)
{
    int errind, flag[13];

    GQASF (&errind, flag);

    asflist->ln_type    = (Gasf)flag[0];
    asflist->ln_width   = (Gasf)flag[1];
    asflist->ln_colour  = (Gasf)flag[2];
    asflist->mk_type    = (Gasf)flag[3];
    asflist->mk_size    = (Gasf)flag[4];
    asflist->mk_colour  = (Gasf)flag[5];
    asflist->tx_fp      = (Gasf)flag[6];
    asflist->tx_exp     = (Gasf)flag[7];
    asflist->tx_space   = (Gasf)flag[8];
    asflist->tx_colour  = (Gasf)flag[9];
    asflist->fl_inter   = (Gasf)flag[10];
    asflist->fl_style   = (Gasf)flag[11];
    asflist->fl_colour  = (Gasf)flag[12];
    *error_status = errind;

    return OK;
}
 

int ginqlineind (Gint *index, Gint *error_status)
{
    int errind, pli;

    GQPLI (&errind, &pli);

    *index = pli;
    *error_status = errind;

    return OK;
}


int ginqlinetype (Gint *type, Gint *error_status)
{
    int errind, ltype;

    GQLN (&errind, &ltype);

    *type = ltype;
    *error_status = errind;

    return OK;
}


int ginqlinewidth (Gfloat *width, Gint *error_status)
{
    int errind;
    float lwidth;

    GQLWSC (&errind, &lwidth);

    *width = lwidth;
    *error_status = errind;

    return OK;
}


int ginqlinecolourind (Gint *index, Gint *error_status)
{
    int errind, coli;

    GQPLCI (&errind, &coli);

    *index = coli;
    *error_status = errind;

    return OK;
}
 

int ginqmarkerind (Gint *index, Gint *error_status)
{
    int errind, pmi;

    GQPMI (&errind, &pmi);

    *index = pmi;
    *error_status = errind;

    return OK;
}


int ginqmarkertype (Gint *level, Gint *error_status)
{
    int errind, mtype;

    GQMK (&errind, &mtype);

    *level = mtype;
    *error_status = errind;

    return OK;
}
 

int ginqmarkersize (Gfloat *mksize, Gint *error_status)
{
    int errind;
    float mszsc;

    GQMKSC (&errind, &mszsc);

    *mksize = mszsc;
    *error_status = errind;

    return OK;
}
 

int ginqmarkercolourind (Gint *index, Gint *error_status)
{
    int errind, coli;

    GQPMCI (&errind, &coli);

    *index = coli;
    *error_status = errind;

    return OK;
}
 

int ginqtextind (Gint *index, Gint *error_status)
{
    int errind, txi;

    GQTXI (&errind, &txi);

    *index = txi;
    *error_status = errind;

    return OK;
}


int ginqtextfontprec (Gtxfp *txfp, Gint *error_status)
{
    int errind, font, prec;

    GQTXFP (&errind, &font, &prec);

    txfp->font = font;
    txfp->prec = (Gtxprec)prec;
    *error_status = errind;

    return OK;
}
 

int ginqcharexpan (Gfloat *chexp, Gint *error_status)
{
    int errind;
    float chxp;

    GQCHXP (&errind, &chxp);

    *chexp = chxp;
    *error_status = errind;

    return OK;
}
 

int ginqcharspace (Gfloat *chspc, Gint *error_status)
{
    int errind;
    float chsp;

    GQCHSP (&errind, &chsp);

    *chspc = chsp;
    *error_status = errind;

    return OK;
}
 

int ginqtextcolourind (Gint *index, Gint *error_status)
{
    int errind, coli;

    GQTXCI (&errind, &coli);

    *index = coli;
    *error_status = errind;

    return OK;
}
 

int ginqfillind (Gint *index, Gint *error_status)
{
    int errind, fai;

    GQFAI (&errind, &fai);

    *index = fai;
    *error_status = errind;

    return OK;
}
 

int ginqfillintstyle (Gint *style, Gint *error_status)
{
    int ints, errind;

    GQFAIS (&errind, &ints);

    *style = ints;
    *error_status = errind;

    return OK;
}
 

int ginqfillstyle (Gint *index, Gint *error_status)
{
    int styli, errind;

    GQFASI (&errind, &styli);

    *index = styli;
    *error_status = errind;

    return OK;
}
 

int ginqfillcolourind (Gint *index, Gint *error_status)
{
    int errind, coli;

    GQFACI (&errind, &coli);

    *index = coli;
    *error_status = errind;

    return OK;
}


int ginqcurntrannum (Gint *tran, Gint *error_status)
{
    int errind, tnr;

    GQCNTN (&errind, &tnr);

    *tran = tnr;
    *error_status = errind;

    return OK;
}
 

int ginqntran (Gint num, Gtran *tran, Gint *error_status)
{
    int tnr = num, errind;
    float wn[4], vp[4];

    GQNT (&tnr, &errind, wn, vp);

    tran->w.xmin = wn[0];
    tran->w.xmax = wn[1];
    tran->w.ymin = wn[2];
    tran->w.ymax = wn[3];
    tran->v.xmin = vp[0];
    tran->v.xmax = vp[1];
    tran->v.ymin = vp[2];
    tran->v.ymax = vp[3];
    *error_status = errind;

    return OK;
}
 

int ginqclip (Gcliprect *clipping, Gint *error_status)
{
    int errind, clsw;
    float clrt[4];

    GQCLIP (&errind, &clsw, clrt);

    clipping->rec.xmin = clrt[0];
    clipping->rec.xmax = clrt[1];
    clipping->rec.ymin = clrt[2];
    clipping->rec.ymax = clrt[3];
    clipping->ind = (Gclip)clsw;
    *error_status = errind;

    return OK;
}
 

/*
int ginqwsconntype (void)
{
} */


int ginqwsst (Gint workstation_id, Gint *state, Gint *error_status)
{
    int wkid = workstation_id, errind, wsstate;

    GQWKS (&wkid, &errind, &wsstate);

    *state = wsstate;
    *error_status = errind;

    return OK;
}
 

int ginqwscategory (Gwstype *workstation_type, Gint *cat, Gint *error_status)
{
    int wstype = *workstation_type;
    int errind, category;

    GQWKCA (&wstype, &errind, &category);

    *cat = category;
    *error_status = errind;

    return OK;
}
 

int ginqwscf (Gwstype *workstation_type, Gint *ncolintensities, Gint *colavail,
    Gint *npcinds, Gint *error_status)
{
    int wstype = *workstation_type;
    int errind, ncoli, cola, npci;

    GQCF (&wstype, &errind, &ncoli, &cola, &npci);

    *ncolintensities = ncoli;
    *colavail = cola;
    *npcinds = npci;
    *error_status = errind;

    return OK;
}
 

int ginqdisplaysize (Gwstype *workstation_type, Gdspsize *dspsz,
    Gint *error_status)
{
    int wstype = *workstation_type;
    int errind, units, ras_x, ras_y;
    float px, py;
 
    GQDSP (&wstype, &errind, &units, &px, &py, &ras_x, &ras_y);

    dspsz->units = (Gdevunits)units;
    dspsz->device.x = px;
    dspsz->device.y = py;
    dspsz->raster.x = ras_x;
    dspsz->raster.y = ras_y;
    *error_status = errind;

    return OK;
}
 

int ginqtextextent (Gint workstation_id, Gpoint *position, Gchar *string,
    Gextent *extent, Gint *error_status)
{
    int wkid = workstation_id, errind, nchars;
    unsigned short chars_len;
    float qx, qy, cpx, cpy, tx[4], ty[4];
    char *chars;
#ifdef VMS
    struct dsc$descriptor_s text;
#endif
#ifdef cray
    _fcd text;
#endif

    qx = position->x;
    qy = position->y;
    chars = string;
    nchars = strlen(chars);
    chars_len = (unsigned short)nchars;
#ifdef VMS
    text.dsc$b_dtype = DSC$K_DTYPE_T;
    text.dsc$b_class = DSC$K_CLASS_S;
    text.dsc$w_length = nchars;
    text.dsc$a_pointer = chars;

    GQTXXS (&wkid, &qx, &qy, &nchars, &text, &errind, &cpx, &cpy, tx, ty);
#else
#ifdef cray
    text = _cptofcd(chars, nchars);

    GQTXXS (&wkid, &qx, &qy, &nchars, text, &errind, &cpx, &cpy, tx, ty);
#else
    GQTXXS (&wkid, &qx, &qy, &nchars, chars, &errind, &cpx, &cpy, tx, ty,
	chars_len);
#endif /* cray */
#endif /* VMS */

    extent->concat.x = cpx;
    extent->concat.y = cpy;
    extent->corner_1.x = tx[0];
    extent->corner_1.y = ty[0];
    extent->corner_2.x = tx[1];
    extent->corner_2.y = ty[1];
    extent->corner_3.x = tx[2];
    extent->corner_3.y = ty[2];
    extent->corner_4.x = tx[3];
    extent->corner_4.y = ty[3];

    return OK;
}
 

int ginqnameopenseg (Gint *segment_name, Gint *error_status)
{
    int errind, segn;

    GQOPSG (&errind, &segn);

    *segment_name = segn;
    *error_status = errind;

    return OK;
}


int gemergencyclosegks (void)
{
    GECLKS ();

    return OK;
}
