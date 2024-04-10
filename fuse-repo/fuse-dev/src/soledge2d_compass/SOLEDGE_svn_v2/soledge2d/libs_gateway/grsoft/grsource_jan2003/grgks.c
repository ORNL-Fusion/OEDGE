/*
 * Copyright @ 1993, 1994   Josef Heinen
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
 *      GR/GR3 Software
 *
 * ABSTRACT:
 *
 *      This module contains the GLIGKS layer for the GR/GR3 software.
 *
 * AUTHOR:
 *
 *      Josef Heinen
 *
 * VERSION:
 *
 *      V1.0-00
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include <sys/types.h>

#if defined (cray) || defined (__SVR4)
#include <fcntl.h>
#else
#include <sys/file.h>
#endif

#ifdef VMS
#include <descrip.h>
#endif

#ifdef cray
#include <fortran.h>
#endif

/* Entry point definitions */

#ifndef cray
#if defined (VMS) || defined (hpux) || defined (aix)

#define GRCONF grconf
#define GRPAN grpan
#define GRGIF grgif
#define GRQSCL grqscl

#define GACWK gacwk
#define GCLRWK gclrwk
#define GCLSG gclsg
#define GCLWK gclwk
#define GCRSG gcrsg
#define GCSGWK gcsgwk
#define GDAWK gdawk
#define GDSG gdsg
#define GECLKS geclks
#define GOPWK gopwk
#define GQCNTN gqcntn
#define GQDSP gqdsp
#define GQNT gqnt
#define GQOPS gqops
#define GQOPWK gqopwk
#define GQWKC gqwkc
#define GQWKCA gqwkca
#define GRQLC grqlc
#define GSCR gscr
#define GSWKVP gswkvp
#define GSWKWN gswkwn
#define GUWK guwk

#else

#define GRCONF grconf_
#define GRPAN grpan_
#define GRGIF grgif_
#define GRQSCL grqscl_

#define GACWK gacwk_
#define GCLRWK gclrwk_
#define GCLSG gclsg_
#define GCLWK gclwk_
#define GCRSG gcrsg_
#define GCSGWK gcsgwk_
#define GDAWK gdawk_
#define GDSG gdsg_
#define GECLKS geclks_
#define GOPWK gopwk_
#define GQCNTN gqcntn_
#define GQDSP gqdsp_
#define GQNT gqnt_
#define GQOPS gqops_
#define GQOPWK gqopwk_
#define GQWKC gqwkc_
#define GQWKCA gqwkca_
#define GRQLC grqlc_
#define GSCR gscr_
#define GSWKVP gswkvp_
#define GSWKWN gswkwn_
#define GUWK guwk_

#endif
#endif /* cray */

#ifdef VMS
#include <descrip.h>
#endif /* VMS */

#define WSAC 3                      /* at least one workstation active */
#define OUTPT 0                     /* output workstation */
#define OUTIN 2                     /* output/input workstation */
#define MO 4                        /* metafile output */

#define Bool int
#define True 1
#define False 0

#define NoOutput -1
#define Segment 1

#define cl "\33[H\33[J"             /* clear screen and home cursor */

#define odd(status) ((status) & 01)
#define max(a,b) (((a) > (b)) ? (a) : (b))

typedef char String[80];
typedef unsigned char byte;

typedef struct {
    int wkid, wktyp;
    char *env, *name, *type;
    FILE *stream;
} ws_cap;

static
ws_cap ws_info[] = {
    { 0,  NoOutput, "", "", "no output", NULL },
    { 0,  16, "", "", "VT330", NULL },
    { 0,  17, "", "", "VT340", NULL },
    { 0,  72, "", "", "TEK 401x", NULL },
    { 0,  82, "", "", "TEK 42xx", NULL },
    { 0, 201, "", "", "TAB 132/15-G", NULL },
    { 0, 204, "", "", "Monterey", NULL },
    { 0, 207, "", "", "IBM/PC", NULL },
    { 0, 210, "", "", "X display (output only)", NULL },
    { 0, 211, "", "", "X display", NULL },
    { 0, 217, "", "", "X display w\\ frame grabber\n", NULL },
    { 0,   7, "GLI_CGM", "gli%2d.cgm", "CGM Binary", NULL },
    { 0,   8, "GLI_CGM", "gli%2d.cgm", "CGM Clear Text", NULL },
    { 0,  38, "GLI_DEC", "gli%2d.dec", "DEC LN03+", NULL },
    { 0,  53, "GLI_HP", "gli%2d.hp", "HP-GL", NULL },
    { 0,  61, "GLI_EPS", "gli%2d.eps", "PostScript (b/w)", NULL },
    { 0,  62, "GLI_EPS", "gli%2d.eps", "Color PostScript", NULL },
    { 0,  63, "", "", "Display PostScript (b/w) w\\ Compuserve GIF dump",
	NULL },
    { 0,  64, "", "", "Color Display PostScript w\\ Compuserve GIF dump",
	NULL },
    { 0, 104, "GLI_PBM", "gli%2d.pbm", "PBM (Portable BitMap)", NULL },
    { 0, 214, "GLI_RF", "gli%2d.rf", "X display w\\ Sun rle rasterfile dump",
	NULL },
    { 0, 215, "GLI_GIF", "gli%2d.gif", "X display w\\ Compuserve GIF dump",
	NULL }
};
static
int n_ws_infos = sizeof (ws_info) / sizeof (ws_info[1]);

static
int n_frame = 0, def_wsid = 0, def_wstype = 0, def_devtype = NoOutput;

static
ws_cap *ws, *dev;

#ifdef VAXC
noshare long xtConId = 0;
#else
#if defined (__ALPHA) || defined (__alpha)
long xtConId = 0;
#else
int xtConId = 0;
#endif
#endif

typedef struct {
    float r, g, b;
} rgb_color;

static
rgb_color color[] = {
    { 1, 1, 1 },
    { 0, 0, 0 },
    { 1, 0, 0 },
    { 0, 0, 1 },
    { 0, 1, 0 },
    { 1, 0, 1 },
    { 1, 1, 0 },
    { 0, 1, 1 },
};

#ifdef VMS

static
int str_locate (str, ch)
    char *str;
    char ch;
{
    int i = 0;
    while (*(str + i) && str[i] != ch)
        i++;
    return i;
}

static
void putenv (env)
    char *env;
{
    char log_name[255];
    int pos, stat;
    char equ_name[255];
    struct dsc$descriptor_s logical_name, value_string;

    strcpy (log_name, env);
    pos = str_locate(env, '=');
    log_name[pos++] = '\0';

    logical_name.dsc$b_dtype = DSC$K_DTYPE_T;
    logical_name.dsc$b_class = DSC$K_CLASS_S;
    logical_name.dsc$w_length = strlen(log_name);
    logical_name.dsc$a_pointer = log_name;

    if (pos < strlen(env))
	strcpy (equ_name, env+pos);
    else
	strcpy (equ_name, " ");

    value_string.dsc$b_dtype = DSC$K_DTYPE_T;
    value_string.dsc$b_class = DSC$K_CLASS_S;
    value_string.dsc$w_length = strlen(equ_name);
    value_string.dsc$a_pointer = equ_name;

    stat = LIB$SET_LOGICAL (&logical_name, &value_string, NULL, NULL, NULL);
    if (!odd(stat))
        LIB$SIGNAL (stat);
}

#endif /* VMS */

static
ws_cap *id(wstype)
    int wstype;
{
    int i;

    for (i=0; i<n_ws_infos; i++)
        if (wstype == ws_info[i].wktyp)
            return &ws_info[i];

    return NULL;
}

static
void panel ()
{
    int i, wsid = def_wsid, wstype = def_wstype, devtype = def_devtype;
    String s;

    while (!(ws = id(wstype)) || !(dev = id(devtype)))
    {
	printf ("%s", cl);
	for (i=0; i<n_ws_infos; i++)
	    printf ("%4d:    %s\n", ws_info[i].wktyp, ws_info[i].type);
	printf ("? ");
	gets (s);

	switch (sscanf(s, "%d %d", &wstype, &devtype))
	{
	    default:
	    case 0 : wstype = NoOutput;
	    case 1 : devtype = NoOutput;
	    case 2 : break;
	}
	if (devtype > 0 && devtype == wstype)
	    devtype = NoOutput;
    }
    if (wstype > 0)
        ws->wkid = ++wsid;
    if (devtype > 0)
        dev->wkid = ++wsid;

    if (wstype > 0 || devtype > 0)
	n_frame++;
}

static
void req_loc (ws)
    ws_cap *ws;
{
    int errind, conid, wtype;
    int wkcat, lcdnr = 1, inp_dev_stat, inp_tnr;
    float px, py;

    if (!*(ws->env))
    {
        GQWKC (&ws->wkid, &errind, &conid, &wtype);
        GQWKCA (&wtype, &errind, &wkcat);

        if (wkcat == OUTIN && isatty(0))
	{
            GRQLC (&ws->wkid, &lcdnr, &inp_dev_stat, &inp_tnr, &px, &py);
            if (inp_dev_stat == 0)
	    {
		GECLKS ();
                exit (0);
	    }
        }
    }
}

static
void fit_ws_viewport (xcm, ycm, rx, ry)
    float xcm, ycm, *rx, *ry;
{
    float xratio, yratio;

    if (*rx > xcm && *ry > ycm)
    {
	*rx = xcm;
	*ry = ycm;
    }
    else
    {
        xratio = xcm / *rx;
	yratio = ycm / *ry;

	if (yratio > xratio) {
	    *rx = xcm / yratio;
	    *ry = ycm / yratio;
	} else {
	    *rx = xcm / xratio;
	    *ry = ycm / xratio;
	}
    }
}

static
void open_wk (ws)
    ws_cap *ws;
{
#if defined (__ALPHA) || defined (__alpha)
    long conid;
#else
    int conid;
#endif
    int coli, err, devunits, lx, ly;
    float rx, ry, xcm, ycm, aspect_ratio;
    float xmin, xmax, ymin, ymax;
    char *cp, *ep, path[100];

    if (ws->wktyp == 211 && xtConId > 0)
    {
	conid = xtConId;
    }
    else if (*(ws->env))
    {
	if ((cp = getenv(ws->env)) == NULL)
	    cp = ws->name;

        sprintf (path, cp, n_frame);
	cp = path;
	while (isspace(*cp))
	    cp++;

#ifndef VMS
	if (*cp != '|')
	{
#endif
	    for (cp = path; *cp; cp++)
		if (*cp == ' ')
		    *cp = '0';

	    if (n_frame == 0)
	    {
		ep = path;
		for (cp = path; *cp; cp++)
		    if (!isdigit(*cp))
			*ep++ = *cp;
		*ep = '\0';
	    }
#ifndef VMS
        }

	if (*cp == '|')
	{
	    if ((ws->stream = popen (cp+1, "w")) != NULL)
		conid = fileno(ws->stream);
	    else
		conid = -1;
	}
	else
#endif
	    conid = open (path, O_CREAT | O_TRUNC | O_WRONLY, 0644);

        if (conid <= 0)
	{
            fprintf (stderr, "grgks: can't create output file: %s\n", path);
	    conid = 0;
	}
    }
    else
        conid = 0;

    GOPWK (&ws->wkid, &conid, &ws->wktyp);
    GACWK (&ws->wkid);

    GQDSP (&ws->wktyp, &err, &devunits, &rx, &ry, &lx, &ly);

    GRQSCL (&xcm, &ycm);
    aspect_ratio = ycm / xcm;

    if (ws->wktyp == 217)
    {
	if (aspect_ratio < 1) {
            rx = rx * (640.0 / lx / aspect_ratio);
	    ry = rx * aspect_ratio;
	} else {
            ry = ry * (640.0 / ly * aspect_ratio);
	    rx = ry / aspect_ratio;
	}
    }
    else
    {
	if (ws->wktyp >= 210 && ws->wktyp <= 215)
	{
	    if (rx > ry)
		rx *= 0.95;
	    else
        	ry *= 0.95;
	}
	fit_ws_viewport (xcm, ycm, &rx, &ry);
    }

    xmin = ymin = 0.0;
    xmax = rx; ymax = ry;
    GSWKVP (&ws->wkid, &xmin, &xmax, &ymin, &ymax);

    if (aspect_ratio < 1) {
        xmax = 1.0; ymax = aspect_ratio;
    } else {
        xmax = 1.0 / aspect_ratio; ymax = 1.0;
    }
    GSWKWN (&ws->wkid, &xmin, &xmax, &ymin, &ymax);

    /* set colors */
    for (coli=0; coli<8; coli++)
        GSCR (&ws->wkid, &coli, &color[coli].r, &color[coli].g, &color[coli].b);
}

static
void clear_wk ()
{
    int count, n = 1, errind, ol, wkid;
    int state, wkcat, conid, wtype, always = 1;

    GQOPS (&state);
    if (state == WSAC)
    {
        GQOPWK (&n, &errind, &ol, &wkid);
        for (count = 1; count <= ol; count++)
        {
            n = count;
            GQOPWK (&n, &errind, &ol, &wkid);

            GQWKC (&wkid, &errind, &conid, &wtype);
            GQWKCA (&wtype, &errind, &wkcat);

            if (wkcat == OUTPT || wkcat == OUTIN || wkcat == MO)
                GCLRWK (&wkid, &always);
        }
    }
}

static
void update_wk ()
{
    int count, n = 1, errind, ol, wkid;
    int conid, wtype, wkcat, perform = 1;

    GQOPWK (&n, &errind, &ol, &wkid);
    for (count = 1; count <= ol; count++)
    {
        n = count;
        GQOPWK (&n, &errind, &ol, &wkid);

        GQWKC (&wkid, &errind, &conid, &wtype);
        GQWKCA (&wtype, &errind, &wkcat);

        if (wkcat == OUTPT || wkcat == OUTIN)
	    GUWK (&wkid, &perform);
    }
}

static
void close_wk (ws)
    ws_cap *ws;
{
    int err, conid, wstype, always = 1;

    GQWKC (&ws->wkid, &err, &conid, &wstype);
    if (wstype == 217)
        GCLRWK (&ws->wkid, &always);

    GDAWK (&ws->wkid);
    GCLWK (&ws->wkid);

#ifndef VMS
    if (ws->stream != NULL)
    {
	pclose (ws->stream);
	ws->stream = NULL;
    }
    else if (conid > 2)
        close (conid);
#else
    if (conid > 2)
        close (conid);
#endif
}

static
void copy_sg_to_wk (ws)
    ws_cap *ws;
{
    int segnr = Segment;

    open_wk (ws);

    /* copy segment */
    GCSGWK (&ws->wkid, &segnr);

    req_loc (ws);
    close_wk (ws);
}

static
void open_wiss ()
{
    int wiss = ++def_wsid, conid = 88, wstype = 5;
    int segnr = Segment;

    GOPWK (&wiss, &conid, &wstype);
    GACWK (&wiss);
    GCRSG (&segnr);
}

static
ws_cap *configure (wstype)
    int *wstype;
{
    ws_cap *ws;

    if (*wstype)
    {
	if (ws = id(*wstype))
	{
	    *wstype = ws->wktyp;
	    if (*wstype > 0)
	    {
		ws->wkid = ++def_wsid;
		open_wk (ws);
	    }
	}
	else
	{
	    fprintf (stderr, "grgks: invalid workstation type: %d\n", *wstype);
	    ws = id(*wstype = NoOutput);
	}
    }
    else
	ws = id(NoOutput);

    return (ws);
}

void GRCONF (wstype)
    int *wstype;
{
    char *env, path[100];
    FILE *file;

    n_frame = 0;

    def_wsid = 0;
    def_wstype = 0;
    def_devtype = NoOutput;

    if (*wstype != 35)
    {
        if (env = getenv("GLI_WSTYPE"))
        {
            def_wstype = atoi(env);
            def_devtype = *wstype;
        }
    }
    else if (env = getenv("GRSOFT_DEVICE"))
    {
        switch (sscanf(env, "%d %d", &def_wstype, &def_devtype))
        {
            default:
            case 0 : def_wstype = NoOutput;
            case 1 : def_devtype = NoOutput;
            case 2 : break;
        }
    }
    else
    {
        if (env = getenv("HOME"))
        {
            strcpy (path, env);
#ifndef VMS
            strcat (path, "/");
#endif
        }
        else
            strcpy (path, "");

        strcat (path, "grsoft.device");
        if (file = fopen(path, "r"))
        {
            switch (fscanf(file, "%d %d", &def_wstype, &def_devtype))
            {
                default:
                case 0 : def_wstype = NoOutput;
                case 1 : def_devtype = NoOutput;
                case 2 : break;
            }
            fclose (file);
        }
    }

    putenv ("GLI_GKS=GKSGRAL");
    putenv ("GLI_GKS_CMAP_EXTENT=84");

    ws = configure (&def_wstype);
    dev = configure (&def_devtype);

    if (!def_wstype || !def_devtype || (def_wstype == -1 && def_devtype == -1))
	open_wiss ();
}

void GRPAN (closewk)
    int *closewk;
{
    int wiss = def_wsid, segnr = Segment;

    if (!def_wstype || !def_devtype || (def_wstype == -1 && def_devtype == -1))
    {
	/* close segment */
	GCLSG ();

	/* deactivate WISS */
	GDAWK (&wiss);

	update_wk ();
	panel ();
        clear_wk ();

	if (dev->wktyp > 0)
	    copy_sg_to_wk (dev);
	if (ws->wktyp > 0)
	    copy_sg_to_wk (ws);

	/* activate WISS */
	GACWK (&wiss);

	GDSG (&segnr);
        GCRSG (&segnr);
    }
    else
    {
	if (ws && xtConId == 0)
	    req_loc (ws);

	if (ws && *closewk)
	    close_wk (ws);
	else
	{
	    update_wk ();
	    if (xtConId == 0)
		clear_wk ();
	}
    }
}

static write_gif_word (fd, w)
    int fd, w;
{
    byte c;

    c = (w & 0xff);
    write (fd, &c, 1);
    c = ((w >> 8) & 0xff);
    write (fd, &c, 1);
}

void GRGIF (name, ncolors, rgb, width, height, colia, len_name)
#ifdef VMS
    struct dsc$descriptor *name;
#else
#ifdef cray
    _fcd name;
#else
    char *name;
#endif /* cray */
#endif /* VMS */
    int *ncolors;
    float *rgb;
    int *width, *height;
#ifdef VMS
    struct dsc$descriptor *colia;
#else
#ifdef cray
    _fcd colia;
#else
    char *colia;
#endif /* cray */
#endif /* VMS */
    int len_name;
{
    char gif_name[255];
    int fd, size, besize;
    byte c, r, g, b, *pix, *beimage;
    register int i, j;
    int BitsPerPixel, ColorMapSize, InitCodeSize;

#ifdef VMS
    strncpy (gif_name, name->dsc$a_pointer, name->dsc$w_length);
    gif_name[name->dsc$w_length] = '\0';
    pix = (byte *) colia->dsc$a_pointer;
#else
#ifdef cray
    strncpy (gif_name, _fcdtocp(name), _fcdlen(name));
    gif_name[_fcdlen(name)] = '\0';
    pix = (byte *) _fcdtocp(colia);
#else
    strncpy (gif_name, name, len_name);
    gif_name[len_name] = '\0';
    pix = (byte *) colia;
#endif /* cray */
#endif /* VMS */

    fd = open(gif_name, O_CREAT | O_TRUNC | O_WRONLY, 0644);
    if (fd >= 0)
    {
	for (BitsPerPixel = 1; BitsPerPixel < 8; BitsPerPixel++)
	    if ((1 << BitsPerPixel) > *ncolors) break;

	/* write the GIF header */

	write (fd, "GIF87a", 6);

	write_gif_word (fd, *width);	/* screen descriptor */
	write_gif_word (fd, *height);

	c = 0x80;			/* yes, there is a color map */
	c |= (8 - 1) << 4;		/* OR in the color resolution */
	c |= (BitsPerPixel - 1);	/* OR in the # of bits per pixel */
	write (fd, &c, 1);

	c = 0x0;
	write (fd, &c, 1);		/* background color */
	write (fd, &c, 1);		/* future expansion byte */

	/* write colormap */

	ColorMapSize = 1 << BitsPerPixel;

	for (i = 0, j = 0; i < *ncolors; i++) {
	    r = (byte) (255 * rgb[j++]);
	    g = (byte) (255 * rgb[j++]);
	    b = (byte) (255 * rgb[j++]);
	    write (fd, &r, 1);
	    write (fd, &g, 1);
	    write (fd, &b, 1);
	}
	for ( ; i < ColorMapSize; i++) {
	    r = g = b = 0;
	    write (fd, &r, 1);
	    write (fd, &g, 1);
	    write (fd, &b, 1);
	}

	write (fd, ",", 1);		/* image separator */

	write_gif_word (fd, 0);		/* image header */
	write_gif_word (fd, 0);
	write_gif_word (fd, *width);
	write_gif_word (fd, *height);

	write (fd, &c, 1);

	size = *width * *height;
	beimage = (byte *) malloc (sizeof(byte) * size * 3 / 2);

	if (beimage != NULL)
        {
	    InitCodeSize = max(BitsPerPixel, 2);
	    grcompress (InitCodeSize + 1, pix, size, beimage, &besize);

	    c = InitCodeSize;
	    write (fd, &c, 1);
	    if (write (fd, beimage, besize) != besize) {
		fprintf (stderr, "GKS: can't write GIF file\n");
		perror ("write");
            }

	    free (beimage);
        }
	else
	    fprintf (stderr, "GKS: can't allocate temporary storage\n");

	write (fd, "\0", 1);		/* write out a zero-length packet */
	write (fd, ";", 1);		/* terminator */

	close (fd);
    }
    else
	fprintf (stderr, "grgks: can't open GIF file\n");
}
