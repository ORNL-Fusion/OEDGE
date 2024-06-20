/*
 * Copyright @ 1984 - 1994   Josef Heinen
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
 *	Graphical Kernel System (GKS)
 *
 * ABSTRACT:
 *
 *	This module contains some I/O routines for GLI GKS.
 *
 * AUTHOR(S):
 *
 *	J. Heinen
 *
 * VERSION:
 *
 *	V1.0
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>

#ifndef MSDOS
#include <sys/types.h>
#endif

#if !defined(VMS) && !defined(MSDOS) && !defined(_WIN32)
#include <unistd.h>
#endif

#ifdef hpux
#include <sys/utsname.h>
#endif

#if defined (cray) || defined (__SVR4) || defined(MSDOS) || defined(_WIN32)
#include <fcntl.h>
#else
#include <sys/file.h>
#endif

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef VMS
#include <descrip.h>
#include <iodef.h>
#include <ssdef.h>
#else
#include <sys/stat.h>
#ifndef _WIN32
#ifdef NO_TERMIO
#include <sys/ioctl.h>
#include <sgtty.h>
#else
#include <termios.h>
#endif /* NO_TERMIO */
#endif /* _WIN32 */
#endif /* VMS */

#ifdef cray
#include <fortran.h>
#endif

#ifdef DPS
#include <X11/Xlib.h>

static char *options[] = {
    "-sixel", "-gif"
};
#endif /* DPS */

#ifdef sun
extern char *sys_errlist[];
extern int sys_nerr;
#define strerror(errno) (errno < sys_nerr) ? sys_errlist[errno] : NULL
#endif

#ifdef VMS
#define CHARARG(a) struct dsc$descriptor *a
#define CHARPAR(a) a
#else
#ifdef cray
#define CHARARG(a) _fcd a
#define CHARPAR(a) a
#else
#if defined(_WIN32) && !defined(__GNUC__)
#define CHARARG(a) char *(a), unsigned short a##_len
#define CHARPAR(a) a, (unsigned short) a##_len
#else
#define CHARARG(a) char *a
#define CHARPAR(a) a
#endif /* _WIN32 */
#endif /* cray */
#endif /* VMS */

#if defined(_WIN32) && !defined(__GNUC__)
#define STDCALL __stdcall
#else
#ifdef STDCALL
#undef STDCALL
#endif
#define STDCALL
#endif

#define __gks 1

#include "gksdefs.h"
#include "system.h"


#define Bool	int
#define True    1
#define False   0
#define Nil 	0

#define odd(status) ((status) & 01)


typedef struct {
    char *name;
    unsigned number;
} fac_struct;

typedef struct {
    char *name;
    unsigned code;
    char *text;
} msg_struct;


char *gks_a_error_info = Nil;
int (*gks_a_fprintf_routine)(FILE *, char *) = Nil;


static fac_struct facility[] = {

{"GKS",		0x0000085A}
};

static int n_facilities = sizeof(facility)/sizeof(facility[0]);

static char *severity_code[] = {
    "W", "S", "E", "I", "F", "?", "?", "?"
};

static msg_struct message[] = {

{"NORMAL",	0x085A8001, "normal successful completion"},
{"ERROR_1",	0x085A800A, "GKS not in proper state. GKS must be in the state\
 GKCL in routine %s"},
{"ERROR_2",	0x085A8012, "GKS not in proper state. GKS must be in the state\
 GKOP in routine %s"},
{"ERROR_3",	0x085A801A, "GKS not in proper state. GKS must be in the state\
 WSAC in routine %s"},
{"ERROR_4",	0x085A8022, "GKS not in proper state. GKS must be in the state\
 SGOP in routine %s"},
{"ERROR_5",	0x085A802A, "GKS not in proper state. GKS must be either in the\
 state WSAC or SGOP in routine %s"},
{"ERROR_6",	0x085A8032, "GKS not in proper state. GKS must be either in the\
 state WSOP or WSAC in routine %s"},
{"ERROR_7",	0x085A803A, "GKS not in proper state. GKS must be in one of the\
 states WSOP,WSAC,SGOP in routine %s"},
{"ERROR_8",	0x085A8042, "GKS not in proper state. GKS must be in one of the\
 states GKOP,WSOP,WSAC,SGOP in routine %s"},
{"ERROR_20",	0x085A80A2, "Specified workstation identifier is invalid in\
 routine %s"},
{"ERROR_21",	0x085A80AA, "Specified connection identifier is invalid in\
 routine %s"},
{"ERROR_22",	0x085A80B2, "Specified workstation type is invalid in routine\
 %s"},
{"ERROR_24",	0x085A80C2, "Specified workstation is open in routine %s"},
{"ERROR_25",	0x085A80CA, "Specified workstation is not open in routine %s"},
{"ERROR_26",	0x085A80D2, "Specified workstation cannot be opened in routine\
 %s"},
{"ERROR_27",	0x085A80DA, "Workstation Independent Segment Storage is not\
 open in routine %s"},
{"ERROR_28",	0x085A80E2, "Workstation Independent Segment Storage is already\
 open in routine %s"},
{"ERROR_29",	0x085A80EA, "Specified workstation is active in routine %s"},
{"ERROR_30",	0x085A80F2, "Specified workstation is not active in routine\
 %s"},
{"ERROR_34",	0x085A8110, "Specified workstation is not of category MI in\
 routine %s"},
{"ERROR_50",	0x085A8192, "Transformation number is invalid in routine %s"},
{"ERROR_51",	0x085A819A, "Rectangle definition is invalid in routine %s"},
{"ERROR_52",	0x085A81A2, "Viewport is not within the NDC unit square in\
 routine %s"},
{"ERROR_53",	0x085A81AA, "Workstation window is not within the NDC unit\
 square in routine %s"},
{"ERROR_60",	0x085A81E2, "Polyline index is invalid in routine %s"},
{"ERROR_62",	0x085A81F2, "Linetype is invalid in routine %s"},
{"ERROR_64",	0x085A8202, "Polymarker index is invalid in routine %s"},
{"ERROR_65",	0x085A820A, "Colour index is invalid in routine %s"},
{"ERROR_66",	0x085A8212, "Marker type is invalid in routine %s"},
{"ERROR_68",	0x085A8222, "Text index is invalid in routine %s"},
{"ERROR_70",	0x085A8232, "Text font is invalid in routine %s"},
{"ERROR_72",	0x085A8242, "Character expansion factor is invalid in routine\
 %s"},
{"ERROR_73",	0x085A824A, "Character height is invalid in routine %s"},
{"ERROR_74",	0x085A8252, "Character up vector is invalid in routine %s"},
{"ERROR_75",	0x085A825A, "Fill area index is invalid in routine %s"},
{"ERROR_78",	0x085A8272, "Style index is invalid in routine %s"},
{"ERROR_81",	0x085A828A, "Pattern size value is invalid in routine %s"},
{"ERROR_84",	0x085A82A2, "Dimensions of colour index array are invalid in\
 routine %s"},
{"ERROR_85",	0x085A82AA, "Colour index is invalid in routine %s"},
{"ERROR_88",	0x085A82C2, "Colour is invalid in routine %s"},
{"ERROR_100",	0x085A8322, "Number of points is invalid in routine %s"},
{"ERROR_161",	0x085A850A, "Item length is invalid in routine %s"},
{"ERROR_163",	0x085A851A, "Metafile item is invalid in routine %s"},
{"ERROR_164",	0x085A8522, "Item type is not a valid GKS item in routine %s"},
{"DEBUG",	0x085A95E3, "stepped to routine %s"},
{"ERROR_901",	0x085A9C2A, "Invalid connection identifier in routine %s"},
{"ERROR_902",	0x085A9C32, "Invalid workstation type in routine %s"},
{"ERROR_903",	0x085A9C3A, "A workstation of this type is already open in\
 routine %s"},
{"ERROR_904",	0x085A9C42, "No logical unit number available in routine %s"},
{"ERROR_905",	0x085A9C4A, "Open failed in routine %s"},
{"ERROR_906",	0x085A9C52, "Assign failed in routine %s"},
{"ERROR_907",	0x085A9C5A, "Deassign failed in routine %s"},
{"ERROR_908",	0x085A9C62, "VWS not accessible in routine %s"},
{"ERROR_909",	0x085A9C6A, "X display not accessible in routine %s"}
};

static int n_messages = sizeof(message)/sizeof(message[0]);

static int fontfile = -1;

static int hash_table[95], bufcache[95][256];

static enum {
    ANY, GLIGKS, GKSGRAL, GKGKS
} gks;

static int gks_errfile = 0;

#ifdef VMS
static int chan = 0;
#endif

#ifdef VAXC
noshare int gks_errno = 0, gks_status = -1;
#else
int gks_errno = 0, gks_status = -1;
#endif


void gks_fprintf (FILE *file, char *format, ...)
{
    va_list args;
    char s[BUFSIZ];

    va_start (args, format);
    vsprintf (s, format, args);
    va_end (args);

    if (gks_a_fprintf_routine != Nil)
        (*gks_a_fprintf_routine) (file, s);
    else
	fprintf (file, s);
}


void STDCALL OPFONT (int *gr_bc)
{
    char *path, fontdb[80];
    FILE *file;
    int i, nitems, pattern, pa[33], key;

    if ((char *) getenv ("GR_BC") != NULL)
	*gr_bc = 1;
    else
	*gr_bc = 0;

    path = (char *) getenv ("GLI_HOME");
    if (path)
    {
	strcpy (fontdb, path);
#ifndef VMS
#ifdef _WIN32
	strcat (fontdb, "\\");
#else
	strcat (fontdb, "/");
#endif
#endif
	strcat (fontdb, "gksfont.dat");
    }
    else
#ifdef VMS
	strcpy (fontdb, "sys$sysdevice:[gli]gksfont.dat");
#else
#ifdef _WIN32
	strcpy (fontdb, "c:\\gli\\gksfont.dat");
#else
	strcpy (fontdb, "/usr/local/gli/gksfont.dat");
#endif
#endif

    fontfile = open (fontdb,
#ifdef _WIN32
	O_RDONLY | O_BINARY, 0);
#else
	O_RDONLY, 0);
#endif
    if (fontfile == -1)
#ifdef VMS
        fontfile = open ("sys$library:gksfont.dat", O_RDONLY, 0);
#else
#ifdef _WIN32
        fontfile = open ("c:\\gli\\gksfont.dat", O_RDONLY | O_BINARY, 0);
#else
        fontfile = open ("/usr/local/lib/gksfont.dat", O_RDONLY, 0);
#endif
#endif
    gks = ANY;

    path = (char *) getenv ("GLI_GKS_PATTERN");
    if (path)
    {
        file = fopen (path, "r");
        if (file)
        {
            while (!feof(file))
            {
		nitems = 0;
		*pa = -1;

                if (fscanf (file, "%d %d", &pattern, pa) == 2)
		{
		    if (*pa == 4 || *pa == 8 || *pa == 32)
		    {
			for (i = 1; i <= *pa; i++)
			{
			    if (fscanf (file, "%d", pa + i) != 1)
				break;
			    nitems++;
			}
			fscanf (file, "\n");
		    }
                }
		if (nitems != *pa)
		{
                    gks_fprintf (stderr,
			"GKS (gksio): can't read pattern file\n");
                    break;
                }
                else
                    GKSPA (&pattern, pa);
            }            
            fclose (file);
        }
	else
            gks_fprintf (stderr,
		"GKS (gksio): can't access pattern database\n");
    }

    for (key = 0; key < 95; key++)
        hash_table[key] = -1;

    gks_status = -1;
}


void STDCALL LOOKUP (int *font, int *ch, int *buffer)
{
   /*  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 */
    static int map[] = {
       1,18, 1, 6,12, 3, 8,11, 4, 7,10, 2,13,14, 5, 9,15,16,17,20,21,19,22,23 };
    static int gksgralmap[] = {
       1,12, 6, 9, 8,11, 5,13,18,17,19, 1, 4, 7,24, 1, 1, 1, 1, 1, 1, 1,23,24 };
    static int gkgksmap[] = {
       1, 2, 3, 4, 5, 6, 7,13,11,14,12,15,16,17,18,19,20,21, 1, 1, 1, 1,23,24 };
    static int s_map[] = {
       4, 4, 4, 4, 4, 7, 7, 7,10,10,10, 7, 7, 7, 4, 4, 7, 7, 7, 4, 4, 4, 4, 4 };

    static int german[] = {
	196, 214, 220, 228, 246, 252, 223, 171, 187, 183, 169};
    static char ansi[] = {
	'A', 'O', 'U', 'a', 'o', 'u', 'b', '<', '>', '.', '@'};
    static char greek[] = {
	'j', 'o', 'q', 'u', 'v', 'w', 'y', 'J', 'O', 'Q', 'U', 'V', 'W', 'Y' };
    static char g_map[] = {
	' ', 'w', ' ', 'o', 'y', 'v', 'q', ' ', 'W', ' ', 'O', 'Y', 'V', 'Q' };
 
    char *env, buf[256];
    int fn, chr, umlaut, sharp_s, offset;
    register int i, *elptr;
    register char *ebptr;
    
    if (gks == ANY)
    {
        gks = GLIGKS;
        if (env = (char *) getenv ("GLI_GKS"))
        {
	    if (strcmp (env, "GKSGRAL") == 0)
                gks = GKSGRAL;
            else if (strcmp (env, "GKGKS") == 0)
                gks = GKGKS;
        }
    }

    if (fontfile != -1)
    {
	chr = *ch;
	umlaut = sharp_s = False;

	if (chr >= 127)
        {
	    for (i=0; i<=10; i++)
            {
		if (chr == german[i]) {
		    chr = ansi[i];
		    if (i < 6)
			umlaut = True;
		    else if (i == 6)
			sharp_s = True;
		}
            }
        }
	if (chr < ' ' || chr >= 127) chr = ' ';

	fn = abs(*font) % 100;
        if (fn == 51)
            fn = 23;            /* fill font */
        else if (fn > 23)
            fn = 1;

	if (chr == '_' ) {
	    if (fn < 20)
		fn = 23;
	} else if (sharp_s) {
            if (fn != 23)
                fn = s_map[fn-1];
            else
                chr = 126;	/* ~ */
        }
	else if (gks == GKSGRAL)
        {
            if (fn == 13 || fn == 14)
            {
                for (i=0; i<14; i++) {
                    if (chr == greek[i]) {
                        chr = g_map[i];
                        break;
                    }
                }
            }
            fn = gksgralmap[fn-1];
        }
        else if (gks == GKGKS)
            fn = gkgksmap[fn-1];

        chr -= ' ';
	offset = ((map[fn-1]-1)*95 + chr)*256;

        if (hash_table[chr] != offset)
        {
            if (lseek (fontfile, offset, 0) != -1)
            {
                if (read (fontfile, buf, 256) != -1)
                {
                    hash_table[chr] = offset;

                    elptr = bufcache[chr];
                    ebptr = buf;
                    for (i = 0; i < 256; i++)
                        *elptr++ = *ebptr++;
                }
                else
                    gks_fprintf (stderr, "GKS (gksio): file read error\n");
            }
            else
                gks_fprintf (stderr, "GKS (gksio): file position error\n");
        }
        memcpy ((void *) buffer, (void *) bufcache[chr], 256 * sizeof(int));

        if (umlaut && (buffer[7] < 120-20))
            buffer[7] = buffer[7] + 10;
    }
    else
    {
	gks_fprintf (stderr, "GKS (gksio): can't access font database\n");
	exit (-1);
    }
}


void STDCALL CLFONT (void)
{
#ifdef VMS
    int ret_status;
#endif

    if (fontfile != -1)
	close (fontfile);

#ifdef VMS
    if (chan != 0)
    {
	ret_status = SYS$DASSGN (chan);
	if (odd(ret_status))
	    chan = 0;
	else
	    LIB$STOP (ret_status);
    }
#endif
}


void STDCALL BUFIN (
    int *conid, int *status, int *terminator, int *nchars, CHARARG(buffer))
{
    char *cp;
#if defined (MSDOS) || defined (_WIN32)
    int count;
#else
#ifdef VMS
    $DESCRIPTOR (terminal, "SYS$COMMAND");
    int ret_status, mask[] = {0, 0};
    short iosb[4];
#else
#ifdef NO_TERMIO
    struct sgttyb saved_term, term;
    int saved_term2, term2;
#else
    struct termios saved_term, term;
#endif /* NO_TERMIO */
#endif /* VMS */
#endif /* MSDOS, _WIN32 */

#ifdef VMS
    cp = buffer->dsc$a_pointer;
#else
#ifdef cray
    cp = _fcdtocp(buffer);
#else
    cp = buffer;
#endif /* cray */
#endif /* VMS */

    if (*conid == 0)
    {
#if defined (MSDOS) || defined (_WIN32)
        *status = 1;
	for (count = 0; count < *nchars; count++)
	    cp[count] = getchar ();
#else
#ifdef VMS
	if (chan == 0)
	{
	    ret_status = SYS$ASSIGN (&terminal, &chan, 0, 0);
	    if (!odd(ret_status))
		LIB$STOP (ret_status);
	}

	mask[1] = *terminator;
	ret_status = SYS$QIOW (0, chan, IO$_READVBLK | IO$M_NOECHO |
	    IO$M_TRMNOECHO | IO$M_NOFILTR, iosb, 0, 0, cp, *nchars, 0, mask,
	    0, 0);
	if (!odd(ret_status))
	    LIB$STOP (ret_status);

	*status = odd(iosb[0]) ? 1 : 0;
#else
#ifdef NO_TERMIO
        ioctl (0, TIOCGETP, &term);
        ioctl (0, TIOCLGET, &term2);
        saved_term = term;
        saved_term2 = term2;

        /* Turn off CR mapping and echo, turn on full raw mode */

        term.sg_flags &= ~(CRMOD | ECHO);
        term.sg_flags |= (CBREAK | RAW);

        ioctl (0, TIOCSETN, &term);
        *status = (read (0, cp, *nchars) == *nchars) ? 1 : 0;

        ioctl (0, TIOCSETN, &saved_term);
        ioctl (0, TIOCLSET, &saved_term2);
#else
	if (tcgetattr (0, &term) < 0)
            perror ("tcgetattr");	
        saved_term = term;

        term.c_iflag &= ~(IGNCR | INLCR | ICRNL);
        term.c_lflag &= ~(ICANON | ISIG | ECHO);
        term.c_cc[VMIN] = *nchars;
        term.c_cc[VTIME] = 0;

	if (tcsetattr (0, TCSADRAIN, &term) < 0)
            perror ("tcsetattr");
        *status = (read (0, cp, *nchars) == *nchars) ? 1 : 0;

	if (tcsetattr (0, TCSADRAIN, &saved_term) < 0) 
            perror ("tcsetattr");
#endif /* NO_TERMIO */
#endif /* VMS */
#endif /* MSDOS, _WIN32 */
    }
    else
        *status = 0; /* GNONE */
}


void STDCALL BUFOUT (int *conid, int *nchars, CHARARG(chars))
{
    char *cp;
    char buf[256];
    static char *err = NULL;
    int fd, n;
#ifdef VMS
    int status = -1;

    cp = chars->dsc$a_pointer;
#else
#ifdef cray
    cp = _fcdtocp(chars);
#else
    cp = chars;
#endif /* cray */
#endif /* VMS */

    if (*nchars > 0)
    {
	if (*conid != 0)
	{
	    if (*conid > 100)
            {
                fd = *conid - 100;
                n = write (fd, cp, *nchars);
            }
            else
            {
#if defined(_WIN32) && !defined(__GNUC__)
                unsigned short chars_len = *nchars;

		FORTWR (conid, nchars, chars, chars_len, &n);
#else
		FORTWR (conid, nchars, chars, &n, *nchars);
#endif
            }
	}
	else
	    n = write (1, cp, *nchars);

	if (n < 0 && !err)
	{
	    if (n != -1)
#ifdef VMS
	    {
		sprintf (buf, "FORTRAN I/O error number %d", n);
		err = buf;
	    }
#else
	    {
		sprintf (buf, "Unknown I/O error (%d)", n);
		err = buf;
	    }
#endif

	    else if (!(err = (char *) strerror (errno)))
#ifdef VMS
	    {
		err = "VMS message:";
		status = vaxc$errno;
	    }
#else
	    {
		sprintf (buf, "Unknown error (%d)", errno);
		err = buf;
	    }
#endif
	    gks_fprintf (stderr, "GKS (gksio): file write error. %s\n", err);
#ifdef VMS
	    if (status != -1)
		lib$signal (status);
#endif
	    exit (-1);
	}
    }
}


void STDCALL BINOUT (int *conid, int *nchars, CHARARG(chars))
{
    BUFOUT (conid, nchars, CHARPAR(chars));
}


#ifdef DPS

static FILE *stream = NULL;
static int fd = 0;

#endif


void STDCALL DPSOP (float *sizex, float *sizey, int *format)
{
#ifdef DPS
    char *path, *delim, command[255];
    Display *dpy;
    Screen *screen;
    int dpi, width, height;

    if ((dpy = XOpenDisplay(NULL)) != NULL)
    {
	screen = XDefaultScreenOfDisplay(dpy);
	dpi = (int)(XWidthOfScreen(screen) / (XWidthMMOfScreen(screen) / 25.4));
	XCloseDisplay(dpy);
    }
    else
    {
	gks_fprintf(stderr, "GKS: can't open display\n");
	exit (-1);
    }

    path = (char *) getenv ("GLI_HOME");
    if (path == NULL)
#ifdef VMS
	path = "sys$sysdevice:[gli]";
#else
#ifdef _WIN32
	path = "c:\\gli";
#else
	path = "/usr/local/gli";
#endif
#endif

    width = (int)(*sizex * dpi / 600);
    height = (int)(*sizey * dpi / 600);
#ifdef VMS
    delim = "";
#else
#ifdef _WIN32
    delim = "\\";
#else
    delim = "/";
#endif
#endif
    sprintf(command, "%s%spsprint -width %d -height %d %s",
	path, delim, width, height, options[*format]);

    if ((stream = popen(command, "w")) == NULL)
    {
	gks_fprintf(stderr, "GKS: can't initiate pipe to DPS system\n");
	exit (-1);
    }
    else
	fd = fileno(stream);
#else
    gks_fprintf(stderr,
	"GKS: X client does not have Display PostScript system\n");
#endif /* DPS */
}


void STDCALL DPSWR (int *nchars, CHARARG(chars))
{
#ifdef DPS
    char *cp;

#ifdef VMS
    cp = chars->dsc$a_pointer;
#else
#ifdef cray
    cp = _fcdtocp(chars);
#else
    cp = chars;
#endif /* cray */
#endif /* VMS */

    if (*nchars > 0)
    {
	if (write(fd, cp, *nchars) != *nchars) {
	    gks_fprintf(stderr, "GKS: can't write to DPS pipe\n");
	    exit (-1);
	}
    }
#endif /* DPS */
}


void STDCALL DPSFL (void)
{
#ifdef DPS
    if (stream)
	fflush(stream);
#endif /* DPS */
}


void STDCALL DPSCL (void)
{
#ifdef DPS
    if (stream) {
	fflush(stream);
	pclose(stream);
	stream = NULL;
    }
#endif /* DPS */
}


void STDCALL GKINFO (int *nchars, CHARARG(chars))
{
    char *cp, *date, *user, host[100];
    time_t elapsed_time;
#ifdef hpux
    struct utsname utsname;
#endif
#ifdef _WIN32
    char lpBuffer[100];
    DWORD nSize = 100;
#endif

#ifdef VMS
    cp = chars->dsc$a_pointer;
#else
#ifdef cray
    cp = _fcdtocp(chars);
#else
    cp = chars;
#endif /* cray */
#endif /* VMS */

    time (&elapsed_time);
    date = ctime (&elapsed_time);

#ifndef _WIN32
    user = (char *)getenv("USER");
#else
    if (GetUserName(lpBuffer, &nSize) != 0)
    {
        user = lpBuffer;
        lpBuffer[nSize] = '\0';
    }
    else
        user = NULL;
#endif
    if (user == NULL)
	user = "(?)";

#ifdef VMS
    strcpy (host, (char *) getenv ("SYS$NODE"));
#else
#ifdef hpux
    uname (&utsname);
    strcpy (host, utsname.nodename);
#else
#if defined(OS2) || (defined(_WIN32) && defined(__GNUC__))
    strcpy (host, "(unknown)"); /* FIXME */
#else
    gethostname (host, 100);
#endif /* _WIN32 */
#endif /* hpux */
#endif /* VMS */

    strtok (date, "\n");
    strtok (host, ".");
	
    sprintf (cp, "%s  by user  %s @ %s", date, user, host);
    *nchars = strlen(cp);
}


static
char *base_name(char *name)
{
    char *base = name;

    while (*name)
    {
        if (*name++ == '/')
            base = name;
    }
    return (base);
}


void STDCALL GKTMP (int *nchars, CHARARG(chars))
{
    char *cp, s[100], *env;

#ifdef VMS
    cp = chars->dsc$a_pointer;
#else
#ifdef cray
    cp = _fcdtocp(chars);
#else
    cp = chars;
#endif /* cray */
#endif /* VMS */

#ifdef __NetBSD__
    sprintf(s, "/tmp/gks%d.tmp\n", getpid());
#else
    tmpnam(s);
#endif
    if ((env = (char *) getenv("GLI_GKS_TMP")) != NULL)
	sprintf(cp, "%s/%s", env, base_name(s));
    else
	strcpy (cp, s);

    *nchars = strlen(cp);
}


void STDCALL GKMAGS (float *magstep, int *dpi)
{
    char *env;
#ifdef DPS
    Display *dpy;
    Screen *screen;
#endif /* DPS */

    if (env = (char *) getenv ("GLI_GKS_MAGSTEP"))
        *magstep = atof(env);
    else
        *magstep = 0;

    *dpi = 75;
#ifdef DPS
    if ((dpy = XOpenDisplay(NULL)) != NULL)
    {
	screen = XDefaultScreenOfDisplay(dpy);
	*dpi = (int)(XWidthOfScreen(screen) /
	    (XWidthMMOfScreen(screen) / 25.4));
	XCloseDisplay(dpy);
    }
#endif /* DPS */
}


static int gks_get_status_text (int msgid, unsigned int flags, char *bufadr)
{
    static char fac[32], sev[2], id[32], text[256];
    int i;

    strcpy (fac, "NONAME");
    for (i = 0; i < n_facilities; i++)
	if (STATUS_FAC_NO (msgid) == facility[i].number) {
	    strcpy (fac, facility[i].name);
	    break;
	}
 
    strcpy (sev, severity_code [STATUS_SEVERITY (msgid)]);
    strcpy (id, "NOMSG");
    strcpy (text, "");

    for (i = 0; i < n_messages; i++)
	if (STATUS_COND_ID (msgid) == STATUS_COND_ID (message[i].code)) {
	    strcpy (id, message[i].name);
	    strcpy (text, message[i].text);
	    break;
	}

    if (strlen(text) == 0)
	sprintf (text, "Message number %x", msgid);

    if (flags == 0)
	flags = STS_M_MSG_FAC | STS_M_MSG_SEV | STS_M_MSG_ID | STS_M_MSG_TEXT;

    strcpy (bufadr, "");
    if (flags & STS_M_MSG_FAC) {
	strcpy (bufadr, "%%");
	strcat (bufadr, fac);
    }

    if (flags & STS_M_MSG_SEV) {
	if (strlen(bufadr) == 0)
	    strcpy (bufadr, "%%");
	else
	    strcat (bufadr, "-");
	strcat (bufadr, sev);
    }

    if (flags & STS_M_MSG_ID) {
	if (strlen(bufadr) == 0)
	    strcpy (bufadr, "%%");
	else
	    strcat (bufadr, "-");
	strcat (bufadr, id);
    }

    if (flags & STS_M_MSG_TEXT) {
	if (strlen(bufadr) != 0) strcat (bufadr, ", ");
	strcat (bufadr, text);
    }

    return SS__NORMAL;
}


void STDCALL LIBSIG (int *status, CHARARG(arg))

/*
 *  LIBSIG - signal exception condition
 */

{
    char message[255], chars[255], name[7];
#ifdef VMS
    struct dsc$descriptor_s text;
#endif
#ifdef cray
    _fcd text;
#endif
    int nchars;
#if defined(_WIN32) && !defined(__GNUC__)
    int chars_len = arg_len;
#endif

#ifdef VMS
    strncpy(name, arg->dsc$a_pointer, 6);
#else
#ifdef cray
    strncpy(name, _fcdtocp(arg), 6);
#else
    strncpy(name, arg, 6);
#endif /* cray */
#endif /* VMS */
    name[6] = '\0';

    if (*status != gks_status)
    {
	if (gks_a_error_info == Nil)
	{
	    gks_get_status_text (*status, 0, message);
	    sprintf (chars, message, name);

	    strcat (chars, "\r\n");
	    nchars = strlen(chars);
#ifdef VMS
	    text.dsc$b_dtype = DSC$K_DTYPE_T;
	    text.dsc$b_class = DSC$K_CLASS_S;
	    text.dsc$w_length = nchars;
	    text.dsc$a_pointer = chars;

	    BUFOUT (&gks_errfile, &nchars, &text);
#else
#ifdef cray
	    text = _cptofcd(chars, nchars);

	    BUFOUT (&gks_errfile, &nchars, text);
#else
	    BUFOUT (&gks_errfile, &nchars, CHARPAR(chars));
#endif /* cray */
#endif /* VMS */
	}
	else
	{
	    gks_get_status_text (*status, STS_M_MSG_TEXT, message);
	    sprintf (gks_a_error_info, message, name);
	}

	gks_status = *status;
    }
}


void STDCALL GKFD (int *conid)

/*
 *  GKFD - check for a valid file descriptor
 */

{
#if !defined (VMS)
    struct stat buf;
#else
    char buf[256];
#endif

    if (*conid != 0)
    {
#if !defined (VMS)
        if (isatty(*conid) == 1 || (fstat (*conid, &buf) == 0))
#else
        if (isatty(*conid) == 1 || getname(*conid, &buf) > 0)
#endif /* VMS */
            *conid += 100;
    }
}


void STDCALL GTWSTY (int *wstype)

/*
 *  GTWSTY - get workstation type
 */

{
    char *env;

    env = (char *) getenv ("GLI_WSTYPE");
#ifdef VMS
    if (!env)
	env = (char *) getenv ("GKS3D$WSTYPE");
    if (!env)
	env = (char *) getenv ("GKS$WSTYPE");
#endif
#if defined (ultrix) || defined(__osf__)
    if (!env)
	env = (char *) getenv ("GKS3Dwstype");
    if (!env)
	env = (char *) getenv ("GKSwstype");
#endif

    if (env)
	*wstype = atoi(env);
    else
#ifndef _WIN32
	*wstype = 211;
#else
	*wstype = 41;
#endif
}


void STDCALL GKLTOI (unsigned long *l, unsigned int *i)
{
    unsigned int *p;

    if (sizeof(int) < sizeof(long))
    {
	p      = (unsigned int *) l;
	*i     = *p;
	*(i+1) = *(p+1);
    }
}


void STDCALL GERSET (int *errnum, int *errfil)
{
    gks_errno = *errnum;
    if (*errfil >= 0)
        gks_errfile = *errfil;
}


#ifdef ultrix
int s_abort ()
{
    return 0;
}
#endif /* ultrix */


#if defined(linux) || defined(F2C)
unsigned int iand_ (unsigned int *a, unsigned int *b)
{
    return *a & *b;
}

unsigned int ior_ (unsigned int *a, unsigned int *b)
{
    return *a | *b;
}

#endif /* linux, F2C */


void STDCALL GKOPEN (int *conid, CHARARG(name))
{
    char *path;

#ifdef VMS
    path = name->dsc$a_pointer;
#else
#ifdef cray
    path = _fcdtocp(name);
#else
    path = name;
#endif /* cray */
#endif /* VMS */

    *conid = open(path, O_CREAT | O_TRUNC | O_WRONLY, 0644);
    if (*conid == -1)
        gks_fprintf (stderr, "GKS (gksio): can't open file %s\n", path);
}


void STDCALL GKCLOS (int *conid)
{
    close(*conid);
}
