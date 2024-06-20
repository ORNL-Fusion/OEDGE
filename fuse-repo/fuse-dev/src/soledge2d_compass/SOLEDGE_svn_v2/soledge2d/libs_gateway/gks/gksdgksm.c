/*
 * Copyright @ 1995 - 1999   Josef Heinen
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
 *	This module contains a logical device driver for GKSM metafiles.
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


#if defined(RPC) || (defined(_WIN32) && !defined(__GNUC__))
#define HAVE_SOCKETS
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#if !defined(VMS) && !defined(MSDOS) && !defined(_WIN32)
#include <unistd.h>
#endif

#if defined (cray) || defined (__SVR4) || defined (MSDOS) || defined(_WIN32)
#include <fcntl.h>
#else
#include <sys/types.h>
#include <sys/file.h>
#endif

#ifdef HAVE_SOCKETS
#ifndef _WIN32
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#else
#include <windows.h>
#endif
#endif

#ifndef MSDOS
#include <sys/stat.h>
#endif

#ifdef VMS
#include <descrip.h>
#endif

#ifdef cray
#include <fortran.h>
#endif


#include "gksdefs.h"


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

#ifndef STDCALL
#if defined(_WIN32) && !defined(__GNUC__)
#define STDCALL __stdcall
#else
#define STDCALL
#endif
#endif

#ifndef MSDOS
#define SEGM_SIZE       262144			/* 256K */
#else
#define SEGM_SIZE       32767			/* 32K */
#endif

#define COPY(s, n) memcpy(p->buffer + p->nbytes, (void *) s, n); p->nbytes += n

#define RESOLVE(arg, type, nbytes) arg = (type *)(s + sp); sp += nbytes


typedef struct ws_state_list_struct {
    int conid, state;
    int socket;
    int empty;
    char *buffer;
    int size, nbytes, position;
} ws_state_list;


static ws_state_list *p;
static gks_state_list *gksl;
static int wkid = 1;


static
void reallocate(int len)
{
    while (p->nbytes + len > p->size)
	p->size += SEGM_SIZE;

    p->buffer = (char *) realloc(p->buffer, p->size + 1);
}

static
void write_item(int *fctid, int *dx, int *dy, int *dimx, int *ia,
    int *lr1, float *r1, int *lr2, float *r2,
    int *lc,
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
    char s[132];
    int len = -1, slen;

    switch (*fctid)
    {
      case 12 :		/* polyline */
      case 13 :		/* polymarker */
      case 15 :		/* fill area */

	len = 3 * sizeof(int) + 2 * *ia * sizeof(float);
	if (p->nbytes + len > p->size)
	    reallocate(len);

	COPY(&len, sizeof(int));
	COPY(fctid, sizeof(int));
	COPY(ia, sizeof(int));
	COPY(r1, *ia * sizeof(float));
	COPY(r2, *ia * sizeof(float));
	break;

      case 14 :		/* text */

	len = 3 * sizeof(int) + 2 * sizeof(float) + 132;
	if (p->nbytes + len > p->size)
	    reallocate(len);

	memset((void *) s, 0, 132);
#ifdef VMS
	slen = chars->dsc$w_length;
	strncpy(s, chars->dsc$a_pointer, slen);
#else
#ifdef cray
	slen = _fcdlen(chars);
	strncpy(s, _fcdtocp(chars), slen);
#else
	slen = strlen(chars);
	strncpy(s, chars, slen);
#endif /* cray */
#endif /* VMS */

	COPY(&len, sizeof(int));
	COPY(fctid, sizeof(int));
	COPY(r1, sizeof(float));
	COPY(r2, sizeof(float));
	COPY(&slen, sizeof(int));
	COPY(s, 132);
	break;

      case 16 :		/* cell array */

	len = (5 + *dimx * *dy) * sizeof(int) + 4 * sizeof(float);
	if (p->nbytes + len > p->size)
	    reallocate(len);

	COPY(&len, sizeof(int));
	COPY(fctid, sizeof(int));
	COPY(r1, 2 * sizeof(float));
	COPY(r2, 2 * sizeof(float));
	COPY(dx, sizeof(int));
	COPY(dy, sizeof(int));
	COPY(dimx, sizeof(int));
	COPY(ia, *dimx * *dy * sizeof(int));
	break;

      case 19 :		/* set linetype */
      case 21 :		/* set polyline color index */
      case 23 :		/* set markertype */
      case 25 :		/* set polymarker color index */
      case 30 :		/* set text color index */
      case 33 :		/* set text path */
      case 36 :		/* set fillarea interior style */
      case 37 :		/* set fillarea style index */
      case 38 :		/* set fillarea color index */
      case 52 :		/* select normalization transformation */
      case 53 :		/* set clipping indicator */

	len = 3 * sizeof(int);
	if (p->nbytes + len > p->size)
	    reallocate(len);

	COPY(&len, sizeof(int));
	COPY(fctid, sizeof(int));
	COPY(ia, sizeof(int));
	break;

      case 27 :		/* set text font and precision */
      case 34 :		/* set text alignment */

	len = 4 * sizeof(int);
	if (p->nbytes + len > p->size)
	    reallocate(len);

	COPY(&len, sizeof(int));
	COPY(fctid, sizeof(int));
	COPY(ia, 2 * sizeof(int));
	break;

      case 20 :		/* set linewidth scale factor */
      case 24 :		/* set marker size scale factor */
      case 28 :		/* set character expansion factor */
      case 29 :		/* set character spacing */
      case 31 :		/* set character height */

	len = 2 * sizeof(int) + sizeof(float);
	if (p->nbytes + len > p->size)
	    reallocate(len);

	COPY(&len, sizeof(int));
	COPY(fctid, sizeof(int));
	COPY(r1, sizeof(float));
	break;

      case 32 :		/* set character up vector */

	len = 2 * sizeof(int) + 2 * sizeof(float);
	if (p->nbytes + len > p->size)
	    reallocate(len);

	COPY(&len, sizeof(int));
	COPY(fctid, sizeof(int));
	COPY(r1, sizeof(float));
	COPY(r2, sizeof(float));
	break;

      case 48 :		/* set color representation */

	len = 3 * sizeof(int) + 3 * sizeof(float);
	if (p->nbytes + len > p->size)
	    reallocate(len);

	COPY(&len, sizeof(int));
	COPY(fctid, sizeof(int));
	COPY(&ia[1], sizeof(int));
	COPY(r1, 3 * sizeof(float));
	break;

      case 49 :		/* set window */
      case 50 :		/* set viewport */

	len = 3 * sizeof(int) + 4 * sizeof(float);
	if (p->nbytes + len > p->size)
	    reallocate(len);

	COPY(&len, sizeof(int));
	COPY(fctid, sizeof(int));
	COPY(ia, sizeof(int));
	COPY(r1, 2 * sizeof(float));
	COPY(r2, 2 * sizeof(float));
	break;
    }
}

#ifdef HAVE_SOCKETS

static
int open_socket(char *name, char flag)
{
    static int s = -1;
    struct hostent *hp;
    struct sockaddr_in sin;
    char *hostname, *port, tmp[BUFSIZ];
    int sd;

#if defined(_WIN32) && !defined(__GNUC__)
    WORD wVersionRequested = MAKEWORD(2, 0);
    WSADATA wsaData;

    if (WSAStartup(wVersionRequested, &wsaData) != 0)
    {
        gks_fprintf(stderr, "Can't find a usable WinSock DLL\n");
        exit(1);
    }
#endif
    if (s == -1 || flag == 'w')
    {
	strcpy(tmp, name);
	hostname = strtok(tmp, ":");
	port = strtok(NULL, ":");
	if (port == NULL) {
	    port = hostname;
	    hostname = "localhost";
	}

	if (flag == 'w') {
	    if ((hp = gethostbyname(hostname)) == 0) {
		perror("gethostbyname");
		exit(1);
	    }
	}

	memset(&sin, 0, sizeof(sin));
	sin.sin_family = AF_INET;
	sin.sin_addr.s_addr = flag == 'w' ?
	    ((struct in_addr *)(hp->h_addr))->s_addr : INADDR_ANY;
	sin.sin_port = htons((unsigned short)atoi(port));
    }

    if (s == -1)
    {
	s = socket(AF_INET, SOCK_STREAM, 0);
	if (s == -1) {
	    perror("socket");
	    exit(1);
	}

	if (flag == 'r') {
	    if (bind(s, (struct sockaddr *)&sin, sizeof(sin)) == -1) {
		perror("bind");
		exit(1);
	    }
	    if (listen(s, 1) == -1) {
		perror("listen");
		exit(1);
	    }
	}
    }

    if (flag == 'r') {
	sd = accept(s, NULL, NULL);
	if (sd == -1) {
	    perror("accept");
	    exit(1);
	}
    }
    else {
	if (connect(s, (struct sockaddr *)&sin, sizeof(sin)) == -1) {
	    perror("connect");
	    exit(1);
	}
	sd = s;
    }

    return sd;
}

#endif

static
void write_gksm(int stream)
{
    char *gksm_name;
#ifndef MSDOS
    struct stat buf;
#endif
    int fd, oflag, nbytes;
    char *buffer;
    int opened = 0;

    if (p->socket == -1)
	fd = stream > 100 ? stream - 100 : stream;
    else
	fd = p->socket;

    buffer = p->buffer;
    nbytes = p->nbytes;
    if (p->position != 0)
    {
	buffer += p->position;
	nbytes -= p->position;
    }

#ifndef MSDOS
    if (isatty(fd) || fstat(fd, &buf) != 0)
#else
    if (isatty(fd))
#endif
    {
        if ((gksm_name = (char *) getenv("GLI_GKSM")) == NULL)
            gksm_name = "gli.gksm";

#ifdef HAVE_SOCKETS
	if (strchr(gksm_name, ':') == NULL)
	{
#endif
	    oflag = O_CREAT | O_WRONLY | (p->position != 0 ?
		O_APPEND : O_TRUNC);
#ifdef _WIN32
            oflag |= O_BINARY;
#endif
	    if ((fd = open(gksm_name, oflag, 0644)) < 0) {
		gks_fprintf(stderr, "GKS: can't open GKSM metafile\n");
		perror("open");
	    }
	    else
		opened = 1;
#ifdef HAVE_SOCKETS
	}
	else if (p->socket == -1)
	{
	    fd = open_socket(gksm_name, 'w');
	    if (fd >= 0)
		p->socket = fd;
	}
#endif
    }

    if (fd >= 0)
    {
	int offset = 0, bufsiz, cc;

	while (offset < nbytes)
	{
	    bufsiz = (nbytes - offset <= BUFSIZ) ? nbytes - offset : BUFSIZ;
	    if ((cc = 
#ifdef HAVE_SOCKETS
		p->socket != -1 ? send(fd, buffer + offset, bufsiz, 0) :
#endif
		write(fd, buffer + offset, bufsiz)) <= 0)
	    {
		gks_fprintf(stderr, "GKS: can't write GKSM metafile\n");
		perror(p->socket != -1 ? "send" : "write");
		break;
	    }
	    offset += cc;
	}
	if (opened && p->socket == -1)
	    close(fd);
    }
}

void STDCALL GKDMFO(
    int *fctid, int *dx, int *dy, int *dimx, int *ia, int *lr1, float *r1,
    int *lr2, float *r2, int *lc, CHARARG(chars), ws_state_list **ptr)
{
    int gksm = 2;
    p = *ptr;

    switch (*fctid)
    {
      case 2 :		/* open workstation */

	p = (ws_state_list *) malloc(sizeof(ws_state_list));

	p->conid = ia[1];
	p->state = GINACT;
	p->socket = -1;
	p->empty = 1;

	p->buffer = (char *) malloc(SEGM_SIZE + 1);
	p->size = SEGM_SIZE;
	p->nbytes = p->position = 0;

        gksl = (gks_state_list *) (ia + 4);

	if (sizeof(char *) > sizeof(int)) {
	    long *la = (long *) ia;
	    *la = (long) p;
	} else
	    *ia = (long) p;

	break;

      case 3 :		/* close workstation */

	if (p->position < p->nbytes && !p->empty)
	    write_gksm(p->conid);

	if (p->socket != -1)
#if defined (_WIN32) && !defined (__GNUC__)
	{
	    closesocket(p->socket);
	    WSACleanup();
	}
#else
	    close(p->socket);
#endif
	free(p->buffer);
	free(p);

	p = NULL;
	break;

      case 4 :		/* activate workstation */

	p->state = GACTIV;
	break;

      case 5 :		/* deactivate workstation */

	p->state = GINACT;
	break;

      case 6 :		/* clear workstation */

	p->nbytes = p->position = 0;
	p->empty = 1;
	break;

      case 8 :		/* update workstation */

	if (ia[1] == GPERFO)
	{
	    if (p->position < p->nbytes && !p->empty) {
		write_gksm(p->conid);
		p->position = p->nbytes;
	    }
	}
	break;

      case 12 : case 13 : case 14 : case 15 : case 16 :	p->empty = 0;
      case 19 : case 20 : case 21 : case 23 : case 24 : case 25 : case 27 :
      case 28 : case 29 : case 30 : case 31 : case 32 : case 33 : case 34 :
      case 36 : case 37 : case 38 : case 48 : case 49 : case 50 : case 52 :
      case 53 :

	if (p->state == GACTIV)
	{
	    int len = 2 * sizeof(int) + sizeof(gks_state_list);

	    if (p->nbytes == 0) {
		COPY(&len, sizeof(int));
		COPY(&gksm, sizeof(int));
		COPY(gksl, sizeof(gks_state_list));
	    }
	    write_item(fctid, dx, dy, dimx, ia, lr1, r1, lr2, r2, lc, chars);
	}
	break;
    }
}

static
char *readfile(int fd)
{
    char *gksm_name;
    int cc;
#ifndef MSDOS
    struct stat buf;
#endif
    char *s = NULL;
    int size;
#ifdef HAVE_SOCKETS
    int is_socket = 0;
#endif

    if (fd != -1)
    {
#ifndef MSDOS
        if (fstat(fd, &buf) != 0)
	{
	    if ((gksm_name = (char *) getenv("GLI_GKSM")) == NULL)
		gksm_name = "gli.gksm";

#ifdef HAVE_SOCKETS
	    if (strchr(gksm_name, ':') == NULL)
#endif
		fd = open(gksm_name,
#ifdef _WIN32
		    O_RDONLY | O_BINARY, 0);
#else
		    O_RDONLY, 0);
#endif
#ifdef HAVE_SOCKETS
	    else if ((fd = open_socket(gksm_name, 'r')) != -1)
		is_socket = 1;
#endif
	}
        fstat(fd, &buf);
	size = (buf.st_size > 0) ? buf.st_size : 1000000;
#else
        size = 32767;
#endif
	s = (char *) malloc(size + 1);
	if ((cc =
#ifdef HAVE_SOCKETS
	    is_socket ? recv(fd, s, size, 0) :
#endif
	    read(fd, s, size)) != -1)
	    s[cc] = '\0';
    }
    else
        gks_fprintf(stderr, "GKS: invalid file descriptor (%d)\n", fd);

    return s;
}

static
void gksinit(gks_state_list *gksl)
{
    int i, tnr;
    float *wn, *vp;

    GSPLI(&gksl->lindex);
    GSLN(&gksl->ltype);
    GSLWSC(&gksl->lwidth);
    GSPLCI(&gksl->plcoli);

    GSPMI(&gksl->mindex);
    GSMK(&gksl->mtype);
    GSMKSC(&gksl->mszsc);
    GSPMCI(&gksl->pmcoli);

    GSTXI(&gksl->tindex);
    GSTXFP(&gksl->txfont, &gksl->txprec);
    GSCHXP(&gksl->chxp);
    GSCHSP(&gksl->chsp);
    GSTXCI(&gksl->txcoli);
    GSCHH(&gksl->chh);
    GSCHUP(&gksl->chup[0], &gksl->chup[1]);
    GSTXP(&gksl->txp);
    GSTXAL(&gksl->txal[0], &gksl->txal[1]);

    GSFAI(&gksl->findex);
    GSFAIS(&gksl->ints);
    GSFASI(&gksl->styli);
    GSFACI(&gksl->facoli);

    for (i = 1; i <= 8; i++)
    {
	tnr = i;
	wn = gksl->window[tnr];
	vp = gksl->viewport[tnr];
	GSWN(&tnr, &wn[0], &wn[1], &wn[2], &wn[3]);
	GSVP(&tnr, &vp[0], &vp[1], &vp[2], &vp[3]);
    }

    GSELNT(&gksl->cntnr);
    GSCLIP(&gksl->clip);

    GSSGT(&gksl->opsg, gksl->mat);

    GSASF(gksl->asf);
}

static
void interp(
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
    char *s;
    gks_state_list *gksl;
    int sp = 0, *len, *f, *ia, *dx, *dy, *dimx, sx = 1, sy = 1, *lc;
    float *r1, *r2;
    char *c;
#ifdef VMS
    struct dsc$descriptor_s text;
#endif
#ifdef cray
    _fcd text;
#endif

#ifdef VMS
    s = chars->dsc$a_pointer;
#else
#ifdef cray
    s = _fcdtocp(chars);
#else
    s = chars;
#endif /* cray */
#endif /* VMS */

    while (s[sp])
    {
	RESOLVE(len, int, sizeof(int));
	RESOLVE(f, int, sizeof(int));

        switch (*f)
	{
	    case  2 :

		RESOLVE(gksl, gks_state_list, sizeof(gks_state_list));
		break;

	    case 12 :		/* polyline */
	    case 13 :		/* polymarker */
	    case 15 :		/* fill area */

		RESOLVE(ia, int, sizeof(int));
		RESOLVE(r1, float, *ia * sizeof(float));
		RESOLVE(r2, float, *ia * sizeof(float));
		break;

	    case 14 :                       /* text */

	        RESOLVE(r1, float, sizeof(float));
	        RESOLVE(r2, float, sizeof(float));
	        RESOLVE(lc, int, sizeof(int));
	        RESOLVE(c, char, 132);
	        break;

	    case 16 :               /* cell array */

	        RESOLVE(r1, float, 2 * sizeof(float));
	        RESOLVE(r2, float, 2 * sizeof(float));
	        RESOLVE(dx, int, sizeof(int));
	        RESOLVE(dy, int, sizeof(int));
	        RESOLVE(dimx, int, sizeof(int));
	        RESOLVE(ia, int, *dimx * *dy * sizeof(int));
	        break;

	    case 19 :		/* set linetype */
	    case 21 :		/* set polyline color index */
	    case 23 :		/* set markertype */
	    case 25 :           /* set polymarker color index */
	    case 30 :           /* set text color index */
	    case 33 :           /* set text path */
	    case 36 :           /* set fillarea interior style */
	    case 37 :           /* set fillarea style index */
	    case 38 :           /* set fillarea color index */
	    case 52 :           /* select normalization transformation */
	    case 53 :           /* set clipping indicator */

		RESOLVE(ia, int, sizeof(int));
	        break;

	    case 27 :		/* set text font and precision */
	    case 34 :		/* set text alignment */

		RESOLVE(ia, int, 2 * sizeof(int));
	        break;

	    case 20 :           /* set linewidth scale factor */
	    case 24 :           /* set marker size scale factor */
	    case 28 :           /* set character expansion factor */
	    case 29 :           /* set character spacing */
	    case 31 :           /* set character height */

		RESOLVE(r1, float, sizeof(float));
	        break;

	    case 32 :           /* set character up vector */

		RESOLVE(r1, float, sizeof(float));
		RESOLVE(r2, float, sizeof(float));
	        break;

	    case 48 :           /* set color representation */

	        RESOLVE(ia, int, sizeof(int));
	        RESOLVE(r1, float, 3 * sizeof(float));
		break;

	    case 49 :           /* set window */
	    case 50 :           /* set viewport */

	        RESOLVE(ia, int, sizeof(int));
	        RESOLVE(r1, float, 2 * sizeof(float));
	        RESOLVE(r2, float, 2 * sizeof(float));
	        break;

	    default:
		gks_fprintf(stderr, "GKS: metafile is corrupted ",
		    "(len=%d, fctid=%d)\n", *len, *f);
		exit(1);
	}

        switch (*f)
	{
	    case  2 : gksinit(gksl); break;

	    case 12 : GPL(ia, r1, r2); break;
	    case 13 : GPM(ia, r1, r2); break;

	    case 14 :
#ifdef VMS
		text.dsc$b_dtype = DSC$K_DTYPE_T;
		text.dsc$b_class = DSC$K_CLASS_S;
		text.dsc$w_length = *lc;
		text.dsc$a_pointer = c;

		GTXS(r1, r2, lc, &text);
#else
#ifdef cray
		text = _cptofcd(c, *lc);

		GTXS(r1, r2, lc, text);
#else
#if defined(_WIN32) && !defined(__GNUC__)
		{
		    unsigned short chars_len = *lc;

		    GTXS(r1, r2, lc, c, chars_len);
		}
#else
		GTXS(r1, r2, lc, c, *lc);
#endif /* _WIN32 */
#endif /* cray */
#endif /* VMS */
	        break;

	    case 15 : GFA(ia, r1, r2); break;
	    case 16 : GCA(&r1[0], &r2[0], &r1[1], &r2[1], dx, dy, &sx, &sy,
			  dimx, dy, ia); break;

	    case 19 : GSLN(ia); break;
	    case 20 : GSLWSC(r1); break;
	    case 21 : GSPLCI(ia); break;
	    case 23 : GSMK(ia); break;
	    case 24 : GSMKSC(r1); break;
	    case 25 : GSPMCI(ia); break;
	    case 27 : GSTXFP(&ia[0], &ia[1]); break;
	    case 28 : GSCHXP(r1); break;
	    case 29 : GSCHSP(r1); break;
	    case 30 : GSTXCI(ia); break;
	    case 31 : GSCHH(r1); break;
	    case 32 : GSCHUP(r1, r2); break;
	    case 33 : GSTXP(ia); break;
	    case 34 : GSTXAL(&ia[0], &ia[1]); break;
	    case 36 : GSFAIS(ia); break;
	    case 37 : GSFASI(ia); break;
	    case 38 : GSFACI(ia); break;

	    case 48 : GSCR(&wkid, ia, &r1[0], &r1[1], &r1[2]); break;

	    case 49 : GSWN(ia, &r1[0], &r1[1], &r2[0], &r2[1]); break;
	    case 50 : GSVP(ia, &r1[0], &r1[1], &r2[0], &r2[1]); break;
	    case 52 : GSELNT(ia); break;
	    case 53 : GSCLIP(ia); break;
	}
    }
}

void STDCALL GKDMFI(
    int *fctid, int *dx, int *dy, int *dimx, int *ia, int *lr1, float *r1,
    int *lr2, float *r2, int *lc, CHARARG(chars), ws_state_list **ptr)
{
    char *s;
    int len;

    p = *ptr;

    switch (*fctid)
    {
      case 2 :		/* open workstation */

	p = (ws_state_list *) malloc(sizeof(ws_state_list));

	p->conid = ia[1];
	p->state = GINACT;

	p->buffer = (char *) readfile(p->conid);
	p->position = 0;

	if (sizeof(char *) > sizeof(int)) {
	    long *la = (long *) ia;
	    *la = (long) p;
	} else
	    *ia = (long) p;

	break;

      case 3 :		/* close workstation */

	if (p->buffer != NULL)
	    free(p->buffer);
	free(p);

	p = NULL;
	break;

      case 102 :	/* get item */

	if (p->buffer != NULL)
	{
	    ia[0] = *(int *)(p->buffer + p->position + sizeof(int));
	    ia[1] = *(int *)(p->buffer + p->position);
	}
	else
	{
	    ia[0] = ia[1] = 0;
	}
	break;

      case 103 :	/* read item */

	if (p->buffer == NULL)
	    break;
#ifdef VMS
	s = chars->dsc$a_pointer;
#else
#ifdef cray
	s = _fcdtocp(chars);
#else
	s = chars;
#endif /* cray */
#endif /* VMS */

	len = *(int *)(p->buffer + p->position);
	if (len < ia[1] * 80)
	{
	    memcpy(s, p->buffer + p->position, len);
#ifdef VMS
	    chars.dsc$a_pointer[len] = '\0';
#endif
	    s[len] = '\0';
	}
	else
	{
	    memset(s, 0, ia[1] * 80);
	    gks_fprintf(stderr, "GKS: item data record is too long\n");
	}
	p->position += len;
	break;

      case 104 :	/* interpret item */

	if (p->buffer != NULL)
	    interp(chars);
	break;
    }
}
