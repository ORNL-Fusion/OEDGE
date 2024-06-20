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
 *	GLIGKS
 *
 * ABSTRACT:
 *
 *	This module contains some definitions for the Graphical
 *	Kernel System (GKS)
 *
 * AUTHOR:
 *
 *	Josef Heinen
 *
 * VERSION:
 *
 *	V4.5
 *
 */


/* Constant definitions */

static int GCONID = 0;	    /* default connection identifier */

/* Workstation types */

static int GWSDEF = 0;	    /* default workstation */
static int GMOUTP = 2;	    /* GKS output metafile */
static int GMINPT = 3;	    /* GKS input metafile */
static int GWSWIS = 5;	    /* Workstation Independent Storage */
static int GCGMO  = 7;	    /* CGM output metafile */
static int GCGMC  = 0x20000;/* CGM character encoding */
static int GCGMB  = 0x30000;/* CGM binary encoding */
static int GCGMT  = 0x40000;/* CGM clear text encoding */
static int GCGMGK = 0x50000;/* CGM GRAFkit compatible encoding */
static int GV330  = 16;     /* DIGITAL VT330 b&w */
static int GV340  = 17;     /* DIGITAL VT340 color */
static int GLN03P = 38;     /* LN03 Plus printer */
static int GVSI   = 42;     /* DIGITAL uVAX I */
static int GHP747 = 51;     /* Hewlett Packard HP7475 */
static int GHP755 = 53;     /* Hewlett Packard HP7550's */
static int GPS    = 61;     /* Postscript */
static int GT4114 = 72;     /* Tektronix 4014 */
static int GT4107 = 82;     /* Tektronix 4107 */
 
/* Aspect source */

static int GBUNDL = 0;	    /* aspect bundled */
static int GINDIV = 1;	    /* aspect individual */

/* Clear control flag */

static int GCONDI = 0;	    /* clear conditionally */
static int GALWAY = 1;	    /* clear always */
 
/* Clipping indicator */

static int GNCLIP = 0;	    /* no clipping */
static int GCLIP  = 1;	    /* clipping */
 
/* Color available */

static int GMONOC = 0;	    /* monochrome display */
static int GCOLOR = 1;	    /* color display */

/* Coordinate switch */

static int GWC    = 0;	    /* world coordinates */
static int GNDC   = 1;	    /* NDC */

/* Deferral mode */

static int GASAP  = 0;	    /* as soon as possible */
static int GBNIG  = 1;	    /* before next global interaction */
static int GBNIL  = 2;	    /* before next interaction */
static int GASTI  = 3;	    /* at some time */

/* Detectability */

static int GUNDET = 0;	    /* undetectable */
static int GDETEC = 1;	    /* detectable */

/* Device coordinate units */

static int GMETRE = 0;	    /* meters */
static int GOTHU  = 1;	    /* other units */
 
/* Display surface empty */

static int GNEMPT = 0;	    /* display surface not empty */
static int GEMPTY = 1;	    /* display surface empty */

/* Dynamic modification */

static int GIRG   = 0;	    /* implicit regeneration */
static int GIMN   = 1;	    /* immediate */

/* Echo switch */

static int GNECHO = 0;	    /* echo disabled */
static int GECHO  = 1;	    /* echo enabled */

/* Fill area interior style */

static int GHOLLO = 0;	    /* interior style hollow */
static int GSOLID = 1;	    /* interior style solid */
static int GPATTR = 2;      /* interior style pattern */
static int GHATCH = 3;      /* interior style hatch */

/* Highlighting */

static int GNORML = 0;	    /* unhighlighted */
static int GHILIT = 1;	    /* highlighted */

/* Input device status */

static int GNONE  = 0;	    /* no input obtained */
static int GOK    = 1;	    /* input obtained */
static int GNCHOI = 2;	    /* no choice input */
static int GNPICK = 2;

/* Input class */

static int GNCLAS = 0;	    /* none */
static int GLOCAT = 1;	    /* locator */
static int GSTROK = 2;	    /* stroke */
static int GVALUA = 3;	    /* valuator */
static int GCHOIC = 4;	    /* choice */
static int GPICK  = 5;	    /* pick */
static int GSTRIN = 6;	    /* string */

/* Implicit regeneration mode */

static int GSUPPD = 0;	    /* implicit regeneration suppressed */
static int GALLOW = 1;	    /* implicit regeneration allowed */

/* Level of GKS */

static int GL0A   = 0;	    /* level 0a */
static int GL0B   = 1;	    /* level 0b */
static int GL0C   = 2;	    /* level 0c */
static int GL1A   = 3;	    /* level 1a */
static int GL1B   = 4;	    /* level 1b */
static int GL1C   = 5;	    /* level 1c */
static int GL2A   = 6;	    /* level 2a */
static int GL2B   = 7;	    /* level 2b */
static int GL2C   = 8;	    /* level 2c */

/* New frame action necessary */

static int GNO    = 0;	    /* new frame action not necessary on update */
static int GYES   = 1;	    /* new frame necessary on update */
 
/* Operating mode */

static int GREQU  = 0;	    /* request mode */
static int GSAMPL = 1;      /* sample mode */
static int GEVENT = 2;	    /* event mode */

/* Operating state value */

static int GGKCL  = 0;	    /* GKS closed */
static int GGKOP  = 1;	    /* GKS open */
static int GWSOP  = 2;	    /* at least one workstation open */
static int GWSAC  = 3;	    /* at least one workstation active */
static int GSGOP  = 4;	    /* at least one segment open */

/* Presence of invalid values */

static int GABSNT = 0;	    /* invalid values absent */
static int GPRSNT = 1;	    /* invalid values present */

/* Regeneration flag */

static int GPOSTP = 0;	    /* regeneration flag suppress */
static int GPERFO = 1;	    /* regeneration flag perform */

/* Relative input priority */

static int GHIGHR = 0;	    /* relative input priority higher */
static int GLOWER = 1;	    /* relative input priority lower */

/* Simultaneous events flag */

static int GNMORE = 0;	    /* nomore */
static int GMORE  = 1;	    /* more */

/* Text alignment horizontal */

static int GAHNOR = 0;	    /* horizontal align normal */
static int GALEFT = 1;	    /* horizontal align left */
static int GACENT = 2;	    /* horizontal align center */
static int GARITE = 3;	    /* horizontal align right */
 
/* Text alignment vertical */

static int GAVNOR = 0;	    /* vertical align normal */
static int GATOP  = 1;	    /* vertical align top */
static int GACAP  = 2;	    /* vertical align cap */
static int GAHALF = 3;	    /* vertical align half */
static int GABASE = 4;	    /* vertical align base */
static int GABOTT = 5;	    /* vertical align bottom */

/* Text path */

static int GRIGHT = 0;	    /* path right */
static int GLEFT  = 1;	    /* path left */
static int GUP    = 2;	    /* path up */
static int GDOWN  = 3;	    /* path down */

/* Text precision */

static int GSTRP  = 0;	    /* text precision string */
static int GCHARP = 1;	    /* text precision character */
static int GSTRKP = 2;	    /* text precision stroke */

/* Type of returned values */

static int GSET   = 0;      /* returned value is set */
static int GREALI = 1;	    /* returned value is realized */

/* Update state */

static int GNPEND = 0;	    /* not pending */
static int GPEND  = 1;	    /* pending */

/* Vector/raster/other type */

static int GVECTR = 0;	    /* vector workstation */
static int GRASTR = 1;	    /* raster workstation */
static int GOTHWK = 2;	    /* other device */

/* Visibility */

static int GINVIS = 0;	    /* invisible */
static int GVISI  = 1;	    /* visible */

/* Workstation category */

static int GOUTPT = 0;	    /* output workstation */
static int GINPUT = 1;	    /* input workstation */
static int GOUTIN = 2;	    /* output/input workstation */
static int GWISS  = 3;	    /* workstation independent segment storage */
static int GMO    = 4;	    /* metafile output */
static int GMI    = 5;	    /* metafile input */

/* Workstation state */

static int GINACT = 0;	    /* work station active */
static int GACTIV = 1;	    /* work station inactive */

/* List of GDP attributes */

static int GPLATT = 0;	    /* GDP polyline bundle */
static int GPMATT = 1;	    /* GDP polymarker bundled */
static int GTXATT = 2;	    /* GDP text bundle */
static int GFAATT = 3;	    /* GDP fill area bundle */

/* Standard linetypes */

static int GLSOLI = 1;	    /* linetype solid */
static int GLDASH = 2;	    /* linetype dashed */
static int GLDOT  = 3;	    /* linetype dotted */
static int GLDASD = 4;	    /* linetype dashed-dotted */
 
/* GKS specific linetypes */

static int GLDS2D = -1;     /* linetype dash-2-dots */
static int GLDS3D = -2;     /* linetype dash-3-dots */
static int GLLGDS = -3;     /* linetype long-dash */
static int GLLSDS = -4;     /* linetype long-short-dash */
static int GLSPDS = -5;     /* linetype spaced-dash */
static int GLSPDT = -6;     /* linetype spaced-dot */
static int GLDBDT = -7;     /* linetype double dots */
static int GLTPDT = -8;     /* linetype triple dots */
 
/* Standard markertypes */

static int GPOINT = 1;	    /* markertype dot */
static int GPLUS  = 2;	    /* markertype plus */
static int GAST   = 3;	    /* markertype asterisk */
static int GOMARK = 4;	    /* markertype circle */
static int GXMARK = 5;	    /* markertype diagonal cross */
 
/* GKS specific markertypes */

static int GMSCIR = -1;     /* markertype solid circle */
static int GMTRU  = -2;     /* markertype hollow up triangle */
static int GMSTRU = -3;     /* markertype solid up triangle */
static int GMTRD  = -4;     /* markertype hollow down triangle */
static int GMSTRD = -5;     /* markertype solid down triangle */
static int GMSQ   = -6;     /* markertype hollow square */
static int GMSSQ  = -7;     /* markertype solid square */
static int GMBT   = -8;     /* markertype hollow bow tie */
static int GMSBT  = -9;     /* markertype solid bow tie */
static int GMHG   = -10;    /* markertype hollow hour glass */
static int GMSHG  = -11;    /* markertype solid hour glass */
static int GMDIA  = -12;    /* markertype hollow diamond */
static int GMSDIA = -13;    /* markertype solid diamond */
 

/* Status code definitions */

#define gks__facility 0x0000085A

#define gks__normal	0x085A8001  /* normal successful completion */
#define gks__error_1	0x085A800A  /* GKS not in proper state. GKS must be in
				       the state GKCL in routine !AS */
#define gks__error_2	0x085A8012  /* GKS not in proper state. GKS must be in
				       the state GKOP in routine !AS */
#define gks__error_3	0x085A801A  /* GKS not in proper state. GKS must be in
				       the state WSAC in routine !AS */
#define gks__error_4	0x085A8022  /* GKS not in proper state. GKS must be in
				       the state SGOP in routine !AS */
#define gks__error_5	0x085A802A  /* GKS not in proper state. GKS must be
				       either in the state WSAC or SGOP in
				       routine !AS */
#define gks__error_6	0x085A8032  /* GKS not in proper state. GKS must be
				       either in the state WSOP or WSAC in
				       routine !AS */
#define gks__error_7	0x085A803A  /* GKS not in proper state. GKS must be in
				       one of the states WSOP,WSAC,SGOP in
				       routine !AS */
#define gks__error_8	0x085A8042  /* GKS not in proper state. GKS must be in
				       one of the states GKOP,WSOP,WSAC,SGOP in
				       routine !AS */
#define gks__error_20	0x085A80A2  /* Specified workstation identifier is
				       invalid in routine !AS */
#define gks__error_21	0x085A80AA  /* Specified connection identifier is
				       invalid in routine !AS */
#define gks__error_22	0x085A80B2  /* Specified workstation type is invalid in
				       routine !AS */
#define gks__error_24	0x085A80C2  /* Specified workstation is open in routine
				       !AS */
#define gks__error_25	0x085A80CA  /* Specified workstation is not open in
				       routine !AS */
#define gks__error_26	0x085A80D2  /* Specified workstation cannot be opened
				       in routine !AS */
#define gks__error_29	0x085A80EA  /* Specified workstation is active in
				       routine !AS */
#define gks__error_30	0x085A80F2  /* Specified workstation is not active in
				       routine !AS */
#define gks__error_50	0x085A8192  /* Transformation number is invalid in
				       routine !AS */
#define gks__error_51	0x085A819A  /* Rectangle definition is invalid in
				       routine !AS */
#define gks__error_52	0x085A81A2  /* Viewport is not within the NDC unit
				       square in routine !AS */
#define gks__error_53	0x085A81AA  /* Workstation window is not within the
				       NDC unit square in routine !AS */
#define gks__error_60	0x085A81E2  /* Polyline index is invalid in routine !AS
				     */
#define gks__error_62	0x085A81F2  /* Linetype is invalid in routine !AS */
#define gks__error_64	0x085A8202  /* Polymarker index is invalid in routine
				       !AS */
#define gks__error_65	0x085A820A  /* Colour index is invalid in routine !AS */
#define gks__error_66	0x085A8212  /* Marker type is invalid in routine !AS */
#define gks__error_68	0x085A8222  /* Text index is invalid in routine !AS */
#define gks__error_70	0x085A8232  /* Text font is invalid in routine !AS */
#define gks__error_72	0x085A8242  /* Character expansion factor is invalid in
				       routine !AS */
#define gks__error_73	0x085A824A  /* Character height is invalid in routine
				       !AS */
#define gks__error_74	0x085A8252  /* Character up vector is invalid in routine
				       !AS */
#define gks__error_75	0x085A825A  /* Fill area index is invalid in routine !AS
				     */
#define gks__error_78	0x085A8272  /* Style index is invalid in routine !AS */
#define gks__error_81	0x085A828A  /* Pattern size value is invalid in routine
				       !AS */
#define gks__error_84	0x085A82A2  /* Dimensions of colour index array are
				       invalid in routine !AS */
#define gks__error_85	0x085A82AA  /* Colour index is invalid in routine !AS */
#define gks__error_88	0x085A82C2  /* Colour is invalid in routine !AS */
#define gks__error_100	0x085A8322  /* Number of points is invalid in routine
				       !AS */
#define gks__error_161	0x085A850A  /* Item length is invalid in routine !AS */
#define gks__debug	0x085A95E3  /* stepped to routine !AS */
#define gks__error_901	0x085A9C2A  /* Invalid connection identifier in routine
				       !AS */
#define gks__error_902	0x085A9C32  /* Invalid workstation type in routine !AS
				     */
#define gks__error_903	0x085A9C3A  /* A workstation of this type is already
				       open in routine !AS */
#define gks__error_904	0x085A9C42  /* No logical unit number available in
				       routine !AS */
#define gks__error_905	0x085A9C4A  /* Open failed in routine !AS */
#define gks__error_906	0x085A9C52  /* Assign failed in routine !AS */
#define gks__error_907	0x085A9C5A  /* Deassign failed in routine !AS */
#define gks__error_908	0x085A9C62  /* VWS not accessible in routine !AS */
#define gks__error_909	0x085A9C6A  /* X display not accessible in routine !AS
				     */

/* Entry point definitions */

#if !defined (cray) && !(defined (_WIN32) && !defined (__GNUC__))
#if defined (VMS) || ((defined (hpux) || defined (aix)) && !defined(NAGware))

#define GOPKS	gopks
#define GCLKS	gclks
#define GOPWK	gopwk
#define GRSGWK  grsgwk
#define GCLWK	gclwk
#define GACWK	gacwk
#define GDAWK	gdawk
#define GCLRWK	gclrwk
#define GUWK	guwk
#define GECLKS	geclks
#define GESC	gesc
#define GMSG	gmsg
#define GMSGS	gmsgs

#define GPL	gpl
#define GPM	gpm
#define GTX	gtx
#define GTXS	gtxs
#define GFA	gfa
#define GCA	gca

#define GSASF	gsasf
#define GSPLI	gspli
#define GSLN	gsln
#define GSLWSC	gslwsc
#define GSPLCI	gsplci
#define GSPMI	gspmi
#define GSMK	gsmk
#define GSMKSC	gsmksc
#define GSPMCI	gspmci
#define GSCR	gscr
#define GSTXI	gstxi
#define GSTXFP	gstxfp
#define GSCHXP	gschxp
#define GSCHH	gschh
#define GSCHUP	gschup
#define GSTXP	gstxp
#define GSCHSP	gschsp
#define GSTXAL	gstxal
#define GSTXCI	gstxci
#define GSFAI	gsfai
#define GSFAIS	gsfais
#define GSFASI	gsfasi
#define GSFACI	gsfaci

#define GSWN	gswn
#define GSVP	gsvp
#define GSELNT	gselnt
#define GSCLIP	gsclip
#define GSWKWN	gswkwn
#define GSWKVP	gswkvp

#define GCRSG   gcrsg
#define GCLSG   gclsg
#define GDSG    gdsg
#define GCSGWK  gcsgwk
#define GASGWK  gasgwk
#define GSSGT   gssgt
#define GSDS    gsds
#define GEVTM   gevtm

#define GINLC	ginlc
#define GRQLC	grqlc
#define GRQSK	grqsk
#define GRQST	grqst
#define GRQCH	grqch

#define GQOPS	gqops
#define GQLVKS	gqlvks
#define GQEWK	gqewk
#define GQMNTN	gqmntn
#define GQOPWK	gqopwk
#define GQACWK	gqacwk
#define GQSGWK  gqsqwk
#define GQOPSG  gqopsg

#define GQASF	gqasf
#define GQPLI	gqpli
#define GQLN	gqln
#define GQLWSC	gqlwsc
#define GQPLCI	gqplci
#define GQPMI	gqpmi
#define GQMK	gqmk
#define GQMKSC	gqmksc
#define GQPMCI	gqpmci
#define GQTXI	gqtxi
#define GQTXFP	gqtxfp
#define GQCHXP	gqchxp
#define GQCHSP	gqchsp
#define GQTXCI	gqtxci
#define GQCHH	gqchh
#define GQCHUP	gqchup
#define GQTXP	gqtxp
#define GQTXAL	gqtxal
#define GQFAI	gqfai
#define GQFAIS	gqfais
#define GQFASI	gqfasi
#define GQFACI	gqfaci

#define GQCNTN	gqcntn
#define GQNT	gqnt
#define GQCLIP	gqclip
#define GQWKC	gqwkc
#define GQWKS	gqwks
#define GQWKCA	gqwkca
#define GQCF	gqcf
#define GQTXX	gqtxx
#define GQTXXS	gqtxxs
#define GQDSP	gqdsp
#define GQMDS	gqmds

#define GKSTST	gkstst

#define GNT	gnt
#define GDNT	gdnt
#define GCNT	gcnt
#define GSDT	gsdt
#define GST	gst
#define GCST	gcst
#define GSCT    gsct
#define GCHH    gchh
#define GPOLIN	gpolin
#define GPOLMK	gpolmk
#define GFILLA	gfilla
#define GKSFA	gksfa
#define GSIMTX	gsimtx
#define GMOVE	gmove
#define GDASH	gdash
#define GSIMPM	gsimpm
#define GTEXTS	gtexts
#define GQRGB	gqrgb
#define GSRGB	gsrgb
#define GQPIX	gqpix
#define GSPIX	gspix
#define GKAFM	gkafm

#define GKDXW	gkdxw
#define GKDGA	gkdga
#define GKDCGM	gkdcgm
#define GKDBM	gkdbm
#define GKDMFO	gkdmfo
#define GKDMFI	gkdmfi
#define GKDPDF	gkdpdf
#define GKDCSG	gkdcsg
#define GKCSG	gkcsg
#define GKFVP	gkfvp
#define GKQPA	gkqpa
#define GKSPA	gkspa

#define OPFONT	opfont
#define LOOKUP	lookup
#define CLFONT	clfont
#define BUFIN	bufin
#define BUFOUT	bufout
#define BINOUT	binout
#define DPSOP	dpsop
#define DPSWR	dpswr
#define DPSFL	dpsfl
#define DPSCL	dpscl
#define GKINFO  gkinfo
#define GKTMP   gktmp
#define GKMAGS  gkmags
#define LIBSIG	libsig
#define GKFD    gkfd
#define GTWSTY	gtwsty
#define GKLTOI	gkltoi
#define GERSET  gerset
#define FORTWR	fortwr
#define GKOPEN  gkopen
#define GKCLOS  gkclos

#else

#define GOPKS	gopks_
#define GCLKS	gclks_
#define GOPWK	gopwk_
#define GRSGWK  grsgwk_
#define GCLWK	gclwk_
#define GACWK	gacwk_
#define GDAWK	gdawk_
#define GCLRWK	gclrwk_
#define GUWK	guwk_
#define GECLKS	geclks_
#define GESC	gesc_
#define GMSG	gmsg_
#define GMSGS	gmsgs_

#define GPL	gpl_
#define GPM	gpm_
#define GTX	gtx_
#define GTXS	gtxs_
#define GFA	gfa_
#define GCA	gca_

#define GSASF	gsasf_
#define GSPLI	gspli_
#define GSLN	gsln_
#define GSLWSC	gslwsc_
#define GSPLCI	gsplci_
#define GSPMI	gspmi_
#define GSMK	gsmk_
#define GSMKSC	gsmksc_
#define GSPMCI	gspmci_
#define GSCR	gscr_
#define GSTXI	gstxi_
#define GSTXFP	gstxfp_
#define GSCHXP	gschxp_
#define GSCHH	gschh_
#define GSCHUP	gschup_
#define GSTXP	gstxp_
#define GSCHSP	gschsp_
#define GSTXAL	gstxal_
#define GSTXCI	gstxci_
#define GSFAI	gsfai_
#define GSFAIS	gsfais_
#define GSFASI	gsfasi_
#define GSFACI	gsfaci_

#define GSWN	gswn_
#define GSVP	gsvp_
#define GSELNT	gselnt_
#define GSCLIP	gsclip_
#define GSWKWN	gswkwn_
#define GSWKVP	gswkvp_

#define GCRSG   gcrsg_
#define GCLSG   gclsg_
#define GDSG    gdsg_
#define GCSGWK  gcsgwk_
#define GASGWK  gasgwk_
#define GSSGT   gssgt_
#define GSDS    gsds_
#define GEVTM   gevtm_

#define GINLC	ginlc_
#define GRQLC	grqlc_
#define GRQSK	grqsk_
#define GRQST	grqst_
#define GRQCH	grqch_

#define GQOPS	gqops_
#define GQLVKS	gqlvks_
#define GQEWK	gqewk_
#define GQMNTN	gqmntn_
#define GQOPWK	gqopwk_
#define GQACWK	gqacwk_
#define GQSGWK  gqsqwk_
#define GQOPSG  gqopsg_

#define GQASF	gqasf_
#define GQPLI	gqpli_
#define GQLN	gqln_
#define GQLWSC	gqlwsc_
#define GQPLCI	gqplci_
#define GQPMI	gqpmi_
#define GQMK	gqmk_
#define GQMKSC	gqmksc_
#define GQPMCI	gqpmci_
#define GQTXFP	gqtxfp_
#define GQCHXP	gqchxp_
#define GQCHSP	gqchsp_
#define GQTXI	gqtxi_
#define GQTXCI	gqtxci_
#define GQCHH	gqchh_
#define GQCHUP	gqchup_
#define GQTXP	gqtxp_
#define GQTXAL	gqtxal_
#define GQFAI	gqfai_
#define GQFAIS	gqfais_
#define GQFASI	gqfasi_
#define GQFACI	gqfaci_

#define GQCNTN	gqcntn_
#define GQNT	gqnt_
#define GQCLIP	gqclip_
#define GQWKC	gqwkc_
#define GQWKS	gqwks_
#define GQWKCA	gqwkca_
#define GQCF	gqcf_
#define GQTXX	gqtxx_
#define GQTXXS	gqtxxs_
#define GQDSP	gqdsp_
#define GQMDS	gqmds_

#define GKSTST	gkstst_

#define GNT	gnt_
#define GDNT	gdnt_
#define GCNT	gcnt_
#define GSDT	gsdt_
#define GST	gst_
#define GCST	gcst_
#define GSCT    gsct_
#define GCHH    gchh_
#define GPOLIN	gpolin_
#define GPOLMK	gpolmk_
#define GFILLA	gfilla_
#define GKSFA	gksfa_
#define GSIMTX	gsimtx_
#define GMOVE	gmove_
#define GDASH	gdash_
#define GSIMPM	gsimpm_
#define GTEXTS	gtexts_
#define GQRGB	gqrgb_
#define GSRGB	gsrgb_
#define GQPIX	gqpix_
#define GSPIX	gspix_
#define GKAFM	gkafm_

#define GKDXW	gkdxw_
#define GKDGA	gkdga_
#define GKDCGM	gkdcgm_
#define GKDBM	gkdbm_
#define GKDMFO	gkdmfo_
#define GKDMFI	gkdmfi_
#define GKDPDF	gkdpdf_
#define GKDCSG	gkdcsg_
#define GKCSG	gkcsg_
#define GKFVP	gkfvp_
#define GKQPA	gkqpa_
#define GKSPA	gkspa_

#define OPFONT	opfont_
#define LOOKUP	lookup_
#define CLFONT	clfont_
#define BUFIN	bufin_
#define BUFOUT	bufout_
#define BINOUT	binout_
#define DPSOP	dpsop_
#define DPSWR	dpswr_
#define DPSFL	dpsfl_
#define DPSCL	dpscl_
#define GKINFO  gkinfo_
#define GKTMP   gktmp_
#define GKMAGS  gkmags_
#define LIBSIG	libsig_
#define GKFD    gkfd_
#define GTWSTY	gtwsty_
#define GKLTOI	gkltoi_
#define GERSET  gerset_
#define FORTWR	fortwr_
#define GKOPEN  gkopen_
#define GKCLOS  gkclos_

#endif
#endif /* cray, _WIN32 */

#define FPRINTF (*gks_a_fprintf_routine)

#if defined (_WIN32) && !defined (__GNUC__)

extern void __stdcall GOPKS(int *, int *);
extern void __stdcall GCLKS();
extern void __stdcall GOPWK(int *, int *, int *);
extern void __stdcall GRSGWK(int *);
extern void __stdcall GCLWK(int *);
extern void __stdcall GACWK(int *);
extern void __stdcall GDAWK(int *);
extern void __stdcall GCLRWK(int *, int *);
extern void __stdcall GUWK(int *, int *);
extern void __stdcall GECLKS();
extern void __stdcall GESC(int *, int *, char *, unsigned short, int *, int *,
    char *, int);
extern void __stdcall GMSG(int *, char *, unsigned short);
extern void __stdcall GMSGS(int *, int *, char *, unsigned short);

extern void __stdcall GPL(int *, float *, float *);
extern void __stdcall GPM(int *, float *, float *);
extern void __stdcall GTX(float *, float *, char *, unsigned short);
extern void __stdcall GTXS(float *, float *, int *, char *, unsigned short);
extern void __stdcall GFA(int *, float *, float *);
extern void __stdcall GCA(float *, float *, float *, float *, int *, int *,
    int *, int *, int *, int *, int *);

extern void __stdcall GSASF(int *);
extern void __stdcall GSPLI(int *);
extern void __stdcall GSLN(int *);
extern void __stdcall GSLWSC(float *);
extern void __stdcall GSPLCI(int *);
extern void __stdcall GSPMI(int *);
extern void __stdcall GSMK(int *);
extern void __stdcall GSMKSC(float *);
extern void __stdcall GSPMCI(int *);
extern void __stdcall GSCR(int *, int *, float *, float *, float *);
extern void __stdcall GSTXI(int *);
extern void __stdcall GSTXFP(int *, int *);
extern void __stdcall GSCHXP(float *);
extern void __stdcall GSCHH(float *);
extern void __stdcall GSCHUP(float *, float *);
extern void __stdcall GSTXP(int *);
extern void __stdcall GSCHSP(float *);
extern void __stdcall GSTXAL(int *, int *);
extern void __stdcall GSTXCI(int *);
extern void __stdcall GSFAI(int *);
extern void __stdcall GSFAIS(int *);
extern void __stdcall GSFASI(int *);
extern void __stdcall GSFACI(int *);

extern void __stdcall GSWN(int *, float *, float *, float *, float *);
extern void __stdcall GSVP(int *, float *, float *, float *, float *);
extern void __stdcall GSELNT(int *);
extern void __stdcall GSCLIP(int *);
extern void __stdcall GSWKWN(int *, float *, float *, float *, float *);
extern void __stdcall GSWKVP(int *, float *, float *, float *, float *);

extern void __stdcall GCRSG(int *);
extern void __stdcall GCLSG();
extern void __stdcall GDSG(int *);
extern void __stdcall GCSGWK(int *, int *);
extern void __stdcall GASGWK(int *, int *);
extern void __stdcall GSSGT(int *, float [3][2]);
extern void __stdcall GSDS(int *, int *, int *);
extern void __stdcall GEVTM(float *, float *, float *, float *,
    float *, float *, float *, int *, float [3][2]);

extern void __stdcall GINLC(int *, int *, int *, float *, float *, int *,
    float *, float *, float *, float *, int *, char *, unsigned short);
extern void __stdcall GRQLC(int *, int *, int *, int *, float *, float *);
extern void __stdcall GRQSK(int *, int *, int *, int *, int *, int *,
    float *, float *);
extern void __stdcall GRQST(int *, int *, int *, int *, char *, unsigned short);
extern void __stdcall GRQCH(int *, int *, int *, int *);

extern void __stdcall GQOPS(int *);
extern void __stdcall GQLVKS(int *, int *);
extern void __stdcall GQEWK(int *, int *, int *, int *);
extern void __stdcall GQMNTN(int *, int *);
extern void __stdcall GQOPWK(int *, int *, int *, int *);
extern void __stdcall GQACWK(int *, int *, int *, int *);
extern void __stdcall GQSGWK(int *, int *, int *, int *, int *);
extern void __stdcall GQOPSG(int *, int *);

extern void __stdcall GQASF(int *, int *);
extern void __stdcall GQPLI(int *, int *);
extern void __stdcall GQLN(int *, int *);
extern void __stdcall GQLWSC(int *, float *);
extern void __stdcall GQPLCI(int *, int *);
extern void __stdcall GQPMI(int *, int *);
extern void __stdcall GQMK(int *, int *);
extern void __stdcall GQMKSC(int *, float *);
extern void __stdcall GQPMCI(int *, int *);
extern void __stdcall GQTXI(int *, int *);
extern void __stdcall GQTXFP(int *, int *, int *);
extern void __stdcall GQCHXP(int *, float *);
extern void __stdcall GQCHSP(int *, float *);
extern void __stdcall GQTXCI(int *, int *);
extern void __stdcall GQCHH(int *, float *);
extern void __stdcall GQCHUP(int *, float *, float *);
extern void __stdcall GQTXP(int *, int *);
extern void __stdcall GQTXAL(int *, int *, int *);
extern void __stdcall GQFAI(int *, int *);
extern void __stdcall GQFAIS(int *, int *);
extern void __stdcall GQFASI(int *, int *);
extern void __stdcall GQFACI(int *, int *);

extern void __stdcall GQCNTN(int *, int *);
extern void __stdcall GQNT(int *, int *, float *, float *);
extern void __stdcall GQCLIP(int *, int *, float *);
extern void __stdcall GQWKC(int *, int *, int *, int *);
extern void __stdcall GQWKS(int *, int *, int *);
extern void __stdcall GQWKCA(int *, int *, int *);
extern void __stdcall GQCF(int *, int *, int *, int *, int *);
extern void __stdcall GQTXX(
    int *, float *, float *, char *, unsigned short,
    int *, float *, float *, float *, float *);
extern void __stdcall GQTXXS(
    int *, float *, float *, int *, char *, unsigned short,
    int *, float *, float *, float *, float *);
extern void __stdcall GQDSP(int *, int *, int *, float *, float *,
    int *, int *);
extern void __stdcall GQMDS(int *, int *, int *, float *, float *,
    int *, int *);

extern void __stdcall GKSTST(void);

extern void __stdcall GCST(float *, float *);
extern void __stdcall GCHH(float *);
extern void __stdcall GKDCSG(int *, int *, int *, float *, float *,
    char *, unsigned short);
extern void __stdcall GKFVP(float *, float *, float *);
extern void __stdcall GKQPA(int *, int [33]);
extern void __stdcall GSDT(float [4], float [4]);
extern void __stdcall GPOLIN(int *, float *, float *, int *, int *,
    void (*)(), void (*)());
extern void __stdcall GMOVE(float *, float *, void (*)());
extern void __stdcall GDASH(float *, float *, void (*)(), void (*)());
extern void __stdcall GSIMPM(int *, float *, float *, void (*)());
extern void __stdcall GSIMTX(float *, float *, int *, char *, unsigned short,
    void (*)(), void (*)());
extern void __stdcall GQRGB(int *, float *, float *, float *);
extern void __stdcall GSCT();
extern void __stdcall GSRGB(int *, float *, float *, float *);
extern void __stdcall GST(float *, float *);
extern void __stdcall GKAFM(int *, int *, int *);

#ifndef __gks
extern void __stdcall BUFOUT(int *, int *, char *, unsigned short);
#endif

#else

#if !defined(VMS) && !defined(cray)

void GOPKS(int *, int *);
void GCLKS();
void GOPWK(int *, int *, int *);
void GRSGWK(int *);
void GCLWK(int *);
void GACWK(int *);
void GDAWK(int *);
void GCLRWK(int *, int *);
void GUWK(int *, int *);
void GECLKS();
void GESC(int *, int *, char *, int *, int *, char *, unsigned short,
    unsigned short);
void GMSG(int *, char *, unsigned short);
void GMSGS(int *, int *, char *, unsigned short);

void GPL(int *, float *, float *);
void GPM(int *, float *, float *);
void GTX(float *, float *, char *, unsigned short);
void GTXS(float *, float *, int *, char *, unsigned short);
void GFA(int *, float *, float *);
void GCA(float *, float *, float *, float *, int *, int *, int *, int *, int *,
    int *, int *);

void GSASF(int *);
void GSPLI(int *);
void GSLN(int *);
void GSLWSC(float *);
void GSPLCI(int *);
void GSPMI(int *);
void GSMK(int *);
void GSMKSC(float *);
void GSPMCI(int *);
void GSCR(int *, int *, float *, float *, float *);
void GSTXI(int *);
void GSTXFP(int *, int *);
void GSCHXP(float *);
void GSCHH(float *);
void GSCHUP(float *, float *);
void GSTXP(int *);
void GSCHSP(float *);
void GSTXAL(int *, int *);
void GSTXCI(int *);
void GSFAI(int *);
void GSFAIS(int *);
void GSFASI(int *);
void GSFACI(int *);

void GSWN(int *, float *, float *, float *, float *);
void GSVP(int *, float *, float *, float *, float *);
void GSELNT(int *);
void GSCLIP(int *);
void GSWKWN(int *, float *, float *, float *, float *);
void GSWKVP(int *, float *, float *, float *, float *);

void GCRSG(int *);
void GCLSG();
void GDSG(int *);
void GCSGWK(int *, int *);
void GASGWK(int *, int *);
void GSSGT(int *, float [3][2]);
void GSDS(int *, int *, int *);
void GEVTM(float *, float *, float *, float *, float *, float *, float *,
    int *, float [3][2]);

void GINLC(int *, int *, int *, float *, float *, int *, float *, float *,
    float *, float *, int *, char *, unsigned short);
void GRQLC(int *, int *, int *, int *, float *, float *);
void GRQSK(int *, int *, int *, int *, int *, int *, float *, float *);
void GRQST(int *, int *, int *, int *, char *, unsigned short);
void GRQCH(int *, int *, int *, int *);

void GQOPS(int *);
void GQLVKS(int *, int *);
void GQEWK(int *, int *, int *, int *);
void GQMNTN(int *, int *);
void GQOPWK(int *, int *, int *, int *);
void GQACWK(int *, int *, int *, int *);
void GQSGWK(int *, int *, int *, int *, int *);
void GQOPSG(int *, int *);

void GQASF(int *, int *);
void GQPLI(int *, int *);
void GQLN(int *, int *);
void GQLWSC(int *, float *);
void GQPLCI(int *, int *);
void GQPMI(int *, int *);
void GQMK(int *, int *);
void GQMKSC(int *, float *);
void GQPMCI(int *, int *);
void GQTXI(int *, int *);
void GQTXFP(int *, int *, int *);
void GQCHXP(int *, float *);
void GQCHSP(int *, float *);
void GQTXCI(int *, int *);
void GQCHH(int *, float *);
void GQCHUP(int *, float *, float *);
void GQTXP(int *, int *);
void GQTXAL(int *, int *, int *);
void GQFAI(int *, int *);
void GQFAIS(int *, int *);
void GQFASI(int *, int *);
void GQFACI(int *, int *);

void GQCNTN(int *, int *);
void GQNT(int *, int *, float *, float *);
void GQCLIP(int *, int *, float *);
void GQWKC(int *, int *, int *, int *);
void GQWKS(int *, int *, int *);
void GQWKCA(int *, int *, int *);
void GQCF(int *, int *, int *, int *, int *);
void GQTXX(int *, float *, float *, char *, int *, float *, float *, float *,
    float *, unsigned short);
void GQTXXS(int *, float *, float *, int *, char *, int *, float *, float *,
    float *, float *, unsigned short);
void GQDSP(int *, int *, int *, float *, float *, int *, int *);
void GQMDS(int *, int *, int *, float *, float *, int *, int *);

void GKSTST(void);

void GCST(float *, float *);
void GCHH(float *);
void GKDCSG(int *, int *, int *, float *, float *, char *, unsigned short);
void GKFVP(float *, float *, float *);
void GKQPA(int *, int [33]);
void GSDT(float [4], float [4]);
void GPOLIN(int *, float *, float *, int *, int *, void (*)(float *, float *),
    void (*)(float *, float *));
void GMOVE(float *, float *, void (*)(float *, float *));
void GDASH(float *, float *, void (*)(float *, float *),
    void (*)(float *, float *));
void GSIMPM(int *, float *, float *, void (*)(float *, float *));
void GSIMTX(float *, float *, int *, char *,
    void (*)(int *, float *, float *, int *, int *),
    void (*)(int *, float *, float *, int *), unsigned short);
void GQRGB(int *, float *, float *, float *);
void GSCT(void);
void GSRGB(int *, float *, float *, float *);
void GST(float *, float *);
void GSPIX(int *, int *);
void GKAFM(int *, int *, int *);
void GTEXTS(float *, float *, int *, char *, int *,
    void (*)(int *, float *, float *, int *, int *),
    void (*)(int *, float *, float *, int *),
    void (*)(float *, float *, int *, char *));
void GKSPA(int *, int *);
void FORTWR(int *, int *, char *, int *, unsigned short);

#ifndef __gks
void BUFOUT(int *, int *, char *, unsigned short);
#endif

void gks_fprintf (FILE *file, char *format, ...);
void gkscompress(int, unsigned char *, int, unsigned char  *, int *);

#endif /* VMS */

#endif

typedef struct gks_state_list_struct {
    int lindex;
    int ltype;
    float lwidth;
    int plcoli;
    int mindex;
    int mtype;
    float mszsc;
    int pmcoli;
    int tindex;
    int txfont, txprec;
    float chxp;
    float chsp;
    int txcoli;
    float chh;
    float chup[2];
    int txp;
    int txal[2];
    int findex;
    int ints;
    int styli;
    int facoli;
    float window[9][4], viewport[9][4];
    int cntnr, clip, opsg;
    float mat[3][2];
    int asf[13];
    } gks_state_list;


/* Global variables */

#ifndef __gks

extern char *gks_a_error_info;
extern void (*gks_a_fprintf_routine)();

#endif /* __gks */
