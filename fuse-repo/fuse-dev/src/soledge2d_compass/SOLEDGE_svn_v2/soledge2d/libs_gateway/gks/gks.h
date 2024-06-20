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
 */

#include <stdio.h>

/* Simple types */
 
typedef FILE            Gfile;
typedef char            Gchar;
typedef char            Gconn;
typedef float           Gfloat ;
typedef int             Gwstype;
typedef int             Gint;
typedef unsigned int    Guint;
typedef long            Glong;

#define GWC_DEF         NULL
#define GCONID_DEFAULT  NULL

#define GWS_DEFAULT     NULL		/* Default Workstation type */
#define GWS_DEF         NULL

/* Enumerations */

typedef enum {          /* aspect control flag */
        GBUNDLED,
        GINDIVIDUAL
} Gasf;

typedef enum {          /* clipping indicator */
        GCLIP,
        GNOCLIP
} Gclip;

typedef enum {          /* clear control flag */
        GCONDITIONALLY,
        GALWAYS
} Gclrflag;

typedef enum {          /* co-ord switch */
        GWC,
        GNDC
} Gcsw;

typedef enum {          /* device co-ord units */
        GDC_METRES,
        GDC_OTHER
} Gdevunits;

typedef enum {          /* fill area interior style */
        GHOLLOW,
        GSOLID,
        GPATTERN,
        GHATCH
} Gflinter;

typedef enum {          /* request status */
        GOK,
        GNONE
} Gistat;

typedef enum {          /* line type - not standard */
        GLN_SOLID = 1,
        GLN_DASHED,
        GLN_DOTTED,
        GLN_DASHDOT,
        GLN_TRIPLE_DOT = -8,
        GLN_DOUBLE_DOT,
        GLN_SPACED_DOT,
        GLN_SPACED_DASH,
        GLN_LONG_SHORT_DASH,
        GLN_LONG_DASH,
        GLN_DASH_3_DOT,
        GLN_DASH_2_DOT
} Glntype;

typedef enum {          /* marker type - not standard */
        GMK_POINT = 1,
        GMK_PLUS,
        GMK_STAR,
        GMK_O,
        GMK_X,
        GMK_SOLID_DIAMOND = -13,
        GMK_DIAMOND,
        GMK_SOLID_HGLASS,
        GMK_HOURGLASS,
        GMK_SOLID_BOWTIE,
        GMK_BOWTIE,
        GMK_SOLID_SQUARE,
        GMK_SQUARE,
        GMK_SOLID_TRI_DOWN,
        GMK_TRIANGLE_DOWN,
        GMK_SOLID_TRI_UP,
        GMK_TRIANGLE_UP,
        GMK_SOLID_CIRCLE
} Gmktype;

typedef enum {          /* GKS operating state */
        GGKCL,          /* closed */
        GGKOP,          /* open */
        GWSOP,          /* workstation open */
        GWSAC,          /* workstation active */
        GSGOP           /* segment open */
} Gopst;

typedef enum {          /* regeneration flag */
        GPERFORM,
        GPOSTPONE
} Gregen;

typedef enum {          /* horiz text alignment component */
        GAH_NORMAL,
        GAH_LEFT,
        GAH_CENTRE,
        GAH_RIGHT
} Gtxhor;

typedef enum {          /* text path */
        GTP_RIGHT,
        GTP_LEFT,
        GTP_UP,
        GTP_DOWN
} Gtxpath;

typedef enum {          /* text precision */
        GP_STRING,
        GP_CHAR,
        GP_STROKE
} Gtxprec;

typedef enum {          /* vert text alignment component */
        GAV_NORMAL,
        GAV_TOP,
        GAV_CAP,
        GAV_HALF,
        GAV_BASE,
        GAV_BOTTOM
} Gtxver;

typedef enum {          /* WS category */
        GOUTPUT,
        GINPUT,
        GOUTIN,
        GWISS,
        GMO,
        GMI
} Gwscat;
 
/* Forward type definitions */

typedef struct {        /* integer point */
        Gint    x;              /* x coordinate */
        Gint    y;              /* y coordinate */
} Gipoint;

typedef struct {        /* coordinate point */
        Gfloat  x;              /* X coordinate */
        Gfloat  y;              /* Y coordinate */
} Gpoint;

typedef struct {        /* coordinate limits */
        Gfloat  xmin;           /* x minimum limit */
        Gfloat  xmax;           /* x maximum limit */
        Gfloat  ymin;           /* y minimum limit */
        Gfloat  ymax;           /* y maximum limit */
} Glimit;

typedef struct {        /* text facilities */
        Gint    font;           /* text font */
        Gtxprec prec;           /* text precision */
} Gtxfp;

typedef struct {        /* text alignment */
        Gtxhor  hor;            /* horizontal component */
        Gtxver  ver;            /* vertical component */
} Gtxalign;

typedef struct {        /* coordinate rectangle pointer */
        Gpoint  ul;             /* upper left point */
        Gpoint  lr;             /* lower right point */
} Grect;

typedef struct {        /* dimension in integer values */
        Guint   x_dim;          /* X dimension */
        Guint   y_dim;          /* Y dimension */
} Gidim;

/* Structs */
 
typedef struct {        /* aspect source flags */
        Gasf    ln_type;        /* line type */
        Gasf    ln_width;       /* line width */
        Gasf    ln_colour;      /* line colour */
        Gasf    mk_type;        /* marker type */
        Gasf    mk_size;        /* marker size */
        Gasf    mk_colour;      /* marker colour */
        Gasf    tx_fp;          /* text font and precision */
        Gasf    tx_exp;         /* text expansion */
        Gasf    tx_space;       /* text character spacing */
        Gasf    tx_colour;      /* text colour */
        Gasf    fl_inter;       /* fill area interior style */
        Gasf    fl_style;       /* fill area style index */
        Gasf    fl_colour;      /* fill area colour */
} Gasfs;

typedef struct {        /* colour bundle */
        Gfloat  red;            /* red intensity */
        Gfloat  green;          /* green intensity */
        Gfloat  blue;           /* blue intensity */
} Gcobundl;

typedef struct {        /* clipping rectangle */
        Gclip   ind;    /* clipping indicator */
        Glimit  rec;    /* clipping rectangle */
} Gcliprect;

typedef struct {        /* display size */
        Gdevunits units;        /* device coordinate units */
        Gpoint  device;         /* device coordinate unit size */
        Gipoint raster;         /* raster unit size */
} Gdspsize;

typedef struct {        /* text extent */
        Gpoint  concat;         /* concatenation point */
        Gpoint  corner_1;       /* corner 1 */
        Gpoint  corner_2;       /* corner 2 */
        Gpoint  corner_3;       /* corner 3 */
        Gpoint  corner_4;       /* corner 4 */
} Gextent;

typedef struct {        /* locator data */
        Gint    transform;      /* normalization transformation number */
        Gpoint  position;       /* locator position */
} Gloc;

typedef struct {        /* request locator */
        Gistat  status;         /* request status */
        Gloc    loc;            /* locator data */
} Gqloc;

typedef struct {        /* request string */
        Gistat  status;         /* request status */
        Gchar   *string;        /* string data */
} Gqstring;
 
typedef struct {        /* scale vector */
        Gfloat x_scale;
        Gfloat y_scale;
} Gscale;

typedef struct {        /* transformation */
        Glimit  w;              /* window */
        Glimit  v;              /* viewport */
} Gtran;

typedef struct {        /* metafile item information */
        Gint    type;           /* item type */
        Gint    length;         /* item length */
} Ggksmit;

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

#endif
#endif /* cray, _WIN32 */

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

#endif

#define gsetlinecolorind gsetlinecolourind
#define gsetmarkercolorind gsetmarkercolourind
#define gsettextcolorind gsettextcolourind
#define gsetfillcolorind gsetfillcolourind
#define gsetcolorrep gsetcolourrep
#define ginqlinecolorind ginqlinecolourind
#define ginqmarkercolorind ginqmarkercolourind
#define ginqtextcolorind ginqtextcolourind
#define ginqfillcolorind ginqfillcolourind

int gopengks(Gfile *, Glong);
int gclosegks(void);
int gopenws(Gint, Gconn *, Gwstype *);
int gclosews(Gint);
int gactivatews(Gint);
int gdeactivatews(Gint);
int gclearws(Gint, Gclrflag);
int gupdatews(Gint, Gregen);
int gmessage(Gint, Gchar *);
int gpolyline(Gint, Gpoint *);
int gpolymarker(Gint, Gpoint *);
int gtext(Gpoint *, Gchar *);
int gfillarea(Gint, Gpoint *);
int gcellarray(Grect *, Gidim *, Gint *);
int gsetasf(Gasfs *);
int gsetlineind(Gint);
int gsetlinetype(Gint);
int gsetlinewidth(Gfloat);
int gsetlinecolourind(Gint);
int gsetmarkerind(Gint);
int gsetmarkertype(Gint);
int gsetmarkersize(Gfloat);
int gsetmarkercolourind(Gint);
int gsettextind(Gint);
int gsettextfontprec(Gtxfp *);
int gsetcharexpan(Gfloat);
int gsetcharspace(Gfloat);
int gsettextcolourind(Gint);
int gsetcharheight(Gfloat);
int gsetcharup(Gpoint *);
int gsettextpath(Gtxpath);
int gsettextalign(Gtxalign *);
int gsetfillind(Gint);
int gsetfillintstyle(Gflinter);
int gsetfillstyle(Gint);
int gsetfillcolourind(Gint);
int gsetcolourrep(Gint, Gint, Gcobundl *);
int gsetwindow(Gint, Glimit *);
int gsetviewport(Gint, Glimit *);
int gselntran(Gint);
int gsetclip(Gclip);
int gsetwswindow(Gint, Glimit *);
int gsetwsviewport(Gint, Glimit *);
int greqloc(Gint, Gint, Gqloc *);
int greqstring(Gint, Gint, Gqstring *);
int gcreateseg(Gint);
int gcopysegws(Gint, Gint);
int gredrawsegws(Gint);
int gcloseseg(void);
int gevaltran(Gpoint *, Gpoint *, Gfloat, Gscale *, Gcsw, Gfloat [3][2]);
int gsetsegtran(Gint, Gfloat [3][2]);
int ginqopst(Gint *);
int ginqlevelgks(Gint *, Gint *);
int ginqmaxntrannum( Gint *, Gint *);
int ginqcharheight(Gfloat *, Gint *);
int ginqcharup(Gpoint *, Gint *);
int ginqtextpath(Gtxpath *, Gint *);
int ginqtextalign(Gtxalign *, Gint *);
int ginqasf(Gasfs *, Gint *);
int ginqlineind(Gint *, Gint *);
int ginqlinetype(Gint *, Gint *);
int ginqlinewidth(Gfloat *, Gint *);
int ginqlinecolourind(Gint *, Gint *);
int ginqmarkerind(Gint *, Gint *);
int ginqmarkertype(Gint *, Gint *);
int ginqmarkersize(Gfloat *, Gint *);
int ginqmarkercolourind(Gint *, Gint *);
int ginqtextind(Gint *, Gint *);
int ginqtextfontprec(Gtxfp *, Gint *);
int ginqcharexpan(Gfloat *, Gint *);
int ginqcharspace(Gfloat *, Gint *);
int ginqtextcolourind(Gint *, Gint *);
int ginqfillind(Gint *, Gint *);
int ginqfillintstyle(Gint *, Gint *);
int ginqfillstyle(Gint *, Gint *);
int ginqfillcolourind(Gint *, Gint *);
int ginqcurntrannum(Gint *, Gint *);
int ginqntran(Gint, Gtran *, Gint *);
int ginqclip(Gcliprect *, Gint *);
int ginqwsst(Gint, Gint *, Gint *);
int ginqwscategory(Gwstype *, Gint *, Gint *);
int ginqwscf(Gwstype *, Gint *, Gint *, Gint *, Gint *);
int ginqdisplaysize(Gwstype *, Gdspsize *, Gint *);
int ginqtextextent(Gint, Gpoint *, Gchar *, Gextent *, Gint *);
int ginqnameopenseg(Gint *, Gint *);
int gemergencyclosegks(void);
