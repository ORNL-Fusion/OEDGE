C*
C* Copyright @ 1984 - 1993   Josef Heinen
C*
C* Permission to use, copy, and distribute this software and its
C* documentation for any purpose with or without fee is hereby granted,
C* provided that the above copyright notice appear in all copies and
C* that both that copyright notice and this permission notice appear
C* in supporting documentation.
C*
C* Permission to modify the software is granted, but not the right to
C* distribute the modified code.  Modifications are to be distributed
C* as patches to released version.
C*
C* This software is provided "as is" without express or implied warranty.
C*
C* Send your comments or suggestions to
C*  J.Heinen@kfa-juelich.de.
C*
C*
C*      GKS Symbol Table

C  Mnemonic FORTRAN names and their values for GKS ENUMERATION type values

C    aspect source
        INTEGER GBUNDL,GINDIV
        PARAMETER (GBUNDL=0,GINDIV=1)

C    clear control flag
        INTEGER GCONDI,GALWAY
        PARAMETER (GCONDI=0,GALWAY=1)

C    clipping indicator
        INTEGER GNCLIP,GCLIP
        PARAMETER (GNCLIP=0,GCLIP=1)

C    colour available
        INTEGER GMONOC,GCOLOR
        PARAMETER (GMONOC=0,GCOLOR=1)

C    coordinate switch
        INTEGER GWC,GNDC
        PARAMETER (GWC=0,GNDC=1)

C    deferral mode
        INTEGER GASAP,GBNIL,GBNIG,GASTI
        PARAMETER (GASAP=0,GBNIL=1,GBNIG=2,GASTI=3)

C    detectability
        INTEGER GUNDET,GDETEC
        PARAMETER (GUNDET=0,GDETEC=1)

C    device coordinate units
        INTEGER GMETRE,GOTHER
        PARAMETER (GMETRE=0,GOTHER=1)

C    display surface empty
        INTEGER GNEMPT,GEMPTY
        PARAMETER (GNEMPT=0,GEMPTY=1)

C    dynamic modification
        INTEGER GIRG,GIMM
        PARAMETER (GIRG=0,GIMM=1)

C    echo switch
        INTEGER GNECHO,GECHO
        PARAMETER (GNECHO=0,GECHO=1)

C    fill area interior style
        INTEGER GHOLLO,GSOLID,GPATTR,GHATCH
        PARAMETER (GHOLLO=0,GSOLID=1,GPATTR=2,GHATCH=3)

C    highlighting
        INTEGER GNORML,GHILIT
        PARAMETER (GNORML=0,GHILIT=1)

C    input device status
        INTEGER GNONE,GOK,GNPICK
        PARAMETER (GNONE=0,GOK=1,GNPICK=2)

C    input class
        INTEGER GNCLAS,GLOCAT,GSTROK,GVALUA,GCHOIC,GPICK,GSTRIN
        PARAMETER (GNCLAS=0,GLOCAT=1,GSTROK=2,GVALUA=3,GCHOIC=4)
        PARAMETER (GPICK=5,GSTRIN=6)

C    implicit regeneration mode
        INTEGER GSUPPD,GALLOW
        PARAMETER (GSUPPD=0,GALLOW=1)

C    level of GKS
        INTEGER GL0A,GL0B,GL0C,GL1A,GL1B,GL1C,GL2A,GL2B,GL2C
        PARAMETER (GL0A=0,GL0B=1,GL0C=2,GL1A=3,GL1B=4,GL1C=5)
        PARAMETER (GL2A=6,GL2B=7,GL2C=8)

C    new frame action necessary
        INTEGER GNO,GYES
        PARAMETER (GNO=0,GYES=1)

C    operating mode
        INTEGER GREQU,GSAMPL,GEVENT
        PARAMETER (GREQU=0,GSAMPL=1,GEVENT=2)

C    operating state value
        INTEGER GGKCL,GGKOP,GWSOP,GWSAC,GSGOP
        PARAMETER (GGKCL=0,GGKOP=1,GWSOP=2,GWSAC=3,GSGOP=4)

C    presence of invalid values
        INTEGER GABSNT,GPRSNT
        PARAMETER (GABSNT=0,GPRSNT=1)

C     regeneration flag
        INTEGER GSUPP,GPERFO
        PARAMETER (GSUPP=0,GPERFO=1)

C    relative input priority
        INTEGER GHIGHR,GLOWER
        PARAMETER (GHIGHR=0,GLOWER=1)

C    simultaneous events flag
        INTEGER GNMORE,GMORE
        PARAMETER (GNMORE=0,GMORE=1)

C    text alignment horizontal
        INTEGER GAHNOR,GALEFT,GACENT,GARITE
        PARAMETER (GAHNOR=0,GALEFT=1,GACENT=2,GARITE=3)

C    text alignment vertical
        INTEGER GAVNOR,GATOP,GACAP,GAHALF,GABASE,GABOTT
        PARAMETER (GAVNOR=0,GATOP=1,GACAP=2,GAHALF=3,GABASE=4,GABOTT=5)

C    text path
        INTEGER GRIGHT,GLEFT,GUP,GDOWN
        PARAMETER (GRIGHT=0,GLEFT=1,GUP=2,GDOWN=3)

C    text precision
        INTEGER GSTRP,GCHARP,GSTRKP
        PARAMETER (GSTRP=0,GCHARP=1,GSTRKP=2)

C    type of returned values
        INTEGER GSET,GREALI
        PARAMETER (GSET=0,GREALI=1)

C    update state
        INTEGER GNPEND,GPEND
        PARAMETER (GNPEND=0,GPEND=1)

C    vector/raster/other type
        INTEGER GVECTR,GRASTR,GOTHWK
        PARAMETER (GVECTR=0,GRASTR=1,GOTHWK=2)

C    visibility
        INTEGER GINVIS,GVISI
        PARAMETER (GINVIS=0,GVISI=1)

C    workstation category
        INTEGER GOUTPT,GINPUT,GOUTIN,GWISS,GMO,GMI
        PARAMETER (GOUTPT=0,GINPUT=1,GOUTIN=2,GWISS=3,GMO=4,GMI=5)

C    workstation state
        INTEGER GINACT,GACTIV
        PARAMETER (GINACT=0,GACTIV=1)

C    list of GDP attributes
        INTEGER GPLBND,GPMBND,GTXBND,GFABND
        PARAMETER (GPLBND=0,GPMBND=1,GTXBND=2,GFABND=3)

C    line type
        INTEGER GLSOLI,GLDASH,GLDOT,GLDASD
        PARAMETER (GLSOLI=1,GLDASH=2,GLDOT=3,GLDASD=4)

C    marker type
        INTEGER GPOINT,GPLUS,GAST,GOMARK,GXMARK
        PARAMETER (GPOINT=1,GPLUS=2,GAST=3,GOMARK=4,GXMARK=5)

C
C   GKS functions - These names are used for error handling. The names
C                   are the same as the GKS function names except that
C                   the sentinel character 'G' is replaced by 'E'.

        INTEGER EOPKS,ECLKS,EOPWK,ECLWK
        PARAMETER (EOPKS =  0,ECLKS =  1,EOPWK =  2,ECLWK =  3)
        INTEGER EACWK,EDAWK,ECLRWK
        PARAMETER (EACWK =  4,EDAWK =  5,ECLRWK=  6)

        INTEGER ERSGWK,EUWK,ESDS,EMSG
        PARAMETER (ERSGWK=  7,EUWK  =  8,ESDS  =  9,EMSG  = 10)
        INTEGER EESC,EPL,EPM
        PARAMETER (EESC  = 11,EPL   = 12,EPM   = 13)

        INTEGER ETX,EFA,ECA,EGDP
        PARAMETER (ETX   = 14,EFA   = 15,ECA   = 16,EGDP   = 17)
        INTEGER ESPLI,ESLN,ESLWSC
        PARAMETER (ESPLI = 18,ESLN  = 19,ESLWSC= 20)

        INTEGER ESPLCI,ESPMI,ESMK,ESMKSC
        PARAMETER (ESPLCI= 21,ESPMI = 22,ESMK  = 23,ESMKSC= 24)
        INTEGER ESPMCI,ESTXI,ESTXFP
        PARAMETER (ESPMCI= 25,ESTXI = 26,ESTXFP= 27)

        INTEGER ESCHXP,ESCHSP,ESTXCI,ESCHH
        PARAMETER (ESCHXP= 28,ESCHSP= 29,ESTXCI= 30,ESCHH = 31)
        INTEGER ESCHUP,ESTXP,ESTXAL
        PARAMETER (ESCHUP= 32,ESTXP = 33,ESTXAL= 34)

        INTEGER ESFAI,ESFAIS,ESFASI,ESFACI
        PARAMETER (ESFAI = 35,ESFAIS= 36,ESFASI= 37,ESFACI= 38)
        INTEGER ESPA,ESPARF,ESASF
        PARAMETER (ESPA  = 39,ESPARF= 40,ESASF = 41)

        INTEGER ESPKID,ESPLR,ESPMR,ESTXR
        PARAMETER (ESPKID= 42,ESPLR = 43,ESPMR = 44,ESTXR = 45)
        INTEGER ESFAR,ESPAR,ESCR
        PARAMETER (ESFAR = 46,ESPAR = 47,ESCR  = 48)

        INTEGER ESWN,ESVP,ESVPIP,ESELNT
        PARAMETER (ESWN  = 49,ESVP  = 50,ESVPIP= 51,ESELNT= 52)
        INTEGER ESCLIP,ESWKWN,ESWKVP
        PARAMETER (ESCLIP= 53,ESWKWN= 54,ESWKVP= 55)

        INTEGER ECRSG,ECLSG,ERENSG,EDSG
        PARAMETER (ECRSG = 56,ECLSG = 57,ERENSG= 58,EDSG  = 59)
        INTEGER EDSGWK,EASGWK,ECSGWK
        PARAMETER (EDSGWK= 60,EASGWK= 61,ECSGWK= 62)

        INTEGER EINSG,ESSGT,ESVIS,ESHLIT
        PARAMETER (EINSG = 63,ESSGT = 64,ESVIS = 65,ESHLIT= 66)
        INTEGER ESSGP,ESDTEC,EINLC
        PARAMETER (ESSGP = 67,ESDTEC= 68,EINLC = 69)

        INTEGER EINSK,EINVL,EINCH,EINPK
        PARAMETER (EINSK = 70,EINVL = 71,EINCH = 72,EINPK = 73)
        INTEGER EINST,ESLCM,ESSKM
        PARAMETER (EINST = 74,ESLCM = 75,ESSKM = 76)

        INTEGER ESVLM,ESCHM,ESPKM,ESSTM
        PARAMETER (ESVLM = 77,ESCHM = 78,ESPKM = 79,ESSTM = 80)
        INTEGER ERQLC,ERQSK,ERQVL
        PARAMETER (ERQLC = 81,ERQSK = 82,ERQVL = 83)

        INTEGER ERQCH,ERQPK,ERQST,ESMLC
        PARAMETER (ERQCH = 84,ERQPK = 85,ERQST = 86,ESMLC = 87)
        INTEGER ESMSK,ESMVL,ESMCH
        PARAMETER (ESMSK = 88,ESMVL = 89,ESMCH = 90)

        INTEGER ESMPK,ESMST,EWAIT,EFLUSH
        PARAMETER (ESMPK = 91,ESMST = 92,EWAIT = 93,EFLUSH= 94)
        INTEGER EGTLC,EGTSK,EGTVL
        PARAMETER (EGTLC = 95,EGTSK = 96,EGTVL = 97)

        INTEGER EGTCH,EGTPK,EGTST,EWITM
        PARAMETER (EGTCH = 98,EGTPK = 99,EGTST =100,EWITM =101)
        INTEGER EGTITM,ERDITM,EIITM
        PARAMETER (EGTITM=102,ERDITM=103,EIITM =104)

        INTEGER EEVTM,EACTM
        PARAMETER (EEVTM =105,EACTM=106)
