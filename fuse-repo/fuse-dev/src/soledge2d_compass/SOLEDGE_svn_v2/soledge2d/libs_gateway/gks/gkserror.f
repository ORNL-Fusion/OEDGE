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

        SUBROUTINE GECLKS
C               emergency close GKS

        INTEGER STATE
        INTEGER I, IDUM, IWKID, JERR, NMAX

C               include GKS symbol definitions
        INCLUDE 'gksdefs.i'

C               first check GKS operating state
        CALL GQOPS(STATE)
        IF (STATE.EQ.0) GOTO 5

C               actions necessary only when GKS open
        GOTO (1,2,3,4),STATE

C               STATE=GSGOP : there are open segments
C
  4     CALL GCLSG
        CONTINUE

C               STATE=GWSAC : there are active workstations.
C                             deactivate them
  3     CONTINUE
        CALL GQACWK(1,JERR,NMAX,IWKID)
        IF(NMAX.LE.0) GOTO 2
        DO 300 I=1,NMAX
          CALL GQACWK(1,JERR,IDUM,IWKID)
          CALL GDAWK(IWKID)
300     CONTINUE

C               STATE=GWSOP : there are open workstations.
C                             close them
  2     CONTINUE
        CALL GQOPWK(1,JERR,NMAX,IWKID)
        IF(NMAX.LE.0) GOTO 1
        DO 200 I=1,NMAX
          CALL GQOPWK(1,JERR,IDUM,IWKID)
          CALL GCLWK(IWKID)
200     CONTINUE

C               STATE=GGKOP : GKS open
C                             close it
  1     CONTINUE
        CALL GCLKS

  5     CONTINUE

        RETURN
        END



        SUBROUTINE GERHND (ERRNO,FCTID,ERRFIL)
C               error handler

        INTEGER ERRNO,FCTID,ERRFIL

        CALL GERLOG (ERRNO,FCTID,ERRFIL)

        RETURN
        END



        SUBROUTINE GERLOG (ERRNO,FCTID,ERRFIL)
C               error logging
 
        INTEGER ERRNO,FCTID,ERRFIL

        INTEGER BASE,STATUS

        INCLUDE 'gksdefs.i'

        CHARACTER*6 NMS(EOPKS:EACTM)

        DATA NMS(EOPKS),NMS(ECLKS),NMS(EOPWK),NMS(ECLWK),NMS(EACWK)
     &        /'GOPKS ','GCLKS ','GOPWK ','GCLWK ','GACWK '/
        DATA NMS(EDAWK),NMS(ECLRWK)
     &        /'GDAWK ','GCLRWK'/
        DATA NMS(ERSGWK),NMS(EUWK),NMS(ESDS),NMS(EMSG),NMS(EESC)
     &        /'GRSGWK','GUWK  ','GSDS  ','GMSG  ','GESC  '/
        DATA NMS(EPL),NMS(EPM),NMS(ETX),NMS(EFA),NMS(ECA)
     &        /'GPL   ','GPM   ','GTX   ','GFA   ','GCA   '/
        DATA NMS(EGDP),NMS(ESPLI),NMS(ESLN),NMS(ESLWSC),NMS(ESPLCI)
     &        /'GGDP  ','GSPLI ','GSLN  ','GSLWSC','GSPLCI'/
        DATA NMS(ESPMI),NMS(ESMK),NMS(ESMKSC),NMS(ESPMCI),NMS(ESTXI)
     &        /'GSPMI ','GSMK  ','GSMKSC','GSPMCI','GSTXI '/
        DATA NMS(ESTXFP)
     &        /'GSTXFP'/
        DATA NMS(ESCHXP),NMS(ESCHSP),NMS(ESTXCI),NMS(ESCHH),NMS(ESCHUP)
     &        /'GSCHXP','GSCHSP','GSTXCI','GSCHH ','GSCHUP'/
        DATA NMS(ESTXP),NMS(ESTXAL)
     &        /'GSTXP ','GSTXAL'/
        DATA NMS(ESFAI),NMS(ESFAIS),NMS(ESFASI),NMS(ESFACI),NMS(ESPA)
     &        /'GSFAI ','GSFAIS','GSFASI','GSFACI','GSPA  '/
        DATA NMS(ESPARF),NMS(ESASF)
     &        /'GSPARF','GSASF '/
        DATA NMS(ESPKID),NMS(ESPLR),NMS(ESPMR),NMS(ESTXR),NMS(ESFAR)
     &        /'GSPKID','GSPLR ','GSPMR ','GSTXR ','GSFAR '/
        DATA NMS(ESPAR),NMS(ESCR)
     &        /'GSPAR ','GSCR  '/
        DATA NMS(ESWN),NMS(ESVP),NMS(ESVPIP),NMS(ESELNT),NMS(ESCLIP)
     &        /'GSWN  ','GSVP  ','GSVPIP','GSELNT','GSCLIP'/
        DATA NMS(ESWKWN),NMS(ESWKVP)
     &        /'GSWKWN','GSWKVP'/
        DATA NMS(ECRSG),NMS(ECLSG),NMS(ERENSG),NMS(EDSG),NMS(EDSGWK)
     &        /'GCRSG ','GCLSG ','GRENSG','GDSG  ','GDSGWK'/
        DATA NMS(EASGWK),NMS(ECSGWK)
     &        /'GASGWK','GCSGWK'/
        DATA NMS(EINSG),NMS(ESSGT),NMS(ESVIS),NMS(ESHLIT),NMS(ESSGP)
     &        /'GINSG ','GSSGT ','GSVIS ','GSHLIT','GSSGP '/
        DATA NMS(ESDTEC),NMS(EINLC)
     &        /'GSDTEC','GINLC '/
        DATA NMS(EINSK),NMS(EINVL),NMS(EINCH),NMS(EINPK),NMS(EINST)
     &        /'GINSK ','GINVL ','GINCH ','GINPK ','GINST '/
        DATA NMS(ESLCM),NMS(ESSKM)
     &        /'GSLCM ','GSSKM '/
        DATA NMS(ESVLM),NMS(ESCHM),NMS(ESPKM),NMS(ESSTM),NMS(ERQLC)
     &        /'GSVLM ','GSCHM ','GSPKM ','GSSTM ','GRQLC '/
        DATA NMS(ERQSK),NMS(ERQVL)
     &        /'GRQSK ','GRQVL '/
        DATA NMS(ERQCH),NMS(ERQPK),NMS(ERQST),NMS(ESMLC),NMS(ESMSK)
     &        /'GRQCH ','GRQPK ','GRQST ','GSMLC ','GSMSK '/
        DATA NMS(ESMVL),NMS(ESMCH)
     &        /'GSMVL ','GSMCH '/
        DATA NMS(ESMPK),NMS(ESMST),NMS(EWAIT),NMS(EFLUSH),NMS(EGTLC)
     &        /'GSMPK ','GSMST ','GWAIT ','GFLUSH','GGTLC '/
        DATA NMS(EGTSK),NMS(EGTVL)
     &        /'GGTSK ','GGTVL '/
        DATA NMS(EGTCH),NMS(EGTPK),NMS(EGTST),NMS(EWITM),NMS(EGTITM)
     &        /'GGTCH ','GGTPK ','GGTST ','GWITM ','GGTITM'/
        DATA NMS(ERDITM),NMS(EIITM)
     &        /'GRDITM','GIITM '/
        DATA NMS(EEVTM),NMS(EACTM)
     &        /'GEVTM ','GACTM '/

        DATA BASE /140148736/

        CALL GERSET (ERRNO,ERRFIL)

        STATUS = BASE + ERRNO*8
        CALL LIBSIG (STATUS, NMS(FCTID)(1:6))

        RETURN
        END
