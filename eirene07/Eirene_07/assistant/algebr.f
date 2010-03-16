C
C
C-----------------------------------------------------------------------
      SUBROUTINE ALGEBR (TERM,OPER,IZIF,CONST,NOP)
C-----------------------------------------------------------------------
C
C     AUTOR:           ST. HUBER
C     IHK-KENNZIFFER:  121
C     DATUM:           25-NOV-1988
C
C     FUNKTION:
C
C     DAS PROGRAMM LIEST ARITHMETISCHE AUSDRUECKE EIN,
C     RUFT DAS UNTERPROGRAMM ZERLEG AUF
C     UND GIBT ENTWEDER DIE EINZELNEN ZERLEGUNGEN ODER DIE
C     REGELVERLETZUNG AUS
C
C-----------------------------------------------------------------------
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
C
C     ARGUMENT:
C
         CHARACTER(*), INTENT(INOUT) :: TERM
C           : EINZULESENDER AUSDRUCK

         CHARACTER(1), INTENT(OUT) :: OPER(*)
         INTEGER, INTENT(OUT) :: IZIF(4,*)
         REAL(DP), INTENT(INOUT) :: CONST(*)
         INTEGER, INTENT(OUT) :: NOP
C
C     KONSTANTENDEKLARATION :
C
         INTEGER, PARAMETER :: ZMAX = 20
C           : ANZAHL DER MAXIMALEN ZERLEGUNGEN

         INTEGER, PARAMETER :: MAXLEN = 72
C           : MAXIMALE STRINGLAENGE

C
C     LOKALE VARIABLEN :
C
         INTEGER ::  LAENGE
C           : AKTUELLE LAENGE VON TERM

         CHARACTER(MAXLEN) :: HLFTERM
C           : HILFSSTRING ZUM UMSPEICHERN

         CHARACTER(MAXLEN+2):: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

         INTEGER :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         INTEGER :: TEIL
C           : AKTUELLE ANZAHL DER ZERLEGUNGEN

         CHARACTER(MAXLEN) :: PART(ZMAX)
C           : FELD VON STRINGS, AUF DENEN DIE EINZELNEN
C             ELEMENTARZERLEGUNGEN FESTGEHALTEN WERDEN

         INTEGER :: IPART(ZMAX)
C           : AKTUELLE LAENGEN VON PART(ZMAX)

         CHARACTER(MAXLEN) :: ARITH(ZMAX)
C           : FELD VON STRINGS, AUF DENEN DIE TEIL-TE GENERATION
C             VON AUSDRU FESTGEHALTEN WIRD

         INTEGER :: IARITH(ZMAX)
C           : AKTUELLE LAENGEN VON ARITH(ZMAX)

         CHARACTER(MAXLEN) :: HILFE
C           : ARBEITSSPEICHER FUER UNTERPROGRAMM ZERLEG

         INTEGER :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN

CHR
CHR      VARIABLEN ZUR MODIFIKATION DES PROGRAMMES
         INTEGER :: NR, ANFANG, ENDE, FELDIND, IK, IKM, IKP
         CHARACTER(4) :: FO
         CHARACTER(12) :: ERSETZ(10)
         character(10) :: buchst
chr
C
C     HILFSVARIABLEN :
C
         INTEGER :: MAXI, I
chr
chr   string, der die neuen variablennamen enthaelt
      buchst='ABCDEFGHIJ'
chr
C
C        LESE TERM UND WERTE AUS
C
         IF (TERM .NE. ' ') THEN
C
C           VERARBEITUNG, FUER DEN FALL, DASS KEINE LEERZEILE
C           EINGELESEN WURDE
C
chr         vorbereiten des terms fuer die weitere verarbeitung, d.h.
chr         bringen der operanden in die vom programm verlangte
chr         zweistellige alphabetische form
            nr=0
101         anfang=index(term,'<')
            if (anfang.ne.0) then
               nr=nr+1
               ende=index(term,'>')
chr            abspeichern der ersetzten operanden
chr            der i-te operand wird hierbei in der i-ten dimension
chr            des feldes ersetz abgelegt und durch den i-ten
chr            buchstaben im alphabet ersetzt. gleichheit von
chr            operanden wird hierbei nicht beruecksichtigt
               ersetz(nr)=term(anfang:ende)
               if (anfang.gt.1) then
               hlfterm=term(1:anfang-1)//buchst(nr:nr)//buchst(nr:nr)//
     .              term(ende+1:)
               else
               hlfterm=buchst(nr:nr)//buchst(nr:nr)//
     .              term(ende+1:)
               endif
               term=hlfterm
               goto 101
            endif
chr
C
C           ERMITTELN DER LAENGE VON TERM
C
            LAENGE=LEN(TERM)
20          IF (TERM(LAENGE:LAENGE) .EQ. ' ') THEN
               LAENGE=LAENGE-1
               GOTO 20
            ENDIF

            AUSDRU=TERM
            AKTLEN=LAENGE

            CALL ZERLEG(AUSDRU, AKTLEN, IPART, PART, IARITH, ARITH,
     >                  TEIL, HILFE, ERROR)

            IF (ERROR .EQ. 0) THEN
C
C              AUSGABE DER ZERLEGUNG
C
               MAXI=0
               DO 45, I=1,TEIL
                   MAXI=MAX(IPART(I),MAXI)
45                 CONTINUE

               NOP=TEIL
               DO 30, I=1,TEIL
chr               ausgabe der zerlegung in der form:
chr                    operator ziffer1 ziffer2 ziffer3 ziffer4
chr               wobei jeweils die 1. und 2. sowie die 3. und 4.
chr               ziffer eine einheit bilden.
chr               entweder stellt eine solche einheit einen operanden
chr               dar oder ein zwischenergebnis mit der 1. ziffer als
chr               nummer und der 2. als 0 zur kennzeichnung des paares
chr               als zwischenergebnis
                  if (part(i)(7:7).ne.'Z') then
                     OPER(I)=PART(I)(9:9)
                     feldind=index(buchst,part(i)(7:7))
                     IK=INDEX(ERSETZ(FELDIND),',')
                     IF (IK.EQ.0) THEN
                       IZIF(1,I)=-I
                       IZIF(2,I)=0
                       CALL RDCN (ERSETZ(FELDIND),CONST(I))
                     ELSE
                       IKM=IK-1
                       IKP=IK+1
                       FO(1:4)='(I )'
                       WRITE (FO(3:3),'(I1)') IK-2
                       READ(ERSETZ(FELDIND)(2:IKM),FO) IZIF(1,I)
                       IF (ERSETZ(FELDIND)(IK+2:IK+2).EQ.'>') THEN
                          READ(ERSETZ(FELDIND)(IKP:IKP),'(I1)')
     .                         IZIF(2,I)
                       ELSEIF (ERSETZ(FELDIND)(IK+3:IK+3).EQ.'>') THEN
                          READ(ERSETZ(FELDIND)(IKP:IKP+1),'(I2)')
     .                         IZIF(2,I)
                       ELSEIF (ERSETZ(FELDIND)(IK+4:IK+4).EQ.'>') THEN
                          READ(ERSETZ(FELDIND)(IKP:IKP+2),'(I3)')
     .                         IZIF(2,I)
                       ENDIF
                     ENDIF

                     IF (PART(I)(10:10).NE.'Z') THEN
                       FELDIND=INDEX(BUCHST,PART(I)(10:10))
                       IK=INDEX(ERSETZ(FELDIND),',')
                       IF (IK.EQ.0) THEN
                         IZIF(3,I)=-I
                         IZIF(4,I)=0
                         CALL RDCN (ERSETZ(FELDIND),CONST(I))
                       ELSE
                         IKM=IK-1
                         IKP=IK+1
                         FO(1:4)='(I )'
                         WRITE (FO(3:3),'(I1)') IK-2
                         READ(ERSETZ(FELDIND)(2:IKM),FO) IZIF(3,I)
                         IF (ERSETZ(FELDIND)(IK+2:IK+2).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP),'(I1)')
     .                          IZIF(4,I)
                         ELSEIF (ERSETZ(FELDIND)(IK+3:IK+3).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP+1),'(I2)')
     .                          IZIF(4,I)
                         ELSEIF (ERSETZ(FELDIND)(IK+4:IK+4).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP+2),'(I3)')
     .                          IZIF(4,I)
                         ENDIF
                       ENDIF
                     ELSE
                        READ(PART(I)(12:12),'(I1)') IZIF(3,I)
                        IZIF(4,I)=0
                     ENDIF
                  ELSE
                     OPER(I)=PART(I)(10:10)
                     READ(PART(I)(9:9),'(I1)') IZIF(1,I)
                     IZIF(2,I)=0
                     if (part(i)(11:11).ne.'Z') then
                       FELDIND=INDEX(BUCHST,PART(I)(11:11))
                       IK=INDEX(ERSETZ(FELDIND),',')
                       IF (IK.EQ.0) THEN
                         IZIF(3,I)=-I
                         IZIF(4,I)=0
                         CALL RDCN (ERSETZ(FELDIND),CONST(I))
                       ELSE
                         IKM=IK-1
                         IKP=IK+1
                         FO(1:4)='(I )'
                         WRITE (FO(3:3),'(I1)') IK-2
                         READ(ERSETZ(FELDIND)(2:IKM),FO) IZIF(3,I)
                         IF (ERSETZ(FELDIND)(IK+2:IK+2).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP),'(I1)')
     .                          IZIF(4,I)
                         ELSEIF (ERSETZ(FELDIND)(IK+3:IK+3).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP+1),'(I2)')
     .                          IZIF(4,I)
                         ELSEIF (ERSETZ(FELDIND)(IK+4:IK+4).EQ.'>') THEN
                           READ(ERSETZ(FELDIND)(IKP:IKP+2),'(I3)')
     .                          IZIF(4,I)
                         ENDIF
                       ENDIF
                     else
                        READ(PART(I)(13:13),'(I1)') IZIF(3,I)
                        IZIF(4,I)=0
                     endif
                  endif
30                CONTINUE
            ELSE
C
C              AUSGABE DER FEHLERMELDUNG
C
               WRITE(iunout,'(2A)') '0FOLGENDE REGELVERLETZUNG ',
     >                        'WURDE ERKANNT:'
               CALL MECKER(ERROR)
               NOP=0
            ENDIF
         ENDIF
C
C     ENDE VON ALGEBR
C
      RETURN
      END
