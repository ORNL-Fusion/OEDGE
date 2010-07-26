c     -*-Fortran-*-
c
      SUBROUTINE RDG (IREF,GRAPH,IOPT,IERR)
      IMPLICIT  NONE
      INTEGER   IREF,IOPT,IERR
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG  : READ IN A ROW OF GRAPH DETAILS                            *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
      integer len,lenstr
      external lenstr 

C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
      if (buffer(2:4).eq.'000') goto 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
C
      MESAGE = 'EXPECTING CHARACTER STRING AND 1 INTEGER'
      if (buffer(1:3).eq.'000') goto 100
c
c     Make sure buffer isn't a blank line or single character
c
      len = lenstr(buffer)  
      if (len.le.2) goto 100
c
c      write(0,'(a)') 'BUFFER:',buffer(1:len),':'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,IOPT
      READ (GRAPH(1:3),'(I3)',ERR=9999,END=9999) IREF
      RETURN
C
 9998 IERR = 1
      RETURN
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
      END
C
C
C
      SUBROUTINE RDFN(GRAPH,STRING,iseld,IERR)
      use error_handling
      IMPLICIT  NONE
      INTEGER   Iseld,IERR
      CHARACTER GRAPH*(*),STRING*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDFN  : READS IN A FILENAME AND AN INTEGER OPTION                *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDFN'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 2 CHARACTER STRINGS AND AN INTEGER'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,string,iseld
      RETURN
C
 9998 IERR = 1
      RETURN
C
 9999 IERR = 1
      call errmsg('IOOUT:RDFN:','ERROR READING '//
     >    graph(1:len_trim(graph))//' '//
     >    mesage(1:len_trim(mesage))//' '//
     >    ' :LAST LINE READ:'//buffer(1:len_trim(buffer)))
c    
c      WRITE (7,'(1X,2A,3(/1X,A))')
c     >  'RDFN: ERROR READING ',GRAPH,MESAGE,
c     >  'LAST LINE READ :-',BUFFER
c
      RETURN
      END
c
c
c
      SUBROUTINE RDFN_MULTICASE(GRAPH,CMD,NAME,start,increment,
     >                          ncases,IERR)
      IMPLICIT  NONE
      INTEGER   start,increment,ncases,IERR
      CHARACTER GRAPH*(*),CMD*(*),name*(*)
    
C
C  *********************************************************************
C  *                                                                   *
C  *  RDFN_MULTICASE  : READS IN A FILENAME DEFINITION FOR MULTI-CASE  *
C  *                    PLOTS                                          *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDFNMC'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 3 CHARACTER STRINGS AND 3 INTEGERS'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,cmd,name,
     >                                start,increment,ncases
      RETURN
C
 9998 IERR = 1
      RETURN
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDFN_MULTICASE: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDFN_CMD(GRAPH,CMD,NAME,iopt,IERR)
      IMPLICIT  NONE
      INTEGER   iopt,IERR
      CHARACTER GRAPH*(*),CMD*(*),name*(*)
    
C
C  *********************************************************************
C  *                                                                   *
C  *  RDFN_CMD  : READS IN A FILENAME, COMMAND AND INTEGER OPTION      *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDFNRW'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 3 CHARACTER STRINGS AND 1 INTEGER'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,cmd,name,iopt
      RETURN
C
 9998 IERR = 1
      RETURN
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDFN_CMD: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
C
C
C
      SUBROUTINE RDFN_BOLO(GRAPH,STRING,iseld,iflag,IERR)
      IMPLICIT  NONE
      INTEGER   Iseld,iflag,IERR
      CHARACTER GRAPH*(*),STRING*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDFN_BOLO  : READS IN A FILENAME AND AND 2 INTEGER OPTIONS       *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDFN'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 2 CHARACTER STRINGS AND 2 INTEGERS'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,string,iseld,iflag
      RETURN
C
 9998 IERR = 1
      RETURN
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDFN_BOLO: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG_XSECTION(GRAPH,r1p,z1p,r2p,z2p,npts,iselect,
     >                        istate,iexpt,iavg,ierr)
      IMPLICIT  NONE
      real      r1p,z1p,r2p,z2p 
      INTEGER   Iselect,IERR,istate,npts,iexpt,iavg
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG_XSECTION  : 2D CONTOUR CROSS-SECTION PLOT INPUT              *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG_XSECTION'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRINGS 4 REALS AND 5 INTEGERS'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,r1p,z1p,r2p,z2p,
     >                               npts,iselect,istate,iexpt,iavg
      RETURN

 9998 ierr = 1
       
      return 
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG_XSECTION: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      WRITE (6,'(1X,2A,3(/1X,A))')
     >  'RDG_XSECTION: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      WRITE (0,'(1X,2A,3(/1X,A))')
     >  'RDG_XSECTION: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG_contopts(GRAPH,icntr,ncntr,uconts,maxpts,
     >                        xcen,ycen,
     >                        xnear,ynear,ierr)
      IMPLICIT  NONE
      INTEGER   icntr,ncntr,ierr,in,maxpts
      real xcen,ycen,xnear,ynear,uconts(maxpts)
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG_CONTOPTS  : LOAD CONTOUR OPTIONS AND PLOT RANGE              *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A)') BUFFER,'RDG_CONTOPTS'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING + 4 REALS +'//
     >         ' 2 INTS + N REALS'
c
      READ (BUFFER,*,ERR=9997,END=9997) GRAPH,xcen,ycen,
     >                         xnear,ynear,icntr,ncntr,
     >                         (uconts(in),in=1,ncntr)

 9997 continue

      RETURN

 9998 ierr = 1
       
      return 
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG_CONTOPTS: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG_CONTOUR(GRAPH,iselect,
     >                   istate,iexpt,optval,ierr)
      IMPLICIT  NONE
      INTEGER   Iselect,IERR,istate,ifact,iexpt,icntr,ncntr
      real optval  
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG_CONTOUR  : LOAD BASE PARAMETERS FOR GENERALIZED CONTOUR PLOTS*
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A)') BUFFER,'RDG_CONTOUR'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING AND 3 INTEGERS AND 1 REAL'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,iselect,istate,
     >                         iexpt,optval
      RETURN

 9998 ierr = 1
       
      return 
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG_CONTOUR: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG_3I1R(GRAPH,iselect,
     >                   istate,iexpt,optval,ierr)
      IMPLICIT  NONE
      INTEGER   Iselect,IERR,istate,ifact,iexpt,icntr,ncntr
      real optval  
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG_3I1R  :   READS 3 INTEGER PARAMETERS AND 1 REAL              *
C  *                -LOAD PARAMETERS FOR COMBINED ERO/DEP PLOTS        * 
C  *                -LOAD BASE PARAMETERS FOR GENERALIZED CONTOUR PLOTS*
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A)') BUFFER,'RDG_CONTOUR'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING AND 3 INTEGERS AND 1 REAL'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,iselect,istate,
     >                         iexpt,optval
      RETURN

 9998 ierr = 1
       
      return 
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG_3I1R: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG_ring(GRAPH,iselect,
     >                   istate,iexpt,minfrac,maxfrac,
     >                   axis_type,plot_type,ierr)
      IMPLICIT  NONE
      INTEGER   Iselect,IERR,istate,iexpt
      integer   axis_type,plot_type
      real      minfrac,maxfrac
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG_RING : LOAD BASE PARAMETERS FOR GENERALIZED ALONG RING PLOTS *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A)') BUFFER,'RDG_CONTOUR'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING AND 3 INTEGERS, 2 REALS '
     >         //'AND 3 INTEGERS'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,iselect,istate,
     >                    iexpt,minfrac,maxfrac,axis_type,plot_type
      RETURN

 9998 ierr = 1
       
      return 
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG_RING: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG_LOS(GRAPH,npts,nlines,iselect,
     >                   istate,iexpt,iaxis,iavg,ifact,optval,ierr)
      IMPLICIT  NONE
      INTEGER   Iselect,IERR,istate,npts,nlines,iaxis,iavg,ifact,iexpt
      real optval  
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG_LOS  : LOAD BASE PARAMETERS FOR GENERALIZED LOS PLOTS        *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG_LOS'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING, 8 INTEGERS AND 1 REAL'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,iselect,istate,
     >                         npts,nlines,iaxis,
     >                         iexpt,iavg,ifact,optval
      RETURN

 9998 ierr = 1
       
      return 
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG_LOS: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG_819(GRAPH,mindist,maxdist,shift_dist,
     >                   scale_min,scale_max,scale_factor,iexpt,ierr)
      IMPLICIT  NONE
      INTEGER   iexpt,ierr
      real mindist,maxdist,shift_dist,scale_min,scale_max,scale_factor
      
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG_LOS  : LOAD BASE PARAMETERS FOR GENERALIZED LOS PLOTS        *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG_819'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING, 8 INTEGERS AND 1 REAL'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,
     >   mindist,maxdist,shift_dist,scale_min,scale_max,scale_factor,
     >   iexpt
      RETURN

 9998 ierr = 1
       
      return 
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG_819: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG_EXPT(GRAPH,plotid,nexpt,maxexpt,
     >                    expt_ds,expt_col,ierr)
      IMPLICIT  NONE
      INTEGER   plotid,nexpt,maxexpt,expt_ds(*),expt_col(*),ierr
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG_EXPT  : LOADS THE INDICES FOR NEXPT EXPERIMENTAL DATA SETS   *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
      integer   in
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG_LOS'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING, 2 INTEGERS AND'//
     >         ' NEXPT INTEGER PAIRS'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,plotid,nexpt,
     >       ((expt_ds(in),expt_col(in)),in=1,min(nexpt,maxexpt))

      nexpt = min(nexpt,maxexpt)  
c
      RETURN

 9998 ierr = 1
       
      return 
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG_LOS: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG_LP_PLOTDATA(GRAPH,iexpt,lp_plot_type,lp_plot_avg,
     >                    expt_axis_offset,
     >                    lp_axis_offset,axis_callibration,ierr)
      IMPLICIT  NONE
      INTEGER   iexpt,lp_plot_type,lp_plot_avg,ierr
      real      expt_axis_offset,
     >          lp_axis_offset,axis_callibration
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG_EXPT  : LOADS THE INDICES FOR NEXPT EXPERIMENTAL DATA SETS   *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
      integer   in
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A20)') BUFFER,'RDG_LP_PLOTDATA'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 2
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING, 3 INTEGERS AND'//
     >         ' 3 REALS'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,iexpt,lp_plot_type,
     >       lp_plot_avg,expt_axis_offset,
     >       lp_axis_offset,axis_callibration

c
      RETURN

 9998 ierr = 1
       
      return 
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG_LP_PLOTDATA: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG_LOS3D(GRAPH,iselect,istate,iexpt,iaxis,
     >                     minsteps,stepsize,optval,ierr)
      IMPLICIT  NONE
      INTEGER   Iselect,IERR,istate,iaxis,minsteps,iexpt
      real optval,stepsize  
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG_LOS3D  : LOAD BASE PARAMETERS FOR 3D IMAGE LOS PLOTS         *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A)') BUFFER,'RDG_LOS3D'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING, 5 INTEGERS AND 1 REAL'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,iselect,istate,
     >                         iaxis,
     >                         iexpt,minsteps,stepsize,optval
      RETURN

 9998 ierr = 1
       
      return 
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG_LOS3D: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG_CAMERA(GRAPH,xres,yres,position,direction,
     >                      upvec,rightvec,lookat,camera,ierr)
      IMPLICIT  NONE
      integer xres,yres,ierr
      real*8 position(3),direction(3),upvec(3),rightvec(3),lookat(3)
      CHARACTER GRAPH*(*),camera*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG_CAMERA  : LOADS A CAMERA VIEW IN POVRAY STYLE SPECIFICATION  *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
      integer i,len,lenstr
      external lenstr
c
      ierr = 0
C
c-----------------------------------------------------------------------
c     Read camera pixel resolution 
c-----------------------------------------------------------------------
c
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A)') BUFFER,'RDG_CAM1'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING AND 2 INTEGERS'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,xres,yres
      len = lenstr(graph) 
      CAMERA = graph(5:)
C
c----------------------------------------------------------------
c     Read camera position
c----------------------------------------------------------------
c
  200 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A)') BUFFER,'RDG_CAM2'
      IF (BUFFER(1:1).EQ.'$') GOTO 200
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING AND 3 REALS'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,(position(i),i=1,3)
C
c----------------------------------------------------------------
c     Read camera direction vector 
c     - magnitude is used to define view area
c----------------------------------------------------------------
c
  300 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A)') BUFFER,'RDG_CAM3'
      IF (BUFFER(1:1).EQ.'$') GOTO 300
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING AND 3 REALS'
c
      READ (BUFFER,*,ERR=9999,END=9999)GRAPH,(direction(i),i=1,3)
C
c----------------------------------------------------------------
c     Read camera UP vector 
c     - used to define vertical range of viewing region
c----------------------------------------------------------------
c
  400 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A)') BUFFER,'RDG_CAM4'
      IF (BUFFER(1:1).EQ.'$') GOTO 400
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING AND 3 REALS'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,(upvec(i),i=1,3)
C
c----------------------------------------------------------------
c     Read camera RIGHT vector 
c     - used to define horizontal range of viewing region
c----------------------------------------------------------------
c
  500 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A)') BUFFER,'RDG_CAM5'
      IF (BUFFER(1:1).EQ.'$') GOTO 500
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING AND 3 REALS'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,(rightvec(i),i=1,3)
C
c----------------------------------------------------------------
c     Read camera LOOKAT vector 
c     - used to define actual location where the camera looks
c----------------------------------------------------------------
c
  600 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A)') BUFFER,'RDG_CAM6'
      IF (BUFFER(1:1).EQ.'$') GOTO 600
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING AND 3 REALS'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,(lookat(i),i=1,3)
c
c     Return
c
      RETURN

C
 9998 continue
c
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG_CAMERA: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG_REAL_ARRAY(GRAPH,data,maxpts,npts,ierr)
      IMPLICIT  NONE
      INTEGER   npts,ierr,maxpts
      real data(maxpts) 
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG_REAL_ARRAY  : READS A 1D ARRAY OF REAL DATA ITEMS            *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
      integer in,ios
c
      ierr = 0
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9997,END=9998,IOSTAT=ios) 
     >                  BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG_RA'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c     If it is not a '000' designation 
c
      if (buffer(2:4).ne.'000') then
         ierr = 1
         rewind(5)
         return
      endif
C
      MESAGE = 'EXPECTING 1 CHARACTER STRING AND A REAL DATA ARRAY'
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,npts
c
      if (npts.le.maxpts) then 
c
         READ (BUFFER,*,ERR=9999,END=9999) GRAPH,npts,
     >                (data(in),in=1,npts)
c
      else
c 
         MESAGE = 'NUMBER OF DATA ITEMS TOO LARGE FOR ARRAY'
         goto 9999
c
      endif
c
      RETURN

 9997 ierr = 1
      WRITE (6,'(1X,3A,i6,3(/1X,A))')
     >  'RDG_REAL_ARRAY: ERROR LOADING ',GRAPH,' ERRNUM = ',
     >  ios,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      WRITE (7,'(1X,3A,i5,3(/1X,A))')
     >  'RDG_REAL_ARRAY: ERROR LOADING ',GRAPH,' ERRNUM = ',
     >  ios,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      return

 9998 ierr = 1
      WRITE (6,'(1X,3A,i6,3(/1X,A))')
     >  'RDG_REAL_ARRAY: EOF READING ',GRAPH,' ERRNUM = ',
     >  ios,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      WRITE (7,'(1X,3A,i5,3(/1X,A))')
     >  'RDG_REAL_ARRAY: EOF READING ',GRAPH,' ERRNUM = ',
     >  ios,MESAGE,
     >  'LAST LINE READ :-',BUFFER
       
      return 
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG_REAL_ARRAY: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      WRITE (6,'(1X,2A,3(/1X,A))')
     >  'RDG_REAL_ARRAY: ERROR READING ',GRAPH,MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG2 (GRAPH,ROBS,ZOBS,THEMIN,DTHE,drtheres,NUMTHE,
     >                 IZMIN,IZMAX,AVPTS,NUMSMOOTH,ATYPE,IERR)
      implicit none
      INTEGER   NUMTHE,IZMIN,IZMAX,AVPTS,NUMSMOOTH,ATYPE,IERR
      REAL      ROBS,ZOBS,THEMIN,DTHE,DRtheres
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG2 : READ IN A SECOND ROW OF GRAPH DETAILS FOR LOS PLOTS       *
C  *         This routine can read in plot details to support both     *
c  *         kinds of LOS plots. The reason this is required is that   *
c  *         JET grid format is not available for all of the           *
c  *         geometries that can be run using DIVIMP and thus support  *
c  *         for the other methods is still required.                  *
c  *         The variable drtheres represents the deltaR value for     *
c  *         the INTLOS type of plot and the theres value for the      *
c  *         LOSINT type of plot.                                      *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG2'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
C
      MESAGE = 'EXPECTING 1 CHAR, 5 REALS AND 6 INTEGERS'
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,ROBS,ZOBS,THEMIN,DTHE,
     >      DRtheres,NUMTHE,AVPTS,NUMSMOOTH,IZMIN,IZMAX,ATYPE
      RETURN
C
 9998 IERR = 1
      RETURN
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
      END
C
C
C
      SUBROUTINE RDG2A (GRAPH,ROBS,ZOBS,THEMIN,DTHE,NUMTHE,THERES,
     >                 IZMIN,IZMAX,AVPTS,NUMSMOOTH,ATYPE,IERR)
      implicit none
      INTEGER   NUMTHE,IZMIN,IZMAX,AVPTS,NUMSMOOTH,ATYPE,IERR
      REAL      ROBS,ZOBS,THEMIN,DTHE,THERES
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG2A: READ IN A SECOND ROW OF GRAPH DETAILS FOR LOS PLOTS       *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG2'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
C
      MESAGE = 'EXPECTING 1 CHAR, 5 REALS AND 6 INTEGERS'
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,ROBS,ZOBS,THEMIN,DTHE,
     >      NUMTHE,THERES,AVPTS,NUMSMOOTH,IZMIN,IZMAX,ATYPE
      RETURN
C
 9998 IERR = 1
      RETURN
C
 9999 IERR = 1
      WRITE (7,'(1X,A,4(/1X,A))')
     >  'RDG2: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
      END
C
C
C
      SUBROUTINE RDG3 (GRAPH,ROBS,ZOBS,CX1,CY1,CZ1,CX2,CY2,CZ2,
     > NUMTHE,THERES,IZMIN,IZMAX,AVPTS,NUMSMOOTH,ATYPE,IERR)
      implicit none
      INTEGER   NUMTHE,IZMIN,IZMAX,AVPTS,NUMSMOOTH,ATYPE,IERR
      REAL      ROBS,ZOBS,CX1,CY1,CZ1,CX2,CY2,CZ2,THERES
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG3 : READ IN 2 ROWS OF GRAPH DETAILS FOR 3D ARRAYS             *
C  *                                                                   *
C  *********************************************************************
C
      CHARACTER MESAGE*72,BUF1*72,BUF2*72
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 READ (5,'(A72)',ERR=9998,END=9998) BUF1
      WRITE (9,'(1X,A72,1X,A6)') BUF1,'RDG3'
      IF (BUF1(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUF1(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUF1)
        GOTO 100
      ENDIF
c
  101 READ (5,'(A72)',ERR=9998,END=9998) BUF2
      WRITE (9,'(1X,A72,1X,A6)') BUF2,'RDG3'
      IF (BUF2(1:1).EQ.'$') GOTO 101
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUF2(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUF2)
        GOTO 101
      ENDIF
C
      MESAGE = 'EXPECTING 1 CHAR AND 9 REALS'
      READ (BUF1,*,ERR=9996,END=9996) GRAPH,ROBS,ZOBS,CX1,CY1,CZ1,
     >      CX2,CY2,CZ2,THERES
      MESAGE = 'EXPECTING 1 CHAR AND 6 INTEGERS'
      READ (BUF2,*,ERR=9997,END=9997) GRAPH,
     >      NUMTHE,AVPTS,NUMSMOOTH,IZMIN,IZMAX,ATYPE
      RETURN
C
 9998 IERR = 1
      RETURN
C
 9996 IERR = 1
      WRITE (7,'(1X,A,4(/1X,A))')
     >  'RDG3: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUF1
      RETURN
C
 9997 IERR = 1
      WRITE (7,'(1X,A,4(/1X,A))')
     >  'RDG3: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUF2
      RETURN
      END
C
C
C
      SUBROUTINE RDG4 (graph,ngrm,nplts,ringnos,maxplts,pltfact,
     >                 ierr)
      IMPLICIT  none
c
      INTEGER   ngrm,nplts,maxplts,ringnos(maxplts),ierr
      real      pltfact
      character*(*) graph
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG4 : READ IN A SECOND ROW OF GRAPH DETAILS FOR MULTI PLOTS     *
c  *         This reads in two integers followed by an array           *
c  *         of integers. The first integer is the number of plots     *
c  *         on the page, the second the total number of plots and     *
c  *         theother integers are a list fo the ring numbers for      *
c  *         which plots are to be produced.			       *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
      integer i,npltstmp
c
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG4'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
C
c     NEED TO PRE-READ to check the validity of NPLTS 
c     before reading the list of ring numbers. 
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,ngrm,nplts
c
      if (nplts.gt.maxplts) then 
         nplts = maxplts
         MESAGE=' NPLTS > MAXPLTS '
         WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG4: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
         write (7,'(a)') ' NPLTS SET EQUAL TO MAXPLTS '
      endif 
c
      MESAGE = 'EXPECTING 2 Integers + 1 REAL + Nt integers'
c
c     Now READ in the list of ring numbers to the now 
c     possibly adjusted value of NPLTS
c
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,ngrm,npltstmp,
     >                        pltfact, (ringnos(i),i=1,nplts)
c
      RETURN
C
 9998 IERR = 1
      RETURN
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG4: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDG5 (GRAPH,ADASID,ADASYR,ADASEX,
     >                 ISELE,ISELR,ISELX,ISELD,IZ,ZION,IERR)
      implicit none
      INTEGER   ISELE,ISELR,ISELX,ISELD,IERR,ADASYR,IZ,ZION
      CHARACTER GRAPH*(*), ADASID*(*),ADASEX*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG5 : READ IN SELECTOR SWITCHES FOR ADAS PLRP CALCULATIONS      *
c  *         INCLUDING CHARGE STATE AND ATOMIC NUMBER                  *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
C
      IERR = 0
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG5'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
C
      MESAGE = 'EXPECTING 2 CHAR, 1 INT, 1 CHAR  AND 6 INTEGERS'
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,ADASID,ADASYR,ADASEX,
     >                                  ISELE,ISELR,ISELX,ISELD,
     >                                  IZ,ZION
      RETURN
C
 9998 IERR = 1
      RETURN
C
 9999 IERR = 1
      WRITE (7,'(1X,A,4(/1X,A))')
     >  'RDG1: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
      END
C
C
C
      SUBROUTINE RDG6 (GRAPH,interp_opt,iseld,axis_offset_r,
     >                 axis_offset_z,nsets,plotpts,maxnpts,
     >                 maxndata,ierr)
      implicit none
      INTEGER   interp_opt,iseld,nsets,in,maxnpts,maxndata,ierr
      REAL      plotpts(maxnpts,maxndata)
      real      axis_offset_r,axis_offset_z
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG6 : READ IN N ROWS OF GRAPH DETAILS FOR THOMPSON PLOTS        *
C  *                                                                   *
C  *********************************************************************
C
      CHARACTER MESAGE*72,buf1*132,buf2*132
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 READ (5,'(A132)',ERR=9998,END=9998) BUF1
      WRITE (9,'(1X,A72,1X,A6)') BUF1,'RDG6'
      IF (BUF1(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUF1(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUF1)
        GOTO 100
      ENDIF
C
      MESAGE = 'EXPECTING 1 CHAR 2 INT 2 REAL 1 INT'
      READ (BUF1,*,ERR=9996,END=9996) GRAPH,interp_opt,iseld,
     >      axis_offset_r,axis_offset_z,nsets
c
      write(mesage,'(''EXPECTING '',i6,'' SETS OF R,Z COORDINATE'',
     >    '' PAIRS (ONE SET/LINE)'')') nsets
c
      do in = 1, nsets 
c
  101    READ (5,'(A132)',ERR=9998,END=9998) BUF2
         WRITE (9,'(1X,A72,1X,A6)') BUF2,'RDG6'
         IF (BUF2(1:1).EQ.'$') GOTO 101
c
         if (buf2(2:4).ne.'000') goto 9997 
c       
         READ (BUF2,*,ERR=9997,END=9997) graph,plotpts(in,1),
     >                                         plotpts(in,2)
c
      end do         

      RETURN
C
 9998 IERR = 1
      RETURN
C
 9996 IERR = 1
      WRITE (7,'(1X,A,4(/1X,A))')
     >  'RDG6: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUF1
      RETURN
C
 9997 IERR = 1
      WRITE (7,'(1X,A,4(/1X,A))')
     >  'RDG6: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUF2
      RETURN
      END
c
c     Krieger IPP/98  
c
      SUBROUTINE rdg7(graph, minscale,maxscale,localcngs,ierr)
      IMPLICIT  none
      real   minscale,maxscale
      INTEGER   localcngs,ierr, iref
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG7: read in extra line with contour/false color plot details   *
C  *                                                                   *
C  *********************************************************************
C
      include 'reader'
      CHARACTER MESAGE*72
      character*(*) graph
c
      localcngs = 0
      minscale = 0.
      maxscale = 0.
C
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(A72)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG7'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
      if (buffer(2:4).eq.'000') then
        backspace 5
        return
      endif
C
      MESAGE = 'EXPECTING 2 REALS AND 1 INTEGER'
      if (buffer(2:4).eq.'001') then
        READ (BUFFER,*,ERR=9999,END=9999) graph,
     >                                    minscale,maxscale,localcngs
      else
        backspace 5
      endif
      RETURN
C
 9998 IERR = 1
      RETURN
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDG7: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c
c
      SUBROUTINE RDGRCP(graph,r1p,z1p,r2p,z2p,int_type,
     >                  exp_ds,exp_offsets,exp_dataopt,
     >                  exp_vcalcopt,exp_tcalcopt,exp_param,
     >                  ierr)
      implicit none
      real r1p,z1p,r2p,z2p,exp_param,exp_offsets(4)
      integer int_type,exp_ds,exp_vcalcopt,exp_tcalcopt,ierr
      integer exp_dataopt
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG2 : READ IN A SECOND ROW OF GRAPH DETAILS FOR LOS PLOTS       *
C  *         This routine can read in plot details to support both     *
c  *         kinds of LOS plots. The reason this is required is that   *
c  *         JET grid format is not available for all of the           *
c  *         geometries that can be run using DIVIMP and thus support  *
c  *         for the other methods is still required.                  *
c  *         The variable drtheres represents the deltaR value for     *
c  *         the INTLOS type of plot and the theres value for the      *
c  *         LOSINT type of plot.                                      *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
      integer in
c
      ierr = 0
C
c     Read first line - osm probe parameters
c
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDGRCP'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
C
      MESAGE = 'EXPECTING 1 CHAR, 4 REALS, 1 INTEGER ON FIRST LINE'
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,r1p,z1p,r2p,z2p,
     >                                  int_type
c
c     Read second line - rcp probe parameters
c
      MESAGE = 'END OF FILE ON UNIT 5'
  200 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDGRCP'
      IF (BUFFER(1:1).EQ.'$') GOTO 200
c
c     read in as many offsets as are present - no errors allowed
c
      MESAGE = 'EXPECTING 1 CHAR, 4 INTEGERS, 2 REAL ON SECOND LINE'
      READ (BUFFER,*,ERR=9997,END=9997) GRAPH,exp_ds,exp_dataopt,
     >                                  exp_vcalcopt,exp_tcalcopt,
     >                                  exp_param,
     >                                  (exp_offsets(in),in=1,4)

9997  continue 
c
      RETURN
C
 9998 IERR = 1
      RETURN
C
 9999 IERR = 1
      WRITE (7,'(1X,2A,3(/1X,A))')
     >  'RDGRCP: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',
     >               BUFFER
      RETURN
      END
c
c psmod
c
      SUBROUTINE RDG579(GRAPH,ZB,Z,MB,M,NB,TBG,TB,VBG,XPER,FILED,IERR)
 
      implicit none
      INTEGER   ZB,Z,MB,M,XPER,IERR,FILED
      REAL      TBG,VBG,NB,TB
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG579: READ IN DATA FOR REISER/DIVIMP PLOTS OF THE FRICTIONAL   *
C  *          AND THERMAL FORCES PLOTTED VERSUS CHI                    *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
C
      IERR = 0
*      WRITE(0,*)'HI AGAIN'
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(A132)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG579'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
C
      MESAGE = 'EXPECTING 1 CHAR, 5 INT, AND 3 REALS'
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,ZB,Z,MB,M,NB,TBG,TB,VBG,
     >                                  XPER,FILED
      RETURN
C
 9998 IERR = 1
      RETURN
C
 9999 IERR = 1
      WRITE (7,'(1X,A,4(/1X,A))')
     >  'RDG1: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
      END
c
c psmod
c
c
      SUBROUTINE GET (TITLE,NIZS,JOB,equil,
     >                FACTA,FACTB,ITER,NITERS)
      use subgrid
      IMPLICIT  NONE
C     INCLUDE   "PARAMS"
      include 'params'
      CHARACTER*(*) TITLE,JOB,equil
      INTEGER   NIZS,ITER,NITERS
      REAL      FACTA(-1:MAXIZS),FACTB(-1:MAXIZS)
C
C  *********************************************************************
C  *                                                                   *
C  *   GET:  READ  RESULTS OF DIVIMP RUN ON UNFORMATTED UNIT 8         *
C  *         GET MAY BE CALLED SEVERAL TIMES, ONCE FOR                 *
C  *         EACH ITERATION.                                           *
C  *                                   CHRIS FARRELL    MARCH 1989     *
C  *********************************************************************
C
C     INCLUDE   "COMTOR"
      include 'comtor'
C     INCLUDE   "CNEUT2"
      include 'cneut2'
C     INCLUDE   "CGEOM"
      include 'cgeom'
C     INCLUDE   "DYNAM2"
      include 'dynam2'
C     INCLUDE   "DYNAM3"
      include 'dynam3'
C     INCLUDE   "DYNAM4"
      include 'dynam4'
C     INCLUDE   "PINDATA"
      include 'pindata'
c
      include 'cadas'
c
      include 'outxy'
C
      include 'grbound'
c
      include 'cedge2d' 
C
      include 'transcoef'
c
      include 'cioniz'
c
      include 'promptdep'
c
      include 'reiser_com' 
      include 'line_profile'
c
      include 'hc_global_opts'
c
      CHARACTER VERSE*5
C
C     Temporary title strings
C
      character*80  tmptitle1  
      character*174 tmptitle2
c
c     These should be placed in a common block ideally so 
c     that changes to these variables in DIVIMP are automatically
c     included in OUT.
c
      character*174 tmpdesc
      character*1024 tmpdesc2
c   
      INTEGER   IR,IZ,IT,in,ik
c
c     Version control - calculate a unique and always increasing
c     version code.
c     Maximum revison number for a given version number is maxrev-1 
c
      integer   maxrev,version_code
      parameter (maxrev=100)
c slmod begin 
      INCLUDE 'diagvel'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      INTEGER i1,i2,i3,i4,idum1,idum2,idum3,idum4,idum5
      REAL    slver

c...  Crap, but needed for backward compatability:
      INTEGER s1,v1,maxasc3

      INTEGER    MAXASC2
      PARAMETER (MAXASC2=5000)
      COMMON /ASCCOM/  asc_cell2   ,asc_link2  ,asc_nvp2   ,asc_rvp2   ,
     .                 asc_grid2   ,asc_region2,asc_vol2   ,
     .                 asc_zvp2   ,
     .                 ascvertex2  ,ascnvertex2,ascdata2   ,
     .                 asc_link3D2 ,asc_nvp3D2 ,asc_zmin3D2,asc_zmax3D2,
     .                 asc_xvp3D 2 ,asc_yvp3D2 ,asc_zvp3D2 
      INTEGER
     .  asc_cell2  (MAXASC2),
     .  asc_region2(MAXASC2),
     .  asc_link2(4,MAXASC2),
     .  asc_grid2(2,MAXASC2),
     .  asc_nvp2   (MAXASC2),
     .  ascnvertex2(MAXASC2),
     .  asc_link3D2(6,MAXASC2),
     .  asc_nvp3D2 (  MAXASC2)
      REAL
     .  asc_rvp2  ( 8,MAXASC2),
     .  asc_zvp2  ( 8,MAXASC2),
     .  asc_vol2     (MAXASC2),
     .  ascvertex2(40,MAXASC2),
     .  ascdata2  (MAXASC2,5),
     .  asc_zmin3D2(   MAXASC2),
     .  asc_zmax3D2(   MAXASC2),
     .  asc_xvp3D2(6,8,MAXASC2),
     .  asc_yvp3D2(6,8,MAXASC2),
     .  asc_zvp3D2(6,8,MAXASC2)

      REAL pinasd2 (2000 ,MAXASD2,MAXASS,3),eirpuff2(20,MAXPUFF),
     .     pinasd3 (20000,MAXASD2,MAXASS,2),eirpuff3(MAXNAS,10),
     .     ascdata3(20000,5),pinstrata2(MAXNKS,MAXNRS,3,6)


      INTEGER MAXNDS_,MAXPTS_

c...TEMP
      pincode = -1

      eirzaa = -99.0

c slmod end
C
c
c     Add Steve's requested rewind
c
      rewind(8) 
c
c     Read in version string to define output file version - then rewind
c     if necessary and read the rest. 
c
      read(8) verse
      read(verse,'(i1,2x,i2)') vernum,revnum
c
      version_code = vernum * maxrev + revnum
c
      write(6,*) 'DIVIMP VERSION = ', vernum,
     >           ' REVISION = ',revnum,
     >           ' VERSION CODE = ',version_code
C
c slmod begin
      IF (version_code.GE.6*maxrev+35) THEN
        MAXNDS_ = MAXNDS
      ELSE
c...    Based on MAXNDS value as of January 29, 2004 for version 6A/34:
        MAXNDS_ = 110
        IF (MAXNDS.NE.110) 
     .    WRITE(0,*) 'WARNING: MAXNDS problems for DEPS, NEROS, '//
     .               'ROMPRDEPS, TARGFLUXDATA, and KPREDBAR'
      ENDIF
      IF (version_code.GE.6*maxrev+36) THEN
        MAXPTS_ = MAXPTS
      ELSE
        IF (MAXNRS.EQ.60) THEN
c...      TEMPORARY: THIS DETECTS MY THESIS CODE, WHERE MAXPTS WAS NOT CHANGED.
          MAXPTS_ = MAXPTS
        ELSE
c...      Based on MAXPTS value as of February 26, 2004 for version 6A/35:
          MAXPTS_ = 300
          IF (MAXPTS.NE.300) 
     .      WRITE(0,*) 'WARNING: MAXPTS problems in IOOUT'
        ENDIF
      ENDIF
c slmod end
      if (version_code.ge.6*maxrev+6) then 
c
c
c        Basic information - title, job description,
c                            equilibrium, shot number and time slice
c
c        Added description - not used in OUT so only loaded locally
c                            for now - can propagate later if req'd 

         if (version_code.ge.6*maxrev+23) then 
            READ  (8) tmpTITLE2,tmpdesc2,JOB,EQUIL,ISHOT,TSLICE
         elseif (version_code.ge.6*maxrev+10) then 
            READ  (8) tmpTITLE2,tmpdesc,JOB,EQUIL,ISHOT,TSLICE
         else
            READ  (8) tmpTITLE2,JOB,EQUIL,ISHOT,TSLICE
         endif
c
c
c        Simulation values
c

         read(8) ITER,NITERS,NIZS,NTS,CRMB,CRMI,CIZB,CION,IMODE,
     >           NITERSOL
c
c        Geometry values 
c
         read(8) R0,Z0,RXP,ZXP,NRS,
     >           RMIN,RMAX,ZMIN,ZMAX,DR,DZ,NDS,NDSIN,NXS,NYS,
     >           IRCENT,IRSEP,IRWALL,IRTRAP,IKT,IKREF,IKTO,IKTI,
     >           irsep2,irwall2,nrs2,ndsin2,ndsin3,irtrap2,
     >           nves,nvesm,nvesp,inmid,oumid,refct,CIRHR,NPOLYP,
     >           cxnear,cynear
c
         CALL IINOUT ('R NKS    ',nks ,maxnrs)
c
c        Scaling factors  
c
         if (version_code.ge.6*maxrev+32) then 
            read(8) QTIM,FSRATE,ABSFAC,absfac_neut,CBPHI,CK0
         else
            read(8) QTIM,FSRATE,ABSFAC,CBPHI,CK0
         endif 
c
         CALL RINOUT ('R FACTA  ',facta ,maxizs+2)
         CALL RINOUT ('R FACTB  ',factb ,maxizs+2)
         CALL RINOUT ('R DWELTS ',dwelts,maxizs+2)
c slmod begin
         IF (IMODE.EQ.1) THEN
           CALL RINOUT ('R DWELFS ',dwelfs,maxnts)
         ELSE
           CALL RINOUT ('R DWELFS ',dwelfs,1     )       
         ENDIF
c
c         CALL RINOUT ('R DWELFS ',dwelfs,maxnts)
c slmod end
         CALL RINOUT ('R KALPHS ',kalphs,maxizs)
         CALL RINOUT ('R KBETAS ',kbetas,maxizs)
c
c        DIVIMP Options
c
         read(8) CIOPTF,cdatopt,cneur,cgridopt,cre2d,cre2dizs,
     >           xygrid
c
c        Result numbers
c
         read(8) nleakcore,cleaksn
c
c        ADAS       
c
         read(8) useridh,iyearh,useridz,iyearz
c
c
c        The following are the piece-wise wall specifications
c        from the GA15 routines.
c
         READ(8) ioncpts,ionwpts
c
c        Wall definitions
c
         CALL RINOUT ('R RIW    ',riw ,maxpts_)
         CALL RINOUT ('R ZIW    ',ziw ,maxpts_)
         CALL RINOUT ('R RCW    ',rcw ,maxpts_)
         CALL RINOUT ('R ZCW    ',zcw ,maxpts_)
         CALL RINOUT ('R RW     ',rw  ,maxpts_)
         CALL RINOUT ('R ZW     ',zw  ,maxpts_)
         CALL RINOUT ('R RVES   ',rves,maxpts_)
         CALL RINOUT ('R ZVES   ',zves,maxpts_)
c
c        GA15 workspace 
c
         CALL IINOUT ('R IWINDW ',iwindw,2*maxpts_)
         CALL RINOUT ('R IWWORK ',iwwork,4*maxpts_)
         CALL RINOUT ('R IWTDUM ',iwtdum,maxpts_)
         CALL RINOUT ('R IWXDUM ',iwxdum,maxpts_)
         CALL RINOUT ('R IWYDUM ',iwydum,maxpts_)
c
         CALL IINOUT ('R ICINDW ',icindw,2*maxpts_)
         CALL RINOUT ('R ICWORK ',icwork,4*maxpts_)
         CALL RINOUT ('R ICTDUM ',ictdum,maxpts_)
         CALL RINOUT ('R ICXDUM ',icxdum,maxpts_)
         CALL RINOUT ('R ICYDUM ',icydum,maxpts_)
c
c
c        The following relate to the wall definition as 
c        calculated in DIVIMP.
c
         READ(8) wlwall1,wlwall2,wltrap1,wltrap2,wallpts
c
         title = tmptitle2 
c
      else
c
c        For older raw files - rewind and re-read the entire first
c        record.
c
         rewind(8)
c
         READ(8) VERSE,ITER,NITERS,NIZS,tmpTITLE1,JOB,NTS,R0,Z0,RXP,ZXP,
     >          NRS,(NKS(IR),IR=1,MAXNRS),RMIN,RMAX,ZMIN,ZMAX,DR,DZ,NDS,
     >          (FACTA(IZ),FACTB(IZ),DWELTS(IZ),IZ=-1,MAXIZS),
     >          (DWELFS(IT),IT=1,MAXNTS),QTIM,FSRATE,NXS,NYS,ABSFAC,
     >          CRMB,CRMI,CIZB,CION,IMODE,IRCENT,IRSEP,IRWALL,IRTRAP,
     >          IKT,NDSIN,(KALPHS(IZ),KBETAS(IZ),IZ=1,MAXIZS),IKREF,
     >          CK0,CIRHR,IKTO,IKTI,NPOLYP,ISHOT,TSLICE,EQUIL,
     >          NITERSOL,CIOPTF,CBPHI,cdatopt,cneur,nves,nleakcore,
     >          irsep2,irwall2,nrs2,ndsin2,ndsin3,cgridopt,
     >          irtrap2,xygrid,cleaksn,cxnear,cynear,refct,
     >          nvesm,nvesp,inmid,oumid,useridh,iyearh,useridz,iyearz,
     >          cre2d,cre2dizs,
c
c
c               The following are the piece-wise wall specifications
c               from the GA15 routines.
c
     >          ioncpts,ionwpts,(riw(in),ziw(in),rcw(in),zcw(in),
     >          rw(in),zw(in),rves(in),zves(in),in=1,maxpts_),
     >          ((iwindw(it,in),icindw(it,in),it=1,2),in=1,maxpts_),
     >          (iwwork(in),icwork(in),in=1,4*maxpts_),
     >          (iwtdum(in),iwxdum(in),iwydum(in),in=1,maxpts_),
     >          (ictdum(in),icxdum(in),icydum(in),in=1,maxpts_),
c
c               The following relate to the wall definition as 
c               calculated in DIVIMP.
c
     >          wlwall1,wlwall2,wltrap1,wltrap2,wallpts
c
         title = tmptitle1
c
      endif
c
C
c      read(verse,'(i1,2x,i2)') vernum,revnum
c
c      version_code = vernum * maxrev + revnum
c
      WRITE (6,9001) NXS,NYS,NRS,NDS,NIZS,NTS,
     >  MAXNXS,MAXNYS,MAXNRS,MAXNDS,MAXIZS,MAXNTS,
     >  TITLE,JOB,ITER,IMODE,refct,maxseg,nvesm,nvesp
c
      CALL RINOUT ('R POWLS ',POWLS ,MAXNKS*MAXNRS*(MAXIZS+2))
      CALL RINOUT ('R LINES ',LINES ,MAXNKS*MAXNRS*(MAXIZS+2))
      CALL RINOUT ('R HPOWLS',HPOWLS,MAXNKS*MAXNRS*2)
      CALL RINOUT ('R HLINES',HLINES,MAXNKS*MAXNRS*2)
      CALL RINOUT ('R TIZS  ',TIZS  ,MAXNKS*MAXNRS*(MAXIZS+2))
      CALL RINOUT ('R ZEFFS ',ZEFFS ,MAXNKS*MAXNRS*3)
      CALL RINOUT ('R WALLS ',WALLS ,MAXNKS*MAXNRS*(MAXIZS+2))
c...  MAXNDS repair required:
      CALL RINOUT ('R DEPS  ',DEPS  ,MAXNDS_*MAXIZS)
      CALL RINOUT ('R NEROS ',NEROS ,MAXNDS_*5)
c
      if (version_code.ge.(5*maxrev+11)) then 
c...     MAXNDS repair required:
         CALL RINOUT ('R PRDEPS',PROMPTDEPS,MAXNDS_*6)
      elseif (version_code.ge.(5*maxrev+8)) then 
         CALL RINOUT ('R PRDEPS',PROMPTDEPS,MAXNDS_)
      endif  
c
      if (version_code.ge.(5*maxrev+14)) then 
         CALL RINOUT ('R WALLSN',WALLSN,MAXPTS_+1)
         CALL RINOUT ('R WALLSE',WALLSE,MAXPTS_+1)
         if (version_code.ge.(6*maxrev+24)) then 
            CALL RINOUT ('R WALLSEI',WALLSE_I,MAXPTS_+1)
         endif
         CALL RINOUT ('R WALLSI',WALLSI,MAXPTS_+1)
      else 
         CALL RINOUT ('R WALLSN',WALLSN,MAXPTS_)
         CALL RINOUT ('R WALLSE',WALLSE,MAXPTS_)
      endif

      if (version_code.ge.(6*maxrev+38)) then 
        CALL RINOUT ('R WALLPT',WALLPT,MAXPTS_*32)
      elseif (version_code.ge.(6*maxrev+1)) then 
        CALL RINOUT ('R WALLPT',WALLPT,MAXPTS_*25)
      else
        CALL RINOUT ('R WALLPT',WALLPT,MAXPTS_*19)
      endif
c
      CALL RINOUT ('R RS    ',RS    ,MAXNKS*MAXNRS)
      CALL RINOUT ('R ZS    ',ZS    ,MAXNKS*MAXNRS)
c
c     More geometry data
c
      if (version_code.ge.(6*maxrev+3)) then 
c
         CALL RINOUT ('R KSB   ',KSB   ,(MAXNKS+1)*MAXNRS)
         CALL RINOUT ('R KPB   ',KPB   ,(MAXNKS+1)*MAXNRS)
c
      endif
c
c     The storing of these arrays needed to be customized
c     because of a likely size mismatch between the
c     array in DIVIMP and that in OUT.
c
c      CALL IINOUT ('R IKXYS ',IKXYS ,MAXNXS*MAXNYS)
c      CALL IINOUT ('R IRXYS ',IRXYS ,MAXNXS*MAXNYS)
c      CALL IINOUT ('R IFXYS ',IFXYS ,MAXNXS*MAXNYS)
c
      CALL IINOUT2 ('R IKXYS ',IKXYS ,MAXNXS,MAXNYS,MAXIXS,MAXIYS)
      CALL IINOUT2 ('R IRXYS ',IRXYS ,MAXNXS,MAXNYS,MAXIXS,MAXIYS)
      CALL IINOUT2 ('R IFXYS ',IFXYS ,MAXNXS,MAXNYS,MAXIXS,MAXIYS)
c
      CALL IINOUT ('R IKDS  ',IKDS  ,MAXNDS_)
      CALL IINOUT ('R IRDS  ',IRDS  ,MAXNDS_)
C
      CALL IINOUT ('R IKINS ',IKINS ,MAXNKS*MAXNRS)
      CALL IINOUT ('R IKOUTS',IKOUTS,MAXNKS*MAXNRS)
      CALL IINOUT ('R IRINS ',IRINS ,MAXNKS*MAXNRS)
      CALL IINOUT ('R IROUTS',IROUTS,MAXNKS*MAXNRS)
C
      CALL IINOUT ('R KORY  ',KORY  ,MAXNKS*MAXNRS)
      CALL IINOUT ('R KORPG ',KORPG ,MAXNKS*MAXNRS)
      CALL IINOUT ('R NVERTP',NVERTP,MAXNKS*MAXNRS)
      CALL RINOUT ('R RVERTP',RVERTP,5*MAXNKS*MAXNRS)
      CALL RINOUT ('R ZVERTP',ZVERTP,5*MAXNKS*MAXNRS)
C
      CALL RINOUT ('R SDLIMS',SDLIMS,MAXNKS*MAXNRS*(MAXIZS+2))
      CALL RINOUT ('R SDTS  ',SDTS  ,MAXNKS*MAXNRS*(MAXIZS+2))
      CALL RINOUT ('R ELIMS ',ELIMS ,MAXNKS*3*(MAXIZS+2))
      CALL RINOUT ('R WALKS ',WALKS ,MAXNWS*2)
c
      call rinout ('R CHEM D',chemden,maxnks*maxnrs)
      call rinout ('R CHEMIZ',chemizs,maxnks*maxnrs)
C
      CALL RINOUT ('R KKS   ',KKS   ,MAXNRS)
      CALL RINOUT ('R KSS   ',KSS   ,MAXNKS*MAXNRS)
      CALL RINOUT ('R KPS   ',KPS   ,MAXNKS*MAXNRS)
      CALL RINOUT ('R KNORMS',KNORMS,MAXNKS*MAXNRS)
      CALL RINOUT ('R KPERPS',KPERPS,MAXNKS*MAXNRS)
      CALL RINOUT ('R KCURVS',KCURVS,MAXNKS*MAXNRS)
      CALL RINOUT ('R KVOLS ',KVOLS ,MAXNKS*MAXNRS)
      CALL RINOUT ('R KAREAS',KAREAS,MAXNKS*MAXNRS)
      CALL RINOUT ('R KTOTAS',KTOTAS,MAXNRS)
      CALL RINOUT ('R KTOTVS',KTOTVS,MAXNRS)
      CALL RINOUT ('R KVOL2 ',KVOL2 ,MAXNKS*MAXNRS)
      CALL RINOUT ('R KAREA2',KAREA2,MAXNKS*MAXNRS)
      CALL RINOUT ('R KTOTA2',KTOTA2,MAXNRS)
      CALL RINOUT ('R KTOTV2',KTOTV2,MAXNRS)
      CALL RINOUT ('R KBFS  ',KBFS  ,MAXNKS*MAXNRS)
      CALL RINOUT ('R KINS  ',KINS  ,MAXNKS*MAXNRS)
      CALL RINOUT ('R KSMAXS',KSMAXS,MAXNRS)
      CALL RINOUT ('R KPMAXS',KPMAXS,MAXNRS)
      CALL RINOUT ('R KTEBS ',KTEBS ,MAXNKS*MAXNRS)
      CALL RINOUT ('R KTIBS ',KTIBS ,MAXNKS*MAXNRS)
      CALL RINOUT ('R KNBS  ',KNBS  ,MAXNKS*MAXNRS)
c
      CALL RINOUT ('R KFIZS ',KFIZS ,MAXNKS*MAXNRS*(MAXIZS+1))
c
c      CALL RINOUT ('R KFSSMO',KFSSMOD,MAXNKS*MAXNRS)
c
      CALL RINOUT ('R KINDS ',KINDS ,MAXNKS*MAXNRS)
      CALL RINOUT ('R KOUTDS',KOUTDS,MAXNKS*MAXNRS)
      CALL RINOUT ('R KFORDS',KFORDS,MAXNKS*MAXNRS)
      CALL RINOUT ('R KBACDS',KBACDS,MAXNKS*MAXNRS)

c
      if (version_code.ge.(6*maxrev+39)) then 
c
c        Load more geometry data
c
         CALL RINOUT ('R COSALI',COSALI,MAXNKS*MAXNRS)
         CALL RINOUT ('R COSALO',COSALO,MAXNKS*MAXNRS)
         CALL RINOUT ('R DISTIN',DISTIN,MAXNKS*MAXNRS)
         CALL RINOUT ('R DISTOU',DISTOUT,MAXNKS*MAXNRS)
c
      endif

c
      CALL RINOUT ('R OKTEBS',OKTEBS,MAXNKS*MAXNRS)
      CALL RINOUT ('R OKTIBS',OKTIBS,MAXNKS*MAXNRS)
      CALL RINOUT ('R OKNBS ',OKNBS ,MAXNKS*MAXNRS)
      CALL RINOUT ('R OKVHS ',OKVHS ,MAXNKS*MAXNRS)
      CALL RINOUT ('R OKES  ',OKES  ,MAXNKS*MAXNRS)
C
      CALL RINOUT ('R KFEGS ',KFEGS ,MAXNKS*MAXNRS)
      CALL RINOUT ('R KFIGS ',KFIGS ,MAXNKS*MAXNRS)
      CALL RINOUT ('R KES   ',KES   ,MAXNKS*MAXNRS)
      CALL RINOUT ('R KVHS  ',KVHS  ,MAXNKS*MAXNRS)
c
      if (version_code.ge.(6*maxrev+7)) then 
c
c     Average force arrays 
c
         if (version_code.ge.(6*maxrev+18)) then 
c
c           Coulomb logarithm for the Drift-Kinetic Model: CIOPTR
c 
            READ(8) Coulomb_log
c
         endif
c
         CALL RINOUT ('R Fcell ',Fcell ,MAXNKS*MAXNRS*MAXIZS)
         CALL RINOUT ('R Fthi  ',Fthi  ,MAXNKS*MAXNRS*MAXIZS)
         CALL RINOUT ('R Ffi   ',Ffi   ,MAXNKS*MAXNRS*MAXIZS)
         CALL RINOUT ('R Fvbg  ',Fvbg  ,MAXNKS*MAXNRS*MAXIZS)
c
         if (version_code.ge.(6*maxrev+18)) then 
            CALL RINOUT ('R DIFF  ',DIFF  ,MAXNKS*MAXNRS*MAXIZS)
         endif
c
         CALL RINOUT ('R VELavg',VELavg,MAXNKS*MAXNRS*MAXIZS)
      endif
c
c     Background data at plates
c
      CALL RINOUT ('R RP    ',RP    ,MAXNDS_)
      CALL RINOUT ('R ZP    ',ZP    ,MAXNDS_)
      call iinout ('R IDDS  ',idds  ,maxnrs*2)
c
      if (version_code.ge.(6*maxrev+11)) then  

         call rinout ('R PSITAR',psitarg ,maxnrs*2)

      endif
c
      CALL RINOUT ('R KTEDS ',KTEDS ,MAXNDS_)
      CALL RINOUT ('R KTIDS ',KTIDS ,MAXNDS_)
      CALL RINOUT ('R KTI3LS',KTI3LS,MAXNDS_)
      CALL RINOUT ('R KTINJ ',KTINJ ,MAXNDS_)
      CALL RINOUT ('R KNDS  ',KNDS  ,MAXNDS_)
      CALL RINOUT ('R KFEDS ',KFEDS ,MAXNDS_)
      CALL RINOUT ('R KFIDS ',KFIDS ,MAXNDS_)
      CALL RINOUT ('R KEDS  ',KEDS  ,MAXNDS_)
      CALL RINOUT ('R KVDS  ',KVDS  ,MAXNDS_)
c
      if (version_code.ge.(6*maxrev+9)) then 
c...     MAXNDS repair required:
         call rinout ('R HEATF ',targfluxdata,(MAXNDS_+3)*4*4)      
      endif  
c
c...  MAXNDS repair required:
      call rinout ('R KPREDB',kpredbar,MAXNDS_*3*2)
c
      CALL RINOUT ('R KFLUX ',KFLUX ,MAXNDS_)
      CALL RINOUT ('R KENER ',KENER ,MAXNDS_)
      CALL RINOUT ('R KYIELD',KYIELD,MAXNDS_)
      CALL RINOUT ('R KFY   ',KFY   ,MAXNDS_)
      CALL RINOUT ('R KRMAX ',KRMAX ,MAXNDS_)
      CALL RINOUT ('R KCUM  ',KCUM  ,MAXNDS_)
      CALL RINOUT ('R DDS   ',DDS   ,MAXNDS_)
      CALL RINOUT ('R THETAS',THETAS,MAXNDS_)
      CALL RINOUT ('R DDS2  ',DDS2  ,MAXNDS_)
      CALL RINOUT ('R THETA2',THETAS2,MAXNDS_)
      CALL RINOUT ('R COSTET',COSTET,MAXNDS_)
c
      call rinout ('R RHOG  ',rhog  ,maxnrs*maxnks)
      call rinout ('R THETAG',thetag,maxnrs*maxnks)
      call rinout ('R HRO   ',hro   ,maxnrs*maxnks)
      call rinout ('R HTETA ',hteta ,maxnrs*maxnks)
      call rinout ('R BTS   ',bts   ,maxnrs*maxnks)
c
      if (version_code.ge.(5*MAXREV+12)) then 
         call rinout ('R PSIFL ',psifl ,maxnrs*maxnks)
      endif 
c
      call iinout ('R TAGDV ',tagdv ,maxnrs*maxnks)
c
c     Chi Squared data
c
      CALL RINOUT ('R CHISQ1',SCHISQ1,25)
      CALL RINOUT ('R CHISQ2',SCHISQ2,25)
      CALL RINOUT ('R CHISQ3',SCHISQ3,25)
      CALL RINOUT ('R CHISQ4',SCHISQ4,25)
      CALL RINOUT ('R CHISQ5',SCHISQ5,25)
C
      CALL RINOUT ('R PINATO',PINATOM ,MAXNKS*MAXNRS)
      CALL RINOUT ('R PINION',PINION  ,MAXNKS*MAXNRS)
      CALL RINOUT ('R PINALP',PINALPHA,MAXNKS*MAXNRS)
      CALL RINOUT ('R PINMOL',PINMOL  ,MAXNKS*MAXNRS)
      CALL RINOUT ('R PINZ0 ',PINZ0   ,MAXNKS*MAXNRS)
      CALL RINOUT ('R PININZ',PINIONZ ,MAXNKS*MAXNRS)
      CALL RINOUT ('R PINENA',PINENA  ,MAXNKS*MAXNRS)
      CALL RINOUT ('R PINENM',PINENM  ,MAXNKS*MAXNRS)
      CALL RINOUT ('R PINENZ',PINENZ  ,MAXNKS*MAXNRS)
      CALL RINOUT ('R PINQI ',PINQI   ,MAXNKS*MAXNRS)
      CALL RINOUT ('R PINQE ',PINQE   ,MAXNKS*MAXNRS)
      CALL RINOUT ('R PINMP ',PINMP   ,MAXNKS*MAXNRS)
      CALL RINOUT ('R PINVDI',PINVDIST,3*14*MAXNKS*MAXNRS)
      CALL RINOUT ('R PINREC',PINREC  ,MAXNKS*MAXNRS)
c 
      if (version_code.ge.(5*maxrev+13)) then 
c
         CALL RINOUT ('R PININF',PINIZ_INFO,MAXNRS*4)
c
      endif  
C
      call rinout ('R RVESM ',RVESM   ,2*MAXSEG)
      call rinout ('R ZVESM ',ZVESM   ,2*MAXSEG)
      call iinout ('R JVESM ',JVESM   ,MAXSEG)
      call rinout ('R FLUXHW',FLUXHW  ,MAXSEG)
      call rinout ('R FLXHW2',FLXHW2  ,MAXSEG)
      call rinout ('R FLXHW3',FLXHW3  ,MAXSEG)
      call rinout ('R FLXHW4',FLXHW4  ,MAXSEG)
      call rinout ('R FLXHW5',FLXHW5  ,MAXSEG)

      if (version_code.ge.(5*maxrev+12)) then
         call rinout ('R FLXHW6',FLXHW6  ,MAXSEG)
      endif
      if (version_code.ge.(6*maxrev+9)) then
         call rinout ('R FLXHW7',FLXHW7  ,MAXSEG)
      endif
      if (version_code.ge.(6*maxrev+12)) then
         call rinout ('R FLXHW8',FLXHW8  ,MAXSEG)
      endif
c
      CALL RINOUT ('R HWALKS',HWALKS  ,MAXNWS*2)
C
      CALL rINOUT ('R SOLTE ',solte,maxnks*msolpt+msolpt+1)
      CALL rINOUT ('R SOLTI ',solti,maxnks*msolpt+msolpt+1)
      CALL rINOUT ('R SOLNE ',solne,maxnks*msolpt+msolpt+1)
      CALL rINOUT ('R SOLVEL',solvel,maxnks*msolpt+msolpt+1)
      CALL rINOUT ('R SOLCOR',solcor,maxnks*msolpt+msolpt+1)
C
c     Leakage data
c
      CALL RINOUT ('R CLEAKS',cleaks  ,Maxpts_)
      CALL RINOUT ('R CLEAKN',cleakn  ,Maxpts_*(maxizs+1))
c
c      write (0,*) 'MAXIMP:',maximp 
c
      call rinout ('R LEAKPS',cleakpos,maximp*2) 
c
c     More arrays related to leakage results 
c
      call rinout ('R ncore ',ncore,maxnks*maxnrs)
      call rinout ('R nedge ',nedge,maxnks*maxnrs)
      call rinout ('R ntrap ',ntrap,maxnks*maxnrs)
      call rinout ('R ndivt ',ndivert,maxnks*maxnrs)
      call rinout ('R nmsol ',nmsol,maxnks*maxnrs)
      call rinout ('R WTSOU ',wtsource,maxpts_*maxnrs*4*5)
c
      if (version_code.ge.(6*maxrev+24)) then 
         call rinout ('R WTDEP ',wtdep,maxpts_*(maxpts_+1)*3)
      elseif (version_code.ge.(6*maxrev+23)) then 
         call rinout ('R WTDEP ',wtdep,maxpts_*(maxpts_+1))
      endif
c    
      call rinout ('R TSOUR ',targsrc,3*4)
      call rinout ('R TLEAK ',targleak,3*4)
      call rinout ('R WSOUR ',wallsrc,5*3)
      call rinout ('R WLEAK ',wallleak,5*3)

c
c     READ any EDGE2D BG data that has been saved 
c

      if (cre2d.eq.1.or.cre2d.eq.2.or.cre2d.eq.3.or.cre2d.eq.5) then 

        call rinout ('R E2D N ',e2dnbs,maxnks*maxnrs)
        call rinout ('R E2D TE',e2dtebs,maxnks*maxnrs)
        call rinout ('R E2D TI',e2dtibs,maxnks*maxnrs)
        call rinout ('R E2D VB',e2dvhs,maxnks*maxnrs)
        call rinout ('R E2D E ',e2des,maxnks*maxnrs)
        call rinout ('R E2D I ',e2dion,maxnks*maxnrs)
        if (version_code.ge.(5*maxrev+9)) then 
           call rinout ('R E2D A ',e2datom,maxnks*maxnrs)
        endif
        call rinout ('R E2D TA',e2dtarg,maxnrs*8*2)
        call rinout ('R E2D GP',e2dgpara,maxnks*maxnrs)
        call rinout ('R E2D GD',e2dgdown,maxnks*maxnrs)
c
        if (version_code.ge.(5*maxrev+8)) then 
           call rinout ('R E2D G ',e2dflux,(maxnks+1)*maxnrs)
        else
           call rinout ('R E2D G ',e2dflux,maxnks*maxnrs)
        endif
c
        if (version_code.ge.(5*maxrev+8)) then 
           call rinout ('R E2D VE',e2dbvel,(maxnks+1)*maxnrs)
        endif
c
        if (version_code.ge.(5*maxrev+9)) then 
           call rinout ('R E2D Z0',e2dz0,maxnks*maxnrs)
        endif
c
        if (version_code.ge.(5*maxrev+15)) then 
           call rinout ('R E2D RC',e2dhrec,maxnks*maxnrs)
        endif
c
        if (version_code.ge.(5*maxrev+10)) then 
           call rinout ('R E2D RC',e2drec,maxnks*maxnrs)
           call rinout ('R E2D CX',e2dcxrec,maxnks*maxnrs)
        endif
c
        if (cre2dizs.gt.0) then 

          call rinout ('R E2D NZ ',e2dnzs,maxnks*maxnrs*
     >                              (maxe2dizs+1))
          call rinout ('R E2D PW',e2dpowls,maxnks*maxnrs*
     >                              (maxe2dizs+1))
          call rinout ('R E2D LI',e2dlines,maxnks*maxnrs*
     >                              (maxe2dizs+1))

        endif
c
      endif
c
C
c     Read Data related to transport coefficient calculations 
c
      call rinout ('R DPERP ',DPERP,maxnrs)
      call rinout ('R DPERPO',odperp,maxnrs)
      call rinout ('R DPERPI',idperp,maxnrs)
      call rinout ('R XPERP ',xPERPt,maxnrs)
      call rinout ('R XPERPO',oxperpt,maxnrs)
      call rinout ('R XPERPI',ixperpt,maxnrs)
      call rinout ('R XPI   ',chiperpi,maxnrs)
      call rinout ('R XPI  O',ochiperpi,maxnrs)
      call rinout ('R XPI  I',ichiperpi,maxnrs)
      call rinout ('R XPE   ',chiperpe,maxnrs)
      call rinout ('R XPE  O',ochiperpe,maxnrs)
      call rinout ('R XPE  I',ichiperpe,maxnrs)
      call rinout ('R RC OUT',rcouter,maxnrs)
      call rinout ('R RC IN ',rcinner,maxnrs)
      if (version_code.ge.(5*maxrev+13)) then 
         call rinout ('R ZC OUT',zcouter,maxnrs)
         call rinout ('R ZC IN ',zcinner,maxnrs)
         call rinout ('R MIDIST',middist,maxnrs*2)
      endif
      call rinout ('R GRADNE',gradn,maxnks*maxnrs)
      call rinout ('R GRADTE',gradte,maxnks*maxnrs)
      call rinout ('R GRADTI',gradti,maxnks*maxnrs)
      call rinout ('R E2DGNE',e2dgradn,maxnks*maxnrs)
      call rinout ('R E2DGTE',e2dgradte,maxnks*maxnrs)
      call rinout ('R E2DGTI',e2dgradti,maxnks*maxnrs)
c
c     Read in the line profile data if it is present 
c
      if (version_code.ge.(6*maxrev+31)) then 
c
         read(8) line_profile_opt
c
         if (line_profile_opt.ne.0) then 
             read(8) lp_wave,lp_instrument_width,
     >             lp_bin_width,lp_robs,lp_zobs,lp_theta,lp_dtheta
             CALL R8INOUT ('R LP',line_profile,max_lp_bins*2+1)
             CALL R8INOUT ('R MOD_LP',mod_line_profile,max_lp_bins*2+1)
         endif     
c
      endif
c
c     Read the pressure - from SOL option 22 
c
      call rinout ('R KPRESS',kpress,maxnks*maxnrs*2)
c   
      if (version_code.ge.(6*maxrev+2)) then 
c
         call rinout ('R KPRAD',kprad,maxnks*maxnrs)
 
      endif
c
c     Read in HC related data 
c
c slmod begin           
c...  Seems that .raw modifications were made in more than one place for
c     +35, so I am advancing the subversion number to 37 in params and
c     eliminating the following load from my most recent results, 
c     which were for +36 (and this subversion should be skipped by everyone
c     else):
      if (version_code.ge.(6*maxrev+35).and.
     .    version_code.ne.(6*maxrev+36)) then 
c
c      if (version_code.ge.(6*maxrev+35).AND..FALSE.) then 
c slmod end
c
! ammod begin
c
         read(8) global_hc_follow_option 
c
         if (global_hc_follow_option.ne.0) then
c
            call global_hc_read_raw_data 
c
         endif
c
         Call rinout ('R BRATIO',BRATIO,maxnks*maxnrs)
c
! ammod end
c
      endif
c
c     Add code to read subgrid data if it was created from the DIVIMP run
c
      if (version_code.ge.(6*maxrev+40)) then 

         call reload_subgrid(8)

      endif

c
c     Temporarily Add the following
c
      call rinout ('R FLUXES',fluxes,maxnks*maxnrs*16)
c
c
c
      IF (IMODE.EQ.1) THEN
      CALL RINOUT ('R LIMS  ',LIMS  ,MAXNKS*MAXNRS*(MAXIZS+2)*MAXNTS)
      ENDIF

c
c slmod begin - new
c
c
c I have modified the version checking code so that it is compatible
c with both of our current .raw file formats.  I am abandoning the slver
c method from now on.  However, we have to be sure to notify each other
c when we change the DIVIMP revision number and then do an immediate
c synchronization of iodiv/ioout, so that we will still be able to read each
c others input files. - Dec 4, 1999
c
c
c...  Just in case on old .raw file has this version_code, but
c     did not yet include slver -- may be the case for some of
c     my old files:
      READ(8,ERR=9003,END=90) slver

      sldata = slver

      IF ( version_code.GE.(6*maxrev+4).OR.
     .    (version_code.GE.(6*maxrev+3).AND.slver.GE.1.0)) THEN
        READ(8) idum1,idum2,
     .          eirnpgdat,((eirpgdat(i1,i2),i2=1,idum1),i1=1,idum2)
      ENDIF

      IF ( version_code.GE.(6*maxrev+5).OR.
     .    (version_code.GE.(6*maxrev+4).AND.slver.GE.1.1)) THEN

        IF (slver.GE.2.6) THEN
          READ(8) asc_ncell,MAXASC3
          IF (maxasc3.GT.MAXASC) CALL ER('Get','MAXASC too small',*9004)
        ELSE
c...      NOTE: THIS IS *TEMPORARY* AS THE DATA LOADED INTO MULTI-DIMENSIONAL
c               ARRAYS WILL BE BOGUS
          READ(8) asc_ncell
          maxasc3 = 1000
        ENDIF

        IF (slver.GE.3.2) THEN

          CALL IINOUT('R CELL  ',asc_cell  ,MAXASC)
          CALL IINOUT('R REGION',asc_region,MAXASC)
          CALL IINOUT('R LINK  ',asc_link  ,MAXASC*4)
          CALL IINOUT('R GRID  ',asc_grid  ,MAXASC*2)
          CALL IINOUT('R NVP   ',asc_nvp   ,MAXASC)
          CALL RINOUT('R RVP   ',asc_rvp   ,MAXASC*8)
          CALL RINOUT('R ZVP   ',asc_zvp   ,MAXASC*8)
          IF (version_code.LE.(6*maxrev+13)) THEN
c...        MAXASCDAT was 20000 before moving to 6.14:
            CALL RINOUT('R VOL   ',asc_vol   ,20000)
          ENDIF
          CALL IINOUT('R NVERTX',ascnvertex,MAXASC)
          IF (version_code.GE.(6*maxrev+21)) THEN
            CALL RINOUT('R VERTEX',ascvertex ,MAXASC*40)       
          ELSE
            CALL RINOUT('R VERTEX',ascvertex2,MAXASC*40)       
            DO i1 = 1, asc_ncell
              DO i2 = 1, 40
                ascvertex(i2,i1) = ascvertex2(i2,i1)
              ENDDO
            ENDDO
          ENDIF

        ELSE

          WRITE(0,*) 'GET: READING PINASC DATA TO DUMMY ARRAY.'
          WRITE(0,*) '     2D VACUUM GRID PLOTS WILL NOT WORK.'

          CALL IINOUT('R CELL  ',asc_cell2  ,maxasc3)
          CALL IINOUT('R REGION',asc_region2,maxasc3)
          CALL IINOUT('R LINK  ',asc_link2  ,maxasc3*4)
          CALL IINOUT('R GRID  ',asc_grid2  ,maxasc3*2)
          CALL IINOUT('R NVP   ',asc_nvp2   ,maxasc3)
          CALL RINOUT('R RVP   ',asc_rvp2   ,maxasc3*8)
          CALL RINOUT('R ZVP   ',asc_zvp2   ,maxasc3*8)
          CALL RINOUT('R VOL   ',asc_vol2   ,maxasc3)
          IF     (version_code.GE.(6*maxrev+10).AND.slver.GE.2.2) THEN
            CALL IINOUT('R NVERTX',ascnvertex2,maxasc3)
            CALL RINOUT('R VERTEX',ascvertex2 ,maxasc3*40)       
          ELSEIF (version_code.GE.(6*maxrev+10).AND.slver.GE.2.1) THEN
            CALL IINOUT('R NVERTX',ascnvertex2,maxasc3)
            CALL RINOUT('R VERTEX',ascvertex2 ,maxasc3*20)       
          ELSEIF (version_code.GE.(6*maxrev+10).AND.slver.GE.2.0) THEN
c...        Messed this up.  Fixed in with slver=2.1:
            READ(8) idum1
            CALL RINOUT('R VERTEX',ascvertex2 ,idum1*2) 
          ENDIF

        ENDIF

      ELSE
        CALL ER('Get','Invalid MAXASC3 value (old RAW file)',*9004)
      ENDIF

      IF     ( version_code.GE.(6*maxrev+8).AND.slver.LE.1.8) THEN
        READ(8) eirnres
        CALL RINOUT('R EIRRES',eirres,20*MAXPINITER)
      ELSEIF ((version_code.GE.(6*maxrev+5).OR.
     .        (version_code.GE.(6*maxrev+4).AND.slver.GE.1.2)).AND.
     .                                          slver.LE.1.8) THEN
c...    Temporarily avoiding MAXPINITER:
        WRITE(0,*) 'TEMPORARILY AVOIDING MAXPINITER'
        READ(8) eirnres
        CALL RINOUT('R EIRRES',eirres,20*50)
      ENDIF

      IF     (version_code.GE.(6*maxrev+34)) THEN
        CALL RINOUT('R PINLN1',pinline(1,1,1,H_BALPHA),MAXNKS*MAXNRS*6)
        CALL RINOUT('R PINLN2',pinline(1,1,1,H_BBETA ),MAXNKS*MAXNRS*6)
        CALL RINOUT('R PINLN3',pinline(1,1,1,H_BGAMMA),MAXNKS*MAXNRS*6)
      ELSEIF (version_code.GE.(6*maxrev+5).OR.
     .       (version_code.GE.(6*maxrev+4).AND.slver.GE.1.3)) THEN
        CALL RINOUT('R PINLN1',pinline(1,1,1,H_BALPHA),MAXNKS*MAXNRS*6)
        CALL RINOUT('R PINLN2',pinline(1,1,1,H_BGAMMA),MAXNKS*MAXNRS*6)
      ENDIF

      IF (version_code.GE.(6*maxrev+5).AND.slver.GE.1.4) THEN
        READ(8) pincode
        CALL RINOUT ('R PINMOI',pinmoi,MAXNKS*MAXNRS)
      ENDIF

      IF (version_code.GE.(6*maxrev+5).AND.slver.GE.1.5) THEN
        CALL RINOUT('R OSMCDE',osmcde,MAXNKS*MAXNRS)
        CALL RINOUT('R OSMCDI',osmcdi,MAXNKS*MAXNRS)      
        CALL RINOUT('R OSMCVE',osmcve,MAXNKS*MAXNRS)
        CALL RINOUT('R OSMCVI',osmcvi,MAXNKS*MAXNRS)      
      ENDIF

      IF (version_code.GE.(6*maxrev+5).AND.slver.GE.1.6) THEN
        READ(8) tarshift(IKLO),tarshift(IKHI)
        READ(8) te_mult_o,ti_mult_o,n_mult_o,
     .          te_mult_i,ti_mult_i,n_mult_i
      ENDIF

      IF (slver.GE.1.7) THEN
        READ(8) lpdatsw
        CALL RINOUT('R SEPDIS',sepdist ,MAXNDS_)
        CALL RINOUT('R SEPDI2',sepdist2,MAXNDS_)
      ENDIF

      IF (slver.GE.1.8) THEN
        CALL IINOUT('R PRBNUM',prb_num ,NUMPRB) 
        CALL RINOUT('R PRBRHO',prb_rho ,MAXPRB*NUMPRB)
        CALL RINOUT('R PRBTE ',prb_te  ,MAXPRB*NUMPRB)
        CALL RINOUT('R PRBTI ',prb_ti  ,MAXPRB*NUMPRB)
        CALL RINOUT('R PRBNE ',prb_ne  ,MAXPRB*NUMPRB)
        CALL RINOUT('R PRBR  ',prb_r   ,MAXPRB*NUMPRB)
        CALL RINOUT('R PRBZ  ',prb_z   ,MAXPRB*NUMPRB)
      ENDIF

      IF (slver.GE.1.9) THEN
        READ(8) eirnres
        CALL RINOUT('R EIRRES',eirres,6*7*MAXPINITER)
      ENDIF

      IF     (slver.GE.2.6.AND.version_code.LE.(6*maxrev+13)) THEN
        IF (slver.GE.3.2) THEN
          WRITE(0,*) 'GET: LOADING OLD PINASD DATA AND CONVERTING'
          CALL RINOUT('R PINASD',pinasd3,20000*MAXASD2*MAXASS*2)
c...      Copy data from PINASD3 to PINASD:
          DO i1 = 1, 20000
            DO i2 = 1, MAXASD2
              DO i3 = 1, MAXASS
                pinasd(i1,i2,i3,1) = pinasd3(i1,i2,i3,1)
                pinasd(i1,i2,i3,2) = pinasd3(i1,i2,i3,2)
              ENDDO
            ENDDO
          ENDDO
c          CALL RINOUT('R PINASD',pinasd,MAXASCDAT*MAXASD2*MAXASS*2)
        ELSE
          WRITE(6,*) 'GET: LOADING OLD PINASD DATA AND CONVERTING'
          CALL RINOUT('R PINASD',pinasd2,2000*MAXASD2*MAXASS*3)
c...      Copy data from PINASD2 to PINASD:
          DO i1 = 1, 2000
            DO i2 = 1, MAXASD2
              DO i3 = 1, MAXASS
                pinasd(i1,i2,i3,1) = pinasd2(i1,i2,i3,1)
                pinasd(i1,i2,i3,2) = pinasd2(i1,i2,i3,3)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ELSEIF (slver.GE.2.3.AND.version_code.LE.(6*maxrev+13)) THEN
        CALL RINOUT('R PINASD',pinasd,300    *MAXASD2*MAXASS*3)
      ENDIF

      IF (slver.GE.2.4) THEN
        IF (version_code.GE.(6*maxrev+17)) THEN
        ELSE
          CALL RINOUT('R PINSTR',pinstrata2,MAXNKS*MAXNRS*3*6)
        ENDIF
        CALL RINOUT('R PINPLO',pinploss ,MAXNKS*MAXNRS*NMOMCHA)
      ENDIF

      IF (version_code.GE.(6*maxrev+19)) THEN
        READ(8) eirnpuff,eirpmode
        CALL RINOUT('R PUFF  ',eirpuff,MAXNAS*MAXPUFF)
      ELSEIF (slver.GE.2.5) THEN
        READ(8) eirnpuff,eirpmode
        IF (slver.GE.3.4) THEN
          CALL RINOUT('R PUFF  ',eirpuff3,MAXNAS*10 )
          DO i1 = 1, MAXNAS
            DO i2 = 1, 10
              eirpuff(i1,i2) = eirpuff3(i1,i2)
            ENDDO
          ENDDO
        ELSE
          CALL RINOUT('R PUFF  ',eirpuff2,20*MAXPUFF)
          DO i1 = 1, 20
            DO i2 = 1, MAXPUFF
              eirpuff(i1,i2) = eirpuff2(i1,i2)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      IF (slver.GE.2.7) THEN
        CALL RINOUT('R BGK   ',pinbgk   ,MAXNKS*MAXNRS*MAXBGK)
        READ(8) idum1,idum2,
     .          eirnspdat,((eirspdat(i1,i2),i2=1,idum1),i1=1,idum2)
      ENDIF

      IF (slver.GE.3.2) THEN

        IF (version_code.LE.(6*maxrev+13)) THEN
c...      MAXASCDAT was 2000 before 6.14:
          WRITE(0,*) 'GET: LOADING OLD ASCDATA AND CONVERTING'
          CALL RINOUT('R ASCDAT',ascdata3,20000*5)
          DO i1 = 1, 20000
            DO i2 = 1, 5
              ascdata(i1,i2) = ascdata3(i1,i2)
            ENDDO
          ENDDO
c          CALL RINOUT('R ASCDAT',ascdata   ,MAXASCDAT*5)
        ENDIF

        READ (8) asc_3dmode
        CALL IINOUT('R NVP3D ',asc_nvp3D ,MAXASC3D)
        CALL IINOUT('R LIMK3D',asc_link3D,MAXASC3D*6)
        CALL RINOUT('R ZMIN3D',asc_zmin3D,MAXASC3D)
        CALL RINOUT('R ZMAX3D',asc_zmax3D,MAXASC3D)
        CALL RINOUT('R XVP3D ',asc_xvp3D ,MAXASC3D*6*8)
        CALL RINOUT('R YVP3D ',asc_yvp3D ,MAXASC3D*6*8)
        CALL RINOUT('R ZVP3D ',asc_zvp3D ,MAXASC3D*6*8)

      ELSEIF (MAXASC3.EQ.5000) THEN

        IF (slver.GE.2.8) THEN
          CALL RINOUT('R ASCDAT',ascdata2,MAXASC3*5)
        ENDIF
        IF (slver.GE.3.0) THEN
          READ (8) asc_3dmode
          CALL IINOUT('R NVP3D ',asc_nvp3D2 ,MAXASC3)
          CALL IINOUT('R LIMK3D',asc_link3D2,MAXASC3*6)
          CALL RINOUT('R ZMIN3D',asc_zmin3D2,MAXASC3)
          CALL RINOUT('R ZMAX3D',asc_zmax3D2,MAXASC3)
          CALL RINOUT('R XVP3D ',asc_xvp3D2 ,MAXASC3*6*8)
          CALL RINOUT('R YVP3D ',asc_yvp3D2 ,MAXASC3*6*8)
          CALL RINOUT('R ZVP3D ',asc_zvp3D2 ,MAXASC3*6*8)
        ELSE
          asc_3dmode = 0
        ENDIF

      ELSE
        CALL ER('Get','Invalid MAXASC3 value',*9004)
      ENDIF

      IF (slver.GE.3.1) THEN
        READ (8) ascncut,idum1,idum2
      ENDIF

      IF (slver.GE.3.3) THEN
        READ(8) eirnsdtor,(eirsdtor(i1),eirsdvol(i1),i1=1,eirnsdtor)
        IF (eirnsdtor.GT.MAXTOR) THEN
          WRITE(0,*) 'ERROR: EIRNSDTOR.GT.MAXTOR, INCREASE MAXTOR'
          WRITE(0,*) '  EIRNSTDOR=',eirnsdtor
          WRITE(0,*) '  MAXTOR   =',MAXTOR
          WRITE(0,*) 'STOP'
        ENDIF
        DO i1 = 2, eirnsdtor
          i2 = (i1-1)*MAXBGK+1
          CALL RINOUT('R BGKTOR',pinbgk(1,1,i2),MAXNKS*MAXNRS*MAXBGK)
        ENDDO
      ENDIF

      IF (version_code.ge.(6*MAXREV+13)) THEN
        READ(8) eirniontime,idum1,idum2
        CALL RINOUT('R IONTIM',eiriontime,MAXIONTIME*(20+MAXBIN*3))
      ENDIF



      IF     (version_code.GE.(6*maxrev+20)) THEN
        READ (8) idum1,idum2,idum3,idum4
        IF (idum1.GT.MAXASD2) CALL ER('Get','MAXASD2   too small',*9004)
        IF (idum2.GT.MAXASS ) CALL ER('Get','MAXASS    too small',*9004)
        IF (idum3*idum4.GT.MAXASCDAT) 
     .                        CALL ER('Get','MAXASCDAT too small',*9004)
        READ (8) 
     .    (asc_vol(i1),
     .     (ascdata(i1,i2),i2=1,5),
     .     ((pinasd(i1,i2,i3,1),pinasd(i1,i2,i3,2),i2=1,idum1),
     .                                             i3=1,idum2 ), 
     .    i1=1,idum3*idum4+1+eirnpgdat),idum5
        IF (idum5.NE.999999) 
     .    CALL ER('Get','Version 6.14 data flag not found',*9004)
      ELSEIF (version_code.GE.(6*maxrev+14)) THEN
        READ (8) idum1,idum2,idum3,idum4
        IF (idum1.GT.MAXASD2) CALL ER('Get','MAXASD2   too small',*9004)
        IF (idum2.GT.MAXASS ) CALL ER('Get','MAXASS    too small',*9004)
        IF (idum3*idum4.GT.MAXASCDAT) 
     .                        CALL ER('Get','MAXASCDAT too small',*9004)
        READ (8) 
     .    (asc_vol(i1),
     .     (ascdata(i1,i2),i2=1,5),
     .     ((pinasd(i1,i2,i3,1),pinasd(i1,i2,i3,2),i2=1,idum1),
     .                                             i3=1,idum2 ), 
     .    i1=1,idum3*idum4),idum5
        IF (idum5.NE.999999) 
     .    CALL ER('Get','Version 6.14 data flag not found',*9004)
      ENDIF

      IF (version_code.GE.(6*maxrev+15)) THEN
        READ (8) eirntally,idum2,idum3,idum4
        IF (idum2.GT.MAXNKS)   CALL ER('Get','MAXNKS   too small',*9004)
        IF (idum3.GT.MAXNRS)   CALL ER('Get','MAXNRS   too small',*9004)
        IF (idum4.GT.MAXTALLY) CALL ER('Get','MAXTALLY too small',*9004) 
        READ (8) 
     .    ((eirtally(i1,i2),i2=1,4),
     .     ((pinalgv(ik,ir,i1),ik=1,nks(ir)),ir=1,nrs),
     .     i1=1,eirntally),idum5
        IF (idum5.NE.999999) 
     .    CALL ER('Get','Version 6.15 data flag not found',*9004)
      ENDIF
      IF (version_code.GE.(6*maxrev+16)) THEN
        READ(8) eirzaa,eirtorfrac,eirsrcmul
      ENDIF
      IF (version_code.GE.(6*maxrev+22)) THEN
        READ(8) osmns28,idum1,osm_nfnt,idum2,idum3
        READ(8) ((osms28(i1,i2),i2=1,idum1),i1=1,osmns28),
     .          ((osm_ionfnt(i1,i2),i2=1,idum2),i1=1,osm_nfnt)
      ENDIF
      IF (version_code.GE.(6*maxrev+25)) THEN
c...    Data for contributions to pressure for selected additional cells
c       for each stratum:
        READ(8) eirnstrdat,eirnstrai,eirnatmi,eirnmoli,
     .          idum1,idum2,idum3
        IF (idum1.GT.MAXSTRDAT) CALL ER('Get','MAXSTRDAT bad',*9004)
        IF (idum2.GT.MAXSTRATA) CALL ER('Get','MAXSTRATA bad',*9004)
        READ(8) (((eirstrdat(i1,i2,i3),i3=1,idum3),
     .                                 i2=1,idum2),
     .           eirstrlis(i1)        ,i1=1,idum1),idum4
        IF (idum4.NE.999999) CALL ER('Get','6.25 flag lost',*9004)
      ENDIF
      IF (version_code.GE.(6*maxrev+26)) THEN
        CALL IINOUT('R IKBOUN',ikbound,MAXNRS*2)
        CALL IINOUT('R IKFLUI',ikfluid,MAXNRS*2)
      ENDIF

      IF (version_code.GE.(6*maxrev+30)) THEN
c...    This is temporary:
        DO i1 = H_ION1, H_ION1+11
          CALL IINOUT('R PINDAT',pindata(1,1,i1),MAXNKS*MAXNRS)
        ENDDO
      ELSEIF (version_code.GE.(6*maxrev+27)) THEN
        CALL IINOUT('R PINDI1',pindata(1,1,H_ION1),MAXNKS*MAXNRS)
        CALL IINOUT('R PINDI2',pindata(1,1,H_ION2),MAXNKS*MAXNRS)
        CALL IINOUT('R PINDI3',pindata(1,1,H_ION3),MAXNKS*MAXNRS)
      ENDIF

      IF (version_code.GE.(6*maxrev+28)) THEN
        DO i1 = 1, 3
          CALL IINOUT('R PINICP',pinioncomp(1,1,i1),MAXNKS*MAXNRS)
        ENDDO
      ENDIF
      IF (version_code.GE.(6*maxrev+29)) THEN
        READ(8) eirntorseg
      ENDIF

      IF (version_code.GE.(6*maxrev+33)) THEN
        READ(8) ciopte,cxsc,cysc,cxsca,cysca,cxscb,cyscb
      ENDIF

c *TEMP*
      IF (slver.GE.3.5) THEN
        CALL RINOUT('R EIRPH1',eirpho1,MAXNKS*MAXNRS)
        CALL RINOUT('R EIRPH2',eirpho2,MAXNKS*MAXNRS)
      ENDIF

      IF (version_code.GE.(6*maxrev+41)) THEN
        READ (8) debugv,cstepv
        IF (debugv) CALL RINOUT ('R SDVS',sdvs,MAXNKS*MAXNRS*(MAXIZS+2))      
      ENDIF

      IF (version_code.GE.(6*maxrev+14)) THEN
        READ(8) idum1
        IF (idum1.NE.123456789)
     .    CALL ER('Get','End of data flag not found',*9004)
      ENDIF

90    CONTINUE

c...  This is temporary, and is the appropriate value
c     for all C-Mod plenum work prior to Oct 16, 2001:
      IF (eirzaa.EQ.-99.0) eirzaa = 0.47


c...  Generate file with title and description strings for posting to results web
c     page:
      

c *TEMP*
c...  Decide which machine is being modeled:
      IF     (rxp.GT.0.4.AND.rxp.LT.0.6) THEN
        machine = CMOD
      ELSEIF (rxp.GT.1.0.AND.rxp.LT.1.7) THEN
        machine = DIIID
      ELSE
        machine = NULL
      ENDIF

      nopriv = .FALSE.
      IF (ikto.LE.1) nopriv = .TRUE.
c slmod end
c
c
c------------------------------------------------------------------------------------
c
c
c     Temporary work-around - recalculate wallindex array values from wallpt data array since
c     wallindex is not in the .raw file
c
      call rzero(wallindex,maxnds)
c
      do in = 1,wallpts 
c
         if (wallpt(in,18).gt.0.and.wallpt(in,18).le.maxnds) then 
            wallindex(int(wallpt(in,18)))  = in
         endif
c
      end do
c
c------------------------------------------------------------------------------------
c


c
c
      RETURN
C
 9001 FORMAT(//1X,'GET:     NXS    NYS    NRS    NDS   NIZS    NTS',
     >  /6X,6I7,/1X,'(MAX)',6I7,/6X,A,/6X,A,
     >        /1X,'        ITER   MODE   REFCT  MAXSEG NVESM  NVESP',
     >  /6x,6I7,/)
 9002 FORMAT(//1X,'WARNING! PROGRAM HAS BEEN RECOMPILED SINCE ',
     >  'RESULTS FILE CREATED.',/10X,'RESULTS  WRITTEN  BY  DIV',A5,
     >                          /10X,'RESULTS BEING READ BY OUT',A5,/)
c slmod begin - new
 9003 CALL ER('Get','Invald .RAW file format',*9004)
 9004 STOP
c slmod end
      END
C
C  *********************************************************************
C  *  PRDMAT:  PRINT A PLASMA PROFILE                                  *
C  *********************************************************************
C
      SUBROUTINE PRDMAT(A,MAXNKS,NK,NR,IRSEP,IRWALL,JFP,JLP,IWT,TIT)
      implicit none
C
C INPUT
C -------
C BY ARGUMENT-LIST:
C    A        - MATRIX OF DIMENSION A(MAXNKS,NR)
C    MAXNKS   - LEADING DIMENSION OF A
C    NK       - NUMBER OF KNOTS IN SOL RINGS
C    NR       - NUMBER OF RINGS
C    IRSEP    - FIRST SOL RING
C    IRWALL   - LAST SOL RING
C    JFP      - FIRST KNOT IN MAIN SOL (AFTER OUTER DIVERTOR)
C    JLP      - LAST  KNOT IN MAIN SOL (BEFORE INNER DIVERTOR)
C    IWT      - OUTPUT-CHANNEL FOR PRINTOUT
C    TIT      - CHARACTER STRING FOR TITLE
C=======================================================================
      integer  maxnks,nr,nk,irsep,irwall,jfp,flp,iwt
      integer  ii,nrf,nrl,inum,npage,k,l1,l2,i,j,n,jlp
      REAL A(MAXNKS,NR)  
      CHARACTER*124 LINE
      CHARACTER*13 STAR,CORE,BLANK
      DATA STAR /'       U     '/
      DATA CORE /'.............'/
      DATA BLANK/'             '/
      CHARACTER*(*) TIT
C INUM = NUMBER OF COLUMNS ON PAGE
      DATA INUM/9/
C
      WRITE(IWT,50)TIT
C
      NPAGE = (NR-1)/INUM + 1
      DO 40 II=1,NPAGE
         N   = MIN0(NR-(II-1)*INUM,INUM)
         NRF = (II-1)*INUM + 1
         NRL = (II-1)*INUM + N
C
         WRITE(LINE,85) (J,J=NRF,NRL)
         WRITE(IWT,95) LINE
         WRITE(IWT,95)
C
         DO 30 I=1,NK
            WRITE(LINE,60) I
            DO J = NRF, NRL
              K = J - NRF + 1
              L1 = 8 + (K-1)*13
              L2 = L1 + 12
              IF (J.LT.IRSEP) THEN
                IF (I.LT.JFP) THEN
                  LINE(L1:L2) = STAR
                ELSE IF (I.GT.JLP) THEN
                  LINE(L1:L2) = STAR
                ELSE
                  WRITE(LINE(L1:L2),70) A(I-JFP+1,J)
                ENDIF
              ELSE IF (J.GT.IRWALL) THEN
                IF (I.LT.JFP) THEN
                  WRITE(LINE(L1:L2),70) A(I,J)
                ELSE IF (I.LE.JLP) THEN
                  LINE(L1:L2) = STAR
                ELSE
                  WRITE(LINE(L1:L2),70) A(I-JLP+JFP-1,J)
                ENDIF
              ELSE
                WRITE(LINE(L1:L2),70) A(I,J)
              ENDIF
            ENDDO
            IF( I.EQ.JFP .AND. NRF.LT.IRSEP )
     &          WRITE(IWT,80) (CORE,J=NRF,MIN0(NRL,IRSEP-1))
            WRITE(IWT,90) LINE
            IF( I.EQ.JLP .AND. NRF.LT.IRSEP )
     &          WRITE(IWT,80) (CORE,J=NRF,MIN0(NRL,IRSEP-1))
   30    CONTINUE
   40 CONTINUE
C
      RETURN
C
   50 FORMAT(////,1X,A/1X,132('-'))
   60 FORMAT(1X,I4,2X)
   70 FORMAT(1P,E13.5)
   80 FORMAT(7X,9(A))
   85 FORMAT(11X,9(I4,9X))
   90 FORMAT(A)
   95 FORMAT(/A)
      END
c
c
