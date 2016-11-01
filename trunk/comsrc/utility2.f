c -*-Fortran-*- 
c
c ======================================================================
c
c function: CalcWidth
c
c Return an estimate of the width of cell IK,IR.
c
      REAL FUNCTION CalcWidth(ik,ir,s,mode)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

c     Input:
      INTEGER ik,ir,mode
      REAL    s

      REAL    r,z,r1,z1,r2,z2,c1,c2,wid(3)
      REAL    z3,z4,r3,r4
      REAL    widtop,widbot,cbot1,ctop1,ctop2,cbot2
      REAL    deltar,deltaz
      REAL    rside1,zside1,rside2,zside2
      INTEGER id
c
c Hard coded width option (need to restore this when adding to
c Dave's version, since slcom is not included there):
c
c      INTEGER cwcopt

      cwcopt = 2

      id = korpg(ik,ir)
c
c
c
      IF (mode.NE.SIDE14.AND.mode.NE.side23.AND.mode.NE.TOTAL)
     .  STOP 'ERROR (CalcWidth): Invalid MODE value'
c
c
c
      IF (ir.EQ.1.OR.ir.EQ.irwall.OR.ir.EQ.irtrap) THEN
        WRITE(0,*) 'ERROR (CalcWidth): Virtual ring (IR = ',ir,')'
        STOP
      ENDIF
c
c
c
      IF (nvertp(id).NE.4) THEN
        write (0,*) 'ERROR (CalcWidth): Invalid cell',ik,ir,id,
     >   nvertp(id)
        STOP 'ERROR (CalcWidth): Invalid cell'
      ENDIF
c
c
c
      IF (s.EQ.CENTER) THEN
        r = rs(ik,ir)
        z = zs(ik,ir)
      ELSE
        IF (cwcopt.NE.2)
     .    CALL ER('CalcWidth','Cell center must be selected',*99)

c change ID if changing IK
c assign IK to a non-refenced variable if altering its value
        STOP 'ERROR (CalcWidth): Invalid S value'
      ENDIF
c
c
c
      IF (cwcopt.EQ.1) THEN
c
c       Approximate the cell width as the sum of the distances between the
c       cell side mid-points and the cell center:
c
        r1 = 0.5 * (rvertp(1,id) + rvertp(4,id))
        z1 = 0.5 * (zvertp(1,id) + zvertp(4,id))

        r2 = 0.5 * (rvertp(2,id) + rvertp(3,id))
        z2 = 0.5 * (zvertp(2,id) + zvertp(3,id))

        wid(SIDE14) = SQRT((r - r1)**2 + (z - z1)**2)
        wid(SIDE23) = SQRT((r - r2)**2 + (z - z2)**2)

c Debug:
c        WRITE(SLOUT,'(10X,A,2I3,F12.6)')
c     +    'WidDat:',IK,IR,CalcWidth

      ELSEIF (cwcopt.EQ.2) THEN
c
c       Perpendicular distance from point on center line to
c       each side:
c
c       Side 1-4:
        deltar = rvertp(4,id) - rvertp(1,id)
        deltaz = zvertp(4,id) - zvertp(1,id)

        c1 = ((r - rvertp(1,id)) * deltar +
     .        (z - zvertp(1,id)) * deltaz) / (deltar**2 + deltaz**2)

        r1 = rvertp(1,id) + c1 * deltar
        z1 = zvertp(1,id) + c1 * deltaz

c       Side 2-3:
        deltar = rvertp(3,id) - rvertp(2,id)
        deltaz = zvertp(3,id) - zvertp(2,id)

        c2 = ((r - rvertp(2,id)) * deltar +
     .        (z - zvertp(2,id)) * deltaz) / (deltar**2 + deltaz**2)

        r2 = rvertp(2,id) + c2 * deltar
        z2 = zvertp(2,id) + c2 * deltaz
c
c       Calculate distances:
c
        wid(SIDE14) = SQRT((r - r1)**2 + (z - z1)**2)
        wid(SIDE23) = SQRT((r - r2)**2 + (z - z2)**2)

c Debug:
c        WRITE(SLOUT,'(10X,A,2I3,4F10.4)')
c     +    'WidDat:',ik,ir,c1,c2,CalcWidth

      ELSEIF (cwcopt.EQ.3) THEN
        wid(SIDE14) = 0.5
        wid(SIDE23) = 0.5
      ELSE
        STOP 'ERROR (CalcWidth): Invalid option'
      ENDIF

      wid(TOTAL) = wid(SIDE14) + wid(SIDE23)

      CalcWidth = wid(mode)

      RETURN
99    STOP
      END
c
c
c
c
c
c ======================================================================
c
c function: OutsideBreak
c
      LOGICAL FUNCTION OutsideBreak(ir)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER ir

      OutsideBreak = .FALSE.
     
      IF ( ir.GE.irbreak.OR.
     .    (ir.LE.irwall.AND.irbreak.GE.irtrap)) OutsideBreak = .TRUE.

      RETURN
99    STOP
      END



      LOGICAL FUNCTION CHKINT(a1,a2,b1,b2,c1,c2,d1,d2)
      IMPLICIT none

      REAL*8 a1,a2,b1,b2,c1,c2,d1,d2

      REAL*8 TOL
      PARAMETER (TOL = 1.0D-07)

      REAL*8 t0,e1,e2,f1,f2,fact,tab,tcd

      CHKINT = .FALSE.

      IF (DABS(c1-d1).LT.TOL.AND.DABS(c2-d2).LT.TOL) THEN
        WRITE(0,*) 'CHKINT: ZERO LENGTH LINE SEGMENT' 
        WRITE(6,*) 'CHKINT: ZERO LENGTH LINE SEGMENT' 
        CALL EXIT
      ENDIF

c...  Calculate cross-product to see if the lines are parallel:
      IF (DABS((a1 - b1) * (c2 - d2)).LT.1.0D-8.AND.
     .    DABS((a2 - b2) * (c1 - d1)).LT.1.0D-8) RETURN

c...  Find projection of AB onto CD:
      t0 = ((a1 - c1) * (d1 - c1) + (a2 - c2) * (d2 - c2)) /
     .     ((c1 - d1)**2 + (c2 - d2)**2)
      e1 = c1 + t0 * (d1 - c1)
      e2 = c2 + t0 * (d2 - c2)

c...  Calculate F, the intersection point between AB and CD:
      fact = (e1 - a1) * (b1 - a1) + (e2 - a2) * (b2 - a2)

c...  Determine the parametric location of F on AB:
      IF (DABS(fact).LT.1.0D-10) THEN
        tab = 0.0
      ELSE
        tab = ((a1 - e1)**2 + (a2 - e2)**2) / fact
      ENDIF

c...  Determine the parametric location of F on CD:
      f1 = a1 + tab * (b1 - a1)
      f2 = a2 + tab * (b2 - a2)

      IF (DABS(d1-c1).GT.DABS(d2-c2)) THEN
        tcd = (f1 - c1) / (d1 - c1)
      ELSE
        tcd = (f2 - c2) / (d2 - c2)
      ENDIF

c...  Decide if the lines interect between their respective end points:
      IF (tab.GT.0.0D0+1.0D-8.AND.tab.LT.1.0D0-1.0D-8.AND.
     .    tcd.GT.0.0D0+1.0D-8.AND.tcd.LT.1.0D0-1.0D-8)
     .  CHKINT = .TRUE.

      RETURN
      END








      subroutine load_expt_data(dataunit,iseld,expt_axis,axis_type,
     >                    expt_data,maxcols,maxdatx,
     >                    num_expt,ncols,datatitle)
      use error_handling
c slmod begin - new
c
c Input:
c ISELD         - data index number
c MAXCOLS       - maximum number of data columns
c MAXDATX       - maximum number of data items
c
c Output:
c EXPT_AXIS     - independent data
c AXIS_TYPE     - read from file 13 (optional listing), default is 1
c               - 1 = Theta 
c               - 2 = R
c               - 3 = Dist
c               - 4 = R,Z
c               - 5 = Channel or count
c               - 6 = PSIN data   
c  
c
c EXPT_DATA     - MAXDATX,MAXCOLS dependent data
c NUM_EXPT      - number of data items in the data block, default is 0
c NCOLS         - number of data columns in data block, default is 1
c DATATITLE     - title, default is 'NO TITLE'
c
c 
c
c
c slmod end
      implicit none
c      include 'params' 
      integer dataunit,axis_type,num_expt,iseld,maxdatx,colindex
      integer maxcols,ncols
      character*(*) datatitle
      real expt_data(maxdatx,maxcols),expt_axis(maxdatx)
c
c     LOAD_EXPT_DATA: This routine loads the experimental data
c                     into the local arrays and passes
c                     back the relevant information.
c
c
c     Local variables
c
      character*200 buffer
      real xval,yval(maxcols),scalef
      real axis_shift,expt_norm  
      integer extstr,len,start,in,startn,endn,lenstr,colcnt
      external extstr,lenstr
      integer dataset_num,dataset_cnt,total_dataset_cnt
c
c     The following quantities were introduced to allow adjustments
c     to invlaid experimental values without actually changing the
c     contents of the file.
c
c     cutoff_val - this is for the elimination of spikes - any entries
c                  that are greater than this value will be set to zero.
c     offset_val - this is to graphically compensate for non-zero
c                  calibrartion offsets so plot comparison is easier.
c
      integer maxvals_cutoff
      parameter(maxvals_cutoff=25)
      real cutoff_val(maxvals_cutoff), offset_val,minexpt,maxexpt
      character*255 errormsg


      REAL    HI,LO
      PARAMETER( HI=1.E37 ,LO=1.E-37)
c
c
c     Initialize
c
      scalef = 1.0
      axis_shift = 0.0
      expt_norm = 0.0
      startn = 1
      endn = maxdatx
      offset_val = 0.0
      call rinit(cutoff_val,maxvals_cutoff,HI)
      maxexpt = -HI   
      minexpt = HI
      axis_type = 1
      num_expt = 0
      dataset_cnt = 0
      ncols = 1
      datatitle = 'NO TITLE'
c
      write (6,*) 'ISELD:',iseld
c
      rewind (dataunit)
c
c     Scan through data unit looking for INDEX keyword
c     with the appropriate data tag number.
c
 100  read(dataunit,'(a200)',end=500,err=500) buffer
c
c      len = lenstr(buffer)
c      write(6,*) 'BUFF:',buffer(1:len),':'
c
c     Ignore Empty Lines
c
c      if (buffer.eq.'') goto 100
c
c     Ignore comments
c
      if (buffer(1:1).eq.'$') goto 100
c
c     Check number of datasets in file
c
      if (buffer(1:10).eq.'FILECOUNT:') then
         read (buffer(11:),*) total_dataset_cnt
c
         write (6,*) 'DATASETS:',total_dataset_cnt
c
         if (iseld.gt.total_dataset_cnt) then
            write (6,*) 'REQUESTED EXPERIMENTAL DATA SET # ',
     >                  iseld,' DOES NOT EXIST'
            write (6,*) 'ONLY ', total_dataset_cnt,' DATA SETS ARE '//
     >                       'SPECIFIED BY FILECOUNT:'
            write (0,*) 'REQUESTED EXPERIMENTAL DATA SET # ',
     >                  iseld,' DOES NOT EXIST'
            write (0,*) 'ONLY ', total_dataset_cnt,' DATA SETS ARE '//
     >                       'SPECIFIED BY FILECOUNT:'
         endif

      endif
c
c     Look for INDEX and count for datasets
c
      if (buffer(1:6).eq.'INDEX:') then
         dataset_cnt = dataset_cnt+1
c
         read(buffer(7:),*) dataset_num
c
c         write (6,*) 'INDEX:',dataset_cnt,dataset_num,iseld
c
c        Indexing error - write out error message and continue
c
         if (dataset_cnt.ne.dataset_num) then
c
            write (errormsg,'(a,i6,a,i6)') 
     >               'ERROR IN EXPERIMENTAL DATA'//
     >               ' FILE: DATASET # ',dataset_num, ' IS IN'//
     >                  ' FILE POSITION ',DATASET_CNT
            call errmsg('LOAD_EXPT_DATA',errormsg)
C
         endif
c
         if (dataset_num.eq.iseld.or.dataset_cnt.eq.iseld) goto 300
c
      endif
c
c     Loop back and read more data until exit
c
      goto 100
c
c     Index or count for dataset found - continue processing
c
 300  continue
c
c     Read in the rest of the HEADER block until data found - then
c     start processing data until EOF or next INDEX.
c
 350  read(dataunit,'(a200)',end=400,err=400) buffer
c
c      len = lenstr(buffer)
c      write(6,*) 'DATA:',len,':',buffer(1:len),':'
c
c     Ignore Empty lines
c
c      if (buffer.eq.'') goto 350
c
c     Ignore comments
c
      if (buffer(1:1).eq.'$') goto 350
c
c     Exit if Next INDEX is found
c
      if (buffer(1:6).eq.'INDEX:') goto 400
c
c     Extract a Title if one is specified
c
      if (buffer(1:6).eq.'TITLE:') then
         len = extstr(buffer(7:),start)
         datatitle = buffer(8:LEN_TRIM(buffer))
c         datatitle = buffer(6+start:len)
         write(6,*) 'TITLE:',datatitle(1:LEN_TRIM(datatitle)),':'
      endif
c
c     Extract Scaling Factor if specified
c
      if (buffer(1:7).eq.'SCALEF:') read(buffer(8:),*) scalef
c
c     Extract Axis_type if given
c
      if (buffer(1:5).eq.'AXIS:') read(buffer(6:),*) axis_type
c
c     Extract Number of colums of data if given - one assumed
c
      if (buffer(1:6).eq.'NCOLS:') then
c
         read(buffer(7:),*) ncols
c
         write(6,*) 'NCOLS:',ncols
c
         if (ncols.gt.maxcols) then
c
c           Issue error message and only load the first column of
c           data
c
            write(errormsg,'(a,i6,a,i6)') 
     >                  'REQUESTED EXPERIMENTAL DATA SET # ',
     >                  iseld,' HAS MORE COLUMNS THAN MAX =',maxcols
            call errmsg('LOAD_EXPT_DATA',errormsg)
c
            ncols = 1
c
         endif
      endif
c
c     Extract Offset Value if specified
c
      if (buffer(1:7).eq.'OFFSET:') read(buffer(8:),*) offset_val
c
c     Extract Axis Shift if given
c
      if (buffer(1:6).eq.'SHIFT:') read(buffer(7:),*) axis_shift
c
c     Extract Normalization value if given
c
      if (buffer(1:5).eq.'NORM:') read(buffer(6:),*) expt_norm
c
c     Extract Cutoff Value if specified
c
      if (buffer(1:7).eq.'CUTOFF:') read(buffer(8:),*) 
     >                (cutoff_val(in),in=1,min(ncols,maxvals_cutoff))
c
c     Extract data limit counters if present
c
      if (buffer(1:6).eq.'COUNT:') then
         read(buffer(7:),*) startn,endn
         write (6,*) 'COUNT:',startn,endn
         if (endn-startn+1.gt.maxdatx) then 

            write(errormsg,'(a,i6,a,i6,a,i6,a)') 
     >           'NUMBER OF EXPERIMENTAL DATA POINTS IN DATASET ',
     >            iseld,' IS ',endn-startn+1,
     >           ' WHICH EXCEEDS THE AVAILABLE STORAGE OF ',maxdatx,
     >           ' - ONLY LOADING TO ARRAY LIMIT'

            call errmsg('LOAD_EXPT_DATA',errormsg)
            endn = maxdatx+startn-1
         endif

      endif
c
c     Other headers will be ignored - check for data and load it.
c
c     Data lines start with a blank
c
c     Need to add array checking so that the code does not try 
c     to load more data than can be held in the expt_axis and
c     expt_data arrays - checking on ncols is performed above.
c     
c
      if (buffer(1:1).eq.' ') then
c
         if(lenstr(buffer).gt.6) then
c
c            write (0,*) 'BUFFER:'//buffer//':',ncols
c
            read(buffer,*) in,xval,
     >                        (yval(colcnt),colcnt=1,ncols)
c
            if (in.ge.startn.and.in.le.endn) then
c
c              Assign a maximum index equal to the array size - issue 
c              an error message the first time this occurs.  
c
               num_expt = min(in-startn+1,maxdatx)
c
               if ((in-startn+1).eq.maxdatx+1) then 
                  call errmsg('LOAD_EXPT_DATA:',
     >               'EXPERIMENTAL DATA SET HAS TOO MANY'//
     >               ' ELEMENTS FOR STORAGE ARRAY - '//
     >               'NUMBER LIMITED TO MAXDATX')
               endif
c
               expt_axis(num_expt) = xval
c
c               if (axis_type.eq.5) then
c
c                  expt_axis(in-startn+1) = expt_axis(in-startn+1)
c     >                                     /360.0 * 2.0 * PI
c               endif
c
               do colcnt = 1,ncols
                  expt_data(num_expt,colcnt) = yval(colcnt)*scalef
               end do
c
c              Keep track of minimum value in first experimental data.
c
               minexpt = min(minexpt,expt_data(num_expt,1))
c
c              Keep track of maximum value in first experimental data.
c
               maxexpt = max(maxexpt,expt_data(num_expt,1))
c
c               write (6,'(A,2I6,5G12.4)') 'NUM:',
c     >           num_expt,in,xval,yval(1),yval(2),yval(3),scalef
c
            endif
c
         endif
c
      endif
c
      goto 350
c
c     Continue and wrap up file processing - experimental data
c     has been read.
c
 400  continue
c
c     Data cleanup processing ...
c
c     Check experimental data for cutoff and adjust by offset.
c
      if (offset_val.eq.-1.0) offset_val = minexpt
c
      write (6,*) 'OFFSET_VAL:',offset_val
c
      do in = 1,num_expt
c
         if (axis_shift.ne.0.0) expt_axis(in) = expt_axis(in) 
     >                          + axis_shift
c
         do colcnt = 1,ncols

            expt_data(in,colcnt) = expt_data(in,colcnt) - offset_val
c
            if (expt_norm.ne.0.0.and.maxexpt.ne.0.0) then 
               expt_data(in,colcnt) = 
     >                  expt_data(in,colcnt)/maxexpt * expt_norm
            endif
c
            if (expt_data(in,colcnt).gt.
     >          cutoff_val(min(colcnt,maxvals_cutoff)))
     >                                  expt_data(in,colcnt) = 0.0

         end do
c
      end do
c
      return
c
 500  continue
c
c     Error exit condition - DATA SET NOT FOUND
c
      write(errormsg,'(a,i6,a)') 
     >            'ERROR IN EXPERIMENTAL DATA: EOF OR DATASET # ',
     >             iseld,' NOT FOUND'
      call errmsg('LOAD_EXPT_DATA',errormsg)

c
c     Exit
c
      return
      end








c
c ======================================================================
c
c subroutine: MapArray
c
c
      SUBROUTINE MapArray(a1,b1,n11,n12,a2,b2,n21,n22)
      IMPLICIT none

      INCLUDE 'params'

      REAL    a1(*),a2(*),b1(*),b2(*)
      INTEGER n11,n12,n21,n22

      INTEGER i1,i2, select_opt

      DO i2 = n21, n22
        DO i1 = n11, n12
          IF (a1(i1).EQ.a2(i2)) THEN
            IF (b1(i1).NE.0.0.AND.b1(i1).NE.LO) then 
           WRITE(0,*) 'WARNING MapArray: Overwriting existing value'
           write(6,'(a,2i6,4(1x,g18.8))') 
     >        'Maparray overwrite:',i1,i2,a1(i1),a2(i2),b1(i1),b2(i2)
           write(0,'(a,2i6,4(1x,g18.8))') 
     >        'Maparray overwrite:',i1,i2,a1(i1),a2(i2),b1(i1),b2(i2)
            endif
            b1(i1) = b2(i2)
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END
c
c ======================================================================
c
c subroutine: LoadArray
c
c
      SUBROUTINE LoadArray(a1,n1,a2,n21,n22)

      IMPLICIT none

      REAL a1(*),a2(*)
      INTEGER n1,n21,n22

      INCLUDE 'params'

      INTEGER i1,i2,n31,n32
      REAL    curv,curm,a3(MAXGXS)

      n31 = 1
      n32 = n1
      
      DO i1 = 1, n1
        a3(i1) = a1(i1)
      ENDDO

      n1   =  0
      curv =  HI
      curm = -HI

c      WRITE(6,*)
c      DO i1 = n21,n22
c        WRITE(6,'(A,I4,1P,E15.7)')
c     .    '1st = ',i1,a2(i1)
c      ENDDO
c      WRITE(6,*)
c      DO i1 = n31,n32
c        WRITE(6,'(A,I4,1P,E15.7)')
c     .    '2nd = ',i1,a3(i1)
c      ENDDO

      DO i1 = 1, (n22 - n21) + (n32 - n31) + 2
        DO i2 = n21, n22
          IF (a2(i2).GT.curm) curv = MIN(curv,a2(i2))
        ENDDO
        DO i2 = n31, n32
          IF (a3(i2).GT.curm) curv = MIN(curv,a3(i2))
        ENDDO

        IF (curv.GT.curm.AND.curv.LT.HI) THEN
          n1     = n1 + 1
          a1(n1) = curv
          curm   = curv
          curv   = HI
        ENDIF
      ENDDO

      RETURN
      END


c
c ======================================================================
c
c subroutine LoadThomsonData
c
c input:
c
c output:
c   raw - 1 - R
c         2 - Z
c         3 - cell index
c         4 - ring index
c         5 - density
c         6 - temperature
c
      SUBROUTINE LoadThomsonData(nraw,raw,MAXTDAT,MAXCOLS,
     .                           xshift1,yshift1,mode)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER    DATAUNIT   ,MAXBLOCK
      PARAMETER (DATAUNIT=13,MAXBLOCK=100)

      REAL PsinToR

c...  Input:
      INTEGER mode,MAXTDAT,MAXCOLS
c      INTEGER mode,quantity,MAXTDAT,MAXCOLS
      REAL    xshift1,yshift1
c...  Output:
      INTEGER nraw
      REAL    raw(MAXTDAT,MAXCOLS)

      INTEGER   etype,ndata,ncol,i1,i2,ik,ir,i,nblock,
     .          lblock(MAXBLOCK),back,lcount
      LOGICAL   outofgrid,mode2,moredata
      REAL      xdata(MAXTDAT),edata(MAXTDAT,MAXCOLS),
     .          count(MAXNKS,MAXNRS),f,val,xshift2,yshift2
      CHARACTER datatitle*128,dataline*128,cdum1*128,cdum2*128

      DATA lcount /0/
      SAVE

      CALL RZero(count,MAXNKS*MAXNRS)

      moredata = .FALSE. 
      mode2    = .FALSE.
      back     = 0
      nraw     = 0



10    CONTINUE

      back = back + 1

      IF (mode.EQ.1) THEN

        READ(5,'(A128)',END=97,ERR=98) dataline

        IF   (dataline(8:11).EQ.'THOM'.OR.dataline(8:11).EQ.'thom'.OR.
     .        dataline(8:11).EQ.'Thom') THEN
c...      Found listing for location of Thomson data in unit 13 (experimental
c         data file):
          READ(dataline,*) cdum1,cdum2,xshift2,yshift2,nblock
          IF (nblock.LT.0.OR.nblock.GT.MAXBLOCK) THEN
            CALL ER('LoadThomsonData','Invalid index size for Thomson'//
     .                                ' data block list',*99)
          ELSE
            READ(dataline,*) cdum1,cdum2,xshift2,yshift2,nblock,
     .                       (lblock(i1),i1=1,nblock)

c...        MAJOR HACK JOB - can load the next line before leaving this
c           routine so that additional Thomson data can be loaded but a 
c           different shift applied (need to mention this on the plot
c           somewhere):
            IF (lblock(nblock).EQ.-1) THEN
              moredata = .TRUE.
              nblock   = nblock -1
            ENDIF

c...        This is rather selective, obtuse, but what can one do:
            IF (xshift2.EQ.98.0) xshift2 = tarshift(IKLO)
            IF (xshift2.EQ.99.0) xshift2 = tarshift(IKHI)
          ENDIF
        ELSE
c...      Required Thomson data index information not found:
          GOTO 97
        ENDIF
      ELSE
c...    Temp: Need to specify index list elsewhere when using DIVIMP, since
c       can't just read the next input line as in OUT:

        WRITE(0,*) 'LIMITED THOMSON DATA LOADING METHOD - FIX'

        xshift2 = 0.0
        yshift2 = 1.6

        nblock = 1
        lblock(1) = 1
c        nblock  = 17
c        lblock(1) = 1
c        DO i = 5, 20
c           lblock(i-4+1) = i
c        ENDDO
      ENDIF

      DO i = 1, nblock
        CALL Load_Expt_Data(DATAUNIT,lblock(i),xdata,etype,edata,
     .                      MAXCOLS,MAXTDAT,ndata,ncol,datatitle)

        WRITE(6,*) '============================================='
        WRITE(6,*) 'TITLE  = ',datatitle(1:LEN_TRIM(datatitle))
        WRITE(6,*) 'INDEX  = ',lblock(i)
        WRITE(6,*) 'TYPE   = ',etype
        WRITE(6,*) 'NDATA  = ',ndata
        WRITE(6,*) 'NCOL   = ',ncol

        IF (lblock(i).EQ.2) lcount = lcount + 1   ! * TEMP *

        DO i1 = 1, ndata
c
c         jdemod - what sort of hack is this? - should check the axis type of the 
c                  data in the experimental data file before hacking good R values
c                  to something else ... if (.true.) just doesn't cut it ...

c...      Need to recalculate R and rearrange as necessary for alternate DTS 
c         format (from Eric, EPS05):
c
c
c         jdemod - ONLY execute this code if the experimental data loaded is not already flagged as R,Z data
c                - IDEALLY - this code should only run when the input data is actually  flagged as PSIN,Z!   
c
          if (etype.ne.4) then
c          IF (.TRUE.) THEN
c            WIRTE(0,*) 'CONVERTING PSIN TO R'

c            IF (lblock(i).EQ.2.AND.MOD(REAL(lcount),2.0).EQ.0) THEN 
c              edata(i1,1) = (edata(i1,1) - (-1.36600006)) *  
c     .                      (0.226 / 0.281) + (-1.366)
c            ENDIF

            xdata(i1) = PsinToR(xdata(i1),edata(i1,1))
c            xdata(i1) = 1.45
          ENDIF

          xdata(i1)   = xdata(i1)   + xshift1 + xshift2
          edata(i1,1) = edata(i1,1) + yshift1 + yshift2

          ik = 0
          ir = 0
          CALL GridPos(ik,ir,xdata(i1),edata(i1,1),.FALSE.,outofgrid)

          IF (.NOT.outofgrid) THEN

            write(6,'(a,3i6,l4,3(1x,g12.5))') 'LoadThom:',ik,ir,i1,
     >             outofgrid,xdata(i1),edata(i1,1)
             
            IF (edata(i1,2).NE.0.0.AND.edata(i1,3).NE.0.0) THEN
              nraw = nraw + 1
              raw(nraw,1) = xdata(i1)    ! R
              raw(nraw,2) = edata(i1,1)  ! Z
              raw(nraw,3) = REAL(ik)     ! ik index
              raw(nraw,4) = REAL(ir)     ! ir index
              raw(nraw,5) = edata(i1,2)  ! plasma density
              raw(nraw,6) = edata(i1,3)  ! plasma temperature

c...HARDCODED, HACK
              IF (raw(nraw,5).LT.1.E+10) raw(nraw,5) =raw(nraw,5)*1.E+20

c              count(ik,ir) = count(ik,ir) + 1.0
c              f = 1.0 / count(ik,ir)
c              array(ik,ir) =  f * val + (1.0 - f) * array(ik,ir)
            ENDIF
          ENDIF

c          IF (ir.EQ.15) 
c     .     WRITE(6,90) ' TD : ',i1,xdata(i1),edata(i1,1),ik,ir,
c     .                outofgrid,rs(ik,ir),zs(ik,ir),
c     .                xdata(i1),(edata(i1,i2),i2=1,ncol)

c
c     jdemod - removed f from write statement since some warnings
c              were issued because it is not assigned a value before
c              use. 
c     .                xdata(i1),(edata(i1,i2),i2=1,ncol),f
c
90        FORMAT(A,I4,2F8.4,2I4,L2,2F8.4,1P,4E10.2,0P,F6.3)
c          WRITE(6,90) ' TD : ',i1,xdata(i1),edata(i1,1),ik,ir,
c     .                outofgrid,rs(ik,ir),zs(ik,ir),
c     .                xdata(i1),(edata(i1,i2),i2=1,ncol),f,array(ik,ir)
c90        FORMAT(A,I4,2F8.4,2I4,L2,2F8.4,1P,4E10.2,0P,F6.3,1P,E10.2,
c     .           0P)
        ENDDO
      ENDDO
c
c...  Check to see if there is another data line for Thomson data (ugly).  This
c     is done so that upstream Thomson data can be included, but with
c     different spatial shifts:
      IF (mode.EQ.1) THEN
        IF (moredata) GOTO 10
      ENDIF


c...  Return the total shift in XSHIFT1 for plot notes:
      xshift1 = xshift1 + xshift2

      RETURN
97    WRITE(6,*) 'THOMSON DATA BLOCK LIST NOT FOUND'
      WRITE(0,*) 'THOMSON DATA BLOCK LIST NOT FOUND'
      BACKSPACE 5
      RETURN
98    CALL ER('LoadThomsonData','OUT input file access error',*99)
99    STOP
      END


















c
c ======================================================================
c
c subroutine: ReadGeometry
c
      SUBROUTINE ReadGeometry(fp,error)
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER fp,error

      REAL    slver
      INTEGER ik,ir,id,ii,i1

      INTEGER nitersol1,rel_opt1 ,rel_nstep1,rel_niter1,
     .        rel_step1,rel_iter1,rel_count1,adp_opt1,
     .        hold_nvesm,hold_nvesp,hold_jvesm(MAXSEG)
      REAL    rel_frac1,hold_rvesm(MAXSEG,2),hold_zvesm(MAXSEG,2)

c...  Store neutral wall:
      hold_nvesm = nvesm
      hold_nvesp = nvesp
      DO i1 = 1, nvesm+nvesp
        hold_jvesm(i1)   = jvesm(i1)
        hold_rvesm(i1,1) = rvesm(i1,1)
        hold_rvesm(i1,2) = rvesm(i1,2)
        hold_zvesm(i1,1) = zvesm(i1,1)
        hold_zvesm(i1,2) = zvesm(i1,2)
      ENDDO


c.... Open file:
c      fp = 98
c      OPEN(UNIT=fp,FILE='geomty.dat',FORM='UNFORMATTED',
c     .     STATUS='OLD',ERR=98)

      READ(fp,ERR=98,END=99) slver

      READ(fp,ERR=98,END=99)
     .  nitersol1,rel_opt1 ,rel_frac1 ,rel_nstep1,rel_niter1,
     .  rel_step1,rel_iter1,rel_count1,adp_opt1
c
c     In attempt to save disk space the geometry data was only stored
c     when the grid was modified (during grid adaptation).  As a result,
c     some checking is required to verify that the geometry data
c     correctly corresponds to the plasma and PIN data:
c
      IF (rel_step1 .NE.rel_step.OR.rel_iter1.NE.rel_iter.OR.
     .    rel_count1.NE.rel_count) THEN
        WRITE(0,'(A)')
     .    'WARNING ReadGeometry: Sorry mate, wrong geometry data'
        WRITE(0,'(5X,A,3I4)')
     .    'STEP  ITER  COUNT  = ',rel_step ,rel_iter ,rel_count
        WRITE(0,'(5X,A,3I4)')
     .    'STEP1 ITER1 COUNT1 = ',rel_step1,rel_iter1,rel_count1
      ENDIF

      IF (slver.GE.1.02) THEN
        READ(fp,ERR=98,END=99)
     .    nrs,nds,ndsin,irwall,irtrap,irbreak,nbr,ikto,ikti,npolyp,
     .    vpolyp,vpolmin,r0,z0,rxp,zxp,dthetg
      ELSE
        READ(fp,ERR=98,END=99)
     .    nrs,nds,ndsin,irwall,irtrap,irbreak,nbr
      ENDIF

      IF     (slver.GE.1.08) THEN

        nvesm = 0
        nvesp = 0
        CALL RZero(rvesm,MAXSEG*2)
        CALL RZero(zvesm,MAXSEG*2)

        READ(fp,ERR=98,END=99) (
     .      nks(ir),ikmids(ir),ksmaxs(ir),idds(ir,1),idds(ir,2),
     .      idring(ir),ikto2(ir),ikti2(ir),
     .      rho(ir,IN14),rho(ir,CELL1),rho(ir,OUT23),
     .    ir=1,nrs),(
     .      thetat(id),ikds(id),irds(id),sepdist(id),sepdist2(id),
     .    id=1,nds),((
     .      rs    (ik,ir),zs    (ik,ir),kss   (ik,ir),kps  (ik,ir),
     .      kbacds(ik,ir),kfords(ik,ir),thetag(ik,ir),korpg(ik,ir),
     .      bratio(ik,ir),kbfs  (ik,ir),
     .      nvertp(MAX(1,korpg(ik,ir))),
     .      (rvertp(ii,MAX(1,korpg(ik,ir))),
     .       zvertp(ii,MAX(1,korpg(ik,ir))),
     .         ii=1,nvertp(MAX(1,korpg(ik,ir)))),
     .    ik=1,nks(ir)),ir=1,nrs),((
     .      ksb(ik,ir),kpb(ik,ir),
     .    ik=0,nks(ir)),ir=1,nrs),
     .    nvesm,nvesp,(rvesm(i1,1),rvesm(i1,2),zvesm(i1,1),zvesm(i1,2),
     .                 jvesm(i1),i1=1,nvesm+nvesp)       

      ELSEIF (slver.GE.1.07) THEN

        nvesm = 0
        nvesp = 0
        CALL RZero(rvesm,MAXSEG*2)
        CALL RZero(zvesm,MAXSEG*2)

        READ(fp,ERR=98,END=99) (
     .      nks(ir),ikmids(ir),ksmaxs(ir),idds(ir,1),idds(ir,2),
     .      idring(ir),ikto2(ir),ikti2(ir),
     .    ir=1,nrs),(
     .      thetat(id),ikds(id),irds(id),sepdist(id),sepdist2(id),
     .    id=1,nds),((
     .      rs    (ik,ir),zs    (ik,ir),kss   (ik,ir),kps  (ik,ir),
     .      kbacds(ik,ir),kfords(ik,ir),thetag(ik,ir),korpg(ik,ir),
     .      bratio(ik,ir),kbfs  (ik,ir),
     .      nvertp(MAX(1,korpg(ik,ir))),
     .      (rvertp(ii,MAX(1,korpg(ik,ir))),
     .       zvertp(ii,MAX(1,korpg(ik,ir))),
     .         ii=1,nvertp(MAX(1,korpg(ik,ir)))),
     .    ik=1,nks(ir)),ir=1,nrs),((
     .      ksb(ik,ir),kpb(ik,ir),
     .    ik=0,nks(ir)),ir=1,nrs),
     .    nvesm,nvesp,(rvesm(i1,1),rvesm(i1,2),zvesm(i1,1),zvesm(i1,2),
     .                 jvesm(i1),i1=1,nvesm+nvesp)       

      ELSEIF (slver.GE.1.06) THEN

        nvesm = 0
        nvesp = 0
        CALL RZero(rvesm,MAXSEG*2)
        CALL RZero(zvesm,MAXSEG*2)

        READ(fp,ERR=98,END=99) (
     .      nks(ir),ikmids(ir),ksmaxs(ir),idds(ir,1),idds(ir,2),
     .      idring(ir),
     .    ir=1,nrs),(
     .      thetat(id),ikds(id),irds(id),sepdist(id),sepdist2(id),
     .    id=1,nds),((
     .      rs    (ik,ir),zs    (ik,ir),kss   (ik,ir),kps  (ik,ir),
     .      kbacds(ik,ir),kfords(ik,ir),thetag(ik,ir),korpg(ik,ir),
     .      bratio(ik,ir),kbfs  (ik,ir),
     .      nvertp(MAX(1,korpg(ik,ir))),
     .      (rvertp(ii,MAX(1,korpg(ik,ir))),
     .       zvertp(ii,MAX(1,korpg(ik,ir))),
     .         ii=1,nvertp(MAX(1,korpg(ik,ir)))),
     .    ik=1,nks(ir)),ir=1,nrs),((
     .      ksb(ik,ir),kpb(ik,ir),
     .    ik=0,nks(ir)),ir=1,nrs),
     .    nvesm,nvesp,(rvesm(i1,1),rvesm(i1,2),zvesm(i1,1),zvesm(i1,2),
     .                 jvesm(i1),i1=1,nvesm+nvesp)

      ELSEIF (slver.GE.1.05) THEN

        nvesm = 0
        nvesp = 0
        CALL RZero(rvesm,MAXSEG*2)
        CALL RZero(zvesm,MAXSEG*2)

        READ(fp,ERR=98,END=99) (
     .      nks(ir),ikmids(ir),ksmaxs(ir),idds(ir,1),idds(ir,2),
     .    ir=1,nrs),(
     .      thetat(id),ikds(id),irds(id),sepdist(id),sepdist2(id),
     .    id=1,nds),((
     .      rs    (ik,ir),zs    (ik,ir),kss   (ik,ir),kps  (ik,ir),
     .      kbacds(ik,ir),kfords(ik,ir),thetag(ik,ir),korpg(ik,ir),
     .      bratio(ik,ir),kbfs  (ik,ir),
     .      nvertp(MAX(1,korpg(ik,ir))),
     .      (rvertp(ii,MAX(1,korpg(ik,ir))),
     .       zvertp(ii,MAX(1,korpg(ik,ir))),
     .         ii=1,nvertp(MAX(1,korpg(ik,ir)))),
     .    ik=1,nks(ir)),ir=1,nrs),((
     .      ksb(ik,ir),kpb(ik,ir),
     .    ik=0,nks(ir)),ir=1,nrs),
     .    nvesm,nvesp,(rvesm(i1,1),rvesm(i1,2),zvesm(i1,1),zvesm(i1,2),
     .                 jvesm(i1),i1=1,nvesm+nvesp)

      ELSEIF (slver.GE.1.04) THEN
        READ(fp,ERR=98,END=99) (
     .      nks(ir),ikmids(ir),ksmaxs(ir),idds(ir,1),idds(ir,2),
     .    ir=1,nrs),(
     .      thetat(id),ikds(id),irds(id),sepdist(id),sepdist2(id),
     .    id=1,nds),((
     .      rs    (ik,ir),zs    (ik,ir),kss   (ik,ir),kps  (ik,ir),
     .      kbacds(ik,ir),kfords(ik,ir),thetag(ik,ir),korpg(ik,ir),
     .      bratio(ik,ir),kbfs  (ik,ir),
     .      nvertp(MAX(1,korpg(ik,ir))),
     .      (rvertp(ii,MAX(1,korpg(ik,ir))),
     .       zvertp(ii,MAX(1,korpg(ik,ir))),
     .         ii=1,nvertp(MAX(1,korpg(ik,ir)))),
     .    ik=1,nks(ir)),ir=1,nrs),((
     .      ksb(ik,ir),kpb(ik,ir),
     .    ik=0,nks(ir)),ir=1,nrs),
     .    nvesm,(rvesm(i1,1),rvesm(i1,2),zvesm(i1,1),zvesm(i1,2),
     .           jvesm(i1),i1=1,nvesm)
      ELSEIF (slver.GE.1.03) THEN
        READ(fp,ERR=98,END=99) (
     .      nks(ir),ikmids(ir),ksmaxs(ir),idds(ir,1),idds(ir,2),
     .    ir=1,nrs),(
     .      thetat(id),ikds(id),irds(id),
     .    id=1,nds),((
     .      rs    (ik,ir),zs    (ik,ir),kss   (ik,ir),kps  (ik,ir),
     .      kbacds(ik,ir),kfords(ik,ir),thetag(ik,ir),korpg(ik,ir),
     .      bratio(ik,ir),kbfs  (ik,ir),
     .      nvertp(MAX(1,korpg(ik,ir))),
     .      (rvertp(ii,MAX(1,korpg(ik,ir))),
     .       zvertp(ii,MAX(1,korpg(ik,ir))),
     .         ii=1,nvertp(MAX(1,korpg(ik,ir)))),
     .    ik=1,nks(ir)),ir=1,nrs),((
     .      ksb(ik,ir),kpb(ik,ir),
     .    ik=0,nks(ir)),ir=1,nrs),
     .    nvesm,(rvesm(i1,1),rvesm(i1,2),zvesm(i1,1),zvesm(i1,2),
     .           jvesm(i1),i1=1,nvesm)

      ELSEIF (slver.GE.1.02) THEN
        READ(fp,ERR=98,END=99) (
     .      nks(ir),ikmids(ir),ksmaxs(ir),idds(ir,1),idds(ir,2),
     .    ir=1,nrs),(
     .      thetat(id),ikds(id),irds(id),
     .    id=1,nds),((
     .      rs    (ik,ir),zs    (ik,ir),kss   (ik,ir),kps  (ik,ir),
     .      kbacds(ik,ir),kfords(ik,ir),thetag(ik,ir),korpg(ik,ir),
     .      bratio(ik,ir),kbfs  (ik,ir),
     .      nvertp(korpg(ik,ir)),
     .      (rvertp(ii,korpg(ik,ir)),
     .       zvertp(ii,korpg(ik,ir)),ii=1,nvertp(korpg(ik,ir))),
     .    ik=1,nks(ir)),ir=1,nrs),((
     .      ksb(ik,ir),kpb(ik,ir),
     .    ik=0,nks(ir)),ir=1,nrs)

      ELSEIF (slver.GE.1.01) THEN
        READ(fp,ERR=98,END=99) (
     .      nks(ir),ikmids(ir),ksmaxs(ir),idds(ir,1),idds(ir,2),
     .    ir=1,nrs),(
     .      thetat(id),ikds(id),irds(id),
     .    id=1,nds),((
     .      rs    (ik,ir),zs    (ik,ir),kss   (ik,ir),kps  (ik,ir),
     .      kbacds(ik,ir),kfords(ik,ir),thetag(ik,ir),korpg(ik,ir),
     .      nvertp(korpg(ik,ir)),
     .      (rvertp(ii,korpg(ik,ir)),
     .       zvertp(ii,korpg(ik,ir)),ii=1,nvertp(korpg(ik,ir))),
     .    ik=1,nks(ir)),ir=1,nrs),((
     .      ksb(ik,ir),kpb(ik,ir),
     .    ik=0,nks(ir)),ir=1,nrs)
      ELSE
        READ(fp,ERR=98,END=99) (
     .    nks(ir),ikmids(ir),
     .    ir=1,nrs)

        READ(fp,ERR=98,END=99) ((
     .    rs(ik,ir),zs(ik,ir),
     .    korpg(ik,ir),nvertp(korpg(ik,ir)),
     .    (rvertp(ii,korpg(ik,ir)),
     .     zvertp(ii,korpg(ik,ir)),ii=1,nvertp(korpg(ik,ir))),
     .    ik=1,nks(ir)),ir=1,nrs)
      ENDIF

c...  Restore saved neutral wall:
      IF (.FALSE.) THEN
        nvesm = hold_nvesm
        nvesp = hold_nvesp
        DO i1 = 1, nvesm+nvesp
          jvesm(i1)   = hold_jvesm(i1)
          rvesm(i1,1) = hold_rvesm(i1,1)
          rvesm(i1,2) = hold_rvesm(i1,2)
          zvesm(i1,1) = hold_zvesm(i1,1)
          zvesm(i1,2) = hold_zvesm(i1,2)
        ENDDO
      ENDIF

c      CLOSE(fp)

      RETURN
98    error = 1
      RETURN
99    error = 2
      RETURN
      END
c
c ======================================================================
c
c subroutine: ReadPlasma
c
      SUBROUTINE ReadPlasma(fp,error)
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER fp,error

      REAL    slver
      INTEGER ik,ir,id,i1

      INTEGER s21_idumi,          s21_idumo,idum1
      REAL    s21_rdumi(MAXNRS,9),s21_rdumo(MAXNRS,9)


c...  Open file:
c      fp = 98
c      OPEN(UNIT=fp,FILE='plasma.dat',FORM='UNFORMATTED',
c     .     STATUS='OLD',ERR=98)


      READ (fp,ERR=98,END=99) slver

      IF (slver.GE.1.01) THEN
        READ (fp,ERR=98,END=99)
     .    nitersol,rel_opt ,rel_frac ,rel_nstep,rel_niter,
     .    rel_step,rel_iter,rel_count,adp_opt
      ELSE
        READ (fp,ERR=98,END=99)
     .    nitersol,rel_opt,rel_frac,rel_nstep,rel_niter,rel_step,
     .    rel_iter,rel_count
      ENDIF

      IF (slver.LE.1.02) THEN
        READ (fp,ERR=98,END=99)
     .    nrs,nds,ndsin,irwall,irtrap,irbreak,nbr

        READ (fp,ERR=98,END=99) (
     .    nks(ir),ikmids(ir),idds(ir,1),idds(ir,2),
     .    ksb(0,ir),ksb(nks(ir),ir),ksmaxs(ir),
     .    ir=1,nrs),(
     .    kteds(id),ktids(id),knds(id),kvds(id),thetat(id),
     .    id=1,nds)
      ENDIF

      IF     (slver.GE.1.15) THEN
        READ (fp,ERR=98,END=99) nds,(
     .      kteds(id),ktids(id),knds(id),kvds(id),keds(id),
     .    id=1,nds),
     .    nrs,(nks(ir),cmachno(ir,1),cmachno(ir,2),
     .      s28ionfrac(IKLO,ir),s28recfrac(IKLO,ir),s28momfrac(IKLO,ir),
     .      s28ionfrac(IKHI,ir),s28recfrac(IKHI,ir),s28momfrac(IKHI,ir),
     .      rel_hte  (ir),rel_hti  (ir),rel_hne  (ir)  ,rel_deltati(ir),
     .      rel_dirtg(ir),rel_symfr(ir),rel_prbfr(1,ir),rel_prbfr(2,ir),
     .      osm_sympt(ir),osm_peimul(1,ir),osm_peimul(2,ir),
     .      (rel_hproe(i1,ir),rel_hproi(i1,ir),i1=1,3),
     .      rel_qemul(IKLO,ir),rel_qemul(IKHI,ir),
     .      osm_model(IKLO,ir),osm_model(IKHI,ir),
     .      osm_code (IKLO,ir),osm_code (IKHI,ir),     
     .      ikbound  (ir,IKLO),ikbound  (ir,IKHI),
     .     (kpress (ik,ir,1),kpress(ik,ir,2),
     .      ktebs  (ik,ir)  ,ktibs (ik,ir)  ,knbs(ik,ir),kvhs(ik,ir),
     .      kes    (ik,ir)  ,
     .      osm_dp6(ik,ir)  ,
     .    ik=1,nks(ir)),ir=1,nrs),
     .    idum1,s21_idumi,(s21_rdumi(i1,7),i1=1,s21_idumi),
     .          s21_idumo,(s21_rdumo(i1,7),i1=1,s21_idumo),
     .    tarshift(IKLO),tarshift(IKHI),
     .    te_mult_o,ti_mult_o,n_mult_o,te_mult_i,ti_mult_i,n_mult_i,
     .    s28ionset,s28recset,s28momset,
     .    s28ionsetpfz,s28recsetpfz,s28momsetpfz

c... Still necessary?
        IF (idum1.EQ.1) THEN
          WRITE(0,*) 'OVERWRITING S21_DATAx7'
          DO i1 = 1, s21_idumo
            WRITE(PINOUT,'(A,I4,2F10.4)') 'I S21_7: ',
     .        i1,s21_datao(i1,7),s21_rdumo(i1,7)
            s21_datao(i1,7) = s21_rdumo(i1,7)
          ENDDO
          DO i1 = 1, s21_idumi
            WRITE(PINOUT,'(A,I4,2F10.4)') 'O S21_7: ',
     .        i1,s21_datai(i1,7),s21_rdumi(i1,7)
            s21_datai(i1,7) = s21_rdumi(i1,7)
          ENDDO
        ENDIF

      ELSEIF (slver.GE.1.14) THEN
        READ (fp,ERR=98,END=99) nds,(
     .      kteds(id),ktids(id),knds(id),kvds(id),keds(id),
     .    id=1,nds),nrs,(nks(ir),cmachno(ir,1),cmachno(ir,2),
     .      rel_hte  (ir),rel_hti  (ir),rel_hne  (ir)  ,rel_deltati(ir),
     .      rel_dirtg(ir),rel_symfr(ir),rel_prbfr(1,ir),rel_prbfr(2,ir),
     .      osm_sympt(ir),osm_peimul(1,ir),osm_peimul(2,ir),
     .      (rel_hproe(i1,ir),rel_hproi(i1,ir),i1=1,3),
     .      rel_qemul(IKLO,ir),rel_qemul(IKHI,ir),
     .      osm_model(IKLO,ir),osm_model(IKHI,ir),
     .      osm_code (IKLO,ir),osm_code (IKHI,ir),     
     .      ikbound  (ir,IKLO),ikbound  (ir,IKHI),
     .     (kpress (ik,ir,1),kpress(ik,ir,2),
     .      ktebs  (ik,ir)  ,ktibs (ik,ir)  ,knbs(ik,ir),kvhs(ik,ir),
     .      kes    (ik,ir)  ,
     .      osm_dp6(ik,ir)  ,
     .    ik=1,nks(ir)),ir=1,nrs),
     .    idum1,s21_idumi,(s21_rdumi(i1,7),i1=1,s21_idumi),
     .          s21_idumo,(s21_rdumo(i1,7),i1=1,s21_idumo),
     .    tarshift(IKLO),tarshift(IKHI),
     .    te_mult_o,ti_mult_o,n_mult_o,te_mult_i,ti_mult_i,n_mult_i

c... Still necessary?
        IF (idum1.EQ.1) THEN
          WRITE(0,*) 'OVERWRITING S21_DATAx7'
          DO i1 = 1, s21_idumo
            WRITE(PINOUT,'(A,I4,2F10.4)') 'I S21_7: ',
     .        i1,s21_datao(i1,7),s21_rdumo(i1,7)
            s21_datao(i1,7) = s21_rdumo(i1,7)
          ENDDO
          DO i1 = 1, s21_idumi
            WRITE(PINOUT,'(A,I4,2F10.4)') 'O S21_7: ',
     .        i1,s21_datai(i1,7),s21_rdumi(i1,7)
            s21_datai(i1,7) = s21_rdumi(i1,7)
          ENDDO
        ENDIF

      ELSEIF (slver.GE.1.13) THEN
        READ (fp,ERR=98,END=99) nds,(
     .      kteds(id),ktids(id),knds(id),kvds(id),keds(id),
     .    id=1,nds),nrs,(nks(ir),
     .      rel_hte  (ir),rel_hti  (ir),rel_hne  (ir)  ,rel_deltati(ir),
     .      rel_dirtg(ir),rel_symfr(ir),rel_prbfr(1,ir),rel_prbfr(2,ir),
     .      osm_sympt(ir),osm_peimul(1,ir),osm_peimul(2,ir),
     .      (rel_hproe(i1,ir),rel_hproi(i1,ir),i1=1,3),
     .      rel_qemul(IKLO,ir),rel_qemul(IKHI,ir),
     .      osm_model(IKLO,ir),osm_model(IKHI,ir),
     .      osm_code (IKLO,ir),osm_code (IKHI,ir),     
     .      ikbound  (ir,IKLO),ikbound  (ir,IKHI),
     .     (kpress (ik,ir,1),kpress(ik,ir,2),
     .      ktebs  (ik,ir)  ,ktibs (ik,ir)  ,knbs(ik,ir),kvhs(ik,ir),
     .      kes    (ik,ir)  ,
     .      osm_dp6(ik,ir)  ,
     .    ik=1,nks(ir)),ir=1,nrs),
     .    idum1,s21_idumi,(s21_rdumi(i1,7),i1=1,s21_idumi),
     .          s21_idumo,(s21_rdumo(i1,7),i1=1,s21_idumo),
     .    tarshift(IKLO),tarshift(IKHI),
     .    te_mult_o,ti_mult_o,n_mult_o,te_mult_i,ti_mult_i,n_mult_i

        IF (idum1.EQ.1) THEN
          WRITE(0,*) 'OVERWRITING S21_DATAx7'
          DO i1 = 1, s21_idumo
            WRITE(PINOUT,'(A,I4,2F10.4)') 'I S21_7: ',
     .        i1,s21_datao(i1,7),s21_rdumo(i1,7)
            s21_datao(i1,7) = s21_rdumo(i1,7)
          ENDDO
          DO i1 = 1, s21_idumi
            WRITE(PINOUT,'(A,I4,2F10.4)') 'O S21_7: ',
     .        i1,s21_datai(i1,7),s21_rdumi(i1,7)
            s21_datai(i1,7) = s21_rdumi(i1,7)
          ENDDO
        ENDIF
      ELSEIF (slver.GE.1.12) THEN
        READ (fp,ERR=98,END=99) nds,(
     .      kteds(id),ktids(id),knds(id),kvds(id),keds(id),
     .    id=1,nds),nrs,(nks(ir),
     .      rel_hte  (ir),rel_hti  (ir),rel_hne  (ir)  ,rel_deltati(ir),
     .      rel_dirtg(ir),rel_symfr(ir),rel_prbfr(1,ir),rel_prbfr(2,ir),
     .      osm_sympt(ir),osm_peimul(1,ir),osm_peimul(2,ir),
     .      (rel_hproe(i1,ir),rel_hproi(i1,ir),i1=1,3),
     .      rel_qemul(IKLO,ir),rel_qemul(IKHI,ir),
     .      osm_model(IKLO,ir),osm_model(IKHI,ir),
     .      osm_code (IKLO,ir),osm_code (IKHI,ir),     
     .      ikbound  (ir,IKLO),ikbound  (ir,IKHI),
     .     (kpress (ik,ir,1),kpress(ik,ir,2),
     .      ktebs  (ik,ir)  ,ktibs (ik,ir)  ,knbs(ik,ir),kvhs(ik,ir),
     .      kes    (ik,ir)  ,
     .      osm_dp6(ik,ir)  ,
     .    ik=1,nks(ir)),ir=1,nrs),
     .    idum1,s21_idumi,(s21_rdumi(i1,7),i1=1,s21_idumi),
     .          s21_idumo,(s21_rdumo(i1,7),i1=1,s21_idumo)

        IF (idum1.EQ.1) THEN
          WRITE(0,*) 'OVERWRITING S21_DATAx7'
          DO i1 = 1, s21_idumo
            WRITE(PINOUT,'(A,I4,2F10.4)') 'I S21_7: ',
     .        i1,s21_datao(i1,7),s21_rdumo(i1,7)
            s21_datao(i1,7) = s21_rdumo(i1,7)
          ENDDO
          DO i1 = 1, s21_idumi
            WRITE(PINOUT,'(A,I4,2F10.4)') 'O S21_7: ',
     .        i1,s21_datai(i1,7),s21_rdumi(i1,7)
            s21_datai(i1,7) = s21_rdumi(i1,7)
          ENDDO
        ENDIF

      ELSEIF (slver.GE.1.11) THEN
        READ (fp,ERR=98,END=99) nds,(
     .      kteds(id),ktids(id),knds(id),kvds(id),keds(id),
     .    id=1,nds),nrs,(nks(ir),
     .      rel_hte  (ir),rel_hti  (ir),rel_hne  (ir)  ,rel_deltati(ir),
     .      rel_dirtg(ir),rel_symfr(ir),rel_prbfr(1,ir),rel_prbfr(2,ir),
     .      osm_sympt(ir),osm_peimul(1,ir),osm_peimul(2,ir),
     .      (rel_hproe(i1,ir),rel_hproi(i1,ir),i1=1,3),
     .      rel_qemul(IKLO,ir),rel_qemul(IKHI,ir),
     .      osm_model(IKLO,ir),osm_model(IKHI,ir),
     .      ikbound  (ir,IKLO),ikbound  (ir,IKHI),
     .     (kpress (ik,ir,1),kpress(ik,ir,2),
     .      ktebs  (ik,ir)  ,ktibs (ik,ir)  ,knbs(ik,ir),kvhs(ik,ir),
     .      kes    (ik,ir)  ,
     .      osm_dp6(ik,ir)  ,
     .    ik=1,nks(ir)),ir=1,nrs),
     .    idum1,s21_idumi,(s21_rdumi(i1,7),i1=1,s21_idumi),
     .          s21_idumo,(s21_rdumo(i1,7),i1=1,s21_idumo)

        IF (idum1.EQ.1) THEN
          WRITE(0,*) 'OVERWRITING S21_DATAx7'
          DO i1 = 1, s21_idumo
            WRITE(PINOUT,'(A,I4,2F10.4)') 'I S21_7: ',
     .        i1,s21_datao(i1,7),s21_rdumo(i1,7)
            s21_datao(i1,7) = s21_rdumo(i1,7)
          ENDDO
          DO i1 = 1, s21_idumi
            WRITE(PINOUT,'(A,I4,2F10.4)') 'O S21_7: ',
     .        i1,s21_datai(i1,7),s21_rdumi(i1,7)
            s21_datai(i1,7) = s21_rdumi(i1,7)
          ENDDO
        ENDIF

      ELSEIF (slver.GE.1.10) THEN
        READ (fp,ERR=98,END=99) nds,(
     .      kteds(id),ktids(id),knds(id),kvds(id),keds(id),
     .    id=1,nds),nrs,(nks(ir),
     .      rel_hte  (ir),rel_hti  (ir),rel_hne  (ir)  ,rel_deltati(ir),
     .      rel_dirtg(ir),rel_symfr(ir),rel_prbfr(1,ir),rel_prbfr(2,ir),
     .      osm_sympt(ir),osm_peimul(1,ir),osm_peimul(2,ir),
     .      (rel_hproe(i1,ir),rel_hproi(i1,ir),i1=1,3),
     .      rel_qemul(IKLO,ir),rel_qemul(IKHI,ir),
     .      osm_model(IKLO,ir),osm_model(IKHI,ir),
     .      ikbound  (ir,IKLO),ikbound  (ir,IKHI),
     .     (kpress (ik,ir,1),kpress(ik,ir,2),
     .      ktebs  (ik,ir)  ,ktibs (ik,ir)  ,knbs(ik,ir),kvhs(ik,ir),
     .      kes    (ik,ir)  ,
     .      osm_dp6(ik,ir)  ,
     .    ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.GE.1.09) THEN
        READ (fp,ERR=98,END=99) nds,(
     .      kteds(id),ktids(id),knds(id),kvds(id),keds(id),
     .    id=1,nds),nrs,(nks(ir),
     .      rel_hte  (ir),rel_hti  (ir),rel_hne  (ir)  ,rel_deltati(ir),
     .      rel_dirtg(ir),rel_symfr(ir),rel_prbfr(1,ir),rel_prbfr(2,ir),
     .      osm_sympt(ir),osm_peimul(1,ir),osm_peimul(2,ir),
     .      (rel_hproe(i1,ir),rel_hproi(i1,ir),i1=1,3),
     .      rel_qemul(IKLO,ir),rel_qemul(IKHI,ir),
     .     (kpress (ik,ir,1),kpress(ik,ir,2),
     .      ktebs  (ik,ir)  ,ktibs (ik,ir)  ,knbs(ik,ir),kvhs(ik,ir),
     .      kes    (ik,ir)  ,
     .      osm_dp6(ik,ir)  ,
     .    ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.GE.1.08) THEN
        READ (fp,ERR=98,END=99) nds,(
     .      kteds(id),ktids(id),knds(id),kvds(id),keds(id),
     .    id=1,nds),nrs,(nks(ir),
     .      rel_hte  (ir),rel_hti  (ir),rel_hne  (ir)  ,rel_deltati(ir),
     .      rel_dirtg(ir),rel_symfr(ir),rel_prbfr(1,ir),rel_prbfr(2,ir),
     .      osm_sympt(ir),osm_peimul(1,ir),osm_peimul(2,ir),
     .      (rel_hproe(i1,ir),rel_hproi(i1,ir),i1=1,3),
     .      rel_qemul(IKLO,ir),rel_qemul(IKHI,ir),
     .     (kpress(ik,ir,1),kpress(ik,ir,2),
     .      ktebs (ik,ir)  ,ktibs (ik,ir)  ,knbs(ik,ir),kvhs(ik,ir),
     .      kes   (ik,ir)  ,
     .    ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.GE.1.07) THEN
        READ (fp,ERR=98,END=99) nds,(
     .      kteds(id),ktids(id),knds(id),kvds(id),keds(id),
     .    id=1,nds),nrs,(nks(ir),
     .      rel_hte  (ir),rel_hti  (ir),rel_hne  (ir)  ,rel_deltati(ir),
     .      rel_dirtg(ir),rel_symfr(ir),rel_prbfr(1,ir),rel_prbfr(2,ir),
     .      osm_sympt(ir),
     .      (rel_hproe(i1,ir),rel_hproi(i1,ir),i1=1,3),
     .      rel_qemul(IKLO,ir),rel_qemul(IKHI,ir),
     .     (kpress(ik,ir,1),kpress(ik,ir,2),
     .      ktebs (ik,ir)  ,ktibs (ik,ir)  ,knbs(ik,ir),kvhs(ik,ir),
     .      kes   (ik,ir)  ,
     .    ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.GE.1.06) THEN
        READ (fp,ERR=98,END=99) nds,(
     .    kteds(id),ktids(id),knds(id),kvds(id),keds(id),
     .  id=1,nds),nrs,(nks(ir),
     .    rel_hte  (ir),rel_hti  (ir),rel_hne(ir),rel_deltati(ir),
     .    rel_dirtg(ir),rel_symfr(ir),
     .    (rel_hproe(i1,ir),rel_hproi(i1,ir),i1=1,3),
     .    rel_qemul(IKLO,ir),rel_qemul(IKHI,ir),
     .   (kpress(ik,ir,1),kpress(ik,ir,2),
     .    ktebs (ik,ir)  ,ktibs (ik,ir)  ,knbs(ik,ir),kvhs(ik,ir),
     .    kes   (ik,ir)  ,
     .  ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.GE.1.05) THEN
        READ (fp,ERR=98,END=99) nds,(
     .    kteds(id),ktids(id),knds(id),kvds(id),keds(id),
     .  id=1,nds),nrs,(nks(ir),
     .    rel_hte  (ir),rel_hti  (ir),rel_hne(ir),rel_deltati(ir),
     .    rel_dirtg(ir),rel_symfr(ir),
     .    (rel_hproe(i1,ir),rel_hproi(i1,ir),i1=1,3),(
     .    kpress(ik,ir,1),kpress(ik,ir,2),
     .    ktebs (ik,ir)  ,ktibs (ik,ir)  ,knbs(ik,ir),kvhs(ik,ir),
     .    kes   (ik,ir)  ,
     .  ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.GE.1.04) THEN
        READ (fp,ERR=98,END=99) nds,(
     .    kteds(id),ktids(id),knds(id),kvds(id),keds(id),
     .  id=1,nds),nrs,(nks(ir),(
     .    kpress(ik,ir,1),kpress(ik,ir,2),
     .    ktebs (ik,ir)  ,ktibs (ik,ir)  ,knbs(ik,ir),kvhs(ik,ir),
     .    kes   (ik,ir)  ,
     .  ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.GE.1.03) THEN
        READ (fp,ERR=98,END=99) nds,(
     .      kteds(id),ktids(id),knds(id),kvds(id),
     .    id=1,nds),nrs,(nks(ir),(
     .      kpress(ik,ir,1),kpress(ik,ir,2),
     .      ktebs (ik,ir)  ,ktibs (ik,ir)  ,knbs(ik,ir),kvhs(ik,ir),
     .    ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.GE.1.02) THEN
        READ (fp,ERR=98,END=99) ((
     .    kss   (ik,ir),kps   (ik,ir),thetag(ik,ir),
     .    kbacds(ik,ir),kfords(ik,ir),kpress(ik,ir,1),kpress(ik,ir,2),
     .    ktebs (ik,ir),ktibs (ik,ir),knbs  (ik,ir)  ,kvhs  (ik,ir)  ,
     .    ik=1,nks(ir)),ir=1,nrs)
      ELSE
        READ (fp,ERR=98,END=99) ((
     .    kss   (ik,ir),kps   (ik,ir),thetag(ik,ir),
     .    kbacds(ik,ir),kfords(ik,ir),
     .    ktebs (ik,ir),ktibs (ik,ir),knbs(ik,ir),kvhs(ik,ir),
     .    ik=1,nks(ir)),ir=1,nrs)
      ENDIF


c      CLOSE(fp)


      RETURN
98    error = 1
      RETURN
99    error = 2
      RETURN
      END
c
c ======================================================================
c
c subroutine: ReadSources
c
      SUBROUTINE ReadSources(fp,error)
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER fp,error

      REAL    slver
      INTEGER ik,ir,id,i1,i2,i3,i4,idum1,idum2,idum3,idum4

c...  Sucks, but for backward compatability (MAKE DYNAMIC!):
      REAL pinasd2(2000,MAXASD2,MAXASS,3) 

      error = 0

c.... Open file:
c      fp = 98
c      OPEN(UNIT=fp,FILE='source.dat',FORM='UNFORMATTED',
c     .     STATUS='OLD',ERR=98)

      READ (fp,ERR=98,END=99) slver

      IF (slver.GE.1.01) THEN
        READ (fp,ERR=98,END=99)
     .    nitersol,rel_opt ,rel_frac ,rel_nstep,rel_niter,
     .    rel_step,rel_iter,rel_count,adp_opt
      ELSE
        READ (fp,ERR=98,END=99)
     .    nitersol,rel_opt,rel_frac,rel_nstep,rel_niter,rel_step,
     .    rel_iter
      ENDIF

      IF (slver.LE.1.02) THEN
        READ (fp,ERR=98,END=99)
     .    nrs,nds,ndsin,irwall,irtrap,irbreak,nbr

        READ (fp,ERR=98,END=99) (
     .    nks(ir),ikmids(ir),idds(ir,1),idds(ir,2),
     .    ksb(0,ir),ksb(nks(ir),ir),ksmaxs(ir),
     .    ir=1,nrs),(
     .    kteds(id),ktids(id),knds(id),kvds(id),thetat(id),
     .    id=1,nds)
      ENDIF

      IF     (slver.GE.1.21) THEN

        tagpinatom = .TRUE.

        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .      pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .      pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .      pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .      osmcfe (ik,ir),osmmp (ik,ir),osmqe   (ik,ir),
     .      osmcve (ik,ir),osmcvi(ik,ir),osmcde  (ik,ir),osmcdi(ik,ir),
     .      osmion (ik,ir),osmrec(ik,ir),
     .      (osmcfpflx(ik,ir,i1),i1=1,5),
     .      pinqir (ik,ir),pinqer(ik,ir),
     .      idum1  ,(pinline(ik,ir,i1,H_BALPHA),
     .               pinline(ik,ir,i1,H_BGAMMA)   ,i1=1,idum1),
     .      idum2,idum1,
     .        ((pinstrata(ik,ir,i2,i1),i2=1,idum2),i1=1,idum1),
     .      idum1,(pinploss(ik,ir,i1)             ,i1=1,idum1),
     .      idum1,(pinbgk  (ik,ir,i1)             ,i1=1,idum1),
     .    ik=1,nks(ir)),
     .                (nks(ir),osmpmk(ik,ir),
     .    ik=1,nks(ir)+1),
     .    ir=1,nrs),
     .    idum1,idum2,((eirpgdat(i1,i2),i1=1,idum1),i2=1,idum2),
     .    indasd,idum1,idum2,idum3,idum4,
     .    ((((pinasd(i1,i2,i3,i4),i1=1,idum1),i2=1,idum2),
     .                            i3=1,idum3),i4=1,idum4)

      ELSEIF (slver.GE.1.20) THEN

        tagpinatom = .TRUE.

        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .      pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .      pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .      pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .      osmcfe (ik,ir),osmmp (ik,ir),osmqe   (ik,ir),
     .      osmcve (ik,ir),osmcvi(ik,ir),osmcde  (ik,ir),osmcdi(ik,ir),
     .      osmion (ik,ir),osmrec(ik,ir),
     .      pinqir (ik,ir),pinqer(ik,ir),
     .      idum1  ,(pinline(ik,ir,i1,H_BALPHA),
     .               pinline(ik,ir,i1,H_BGAMMA)   ,i1=1,idum1),
     .      idum2,idum1,
     .        ((pinstrata(ik,ir,i2,i1),i2=1,idum2),i1=1,idum1),
     .      idum1,(pinploss(ik,ir,i1)             ,i1=1,idum1),
     .      idum1,(pinbgk  (ik,ir,i1)             ,i1=1,idum1),
     .    ik=1,nks(ir)),
     .                (nks(ir),osmpmk(ik,ir),
     .    ik=1,nks(ir)+1),
     .    ir=1,nrs),
     .    idum1,idum2,((eirpgdat(i1,i2),i1=1,idum1),i2=1,idum2),
     .    indasd,idum1,idum2,idum3,idum4,
     .    ((((pinasd(i1,i2,i3,i4),i1=1,idum1),i2=1,idum2),
     .                            i3=1,idum3),i4=1,idum4)

      ELSEIF (slver.GE.1.19) THEN

        tagpinatom = .TRUE.

        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .      pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .      pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .      pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .      osmcfe (ik,ir),osmmp (ik,ir),osmqe   (ik,ir),
     .      osmcve (ik,ir),osmcvi(ik,ir),osmcde  (ik,ir),osmcdi(ik,ir),
     .      pinqir (ik,ir),pinqer(ik,ir),
     .      idum1  ,(pinline(ik,ir,i1,H_BALPHA),
     .               pinline(ik,ir,i1,H_BGAMMA)   ,i1=1,idum1),
     .      idum2,idum1,
     .        ((pinstrata(ik,ir,i2,i1),i2=1,idum2),i1=1,idum1),
     .      idum1,(pinploss(ik,ir,i1)             ,i1=1,idum1),
     .      idum1,(pinbgk  (ik,ir,i1)             ,i1=1,idum1),
     .    ik=1,nks(ir)),
     .                (nks(ir),osmpmk(ik,ir),
     .    ik=1,nks(ir)+1),
     .    ir=1,nrs),
     .    idum1,idum2,((eirpgdat(i1,i2),i1=1,idum1),i2=1,idum2),
     .    indasd,idum1,idum2,idum3,idum4,
     .    ((((pinasd(i1,i2,i3,i4),i1=1,idum1),i2=1,idum2),
     .                            i3=1,idum3),i4=1,idum4)

      ELSEIF (slver.GE.1.18) THEN

        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .      pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .      pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .      pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .      osmcfe (ik,ir),osmmp (ik,ir),osmqe   (ik,ir),
     .      osmcve (ik,ir),osmcvi(ik,ir),osmcde  (ik,ir),osmcdi(ik,ir),
     .      pinqir (ik,ir),pinqer(ik,ir),
     .      idum1  ,(pinline(ik,ir,i1,H_BALPHA),
     .               pinline(ik,ir,i1,H_BGAMMA)   ,i1=1,idum1),
     .      idum2,idum1,
     .        ((pinstrata(ik,ir,i2,i1),i2=1,idum2),i1=1,idum1),
     .      idum1,(pinploss(ik,ir,i1)             ,i1=1,idum1),
     .      idum1,(pinbgk  (ik,ir,i1)             ,i1=1,idum1),
     .    ik=1,nks(ir)),
     .                (nks(ir),osmpmk(ik,ir),
     .    ik=1,nks(ir)+1),
     .    ir=1,nrs),
     .    idum1,idum2,((eirpgdat(i1,i2),i1=1,idum1),i2=1,idum2),
     .    indasd,idum1,idum2,idum3,idum4,
     .    ((((pinasd2(i1,i2,i3,i4),i1=1,idum1),i2=1,idum2),
     .                             i3=1,idum3),i4=1,idum4)

          WRITE(0,*) 'IDUM1=',idum1


      ELSEIF (slver.GE.1.17) THEN
        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .      pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .      pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .      pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .      osmcfe (ik,ir),osmmp (ik,ir),osmqe   (ik,ir),
     .      osmcve (ik,ir),osmcvi(ik,ir),osmcde  (ik,ir),osmcdi(ik,ir),
     .      pinqir (ik,ir),pinqer(ik,ir),
     .      idum1  ,(pinline(ik,ir,i1,H_BALPHA),
     .               pinline(ik,ir,i1,H_BGAMMA)   ,i1=1,idum1),
     .      idum2,idum1,
     .        ((pinstrata(ik,ir,i2,i1),i2=1,idum2),i1=1,idum1),
     .      idum1,(pinploss(ik,ir,i1)             ,i1=1,idum1),
     .      idum1,(pinbgk  (ik,ir,i1)             ,i1=1,idum1),
     .    ik=1,nks(ir)),ir=1,nrs),
     .    idum1,idum2,((eirpgdat(i1,i2),i1=1,idum1),i2=1,idum2),
     .    indasd,idum1,idum2,idum3,idum4,
     .    ((((pinasd2(i1,i2,i3,i4),i1=1,idum1),i2=1,idum2),
     .                             i3=1,idum3),i4=1,idum4)

      ELSEIF (slver.GE.1.16) THEN

        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .      pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .      pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .      pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .      osmcfe (ik,ir),osmmp (ik,ir),osmqe   (ik,ir),
     .      osmcve (ik,ir),osmcvi(ik,ir),osmcde  (ik,ir),osmcdi(ik,ir),
     .      pinqir (ik,ir),pinqer(ik,ir),
     .      idum1,(pinline(ik,ir,i1,H_BALPHA),
     .             pinline(ik,ir,i1,H_BGAMMA),i1=1,idum1),
     .      idum1,(pindata(ik,ir,i1)         ,i1=1,idum1),
     .      idum1,(pinbgk (ik,ir,i1)         ,i1=1,idum1),
     .    ik=1,nks(ir)),ir=1,nrs),
     .    idum1,idum2,((eirpgdat(i1,i2),i1=1,idum1),i2=1,idum2),
     .    indasd,idum1,idum2,idum3,idum4,
     .    ((((pinasd2(i1,i2,i3,i4),i1=1,idum1),i2=1,idum2),
     .                             i3=1,idum3),i4=1,idum4)

      ELSEIF (slver.GE.1.15) THEN
        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .      pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .      pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .      pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .      osmcfe (ik,ir),osmmp (ik,ir),osmqe   (ik,ir),
     .      osmcve (ik,ir),osmcvi(ik,ir),osmcde  (ik,ir),osmcdi(ik,ir),
     .      idum1,(pinline(ik,ir,i1,H_BALPHA),
     .             pinline(ik,ir,i1,H_BGAMMA),i1=1,idum1),
     .      idum1,(pindata(ik,ir,i1)         ,i1=1,idum1),
     .      idum1,(pinbgk (ik,ir,i1)         ,i1=1,idum1),
     .    ik=1,nks(ir)),ir=1,nrs),
     .    idum1,idum2,((eirpgdat(i1,i2),i1=1,idum1),i2=1,idum2),
     .    indasd,idum1,idum2,idum3,idum4,
     .    ((((pinasd2(i1,i2,i3,i4),i1=1,idum1),i2=1,idum2),
     .                             i3=1,idum3),i4=1,idum4)
      ELSEIF (slver.GE.1.14) THEN
        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .    pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .    pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .    pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .    osmcfe (ik,ir),osmmp (ik,ir),osmqe   (ik,ir),
     .    osmcve (ik,ir),osmcvi(ik,ir),osmcde  (ik,ir),osmcdi(ik,ir),
     .    (pinline(ik,ir,i1,H_BALPHA),
     .     pinline(ik,ir,i1,H_BGAMMA),i1=1,6      ),
     .    idum1,(pindata(ik,ir,i1)         ,i1=1,idum1),
     .    ik=1,nks(ir)),ir=1,nrs),
     .    idum1,idum2,((eirpgdat(i1,i2),i1=1,idum1),i2=1,idum2),
     .    indasd,idum1,idum2,idum3,idum4,
     .    ((((pinasd2(i1,i2,i3,i4),i1=1,idum1),i2=1,idum2),
     .                             i3=1,idum3),i4=1,idum4)
      ELSEIF (slver.GE.1.13) THEN
        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .    pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .    pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .    pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .    osmcfe (ik,ir),osmmp (ik,ir),osmqe   (ik,ir),
     .    osmcve (ik,ir),osmcvi(ik,ir),osmcde  (ik,ir),osmcdi(ik,ir),
     .    (pinline(ik,ir,i1,H_BALPHA),
     .     pinline(ik,ir,i1,H_BGAMMA),i1=1,6      ),
     .    (pindata(ik,ir,i1)         ,i1=1,MAXDATA),
     .    ik=1,nks(ir)),ir=1,nrs),
     .    ((eirpgdat(i1,i2),i1=1,MAXNAS),i2=1,MAXASD)
      ELSEIF (slver.GE.1.12) THEN
        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .    pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .    pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .    pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .    osmcfe (ik,ir),osmmp (ik,ir),osmqe   (ik,ir),
     .    osmcve (ik,ir),osmcvi(ik,ir),osmcde  (ik,ir),osmcdi(ik,ir),
     .    (pinline(ik,ir,i1,H_BALPHA),
     .     pinline(ik,ir,i1,H_BGAMMA),i1=1,6      ),
     .    (pindata(ik,ir,H_MP1+i1-1) ,i1=1,NMOMCHA),
     .    ik=1,nks(ir)),ir=1,nrs),
     .    ((eirpgdat(i1,i2),i1=1,MAXNAS),i2=1,MAXASD)
      ELSEIF (slver.GE.1.11) THEN
        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .    pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .    pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .    pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .    osmcfe (ik,ir),osmmp (ik,ir),osmqe   (ik,ir),
     .    osmcve (ik,ir),osmcvi(ik,ir),osmcde  (ik,ir),osmcdi(ik,ir),
     .    (pinline(ik,ir,i1,H_BALPHA),
     .     pinline(ik,ir,i1,H_BGAMMA),i1=1,6      ),
     .    (pindata(ik,ir,H_MP1+i1-1) ,i1=1,NMOMCHA),
     .    ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.GE.1.10) THEN
        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .    pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .    pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .    pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .    osmcfe (ik,ir),osmmp (ik,ir),osmqe   (ik,ir),
     .    osmcve (ik,ir),osmcvi(ik,ir),osmcde  (ik,ir),osmcdi(ik,ir),
     .    (pinline(ik,ir,i1,H_BALPHA),
     .     pinline(ik,ir,i1,H_BGAMMA),i1=1,6),
     .    ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.GE.1.09) THEN
        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .    pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .    pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .    pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .    osmcfe (ik,ir),osmmp (ik,ir),
     .    osmcve (ik,ir),osmcvi(ik,ir),osmcde  (ik,ir),osmcdi(ik,ir),
     .    (pinline(ik,ir,i1,H_BALPHA),
     .     pinline(ik,ir,i1,H_BGAMMA),i1=1,6),
     .    ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.GE.1.08) THEN
        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .    pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .    pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .    pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .    osmcfe (ik,ir),osmmp (ik,ir),
     .    osmcve (ik,ir),osmcvi(ik,ir),osmcde  (ik,ir),osmcdi(ik,ir),
     .    pinline(ik,ir,6,H_BALPHA),pinline(ik,ir,6,H_BGAMMA),
     .    ik=1,nks(ir)),ir=1,nrs)

      ELSEIF (slver.GE.1.07) THEN
        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .    pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .    pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .    pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .    osmcfe (ik,ir),
     .    pinline(ik,ir,6,H_BALPHA),pinline(ik,ir,6,H_BGAMMA),
     .    ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.EQ.1.06) THEN
        CALL RZero(pinline(1,1,6,H_BALPHA),MAXNKS*MAXNRS)
        CALL RZero(pinline(1,1,6,H_BGAMMA),MAXNKS*MAXNRS)

        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .    pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi (ik,ir),
     .    pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp (ik,ir),
     .    pinena (ik,ir),osmpei(ik,ir),osmcfp  (ik,ir),osmcfi(ik,ir),
     .    osmcfe (ik,ir),
     .    ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.EQ.1.05) THEN
        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .    pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi(ik,ir),
     .    pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp(ik,ir),
     .    pinena (ik,ir),osmpei(ik,ir),
     .    ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.GE.1.04) THEN
        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .    pinion (ik,ir),pinrec(ik,ir)  ,pinqe   (ik,ir),pinqi(ik,ir),
     .    pinatom(ik,ir),pinmol(ik,ir)  ,pinalpha(ik,ir),pinmp(ik,ir),
     .    pinena (ik,ir),
     .    ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.GE.1.03) THEN
        READ (fp,ERR=98,END=99) nrs,(nks(ir),(
     .    pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi(ik,ir),
     .    pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp(ik,ir),
     .    ik=1,nks(ir)),ir=1,nrs)
      ELSEIF (slver.GE.1.02) THEN
        READ (fp,ERR=98,END=99) ((
     .    kss    (ik,ir),kps   (ik,ir),thetag(ik,ir),
     .    kbacds (ik,ir),kfords(ik,ir),
     .    pinion (ik,ir),pinrec(ik,ir),pinqe   (ik,ir),pinqi(ik,ir),
     .    pinatom(ik,ir),pinmol(ik,ir),pinalpha(ik,ir),pinmp(ik,ir),
     .    ik=1,nks(ir)),ir=1,nrs)
      ELSE
        READ (fp,ERR=98,END=99) ((
     .    kss    (ik,ir),kps   (ik,ir),thetag(ik,ir),
     .    kbacds (ik,ir),kfords(ik,ir),
     .    pinion (ik,ir),pinrec(ik,ir),pinqe(ik,ir),pinqi(ik,ir),
     .    pinatom(ik,ir),pinmol(ik,ir),
     .    ik=1,nks(ir)),ir=1,nrs)
      ENDIF

c...  Sucks:
      IF (slver.GE.1.14.AND.slver.LE.1.18) THEN
        DO i1 = 1, idum1
          DO i2 = 1, MAXASD2
            DO i3 = 1, MAXASS         
              pinasd(i1,i2,i3,1) = pinasd2(i1,i2,i3,1)  
              pinasd(i1,i2,i3,2) = pinasd2(i1,i2,i3,3)  
            ENDDO
          ENDDO
        ENDDO  
      ENDIF

c      CLOSE(fp)

      RETURN
98    error = 1
      RETURN
99    error = 2
      RETURN
      END
c
c ======================================================================
c
c function: CalcPoint
c
c Calculate the perpendicular distance from a point to a line.
c
      INTEGER FUNCTION CalcPoint(ar,az,br,bz,cr,cz,t)

      IMPLICIT none

c     Input:
      REAL*8 ar,az,br,bz,cr,cz

c     Output:
      REAL*8 t

      REAL*8     TOL
      PARAMETER (TOL = 3.0E-7)

      REAL*8 r,z,deltar,deltaz

      CalcPoint = 0

      deltar = br - ar
      deltaz = bz - az

      IF (ABS(deltar).GT.TOL.OR.ABS(deltaz).GT.TOL) THEN

        t = ((cr - ar) * deltar + (cz - az) * deltaz) /
     .      (deltar**2 + deltaz**2)

        r = ar + t * deltar
        z = az + t * deltaz

c...NEW
        IF     (DSQRT((ar-cr)**2+(az-cz)**2).LT.TOL) THEN
           CalcPoint = 1
           t = 0.0D0
        ELSEIF (DSQRT((br-cr)**2+(bz-cz)**2).LT.TOL) THEN
           CalcPoint = 1
           t = 1.0D0
        ELSEIF ((t+TOL).GE.0.0.AND.(t-TOL).LE.1.0) THEN
c...OLD
c        IF ((t+TOL).GE.0.0.AND.(t-TOL).LE.1.0) THEN
          IF (SQRT((r - cr)**2 + (z - cz)**2).LT.TOL*10.0) THEN
c
c           Point C is on the line AB:
c
            CalcPoint = 1
          ELSE
c
c           Point of perpendicular intersection is displaced from
c           the point C:
c
c           WRITE(50,*) '  DIST = ',SQRT((r - cr)**2 + (z - cz)**2)
c           WRITE(50,*) '  R,Z  = ',r,z

            CalcPoint = 2
          ENDIF
        ELSEIF ((t-TOL).LE.1.0) THEN
          CalcPoint = 3
        ENDIF
      ELSE
c
c       If the points are all identical, then return a positive result,
c       otherwise indicate that the problem was ill-posed:
c
        IF (ABS(ar-cr).LT.TOL.AND.ABS(az-cz).LT.TOL) THEN
          CalcPoint =  1
        ELSE
          CalcPoint = -1
        ENDIF
      ENDIF

      RETURN
      END
c
c ======================================================================
c
c function: CalcPressure
c
c     Returns pressure in units of [eV m-3].
c
      REAL FUNCTION CalcPressure(n,te,ti,v)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'

      REAL n,te,ti,v

C     IPP/01 - Krieger: incredibly enough, the FUJI f90 compiler
C     chokes over v**2.0 if v is negative :-(
C     fixed by v**2.0 -> v*v or v**2

      CalcPressure = n * (te + ti + (crmb * AMU / ECH) * v**2)
c      CalcPressure = n * (te + ti + (crmb * AMU / ECH) * v**2.0)

      RETURN
99    STOP
      END














c
c ======================================================================
c
c subroutine: InterpolateProbeData
c
c
c Need to check LPDATSW, and adjust density accordingly...
c
c
      SUBROUTINE InterpolateProbeData(probe)
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER probe

      INTEGER ir

      IF (prb_num(probe).EQ.0) THEN
        CALL WN('InterpolateProbeData','Probe data not found')
        RETURN
      ENDIF

      IF (probe.NE.FSP1) THEN
c
c Currently, dip_v is only calculated for FSP1, so there will be an error in
c the SOLEDGE routine when calculating the upstream pressure.  Development
c is needed in PROBEPATH so that dip_v is available for other probes.
c
        CALL ER('InterpolateProbeData','Development required',*99)
      ENDIF

      CALL ProbePath

      WRITE(SLOUT,*)
      WRITE(SLOUT,*) 'INTERPOLATE FSP1 PROBE DATA TO GRID'

      DO ir = irsep, irwall-1
        CALL Fitter(prb_num(probe),prb_rho(1,probe),prb_te(1,probe),
     .              1,rho(ir,CELL1),prp_te(ir,probe),'LINEAR')
        CALL Fitter(prb_num(probe),prb_rho(1,probe),prb_ti(1,probe),
     .              1,rho(ir,CELL1),prp_ti(ir,probe),'LINEAR')
        CALL Fitter(prb_num(probe),prb_rho(1,probe),prb_ne(1,probe),
     .              1,rho(ir,CELL1),prp_ne(ir,probe),'LINEAR')

        IF (osm_recopt.EQ.6) THEN

c          prp_te(ir,probe) = prp_te(ir,probe) * rflexopt(3)
c          prp_ti(ir,probe) = prp_ti(ir,probe) * rflexopt(3)
          prp_ne(ir,probe) = prp_ne(ir,probe) * rflexopt(3)

        ENDIF




        WRITE(SLOUT,'(A,I4,I4,F15.7,1X,2F12.4,1P,E12.4)')
     .    'probe rho PRB Te,Ti,ne = ',
     .    probe,ir,rho(ir,CELL1),
     .    prp_te(ir,probe),prp_ti(ir,probe),prp_ne(ir,probe)

      ENDDO

      RETURN
c     Error code:
99    CONTINUE
      STOP
      END
c
c ======================================================================
c
c subroutine: ProbePath
c
c It is possible that there be two intesections of the probe path
c near each other on a given ring... this is unlikely... especially
c for the FSP... but an error check should be made somehow...
c
      SUBROUTINE ProbePath
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER mark,i1,ik,ir,ikprb
      DOUBLE PRECISION a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd,tmin,tprb


c... WARN THAT FSP1 is assumed!

      IF (prb_num(FSP1).EQ.0) THEN
        CALL WN('ProbePath','Probe data not available')
        RETURN
      ENDIF
c
c     Find starting point point and ending point for probe trajectory,
c     and assume that the probe moves in a straight line:
c
      mark = 0
      DO i1 = 1, prb_num(FSP1)
        IF (prb_r(i1,FSP1).LT.HI.AND.prb_z(i1,FSP1).LT.HI) THEN
          IF (mark.EQ.0) THEN
            a1 = prb_r(i1,FSP1)
            a2 = prb_z(i1,FSP1)
            mark = 1
c            WRITE(SLOUT,'(A,1P,2E15.7)') 'A1,A2 = ',a1,a2
         ELSE
            b1 = prb_r(i1,FSP1)
            b2 = prb_z(i1,FSP1)
c            WRITE(SLOUT,'(A,1P,2E15.7)') 'B1,B2 = ',b1,b2
          ENDIF
        ENDIF
      ENDDO

      DO ir = irsep, irwall-1
c        WRITE(SLOUT,*) ' '
        tmin  = HI
        ikprb = NULL
        DO ik = 1, nks(ir)-1
          c1 = rs(ik,ir)
          c2 = zs(ik,ir)
          d1 = rs(ik+1,ir)
          d2 = zs(ik+1,ir)

          CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)

c          WRITE(SLOUT,'(A,2I4,1P,2E15.7,::)')
c     .      'path = ',ik,ir,tab,tcd

          IF (tcd.GT.0.0.AND.tcd.LT.1.0.AND.ABS(tab).LT.tmin) THEN
            ikprb = ik
            tprb  = tcd
            tmin  = ABS(tab)

c            WRITE(SLOUT,'(A,::)') ' *'
          ENDIF
c            WRITE(SLOUT,*) ' '
        ENDDO

        IF (ikprb.GT.NULL) THEN
          dip_ik(ir,FSP1) = ikprb
          dip_te(ir,FSP1) = ktebs(ikprb,ir) +
     .                      tprb * (ktebs(ikprb+1,ir) - ktebs(ikprb,ir))
          dip_ti(ir,FSP1) = ktibs(ikprb,ir) +
     .                      tprb * (ktibs(ikprb+1,ir) - ktibs(ikprb,ir))
          dip_ne(ir,FSP1) = knbs(ikprb,ir) +
     .                      tprb * (knbs(ikprb+1,ir) - knbs(ikprb,ir))
          dip_v (ir,FSP1) = kvhs(ikprb,ir) +
     .                      tprb * (kvhs(ikprb+1,ir) - kvhs(ikprb,ir))
          dip_s (ir,FSP1) = kss(ikprb,ir) +
     .                      tprb * (kss(ikprb+1,ir) - kss(ikprb,ir))
        ELSE
          CALL WN('ProbePath','Intersection not found')
c          CALL ER('ProbePath','Intersection not found',*99)
        ENDIF


        IF (stopopt3.NE.3) dip_v(ir,FSP1) = 0.0


      ENDDO


      IF (stopopt3.NE.3) WRITE(0,*) 'SETTING DIP_V TO ZERO'


      WRITE(SLOUT,*)
      WRITE(SLOUT,*) 'DIVIMP data at FSP location'
      WRITE(SLOUT,*)
      DO ir = irsep, irwall-1
        WRITE(SLOUT,'(A,I4,F8.3,2X,2F8.3,1P,2E10.2)')
     .    'IR dip_s,te,ti,ne,v  = ',
     .    ir,dip_s (ir,FSP1),dip_te(ir,FSP1),
     .       dip_ti(ir,FSP1),dip_ne(ir,FSP1),
     .       dip_v (ir,FSP1)
      ENDDO

      RETURN
99    WRITE(EROUT,'(5X,A,I4)') 'IR = ',ir
      WRITE(EROUT,'(5X,A,1P,4E15.7)') 'A1,A2 B1,B2 = ',a1,a2,b1,b2
      STOP
      END
c
c ======================================================================
c
c subroutine: ReadProbeData
c
      SUBROUTINE ReadProbeData(status)

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      LOGICAL status

      REAL GetNe,GetJsat

      CHARACTER buffer*200
      INTEGER   fp,i1,i2,i3,ir,npts
      REAL             dum (5),jsat
      DOUBLE PRECISION dum2(5)

      DOUBLE PRECISION DLO,DHI
c...  BUG: Low range should have been negative. SL June 28, 2000
c      PARAMETER (DLO=1.0E-37,DHI=1.0E+37)
      PARAMETER (DLO=-1.0E+37,DHI=1.0E+37)

      INTEGER CalcPoint
      INTEGER result,idum1,id
      REAL*8  r1,z1,r2,z2,r3,z3,r4,z4,t,dist,mindist

      status = .FALSE.

      IF (.TRUE..OR.cmodopt.EQ.1) THEN

        IF (cmodopt.NE.1) WRITE(0,*) 'USING THESIS PROBE DATA '//
     .                               'LOADING CODE INSTEAD OF MAIN '//
     .                               'VERSION CODE'

        fp = 91

        REWIND(fp)

        WRITE(SLOUT,*)
        WRITE(SLOUT,'(A)') 'FMP AND FSP DATA'

        READ(fp,'(A200)',END=98,ERR=98) buffer
        BACKSPACE fp

        CALL TransferLine(fp,SLOUT,buffer,2)
        READ(buffer,*) shot,slice_time

        DO i1 = 1, NUMPRB
          CALL TransferLine(fp,SLOUT,buffer,2)
          READ(buffer,*) npts
          CALL TransferLine(fp,SLOUT,buffer,1)

          prb_num(i1) = 0
          DO i2 = 1, npts
            CALL TransferLine(fp,SLOUT,buffer,1)
            READ(buffer,*) (dum2(i3),i3=1,5)

            DO i3 = 1, 5
              dum(i3) = REAL(DMAX1(DLO,DMIN1(DHI,dum2(i3))))
            ENDDO

            IF (dum(2).LT.HI.AND.dum(3).LT.HI) THEN
              i3          = prb_num(i1) + 1
              prb_num(i1) = i3
              prb_rho(i3,i1) = dum(1)
              prb_te (i3,i1) = dum(2)
              prb_ti (i3,i1) = dum(2)
              prb_ne (i3,i1) = dum(3)
              prb_r  (i3,i1) = dum(4)
              prb_z  (i3,i1) = dum(5)
            ENDIF
          ENDDO


          WRITE(SLOUT,*)
          DO i2 = 1, prb_num(i1)
            WRITE(SLOUT,'(A,I2,A,I4,1P,3E15.7)')
     .        'PROBE ',i1,': I RHO,TE,NE = ',
     .        i2,prb_rho(i2,i1),prb_te(i2,i1),prb_ne(i2,i1)
          ENDDO
        ENDDO


        IF (prb_shift.NE.99.0) CALL AlignProbeData

c...temp!
        i1 = FSP1


        IF (tarsource.GE.5.AND.tarsource.LE.8) THEN
c...      Mapping targets to probe data using R,Z data:
	 
c...      Outer target:
          DO i1 = 1, prb_num(OFMP)
            mindist = HI
            DO ir = irsep, nrs
              IF (idring(ir).EQ.-1) CYCLE

              id = korpg(nks(ir),ir)
              r1 = rvertp(3,id)
              z1 = zvertp(3,id)
              r2 = rvertp(4,id)
              z2 = zvertp(4,id)
              r3 = prb_r(i1,OFMP)
              z3 = prb_z(i1,OFMP)
	 
              result = CalcPoint(r1,z1,r2,z2,r3,z3,t)
	 
              WRITE(SLOUT,'(A,2I6,3(2F8.4,2X),F10.2,I6)') 'OUT FMP: ',
     .          ir,i1,r1,z1,r2,z2,r3,z3,t,result

              IF (result.EQ.1.OR.result.EQ.2) THEN
                r4 = r1 + t * (r2 - r1)
                z4 = z1 + t * (z2 - z1)
                dist = DSQRT((r3 - r4)**2.0D0 + (z3 - z4)**2.0D0)
                IF (dist.LT.mindist) THEN
c...BUG: Sep 25, 2000 - SL
                  prb_rho(i1,OFMP) =  rho(ir,OUT23) + t * 
     .                               (rho(ir,IN14) - rho(ir,OUT23))
c                  prb_rho(i1,IFMP) =  rho(ir,OUT23) + t * 
c     .                               (rho(ir,IN14) - rho(ir,OUT23))
                  mindist = dist
                  WRITE(SLOUT,*) '*',t
                ENDIF
              ENDIF
            ENDDO

            IF (mindist.EQ.HI)
     .        CALL ER('ReadProbeData','Outer probe not mapped to '//
     .                                'grid',*99)
            
          ENDDO

          i1 = OFMP
          WRITE(SLOUT,*)
          DO i2 = 1, prb_num(i1)
            WRITE(SLOUT,'(A,I2,A,I4,1P,3E15.7)')
     .        'PROBE RHO ',i1,': I RHO,TE,NE = ',
     .        i2,prb_rho(i2,i1),prb_te(i2,i1),prb_ne(i2,i1)
          ENDDO


c...      Inner target:
          DO i1 = 1, prb_num(IFMP)
            mindist = HI
            DO ir = irsep, nrs
              IF (idring(ir).EQ.-1) CYCLE

              id = korpg(1,ir)
              r1 = rvertp(2,id)
              z1 = zvertp(2,id)
              r2 = rvertp(1,id)
              z2 = zvertp(1,id)
              r3 = prb_r(i1,IFMP)
              z3 = prb_z(i1,IFMP)
	 
              result = CalcPoint(r1,z1,r2,z2,r3,z3,t)
	 
              WRITE(SLOUT,'(A,2I6,3(2F8.4,2X),F10.2,I6)') 'IN  FMP: ',
     .          ir,i1,r1,z1,r2,z2,r3,z3,t,result

              IF (result.EQ.1.OR.result.EQ.2) THEN
                r4 = r1 + t * (r2 - r1)
                z4 = z1 + t * (z2 - z1)
                dist = DSQRT((r3 - r4)**2.0D0 + (z3 - z4)**2.0D0)
                IF (dist.LT.mindist) THEN
                  prb_rho(i1,IFMP) =  rho(ir,OUT23) + t * 
     .                               (rho(ir,IN14) - rho(ir,OUT23))
                  mindist = dist
                  WRITE(SLOUT,*) '*',t
                ENDIF
              ENDIF
            ENDDO

            IF (mindist.EQ.HI)
     .        CALL ER('ReadProbeData','Inner probe not mapped to '//
     .                                'grid',*99)

          ENDDO

          i1 = IFMP
          WRITE(SLOUT,*)
          DO i2 = 1, prb_num(i1)
            WRITE(SLOUT,'(A,I2,A,I4,1P,3E15.7)')
     .        'PROBE RHO ',i1,': I RHO,TE,NE = ',
     .        i2,prb_rho(i2,i1),prb_te(i2,i1),prb_ne(i2,i1)
          ENDDO

        ENDIF

      ELSE

c
c       Gets probe data from .13 experimental data file:
c

      ENDIF

      status = .TRUE.

      RETURN
98    WRITE(0,*) 'PROBE DATA NOT FOUND'
      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: AlignProbeData
c
      SUBROUTINE AlignProbeData

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER   ii,nfmp,nfsp
      REAL      n,r,r1,r2,rnum,rstep,rshift,shift,sum,minsum,vfmp,vfsp,
     .          bestshift,weight,
     .          dum1(MAXPRB),dum2(MAXPRB),dum3(MAXPRB),dum4(MAXPRB)
      CHARACTER type*6

c      prb_align = 1

      WRITE(SLOUT,*)
      WRITE(SLOUT,'(A)') 'FSP AND FMP PRESSURE ALIGNMENT'
      WRITE(SLOUT,*)

      type   = 'LINEAR'
      rshift =  0.0020
      rstep  =  0.0001
      rnum   =  50.0
      minsum =  HI
      bestshift = 0.0

c...set this to use OSM_PROBE...
      nfmp = prb_num(OFMP)
      nfsp = prb_num(FSP1)

      r1 = MAX(prb_rho(1   ,OFMP),prb_rho(1   ,FSP1)+rshift)
      r2 = MIN(prb_rho(nfmp,OFMP),prb_rho(nfsp,FSP1)-rshift)

      r1 = MAX(r1,0.20*r2)
c
c
c
      DO ii = 1, MAXPRB
        dum1(ii) = prb_rho(ii,OFMP)

        dum3(ii) = 2.0 * 0.5 * prb_ne(ii,OFMP) * 2.0 * prb_te(ii,OFMP)
        dum4(ii) =             prb_ne(ii,FSP1) * 2.0 * prb_te(ii,FSP1)
      ENDDO

      DO shift = -rshift, rshift, rstep
        sum = 0.0

        DO ii = 1, MAXPRB
          dum2(ii) = prb_rho(ii,FSP1) + shift
        ENDDO

c        WRITE(SLOUT,'(A,3F12.6)') 'SHIFTS = ',-rshift,shift,rshift

        DO r = r1, r2, (r2 - r1) / rnum
          CALL Fitter(nfmp,dum1,dum3,1,r,vfmp,type)
          CALL Fitter(nfsp,dum2,dum4,1,r,vfsp,type)

          weight = vfmp / 1.0E+21

          sum = sum + ABS(vfmp - vfsp) * weight

c          WRITE(SLOUT,'(3X,A,1P,3E11.3,2X,2E11.3,2X,E11.3)')
c     .      'Stat R1,R,R2 VFMP,VFST SUM = ',
c     .      r1,r,r2,vfmp,vfsp,sum
         ENDDO

         WRITE(SLOUT,'(A,1P,2E15.7)') 'SHIFT SUM = ',shift,sum

         IF (sum.LT.minsum) THEN
           minsum    = sum
           bestshift = shift
         ENDIF
      ENDDO


      IF (prb_shift.NE.0.0) THEN
        WRITE(0,*) 'FORCING PROBE SHIFT'
        WRITE(SLOUT,'(2(A,F10.5,3X))') 'SUGGESTED SHIFT = ',bestshift,
     .                                 'FORCED SHIFT    = ',prb_shift
      ELSE
        WRITE(SLOUT,'(A,F10.5)') 'SUGGESTED SHIFT = ',bestshift
      ENDIF

      WRITE(SLOUT,*) 'PRB_ALIGN= ',prb_align

      IF (prb_align.EQ.1) THEN
        IF (prb_shift.NE.99.0) THEN
          WRITE(SLOUT,'(A)') 'FORCED SHIFT PERFORMED'
          DO ii = 1, prb_num(FSP1)
            prb_rho(ii,FSP1) = prb_rho(ii,FSP1) + prb_shift
          ENDDO
        ELSE
          WRITE(SLOUT,'(A)') 'CALCULATED SHIFT PERFORMED'
          DO ii = 1, prb_num(FSP1)
            prb_rho(ii,FSP1) = prb_rho(ii,FSP1) + bestshift
          ENDDO
          prb_shift = bestshift
        ENDIF
      ENDIF

      IF (tarsource.EQ.0) tarsource = 1

      RETURN
99    STOP
      END

c
c ======================================================================
c
c
c ======================================================================
c
c function:FindPeak
c
c ...will crash on uniform QUANT values (or two equal values anyway)!
c
      REAL FUNCTION FindPeak(ik,ir,quant)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      INTEGER ik,ir
      REAL    quant(MAXNKS,MAXNRS)
c...note: If the PEAKSPAN is not representative of the number of cells on a
c         given ring, then one noisy broad peak could be identified as two
c         peaks.
      INTEGER    PEAKSPAN
c
c...THIS IS TERRIBLY GRID DEPENDENT!
c
      PARAMETER (PEAKSPAN = 8)
c      PARAMETER (PEAKSPAN = 10)      

      INTEGER SymmetryPoint

      INTEGER iks,ike,ik1,ikp1,ikp2,ikp3,ikm
      REAL    pk,pk2

      FindPeak = 1.0

      pk  = -1.0
      pk2 = -1.0

      ikm = SymmetryPoint(ir)

      IF (ik.LE.ikm) THEN
        iks  = ikbound(ir,IKLO)
        ikp1 = iks
        ikp2 = iks
        ikp3 = iks

        DO ik1 = iks, ik
          IF (ABS(quant(ik1,ir)).GT.ABS(quant(ikp2,ir))) ikp2 = ik1
        ENDDO

        ikp1 = ikp2
        DO WHILE (ABS(quant(ikp1,ir)).GT.0.10*ABS(quant(ikp2,ir)).AND.
     .            ikp1.GT.iks)
          ikp1 = ikp1 - 1
        ENDDO

        ikp3 = ikp2

c        WRITE(0,*) '-->',ik,IR,IKP1,ikp2,ikp3
        DO WHILE (ABS(quant(ikp3,ir)).GT.0.10*ABS(quant(ikp2,ir)).AND.
     .            ikp3.LT.ikm)
          ikp3 = ikp3 + 1
        ENDDO

        CALL CalcIntegral4(quant,ikp1,ikp3,ir,pk,2)
        
c        IF (ik.EQ.43.OR.ik.EQ.45)
c     .    WRITE(0,'(3X,A,4I6,1P,E12.4,0P)') 
c     .      'FIND PEAK  IKP1-3,IR PK  = ',ikp1,ikp2,ikp3,ir,pk
      ELSE

        iks  = ikbound(ir,IKHI)
        ikp1 = iks
        ikp2 = iks
        ikp3 = iks

        DO ik1 = iks,  ik, -1
          IF (ABS(quant(ik1,ir)).GT.ABS(quant(ikp2,ir))) ikp2 = ik1
        ENDDO

        pk = ABS(quant(ikp2,ir))

c        WRITE(0,*) 'FIND PEAK  IKP3  ,IR PK = ',ikp2,ir,pk
      ENDIF

      DO ik1 = MAX(ikbound(ir,IKLO),ik-PEAKSPAN),
     .         MIN(ikbound(ir,IKHI),ik+PEAKSPAN)
c        WRITE(0,*) '-->',ik1,ik
        IF (ABS(quant(ik1,ir)).GT.ABS(quant(ik,ir))) FindPeak = 0.0

c        IF (ik.EQ.43.OR.ik.EQ.45)
c     .    WRITE(0,'(3X,A,2(I6,1PE12.4,0P),F4.1)')
c     .      ' PEAK CRAP ',ik1,quant(ik1,ir),ik,quant(ik,ir),FindPeak

      ENDDO
c
c     Calculate integral of local peak:
c
      IF (FindPeak.EQ.1.0.AND.ik.NE.ikp2) THEN

        IF (ik.LE.ikm) THEN
          iks  = ikbound(ir,IKLO)
          ikp1 = ik
          ikp2 = ik
          ikp3 = ik

          DO WHILE (ABS(quant(ikp1,ir)).GT.0.10*ABS(quant(ik,ir)).AND.
     .              ikp1.GT.iks)
            ikp1 = ikp1 - 1
          ENDDO

          DO WHILE (ABS(quant(ikp3,ir)).GT.0.10*ABS(quant(ik,ir)).AND.
     .              ikp3.LT.ikm)
            ikp3 = ikp3 + 1
          ENDDO

          ikp1 = MAX(ikp1,ik-PEAKSPAN)
          ikp3 = MIN(ikp3,ik+PEAKSPAN)

          CALL CalcIntegral4(quant,ikp1,ikp3,ir,pk2,2)

          IF (ik.EQ.43.OR.ik.EQ.45)
     .      WRITE(0,'(3X,A,4I6,1P,E12.4,0P)')  
     .      '  TEST     IKP1-3,IR PK2 = ',
     .        ikp1,ikp2,ikp3,ir,pk2

          IF (pk2.LT.0.10*pk) THEN
            FindPeak = 0.0
          ELSE
            FindPeak = pk2 / pk
          ENDIF

        ELSE

          IF (ABS(quant(ik,ir)).LT.0.10*pk) FindPeak = 0.0

        ENDIF

      ENDIF

c          IF (ik.EQ.43.OR.ik.EQ.45)
c     .      WRITE(0,*) '   FINDPEAK= ',FindPeak


      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c
      REAL FUNCTION GetMach(v,te,ti)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'

      REAL te,ti,v

      REAL GetCs

      GetMach = ABS(v / GetCs(te,ti))

      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c
      REAL FUNCTION GetNe(te,ti,jsat,machno)
      IMPLICIT none
      INCLUDE 'params'
      REAL te,ti,jsat,machno

      REAL GetCs

      GetNe = ABS(jsat / (ECH * machno * GetCs(te,ti)))

      RETURN
99    STOP
      END
c
c
c ======================================================================
c
c
c
      REAL FUNCTION GetCs(te,ti)
      IMPLICIT none
      INCLUDE 'params'
      INCLUDE 'comtor'
      REAL, INTENT(IN) :: te,ti

c     Te,i in eV
c     a    in amu
c     result is m s-1

      GetCs = 9.78817E+03 * SQRT(0.5 * (1.0 + rizb) * (te + ti) / crmb)

      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c
      REAL FUNCTION GetJsat(te,ti,ne,v)
      IMPLICIT none

      REAL te,ti,ne,v,vb

      REAL GetCs

      IF (v.EQ.1.0) THEN
        vb = GetCs(te,ti)
      ELSE
        vb = v
      ENDIF

      GetJsat = ne * 1.602E-19 * vb

c      WRITE(0,*) 'GETJSAT:',ne,ECH,vb,getjsat,te,ti

      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c
      REAL FUNCTION GetFlux(region,ring)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      REAL GetCs

      INTEGER region,ring
      INTEGER id,ir
      REAL    brat
c      REAL    brat,mn

      ir = ring

      IF (region.EQ.IKLO) THEN
        id = idds(ir,2)
        brat = 1.0 / kbfs(1,ir)
c        brat = 1.0 / kbfst(ir,2)
c        brat = bratio (1,ir)
c        mn   = cmachno(ir,2)
      ELSEIF (region.EQ.IKHI) THEN
        id = idds(ir,1)
        brat = 1.0 / kbfs(nks(ir),ir)
c        brat = 1.0 / kbfst(ir,1)
c        brat = bratio (nks(ir),ir)
c        mn   = cmachno(ir,1)
      ELSE
        CALL ER('GetFlux','Invalid region',*99)
      ENDIF

      IF (ir.LT.irsep.OR.idring(ir).EQ.BOUNDARY) THEN 
        GetFlux = 0.0
      ELSE
c...fix? dds to dds2?  kbfst, above
        GetFlux = knds(id) * kvds(id) * dds2(id) * brat *
     .            2.0 * PI * rp(id) * costet(id) * eirsrcmul * 
     .            eirtorfrac

c         write(88,'(A,1P,2I6,9E10.2,0P)')  'calculating',
c     .     ikds(id),irds(id),knds(id),kvds(id),dds2(id),brat,
c     .     rp(id),costet(id),eirsrcmul,
c     .     eirtorfrac,GetFlux

      ENDIF

c      GetFlux = knds(id) * GetCs(kteds(id),ktids(id),) * 
c                mn * dds(id) *
c     .          brat * 2.0 * PI * rp(id) * costet(id)


      IF (supflx(region,ir).EQ.1) GetFlux = GetFlux * 1.0E-15


      RETURN
99    WRITE(0,*) '  RING  =',ir
      WRITE(0,*) '  REGION=',region
      STOP
      END
c
c ======================================================================
c
c
c
      REAL FUNCTION GetGamma(region,ring)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER region,ring
      INTEGER id,ir
      REAL    brat,isat,t_ratio,delta_e,m_ratio

      ir = ring

      IF (region.EQ.IKLO) THEN
        id = idds(ir,2)
      ELSEIF (region.EQ.IKHI) THEN
        id = idds(ir,1)
      ELSE
        CALL ER('GetGamma','Invalid region',*99)
      ENDIF

      IF (ir.LT.irsep.OR.idring(ir).EQ.BOUNDARY) THEN 
        GetGamma = 0.0
      ELSE
c...    Taken from Stangeby 1st edition, equation 25.46, pg 649.  Note that many effects
c       are missing, i.e. realistic secondary electron emission (0.0 here), e-i recombination
c       energy, atom-atom recombination energy, low collisionality effects, space charge
c       effects, etc. see the discussion by Stangeby pp 646-654. -SL, 29.03.2010
        t_ratio = ktids(id) / kteds(id)
        delta_e = 0.0
        m_ratio = 9.11E-31 / (crmb * 1.67E-27)
        GetGamma = 2.5 * t_ratio + 2.0 / (1.0 - delta_e) - 
     .             0.5 * LOG( (2.0 * PI * m_ratio) * (1 + t_ratio) * 
     .                        (1 - delta_e)**-2 )
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c
      REAL FUNCTION GetHeatFlux(region,ring)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      REAL GetFlux,GetGamma

      INTEGER region,ring
      INTEGER id,ir
      REAL    brat,isat,gamma

      ir = ring

      IF (region.EQ.IKLO) THEN
        id = idds(ir,2)
      ELSEIF (region.EQ.IKHI) THEN
        id = idds(ir,1)
      ELSE
        CALL ER('GetHeatFlux','Invalid region',*99)
      ENDIF

      IF (ir.LT.irsep.OR.idring(ir).EQ.BOUNDARY) THEN 
        GetHeatFlux = 0.0
      ELSE
        isat  = GetFlux(region ,ring)
        gamma = GetGamma(region,ring) 
        GetHeatFlux = gamma * ABS(isat) * ECH * kteds(id)
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c
      REAL FUNCTION GetIonRelTime(ni,ti,mi,Zi)
      IMPLICIT none

      REAL ni,ti,mi,Zi

      REAL lambda,coeff,time

      lambda = 30.0 - 0.5 * LOG(ni) + 1.5 * LOG(ti)
  
      coeff = 2.502D+26

      GetIonRelTime = (coeff * SQRT(mi) * ti**1.5) / 
     .                (ni * Zi**4.0 * lambda)


      RETURN
99    STOP
      END
c
c ======================================================================
c
c ...work in progress... 
c
      REAL FUNCTION IonViscosity(ik,ir,ni,ti,mi,Zi)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'

      REAL ni,ti,mi,Zi

      REAL GetIonRelTime,GetMach

      INTEGER ik,ir
      REAL    ionreltime,ionpressure,velgradient,frac1,frac2,vb1,vb2

      LOGICAL, SAVE :: message = .FALSE.

c      REAL       ECH
c      PARAMETER (ECH = 1.602E-19)

c...  
      ionpressure = ni * ti * ECH

c...
      ionreltime = GetIonRelTime(ni,ti,mi,Zi)

c...
      IF (ir.LT.irsep) THEN
c...    
        ionviscosity = 0.0       
        IF (.NOT.message) THEN
          WRITE(0,*)
          WRITE(0,*) 'VISCOSITY NOT READY IN CORE'
          WRITE(0,*)
          message = .TRUE.
        ENDIF
        RETURN

        
c        STOP 'VISCOSITY NOT READY IN CORE'
      ELSE
        IF (ik.EQ.1) THEN
          vb1 = 1.0
          vb2 = 1.0

        ELSEIF (ik.EQ.nks(ir)) THEN
          vb1 = 1.0
          vb2 = 1.0

        ELSE
          frac1 = (ksb(ik-1,ir) - kss(ik-1,ir)) / 
     .            (kss(ik  ,ir) - kss(ik-1,ir))

          frac2 = (ksb(ik  ,ir) - kss(ik  ,ir)) / 
     .            (kss(ik+1,ir) - kss(ik  ,ir))

          vb1 = (1.0 - frac1) * kvhs(ik-1,ir) + frac1 * kvhs(ik  ,ir)
          vb2 = (1.0 - frac2) * kvhs(ik  ,ir) + frac2 * kvhs(ik+1,ir)



        ENDIF
      ENDIF

      velgradient = (vb2 - vb1) / (ksb(ik,ir) - ksb(ik-1,ir))

c      WRITE(0,'(A,2I6,1P,3E10.2,4X,3E10.2,0P,4X,F10.2,4X,F10.2)') 
c     .  'VIS:',ik,ir,vb1,kvhs(ik,ir),vb2,
c     .         velgradient,ionreltime,ionpressure,
c     .         -ionpressure*ionreltime*velgradient,
c     .         GetMach(kvhs(ik,ir),ktebs(ik,ir),ktibs(ik,ir))

      ionviscosity = 0.0

      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c
      REAL FUNCTION GetIonSrc(region,ir,stratum,mode)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      REAL GetFlux

      INTEGER region,ir,stratum,mode,ik,iks,ike
      REAL    sumion,sumrec,valion,valrec,flxtar


c...  Add a check on STRATUM:


      IF     (region.EQ.IKLO) THEN
        iks = 1
        ike = nks(ir) / 2
        flxtar = ABS(GetFlux(IKLO,ir)) / eirtorfrac
      ELSEIF (region.EQ.IKHI) THEN 
        iks = nks(ir) / 2 + 1
        ike = nks(ir)
        flxtar = ABS(GetFlux(IKHI,ir)) / eirtorfrac
      ELSE
        CALL ER('GetIonSrc','Invalid REGION',*99)
      ENDIF

      sumion = 0.0
      sumrec = 0.0
      DO ik = iks, ike
        IF (stratum.EQ.0) THEN
          valion = pinion(ik,ir)
        ELSE
          valion = pindata(ik,ir,H_ION1+stratum-1)
        ENDIF

        valrec = pinrec(ik,ir)

c...    Integrate.  Get values for the entire torus by dividing
c       by EIRTORFRAC:
        sumion = sumion + kvols(ik,ir) * valion / eirtorfrac
        sumrec = sumrec + kvols(ik,ir) * valrec / eirtorfrac  
      ENDDO

      IF     (mode.EQ.0) THEN
c...    No normalization:       
      ELSEIF (mode.EQ.1) THEN
c...    Normalize ionisation source by total ion sink:
        sumion = sumion / (flxtar + sumrec)
      ELSE
        CALL ER('GetIonSrc','Invalid MODE',*99)
      ENDIF

      GetIonSrc = sumion



      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: CalcIntegral4
c
c Calculate integral based on the absolute value of array elements.
c
c
      SUBROUTINE CalcIntegral4(quant,ik1,ik2,ir,integral,mode)
      IMPLICIT none

      INCLUDE 'params'

c     Input:
      INTEGER ik1,ik2,ir,mode
      REAL    quant(MAXNKS,MAXNRS),integral

      INTEGER          ik
      REAL             quantabs(MAXNKS,MAXNRS)
      DOUBLE PRECISION dpdum1

      dpdum1 = 0.0

      CALL RZero(quantabs,MAXNKS*MAXNRS)

      DO ik = ik1, ik2
        quantabs(ik,ir) = ABS(quant(ik,ir))
      ENDDO

      CALL CalcIntegral2(quantabs(1,ir),ik1,ik2,ir,dpdum1,mode)

      integral = SNGL(dpdum1)

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: CalcIntegral3
c
      SUBROUTINE CalcIntegral3(quant,ik1,ik2,ir,integral,mode)
      IMPLICIT none

      INCLUDE 'params'

c     Input:
      INTEGER ik1,ik2,ir,mode
      REAL    quant(MAXNKS,MAXNRS),integral

      DOUBLE PRECISION dpdum1

      dpdum1 = 0.0

      CALL CalcIntegral2(quant(1,ir),ik1,ik2,ir,dpdum1,mode)

c      WRITE(0,*) '????:',dpdum1

      integral = SNGL(dpdum1)

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: CalcIntegral2
c
      SUBROUTINE CalcIntegral2(quant,ik1,ik2,ir,integral,mode)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

c     Input:
      REAL    quant(*)
      REAL*8  integral
      INTEGER mode,ir,ik1,ik2

c     Local:
      REAL*8  val
      INTEGER ik,iks,ike

      integral = 0.0D0

      IF (mode.EQ.5) THEN
        iks = ik1
        ike = ik2
      ELSE
        iks = ik1
        ike = ik2 
      ENDIF

      DO ik = iks, ike

        val = 0.0D0

        IF (mode.EQ.1) THEN

          IF     (ik.EQ.1        ) THEN
            val = DBLE(quant(1)) * DBLE(kss(1,ir) - ksb(0,ir))
          ELSEIF (ik.EQ.nks(ir)+1) THEN
            val = DBLE(quant(nks(ir))) * DBLE(ksmaxs(ir) - 
     .                                        kss(nks(ir),ir))
          ELSE
            val = DBLE(quant(ik-1)) * DBLE(ksb(ik-1,ir) - kss(ik-1,ir))+
     .            DBLE(quant(ik  )) * DBLE(kss(ik  ,ir) - ksb(ik-1,ir))
          ENDIF
        ELSEIF (mode.EQ.2) THEN
          IF     (iks.EQ.1.AND.iks.EQ.ike) THEN
            val = DBLE(quant(ik)) * DBLE(kss(ik,ir) - ksb(ik-1,ir))
          ELSEIF (iks.EQ.nks(ir).AND.iks.EQ.ike) THEN
            val = DBLE(quant(ik)) * DBLE(ksb(ik,ir) - kss(ik,ir))
          ELSEIF (ik.EQ.iks.AND.iks.GT.1      ) THEN
            val = DBLE(quant(ik)) * DBLE(ksb(ik,ir) - kss(ik,ir))
          ELSEIF (ik.EQ.ike.AND.ike.LT.nks(ir)) THEN
            val = DBLE(quant(ik)) * DBLE(kss(ik,ir) - ksb(ik-1,ir))
          ELSE
            val = DBLE(quant(ik)) * DBLE(ksb(ik,ir) - ksb(ik-1,ir))
          ENDIF
        ELSEIF (mode.EQ.3) THEN
          IF     (ik.EQ.iks) THEN
            val = DBLE(quant(ik)) * DBLE(ksb(ik,ir) - kss(ik,ir))
          ELSEIF (ik.EQ.ike) THEN
            val = DBLE(quant(ik)) * DBLE(kss(ik,ir) - ksb(ik-1,ir))
          ELSE
            val = DBLE(quant(ik)) * DBLE(ksb(ik,ir) - ksb(ik-1,ir))
          ENDIF
        ELSEIF (mode.EQ.4) THEN
          val = DBLE(quant(ik)) * DBLE(ksb(ik,ir) - ksb(ik-1,ir))
        ELSEIF (mode.EQ.5) THEN

          val = DBLE(quant(ik+1)) * DBLE(kss(ik+1,ir) - kss(ik,ir))

          IF (ik.EQ.ike.AND.ik.EQ.nks(ir))
     .      val = val + DBLE(quant(ik+1)) * DBLE(ksb(ik,ir) - 
     .                                           kss(ik,ir))

          IF (ik.EQ.iks.AND.ik.EQ.1) 
     .      val = val + DBLE(quant(1)) * DBLE(kss(1,ir) - ksb(0,ir))

        ELSEIF (mode.EQ.6) THEN
          IF (ik.EQ.iks.AND.iks.GT.1      ) THEN
            val = DBLE(quant(ik)) * DBLE(ksb(ik,ir) - kss(ik,ir))
          ELSE
            val = DBLE(quant(ik)) * DBLE(ksb(ik,ir) - ksb(ik-1,ir))
          ENDIF
        ELSEIF (mode.EQ.7) THEN
          IF (ik.EQ.ike.AND.ike.LT.nks(ir)) THEN
            val = DBLE(quant(ik)) * DBLE(kss(ik,ir) - ksb(ik-1,ir))
          ELSE
            val = DBLE(quant(ik)) * DBLE(ksb(ik,ir) - ksb(ik-1,ir))
          ENDIF
        ELSEIF (mode.EQ.8) THEN
c...      Integration from start of ring to center of cell IK,IR:
          IF (ik.EQ.ike) THEN
            val = DBLE(quant(ik)) * DBLE(kss(ik,ir) - ksb(ik-1,ir))
          ELSE
            val = DBLE(quant(ik)) * DBLE(ksb(ik,ir) - ksb(ik-1,ir))
          ENDIF
        ELSEIF (mode.EQ.10) THEN
c...      Integration from center of cell IK,IR to end of ring:
          IF (ik.EQ.iks) THEN
            val = DBLE(quant(ik)) * DBLE(ksb(ik,ir) - kss(ik  ,ir))
          ELSE
            val = DBLE(quant(ik)) * DBLE(ksb(ik,ir) - ksb(ik-1,ir))
          ENDIF
        ELSEIF (mode.EQ.9) THEN
c...      Cell centre to centre:
          IF (ik.GE.1.AND.ik.LE.nks(ir)) THEN
            IF     (ik.EQ.iks) THEN
              val = DBLE(quant(ik)) * DBLE(ksb(ik,ir) - kss(ik,ir))
c              WRITE(0,*) 'QHANT 1:',ik,iks,ike,val
            ELSEIF (ik.EQ.ike.AND.ik.LE.nks(ir)) THEN
              val = DBLE(quant(ik)) * DBLE(kss(ik,ir) - ksb(ik-1,ir))
c              WRITE(0,*) 'QHANT 2:',ik,iks,ike,val
            ELSE
              val = DBLE(quant(ik)) * DBLE(ksb(ik,ir) - ksb(ik-1,ir))
c              WRITE(0,*) 'QHANT  :',ik,iks,ike,val
            ENDIF
          ENDIF
        ELSEIF (mode.EQ.99) THEN
          val = DBLE(quant(ik)) * DBLE(kss(ik+1,ir) - kss(ik,ir))
        ELSE
          CALL ER('CalcIntergral2','Invalid MODE specified',*99)
        ENDIF

c...DEV:
c        WRITE(88,*) 'CHECK: ',ik,ir,val

        integral = integral + val
c        WRITE(0,*) 'CHECK: ',integral,mode,ik
      ENDDO



      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: CalcIntegral1
c
      SUBROUTINE CalcIntegral1(quant,ik1,ik2,ir,integral,mode)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

c     Input:
      REAL*8           quant(*)
      DOUBLE PRECISION integral
      INTEGER mode,ir,ik1,ik2

c     Local:
      REAL*8  val,dksb(0:MAXNKS),dkss(MAXNKS),dksmaxs
      INTEGER ik,iks,ike

      integral = 0.0D0

      dksb(0) = DBLE(ksb(0,ir))
      DO ik = 1, nks(ir)
        dksb(ik) = DBLE(ksb(ik,ir))
        dkss(ik) = DBLE(kss(ik,ir))
      ENDDO
      dksmaxs = DBLE(ksmaxs(ir))

      IF (mode.EQ.5) THEN
        iks = ik1
        ike = ik2
      ELSE
        iks = ik1
        ike = ik2 
      ENDIF

      DO ik = iks, ike
        IF (mode.EQ.1) THEN
          IF     (ik.EQ.1        ) THEN
            val = quant(1) * (dkss(1) - dksb(0))
          ELSEIF (ik.EQ.nks(ir)+1) THEN
            val = quant(nks(ir)) * (dksmaxs - dkss(nks(ir)))
          ELSE
            val = quant(ik-1) * (dksb(ik-1) - dkss(ik-1)) +
     .            quant(ik  ) * (dkss(ik  ) - dksb(ik-1))
          ENDIF
        ELSEIF (mode.EQ.2) THEN
          IF     (ik.EQ.iks.AND.iks.GT.1      ) THEN
            val = quant(ik) * (dksb(ik) - dkss(ik))
          ELSEIF (ik.EQ.ike.AND.ike.LT.nks(ir)) THEN
            val = quant(ik) * (dkss(ik) - dksb(ik-1))
          ELSE
            val = quant(ik) * (dksb(ik) - dksb(ik-1))
          ENDIF
        ELSEIF (mode.EQ.3) THEN
          IF     (ik.EQ.iks) THEN
            val = quant(ik) * (dksb(ik) - dkss(ik))
          ELSEIF (ik.EQ.ike) THEN
            val = quant(ik) * (dkss(ik) - dksb(ik-1))
          ELSE		      
            val = quant(ik) * (dksb(ik) - dksb(ik-1))
          ENDIF
        ELSEIF (mode.EQ.4) THEN
          val = quant(ik) * (dksb(ik) - dksb(ik-1))
        ELSEIF (mode.EQ.5) THEN

          val = quant(ik+1) * (dkss(ik+1) - dkss(ik))

          IF (ik.EQ.ike.AND.ik.EQ.nks(ir))
     .      val = val + quant(ik+1) * (dksb(ik) - dkss(ik))

          IF (ik.EQ.iks.AND.ik.EQ.1) 
     .      val = val + quant(1) * (dkss(1) - dksb(0))

        ELSEIF (mode.EQ.6) THEN
          IF (ik.EQ.iks.AND.iks.GT.1) THEN
            val = quant(ik) * (dksb(ik) - dkss(ik))
          ELSE
            val = quant(ik) * (dksb(ik) - dksb(ik-1))
          ENDIF
        ELSEIF (mode.EQ.7) THEN
          IF (ik.EQ.ike.AND.ike.LT.nks(ir)) THEN
            val = quant(ik) * (dkss(ik) - dksb(ik-1))
          ELSE
            val = quant(ik) * (dksb(ik) - dksb(ik-1))
          ENDIF
        ELSEIF (mode.EQ.99) THEN
          val = quant(ik) * (dkss(ik+1) - dkss(ik))
        ENDIF

c...DEV:
c        WRITE(88,*) 'CHECK: ',ik,ir,val

        integral = integral + val
      ENDDO

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: CalcIntegral
c
      REAL FUNCTION CalcIntegral(source,mode,ir)

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      REAL    integral,source(MAXNKS,MAXNRS),val
      INTEGER mode,ir

      INTEGER ik,iks,ike

      integral = 0.0

      IF     (mode.EQ.1.OR.mode.EQ.4) THEN
        iks = 1
        ike = ikmids(ir)
      ELSEIF (mode.EQ.2.OR.mode.EQ.5) THEN
        iks = ikmids(ir) + 1
        ike = nks(ir)
      ELSE
        iks = 1
        ike = nks(ir)
      ENDIF

      DO ik = iks, ike
        IF (mode.LE.3) THEN
          val = karea2(ik,ir)
        ELSE
          val = ksb(ik,ir) - ksb(ik-1,ir)
        ENDIF

        integral = integral + source(ik,ir) * val
      ENDDO

      CalcIntegral = integral

      RETURN
99    STOP
      END











c
c ======================================================================
c
      SUBROUTINE UpdateLine1I(fp1,fp2,buffer,loc,int)

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

c     Input:
      INTEGER   fp1,fp2,loc,int
      CHARACTER buffer*(*)

      CALL ReadLine (fp1,buffer,1,*97,*98)
      WRITE(buffer((loc-1)*6+1:loc*6),'(I6)') int
      CALL WriteLine(fp2,buffer)

      RETURN
c
c     Error code:
c
97    CONTINUE
      CALL ER('UpdateLine1I','Unexpected end of file',*99)
98    CONTINUE
      CALL ER('UpdateLine1I','Problems reading template file',*99)
99    CONTINUE
      WRITE(EROUT,*) '  Last line read: '
      WRITE(EROUT,*) '  "',buffer,'"'
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE UpdateLine2I(fp1,fp2,buffer,loc1,loc2,int1,int2)

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

c     Input:
      INTEGER   fp1,fp2,loc1,loc2,int1,int2
      CHARACTER buffer*(*)

      CALL ReadLine (fp1,buffer,1,*97,*98)
      WRITE(buffer((loc1-1)*6+1:loc1*6),'(I6)') int1
      WRITE(buffer((loc2-1)*6+1:loc2*6),'(I6)') int2
      CALL WriteLine(fp2,buffer)

      RETURN
c
c     Error code:
c
97    CONTINUE
      CALL ER('UpdateLine2I','Unexpected end of file',*99)
98    CONTINUE
      CALL ER('UpdateLine2I','Problems reading template file',*99)
99    CONTINUE
      WRITE(EROUT,*) '  Last line read: '
      WRITE(EROUT,*) '  "',buffer,'"'
      STOP
      END
c


c
c ======================================================================
c
      SUBROUTINE TransferLine(fp1,fp2,buffer,numline)

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

c     Input:
      INTEGER   fp1,fp2,numline
      CHARACTER buffer*(*)

      INTEGER   ii

      DO ii = 1, numline
10      CALL ReadLine (fp1,buffer,1,*97,*98)
        IF (buffer(1:2).EQ.'* ') GOTO 10
        CALL WriteLine(fp2,buffer)
      ENDDO

      RETURN
c
c     Error code:
c
97    CONTINUE
      CALL ER('TransferLine','Unexpected end of file',*99)
98    CONTINUE
      CALL ER('TransferLine','Problems reading file',*99)
99    CONTINUE
      WRITE(EROUT,*) '  Last line read: '
      WRITE(EROUT,*) '  "',buffer,'"'
      STOP
      END
c
c ======================================================================
c
c Same as above but doesn't skip lines starting with '*':
c
      SUBROUTINE TransferLine2(fp1,fp2,buffer,numline)

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

c     Input:
      INTEGER   fp1,fp2,numline
      CHARACTER buffer*(*)

      INTEGER   ii

      DO ii = 1, numline
10      CALL ReadLine (fp1,buffer,1,*97,*98)
        CALL WriteLine(fp2,buffer)
      ENDDO

      RETURN
97    CONTINUE
      CALL ER('TransferLine','Unexpected end of file',*99)
98    CONTINUE
      CALL ER('TransferLine','Problems reading file',*99)
99    CONTINUE
      WRITE(EROUT,*) '  Last line read: '
      WRITE(EROUT,*) '  "',buffer,'"'
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE ReadLine(fp,buffer,numline,*,*)

      IMPLICIT none

c     Input:
      INTEGER   fp,numline
      CHARACTER buffer*(*)

      INTEGER ii

      DO ii = 1, numline
        WRITE(buffer,'(200X)')
        READ(fp,'(A200)',END=98,ERR=99) buffer
      ENDDO

      RETURN
98    RETURN 1
99    RETURN 2
      END
c
c ======================================================================
c
      SUBROUTINE WriteLine(fp,buffer)

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

c     Input:
      INTEGER   fp
      CHARACTER buffer*(*)

      INTEGER i1,i2

      i1 = 0

      WRITE(fp,'(A)',ERR=98) buffer(1:LEN_TRIM(buffer))

      RETURN
c
c     Error code:
c
98    CONTINUE
      WRITE(EROUT,*) 'ERROR (WriteLine): Problem writing to file'
      STOP
      END
c














c
c
c
      LOGICAL FUNCTION SameString(string1,string2)
      IMPLICIT         none

      CHARACTER string1*(*),string2*(*)
      INTEGER   len1,len2,i1

      SameString = .TRUE.

      len1 = LEN_TRIM(string1)
      len2 = LEN_TRIM(string2)

      IF (len1.NE.len2) THEN
        SameString = .FALSE.
      ELSE
        DO i1 = 1, len1
          IF (string1(i1:i1).NE.string2(i1:i1)) SameString = .FALSE.
        ENDDO
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c
c
      SUBROUTINE DB(message)
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'slcom'

      CHARACTER message*(*)

      IF (outmode.GE.3) THEN
        WRITE(0    ,'(A)') 'DEBUG: '//message
        WRITE(SLOUT,'(A)') 'DEBUG: '//message
      ENDIF

      RETURN
      END
c
c ======================================================================
c
c
c
c
      SUBROUTINE MS(routine,message)

      IMPLICIT none

      CHARACTER routine*(*),message*(*)

      INCLUDE 'params'
      INCLUDE 'slcom'

      LOGICAL SameString

      INTEGER MAXMESS
      PARAMETER (MAXMESS=128)

      CHARACTER*128 store(MAXMESS),cdum1
      INTEGER       count(MAXMESS),i1,ncount,fp,nmess
      DATA ncount,nmess /0, 0/
      SAVE

      IF (routine.EQ.'dump') THEN
        READ(message,*) fp

        WRITE(fp,'(I6,A)')
     .    nmess,' message(s) posted'

        IF (outmode.EQ.2) THEN
          DO i1 = 1, ncount
            cdum1 = store(i1)

            WRITE(fp,'(5X,I4,1X,A)') count(i1),cdum1(1:LEN_TRIM(cdum1))
          ENDDO
        ENDIF
      ELSE

        IF (outmode.GE.3)
     .    WRITE(0    ,'(4A)') ' MESSAGE ',routine,': ',message

        WRITE(SLOUT,'(4A)') ' MESSAGE ',routine,': ',message

        DO i1 = 1, ncount
          IF (SameString(message,store(i1))) THEN
            count(i1) = count(i1) + 1
            GOTO 10
          ENDIF
        ENDDO

10      IF (i1.GT.ncount.AND.ncount.LT.MAXMESS) THEN
          WRITE(store(i1),'(128X)')

          store(i1) = message
          count(i1) = 1
          ncount    = ncount + 1
        ENDIF

        nmess = nmess + 1

      ENDIF

      RETURN
      END
c
c ======================================================================
c
c
c
c
      SUBROUTINE WN(routine,message)
      IMPLICIT   none

      CHARACTER routine*(*),message*(*)

      INCLUDE 'params'
      INCLUDE 'slcom'

      LOGICAL SameString

      INTEGER MAXMESS
      PARAMETER (MAXMESS=128)

      CHARACTER*128 store(MAXMESS),cdum1
      INTEGER       count(MAXMESS),i1,ncount,fp,nwarn
      DATA ncount,nwarn /0, 0/
      SAVE

      IF (routine.EQ.'dump') THEN
        READ(message,*) fp

        WRITE(fp,'(I6,A)')
     .    nwarn,' warning(s) posted'

        IF (outmode.EQ.0) THEN
          DO i1 = 1, ncount
            cdum1 = store(i1)

            WRITE(fp,'(5X,I4,1X,A)') count(i1),cdum1(1:LEN_TRIM(cdum1))
          ENDDO
        ENDIF

      ELSE

        IF (outmode.GE.1) THEN
          WRITE(0     ,'(4A)') ' WARNING ',routine,': ',message
          WRITE(PINOUT,'(4A)') ' WARNING ',routine,': ',message
          WRITE(6     ,'(4A)') ' WARNING ',routine,': ',message
        ENDIF

        WRITE(SLOUT,'(4A)') 'WARNING ',routine,': ',message

        DO i1 = 1, ncount
          IF (SameString(message,store(i1))) THEN
            count(i1) = count(i1) + 1
            GOTO 10
          ENDIF
        ENDDO

10      IF (i1.GT.ncount.AND.ncount.LT.MAXMESS) THEN
          WRITE(store(i1),'(128X)')

          store(i1) = message
          count(i1) = 1
          ncount    = ncount + 1
        ENDIF

        nwarn = nwarn + 1

      ENDIF

      RETURN
      END

c
c ======================================================================
c
c
c
c
      SUBROUTINE HD(fp,message,code,start,end)
      IMPLICIT   none

      INTEGER   fp,start,end
      CHARACTER message*(*),code*(*)

      INCLUDE 'params'
      INCLUDE 'slcom'

      INTEGER fpin

      fpin = 97

      IF (fp.EQ.fpin) 
     .  CALL ER('HD','Input and output streams are the same',*99)

      WRITE(fp,*)
      WRITE(fp,*) message(1:LEN_TRIM(message))//':'

      OPEN(UNIT=fpin,FILE='info.dat',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)

      CALL GetInfo(fpin,fp,code,start,end)

      CLOSE(fpin)

      RETURN
98    CALL WN('HD','Unable to access info.dat data file')
      WRITE(fp,*) '  ERROR: UNABLE TO OPEN INFO.DAT FILE'
      RETURN
99    STOP
      END



c
c ======================================================================
c
c subroutine: GetInfo
c
      SUBROUTINE GetInfo(fpin,fpout,tag,indent,maxlen)

      IMPLICIT none

      INTEGER   fpin,fpout,indent,maxlen
      CHARACTER tag*(*),line*1024

      INCLUDE 'params'
      INCLUDE 'slcom'

      CHARACTER comment*80,buffer*2048,word*50,sp*40
      INTEGER   lentag,lenbuf,i1,i2,i3,linepos

      WRITE(sp  ,'(  40X)')
      WRITE(line,'(1024X)')
      lentag = LEN_TRIM(tag)

      REWIND(fpin)

10    CALL ReadLine(fpin,buffer,1,*96,*98)
      IF (buffer(1:lentag+1).NE.'['//tag(1:lentag)) GOTO 10

20    READ (fpin,'(A2048)',END=97,ERR=98) buffer
      IF (buffer(1:1).NE.'[') THEN
        IF (buffer(1:1).EQ.'$') GOTO 20

        lenbuf = LEN_TRIM(buffer)

        buffer(lenbuf+1:lenbuf+1) = ' '

        lenbuf  = lenbuf + 1
        i3      = 0
        linepos = indent

        WRITE(line(1:1024),'(A)') sp(1:indent)

        DO i2 = 1, lenbuf
          IF (buffer(i2:i2).EQ.' ') THEN
            IF (linepos+i3+1.GT.maxlen) THEN
              IF (i3.LT.maxlen-indent) THEN
                WRITE(fpout,'(A)') line(1:LEN_TRIM(line))
                WRITE(line,'(1024X)')
                WRITE(line,'(A)') sp(1:indent)
                linepos = indent
              ELSE
                CALL ER('GetInfo','Word too long',*99)
              ENDIF
            ENDIF
            WRITE(line(linepos+1:1024),'(A)') word(1:i3)//' '
            linepos = linepos + i3 + 1
            i3 = 0
          ELSE
            i3 = i3 + 1
            word(i3:i3) = buffer(i2:i2)
          ENDIF
        ENDDO

c        WRITE(fpout,*)
        WRITE(fpout,'(A)') line(1:LEN_TRIM(line))
        WRITE(line  ,'(1024X)')
        GOTO 20
      ENDIF

      RETURN
96    WRITE(fpout,'(A)')
     .  sp(1:indent)//'[SORRY, TAG "'//tag(1:lentag)//'" NOT FOUND]'
      CALL MS('GetInfo','Tag "'//tag(1:lentag)//'" not found')
97    RETURN
98    CALL ER('GetInfo','Problem reading data file',*99)
99    WRITE(EROUT,'(5X,A,2I4)') 'FPIN,FPOUT = ',fpin,fpout
      WRITE(EROUT,'(5X,A)')     'TAG        = "'//tag(1:lentag)//'"'
      END








c
c
c ======================================================================
c
c
c
      SUBROUTINE ReadIR(line,ival,rval,imin,imax,tag)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      CHARACTER line*72,tag*(*)
      INTEGER fp,ival,imin,imax
      REAL    rval

      INTEGER i
      REAL    r
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,i,r

      IF (i.LT.imin.OR.i.GT.imax)
     .  CALL ER('ReadI','Out of bounds: '//line,*99)

      ival = i
      rval = r

      WRITE(SLOUT,'(A)') line
      WRITE(SLOUT,'(5X,2A,I4,1P,E10.2)') tag,' = ',ival,rval

      RETURN
98    WRITE(EROUT,*) 'Problem reading unstructured input'
99    WRITE(EROUT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(EROUT,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(EROUT,'(5X,A,3I4)') 'I,IVAL,IMIN,IMAX = ',i,ival,imin,imax
      STOP
      END
c
c
c ======================================================================
c
c
c
      SUBROUTINE ReadI(line,ival,imin,imax,tag)

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      CHARACTER line*72,tag*(*)
      INTEGER fp,ival,imin,imax

      INTEGER i
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,i

      IF (i.LT.imin.OR.i.GT.imax) then 

        write (0,*)  'READI:ERROR:',i,imin,imax 
        CALL ER('ReadI','Out of bounds: '//line,*99)

      endif

      ival = i

      WRITE(SLOUT,'(A)')        line
      WRITE(SLOUT,'(5X,2A,I4)') tag,' = ',ival

      RETURN
98    WRITE(EROUT,*) 'Problem reading unstructured input'
99    WRITE(EROUT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(EROUT,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(EROUT,'(5X,A,3I4)') 'I,IVAL,IMIN,IMAX = ',i,ival,imin,imax
      STOP
      END
c
c
c ======================================================================
c
c
c
      SUBROUTINE ReadC(line,cval,tag)

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      CHARACTER line*(*),tag*(*),cval*(*)
      INTEGER fp,ival,imin,imax

      INTEGER i
      CHARACTER comment*72

      integer erout1

      erout1 = 0

      READ (line,*,ERR=98,END=98) comment,cval

      WRITE(SLOUT,'(A)')        line
      WRITE(SLOUT,'(5X,2A,A)') tag,' = ',cval

      RETURN
98    WRITE(EROUT1,*) 'Problem reading unstructured input'
99    WRITE(EROUT1,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(EROUT1,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(EROUT1,'(5X,2A)')    'CVAL = ''',cval,''''
      STOP 'READC'
      END
c
c ======================================================================
c
c
c
      SUBROUTINE Read2I(line,ival1,ival2,imin,imax,tag)

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      CHARACTER line*72,tag*(*)
      INTEGER fp,ival1,ival2,imin,imax

      INTEGER i1,i2
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,i1,i2

      IF (i1.LT.imin.OR.i1.GT.imax.OR.
     .    i2.LT.imin.OR.i2.GT.imax)
     .  CALL ER('Read2I','Out of bounds: '//line,*99)

      ival1 = i1
      ival2 = i2

      WRITE(SLOUT,'(A)')        line
      WRITE(SLOUT,'(5X,2A,I4)') tag,' = ',ival1
      WRITE(SLOUT,'(5X,2A,I4)') tag,' = ',ival2

      RETURN
98    WRITE(EROUT,*) 'Problem reading unstructured input'
99    WRITE(EROUT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(EROUT,'(5X,2A)')    'TAG  = ''',tag,''''
      STOP
      END
c
c ======================================================================
c
c
c
      SUBROUTINE ReadR(line,rval,rmin,rmax,tag)

      IMPLICIT none

      CHARACTER line*72,tag*(*)
      REAL rval,rmin,rmax

      INCLUDE 'params'
      INCLUDE 'slcom'

      REAL r
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,r

      IF (r.LT.rmin.OR.r.GT.rmax)
     .  CALL ER('ReadR','Out of bounds: '//line,*99)

      rval = r

      WRITE(SLOUT,'(A)')        line
      WRITE(SLOUT,'(2A,G10.3)') tag,' = ',rval

      RETURN
98    WRITE(EROUT,*) 'Problem reading unstructured input'
99    WRITE(EROUT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(EROUT,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(EROUT,'(5X,A,3G10.3)')
     .  'R,RVAL,RMIN,RMAX = ',r,rval,rmin,rmax
      STOP
      END
c
c
c ======================================================================
c
c
c
      SUBROUTINE Read2R(line,rval1,rval2,rmin,rmax,tag)

      IMPLICIT none

      CHARACTER line*72,tag*(*)
      REAL rval1,rval2,rmin,rmax

      INCLUDE 'params'
      INCLUDE 'slcom'

      REAL r1,r2
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,r1,r2

      IF (r1.LT.rmin.OR.r1.GT.rmax.OR.
     .    r2.LT.rmin.OR.r2.GT.rmax)
     .  CALL ER('ReadR','Out of bounds: '//line,*99)

      rval1 = r1
      rval2 = r2

      WRITE(SLOUT,'(A)')        line
      WRITE(SLOUT,'(2A,2G10.3)') tag,' = ',rval1,rval2

      RETURN
98    WRITE(EROUT,*) 'Problem reading unstructured input'
99    WRITE(EROUT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(EROUT,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(EROUT,'(5X,A,6G10.3)')
     .  'R,RVAL,RMIN,RMAX = ',r1,r2,rval1,rval2,rmin,rmax
      STOP
      END














c
c ======================================================================
c
c subroutine: ISet
c
c
      SUBROUTINE ISet(array,ndim,val)
      implicit none   

      INTEGER ndim
      INTEGER array(ndim),val

      INTEGER ii

      DO ii = 1, ndim
        array(ii) = val
      ENDDO

      RETURN
      END

c
c ======================================================================
c
c subroutine: RSet
c
c
      SUBROUTINE RSet(array,ndim,val)
      implicit none 
      INTEGER ndim
      REAL array(ndim),val

      INTEGER ii

      DO ii = 1, ndim
        array(ii) = val
      ENDDO

      RETURN
      END











c
c ======================================================================
c
c subroutine:
c
c
c
c ... CHECK OPTIONS!
c
      INTEGER FUNCTION SymmetryPoint(ir)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'


      INTEGER ir,ik,ik1,ik2


      IF     (osm_symopt.EQ.0.OR.
     .        (osm_symopt.EQ.1.AND.
     .         ir.GT.irwall.AND.ir.LE.nrs)) THEN

        SymmetryPoint = ikmids(ir)

      ELSEIF (osm_symopt.EQ.1) THEN
        IF (osm_probe.EQ.1) THEN
          IF (prb_num(osm_probe).GT.0) THEN
            SymmetryPoint = dip_ik(ir,osm_probe)
          ELSE
            CALL ER('Interface','Probe data not found',*99)
          ENDIF
        ELSE
          CALL ER('Interface','Unsupported reference probe',*99)
        ENDIF

      ELSEIF (osm_symopt.EQ.2) THEN
 

c... BUG:THIS OPTION FAILE WITHOUT THIS!
        SymmetryPoint = ikmids(ir)

        DO ik = nks(ir), 1, -1
          IF (ksb(ik-1,ir)/ksmaxs(ir).GT.0.667) THEN
           SymmetryPoint = ik
c          ELSE
c            EXIT
          ENDIF
        ENDDO

c        WRITE(0,*) 'MARK: SYMPT',ir,nks(ir),ik,SymmetryPoint

c        IF (symmetrypoint.EQ.2079) CALL OUTPUTDATA(85,'crap')
c        IF (symmetrypoint.EQ.2079) STOP 'crap'

      ELSEIF (osm_symopt.EQ.3) THEN
 
        SymmetryPoint = ikmids(ir)

        ik1 = ikbound(ir,IKLO)
        ik2 = ikbound(ir,IKHI)

        DO ik = ik2, ik1, -1
          IF (ksb(ik-1,ir).GT.(0.5*kss(ik1,ir)+0.5*kss(ik2,ir))) THEN
            SymmetryPoint = ik
          ENDIF
        ENDDO

      ELSEIF (osm_symopt.EQ.4) THEN
 
        SymmetryPoint = ikmids(ir)

        ik1 = ikfluid(IKLO,ir)
        ik2 = ikfluid(IKHI,ir)

        DO ik = ik2, ik1, -1
          IF (ksb(ik-1,ir).GT.(0.5*kss(ik1,ir)+0.5*kss(ik2,ir))) THEN
            SymmetryPoint = ik
          ENDIF
        ENDDO

      ELSEIF (osm_symopt.EQ.5) THEN
 
        SymmetryPoint = ikmids(ir)

        ik1 = ikfluid(IKLO,ir)
        ik2 = ikfluid(IKHI,ir)

        DO ik = ik2, ik1, -1
          IF (kpb(ik-1,ir).GT.(0.5*kps(ik1,ir)+0.5*kps(ik2,ir))) THEN
            SymmetryPoint = ik
          ENDIF
        ENDDO

      ELSE
        CALL ER('Interface','Unsupported symmetry point option',*99)
      ENDIF



      RETURN
99    STOP 'SymmetryPoint'
      END


c
c ======================================================================
c
c function: GetModel
c
      INTEGER FUNCTION GetModel(region,ir)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER region,ir

      INTEGER in,reg,ir1,ir2

      GetModel = -1

      IF (orgcioptg.EQ.90.OR.orgcioptg.EQ.91.OR.
     .    orgcioptg.EQ.92) THEN
c...    GetModel will be assigned to the first instance of the half-ring 
c       in the list:
        DO in = 1, nbgplas
          ir1 = bgplasopt(in,1)
          ir2 = bgplasopt(in,2)
          reg = bgplasopt(in,3)

          IF (ir.GE.ir1.AND.ir.LE.ir2.AND.
     .        (reg.EQ.region.OR.reg.EQ.3).AND.
     .        GetModel.EQ.-1)
     .      GetModel = INT(bgplasopt(in,5))
        ENDDO
      ENDIF

c...  No SOL model found, apply default:
      IF (GetModel.EQ.-1) GetModel = orgcioptf

      IF (GetModel.EQ.22.AND.osm_preopt.GT.0) THEN
c
c       If SOL22 prescription mode is turned on, then check if
c       SOL21 data is specified for the ring.  If so, indicate that
c       the prescription mode by returning model number 24:
c
        IF (region.EQ.IKLO) THEN
          DO in = 1, s21_ndatao
            IF (ir.EQ.INT(s21_datao(in,1))) GetModel = 24
          ENDDO
        ELSE
          DO in = 1, s21_ndatai
            IF (ir.EQ.INT(s21_datai(in,1))) GetModel = 24
          ENDDO
        ENDIF
      ENDIF

c...  SOL27:
      IF (GetModel.EQ.22.AND.osm_preopt.LT.0) THEN
c...    If SOL27 source prescription mode is turned on, then check if
c       the prescription data is specified for the ring.  If so, 
c       return model number 27:
        IF (region.EQ.IKLO) THEN
          DO in = 1, s21_ndatao
            IF (ir.EQ.INT(s21_datao(in,1))) GetModel = 27
          ENDDO
        ELSE
          DO in = 1, s21_ndatai
            IF (ir.EQ.INT(s21_datai(in,1))) GetModel = 27
          ENDDO
        ENDIF
      ENDIF

c      WRITE(0,*) 'GETMODEL B:',region,ir,getmodel

      RETURN
99    WRITE(EROUT,'(A,2I4)') 'IR REGION = ',ir,region
      STOP 
      END





c
c subroutine: VolInteg
c
c
      SUBROUTINE VolInteg(quant,region,ir1,ir2,integ)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER region,ir1,ir2
      REAL    quant(MAXNKS,MAXNRS),integ

      INTEGER ik,ik1,ik2,ir

      integ = 0.0

c REGION:
c
c 1 - IKLO
c 2 - IKHI
c 3 - ENTIRE RING
c 4 - INNER HALF-RING, BELOW X-POINT
c 5 - INNER HALF-RING, ABOVE X-POINT
c 6 - OUTER HALF-RING, BELOW X-POINT
c 7 - OUTER HALF-RING, ABOVE X-POINT
c

      DO ir = ir1, ir2
        IF (idring(ir).EQ.BOUNDARY) CYCLE

        IF     (region.EQ.IKLO.OR.region.EQ.4.OR.
     .                            region.EQ.5) THEN
          ik1 = 1
          ik2 = ikmids(ir)
        ELSEIF (region.EQ.IKHI.OR.region.EQ.6.OR.
     .                            region.EQ.7) THEN
          ik1 = ikmids(ir) + 1
          ik2 = nks(ir)
        ELSEIF (region.EQ.3) THEN
          ik1 = 1
          ik2 = nks(ir)
        ELSE
          CALL ER('VolInteg','Invlid REGION',*99)
        ENDIF

        DO ik = ik1, ik2
          IF ((region.EQ.4.AND.zs(ik,ir).GT.zxp).OR.
     .        (region.EQ.5.AND.zs(ik,ir).LT.zxp).OR.
     .        (region.EQ.6.AND.zs(ik,ir).GT.zxp).OR.
     .        (region.EQ.7.AND.zs(ik,ir).LT.zxp)) CYCLE

          integ = integ + quant(ik,ir) * kvols(ik,ir)
        ENDDO
      ENDDO


      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: ShiftTargetData
c
c Target data (LPDATO and LPDATI) is shifted in "along the target"
c space.  A positive shift moves the data outward (generally, in the
c direction of increasing major radius).
c
c
      SUBROUTINE ShiftTargetData
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INCLUDE 'solparams'
      INCLUDE 'solswitch'

      INTEGER RingNo,GetJsat,getne

      INTEGER fp,in,ir,ir1,ir2,id1,id2,i1,idat,ierr,ik1
      REAL    sddat(MAXINS),tedat(MAXINS),tidat(MAXINS),nedat(MAXINS),
     .        sd,sgn(MAXNRS),te,ti,ne,cs,rdum1

      fp = PINOUT

      CALL RSet(sgn,MAXNRS,1.0)
c
c     jdemod - sepdist2 is now adjusted to be negative in PFZ so the
c              sgn array should no longer be needed
c
c      DO ir = irtrap, nrs
c        sgn(ir) = -1.0
c      ENDDO
c
      CALL RZero(sddat,MAXINS)
      CALL RZero(tedat,MAXINS)
      CALL RZero(tidat,MAXINS)
      CALL RZero(nedat,MAXINS)
c
c...  Output:
      CALL HD(fp,'Shifting target data','TARSHIFT-HD',5,67)
      WRITE(fp,'(4X,2F10.4)') tarshift(IKLO),tarshift(IKHI)
      WRITE(fp,*) '  Before:'
      ir = irtrap + 1
      DO WHILE (ir.NE.irwall)
        id1 = RingNo(ir,lpdato,nlpdato,MAXINS,4,ierr)
        id2 = RingNo(ir,lpdati,nlpdati,MAXINS,4,ierr)
        WRITE(fp,'(3X,I4,2(2F8.2,1P,E10.2,0P,2X))')
     .    ir,lpdato(id1,2),lpdato(id1,3),lpdato(id1,4),
     .       lpdati(id2,2),lpdati(id2,3),lpdati(id2,4)
        ir = irouts(1,ir)
      ENDDO

c
c...  Adjust low index target data (inner target on JET; outer target on
c     C-Mod, DIII-D and ASDEX-U):
      IF (tarshift(IKLO).NE.0.0) THEN
        ir1 = irtrap + 1
        DO WHILE (ir1.NE.irwall)
          idat = 0
          ir2  = ir1
          ik1  = 1
          id1  = korpg(ik1,ir2)
          id2  = korpg(ik1,irouts(ik1,ir2))
          DO WHILE (rvertp(2,id1).EQ.rvertp(1,id2).AND.
     .              zvertp(2,id1).EQ.zvertp(1,id2).AND.ir2.NE.irwall)
            in = RingNo(ir2,lpdato,nlpdato,MAXINS,4,ierr)
            IF (ierr.EQ.0) THEN
              idat = idat + 1
              sddat(idat) = sepdist2(idds(ir2,2)) * sgn(ir2)
              tedat(idat) = lpdato(in,2)
              tidat(idat) = lpdato(in,3)
              nedat(idat) = lpdato(in,4)

c...          Patch that forces the data listed in unstructured input block
c             079 when modifying LODATO:
              IF (switch(SWIONP).EQ.-6.0) THEN
                DO i1 = 1, osmnppv
                  IF ( osmppv(i1,1).EQ.REAL(ir2).AND.
     .                (osmppv(i1,2).EQ.1.0.OR.osmppv(i1,2).EQ.3.0)) THEN
                    WRITE(0,*) 'USING 079 DATA FOR INNER TARGET ',ir2
                    tedat(idat) = osmppv(i1,3)
                    tidat(idat) = osmppv(i1,4)                    
                    IF (lpdatsw.EQ.1) THEN
c...                  Any attempt to call GETJSAT, which I would have preferred since
c                     keeping the Isat calculation in one place helps to
c                     avoid inconsistencies, crashes the case under RedHat Linux 6.0 and
c                     PGI compiler version 3.0.  At the moment it looks to be
c                     a compiler bug. (SL, Aug 14, 2000)
                      te = osmppv(i1,3)
                      ti = osmppv(i1,4)
                      ne = osmppv(i1,5)
                      cs = 9.78817E+03 * 
     .                     SQRT(0.5 * (1.0 + rizb) * (te + ti) / crmb)
                      nedat(idat) =  ne * ECH * cs
                    ELSE
                      nedat(idat) = osmppv(i1,5)
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF

              ir2 = irouts(ik1,ir2)
              ik1 = 1
              id1 = korpg(ik1,ir2)
              id2 = MAX(korpg(ik1,irouts(1,ir2)),1)
            ELSE
              CALL ER('ShiftTargetData','Low index data not found',*99)
            ENDIF
          ENDDO
          IF (idat.GT.0) THEN
            in = RingNo(ir2,lpdato,nlpdato,MAXINS,4,ierr)
            IF (ierr.EQ.0) THEN
              idat = idat + 1
              sddat(idat) = sepdist2(idds(ir2,2)) * sgn(ir2)
              tedat(idat) = lpdato(in,2)
              tidat(idat) = lpdato(in,3)
              nedat(idat) = lpdato(in,4)
            ELSE
              CALL ER('ShiftTargetData','Low index data not found',*99)
            ENDIF
c
c...        Shift target data by interpolating target data stored in local
c           xxdat arrays:
            ir = ir1
            DO i1 = 1, idat
              in = RingNo(ir,lpdato,nlpdato,MAXINS,4,ierr)
              sd = sepdist2(idds(ir,2)) * sgn(ir) + tarshift(IKLO)

              IF ( tarshiftopt.EQ.0.OR.
     .            (tarshiftopt.EQ.1.AND.ir.LT.irwall)) THEN
                CALL Fitter(idat,sddat,tedat,1,sd,lpdato(in,2),'LINEAR')
                CALL Fitter(idat,sddat,tidat,1,sd,lpdato(in,3),'LINEAR')
                CALL Fitter(idat,sddat,nedat,1,sd,lpdato(in,4),'LINEAR')
              ENDIF

              ik1 = 1
              ir  = irouts(ik1,ir)
            ENDDO
          ENDIF
          ir1 = MIN(ir2+1,irwall)
        ENDDO
      ENDIF

c
c...  Adjust high index target data:
      IF (tarshift(IKHI).NE.0.0) THEN
        ir1 = irtrap + 1
        DO WHILE (ir1.NE.irwall)
          idat = 0
          ir2  = ir1
          ik1  = nks(ir2)
          id1  = korpg(ik1,ir2)
          id2  = korpg(nks(irouts(ik1,ir2)),irouts(ik1,ir2))
          DO WHILE (rvertp(3,id1).EQ.rvertp(4,id2).AND.
     .              zvertp(3,id1).EQ.zvertp(4,id2).AND.ir2.NE.irwall)
            in = RingNo(ir2,lpdati,nlpdati,MAXINS,4,ierr)
            IF (ierr.EQ.0) THEN
              idat = idat + 1
              sddat(idat) = sepdist2(idds(ir2,1)) * sgn(ir2)
              tedat(idat) = lpdati(in,2)
              tidat(idat) = lpdati(in,3)
              nedat(idat) = lpdati(in,4)

c...          Patch that forces the data listed in unstructured input block
c             079 when modifying LODATI:
              IF (switch(SWIONP).EQ.-6.0) THEN
                DO i1 = 1, osmnppv
                  IF ( osmppv(i1,1).EQ.REAL(ir2).AND.
     .                (osmppv(i1,2).EQ.2.0.OR.osmppv(i1,2).EQ.3.0)) THEN
                    WRITE(0,*) 'USING 079 DATA FOR OUTER TARGET ',ir2
                    tedat(idat) = osmppv(i1,3)
                    tidat(idat) = osmppv(i1,4)                    
                    IF (lpdatsw.EQ.1) THEN
                      te = osmppv(i1,3)
                      ti = osmppv(i1,4)
                      ne = osmppv(i1,5)
                      cs = 9.78817E+03 * 
     .                     SQRT(0.5 * (1.0 + rizb) * (te + ti) / crmb)
                      nedat(idat) =  ne * ECH * cs
                    ELSE
                      nedat(idat) = osmppv(i1,5)
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF

              ir2 = irouts(ik1,ir2)
              ik1 = nks(ir2)
              id1 = korpg(ik1,ir2)
              id2 = MAX(korpg(nks(irouts(ik1,ir2)),irouts(ik1,ir2)),1)
            ELSE
              CALL ER('ShiftTargetData','High index data not found',*99)
            ENDIF
          ENDDO
          IF (idat.GT.0) THEN
            in = RingNo(ir2,lpdati,nlpdati,MAXINS,4,ierr)
            IF (ierr.EQ.0) THEN
              idat = idat + 1
              sddat(idat) = sepdist2(idds(ir2,1)) * sgn(ir2)
              tedat(idat) = lpdati(in,2)
              tidat(idat) = lpdati(in,3)
              nedat(idat) = lpdati(in,4)
            ELSE
              CALL ER('ShiftTargetData','High index data not found',*99)
            ENDIF
c
c...        Shift target data:
            ir = ir1
            DO i1 = 1, idat
              in = RingNo(ir,lpdati,nlpdati,MAXINS,4,ierr)
              sd = sepdist2(idds(ir,1)) * sgn(ir) - tarshift(IKHI)

              IF ( tarshiftopt.EQ.0.OR.
     .            (tarshiftopt.EQ.1.AND.ir.LT.irwall)) THEN
                CALL Fitter(idat,sddat,tedat,1,sd,lpdati(in,2),'LINEAR')
                CALL Fitter(idat,sddat,tidat,1,sd,lpdati(in,3),'LINEAR')
                CALL Fitter(idat,sddat,nedat,1,sd,lpdati(in,4),'LINEAR')
              ENDIF

              ik1 = nks(ir)
              ir  = irouts(ik1,ir)
            ENDDO
          ENDIF
          ir1 = MIN(ir2+1,irwall)
        ENDDO
      ENDIF

c
c...  Output:
      WRITE(fp,*) '  After:'
      ir = irtrap + 1
      DO WHILE (ir.NE.irwall)
        id1 = RingNo(ir,lpdato,nlpdato,MAXINS,4,ierr)
        id2 = RingNo(ir,lpdati,nlpdati,MAXINS,4,ierr)
        WRITE(fp,'(3X,I4,2(2F8.2,1P,E10.2,0P,2X))')
     .    ir,lpdato(id1,2),lpdato(id1,3),lpdato(id1,4),
     .       lpdati(id2,2),lpdati(id2,3),lpdati(id2,4)
        ir = irouts(1,ir)
      ENDDO

      RETURN
99    STOP
      END

c
c ======================================================================
c
c subroutine: OutputMomentumData
c
      SUBROUTINE OutputMomentumData
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER   fp,ik,ir,i1,i2
      REAL      sum
      CHARACTER note*20

      fp = 6

c      OPEN(UNIT=fp,FILE='rec.dat',ACCESS='SEQUENTIAL',STATUS='REPLACE')

      WRITE(fp,*) 'Momentum loss data:'

      WRITE(fp,*) 'NOT IN USE AT THE MOMENT'
      RETURN

      DO ir = osm_watch1, osm_watch2
c      DO ir = irsep, nrs
        IF (ir.EQ.0.OR.idring(ir).EQ.-1) CYCLE

        WRITE(fp,'(2A4,14A9,2A12)') 'ik','ir','1','2','3','4','5','6',
     .                              '7','8','9','10','11','12','13',
     .                              '14','sum','pinmp'

        DO ik = 1, nks(ir)

          note(1:20) = '                   '
          IF (ik.EQ.ikto2 (ir)) note = note(1:LEN_TRIM(note))//' IKTO2'
          IF (ik.EQ.ikti2 (ir)) note = note(1:LEN_TRIM(note))//' IKTI2'
          IF (ik.EQ.ikmids(ir)) note = note(1:LEN_TRIM(note))//' IKMIDS'
          IF (ik.EQ.ikbound(ir,IKLO))
     .      note = note(1:LEN_TRIM(note))//' IK1'
          IF (ik.EQ.ikbound(ir,IKHI))
     .      note = note(1:LEN_TRIM(note))//' IK2'

          sum = 0.0
          DO i1 = H_MP1, MAXDATA-1
            sum = sum + pindata(ik,ir,i1)
          ENDDO

          WRITE(fp,'(2I4,1P,14E9.1,2E12.4,0P,A)')
     .      ik,ir,(pindata(ik,ir,i1),i1=H_MP1,MAXDATA),sum,pinmp(ik,ir),
     .      note(1:LEN_TRIM(note))
        ENDDO
      ENDDO

c      CLOSE(fp)

      RETURN
99    STOP
      END
c slnote - moved here from rundiv.d6a
c
      integer function calc_random_seed(n)
      implicit none
      integer n
c
c     Calculates a random number seed based on the current date and 
c     time. The seed number will be limited to n digits if n is specified
c     greater than zero.
c 
c     Local varaiables
c
      integer j,ciseed
      character*8 systim,sysdat
c
      CALL ZA08AS (SYSTIM)                                              
      CALL ZA09AS (SYSDAT)                                              
c
      CISEED = 1                                                      
      DO 300 J = 1, 8                                                 
         IF (J.EQ.3.OR.J.EQ.6) GOTO 300                                
         CISEED = CISEED * (ICHAR(SYSTIM(J:J)) - ICHAR('0') + 1) *     
     >                      (ICHAR(SYSDAT(J:J)) - ICHAR('0') + 1)       
 300  CONTINUE                                                        
      CISEED = CISEED - 1                                             
c
      if (n.gt.0) then 
c
         do while (ciseed.ge.10**n)
            ciseed = ciseed /10
         end do  
c
      endif 
c
      calc_random_seed = ciseed
c
      return
      end 
c
c
c
c ======================================================================
c
c function: atan3c
c
      REAL FUNCTION atan3c(deltaz,deltar)

      INCLUDE 'params'

c     Input:
      REAL deltaz,deltar

      REAL atan2c

      atan3c = atan2c(ABS(deltaz),ABS(deltar)) * 180.0 / PI

      IF (deltar.LT.0.0.AND.deltaz.GT.0.0) THEN
        atan3c = 180.0 - atan3c
      ELSEIF (deltar.LT.0.0.AND.deltaz.LT.0.0) THEN
        atan3c = 180.0 + atan3c
      ELSEIF (deltar.GT.0.0.AND.deltaz.LT.0.0) THEN
        atan3c = 360.0 - atan3c
      ENDIF

      RETURN
      END
c
c
c
c
c ======================================================================
c Routines to calcuation ring distributions (from grad.d5a):
c
c
c

c
c ======================================================================
c
c subroutine: GetDist
c
c
      SUBROUTINE GetDist(iks,ike,ir,oldquant,tquant,dist,mode1,
     .                           pinavail,tag)
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER SymmetryPoint,GetModel
      REAL    CalcDist

      INTEGER   iks,ike,ir,mode,mode1
      LOGICAL   pinavail
      REAL      oldquant(MAXNKS),tquant(MAXNDS),dist(*)
c      REAL      oldquant(MAXNKS),tquant(MAXNDS),dist(MAXNKS)
      CHARACTER tag*(*)

      INTEGER ik,ikm
      REAL    deltas,sum,sval

      CALL DB('Calculating cross-field distributions')

      CALL SetBounds

c      CALL RZero(dist,MAXNKS)

c      WRITE(0,*) 'DIST: ring = ',ir

      sum = 0.0
      mode = mode1

c      WRITE(0,*) tag

      IF (.NOT.pinavail.AND.
     .    ((mode.GT.9.AND.mode.LT.17).OR.mode.EQ.20)) THEN
c
c       If an old background plasma is unavailable, make sure an
c       'old background independent' method of calculating the
c       distribution is used:
c
        mode = 1

        CALL WN('GetDist','Clipping distribution mode')
      ENDIF

      ikm = SymmetryPoint(ir)

       WRITE(88,*) 'MODE:',mode

      DO ik = iks, ike

        IF     (tag.EQ.'MOCK P2') THEN
          IF (iks.EQ.1.AND.ik.EQ.iks) THEN
            deltas = kss(1,ir) - ksb(0,ir)
            sval   = 0.5 * (kss(1,ir) + ksb(0,ir))
            dist(1) = deltas * CalcDist(sval,ir,mode)
            sum = sum + dist(1)
          ENDIF

          IF (ike.EQ.nks(ir)) 
     .      STOP 'MOCK POWER NOT DONE FOR NON-SOL24 RINGS'

c...      THERE WILL BE SOME ISSUES HERE, SINCE THE INNER
c         AND OUTER RINGS WILL BE FIGHTING OVER THE MIDDLE CELL

          deltas = kss(ik+1,ir) - kss(ik,ir)
          sval   = 0.5 * (kss(ik,ir) + ksb(ik,ir))

          dist(ik+1) = deltas * CalcDist(sval,ir,mode)

c          WRITE(88,*) '--? DIST',ik,ir,sval
c          WRITE(88,*) '--? DIST',ik,iks,ike,dist(ik+1),deltas

          sum = sum + dist(ik+1)
        ELSEIF (tag.EQ.'MOCK P3') THEN
          STOP 'MOCK POWER FOR OUTER TARGET NOT DONE'
        ELSE
          deltas = ksb(ik,ir) - ksb(ik-1,ir)
          sval   = kss(ik,ir)

          dist(ik) = deltas * CalcDist(sval,ir,mode)

          IF (((ik.LT.ikbound(ir,IKLO).AND.GetModel(IKLO,ir).EQ.24).OR.
     .         (ik.GT.ikbound(ir,IKHI).AND.GetModel(IKHI,ir).EQ.24))
     .        .AND.
     .        (tag.EQ.'PARTICLE'.OR.tag.EQ.'ION POW'.OR.
     .         tag.EQ.'ELE POW')) dist(ik) = 0.0

          sum = sum + dist(ik)
        ENDIF


      ENDDO



c
c
c



c
c     Normalize distribution:
c
      DO ik = iks, ike

        IF (tag.EQ.'MOCK P2'.OR.tag.EQ.'MOCK P3') THEN

          IF (ike.EQ.nks(ir)) 
     .      STOP 'NOT READY FOR TRUE MOCK YET'

          IF (iks.EQ.1.AND.ik.EQ.iks) THEN
            deltas = kss(1,ir) - ksb(0,ir)
            dist(1) = dist(1) / sum / deltas
          ENDIF

          deltas = kss(ik+1,ir) - kss(ik,ir)
          dist(ik+1) = dist(ik+1) / sum / deltas

c          WRITE(88,*) '--? DISTB ',ik,dist(ik+1),deltas
        ELSE
          deltas = ksb(ik,ir) - ksb(ik-1,ir)
          dist(ik) = dist(ik) / sum / deltas
        ENDIF

c...DEV:
c        WRITE(88,'(A,2I6,2F16.8,E14.6)') 'DIST-> '//tag//' ',ik,ir,
c     .                  deltas,dist(ik),pinqe(ik,ir)



      ENDDO

      RETURN
99    STOP 'GetDist'
      END
c
c ======================================================================
c
c function: CalcDist
c
c Should I be working with P or S?  Does it really matter?
c
      REAL FUNCTION CalcDist(sval1,irval1,mode)
c      REAL FUNCTION CalcDist(sval1,irval1,quantv,tquantv,mode)
      IMPLICIT      none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

c     Input:
      INTEGER irval1,mode
      LOGICAL initialize
      REAL    sval1,quant(MAXNKS,MAXNRS,4),tquant(MAXNDS,4)

      INTEGER SymmetryPoint,GetModel
      REAL    StoT,TtoS,TtoP,CalcWidth,GetL1
      DOUBLE PRECISION LnLam


      INTEGER ik,iki,iko,ir,iri,iro,ini,ino,ik1,irval,midpt,qopt,
     .        ik2,ikplane(MAXNRS),id,iks
      REAL    s,thetav,val,vali,valo,gradi,grado,tolflow,
     .        grad2,sval,thetai,thetao,tgrado,tgradi,tgrad,
     .        ds1,ds2,tarflx1,tarflx2,ionsrc,recsrc,intcfp,intcfptot
      REAL    widi,wido,wid,l1,l2,s1,sdist,sfrac,s2,z1,z2
      REAL    tmini,tmino,tmaxi,tmaxo,t1i,t2i,t1o,t2o,p1i,p2i,
     .        p1o,p2o,pdum,deltapi,deltapo,deltapt,p1,p2,pmid


      COMMON /OLDPLASMA/ oldknbs ,oldktebs ,oldktibs ,oldkvhs ,
     .                   oldknbs2,oldktebs2,oldktibs2,oldkvhs2
      REAL
     .     oldknbs  (MAXNKS,MAXNRS),oldktebs (MAXNKS,MAXNRS),
     .     oldktibs (MAXNKS,MAXNRS),oldkvhs  (MAXNKS,MAXNRS),
     .     oldktebs2(MAXNKS,MAXNRS),oldktibs2(MAXNKS,MAXNRS),
     .     oldknbs2 (MAXNKS,MAXNRS),oldkvhs2 (MAXNKS,MAXNRS)


      EQUIVALENCE (quant(1,1,1),oldknbs (1,1)),
     .            (quant(1,1,2),oldktebs(1,1)),
     .            (quant(1,1,3),oldktibs(1,1)),
     .            (quant(1,1,4),oldkvhs (1,1))
c...this isn't completely consistent with using oldxxxx for calculating
c   gradients:
c      EQUIVALENCE (tquant(1,1),knds (1))
c      EQUIVALENCE (tquant(1,2),kteds(1))
c      EQUIVALENCE (tquant(1,3),ktids(1))
c      EQUIVALENCE (tquant(1,4),kvds (1))

      DATA initialize /.TRUE./

      SAVE

      tquant(1:MAXNDS,1) = knds (1:MAXNDS)
      tquant(1:MAXNDS,2) = kteds(1:MAXNDS)
      tquant(1:MAXNDS,3) = ktids(1:MAXNDS)
      tquant(1:MAXNDS,4) = kvds (1:MAXNDS)


      IF (initialize.AND.mode.EQ.17) THEN
        initialize = .FALSE.
c...    For mode=17:
        DO ir = irsep, irwall-1
          DO ik = 1, nks(ir)
            id = korpg(ik,ir)
            z1 = 0.5 * (zvertp(1,id) + zvertp(2,id))
            z2 = 0.5 * (zvertp(3,id) + zvertp(4,id))
            IF (((z1.GT.0.0.AND.z2.LT.0.0).OR.
     .           (z1.LT.0.0.AND.z2.GT.0.0)).AND.
     .          rs(ik,ir).GT.r0) THEN
              ikplane(ir) = ik
              WRITE(0,*) 'IKPLANE:',ir,ikplane(ir)
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      s       = 0.0
      thetav  = 0.0
      ik      = 0
      ir      = 0
      vali  = 1.0
      val   = 1.0
      valo  = 1.0
      widi    = 1.0
      wid     = 1.0
      wido    = 1.0
      deltapi = 1.0
      deltapt = 1.0
      deltapo = 1.0


      ik1 = 0

      CalcDist = 1.0
      IF (stopopt.EQ.60) RETURN


      irval = irval1
      sval  = sval1

10    CONTINUE

      ir   = irval
      s    = sval
c
c     Check for valid ring number:
c
      IF (ir.LT.irsep.OR.idring(ir).EQ.-1)
     .  CALL ER('CalcDist','Invalid ring number',*99)




      thetav = StoT(ir,s)



c
c  ... make sure that this is fool proof for broken grids...
c
      pdum = TtoP(ir,thetav,ik)








      IF (mode.GT.10.AND.mode.LE.16) THEN

        iri = irins (ik,ir)
        iro = irouts(ik,ir)

        IF     (ir.EQ.irsep.AND.ik.GE.ikti     ) THEN
          thetai = thetav - dthetg
          thetao = thetav
        ELSEIF (ir.EQ.nrs  .AND.ik.GE.ikti2(ir)) THEN
          thetai = thetav
          thetao = thetav + dthetg
        ELSE
          thetai = thetav
          thetao = thetav
        ENDIF

c
c This may not work so well for generalized grids.  Specifically, if
c there was only one wall to target ring, in which case the S mapping
c to the next inner ring would be way off... check for this...
c
c


        IF (iri.EQ.1.OR.iri.EQ.irtrap.OR.
     .      (ir.EQ.irsep.AND.ik.GT.ikto.AND.ik.LT.ikti)) THEN
          sval  = (s - ksb(0,ir)) * ksmaxs(iro) / ksmaxs(ir)
          irval = iro
          GOTO 10
        ELSEIF (iro.EQ.irwall) THEN
          sval  = (s - ksb(0,ir)) * ksmaxs(iri) / ksmaxs(ir)
          irval = iri
          GOTO 10
        ENDIF
c
c       Make sure that the denisty is not taken across the symmetry points
c       of neighbouring rings:
c
c        IF (ik.LE.ikmids(ir)) THEN
c          IF ((thetav.GT.thetag(ikmids(ir ),ir )).OR.
c     .        (thetav.GT.thetag(ikmids(iri),iri).AND.iri.NE.nrs  ).OR.
c     .        (thetav.GT.thetag(ikmids(iro),iro).AND.iro.NE.irsep)) THEN
c
c            sval = s - 1.0
c            GOTO 10
c           ENDIF
c        ELSE
c          IF ((thetav.LT.thetag(ikmids(ir ),ir )).OR.
c     .        (thetav.LT.thetag(ikmids(iri),iri).AND.iri.NE.nrs  ).OR.
c     .        (thetav.LT.thetag(ikmids(iro),iro).AND.iro.NE.irsep)) THEN
c
c            sval = s + 1.0
c            GOTO 10
c          ENDIF
c        ENDIF


        IF     (mode.EQ.11.OR.mode.EQ.12) THEN
c...denisty
          qopt = 1
c        ELSEIF (mode.EQ.12) THEN
c...electron temperature
c          qopt = 2
        ELSE
          STOP 'KSJFSKJF'
        ENDIF


c
c       Get quantity:
c
        CALL GetQuant(ir ,thetav,val ,NULL,NULL,quant(1,1,qopt),
     .                tquant(1,qopt),2)
        CALL GetQuant(iri,thetai,vali,iki ,ir  ,quant(1,1,qopt),
     .                tquant(1,qopt),2)
        CALL GetQuant(iro,thetao,valo,iko ,ir  ,quant(1,1,qopt),
     .                tquant(1,qopt),2)

c
c       Determine cell widths:
c
c I should really update this width estimate to be based on S, instead of
c always the cell center...
c
        wid  = 0.5 * CalcWidth(ik ,ir ,CENTER,TOTAL )
        widi =       CalcWidth(iki,iri,CENTER,SIDE23)
        wido =       CalcWidth(iko,iro,CENTER,SIDE14)

c
c       Use the cell containing the S point to *estimate* the grid
c       curvature at that point:
c

        CALL CalcCurvature(ik,ir,wid,widi,wido,deltapi,deltapo)

        deltapt = deltapi + deltapo
      ELSEIF (mode.EQ.9) THEN

c        WRITE(0,*) 'PROBLEM!'
c        CALL GetQuant(ir,thetav,val,NULL,NULL,quant,tquant,2)

        wid = 0.5 * CalcWidth(ik,ir,CENTER,TOTAL)
      ENDIF






















c        WRITE(0,*) 'MODE:',mode


      IF     (mode.EQ.1.OR.mode.EQ.21) THEN
        gradi    = 0.5 * 1.0
        grado    = 0.5 * 1.0
        CalcDist = gradi + grado

c        CalcDist = CalcDist * osm_dp4(ik,ir)

      ELSEIF (mode.EQ.2) THEN
        IF (kss(ik,ir).LT.0.5*ksmaxs(ir)) THEN
          wid = ABS(kss(ik,ir) - 0.25 * ksmaxs(ir))
        ELSE
          wid = ABS(kss(ik,ir) - 0.75 * ksmaxs(ir))
        ENDIF

        wid = wid / (0.25 * ksmaxs(ir))
        wid = EXP( wid / 0.2)

        gradi    = 0.5 / wid
        grado    = 0.5 / wid
        CalcDist = gradi + grado

      ELSEIF (mode.EQ.3) THEN

        CalcDist = MAX(1.0,ABS(pinion(ik,ir)))

c        wid = ABS(kss(ik,ir) - 0.5 * ksmaxs(ir)) / (0.5 * ksmaxs(ir))
c        wid = EXP(wid / 0.2)
c        gradi    = 0.5 / wid
c        grado    = 0.5 / wid
c        CalcDist = gradi + grado

      ELSEIF (mode.EQ.4) THEN
c
c       Peaked near target:
c
        midpt = SymmetryPoint(ir)

c        IF (ir.EQ.10) WRITE(0,*) '--> ',ik,midpt

        IF (ik.LE.midpt)THEN
          l1 = osm_dp2(IKLO,ir) * ksmaxs(ir)
          l2 = osm_dp1(IKLO,ir) * ksmaxs(ir)
          s1 = s

          CalcDist = MAX(-ABS(s1 - l1) / l2 + 1.0,0.0)

c...dp3
c          CalcDist = CalcDist * osm_dp4(ik,ir)
c          IF (s1.LT.l1) CalcDist = CalcDist * osm_dp3(IKLO,ir)


c          IF (ir.EQ.9) WRITE(0,*) ' -> ',l1,l2,s1,CalcDist
c          IF (ir.EQ.10) WRITE(0,*) ' -> ',CalcDist
        ELSE
          l1 = osm_dp2(IKHI,ir) * ksmaxs(ir)
          l2 = osm_dp1(IKHI,ir) * ksmaxs(ir)
          s1 = ksmaxs(ir) - s

c          IF (ir.EQ.10) WRITE(0,*) ' -> ',l1,l2,s1

          CalcDist = MAX(-ABS(s1 - l1) / l2 + 1.0,0.0)

c...dp3
c          CalcDist = CalcDist * osm_dp4(ik,ir)
c          IF (s1.LT.l1) CalcDist = CalcDist * osm_dp3(IKHI,ir)


c          IF (ir.EQ.10) WRITE(0,*) ' -> ',CalcDist
        ENDIF


c        IF (ik.LE.midpt)THEN
c          l1 = GetL1  (IKLO,ir) * ksmaxs(ir)
c          l2 = osm_dp1(IKLO,ir) * ksmaxs(ir)
c          s1 = s
c
c          CalcDist = -MAX(l2 - (s1 - l1),0.0)
c        ELSE
c          l1 = GetL1  (IKHI,ir) * ksmaxs(ir)
c          l2 = osm_dp1(IKHI,ir) * ksmaxs(ir)
c          s1 = ksmaxs(ir) - s
c
c          CalcDist = -MAX(l2 - (s1 - l1),0.0)
c        ENDIF


      ELSEIF (mode.EQ.5) THEN
c
c       Gaussian:
c
        midpt = SymmetryPoint(ir)

c        IF (ir.EQ.10) WRITE(0,*) '--> ',ik,midpt

        IF (ik.LE.midpt)THEN
c          l1 = (GetL1(IKLO,ir) + osm_dp2(IKLO,ir)) * ksmaxs(ir)
c          l2 =                   osm_dp1(IKLO,ir)  * ksmaxs(ir)

          l1 = osm_dp2(IKLO,ir) * ksmaxs(ir)
          l2 = osm_dp1(IKLO,ir) * ksmaxs(ir)
          s1 = s

c          IF (ir.EQ.10) WRITE(0,*) ' -> ',l1,l2,s1

          CalcDist = EXP(-0.5 * MIN(100.0,((s1 - l1) / l2)**2.0))

          IF (s1.LT.l1) CalcDist = CalcDist * osm_dp3(IKLO,ir)

c          IF (ir.EQ.10) WRITE(0,*) ' -> ',CalcDist
        ELSE
c          l1 = (GetL1(IKHI,ir) + osm_dp2(IKHI,ir)) * ksmaxs(ir)
c          l2 =                   osm_dp1(IKHI,ir)  * ksmaxs(ir)
          l1 = osm_dp2(IKHI,ir) * ksmaxs(ir)
          l2 = osm_dp1(IKHI,ir) * ksmaxs(ir)
          s1 = ksmaxs(ir) - s

c          IF (ir.EQ.10) WRITE(0,*) ' -> ',l1,l2,s1

          CalcDist = EXP(-0.5 * MIN(100.0,((s1 - l1) / l2)**2.0))

          IF (s1.LT.l1) CalcDist = CalcDist * osm_dp3(IKHI,ir)

c          IF (ir.EQ.10) WRITE(0,*) ' -> ',CalcDist
        ENDIF

      ELSEIF (mode.EQ.6) THEN
c
c       PINQe distribution:
c
c...new (wrong place)
         CalcDist = MAX(1.0,ABS(pinqe(ik,ir) + osmqe(ik,ir)))

      ELSEIF (mode.EQ.7.OR.mode.EQ.8) THEN

        IF     (ik.EQ.1      .AND.s.LT.kss(ik,ir)) THEN
          CalcDist = ABS((pinqe(1,ir) + osmqe(1,ir)) * 
     .                   (kss  (1,ir) - ksb  (0,ir)))
        ELSEIF (ik.EQ.nks(ir).AND.s.GT.kss(ik,ir))THEN
          STOP 'TRUE MOCK POWER NOT READY HERE'
          WRITE(88,*) 'HERE->',ik,ir,pinqe(ik,ir)
          CalcDist = ABS(pinqe(ik,ir) + osmqe(ik,ir))
        ELSEIF (s.LT.kss(ik,ir)) THEN
          CalcDist = 0.5 * ABS(pinqe(ik-1,ir) + osmqe(ik-1,ir) + 
     .                         pinqe(ik  ,ir) + osmqe(ik  ,ir))
          STOP 'SHOULD NOT BE HERE FOR MODE 7'
        ELSEIF (s.GT.kss(ik,ir)) THEN
          CalcDist = 0.5 * ABS(pinqe(ik  ,ir) + osmqe(ik  ,ir) + 
     .                         pinqe(ik+1,ir) + osmqe(ik+1,ir))
c          WRITE(88,*) 'HERE A->',ik,ir,s,pinqe(ik,ir)
        ELSE
          CALL ER('CalcDist','Something isn''t right',*99)
        ENDIF

        IF ((mode.EQ.8.AND.stopopt3.EQ.18.AND.ir.EQ.27).OR.
     .      (mode.EQ.8.AND.(stopopt3.EQ.19.OR.stopopt3.eq.20))) THEN

          STOP 'DEVELOPMENT 000A'

          sfrac = 1.0 
          IF     (ik.LE.osm_sympt(ir).AND.ik.GE.ikfluid(IKLO,ir)) THEN
            sfrac = (ksb(osm_sympt(ir),ir) - s) / 
     .              (ksb(osm_sympt(ir),ir) - kss(ikfluid(IKLO,ir),ir))
          ELSEIF (ik.GT.osm_sympt(ir).AND.ik.LE.ikfluid(IKHI,ir)) THEN
            sfrac = (s - ksb(osm_sympt(ir),ir)) / 
     .              (kss(ikfluid(IKHI,ir),ir) - ksb(osm_sympt(ir),ir))
          ENDIF

          sfrac  = MIN(MAX(0.0,sfrac),1.0)
c          sfrac  = MIN((2.0*sfrac)**2,1.0)
            
c          IF (ik.EQ.osm_sympt(ir)) sfrac = 0.0

          CalcDist = CalcDist * sfrac

c          WRITE(0,'(A,5I6,F10.4)') '-->',ik,ir,
c     .      ikfluid(IKLO,ir),ikfluid(IKHI,ir),
c     .      osm_sympt(ir),sfrac

        ENDIF


c      ELSEIF (mode.EQ.8) THEN
c
c       Inverse field line separation:
c
c        gradi    = 0.5 / wid
c        grado    = 0.5 / wid
c        CalcDist = gradi + grado



      ELSEIF (mode.EQ.9) THEN
c
c       Major radius effect:
c
        gradi    = 0.5
        grado    = 0.5

        CalcDist = (gradi + grado) * rs(ik,ir)

cc
cc       Inverse field line separation and major radius effect:
cc
c        gradi    = 0.5 / wid
c        grado    = 0.5 / wid
c
c        CalcDist = (gradi + grado) * rs(ik,ir)


c        WRITE(0,*) ik,ir,gradi+grado,CalcDist
c        WRITE(PINOUT,*) ik,ir,gradi+grado,rs(ik,ir),CalcDist

      ELSEIF (mode.EQ.10) THEN
c
c       Direct with denisty:


        CalcDist = oldknbs2(ik,ir)

      ELSEIF (mode.EQ.11) THEN
        gradi    = 0.5 * (vali - val)
        grado    = 0.5 * (valo - val)
        CalcDist = gradi + grado

      ELSEIF (mode.EQ.12) THEN
        gradi    = (vali - val) / widi
        grado    = (valo - val) / wido

        CalcDist = (gradi + grado) / wid

      ELSEIF (mode.EQ.13) THEN
        gradi    = deltapi / deltapt * (vali - val)
        grado    = deltapo / deltapt * (valo - val)
        CalcDist = gradi + grado

      ELSEIF (mode.EQ.14) THEN
c        gradi    = deltapi / deltapt * (vali - val) / wid
c        grado    = deltapo / deltapt * (valo - val) / wid
c        CalcDist = gradi + grado
        CalcDist = 1.0

      ELSEIF (mode.EQ.15) THEN
        gradi    = 0.5 * (vali - val) / (widi + wid)
        grado    = 0.5 * (valo - val) / (wido + wid)
        CalcDist = gradi + grado

      ELSEIF (mode.EQ.16) THEN
        gradi    = deltapi / deltapt * (vali - val) / (widi + wid)
        grado    = deltapo / deltapt * (valo - val) / (wido + wid)
        CalcDist = gradi + grado

      ELSEIF (mode.EQ.17) THEN
c...    Distribution is a maximum at the outer midplane, and goes to zero at the
c       IKBOUND:
        gradi    = 0.5
        grado    = 0.5

        IF (ir.GE.irsep.AND.ir.LT.irwall) THEN
          p1 = kps(ikbound(ir,IKLO),ir) 
          p2 = kps(ikbound(ir,IKHI),ir) 
          pmid = kps(ikplane(ir),ir)
          
          IF (kps(ik,ir).LT.pmid) THEN
            gradi = MAX(0.0,kps(ik,ir) - p1) / ABS(pmid - p1)
c            gradi = MAX(0.05,gradi)
          ELSE
            gradi = MAX(0.0,p2 - kps(ik,ir)) / ABS(pmid - p2)
          ENDIF

          gradi = gradi**1.5

          CalcDist = gradi
        ELSE 
          CALL ER('CalcDist','Option 17 not ready for the PFZ',*99)
        ENDIF

      ELSEIF (mode.EQ.18) THEN
c...    Power between cutpoints:
        gradi    = 0.5
        grado    = 0.5

        IF     ((iflexopt(6).EQ.20.OR.iflexopt(6).EQ.24).AND.
     .          ikfluid(IKHI,ir).NE.nks(ir)) THEN
          IF (ik.GT.ikbound(ir,IKLO).AND.ik.LT.ikfluid(IKHI,ir)) THEN
c        IF     ((iflexopt(6).EQ.20.OR.iflexopt(6).EQ.24).AND.
c     .          ir.GE.14.AND.ir.LE.16) THEN
c          IF (ik.GT.ikbound(ir,IKLO).AND.ik.LT.ikti2(ir)) THEN
            CalcDist = 1.0
          ELSE
            CalcDist = 0.0
          ENDIF
        ELSEIF (iflexopt(6).EQ.21.AND.ir.GE.19.AND.ir.LE.24) THEN
          IF (ik.GE.ikfluid(IKLO,ir).AND.ik.LE.nks(ir)) THEN
            CalcDist = 1.0
          ELSE
            CalcDist = 0.0
          ENDIF
        ELSE 
          CalcDist = 1.0
        ENDIF

      ELSEIF (mode.EQ.19) THEN

        IF     (ik.EQ.1      .AND.s.LT.kss(ik,ir)) THEN
          CalcDist = ABS((pinqe(1,ir)) * 
     .                   (kss  (1,ir) - ksb  (0,ir)))
        ELSEIF (ik.EQ.nks(ir).AND.s.GT.kss(ik,ir))THEN
          STOP 'TRUE MOCK POWER NOT READY HERE'
        ELSEIF (s.LT.kss(ik,ir)) THEN
          CalcDist = 0.5 * ABS(pinqe(ik-1,ir) + pinqe(ik  ,ir))
          STOP 'SHOULD NOT BE HERE FOR MODE 7'
        ELSEIF (s.GT.kss(ik,ir)) THEN
          CalcDist = 0.5 * ABS(pinqe(ik  ,ir) + pinqe(ik+1,ir))
        ELSE
          CALL ER('CalcDist','Something isn''t right',*99)
        ENDIF

      ELSEIF (mode.EQ.20) THEN
c...    Uniform distribution of cross-field particle source, unless there is
c       flow reversal on the outside, in which case redistribute
c       the cross-field flux between the x-point and the outer target so that 
c       the flow reversal source fraction is limited according to ABS(rflexopt(1)).  
c       Note: there is the implicit assumption that the majority of the ionisation
c       and recombination on the outside occurs between the target and the
c       x-point:

        IF     (disindex(ir).EQ.0.AND.
     .          osm_model(IKLO,ir).EQ.24.AND.
     .          osm_model(IKHI,ir).EQ.22.AND.
c... separatrix only for now!
     .          ir.EQ.irsep) THEN

          WRITE(0     ,*) '(((((((((((((CHECKING!))))))))))))',ir
          WRITE(PINOUT,*) '(((((((((((((CHECKING!))))))))))))',ir

          iks = osm_sympt(ir)

c...      Assign uniform distribution as the default:
          disindex(ir) = -1

c...      Check the ring to see if there is significant flow reversal
c         on the outside:
          tarflx1 = ABS(knds(idds(ir,2)) * kvds(idds(ir,2)))
          tarflx2 = ABS(knds(idds(ir,1)) * kvds(idds(ir,1)))

          CALL CalcIntegral4(pinrec,iks+1,nks(ir),ir,recsrc,4)
          CALL CalcIntegral3(pinion,iks+1,nks(ir),ir,ionsrc,4)

c          tolflow = ABS(rflexopt(1))
          tolflow = 0.0

          IF (ABS(ionsrc)/ABS(recsrc+tarflx2).GT.
     .        1.0+tolflow) THEN
c...        Excessive flow reversal detected:

c...        Define cell range between x-point and the target (not sure if this
c           will work for JET grids):
            disindex(ir) = 0
            DO ik1 = iks+1, nks(ir)
              IF (zs(ik1,ir).GT.zxp) disindex(ir) = ik1 + 1
            ENDDO
           
            IF (disindex(ir).GT.0) THEN
              ds1 = ksb(disindex(ir)-1,ir) - kss(ikbound(ir,IKLO),ir)
              ds2 = ksmaxs(ir) - ksb(disindex(ir)-1,ir)

              intcfp = ionsrc - (1.0 + tolflow) * (recsrc + tarflx2)

              CALL CalcIntegral4(pinrec,1,nks(ir),ir,recsrc,2)
              CALL CalcIntegral3(pinion,1,nks(ir),ir,ionsrc,2)

              intcfptot = ionsrc - recsrc - tarflx1 - tarflx2

              dislev(ir) = (ds2 * (intcfptot - intcfp)) / (ds1 * intcfp)

c              WRITE(0,*) 'CALCDIST:',ir,dislev(ir),disindex(ir)
c              WRITE(0,*) 'CALCDIST:',intcfp,intcfptot
c              WRITE(0,*) 'CALCDIST:',ds1,ds2

c          CALL CalcIntegral4(pinrec,1,iks,ir,recsrc,4)
c          CALL CalcIntegral3(pinion,1,iks,ir,ionsrc,4)
c              WRITE(0,*) 'CALCDIST:',ionsrc,recsrc,tarflx1
c              WRITE(0,*) 'CALCDIST:',ionsrc-recsrc-tarflx1

c              STOP  'testiog'
            ENDIF

          ENDIF

        ELSEIF (disindex(ir).EQ.0) THEN
c...      Assign uniform distribution:
          disindex(ir) = -1
        ENDIF

        IF (disindex(ir).EQ.-1) THEN
c...      Uniform distribution:
          CalcDist = 1.0
        ELSE
c...      Step distribution:
          IF (ik.LT.disindex(ir)) THEN
            CalcDist = dislev(ir)
          ELSE
            CalcDist = 1.0
          ENDIF
        ENDIF

      ELSEIF (mode.EQ.22) THEN
c...    Uniform power above the x-point only:

c...    Problem with JET grids?
        IF (zs(ik,ir).GT.zxp) THEN
          CalcDist = 1.0
        ELSE
          CalcDist = 0.0
        ENDIF
      ELSEIF (mode.EQ.-1) THEN
        gradi    = (vali - val) / (widi + wid)
        grado    = (valo - val) / (wido + wid)
        grad2    = (gradi + grado) / (2.0 * wid)
        CalcDist = grad2

      ELSE
        CALL ER('CalcDist','Invalid mode',*99)
      ENDIF




      IF (mode.GT.0) THEN


c        WRITE(0,*) '============================='
c        WRITE(0,*) s
c        WRITE(0,*) thetav
c        WRITE(0,*) ik
c        WRITE(0,*) ir
c        WRITE(0,*) CalcDist
c        WRITE(0,*) vali
c        WRITE(0,*) val
c        WRITE(0,*) valo
c        WRITE(0,*) valo
c        WRITE(0,*) widi
c        WRITE(0,*) wid
c        WRITE(0,*) wido
c        WRITE(0,*) ik1
c        WRITE(0,*) deltapi
c        WRITE(0,*) deltapo
c        WRITE(0,*) deltapt



c        WRITE(SLOUT,'(A,2F8.3,1X,2I4,1X,1P,E10.2,
c     .                1X,3E10.2,0P,1X,3F7.5,
c     .                1X,I4,1X,2F4.2)')
c     .    'G: ',
c     .    s,thetav,ik,ir,CalcDist,
c     .    vali,val,valo,widi,wid,wido,
c     .    ik1,deltapi/deltapt,deltapo/deltapt


      ENDIF

      RETURN
c
c     Error code:
99    WRITE(EROUT,'(5X,A,I4)') 'IR = ',ir
      STOP 'End of CalcDist'
      END
c
c ======================================================================
c
c
c ======================================================================
c
c subroutine: CalcCurvature
c
c
      SUBROUTINE CalcCurvature(ik,ir,wid,widi,wido,deltapi,deltapo)
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER ik,ir
      REAL    wid,widi,wido

      REAL TtoP,StoT

      INTEGER ik1,iro,iri
      REAL    tmini,tmino,tmaxi,tmaxo,t1i,t2i,t1o,t2o,p1i,p2i,
     .        p1o,p2o,pdum,deltapi,deltapo,p1,p2

      iri = irins (ik,ir)
      iro = irouts(ik,ir)

      ik1 = ik
10    CONTINUE

      t1i = StoT(ir,ksb(ik1-1,ir))
      t2i = StoT(ir,ksb(ik1  ,ir))
      t1o = t1i
      t2o = t2i

      IF     (ir.EQ.irsep.AND.ik1.GE.ikti2(ir)) THEN
        t1i = t1i - dthetg
        t2i = t2i - dthetg
      ELSEIF (ir.EQ.nrs  .AND.ik1.GE.ikti2(ir)) THEN
        t1o = t1o + dthetg
        t2o = t2o + dthetg
      ENDIF

      tmini = thetat(idds(iri,2))
      tmino = thetat(idds(iro,2))
      tmaxi = thetat(idds(iri,1))
      tmaxo = thetat(idds(iro,1))

      p1      = kpb(ik1-1,ir)
      p2      = kpb(ik1  ,ir)
      p1i     = 0.0
      p2i     = 0.0
      p1o     = 0.0
      p2o     = 0.0
      deltapi = NIL
      deltapo = NIL

      IF (t1i.GE.tmini.AND.t2i.LE.tmaxi) THEN
        p1i = TtoP(iri,t1i,NULL)
        p2i = TtoP(iri,t2i,NULL)
        deltapi = (wid * (p2i - p1i) + widi * (p2 - p1)) / (wid + widi)
      ENDIF

      IF (t1o.GE.tmino.AND.t2o.LE.tmaxo) THEN
        p1o = TtoP(iro,t1o,NULL)
        p2o = TtoP(iro,t2o,NULL)
        deltapo = (wid * (p2o - p1o) + wido * (p2 - p1)) / (wid + wido)
      ENDIF

c
c     If either side is invalid then move to the nearest place along
c     the ring with a valid estimate:
c
      IF (deltapo.EQ.NIL.OR.deltapi.EQ.NIL) THEN
        IF     (ik1.LT.ikmids(ir)) THEN
          ik1 = ik1 + 1
        ELSEIF (ik1.GT.ikmids(ir)) THEN
          ik1 = ik1 - 1
        ELSE
          CALL ER('CalcDist','Cannot find curvature estimate',*99)
        ENDIF
c        WRITE(SLOUT,'(5X,A,2I4,2X,I4)') 'IR,IK IK1 = ',ir,ik,ik1
        GOTO 10
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c
c ======================================================================
c
c function: StoT
c
c Returns the THETA value corresponding to a given S value for the
c ring IR.
c
c Aug 14, 97 - The current method of handling the x-point
c region in the core will be inaccurate if the method of determining
c the distance between points on neighbouring rings differs from
c from the current implimentation, which uses an estimate of
c the cell widths (sorry about the run-on sentence).
c
      REAL FUNCTION StoT(ir,s)

      REAL    s
      INTEGER ir

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER ik
c
c     Boundary rings not allowed:
c
      IF (ir.EQ.1.OR.ir.EQ.irtrap.OR.ir.EQ.irwall)
     .  CALL ER('StoT','Invalid ring',*99)
c
c     Find cell index:
c
      DO ik = 1, nks(ir)
        IF (s.LE.kss(ik,ir)) GOTO 10
      ENDDO
10    CONTINUE
c
c     Interpolate linearly to find StoT:
c
      IF (ik.EQ.1) THEN
        StoT = thetat(idds(ir,2)) +
     .               (s             - ksb   (0,ir)      ) *
     .               (thetag(ik,ir) - thetat(idds(ir,2))) /
     .               (kss   (ik,ir) - ksb   (0,ir)      )
      ELSEIF (ik.EQ.nks(ir)+1) THEN
        ik = nks(ir)
        StoT = thetag(ik,ir) +
     .               (s                  - kss   (ik,ir)) *
     .               (thetat(idds(ir,1)) - thetag(ik,ir)) *
     .               (ksb   (ik,ir)      - kss   (ik,ir))
      ELSE
        StoT = thetag(ik-1,ir) +
     .               (s             - kss   (ik-1,ir)) *
     .               (thetag(ik,ir) - thetag(ik-1,ir)) /
     .               (kss   (ik,ir) - kss   (ik-1,ir))
      ENDIF
c
c     Check that StoT is not out of bounds:
c
      IF ((ir.GE.irsep.AND.(StoT.LT.thetat(idds(ir,2)).OR.
     .     StoT.GT.thetat(idds(ir,1)))).OR.
     .    (ir.LT.irsep.AND.(StoT.LT.thetag(1,ir).OR.
     .     StoT.GT.thetag(nks(ir),ir))))
     .  CALL ER('StoT','Calculated value out of bounds',*99)

      RETURN
99    WRITE(EROUT,'(5X,A,3I4)')    'IK,IR NKS = '   ,ik,ir,nks(ir)
      WRITE(EROUT,'(5X,A,2F15.7)') 'S,theta = ',s,StoT
      WRITE(EROUT,'(5X,A,2F15.7)') 'thetat2,thetat1 = ',
     .  thetat(idds(ir,2)),thetat(idds(ir,1))
      WRITE(EROUT,'(5X,A,2F15.7)') 'thetag-1,thetag = ',
     .  thetag(MAX(1,ik-1),ir),thetag(ik,ir)
      WRITE(EROUT,'(5X,A,2F15.7)') 'ksb0,ksbmax = ',
     .  ksb(0,ir),ksb(nks(ir),ir)
      WRITE(EROUT,'(5X,A,2F15.7)') 'kss-1,kss = ',
     .  kss(MAX(1,ik-1),ir),kss(ik,ir)

      WRITE(6    ,'(5X,A,3I4)')    'IK,IR NKS = '   ,ik,ir,nks(ir)
      WRITE(6    ,'(5X,A,2F15.7)') 'S,theta = ',s,StoT
      WRITE(6    ,'(5X,A,2F15.7)') 'thetat2,thetat1 = ',
     .  thetat(idds(ir,2)),thetat(idds(ir,1))
      WRITE(6    ,'(5X,A,2F15.7)') 'thetag-1,thetag = ',
     .  thetag(MAX(1,ik-1),ir),thetag(ik,ir)
      WRITE(6    ,'(5X,A,2F15.7)') 'ksb0,ksbmax = ',
     .  ksb(0,ir),ksb(nks(ir),ir)
      WRITE(6    ,'(5X,A,2F15.7)') 'kss-1,kss = ',
     .  kss(MAX(1,ik-1),ir),kss(ik,ir)

      CALL OutputData(87,'Failure in StoT')
      STOP
      END
c
c ======================================================================
c
c
c ======================================================================
c
c function: TtoS
c
c Returns the S value corresponding to a given THETA value for the
c ring IR.
c
c some liberties around the x-point
c check that ksb is defined for core rings
c
      REAL FUNCTION TtoS(ir,thetav)

      REAL    thetav
      INTEGER ir

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER ik
c
c     Check for valid IR value (virtual rings not allowed):
c
      IF (idring(ir).EQ.-1) CALL ER('StoT','Boundary ring',*99)
c
c     Determine if THETA is out of bounds:
c
      IF (ir.GE.irsep) THEN
        IF (thetav.LT.thetat(idds(ir,2))) THEN
          TtoS = ksb(0,ir)
          RETURN
        ELSEIF (thetav.GT.thetat(idds(ir,1))) THEN
          TtoS = ksb(nks(ir),ir)
          RETURN
        ENDIF
      ELSE
        IF (thetav.LE.thetag(1,ir)) THEN
          TtoS = kss(1,ir)
          RETURN
        ELSEIF (thetav.GE.thetag(nks(ir),ir)) THEN
          TtoS = kss(nks(ir),ir)
          RETURN
        ENDIF
      ENDIF
c
c     Find cell index:
c
      DO ik = 1, nks(ir)
        IF (thetav.LE.thetag(ik,ir)) GOTO 10
      ENDDO
10    CONTINUE
c
c     Interpolate linearly to find S value:
c
      IF     (ik.EQ.1        ) THEN
        TtoS = ksb(0,ir) +
     .         (thetav        - thetat(idds(ir,2))) *
     .         (kss   (ik,ir) - ksb   (0,ir)      ) /
     .         (thetag(ik,ir) - thetat(idds(ir,2)))
      ELSEIF (ik.EQ.nks(ir)+1) THEN
        ik = nks(ir)
        TtoS = kss(ik,ir) +
     .         (thetav             - thetag(ik,ir)) *
     .         (ksb   (nks(ir),ir) - kss   (ik,ir)) *
     .         (thetat(idds(ir,1)) - thetag(ik,ir))
      ELSE
        TtoS = kss(ik-1,ir) +
     .         (thetav        - thetag(ik-1,ir)) *
     .         (kss   (ik,ir) - kss   (ik-1,ir)) /
     .         (thetag(ik,ir) - thetag(ik-1,ir))
      ENDIF

c Debug:
c   write(0,*) 's stuff: ',ir,ksb(ik-1,ir),TtoS,ksb(ik,ir)

c
c Check that TtoS is not out of bounds:
c
      IF (TtoS.LT.ksb(0,ir).OR.TtoS.GT.ksb(nks(ir),ir)) THEN
        WRITE(0,*) 'ERROR (TtoS): Calculated S value is out '//
     .             'of bounds'
        WRITE(0,*) '      IR = ',ir
        WRITE(0,*) '       S = ',TtoS
        WRITE(0,*) '  KSBmin = ',ksb(0,ir)
        WRITE(0,*) '  KSBmax = ',ksb(nks(ir),ir)
        STOP
      ENDIF

      RETURN
99    WRITE(EROUT,'(5X,A,I4)') 'IR = ',ir
      STOP
      END
c
c ======================================================================
c
c subroutine: TtoP
c
c Determine P, the poloidal co-ordinate, and IK, the cell index,
c for ring IR given THETA.
c
      REAL FUNCTION TtoP(ir,thetav,ik)

      IMPLICIT none

c     Input:
      REAL    thetav
      INTEGER ir

c     Output:
      INTEGER ik

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      REAL p
c
c     Check for valid IR value:
c
      IF (idring(ir).EQ.-1) CALL ER('StoT','Boundary ring',*99)
c
c     Check that THETAV is valid for the ring, otherwise return
c     target location:
c
      IF (ir.GE.irsep) THEN
        IF (thetav.LT.thetat(idds(ir,2))) THEN
          ik   = 1
          TtoP = kpb(0,ir)
          RETURN
        ELSEIF (thetav.GT.thetat(idds(ir,1))) THEN
          ik   = nks(ir)
          TtoP = kpb(ik,ir)
          RETURN
        ENDIF
      ELSE
        IF (thetav.LE.thetag(1,ir)) THEN
          ik   = 1
          TtoP = kps(ik,ir)
          RETURN
        ELSEIF (thetav.GT.thetag(nks(ir),ir)) THEN
          ik   = nks(ir)
          TtoP = kps(ik,ir)
          RETURN
        ENDIF
      ENDIF
c
c     Find cell index:
c
      DO ik = 1, nks(ir)
        IF (thetav.LE.thetag(ik,ir)) GOTO 10
      ENDDO
10    CONTINUE
c
c     Linear interpolation:
c
      IF (ik.EQ.1) THEN
        p = kpb(0,ir) + (thetav        - thetat(idds(ir,2))) *
     .                  (kps   (ik,ir) - kpb   (0,ir)      ) /
     .                  (thetag(ik,ir) - thetat(idds(ir,2)))
      ELSEIF (ik.EQ.nks(ir)+1) THEN
        ik = nks(ir)
        p  = kps(ik,ir) + (thetav             - thetag(ik,ir)) *
     .                    (kpb   (ik,ir)      - kps   (ik,ir)) /
     .                    (thetat(idds(ir,1)) - thetag(ik,ir))
      ELSE
        p = kps(ik-1,ir) + (thetav        - thetag(ik-1,ir)) *
     .                     (kps   (ik,ir) - kps   (ik-1,ir)) /
     .                     (thetag(ik,ir) - thetag(ik-1,ir))
c
c       Adjust IK based on P if necessary:
c
        IF (p.LT.kpb(ik-1,ir)) ik = ik - 1
      ENDIF
c
c     Check that P and IK are not out of bounds:
c
      IF (p.LT.kpb(0,ir).OR.p.GT.kpb(nks(ir),ir))
     .  CALL ER('TtoP','Calculated value out of bounds',*99)

      IF (ik.LT.1.OR.ik.GT.nks(ir))
     .  CALL ER('TtoP','Cell index out of bounds',*99)

      TtoP = p

      RETURN
99    WRITE(EROUT,'(5X,A,I4)') 'IR = ',ir
      STOP
      END
c
c ======================================================================
c
c
c ======================================================================
c
c subroutine: GetQuant
c
c
      SUBROUTINE GetQuant(ir,thetav,density,ik,ir1,quant,tquant,vopt)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

c     Input:
      REAL    thetav,quant(MAXNKS,MAXNRS),tquant(MAXNDS)
      INTEGER ir,ir1,vopt

c     Output:
      REAL    density
      INTEGER ik

      REAL    CalcWidth,TtoP,StoT

      INTEGER valopt,idum1,idum2,idum3
      REAL    val1,val2,val3,val4,p,p1,p2,theta1

      DATA valopt,idum1,idum2,idum3 /0,0,0,0/
      DATA val1,val2,val3,val4,p,p1,p2,theta1
     .     /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/

      ik      = 0
      density = 0.0

      theta1 = thetav
      valopt = 1
c      valopt = vopt
c
c     Check for valid IR value (virtual rings not allowed):
c
      IF (ir.EQ.1.OR.ir.EQ.irtrap.OR.ir.EQ.irwall)
     .  CALL ER('StoT','Boundary rings not allowed',*99)
c
c     Make sure that THETA1 is not out of bounds:
c

c
c     Find quantity and cell index:
c

      IF (valopt.EQ.1) THEN
c
c       Assign value at cell center:
c
        p       = TtoP(ir,theta1,ik)
        density = quant(ik,ir)

      ELSEIF (valopt.EQ.2) THEN
c
c       Linear interpolation between cell centers:
c
        IF     (theta1.LT.thetat(idds(ir,2))) THEN
          ik = 1

          IF (ir1.EQ.NULL)
     .      CALL ER('CalcQuant','Problems with low target value',*99)

          val1 = tquant(idds(ir1,1))
          val2 = tquant(idds(ir ,1))
          p1   = thetat(idds(ir1,1))
          p2   = thetat(idds(ir ,1))
          p    = theta1


c          WRITE(SLOUT,'(5X,A,2I4,3F10.5,1P,2E15.7)')
c     .      'Fancy low : ',ir,ir1,p,p1,p2,val1,val2


        ELSEIF (theta1.GT.thetat(idds(ir,1))) THEN
          ik = nks(ir)

          IF (ir1.Eq.NULL)
     .      CALL ER('CalcQuant','Problems with high target value',*99)


          val1 = tquant(idds(ir1,1))
          val2 = tquant(idds(ir ,1))
          p1   = thetat(idds(ir1,1))
          p2   = thetat(idds(ir ,1))
          p    = theta1



c          WRITE(SLOUT,'(5X,A,2I4,3F10.5,1P,2E15.7)')
c     .      'Fancy high: ',ir,ir1,p,p1,p2,val1,val2

        ELSE
          p = TtoP(ir,theta1,ik)

          IF (p.LT.kps(ik,ir).AND.ik.EQ.1) THEN
            p1   = kps (ik,ir)
            p2   = kpb (0 ,ir)
            val1 = quant(ik,ir)
            if (ir.ge.irsep) then
               val2 = tquant(idds(ir,2))
            else
               val2 = (quant(nks(ir)-1,ir)+quant(ik,ir))/2.0
            endif
          ELSEIF (p.GT.kps(ik,ir).AND.ik.EQ.nks(ir)) THEN
            p1   = kpb (ik,ir)
            p2   = kps (ik,ir)
            if (ir.ge.irsep) then
              val1 = tquant(idds(ir,1))
            else
              val1 = (quant(2,ir)+quant(ik,ir))/2.0
            endif
            val2 = quant(ik,ir)
          ELSE
            IF (p.GE.kps(ik,ir)) THEN
              p1   = kps (ik+1,ir)
              p2   = kps (ik  ,ir)
              val1 = quant(ik+1,ir)
              val2 = quant(ik  ,ir)
            ELSE
              p1   = kps (ik  ,ir)
              p2   = kps (ik-1,ir)
              val1 = quant(ik  ,ir)
              val2 = quant(ik-1,ir)
            ENDIF
          ENDIF
        ENDIF

        density = val2 + (p - p2) / (p1 - p2) * (val1 - val2)

      ELSEIF (valopt.EQ.3) THEN
        density = quant(ik,ir)

      ELSE
        CALL ER('GetQuant','Invalid option',*99)
      ENDIF

      RETURN
99    WRITE(EROUT,'(5X,A,3I5)')
     .  'IK IR IR1 = ',ik,ir,ir1
      WRITE(EROUT,'(5X,A,F10.5)')
     .  'THETA1 = ',theta1
      STOP
      END
c
c ======================================================================
c
c subroutine: GetPotential
c
c Estimate the plasma potential from (an estimate) of the parallel 
c electric field (KES).
c
      SUBROUTINE CalcPotential
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER mode,ik,ir
      REAL    deltas,efield

 

c      STOP 'CALCULATING POTENTIAL: NEED SHEATH DROP INCLUDED!'


      mode = 1

      DO ir = irsep, irwall
        IF (idring(ir).EQ.BOUNDARY) CYCLE

c...      Assume          
        IF     (mode.EQ.1) THEN
c...      Calculate from low IK target to high IK target:
          DO ik = 1, nks(ir)
            IF (ik.EQ.1) THEN
              deltas = kss(ik,ir) - ksb(0 ,ir)
              efield = kes(ik,ir) / (qtim * qtim * emi / crmi)
              osmpot(ik,ir) = -efield * deltas
            ELSE
              deltas = kss(ik,ir) - kss(ik-1,ir)
              efield = kes(ik,ir) / (qtim * qtim * emi / crmi)
              osmpot(ik,ir) = osmpot(ik-1,ir) - efield * deltas

              efield = -(ktebs(ik,ir) - ktebs(ik-1,ir)) / deltas -
     .                (1 / knbs(ik,ir)) *
     .         (ktebs(ik,ir)*knbs(ik,ir) -
     .          ktebs(ik-1,ir)*knbs(ik-1,ir)) / deltas

            ENDIF

      
c          DO ik = 1, nks(ir)/2
c            IF (ik.EQ.1) THEN
c              deltas = kss(ik,ir) - ksb(0 ,ir)
c              osmpot(ik,ir+1) = -kes(ik,ir) * deltas
c            ELSE
c              deltas = kss(ik,ir) - kss(ik-1,ir)
c              osmpot(ik,ir+1) = osmpot(ik-1,ir+1) -kes(ik,ir) * deltas
c            ENDIF
c          ENDDO
c          DO ik = nks(ir), nks(ir)/2+1, -1
c            IF (ik.EQ.1) THEN
c              deltas = kss(ik,ir) - ksb(0 ,ir)
c              osmpot(ik,ir+1) = -kes(ik,ir) * deltas
c            ELSE
c              deltas = kss(ik,ir) - kss(ik-1,ir)
c              osmpot(ik,ir+1) = osmpot(ik-1,ir+1) -kes(ik,ir) * deltas
c            ENDIF
c          ENDDO


            WRITE(0,*) 'POT:',ik,ir,osmpot(ik,ir),
     .        kes(ik,ir)/(QTIM * QTIM * EMI / CRMI), 
     .        efield
          ENDDO



        ENDIF

      ENDDO

      RETURN
 99   STOP
      END

c subroutine: CalcPotential2
c
c
c jdemod - same purpose - include target potentials
c                       - calculate independently from each target since
c                         the purpose is primarily to calculate drifts
c                       - will need to look at discrepancies at mid-point
c                         might be some useful physics or at least 
c                         inner/outer alignment that can be pulled out
c
c Estimate the plasma potential from (an estimate) of the parallel 
c electric field (KES).
c
      SUBROUTINE CalcPotential2
      use mod_interpolate
      use plasma_overlay
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'
c     
c     jdemod - put the potential and drift variables in the common block driftvel
c     
c     this routine is for use in DIVIMP for calculating ExB related drifts
c     Steve's original calcpotential is used in OUT
c     
      include 'driftvel'
c     
      INTEGER mode,ik,ir
      REAL    deltas1,deltas2,efield

      real    vpot
      real*8 :: rv(4,4),zv(4,4),phiv(4,4),drn,dzn,e_rad_tmp

      real ormin,ormax,ozmin,ozmax
      real phi1,phi2,ds1,ds2,fact
      integer irin,ikin,irout,ikout
      integer ikinu,irinu,ikoutu,iroutu
      integer ikind,irind,ikoutd,iroutd
      integer iku,ikd
      integer id,nc,nv,ic,iv

      
c
c     Initialization
c      
c     Number of cells and number of vertices
c
      nc = 4
      nv = 4
      osmpot2 = 0.0
      e_rad = 0.0
      e_pol = 0.0
      exb_rad_drft = 0.0
      exb_pol_drft = 0.0

c     
c     jdemod - depending on option chosen the sheath potential is either approximated as 3Te or can be loaded as
c     a profile from probe data
c
c     NOTE: The code defines the potential as ZERO at the target surface and thus the potential increases going along the ring ... the 
c           potential at the sheath edge is thus +3Te. 
c
c     
c     - assumes below that efield is zero in core - not sure what to do about core plasma potential since 
c     one would expect a potential discrepancy across the separatrix and pedestal. 
c     
c     
c     Perform calculation for all SOL rings including PFZ
c     

      DO ir = irsep, nrs
         IF (idring(ir).EQ.BOUNDARY) CYCLE

c     
c     jdemod - calculate from each target independently
c     
c     
c     jdemod - code to calculate the initial potential from LP data would go here
c     
!         write (6,'(a,i8,a)') '---------------------  ',ir,
!     >                        '  ---------------------' 

         osmpot2(0,ir) =  3.0 * kteds(idds(ir,2))
         vpot = osmpot2(0,ir)
c     
         DO ik = 1, ikmids(ir)
c     
c     jdemod - need to increment vpot for both halves of cell but the value 
c     assigned to the cell is the centerpoint
c     - could interpolate the efield along the field line but I don't think
c     it would change the integral between the two cell centers
c     
            deltas1 = kss(ik,ir) - ksb(ik-1 ,ir)
            deltas2 = ksb(ik,ir) - kss(ik,ir)

            efield = kes(ik,ir) / (qtim * qtim * emi / crmi)
c     
c     first half
c     
            vpot = vpot - efield * deltas1 
c     
c     assign
c     
            osmpot2(ik,ir) = vpot
c     
c     second half
c     
            vpot = vpot - efield * deltas2
c     
!            WRITE(6,'(a,2i8,10(1x,g18.8))') 'POT:',ik,ir,kss(ik,ir),
!     .           osmpot2(ik,ir),
!     .           kes(ik,ir)/(QTIM * QTIM * EMI / CRMI), 
!     .           efield,vpot
!

         end do
c     
c     jdemod - code to calculate the initial potential from LP data would go here
c     
         write (6,*) 

         osmpot2(nks(ir)+1,ir) = 3.0 * kteds(idds(ir,1))
         vpot = osmpot2(nks(ir)+1,ir)

         DO ik = nks(ir),ikmids(ir)+1,-1

c     
c     jdemod - need to increment vpot for both halves of cell but the value 
c     assigned to the cell is the centerpoint
c     - could interpolate the efield along the field line but I don't think
c     it would change the integral between the two cell centers
c     
            deltas1 = ksb(ik,ir) - kss(ik,ir)
            deltas2 = kss(ik,ir) - ksb(ik-1,ir)

            efield = kes(ik,ir) / (qtim * qtim * emi / crmi)

c     
c     first half
c     
c     Note: Efield has opposite sign on second half ring in order to get 
c           transport toward the target. Change the sign in the calculation
c           of potential.
c
            vpot = vpot + efield * deltas1 
c     
c     assign
c     
            osmpot2(ik,ir) = vpot
c     
c     second half
c
c     Note: Efield has opposite sign on second half ring in order to get 
c           transport toward the target. Change the sign in the calculation
c           of potential.
c     
            vpot = vpot + efield * deltas2
c     

!            WRITE(6,'(a,2i8,10(1x,g18.8))') 'POT:',ik,ir,kss(ik,ir),
!     .           osmpot2(ik,ir),
!     .           kes(ik,ir)/(QTIM * QTIM * EMI / CRMI), 
!     .           efield,vpot

         end do

      ENDDO

c     
c     jdemod - run through and calculate the potential at the cell corners
c     so that the derivatives can be calculated on a cell by cell
c     basis without referring to other cells
c     

c     
c     jdemod - or just calculate the gradients and resulting drifts now
c     obtaining the data for each corner might still be desirable(?)
c     
c     - now that osmpot2 is available ... can I use existing interpolation code? 
c     
c     These calculations need to go from one end of the ring to the other since the sign
c     of the gradients is important. 
c     - for radial drifts (a negative sign is an "inward" toward core or PFZ drift while a
c     positive value is toward the plasma edge.
c     - for poloidal drifts ... a negative sign is toward the "inner" target or target 1 while
c     a postive value is toward the "outer" target or target 2.
c     
c     
c     
c     1) Calculate grad_phi_poloidal
c     - take E-field interpolated to upper and lower cell boundaries
c     - calculate the distance between them using the formula
c     Need to calculate Btor/Bpol from KBFS which is Btot/Bpol
c     Btot = sqrt(Bpol**2 + Btor**2)
c     KBFS**2 = (Bpol**2 + Btor**2) / Bpol**2 = 1 + (Btor/Bpol)**2
c     Btor/Bpol = sqrt(kbfs**2 -1)
c     Note: Bpol = sqrt(Br**2 + Bz**2)
c     
c     fact = sqrt(kbfs(ik,ir)**2-1.0)
c     d_pol = deltas / fact = delta_s * Bpol/Btor
c     
c     e_pol = -grad_phi_poloidal = (Phi_top_cell - Phi_bot_cell) / d_pol
c     deltas = ksb(ik,ir)-ksb(ik-1,ir)
c     
c     
      do ir = irsep, nrs
         ! Do not calculate E fields for boundary rings
         if (idring(ir).eq.BOUNDARY) cycle

!         write (6,'(a,i8,a)') '---------------------  ',ir,
!     >                        '  ---------------------'          
         do ik = 1,nks(ir)
            if (ik.eq.1) then 
               phi1 = osmpot2(0,ir)
               ds1 = kss(ik,ir)-ksb(0,ir)
               
               ds2 = ksb(ik,ir)-kss(ik,ir)
               phi2 = osmpot2(ik,ir) + 
     >              (osmpot2(ik+1,ir)-osmpot2(ik,ir)) * 
     >              ds2/(kss(ik+1,ir)-kss(ik,ir)) 


            elseif (ik.eq.nks(ir)) then 

               ds1 =  kss(ik,ir)-ksb(ik-1,ir)
               phi1 = osmpot2(ik,ir) +
     >              (osmpot2(ik-1,ir)-osmpot2(ik,ir)) * 
     >              ds1/(kss(ik,ir)-kss(ik-1,ir)) 

               ds2 = ksb(ik,ir)-kss(ik,ir)
               phi2 = osmpot2(nks(ir)+1,ir)

            else
               ds1 =  kss(ik,ir)-ksb(ik-1,ir)
               phi1 = osmpot2(ik,ir) +
     >              (osmpot2(ik-1,ir)-osmpot2(ik,ir)) * 
     >              ds1/(kss(ik,ir)-kss(ik-1,ir)) 
               
               ds2 = ksb(ik,ir)-kss(ik,ir)
               phi2 = osmpot2(ik,ir) + 
     >              (osmpot2(ik+1,ir)-osmpot2(ik,ir)) * 
     >              ds2/(kss(ik+1,ir)-kss(ik,ir)) 
            endif
            
            fact = sqrt(kbfs(ik,ir)**2-1.0)
            if (fact.le.0.0) then 
               e_pol(ik,ir) = 0.0
            else
               e_pol(ik,ir) = -(phi2-phi1)/(ds1+ds2) * fact
            endif

!            write(6,'(a,2i8,20(1x,g18.8))') 'E_POL:',ik,ir,e_pol(ik,ir),
!     >           phi2,phi1,ds2,ds1,(ds2+ds1)/fact,fact,ksb(ik-1,ir),
!     >           kss(ik,ir),ksb(ik,ir),osmpot2(ik,ir)

         enddo

         ! Zero out midpoint where potentials do not match
         e_pol(ikmids(ir),ir) = 0.0
         e_pol(ikmids(ir)+1,ir) = 0.0

      enddo

c     
c     
c     2) Calculate e_rad = -grad_phi_radial
c     - This is more complicated
c     - extract 9 phi points centered on the current cell
c     - Load the normalize vector for the cell center line [drn,dzn]
c     - Take the normal vector to this (90 degrees CCW ... [-dzn,drn] ... this defines the direction for interpolation of
c     grad_phi_radial since the normal vector defines the radial direction to the field line at the cell center.
c     - Go a small distance along the normal in each direction ... interpolate to get Phi at these locations ... calculate 
c     the gradient ... grad_phi_radial = (ph1-ph2)/(dx1+dx2) where dx1 and dx2 are the distances from the cell center
c     along the normal line to the phi interpolation locations. 
c     - the one requirement is that the interpolation locations must be within the region defined by the 8 adjacent cell centers
c     no matter what the direction of the cell center normal vector. This means that the distance needs to be resampled
c     until the test locations lie within the defined cells for interpolation (recursive code?).
c     - also need to test for degenerate cases or boundary rings where there could be geometry issues. 
c     - alternatively, could implement a routine to find the cell wall intersections and do the interpolation to that location,
c     however, in this case the contributions from the other corners of the cell would be left out so it might not be as 
c     accurate (?) ... should probably use short distances relative to the cells size and rely on the interpolation routines
c     

      do ir = irsep,nrs
         ! Do not calculate E fields for boundary rings
         if (idring(ir).eq.BOUNDARY) cycle

         do ik = 1,nks(ir)
            iku = ik+1
            ikd = ik-1
            if (ik.eq.1) then 
               ikd = ik
            elseif (ik.eq.nks(ir)) then 
               iku = ik
            endif
c
c           center cell
c
            irin = irins(ik,ir)
            ikin = ikins(ik,ir)
            irout = irouts(ik,ir)
            ikout = ikouts(ik,ir)
c
c           Next cell up ring  
c
            irinu = irins(iku,ir)
            ikinu = ikins(iku,ir)
            iroutu = irouts(iku,ir)
            ikoutu = ikouts(iku,ir)
c
c           Last cell down ring
c
            irind = irins(ikd,ir)
            ikind = ikins(ikd,ir)
            iroutd = irouts(ikd,ir)
            ikoutd = ikouts(ikd,ir)
c
            call get_cell_norm(ik,ir,drn,dzn)

c
c           Deal with cells at targets since rs and zs are
c           replaced by krb and kzb (cell bounds for first and
c           last cells). This code should not need modificaiton
c           for dealing with cells adjacent to core plasma
c
            if (ik.eq.1) then 

c     
c     Load up an array containing the corner points and potential values for 
c     each quadrant from the cell center
c     

! first polygon
               rv(1,1) = krb(ik-1,ir)
               zv(1,1) = kzb(ik-1,ir)
               phiv(1,1) = osmpot2(ik-1,ir)
               
               rv(2,1) = krb(ikout-1,irout)
               zv(2,1) = kzb(ikout-1,irout)
               phiv(2,1) = osmpot2(ikout-1,irout)
               
               rv(3,1) = rs(ikout,irout)
               zv(3,1) = zs(ikout,irout)
               phiv(3,1) = osmpot2(ikout,irout)

               rv(4,1) = rs(ik,ir)
               zv(4,1) = zs(ik,ir)
               phiv(4,1) = osmpot2(ik,ir)

! second polygon

               rv(1,2) = rs(ik,ir)
               zv(1,2) = zs(ik,ir)
               phiv(1,2) = osmpot2(ik,ir)

               rv(2,2) = rs(ikout,irout)
               zv(2,2) = zs(ikout,irout)
               phiv(2,2) = osmpot2(ikout,irout)

               rv(3,2) = rs(ikout+1,irout)
               zv(3,2) = zs(ikout+1,irout)
               phiv(3,2) = osmpot2(ikout+1,irout)

               rv(4,2) = rs(ik+1,ir)
               zv(4,2) = zs(ik+1,ir)
               phiv(4,2) = osmpot2(ik+1,ir)

! third polygon

               rv(1,3) = rs(ikin,irin)
               zv(1,3) = zs(ikin,irin)
               phiv(1,3) = osmpot2(ikin,irin)

               rv(2,3) = rs(ik,ir)
               zv(2,3) = zs(ik,ir)
               phiv(2,3) = osmpot2(ik,ir)

               rv(3,3) = rs(ik+1,ir)
               zv(3,3) = zs(ik+1,ir)
               phiv(3,3) = osmpot2(ik+1,ir)

               rv(4,3) = rs(ikin+1,irin)
               zv(4,3) = zs(ikin+1,irin)
               phiv(4,3) = osmpot2(ikin+1,irin)

! fourth polygon

               rv(1,4) = krb(ikin-1,irin)
               zv(1,4) = kzb(ikin-1,irin)
               phiv(1,4) = osmpot2(ikin-1,irin)
               
               rv(2,4) = krb(ik-1,ir)
               zv(2,4) = kzb(ik-1,ir)
               phiv(2,4) = osmpot2(ik-1,ir)
               
               rv(3,4) = rs(ik,ir)
               zv(3,4) = zs(ik,ir)
               phiv(3,4) = osmpot2(ik,ir)

               rv(4,4) = rs(ikin,irin)
               zv(4,4) = zs(ikin,irin)
               phiv(4,4) = osmpot2(ikin,irin)
c     
            elseif (ik.eq.nks(ir)) then 
c     
c     Load up an array containing the corner points and potential values for 
c     each quadrant from the cell center
c     

! first polygon
               rv(1,1) = rs(ik-1,ir)
               zv(1,1) = zs(ik-1,ir)
               phiv(1,1) = osmpot2(ik-1,ir)
               
               rv(2,1) = rs(ikout-1,irout)
               zv(2,1) = zs(ikout-1,irout)
               phiv(2,1) = osmpot2(ikout-1,irout)
               
               rv(3,1) = rs(ikout,irout)
               zv(3,1) = zs(ikout,irout)
               phiv(3,1) = osmpot2(ikout,irout)

               rv(4,1) = rs(ik,ir)
               zv(4,1) = zs(ik,ir)
               phiv(4,1) = osmpot2(ik,ir)

! second polygon

               rv(1,2) = rs(ik,ir)
               zv(1,2) = zs(ik,ir)
               phiv(1,2) = osmpot2(ik,ir)

               rv(2,2) = rs(ikout,irout)
               zv(2,2) = zs(ikout,irout)
               phiv(2,2) = osmpot2(ikout,irout)

               rv(3,2) = krb(ikout,irout)
               zv(3,2) = kzb(ikout,irout)
               phiv(3,2) = osmpot2(ikout+1,irout)

               rv(4,2) = krb(ik,ir)
               zv(4,2) = kzb(ik,ir)
               phiv(4,2) = osmpot2(ik+1,ir)

! third polygon

               rv(1,3) = rs(ikin,irin)
               zv(1,3) = zs(ikin,irin)
               phiv(1,3) = osmpot2(ikin,irin)

               rv(2,3) = rs(ik,ir)
               zv(2,3) = zs(ik,ir)
               phiv(2,3) = osmpot2(ik,ir)

               rv(3,3) = krb(ik,ir)
               zv(3,3) = kzb(ik,ir)
               phiv(3,3) = osmpot2(ik+1,ir)

               rv(4,3) = krb(ikin,irin)
               zv(4,3) = kzb(ikin,irin)
               phiv(4,3) = osmpot2(ikin+1,irin)

! fourth polygon

               rv(1,4) = rs(ikin-1,irin)
               zv(1,4) = zs(ikin-1,irin)
               phiv(1,4) = osmpot2(ikin-1,irin)
               
               rv(2,4) = rs(ik-1,ir)
               zv(2,4) = zs(ik-1,ir)
               phiv(2,4) = osmpot2(ik-1,ir)
               
               rv(3,4) = rs(ik,ir)
               zv(3,4) = zs(ik,ir)
               phiv(3,4) = osmpot2(ik,ir)

               rv(4,4) = rs(ikin,irin)
               zv(4,4) = zs(ikin,irin)
               phiv(4,4) = osmpot2(ikin,irin)


            else

               ! deal with standard cells on the grid
               ! first polygon - down and out :) 
               rv(1,1) = rs(ikd,ir)
               zv(1,1) = zs(ikd,ir)
               phiv(1,1) = osmpot2(ikd,ir)
               
               rv(2,1) = rs(ikoutd,iroutd)
               zv(2,1) = zs(ikoutd,iroutd)
               phiv(2,1) = osmpot2(ikoutd,iroutd)
               
               rv(3,1) = rs(ikout,irout)
               zv(3,1) = zs(ikout,irout)
               phiv(3,1) = osmpot2(ikout,irout)

               rv(4,1) = rs(ik,ir)
               zv(4,1) = zs(ik,ir)
               phiv(4,1) = osmpot2(ik,ir)

               ! second polygon
               ! up and out

               rv(1,2) = rs(ik,ir)
               zv(1,2) = zs(ik,ir)
               phiv(1,2) = osmpot2(ik,ir)

               rv(2,2) = rs(ikout,irout)
               zv(2,2) = zs(ikout,irout)
               phiv(2,2) = osmpot2(ikout,irout)

               rv(3,2) = rs(ikoutu,iroutu)
               zv(3,2) = zs(ikoutu,iroutu)
               phiv(3,2) = osmpot2(ikoutu,iroutu)

               rv(4,2) = rs(iku,ir)
               zv(4,2) = zs(iku,ir)
               phiv(4,2) = osmpot2(iku,ir)

               ! third polygon
               ! up and in 

               rv(1,3) = rs(ikin,irin)
               zv(1,3) = zs(ikin,irin)
               phiv(1,3) = osmpot2(ikin,irin)

               rv(2,3) = rs(ik,ir)
               zv(2,3) = zs(ik,ir)
               phiv(2,3) = osmpot2(ik,ir)

               rv(3,3) = rs(iku,ir)
               zv(3,3) = zs(iku,ir)
               phiv(3,3) = osmpot2(iku,ir)

               rv(4,3) = rs(ikinu,irinu)
               zv(4,3) = zs(ikinu,irinu)
               phiv(4,3) = osmpot2(ikinu,irinu)

               ! fourth polygon
               ! down and in

               rv(1,4) = rs(ikind,irind)
               zv(1,4) = zs(ikind,irind)
               phiv(1,4) = osmpot2(ikind,irind)
               
               rv(2,4) = rs(ikd,ir)
               zv(2,4) = zs(ikd,ir)
               phiv(2,4) = osmpot2(ikd,ir)
               
               rv(3,4) = rs(ik,ir)
               zv(3,4) = zs(ik,ir)
               phiv(3,4) = osmpot2(ik,ir)

               rv(4,4) = rs(ikin,irin)
               zv(4,4) = zs(ikin,irin)
               phiv(4,4) = osmpot2(ikin,irin)

            endif

            ! make adjustments to phi values for cells that are adjacent to the core
            ! Need to do more to detect when adjacent to the Xpoint. The cells that 
            ! are just below and just above the Xpoint need adjustment since the
            ! code tries to use the center points of inappropriate non-adjacent cells. 
            ! 
            ! Adjust osmpot2 values if center cell is adjacent to core
            !
            if (irin.lt.irsep) then 
               ! Need to change the phi values for the "core" elements since they will be zero
               ! at this point

               ! third polygon
               ! Set potentials equal to separatrix value

               phiv(1,3) = osmpot2(ik,ir)

               ! fourth polygon
               ! Set potentials equal to separatrix value

               phiv(4,4) = osmpot2(ik,ir)

            endif

            ! Adjust osmpot2 values if up cell is adjacent to core
            if (irinu.lt.irsep) then
               phiv(4,3) = osmpot2(iku,ir)
            endif

            ! Adjust osmpot2 values if down cell is adjacent to core
            if (irind.lt.irsep) then 
               phiv(1,4) = osmpot2(ikd,ir)               
            endif
            
            ! make adjustments for phi values at the grid edge ... 
            ! irout = irwall 
            ! irin  = irtrap 
            !
            ! idring(ir) = BOUNDARY 
            ! Assign the same phi value as cell center when at a boundary 
            ! ... which gives a zero gradient for boundary side of cell. 

            if (idring(irout).eq.BOUNDARY) then 
               phiv(3,1) = phiv(4,1)
               phiv(2,2) = phiv(1,2)
            endif
            
            if (idring(iroutu).eq.BOUNDARY) then 
               phiv(3,2) = phiv(4,2)
            endif

            if (idring(iroutd).eq.BOUNDARY) then 
               phiv(2,1) = phiv(1,1)

            endif

            if (idring(irin).eq.BOUNDARY) then 
               phiv(1,3) = phiv(2,3)
               phiv(4,4) = phiv(3,4)
            endif
            
            if (idring(irinu).eq.BOUNDARY) then 
               phiv(4,3) = phiv(3,3)
            endif


            if (idring(irind).eq.BOUNDARY) then 
               phiv(1,4) = phiv(2,4)
            endif



!               do ic = 1,nc
!                  write(6,'(a,7i8,20(1x,g15.8))') 'CELLS:',
!     >                 ic,ik,ir,ikin,irin,ikout,irout,
!     >                 (rv(iv,ic),zv(iv,ic),phiv(iv,ic),iv=1,nv),
!     >                 rs(ik,ir),zs(ik,ir),drn,dzn,sqrt(drn**2+dzn**2)
!               end do


            call get_e_rad(rv,zv,phiv,drn,dzn,nc,nv,ik,ir,e_rad_tmp,
     >           dble(rs(ik,ir)),dble(zs(ik,ir)))

            e_rad(ik,ir) = e_rad_tmp

!            write(6,'(a,2i8,10(1x,g15.8))') 'E_RAD:',ik,ir,e_rad(ik,ir),
!     >           drn,dzn,rs(ik,ir),zs(ik,ir)
         end do
         ! Zero out midpoint where potentials do not match
         e_rad(ikmids(ir),ir) = 0.0
         e_rad(ikmids(ir)+1,ir) = 0.0

      end do

c     
c     Now that the radial and poloidal electric fields have been calculated
c     Calculate the drifts ... scale the drifts to distances by 
c     multiplying by the ion time step and applying geometric factors
c     to the poloidal drift
c     
      do ir = irsep,nrs
         do ik = 1,nks(ir)
            if (bts(ik,ir).ne.0.0) then 
! include timestep in radial drift

! project the poloidal drift to S and include timestep
               fact = sqrt(kbfs(ik,ir)**2-1.0)
               exb_pol_drft(ik,ir)= e_rad(ik,ir)/bts(ik,ir) 
     >                                 * fact *qtim * exb_scale
               exb_rad_drft(ik,ir) = -e_pol(ik,ir)/bts(ik,ir)*qtim 
     >                                 * exb_scale
            endif
         end do

      end do
c
c     If the plasma overlay option is turned on then 
c     MASK out all of the results calculated except for values
c     inside the overlay region. 
c
      if (external_plasma_overlay.eq.1) then 
         call get_overlay_limits(ormin,ormax,ozmin,ozmax)

         do ir = 1,nrs
            do ik = 1,nks(ir)
               if (rs(ik,ir).lt.ormin.or.rs(ik,ir).gt.ormax.or.
     >             zs(ik,ir).lt.ozmin.or.zs(ik,ir).gt.ozmax) then 
                  osmpot2(ik,ir) = 0.0
                  e_rad(ik,ir) = 0.0
                  e_pol(ik,ir) = 0.0
                  exb_pol_drft(ik,ir) = 0.0
                  exb_rad_drft(ik,ir) = 0.0
               end if
            end do
         end do
      endif 
c     
c     Lets print out everything 
c     
      if (cprint.ge.9) then 
         write(6,'(a)') 'DRIFT INFORMATION:',qtim
         do ir = irsep,nrs
            do ik = 1,nks(ir)
               fact = sqrt(kbfs(ik,ir)**2-1.0)

               write(6,'(2i8,20(1x,g14.6))') ik,ir,
     >              osmpot2(ik,ir),e_rad(ik,ir),e_pol(ik,ir),
     >              exb_rad_drft(ik,ir),exb_pol_drft(ik,ir),bts(ik,ir),
     >              fact

            end do
         end do
      endif


      RETURN
 99   STOP
      END
c
c
c
      subroutine get_e_rad(rv,zv,phi,drn,dzn,nc,nv,ik,ir,e_rad_tmp,r,z)
      use mod_interpolate
      implicit none
      integer :: nv,nc,ik,ir
      real*8 :: rv(nv,nc),zv(nv,nc),phi(nv,nc),drn,dzn,r,z,e_rad_tmp

      real :: areas(4)

      real*8 :: rt(2),zt(2),dt(2)
      integer :: ct(2),ic
      real*8 :: dist
      real*8 :: phi1,phi2

      ! calculate the radial electric field for the cell

      e_rad_tmp = 0.0

      ! calculate areas of the cells to detect boundary conditions
      ! If any of the cell areas are zero then set e_rad = 0.0 since 
      ! we are on a boundary ring or grid edge
       do ic = 1,nc
          areas(ic) = cell_area(rv(:,ic),zv(:,ic),nv)
          if (areas(ic) .le. 0.0) then
             write(6,'(a,2i8)') 'Debug: Zero area cells in e_rad'//
     >                  ' calculation - boundary?', ik,ir
             e_rad_tmp = 0.0
             return
          endif
       end do

       ! Determine test points ... they must lie along normal and be in one of the 4 cells
       ! around the point where we are determining the value of -grad_phi_radial
       
       ! Note that the vertices 2-3 and 1-4 are parallel to the field lines in the 
       ! cells passed into this routine. Sides 1-2 and 3-4 go across the field lines
       ! but may be at almost any angle to the field lines. Cells are non-orthogonal in general.

       ! Perform calculations in double precision. 

       ! Start off using +/- 0.00001m or 0.01mm from the cell center - this should pretty much always 
       ! lie in the cell. 
       ! dist = 1.0d-3

       call get_test_points(rv,zv,r,z,drn,dzn,rt,zt,ct,dt,nv,nc)

       ! now that we have test points that are guaranteed to be within cells and we have identified
       ! the cells and distances from the cell center ... calculate the gradient. 
       ! Note: cell vertices are listed clockwise .. sides 1->2, 3->4 are across the field lines while
       !       sides 2->3 and 1->4 are parallel to the field lines ... this must match the convention
       !        used in the interpolation routine. 
       ! Need to define the sign used. Test point 1 will be along the "positive" direction of the normal
       ! or "outward" from the cell center while test point 2 is "inward" or along the negative normal direction. 
       ! e = -grad(phi) ... for now lets use phi2 = phi(test_point_2) ... phi1 = phi(test_point_1)
       ! e = -grad(phi) = -(phi2-phi1)/(dist1+dist2)
       ! NOTE: if you consider the positive R axis as the reference line in the outer divertor ... then 
       !       outward radially in the outer divertor the definition would be ... 
       !       E-rad = -(phi_outer-phi_inner)/dRad = - (phi1 - phi2)/(dist1+dist2)
       !       Switch to this definition for now. 

       call cell_interpolate(rt(1),zt(1),phi1,
     >                       rv(:,ct(1)),zv(:,ct(1)),phi(:,ct(1)))
       
       call cell_interpolate(rt(2),zt(2),phi2,
     >                       rv(:,ct(2)),zv(:,ct(2)),phi(:,ct(2)))
       
       dist = dt(1) + dt(2)

       if (dist.le.0.0) then 
          write(0,*) 'ERROR in get_e_rad: dist = 0.0'
          e_rad_tmp = 0.0
       else
!
!         Set definiton of E-rad = -(phi1-phi2)/dist
!         
!          e_rad_tmp = -(phi2-phi1)/dist
!
          e_rad_tmp = -(phi1-phi2)/dist
       endif

!       write(6,'(a,16x,4(1x,g15.8),i8,4(1x,g15.8),i8,g15.8)') 'E_RAD:',
!     >      phi1,rt(1),zt(1),dt(1),ct(1),
!     >      phi2,rt(2),zt(2),dt(2),ct(2),
!     >      e_rad_tmp


       return
       end


      subroutine get_test_points(rv,zv,r,z,drn,dzn,rt,zt,ct,dt,nv,nc)
      implicit none
      integer nv,nc
      real*8 :: rv(nv,nc),zv(nv,nc),drn,dzn,r,z
      real*8 :: rt(2),zt(2),dt(2)
      integer :: ct(2),found
      logical :: res            

      real*8 :: dist
      integer :: npt,ic
      logical,external :: inpolydp
      integer :: iv

      found = 0
      !dist = 0.001d0            ! start with +/- 1mm from center
      !dist = 0.0001d0            ! start with +/- 0.1mm from center
      dist = 0.00001d0            ! start with +/- 0.01mm from center
      ct = 0.0
      dt = dist
!      write(6,'(a,20(1x,g14.7))') 'TESTA:',r,z,drn,dzn


      do while (found.ne.2)
         !     need to get two points found within the cells
         do npt = 1,2
            rt(npt) = r - (-1)**npt * drn * dist
            zt(npt) = z - (-1)**npt * dzn * dist
!            write(6,'(a,2i8,20(1x,g14.7))') 'TESTB:',npt,found,
!     >           r,z, rt(npt),zt(npt),drn,dzn,dist,(-1)**npt

            do ic = 1,nc
               res = inpolydp(rt(npt),zt(npt),nv,rv(:,ic),zv(:,ic)) 
               if (res) then
                  !     point has been found in polygon
                  found = found + 1
                  ct(npt) = ic
                  dt(npt) = dist                   
!                  write(6,'(a,2i8,20(1x,g14.7))') 'TESTC:',npt,found,
!     >                 r,z, rt(npt),zt(npt),ct(npt),
!     >                 drn,dzn,dist,(-1)**npt
                  exit
               endif
            end do
         end do

         if (found .ne.2) then
            found = 0
            ! reduce dist and try again
            ! check to see if dist is too small and issue an error message
            dist = dist * 0.5d0
            if (dist.lt.1.0d-6) then 
               write(0,*) 'ERROR finding test points'//
     >              ' in e_rad calculation: DIST TOO SMALL =',dist
               write(6,*) 'ERROR finding test points'//
     >              ' in e_rad calculation: DIST TOO SMALL =',dist
               do ic = 1,nc
                  write(0,'(a,i8,20(1x,g14.6))') 'GET_TEST_POINTS:',
     >                 ic,(rv(iv,ic),zv(iv,ic),iv=1,nv),r,z,dist,
     >                 rt(1),zt(1),rt(2),zt(2)
                  write(6,'(a,i8,20(1x,g14.6))') 'GET_TEST_POINTS:',
     >                 ic,(rv(iv,ic),zv(iv,ic),iv=1,nv),r,z,dist,
     >                 rt(1),zt(1),rt(2),zt(2)
               end do
               stop 'ERROR in GET_TEST_POINTS'
            endif
         endif
      end do

      return 
      end


      subroutine get_cell_norm(ik,ir,drn,dzn)
      implicit none
      include 'params'
      include 'cgeom'
      integer ik,ir
      real*8 :: dn
      real*8 ::  drn,dzn
      real*8 :: drt, dzt

      ! Find the unit vector normal to the axis of the cell

      drt = krb(ik,ir)-krb(ik-1,ir)
      dzt = kzb(ik,ir)-kzb(ik-1,ir)

      dn = sqrt(drt**2+dzt**2)


      ! Outward normal unit vector from cell center

      drn = -dzt/dn
      dzn =  drt/dn

!      write(6,'(a,2i4,20(1x,g14.7))') 'Norm:',ik,ir,krb(ik,ir),
!     >       kzb(ik,ir),krb(ik-1,ir),kzb(ik-1,ir),drn,dzn,drt,dzt,dn


      return
      end

c
c ======================================================================
c
c
c ======================================================================
c
      SUBROUTINE LoadTriangles
      USE mod_eirene04
      IMPLICIT none

      INTEGER fp,i1,i2
      REAL    version

      fp = 99
      OPEN(UNIT=fp,FILE='triangles.raw',ACCESS='SEQUENTIAL',
     .     FORM='UNFORMATTED',STATUS='OLD',ERR=98)            
      READ(fp,ERR=98) version,ntri,nver,nadd
c      WRITE(0,*) 'TRI:',version,ntri,nver,nadd

      IF (version.NE.1.0)
     .  CALL ER('LoadTriangles','Unsupporting version',*99)

      CALL ALLOC_VERTEX(nver)
      CALL ALLOC_SURFACE(nadd)
      CALL ALLOC_TRIANGLE(ntri)
      READ(fp,ERR=98) (tri(i1),i1=1,ntri)
      READ(fp,ERR=98) ((ver(i1,i2),i2=1,3),i1=1,nver)
      READ(fp,ERR=98) (add(i1),i1=1,nadd)

c      READ(fp,ERR=98) tri,ver,add
      CLOSE (fp)
      
      RETURN
 98   CALL ER('LoadTriangles','Problems reading data file',*99)
 99   STOP
      END

