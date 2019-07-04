c     -*-Fortran-*-
c
c ======================================================================
c
c ray: program init
c
c
c
      PROGRAM RAY
      IMPLICIT none

c...  Input:
      INTEGER iopt,i,j,n
      CHARACTER buffer*1024,cdum1*1024

c...  Postscript initialization:      
      CALL GPSTOP (100)
      CALL PAPER  (1)
c
c     System dependent printer initialization
c
c      call printerinit
c
c      CALL XUFLOW (0)
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      CALL FULL
      CALL LINCOL(1)	
      CALL CTRMAG(10)
      CALL PLOTST(0.0,0.0,'[i]')

      OPEN(5,FILE='ray.input',FORM='FORMATTED',STATUS='OLD',ERR=97)     

      DO WHILE(.TRUE.) 

        READ(5,'(A1024)',ERR=98,END=10) buffer 

c        WRITE(0,*) 'buffer:',buffer(1:50)

c...    Legacy input (OUT):
        SELECTCASE (buffer(2:4))        
          CASE('985')
            READ(buffer,*) cdum1,iopt
            IF (iopt.NE.0) CALL Main985(iopt,'no title')
          CASE('989')
            READ(buffer,*) cdum1,iopt
            IF (iopt.NE.0) CALL Main989(iopt)
        ENDSELECT

c...    Isolate tag string:
        n = LEN_TRIM(buffer)
        DO i = 1, n
          IF (buffer(i:i).EQ.'{') EXIT
        ENDDO
        DO j = i+1, n
          IF (buffer(j:j).EQ.'}') EXIT
        ENDDO
c...    Process proper tagged input:
        IF (i.LT.n.AND.j.LE.n) THEN
          SELECTCASE (buffer(i+1:j-1))        
            CASE('RAY TRACE')
           write(0,*) 'buffer '//buffer
              READ(buffer,*) cdum1,iopt
              IF (iopt.NE.0) CALL Main985(iopt)
            CASE('INVERSION')
              READ(buffer,*) cdum1,iopt
              IF (iopt.NE.0) CALL Main989(iopt)
          ENDSELECT
        ENDIF

      ENDDO
 10   CONTINUE

      CLOSE(5)

      STOP
 97   CALL ER('Ray','Unable to open ray.input',*99)
 98   CALL ER('Ray','Error reading input file',*99)
 99   STOP
      END
c
c ======================================================================
c
