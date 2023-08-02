c
c     The following routine reads in the unstructured input (aka optional
c     input) in the OEDGE/DIVIMP input file. The optional inputs are identified
c     by a tag in the input file. These tags are a letter followed by a 2-digit
c     number - leaving room for up to 2600 options in the current configuration
c     These options are usually arranged into letter groupings that are related
c     to a common theme. In addition, "tag series" which have a large number of
c     options are read in by tag specific subroutines:
c
c     e.g.
c
c        SUBROUTINE ReadTagSeries_G
c        SUBROUTINE ReadTagSeries_H  
c        SUBROUTINE ReadTagSeries_I  
c
c     In addition to the letter indexed tags the original implementation used
c     three digit numbers as the tags. These options are somewhat less well
c     documented. 
c
c     Finally, the OEDGE/DIVIMP input file has a tag assigned for every input
c     value whether optional or not. This indexing will allow for all of the input
c     to be made optional when this is desirable. Any additional optional input
c     must be assigned a new and unique tag.  
c
c
      SUBROUTINE ReadUnstructuredInput(line2)
      USE mod_osm_input
      use allocatable_input_data
      use mod_sol22_input
      use ero_interface
      use mod_params
      use mod_slcom
      use mod_cadas
      use mod_comtor
      use mod_cgeom
      use mod_cedge2d
      use mod_solparams
      use mod_solswitch
      !use mod_solcommon
      use mod_reiser_com
      use mod_line_profile
      use mod_out_unstruc
      use mod_driftvel
      use mod_diagvel
      use mod_dperpz
      use mod_lambda
      use mod_sol29_input
      use mod_io
      IMPLICIT none

c
c     jdemod - cleaning up some of the line length constraints since the inputs really don't need them
c
c      CHARACTER line2*(*),LINE*72,TAG*3,coment*72,cdum1*1024
      CHARACTER line2*(*),LINE*128,TAG*3,cdum1*1024
      REAL      R,vol,z1,version
      INTEGER   I,ir,ierr,i1,i2

      REAL GetInputVersion
c
c jdemod - added variable to hold line read when calling RDG1 to get 
c          ADAS data.
c
      character line3*512
c
c     INCLUDE 'params'
c     INCLUDE 'slcom'
c     include 'cadas'
c     include 'comtor' 
c     include 'cgeom'  
c     include 'cedge2d'

c     INCLUDE 'solparams'
c     INCLUDE 'solswitch'
c     INCLUDE 'solcommon'

c     include 'reiser_com'
c     include 'line_profile'
c
c     include 'out_unstruc'
c     include 'driftvel'
c     include 'dperpz'  
c
c      include 'cyield'
c
      INTEGER tagnum,fp
      DATA    tagnum /0/

      COMMON /OPTTEMP/ osm_matcht,forcet1
      INTEGER          osm_matcht,forcet1

      COMMON /INPUTCHK/ inputflag
      INTEGER           inputflag(100)

      COMMON /MACHCOM/ machine2
      CHARACTER*64     machine2


      INTEGER    MAXTAG
      PARAMETER (MAXTAG=1000)
      COMMON /INPCHK/ ntaglist,taglist
      INTEGER     ntaglist,idum1
      REAL        rdum1
      CHARACTER*3 taglist(MAXTAG)

      INTEGER    itag
      LOGICAL    status,sol28_first
      CHARACTER  local_buffer*1024

      DATA sol28_first /.TRUE./
      SAVE
c
c     Function declaration for TAG T29 
c
c      real vtest,res,vr_pdf_int
c      external vr_pdf_int
      integer in
      real deltav1, deltav2
c
      WRITE(line,'(A128)') line2

      WRITE(TAG,'(A3)') LINE(3:5)

      ierr = 0

      fp = PINOUT

      ntaglist = ntaglist + 1
      taglist(ntaglist) = tag


      WRITE(SLOUT,*) 'TAG:',tag

c
c     *** NOTE ***
c
c...  This routine is getting unruly so it is being divided into 
c     a series of routines that process the various input categories
c     individually:
c -----------------------------------------------------------------------
      IF     (line(2:2).EQ.'{') THEN
        IF (sol28_first) THEN 
c          WRITE(0,*) 'INITIALIZING SOL28 OPTIONS'
c          CALL InitializeOptions
          sol28_first = .FALSE.
        ENDIF
c...    SOL28/OSM input options:
        IF (line(3:5).EQ.'999'.OR.line(3:6).EQ.'EXIT') THEN   ! Not great...
          CALL ProcessIterationBlocks
        ELSE
          s28mode = 4.0
          WRITE(local_buffer,'(A)') line2
c          WRITE(0,*) 'local_buffer:',local_buffer(1:100)
c...      Isolate tag string:
          DO itag = 2, LEN_TRIM(local_buffer)
            IF (local_buffer(itag:itag).EQ.'}') EXIT
          ENDDO
c...      Remove the portion of the data line that comes after a comment
c         character:
          DO i1 = 1, LEN_TRIM(local_buffer)
             IF (local_buffer(i1:i1).EQ.'$'.OR.
     >           local_buffer(i1:i1).EQ.'*') EXIT
          ENDDO 
          local_buffer(i1:LEN_TRIM(local_buffer)) = ' '
          CALL ProcessInputTag(5,itag,local_buffer,status)
        ENDIF
        ntaglist = ntaglist - 1
C
C Series G
C
      ELSEIF (tag(1:1).EQ.'G') THEN
        CALL ReadTagSeries_G(tag,line,fp)
c
c Series H
c
      elseif (tag(1:1).eq.'H') then 
        CALL ReadTagSeries_H(line,tag,fp)
c
c Series I
c
      ELSEIF (tag(1:1).EQ.'I') THEN
        CALL ReadTagSeries_I(line,tag,fp)

c
c Series K
c     - this tag series is assigned to ERO interface related 
c       quantities
c
c
      ELSEIF (tag(1:1).EQ.'K') THEN
        CALL read_ero_unstructured_input(line,tag,fp)
c
c Series O
c     - this tag refers to unstructured input related to OUT
c     - OUT shares this code with DIVIMP but the initialization routines 
c       are separate. This does cause some overhead storing variables not 
c       used in the specific codes but for now this is only 3 reals in 
c       OUT or about 12 bytes ... revisit later if it becomes an issue
c
      elseif (tag(1:1).eq.'O') then 
         call read_out_unstructured_input(line,tag,fp)

c
c
c
c
c
c -----------------------------------------------------------------------
c      ELSEIF (TAG(1:3).EQ.'001') THEN
c -----------------------------------------------------------------------
c      ELSEIF (TAG(1:3).EQ.'002') THEN
c -----------------------------------------------------------------------
c      ELSEIF (TAG(1:3).EQ.'003') THEN
c -----------------------------------------------------------------------
c      ELSEIF (TAG(1:3).EQ.'004') THEN
c -----------------------------------------------------------------------
c      ELSEIF (TAG(1:3).EQ.'005') THEN
c -----------------------------------------------------------------------
c      ELSEIF (TAG(1:3).EQ.'006') THEN
c -----------------------------------------------------------------------
c      ELSEIF (TAG(1:3).EQ.'007') THEN
c -----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'008') THEN
        CALL ReadI(line,outmode,0,3,'Output mode for user information')
c -----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'009') THEN
        CALL ReadI(line,cgridst,0,1,'Stuff grid')
c -----------------------------------------------------------------------
c
c     EIRENE related options:
c
      ELSEIF (TAG(1:3).EQ.'010') THEN
        CALL ReadI(line,eirgeom,0,3,'EIRENE geometry')
      ELSEIF (TAG(1:3).EQ.'011') THEN
        CALL ReadI(line,eirgrid,0,1,'EIRENE grid type')
      ELSEIF (TAG(1:3).EQ.'014') THEN
        CALL ReadI(line,eiradd,0,100,'EIRENE AddUsr option')
      ELSEIF (TAG(1:3).EQ.'018') THEN
        CALL ReadI(line,eirneut,0,1,'EIRENE neutral data option')
      ELSEIF (TAG(1:3).EQ.'019') THEN
        CALL ReadI(line,eirdebug,-100,10000,'EIRENE debugging option')
      ELSEIF (TAG(1:3).EQ.'020') THEN
        CALL ReadI(line,eirtime,0,600000,'EIRENE execution time')

        CALL GetEnv('DIVNAME',machine2)
c        WRITE(0,*) 'MARK: MACHINE2= '//machine2(1:LEN_TRIM(machine2))
c
C       IPP/01 - Krieger: SUN Fortran chokes if machine has zero chars
c
        IF (.FALSE..AND.len_trim(machine2).gt.0) THEN
          IF     (machine2(1:LEN_TRIM(machine2)).EQ.'hannah') THEN
            CALL WN('GetInput','Increasing EIRENE runtime for HANNAH')
            eirtime = eirtime * 1.1
          ELSEIF (machine2(1:LEN_TRIM(machine2)).EQ.'jonah') THEN
c            CALL WN('GetInput','Decreasing EIRENE runtime for JONAH')
c            eirtime = eirtime * 0.25
          ELSEIF (machine2(1:LEN_TRIM(machine2)).EQ.'claire') THEN
            CALL WN('GetInput','Increasing EIRENE runtime for CLAIRE')
            eirtime = eirtime * 2.0
          ELSEIF (machine2(1:LEN_TRIM(machine2)).EQ.'juelich') THEN
            CALL WN('GetInput','Increasing EIRENE runtime for JUELICH')
            eirtime = eirtime * 1.2
          ELSEIF (machine2(1:LEN_TRIM(machine2)).EQ.'joshua') THEN
          ELSE
            CALL WN('GetInput','Unidentified computer: '//
     .                       machine2(1:LEN_TRIM(machine2))//
     .                       ', EIRENE runtime not modified')
          ENDIF
        ENDIF
      ELSEIF (TAG(1:3).EQ.'021') THEN
        CALL ReadI(line,eirdata,0,1,'EIRENE input file')
      ELSEIF (TAG(1:3).EQ.'022') THEN
        CALL ReadI(line,eirmat1,1,5,'EIRENE target material')
c      ELSEIF (TAG(1:3).EQ.'023') THEN
      ELSEIF (TAG(1:3).EQ.'024') THEN
        CALL ReadI(line,eirmat2,1,5,'EIRENE wall material')
      ELSEIF (TAG(1:3).EQ.'058') THEN
c...    Load data for EIRENE pressure gauge specifications (additional
c       surface data):
        IF (line(7:9).EQ.'1.1') THEN
          CALL RdRarn(eirpgdat,eirnpgdat,MAXNAS,-MACHHI,MACHHI,.FALSE.,
     .              -MACHHI,MACHHI,6,'EIRENE pressure gauge specs',ierr)
          IF (ierr.NE.0) CALL ER('GetInput','eirpgdat',*99)
c...      Calculate volume (only toroidal volume at the moment):
          WRITE(fp,*)
          WRITE(fp,'(A)') 'Pressure gauge data:'
          WRITE(fp,'(5X,A4,5A8,2A10)')
     .      'No.','x','y','T (deg)','r (m)','Len (m)','P (mTorr)',
     .      'Vol (m3)'
          DO i1 = 1, eirnpgdat
c            eirpgdat(i1,8) = PI * eirpgdat(i1,5)**2.0 * 2.0 * PI *
c     .                       eirpgdat(i1,2)
            WRITE(fp,'(5X,I4,5F8.3,1P,2E10.2,0P)') INT(eirpgdat(i1,1)),
     .        (eirpgdat(i1,i2),i2=2,8)
          ENDDO
        ELSE
          CALL RdRarn(eirpgdat,eirnpgdat,MAXNAS,-MACHHI,MACHHI,.FALSE.,
     .              -MACHHI,MACHHI,4,'EIRENE pressure gauge specs',ierr)
          IF (ierr.NE.0) CALL ER('GetInput','EIRPGDAT',*99)

          DO i1 = 1, eirnpgdat
            eirpgdat(i1,7) = eirpgdat(i1,5)
            eirpgdat(i1,5) = eirpgdat(i1,4)
            eirpgdat(i1,4) = -1.0
            eirpgdat(i1,6) = -1.0
            eirpgdat(i1,8) =  1.0
          ENDDO
c...      Output gauge specifications:
          WRITE(fp,*)
          WRITE(fp,'(A)') 'Pressure gauge data (no volume):'
          WRITE(fp,'(5X,A4,5A8,2A10)')
     .      'No.','x','y','T (deg)','r (m)','Len (m)','P (mTorr)',
     .      'Vol (m3)'
          DO i1 = 1, eirnpgdat
            WRITE(fp,'(5X,I4,5F8.3,1P,2E10.2,0P)') INT(eirpgdat(i1,1)),
     .        (eirpgdat(i1,i2),i2=2,8)
          ENDDO

        ENDIF

      ELSEIF (TAG(1:3).EQ.'E11') THEN
        CALL ReadI(line,eirbgk    ,0,5,'EIRENE n-n collision option')
      ELSEIF (TAG(1:3).EQ.'E12') THEN
        CALL ReadI(line,eiropacity,-5,6,'EIRENE Lyman alpha opacity')
      ELSEIF (TAG(1:3).EQ.'E13') THEN
        CALL ReadI(line,eirnstrata,0,0,'EIRENE strata specification')
      ELSEIF (TAG(1:3).EQ.'076') THEN
c...    EIRENE surface properties data (default override):
        IF     (line(7:9).EQ.'1.0') THEN

          CALL RdRarn(eirspdat,eirnspdat,MAXNAS3,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,6,'EIRENE surface properties',ierr)
          IF (ierr.NE.0) CALL ER('GetInput','EIRSPDAT',*99)

          WRITE(fp,*) 
          WRITE(fp,*) 'SURFACE PROPERTIES:'
          DO i1 = 1, eirnspdat
c...        Assign default value to the surface recycling coefficient:
            eirspdat(i1,8) = 1.0

           WRITE(fp,'(A,I4,8F12.6)') '> ',i1,(eirspdat(i1,i2),i2=1,8)
          ENDDO

          IF (eirnspdat.GT.0) eirspdatmode = 1

        ELSEIF (line(7:9).EQ.'2.0'.OR.line(7:9).EQ.'2.1') THEN

          CALL RdRarn(eirspdat,eirnspdat,MAXNAS3,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,7,'EIRENE surface properties',ierr)
          IF (ierr.NE.0) CALL ER('GetInput','EIRSPDAT',*99)

          IF (eirnspdat.GT.0) THEN
            IF (line(7:9).EQ.'2.0') THEN
              eirspdatmode = 2
            ELSE
              eirspdatmode = 3 
            ENDIF
          ENDIF

          WRITE(fp,*) 
          WRITE(fp,*) 'SURFACE PROPERTIES:'
          DO i1 = 1, eirnspdat
           WRITE(fp,'(A,I4,8F12.6)') '> ',i1,(eirspdat(i1,i2),i2=1,8)
          ENDDO

        ELSEIF (line(7:9).EQ.'2.2') THEN
          CALL RdRarn(eirspdat,eirnspdat,MAXNAS3,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,8,'EIRENE surface properties',ierr)
          IF (ierr.NE.0) CALL ER('GetInput','EIRSPDAT',*99)
          eirspdatmode = 3 
          WRITE(fp,*) 
          WRITE(fp,*) 'SURFACE PROPERTIES:'
          DO i1 = 1, eirnspdat
c...       Move data to accomodate 2.3:
           eirspdat(i1,10) = eirspdat(i1,9)
           eirspdat(i1,9) = 1.0
           WRITE(fp,'(A,I4,9F12.6)') '> ',i1,(eirspdat(i1,i2),i2=1,9)
          ENDDO
        ELSEIF (line(7:9).EQ.'2.3') THEN
          CALL RdRarn(eirspdat,eirnspdat,MAXNAS3,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,8,'EIRENE surface properties',ierr)
          IF (ierr.NE.0) CALL ER('GetInput','EIRSPDAT',*99)
          eirspdatmode = 3 
          WRITE(fp,*) 
          WRITE(fp,*) 'SURFACE PROPERTIES:'
          DO i1 = 1, eirnspdat
            WRITE(fp,'(A,I4,9F12.6)') '> ',i1,(eirspdat(i1,i2),i2=1,9)
          ENDDO
        ELSE
          CALL ER('RUI','Unsupported version for 076 EIRSPDAT',*99)
        ENDIF

      ELSEIF (TAG(1:3).EQ.'077') THEN
c...    Additional surfaces for EIRENE
        IF     (line(7:9).EQ.'2.0') THEN
          CALL RdRarn(eirasdat,eirnasdat,MAXNAS2,-MACHHI,MACHHI,.FALSE.,
     .              -MACHHI,MACHHI,9,'EIRENE additional surfaces',ierr)
          IF (ierr.NE.0) CALL ER('GetInput','EIRASDAT',*99)

          eirnasdat = MAX(1,eirnasdat)  ! To preserve the default settings in InitializeUnstructuredInput
                                        ! in unstructured_input.f

c...      Shift data around so that it is compatible with the existing code.  The 
c         z-coordinate data is now in elements 8 and 9:
          DO i1 = 1, eirnasdat
            IF (eirasdat(i1,1).EQ.1.0.OR.
     .          eirasdat(i1,1).EQ.98.0.OR.
     .          eirasdat(i1,1).LT.0.0) THEN
              z1 = eirasdat(i1,5)
              DO i2 = 5, 7
                eirasdat(i1,i2) = eirasdat(i1,i2+1)
              ENDDO
              eirasdat(i1,8) = z1
            ENDIF
          ENDDO
          WRITE(SLOUT,*) 'Additional surface data 2.0:'
          DO i1 = 1, eirnasdat
            WRITE(SLOUT,'(I6,10F10.4)') i1,(eirasdat(i1,i2),i2=1,10)
          ENDDO
c For deletion...
        ELSEIF (line(7:9).EQ.'1.1') THEN
          CALL WN('RUI','Additional wall data format 1.1 soon to be '//
     .                  'obsolete (unstructured input tag *077)')
          CALL RdRarn(eirasdat,eirnasdat,MAXNAS2,-MACHHI,MACHHI,.FALSE.,
     .              -MACHHI,MACHHI,8,'EIRENE additional surfaces',ierr)
          IF (ierr.NE.0) CALL ER('GetInput','EIRASDAT',*99)
c...      Shift data around so that it is compatible with the existing code.  The 
c         z-coordinate data is now in elements 8 and 9:
          DO i1 = 1, eirnasdat
            z1 = eirasdat(i1,5)
            DO i2 = 5, 7
              eirasdat(i1,i2) = eirasdat(i1,i2+1)
            ENDDO
            eirasdat(i1,8) = z1
          ENDDO
          WRITE(SLOUT,*) 'Additional surface data:'
          DO i1 = 1, eirnasdat
            WRITE(SLOUT,'(I6,9F10.4)') i1,(eirasdat(i1,i2),i2=1,9)
          ENDDO
c        ELSEIF (line(7:9).EQ.'1.0') THEN
c          CALL RdRarn(eirasdat,eirnasdat,MAXNAS2,-MACHHI,MACHHI,.FALSE.,
c     .              -MACHHI,MACHHI,6,'EIRENE additional surfaces',ierr)
c          IF (ierr.NE.0) CALL ER('GetInput','EIRASDAT',*99)
        ELSE
          CALL ER('RUI','Unsupported version for 077 EIRASDAT',*99)
        ENDIF
      ELSEIF (TAG(1:3).EQ.'E14') THEN
        CALL ReadI(line,eircxd2   ,0,2,'EIRENE CX D2+ production')
      ELSEIF (TAG(1:3).EQ.'E15') THEN
        CALL ReadI(line,eirph2    ,0,1,'EIRENE p-H2 collisions inc.')
      ELSEIF (TAG(1:3).EQ.'E16') THEN
c...    EIRENE puffing surface data:
        IF     (line(7:9).EQ.'1.0') THEN
          CALL RdRarn(eirpuff,eirnpuff,MAXNAS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,7,'EIRENE puffing surfaces',ierr)
          IF (ierr.NE.0) CALL ER('GetInput','E16 EIRPUFF 1',*99)
          IF (eirnpuff.GT.0) eirpmode = 1
          DO i1 = 1, eirnpuff
           WRITE(fp,'(A,I4,6F10.5)') '> ',i1,(eirpuff(i1,i2),i2=1,6)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.1') THEN
          CALL RdRarn(eirpuff,eirnpuff,MAXNAS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,5,'EIRENE puffing surfaces',ierr)
          IF (eirnpuff.GT.0) eirpmode = 2
          IF (ierr.NE.0) CALL ER('GetInput','E16 EIRPUFF 2',*99)
          DO i1 = 1, eirnpuff
           WRITE(fp,'(A,I4,6F10.5)') '> ',i1,(eirpuff(i1,i2),i2=1,6)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.2') THEN
          CALL RdRarn(eirpuff,eirnpuff,MAXNAS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,7,'EIRENE puffing surfaces',ierr)
          IF (eirnpuff.GT.0) eirpmode = 3
          IF (ierr.NE.0) CALL ER('GetInput','E16 EIRPUFF 2',*99)
          DO i1 = 1, eirnpuff
           WRITE(fp,'(A,I4,8F10.5)') '> ',i1,(eirpuff(i1,i2),i2=1,8)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.3') THEN
          CALL RdRarn(eirpuff,eirnpuff,MAXNAS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,7,'EIRENE puffing surfaces',ierr)
          IF (eirnpuff.GT.0) eirpmode = 4
          IF (ierr.NE.0) CALL ER('GetInput','E16 EIRPUFF 2',*99)
          DO i1 = 1, eirnpuff
           WRITE(fp,'(A,I4,8F10.5)') '> ',i1,(eirpuff(i1,i2),i2=1,8)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.4') THEN
          CALL RdRarn(eirpuff,eirnpuff,MAXNAS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,9,'EIRENE puffing surfaces',ierr)
          IF (eirnpuff.GT.0) eirpmode = 4
          IF (ierr.NE.0) CALL ER('GetInput','E16 EIRPUFF 4',*99)
          DO i1 = 1, eirnpuff
           eirpuff (i1,11) = 0.0
           WRITE(fp,'(A,I4,10F10.5)') '> ',i1,(eirpuff(i1,i2),i2=1,10)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.5') THEN
          CALL RdRarn(eirpuff,eirnpuff,MAXNAS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,10,'EIRENE puffing surfaces',ierr)
          IF (eirnpuff.GT.0) eirpmode = 4
          IF (ierr.NE.0) CALL ER('GetInput','E16 EIRPUFF 1.5',*99)
          DO i1 = 1, eirnpuff
           WRITE(fp,'(A,I4,11F10.5)') '> ',i1,(eirpuff(i1,i2),i2=1,11)
          ENDDO
        ELSEIF (line(7:9).EQ.'2.0') THEN
          READ(5,*) cdum1,eirnpuff
          DO i1 = 1, eirnpuff
            READ(5,*) eirpuff(1:14,i1),eircpuff(i1)
            WRITE(fp,'(A,I4,14E10.2,1X,A10)') 'puff: ',i1,
     .        (eirpuff(i2,i1),i2=1,14),eircpuff(i1)
          ENDDO
        ELSE
          CALL ER('RUI','Unsupported version for E16 EIRPUFF',*99)
        ENDIF
      ELSEIF (TAG(1:3).EQ.'E17') THEN
        CALL ReadI(line,eirniter,0,100,'EIRENE self-iterations')
      ELSEIF (TAG(1:3).EQ.'E18') THEN
        CALL ReadR(line,eirrinteg,-1.0,1.0,'EIRENE surface par reflec')
      ELSEIF (TAG(1:3).EQ.'E19') THEN
        CALL ReadR(line,eireinteg,-1.0,1.0,'EIRENE surface eng reflec')
      ELSEIF (TAG(1:3).EQ.'E20') THEN
        CALL ReadR(line,eirermin,-1.0,10.0,'EIRENE thermal cut off')
      ELSEIF (TAG(1:3).EQ.'E21') THEN
        CALL Read2I(line,eirtrc1,eirtrc2,0,999999,'Particle trck range')
      ELSEIF (TAG(1:3).EQ.'E22') THEN
        CALL ReadR(line,eirzaa,-10.0,1000.0,'toroidal circumference')
      ELSEIF (TAG(1:3).EQ.'E23') THEN
        IF     (line(7:9).EQ.'1.1') THEN
          CALL RdIarn(eiraout,eirnaout,MAXASD,1,2,.TRUE.,
     .               -99,MAXASCDAT,5,'EIRENE add cell output',ierr)
          IF (ierr.NE.0) CALL ER('GetInput','E23 EIRAOUT 1',*99)

          IF (eirnaout.EQ.0) THEN
            WRITE(fp,*) 'SETTING DEFAULT FOR EIRAOUT'
            eirnaout = 1
            eiraout(1,1) =  2
            eiraout(1,2) =  0
            eiraout(1,3) = -MAXASCDAT
            eiraout(1,4) =  MAXASCDAT
            eiraout(1,5) =  0
            eiraout(1,6) =  0
          ENDIF

          WRITE(fp,*) 'ADDITIONAL CELL OUTPUT LIST 1.1:'
          DO i1 = 1, eirnaout
            WRITE(fp,'(A,I4,6I10)') '> ',i1,(eiraout(i1,i2),i2=1,6)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.0') THEN
          CALL RdIarn(eiraout,eirnaout,MAXASD,1,2,.TRUE.,
     .               -99,MAXASCDAT,4,'EIRENE add cell output',ierr)
          IF (ierr.NE.0) CALL ER('GetInput','E23 EIRAOUT 1',*99)

c...      Shift data to account for version 1.1 format:
          DO i1 = 1, eirnaout
            DO i2 = 6,3,-1
              eiraout(i1,i2) = eiraout(i1,i2-1)
            ENDDO
            eiraout(i1,2) = 1
          ENDDO

          IF (eirnaout.EQ.0) THEN
            WRITE(fp,*) 'SETTING DEFAULT FOR EIRAOUT'
            eirnaout = 1
            eiraout(1,1) =  2
            eiraout(1,2) = -MAXASCDAT
            eiraout(1,3) =  MAXASCDAT
            eiraout(1,4) =  0
            eiraout(1,5) =  0
          ENDIF

          WRITE(fp,*) 'ADDITIONAL CELL OUTPUT LIST 1.0:'
          DO i1 = 1, eirnaout
            WRITE(fp,'(A,I4,6I10)') '> ',i1,(eiraout(i1,i2),i2=1,6)
          ENDDO
        ELSE
          CALL ER('ReadUnstructuredInput','Invalid version',*99)
        ENDIF
      ELSEIF (TAG(1:3).EQ.'E24') THEN
        CALL ReadR(line,eirtorfrac,0.0001,1.0,'toroidal fraction')
      ELSEIF (TAG(1:3).EQ.'E25') THEN
        CALL ReadR(line,eirsrcmul,0.0,1.0E+6,'source multiplication')
      ELSEIF (TAG(1:3).EQ.'E26') THEN
        CALL ReadI(line,eirfuji,-1,1,'Fujimoto D2+ rates')
      ELSEIF (TAG(1:3).EQ.'E27') THEN
        CALL ReadR(line,eiralloc,0.0,1.0,'CPU-time weighting')
      ELSEIF (TAG(1:3).EQ.'E28') THEN
        IF     (line(7:9).EQ.'1.0') THEN
          CALL RdRarn(eirstrata,eirnstrata,MAXSTRATA,-MACHHI,MACHHI,
     .            .FALSE.,-MACHHI,MACHHI,4,'EIRENE stratum data',ierr)
          IF (ierr.NE.0) CALL ER('GetInput','E28 EIRSTRATA 1',*99)
          WRITE(fp,*) 'EIRENE STRATUM DATA 1.0:'
          DO i1 = 1, eirnstrata
           WRITE(fp,'(A,I4,5F10.5)') '> ',i1,(eirstrata(i1,i2),i2=1,5)
          ENDDO
        ELSE
          CALL ER('ReadUnstructuredInput','Invalid E28 version',*99)
        ENDIF
      ELSEIF (TAG(1:3).EQ.'E29') THEN
        CALL ReadI(line,eiracx    ,0,1,'EIRENE CX reaction included')
      ELSEIF (TAG(1:3).EQ.'E30') THEN
        CALL ReadI(line,eird2ion  ,0,1,'EIRENE D2 ionisation')
      ELSEIF (TAG(1:3).EQ.'E31') THEN
        CALL ReadI(line,eird2dis  ,0,1,'EIRENE D2 dissociation')
      ELSEIF (TAG(1:3).EQ.'E32') THEN
        IF     (line(7:9).EQ.'1.0') THEN
          CALL RdRarn(eiriontime,eirniontime,MAXIONTIME,-MACHHI,MACHHI,
     .            .FALSE.,-MACHHI,MACHHI,7,'EIRENE time-to-ion',ierr)
          IF (ierr.NE.0) CALL ER('GetInput','E32 EIRIONTIME',*99)
          WRITE(fp,*) 'EIRENE TIME-TO-IONISATION PARAMETERS 1.0:'
          DO i1 = 1, eirniontime
            WRITE(fp,'(A,I4,8F10.5)') '> ',i1,(eiriontime(i1,i2),i2=1,8)
            IF (eiriontime(i1,5).GT.MAXBIN-2) 
     .        CALL ER('ReadUnstructuredInput','Number of time bins '//
     .                'requested in E32 exceeds maximum',*99)
          ENDDO
        ELSE
          CALL ER('ReadUnstructuredInput','Invalid E32 version',*99)
        ENDIF
      ELSEIF (TAG(1:3).EQ.'E33') THEN
c...    Some custom input here:
        IF     (line(7:9).EQ.'1.0') THEN
          READ(5,*) cdum1,eirntally
          IF (eirntally.GT.MAXTALLY) 
     .      CALL ER('ReadUnstructuredInput','Too many volume tallies '//
     .              'specified.  Increase MAXTALLY.',*99)
          DO i1 = 1, eirntally
            READ(5,*) (eirtally(i1,i2),i2=1,5) 
            WRITE(PINOUT,'(I2,4(1X,A))') i1,
     .        (eirtally(i1,i2)(1:LEN_TRIM(eirtally(i1,i2))),i2=1,5)
          ENDDO
        ELSE
          CALL ER('ReadUnstructuredInput','Invalid E32 version',*99)
        ENDIF   
      ELSEIF (TAG(1:3).EQ.'E34') THEN
c...    Load EIRENE atomic data multipliers:
c         1 - D ionisation
c         2 - D CX
c         3 - D-D collisions
c         4 - D-D2 collisions
c         5 - D2 ionisation
c         6 - D2 dissociation 1
c         7 - D2 dissociation 2
c         8 - D2-D collisions
c         9 - D2-D2 collisions
c        10 - charge exchange production of D2
c        11 - D+ recombination rate
c        12 - p-D2 collision rate
c
        CALL ReadIR(line,idum1,rdum1,1,12,'EIRENE atomic data mul.')
        eirscale(idum1) = rdum1
      ELSEIF (TAG(1:3).EQ.'E35') THEN
        CALL ReadI(line,eirnsection,1,100,'Number of toroidal sections')
      ELSEIF (TAG(1:3).EQ.'E36') THEN
c...    Specify regions where EIRENE non-standard default surfaces are 
c       transparent:
        IF     (line(7:9).EQ.'1.0') THEN
          CALL RdRarn(eirtrans,eirntrans,MAXTOR,-MACHHI,MACHHI,
     .            .FALSE.,-MACHHI,MACHHI,2,'EIRENE time-to-ion',ierr)
          IF (ierr.NE.0) CALL ER('GetInput','E36 EIRTRANS',*99)
          WRITE(fp,*)
          WRITE(fp,*) 'EIRENE TRANSPARENT SURFACE SPAN 1.0:'
          DO i1 = 1, eirntrans
            WRITE(fp,'(A,I4,3F10.5)') '> ',i1,(eirtrans(i1,i2),i2=1,3)
          ENDDO
        ELSE
          CALL ER('ReadUnstructuredInput','Invalid E36 version',*99)
        ENDIF   
      ELSEIF (TAG(1:3).EQ.'E37') THEN
        CALL ReadI(line,eirntorseg,0,500,'Num sections in EIRENE appro')
      ELSEIF (TAG(1:3).EQ.'E38') THEN
        CALL ReadR(line,eirdtimv,0.0,100.0,'Time dependent mode int.')
      ELSEIF (TAG(1:3).EQ.'E39') THEN
        CALL ReadI(line,eirtrim,0,1,'EIRENE TRIM database option')
      ELSEIF (TAG(1:3).EQ.'E40') THEN
        IF (.TRUE.) THEN
          CALL RDRARN(eirtri,eirntri,MAXNRS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,4,'Additional triangles',IERR)
          DO i1 = 1, eirntri
           WRITE(fp,'(A,1P,1P,8E10.2,0P)') 'TRI:',(eirtri(i1,i2),i2=1,8)
          ENDDO
        ENDIF 
      ELSEIF (TAG(1:3).EQ.'E41') THEN
        CALL ReadI(line,eirphoton,-1,2,'photons being followed')
c ----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'012') THEN					
        CALL ReadI(line,fradd,0,100,'Cells to add to front of ring')	
c ----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'013') THEN					
        CALL ReadI(line,bkadd,0,100,'Cells to add to end of ring')	
c ----------------------------------------------------------------------
c ----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'015') THEN
        CALL ReadI(line,cmodopt,0,1,'Grid source option')

        IF (cmodopt.EQ.1) thesis = .TRUE.
c -----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'016') THEN
        CALL ReadI(line,stagopt,0,3,'Stagger grid')
c -----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'017') THEN
        CALL ReadI(line,stopopt ,0,999,'Special function setting')
      ELSEIF (TAG(1:3).EQ.'063') THEN
        CALL ReadI(line,stopopt2,0,999,'Special function setting 2')
      ELSEIF (TAG(1:3).EQ.'064') THEN
        CALL ReadI(line,stopopt3,0,999,'Special function setting 3')
      ELSEIF (TAG(1:3).EQ.'066') THEN
        CALL ReadI(line,iflexopt(4),0,999,'Integer flex option 4')
        IF (iflexopt(4).EQ.30) osm_mode = 2
      ELSEIF (TAG(1:3).EQ.'071') THEN
        CALL ReadI(line,iflexopt(5),0,999,'Integer flex option 5')
      ELSEIF (TAG(1:3).EQ.'073') THEN
        CALL ReadI(line,iflexopt(6),0,999,'Integer flex option 6')
      ELSEIF (TAG(1:3).EQ.'074') THEN
        CALL ReadI(line,iflexopt(7),0,999,'Integer flex option 7')
      ELSEIF (TAG(1:3).EQ.'075') THEN
        CALL ReadI(line,iflexopt(8),0,999,'Integer flex option 8')
      ELSEIF (TAG(1:3).EQ.'065') THEN
        CALL ReadR(line,rflexopt(1),-99.0,1.0,'Real flex option 1')
      ELSEIF (TAG(1:3).EQ.'067') THEN
        CALL ReadR(line,rflexopt(2),-1.0,99.0,'Real flex option 2')
      ELSEIF (TAG(1:3).EQ.'068') THEN
        CALL ReadR(line,rflexopt(3),0.0,1.0E+8,'Real flex option 3')
      ELSEIF (TAG(1:3).EQ.'069') THEN
        CALL ReadR(line,rflexopt(4),-1.0,100.0,'Real flex option 4')
      ELSEIF (TAG(1:3).EQ.'070') THEN
        CALL ReadR(line,rflexopt(5),0.0,10.0,'Real flex option 5')
      ELSEIF (TAG(1:3).EQ.'072') THEN
        CALL ReadR(line,rflexopt(6),0.0,1.0E+25,'Real flex option 6')
c -----------------------------------------------------------------------

c ----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'025') THEN
        CALL ReadI(line,haldata,0,1,'CMOD Halpha data')
c -----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'026') THEN
        CALL ReadI(line,pincode,0,5,'PIN neutral code option')
c -----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'027') THEN
        CALL ReadI(line,tarsource,0,8,'Target data source')
      ELSEIF (TAG(1:3).EQ.'078') THEN
        CALL Read2R(line,tarshift(IKHI),tarshift(IKLO),-HI,HI,'T-shift')
      ELSEIF (TAG(1:3).EQ.'085') THEN
        CALL ReadI(line,tarshiftopt,0,1,'Target shift option')
      ELSEIF (TAG(1:3).EQ.'081') THEN
c...    High index target data for boundary condition relaxation:
        IF     (line(7:9).EQ.'1.1') THEN
          CALL RDRARN(LPDATI2,NLPDATI2,MAXINS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,6,'INNER RELAX. DATA',IERR)
          IF (ierr.NE.0) CALL ER('ReadUnstructuredInput',
     .                           '081 error',*99)          
        ELSE
          CALL ER('ReadUnstructuredInput','Unknown 081 version',*99)        
        ENDIF
      ELSEIF (TAG(1:3).EQ.'082') THEN
c...    Low index target data for boundary condition relaxation:
        IF     (line(7:9).EQ.'1.1') THEN
          CALL RDRARN(LPDATO2,NLPDATO2,MAXINS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,6,'OUTER RELAX. DATA',IERR)
          IF (ierr.NE.0) CALL ER('ReadUnstructuredInput',
     .                           '082 error',*99)
        ELSE
          CALL ER('ReadUnstructuredInput','Unknown 082 version',*99)
        ENDIF

      ELSEIF (TAG(1:3).EQ.'083') THEN
c...    Low index target data for SOL21,24 over-ride:
        IF     (line(7:9).EQ.'1.1') THEN
          CALL RDRARN(s21_datao,s21_ndatao,MAXNRS,-MACHHI,MACHHI,
     >                .FALSE.,-MACHHI,MACHHI,9,'OUTER SOL21 DATA',IERR)
          IF (ierr.NE.0) CALL ER('ReadUnstructuredInput',
     .                           '083 error',*99)
        ELSEIF (line(7:9).EQ.'1.0') THEN
          CALL RDRARN(s21_datao,s21_ndatao,MAXNRS,-MACHHI,MACHHI,
     >                .FALSE.,-MACHHI,MACHHI,7,'OUTER SOL21 DATA',IERR)
          IF (ierr.NE.0) CALL ER('ReadUnstructuredInput',
     .                           '083 error',*99)
          DO i1 = 1, s21_ndatao
            s21_datao(i1,9 ) = -1.0
            s21_datao(i1,10) = -1.0
          ENDDO
        ELSE
          CALL ER('ReadUnstructuredInput','Unknown 083 version',*99)
        ENDIF
      ELSEIF (TAG(1:3).EQ.'084') THEN
c...    High index target data for SOL21,24 over-ride:
        IF     (line(7:9).EQ.'1.1') THEN
          CALL RDRARN(s21_datai,s21_ndatai,MAXNRS,-MACHHI,MACHHI,
     >                .FALSE.,-MACHHI,MACHHI,9,'INNER SOL21 DATA',IERR)
          IF (ierr.NE.0) CALL ER('ReadUnstructuredInput',
     .                           '084 error',*99)
        ELSEIF (line(7:9).EQ.'1.0') THEN
          CALL RDRARN(s21_datai,s21_ndatai,MAXNRS,-MACHHI,MACHHI,
     >                .FALSE.,-MACHHI,MACHHI,7,'INNER SOL21 DATA',IERR)
          IF (ierr.NE.0) CALL ER('ReadUnstructuredInput',
     .                           '084 error',*99)
          DO i1 = 1, s21_ndatai
            s21_datai(i1,9 ) = -1.0
            s21_datai(i1,10) = -1.0
          ENDDO
        ELSE
          CALL ER('ReadUnstructuredInput','Unknown 084 version',*99)
        ENDIF

      ELSEIF (TAG(1:3).EQ.'088') THEN
c...    High index target data for interpolation:
        version = GetInputVersion(line)
        IF     (version.EQ.1.0) THEN
          tarintermode(IKHI) = 1.0
          CALL RDRARN(tarinter(1,1,IKHI),tarninter(IKHI),
     >                MAXNRS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,6,
     .                'HIGH INDEX INTER DATA',IERR)
          IF (ierr.NE.0) CALL ER('ReadUnstructuredInput',
     .                           '088 error',*99)
          DO i1 = 1, tarninter(IKHI)
            WRITE(SLOUT,'(A,I6,1P,4E10.2,0P)') 
     .        '088:',i1,(tarinter(i1,i2,IKHI),i2=1,4)
          ENDDO
        ELSE
          tarintermode(IKHI) = 0.0
          CALL RDRARN(tarinter(1,1,IKHI),tarninter(IKHI),
     >                MAXNRS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,3,
     .                'HIGH INDEX INTER DATA',IERR)
          IF (ierr.NE.0) CALL ER('ReadUnstructuredInput',
     .                           '088 error',*99)

          DO i1 = 1, tarninter(IKHI)
            WRITE(SLOUT,'(A,I6,1P,4E10.2,0P)') 
     .        '088:',i1,(tarinter(i1,i2,IKHI),i2=1,4)
          ENDDO
        ENDIF

      ELSEIF (TAG(1:3).EQ.'089') THEN
c...    Low index target data for interpolation:
        IF     (line(7:9).EQ.'1.0') THEN
          tarintermode(IKLO) = 1.0
          CALL RDRARN(tarinter(1,1,IKLO),tarninter(IKLO),
     .                MAXNRS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,6,
     .                'LOW INDEX INTER DATA',IERR)
          IF (ierr.NE.0) CALL ER('ReadUnstructuredInput',
     .                           '089 error',*99)
        ELSE
          tarintermode(IKLO) = 0.0
          CALL RDRARN(tarinter(1,1,IKLO),tarninter(IKLO),
     .                MAXNRS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,3,
     .                'LOW INDEX INTER DATA',IERR)
          IF (ierr.NE.0) CALL ER('ReadUnstructuredInput',
     .                           '089 error',*99)
        ENDIF
c -----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'028') THEN
        CALL ReadR(line,grd_minpl,LO,HI,'Min poloidal side length')
c -----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'029') THEN
        CALL ReadI(line,grd_refine,0,13,'Grid refinement option')
c -----------------------------------------------------------------------

c -----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'032') THEN
        CALL ReadR(line,grd_range,-HI,HI,'Reginement region')
c      ELSEIF (TAG(1:3).EQ.'G01') THEN
c...    Not sure if Dave has reserved G01 already, and the web server is down:
c        CALL ReadR(line,grd_thresh,0.0,100.0,'Refinement threshold')
c -----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'033') THEN
        CALL ReadR(line,osm_range,0.0,1.0,'SOL 22 error region')
      ELSEIF (TAG(1:3).EQ.'042') THEN
        CALL ReadR(line,switch(SWPOW2),-1.0,22.0,'SOL 22 pow dist')
      ELSEIF (TAG(1:3).EQ.'059') THEN
        CALL ReadR(line,switch(SWPOW3),-1.0,19.0,'SOL 22 mock pow dist')
      ELSEIF (TAG(1:3).EQ.'043') THEN
        CALL ReadR(line,switch(SWION2),-1.0,21.0,'SOL 22 ion dist')
      ELSEIF (TAG(1:3).EQ.'044') THEN
        CALL ReadI(line,osm_store,-1,100,'load stored PIN sources')
      ELSEIF (TAG(1:3).EQ.'046') THEN
        CALL ReadI(line,osm_probe,-2,5,'Probe for pressure reference')
      ELSEIF (TAG(1:3).EQ.'048') THEN
        CALL ReadI(line,s21_mode,0,5,'PIN source use option')
      ELSEIF (TAG(1:3).EQ.'050') THEN
        CALL ReadI(line,osm_matchs,0,2,'Match plasma at symmetry point')
      ELSEIF (TAG(1:3).EQ.'051') THEN
        CALL ReadI(line,osm_matchp,0,4,'Match pressure at probe loc.')
      ELSEIF (TAG(1:3).EQ.'053') THEN
        CALL ReadI(line,osm_symopt,0,5,'SOL 22 symmetry point option')
      ELSEIF (TAG(1:3).EQ.'055') THEN
        CALL ReadI(line,osm_matcht,0,6,'Match temp. at probe loc.')
      ELSEIF (TAG(1:3).EQ.'060') THEN
        CALL ReadI(line,osm_preopt,-1,9,'Use SOL22p prescription')
      ELSEIF (TAG(1:3).EQ.'061') THEN
        CALL ReadI(line,osm_recopt,0,7,'Calculate recombination adj')
      ELSEIF (TAG(1:3).EQ.'087') THEN
        CALL ReadI(line,osmmock,0,1,'Mock power term option')
      ELSEIF (TAG(1:3).EQ.'079') THEN
c...    Plasma data for uniform private flux zone:
        IF (line(7:9).EQ.'1.0') THEN
          CALL RdRarn(osmppv,osmnppv,MAXNRS,-MACHHI,MACHHI,.FALSE.,
     .              -MACHHI,MACHHI,5,'Uniform PP data',ierr)
          IF (ierr.NE.0) CALL ER('GetInput','OSMPPV',*99)
c           DO i1 = 1, osmnppv
c            WRITE(0,'(A,1P,6E10.2,0P)') '> ',(osmppv(i1,i2),i2=1,6)
c          ENDDO
        ELSE
          CALL ER('RUI','Unsupported version for 070 OSMPPV',*99)
        ENDIF

c...    Using 'S' prefix:
      ELSEIF (TAG(1:3).EQ.'S70') THEN
        CALL ReadR(line,osmbulkn,-HI,HI,'bluk n over-ride')
      ELSEIF (TAG(1:3).EQ.'S71') THEN
        CALL ReadR(line,osmbulkte,-HI,HI,'bluk electron T  over-ride')
      ELSEIF (TAG(1:3).EQ.'S72') THEN
        CALL ReadR(line,osmbulkti,-HI,HI,'bluk ion T over-ride')
      ELSEIF (TAG(1:3).EQ.'S73') THEN
        CALL ReadR(line,osmbulkv,-HI,HI,'bluk plasma v over-ride')
      ELSEIF (TAG(1:3).EQ.'S74') THEN
c...
c
c
        READ(line(7:9),*) s28mode
        IF (s28mode.EQ.1.0.OR.s28mode.EQ.1.1) THEN
          CALL RDRARN(osms28,osmns28,MAXNRS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,7,'SOL28 specifications 1.x',IERR)
          DO i1 = 1, osmns28
           osms28(i1,9)  = 0.0
           osms28(i1,10) = 0.0
           WRITE(fp,'(A,1P,1P,8E10.2,0P)') 'S28:',(osms28(i1,i2),i2=1,8)
          ENDDO
        ELSEIF (s28mode.EQ.2.0) THEN
          CALL RDRARN(osms28,osmns28,MAXNRS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,7,'SOL28 specifications 2.0',IERR)
          DO i1 = 1, osmns28
           osms28(i1,9)  = 0.0
           osms28(i1,10) = 0.0
           WRITE(fp,'(A,1P,1P,8E10.2,0P)') 'S28:',(osms28(i1,i2),i2=1,8)
          ENDDO
        ELSEIF (s28mode.EQ.2.1) THEN
c...      Added ability to limit application of SOL28 parameters to specific
c         rings:
          CALL RDRARN(osms28,osmns28,MAXNRS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,9,'SOL28 specifications 2.1',IERR)
          DO i1 = 1, osmns28
           WRITE(fp,'(A,1P,1P,10E10.2,0P)') 
     .       'S28:',(osms28(i1,i2),i2=1,10)
          ENDDO
        ELSEIF (s28mode.GE.3.0) THEN
c...      Improved:
          CALL RDRARN(osms28,osmns28,MAXNRS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,11,'SOL28 specifications 3.0',IERR)
          DO i1 = 1, osmns28
           WRITE(fp,'(A,1P,1P,12E10.2,0P)') 
     .       'S28:',(osms28(i1,i2),i2=1,12)
          ENDDO
        ELSE
          CALL ER('RUI','Unsupported version for S74 OSMS28',*99)
        ENDIF
      ELSEIF (TAG(1:3).EQ.'S75') THEN
c...
        CALL Read2I(line,s28ion,s28ionpfz,-5,20,'SOL28 ion src option')
      ELSEIF (TAG(1:3).EQ.'S83') THEN
c...
        CALL ReadR(line,s28ionbound,-99.0,5.0,'SOL28 ion limit')
      ELSEIF (TAG(1:3).EQ.'S84') THEN
c...
        CALL Read2I(line,s28cfp,s28cfppfz,0,25,'SOL28 cross-f ion opt')
      ELSEIF (TAG(1:3).EQ.'S76') THEN
c...
        CALL Read2I(line,s28mom,s28mompfz,-5,20,'SOL28 momentum option')
      ELSEIF (TAG(1:3).EQ.'S85') THEN
c...
        CALL Read2I(line,s28cfm,s28cfmpfz,0,1,'SOL28 cross-f mom opt')
      ELSEIF (TAG(1:3).EQ.'S98') THEN
c...
        CALL Read2I(line,s28te,s28tepfz,0,21,'SOL28 Te option')
      ELSEIF (TAG(1:3).EQ.'S86') THEN
c...
        CALL Read2I(line,s28ti,s28tipfz,0,20,'SOL28 Ti option')
      ELSEIF (TAG(1:3).EQ.'S87') THEN
c...
        CALL ReadR(line,s28tiratio,0.1,10.0,'SOL28 Ti/Te ratio')
      ELSEIF (TAG(1:3).EQ.'S88') THEN
c...
        CALL Read2I(line,s28superdet,s28superdetpfz,-5,5,'SOL28 M opt')
      ELSEIF (TAG(1:3).EQ.'S77') THEN
c...
        CALL ReadR(line,s28momfr,0.0,HI,'SOL28 momentum balance')
      ELSEIF (TAG(1:3).EQ.'S78') THEN
c...
        CALL ReadI(line,s22pfzion,0,1,'feeding PFZ ionisation')
      ELSEIF (TAG(1:3).EQ.'S79') THEN
c...
        CALL ReadR(line,setmach,0.0,1.0,'setting mach number in SOL24')
      ELSEIF (TAG(1:3).EQ.'S80') THEN
c...
        CALL ReadR(line,muldensity,-HI,HI,'SOL24,28 density multi')
      ELSEIF (TAG(1:3).EQ.'S81') THEN
c...
        CALL ReadR(line,addtemp,-10.0,99.0,'SOL24,28 density temp add')
      ELSEIF (TAG(1:3).EQ.'S82') THEN
c...
        CALL ReadR(line,s28b34,0.0,1.0,'SOL28 Te profile exponent')
      ELSEIF (TAG(1:3).EQ.'S89') THEN
c...
        CALL ReadR(line,s28b34det,0.0,1.0,'SOL28 det Te profile exp')
      ELSEIF (TAG(1:3).EQ.'S90') THEN
c...
        CALL Read2I(line,s28probe,s28probepfz,0,999,'SOL28 probe data')
      ELSEIF (TAG(1:3).EQ.'S91') THEN
c...
        CALL Read2I(line,s28rec,s28recpfz,0,20,'SOL28 recomb. opt')
      ELSEIF (TAG(1:3).EQ.'S92') THEN
c...
        CALL ReadR(line,s28momTe,0.0,30.0,'SOL28 mom Te cutoff')
      ELSEIF (TAG(1:3).EQ.'S93') THEN
c...
c
c
        READ(line(7:9),*) version
        IF (version.EQ.1.0) THEN
          CALL RDRARN(osms28over,osmns28over,MAXNRS,-MACHHI,MACHHI,
     .                .FALSE.,-MACHHI,MACHHI,8,'SOL28 over-ride 1.0',
     .                IERR)
          DO i1 = 1, osmns28over
           WRITE(fp,'(A,1P,9E10.2)') 'OVER:',(osms28over(i1,i2),i2=1,9)
          ENDDO
        ELSE
          CALL ER('RUI','Unsupported version for S93 OSMS28OVER',*99)
        ENDIF
      ELSEIF (TAG(1:3).EQ.'S94') THEN
c...
        CALL ReadI(line,s28sym,0,1,'Symmeterize non-SOL28 rings')
      ELSEIF (TAG(1:3).EQ.'S95') THEN
c...
        CALL Read2I(line,s28cfpnt,s28cfppfznt,0,7,'Near target cf src')
      ELSEIF (TAG(1:3).EQ.'S96') THEN
c...
        CALL Read2I(line,s28cfpdrft,s28cfppfzdrft,-3,6,'Rad. drift')
      ELSEIF (TAG(1:3).EQ.'S97') THEN
c...
        CALL Read2I(line,s28nemode,s28nemodepfz,-3,2,'upstream ne/jsat')

c -----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'030') THEN
        CALL ReadI(line,rel_opt,0,3,'Relaxation option')
      ELSEIF (TAG(1:3).EQ.'031') THEN
        CALL ReadR(line,rel_frac,0.0,1.0,'PIN relaxation fraction')
      ELSEIF (TAG(1:3).EQ.'086') THEN
        CALL ReadI(line,relreset,0,2,'Reset PIN sources with each step')
      ELSEIF (TAG(1:3).EQ.'034') THEN
        CALL ReadI(line,rel_nstep,1,100,'Number of steps')
      ELSEIF (TAG(1:3).EQ.'035') THEN
        CALL ReadI(line,rel_niter,1,100,'Iteratons per step')
      ELSEIF (TAG(1:3).EQ.'036') THEN
        CALL ReadR(line,rel_pace,-1.0,9.0,'Target relaxation pace')
      ELSEIF (TAG(1:3).EQ.'045') THEN
        CALL Read2R(line,rel_bound1,rel_bound2,0.0,1.0,'BC relax. bnd')
      ELSEIF (TAG(1:3).EQ.'047') THEN
        CALL ReadR(line,rel_tol,0.01,100.0,'Pressure balance tolerance')
      ELSEIF (TAG(1:3).EQ.'052') THEN
        CALL ReadI(line,osm_powopt,0,2,'Mock power option')
      ELSEIF (TAG(1:3).EQ.'041') THEN
        CALL Read2I(line,osm_watch1,osm_watch2,0,MAXNRS,'SOL 22 watch')
      ELSEIF (TAG(1:3).EQ.'062') THEN
c
c       Load data for ionisation front specification for SOL22p:
c
        CALL RdRarn(osm_ionfnt,osm_nfnt,MAXFNT,-MACHHI,MACHHI,
     .      .FALSE.,-MACHHI,MACHHI,2,'Ionisation front spec.     ',ierr)

        IF (ierr.NE.0) CALL ER('GetInput','Reading OSM_IONFNT',*99)
      ELSEIF (TAG(1:3).EQ.'080') THEN
c...    Relaxation over-ride:
        IF     (line(7:9).EQ.'1.0') THEN
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE.,
     .                0.0,MACHHI,6,'Relaxation data 1.0',IERR)
          relmode = 0
          DO i1 = 1, rel_ndata
           WRITE(fp,'(A,1P,7E10.2)') 'RD: ',(rel_data(i1,i2),i2=1,7)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.1') THEN
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,8,'Relaxation data 1.1',IERR)
          relmode = 1
          DO i1 = 1, rel_ndata
           WRITE(fp,'(A,1P,9E10.2)') 'RD: ',(rel_data(i1,i2),i2=1,9)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.2') THEN
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,3,'Relaxation data 1.2',IERR)
          relmode = 2
          DO i1 = 1, rel_ndata
           WRITE(fp,'(A,1P,9E10.2)') 'RD: ',(rel_data(i1,i2),i2=1,4)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.3') THEN
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,8,'Relaxation data 1.3',IERR)
          relmode = 3
          DO i1 = 1, rel_ndata
           WRITE(fp,'(A,1P,9E10.2)') 'RD: ',(rel_data(i1,i2),i2=1,9)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.4') THEN
c...      Over-riding the additional cell plasma parameters:
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,2,'Relaxation data 1.4',IERR)
          relmode = 4
          DO i1 = 1, rel_ndata
           WRITE(fp,'(A,1P,9E10.2)') 'RD: ',(rel_data(i1,i2),i2=1,2)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.5') THEN
c...      Over-riding the additional cell plasma parameters and bulk plasma values:
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,6,'Relaxation data 1.5',IERR)
          relmode = 5
          DO i1 = 1, rel_ndata
           WRITE(fp,'(A,1P,9E10.2)') 'RD: ',(rel_data(i1,i2),i2=1,7)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.6') THEN
c...      Over-riding puff strength:
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,1,'Relaxation data 1.6',IERR)
          relmode = 6
          DO i1 = 1, rel_ndata
           WRITE(fp,'(A,1P,2E10.2)') 'RD: ',(rel_data(i1,i2),i2=1,2)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.7') THEN
c...      Over-riding bulk plasma values:
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,4,'Relaxation data 1.7',IERR)
          relmode = 7
          DO i1 = 1, rel_ndata
           WRITE(fp,'(A,1P,5E10.2)') 'RD: ',(rel_data(i1,i2),i2=1,5)
          ENDDO
        ELSEIF (line(7:9).EQ.'2.0') THEN
c...      Flexible over-ride option:
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,6,'Relaxation data 2.0',IERR)
          relmode = 20
          WRITE(fp,*) 'RELMODE:',relmode
          DO i1 = 1, rel_ndata
           WRITE(fp,'(A,1P,7E10.2)') 'RD: ',(rel_data(i1,i2),i2=1,7)
          ENDDO
        ELSE
          CALL ER('RUI','Unsupported version for 080 REL_DATA',*99)
        ENDIF
c...    Set RELMODE to zero if no data is listed in the input file:       
        IF (rel_ndata.EQ.0) relmode = 0
c ----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'037') THEN
        CALL ReadI(line,adp_opt,0,2,'Grid adaptation option')
      ELSEIF (TAG(1:3).EQ.'038') THEN
        CALL ReadI(line,adp_region,1,3,'Adaptation region')
      ELSEIF (TAG(1:3).EQ.'039') THEN
        CALL ReadR(line,adp_upper,1.0,HI,'Bound for refinement')
      ELSEIF (TAG(1:3).EQ.'040') THEN
        CALL ReadR(line,adp_lower,1.0,HI,'Bound for coarsification')
c -----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'056') THEN
        CALL ReadI(line,prb_align,0,1,'Adjust position of FSP data')
      ELSEIF (TAG(1:3).EQ.'057') THEN
        CALL ReadR(line,prb_shift,-5.0,99.0,'Forced FSP shift')
        IF (prb_shift.NE.99.0) prb_shift = prb_shift / 1.0E+03

c -----------------------------------------------------------------------
c...  Vacuum grid options:
      ELSEIF (TAG(1:3).EQ.'V01') THEN
        IF (line(7:9).EQ.'1.0') THEN
          CALL RDRARN(vacseg,vacnseg,MAXNKS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,3,'Vacuum grid data 1.0',IERR)

          DO i1 = 1, vacnseg
           ! jdemod - changed to f12 due to some formatting errors in INTEL compiler
           WRITE(fp,'(A,1P,4G12.4)') 'VS: ',(vacseg(i1,i2),i2=1,4)
          ENDDO

c          iflexopt(2) = 20
        ELSE
          CALL ER('RUI','Unsupported version for V01 VACSEG',*99)
        ENDIF

      ELSEIF (TAG(1:3).EQ.'V02') THEN
        IF     (line(7:9).EQ.'1.3') THEN
          CALL RDRARN(vacpla,vacnpla,MAXNKS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,7,'Vacuum grid plasma 1.3',IERR)
          DO i1 = 1, vacnpla
           IF (vacpla(i1,1).NE.-3.0.AND.vacpla(i1,1).NE.-4.0.AND.
     .         vacpla(i1,1).NE. 0.0 )
     .       CALL ER('UnstrInput','V02 1.3 must have type -3.0, -4.0 '//
     .                            'or 0.0',*99)
          ENDDO
          DO i1 = 1, vacnpla
           ! jdemod - changed to g12 due to some formatting errors in INTEL compiler
           WRITE(fp,'(A,1P,8G12.4)') 'VS: ',(vacpla(i1,i2),i2=1,8)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.2') THEN
          CALL RDRARN(vacpla,vacnpla,MAXNKS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,6,'Vacuum grid plasma 1.2',IERR)
          DO i1 = 1, vacnpla
           ! jdemod - changed to g12 due to some formatting errors in INTEL compiler
           WRITE(fp,'(A,1P,7g12.4)') 'VS: ',(vacpla(i1,i2),i2=1,7)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.1') THEN
          CALL RDRARN(vacpla,vacnpla,MAXNKS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,7,'Vacuum grid plasma 1.1',IERR)
          DO i1 = 1, vacnpla
           ! jdemod - changed to g12 due to some formatting errors in INTEL compiler
           WRITE(fp,'(A,1P,8g12.4)') 'VS: ',(vacpla(i1,i2),i2=1,8)
          ENDDO
        ELSEIF (line(7:9).EQ.'1.0') THEN
          CALL RDRARN(vacpla,vacnpla,MAXNKS,-MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,6,'Vacuum grid plasma 1.0',IERR)
c...      Adjust input to make compatible with version 1.1:
          DO i1 = 1, vacnpla
            DO i2 = 8, 3, -1
              vacpla(i1,i2) = vacpla(i1,i2-1)
            ENDDO
            vacpla(i1,2) = vacpla(i1,1)
          ENDDO
          DO i1 = 1, vacnpla
           ! jdemod - changed to g12 due to some formatting errors in INTEL compiler
           WRITE(fp,'(A,1P,8g12.4)') 'VS: ',(vacpla(i1,i2),i2=1,8)
          ENDDO
        ELSE
          CALL ER('RUI','Unsupported revision for V02 VACPLA',*99)
        ENDIF
c -----------------------------------------------------------------------

c ...next...090 -- but it is in use!


c -----------------------------------------------------------------------
      ELSEIF (TAG(1:3).EQ.'090') THEN
        CALL ReadI(line,outtarget,0,1,'Output target data')
      ELSEIF (TAG(1:3).EQ.'091') THEN
        CALL ReadI(line,outwall,0,1,'Output wall data')
      ELSEIF (TAG(1:3).EQ.'092') THEN
        CALL ReadI(line,outtarget,0,1,'Output target data')
      ELSEIF (TAG(1:3).EQ.'093') THEN
        CALL ReadI(line,outcell,0,1,'Output cell data')
      ELSEIF (TAG(1:3).EQ.'094') THEN
        CALL ReadI(line,outgeom,0,1,'Output geometry data')
      ELSEIF (TAG(1:3).EQ.'095') THEN
        CALL ReadI(line,outpoly,0,1,'Output polygon data')
      ELSEIF (TAG(1:3).EQ.'096') THEN
        CALL ReadI(line,outpin,0,1,'Output PIN data')
      ELSEIF (TAG(1:3).EQ.'097') THEN
        CALL ReadI(line,out_source,0,1,'PIN sources')
      ELSEIF (TAG(1:3).EQ.'098') THEN
        CALL ReadI(line,out_plasma,0,1,'PIN plasma')
      ELSEIF (TAG(1:3).EQ.'099') THEN
        CALL ReadI(line,out_geom  ,0,1,'PIN geometry')
c     
c -----------------------------------------------------------------------
c
c     TAG 2??: 200 series tags are related to SOL22. 
c
c     200 series options are unstructured input for SOL22 - the routine to process
c         these is in mod_sol22_input
c
      elseif (tag(1:1).eq.'2') then       
         call sol22_unstructured_input(tag,line,ierr)
         
      ! TAG X: Options related to SOL29.
      elseif (tag(1:1).eq.'X') then
        call sol29_unstructured_input(tag, line)
c     
c -----------------------------------------------------------------------
c
c     TAG A05: 
c
      ELSEIF (tag(1:3).EQ.'A05') THEN
c
c       ne_opt - selects how electron density is defined 
c       0 -> ne = nb    1 -> ne = nb + sigma nz (from fluid code)
c
        CALL ReadI(line,ne_opt,0,1,'Electron density option')
c
c -----------------------------------------------------------------------
c
c
c     TAG A06: 
c
      ELSEIF (tag(1:3).EQ.'A06') THEN
c
c       Option to write a JET TRAN file for POST-PROCESSOR use
c       Written to Unit 41 - 0 = off - 1 = 0n .. default ON
c
c
        CALL ReadI(line,write_tran,0,1,'TRAN FILE PRINT OPTION')
c
c -----------------------------------------------------------------------
c
c
c     TAG A07: 
c
      ELSEIF (tag(1:3).EQ.'A07') THEN
c
c       Option to write a netcdf version of the raw data file 
c       0 = off   1 = 0n .. default OFF
c
c
        CALL ReadI(line,netcdf_opt,0,1,'WRITE NETCDF FORMAT RAW OUTPUT')
c
c
c -----------------------------------------------------------------------
c
c     TAG C21 - read in a myultiplier for PINQE as used in the 
c                  Dperp/Xperp extractor
c
      ELSEIF (tag(1:3).EQ.'C21') THEN
        CALL ReadR(line,dp_pinqe_mult,0.0,HI,
     >             'PINQE Multiplier for Extractor')
c -----------------------------------------------------------------------
c
c     TAG C22 - Calculate a velocity shifted line profile for the 
c               specified neutral impurity line.
c
      ELSEIF (tag(1:3).EQ.'C22') THEN
        CALL ReadI(line,line_profile_opt,0,1,
     >             'Line Profile Calculation Switch')
c
c       Read in the additional required data. 
c
c       Read in ADAS Selector data
c
c        write(0,'(a,4i5)') 'ADAS:',lp_isele,
c     >          len(lp_adasid),len(lp_adasex),
c     >          lp_adasyr
c
        call rdg1(line3,lp_ADASID,lp_ADASYR,lp_ADASEX,
     >                 lp_ISELE,lp_ISELR,lp_ISELX,lp_ISELD,IERR)
c
c        write(0,'(a,4i5)') 'ADAS:',lp_isele,
c     >          len(lp_adasid),len(lp_adasex),
c     >          lp_adasyr
c        write(0,'(a,a,a))') 'ADAS:',lp_adasid,':'
c        write(0,'(a,a,a))') 'ADAS:',lp_adasex,':'
c
c
c       Read in LOS data. 
c         
        call rd_lp_los(line3,lp_robs,lp_zobs,lp_theta,lp_dtheta,
     >                 lp_instrument_width,lp_bin_width,ierr)
c
c       Convert angles in degrees to radians
c        
        lp_theta = lp_theta * degrad
        lp_dtheta = lp_dtheta * degrad
c
c -----------------------------------------------------------------------
c
c
c -----------------------------------------------------------------------
c
c     TAG D37 and D38 
c     ADAS IONIZATION AND RECOMBINATION RATE MODIFIERS
c
c     These values can be used to modify the rates read in from ADAS
c     The default values should always be set to 1.0 and these
c     should be used with care if used at all.  
c
c     TAG D37 - adas_iz_rate_mult
c
      ELSEIF (tag(1:3).EQ.'D37') THEN
        CALL ReadR(line,adas_iz_rate_mult,0.0,HI,
     >             'ADAS Ionization rate multiplier')
c
c     TAG D38 - adas_rec_rate_mult
c
      ELSEIF (tag(1:3).EQ.'D38') THEN
        CALL ReadR(line,adas_rec_rate_mult,0.0,HI,
     >             'ADAS Recombination rate multiplier')
c
c -----------------------------------------------------------------------
c
c     TAG D39 : Alternate Sputter data specifier - used to select one of 
c               several custom sputter datasets - usually based
c               on different impact angles
c 
      elseif (tag(1:3).EQ.'D39') THEN
c
c     Secondary sputter data specifier
c
        CALL ReadR(line,extra_sputter_angle,-10.0,90.0,
     >             'Extra Sputter Angle Opt')
c
c -----------------------------------------------------------------------
c
c     TAG D40 : Flux fraction for alternate bombarding ion sputter 
c               calculations. Goes with CBOMBF and CBOMBZ to allow for 
c               specification of trace impurity sputtering in cases
c               where hydrogenic sputtering is expected to be negligible.
c
c               Ideally this should be expanded to allow for 
c               hydrogen+a specifiable distribution of impurity charge states                              
c 
      elseif (tag(1:3).EQ.'D40') THEN
c
c     Bombarding ion flux fraction
c
        CALL ReadR(line,cbomb_frac,0.0,1.0,
     >             'Bombarding ion flux fraction')
     
!     TAG D41-43
!     These options are related to sputtering from SiC using the SiC 
!     mixed-material model in Abrams NF 2022. Note this requires 
!     csputopt = 8.
!     D41: Usage switch. Useful for comparisons of SiC to its 
!          constituents using the same yield data that is hardocded in 
!          for SiC.
!          1 = Normal. Target is SiC.
!          2 = Target is graphite.
!          3 = Target is silicon.
!     D42: Flux fraction of the carbon in the incoming plasma flux.
!     D43: Flux fraction of the silicon in the incoming plasma flux.    

      elseif (tag(1:3).eq.'D41') then
        call readi(line, mm_usage, 0, 2, 
     >    'SiC mixed-material model usage switch')
      elseif (tag(1:3).eq.'D42') then
        call readr(line, frac_c, 0.0, 1.0, 
     >   'Fraction of C in SiC mixed-material model')
      elseif (tag(1:3).eq.'D43') then
        call readr(line, frac_si, 0.0, 1.0, 
     >    'Fraction of Si in SiC mixed-material model')
     
!     If using this model, then we must ensure CENUTD = 0, it doesn't
!     make sense otherwise.
      if (cneutd.ne.0.and.csputopt.eq.8) then
        write(0,*) 'Warning! SiC model in use, CNEUTD forced to 0.'
        write(0,*) 'Change CNEUTD to 0 to avoid this warning.'
        cneutd = 0
      endif
      
!     TAG D44
!     Usage switch for TiB2 (0) or ZrB2 (1). 
      elseif (tag(1:3).eq.'D44') then
        call readi(line, tib2_or_zrb2, 0, 1, 
     >    'TiB2 or ZrB2 usage switch')

c
c
c -----------------------------------------------------------------------
c 
c     TAG F11 - Uedge background option 
c
      ELSEIF (tag(1:3).EQ.'F11') THEN
        CALL ReadI(line,uedge_bg,0,1,'UEDGE background option')
c
c -----------------------------------------------------------------------
c
c     Options affecting the calculation/interpretation of fluid code
c     target conditions.
c
c     Option F12 sets groups of these flags so that input can be 
c     simplified. Any individual flag options appearing later in the 
c     input file will overwrite the values specified by the overall 
c     flag option. 
c
      ELSEIF (tag(1:3).EQ.'F12') THEN
        CALL ReadI(line,fc_target_calc_option,0,2,
     >       'Fluid Code Target Data Calculation Option')
c
c       Assign values to related sub-options
c
c       Base JET standard
c
        if (fc_target_calc_option.eq.0) then    
           fc_v_calc_opt  = 0
           fc_te_calc_opt = 1
           fc_ti_calc_opt = 2
           fc_ne_calc_opt = 2
c
c       Base UEDGE standard 
c
        elseif (fc_target_calc_option.eq.1) then
           fc_v_calc_opt  = 0
           fc_te_calc_opt = 0
           fc_ti_calc_opt = 0
           fc_ne_calc_opt = 2
c
c       Base Old Divimp/B2 standard
c
        elseif (fc_target_calc_option.eq.2) then
           fc_v_calc_opt  = 1
           fc_te_calc_opt = 2
           fc_ti_calc_opt = 2
           fc_ne_calc_opt = 2
c
c       Alternate B2/B2.5/Eirene
c
        elseif (fc_target_calc_option.eq.3) then
           fc_v_calc_opt  = 1
           fc_te_calc_opt = 1
           fc_ti_calc_opt = 1
           fc_ne_calc_opt = 2
        endif  
c
      ELSEIF (tag(1:3).EQ.'F13') THEN
        CALL ReadI(line,fc_ne_calc_opt,0,2,
     >       'Fluid Code Target Ne Data Calculation Option')
      ELSEIF (tag(1:3).EQ.'F14') THEN
        CALL ReadI(line,fc_te_calc_opt,0,2,
     >       'Fluid Code Target Te Data Calculation Option')
      ELSEIF (tag(1:3).EQ.'F15') THEN
        CALL ReadI(line,fc_ti_calc_opt,0,2,
     >       'Fluid Code Target Ti Data Calculation Option')
      ELSEIF (tag(1:3).EQ.'F16') THEN
        CALL ReadI(line,fc_v_calc_opt,0,1,
     >       'Fluid Code Target Vb Data Calculation Option')
      ELSEIF (tag(1:3).EQ.'F17') THEN
        CALL ReadI(line,fc_v_interp_opt,0,1,
     >       'Fluid Code Cell Edge Value Interpretation Option')
      ELSEIF (tag(1:3).EQ.'F18') THEN
        CALL ReadI(line,e2dformopt,0,2,
     >       'Fluid Code Plasma File Format Specifier')
c
c     For readaux = 3 (read auxiliary fluid code data option 3 - SOLPS4.3)
c     fort.44 file containing neutral densities.         
c
      ELSEIF (tag(1:3).EQ.'F19') THEN
        CALL ReadI(line,e2dneut_select,1,3,
     >       'Specify which block of impurity neutral data to read')
c
c
c     F20:
c
c     e2dion_select = 1 
c
c     The e2dnzs data stored in fort.31 can contain multiple fluid species
c     in addition to H+. This quantity specifies an offset into this data
c     so that the code can start reading the correct impurity into the e2d
c     fluid code arrays like e2dnzs. This is required in DIVIMP so it can
c     include meaningful comparisons between the DIVIMP and fluid code results.
c
      ELSEIF (tag(1:3).EQ.'F20') THEN
        CALL ReadI(line,e2dion_select,1,100,
     >       'Specify where to start reading the fluid code ion data')
c     
c
c -----------------------------------------------------------------------
c
c     TAG N21 : External sputtering flux data source
c               0 = geier file format for Ar
c               1 = import divimp charge resolved flux and energy
c                   data from a previous divimp run
c
      elseif(tag(1:3).eq.'N21') then 
        CALL ReadI(line,ext_flx_data_src,0,1,
     >       'External sputter flux data source option')
c
c -----------------------------------------------------------------------
c
c     TAG P60 : Density Gradient Option
c
      ELSEIF (tag(1:3).EQ.'P60') THEN
        CALL ReadI(line,ngradopt,0,1,'Density Gradient option')
c
c -----------------------------------------------------------------------
c
c     TAG P61 : Background Flow Velocity Override Option
c
      ELSEIF (tag(1:3).EQ.'P61') THEN
        CALL ReadI(line,override_bg_velocity_opt,0,14,
     >                           'Background Velocity Override Option')
c
c -----------------------------------------------------------------------
c
c     TAG P62 : Over-ride Efield target E-field calculation
c
      ELSEIF (tag(1:3).EQ.'P62') THEN
        CALL ReadI(line,ofield_targ,1,3,
     >                           'Override E-field target')
c
c -----------------------------------------------------------------------
c
c     TAG P63 : External plasma overlay option
c               0 = off
c               1 = on
c
      ELSEIF (tag(1:3).EQ.'P63') THEN
        CALL ReadI(line,external_plasma_overlay,0,1,
     >                           'External plasma overlay')
c
c -----------------------------------------------------------------------
c
c     TAG P64 : External plasma overlay file name
c               - specifies the name of the file to be loaded
c               - full path required unless rundiv script is modified
c
      ELSEIF (tag(1:3).EQ.'P64') THEN
        CALL ReadC(line,external_plasma_file,
     >                          'EXTERNAL PLASMA OVERLAY FILE NAME')
c
c -----------------------------------------------------------------------
c
c     TAG P65 : SOL option 12/13 etc - additional pressure option
c               PMULT - Adds additional pressure to PMULT * PINF
c               Over a distance of PDIST * SMAX
c
      ELSEIF (tag(1:3).EQ.'P65') THEN
        CALL ReadR(line,sol13_padd,0.0,HI,
     >             'SOL13+ ADDITIONAL PRESSURE')

c
c-----------------------------------------------------------------------
c     TAG P66 : SOL option 12/13 etc - additional pressure option
c               PDIST - Adds additional pressure to PMULT * PINF
c               Over a distance of PDIST * SMAX
c
      ELSEIF (tag(1:3).EQ.'P66') THEN
        CALL ReadR(line,sol13_pdist,0.0,HI,
     >             'SOL13+ ADDITIONAL PRESSURE DISTANCE')

c
c -----------------------------------------------------------------------
c
c     TAG Q42 : SHEATH TEMPERATURE - INNER JET/OUTER SONNET
c
      ELSEIF (tag(1:3).EQ.'Q42') THEN
c
c     Note: the tag line precedes a standard DIVIMP array input of
c           three lines.  
c
c     - specifies a temperature value to be used instead of the 
c       target value for the sheath energy calculations. 
c
c     INPUT IS:  IR  TE
c
        CALL RDRARN(sheath_vali,nsheath_vali,MAXNRS,
     >          -MACHHI,MACHHI,.FALSE.,
     >          -machhi,MACHHI,1,'INNER/OUTER SHEATH POTENTIAL',
     >          IERR)

c
c -----------------------------------------------------------------------
c
c     TAG Q43 : SHEATH TEMPERATURE - OUTER JET/INNER SONNET
c
      ELSEIF (tag(1:3).EQ.'Q43') THEN
c
c     Note: the tag line precedes a standard DIVIMP array input of
c           three lines.  
c
c     - specifies a temperature value to be used instead of the 
c       target value for the sheath energy calculations. 
c
c     INPUT IS:  IR  TE
c
        CALL RDRARN(sheath_valo,nsheath_valo,MAXNRS,
     >          -MACHHI,MACHHI,.FALSE.,
     >          -machhi,MACHHI,1,'OUTER/INNER SHEATH POTENTIAL',
     >          IERR)

c
c -----------------------------------------------------------------------
c
c     TAG Q44 : CORE PLASMA PROFILES AS A FUNCTION OF PSIN
c
      ELSEIF (tag(1:3).EQ.'Q44') THEN
c
c     Note: the tag line precedes a standard DIVIMP array input of
c           three lines.  
c
c     - specifies the core plasma data as a function of PSIN 
c       input line should be:  PSIN   Te    Ti    Ne    Vb 
c       which will be linearly interpolated and overwrite data in Q37 if any
c      
c     INPUT IS:  PSIN TE TI NE VB
c

        CALL RDRARN_ALLOC(coreprofile,ncoreprofile,
     >          -MACHHI,MACHHI,.FALSE.,
     >          -machhi,MACHHI,4,'CORE PLASMA PROFILES BY PSIN',
     >          IERR)

c -----------------------------------------------------------------------
c
c     TAG Q45 : DELTA PSI VALUE TO SHIFT THE INPUT CORE PROFILE
c
      ELSEIF (tag(1:3).EQ.'Q45') THEN
c
c     INPUT IS:  DELTA_PSIN_CORE
c
        CALL ReadR(line,delta_psin_core,-machhi,machhi,
     >             'PSIN SHIFT FOR CORE PROFILES')
c
c -----------------------------------------------------------------------
c
c     TAG R13 : SOL21 Additonal parameters - INNER JET/OUTER SONNET
c
      ELSEIF (tag(1:3).EQ.'R13') THEN
c
c     Note: the tag line precedes a standard DIVIMP array input of
c           three lines.  
c
c     - defines two additional points for linear fitting of N within 
c       region A. 
c
c     INPUT IS: IR L1A L1B NR1A NR1B TER1A TER1B TIR1A TIR1B
c
        CALL RDRARN(aux_s21parmi,aux_ns21i,MAXNRS,
     >          -MACHHI,MACHHI,.FALSE.,
     >          -machhi,MACHHI,8,'EXTRA SOL 21 INNER PARAMETERS',
     >          IERR)

c
c -----------------------------------------------------------------------
c
c     TAG R14 : SOL21 Additonal parameters - OUTER JET/INNER SONNET
c
      ELSEIF (tag(1:3).EQ.'R14') THEN
c
c     Note: the tag line precedes a standard DIVIMP array input of
c           three lines.  
c
c     - defines two additional points for linear fitting of N within 
c       region A. 
c
c     INPUT IS: IR L1A L1B NR1A NR1B TER1A TER1B TIR1A TIR1B
c
        CALL RDRARN(aux_s21parmo,aux_ns21o,MAXNRS,
     >          -MACHHI,MACHHI,.FALSE.,
     >          -machhi,MACHHI,8,'EXTRA SOL 21 OUTER PARAMETERS',
     >          IERR)
c
c -----------------------------------------------------------------------
c
c     TAG R15 : SOL21 extra velocity factor ... used to reduce flow
c               velocity near target. 
c              
c
      ELSEIF (tag(1:3).EQ.'R15') THEN
        CALL ReadR(line,aux_vel21,0.0,HI,
     >             'Modify near target velocity in SOL21')
c
c -----------------------------------------------------------------------
c
c     TAG S22 : Option to turn on V/A flag debugging in NEUT
c               OFF =0 = default
c               ON  >0
c
      ELSEIF (tag(1:3).EQ.'S22') THEN
        CALL ReadI(line,debug_neutv,0,1,'Neutral velocity debugging')
c
c -----------------------------------------------------------------------
c
c     TAG S23 : Maximum energy to be used for binning the
c               velocity/energy in the debugging option code
c
      ELSEIF (tag(1:3).EQ.'S23') THEN
        CALL ReadR(line,debug_neutv_einmax,0.0,HI,
     >             'Max EIN for NEUT V/A binning for debugging')
c -----------------------------------------------------------------------
c
c     TAG S24 : Number of bins to be used in Velocity debugging code
c               Maximum of 20000.
c
      ELSEIF (tag(1:3).EQ.'S24') THEN
        CALL ReadI(line,debug_neutv_nbins,1,20000,
     >             'Neutral velocity debugging - number of bins')
c
c -----------------------------------------------------------------------
c
c     TAGS T19 to T27
c
c     The following are all variables related to the REISER
c     force implementation. T10 to T16 allow for detailed control of the
c     active forces while T17 is the coulomb parameter and T18 is 
c     a debug option to linearize the gradients at the midplane in the
c     case of SOL option 7. 
c
c
      ELSEIF (tag(1:3).EQ.'T19') THEN
        CALL ReadI(line,aswitch,0,1,'Master Switch for force switches')
      ELSEIF (tag(1:3).EQ.'T20') THEN
        CALL ReadI(line,sk11,0,1,'FF switch')
      ELSEIF (tag(1:3).EQ.'T21') THEN
        CALL ReadI(line,sk12,0,1,'FIG switch')
      ELSEIF (tag(1:3).EQ.'T22') THEN
        CALL ReadI(line,sk13,0,1,'FVG switch')
      ELSEIF (tag(1:3).EQ.'T23') THEN
        CALL ReadI(line,sd11,0,1,'Maxwellian Velocity diffusion switch')
      ELSEIF (tag(1:3).EQ.'T24') THEN
        CALL ReadI(line,sd12,0,1,'TiGradient Velocity diffusion switch')
      ELSEIF (tag(1:3).EQ.'T25') THEN
        CALL ReadI(line,sd13,0,1,'V-Gradient Velocity diffusion switch')
      ELSEIF (tag(1:3).EQ.'T26') THEN
        CALL ReadR(line,coulomb_log,0.0,HI,'Coulomb Logarithm Value')
      ELSEIF (tag(1:3).EQ.'T27') THEN
        CALL ReadI(line,linearpeak,0,1,'Linear Peak Debug switch')
c
c
c -----------------------------------------------------------------------
c
c     TAG T28
c
c     T28 - pinch_loc_opt - this option specifies the region of the 
c                           grid where the radial velocity in
c                           pinchopt 4 will be applied.
c
c         = 0 = Entire grid excluding PFZ (default)
c         = 1 = main SOL only
c         = 2 = Entire grid including PFZ 
c               - this requires a sign change to the value assigned to 
c                 Vr (or Vpinch depending on terminology)
c         = 3 = main SOL past X-point region ONLY
c         = 4 = main SOL past X-point + core 
c                                 
      ELSEIF (tag(1:3).EQ.'T28') THEN
c
        CALL ReadI(line,pinch_loc_opt,0,4,'Grid region for pinch')
c
c -----------------------------------------------------------------------
c
c     TAG T29
c
c     T29 - pinch_npdf, pinch_pdf - this option loads the probability
c           distribution function to be used when randomly 
c           determining the value of the pinch/radial velocity at 
c           each time step. 
c
c     Note: the tag line precedes a standard DIVIMP array input of
c           three lines.  
c
      ELSEIF (tag(1:3).EQ.'T29') THEN
c
        CALL RDRARN(pinch_pdf,pinch_npdf,MAXPTS,
     >          -MACHHI,MACHHI,.TRUE.,
     >          -machhi,MACHHI,1,'RADIAL/PINCH V - PDF',
     >          ierr)
c
c       Calculate the integral value of the PDF and assign 
c       this to the normalization value. 
c
c       Note: the integration assumes that P=0 for V<Vmin-1/2*dVbin 
c             and V>Vmax+1/2*dVbin
c             so the implicit integration range is 
c             [Vmin-1/2dV,Vmax+1/2dV] of the 
c             input pdf values.
c       
c       The PDF is input in the form:
c
c           Vel    PDF-value
c
c       and the values are linearly interpolated between the points on 
c       the function to generate intermediate PRPBABILITIES.  
c
c       Calculate the integral and store it at each point in the data
c       - take the total and assign it to pdf_norm_val  
c
c       Set maximum and minimum velocities
c
        npdf_data = pinch_npdf+2
c
        pinch_pdf_data(1,1) = pinch_pdf(1,1) -
     >               (pinch_pdf(2,1)-pinch_pdf(1,1))/2.0
        pinch_pdf_data(1,2) = 0.0 
c
        pinch_pdf_data(npdf_data,1) = pinch_pdf(pinch_npdf,1)+ 
     >    (pinch_pdf(pinch_npdf,1)-pinch_pdf(pinch_npdf-1,1))/2.0
        pinch_pdf_data(npdf_data,2) = 0.0
c
c       Copy probability information to data array. 
c 
        do in = 1,pinch_npdf       
           pinch_pdf_data(in+1,1) = pinch_pdf(in,1)
           pinch_pdf_data(in+1,2) = pinch_pdf(in,2)
        end do 
c
c       Calculate integral
c 
c       Integration starts at zero 
c
        pdf_norm_val = 0.0
        pinch_pdf_data(1,3) = pdf_norm_val
c
        do in = 1,npdf_data-1
c
c          Define the velocity bin sizes based on input
c
c          Calculate the integral by assuming the probability
c          at each point follows the specified input rather 
c          than being constant for each dVbin. 
c
c           if (in.eq.1) then 
c              deltav2 = (pinch_pdf(in+1,1)-pinch_pdf(in,1))/2.0
c              deltav1 = deltav2 
c           elseif (in.eq.pinch_npdf) then 
c              deltav1 = (pinch_pdf(in,1)-pinch_pdf(in-1,1))/2.0
c              deltav2 = deltav1
c           else
c              deltav1 = (pinch_pdf(in,1)-pinch_pdf(in-1,1))/2.0
c              deltav2 = (pinch_pdf(in+1,1)-pinch_pdf(in,1))/2.0
c           endif
c
           deltav1 = pinch_pdf_data(in+1,1) - pinch_pdf_data(in,1) 
c
           pdf_norm_val = pdf_norm_val 
     >         + (pinch_pdf_data(in+1,2)+pinch_pdf_data(in,2))/2.0
     >            * deltav1
c
c          Accumulate integral at the bin centers
c          Integrated probability is zero in first cell. 
c
           pinch_pdf_data(in+1,3) = pdf_norm_val
c
        end do  
c
c       Normalize the integral of the PDF - used in random number selection
c
        if (pdf_norm_val.ne.0.0) then 
c
           do in = 1,npdf_data
c
              pinch_pdf_data(in,3) = pinch_pdf_data(in,3)/pdf_norm_val
c
              write(6,'(a,i4,4(1x,f14.8))') 'PINCH INT:',in,
     >          pinch_pdf_data(in,1),pinch_pdf_data(in,2),
     >          pinch_pdf_data(in,3), pdf_norm_val
c 
           end do
c
       endif  
c
c      Print out vr_pdf_int
c
c       do in = -40,50
c          vtest = 10.0 * in
c          res   = vr_pdf_int(vtest,-1)
c          write(6,'(a,i4,3(1x,f14.8))') 'PDF_INT:',in,
c     >       vtest,res,res/pdf_norm_val
c       end do   
c
c
c -----------------------------------------------------------------------
c
c     TAG T30
c
c     T30 - pinch correlation time. A new radial velocity value will be 
c           chosen periodically based on the value of this quantity.
c           The default value of 0.0 will result in a new velocity 
c           being selected every time step. This value is specified 
c           in second. On JET it is typically 5 to 20 microseconds.
c                                 
      ELSEIF (tag(1:3).EQ.'T30') THEN
c
        CALL ReadR(line,pinch_correlation_time,0.0,HI,
     >                                 'Pinch Correlation Time')      


c 
c -----------------------------------------------------------------------
c
c     TAG T31
c
c     T31 - Drift region - specifies the region to which poloidal 
c           drifts should be applied. 
c           1 - SOL + PFZ
c           2 - SOL only
c           3 - PFZ only
c           4 - CORE only
c
c           Other options can easily be added as needed - the default is
c           option 1. 
c
      ELSEIF (tag(1:3).EQ.'T31') THEN
c
        CALL ReadI(line,drft_region,0,4,'Region over which to apply'//
     >                                   'poloidal drifts')
c
c -----------------------------------------------------------------------
c
c     TAG T32
c
c     T32 - Drift Mach Option - Detailed drift velocity input on a ring
c           ring basis is specified as a mach number to be multiplied
c           by the sound speed at the top of the torus for each ring
c
c           Option 0 : OFF
c           Option 1 : CS calculated from 2*Te
c           Option 2 : CS calculated from Te+Ti
c         
c           Default value is 0 - OFF - data is specified in terms of 
c                                      velocity
c
      ELSEIF (tag(1:3).EQ.'T32') THEN

        CALL ReadI(line,drftvel_machopt,0,2,'Drift velocity specified'//
     >                                   'by mach values')
c
c -----------------------------------------------------------------------
c
c     TAG T33
c
c     T33 - Detailed specifications of data ring by ring - this is 
c           an array listing 
c                      ring number       velocity/mach
c           Data does not need to be specified for each ring - the 
c           default value will be applied instead.
c
c           The number of array elements is initialized to zero  
c
      ELSEIF (tag(1:3).EQ.'T33') THEN
c
c        Read in array data
c
        CALL RDRARN(ringdrftvel,ndrftvel,MAXNRS,
     >          real(1),real(maxnrs),.TRUE.,
     >          -machhi,MACHHI,1,'DETAILED RING DRIFT VELOCITY',
     >          ierr)

c
c -----------------------------------------------------------------------
c
c     TAG T34
c
c     T34 - S displacement in 2D resulting from a perpendicular step in
c           paramagnetic "Z" direction in 3D - this actually moves the 
c           particle onto an adjacent flux tube - however, since 
c           DIVIMP is 2D - the effect is to actually move the particle
c           onto an adjacent identical flux tube at a different S location
c           - thus effectively giving a net S displacement.
c
c           The first approximation to this is to use
c
c           ds = cross_step * Btor/Bpol
c
c           A value of 0 for this option is OFF 
c                      1 is ON

      ELSEIF (tag(1:3).EQ.'T34') THEN
        CALL ReadI(line,dperpz_opt,0,1,'3D Dperp Delts S Option')
c
c
c -----------------------------------------------------------------------
c
c     TAG T35 - related to poloidal drift options - T31,T32,T33
c
c     T35 - Drift region calculation option
c
c     Option 0: Input values are specified in terms of S
c     Option 1: Input values are specified in terms of P (poloidal distance)
c     Option 2: Input values are specified in terms of Z (single null only)
c
c     Default value is S para (option 0)
c
      ELSEIF (tag(1:3).EQ.'T35') THEN
c
c     jdemod - I don't really understand commenting this out 
c            - it defaults to a value of 0 so this will only stop the code 
c              if the T35 line is in the input file
c            - other than that I think it works as intended so I am uncommenting
c              it .. its only purpose is to change the interpretation of other 
c              inputs
c
c slmod begin - *** TEMP ***
        CALL ReadI(line,drft_distopt,0,2,'Drift velocity range'//
     >                                   ' specification (S,P or Z)')
c
c        STOP 'OPTION TURNED OFF FOR NOW...'
c slmod end
c
c        write(0,*) 'READIN: drft_distopt:',drft_distopt
c
c
c
c -----------------------------------------------------------------------
c
c     TAG T36 to T39 - options related to the implementation of 
c                      impurity exb drifts 
c
c     TAG T36 - potopt
c
c             - This option is used to determine the method of calculatng
c               plasma potential
c             - potopt = 0    Use 3xTe(0) at each target as the floating 
c                              potential start start point
c             - potopt = 1    Import LP data listing the measured floating 
c                              potential ... if imported data not available it 
c                              defaults to option 0.(not yet implemented) 
c
      ELSEIF (tag(1:3).EQ.'T36') THEN
        CALL ReadI(line,potopt,0,0,'Option for calculating the floating'
     >             //' potential')
c
c
c     TAG 37
c
c     exb_rad_opt = 0 ... no exb radial drift is applied
c                 = 1 ... exb radial drift is turned on
c
      ELSEIF (tag(1:3).EQ.'T37') THEN
        CALL ReadI(line,exb_rad_opt,0,1,'ExB radial drift option')
c
c
c     TAG 38
c
c     exb_pol_opt = 0 ... no exb poloidal drift is applied
c                 = 1 ... exb poloidal drift is turned on
c
      ELSEIF (tag(1:3).EQ.'T38') THEN
        CALL ReadI(line,exb_pol_opt,0,1,'ExB poloidal drift option')
c
c     TAG 39
c
c     exb_scale - real number
c               - the basic function of this is to switch the sign of 
c                 the ExB drift for cases of forward (+1.0) and reverse (-1.0)
c                 B-field orientation. However, it can also be used as a scaling
c                 factor if the drifts are found to be either too large or too 
c                 small. 
c
      ELSEIF (tag(1:3).EQ.'T39') THEN
        CALL ReadR(line,exb_scale,-HI,HI,
     >                    'ExB scaling factor (usually +/-1.0)')


c
c -----------------------------------------------------------------------
c
c     Force scaling factors - all default to 1.0
c     T40 to T44
c     T40 = Friction force scaling factor (SF_FRIC)
c     T41 = Ion temperature force scaling factor (SF_TI)
c     T42 = Electron temperature force scaling factor (SF_TE)   
c     T43 = Electric field force scaling factor (SF_EF)
c     T44 = Velocity diffusion scaling factor (SF_VDIFF)
c     T45 = Scaling factor for TAU (sf_tau)
c     
c     Defaults:
c     sf_fric = 1.0
c     sf_ti   = 1.0
c     sf_te   = 1.0
c     sf_ef   = 1.0
c     sf_vdiff= 1.0
c     sf_tau  = 1.0        
c     
c-----------------------------------------------------------------------
c       
c
      ELSEIF (tag(1:3).EQ.'T40') THEN
        CALL ReadR(line,sf_fric,-HI,HI,
     >                    'Friction force scaling factor (def=1.0)')
      ELSEIF (tag(1:3).EQ.'T41') THEN
        CALL ReadR(line,sf_ti,-HI,HI,
     >                    'Ti gradient force scaling factor (def=1.0)')
      ELSEIF (tag(1:3).EQ.'T42') THEN
        CALL ReadR(line,sf_te,-HI,HI,
     >                    'Te gradient force scaling factor (def=1.0)')
      ELSEIF (tag(1:3).EQ.'T43') THEN
        CALL ReadR(line,sf_ef,-HI,HI,
     >                 'Electric field force scaling factor (def=1.0)')
      ELSEIF (tag(1:3).EQ.'T44') THEN
        CALL ReadR(line,sf_vdiff,-HI,HI,
     >                   'Velocity diffusion scaling factor (def=1.0)') 
      ELSEIF (tag(1:3).EQ.'T45') THEN
        CALL ReadR(line,sf_tau,-HI,HI,
     >                    'TAU scaling factor (def=1.0)')

c        
c -----------------------------------------------------------------------
c
c    T46 Velocity based temperature calculation option
c
c    0  = default
c    1+ = other options
c     
      ELSEIF (tag(1:3).EQ.'T46') THEN
        CALL ReadI(line,ti_calc_opt,0,3,'Impurity Ti Calculation Opt')
c        
c -----------------------------------------------------------------------
c
c    T47 Coulomb logarithm calculation options
c 
c     0  = default = constant (default value = 15.0)
c     1  = 30.0 - 0.5 * LOG(ni) + 1.5 * LOG(ti)  [HC code - ]
c          Originally in Sivukhin, D.V., Coulomb collisions in a fully ionized plasma in
c          Review of Plasma Physics (Consultation Bureau, New York, 1966) Vol. 4, p.88.
c
c     2  = 17.3 - 0.5*LOG(n/1.0E20) + 1.5*LOG(t/1000.0)  [LIM code]
c     3  = log(1.5e13 * t**(1.5) / sqrt(n))   [SOL22 PEI term]
c    
      ELSEIF (tag(1:3).EQ.'T47') THEN
        CALL ReadI(line,lambda_opt,0,3,'Coulomb logarithm calc opt')
c
c    T48 Coulomb logarithm calculation options
c     Coulomb logarithm constant value - default value is 15.0 - this allows
c     specification of alternate constant values for option 0.         
c
c        
      ELSEIF (tag(1:3).EQ.'T48') THEN
        CALL ReadR(line,lambda_val,0.0,HI,'Coulomb logarithm const val')
        
    
      ! T49 Blob frequency for pinch velocity.
      !
      ! T29 is the PDF of the radial velocities. An additional option
      ! is to specify the rate of blobs for the entire plasma. This is 
      ! to effect that the transport will be intermittent. I.e., instead
      ! of always sampling the velocity PDF at every step, sample it
      ! proportional to the amount of blobs it would see. For instance,
      ! if fblob = 1e6 and qtim = 1e-8, then then probability of choosing
      ! a velocity from the PDF is:
      ! prob_choosing = fblob * qtim = 1e-2. 
      ! In words, the ion sees 0.01 blobs every time step. So the 
      ! probability of choosing from the PDF is 1%, which is easily done.
      ! If fblob = 0.0 this option has no effect.
      
      elseif (tag(1:3).eq.'T49') then
        call readr(line, fblob, -1.0, HI, 'Blob frequency')
        
      ! T50 Core diffusion coefficient. 
      ! -1.0 = Same as CDPERP.
      elseif (tag(1:3).eq.'T50') then
        call readr(line, cdperpc,-1.0,HI, 'Core diffusion coefficient')
        
      ! T51-53 are free to use. They were old options I scrapped.
     
      elseif (tag(1:3).eq.'T54') then
        call readi(line, balloon_opt, 0, 1, 
     >   'Ballooning transport approx. switch')
     
      ! T55  Divertor radial velocity factor for PDF option. If the ion
      ! is below/above the X point (LSN/USN), then multiply the radial
      ! velocity when chosen by this value. This is to simulate lower
      ! (<1) or higher (>1) radial transport in the divertor. Experiments
      ! have shown lower before.
      elseif (tag(1:3).eq.'T55') then
        call readr(line, div_vr_fact, -HI, HI, 
     >   'Divertor radial velocity factor')
      
      elseif (tag(1:3).eq.'T56') then
        call readi(line, in_blob_switch, 0, 1, 
     >    'Turn off parallel transport when in blob switch')
     
      ! T57 T58
      ! Options for checking inward moving hole-like transport near 
      ! separatrix.
      elseif (tag(1:3).eq.'T57') then
        call readi(line, hole_switch, 0, 1, 
     >    'Check for holes near separatrix switch')
      elseif (tag(1:3).eq.'T58') then
        call readr(line, hole_lambda, 0.0, HI, 
     >   'Hole frequency decay length')
     
      ! T59 Core pinch value.
      elseif (tag(1:3).eq.'T59') then
        call readr(line, core_pinch, -HI, HI, 
     >   'Core pinch value')
     
      ! T60 Minium psin for blob-like impurity transport model.
      elseif (tag(1:3).eq.'T60') then
        call readr(line, blob_min_rmrsomp, -HI, HI, 
     >   'Minimum R-Rsep @ OMP for blob-like transport model')
     
      ! T61 Birth location of blob/holes.
      elseif (tag(1:3).eq.'T61') then
        call readr(line, blob_birth_rmrsomp, -HI, HI, 
     >   'R-Rsep @ OMP for blob/hole birth')
     
      ! T62 exponential decay length for blob frequency.
      elseif (tag(1:3).eq.'T62') then
        call readr(line, blob_lambda, -HI, HI, 
     >   'Blob frequency decay length')
      
       
c        
c
c -----------------------------------------------------------------------
c
c
c -----------------------------------------------------------------------
c 
c     TAGS W01 and W02
c
c     Options for walls
c
c     wall_plasma_opt - option to select algorithm to calculate 
c     plasma conditions associated with wall elements
c
      ELSEIF (tag(1:3).EQ.'W01') THEN
        CALL ReadI(line,wall_plasma_opt,0,2,'Wall Plasma Option')
c
c     wall_plasma_fact - scaling factor used by wall_plasma_opt 
c     algorithms
c
      ELSEIF (tag(1:3).EQ.'W02') THEN
        CALL ReadR(line,wall_plasma_fact,0.0,HI,
     >                   'Wall Plasma Scaling Factor')
c         
      ELSE
          CALL ER('ReadUnstructuredInput','Unrecognized tag',*99)
      ENDIF


c...  Check tag list:
      DO i1 = 1, ntaglist-1
        IF (tag.EQ.taglist(i1)) THEN
          CALL ER('ReadUnstructuredInput','Duplicate tag detected',*10)
 10       WRITE(0,*) 'TAG: "'//tag//'"'
          WRITE(0,*) 'PROGRAM STOP'
          STOP
        ENDIF
      ENDDO

      RETURN
c
c There is an error:
c
99    WRITE(SLOUT,'(5X,3A)') 'LINE = "',line,'"'
      WRITE(SLOUT,'(5X,3A)') 'TAG  = "',tag ,'"'
      WRITE(0    ,'(5X,3A)') 'LINE = "',line(1:LEN_TRIM(line)),'"'
      WRITE(0    ,'(5X,3A)') 'TAG  = "',tag ,'"'
      WRITE(0,*) '    DIVIMP HALTED'
      STOP
      END
c
c ======================================================================
c
      REAL FUNCTION GetInputVersion(line)
      IMPLICIT none

      CHARACTER*(*) line
     
      INTEGER i,j,k

      GetInputVersion = -1.0

c...  Find location of first 2 spaces on the line:
      i = 0
      j = 0
      DO k = 1, LEN_TRIM(line)
        IF (line(k:k).EQ.' '.AND.i.NE.0) j = k
        IF (line(k:k).EQ.' '.AND.i.EQ.0) i = k
        IF (j.NE.0) EXIT
      ENDDO
      
      IF (i.NE.0.AND.k.NE.0.AND.j-1.GT.1) 
     .  READ(line(i+1:j-1),*) GetInputVersion
 
      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: ReadTagSeries_G
c
c Unstructure input related to the magnetic grid are loaded.  The input
c tag starts with 'G'.
c
      SUBROUTINE ReadTagSeries_G(tag,line,fp)
      use subgrid_options
      use ribbon_grid_options
      use mod_params
      use mod_slcom
      use mod_comtor
      use mod_cgeom
      use mod_io
      IMPLICIT none
c
c     READ "G" Series Unstructured input
c

      INTEGER   fp
      CHARACTER line*(*),tag*3

c     INCLUDE 'params'
c     INCLUDE 'slcom'
c     INCLUDE 'comtor' 
c     INCLUDE 'cgeom'  

      INTEGER ierr,i1,i2,ir

      IF     (TAG(1:3).EQ.'G01') THEN
c...    Not sure if Dave has reserved G01 already, and the web server is down:
        CALL ReadR(line,grd_thresh,0.0,100.0,'Refinement threshold')
c
c -----------------------------------------------------------------------
c
c     TAG G23: SONNET grid SUB type option
c
c     This option is used to tag a SONNET style grid as being an 
c     FRC (Field Reversed Configuration) custom grid. At the moment 
c     only one option is supported but further subtypes could be 
c     added for Sonnet grids requiring special processing. 
c
c     0 = stanard Sonnet grid = Default 
c     1 = FRC version 1 - type of Sonnet grid 
c         - used to set various FRC related options
c     2 = sonnet grid without boundary cells - boundary cells are added
c         - useful for carre grids
c
      ELSEIF (tag(1:3).EQ.'G23') THEN
        CALL ReadI(line,sonnet_grid_sub_type,0,2,'SONNET Grid '//
     >                                     'SUB-Type Option')
c
c -----------------------------------------------------------------------
c
c     TAG G34 : Machine Type
c
      ELSEIF (tag(1:3).EQ.'G34') THEN
c
c       0 = jet
c       1 = diiid
c       2 = alcator (cmod)
c       3 = aug (asdex upgrade)
c       4 = iter
c       5 = ignitor
c       6 = fire
c
        CALL ReadI(line,tmachine_opt,0,6,'MACHINE TYPE option')
c
c
c -----------------------------------------------------------------------
c
c     TAG G35 : Shot Identification string for JET cataloguing system 
c
      ELSEIF (tag(1:3).EQ.'G35') THEN
c
c       Shot identification string for cataloguing on JET
c
        CALL ReadC(line,divshotid,'OPTIONAL SHOT ID')
c
c -----------------------------------------------------------------------
c
c     TAG G36 : Parallel Ion Reflection option
c
      ELSEIF (tag(1:3).EQ.'G36') THEN
c
c       Option to activate parallel ion reflection
c
        CALL ReadI(line,s_reflect_opt,0,1,'PARALLEL ION REFLECTION OPT')
c
c -----------------------------------------------------------------------
c
      ELSEIF (TAG(1:3).EQ.'G37') THEN
c...    Data for modifying (cutting/extending) the magnetic grid:    
c       
        IF (line(7:9).EQ.'1.0') THEN
c...      Inputs:
c     
c         Type - 1: Cut ring to specified line segment.
c                2: Extend end of ring to specified line segment.
c         Mode - 1: Work from low index target (and remove portion from cut to outer target if
c                   necessary).
c                2: Work from high index target (and remove portion from cut to outer target if
c                   necessary)
c         Refinement - 0: none
c         Range1 - Start of IR index range of rings to be cut (inclusive).
c         Range2 - End of IR index range of rings to be cut.
c
          CALL RdRarn(grdmod,grdnmod,1000,-MACHHI,MACHHI,.FALSE.,
     .              -MACHHI,MACHHI,4,'Magnetic grid modification',ierr)
          IF (ierr.NE.0) CALL ER('ReadTagSeriesG','Error GRDMOD',*99)

          WRITE(fp,*)
          WRITE(fp,'(A)') 'Grid modification data:'
          DO i1 = 1, grdnmod
            WRITE(fp,'(5X,5F9.5)') (grdmod(i1,i2),i2=1,5)
          ENDDO

        ENDIF
c
c -----------------------------------------------------------------------
c
c     The following tags are related to the subgrid option for recording 
c     more detailed data on a finer grid
c
c     G38: Base subgrid option ON/OFF
c     G39: R,Z dimensions of the grid region
c     G40: RMIN,RMAX of the subgrid region
c     G41: ZMIN,ZMAX of the subgrid region
c
c -----------------------------------------------------------------------
c
c     TAG G38 : Subgrid option
c
      ELSEIF (tag(1:3).EQ.'G38') THEN
c
c     G38: Option to activate the subgrid data collection
c
        CALL ReadI(line,subgrid_opt,0,1,'BASE SUBGRID OPTION')
c
      ELSEIF (tag(1:3).EQ.'G39') THEN
c
c     G39: R,Z dimensions of the grid region
c
        CALL Read2I(line,sg_rdim,sg_zdim,1,500,'SUBGRID RDIM,ZDIM')
c
      ELSEIF (tag(1:3).EQ.'G40') THEN
c
c     G40: RMIN,RMAX of the subgrid region
c
        CALL Read2R(line,sg_rmin,sg_rmax,-HI,HI,'SUBGRID RMIN,RMAX')
c
      ELSEIF (tag(1:3).EQ.'G41') THEN
c
c     G41: ZMIN,ZMAX of the subgrid region
c
        CALL Read2R(line,sg_zmin,sg_zmax,-HI,HI,'SUBGRID ZMIN,ZMAX')
c
c -----------------------------------------------------------------------
c
c     Options related to ribbon grids
c     G42 - grid generation option - <i4>
c     G43 - intersection point averaging option - opt_block_av - <r4>
c     G44 - maximum R separation in grid generator - max_r_sep - <r4>
c     G45 - maximum S/Z separation in grid generator - max_s_sep - <r4>
c     G46 - min number of cells on ring - min_cells - <i4>
c     G47 - castem output identifier - <string>
c     G48 - min and max S for selecting intersection subset  2 x <r4>
c     G49 - min and max R for intersection subset grid generation 2 x <r4>
c     G50 - min and max S for intersection subset grid generation 2 x <r4>
c
      ELSEIF (tag(1:3).EQ.'G42') THEN
c
c     G42: Ribbon grid option  0=unstructured  1=structured
c
        CALL ReadI(line,rg_grid_opt,0,1,'RIBBON GRID OPTION')
c
c     G43: Ribbon grid option  1=unstructured  2=structured
c
      ELSEIF (tag(1:3).EQ.'G43') THEN
        CALL ReadI(line,rg_block_av,0,1,'BLOCK AVERAGE OPTION')
c
c     G44: Maximum row separation
c
      ElseIf (tag(1:3).eq.'G44') Then
        Call ReadR(line,rg_max_r_sep,0.0,HI,
     >                       'Maximum Row separation (m)')
c
c     G45: Maximum cell separation
c
      ElseIf (tag(1:3).eq.'G45') Then
        Call ReadR(line,rg_max_s_sep,0.0,HI,
     >                       'Maximum cell separation (m)')
c
c     G46: Minimum number of cells in a row
c
      ElseIf (tag(1:3).eq.'G46') Then
        CALL ReadI(line,rg_min_cells,1,9999,'Minimum cells in a row')

c
c     G47: Castem data set identifier
c
      ElseIf (tag(1:3).eq.'G47') Then
        CALL ReadC(line,rg_castem_data,'CASTEM DATA SET IDENTIFIER')
        write(0,*) 'RG_CASTEM_DATA:',trim(rg_castem_data),':'
c
c     G48 - min and max S for selecting intersection subset  2 x <r4>
c
      ELSEIF (tag(1:3).EQ.'G48') THEN
        CALL Read2R(line,rg_int_win_mins,rg_int_win_maxs,-HI,HI,
     >             'RIBBON GRID INTERSECTION SUBSET RANGE [S1,S2]')
c
c     G49 - min and max S for selecting intersection subset  2 x <r4>
c
      ELSEIF (tag(1:3).EQ.'G49') THEN
        CALL Read2R(line,rg_minr,rg_maxr,-HI,HI,
     >             'RIBBON GRID SUBSET R RANGE [R1,R2]')
c
c     G50 - min and max S for intersection subset grid generation 2 x <r4>
c
      ELSEIF (tag(1:3).EQ.'G50') THEN
        CALL Read2R(line,rg_mins,rg_maxs,-HI,HI,
     >             'RIBBON GRID SUBSET S RANGE [S1,S2]')
c
c     G51 - cutoff factor for ring generation - rings
c           in a ribbon grid with a length factor smaller than 
c           this value will not be generated. 
c
      ELSEIF (tag(1:3).EQ.'G51') THEN
        CALL ReadR(line,lcutoff,-HI,HI,
     >             'RING CUTOFF LENGTH FACTOR')
c
c     G52 - Cell spacing option ... option to calculate the cell
c           boundary spacing along a ring. Only option 0 is currently
c           available which uses an exponential factor given in G53. 
c           A cell_spacing_factor of 1 gives a linear spacing
c
      ELSEIF (tag(1:3).EQ.'G52') THEN
        CALL ReadI(line,cell_spacing_option,0,0,'CELL SPACING OPTION')

c
c     G53 - cell spacing factor
c           Used to determine cell boundary spacing along the rings
c
      ELSEIF (tag(1:3).EQ.'G53') THEN
       CALL ReadR(line,cell_spacing_factor,-HI,HI,
     >             'CELL SPACING FACTOR')
c
c     G54 - Cell spacing option ... option to calculate the cell
c           boundary spacing along a ring. Only option 0 is currently
c           available which uses an exponential factor given in G53. 
c           A cell_spacing_factor of 1 gives a linear spacing
c
      ELSEIF (tag(1:3).EQ.'G54') THEN
        CALL ReadI(line,ribbon_input_format_opt,0,1,
     >     'INTERSECTION DATA INPUT DATA FORMAT OPTION: 0=CASTEM 1=RAY')
c
c
c -----------------------------------------------------------------------
      ELSEIF (tag(1:3).EQ.'G55') THEN
c
c     G55 - read in a list of ikoffsets to move the center of the 
c           ring for background plasma calculation
c
c           n_ikoffsets
c
c           option   ir1   ir2    offset_value
c
c           option=0 offset = ik_mids + offset_value (index offset)
c           option=1 offset = offset_value * SMAX    (fractional offset)
c           option=2 offset = offset_value * PMAX    (fractional offset)
c  
c
          CALL RDRARN(ik_offset_data,n_ik_offsets,MAXNRS,
     .                -MACHHI,MACHHI,.FALSE.,
     .                -MACHHI,MACHHI,3,'IK OFFSET DATA',IERR)
c
c
      ELSEIF (tag(1:3).EQ.'G56') THEN
c
c     jdemod - sol22_halfringlen_opt = 0 = ringlen/2
c                                    = 1 = ksb(ikmid,ir)
c            - the ikmid value should be at ~ 1/2 the ring length     
c            - note: using option 1 will support calculations that 
c                    move the midpoint
c
        CALL ReadI(line,sol22_halfringlen_opt,0,1,
     >     'SOL22 half ring length determination option')
c
      ELSE
        CALL ER('ReadTagSeriesG','Unrecognized tag',*99)
      ENDIF

      RETURN
c.... Error:
99    WRITE(SLOUT,'(5X,3A)') 'LINE = "',line,'"'
      WRITE(SLOUT,'(5X,3A)') 'TAG  = "',tag ,'"'
      WRITE(0    ,'(5X,3A)') 'LINE = "',line(1:LEN_TRIM(line)),'"'
      WRITE(0    ,'(5X,3A)') 'TAG  = "',tag ,'"'
      WRITE(0,*) '    DIVIMP HALTED'
      STOP
      END
c
c ======================================================================
c
c subroutine: ReadTagSeries_H
c
c Unstructured input data tagged with an H are loaded. The majority of 
c these refer to hydrocarbon code modeling options.
c
      SUBROUTINE ReadTagSeries_H(line,tag,fp)
! ammod begin.
      Use comhc            ! Assign values to Hydrocarbon following common block.
      use hc_kinetics_options
c
c      Use HC_Init_DIV_Data ! HC setup routines.
c
c     Params is included in HC_com and so is not required separately 
c
c      INCLUDE 'params'
! ammod end
      use mod_comtor
      use mod_slcom
      use mod_hc_global_opts
      use mod_io
      IMPLICIT none
c
c     READ "H" Series Unstructured input
c

      INTEGER   fp
      CHARACTER line*(*),tag*3

c     include 'comtor'
c     include 'slcom'
c     include 'hc_global_opts'


! ammod begin.
! ----------------------------
!
! Options added for hydrocarbon following.
!
! ----------------------------

      If (Tag(1:3).eq.'H15') Then
        Call ReadI(line,hc_follow_option,0,1,
     >     'Hydrocarbon following option, 0-off, 1-on')
c
c       Set global_hc_follow_option to match hc value read in
c       Initialization of global_hc_follow_option occurs in the setup 
c       routine where initialization of the hc_follow_option also takes place
c
        global_hc_follow_option = hc_follow_option 
c
      ElseIf (Tag(1:3).eq.'H16') Then
        Call ReadI(line,hc_higher_hcs_option,0,1,
     >     'Follow higher hydrocarbon (C2+) option, 0-off, 1-on')
      ElseIf (Tag(1:3).eq.'H17') Then
        Call ReadI(line,hc_wbc_comp_option,0,1,
     >     'WBC comparison case, 0-off, 1-on')

      ElseIf (Tag(1:3).eq.'H20') Then
        Call ReadI(line,hc_sputtering_model,0,1,
     >     'Model for sputtered species release, 0-preset, 1-Mech')
      ElseIf (Tag(1:3).eq.'H21') Then
        Call ReadI(line,hc_sputtered_hc_species,0,58,
     >     'Preset sputtered hydrocarbon species, 10-Methane(CH4)')
      ElseIf (Tag(1:3).eq.'H22') Then
        Call ReadI(line,hc_evolution_model_primary,1,3,
     >     'Model for HC data primary, 1-E&L, 2-Brooks, 3-Janev')
      ElseIf (Tag(1:3).eq.'H23') Then
        Call ReadI(line,hc_evolution_model_secondary,0,3,
     >     'Model for HC data secondary, 0-none, 1-E&L, 2-Brooks,
     > 3-Janev')
c     Available libraries 1) Ehrhardt and Langer (PPPL, 1987)
c                         2) Alman, Ruzic, Brooks (Phys. Plasmas, 2000)
c                         3) Janev, et al. (NIFS, 2001)
      ElseIf (Tag(1:3).eq.'H24') Then
        Call ReadI(line,hc_launch_location,-1,6,
     >     'Model for HC launch location (same options as CNEUTB)')
c
c        Values assigned in global_hc_assign_inputs after the 
c        entire input file has been read in
c
c	If (hc_launch_location .eq. -1) Then
c		hc_launch_location = CNEUTB
c	End If
      ElseIf (Tag(1:3).eq.'H25') Then
        Call ReadI(line,hc_launch_angle_velocity,-1,20,
     >    'Model for HC launch angle/velocity (same options as CNEUTC)')
c
c        Values assigned in global_hc_assign_inputs after the 
c        entire input file has been read in
c
c	If (hc_launch_angle_velocity .eq. -1) Then
c		hc_launch_angle_velocity = CNEUTC
c	End If
      ElseIf (Tag(1:3).eq.'H26') Then
        Call ReadI(line,hc_launch_velocity_model,0,2,
     >     'Launch velocity model, 0-const, 1-MB dist., 2-dual MB')     
      ElseIf (Tag(1:3).eq.'H27') Then
        Call ReadR(line,hc_dual_mb_pri_vel_flux,0.0,1.0,
     >     'Dual MB velocity flux primary MB contrib., 0.0-1.0')
      ElseIf (Tag(1:3).eq.'H28') Then
        Call ReadR(line,hc_dual_mb_sec_mean_temp,0.0,2000.0,
     >     'Dual MB launch velocity T2 (deg K), 0.0-2000.0')
	
      ElseIf (Tag(1:3).eq.'H30') Then
        Call ReadI(line,hc_neut_ion_velocity,-1,3,
     >     'Neutral->Ion initial velocity (same options as CNEUTG)')
c
c        Values assigned in global_hc_assign_inputs after the 
c        entire input file has been read in
c
c       If (hc_neut_ion_velocity .eq. -1) Then
c		hc_neut_ion_velocity = CNEUTG
c	End If
      ElseIf (Tag(1:3).eq.'H31') Then
        Call ReadI(line,hc_ion_neut_angle,0,2,
     >     'Ion->neutral angle emission, 0-isotropic, 1-sine, 2-S dir')
      ElseIf (Tag(1:3).eq.'H32') Then
        Call ReadI(line,hc_ion_neut_velocity,0,1,
     >     'Ion->neutral velocity, 0-ion energy, 1-')
      ElseIf (Tag(1:3).eq.'H33') Then
        Call ReadI(line,hc_lambda_calc,0,1,
     >     'Improved calculation for lambda (Sivukhin),0-off,1-on')
      ElseIf (Tag(1:3).eq.'H34') Then
        Call ReadI(line,hc_disable_transitions,0,1,
     >     'Disable HC transitions, 0-off,1-on')
      ElseIf (Tag(1:3).eq.'H35') Then
        Call ReadI(line,hc_presheath_efield,0,1,
     >     'Improved model for electric field force,0-off,1-on')
      ElseIf (Tag(1:3).eq.'H36') Then
        Call ReadR(line,hc_efield_drop_fraction,0.0,1.0,
     >     'Fraction of potential drop in Debye region,0.0-1.0')
      ElseIf (Tag(1:3).eq.'H37') Then
        Call ReadI(line,hc_efield_cells,0,5,
     >     'Cells from target to apply improved e-field, 0-5')

      ElseIf (Tag(1:3).eq.'H40') Then
        Call ReadI(line,hc_neutral_reflection_option,0,1,
     >     'Neutral HC reflection switch')
      ElseIf (Tag(1:3).eq.'H41') Then
        Call ReadI(line,hc_ion_reflection_option,0,1,
     >     'Ion HC reflection switch')
      ElseIf (Tag(1:3).eq.'H42') Then
        Call ReadI(line,hc_reflection_coef_model,0,4,
     >     'Reflection model, 0-preset, 1-Janev, 2-A&R, 3-CH4=1.0')
      ElseIf (Tag(1:3).eq.'H43') Then
        Call ReadR(line,hc_reflection_coef_preset,0.0,1.0,
     >     'Preset reflection coef')
      ElseIf (Tag(1:3).eq.'H44') Then
        Call ReadI(line,hc_reflection_species_model,0,1,
     >     'Reflected species model, 0-preset, 1-Alman&Ruzic')
      ElseIf (Tag(1:3).eq.'H45') Then
        Call ReadI(line,hc_reflection_energy_model,0,3,
     >     'Reflection energy model,0-set,1-impact,2-thermal,3-AR')
      ElseIf (Tag(1:3).eq.'H46') Then
        Call ReadR(line,hc_refl_energy_neutral_preset,0.0,1E4,
     >     'Preset reflecting particle energy, neutral impact (eV)')
      ElseIf (Tag(1:3).eq.'H47') Then
        Call ReadR(line,hc_refl_energy_ion_preset,0.0,1E4,
     >     'Preset reflecting particle energy, ion impact (eV)')
      ElseIf (Tag(1:3).eq.'H48') Then
        Call ReadI(line,hc_reflection_angle_model,-1,11,
     >     'Refl. angle model,-1=NRFOPT,0-off,1-spec,2-isotr,'
     >     // '3-norm,4-A&R ')
c
c        Values assigned in global_hc_assign_inputs after the 
c        entire input file has been read in
c
c        If (hc_reflection_angle_model .eq. -1) Then
c		hc_reflection_angle_model = NRFOPT
c        End If

      ElseIf (Tag(1:3).eq.'H50') Then
        Call ReadI(line,hc_sputtering_option,0,1,
     >     'HC sputtering switch')
      ElseIf (Tag(1:3).eq.'H51') Then
        Call ReadI(line,hc_sticking_coef_model,0,2,
     >     'Sticking model, 0-preset, 1-Janev, 2-Alman&Ruzic')
      ElseIf (Tag(1:3).eq.'H52') Then
        Call ReadR(line,hc_sticking_coef_preset,-1.0,1.0,
     >     'Preset sticking coefficient, -1.0=CTRESH, 0.0-1.0 '
     >     // 'stuck')
      ElseIf (Tag(1:3).eq.'H53') Then
        Call ReadI(line,hc_sputtering_species_model,0,2,
     >     'Sputtered species model, 0-preset, 1-A&R, 2-same')
      ElseIf (Tag(1:3).eq.'H54') Then
        Call ReadI(line,hc_sputtering_energy_model,0,3,
     >     'Sputtering energy model,0-set,1-impact,2-thermal,3-AR')
      ElseIf (Tag(1:3).eq.'H55') Then
        Call ReadR(line,hc_sput_energy_neutral_preset,0.0,1E4,
     >     'Preset sputtering particle energy, neutral impact (eV)')
      ElseIf (Tag(1:3).eq.'H56') Then
        Call ReadR(line,hc_sput_energy_ion_preset,0.0,1E4,
     >     'Preset sputtering particle energy, ion impact (eV)')
      ElseIf (Tag(1:3).eq.'H57') Then
        Call ReadI(line,hc_sputtering_angle_model,-1,11,
     >     'Sput. angle model,-1=NRFOPT,0-off,1-spec,2-isosin,'
     >     // '3-isocos,4-sqrtcos,5-proj-sqrtsin,10-norm,11-A&R ')
c
c---------------------------
c
c jdemod
c     Wall segment index for HC_launch_location option 6
c 
      ElseIf (Tag(1:3).eq.'H60') Then
        Call ReadI(line,hc_launch6_wall_index,1,maxpts,
     >     'Wall segment index for HC launch option 6')
c
c       Set flag indicating that the value has been read
c
        hc_launch6_wall_index_set=.true.
c
c       Set flag indicating that the value has been read
c
        hc_launch6_wall_index_set=.true.
c
c---------------------------
c
c     HC reaction kinetics option:
c     0 = off
c     1 = on
c     2+= advanced options (to be implemented)
c 
      ElseIf (Tag(1:3).eq.'H61') Then
        Call ReadI(line,hc_reaction_kinetics,0,4,
     >     'HC reaction kinetics option')
c---------------------------
c jdemod
c
c     HC reaction kinetics option:
c     0 = original code
c     1 = new implementation of reaction kinetics
c 
      ElseIf (Tag(1:3).eq.'H62') Then
        Call ReadI(line,hc_kinetics_opt,0,1,
     >     'NEW HC reaction kinetics option')
c---------------------------
c jdemod
c
c     HC kinetic suboption - Vperp evolution:
c     0 = Vperp remains constant between reactions
c     1 = Vperp = Vpara assigned at reaction start
c     2 = Vperp = f(Tperp) - Vperp calculated from temperature
c     3 = Vperp diffuses independent of Vpara at same rate
c     
c 
      ElseIf (Tag(1:3).eq.'H63') Then
        Call ReadI(line,hc_vperp_opt,0,3,
     >     'HC Perpendicular velocity option')

c
c
c     jdemod - added an input parameter to define the H isotope
c              in the hydrocarbon molecules so it can be different
c              from the background plasma 
c
      elseif (Tag(1:3).eq.'H64') then
c slmod begin - *** TEMP ***
        Call ReadR(line,input_HC_H_mass,1.0,3.0,
     >     'Mass of the H isotope in HC molecules 1.0->3.0')
c        STOP 'TURNING OFF FOR NOW...'
c slmod end
c
c
c jdemod
c---------------------------
c
c        Values assigned in global_hc_assign_inputs after the 
c        entire input file has been read in
c
c        If (hc_sputtering_angle_model .eq. -1) Then
c		hc_sputtering_angle_model = NRFOPT
c	End If

! Hydrocarbon output options.

      ElseIf (Tag(1:3).eq.'H90') Then
        Call ReadI(line,hc_coord_print_option,0,1,
     >     'Print r,z position data at each timestep, 0-off, 1-on')
      ElseIf (Tag(1:3).eq.'H91') Then
        Call ReadI(line,hc_evolve_print_option,0,1,
     >     'Print r,z position data at each transition, 0-off, 1-on')


! ----------------------------
!
! End hydrocarbon following options.
!
! ----------------------------
! ammod end.

      ELSE
        CALL ER('ReadTagSeriesH','Unrecognized tag',*99)
      ENDIF
c
c
      RETURN
c
c
99    WRITE(SLOUT,'(5X,3A)') 'LINE = "',line,'"'
      WRITE(SLOUT,'(5X,3A)') 'TAG  = "',tag ,'"'
      WRITE(0    ,'(5X,3A)') 'LINE = "',line(1:LEN_TRIM(line)),'"'
      WRITE(0    ,'(5X,3A)') 'TAG  = "',tag ,'"'
      WRITE(0,*) '    DIVIMP HALTED'
      STOP
      END
c
c ======================================================================
c
c subroutine: ReadTagSeries_I
c
c Unstructure input related to the magnetic grid are loaded.  The input
c tag starts with 'I'.
c
c
      SUBROUTINE ReadTagSeries_I(line,tag,fp)
      use mod_params
      use mod_comtor
      use mod_slcom
      use mod_fperiph_com
      use mod_cedge2d
      use mod_promptdep
      use mod_io
      IMPLICIT none
c
c     READ "I" Series Unstructured input
c
      INTEGER   fp
      CHARACTER line*(*),tag*3

c     INCLUDE 'params'
c     INCLUDE 'comtor'
c     INCLUDE 'slcom'
c     include 'fperiph_com'
c
c -----------------------------------------------------------------------
c
c     TAG I24 : Ion initial positon in cell 
c
c     Option to turn on/off more exact calculation of particle initial 
c     positions. (default = 1 = on)
c
      IF     (tag(1:3).EQ.'I24') THEN
        CALL ReadI(line,init_pos_opt,0,1,
     >             'Refined init position options')
c
c -----------------------------------------------------------------------
c
c     TAG I25 : FP Neutral ionization option
c
      ELSEIF (tag(1:3).EQ.'I25') THEN
        CALL ReadI(line,fp_neut_opt,0,1,
     >             'FP NEUTRAL IONIZATION OPTION')
c
c -----------------------------------------------------------------------
c
c     TAG I26 : FP plasma option
c
      ELSEIF (tag(1:3).EQ.'I26') THEN
        CALL ReadI(line,fp_plasma_opt,0,4,
     >             'FP PLASMA OPTION')
c
c -----------------------------------------------------------------------
c
c     TAG I27 : FP temperature
c
      ELSEIF (tag(1:3).EQ.'I27') THEN
        CALL ReadR(line,fp_te,0.0,HI,
     >             'FP Temperature in eV')
c
c -----------------------------------------------------------------------
c
c     TAG I28: FP density
c
      ELSEIF (tag(1:3).EQ.'I28') THEN
        CALL ReadR(line,fp_ne,0.0,HI,
     >             'FP Density in m-3')
c
c -----------------------------------------------------------------------
c
c     TAG I31: FP flow option
c     Far periphery transport flow option
c     0 - no flow in far periphery
c     1 - flow in far periphery is the same as associated ring
c     2 - flow in far periphery is specified as input
c     3 - flow pattern in far periphery is based on associated virtual ring
c
      ELSEIF (tag(1:3).EQ.'I31') THEN
        CALL ReadI(line,fp_flow_opt,0,3,
     >             'FP FLOW OPTION')
c
c -----------------------------------------------------------------------
c
c     TAG I32: FP flow velocity for option 2

      ELSEIF (tag(1:3).EQ.'I32') THEN
        CALL ReadR(line,fp_flow_velocity_input,-HI,HI,
     >             'FP Flow velocity for FP flow option 2')
c
c -----------------------------------------------------------------------
c
c     TAG I33: FP Number of radial bins in FP grid
c
      ELSEIF (tag(1:3).EQ.'I33') THEN
        CALL ReadI(line,fp_n_bins,1,maxnrs,
     >             'Number of radial bins for FP grid')
c
c -----------------------------------------------------------------------
c
c     TAG I34: FP Grid width option
c     
c     Defines the option used to choose the width of simple crude
c     FP mesh
c     Option 0 = maximum distance from edge cell to wall from 
c               fp_walldist values
c     Option 1 = width of grid is specified using fpxmaxo for MAIN and 
c                fpxmaxi for PFZ
c
      ELSEIF (tag(1:3).EQ.'I34') THEN
        CALL ReadI(line,fp_grid_width_opt,0,1,
     >             'Number of radial bins for FP grid')
c
c -----------------------------------------------------------------------
c
c     TAG I35 : Main Chamber Ion reflection coefficient for fpopt 1
c
      ELSEIF (tag(1:3).EQ.'I35') THEN
        CALL ReadR(line,mc_recyc,0.0,1.0,
     >             'MC Reflection coefficient')
c -----------------------------------------------------------------------
c
c     TAG I36 : Private Fluz Zone Ion reflection coefficient for fpopt 1
c
      ELSEIF (tag(1:3).EQ.'I36') THEN
        CALL ReadR(line,pfz_recyc,0.0,1.0,
     >             'PFZ Reflection coefficient')
c
c
c
c -----------------------------------------------------------------------
c
c     TAG I37 : Fluid code charge state index to use for ion injection
c               options 12 and 13        
c
      ELSEIF (tag(1:3).EQ.'I37') THEN
        CALL ReadI(line,e2diz_inj,1,maxe2dizs,
     >           'Fluid code impurity charge state index for injection')
     
      ! Tag I38 for prompt deposition option 4, the average charge state 
      ! of the ion near the target. This value is used in the calculation
      ! of the gyroradius instead of the calculated value.
      elseif (tag(1:3).eq.'I38') then
        call readr(line, prompt_dep_avg_z, 0.0, hi, 
     >    'Average charge near target for prompt dep option 4')
c
c
c -----------------------------------------------------------------------
        
        

c...  REPLACE! (?) - jdemod - not sure what is up here - why does it need
c                             replacement        
c
c     TAG I29 and I30 - corner points for line/box injections 
c
      ELSEIF (tag.EQ.'I29') THEN
        CALL Read2R(line,cxscA,cyscA,-HI,HI,
     >              'line/box injection point A')

      ELSEIF (tag.EQ.'I30') THEN
        CALL Read2R(line,cxscB,cyscB,-HI,HI,
     >              'line/box injection point B')

      ELSE
        CALL ER('ReadTagSeriesI','Unrecognized tag',*99)
      ENDIF
c
c     jdemod - At least up to I36 - see above
c     

      
      RETURN
99    WRITE(SLOUT,'(5X,3A)') 'LINE = "',line,'"'
      WRITE(SLOUT,'(5X,3A)') 'TAG  = "',tag ,'"'
      WRITE(0    ,'(5X,3A)') 'LINE = "',line(1:LEN_TRIM(line)),'"'
      WRITE(0    ,'(5X,3A)') 'TAG  = "',tag ,'"'
      WRITE(0,*) '    DIVIMP HALTED'
      STOP
      END
c
c ======================================================================
c

      subroutine read_out_unstructured_input(line,tag,fp)
      use mod_params
      use mod_slcom
      use mod_out_unstruc
      use mod_io
      implicit none
c
c     READ "O" Series Unstructured input
c
c     The Oh is used for OUT related tagged input
c
      INTEGER   fp
      CHARACTER line*(*),tag*3

c     INCLUDE 'params'
c     include 'slcom'
c     include 'out_unstruc'
c
c
c      
c
c -----------------------------------------------------------------------
c
c     ADD TAGS RELATED TO OUT - USING SERIES 'O' oooh :) ... for OUT
c
c -----------------------------------------------------------------------
c

      if (tag(1:3).eq.'O01') then 
c
c     new_absfac - use to change scaling of plots from OUT
c
c     This option allows the absolute scaling factor for the DIVIMP
c     run results to be specified in the OUT routine. It's default
c     value is zero.
c

        CALL ReadR(line,new_absfac,0.0,HI,
     >                   'Imposed ABSFAC in OUT')

c        
c     Core fueling code calculates integrated ionization
c     profiles in the core ... these parameters allow the 
c     PSIN inner bound of the integration regions to be set
c     This is used in the pr_eirene_analysis routine
c     
      elseif (tag(1:3).eq.'O02') then 
c
c        O02 - PSIN bound for calculating core ionization 
c              profile 1 (psi1_reg)
c
        CALL ReadR(line,psi1_reg,0.0,HI,
     >                   'PSIN bound for core ionization profile 1')
c
      elseif (tag(1:3).eq.'O03') then 
c
c        O03 - PSIN bound for calculating core ionization 
c              profile 2 (psi2_reg)
c
        CALL ReadR(line,psi2_reg,0.0,HI,
     >                   'PSIN bound for core ionization profile 2')
c
c
      elseif (tag(1:3).eq.'O04') then 
c
c     absfac_opt - options specifying how to calculate scaling factor for DIVIMP
c                  results in OUT         
c
c     This option allows the absolute scaling factor for the DIVIMP
c     run results to be specified in the OUT routine. It's default
c     value is zero.
c

         CALL ReadI(line,absfac_opt,0,3,
     >              'ABSFAC calculation option in OUT')

      elseif (tag(1:3).eq.'O05') then 
c
c     e2dizs_offset - offset to match fluid code impurity charge state results data to 
c                     DIVIMP results (needed sometimes since the fluid code results
c                     sometimes contain multiple fluids). 
c
        CALL ReadI(line,e2dizs_offset,0,100,'FC Impurity offset index')

      elseif (tag(1:3).eq.'O06') then 
c
c     e2dizs_offset - offset to match fluid code impurity charge state results data to 
c                     DIVIMP results (needed sometimes since the fluid code results
c                     sometimes contain multiple fluids). 
c
        CALL ReadI(line,absfac_iz,0,maxizs+1,'FC ABSFAC IZ (no offset)')

      elseif (tag(1:3).eq.'O07') then 
c
c     e2dizs_offset - offset to match fluid code impurity charge state results data to 
c                     DIVIMP results (needed sometimes since the fluid code results
c                     sometimes contain multiple fluids). 
c
        CALL ReadI(line,absfac_ir,1,maxnrs,'FC ABSFAC IR')

      elseif (tag(1:3).eq.'O08') then 
c
c     e2dizs_offset - offset to match fluid code impurity charge state results data to 
c                     DIVIMP results (needed sometimes since the fluid code results
c                     sometimes contain multiple fluids). 
c
        CALL ReadI(line,absfac_ikstart,1,maxnks,'FC ABSFAC IK-START')

      elseif (tag(1:3).eq.'O09') then 
c
c     e2dizs_offset - offset to match fluid code impurity charge state results data to 
c                     DIVIMP results (needed sometimes since the fluid code results
c                     sometimes contain multiple fluids). 
c
        CALL ReadI(line,absfac_ikend,1,maxnks,'FC ABSFAC IK-END')

      elseif (tag(1:3).eq.'O10') then 
c
c     scale_1d: Scale generalized results plots from load_divdata by using the
c               cell length along the field line instead of the cell volume.
c               This applies to code density results and others normalized
c               by cell area.          
c
        CALL ReadI(line,scale_1d,0,2,'SCALE 1D')
c
      ELSE
          CALL ER('ReadUnstructuredInput','Unrecognized tag',*99)
      ENDIF

      return

 99   WRITE(SLOUT,'(5X,3A)') 'LINE = "',line,'"'
      WRITE(SLOUT,'(5X,3A)') 'TAG  = "',tag ,'"'
      WRITE(0    ,'(5X,3A)') 'LINE = "',line(1:LEN_TRIM(line)),'"'
      WRITE(0    ,'(5X,3A)') 'TAG  = "',tag ,'"'
      WRITE(0,*) '    DIVIMP HALTED'
      STOP
      END
c
c ======================================================================
c
