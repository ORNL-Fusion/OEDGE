c
c ======================================================================
c
c subroutine: Plot985
c
c 3D analysis
c
      SUBROUTINE Plot989(job,graph,ref,title,iopt,
     .                   xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .                   ismoth,ignors,itec,avs,navs)
      USE MOD_OUT989
      IMPLICIT none
      INCLUDE 'params'

c...  Input:
      INTEGER IOPT,ISMOTH,IGNORS(MAXNGS),ITEC,NAVS
      REAL    XXMIN,XXMAX,YYMIN,YYMAX,FP,FT,zadj,AVS(0:100)
      CHARACTER TITLE*(*),JOB*72,GRAPH*80
      CHARACTER*36 REF

c...  
      CALL Main989(iopt)

      RETURN
 99   STOP
      END

