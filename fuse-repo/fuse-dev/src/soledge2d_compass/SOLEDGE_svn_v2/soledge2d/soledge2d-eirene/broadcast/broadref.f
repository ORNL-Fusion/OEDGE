 
 
      SUBROUTINE EIRENE_BROADREF
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
!pb      USE EIRMOD_CREFMOD
      USE EIRMOD_CREF
      USE EIRMOD_CPES
      IMPLICIT NONE
 
      INCLUDE 'mpif.h'
      INTEGER :: IER
 
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
 
      CALL MPI_BCAST (RINTEG,NLIMPS+1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EINTEG,NLIMPS+1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (AINTEG,NLIMPS+1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RCREF,NCREF,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ICREF,MCREF,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LTRMOL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
 
      IF (MY_PE /= 0) CALL EIRENE_ALLOC_CREF
      CALL MPI_BCAST (HFTR0,NHD1*NHD2*NHD6,
     .                MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (HFTR1,NHD1*NHD2*NHD3*NHD6,
     .                MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (HFTR2,NHD1*NHD2*NHD3*NHD4*NHD6,
     .                MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (HFTR3,NHD1*NHD2*NHD3*NHD4*NHD5*NHD6,
     .                MPI_REAL8,0,MPI_COMM_WORLD,ier)
 
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
 
      RETURN
      END
