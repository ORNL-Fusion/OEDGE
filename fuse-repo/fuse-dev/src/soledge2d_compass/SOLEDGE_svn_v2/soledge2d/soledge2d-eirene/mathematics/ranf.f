      function ranf_eirene ()
 
C
C RANDOM NUMBER GENERATOR FROM
C  http://www.srcc.msu.su/num_anal/lib_na/cat/g/gsu1r.htm (in russian)
C
C SOURCE:  Knuth, D.E. 1981, Seminumerical Algorithms, 2nd ed., vol. 2 of The Art
C          of Computer Programming (Reading, MA: Addison-Wesley)
C
C ISEED IS THE INTEGER FROM 1 TO  2147483646, AFTER FINISHING ITS VALUE IS
C (2**31) * R (N) AND CAN BE USED FOR THE FUTURE CALLS
C RETURNS ONE RANDOM NUMBER FROM 0 TO 1
C
      USE EIRMOD_PRECISION
      implicit none
      integer :: iseed
      common /cmem/ iseed
      integer, save :: ifirst=0
      integer :: ise
      real(dp) :: ra, dummy, ranf_eirene, ranset_eirene,
     &            ranf_eirene_reinit
 
      INTEGER D2P32M
      DOUBLE PRECISION Z,D2P31M,D2PN31,DMOD,DFLOAT
      DATA  D2PN31/4.656612873077393D-10/,D2P31M/
     .             2147483647.D0/,D2P32M/16807/
 
      if (ifirst == 0) then
         ise = -1
         dummy = ranset_eirene(ise)
         ifirst = 1
      end if
 
!pb      Z=DFLOAT(ISEED)
      Z=REAL(ISEED,KIND=DP)
      Z=DMOD(D2P32M*Z,D2P31M)
      RA=Z*D2PN31
      ISEED=Z
 
      ranf_eirene=ra
      return
 
C     The following ENTRY is for reinitialization of EIRENE
 
      ENTRY ranf_eirene_reinit ()
      ifirst = 0
      ranf_eirene_reinit = 0.D0
      return
      end
