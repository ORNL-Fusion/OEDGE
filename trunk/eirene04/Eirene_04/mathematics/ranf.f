      function ranf_eirene ()
!pb      REAL FUNCTION GSU2R(ISEED)
C
C RUNDOM NUMBER GENERATOR FROM  http://www.srcc.msu.su/num_anal/lib_na/cat/g/gsu1r.htm
C SOURCE [D. Knuth, The Art of  Computer Programmig, vol. 2] 
C ISEED IS THE INTEGER FROM 1 TO  2147483646, AFTER FINISHING ITS VALUE IS
C (2**31) * R (N) AN DCAN BE USED FOR THE FUTURE RUNS
C RETURNS ONE RUNDOM NUMBER FROM 0 TO 1
C
      USE PRECISION
      implicit none
      integer :: iseed
      common /cmem/ iseed
      integer, save :: ifirst=0
      integer :: ise
      real(dp) :: ra, dummy, ranf_eirene, ranset_eirene

      INTEGER D2P32M
      DOUBLE PRECISION Z,D2P31M,D2PN31,DMOD,DFLOAT
      DATA  D2PN31/4.656612873077393D-10/,D2P31M/
     12147483647.D0/,D2P32M/16807/

      if (ifirst == 0) then
         ise = -1
         dummy = ranset_eirene(ise)
         ifirst = 1
      end if

      Z=DFLOAT(ISEED)
      Z=DMOD(D2P32M*Z,D2P31M)
      RA=Z*D2PN31
      ISEED=Z

      ranf_eirene=ra
      return
      end
