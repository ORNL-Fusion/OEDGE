C
C
      SUBROUTINE SCATANG (INENERGY,INRAN,ELTHETA,CTTHETA,SIGM)
c     ***********************************************************
c     *                                                         *
c     *Programm zur Berechnung der Streuwinkel bei vorgegebener *
c     *Energie und Zufallszahl (Monte-Carlo-Simulation) für eine*
c     *atomare Reaktion.                                        *
c     *Die verschiedenen Reaktionen werden mittels der Daten von*
c     *P.S. Krstic und D.R. Schultz (Atomic and Plasma-Material *
c     *Interaction Data for Fusion), Differntielle Querschnitte *
c     *mit zugehörigen Winkeln, berechnet.                      *
c     *Das Programm berechnet Reaktionen zwischen Energien von  *
c     *0.1 eV und 100 eV (CM). Dabei sind 31 Energiewerte aus   *
c     *den o.g.Daten vorgegeben.                                *
c     *Input:    Dateien der jeweiligen Reaktion                *
c     *          Reaktionsenergy: inenergy                      *
c     *          Zufallszahl: inran                             *
c     *Output:   Verwendete Energie: energy                     *
c     *          Streuwinkel: eltheta (elastisch)               *
c     *                       cttheta (ladungsaustausch)        *
c     *Totale Streuquerschnitte + Momente für alle Energien:    *
c     *  sigma(i), gamma(i), moment(i), viscos(i)               *
c     *mit dem Feldindex i der Energie (zugeh. Vorgabe i=m).    *
c     *                                                         *
c     *Mögliche Energien: Feld energy(i)                        *
c     *Winkel der diff.Querschnitte abh.von Energie: theta(k,i) *
c     *mit k=1..768                                             *
c     *                                                         *
c     *(19.9.2000)                          Torsten Haberscheidt*
c     *                                                         *
c     *                                                         *
c     *modifications: 24.3.03:                                  *
c     *   remove "ct"-scattering angle evaluation               *
c     *   remove "el"-scattering angle evaluation for inran.lt.0*
c     *           in this case: only SIGM =SIGM(INENERGY)       *
c     *           is returned                                   *
c     ***********************************************************

c     --Parameter--
      

      USE PRECISION
      IMPLICIT NONE
      integer i,j,k,l,m,n,dat,IFIRST,IUN
      REAL(DP) :: pi,z,w
      REAL(DP) :: energy(31),inenergy
      REAL(DP) :: theta(768,31),el(770),ct(770),dtheta(770)
      REAL(DP) :: sig,sigma(31),gam,gamma(31),SIGM
      REAL(DP) :: mom,moment(31),vis,viscos(31)
      REAL(DP) :: normel(770), normct(770)
      REAL(DP) :: elangle(768,31), ctangle(768,31)
      REAL(DP) :: inran, eltheta, cttheta
      CHARACTER(10) :: FILNAM
      DATA IFIRST /0/
      SAVE
      
c     ***********************************************************
c     -- Einlesen der Daten --

      IF (IFIRST == 0) THEN
        IFIRST = 1

        PI=4.D0*ATAN(1.D0)
        iun=23

        do 10 i=1,31
          IF (I <= 10) THEN
            WRITE (FILNAM,'(a3,i1,a4)') 'el-',i-1,'.dat'
          ELSE
            WRITE (FILNAM,'(a3,i2,a4)') 'el-',i-1,'.dat'
          END IF
        
          open (iun,file=filnam)

          read(iUN,*)
          read(iUN,100) energy(i)
100       format(T19,E9.4)

          do 20 j=1,13
            read(iUN,*)
20        continue
      
          dat=0

          do 30 k=1,768
            read(iun,200,end=5) theta(k,i), el(k), ct(k)
            dat= dat+1
30        continue

200       format(T5, E11.6, T20, E11.6, T35, E11.6)
5         continue
          close(iUN)

c       *******************************************************
c       --Verarbeitung der Daten--

c       Winkelintervalle
      
          dtheta(1) = 0.5*(theta(1,i) + theta(2,i))

          do 40 k=2,dat-1
            dtheta(k) = 0.5*(theta(k+1,i) - theta(k-1,i))
40        continue

          dtheta(dat) = pi - 0.5*(theta(dat-1,i) + theta(dat,i))

c       Streuquerschnitte

          sig=0
          gam=0
          mom=0
          vis=0

          do 50 k=1,dat
            sig = sig + (dtheta(k) * el(k))
            gam = gam + (dtheta(k) * ct(k))
            mom = mom + (dtheta(k) * el(k))*(1-cos(theta(k,i)))
            vis = vis + (dtheta(k) * el(k))*(sin(theta(k,i)))**2
50        continue

          sigma(i) = sig
          gamma(i) = gam
          moment(i)= mom
          viscos(i)= vis

c       Normierung

          do 60 k=1,dat
            normel(k) = (dtheta(k) * el(k)) / sig
            normct(k) = (dtheta(k) * ct(k)) / gam
60        continue

c       Generieren von Zahlen [0;1]

          elangle(1,i) = normel(1)
          ctangle(1,i) = normct(1)

          do 70 k=2,dat
            elangle(k,i) = elangle(k-1,i) + normel(k)
            ctangle(k,i) = ctangle(k-1,i) + normct(k)
70        continue



10      continue

      END IF   ! END OF IFIRST-BLOCK

c     ***********************************************************

c     Berechnung des Monte-Carlo-Winkels mit linearer Interpolation

c  binary search in energy array
      call energyloc(energy,31,inenergy,m)

c     elastic (both with and without charge transfer, IP-model)

      if (inran.ge.0._DP) then

c  binary search in angle array, given the energy energy(m)
      if(inran.lt.elangle(1,m))then
        eltheta = inran*theta(1,m)/elangle(1,m)
      else if(inran.eq.elangle(dat,m))then
        eltheta = pi
      else
        call thetaloc(elangle(1:dat,m),dat,inran,n)
        z=(inran - elangle(n,m))/(elangle(n+1,m) - elangle(n,m))+n
      
        if(z.gt.(dat-1))then
          eltheta = (z-n)*(pi - theta(n,m)) + theta(n,m)
        else
          eltheta = (z-n)*(theta(n+1,m) - theta(n,m))+theta(n,m)
        endif
      endif

      endif

c chargetransfer
c currently not used

c     if(inran.lt.ctangle(1,m))then
c       cttheta = inran*theta(1,m)/ctangle(1,m)
c     else if(inran.eq.ctangle(dat,m))then
c       cttheta = pi
c     else
c       call thetaloc(ctangle(1:dat,m),dat,inran,l)
c       w=(inran - ctangle(l,m))/(ctangle(l+1,m) - ctangle(l,m))+l
c
c       if(w.gt.(dat-1))then
c         cttheta = (w-l)*(pi - theta(l,m)) + theta(l,m)
c       else
c         cttheta = (w-l)*(theta(l+1,m) - theta(l,m))+theta(l,m)
c       endif
c     endif

      SIGM=SIGMA(M)


      RETURN
      end
