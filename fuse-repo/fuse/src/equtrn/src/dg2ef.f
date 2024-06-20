      program dg2ef
c
c  version : 08.01.2002 16:29
c
c=======================================================================
c*** Translation of dg-compatible equilibrium data into efit format
c***
c*** The input and output files must be pre-connected to the units
c*** fort.1 and fort.2, and the field data to the fort.3
c=======================================================================
      implicit none
      integer ngpr,ngpz
      parameter (ngpr=257, ngpz=513)
      real*8 fg(ngpr),pg(ngpr),ffg(ngpr),ppg(ngpr)
      real*8 pfm(ngpr,ngpz),rgr(ngpr),zgr(ngpz)
      real*8 rdim,zdim,rcntc,redge,zmsmid,rma,zma,psimin,psilim,btorc
      integer i,j,ia,ja,nr,nz,iia,iret,ipestg
      real*8 u,v
      character title*40, date*9
      data title/'Conversion from the dg format'/
c=======================================================================
c
c      call date2(date)
      call date_and_time(date)
      call rdeqdg(1,ngpr,ngpz,iret, nr,nz,btorc,rcntc,rgr,zgr,pfm)
c      print *,'After rdeqdg: nr,nz=',nr,nz  !###
      if(iret.ne.0) then
          print *,'==== dg2ef: error in rdeqdg. iret =',iret
          stop
      end if
c
      do i=1,nr !{
        fg(i)=btorc*rcntc
        pg(i)=1.
        ffg(i)=0.
        ppg(i)=0.
      end do !}
      ipestg=3
      psilim=0.
      psimin=1.e30
      do j=nz,1,-1 !{
c        u=vmin(pfm(1,j),nr)
        u=1.e30
        iia=-1
        do i=1,nr !{
          if(pfm(i,j).lt.u) then !{
            u=pfm(i,j)
            iia=i
          end if !}
        end do !}
        if(u.lt.psimin) then !{
          psimin=u
          ja=j
c          ia=lvdmin(pfm(1,j),nr)
          ia=iia
        end if !}
      end do !}
c      print *,'ia,ja,psimin=',ia,ja,psimin  !###
      redge=rgr(1)
      rdim=rgr(nr)-rgr(1)
      zdim=abs(zgr(nz)-zgr(1))
      zmsmid=0.5*(zgr(nz)+zgr(1))
      rma=float(ia)/(nr-1)*rdim+redge
      zma=(float(ja)-(nz+1)/2)/(nz-1)*zdim

      call wrefit(2,ngpr,iret, title,date,ipestg,nr,nz,
     ,           rdim,zdim,zmsmid,rcntc,redge,rma,zma,psimin,psilim,
     ,           btorc,fg,pg,ffg,ppg,pfm)
      if(iret.ne.0) then !{
        print *,'==== dg2ef: error in wrefit. iret = ',iret
      end if !}
c=======================================================================
      end
