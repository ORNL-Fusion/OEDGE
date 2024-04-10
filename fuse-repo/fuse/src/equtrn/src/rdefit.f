      subroutine rdefit(lun,ngpr,ngpz,iret, title,date,ipestg,nr,nz,
     ,           rdim,zdim,zmsmid,rcntc,redge,rma,zma,psimin,psilim,
     ,           btorc,fg,pg,ffg,ppg,pfm,rgr,zgr)
c=====================================================
c*** where:
c***
c*** i)  nr, nz, rdim, redge, and zdim define the rectangular mesh used
c***     to store the psi values via the functions rr & zz defined at
c***     the top of these code,
c***
c*** ii) psimin is the flux value at the magnetic axis; psilim is the
c***     value at the separatrix,
c***
c*** iii) btorc is the toroidal magnetic field at a radius rcntc,
c***
c*** iv) the poloidal flux is pfm,
c***
c*** v) fg is the flux function R*Btor; ffg is its derivative with
c***    respect to psi,
c***
c*** vi) pg is the pressure & ppg is its derivative.
c=====================================================
c
c  version : 18.12.94 18:33
c
      real*8 fg(*),pg(*),ffg(*),ppg(*),pfm(ngpr,*),rgr(*),zgr(*)
      real*8 rdim,zdim,rcntc,redge,zmsmid,rma,zma,psimin,psilim,btorc
      character title*40, date*8
c=====================================================
      rr(r)=r/(nr-1)*rdim+redge
      zz(z)=(z-(nz+1)/2)/(nz-1)*zdim+zmsmid
c=====================================================
c
      iret=0
      rewind lun
      read(lun,'(a40,a8,3i4)') title,date,ipestg,nr,nz
      if(nr.gt.ngpr) then
          print *,'=== rdefit: nr > ngpr'
          iret=2
      end if
      if(nz.gt.ngpz) then
          print *,'=== rdefit: nz > ngpz'
          iret=2
      end if
      if(nr.le.0) then
          print *,'=== rdefit: nr < 1'
          iret=4
      end if
      if(nz.le.0) then
          print *,'=== rdefit: nz < 1'
          iret=4
      end if
      if(iret.ne.0) return
c
      read(lun,*) rdim,zdim,rcntc,redge,zmsmid
      read(lun,*) rma,zma,psimin,psilim,btorc
      read(lun,*)
      read(lun,*)
      read(lun,*) (fg(i),i=1,nr)
      read(lun,*) (pg(i),i=1,nr)
      read(lun,*) (ffg(i),i=1,nr)
      read(lun,*) (ppg(i),i=1,nr)
      read(lun,*) ((pfm(i,j),i=1,nr),j=1,nz)
      do i=1,nr
        rgr(i)=rr(float(i-1))
      enddo
      do i=1,nz
        zgr(i)=zz(float(i))
      enddo
c
      end
