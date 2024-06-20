      subroutine wrefit(lun,ngpr,iret, title,date,ipestg,nr,nz,
     ,           rdim,zdim,zmsmid,rcntc,redge,rma,zma,psimin,psilim,
     ,           btorc,fg,pg,ffg,ppg,pfm)
c=====================================================
c*** where:
c***
c*** i) nr, nz, rdim, redge, zdim, and zmsmid define the rectangular
c***    mesh used to store the psi values via the functions rr & zz
c***    defined at the top of these code,
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
c  version : 18.12.94 18:29
c
      real*8 fg(*),pg(*),ffg(*),ppg(*),pfm(ngpr,*)
      real*8 rdim,zdim,rcntc,redge,zmsmid,rma,zma,psimin,psilim,btorc
      character title*40, date*8
c=====================================================
c
      iret=0
      if(nr.le.0) then
          print *,'=== wrefit: nr < 1'
          iret=4
      end if
      if(nz.le.0) then
          print *,'=== wrefit: nz < 1'
          iret=4
      end if
c
      rewind lun
      write(lun,'(a40,a8,3i4)') title,date,ipestg,nr,nz
      write(lun,'(5e16.9)') rdim,zdim,rcntc,redge,zmsmid
      write(lun,'(5e16.9)') rma,zma,psimin,psilim,btorc
      write(lun,'(5e16.9)')
      write(lun,'(5e16.9)')
      write(lun,'(5e16.9)') (fg(i),i=1,nr)
      write(lun,'(5e16.9)') (pg(i),i=1,nr)
      write(lun,'(5e16.9)') (ffg(i),i=1,nr)
      write(lun,'(5e16.9)') (ppg(i),i=1,nr)
      write(lun,'(5e16.9)') ((pfm(i,j),i=1,nr),j=1,nz)
c
      end
