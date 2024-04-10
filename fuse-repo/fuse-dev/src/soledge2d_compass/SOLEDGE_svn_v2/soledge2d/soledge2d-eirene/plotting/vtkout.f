      subroutine eirene_vtkout

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_cadgeo
      use eirmod_clgin
      use eirmod_cplot
      use eirmod_ctext

      implicit none

      integer :: ilim,j
      integer :: rgbcol(3,7) = reshape( 
     .                         (/ 0, 0, 0,         ! black
     .                            1, 0, 0,         ! red
     .                            0, 0, 1,         ! blue
     .                            0, 1, 0,         ! green
     .                            1, 1, 0,         ! yellow
     .                            0, 1, 1,         ! cyan
     .                            1, 0, 1 /),      ! magenta
     .                         (/ 3, 7 /))
      
      open (unit=28,file='vtk.out')
      write (28,'(a)') 
     .  '<?xml version="1.0" encoding="ISO-8859-1" ?>'

      write (28,'(a)')  
      write (28,'(a)') '<input>'

      write (28,'(a,ss,6es12.4,a)') 
     .                 '<bound>', ch3x0-ch3mx, ch3x0+ch3mx,
     .                            ch3y0-ch3my, ch3y0+ch3my,
     .                            ch3z0-ch3mz, ch3z0+ch3mz,
     .                 ' </bound>'

! camera
      write (28,'(a)') '<camera>'
      write (28,'(a)') '<background color="1,0,1" /> '
      write (28,'(a,ss,3(es12.4,a),a)')
     .                 '<focalpoint coordinate="',
     .                  ch3x0,',',ch3y0,',',ch3z0, '"/>'
      write (28,'(a,ss,3(es12.4,a),a)')
     .                 '<position coordinate="',
     .                  ch3x0+2*ch3mx,',',ch3y0,',',ch3z0+2*ch3mz,
     .                 '"/>'
      write (28,'(a)') '<viewup coordinate="0,0,0"/>'
      write (28,'(a)') '</camera>'


! geometry
      write (28,'(a)') '<geometry resolution="10">'

      do ilim=1,nlimi

        IF (IGJUM0(ilim).NE.0) cycle
        
        write (28,'(3a)') '<!-- ',txtsfl(ilim),' -->'
        if (rlb(ilim) < 7._dp) then

           write (28,'(a,3(i6,a),a)')
     .                 '<quadric color="',
     .                 rgbcol(1,ilcol(ilim)),',',
     .                 rgbcol(2,ilcol(ilim)),',',
     .                 rgbcol(3,ilcol(ilim)),
     .                 ' " >'
           write (28,'(a,ss,10es12.4,a)') 
     .                 '<coefficients>',
     .                 a0lm(ilim),a1lm(ilim),a2lm(ilim),a3lm(ilim),
     .                            a4lm(ilim),a5lm(ilim),a6lm(ilim),
     .                            a7lm(ilim),a8lm(ilim),a9lm(ilim),
     .                 ' </coefficients>'
           if (rlb(ilim) < 0.) then

             do j=1,ilin(ilim)
               write (28,'(a,ss,4es12.4,a)') 
     .                 '<limquadric>',
     .                   alims(j,ilim),xlims(j,ilim),
     .                   ylims(j,ilim),zlims(j,ilim),
     .                 ' </limquadric>'
             end do

             do j=1,iscn(ilim)
               write (28,'(a,ss,10es12.4,a)') 
     .                 '<limquadric>',
     .                   alims0(j,ilim),xlims1(j,ilim),ylims1(j,ilim),
     .                   zlims1(j,ilim),xlims2(j,ilim),ylims2(j,ilim),
     .                   zlims2(j,ilim),xlims3(j,ilim),ylims3(j,ilim),
     .                   zlims3(j,ilim),
     .                 ' </limquadric>'
             end do
              
           else if ((rlb(ilim) > 0.) .and. (rlb(ilim) < 3.)) then

              write (28,'(a,ss,6es12.4,a)') 
     .                 '<bound>', xlims1(1,ilim),xlims2(1,ilim),
     .                            ylims1(1,ilim),ylims2(1,ilim),
     .                            zlims1(1,ilim),zlims2(1,ilim),
     .                 ' </bound>'
              
           else if ((rlb(ilim) > 2.) .and. (rlb(ilim) < 4.)) then
              write (6,*) ' 2 < rlb < 4, ilim = ',ilim
           else if ((rlb(ilim) > 4.) .and. (rlb(ilim) < 5.)) then
              write (6,*) ' 4 < rlb < 5, ilim = ',ilim
           else
              write (6,*) ' 5 < rlb < 5, ilim = ',ilim
           end if
           
           write (28,'(a)') '</quadric>'

        else
        end if

      end do
      
      write (28,'(a)') '</geometry>'

      write (28,'(a)') '</input>'

      return

      end subroutine eirene_vtkout


      
