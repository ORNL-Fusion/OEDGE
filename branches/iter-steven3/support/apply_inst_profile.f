      program apply_inst_profile
      implicit none
      ! this code reads in a block of line profile data from a divimp
      ! case and applies an alternate line profile instrument function
      ! and writes the data out again. 



      integer max_lp_bins
      parameter (max_lp_bins=51) 

      !
      ! Load base wavelength
      !
      real*8,parameter :: wave = 9094.83

      real*8 
     >   line_profile(-max_lp_bins:max_lp_bins,4),
     >   mod_line_profile(-max_lp_bins:max_lp_bins,4)
      real*8 lambda(-max_lp_bins:max_lp_bins)

      real*8 lam,lam0,dlam
      real*8 max_const_inst,max_isler_inst
      real*8 max_isler_inst_95,max_isler_inst_105
      real*8 frac1,frac2,parm1,parm2,base_parm1,base_parm2
      
      real*8 const_inst_funct,isler_inst_funct
      external const_inst_funct,isler_inst_funct
      real*8 isler_inst_funct2
      external isler_inst_funct2

      integer in,cnt,ii
      !

      line_profile = 0

      ! Load the LP data
      !
      ! Discard 4 lines
      !
      read(5,*)
      read(5,*)
      read(5,*)
      read(5,*)
      

      do in =  -max_lp_bins, max_lp_bins

           read(5,'(i5,1x,f10.5,4(1x,g20.12))')
     >             cnt,lambda(in),line_profile(in,1),
     >             line_profile(in,2),
     >             line_profile(in,3),
     >             line_profile(in,4)

           if (cnt.ne.in) write(0,*) 'Index mismatch:',in,cnt


      enddo


      !
      ! Apply 2 instrument functions - constant 0.12 and Isler:
      ! 0.95 exp ( -(dlam/(1.201*0.1))**2) + 0.05 * exp (-(dlam/(1.201*0.35))**2)
      !

      frac1 = 0.95
      frac2 = 0.05
      base_parm1 = 1.201 * 0.1
      base_parm2 = 1.201 * 0.35

      mod_line_profile = 0

      do in = -max_lp_bins+1, max_lp_bins-1
c
         do ii =  -max_lp_bins+1, max_lp_bins-1
c
            lam0 = lambda(in) ! in * lp_bin_width
            lam  = lambda(ii) ! ii * lp_bin_width
c 
            dlam = lam-lam0

            mod_line_profile(in,1)=mod_line_profile(in,1) +
     >         line_profile(ii,1) * 
     >         const_inst_funct(dlam)
c
            parm1 = 1.201 * 0.1
            parm2 = 1.201 * 0.35
            mod_line_profile(in,2)=mod_line_profile(in,2) +
     >         line_profile(ii,1) * 
     >         isler_inst_funct2(dlam,frac1,parm1,frac2,parm2)

            parm1 = 1.201 * 0.12
            parm2 = 1.201 * 0.35
            mod_line_profile(in,3)=mod_line_profile(in,3) +
     >         line_profile(ii,1) * 
     >         isler_inst_funct2(dlam,frac1,parm1,frac2,parm2)

            parm1 = 1.201 * 0.095
            parm2 = 1.201 * 0.225
            mod_line_profile(in,4)=mod_line_profile(in,4) +
     >         line_profile(ii,1) * 
     >         isler_inst_funct2(dlam,frac1,parm1,frac2,parm2)

         end do 
c
         !maxraw = max(maxraw,line_profile(in,1))
         max_const_inst = max(max_const_inst,mod_line_profile(in,1))
         max_isler_inst = max(max_isler_inst,mod_line_profile(in,2))
         max_isler_inst_95 = max(max_isler_inst_95,
     >                           mod_line_profile(in,3))
         max_isler_inst_105 = max(max_isler_inst_105,
     >                           mod_line_profile(in,4))
c
      end do

c      write(0,'(a,4(1x,g20.8)') 'MAXES:',max_const_inst,max_isler_inst,
c     >        max_isler_inst_95,max_isler_inst_105


c
c     Print out tabulated results
c

      
        write(6,'(a)')  
        write(6,'(a)')  'LINE PROFILE CALCULATION:'//
     >                  ' VARIOUS INSTRUMENT FUNCTIONS' 
        write(6,'(a)')  
        write(6,'(a,f12.5)')  'WAVELENGTH:',wave
        write(6,'(a)')  
        write(6,'(a5,1x,a10,1x,a12,4(1x,a20))') 'CNT','DLAM','LAM',
     >               ' CONST 0.12A ',' ISLER 0.1,0.35',
     >               ' ISLER 0.12,0.35',
     >               ' ISLER 0.095,0.225'
        do in =  -max_lp_bins, max_lp_bins
           write(6,'(i5,1x,f10.5,1x,f12.5,4(1x,g20.12))')
     >             in,lambda(in),wave+lambda(in),
     >             mod_line_profile(in,1)/max_const_inst,
     >             mod_line_profile(in,2)/max_isler_inst,
     >             mod_line_profile(in,3)/max_isler_inst_95,
     >             mod_line_profile(in,4)/max_isler_inst_105


        end do



      end



      real*8 function const_inst_funct(dlam)
      implicit none
      real*8 dlam
      real*8,parameter :: lp_instrument_width = 0.12

      const_inst_funct = exp(-(dlam/lp_instrument_width)**2)

      return
      end

      real*8 function isler_inst_funct(dlam)
      implicit none
      real*8 dlam,lp_instrument_width
      
      isler_inst_funct = 0.95 * exp(-(dlam/(1.201*0.1))**2) +
     >                   0.05 * exp(-(dlam/(1.201*0.35))**2)

      return
      end


      real*8 function isler_inst_funct2(dlam,frac1,parm1,frac2,parm2)
      implicit none
      real*8 dlam,lp_instrument_width
      real*8 frac1,frac2,parm1,parm2
      isler_inst_funct2 = frac1 * exp(-(dlam/parm1)**2) +
     >                   frac2 * exp(-(dlam/parm2)**2)

      return
      end
