      program proc_line_profiles
      implicit none

      integer, parameter :: maxcases = 10
      integer, parameter :: maxpoints = 51
      integer, parameter :: unum = 6
      integer, parameter :: onum = 7

      character*10 :: cases(maxcases)
      character*10 :: label(maxcases)
      character*20 :: base = 'd-123417-'
      character*10  :: ext = '.lim'
      character*100 :: outname = 'proc_line_profiles.'
      character*100 :: groupname
      
      integer :: extlen,baselen,caselen,searchlen,filelen
      integer :: outlen,iflen,grouplen
      integer :: ncases
      character*100 :: fname
      character*200 :: line
      character*100 :: search
      character*100 :: ifstring

      real*8 :: lp_dlam(-maxpoints:maxpoints)
      real*8 :: lp_data(maxcases,-maxpoints:maxpoints)
      real*8 :: lp_rawdata(maxcases,-maxpoints:maxpoints)
      real*8 :: baselam
      real*8 :: lp,lp_mod,lpn,lpn_mod
      real*8 :: d_lam

      integer :: in,ierr,id,ip ! loop variable




      !
      ! if_opt = instrument function option
      !      0 = 0.12A constant
      !      1 = Isler
      !      2 = 0.95 Isler
      !      3 = 1.05 Isler
      !
      integer :: if_opt 


      !
      ! Load base wavelength
      !
      real*8,parameter :: wave = 9094.83

      real*8 lam,lam0,dlam
      real*8 max_inst
      real*8 frac1,frac2,parm1,parm2,base_parm1,base_parm2
      
      real*8 const_inst_funct,isler_inst_funct
      external const_inst_funct,isler_inst_funct
      real*8 isler_inst_funct2
      external isler_inst_funct2

      integer cnt,ii
      !

      ! Isler IF = 1
      if_opt = 1
      groupname = 'E-dep-h'
      grouplen = len_trim(groupname)

      
      !
      ! ASSIGN case names and labels
      !

      ncases = 5
      cases(1) = 'f73'
      cases(2) = 'f74'
      cases(3) = 'f75'
      cases(4) = 'f76'

      cases(5) = 'f77'

      if (if_opt.eq.0) then 
         ifstring = 'INSTRUMENT FUNCTION: CONST 0.12A'
      elseif (if_opt.eq.1) then 
         ifstring = 'INSTRUMENT FUNCTION: ISLER'
      elseif (if_opt.eq.2) then 
         ifstring = 'INSTRUMENT FUNCTION: 95% * ISLER'
      elseif (if_opt.eq.3) then 
         ifstring = 'INSTRUMENT FUNCTION: 105% * ISLER'
      endif

      iflen = len_trim(ifstring)

      !
      ! Initilialize
      !
      extlen = len_trim(ext)
      baselen = len_trim(base)

      outname = outname(1:len_trim(outname))//groupname(1:grouplen)
      outlen = len_trim(outname)


      baselam=9094.83

      search='LINE PROFILE CALCULATION:'
      searchlen = len_trim(search)



      do in = 1,ncases

         caselen=len_trim(cases(in))

         fname = base(1:baselen)//cases(in)(1:caselen)//ext(1:extlen)
         filelen = len_trim(fname)
      
         open(unum,form='formatted',status='old',file=fname(1:filelen),
     >        iostat=ierr)
         if (ierr.ne.0) then 
            write(0,*) 'ERROR OPENING:',fname(1:filelen),
     >                 ' ERROR # =',ierr
            stop
         endif


         !
         ! Scan for line profile data
         !

 10      read(unum,'(a200)',end=100,err=100) line
         
         if (line(1:searchlen).eq.search(1:searchlen)) goto 20

         goto 10
         
 20      continue
         ! read raw line profile data
         read(unum,'(a200)') line
         read(unum,'(a200)') line
         read(unum,'(a200)') line

         !
         ! First data line
         !
         do id = -maxpoints,maxpoints

           read(unum,'(i5,1x,f10.5,4(1x,g20.12))')
     >             ip,d_lam,lp,lp_mod,lpn,lpn_mod

           if (ip.ne.id) write(0,*) ' ERROR IN.ne.ID:',ip,id

           lp_dlam(id) = d_lam
           lp_rawdata(in,id) = lp

        end do

 
        close(unum)

        cycle


 100    write(0,'(a)') 'STRING:'//search(1:searchlen)//
     >                ': NOT FOUND IN FILE:'//fname(1:filelen)
        close(unum)


        cycle


      enddo

      ! Data accumulated - process it


      frac1 = 0.95
      frac2 = 0.05
      base_parm1 = 1.201 * 0.1
      base_parm2 = 1.201 * 0.35

      lp_data = 0

      max_inst =0
      do in = 1,ncases

         do id = -maxpoints+1, maxpoints-1
c
            do ii =  -maxpoints+1, maxpoints-1
c
            lam0 = lp_dlam(id) ! id * lp_bin_width
            lam  = lp_dlam(ii) ! ii * lp_bin_width
c 
            dlam = lam-lam0
            
            ! Const 0.12A
            if (if_opt.eq.0) then 
               lp_data(in,id)=lp_data(in,id) +
     >            lp_rawdata(in,ii) * 
     >            const_inst_funct(dlam)
c
            ! Islere
            elseif (if_opt.eq.1) then 
               lp_data(in,id)=lp_data(in,id) +
     >            lp_rawdata(in,ii) * 
     >            isler_inst_funct(dlam)

            ! 95% Isler
            elseif (if_opt.eq.2) then 
               

               parm1 = base_parm1 * 0.95
               parm2 = base_parm2 * 0.95
               lp_data(in,id)=lp_data(in,id) +
     >            lp_rawdata(in,ii) * 
     >            isler_inst_funct2(dlam,frac1,parm1,frac2,parm2)

            ! 105% Isler
            elseif (if_opt.eq.4) then 

               parm1 = base_parm1 * 1.05
               parm2 = base_parm2 * 1.05
               lp_data(in,id)=lp_data(in,id) +
     >            lp_rawdata(in,ii) * 
     >            isler_inst_funct2(dlam,frac1,parm1,frac2,parm2)

            endif

         end do 
c
         max_inst = max(max_inst,lp_data(in,id))
c
      end do

      ! Normalize
      do id = -maxpoints,maxpoints
         lp_data(in,id) = lp_data(in,id)/max_inst
      enddo

      max_inst=0

      enddo


      ! Data accumulated - print it out
 

      open(onum,file=outname(1:outlen))

      write(onum,'(a,a)') 'CALCULATED LINE PROFILES: ',ifstring(1:iflen)
      write(onum,'(2a16,10a19)') 'DLAMBDA','LAMBDA',
     >         (cases(in)(1:len_trim(cases(in))),in=1,ncases)

      do id = -maxpoints,maxpoints

         write(onum,'(2(1x,f15.5),10(1x,g18.6)') 
     >      lp_dlam(id),lp_dlam(id)+baselam,
     >      (lp_data(in,id),in=1,ncases)

      end do

      close (onum)


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
