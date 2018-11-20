module mod_colours

  use mod_params

  !c     -*-Fortran-*-
  !c
  !c     Colour definitions for use in OUT
  !c
  !c
  !      integer maxncols,maxmarkers
  !      parameter(maxncols=255,maxmarkers=8)
  !      
  !c
  !      common /colours/ colour,ncols,defcol,icol,start_col
  !      save /colours/
  !c
  !      integer colour(maxncols),ncols,defcol,icol,col,start_col
  !c 
  !c     Set up selected markers for plotting
  !c
  !      integer marker(maxmarkers) 
  !c
  !      data marker /232,250,224,225,227,248,228,229/
  !c


  implicit none
  private


  !     Colour definitions for use in OUT

  integer,public::  maxncols,maxmarkers
  parameter(maxncols=255,maxmarkers=8)
  !
  integer,public:: colour(maxncols),ncols,defcol,icol,col,start_col
  ! 
  !     Set up selected markers for plotting
  !
  integer,public:: marker(maxmarkers) 
  !
  data marker /232,250,224,225,227,248,228,229/
  !

  public :: allocate_mod_colours, deallocate_mod_colours,init_col,setup_col,get_col,next_col

  integer, public :: colour_plot = 0

  
contains

  subroutine allocate_mod_colours
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)

  end subroutine allocate_mod_colours


  subroutine deallocate_mod_colours
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_colours

  integer function init_col()
    implicit none
    !
    !     Colour management routines
    !
    !     INIT_COL: Returns value of the first colour - initializes the
    !     internal counter.
    !
    !      include 'colours'
    !
    icol     = start_col
    !
    init_col = colour(icol)
    !
    return
  end function init_col

  integer function next_col()
    implicit none
    !
    !     NEXT_COL: Changes the colour pointer/counter to indicate
    !     the next colour and returns that colour
    !
    !      include 'colours'
    !
    icol = icol + 1
    if (icol.gt.ncols) icol = start_col
    !
    next_col = colour(icol)
    !
    return
  end function next_col

  integer function get_col()
    implicit none
    !
    !     GET_COL: Returns the colour at the current index/pointer
    !
    !      include 'colours'
    !
    get_col = colour(icol)

    return
  end function get_col



  subroutine setup_col(n_cols,opt)
    implicit none
    integer n_cols,opt
    !
    !      include 'colours'
    !
    !     SETUP_COL:
    !
    !     This routine assigns the colours in the colour map
    !
    !
    real hue, saturation,value
    integer in,init_col
    external init_col
    !
    ! toggle flag for colour plots
    !
    colour_plot = 1

    !
    !     Initialize number of colours and first colour
    !     n_cols = number of colours to initialize
    !     opt    = option to use for colour initialization
    !              1 = use colour set for gray scale - changing intensity
    !              2 = use colour set for colour plots - changing hue
    !              3 = use 15 colour B2/EIRENE setup
    !              4 = use 40 colour B2/EIRENE velocity plot setup 
    !              5 = special colour setup including an extended
    !                  colour table for Steve's plots  
    !
    ncols = n_cols
    !
    !     Set Black
    !     
    call rgb
    call colset(0.0,0.0,0.0,1)
    colour(1) = 1
    !     
    !     Set Background line drawing colour to black
    !     
    defcol = 1
    !
    !     Set start_col to first real colour - whatever index that is.
    !     
    start_col = 2
    !     
    !     Initialize ICOL to start_col
    !     
    icol = start_col
    !
    !     Set actual colour values
    !
    if ((opt.eq.0.or.opt.eq.1.or.opt.eq.2).and.n_cols.ne.16) then 
       !
       !        Select colour system to use - e.g. HSV, HSI, RGB ...
       !
       call hsi
       !
       !        Play with colour parameters until a reasonable selection is
       !        found.
       !
       do in = 2,ncols
          !
          !           Base colours on opt value
          !
          if (opt.eq.0.or.opt.eq.1) then
             !
             !              Y/G if plotted as colour
             !
             !              HUE:
             !
             !              hue = 1.1 - 0.7/(ncols-1) * (in-1)
             !
             hue = 0.5
             !
             !              SATURATION:
             !
             !              saturation = 1.0/(ncols-1) * (in-1)
             !              saturation = 0.6  + 0.4/(ncols-1) * (in-1)
             !
             saturation = 0.8
             !
             !              VALUE OR INTENSITY:
             !
             !              value = 0.33
             !
             value = 0.1 + 0.75/(ncols-1) * (in-1)

          elseif (opt.eq.2) then
             !
             !              HUE:
             !
             hue = 1.1 - 1.0/(ncols-1) * (in-1)
             !
             !              SATURATION:
             !
             !              saturation = 1.0/(ncols-1) * (in-1)
             !              saturation = 0.6  + 0.4/(ncols-1) * (in-1)
             !
             saturation = 0.8
             !
             !              VALUE OR INTENSITY:
             !
             !              value = 0.1 + 0.75/(ncols-1) * (in-1)
             !
             value = 0.8

          endif
          !
          !           Set colour
          !
          call colset(hue,saturation,value,in)

          colour(in) = in

       end do

    elseif ((opt.eq.0.or.opt.eq.1.or.opt.eq.2).and.n_cols.eq.16) then 
       !
       !        B2/EIRENE 15-color set  Krieger IPP/97
       !
       call rgb

       ncols=16

       do in=2,ncols
          colour(in)=in
       end do
       !
       !        BrightRed
       call colset(1.0, 0.0, 0.1, 16)
       !        Red
       call colset(0.9, 0.25, 0.0, 15)
       !        Orange
       call colset(1.0, 0.65, 0.0, 14)
       !        Golden
       call colset(1.0, 0.85, 0.0, 13)
       !        Yellow
       call colset(1.0, 1.0, 0.0, 12)
       !        GreenYellow
       call colset(0.7, 1.0, 0.2, 11)
       !        Chartreuse
       call colset(0.5, 1.0, 0.0, 10)
       !        Green
       call colset(0.2, 0.9, 0.1, 9)
       !        Aqua
       call colset(0.0, 0.9, 1.0, 8)
       !        DeepSkyBlue
       call colset(0.0, 0.75, 1.0, 7)
       !        RoyalBlue
       call colset(0.25, 0.45, 0.95, 6)
       !        SlateBlue
       call colset(0.4, 0.35, 0.8, 5)
       !        DarkViolet
       call colset(0.6, 0.0, 0.8, 4)
       !        Orchid
       call colset(0.85, 0.45, 0.8, 3)
       !        Lavender
       call colset(0.8, 0.8, 1.0, 2)

    elseif (opt.eq.3) then 
       !
       !        B2/EIRENE 15-color set  Krieger IPP/97
       !        IPP/01 Krieger - changed to 12 colors and reversed order
       !
       call rgb

       ncols=13

       do in=2,ncols
          colour(in)=in
       end do
       !
       !        strong colours
       call colset(1.00, 0.00, 0.00,  2)
       call colset(1.00, 0.40, 0.00,  3)
       call colset(1.00, 0.69, 0.00,  4)
       call colset(1.00, 0.97, 0.00,  5)
       call colset(0.72, 1.00, 0.28,  6)
       call colset(0.44, 1.00, 0.56,  7)
       call colset(0.15, 1.00, 0.84,  8)
       call colset(0.00, 0.89, 1.00,  9)
       call colset(0.00, 0.61, 1.00, 10)
       call colset(0.00, 0.33, 1.00, 11)
       call colset(0.00, 0.00, 1.00, 12)
       call colset(0.33, 0.00, 0.84, 13)
       !
       !        BrightRed
       !         call colset(1.0, 0.0, 0.1, 16)
       !        Red
       !         call colset(0.9, 0.25, 0.0, 15)
       !        Orange
       !         call colset(1.0, 0.65, 0.0, 14)
       !        Golden
       !         call colset(1.0, 0.85, 0.0, 13)
       !        Yellow
       !         call colset(1.0, 1.0, 0.0, 12)
       !        GreenYellow
       !         call colset(0.7, 1.0, 0.2, 11)
       !        Chartreuse
       !         call colset(0.5, 1.0, 0.0, 10)
       !        Green
       !         call colset(0.2, 0.9, 0.1, 9)
       !        Aqua
       !         call colset(0.0, 0.9, 1.0, 8)
       !        DeepSkyBlue
       !         call colset(0.0, 0.75, 1.0, 7)
       !        RoyalBlue
       !         call colset(0.25, 0.45, 0.95, 6)
       !        SlateBlue
       !         call colset(0.4, 0.35, 0.8, 5)
       !        DarkViolet
       !         call colset(0.6, 0.0, 0.8, 4)
       !        Orchid
       !         call colset(0.85, 0.45, 0.8, 3)
       !        Lavender
       !         call colset(0.8, 0.8, 1.0, 2)
       !
    elseif (opt.eq.4) then
       !
       !        B2/EIRENE velocity-color set  Krieger IPP/97
       !
       call rgb

       ncols=41

       do in=2,ncols
          colour(in) = in  
       end do

       do in=1,ncols/2
          call colset( 1.-(1.-.00)*real(in-1)/real(ncols/2-1),0.0, 0.0,1+in)
       end do

       do in=1,ncols/2
          call colset(0.0,0.0+(1.-.00)*real(in-1)/real(ncols/2-1),0.0,1+in+ncols/2)
       end do

    elseif (opt.eq.5) then 


       ! slmod begin - temp
       !...     Temporary colour setup over-ride:
       !
       !          WRITE(0,*) 'COLOUR SETUP BEING OVERWRITTEN'
       !
       CALL HSI
       CALL ColSet(0.0,0.0,0.0,1)
       colour(1) = 1
       DO in = 2, ncols
          colour(in) = in
          hue        = 0.0
          saturation = 0.0
          value      = (REAL(in-1) / REAL(ncols))**0.5
          CALL ColSet(hue,saturation,value,in)
       ENDDO
       !
       !....     "Extended" colour set:
       DO in = ncols+1, (ncols+1) + 12
          colour(in) = in
       ENDDO


       CALL RGB
       CALL ColSet(1.0,0.0,0.0,ncols+1)
       CALL ColSet(1.0,0.0,0.0,ncols+2)
       CALL ColSet(0.0,1.0,0.0,ncols+3)
       CALL ColSet(0.0,0.0,1.0,ncols+4)
       CALL ColSet(0.5,0.7,0.2,ncols+5)
       CALL ColSet(0.5,0.0,0.5,ncols+6)
       CALL ColSet(0.0,0.5,0.5,ncols+7)

       !
       !         Add some more colours to the "extended" colour set
       !

       !         Orange
       call colset(1.0, 0.65, 0.0,ncols+8)

       !         Yellow
       call colset(1.0, 1.0, 0.0,ncols+9)

       !         SlateBlue
       call colset(0.4, 0.35, 0.8,ncols+10)

       !         DeepSkyBlue
       call colset(0.0, 0.75, 1.0,ncols+11)

       !         DarkViolet
       call colset(0.6, 0.0, 0.8,ncols+12)

       CALL HSV
       CALL ColSet(0.33,1.0,1.0,ncols+1)
       CALL ColSet(0.33,1.0,1.0,ncols+2)
       CALL ColSet(0.00,1.0,1.0,ncols+3)
       CALL ColSet(0.67,1.0,1.0,ncols+4)

       CALL ColSet(0.41,1.0,1.0,ncols+5)
       CALL ColSet(0.17,1.0,0.7,ncols+6)
       CALL ColSet(0.83,1.0,1.0,ncols+7)

       CALL ColSet(0.00,0.5,1.0,ncols+8)
       CALL ColSet(0.45,1.0,0.8,ncols+9)
       CALL ColSet(0.67,0.5,0.7,ncols+10)




       !
       ! slmod end
       !
    endif
    !
    return
  end subroutine setup_col


end module mod_colours
