module process_irtv

  use saf_files
  implicit none

  integer :: xpixel,ypixel,nframes_img,data_type,data_len

  character*256 :: reference_string
  character*256 :: frame_string

  integer :: image_count,analysis_count
  character*256 :: input_filename

  integer,allocatable :: frame_set(:)

contains





  subroutine open_irtv(filename,ierr)
    implicit none
    character*(*) :: filename
    integer :: ierr

    ! Open SAF file

    call init_saf_file(filename,ierr)

    if (ierr.ne.0) then 
       write(0,*) 'Error opening SAF file - Exiting'
       return
    endif


    ! Read SAF file header

    call read_saf_file_header(xpixel,ypixel,nframes_img,data_type,data_len,reference_string,ierr)

    ! Load SAF data 

    call read_saf_data(xpixel,ypixel)

    if (ierr.ne.0) then 
       write(0,*) 'Error reading SAF file header - Exiting'
       return
    endif

    image_count = 0
    analysis_count=0
    input_filename = filename

  end subroutine open_irtv


  subroutine process_irtv_commands
    implicit none


    logical :: finished,analysis_finished,background_finished
    integer :: opt,ierr,subopt,start_frame,end_frame
    integer :: pixel_width,start_pixel_x,start_pixel_y,end_pixel_x,end_pixel_y
    real,allocatable :: frame(:,:),bg_frame(:,:)
    real :: calibration, minvalue,maxvalue

    integer :: analysis_opt,background_opt
    integer :: tmp_funit
    integer :: in
    integer :: sample_window(2,2)


    integer :: frame_opt

    ! primitive UI driver loop to execute commands related to processing the IRTV data
    !
    ! Options : 
    !  1 - clear buffer'
    !  2 - Load'
    !    A - single frame
    !    B - multiple summed together
    !  3 - Save'
    !    A - file
    !    B - jpeg
    !  4 - Analyse'
    !    A - Sample Window Average
    !    B - Cross-sections
    !    C - Define Sample Window (default=full frame)
    !  5 - Background
    !    A - recalculate background
    !    B - rescale background
    !  7 - exit'
    !

    !
    ! Initialization and defaults
    !

    finished = .false.
    calibration = 0.0
    minvalue = 0.0
    maxvalue = 0.0


    ! set default sample window
    call set_sample_window(sample_window,1) 


    allocate(frame(xpixel,ypixel),stat=ierr)
    allocate(bg_frame(xpixel,ypixel),stat=ierr)


    ! calculate background

    call calculate_bg_frame(bg_frame,xpixel,ypixel)


    ! main loop

    do while (.not.finished)

       call write_option_list
       read(5,*) opt

       if (opt.eq.1) then 
          ! zero out the frame buffer
          frame = 0.0

       elseif (opt.eq.2) then 

          ! load buffer


          call load_buffer(frame,bg_frame,xpixel,ypixel)


       elseif (opt.eq.3) then 
          ! save buffer

          write(6,'(a)') 'Save frame buffer to file: 1-formatted text   2-jpeg'
          read(5,*) subopt

          ! need to set these every time 
          calibration = 0.0
          minvalue = 0.0
          maxvalue = 0.0

          if (subopt.eq.1) then

             call save_text_image(frame,xpixel,ypixel,calibration,minvalue,maxvalue)


          elseif (subopt.eq.2) then

             call save_jpeg_image(frame,xpixel,ypixel,calibration,minvalue,maxvalue)

          endif


       elseif (opt.eq.4) then 
          ! analyse buffer

          analysis_finished = .false.

          do while (.not.analysis_finished) 
             call write_analyse_list
             read(5,*) analysis_opt

             if (analysis_opt.eq.1) then 

                call characterize_frames(bg_frame,xpixel,ypixel,sample_window)

             elseif (analysis_opt.eq.2) then 

                call calculate_cross_section(frame,xpixel,ypixel,0)

             elseif (analysis_opt.eq.3) then 

                call set_sample_window(sample_window,0)

             elseif (analysis_opt.eq.4) then 

                analysis_finished = .true.
             endif


          end do

       elseif (opt.eq.5) then 

          background_finished = .false.

          do while (.not.background_finished) 

             call write_background_list
             read(5,*) background_opt

             if (background_opt.eq.1) then 


                call calculate_bg_frame(bg_frame,xpixel,ypixel)


             elseif (background_opt.eq.2) then 


                call scale_bg_frame(bg_frame,xpixel,ypixel)


             elseif (background_opt.eq.3) then 
                background_finished = .true.
             endif


          end do



       elseif (opt.eq.7) then 
          ! exit
          finished = .true.
       endif


    end do

    !
    ! Deallocate any storage allocated for SAF data
    !
    call empty_saf_data


  end subroutine process_irtv_commands






  subroutine write_option_list
    implicit none

    write(6,'(a)') '1 - clear buffer'
    write(6,'(a)') '2 - Load'
    write(6,'(a)') '3 - Save'
    write(6,'(a)') '4 - Analyse'
    write(6,'(a)') '5 - Background'
    write(6,'(a)') '7 - exit'


  end subroutine write_option_list


  subroutine write_analyse_list
    implicit none

    write(6,'(a)') '1 - print average and delta for sample window'
    write(6,'(a)') '2 - integrate cross-sections'
    write(6,'(a)') '3 - Define sample window (default=full frame)'
    write(6,'(a)') '4 - exit'


  end subroutine write_analyse_list

  subroutine write_background_list
    implicit none

    write(6,'(a)') '1 - recalculate background'
    write(6,'(a)') '2 - scale background'
    write(6,'(a)') '3 - exit'

  end subroutine write_background_list


  subroutine get_frame_set(frame_opt,nframes)
    implicit none

    integer :: frame_opt,nframes

    integer :: in,start_frame,end_frame


    write(6,'(a,i6)') 'Enter frame numbers: two integers : if identical then only a single frame loaded: Range= 1,',nframes_img
    write(6,'(a)')    'If the first frame number is -1 then the second number is the quantity of discrete frames to enter (you will be prompted)'

    read(5,*) start_frame,end_frame

    frame_opt = 0

    if (start_frame.eq.-1) then

       frame_opt = 3
       nframes = end_frame


       ! allocate frame_set to the required size
       if (allocated(frame_set)) then 
          deallocate(frame_set)
       endif

       allocate(frame_set(nframes))

       write(6,*) 'Enter the list of frames:'
       read(5,*) (frame_set(in),in=1,nframes)

       do in = 1,nframes
          if (frame_set(in).lt.1.or.frame_set(in).gt.nframes_img) then 
             write(6,'(a,i6,a,i6,a)') 'Specified frame with index ',in,' and value ',frame_set(in),' is out of range. Exiting.' 
             frame_opt=0
             deallocate(frame_set)
             frame_string = 'ERROR'
             return
          endif
       end do


    elseif (start_frame.lt.0.or.end_frame.lt.0.or.start_frame.gt.nframes_img.or.end_frame.gt.nframes_img) then 
       frame_opt = 0
       nframes = 0
       write (6,'(2(a,2i10))') 'Frame specified is out of range. Specified = ',start_frame,end_frame,' Available =',1,nframes_img
       frame_string = 'ERROR'
       return
    elseif (start_frame.eq.end_frame) then 

       frame_opt = 1 
       nframes = 1

       ! allocate frame_set to the required size
       if (allocated(frame_set)) then 
          deallocate(frame_set)
       endif

       allocate(frame_set(nframes))


       frame_set(1) = start_frame

    else
       frame_opt = 2 
       nframes = 2

       ! allocate frame_set to the required size
       if (allocated(frame_set)) then 
          deallocate(frame_set)
       endif

       allocate(frame_set(nframes))

       frame_set(1) = start_frame
       frame_set(2) = end_frame


    endif

    if (frame_opt.eq.1) then 
       write(frame_string,'(a,i1,a,i,a,i)') 'FR_',frame_opt,'_',frame_set(1)
    elseif (frame_opt.eq.2) then 
       write(frame_string,'(a,i1,a,i,a,i)') 'FR_',frame_opt,'_',frame_set(1),'-',frame_set(2)
    elseif (frame_opt.eq.3) then 
       write(frame_string,'(a,i1,a,i,a,i,a,i)') 'FR_',frame_opt,'_',nframes,'_',frame_set(1),'+..',frame_set(nframes)       
    endif


  end subroutine get_frame_set



  subroutine load_buffer(frame,bg_frame,xpixel,ypixel)
    implicit none
    integer :: xpixel,ypixel
    real :: frame(xpixel,ypixel),bg_frame(xpixel,ypixel)

    integer :: frame_opt,frame_count,frame_to_load
    integer :: in,subopt

    real, allocatable :: temp_frame(:,:),adj_bg_frame(:,:)
    integer :: bg_cnt, target_frame, nframes


    write(6,'(a)') 'Load general frame=1 , Load single frame with adjacent bg subtraction=2'
    read(5,*) subopt

    if (subopt.eq.1) then 

       write(6,'(a,i6)') 'Enter frame numbers to average and store in buffer :'

       call get_frame_set(frame_opt,nframes)

       if (frame_opt.eq.0) then 
          write(6,*) 'Error in frame specification. Exiting'
          frame = 0.0
          return
       elseif (frame_opt.eq.1.or.frame_opt.eq.3) then 
          ! nframes is assigned in get_frame_set
          frame_count=nframes
       elseif (frame_opt.eq.2) then 
          frame_count = frame_set(2)-frame_set(1) + 1
       endif

       allocate(temp_frame(xpixel,ypixel))

       frame = 0.0
       temp_frame = 0.0

       do in = 1,frame_count

          if (frame_opt.eq.1) then 
             frame_to_load = frame_set(1)
          elseif (frame_opt.eq.2) then
             frame_to_load = frame_set(1) + in -1
          elseif (frame_opt.eq.3) then
             frame_to_load = frame_set(in)
          endif

          call get_saf_frame(temp_frame,xpixel,ypixel,frame_to_load)
          frame = frame + temp_frame - bg_frame

       end do

       frame = frame / (real(frame_count))


    elseif (subopt.eq.2) then 

       write(6,'(a,i6)') 'Enter frame number:'
       read(5,*) target_frame

       write(frame_string,'(a,i)') 'FR_',target_frame

       ! background is the average of the adjacent frames

       allocate(temp_frame(xpixel,ypixel))
       allocate(adj_bg_frame(xpixel,ypixel))

       frame = 0.0
       adj_bg_frame = 0.0
       bg_cnt = 0

       if (target_frame.gt.1.and.target_frame.le.nframes_img) then 
          call get_saf_frame(temp_frame,xpixel,ypixel,target_frame-1)
          bg_cnt = bg_cnt+1
          adj_bg_frame = adj_bg_frame + temp_frame
       endif

       if (target_frame.ge.1.and.target_frame.lt.nframes_img) then 
          call get_saf_frame(temp_frame,xpixel,ypixel,target_frame+1)
          bg_cnt = bg_cnt+1
          adj_bg_frame = adj_bg_frame + temp_frame
       endif

       if (bg_cnt.gt.0) then 
          adj_bg_frame = adj_bg_frame / real(bg_cnt)
       else
          write(6,*) 'ERROR: Specified frame out of range:',target_frame
          return
       endif

       ! load target frame and subtract adjacent background
       call get_saf_frame(temp_frame,xpixel,ypixel,target_frame)
       frame = temp_frame - adj_bg_frame


    endif

    if (allocated(temp_frame)) deallocate(temp_frame)
    if (allocated(adj_bg_frame)) deallocate(adj_bg_frame)



  end subroutine load_buffer


  subroutine calculate_bg_frame(bg_frame,xpixel,ypixel)
    implicit none
    integer :: xpixel,ypixel
    real :: bg_frame(xpixel,ypixel)

    integer :: frame_opt,frame_count,frame_to_load
    integer :: in, nframes

    real, allocatable :: temp_frame(:,:)

    write(6,'(a,i6)') 'Enter frame numbers for background :'
    write(6,'(a,i6)') 'Out of range values assign a background of 0.0'

    call get_frame_set(frame_opt,nframes)

    bg_frame = 0.0

    if (frame_opt.eq.0) then 
       write(6,*) 'Error in BG frame specification. Exiting'
       return
    elseif (frame_opt.eq.1.or.frame_opt.eq.3) then 
       ! nframes is assigned in get_frame_set
       frame_count=nframes
    elseif (frame_opt.eq.2) then 
       frame_count = frame_set(2)-frame_set(1) + 1
    endif

    allocate(temp_frame(xpixel,ypixel))

    do in = 1,frame_count

       if (frame_opt.eq.1) then 
          frame_to_load = frame_set(1)
       elseif (frame_opt.eq.2) then
          frame_to_load = frame_set(1) + in -1
       elseif (frame_opt.eq.3) then
          frame_to_load = frame_set(in)
       endif

       call get_saf_frame(temp_frame,xpixel,ypixel,frame_to_load)

       bg_frame = bg_frame + temp_frame 

    end do

    bg_frame = bg_frame / (real(frame_count))

    if (allocated(temp_frame)) deallocate(temp_frame)


  end subroutine calculate_bg_frame



  subroutine scale_bg_frame(bg_frame,xpixel,ypixel)
    implicit none

    !
    ! Routine is used to apply a uniform factor to the background frame. This value is either specified explicitly or is 
    ! calculated from the ratio of average pixel intensity for two specified frames.
    !
    !

    integer :: xpixel,ypixel
    real :: bg_frame(xpixel,ypixel)

    ! scale background frame
    real :: scalef1, scalef2
    real,allocatable :: temp_frame(:,:)

    real :: avg1, avg2
    integer :: frame1,frame2

    write(6,*) 'Enter scaling factor: (a) Number  -1   (b) frame1   frame 2 '
    write(6,*) 'Scale = Number  OR   avg_intensity_frame1/avg_intensity_frame2'

    read(5,*) scalef1,scalef2 


    if (scalef2.le.0) then 

       bg_frame = bg_frame * scalef1

       write(6,'(a,g12.5)') 'Scale factor applied to background frame = ',scalef1

    else

       allocate(temp_frame(xpixel,ypixel))

       frame1 = int(scalef1)
       frame2 = int(scalef2)


       call get_saf_frame(temp_frame,xpixel,ypixel,frame1)

       avg1 = real(sum(temp_frame)) / (real(xpixel)*real(ypixel))

       call get_saf_frame(temp_frame,xpixel,ypixel,frame2)

       avg2 = real(sum(temp_frame)) / (real(xpixel)*real(ypixel))


       if (allocated(temp_frame)) deallocate(temp_frame)

       scalef1 = avg1/avg2

       bg_frame = bg_frame * scalef1

       write(6,'(a,g12.5)') 'Scale factor applied to background frame = ',scalef1
       write(6,'(2(a,i5,a,1x,g12.5))') 'Frame 1: ',frame1,' Avg = ',avg1,'  Frame 2: ',frame2,' Avg = ',avg2  


    endif


  end subroutine scale_bg_frame


  recursive subroutine set_sample_window(sample_window,opt)
    implicit none
    integer :: sample_window(2,2)
    integer :: x1,x2,y1,y2
    integer :: opt


    if (opt.eq.0) then 
       write(6,'(a)') 'Enter sample window corners: x1  y1  x2  y2  (x2>x1 and y2>y1)'
       write(6,'(a)') 'Enter: -1 0 0 0  for the entire frame'
       write(6,'(a)') 'Enter: -2 0 0 0  for sample lower divertor floor sample region'
       read(5,*) x1,y1,x2,y2

       if (x1.lt.0) then 
          call set_sample_window(sample_window,abs(x1))

       else

          if (x2.le.x1.or.y2.le.y1) then 
             write (6,'(a)') 'Input Error: X2 must be > X1 and Y2 must be > Y1 : Sample window NOT changed'
          else
             sample_window(1,1) = x1
             sample_window(1,2) = y1
             sample_window(2,1) = x2
             sample_window(2,2) = y2
          endif

       endif
    elseif (opt.eq.1) then 
       ! set default sampling window 1 - full image
       sample_window(1,1) = 1
       sample_window(1,2) = 1
       sample_window(2,1) = xpixel
       sample_window(2,2) = ypixel

    elseif (opt.eq.2) then 
       ! default sample window 2 - lower divertor floor section along plotting cross section
       sample_window(1,1) = 330
       sample_window(1,2) = 212
       sample_window(2,1) = 360
       sample_window(2,2) = 220
    endif


  end subroutine set_sample_window


  subroutine characterize_frames(bg_frame,xpixel,ypixel,sample_window)
    implicit none
    integer,intent(in) :: xpixel,ypixel
    real :: bg_frame(xpixel,ypixel)
    integer :: sample_window(2,2)

    ! local variables
    real,allocatable :: temp_frame(:,:)
    integer :: in
    integer :: tmp_funit,slen

    character*(256) fname
    real :: average

    integer :: x1,x2,y1,y2,dx,dy
    integer :: frame_opt

    integer :: frame_count,cnt,frame_to_load, nframes
    real,allocatable :: frame_average(:),frame_delta(:),frame_index(:)

    !real :: previous_average

    write(6,'(a)') 'Enter frame numbers for data characterization:'
    call get_frame_set(frame_opt,nframes)


    ! set easier aliases for the sampling window

    x1 = sample_window(1,1)
    x2 = sample_window(2,1)
    y1 = sample_window(1,2)
    y2 = sample_window(2,2)
    dx = x2-x1+1
    dy = y2-y1+1

    !previous_average = 0.0

    if (frame_opt.eq.0) then
       write(6,'(a)') 'Error in frame specifier entry. Exiting.'
       return
    elseif (frame_opt.eq.1.or.frame_opt.eq.3) then 
       ! nframes is assigned in get_frame_set
       frame_count=nframes
    elseif (frame_opt.eq.2) then 
       frame_count = frame_set(2)-frame_set(1) + 1
    endif

    ! allocate local storage
    allocate(temp_frame(xpixel,ypixel))

    allocate(frame_average(frame_count))
    allocate(frame_delta(frame_count))
    allocate(frame_index(frame_count))

    frame_delta = 0.0


    analysis_count = analysis_count + 1

    slen = len_trim(input_filename)

    ! Create data file name
    call create_file_name(fname,input_filename(1:slen),analysis_count,3)

    ! get free unit number
    call find_free_unit_number(tmp_funit)

    ! open file
    open(tmp_funit,file=fname,form='formatted')

    ! load frame data
    do in = 1,frame_count

       if (frame_opt.eq.1) then 
          frame_to_load = frame_set(1)
       elseif (frame_opt.eq.2) then
          frame_to_load = frame_set(1) + in -1
       elseif (frame_opt.eq.3) then
          frame_to_load = frame_set(in)
       endif

       ! get frame
       call get_saf_frame(temp_frame,xpixel,ypixel,frame_to_load)

       ! background frame has been previously set up
       temp_frame = temp_frame - bg_frame

       average = real(sum(temp_frame(x1:x2,y1:y2))) / (real(dx)*real(dy))

       frame_index(in) = frame_to_load
       frame_average(in) = average

       !frame_delta(in) = average-previous_average
       !previous_average = average

    end do



    ! Analyse the results - calculate the average delta

    do in = 2,frame_count-1

       frame_delta(in) = ((frame_average(in)-frame_average(in-1)) + &
               &         (frame_average(in)-frame_average(in+1)))/2.0

    end do


    ! write out the results

    write(6,'(a)') 'Input filename :',trim(input_filename)
    write(6,'(a)') 'FRAMES :',trim(frame_string)

    write(6,'(a,i6,a,i6,a,i6,a,i6,a)') 'Sample_window:  (',int(x1),',', int(y1),')  to  (',&
         &int(x2),',',int(y2),')'

    write(tmp_funit,'(a)') 'Input filename :',trim(input_filename)
    write(tmp_funit,'(a)') 'FRAMES :',trim(frame_string)
    write(tmp_funit,'(a,i6,a,i6,a,i6,a,i6,a)') 'Sample_window:  (',int(x1),',', int(y1),')  to  (',&
         &int(x2),',',int(y2),')'

    do in = 1,frame_count

       write(6,'(a,a,i6,5(2x,g18.6))') 'NAME,FRAME,AVE,DELTA: ',trim(input_filename),&
            &in,frame_index(in),frame_average(in),frame_delta(in)
       write(tmp_funit,'(a,a,i6,5(2x,g18.6))') 'NAME,FRAME,AVE,DELTA: ',trim(input_filename),&
            &in,frame_index(in),frame_average(in),frame_delta(in)

    end do





    ! deallocate any storage used
    if (allocated(temp_frame)) deallocate (temp_frame)
    if (allocated(frame_index)) deallocate(frame_index)
    if (allocated(frame_average)) deallocate(frame_average)
    if (allocated(frame_delta)) deallocate(frame_delta)

    close(tmp_funit)


  end subroutine characterize_frames


  subroutine calculate_cross_section(frame,xpixel,ypixel,opt)
    implicit none
    integer :: xpixel, ypixel,opt
    real :: frame(xpixel,ypixel)

    !
    ! This routine will calculate cross-sections across an image
    ! Input data required: start pixel x, start pixel y
    !                      end pixel x, end pixel y
    !                      width of strip for averaging in pixels
    !
    !

    real :: start_pixel_x,start_pixel_y,end_pixel_x,end_pixel_y,strip_width_in_pixels

    real,allocatable :: profile_axis(:),profile(:),profile_r_axis(:)

    real :: strip_length

    real :: slope

    real :: delta_x_pixel,delta_y_pixel,strip_delta_x_pixel,strip_delta_y_pixel

    integer :: in,ic,nsteps

    real :: bix,biy,ix,iy

    real :: tmp_val

    integer :: ierr,tmp_funit,slen

    character*(256) fname


    call get_cross_section_input(start_pixel_x,start_pixel_y,end_pixel_x,end_pixel_y,strip_width_in_pixels,opt)

    if (start_pixel_x.lt.1.or.start_pixel_x.gt.xpixel.or. &
         end_pixel_x.lt.1.or.end_pixel_x.gt.xpixel.or.&
         start_pixel_y.lt.1.or.start_pixel_y.gt.xpixel.or. &
         end_pixel_y.lt.1.or.end_pixel_y.gt.xpixel) then 
       write(6,*) 'Invalid pixels specified: exiting'
       return
    endif

    analysis_count = analysis_count + 1


    ! Calculate the profile along the cross-section


    strip_length = sqrt((real(end_pixel_x-start_pixel_x))**2 + (real(end_pixel_y-start_pixel_y))**2)

    ! will step along the strip at unit (1 pixel intervals) ... for motion along an axis this is then easy. 

    delta_x_pixel = (end_pixel_x - start_pixel_x) / strip_length
    delta_y_pixel = (end_pixel_y - start_pixel_y) / strip_length

    strip_delta_x_pixel = delta_y_pixel
    strip_delta_y_pixel = -delta_x_pixel

    nsteps = int(strip_length) +1

    allocate (profile_axis(nsteps),stat=ierr)
    allocate (profile_r_axis(nsteps),stat=ierr)
    allocate (profile(nsteps),stat=ierr)

    do in = 1,nsteps

       bix = start_pixel_x + real((in-1)) * delta_x_pixel
       biy = start_pixel_y + real((in-1)) * delta_y_pixel

       ! sum across strip
       tmp_val = 0.0
       do ic = -strip_width_in_pixels,strip_width_in_pixels
          ix = bix + real(ic) * strip_delta_x_pixel
          iy = biy + real(ic) * strip_delta_y_pixel
          tmp_val = tmp_val + frame(knint(ix),knint(iy))

          !write(6,'(a,2i8,4f10.2,2i8)') 'Pixel indices:',in,ic,bix,biy,ix,iy,knint(ix),knint(iy)
       end do

       ! divide to get average
       tmp_val = tmp_val / (2.0*strip_width_in_pixels+1.0)


       profile_r_axis(in) = calibrated_r_pixel(bix)
       profile_axis(in) = bix
       profile(in) = tmp_val

    end do


    ! output the profile

    slen = len_trim(input_filename)
    ! Create data file name
    call create_file_name(fname,input_filename(1:slen),analysis_count,3)

    ! get free unit number
    call find_free_unit_number(tmp_funit)

    ! open file
    open(tmp_funit,file=fname,form='formatted')



    write(6,'(a)') 'Input filename :',trim(input_filename)
    write(6,'(a)') 'FRAMES :',trim(frame_string)
    write(6,'(a,i6,a,i6,a,i6,a,i6,a)') 'Pixel range selected:  (',int(start_pixel_x),',', int(start_pixel_y),')  to  (',&
         &int(end_pixel_x),',',int(end_pixel_y),')'
    write(6,'(a,i6)') 'Averaging width of strip = ',int(strip_width_in_pixels)


    write(tmp_funit,'(a)') 'Input filename :',trim(input_filename)
    write(tmp_funit,'(a)') 'FRAMES :',trim(frame_string)
    write(tmp_funit,'(a,i6,a,i6,a,i6,a,i6,a)') 'Pixel range selected:  (',int(start_pixel_x),',', int(start_pixel_y),')  to  (',&
         &int(end_pixel_x),',',int(end_pixel_y),')'
    write(tmp_funit,'(a,i6)') 'Averaging width of strip = ',int(strip_width_in_pixels)


    do in = nsteps,1,-1

       ! write data to screen and data file - write in reverse order so that r axis is correct
       write(6,'(a,i6,5(2x,g18.6))') 'CROSS-SECTION:',in,profile_axis(in),profile_r_axis(in),profile(in)
       write(tmp_funit,'(a,i6,5(2x,g18.6))') 'CROSS-SECTION:',in,profile_axis(in),profile_r_axis(in),profile(in)

    end do

    close(tmp_funit)

    if (allocated(profile_r_axis)) deallocate(profile_r_axis)
    if (allocated(profile_axis)) deallocate(profile_axis)
    if (allocated(profile)) deallocate(profile)

  end subroutine calculate_cross_section

  real function calibrated_r_pixel(rix)
    implicit none
    real :: rix
    ! these calibration data are for IRTV camera 6601 for shots greater than 133691
    ! they apply only to rows 212 to 220 of the image

    real,parameter :: offset = 198.46 / 100.0  ! scaling in meters from cm
    real,parameter :: pixel_scale = -0.19621 / 100.0  ! scaling in meters from cm

    calibrated_r_pixel = offset + pixel_scale * rix

  end function calibrated_r_pixel


  recursive subroutine get_cross_section_input(start_pixel_x,start_pixel_y,end_pixel_x,end_pixel_y,strip_width_in_pixels,opt)
    implicit none
    real :: start_pixel_x,start_pixel_y,end_pixel_x,end_pixel_y,strip_width_in_pixels
    integer :: opt

    real,parameter :: left_image_bound=95, right_image_bound=584
    real,parameter :: strip_center = 216.0
    real,parameter :: default_strip_width = 1.0 ! 3 pixel wide strip for averaging

    real,parameter :: left_image_bound_elm=1, right_image_bound_elm=464
    real,parameter :: strip_center_elm_image = 2.0

    if (opt.eq.0) then 

       write(6,'(a)') 'Five input data values required: start_pixel_x start_pixel_y end_pixel_x end_pixel_y  strip_half_width_in_pixels'
       write(6,'(a)') 'Enter: -1 0 0 0 0   for default full image strip'
       write(6,'(a)') 'Enter: -2 0 0 0 0   for default elm  image strip'
       write(6,'(2(a,2i8))') 'X pixel range: ',1,xpixel,'   Y pixel range:',1,ypixel
       read(5,*) start_pixel_x,start_pixel_y,end_pixel_x,end_pixel_y,strip_width_in_pixels

       if (start_pixel_x.lt.0) then
          call get_cross_section_input(start_pixel_x,start_pixel_y,end_pixel_x,end_pixel_y,strip_width_in_pixels,abs(int(start_pixel_x)))
          return
       endif

    elseif (opt.eq.1) then 
       ! define default strip for analysis (standard cross-section of image)
       ! This runs from left_bound to right_bound for the pixels where the calibration applies

       start_pixel_x = left_image_bound
       start_pixel_y = strip_center
       end_pixel_x = right_image_bound
       end_pixel_y = strip_center
       strip_width_in_pixels = default_strip_width

    elseif (opt.eq.2) then 
       ! define default strip for analysis (cross-section of elm image data which is only 4 pixels across)
       ! This runs from left_bound to right_bound for the pixels where the calibration applies

       start_pixel_x = left_image_bound_elm
       start_pixel_y = strip_center_elm_image
       end_pixel_x = right_image_bound_elm
       end_pixel_y = strip_center_elm_image
       strip_width_in_pixels = default_strip_width

    endif

    write(6,'(a,5(1x,g12.5)') 'Strip parameters:',start_pixel_x,start_pixel_y,end_pixel_x,end_pixel_y,strip_width_in_pixels


  end subroutine get_cross_section_input


  subroutine save_jpeg_image(image,xres,yres,calibration,maxvalue,minvalue)
    implicit none
    integer xres,yres
    real :: image(xres,yres),calibration
    real ::  maxvalue,minvalue
    !
    !     Set up and write out a JPEG formatted image
    !
    !
    integer temp_val 
    integer,parameter::SHORT=1
    !
    integer*1,allocatable ::  image_buffer(:,:)  

    integer :: ierr,slen,ic,ir
    real,external :: max

    !integer image_count 
    !data image_count /0/

    character*(256) fname

    real :: local_maxvalue, local_minvalue,local_calibration

    !
    !     Assess inputs
    !
    if (maxvalue.eq.0.0) then 
       local_maxvalue = maxval(image)
    endif

    !
    !---------------------------------------------------------------------- 
    !
    !     Allocate space for image and check to see if space available 
    !
    !      write(0,*) 'Allocating IMAGE_BUFFER:' 
    !
    allocate (image_buffer(xres,yres),STAT=ierr)
    ! 
    if (ierr.ne.0) then 
       write(0,*) 'IMAGE_BUFFER array could not be allocated: ',xres,yres
       write(6,*) 'IMAGE_BUFFER array could not be allocated: ',xres,yres
       return
    endif
    !
    !      write(0,*) 'IMAGE_BUFFER Allocated' 
    !
    !     If a calibration is not given for the camera - scale the 
    !     maximum signal to maximum intensity
    !
    !      calibration = 3.44E+22 / 255.0 / 4.0 / 3.1415
    !
    !
    if (calibration.eq.0.0) then 
       WRITE(0,*) 'CALIBRATION:',local_maxvalue
       local_calibration = local_maxvalue/255.0
    else
       WRITE(0,*) 'CALIBRATION:',calibration/255.0
       local_calibration = calibration/255.0
    endif
    !
    !     Loop through image to generate scaled image_buffer
    !
    do ic = 1,xres
       do ir = 1,yres

          temp_val = int(image(ic,ir)/local_calibration)
          if (temp_val.gt.255) temp_val = 255
          if (temp_val.lt.0) temp_val = 0
          image_buffer(ic,ir) = int(temp_val,SHORT)

          !            WRITE(0,*) 'IMAGE:',ic,ir,image(ic,ir)

       end do
    end do
    !
    !     Determine file name for image
    !
    image_count = image_count+1
    !
    !...  Load base case name from environment variable CASENAME:
    !  CALL CaseName(filename,ierr)
    !        
    slen = len_trim(input_filename)
    !
    call create_file_name(fname,input_filename(1:slen),image_count,1)
    !
    !      len = lenstr(fname)
    !      write(6,*) 'File:',fname(1:len) 
    !
    !     Write out actual JPG file - using quality value of 100 for highest
    !     quality. This routine is written in C and accesses the publicly
    !     available jpegsrc library. The code is based on the example.c
    !     file shipped with that library. 
    !

    slen = len_trim(fname)
    !
    !     For the SUN implementation - the null termination is added in the 
    !     C code since the Fortran compiler seemed to have trouble figuring
    !     out what was required.
    !
    call write_JPEG_file (fname,%VAL(slen),image_buffer,%VAL(yres),%VAL(xres), %VAL(100),%VAL(1))
    !
    !
    !---------------------------------------------------------------------- 
    !
    !     Free space assigned to image_buffer array
    !
    if (allocated(image_buffer)) deallocate(image_buffer)
    !
    !---------------------------------------------------------------------- 

    return 
  end subroutine save_jpeg_image

  subroutine save_text_image(image,xres,yres,calibration,maxvalue,minvalue)
    implicit none
    integer xres,yres
    real :: image(xres,yres),calibration
    real ::  maxvalue,minvalue
    integer :: funit 

    character*(256) fname
    integer ::  image_count 
    integer :: ic,ir,slen

    !
    image_count = image_count+1
    !
    !...  Load base case name from environment variable CASENAME:
    !  CALL CaseName(filename,ierr)
    !        
    slen = len_trim(input_filename)
    !
    call create_file_name(fname,input_filename(1:slen),image_count,2)
    !

    write(6,'(3a)') 'File name:',trim(fname),':'

    call find_free_unit_number(funit)

    open(funit,file=fname,form='formatted')

    do ir = 1,yres
       write(funit,'(1024(1x,g12.5))') (image(ic,ir),ic=1,xres)
    end do

    close(funit)


    return 
  end subroutine save_text_image


  subroutine create_file_name(fname,fbase,image_count,flag)
    implicit none
    character*(*) fname,fbase
    integer image_count,flag
    !
    !     This routine writes system dependent file names.
    !
    integer slen
    !
    if (flag.eq.1) then 
       !
       slen = len_trim(fbase)

       if (image_count.lt.10) then 
          write(fname,'(a,i1,a,Z1)') fbase(1:slen)//'_image0',image_count,'.jpg\0'
       else
          write(fname,'(a,i2,a,Z1)') fbase(1:slen)//'_image',image_count,'.jpg\0'
       endif
       !
    elseif (flag.eq.2) then 
       !
       slen = len_trim(fbase)

       if (image_count.lt.10) then 
          write(fname,'(a,i1,a,Z1)') fbase(1:slen)//'_data0',image_count,'.txt'
       else
          write(fname,'(a,i2,a,Z1)') fbase(1:slen)//'_data',image_count,'.txt'
       endif
       !
    elseif (flag.eq.3) then 
       !
       slen = len_trim(fbase)

       if (image_count.lt.10) then 
          write(fname,'(a,i1,a,Z1)') fbase(1:slen)//'_analysis0',image_count,'.dat'
       else
          write(fname,'(a,i2,a,Z1)') fbase(1:slen)//'_analysis',image_count,'.dat'
       endif
       !
    endif

    return 
  end subroutine create_file_name








end module process_irtv
