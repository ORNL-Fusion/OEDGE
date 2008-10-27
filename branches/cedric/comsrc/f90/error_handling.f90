module error_handling

  interface errmsg

     module procedure rerrmsg,r8errmsg,ierrmsg,cerrmsg,crerrmsg

  end interface


  interface dbgmsg

     module procedure rdbgmsg,r8dbgmsg,idbgmsg,cdbgmsg,crdbgmsg

  end interface

  integer,private :: err1=0,err2=6,err3=-1
  integer,private :: dbg1=6,dbg2=-1,dbg3=-1

contains

  !
  ! Set output channels - allows default output units to be changed
  !

  subroutine set_errmsg_units(u1,u2,u3)
    implicit none
    integer u1,u2,u3

    err1 = u1
    err2 = u2
    err3 = u3

  end subroutine set_errmsg_units

  subroutine set_dbgmsg_units(u1,u2,u3)
    implicit none
    integer u1,u2,u3

    dbg1 = u1
    dbg2 = u2
    dbg3 = u3

  end subroutine set_dbgmsg_units

  !
  ! Error message handling routines
  !

  subroutine rerrmsg(msg,a,unit)
    implicit none
    character*(*) msg
    real a
    integer,optional :: unit

    integer len1,len2

    len1 = len_trim(msg)

    if (present(unit)) then 
       write(unit,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a
    else

       if (err1.ge.0) write(err1,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a
       if (err2.ge.0) write(err2,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a
       if (err3.ge.0) write(err3,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a

    endif

  end subroutine rerrmsg

  subroutine r8errmsg(msg,a,unit)
    implicit none
    character*(*) msg
    real*8 a
    integer,optional :: unit

    integer len1,len2

    len1 = len_trim(msg)

    if (present(unit)) then 
       write(unit,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a
    else

       if (err1.ge.0) write(err1,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a
       if (err2.ge.0) write(err2,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a
       if (err3.ge.0) write(err3,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a

    endif

  end subroutine r8errmsg

  subroutine ierrmsg(msg,a,unit)
    implicit none
    character*(*) msg
    integer a
    integer,optional :: unit

    integer len1,len2

    len1 = len_trim(msg)

    if (present(unit)) then 
       write(unit,'(a,1x,a,1x,a,1x,i10)') 'ERROR:',msg(1:len1),'VALUE =',a
    else
       if (err1.ge.0) write(err1,'(a,1x,a,1x,a,1x,i10)') 'ERROR:',msg(1:len1),'VALUE =',a
       if (err2.ge.0) write(err2,'(a,1x,a,1x,a,1x,i10)') 'ERROR:',msg(1:len1),'VALUE =',a
       if (err3.ge.0) write(err3,'(a,1x,a,1x,a,1x,i10)') 'ERROR:',msg(1:len1),'VALUE =',a
    endif

  end subroutine ierrmsg

  subroutine cerrmsg(msg,a,unit)
    implicit none
    character*(*) msg,a
    integer,optional :: unit

    integer len1,len2

    len1 = len_trim(msg)
    len2 = len_trim(a)

    if (present(unit)) then 
       write(unit,'(a,1x,a,1x,a,1x,a)') 'ERROR:',msg(1:len1),'VALUE =',a(1:len2)
    else
       if (err1.ge.0) write(err1,'(a,1x,a,1x,a,1x,a)') 'ERROR:',msg(1:len1),'VALUE =',a(1:len2)
       if (err2.ge.0) write(err2,'(a,1x,a,1x,a,1x,a)') 'ERROR:',msg(1:len1),'VALUE =',a(1:len2)
       if (err3.ge.0) write(err3,'(a,1x,a,1x,a,1x,a)') 'ERROR:',msg(1:len1),'VALUE =',a(1:len2)
    endif

  end subroutine cerrmsg

  subroutine crerrmsg(msg,a,r,unit)
    implicit none
    character*(*) msg,a
    real :: r
    integer,optional :: unit

    integer len1,len2

    len1 = len_trim(msg)
    len2 = len_trim(a)

    if (present(unit)) then 
       write(unit,'(a,1x,a,1x,a,1x,a,1x,1p,g18.8)') 'ERROR:',msg(1:len1),'MESSAGE =',a(1:len2),r
    else
       if (err1.ge.0) write(err1,'(a,1x,a,1x,a,1x,a,1x,1p,g18.8)') 'ERROR:',msg(1:len1),'MESSAGE =',a(1:len2),r
       if (err2.ge.0) write(err2,'(a,1x,a,1x,a,1x,a,1x,1p,g18.8)') 'ERROR:',msg(1:len1),'MESSAGE =',a(1:len2),r
       if (err3.ge.0) write(err3,'(a,1x,a,1x,a,1x,a,1x,1p,g18.8)') 'ERROR:',msg(1:len1),'MESSAGE =',a(1:len2),r
    endif

  end subroutine crerrmsg


  !
  ! Debug message handling routines
  !


  subroutine rdbgmsg(msg,a,unit)
    implicit none
    character*(*) msg
    real a
    integer,optional :: unit

    integer len1,len2

    len1 = len_trim(msg)

    if (present(unit)) then 
       write(unit,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a
    else
       if (dbg1.ge.0) write(dbg1,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a
       if (dbg2.ge.0) write(dbg2,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a
       if (dbg3.ge.0) write(dbg3,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a
    endif

  end subroutine rdbgmsg

  subroutine r8dbgmsg(msg,a,unit)
    implicit none
    character*(*) msg
    real*8 a
    integer,optional :: unit

    integer len1,len2

    len1 = len_trim(msg)

    if (present(unit)) then 
       write(unit,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a
    else
       if (dbg1.ge.0) write(dbg1,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a
       if (dbg2.ge.0) write(dbg2,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a
       if (dbg3.ge.0) write(dbg3,'(a,1x,a,1x,a,1x,f18.8)') 'ERROR:',msg(1:len1),'VALUE =',a
    endif

  end subroutine r8dbgmsg

  subroutine idbgmsg(msg,a,unit)
    implicit none
    character*(*) msg
    integer a
    integer,optional :: unit

    integer len1,len2

    len1 = len_trim(msg)

    if (present(unit)) then 
       write(unit,'(a,1x,a,1x,a,1x,i10)') 'ERROR:',msg(1:len1),'VALUE =',a
    else
       if (dbg1.ge.0) write(dbg1,'(a,1x,a,1x,a,1x,i10)') 'ERROR:',msg(1:len1),'VALUE =',a
       if (dbg2.ge.0) write(dbg2,'(a,1x,a,1x,a,1x,i10)') 'ERROR:',msg(1:len1),'VALUE =',a
       if (dbg3.ge.0) write(dbg3,'(a,1x,a,1x,a,1x,i10)') 'ERROR:',msg(1:len1),'VALUE =',a
    endif

  end subroutine idbgmsg

  subroutine cdbgmsg(msg,a,unit)
    implicit none
    character*(*) msg,a
    integer,optional :: unit

    integer len1,len2

    len1 = len_trim(msg)
    len2 = len_trim(a)

    if (present(unit)) then 

       write(unit,'(a,1x,a,1x,a,1x,a)') 'ERROR:',msg(1:len1),'VALUE =',a(1:len2)
    else
       if (dbg1.ge.0) write(dbg1,'(a,1x,a,1x,a,1x,a)') 'ERROR:',msg(1:len1),'VALUE =',a(1:len2)
       if (dbg2.ge.0) write(dbg2,'(a,1x,a,1x,a,1x,a)') 'ERROR:',msg(1:len1),'VALUE =',a(1:len2)
       if (dbg3.ge.0) write(dbg3,'(a,1x,a,1x,a,1x,a)') 'ERROR:',msg(1:len1),'VALUE =',a(1:len2)
    endif

  end subroutine cdbgmsg

  subroutine crdbgmsg(msg,a,r,unit)
    implicit none
    character*(*) msg,a
    real :: r
    integer,optional :: unit

    integer len1,len2

    len1 = len_trim(msg)
    len2 = len_trim(a)

    if (present(unit)) then 
       write(unit,'(a,1x,a,1x,a,1x,a,1x,1p,g18.8)') 'ERROR:',msg(1:len1),'MESSAGE =',a(1:len2),r
    else
       if (dbg1.ge.0) write(err1,'(a,1x,a,1x,a,1x,a,1x,1p,g18.8)') 'ERROR:',msg(1:len1),'MESSAGE =',a(1:len2),r
       if (dbg2.ge.0) write(err2,'(a,1x,a,1x,a,1x,a,1x,1p,g18.8)') 'ERROR:',msg(1:len1),'MESSAGE =',a(1:len2),r
       if (dbg3.ge.0) write(err3,'(a,1x,a,1x,a,1x,a,1x,1p,g18.8)') 'ERROR:',msg(1:len1),'MESSAGE =',a(1:len2),r
    endif

  end subroutine crdbgmsg


end module error_handling
