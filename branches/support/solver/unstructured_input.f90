  subroutine ReadUnstructuredInput(buffer)
    use mod_sol22_input
    implicit none
    character*(*) :: buffer
    character*3 :: tag
    character*128 :: line
    integer :: ierr

    ierr= 0

    WRITE(line,'(A128)') buffer

      WRITE(TAG,'(A3)') buffer(3:5)

      call sol22_unstructured_input(tag,line,ierr)


  end subroutine ReadUnstructuredInput

