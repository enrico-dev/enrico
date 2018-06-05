module fobj
  use ISO_C_BINDING
  implicit none

  ! Specifying a `name` in the `bind` atribute is not allowed.  The Fortran 2003 standard says:
  !
  !   "The names of the types and the names of the components are not significant for the purposes 
  !   of determining whether a Fortran derived type is interoperable with a C struct type."
  !
  ! A working example from the standard shows a C struct and Fortran type with differeng names:
  ! 
  !   typedef struct {
  !     int m, n;
  !     float r;
  !   } myctype;
  !   
  !  USE, INTRINSIC :: ISO_C_BINDING
  !  TYPE, BIND(C) :: MYFTYPE
  !    INTEGER(C_INT) :: I, J
  !    REAL(C_FLOAT) :: S
  !  END TYPE MYFTYPE:: Position

  type, bind(C) :: Position
    real(C_DOUBLE) :: x
    real(C_DOUBLE) :: y
    real(C_DOUBLE) :: z
  end type Position



  contains

    subroutine set_position_array(pstns, n) bind(C)
      type(Position), intent(inout), dimension(1:n) :: pstns
      integer(C_INT), intent(in):: n
      integer(C_INT) :: i
      do i = 1, n
        pstns(i)%x = i*10 + 1
        pstns(i)%y = i*10 + 2
        pstns(i)%z = i*10 + 3
      enddo
    end subroutine set_position_array
end module fobj

