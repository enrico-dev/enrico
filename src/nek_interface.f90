module nek_interface
  use, intrinsic :: ISO_C_BINDING
  use geom, only: xm1, ym1, zm1

  implicit none

  ! TODO: Move this to its own header?
  type, bind(C) :: Position
    real(C_DOUBLE) :: x
    real(C_DOUBLE) :: y
    real(C_DOUBLE) :: z
  end type Position

contains

  function nek_get_lelt_centroids(lelts, n_lelts, ctroids) result(err) bind(C)
    integer(C_INT), dimension(n_lelts), intent(in) :: lelts
    integer(C_INT), intent(in), value :: n_lelts
    type(Position), dimension(n_lelts), intent(out) :: ctroids
    integer(C_INT) :: err, i

    do i = 1, n_lelts
      ! TODO: Does not handle GLL indices correctly!  Just demos interface
      ctroids(i)%x = xm1(1,1,1,lelts(i))
      ctroids(i)%y = ym1(1,1,1,lelts(i))
      ctroids(i)%z = zm1(1,1,1,lelts(i))
    end do

    err = 0
  end function nek_get_lelt_centroids

end module nek_interface