module nek_interface
  use, intrinsic :: ISO_C_BINDING
  implicit none

contains

  function nek_get_local_el_centroid(local_el, x, y, z) result(err) bind(C)
    integer(C_INT), intent(in) :: local_el
    real(C_DOUBLE), intent(out) :: x, y, z
    integer(C_INT) :: err

    include 'PARALLEL'
    include 'SIZE'

    if (local_el <= lelt) then
      ! TODO: Does not handle GLL indices correctly!  Just a stub!
      x = xm1(1,1,1,local_el)
      y = ym1(1,1,1,local_el)
      z = zm1(1,1,1,local_el)
      err = 0
    else
      err = -1
    end if
  end function nek_get_local_el_centroid

end module nek_interface