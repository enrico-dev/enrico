module nek_interface
  use, intrinsic :: ISO_C_BINDING
  use nek_geom, only: xm1, ym1, zm1
  use nek_size, only: nelt, nx1, ny1, nz1
  use nek_parallel, only: lglel, gllel, gllnid, nelgt
  use nek_mass, only: bm1
  use nek_size, only: lelg, lelt, lx1

  implicit none

  !> Describes an (x,y,z) coordinate in 3D space
  !!
  !! @var x
  !! @var y
  !! @var z
  !! @todo Move this to its own header?
  type, bind(C) :: Position
    real(C_DOUBLE) :: x
    real(C_DOUBLE) :: y
    real(C_DOUBLE) :: z
  end type Position

contains

  ! TODO: Only works for 3D
  function nek_get_global_elem_centroid(global_elem, centroid) result(ierr) bind(C)
    integer(C_INT), intent(in), value :: global_elem
    type(Position), intent(out) :: centroid
    integer(C_INT) :: ierr
    integer :: i, j, k
    real(C_DOUBLE) :: mass
    integer(C_INT) :: local_elem

    local_elem = gllel(global_elem)

    centroid%x = 0.
    centroid%y = 0.
    centroid%z = 0.
    mass = 0.

    do k = 1, nz1
      do j = 1, ny1
        do i = 1, nx1
          centroid%x = centroid%x + xm1(i,j,k,local_elem)
          centroid%y = centroid%y + ym1(i,j,k,local_elem)
          centroid%z = centroid%z + zm1(i,j,k,local_elem)
          mass = mass + bm1(i,j,k,local_elem)
        end do
      end do
    end do

    centroid%x = centroid%x / mass
    centroid%y = centroid%y / mass
    centroid%z = centroid%z / mass

    ierr = 0
  end function nek_get_global_elem_centroid

  function nek_get_lelg() result(c_lelg) bind(C)
    integer(C_INT) :: c_lelg
    c_lelg = lelg
  end function nek_get_lelg

  function nek_get_lelt() result(c_lelt) bind(C)
    integer(C_INT) :: c_lelt
    c_lelt = lelt
  end function nek_get_lelt

  function nek_get_lx1() result(c_lx1) bind(C)
    integer(C_INT) :: c_lx1
    c_lx1 = lx1
  end function nek_get_lx1

  ! Number of elements in mesh
  function nek_get_nelgt() result(c_nelgt) bind(C)
    integer(C_INT) :: c_nelgt
    c_nelgt = nelgt
  end function

  ! Number of elements on process
  function nek_get_nelt() result(c_nelt) bind(C)
    integer(C_INT) :: c_nelt
    c_nelt = nelt
  end function

end module nek_interface
