module nek_interface
  use, intrinsic :: ISO_C_BINDING
  use nek_geom, only: xm1, ym1, zm1
  use nek_size, only: lelg, lx1
  use nek_parallel, only: lglel, gllel, gllnid

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

!  ROR: 2018-02-06: These aren't part of the current implementations
!
!  !> Retrieves an array of centriods for a given array of local element numbers.
!  !!
!  !! @param[in]  lelts   An array of local element numbers.
!  !! @param[in]  n_lelts The number of local elements in *lelts*.
!  !! @param[out] ctroids An array of centroids corresponding to each local element in *lelts*.
!  function nek_get_lelt_centroids(lelts, n_lelts, ctroids) result(ierr) bind(C)
!    integer(C_INT), dimension(n_lelts), intent(in) :: lelts
!    integer(C_INT), intent(in), value :: n_lelts
!    type(Position), dimension(n_lelts), intent(out) :: ctroids
!    integer(C_INT) :: ierr, i
!
!    do i = 1, n_lelts
!      ! TODO: Does not handle GLL indices correctly!  Just demos interface
!      ctroids(i)%x = xm1(1,1,1,lelts(i))
!      ctroids(i)%y = ym1(1,1,1,lelts(i))
!      ctroids(i)%z = zm1(1,1,1,lelts(i))
!    end do
!
!    ierr = 0
!  end function nek_get_lelt_centroids
!
!  function nek_get_gelt_centroids(gelts, n_gelts, ctroids) result(ierr) bind(C)
!    integer(C_INT), dimension(n_gelts), intent(in) :: gelts
!    integer(C_INT), intent(in), value :: n_gelts
!    type(Position), dimension(n_lelts), intent(out) :: ctroids
!    integer(C_INT) :: ierr, i
!
!    do i = 1, n_gelts
!      ctroids(i)%x = xm1(1,1,1,lelts(i))
!      ctroids(i)%y = ym1(1,1,1,lelts(i))
!      ctroids(i)%z = zm1(1,1,1,lelts(i))
!    end do
!
!  end function nek_get_gelt_centroids

  ! TODO: Only works for 3D
  ! TODO: Does not account for curvature of spectral element.
  function nek_get_global_elem_centroid(global_elem, centroid) result(ierr) bind(C)
    integer(C_INT), intent(in), value :: global_elem
    type(Position), intent(out) :: centroid
    integer(C_INT) :: ierr, i, j, k

    local_elem = gllel(global_elem)

    centroid%x = 0
    centroid%y = 0
    centroid%z = 0
    do k = 1, lz1, lz1-1
      do j = 1, ly1, ly1-1
        do i = 1, lx1, lx1-1
          centroid%x = centroid%x + xm1(i,j,k,local_elem)
          centroid%y = centroid%y + ym1(i,j,k,local_elem)
          centroid%z = centroid%z + zm1(i,j,k,local_elem)
        end do
      end do
    end do

    centroid%x = centroid%x / 8.
    centroid%y = centroid%y / 8.
    centroid%z = centroid%z / 8.

    ierr = 0
  end function nek_get_global_elem_centroid

  function nek_get_lelg() result(c_lelg)
    integer(C_INT) :: c_lelg
    c_lelg = lelg
  end function nek_get_lelg

  function nek_get_lelt() result(c_lelt)
    integer(C_INT) :: c_lelt
    c_lelt = lelt
  end function nek_get_lelt

  function nek_get_lx1() result(c_lx1)
    integer(C_INT) :: c_lx1
    c_lelg = lx1
  end function nek_get_lx1

end module nek_interface