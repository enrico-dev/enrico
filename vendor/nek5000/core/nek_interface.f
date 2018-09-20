      module nek_interface
      use, intrinsic :: ISO_C_BINDING
      use :: nek_interface_types
      implicit none

      include 'SIZE'
      include 'MASS'
      include 'GEOM'
      include 'PARALLEL'

      contains

      !> Get the coordinates of a local element's centroid
      !!
      !! The units of the coordinate are dimensionless and must be interpreted based on the
      !! setup of the Nek5000
      !!
      !! \param[in] local_elem A local element ID
      !! \param[out] centroid The dimensionless coordinates of the local element's centroid
      !! \result Error code
      !! \todo Only works for 3D
      function nek_get_local_elem_centroid(local_elem, centroid)
     &      result(ierr) bind(C)

      integer(C_INT), intent(in), value :: local_elem
         type(Position), intent(out) :: centroid
         integer(C_INT) :: ierr
         integer :: i, j, k
         real(C_DOUBLE) :: mass

         if (local_elem <= nelt) then
            centroid%x = 0.
            centroid%y = 0.
            centroid%z = 0.
            mass = 0.

            do k = 1, nz1
            do j = 1, ny1
            do i = 1, nx1
               centroid%x = centroid%x +
     &                      xm1(i,j,k,local_elem)*bm1(i,j,k,local_elem)
               centroid%y = centroid%y +
     &                      ym1(i,j,k,local_elem)*bm1(i,j,k,local_elem)
               centroid%z = centroid%z +
     &                      zm1(i,j,k,local_elem)*bm1(i,j,k,local_elem)
               mass = mass + bm1(i,j,k,local_elem)
            end do
            end do
            end do

            centroid%x = centroid%x / mass
            centroid%y = centroid%y / mass
            centroid%z = centroid%z / mass

            !print *, local_elem, mass

         ierr = 0
      else
         ierr = 1
      end if
      end function nek_get_local_elem_centroid

      !> Get the coordinates of a global element's centroid
      !!
      !! The units of the coordinate are dimensionless and must be interpreted based on the
      !! setup of the Nek5000
      !!
      !! \param[in] global_elem A global element ID
      !! \param[out] centroid The dimensionless coordinates of the global element's centroid
      !! \result Error code
      !! \todo Only works for 3D
      function nek_get_global_elem_centroid(global_elem, centroid)
     $      result(ierr) bind(C)
         integer(C_INT), intent(in), value :: global_elem
         type(Position), intent(out) :: centroid
         integer(C_INT) :: ierr
         integer :: i, j, k
         real(C_DOUBLE) :: mass
         integer(C_INT) :: local_elem

         if (nek_global_elem_is_in_rank(global_elem, nid) == 0) then
            local_elem = gllel(global_elem)
            ierr = nek_get_local_elem_centroid(local_elem, centroid)
         else
            ierr = 1
      end if
      end function nek_get_global_elem_centroid

      !> Return true if a global element is in a given MPI rank
      !> \param A global element ID
      !> \param An MPI rank
      !> \return True if the global element ID is in the given rank
      function nek_global_elem_is_in_rank(global_elem, rank)
     $      result(result) bind(C)
         integer (C_INT), value :: global_elem, rank
         integer(C_INT) :: result
         if (rank == gllnid(global_elem)) then
           result = 1
         else
           result = 0
         end if
      end function nek_global_elem_is_in_rank

      end module nek_interface
