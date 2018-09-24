      module nek_interface
        use, intrinsic :: ISO_C_BINDING
        use :: nek_interface_types
        implicit none

        include 'SIZE'
        include 'MASS'
        include 'GEOM'
        include 'PARALLEL'
        include 'SOLN'

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

        !> Get the volume of a local element
        !!
        !! The units of the volume are dimensionless and must be interpreted based on the
        !! setup of the Nek5000
        !!
        !! \param[in] local_elem A local element ID
        !! \param[out] volume The dimensionless volume of the local element
        !! \result Error code
        function nek_get_local_elem_volume(local_elem, volume) 
     &      result(ierr) bind(C)
          integer(C_INT), intent(in), value :: local_elem
          real(C_DOUBLE), intent(out) :: volume
          integer(C_INT) :: ierr

          if (local_elem <= nelt) then
            volume = sum(bm1(1:nx1, 1:ny1, 1:nz1, local_elem))
            ierr = 0
          else
            ierr = 1
          end if
        end function nek_get_local_elem_volume

        function nek_get_local_elem_temperature(local_elem, temperature) 
     &      result(ierr) bind(C)
          integer(C_INT), intent(in), value :: local_elem
          real(C_DOUBLE), intent(out) :: temperature
          integer(C_INT) :: ierr

          if (local_elem <= nelt) then
            temperature = sum(t(1:nx1, 1:ny1, 1:nz1, local_elem, 1))
            ierr = 0
          else
            ierr = 1
          end if
        end function nek_get_local_elem_temperature

        !> Get the global element ID for a given local element
        !>
        !> \param[in] local_elem A local element ID
        !> \result The corresponding global element ID
        function nek_get_global_elem(local_elem) 
     &      result(global_elem) bind(C)
          integer(C_INT), value :: local_elem
          integer(C_INT) :: global_elem
          global_elem = lglel(local_elem)
        end function

        !> Get the local element ID for a given global element
        !>
        !> \param[in] global_elem A global element ID
        !> \result The corresponding local element ID
        function nek_get_local_elem(global_elem) 
     &      result(local_elem) bind(C)
          integer(C_INT), value :: global_elem
          integer(C_INT) :: local_elem
          local_elem = gllel(global_elem)
        end function

        !> Get value of lelg (max number of global elements)
        function nek_get_lelg() result(c_lelg) bind(C)
          integer(C_INT) :: c_lelg
          c_lelg = lelg
        end function nek_get_lelg

        !> Get value of lelt (max number of local elements)
        function nek_get_lelt() result(c_lelt) bind(C)
          integer(C_INT) :: c_lelt
          c_lelt = lelt
        end function nek_get_lelt

        !> Get value of lx1 (number of GLL gridpoints in x-dimension)
        function nek_get_lx1() result(c_lx1) bind(C)
          integer(C_INT) :: c_lx1
          c_lx1 = lx1
        end function nek_get_lx1

        !> Get value of nelgt
        function nek_get_nelgt() result(c_nelgt) bind(C)
          integer(C_INT) :: c_nelgt
          c_nelgt = nelgt
        end function

        !> Get value of nelt (number of local elements)
        function nek_get_nelt() result(c_nelt) bind(C)
          integer(C_INT) :: c_nelt
          c_nelt = nelt
        end function

      end module nek_interface
