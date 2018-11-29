      !> Subroutines exposed to C/C++ driver, which use the ISO_C_BINDING module
      module nek_interface
        use, intrinsic :: ISO_C_BINDING
        use :: nek_interface_types
        implicit none

        include 'SIZE'
        include 'MASS'
        include 'GEOM'
        include 'PARALLEL'
        include 'SOLN'
        include 'TSTEP'

      contains

        !> Reset the counters necessary to resume timestepping at the next Picard iteration
        subroutine nek_reset_counters() bind(C)
          time = 0
          lastep = 0
        end subroutine nek_reset_counters

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
        !> \param global_elem A global element ID
        !> \param rank An MPI rank
        !> \return 1 if the global element ID is in the given rank; 0 otherwise
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

        !> Return true if a local element is in the fluid region
        !> \param local_elem  A local element ID
        !> \return 1 if the local element is in fluid; 0 otherwise
        function nek_local_elem_is_in_fluid(local_elem)
     &      result(result) bind(C)
          integer(C_INT), value :: local_elem
          integer(C_INT) :: result, global_elem
          global_elem = lglel(local_elem)
          result = nek_global_elem_is_in_fluid(global_elem)
        end function nek_local_elem_is_in_fluid

        !> Return true if a global element is in the fluid region
        !> \param global_elem  A global element ID
        !> \return 1 if the global element is in fluid; 0 otherwise
        function nek_global_elem_is_in_fluid(global_elem)
     &      result(result) bind(C)
          integer (C_INT), value :: global_elem
          integer (C_INT) :: result
          if (global_elem <= nelgv) then
            result = 1
          else
            result = 0
          end if
        end function nek_global_elem_is_in_fluid

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
          integer :: k

          if (local_elem <= nelt) then
            volume = sum(bm1(1:nx1,1:ny1,1:nz1,local_elem))
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
            temperature = sum(vtrans(1:nx1,1:ny1,1:nz1,local_elem,2)
     &                        * bm1(1:nx1,1:ny1,1:nz1,local_elem)
     &                        * t(1:nx1,1:ny1,1:nz1,local_elem,1))
     &                    / sum(vtrans(1:nx1,1:ny1,1:nz1,local_elem,2)
     &                          * bm1(1:nx1,1:ny1,1:nz1,local_elem))
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

        !> Get value of nelgt (number of global elements)
        function nek_get_nelgt() result(c_nelgt) bind(C)
          integer(C_INT) :: c_nelgt
          c_nelgt = nelgt
        end function nek_get_nelgt

        !> Get value of nelt (number of local elements)
        function nek_get_nelt() result(c_nelt) bind(C)
          integer(C_INT) :: c_nelt
          c_nelt = nelt
        end function nek_get_nelt

        !> Set the heat source term in local heat array
        !!
        !! The units of heat must match on the unit system that was used
        !! to setup the Nek5000 problem. The caller must handle any
        !! necessary conversions.
        !!
        !! \param local_elem A local element ID
        !! \param heat The non-dimensional heat source term
        !! \return Error code
        function nek_set_heat_source(local_elem, heat)
     &      result(ierr) bind(C)
          integer(C_INT), value :: local_elem
          real(C_DOUBLE), value :: heat
          integer(C_INT) :: ierr
          include 'STREAM'
          if (local_elem <= nelt) then
            localq(local_elem) = heat
            ierr = 0
          else
            ierr = 1
          endif
        end function nek_set_heat_source

      end module nek_interface

      !> Get the heat source term for a given gridpoint
      !!
      !! This is called in Nek5000.  It returns the heat source term
      !! from the local heat array, localq
      !!
      !! \param ix x-index of GLL gridpoint
      !! \param iy y-index of GLL gridpoint
      !! \param iz z-index of GLL gridpoint
      !! \param eg A global element ID
      subroutine userq  (ix,iy,iz,eg)
        include 'SIZE'
        include 'TOTAL'
        include 'NEKUSE'
        include 'STREAM'

        integer ix, iy, iz, eg, local_elem

        local_elem = gllel(eg)
        qvol = localq(local_elem)
      end
