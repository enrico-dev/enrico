      !> \file nek_interface.F90
      !> Functions and data types for accessing Nek5000 data
      module nek_interface
        use, intrinsic :: ISO_C_BINDING
        use nek_geom, only: xm1, ym1, zm1
        use nek_size, only: lelg, lelt, lx1, nelt, nx1, ny1, nz1, nid
        use nek_parallel, only: lglel, gllel, gllnid, nelgt
        use nek_mass, only: bm1
        use nek_soln, only: t
      
        implicit none
      
      contains

        !> Get the volume of a local element
        !!
        !! The units of the volume are dimensionless and must be interpreted based on the
        !! setup of the Nek5000
        !!
        !! \param[in] local_elem A local element ID
        !! \param[out] volume The dimensionless volume of the local element
        !! \result Error code
        function nek_get_local_elem_volume(local_elem, volume) result(ierr) bind(C)
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
      
        function nek_get_local_elem_temperature(local_elem, temperature) result(ierr) bind(C)
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
        function nek_get_global_elem(local_elem) result(global_elem) bind(C)
          integer(C_INT), value :: local_elem
          integer(C_INT) :: global_elem
          global_elem = lglel(local_elem)
        end function
      
        !> Get the local element ID for a given global element
        !>
        !> \param[in] global_elem A global element ID
        !> \result The corresponding local element ID
        function nek_get_local_elem(global_elem) result(local_elem) bind(C)
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
