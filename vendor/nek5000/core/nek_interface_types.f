      module nek_interface_types
      use, intrinsic :: ISO_C_BINDING
      implicit none

      !> Describes an (x,y,z) coordinate in 3D space
      !!
      !!
      !! \var x
      !! \var y
      !! \var z
      type, bind(C) :: Position
         real(C_DOUBLE) :: x
         real(C_DOUBLE) :: y
         real(C_DOUBLE) :: z
      end type Position

      end module nek_interface_types
