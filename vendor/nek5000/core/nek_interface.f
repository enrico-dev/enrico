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
        end function nek_get_nelgt

        !> Get value of nelt (number of local elements)
        function nek_get_nelt() result(c_nelt) bind(C)
          integer(C_INT) :: c_nelt
          c_nelt = nelt
        end function nek_get_nelt

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

      !> Depreacted: Do necessary setup at beginning of Picard iteration
      !!
      !! Currently ununsed in STREAM.  Retained from coupling scheme in
      !! MOOSE coupling project
      subroutine nek_init_step()
        include 'SIZE'
        include 'TSTEP'
        include 'INPUT'
        include 'CTIMER'
        real*4 papi_mflops
        integer*8 papi_flops
        integer icall, kstep, i, pstep
        common /cht_coupler/ pstep
        save icall

        if (icall.eq.0) then
        call nekgsync()
        if (instep.eq.0) then
          if(nid.eq.0) write(6,'(/,A,/,A,/)')
     &      ' nsteps=0 -> skip time loop',
     &      ' running solver in post processing mode'
        else
          if(nio.eq.0) write(6,'(/,A,/)') 'Starting time loop ...'
        endif
        isyc  = 0
        itime = 0
        if(ifsync) isyc=1
        itime = 1
        call nek_comm_settings(isyc,itime)
        call nek_comm_startstat()
        istep  = 0
        icall  = 1
        endif

        istep=istep+1

        if (lastep .eq. 1) then
          pstep=2
        else
          call nek_advance
          pstep=2
        endif
      return
      end

      !> Deprecated: Increment pstep counter in Picard iteration
      !!
      !! Currently ununsed in STREAM.  Retained from coupling scheme in
      !! MOOSE coupling project
      subroutine nek_step()
        include 'SIZE'
        include 'TSTEP'
        include 'INPUT'
        include 'CTIMER'
        common /cht_coupler/ pstep
        integer pstep
        pstep=pstep+1
      end subroutine nek_step

      !> Deprecated: Save and dump field data at end of Picard iteration
      !!
      !! Currently ununsed in STREAM.  Retained from coupling scheme in
      !! MOOSE coupling project
      subroutine nek_finalize_step()
        include 'SIZE'
        include 'TSTEP'
        include 'INPUT'
        include 'CTIMER'
        common /cht_coupler/ pstep
        integer pstep

        real*4 papi_mflops
        integer*8 papi_flops
        integer icall, kstep, knstep, i

        if (param(103).gt.0)   call q_filter      (param(103))
        call setup_convect (pstep)  ! Save convective velocity _after_ filter

        call userchk
        call prepost (.false.,'his')
        call in_situ_check()

        if (mod(istep,nsteps).eq.0) lastep=1
        call nek_comm_settings(isyc,0)
        call comment
      end subroutine nek_finalize_step

      !> Deprecated: Expand Fourier coefficients in FET scheme
      !!
      !! Currently ununsed in STREAM.  Retained from FET scheme in
      !! MOOSE coupling project.
      subroutine nek_expansion()
        parameter (nl_max=100)
        parameter (nf_max=100)

        include 'SIZE'
        include 'TOTAL'
        common/expansion_tdata/n_legendre, m_fourier
        common/expansion_tcoef/coeff_tij(nl_max,nf_max)
        integer e,f

        real*8 fmode(lx1,ly1,lz1,lelt), cache(lx1,ly1,lz1,lelt)
        real*8 sint, sint1

        ntot=nx1*ny1*nz1*nelt

        zmin=glmin(zm1,ntot)
        zmax=glmax(zm1,ntot)

        call rzero(fmode,ntot)
        do i0=1,n_legendre
          do j0=1,m_fourier
            call nek_mode(fmode,i0,j0)
            sint1=0.0
            do e=1,nelt
              do f=1,6
                sint=0.0
                if (cbc(f,e,1).eq.'W  ') then
                  call col3(cache,fmode,t,ntot)
                  call surface_int(sint,sarea,cache,e,f)
                  sint1=sint1+sint
                endif
              enddo
            enddo
          call  gop(sint1,wtmp,'+  ',1)
          coeff_tij(i0,j0)=sint1*2.0/(0.5*(zmax-zmin))
          ! Note that R=0.5 here!!!!
          enddo
        enddo
        ! For Testing
        call nek_testp()
      end subroutine nek_expansion

      !> Deprecated: Diagonalization in FET scheme.
      !!
      !! Currently ununsed in STREAM.  Retained from FET scheme in
      !! MOOSE coupling project.
      subroutine nek_diag()
        parameter (nl_max=100)
        parameter (nf_max=100)
        include 'SIZE'
        include 'TOTAL'
        common/expansion_tdata/n_legendre, m_fourier
        common/expansion_tcoef/coeff_tij(nl_max,nf_max)
        common/diags_coeff/diag_c(nl_max,nf_max)
        integer e,f

        real*8 fmode(lx1,ly1,lz1,lelt)
        real*8 cache(lx1,ly1,lz1,lelt)
        real*8 sint, sint1

        ntot=nx1*ny1*nz1*nelt

        zmin=glmin(zm1,ntot)
        zmax=glmax(zm1,ntot)

        call rzero(fmode,ntot)
           do i0=1,n_legendre
             do j0=1,m_fourier
               call nek_mode(fmode,i0,j0)
               sint1=0.0
               do e=1,nelt
                 do f=1,6
                   sint=0.0
                   if (cbc(f,e,1).eq.'W  ') then
                     call col3(cache,fmode,fmode,ntot)
                     call surface_int(sint,sarea,cache,e,f)
                     sint1=sint1+sint
                   endif
                 enddo
               enddo
            call  gop(sint1,wtmp,'+  ',1)
            diag_c(i0,j0)=sint1*2.0/(0.5*(zmax-zmin))
            if (nid.eq.0) write(6,*)i0,j0,diag_c(i0,j0)
          enddo
        enddo
      end subroutine nek_diag

      !> Deprecated: Routine for debugging in FET scheme.
      !!
      !! Currently ununsed in STREAM.  Retained from FET scheme in
      !! MOOSE coupling project.
      subroutine nek_testp()
        parameter (nl_max=100)
        parameter (nf_max=100)
        include 'SIZE'
        include 'TOTAL'
        common/expansion_tdata/n_legendre, m_fourier
        common/expansion_tcoef/coeff_tij(nl_max,nf_max)
        integer e,f
        real*8 fmode(lx1,ly1,lz1,lelt), cache(lx1,ly1,lz1,lelt)
        real*8 fun(lx1,ly1,lz1,lelt)
        ntot=nx1*ny1*nz1*nelt
        call rzero(fmode,ntot)
        call rzero(fun,ntot)
           do i0=1,n_legendre
             do j0=1,m_fourier
               call nek_mode(fmode,i0,j0)
               do i=1,ntot
               fun(i,1,1,1)=fun(i,1,1,1)+fmode(i,1,1,1)*coeff_tij(i0,j0)
               enddo
            enddo
           enddo
        call rzero(cache,ntot)
        call sub3(cache,fun,t,ntot)
        do i=1,ntot
          c1=cache(i,1,1,1)**2
          cache(i,1,1,1)=c1
        enddo

        sint1=0.0
        sarea1=0.0
        do e=1,lelt
          do f=1,6
            call surface_int(sint,sarea,cache,e,f)
            if (cbc(f,e,1).eq.'W  ') then
             sint1=sint1+sint
             sarea1=sarea1+sarea
            endif
           enddo
        enddo
        call  gop(sint1,wtmp,'+  ',1)
        call  gop(sarea1,wtmp,'+  ',1)

        er_avg=sqrt(sint1/sarea1)
        if (nid.eq.0) write(6,*)"Error: ",m_fourier,n_legendre,er_avg
      end subroutine nek_testp

      !> Deprecated: Used in FET scheme
      !!
      !! Currently ununsed in STREAM.  Retained from FET scheme in
      !! MOOSE coupling project.
      subroutine nek_mode(fmode,im,jm)
      include 'SIZE'
      include 'TOTAL'
      integer e,f
      real*8 fmode(lx1,ly1,lz1,lelt)
      ntot=nx1*ny1*nz1*nelt
      zmin=glmin(zm1,ntot)
      zmax=glmax(zm1,ntot)
             do e=1,nelt
               do f=1,6
                 if (cbc(f,e,1).eq.'W  ') then
                   call dsset(nx1,ny1,nz1)
                   iface  = eface1(f)
                   js1    = skpdat(1,iface)
                   jf1    = skpdat(2,iface)
                   jskip1 = skpdat(3,iface)
                   js2    = skpdat(4,iface)
                   jf2    = skpdat(5,iface)
                   jskip2 = skpdat(6,iface)
                   do j2=js2,jf2,jskip2
                     do j1=js1,jf1,jskip1
                       x=xm1(j1,j2,1,e)
                       y=ym1(j1,j2,1,e)
                       z=zm1(j1,j2,1,e)
                       z_leg=2*((z-zmin)/zmax)-1
                       theta=atan2(y,x)
                       fmode(j1,j2,1,e)=
     &                 pl_leg(z_leg,im-1)*fl_four(theta,jm-1)
                     enddo
                   enddo
                 endif
              enddo
            enddo
      return
      end

      !> Deprecated: Reconstruct flux in FET scheme
      !!
      !! Currently ununsed in STREAM.  Retained from FET scheme in
      !! MOOSE coupling project.
      subroutine flux_reconstruction()
        include 'SIZE'
        include 'TOTAL'
        parameter (nl_max=100)
        parameter (nf_max=100)
        common/expansion_tdata/n_legendre, m_fourier
        common/expansion_tcoef/coeff_tij(nl_max,nf_max)
        common/expansion_fcoef/coeff_fij(nl_max,nf_max)
        common/expansion_recon/flux_recon(lx1,ly1,lz1,lelt)
        integer e,f
        real*8 coeff_base
        real*8 fmode(lx1,ly1,lz1,lelt)
        ntot=nx1*ny1*nz1*nelt

        ! --------------------------------------
        ! The flux from MOOSE must have proper sign
        ! --------------------------------------
        ! coeff_base=-1.0
        ! --------------------------------------

        call rzero(flux_recon,ntot)
        do i0=1,n_legendre
          do j0=1,m_fourier
            call rzero(fmode,ntot)
            call nek_mode(fmode,i0,j0)
            do i=1,ntot
              flux_recon(i,1,1,1)= flux_recon(i,1,1,1)
     &                + coeff_base*fmode(i,1,1,1)*coeff_fij(i0,j0)
            enddo
          enddo
        enddo

        ! Below is for testing
        ! do i=1,ntot
        !   flux_recon(i,1,1,1)= 0.0
        ! enddo
        ! return
      end subroutine flux_reconstruction

      !> Deprectaed: Calculates Legendre polynomials Pn(x)
      !!
      !! Calculates Legendre polynomials Pn(x) using the recurrence
      !! relation
      !! if n > 100 the function retuns 0.0
      !!
      !! Currently ununsed in STREAM.  Retained from FET scheme in
      !! MOOSE coupling project.
      function pl_leg(x,n)
      real*8 pl,pl_leg
      real*8 x
      real*8 pln(0:n)
      integer n, k
      pln(0) = 1.0
      pln(1) = x
      if (n.le.1) then
        pl = pln(n)
          else
            do k=1,n-1
              pln(k+1)=
     &  ((2.0*k+1.0)*x*pln(k)-dble(k)*pln(k-1))/(dble(k+1))
            end do
            pl = pln(n)
      end if
      pl_leg=pl*sqrt(dble(2*n+1)/2.0)
      end function pl_leg

      !> Deprecated
      !!
      !! Currently ununsed in STREAM.  Retained from FET scheme in
      !! MOOSE coupling project.
      function fl_four(x,n)
        real*8 fl_four,pi
        real*8 x, A
        integer n, k
        pi=4.0*atan(1.0)
        A=1.0/sqrt(2*pi)
        if (n.gt.0) A=1.0/sqrt(pi)
        fl_four=A*cos(n*x)
c       fl_four=1.0/sqrt(pi)
      end function fl_four

      !> Deprecated
      !!
      !! Currently ununsed in STREAM.  Retained from FET scheme in
      !! MOOSE coupling project.
      subroutine heat_balance(fflux)
        include 'SIZE'
        include 'TOTAL'
        real sint, sint1, sarea, sarea1, wtmp, fflux
        real cache(lx1,ly1,lz1,lelt),dtdz(lx1,ly1,lz1,lelt),
     &   dtdy(lx1,ly1,lz1,lelt),
     &   dtdx(lx1,ly1,lz1,lelt)
        integer e, f

        n1=nelt*lx1*ly1*lz1
        n2=nelt*lx2*ly2*lz2
        nxyz=lx1*ly1*lz1
        call gradm1(dtdx,dtdy,dtdz,t(1,1,1,1,1))

        if (istep.eq.0) then
          sint=0.0
          do e=1,lelt
            do f=1,6
              if (cbc(f,e,1).eq.'O  ') then
               sint=sint+1.0
              endif
            enddo
          enddo
          call  gop(sint,wtmp,'+  ',1)
          if (nid.eq.0) then
            write(6,*)"*** Outlet boundaries: ",sint
          endif
        endif

        sint1=0.0
        sarea1=0.0
        do e=1,lelt
          do i=1,nxyz
            cache(i,1,1,e)=t(i,1,1,e,1)*vz(i,1,1,e)
          enddo
          do f=1,6
            sint=0.0
            call surface_int(sint,sarea,cache,e,f)
            if (cbc(f,e,1).eq.'O  ') then
             sint1=sint1+sint
             sarea1=sarea1+sarea
            endif
           enddo
        enddo
        call  gop(sint1,wtmp,'+  ',1)
        test_flux=sint1

        sint1=0.0
        sarea1=0.0
        do e=1,lelt
          do f=1,6
            sint=0.0
            call surface_int(sint,sarea,dtdz,e,f)
            if (cbc(f,e,1).eq.'v  ') then
             sint1=sint1+sint
             sarea1=sarea1+sarea
            endif
           enddo
        enddo
        call  gop(sint1,wtmp,'+  ',1)
        call  gop(sarea1,wtmp,'+  ',1)
        flux_inlet=sint1

        sint1=0.0
        sarea1=0.0
        do e=1,lelt
          do f=1,6
            sint=0.0
            call surface_int(sint,sarea,vz,e,f)
            if (cbc(f,e,1).eq.'W  ') then
             sarea1=sarea1+sarea
            endif
           enddo
        enddo
        call  gop(sarea1,wtmp,'+  ',1)

        if (nid.eq.0) then
          rr_loss=
     &    flux_inlet*param(8)/(test_flux+flux_inlet*param(8))
          write(6,*)"*** Heat balance: ",
     &    sarea1*fflux,test_flux+flux_inlet*param(8)
          write(6,*)"*** Heat loss: ",rr_loss*100.0," %"
        endif
      end subroutine heat_balance
