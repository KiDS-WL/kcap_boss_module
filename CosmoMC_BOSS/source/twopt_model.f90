!pqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpq
!                             modules                                 !
!bdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbd

module compute_twopt
  use settings, only : dp => mcp, feedback
  use numrec
  use quadpack
  use NLRSD_settings
  use NLRSD_power
  use MPIUtils
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN
  implicit none
  real(kind=dp)                             :: smin, smax, kmin, kmax 
  integer,parameter, private                :: numr = 130, numk = 150  !2016.7
  integer, private                          :: save_ell
  real(kind=dp), private                    :: save_s, save_k
  real(kind=dp), private,dimension(numr)    :: r, xi_0, xi_2, xi_4
  real(kind=dp), private,dimension(numr)    :: xi_0d, xi_2d, xi_4d
  real(kind=dp), private,dimension(numk)    :: kvec, pk_0, pk_2, pk_4
  real(kind=dp), private,dimension(numk)    :: pk_0d, pk_2d, pk_4d
  real(kind=dp), private                    :: save_r
  real(kind=dp), parameter, private         :: linf   = -6.2_dp
  real(kind=dp), parameter, private         :: lsup   = 0._dp
  real(kind=dp), parameter, private         :: inf    = exp(linf)
  real(kind=dp), parameter, private         :: sup    = exp(lsup)
  real(kind=dp), parameter, private         :: kmin_xiell = 1.1e-4_dp
  real(kind=dp), parameter, private         :: kmax_xiell = 1.5_dp
  real(kind=dp), parameter, private         :: epsabs = 1.e-2_dp 
  real(kind=dp), parameter, private         :: epsrel = 1.e-2_dp 
  real(kind=dp), parameter, private         :: epsabs_2 = 1.e-2_dp 
  real(kind=dp), parameter, private         :: epsrel_2 = 1.e-2_dp 
  real(kind=dp), parameter, private         :: epsabs_4 = 2.e-1_dp 
  real(kind=dp), parameter, private         :: epsrel_4 = 1.e-1_dp 
  real(kind=dp), parameter, private         :: kcut_0 = 0.7_dp
  real(kind=dp), parameter, private         :: npow_0 = 2._dp
  real(kind=dp), parameter, private         :: kcut_2 = 0.58_dp
  real(kind=dp), parameter, private         :: npow_2 = 4._dp
  real(kind=dp), parameter, private         :: kcut_4 = 0.6_dp
  real(kind=dp), parameter, private         :: npow_4 = 2._dp

  contains
  
    function cutoff(kh, kcut, npow)
      real(kind=dp)   :: cutoff, kh, kcut, npow
      cutoff = exp(-(kh/kcut)**npow)
    end function cutoff

    function kernel_xi0_sin(kh)
      implicit none
      real(kind=dp)           :: kh, kernel_xi0_sin, P0
      call splint(kvec, pk_0, pk_0d, numk, kh, P0)
      kernel_xi0_sin = P0*kh*cutoff(kh, kcut_0, npow_0)
    end function kernel_xi0_sin

    function kernel_xi0_jell(kh)
      real(kind=dp)   :: kernel_xi0_jell, kr, lk, kh, ci, si
      real(kind=dp)   :: p0, sj, sy, sjp, syp
      kr = kh*save_r
      call sphbes(0, kr, sj, sy, sjp, syp)
      call splint(kvec, pk_0, pk_0d, numk, kh, P0)
      kernel_xi0_jell = p0*sj*kh**2*cutoff(kh, kcut_0, npow_0)
    end function kernel_xi0_jell

    !configuration-space monopole
    function ximonopole(r)
      implicit none
      integer          :: neval, ierr, integr
      real(kind=dp)    :: ximonopole, r, ans, abserr
      integr = 2
      save_r = r
      call qawf(kernel_xi0_sin,inf,r,integr,epsabs,ans,abserr,neval,ierr)
      if(ierr == 0)then
        ximonopole = ans/r/twopisqr
      else
        !call warning('xi0_sin',ierr)  
        !write(*,*)'using alternative integration routine..'
        call qags(kernel_xi0_jell, inf, sup, 1.d-6, 1.d-6, ans, abserr, neval, ierr)
        ximonopole = ans/twopisqr
        !call warning('xi0_jell',ierr)  
      end if
    end function ximonopole
 
    !kernel for the first integral to obtain xi2_nl
    function kernel_xi2_sin(kh)
      implicit none
      real(kind=dp)           :: kh, kr, kernel_xi2_sin, P2
      kr = kh*save_r
      call splint(kvec, pk_2, pk_2d, numk, kh, P2)
      kernel_xi2_sin = P2*(3._dp/kr - kr)*cutoff(kh, kcut_2, npow_2)
    end function kernel_xi2_sin

    !kernel for the second integral to obtain xi2_nl
    function kernel_xi2_cos(kh)
      implicit none
      real(kind=dp)           :: kh, kr, kernel_xi2_cos, P2
      kr = kh*save_r
      call splint(kvec, pk_2, pk_2d, numk, kh, P2)
      kernel_xi2_cos = P2*cutoff(kh, kcut_2, npow_2)
    end function kernel_xi2_cos

    function kernel_xi2_jell(kh)
      real(kind=dp)   :: kernel_xi2_jell, kr, lk, kh, ci, si
      real(kind=dp)   :: p2, sj, sy, sjp, syp
      kr = kh*save_r
      call sphbes(2, kr, sj, sy, sjp, syp)
      call splint(kvec, pk_2, pk_2d, numk, kh, P2)
      kernel_xi2_jell = p2*sj*kh**2*cutoff(kh, kcut_2, npow_2)
    end function kernel_xi2_jell

    !configuration-space quadrupole
    function xiquadrupole(r)
      implicit none
      integer          :: neval, ierr, integr, error
      real(kind=dp)    :: term1, term2, xiquadrupole, r, ans, abserr
      save_r = r
      integr = 2
      error = 0
      call qawf(kernel_xi2_sin,inf,r,integr,epsabs_2,ans,abserr,neval,ierr)
      if(ierr /= 0)  then
        !write(*,*)r
        error = ierr
        !call warning('xi2_sin',ierr)  
      end if
      term1 = ans
      integr = 1
      call qawf(kernel_xi2_cos,inf,r,integr,epsabs_2,ans,abserr,neval,ierr)
      if(ierr /= 0)then
        error = ierr
        !call warning('xi2_cos',ierr)  
      end if
      term2 = -3._dp*ans
      if(error == 0)then
        xiquadrupole = -(term1 + term2)/twopisqr/r**2
      else
        !write(*,*)'r=',r
        call qags(kernel_xi2_jell, inf, sup, 1.d-7, 1.d-7, ans, abserr, neval, ierr)
        !if(ierr /= 0) then
        !  call warning('xi2_jell',ierr)  
        !end if
        xiquadrupole = -ans/twopisqr
      end if
    end function xiquadrupole

    function kernel_xi4_jell(kh)
      real(kind=dp)   :: kernel_xi4_jell, kr, lk, kh, ci, si
      real(kind=dp)   :: p4, sj, sy, sjp, syp
      kr = kh*save_r
      call sphbes(4, kr, sj, sy, sjp, syp)
      call splint(kvec, pk_4, pk_4d, numk, kh, P4)
      kernel_xi4_jell = p4*sj*kh**2*cutoff(kh, kcut_4, npow_4)
    end function kernel_xi4_jell

    !kernel for the first integral to obtain xi4_nl
    function kernel_xi4_sin(kh)
      implicit none
      real(kind=dp)           :: kh, kr, kernel_xi4_sin, P4
      kr = kh*save_r
      call splint(kvec, pk_4, pk_4d, numk, kh, P4)
      kernel_xi4_sin = P4*(105._dp/kr**3 - 45._dp/kr + kr)*cutoff(kh, kcut_4, npow_4)
    end function kernel_xi4_sin

    !kernel for the second integral to obtain xi4_nl
    function kernel_xi4_cos(kh)
      implicit none
      real(kind=dp)           :: kh, kr, kernel_xi4_cos, P4
      kr = kh*save_r
      call splint(kvec, pk_4, pk_4d, numk, kh, P4)
      kernel_xi4_cos = P4*(10._dp - 105._dp/kr**2)*cutoff(kh, kcut_4, npow_4)
    end function kernel_xi4_cos

    !configuration space hexadecapole
    function xihexadecapole(r)
      implicit none
      integer          :: neval, ierr, integr, error
      real(kind=dp)    :: term1, term2, xihexadecapole, r, ans, abserr
      save_r = r
      integr = 2
      error = 0
      call qawf(kernel_xi4_sin,inf,r,integr,epsabs_4,ans,abserr,neval,ierr)
      if(ierr /= 0)then
        error = ierr
        !call warning('xi4_sin',ierr)  
      end if
      term1 = ans
      integr = 1
      call qawf(kernel_xi4_cos,inf,r,integr,epsabs_4,ans,abserr,neval,ierr)
      if(ierr /= 0)then
        error = ierr
        !call warning('xi4_cos',ierr)  
      end if
      term2 = ans
      if(error == 0)then
        xihexadecapole = (term1 + term2)/twopisqr/r**2 
      else
        call qags(kernel_xi4_jell, inf, sup, 1.d-7, 1.d-7, ans, abserr, neval, ierr)
        if(ierr /= 0)  call warning('xi4_jell',ierr)  
        xihexadecapole = ans/twopisqr
      end if
    end function xihexadecapole

    subroutine init_rvec
      implicit none
      integer        :: i
      real(kind=dp)  :: ds
      ds = (smax - smin)/real(numr-1,kind=dp)
      do i = 1,numr
        r(i) = smin + ds*real(i-1,kind=dp)
      end do
    end subroutine init_rvec

    subroutine init_xi_multipoles
      implicit none
      integer                          :: i, erro
      logical                          :: found_nan
      !file for debugging
      call init_rvec
      !setup tables with Fourier-space multipoles for interpolation.
      kmin = kmin_xiell
      kmax = kmax_xiell
      call init_pk_multipoles() 
      found_nan = .false.
      do i = 1,numr
        xi_0(i) = ximonopole(r(i))
        xi_2(i) = xiquadrupole(r(i))
        xi_4(i) = xihexadecapole(r(i))

        found_nan = found_nan .or. (xi_0(i)/=xi_0(i))
        found_nan = found_nan .or. (xi_2(i)/=xi_2(i))
        found_nan = found_nan .or. (xi_4(i)/=xi_4(i))
      end do
      !close(98)
      if(found_nan)then
        write(*,*)'found NaN in multipoles'
        call test_found_nan
      end if
      call spline(r, xi_0, numr, 3d30, 3d30, xi_0d)
      call spline(r, xi_2, numr, 3d30, 3d30, xi_2d)
      call spline(r, xi_4, numr, 3d30, 3d30, xi_4d)
    end subroutine init_xi_multipoles

    function xi2d_at(s,mu)
      implicit none
      real(kind=dp) :: xi2d_at, s, mu
      real(kind=dp) :: ans0, ans2, ans4
      call splint(r, xi_0, xi_0d, numr, s, ans0)
      call splint(r, xi_2, xi_2d, numr, s, ans2)
      call splint(r, xi_4, xi_4d, numr, s, ans4)
      xi2d_at = ans0 + (1.5_dp*mu**2-0.5_dp)*ans2  &
                + 0.125_dp*(35._dp*mu**4 - 30._dp*mu**2 + 3._dp)*ans4
    end function xi2d_at

    function arg_xiell(muf)
      implicit none
      real(kind=dp) :: arg_xiell, muf
      real(kind=dp) :: mu, s, muf2, fac
      real(kind=dp) :: pl 
      muf2 = muf**2
      fac  = sqrt(alpha_lo**2*muf2 + alpha_tr**2*(1._dp - muf2))
      s    = save_s*fac
      mu   = alpha_lo*muf/fac
      pl = plgndr(save_ell, 0, muf)
      arg_xiell = xi2d_at(s,mu)*pl
    end function arg_xiell

    function xi_multipole(s, ell)
      real(kind=dp)  :: xi_multipole, s, ans
      integer        :: ell, neval, ierr
      real(kind=dp)  :: abserr
      save_s   = s
      save_ell = ell
      call qags(arg_xiell, 0._dp, 1._dp, epsabs, epsrel, ans, abserr, neval, ierr)
      if(ierr /= 0)  call warning('xi_multipole',ierr)  
      xi_multipole = ans*(2._dp*real(ell,kind=dp) + 1._dp)
    end function xi_multipole

    function arg_wed(muf)
      implicit none
      real(kind=dp) :: arg_wed,muf
      real(kind=dp) :: mu,s,muf2,fac
      muf2 = muf**2
      fac  = sqrt(alpha_lo**2*muf2+ alpha_tr**2*(1._dp-muf2))
      s    = save_s*fac
      mu   = alpha_lo*muf/fac 
      arg_wed = xi2d_at(s,mu)
    end function arg_wed

    function xi_wedge(s, iwedge, nwedges)
      real(kind=dp)       :: xi_wedge, s
      integer             :: iwedge, nwedges
      integer             :: ierr, neval
      real(kind=dp)       :: mu_min, mu_max, delta_mu
      real(kind=dp)       :: ans, abserr
      delta_mu = 1._dp/real(nwedges, kind=dp) !divide range from 0 to 1 into nwedges bins
      mu_min = real(iwedge-1,kind=dp)*delta_mu  
      mu_max = mu_min + delta_mu
      save_s = s
      call qags(arg_wed, mu_min, mu_max, epsabs, epsrel, ans, abserr, neval, ierr)
      if(ierr /= 0)  call warning('xi_wedges',ierr)  
      xi_wedge = ans/delta_mu
    end function xi_wedge

    subroutine init_pk_kvec
      implicit none
      integer       :: i
      real(kind=dp) :: dk, qmin, qmax
      real(kind=dp) :: kcenter, pow
      qmax = log10(kmax)
      qmin = log10(kmin)
      dk = (qmax - qmin)/real(numk-1,kind=dp)
      do i = 1,numk
        kvec(i) = qmin + dk*real(i-1,kind=dp)
        kvec(i) = 10._dp**kvec(i)
      end do

    end subroutine init_pk_kvec

    subroutine init_pk_multipoles
      implicit none
      integer                          :: i
      !file for debugging
      call init_pk_kvec
      do i = 1, numk
        pk_0(i) = gRPTpower_gal_monopole_at(kvec(i))
        pk_2(i) = gRPTpower_gal_quadrupole_at(kvec(i))
        pk_4(i) = gRPTpower_gal_hexadecapole_at(kvec(i))
      end do
      call spline(kvec, pk_0, numk, 3d30, 3d30, pk_0d)
      call spline(kvec, pk_2, numk, 3d30, 3d30, pk_2d)
      call spline(kvec, pk_4, numk, 3d30, 3d30, pk_4d)
    end subroutine init_pk_multipoles

    function pk2d_at(k, mu)
      implicit none
      real(kind=dp) :: pk2d_at, k, mu
      real(kind=dp) :: ans0, ans2, ans4
      call splint(kvec, pk_0, pk_0d, numk, k, ans0)
      call splint(kvec, pk_2, pk_2d, numk, k, ans2)
      call splint(kvec, pk_4, pk_4d, numk, k, ans4)
      pk2d_at = ans0 + (1.5_dp*mu**2 - 0.5_dp)*ans2  &
                + 0.125_dp*(35._dp*mu**4 - 30._dp*mu**2 + 3._dp)*ans4
    end function pk2d_at

    function arg_pkell(muf)
      implicit none
      real(kind=dp) :: arg_pkell, muf
      real(kind=dp) :: mu, k, muf2, fac
      real(kind=dp) :: pl
      muf2 = muf**2
      fac  = sqrt(muf2/alpha_lo**2 + (1._dp - muf2)/alpha_tr**2)
      k    = save_k*fac
      mu   = muf/alpha_lo/fac
      !the function plgndr gives the Legendre polynomial 
      pl = plgndr(save_ell, 0, muf)
      arg_pkell = pk2d_at(k,mu)*pl
    end function arg_pkell

    function pk_multipole(k, ell)
      real(kind=dp)  :: pk_multipole, k, ans
      integer        :: ell, neval, ierr
      real(kind=dp)  :: abserr
      save_k   = k
      save_ell = ell
      call qags(arg_pkell, 0._dp, 1._dp, epsabs, epsrel, ans, abserr, neval, ierr)
      if(ierr /= 0)  call warning('pk_multipole',ierr)  
      pk_multipole = ans*(2._dp*real(ell,kind=dp) + 1._dp)/alpha**3
    end function pk_multipole

    function arg_pkwed(muf)
      implicit none
      real(kind=dp) :: arg_pkwed, muf
      real(kind=dp) :: mu, k, muf2, fac
      muf2 = muf**2
      fac  = sqrt(muf2/alpha_lo**2 + (1._dp - muf2)/alpha_tr**2)
      k    = save_k*fac
      mu   = muf/alpha_lo/fac
      arg_pkwed = pk2d_at(k,mu)
    end function arg_pkwed

    function pk_wedge(k, iwedge, nwedges)
      real(kind=dp)       :: pk_wedge, k
      integer             :: iwedge, nwedges
      integer             :: ierr, neval
      real(kind=dp)       :: mu_min, mu_max, delta_mu
      real(kind=dp)       :: ans, abserr
      delta_mu = 1._dp/real(nwedges, kind=dp) !divide range from 0 to 1 into nwedges bins
      mu_min = real(iwedge-1,kind=dp)*delta_mu  
      mu_max = mu_min + delta_mu
      save_k = k
      call qags(arg_pkwed, mu_min, mu_max, epsabs, epsrel, ans, abserr, neval, ierr)
      if(ierr /= 0)  call warning('pk_wedges',ierr)  
      pk_wedge = ans/alpha**3/delta_mu
    end function pk_wedge
 
    subroutine get_model_pk_multi(kvec, nell, pvec)
      real(kind=dp), dimension(:), intent(in)  :: kvec
      integer, intent(in)                      :: nell
      real(kind=dp), dimension(:), intent(out) :: pvec
      integer                                  :: numk, i, ik, itot, ell
      numk = size(kvec)
      if(size(pvec) /= nell*numk)then
        call MPIStop('wrong size of pvec array in get_model_pk_multi')
      end if
      kmin = min((1._dp/alpha_lo),(1._dp/alpha_tr))*kvec(1)*0.9_dp
      kmax = max((1._dp/alpha_lo),(1._dp/alpha_tr))*kvec(numk)*1.1_dp
      call init_pk_multipoles
      itot = 0
      do i = 1, nell
        do ik = 1, numk
          itot = itot + 1
          ell = 2*(i-1)
          pvec(itot) = pk_multipole(kvec(ik), ell)
          if(ell ==0 )pvec(itot) = pvec(itot) + shot_noise
        end do
      end do
    end subroutine get_model_pk_multi 

    subroutine get_model_pk_wedges(kvec, nwedges, pvec)
      real(kind=dp), dimension(:), intent(in)  :: kvec
      integer, intent(in)                      :: nwedges
      real(kind=dp), dimension(:), intent(out) :: pvec
      integer                                  :: numk, iw, ik, itot
      numk = size(kvec)
      if(size(pvec) /= nwedges*numk)then
        call MPIStop('wrong size of pvec array in get_model_pk_wedges')
      end if
      kmin = min((1._dp/alpha_lo),(1._dp/alpha_tr))*kvec(1)*0.9_dp
      kmax = max((1._dp/alpha_lo),(1._dp/alpha_tr))*kvec(numk)*1.1_dp
      call init_pk_multipoles
      itot = 0
      do iw = 1, nwedges
        do ik = 1, numk
          itot = itot + 1
          pvec(itot) = pk_wedge(kvec(ik), iw, nwedges) + shot_noise
        end do
      end do
      !open(94,file='pkwedges_model.dat')
      !do ik = 1,numk
      !  write(94,'(9(E16.8,2x))')kvec(ik),(pvec(ik+(iw-1)*numk),iw=1,nwedges)
      !end do 
      !close(94)
      !write(*,*)'model written'
      !read(*,*)

    end subroutine get_model_pk_wedges 

    subroutine get_model_xi_multi(svec, nell, xivec)
      real(kind=dp), dimension(:), intent(in)  :: svec
      integer, intent(in)                      :: nell
      real(kind=dp), dimension(:), intent(out) :: xivec
      integer                                  :: nums, i, is, itot, ell
      nums = size(svec)
      if(size(xivec) /= nell*nums)then
        call MPIStop('wrong size of xivec array in get_model_xi_multi')
      end if
      smin = min(alpha_lo, alpha_tr)*svec(1)*0.9_dp
      smax = max(alpha_lo, alpha_tr)*svec(nums)*1.1_dp
      call init_xi_multipoles
      itot = 0
      do i = 1, nell
        do is = 1, nums
          itot = itot + 1
          ell = 2*(i-1)
          xivec(itot) = xi_multipole(svec(is), ell)
        end do
      end do

    end subroutine get_model_xi_multi 

    subroutine get_model_xi_wedges(svec, nwedges, xivec)
      real(kind=dp), dimension(:), intent(in)  :: svec
      integer, intent(in)                      :: nwedges
      real(kind=dp), dimension(:), intent(out) :: xivec
      integer                                  :: nums, iw, is, itot
      nums = size(svec)
      if(size(xivec) /= nwedges*nums)then
        call MPIStop('wrong size of xivec array in get_model_xi_multi')
      end if
      smin = min(alpha_lo, alpha_tr)*svec(1)*0.9_dp
      smax = max(alpha_lo, alpha_tr)*svec(nums)*1.1_dp
      call init_xi_multipoles
      itot = 0
      do iw = 1, nwedges
        do is = 1, nums
          itot = itot + 1
          xivec(itot) = xi_wedge(svec(is), iw, nwedges)
        end do
      end do
    end subroutine get_model_xi_wedges

    subroutine test_found_nan
      integer      :: i
      write(*,*)'found a NaN!!'
      open(98,file='multipoles_xi_NaN.dat',status='unknown',form='formatted')
      do i = 1,numr
        write(98,'(4(E16.8,2x))')r(i),xi_0(i),xi_2(i),xi_4(i)
      end do
      close(98)

      open(98,file='multipoles_pk_NaN.dat',status='unknown',form='formatted')
      do i = 1,numk
        write(98,'(4(E16.8,2x))')kvec(i),pk_0(i),pk_2(i),pk_4(i)
      end do
      close(98)
      call write_plin
      !stop
    end subroutine test_found_nan


end module compute_twopt

