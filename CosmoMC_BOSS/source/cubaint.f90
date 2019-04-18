!     Cuba Integration Library Header File
! ******************************************************

  module cubaint
    use settings, only: dp => mcp
    use NLRSD_settings, only: pi, pisqr
    implicit none

    !integer gridno   !this is notreally used
    real(kind=dp)             :: xmin, xmax
    !to fix rank missmatch with oneloop.f90
    integer, parameter        :: ncomp = 1

    contains

    subroutine integrate4D(integrand, method, epsrel, epsabs, result, error, prob)
      implicit none
      integer             :: method
      real(kind=dp)       :: epsabs, epsrel
      real(kind=dp)       :: result(1)
      external            :: integrand
! *******************************************************
      integer, parameter  :: ndim = 4    !chosen dimension
      integer, parameter  :: last = 4    !only last set of samples used
      integer, parameter  :: mineval = 0 !minimum integrand evaluations
      integer, parameter  :: maxeval = 110000000 !max evaluations
      !!VEGAS specific flags
!************  Verbosity flag for all integrals **********
      !0*8: sobol quasi random numbers used for sampling
      !1*4: only the last set of smples is used in the final
      !     result
      !1*2: echoes input parameters, verbosity flag
      integer, parameter  :: verbose   = 0*8+0*4+0*2+0*1
      integer, parameter  :: nstart    = 4000    !! 4000
      integer, parameter  :: nincrease = 700
      !! Suave specific flags
      integer, parameter  :: nnew = 200000
      real(kind=dp), parameter :: flatness = 50._dp
      !! Divonne specific flags !!
      !real(kind=dp), parameter    :: epsrel = 1D-2
      !real(kind=dp), parameter    :: epsabs = 1D-4
      !$$$      parameter (verbose = 0*8+0*4+1*2+0*1)
      integer, parameter     :: key1 = 9, key2 = 1, key3 = 1, maxpass = 25
      integer, parameter     :: ngiven = 0, ldxgiven = ndim, nextra = 0
      real(kind=dp)          :: border = 0.D0, maxchisq = 0.10D0, mindeviation = 0.125D0
!      double precision xgiven(ldxgiven,ngiven)
      !!Cuhre specific flags
      integer, parameter     :: key = 0 !default cubature rule for each dimension
      !!Output variables
      integer                        :: nregions, neval, fail
      integer                        :: ic !integer used to access and print the output arrays
      real(kind=dp),dimension(ncomp) ::  integral, error, prob
! *******************************************************

      if(method == 1 .or. method == 0) then
        call vegas(ndim, ncomp, integrand, &
           epsrel, epsabs, verbose, mineval, maxeval, &
           nstart, nincrease, &
           neval, fail, integral, error, prob)

        if (fail /= 0) then
          write(*,*) 'Failure to meet required accuracy'
          write(*,*) "VEGAS: neval =", neval," fail =", fail
          write(*,*) '("VEGAS:   ",E20.12," +- ",E20.12,"   p = ",F8.3)', &
               (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
          if(method /=  0) stop
        endif
!      endif

      elseif(method == 2 .or. method == 0) then

        call suave(ndim, ncomp, integrand, &
           epsrel, epsabs, verbose, mineval, maxeval,nnew,flatness, &
           nregions, neval, fail, integral, error, prob)

        if (fail .ne. 0) then
          write(*,*) 'Failure to meet required accuracy'
          write(*,*) "SUAVE: neval =", neval," fail =", fail
          write(*,*) '("SUAVE:   ",E20.12," +- ",E20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
          if(method /= 0) then
!           stop
          endif
        endif

      elseif(method == 3 .or. method == 0) then

        call cuhre(ndim, ncomp, integrand, &
           epsrel, epsabs, verbose, mineval, maxeval,key, &
           nregions, neval, fail, integral, error, prob)

        if (fail /= 0) then
          write(*,*) 'Failure to meet required accuracy'
          write(*,*) "CUHRE: neval =", neval," fail =", fail
          write(*,*) '("CUHRE:   ",E20.12," +- ",E20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
          if(method /= 0) then
!                stop
          endif
        endif
      elseif (method == 4 .or. method == 0) then
         call divonne(ndim, ncomp, integrand, &
              epsrel, epsabs, verbose, mineval, maxeval, &
              key1, key2, key3, maxpass, &
              border, maxchisq, mindeviation, &
              ngiven, ldxgiven, 0, nextra, 0, &
              nregions, neval, fail, integral, error, prob)
         if (method == 0 .or. fail == 1) then
            write(*,*) "DIVONNE: nregions =", nregions," neval =", neval, &
                 " fail =", fail
            write(*,*) '("DIVONNE: ",E20.12," +- ",E20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
         endif
      endif
      result = integral(1)
    end subroutine integrate4D

! *******************************************************
! *******************************************************
! ******************************************************

! *******************************************************
! *******************************************************
! ******************************************************

    subroutine integrate3D(integrand, method, epsrel, epsabs, result, error,prob)
      implicit none
      integer       :: verbose, ndim,mineval,maxeval,last
      integer       :: nstart,nincrease,method
      real(kind=dp) :: epsabs,epsrel
      real(kind=dp) :: result(1)
      external integrand

! *******************************************************
      parameter (ndim = 3) !chosen dimension
      parameter (last = 4) !only last set of samples used
      parameter (mineval = 0) !minimum integrand evaluations
!      parameter (maxeval = 100000) !max evaluations
      parameter (maxeval = 110000000) !max evaluations

!************  Verbosity flag for all integrals **********
      !0*8: sobol quasi random numbers used for sampling
      !1*4: only the last set of smples is used in the final
      !     result
      !1*2: echoes input parameters, verbosity flag
      parameter (verbose = 0*8+0*4+0*2+0*1)

      !!VEGAS specific flags
      parameter (nstart = 4000)    !! 4000
      parameter (nincrease = 700)


      !! Suave specific flags
      integer nnew
      double precision flatness
      parameter (nnew = 1000)
      parameter (flatness = 50D0)
      !! Divonne specific flags !!

!$$$      parameter (epsrel = 1D-2)
!$$$      parameter (epsabs = 1D-4)       !!-6
!$$$      parameter (verbose = 0*8+0*4+1*2+0*1)

      integer key1, key2, key3, maxpass
      double precision border, maxchisq, mindeviation
      integer ngiven, ldxgiven, nextra
!      parameter (key1 = 47)
      parameter (key1 = 9)
!      parameter (key2 = 100)
      parameter (key2 = 1)
! default key2=67
      parameter (key3 = 1)
      parameter (maxpass = 25) ! was 25
      parameter (border = 1.D-7) ! thickness of border. Important for
                                  ! kernels with angular singularities
      parameter (maxchisq = .10D0)  !!default 10_dp
      parameter (mindeviation = .25D0)  !!default .25_dp
      parameter (ngiven = 0)
      parameter (ldxgiven = ndim)
      parameter (nextra = 0)
!      double precision xgiven(ldxgiven,ngiven)

      !! Cuhre specific flags
      integer key
      parameter (key = 0) ! default cubature rule for each dimension

      !! Output variables
      real(kind=dp) integral(ncomp), error(ncomp), prob(ncomp)
      integer nregions, neval, fail

      integer ic !integer used to access and print the output arrays
! *******************************************************

      if(method .eq. 1 .or. method .eq.0) then

      call vegas(ndim, ncomp, integrand, &
           epsrel, epsabs, verbose, mineval, maxeval, &
           nstart, nincrease, &
           neval, fail, integral, error, prob)

      if (fail .ne. 0) then
            write(*,*) 'Failure to meet required accuracy'
            print *, "VEGAS: neval =", neval," fail =", fail
            print '("VEGAS:   ",E20.12," +- ",E20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
            if(method .ne. 0) then
!                  stop
            endif
            integral(1)=0._dp
      endif
!      endif

      elseif(method .eq. 2 .or. method .eq.0) then

      call suave(ndim, ncomp, integrand, &
           epsrel, epsabs, verbose, mineval, maxeval,nnew,flatness, &
           nregions, neval, fail, integral, error, prob)

      if (fail .ne. 0) then
            write(*,*) 'Failure to meet required accuracy'
            print *, "SUAVE: neval =", neval," fail =", fail
            print '("SUAVE:   ",E20.12," +- ",E20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
            if(method .ne. 0) then
!                  stop
            endif
!      endif
      endif

      elseif(method .eq. 3 .or. method .eq.0) then

      call cuhre(ndim, ncomp, integrand, &
           epsrel, epsabs, verbose, mineval, maxeval,key, &
           nregions, neval, fail, integral, error, prob)

      if (fail .ne. 0) then
            write(*,*) 'Failure to meet required accuracy'
            print *, "CUHRE: neval =", neval," fail =", fail
            print '("CUHRE:   ",E20.12," +- ",E20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
            if(method .ne. 0) then
!                  stop
            endif
            integral(1)=0._dp
      endif

      elseif (method.eq.4 .or. method.eq.0) then
         call divonne(ndim, ncomp, integrand, &
              epsrel, epsabs, verbose, mineval, maxeval, &
              key1, key2, key3, maxpass, &
              border, maxchisq, mindeviation, &
              ngiven, ldxgiven, 0, nextra, 0, &
              nregions, neval, fail, integral, error, prob)
         if (method.eq.0 .or. fail.eq.1) then
            print *, "DIVONNE: nregions =", nregions," neval =", neval, &
                 " fail =", fail
            print '("DIVONNE: ",E20.12," +- ",E20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
         endif
      endif

      result = integral(1)
    end subroutine integrate3D

! *******************************************************
! *******************************************************
! ******************************************************
! *******************************************************
! *******************************************************
! ******************************************************

    subroutine integrate2D(integrand, method, epsrel, epsabs, result, error,prob)
      implicit none
      integer        :: verbose, ndim,mineval,maxeval,last
      integer        :: nstart,nincrease,method
      real(kind=dp)  :: epsabs,epsrel
      real(kind=dp)  :: result(1)
      external       :: integrand

! *******************************************************
      parameter (ndim = 2) !chosen dimension
      parameter (last = 4) !only last set of samples used
      parameter (mineval = 0) !minimum integrand evaluations
      parameter (maxeval = 110000000) !max evaluations

!************  Verbosity flag for all integrals **********
      !0*8: sobol quasi random numbers used for sampling
      !1*4: only the last set of smples is used in the final
      !     result
      !1*2: echoes input parameters, verbosity flag
      parameter (verbose = 0*8+0*4+0*2+0*1)

      !!VEGAS specific flags
      parameter (nstart = 4000)    !! 4000
      parameter (nincrease = 700)

      !! Suave specific flags
      integer nnew
      double precision flatness
      parameter (nnew = 1000)
      parameter (flatness = 50D0)

      !! Divonne specific flags !!

!$$$      parameter (epsrel = 1D-2)
!$$$      parameter (epsabs = 1D-4)       !!-6
!$$$      parameter (verbose = 0*8+0*4+1*2+0*1)

      integer key1, key2, key3, maxpass
      double precision border, maxchisq, mindeviation
      integer ngiven, ldxgiven, nextra
!      parameter (key1 = 47)
      parameter (key1 = 9)
!      parameter (key2 = 100)
      parameter (key2 = 1)
! default key2=67
      parameter (key3 = 1)
      parameter (maxpass = 25)
      parameter (border = 0D0)
      parameter (maxchisq = .10D0)  !!default 10_dp
      parameter (mindeviation = .125D0)  !!default .25_dp
      parameter (ngiven = 0)
      parameter (ldxgiven = ndim)
      parameter (nextra = 0)
!      double precision xgiven(ldxgiven,ngiven)

      !! Cuhre specific flags
      integer key
      parameter (key = 13) ! default cubature rule for each dimension

      !! Output variables
      real(kind=dp) integral(ncomp), error(ncomp), prob(ncomp)
      integer nregions, neval, fail

      integer ic !integer used to access and print the output arrays
! *******************************************************

      if(method .eq. 1 .or. method .eq.0) then

      call vegas(ndim, ncomp, integrand, &
           epsrel, epsabs, verbose, mineval, maxeval, &
           nstart, nincrease, &
           neval, fail, integral, error, prob)

      if (fail .ne. 0) then
            write(*,*) 'Failure to meet required accuracy'
            print *, "VEGAS: neval =", neval," fail =", fail
            print '("VEGAS:   ",E20.12," +- ",E20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
            if(method .ne. 0) then
!                  stop
            endif
      endif
!      endif

      elseif(method .eq. 2 .or. method .eq.0) then

      call suave(ndim, ncomp, integrand, &
           epsrel, epsabs, verbose, mineval, maxeval,nnew,flatness, &
           nregions, neval, fail, integral, error, prob)

      if (fail .ne. 0) then
            write(*,*) 'Failure to meet required accuracy'
            print *, "SUAVE: neval =", neval," fail =", fail
            print '("SUAVE:   ",E20.12," +- ",E20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
            if(method .ne. 0) then
!                  stop
            endif
!      endif
      endif

      elseif(method .eq. 3 .or. method .eq.0) then

      call cuhre(ndim, ncomp, integrand, &
           epsrel, epsabs, verbose, mineval, maxeval,key, &
           nregions, neval, fail, integral, error, prob)

      if (fail .ne. 0) then
            write(*,*) 'Failure to meet required accuracy'
            print *, "CUHRE: neval =", neval," fail =", fail
            print '("CUHRE:   ",E20.12," +- ",E20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
            if(method .ne. 0) then
!                  stop
            endif
            integral(1)=0._dp
      endif
      elseif (method.eq.4 .or. method.eq.0) then
         call divonne(ndim, ncomp, integrand, &
              epsrel, epsabs, verbose, mineval, maxeval, &
              key1, key2, key3, maxpass, &
              border, maxchisq, mindeviation, &
              ngiven, ldxgiven, 0, nextra, 0, &
              nregions, neval, fail, integral, error, prob)
         if (method.eq.0 .or. fail.eq.1) then
            print *, "DIVONNE: nregions =", nregions," neval =", neval, &
                 " fail =", fail
            print '("DIVONNE: ",E20.12," +- ",E20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
         endif
      endif
      result = integral(1)
    end subroutine integrate2D

! *******************************************************
!     Emiliano's 2D integral

    subroutine integrale2D(subint,iint,xepsrel,xepsabs,pout)
      implicit none
      integer xverbose,iint
      real(kind=dp) pout,xepsrel,xepsabs
      external subint

            !! Common arguments !! (last = ??)

      integer ndim,  mineval, maxeval, verbose, last
      double precision epsrel, epsabs
      parameter (ndim = 2)
!$$$      parameter (epsrel = 1D-2)
!$$$      parameter (epsabs = 1D-6)
!$$$      parameter (verbose = 0*8+0*4+1*2+0*1)
      parameter (last = 4)
      parameter (mineval = 0)
      parameter (maxeval = 110000000)

      !! VEGAS !!

      parameter (epsrel = 1D-2)
      parameter (epsabs = 0)
!      parameter (epsrel = 1D-3)
!      parameter (epsabs = 0)

      parameter (verbose = 0*8+0*4+0*2+0*1)

      integer nstart, nincrease
      parameter (nstart = 4000)    !! 4000
      parameter (nincrease = 700)

      !! SUAVE !!

!$$$      parameter (epsrel = 1D-2)
!$$$      parameter (epsabs = 1D-4)
!$$$      parameter (verbose = 0*8+1*4+1*2+0*1)

      integer nnew
      double precision flatness
      parameter (nnew = 1000)
      parameter (flatness = 50D0)

      !! DIVONNE !!

!$$$      parameter (epsrel = 1D-2)
!$$$      parameter (epsabs = 1D-4)       !!-6
!$$$      parameter (verbose = 0*8+0*4+1*2+0*1)


      integer key1, key2, key3, maxpass
      double precision border, maxchisq, mindeviation
      integer ngiven, ldxgiven, nextra
!      parameter (key1 = 47)
      parameter (key1 = 9)
!      parameter (key2 = 100)
      parameter (key2 = 1)
! default key2=67
      parameter (key3 = 1)
      parameter (maxpass = 25)
      parameter (border = 0D0)
      parameter (maxchisq = .10D0)  !!default 10_dp
      parameter (mindeviation = .125D0)  !!default .25_dp
      parameter (ngiven = 0)
      parameter (ldxgiven = ndim)
      parameter (nextra = 0)
!      double precision xgiven(ldxgiven,ngiven)

      !! CUHRE !!

!$$$      parameter (epsrel = 1D-2)
!$$$      parameter (epsabs = 1D-4)       !!-6
!$$$      parameter (verbose = 0*8+1*4+1*2+0*1)

      integer key
      parameter (key = 0)

      !! Output !!

!      external integrand

      real(kind=dp) integral(ncomp), error(ncomp), prob(ncomp)
      integer nregions, neval, fail

      integer ic

      xverbose=0*8+0*4+0*2+0*1

      if (iint.eq.1 .or. iint.eq.0) then

         call vegas(ndim, ncomp, subint,xepsrel, xepsabs, xverbose, &
              mineval,maxeval,nstart, nincrease,neval,fail,integral, &
              error,prob)
         if (iint.eq.0) then
            print *, "VEGAS: neval =", neval," fail =", fail
            print '("VEGAS: ",F20.12," +- ",F20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
         endif

      elseif (iint.eq.2 .or. iint.eq.0) then

            call suave(ndim, ncomp, subint,xepsrel, xepsabs, &
                 xverbose + last, mineval, maxeval,nnew, flatness, &
                 nregions, neval, fail, integral, error, prob)
         if (iint.eq.0 .or. fail.eq.1) then
            print *, "SUAVE: nregions =", nregions," neval =", neval, &
                 " fail =", fail
            print '("SUAVE: ",F20.12," +- ",F20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
         endif

      elseif (iint.eq.3 .or. iint.eq.0) then

         call divonne(ndim,ncomp,subint,xepsrel,xepsabs,xverbose, &
              mineval,maxeval,key1,key2,key3,maxpass,border,maxchisq, &
              mindeviation,ngiven,ldxgiven,0,nextra,0,nregions,neval, &
              fail,integral,error,prob)
         if (iint.eq.0 .or. fail.eq.1) then
            print *, "DIVON: nregions =", nregions," neval =", neval, &
                 " fail =", fail
            print '("DIVON: ",F20.12," +- ",F20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
         endif

      elseif (iint.eq.4 .or. iint.eq.0) then

         call cuhre(ndim,ncomp,subint,xepsrel,xepsabs,xverbose+last, &
              mineval,maxeval,key,nregions,neval,fail,integral,error, &
              prob)
         if (iint.eq.0 .or. fail.eq.1) then
            print *, "CUHRE: nregions =", nregions," neval =", neval, &
                 " fail =", fail
            print '("CUHRE: ",F20.12," +- ",F20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
         endif

      endif
      pout=integral(1)
      return
    end subroutine integrale2D

! *******************************************************
!     Emiliano's 1D integral

    subroutine integrate1D(subint,method,xepsrel,xepsabs,result, error,prob)
      implicit none
      integer       :: xverbose,method
      real(kind=dp) :: xepsrel,xepsabs
      real(kind=dp) :: result(1)
      external      :: subint

            !! Common arguments !! (last = ??)

      integer ndim, mineval, maxeval, verbose, last
      double precision epsrel, epsabs
      parameter (ndim = 1)
!$$$      parameter (epsrel = 1D-2)
!$$$      parameter (epsabs = 1D-6)
!$$$      parameter (verbose = 0*8+0*4+1*2+0*1)
      parameter (last = 4)
      parameter (mineval = 0)
      parameter (maxeval = 110000000)

      !! VEGAS !!

      parameter (epsrel = 1D-2)
      parameter (epsabs = 0)
!      parameter (epsrel = 1D-3)
!      parameter (epsabs = 0)

      parameter (verbose = 0*8+0*4+0*2+0*1)

      integer nstart, nincrease
      parameter (nstart = 4000)    !! 4000
      parameter (nincrease = 700)

      !! SUAVE !!

!$$$      parameter (epsrel = 1D-2)
!$$$      parameter (epsabs = 1D-4)
!$$$      parameter (verbose = 0*8+1*4+1*2+0*1)

      integer nnew
      double precision flatness
      parameter (nnew = 1000)
      parameter (flatness = 50D0)

      !! DIVONNE !!

!$$$      parameter (epsrel = 1D-2)
!$$$      parameter (epsabs = 1D-4)       !!-6
!$$$      parameter (verbose = 0*8+0*4+1*2+0*1)


      integer key1, key2, key3, maxpass
      double precision border, maxchisq, mindeviation
      integer ngiven, ldxgiven, nextra
!      parameter (key1 = 47)
      parameter (key1 = 9)
!      parameter (key2 = 100)
      parameter (key2 = 1)
! default key2=67
      parameter (key3 = 1)
      parameter (maxpass = 25)
      parameter (border = 0D0)
      parameter (maxchisq = .10D0)  !!default 10_dp
      parameter (mindeviation = .125D0)  !!default .25_dp
      parameter (ngiven = 0)
      parameter (ldxgiven = ndim)
      parameter (nextra = 0)
!      double precision xgiven(ldxgiven,ngiven)

      !! CUHRE !!

!$$$      parameter (epsrel = 1D-2)
!$$$      parameter (epsabs = 1D-4)       !!-6
!$$$      parameter (verbose = 0*8+1*4+1*2+0*1)

      integer key
      parameter (key = 0)

      !! Output !!

!      external integrand

      real(kind=dp) integral(ncomp), error(ncomp), prob(ncomp)
      integer nregions, neval, fail

      integer ic

      xverbose=0*8+0*4+0*2+0*1

      if (method.eq.1 .or. method.eq.0) then

         call vegas(ndim, ncomp, subint,xepsrel, xepsabs, xverbose, &
              mineval,maxeval,nstart, nincrease,neval,fail,integral, &
              error,prob)
         if (method.eq.0) then
            print *, "VEGAS: neval =", neval," fail =", fail
            print '("VEGAS: ",F20.12," +- ",F20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
         endif

      elseif (method.eq.2 .or. method.eq.0) then

            call suave(ndim, ncomp, subint,xepsrel, xepsabs, &
                 xverbose + last, mineval, maxeval,nnew, flatness, &
                 nregions, neval, fail, integral, error, prob)
         if (method.eq.0 .or. fail.eq.1) then
            print *, "SUAVE: nregions =", nregions," neval =", neval, &
                 " fail =", fail
            print '("SUAVE: ",F20.12," +- ",F20.12,"   p = ",F8.3)', &
                 (integral(ic), error(ic), prob(ic), ic = 1, ncomp)
         endif

      elseif (method.eq.4 .or. method.eq.3) then
            write(*,*) 'integrate 1-D: method',method,' out of bounds'
            stop
      endif

      result=integral(1)

    end subroutine integrate1D

! **************************** TESTING SUBROUTINE **********************
    subroutine test_cuba()
      IMPLICIT NONE
      real(kind=dp)   :: t1, t2
      real(kind=dp)   :: relacc, absacc
      real(kind=dp)   :: result(1), error(1), prob(1)
      integer         :: method, i
      xmin=0._dp
      xmax=100._dp
      relacc=1.d-3
      absacc=1.d-10
      method=1
      write(*,*) '****************** CUBA TEST *********************'
      write(*,*) '      testing 1-D integrators...'
      do i=1,2
        method=i
        call cpu_time(t1)
        call integrate1D(testc1,method,relacc,absacc,result, error,prob)
        if (abs(result(1)-5000._dp) .gt. 5000*relacc) then
          write(*,*) 'Bad f(x)=x accuracy!'
          write(*,*) 'result',result(1)
          stop
        endif
        call cpu_time(t2)
        write(*,*) '            method:',method,t2-t1,' seconds'
      enddo
      write(*,*) '      testing 2-D integrators...'
      do i=1,4
        method=i
        call cpu_time(t1)
        call integrate2D(testc2,method,relacc,absacc,result, error,prob)
!            write(*,*) 'result=',result(1)
        if (abs(result(1)-4.18879d6) .gt. 4.18879d6*relacc) then
          write(*,*) 'Bad accuracy for method:',method
          stop
        endif
        call cpu_time(t2)
        write(*,*) '            method:',method,t2-t1,' seconds'
      enddo
      write(*,*) '      testing 3-D integrators...'
      do i=1,4
        method=i
        call cpu_time(t1)
        call integrate3D(testc3,method,relacc,absacc,result, error,prob)
!        write(*,*) 'result=',result(1)
        if (abs(result(1)-4.18879d6) .gt. 4.18879d6*relacc) then
          write(*,*) 'Bad accuracy for method:',method
          stop
        endif
        call cpu_time(t2)
        write(*,*) '            method:',method,t2-t1,' seconds'
      enddo
      write(*,*) '      testing 4-D integrators...'
      do i=1,4
        method=i
        call cpu_time(t1)
        call integrate4D(testc4,method,relacc,absacc,result, error,prob)
!        write(*,*) 'result=',result(1)
        if (abs(result(1)-2.09413d10) .gt. 2.09413d10*relacc) then
          write(*,*) 'Bad accuracy for method:',method
          write(*,*) 'result',result(1)
          stop
        endif
        call cpu_time(t2)
        write(*,*) '            method:',method,t2-t1,' seconds'
      enddo

    end subroutine test_cuba

    subroutine testc1(ndim,x,nncomp,f)
      implicit none
      integer ndim,nncomp,i
      real(kind=dp) x(ndim),f(nncomp),q(ndim),j
      do i=1,ndim
        q(i)=xmin+(xmax-xmin)*x(i)
      enddo
      f(1)=q(1)
      j=(xmax-xmin)
      f(1)=f(1)*j
    end subroutine testc1

    subroutine testc2(ndim,x,nncomp,f)
      implicit none
      integer ndim,nncomp
      real(kind=dp) x(ndim),f(nncomp),q(ndim),j
      q(1)=xmin+(xmax-xmin)*x(1)    !  r
      q(2)=pi*x(2)                 ! theta
      f(1)=1._dp
      j=2._dp*pisqr*(xmax-xmin)*q(1)**2*sin(q(2))
      f(1)=f(1)*j
    end subroutine testc2

    subroutine testc3(ndim,x,nncomp,f)
      implicit none
      integer ndim,nncomp
      real(kind=dp) x(ndim),f(nncomp),q(ndim),j
      q(1) = xmin + (xmax - xmin)*x(1)    !  r
      q(2) = pi*x(2)                 ! theta
      q(3) = 2._dp*pi*x(3)             ! phi
      f(1) = 1._dp
      j    = 2._dp*pisqr*(xmax - xmin)*q(1)**2*sin(q(2))
      f(1) = f(1)*j
    end subroutine testc3

    subroutine testc4(ndim,x,nncomp,f)
      implicit none
      integer ndim,nncomp
      real(kind=dp) x(ndim),f(nncomp),q(ndim),j
      q(1)=xmin+(xmax-xmin)*x(1)    !  r
      q(2)=pi*x(2)                 ! theta
      q(3)=2._dp*pi*x(3)             ! phi
      q(4)=xmin+(xmax-xmin)*x(4)
      f(1)=q(4)!1._dp
      j=2._dp*pisqr*(xmax-xmin)**2*q(1)**2*sin(q(2))
      f(1)=f(1)*j
    end subroutine testc4

  end module cubaint










