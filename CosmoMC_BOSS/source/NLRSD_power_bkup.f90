module model_power
  use settings,        only : dp => mcp, feedback
  use numrec
  use quadpack 
  use cubaint
  use NLRSD_settings,  only :  Matterpowerat, fog_model, fac_norm, real_space, tns_full, & 
                              b1, b2, gam2, gam3minus, sigv_P, f, a_roman, sigv_B, A_fast, &
                              do_init_model , init_plin
  use nonlinear_power, only : gRPTpower, init_za, get_root => root, &
                              I0, I2, sphbe0, sphbe1, sphbe2, sphbe3, &
                              sphbe4, sphbe5, sphbe6, sphbe7, sphbe8, sigv
  use bias,            only : Pk_b1b2, Pk_b1g2_prop, Pk_b1g2_mc, Pk_b2b2, Pk_b2g2, &
                              Pk_g2g2, Pk_b1g3m_prop, Pk_b2, Pk_g2_prop, Pk_g2_mc

  implicit none
  real(kind=dp),  private, parameter           :: kmin = 0.0001_dp, kmax = 2.4_dp
  integer, private, parameter                  :: nk_table = 120
  integer, private, parameter                  :: nq_table = 120
  integer, private, parameter                  :: nr_table = 120
  real(kind=dp),  private, dimension(nk_table) :: k_table
  real(kind=dp),  private, dimension(nq_table) :: q_table
  real(kind=dp),  private, dimension(nr_table) :: r_table

  real(kind=dp),  private, dimension(nk_table) :: dP13_dd_table, dP13_dv_table, dP13_vv_table
  real(kind=dp),  private, dimension(nk_table) :: dP13_dd_d_table, dP13_dv_d_table, dP13_vv_d_table
  real(kind=dp),  private, dimension(nk_table) :: P22inv_dd_table, P22inv_dv_table, P22inv_vv_table
  real(kind=dp),  private, dimension(nk_table) :: P22inv_dd_d_table, P22inv_dv_d_table, P22inv_vv_d_table
  real(kind=dp),  private, dimension(nk_table) :: gRPT_root_table, gRPT_root_d_table

  real(kind=dp),  private, dimension(nk_table) :: b1b2_table, b1g2_table, b2b2_table
  real(kind=dp),  private, dimension(nk_table) :: b2g2_table, g2g2_table, b1g3m_table
  real(kind=dp),  private, dimension(nk_table) :: b1b2_d_table, b1g2_d_table, b2b2_d_table
  real(kind=dp),  private, dimension(nk_table) :: b2g2_d_table, g2g2_d_table, b1g3m_d_table

  real(kind=dp),  private, dimension(nk_table) :: b2_table, g2_table, g3m_table
  real(kind=dp),  private, dimension(nk_table) :: b2_d_table, g2_d_table, g3m_d_table

  integer, private  :: ndenvel, ell, nv2, nmu, ellp 
  real(kind=dp),  private  :: k, q, vzdisp, xsq, ysq, zsq, vsq, rpsi, rpw  
  real(kind=dp),  private, parameter :: relacc = 0.01_dp, absacc = 0.01_dp

  real(kind=dp),  private, parameter :: kminRSD = 0.0001_dp, kmaxRSD = 10._dp
  real(kind=dp),  private, parameter :: relaccRSD = 0.001_dp, absaccRSD = 0.01_dp

  real(kind=dp),  private, dimension(nk_table) :: fmu2v1sqb1sq_table, f3mu6v2sq_table, f3mu4v2sq_table
  real(kind=dp),  private, dimension(nk_table) :: fmu2v1sqb1sq_d_table, f3mu6v2sq_d_table, f3mu4v2sq_d_table

  real(kind=dp),  private, dimension(nk_table) :: f2mu4v1v2b1_table, f2mu2v1v2b1_table, fmu2v1sqb1b2_table
  real(kind=dp),  private, dimension(nk_table) :: f2mu4v1v2b1_d_table, f2mu2v1v2b1_d_table, fmu2v1sqb1b2_d_table

  real(kind=dp),  private, dimension(nk_table) :: fmu2v1sqb1g2_table, f2mu4v1v2b2_table, f2mu4v1v2g2_table
  real(kind=dp),  private, dimension(nk_table) :: fmu2v1sqb1g2_d_table, f2mu4v1v2b2_d_table, f2mu4v1v2g2_d_table

  real(kind=dp),  private, dimension(nk_table) :: f2mu2v1v2b2_table, f2mu2v1v2g2_table 
  real(kind=dp),  private, dimension(nk_table) :: f2mu2v1v2b2_d_table, f2mu2v1v2g2_d_table 


  real(kind=dp),  private, dimension(nk_table) :: f4mu4v2sq_table, f4mu4v2sq_d_table 
  real(kind=dp),  private, dimension(nk_table) :: f4mu6v2sq_table, f4mu6v2sq_d_table 
  real(kind=dp),  private, dimension(nk_table) :: f4mu8v2sq_table, f4mu8v2sq_d_table 

  real(kind=dp),  private, dimension(nk_table) :: f3mu4v1v2b1_table, f3mu4v1v2b1_d_table 
  real(kind=dp),  private, dimension(nk_table) :: f3mu6v1v2b1_table, f3mu6v1v2b1_d_table 

  real(kind=dp),  private, dimension(nk_table) :: f2mu2v1sqb1sq_table, f2mu2v1sqb1sq_d_table 
  real(kind=dp),  private, dimension(nk_table) :: f2mu4v1sqb1sq_table, f2mu4v1sqb1sq_d_table 

  real(kind=dp),  private, dimension(nk_table) :: f3mu2v1v2b1_table, f3mu2v1v2b1_d_table 
  real(kind=dp),  private, dimension(nk_table) :: f4mu2v2sq_table, f4mu2v2sq_d_table 

  real(kind=dp),  private, dimension(nq_table) :: PlinDW_table, PlinDW_d_table
  
  real(kind=dp),  private, parameter  :: rmin = 0.01_dp, rmax= 20000._dp

  real(kind=dp),  private, parameter  :: deltaq = 0.001_dp

  real(kind=dp),  private, parameter  :: k1loop = 0.01_dp !kmin at which we start calculating 1loop

  real(kind=dp),  private, dimension(nk_table) :: rsd_convol0_table, rsd_convol0_d_table 
  real(kind=dp),  private, dimension(nk_table) :: rsd_convol2_table, rsd_convol2_d_table 
  real(kind=dp),  private, dimension(nk_table) :: rsd_convol4_table, rsd_convol4_d_table 

  real(kind=dp),  private, dimension(nr_table) :: Psi1m1_table, Psi1m1_d_table 

  integer,  private, dimension(nk_table) :: nsimpson, ibegin
  
  real(kind=dp),  private, parameter  :: safety = 5._dp !safety factor for scale-dep integrals (rel and abs accuracy)

  contains

!************************************************************
! build tables
!************************************************************

    subroutine init_model()
      implicit none 
      real(kind=dp)                :: oldtime, tl, x
      call init_plin
      call init_kvec
      call init_gRPT
      !write(*,*)'in bias_galgal'
      !call cpu_time(oldtime)
      !tl = 0.
      call init_bias_galgal
     ! call cpu_time(x)
     ! tl = x-oldtime
     ! write(*,*)'time=',tl
     ! oldtime = x
      if (f>0._dp) then  
      !write(*,*)'in bias_galvel'
        call init_bias_galvel
      !   call cpu_time(x)
      !   tl = x-oldtime
      !   write(*,*)'time=',tl
      !   oldtime = x
      !write(*,*)'in rsd_velvel (slow)'
        call init_rsd_velvel
      !   call cpu_time(x)
      !   tl = x-oldtime
      !   write(*,*)'time=',tl
      !   oldtime = x
      !write(*,*)'in rsd_galvel (slow)'
        call init_rsd_galvel
      !   call cpu_time(x)
       !  tl = x-oldtime
       !  write(*,*)'time=',tl
       !  oldtime = x
     ! write(*,*)'in rsd_galgal'
         call init_rsd_galgal
      !   call cpu_time(x)
      !   tl = x-oldtime
      !   write(*,*)'time=',tl
      !   oldtime = x
!         call init_VDskew
!         call tns_check
         if(feedback > 1)write(*,'(a)')'> Code Initialized'
!         call init_rsd_convolution
      else
         if(feedback > 1)write(*,*)'running in real-space mode: quadrupole becomes Pgm'
      endif             
       do_init_model = .false.      
    end subroutine init_model

    subroutine init_kvec()
      implicit none
      integer           :: i
      real(kind=dp)     :: qmin, qmax, pow, kcenter

      kcenter = 0.8_dp
      qmin = log10(kmin)
      qmax = log10(kmax)
      pow = 1.5_dp
      qmin = croot(qmin + kcenter, pow)
      qmax = croot(qmax + kcenter, pow)
      do i=1,nk_table
        !linear
 !      k_table(i) = (kmax - kmin)*dble(i-1)/dble(nk_table - 1) + kmin  
        k_table(i) = (qmax - qmin)*dble(i-1)/dble(nk_table - 1) + qmin  
        k_table(i) = sign(abs(k_table(i))**pow,k_table(i)) - kcenter 
      enddo
      k_table = 10._dp**k_table
      do i=1,nq_table
        !linear
        q_table(i) = (kmax - kmin)*dble(i-1)/dble(nq_table - 1) + kmin  
      enddo
      do i=1,nr_table
        !linear
  !      r_table(i) = (rmax - rmin)*dble(i-1)/dble(nr_table - 1) + rmin  
        !log
        r_table(i) = 10._dp**(log10(rmax/rmin)*dble(i-1)/dble(nr_table-1)+log10(rmin))
      enddo
      
      contains
        function croot(x,apow)
          real(kind=dp)  :: croot, a, x, apow
          a = abs(x)**(1._dp/apow)
          croot = sign(a,x)
        end function croot
    end subroutine init_kvec

    subroutine init_gRPT()
      implicit none
      integer        :: i
      real(kind=dp)  :: ak, root, dP13_dd, dP13_dv, dP13_vv,P22inv_dd,P22inv_dv,P22inv_vv

      ak=Matterpowerat(0.1_dp)
      !open(93,file='../bin/fort.43',status='unknown',form='formatted')

      init_za = .true. !initialize Zeldovich

      do i=1,nk_table
        ak = k_table(i)
        if (ak < k1loop) then

          root             = get_root(ak)
          dP13_dd_table(i) = 0._dp
          dP13_dv_table(i) = 0._dp
          dP13_vv_table(i) = 0._dp

          P22inv_dd_table(i) = 0._dp
          P22inv_dv_table(i) = 0._dp
          P22inv_vv_table(i) = 0._dp

          gRPT_root_table(i) = root

        else

          call gRPTpower(ak, root, dP13_dd, dP13_dv, dP13_vv, P22inv_dd, P22inv_dv, P22inv_vv)

          dP13_dd_table(i) = dP13_dd
          dP13_dv_table(i) = dP13_dv
          dP13_vv_table(i) = dP13_vv

          P22inv_dd_table(i) = P22inv_dd
          P22inv_dv_table(i) = P22inv_dv
          P22inv_vv_table(i) = P22inv_vv

          gRPT_root_table(i) = root

        end if
    
      enddo 
      if(feedback > 1)write(*,'(A)')'> Done with init gRPT'
      !close(93)

      call spline(k_table, dP13_dd_table, nk_table, 3d30,3d30, dP13_dd_d_table)
      call spline(k_table, dP13_dv_table, nk_table, 3d30,3d30, dP13_dv_d_table)
      call spline(k_table, dP13_vv_table, nk_table, 3d30,3d30, dP13_vv_d_table)

      call spline(k_table, P22inv_dd_table, nk_table, 3d30,3d30, P22inv_dd_d_table)
      call spline(k_table, P22inv_dv_table, nk_table, 3d30,3d30, P22inv_dv_d_table)
      call spline(k_table, P22inv_vv_table, nk_table, 3d30,3d30, P22inv_vv_d_table)

      call spline(k_table, gRPT_root_table,  nk_table, 3d30,3d30, gRPT_root_d_table)


    end subroutine init_gRPT

    subroutine init_bias_galgal() !checked
      implicit   none
      integer         :: i
      real(kind=dp)   :: ak

      do i = 1, nk_table
        ak = k_table(i)

        if (ak >= k1loop) then !tabulate 1loop integrals
          b1b2_table(i)  = Pk_b1b2(ak)
          b1g2_table(i)  = Pk_b1g2(ak)
          b2b2_table(i)  = Pk_b2b2(ak)
          b2g2_table(i)  = Pk_b2g2(ak)
          g2g2_table(i)  = Pk_g2g2(ak) !was b1g2, bug!!!!
          b1g3m_table(i) = Pk_b1g3m(ak)
        else
          b1b2_table(i)  = 0._dp
          b1g2_table(i)  = 0._dp
          b2b2_table(i)  = 0._dp
          b2g2_table(i)  = 0._dp
          g2g2_table(i)  = 0._dp
          b1g3m_table(i) = 0._dp
        endif        
      enddo 

      call spline(k_table, b1b2_table, nk_table, 3d30, 3d30, b1b2_d_table)
      call spline(k_table, b1g2_table, nk_table, 3d30, 3d30, b1g2_d_table)
      call spline(k_table, b2b2_table, nk_table, 3d30, 3d30, b2b2_d_table)
      call spline(k_table, b2g2_table, nk_table, 3d30, 3d30, b2g2_d_table)
      call spline(k_table, g2g2_table, nk_table, 3d30, 3d30, g2g2_d_table)
      call spline(k_table, b1g3m_table, nk_table, 3d30, 3d30, b1g3m_d_table)

    end subroutine init_bias_galgal
    
    subroutine init_bias_galvel() !checked
      implicit none
      integer       :: i
      real(kind=dp) :: ak

      do i = 1,nk_table
        ak = k_table(i)
        if (ak >= k1loop) then !tabulate 1loop integrals
          b2_table(i)   =  Pk_b2(ak)
          g2_table(i)   =  Pk_g2(ak)
          g3m_table(i)  =  Pk_g3m(ak)
        else
          b2_table(i)  = 0._dp
          g2_table(i)  = 0._dp
          g3m_table(i) = 0._dp
        endif
      enddo 

      call spline(k_table, b2_table, nk_table, 3d30, 3d30, b2_d_table)
      call spline(k_table, g2_table, nk_table, 3d30, 3d30, g2_d_table)
      call spline(k_table, g3m_table, nk_table, 3d30, 3d30, g3m_d_table)

    end subroutine init_bias_galvel

    subroutine init_rsd_velvel()
      implicit   none
      integer         :: i
      real(kind=dp)   :: ak
    

      do i = 1, nk_table
        ak = k_table(i)
        if (ak >= k1loop) then !tabulate 1loop integrals
          f3mu4v2sq_table(i)  = Pk_f3mu4v2sq(ak)
          f3mu6v2sq_table(i)  = Pk_f3mu6v2sq(ak)
        
          f4mu2v2sq_table(i)  = Pk_f4mu2v2sq(ak)
          f4mu4v2sq_table(i)  = Pk_f4mu4v2sq(ak)
          f4mu6v2sq_table(i)  = Pk_f4mu6v2sq(ak)
          f4mu8v2sq_table(i)  = Pk_f4mu8v2sq(ak)
        else
          f3mu4v2sq_table(i)  = 0._dp
          f3mu6v2sq_table(i)  = 0._dp
        
          f4mu2v2sq_table(i)  = 0._dp
          f4mu4v2sq_table(i)  = 0._dp
          f4mu6v2sq_table(i)  = 0._dp
          f4mu8v2sq_table(i)  = 0._dp
        endif
      enddo 

      call spline(k_table, f3mu4v2sq_table, nk_table, 3d30, 3d30, f3mu4v2sq_d_table)
      call spline(k_table, f3mu6v2sq_table, nk_table, 3d30, 3d30, f3mu6v2sq_d_table)

      call spline(k_table, f4mu2v2sq_table, nk_table, 3d30, 3d30, f4mu2v2sq_d_table)
      call spline(k_table, f4mu4v2sq_table, nk_table, 3d30, 3d30, f4mu4v2sq_d_table)
      call spline(k_table, f4mu6v2sq_table, nk_table, 3d30, 3d30, f4mu6v2sq_d_table)
      call spline(k_table, f4mu8v2sq_table, nk_table, 3d30, 3d30, f4mu8v2sq_d_table)
    end subroutine init_rsd_velvel

    subroutine init_rsd_galgal()
      implicit   none
      integer         :: i
      real(kind=dp)   :: ak

      do i = 1, nk_table
        ak = k_table(i)

        if (ak >= k1loop) then !tabulate 1loop integrals
          fmu2v1sqb1sq_table(i)  = Pk_fmu2v1sqb1sq(ak)  
          fmu2v1sqb1b2_table(i)  = Pk_fmu2v1sqb1b2(ak)
          fmu2v1sqb1g2_table(i)  = Pk_fmu2v1sqb1g2(ak)

          f2mu2v1sqb1sq_table(i)  = Pk_f2mu2v1sqb1sq(ak)
          f2mu4v1sqb1sq_table(i)  = Pk_f2mu4v1sqb1sq(ak)
        else
          fmu2v1sqb1sq_table(i)  = 0._dp
          fmu2v1sqb1b2_table(i)  = 0._dp
          fmu2v1sqb1g2_table(i)  = 0._dp

          f2mu2v1sqb1sq_table(i)  = 0._dp
          f2mu4v1sqb1sq_table(i)  = 0._dp
        endif
      enddo 

      call spline(k_table, fmu2v1sqb1sq_table, nk_table, 3d30, 3d30, fmu2v1sqb1sq_d_table)
      call spline(k_table, fmu2v1sqb1b2_table, nk_table, 3d30, 3d30, fmu2v1sqb1b2_d_table)
      call spline(k_table, fmu2v1sqb1g2_table, nk_table, 3d30, 3d30, fmu2v1sqb1g2_d_table)

      call spline(k_table, f2mu2v1sqb1sq_table, nk_table, 3d30, 3d30, f2mu2v1sqb1sq_d_table)
      call spline(k_table, f2mu4v1sqb1sq_table, nk_table, 3d30, 3d30, f2mu4v1sqb1sq_d_table)
    end subroutine init_rsd_galgal

    subroutine init_rsd_galvel()
      implicit   none
      integer         :: i
      real(kind=dp)   :: ak

      do i = 1, nk_table
        ak = k_table(i)

        if (ak >= k1loop) then !tabulate 1loop integrals
          f2mu4v1v2b1_table(i)  = Pk_f2mu4v1v2b1(ak)
          f2mu2v1v2b1_table(i)  = Pk_f2mu2v1v2b1(ak)
          f2mu4v1v2b2_table(i)  = Pk_f2mu4v1v2b2(ak)
          f2mu2v1v2b2_table(i)  = Pk_f2mu2v1v2b2(ak)
          f2mu4v1v2g2_table(i)  = Pk_f2mu4v1v2g2(ak)
          f2mu2v1v2g2_table(i)  = Pk_f2mu2v1v2g2(ak)
  
          f3mu4v1v2b1_table(i)  = Pk_f3mu4v1v2b1(ak)
          f3mu6v1v2b1_table(i)  = Pk_f3mu6v1v2b1(ak)

          f3mu2v1v2b1_table(i)  = Pk_f3mu2v1v2b1(ak)
        else
          f2mu4v1v2b1_table(i)  = 0._dp
          f2mu2v1v2b1_table(i)  = 0._dp
          f2mu4v1v2b2_table(i)  = 0._dp
          f2mu2v1v2b2_table(i)  = 0._dp
          f2mu4v1v2g2_table(i)  = 0._dp
          f2mu2v1v2g2_table(i)  = 0._dp

          f3mu4v1v2b1_table(i)  = 0._dp
          f3mu6v1v2b1_table(i)  = 0._dp

          f3mu2v1v2b1_table(i)  = 0._dp
        endif        
      enddo 

      call spline(k_table, f2mu4v1v2b1_table, nk_table, 3d30, 3d30, f2mu4v1v2b1_d_table)
      call spline(k_table, f2mu2v1v2b1_table, nk_table, 3d30, 3d30, f2mu2v1v2b1_d_table)
      call spline(k_table, f2mu4v1v2b2_table, nk_table, 3d30, 3d30, f2mu4v1v2b2_d_table)
      call spline(k_table, f2mu2v1v2b2_table, nk_table, 3d30, 3d30, f2mu2v1v2b2_d_table)
      call spline(k_table, f2mu4v1v2g2_table, nk_table, 3d30, 3d30, f2mu4v1v2g2_d_table)
      call spline(k_table, f2mu2v1v2g2_table, nk_table, 3d30, 3d30, f2mu2v1v2g2_d_table)

      call spline(k_table, f3mu4v1v2b1_table, nk_table, 3d30, 3d30, f3mu4v1v2b1_d_table)
      call spline(k_table, f3mu6v1v2b1_table, nk_table, 3d30, 3d30, f3mu6v1v2b1_d_table)

      call spline(k_table, f3mu2v1v2b1_table, nk_table, 3d30, 3d30, f3mu2v1v2b1_d_table)
    end subroutine init_rsd_galvel
    
    subroutine init_rsd_convolution()
      implicit   none
      integer         :: i
      real(kind=dp)   :: ak
      real(kind=dp)     :: a000,a002,a004,a020,a022,a024,a040,a042,a044
      real(kind=dp)     :: a200,a202,a204,a220,a222,a224,a240,a242,a244,a260,a262
      real(kind=dp)     :: a400,a402,a404,a420,a422,a424,a440,a442,a444,a460,a462,a464

      if(feedback > 1)write(*,*)'calculating RSD convolution'
      do i = 1,nq_table 
        ak = q_table(i)

        if (ak > k1loop .and. ak < 0.35_dp) then 

          a000 = P000(ak)
          a002 = P002(ak)
          a004 = 0._dp !P004(ak)
          a020 = P020(ak)  
          a022 = P022(ak)
          a024 = 0._dp !P024(ak)
          a040 = 0._dp !P040(ak) 
          a042 = 0._dp !P042(ak) 
          a044 = 0._dp !P044(ak)

          a200 = P200(ak)
          a202 = P202(ak) 
          a204 = 0._dp !P204(ak) 
          a220 = P220(ak)
          a222 = P222(ak)
          a224 = P224(ak)
          a240 = 0._dp !P240(ak)
          a242 = P242(ak)
          a244 = P244(ak)
          
          a260 = 0._dp !P260(ak)
          a262 = 0._dp !P262(ak)

          a400 = P400(ak)
          a402 = P402(ak)  
          a404 = 0._dp !P404(ak) 
          a420 = P420(ak) 
          a422 = P422(ak)
          a424 = 0._dp !P424(ak) 
          a440 = P440(ak) 
          a442 = P442(ak) 
          a444 = P444(ak)

          a460 = 0._dp !P460(ak)
          a462 = 0._dp !P462(ak)
          a464 = 0._dp !P464(ak)

          write(70,'(2x,20e16.6)')ak,a000,a002,a004,a020,a022,a024,a040,a042,a044
          
          write(72,'(2x,20e16.6)')ak,a200,a202,a204,a220,a222,a224,a240,a242,a244,a260,a262
          
          write(74,'(2x,20e16.6)')ak,a400,a402,a404,a420,a422,a424,a440,a442,a444,a460,a462,a464
          
          write(76,'(2x,20e16.6)')ak,power_gal_ell0_novir_at(ak),power_gal_ell2_novir_at(ak) &
          ,power_gal_ell4_novir_at(ak),power_gal_ell6_novir_at(ak),power_gal_ell8_novir_at(ak)
          
          rsd_convol0_table(i) = a000+a002+a004+a020+a022+a024+a040+a042+a044
          rsd_convol2_table(i) = a200+a202+a204+a220+a222+a224+a240+a242+a244+a260+a262
          rsd_convol4_table(i) = a400+a402+a404+a420+a422+a424+a440+a442+a444+a460+a462+a464

        else
          rsd_convol0_table(i) = 0._dp
          rsd_convol2_table(i) = 0._dp
          rsd_convol4_table(i) = 0._dp
        endif
        
      enddo 

      call spline(q_table, rsd_convol0_table, nq_table, 3d30, 3d30, rsd_convol0_d_table)
      call spline(q_table, rsd_convol2_table, nq_table, 3d30, 3d30, rsd_convol2_d_table)
      call spline(q_table, rsd_convol4_table, nq_table, 3d30, 3d30, rsd_convol4_d_table)
      
    end subroutine init_rsd_convolution

    subroutine tns_check()
      implicit none
      integer        :: i, j
      real(kind=dp)  :: ak, Amu2, Amu4, Amu6, Bmu2, Bmu4, Bmu6, Bmu8

      write(*,*)'outputing A and B to fort.77'
      do i = 1, nk_table
        ak = k_table(i)
        if (ak > k1loop) then 
          Amu2 = f*RSDskew_fmu2v1sq(ak) + f**2 *RSDskew_f2mu2v1v2(ak)
          Amu4 = f**2 *RSDskew_f2mu4v1v2(ak)+ f**3 *RSDskew_f3mu4v2sq(ak)
          Amu6 = f**3 *RSDskew_f3mu6v2sq(ak)

          Bmu2 = f**2 *RSDkurt_f2mu2v1sq(ak) + f**3 *RSDkurt_f3mu2v1v2(ak) + f**4 *RSDkurt_f4mu2v2sq(ak)
          Bmu4 = f**2 *RSDkurt_f2mu4v1sq(ak) + f**3 *RSDkurt_f3mu4v1v2(ak) + f**4 *RSDkurt_f4mu4v2sq(ak)
          Bmu6 = f**3 *RSDkurt_f3mu6v1v2(ak) + f**4 *RSDkurt_f4mu6v2sq(ak)
          Bmu8 = f**4 *RSDkurt_f4mu8v2sq(ak)

          write(77,'(2x,8e16.6)')ak,Amu2,Amu4,Amu6,Bmu2,Bmu4,Bmu6,Bmu8
        endif
      enddo 

    end subroutine tns_check

    subroutine init_VDskew()  
      implicit none  
      integer        :: i
      real(kind=dp)  :: p00,p01,p02,p03,p04,ak,pl
      
      write(*,*)'calculating VelDiff skewness'
      do i=1,nr_table
         Psi1m1_table(i) = Psi1m1i(r_table(i))
      enddo

      call spline(r_table, Psi1m1_table, nr_table, 3d30,3d30, Psi1m1_d_table)

!      do i=1,nr_table
!         write(*,*)r_table(i),Psi1m1(r_table(i))
!      enddo
      do i=1,nk_table
         ak =  k_table(i)
         if (0.001_dp<ak .and. ak<0.35_dp) then 
         p00=pw00(ak) 
         p01=pw01(ak) 
         p02=pw02(ak) 
         p03=pw03(ak) 
         p04=pw04(ak) 
         pl=f**2*matterpowerat(ak)/5.d0
         write(*,'(2x,20e16.6)')ak,pl,p00,p01,p02,p03,p04
!         write(*,'(2x,20e16.6)')ak,power_test0(ak)/matterpowerat(ak) &
!                                  ,power_test2(ak)/matterpowerat(ak)
         endif               
      enddo   
    end subroutine init_VDskew

!************************************************************
!  calculate veldiff skewness
!************************************************************

   function Psi1m1i(arg_r)
     implicit none
     real(kind=dp)   :: Psi1m1i,arg_r
      integer, parameter :: limit=10000
      integer            :: neval, ier, iord(limit), last
      real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
      real(kind=dp)      :: epsabs, epsrel, abserr, ans
 
      rpsi=arg_r
      epsabs = 0._dp
      epsrel = 0.0001_dp         !! Relative Error !!
 
      call qage(Psi1m1_dq,kminRSD,kmaxRSD,epsabs,epsrel,30,limit,    &
        ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Psi1m1i   = 4._dp*pi*ans

   end function Psi1m1i
    
   function Psi1m1_dq(arg_q)
     implicit none
     real(kind=dp)   :: Psi1m1_dq,arg_q,x
     
     x = arg_q*rpsi
     Psi1m1_dq   = arg_q*matterpowerat(arg_q)*sphbe1(x)

   end function Psi1m1_dq

   function Psi1m1(arg_r)   
      implicit none
      real(kind=dp)   :: Psi1m1,arg_r

      if (arg_r >= rmin .and. arg_r <= rmax) then 
        call splint(r_table,Psi1m1_table,Psi1m1_d_table,nr_table,arg_r,Psi1m1)
      else 
        Psi1m1 = 0._dp
      endif

    end function Psi1m1

!************************************************************
!  Bias integrals that reqire wiggle exorcism (it turns out, not)
!************************************************************

   function Pk_b1g2(arg_k)
     implicit none
     real(kind=dp)   :: Pk_b1g2,arg_k

      Pk_b1g2   = Pk_b1g2_prop(arg_k)*matterpowerat(arg_k) + Pk_b1g2_mc(arg_k)
!     Pk_b1g2   = Pk_b1g2_prop(arg_k)*gRPT_Pdv(arg_k) + Pk_b1g2_mc(arg_k)

   end function Pk_b1g2
    
   function Pk_b1g3m(arg_k)
     implicit none
     real(kind=dp)   :: Pk_b1g3m,arg_k

     Pk_b1g3m   = Pk_b1g3m_prop(arg_k)*matterpowerat(arg_k)!*gRPT_Pdv(arg_k)

   end function Pk_b1g3m

   function Pk_g2(arg_k)
   implicit none
   real(kind=dp)   :: Pk_g2,arg_k

    Pk_g2   = Pk_g2_prop(arg_k)*matterpowerat(arg_k) + Pk_g2_mc(arg_k)
!    Pk_g2   = Pk_g2_prop(arg_k)*gRPT_Pdv(arg_k) + Pk_g2_mc(arg_k)

   end function Pk_g2

   function Pk_g3m(arg_k)
     implicit none
     real(kind=dp)   :: Pk_g3m,arg_k

     Pk_g3m = Pk_b1g3m(arg_k)/2._dp

   end function Pk_g3m

!************************************************************
!  master 4D integral for generic case 
!************************************************************

    function power_rsd(ell,arg_k)
    implicit none
    real(kind=dp)   :: power_rsd,arg_k
    integer            :: ell,ell1,ell2
    
    power_rsd = 0._dp
!    do ell1=0,8,2
!       do ell2=0,8,2
!          power_rsd = power_rsd + power_gal_novir_at(ell2,arg_q)
!          *(2._dp*ell+1._dp)  *LegendreP(ell,mu)  *LegendreP(ell1,mu)
!          *(2._dp*ell1+1._dp) *LegendreP(ell1,nu) *LegendreP(ell2,nu)  
!          *(2._dp/pi) *sphbessel(ell1,arg_k*ar) *sphbessel(ell2,arg_q*ar) 
!          *( Z0(la,ar,nu) - Z0si(la) )
!       enddo   
!    enddo
    end function power_rsd
    
    
    function power_test0(arg_k)
    implicit none
    real(kind=dp)   :: power_test0,arg_k
    integer, parameter :: limit=10000
    integer            :: neval, ier, iord(limit), last
    real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
    real(kind=dp)      :: epsabs, epsrel, abserr, ans

    k = arg_k
    
    epsabs = 0._dp
    epsrel = 1.e-7_dp         !! Relative Error !!
 
      call qage(ptest0_dr,rmin,rmax,epsabs,epsrel,6,limit,    &
           ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      power_test0 = k**2 *ans*fac_norm/(2._dp*pi**2)
    
    end function power_test0

    function ptest0_dr(arg_r)
    implicit none
     real(kind=dp)   :: ptest0_dr,arg_r,x

     x = k * arg_r
     ptest0_dr = sphbe0(x) * arg_r**2 *I0(arg_r)
    
    end function ptest0_dr

    function power_test2(arg_k)
    implicit none
    real(kind=dp)   :: power_test2,arg_k
    integer, parameter :: limit=10000
    integer            :: neval, ier, iord(limit), last
    real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
    real(kind=dp)      :: epsabs, epsrel, abserr, ans

    k = arg_k
    
    epsabs = 0._dp
    epsrel = 1.e-7_dp         !! Relative Error !!
 
      call qage(ptest2_dr,rmin,rmax,epsabs,epsrel,6,limit,    &
           ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      power_test2 = -k**2 *ans*fac_norm/(2._dp*pi**2)
    
    end function power_test2

    function ptest2_dr(arg_r)
    implicit none
     real(kind=dp)   :: ptest2_dr,arg_r,x

     x = k * arg_r
     ptest2_dr = sphbe2(x) * arg_r**2 *I2(arg_r)
    
    end function ptest2_dr


    function power_w(ell1,ell2,arg_k)
    implicit none
    real(kind=dp)   :: power_w,arg_k
    integer         :: ell1, ell2
    integer, parameter :: limit=10000
    integer            :: neval, ier, iord(limit), last
    real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
    real(kind=dp)      :: epsabs, epsrel, abserr, ans

    k = arg_k
    ell  = ell1
    ellp = ell2
    
    epsabs = 0._dp
    epsrel = 0.0001_dp         !! Relative Error !!
 
    if (mod(ellp,2).eq.0) then 
      call qage(pwint_dr_even,rmin,rmax,epsabs,epsrel,30,limit,    &
           ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      power_w = ans*(-1._dp)**(real(ellp,kind=dp)/2._dp) 
    else
      call qage(pwint_dr_odd,rmin,rmax,epsabs,epsrel,30,limit,    &
           ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
       power_w = ans*(-1._dp)**((real(ellp-1,kind=dp))/2._dp) 
    endif
    power_w = power_w *(2._dp*real(ell ,kind=dp) + 1._dp) &
                      *(2._dp*real(ellp,kind=dp) + 1._dp) 
    power_w = power_w *fac_norm/(2._dp*pi**2)
    
    end function power_w


    function pwint_dr_odd(arg_r)
    implicit none
     real(kind=dp)   :: pwint_dr_odd,arg_r,epsabs,epsrel,x
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     x = k * arg_r
     rpw = arg_r
     
     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(ThreeLegWeff_odd,method,epsrel,epsabs,ans,error,prob)
     pwint_dr_odd = ans(1) * sphbessel(ellp,x) * arg_r**2
    
    end function pwint_dr_odd
    
   subroutine ThreeLegWeff_odd(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: mu,nu,Weff,I0r,I2r,varmin,y,Wns
 
     mu = x(1) 
     nu = x(2) 
     
     I0r = I0(rpw)
     I2r = I2(rpw)
     varmin = 2._dp*( (sigv_P)**2 + f**2  &
              *( sigv**2-I0r/3._dp-2._dp*LegendreP2(nu)*I2r/3._dp) )
     y=(k*mu)**2/(1._dp+(k*mu*a_roman)**2)
     Wns = dexp(-y*varmin/2._dp)/dsqrt(1._dp+(k*mu*a_roman)**2)
     
     Weff = Wns * dsin(-y*varmin*nu*f*k*mu*Psi1m1(rpw))

     h(1) = LegendreP(mu)*LegendrePp(mu)*LegendrePp(nu)*(Weff-1._dp)
              
   end subroutine ThreeLegWeff_odd


    function pwint_dr_even(arg_r)
    implicit none
     real(kind=dp)   :: pwint_dr_even,arg_r,epsabs,epsrel,x
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     x = k * arg_r
     rpw = arg_r

     epsabs    = absacc/safety/10._dp
     epsrel    = relacc/safety/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     
     call   integrate2D(ThreeLegWeff_even,method,epsrel,epsabs,ans,error,prob)
     pwint_dr_even = ans(1)* arg_r**2 * sphbe2(x) !* sphbessel(ellp,x) !
!     write(*,*)ellp
    
    end function pwint_dr_even
    
   subroutine ThreeLegWeff_even(ndim,x,ncomp,h)
     implicit none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: mu,nu,Weff,I0r,I2r,varmin,y,Wns,Weff0,Wns0,varmin0
 
     mu = x(1) 
     nu = x(2) 
     
     I0r = I0(rpw)
     I2r = I2(rpw)
     varmin  = 2._dp*( sigv_P**2 + f**2  &
               *( sigv**2-I0r/3._dp-2._dp*LegendreP2(nu)*I2r/3._dp) )
     varmin0 = 2._dp*( sigv_P**2 + f**2 *sigv**2)
     
     y = (k*mu)**2 !/(1._dp+(k*mu*a_roman)**2)
      
     Wns  = dexp(-y*varmin /2._dp) !/dsqrt(1._dp+(k*mu*a_roman)**2)
     
     Wns0 = dexp(-y*varmin0/2._dp) !/dsqrt(1._dp+(k*mu*a_roman)**2)
     
     Weff = Wns !* dcos(-y*varmin*nu*f*k*mu*Psi1m1(rpw))

     Weff0 = Wns0 
     
     Weff = -y*varmin /2._dp

!     Weff = y*2._dp*LegendreP2(nu)/3._dp/k**2

     h(1) = LegendreP(mu)*LegendrePp(mu)*LegendrePp(nu)*(Weff-Weff0)  
              
   end subroutine ThreeLegWeff_even




    function pw00(arg_k)
    implicit none
    real(kind=dp)   :: pw00,arg_k
    integer, parameter :: limit=10000
    integer            :: neval, ier, iord(limit), last
    real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
    real(kind=dp)      :: epsabs, epsrel, abserr, ans

    k = arg_k    
    epsabs = 0._dp
    epsrel = 0.0001_dp         !! Relative Error !!
    
    if (k.gt.k1loop) then 
       call qage(pw00dr,rmin,rmax,epsabs,epsrel,30,limit,    &
           ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
       ans = ans*fac_norm/(2._dp*pi**2)
       pw00 = ans                !*(-1)^ell/2 (2ell+1)(2ellp+1)
    else
       pw00 = f**2*matterpowerat(k)/9._dp
    endif   
    end function pw00

    function pw01(arg_k)
    implicit none
    real(kind=dp)   :: pw01,arg_k
    integer, parameter :: limit=10000
    integer            :: neval, ier, iord(limit), last
    real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
    real(kind=dp)      :: epsabs, epsrel, abserr, ans

    k = arg_k    
    epsabs = 0._dp
    epsrel = 0.0001_dp         !! Relative Error !!
 
    if (k.gt.k1loop) then 
       call qage(pw01dr,rmin,rmax,epsabs,epsrel,30,limit,    &
           ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
       ans = ans*fac_norm/(2._dp*pi**2)
       pw01 = ans*3._dp               !*(-1)^(ell-1)/2 (2ell+1)(2ellp+1)
    else
       pw01 =0._dp
    endif   
    end function pw01

    function pw02(arg_k)
    implicit none
    real(kind=dp)   :: pw02,arg_k
    integer, parameter :: limit=10000
    integer            :: neval, ier, iord(limit), last
    real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
    real(kind=dp)      :: epsabs, epsrel, abserr, ans

    k = arg_k    
    epsabs = 0._dp
    epsrel = 0.0001_dp         !! Relative Error !!
 
    if (k.gt.k1loop) then 
       call qage(pw02dr,rmin,rmax,epsabs,epsrel,30,limit,    &
           ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
       ans = ans*fac_norm/(2._dp*pi**2)
       pw02 = -ans *5._dp                !*(-1)^ell/2 (2ell+1)(2ellp+1)
    else
       pw02 = 4._dp*f**2*matterpowerat(k)/45._dp
    endif   
    end function pw02

    function pw03(arg_k)
    implicit none
    real(kind=dp)   :: pw03,arg_k
    integer, parameter :: limit=10000
    integer            :: neval, ier, iord(limit), last
    real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
    real(kind=dp)      :: epsabs, epsrel, abserr, ans

    k = arg_k    
    epsabs = 0._dp
    epsrel = 0.0001_dp         !! Relative Error !!
 
    if (k.gt.k1loop) then 
       call qage(pw03dr,rmin,rmax,epsabs,epsrel,30,limit,    &
           ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
       ans = ans*fac_norm/(2._dp*pi**2)
       pw03 = -ans*7._dp               !*(-1)^(ell-1)/2 (2ell+1)(2ellp+1)
    else
       pw03 =0._dp
    endif       
    end function pw03

    function pw04(arg_k)
    implicit none
    real(kind=dp)   :: pw04,arg_k
    integer, parameter :: limit=10000
    integer            :: neval, ier, iord(limit), last
    real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
    real(kind=dp)      :: epsabs, epsrel, abserr, ans

    k = arg_k    
    epsabs = 0._dp
    epsrel = 0.0001_dp         !! Relative Error !!
 
    if (k.gt.k1loop) then 
       call qage(pw04dr,rmin,rmax,epsabs,epsrel,30,limit,    &
           ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
       ans = ans*fac_norm/(2._dp*pi**2)
       pw04 = ans*9._dp                !*(-1)^ell/2 (2ell+1)(2ellp+1)
    else
       pw04 = 0._dp
    endif   
    end function pw04



    function pw00dr(arg_r)
    implicit none
     real(kind=dp)   :: pw00dr,arg_r,epsabs,epsrel,x
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     x = k * arg_r
     rpw = arg_r

     epsabs    = absacc/safety/10._dp
     epsrel    = relacc/safety/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     
     call   integrate2D(ThreeLegWeff00,method,epsrel,epsabs,ans,error,prob)
     pw00dr = ans(1)* arg_r**2 * sphbe0(x) !* sphbessel(ellp,x) !
    
    end function pw00dr
    
   subroutine ThreeLegWeff00(ndim,x,ncomp,h)
     implicit none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: mu,nu,Weff,I0r,I2r,varmin,y,Wns,Wns0,varmin0,y0,a
 
     mu = x(1) 
     nu = x(2) 
     
     I0r = I0(rpw)
     I2r = I2(rpw)
     varmin  = 2._dp*( sigv_P**2 + f**2  &
               *( sigv**2-I0r/3._dp-2._dp*LegendreP2(nu)*I2r/3._dp) )
     varmin0 = 2._dp*( sigv_P**2 + f**2 *sigv**2)
     
     a = dsqrt(a_roman**2 + 0._dp*(f*nu*Psi1m1(rpw))**2)
     
     y = (k*mu)**2/(1._dp+(k*mu*a)**2)
     y0 = (k*mu)**2/(1._dp+(k*mu*a_roman)**2)
      
     Wns  = dexp(-y*varmin /2._dp)/dsqrt(1._dp+(k*mu*a)**2)
     
     Wns0 = dexp(-y0*varmin0/2._dp)/dsqrt(1._dp+(k*mu*a_roman)**2)
     
     Weff = Wns !* dcos(-y*varmin*nu*f*k*mu*Psi1m1(rpw))
     
     h(1) = (Weff-Wns0) !Pell(mu)*Pellp(mu)*Pellp(nu) 
              
   end subroutine ThreeLegWeff00

    function pw01dr(arg_r)
    implicit none
     real(kind=dp)   :: pw01dr,arg_r,epsabs,epsrel,x
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     x = k * arg_r
     rpw = arg_r

     epsabs    = absacc/safety/10._dp
     epsrel    = relacc/safety/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     
     call   integrate2D(ThreeLegWeff01,method,epsrel,epsabs,ans,error,prob)
     pw01dr = ans(1)* arg_r**2 * sphbe1(x) !* sphbessel(ellp,x) !
    
    end function pw01dr
    
   subroutine ThreeLegWeff01(ndim,x,ncomp,h)
     implicit none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: mu,nu,Weff,I0r,I2r,varmin,y,Wns,Wns0,varmin0
 
     mu = x(1)
     nu = x(2)
     
     I0r = I0(rpw)
     I2r = I2(rpw)
     varmin  = 2._dp*( sigv_P**2 + f**2  &
               *( sigv**2-I0r/3._dp-2._dp*LegendreP2(nu)*I2r/3._dp) )
     varmin0 = 2._dp*( sigv_P**2 + f**2 *sigv**2)
     
     y = (k*mu)**2/(1._dp+(k*mu*a_roman)**2)
      
     Wns  = dexp(-y*varmin /2._dp)/dsqrt(1._dp+(k*mu*a_roman)**2)
     
     Wns0 = dexp(-y*varmin0/2._dp)/dsqrt(1._dp+(k*mu*a_roman)**2)
     
     Weff = Wns * dsin(-y*varmin*nu*f*k*mu*Psi1m1(rpw))
     
     h(1) = (Weff)*LegendreP1(mu)*LegendreP1(nu) !Pell(mu)*Pellp(mu)*Pellp(nu) 
              
   end subroutine ThreeLegWeff01

    function pw02dr(arg_r)
    implicit none
     real(kind=dp)   :: pw02dr,arg_r,epsabs,epsrel,x
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     x = k * arg_r
     rpw = arg_r

     epsabs    = absacc/safety/10._dp
     epsrel    = relacc/safety/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     
     call   integrate2D(ThreeLegWeff02,method,epsrel,epsabs,ans,error,prob)
     pw02dr = ans(1)* arg_r**2 * sphbe2(x) !* sphbessel(ellp,x) !
    
    end function pw02dr
    
   subroutine ThreeLegWeff02(ndim,x,ncomp,h)
     implicit none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: mu,nu,Weff,I0r,I2r,varmin,y,Wns,Wns0,varmin0
 
     mu = x(1) 
     nu = x(2) 
     
     I0r = I0(rpw)
     I2r = I2(rpw)
     varmin  = 2._dp*( sigv_P**2 + f**2  &
               *( sigv**2-I0r/3._dp-2._dp*LegendreP2(nu)*I2r/3._dp) )
     varmin0 = 2._dp*( sigv_P**2 + f**2 *sigv**2)
     
     y = (k*mu)**2/(1._dp+(k*mu*a_roman)**2)
      
     Wns  = dexp(-y*varmin /2._dp)/dsqrt(1._dp+(k*mu*a_roman)**2)
     
     Wns0 = dexp(-y*varmin0/2._dp)/dsqrt(1._dp+(k*mu*a_roman)**2)
     
     Weff = Wns !* dcos(-y*varmin*nu*f*k*mu*Psi1m1(rpw))
     
     h(1) = (Weff-Wns0)*LegendreP2(mu)*LegendreP2(nu) !Pell(mu)*Pellp(mu)*Pellp(nu) 
              
   end subroutine ThreeLegWeff02

    function pw03dr(arg_r)
    implicit none
     real(kind=dp)   :: pw03dr,arg_r,epsabs,epsrel,x
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     x = k * arg_r
     rpw = arg_r

     epsabs    = absacc/safety/10._dp
     epsrel    = relacc/safety/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     
     call   integrate2D(ThreeLegWeff03,method,epsrel,epsabs,ans,error,prob)
     pw03dr = ans(1)* arg_r**2 * sphbe3(x) !* sphbessel(ellp,x) !
    
    end function pw03dr
    
   subroutine ThreeLegWeff03(ndim,x,ncomp,h)
     implicit none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: mu,nu,Weff,I0r,I2r,varmin,y,Wns,Wns0,varmin0
 
     mu = x(1) 
     nu = x(2) 
     
     I0r = I0(rpw)
     I2r = I2(rpw)
     varmin  = 2._dp*( sigv_P**2 + f**2  &
               *( sigv**2-I0r/3._dp-2._dp*LegendreP2(nu)*I2r/3._dp) )
     varmin0 = 2._dp*( sigv_P**2 + f**2 *sigv**2)
     
     y = (k*mu)**2/(1._dp+(k*mu*a_roman)**2)
      
     Wns  = dexp(-y*varmin /2._dp)/dsqrt(1._dp+(k*mu*a_roman)**2)
     
     Wns0 = dexp(-y*varmin0/2._dp)/dsqrt(1._dp+(k*mu*a_roman)**2)
     
     Weff = Wns * dsin(-y*varmin*nu*f*k*mu*Psi1m1(rpw))
     
     h(1) = (Weff)*LegendreP3(mu)*LegendreP3(nu) !Pell(mu)*Pellp(mu)*Pellp(nu) 
              
   end subroutine ThreeLegWeff03


    function pw04dr(arg_r)
    implicit none
     real(kind=dp)   :: pw04dr,arg_r,epsabs,epsrel,x
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     x = k * arg_r
     rpw = arg_r

     epsabs    = absacc/safety/10._dp
     epsrel    = relacc/safety/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     
     call   integrate2D(ThreeLegWeff04,method,epsrel,epsabs,ans,error,prob)
     pw04dr = ans(1)* arg_r**2 * sphbe4(x) !* sphbessel(ellp,x) !
    
    end function pw04dr
    
   subroutine ThreeLegWeff04(ndim,x,ncomp,h)
     implicit none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: mu,nu,Weff,I0r,I2r,varmin,y,Wns,Wns0,varmin0
 
     mu = x(1) 
     nu = x(2) 
     
     I0r = I0(rpw)
     I2r = I2(rpw)
     varmin  = 2._dp*( sigv_P**2 + f**2  &
               *( sigv**2-I0r/3._dp-2._dp*LegendreP2(nu)*I2r/3._dp) )
     varmin0 = 2._dp*( sigv_P**2 + f**2 *sigv**2)
     
     y = (k*mu)**2/(1._dp+(k*mu*a_roman)**2)
      
     Wns  = dexp(-y*varmin /2._dp)/dsqrt(1._dp+(k*mu*a_roman)**2)
     
     Wns0 = dexp(-y*varmin0/2._dp)/dsqrt(1._dp+(k*mu*a_roman)**2)
     
     Weff = Wns !* dcos(-y*varmin*nu*f*k*mu*Psi1m1(rpw))
     
     h(1) = (Weff-Wns0)*LegendreP4(mu)*LegendreP4(nu) !Pell(mu)*Pellp(mu)*Pellp(nu) 
              
   end subroutine ThreeLegWeff04




!************************************************************
!  RSD convolution kernels that take into account scale-dependent dispersion
!************************************************************

    function RSD_kernel00(arg_k,arg_q)
      implicit none
      real(kind=dp)   arg_k, arg_q, RSD_kernel00
      integer, parameter :: limit=10000
      integer            :: neval, ier, iord(limit), last
      real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
      real(kind=dp)      :: epsabs, epsrel, abserr, ans
 
      epsabs = 0._dp
      epsrel = 0.0001_dp         !! Relative Error !!
      k=arg_k
      q=arg_q
 
      call qage(RSD_K00_int,rmin,rmax,epsabs,epsrel,30,limit,    &
        ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      RSD_kernel00 = (2._dp/pi)*ans
!      write(*,*)k,q,RSD_kernel00
      
    end function RSD_kernel00

    function RSD_K00_int(arg_r)   
        implicit none
        real(kind=dp)   :: arg_r, psiperp, psipar, psieff, RSD_K00_int
        real(kind=dp)   :: I2r, I0r, vdisp_r, vdisp_infty
 
        I0r = I0(arg_r)
!        I2r = I2(arg_r)
!        psiperp = (I0r-I2r)/3._dp 
!        psipar = (I0r+2._dp*I2r)/3._dp 
        psieff = I0r/3._dp 
!        vdisp_infty = (sigv_P*k)**2  
!        vdisp_r = (sigv_P*k)**2  - (f*k)**2 *psieff
! 
!        RSD_K00_int = (q*arg_r)**2 *sphbe0(k*arg_r)*sphbe0(q*arg_r) &
!                    * (Hkernel00(vdisp_r)-Hkernel00(vdisp_infty))
 
        vdisp_infty = (sigv_P)**2  
        vdisp_r = (sigv_P)**2  - (f)**2 *psieff
 
        RSD_K00_int = (q*arg_r)**2 *sphbe0(k*arg_r)*sphbe0(q*arg_r) &
                    * (RSD_virial4(k,vdisp_r,0,0)-RSD_virial4(k,vdisp_infty,0,0))
 
      end function RSD_K00_int
      

    function RSD_kernel02(arg_k,arg_q)
      implicit none
      real(kind=dp)   arg_k, arg_q, RSD_kernel02
      integer, parameter :: limit=10000
      integer            :: neval, ier, iord(limit), last
      real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
      real(kind=dp)      :: epsabs, epsrel, abserr, ans
 
      epsabs = 0._dp
      epsrel = 0.0001_dp         !! Relative Error !!
      k=arg_k
      q=arg_q
 
      call qage(RSD_K02_int,rmin,rmax,epsabs,epsrel,30,limit,    &
        ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      RSD_kernel02 = (2._dp/pi)*ans
!      write(*,*)k,q,RSD_kernel02
      
    end function RSD_kernel02

    function RSD_K02_int(arg_r)   
        implicit none
        real(kind=dp)   :: arg_r, psiperp, psipar, psieff, RSD_K02_int
        real(kind=dp)   :: I2r, I0r, vdisp_r, vdisp_infty
 
        I0r = I0(arg_r)
!        I2r = I2(arg_r)
!        psiperp = (I0r-I2r)/3._dp 
!        psipar = (I0r+2._dp*I2r)/3._dp 
        psieff = I0r/3._dp 
!        vdisp_infty = (sigv_P*k)**2  
!        vdisp_r = (sigv_P*k)**2  - (f*k)**2 *psieff
! 
!        RSD_K02_int = (q*arg_r)**2 *sphbe2(k*arg_r)*sphbe2(q*arg_r) &
!                    * (Hkernel02(vdisp_r)-Hkernel02(vdisp_infty))
 
        vdisp_infty = (sigv_P)**2  
        vdisp_r = (sigv_P)**2  - (f)**2 *psieff
 
        RSD_K02_int = (q*arg_r)**2 *sphbe2(k*arg_r)*sphbe2(q*arg_r) &
                    * (RSD_virial4(k,vdisp_r,0,2)-RSD_virial4(k,vdisp_infty,0,2))

      end function RSD_K02_int
      

    function RSD_kernel20(arg_k,arg_q)
      implicit none
      real(kind=dp)   arg_k, arg_q, RSD_kernel20
      integer, parameter :: limit=10000
      integer            :: neval, ier, iord(limit), last
      real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
      real(kind=dp)      :: epsabs, epsrel, abserr, ans
 
      epsabs = 0._dp
      epsrel = 0.0001_dp         !! Relative Error !!
      k=arg_k
      q=arg_q
 
      call qage(RSD_K20_int,rmin,rmax,epsabs,epsrel,30,limit,    &
        ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      RSD_kernel20 = (2._dp/pi)*ans *5._dp !(2*ell+1)
!      write(*,*)k,q,RSD_kernel20
      
    end function RSD_kernel20

    function RSD_K20_int(arg_r)   
        implicit none
        real(kind=dp)   :: arg_r, psiperp, psipar, psieff, RSD_K20_int
        real(kind=dp)   :: I2r, I0r, vdisp_r, vdisp_infty
 
        I0r = I0(arg_r)
!        I2r = I2(arg_r)
!        psiperp = (I0r-I2r)/3._dp 
!        psipar = (I0r+2._dp*I2r)/3._dp 
        psieff = I0r/3._dp 
!        vdisp_infty = (sigv_P*k)**2  
!        vdisp_r = (sigv_P*k)**2  - (f*k)**2 *psieff
 
!        RSD_K20_int = (q*arg_r)**2 *sphbe0(k*arg_r)*sphbe0(q*arg_r) &
!                    * (Hkernel02(vdisp_r)-Hkernel02(vdisp_infty))
 
        vdisp_infty = (sigv_P)**2  
        vdisp_r = (sigv_P)**2  - (f)**2 *psieff
 
        RSD_K20_int = (q*arg_r)**2 *sphbe0(k*arg_r)*sphbe0(q*arg_r) &
                    * (RSD_virial4(k,vdisp_r,2,0)-RSD_virial4(k,vdisp_infty,2,0))

      end function RSD_K20_int
          

    function RSD_kernel22(arg_k,arg_q)
      implicit none
      real(kind=dp)   arg_k, arg_q, RSD_kernel22
      integer, parameter :: limit=10000
      integer            :: neval, ier, iord(limit), last
      real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
      real(kind=dp)      :: epsabs, epsrel, abserr, ans
 
      epsabs = 0._dp
      epsrel = 0.0001_dp         !! Relative Error !!
      k=arg_k
      q=arg_q
 
      call qage(RSD_K22_int,rmin,rmax,epsabs,epsrel,30,limit,    &
        ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      RSD_kernel22 = (2._dp/pi)*ans *5._dp !(2ell+1)
!      write(*,*)k,q,RSD_kernel22
      
    end function RSD_kernel22

    function RSD_K22_int(arg_r)   
        implicit none
        real(kind=dp)   :: arg_r, psiperp, psipar, psieff, RSD_K22_int
        real(kind=dp)   :: I2r, I0r, vdisp_r, vdisp_infty
 
        I0r = I0(arg_r)
!        I2r = I2(arg_r)
!        psiperp = (I0r-I2r)/3._dp 
!        psipar = (I0r+2._dp*I2r)/3._dp 
        psieff = I0r/3._dp 
!        vdisp_infty = (sigv_P*k)**2  
!        vdisp_r = (sigv_P*k)**2  - (f*k)**2 *psieff
 
!        RSD_K22_int = (q*arg_r)**2 *sphbe2(k*arg_r)*sphbe2(q*arg_r) &
!                    * (Hkernel22(vdisp_r)-Hkernel22(vdisp_infty))
 
        vdisp_infty = (sigv_P)**2  
        vdisp_r = (sigv_P)**2  - (f)**2 *psieff
 
        RSD_K22_int = (q*arg_r)**2 *sphbe2(k*arg_r)*sphbe2(q*arg_r) &
                    * (RSD_virial4(k,vdisp_r,2,2)-RSD_virial4(k,vdisp_infty,2,2))

      end function RSD_K22_int
      
      
    function RSD_kernel40(arg_k,arg_q)
      implicit none
      real(kind=dp)   arg_k, arg_q, RSD_kernel40
      integer, parameter :: limit=10000
      integer            :: neval, ier, iord(limit), last
      real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
      real(kind=dp)      :: epsabs, epsrel, abserr, ans
 
      epsabs = 0._dp
      epsrel = 0.0001_dp         !! Relative Error !!
      k=arg_k
      q=arg_q
 
      call qage(RSD_K40_int,rmin,rmax,epsabs,epsrel,30,limit,    &
        ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      RSD_kernel40 = (2._dp/pi)*ans *9._dp !(2ell+1)
!      write(*,*)k,q,RSD_kernel40
      
    end function RSD_kernel40

    function RSD_K40_int(arg_r)   
        implicit none
        real(kind=dp)   :: arg_r, psiperp, psipar, psieff, RSD_K40_int
        real(kind=dp)   :: I2r, I0r, vdisp_r, vdisp_infty
 
        I0r = I0(arg_r)
!        I2r = I2(arg_r)
!        psiperp = (I0r-I2r)/3._dp 
!        psipar = (I0r+2._dp*I2r)/3._dp 
        psieff = I0r/3._dp 
!        vdisp_infty = (sigv_P*k)**2  
!        vdisp_r = (sigv_P*k)**2  - (f*k)**2 *psieff
 
!        RSD_K40_int = (q*arg_r)**2 *sphbe0(k*arg_r)*sphbe0(q*arg_r) &
!                    * (Hkernel04(vdisp_r)-Hkernel04(vdisp_infty))
 
        vdisp_infty = (sigv_P)**2  
        vdisp_r = (sigv_P)**2  - (f)**2 *psieff
 
        RSD_K40_int = (q*arg_r)**2 *sphbe0(k*arg_r)*sphbe0(q*arg_r) &
                    * (RSD_virial4(k,vdisp_r,4,0)-RSD_virial4(k,vdisp_infty,4,0))

      end function RSD_K40_int
      

    function RSD_kernel42(arg_k,arg_q)
      implicit none
      real(kind=dp)   arg_k, arg_q, RSD_kernel42
      integer, parameter :: limit=10000
      integer            :: neval, ier, iord(limit), last
      real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
      real(kind=dp)      :: epsabs, epsrel, abserr, ans
 
      epsabs = 0._dp
      epsrel = 0.0001_dp         !! Relative Error !!
      k=arg_k
      q=arg_q
 
      call qage(RSD_K42_int,rmin,rmax,epsabs,epsrel,30,limit,    &
        ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      RSD_kernel42 = (2._dp/pi)*ans *9._dp !(2ell+1)
!      write(*,*)k,q,RSD_kernel42
      
    end function RSD_kernel42

    function RSD_K42_int(arg_r)   
        implicit none
        real(kind=dp)   :: arg_r, psiperp, psipar, psieff, RSD_K42_int
        real(kind=dp)   :: I2r, I0r, vdisp_r, vdisp_infty
 
        I0r = I0(arg_r)
!        I2r = I2(arg_r)
!        psiperp = (I0r-I2r)/3._dp 
!        psipar = (I0r+2._dp*I2r)/3._dp 
        psieff = I0r/3._dp 
!        vdisp_infty = (sigv_P*k)**2  
!        vdisp_r = (sigv_P*k)**2  - (f*k)**2 *psieff
 
!        RSD_K42_int = (q*arg_r)**2 *sphbe2(k*arg_r)*sphbe2(q*arg_r) &
!                    * (Hkernel24(vdisp_r)-Hkernel24(vdisp_infty))
 
        vdisp_infty = (sigv_P)**2  
        vdisp_r = (sigv_P)**2  - (f)**2 *psieff
 
        RSD_K42_int = (q*arg_r)**2 *sphbe2(k*arg_r)*sphbe2(q*arg_r) &
                    * (RSD_virial4(k,vdisp_r,4,2)-RSD_virial4(k,vdisp_infty,4,2))

      end function RSD_K42_int


    function RSD_kernel44(arg_k,arg_q)
      implicit none
      real(kind=dp)   arg_k, arg_q, RSD_kernel44
      integer, parameter :: limit=10000
      integer            :: neval, ier, iord(limit), last
      real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
      real(kind=dp)      :: epsabs, epsrel, abserr, ans
 
      epsabs = 0._dp
      epsrel = 0.0001_dp         !! Relative Error !!
      k=arg_k
      q=arg_q
 
      call qage(RSD_K44_int,rmin,rmax,epsabs,epsrel,30,limit,    &
        ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      RSD_kernel44 = (2._dp/pi)*ans *9._dp !(2ell+1)
!      write(*,*)k,q,RSD_kernel44
      
    end function RSD_kernel44

    function RSD_K44_int(arg_r)   
        implicit none
        real(kind=dp)   :: arg_r, psiperp, psipar, psieff, RSD_K44_int
        real(kind=dp)   :: I2r, I0r, vdisp_r, vdisp_infty
 
        I0r = I0(arg_r)
!        I2r = I2(arg_r)
!        psiperp = (I0r-I2r)/3._dp 
!        psipar = (I0r+2._dp*I2r)/3._dp 
        psieff = I0r/3._dp 
!        vdisp_infty = (sigv_P*k)**2  
!        vdisp_r = (sigv_P*k)**2  - (f*k)**2 *psieff
 
!        RSD_K44_int = (q*arg_r)**2 *sphbe4(k*arg_r)*sphbe4(q*arg_r) &
!                    * (Hkernel44(vdisp_r)-Hkernel44(vdisp_infty))
 
 
        vdisp_infty = (sigv_P)**2  
        vdisp_r = (sigv_P)**2  - (f)**2 *psieff
 
        RSD_K44_int = (q*arg_r)**2 *sphbe4(k*arg_r)*sphbe4(q*arg_r) &
                    * (RSD_virial4(k,vdisp_r,4,4)-RSD_virial4(k,vdisp_infty,4,4))

      end function RSD_K44_int
      
      
!************************************************************
!  2D integrals for scale-dependent and anisotropic dispersion 
!************************************************************

!monopole:
! - aniso
     function P000(arg_k)
     implicit none
     real(kind=dp)   :: P000,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P000_int,method,epsrel,epsabs,ans,error,prob)
     P000 = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P000
   
   subroutine P000_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang000(xsq,ysq,zsq,vsq) *power_gal_ell0_novir_at(q) & !
              *sphbe0(k*r)*sphbe0(q*r) 
              

   end subroutine P000_int

     function P022(arg_k)
     implicit none
     real(kind=dp)   :: P022,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P022_int,method,epsrel,epsabs,ans,error,prob)
     P022 = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P022
   
   subroutine P022_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang022(xsq,ysq,zsq,vsq) *power_gal_ell2_novir_at(q) & !
              *sphbe2(k*r)*sphbe2(q*r) 
              

   end subroutine P022_int

     function P042(arg_k)
     implicit none
     real(kind=dp)   :: P042,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P042_int,method,epsrel,epsabs,ans,error,prob)
     P042 = - ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P042
   
   subroutine P042_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang042(xsq,ysq,zsq) *power_gal_ell2_novir_at(q) & !
              *sphbe4(k*r)*sphbe2(q*r) 
              

   end subroutine P042_int

     function P024(arg_k)
     implicit none
     real(kind=dp)   :: P024,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P024_int,method,epsrel,epsabs,ans,error,prob)
     P024 = - ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P024
   
   subroutine P024_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang024(xsq,ysq,zsq) *power_gal_ell4_novir_at(q) & !
              *sphbe2(k*r)*sphbe4(q*r) 
              

   end subroutine P024_int

     function P044(arg_k)
     implicit none
     real(kind=dp)   :: P044,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P044_int,method,epsrel,epsabs,ans,error,prob)
     P044 = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P044
   
   subroutine P044_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang044(xsq,ysq,zsq,vsq) *power_gal_ell4_novir_at(q) & !
              *sphbe4(k*r)*sphbe4(q*r) 
              

   end subroutine P044_int

     function P002(arg_k)
     implicit none
     real(kind=dp)   :: P002,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P002_int,method,epsrel,epsabs,ans,error,prob)
     P002 = - ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P002
   
   subroutine P002_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang002(xsq,ysq,zsq) *power_gal_ell2_novir_at(q) & !
              *sphbe0(k*r)*sphbe2(q*r) 
              

   end subroutine P002_int

     function P004(arg_k)
     implicit none
     real(kind=dp)   :: P004,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P004_int,method,epsrel,epsabs,ans,error,prob)
     P004 =  ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P004
   
   subroutine P004_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang004(xsq,ysq,zsq) *power_gal_ell4_novir_at(q) & !
              *sphbe0(k*r)*sphbe4(q*r) 
              

   end subroutine P004_int


     function P020(arg_k)
     implicit none
     real(kind=dp)   :: P020,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P020_int,method,epsrel,epsabs,ans,error,prob)
     P020 = - ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P020
   
   subroutine P020_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang020(xsq,ysq,zsq) *power_gal_ell0_novir_at(q) & !
              *sphbe2(k*r)*sphbe0(q*r) 
              

   end subroutine P020_int

     function P040(arg_k)
     implicit none
     real(kind=dp)   :: P040,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P040_int,method,epsrel,epsabs,ans,error,prob)
     P040 =  ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P040
   
   subroutine P040_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang040(xsq,ysq,zsq) *power_gal_ell0_novir_at(q) & !
              *sphbe4(k*r)*sphbe0(q*r) 
              

   end subroutine P040_int


! - iso
     function P000iso(arg_k)
     implicit none
     real(kind=dp)   :: P000iso,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P000iso_int,method,epsrel,epsabs,ans,error,prob)
     P000iso = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P000iso
   
   subroutine P000iso_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     psiperp = I0(r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang000iso(xsq,ysq,vsq) *power_gal_ell0_novir_at(q) & !
              *sphbe0(k*r)*sphbe0(q*r) 
              

   end subroutine P000iso_int

     function P022iso(arg_k)
     implicit none
     real(kind=dp)   :: P022iso,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P022iso_int,method,epsrel,epsabs,ans,error,prob)
     P022iso = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P022iso
   
   subroutine P022iso_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     psiperp = I0(r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang022iso(xsq,ysq,vsq) *power_gal_ell2_novir_at(q) & !
              *sphbe2(k*r)*sphbe2(q*r) 
              

   end subroutine P022iso_int

     function P044iso(arg_k)
     implicit none
     real(kind=dp)   :: P044iso,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P044iso_int,method,epsrel,epsabs,ans,error,prob)
     P044iso = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P044iso
   
   subroutine P044iso_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     psiperp = I0(r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang044iso(xsq,ysq,vsq) *power_gal_ell4_novir_at(q) & !
              *sphbe4(k*r)*sphbe4(q*r) 
              

   end subroutine P044iso_int

!quadrupole:
! - aniso

   function P200(arg_k)
     implicit none
     real(kind=dp)   :: P200,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P200_int,method,epsrel,epsabs,ans,error,prob)
     P200 = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P200
   
   subroutine P200_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang200(xsq,ysq,zsq,vsq) *power_gal_ell0_novir_at(q) & !
              *sphbe0(k*r)*sphbe0(q*r) 

   end subroutine P200_int

     function P222(arg_k)
     implicit none
     real(kind=dp)   :: P222,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P222_int,method,epsrel,epsabs,ans,error,prob)
     P222 = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P222
   
   subroutine P222_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang222(xsq,ysq,zsq,vsq) *power_gal_ell2_novir_at(q) & !
              *sphbe2(k*r)*sphbe2(q*r) 
              

   end subroutine P222_int

   function P244(arg_k)
     implicit none
     real(kind=dp)   :: P244,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P244_int,method,epsrel,epsabs,ans,error,prob)
     P244 = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P244
   
   subroutine P244_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang244(xsq,ysq,zsq,vsq) *power_gal_ell4_novir_at(q) & !
              *sphbe4(k*r)*sphbe4(q*r) 
              

   end subroutine P244_int


   function P202(arg_k)
     implicit none
     real(kind=dp)   :: P202,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P202_int,method,epsrel,epsabs,ans,error,prob)
     P202 = - ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P202
   
   subroutine P202_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang202(xsq,ysq,zsq) *power_gal_ell2_novir_at(q) & !
              *sphbe0(k*r)*sphbe2(q*r) 
              

   end subroutine P202_int


     function P220(arg_k)
     implicit none
     real(kind=dp)   :: P220,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P220_int,method,epsrel,epsabs,ans,error,prob)
     P220 = - ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P220
   
   subroutine P220_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang220(xsq,ysq,zsq) *power_gal_ell0_novir_at(q) & !
              *sphbe2(k*r)*sphbe0(q*r) 
              

   end subroutine P220_int


     function P242(arg_k)
     implicit none
     real(kind=dp)   :: P242,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P242_int,method,epsrel,epsabs,ans,error,prob)
     P242 = - ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P242
   
   subroutine P242_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang242(xsq,ysq,zsq) *power_gal_ell2_novir_at(q) & !
              *sphbe4(k*r)*sphbe2(q*r) 
              

   end subroutine P242_int


     function P224(arg_k)
     implicit none
     real(kind=dp)   :: P224,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P224_int,method,epsrel,epsabs,ans,error,prob)
     P224 = - ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P224
   
   subroutine P224_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang224(xsq,ysq,zsq) *power_gal_ell4_novir_at(q) & !
              *sphbe2(k*r)*sphbe4(q*r) 
              

   end subroutine P224_int


     function P240(arg_k)
     implicit none
     real(kind=dp)   :: P240,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P240_int,method,epsrel,epsabs,ans,error,prob)
     P240 =  ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P240
   
   subroutine P240_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang240(xsq,ysq,zsq) *power_gal_ell0_novir_at(q) & !
              *sphbe4(k*r)*sphbe0(q*r) 
              

   end subroutine P240_int


   function P204(arg_k)
     implicit none
     real(kind=dp)   :: P204,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P204_int,method,epsrel,epsabs,ans,error,prob)
     P204 =  ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P204
   
   subroutine P204_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang204(xsq,ysq,zsq) *power_gal_ell4_novir_at(q) & !
              *sphbe0(k*r)*sphbe4(q*r) 
              

   end subroutine P204_int


   function P260(arg_k)
     implicit none
     real(kind=dp)   :: P260,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P260_int,method,epsrel,epsabs,ans,error,prob)
     P260 = - ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P260
   
   subroutine P260_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang260(xsq,ysq,zsq) *power_gal_ell0_novir_at(q) & !
              *sphbe6(k*r)*sphbe0(q*r) 
              

   end subroutine P260_int

     function P262(arg_k)
     implicit none
     real(kind=dp)   :: P262,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P262_int,method,epsrel,epsabs,ans,error,prob)
     P262 =  ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P262
   
   subroutine P262_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang262(xsq,ysq,zsq) *power_gal_ell2_novir_at(q) & !
              *sphbe6(k*r)*sphbe2(q*r) 
              

   end subroutine P262_int


! - iso
   function P200iso(arg_k)
     implicit none
     real(kind=dp)   :: P200iso,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P200iso_int,method,epsrel,epsabs,ans,error,prob)
     P200iso = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P200iso
   
   subroutine P200iso_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     psiperp = I0(r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang200iso(xsq,ysq,vsq) *power_gal_ell0_novir_at(q) & !
              *sphbe0(k*r)*sphbe0(q*r) 
              

   end subroutine P200iso_int

     function P222iso(arg_k)
     implicit none
     real(kind=dp)   :: P222iso,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P222iso_int,method,epsrel,epsabs,ans,error,prob)
     P222iso = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P222iso
   
   subroutine P222iso_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     psiperp = I0(r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang222iso(xsq,ysq,vsq) *power_gal_ell2_novir_at(q) & !
              *sphbe2(k*r)*sphbe2(q*r) 
              

   end subroutine P222iso_int

     function P244iso(arg_k)
     implicit none
     real(kind=dp)   :: P244iso,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P244iso_int,method,epsrel,epsabs,ans,error,prob)
     P244iso = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P244iso
   
   subroutine P244iso_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     psiperp = I0(r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang244iso(xsq,ysq,vsq) *power_gal_ell4_novir_at(q) & !
              *sphbe4(k*r)*sphbe4(q*r) 
              

   end subroutine P244iso_int


!hexadecapole:
! - aniso

     function P400(arg_k)
     implicit none
     real(kind=dp)   :: P400,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P400_int,method,epsrel,epsabs,ans,error,prob)
     P400 = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
     end function P400
   
   subroutine P400_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang400(xsq,ysq,zsq,vsq) *power_gal_ell0_novir_at(q) & !
              *sphbe0(k*r)*sphbe0(q*r) 
              

   end subroutine P400_int

     function P422(arg_k)
     implicit none
     real(kind=dp)   :: P422,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P422_int,method,epsrel,epsabs,ans,error,prob)
     P422 = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P422
   
   subroutine P422_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang422(xsq,ysq,zsq,vsq) *power_gal_ell2_novir_at(q) & !
              *sphbe2(k*r)*sphbe2(q*r) 
              

   end subroutine P422_int

     function P444(arg_k)
     implicit none
     real(kind=dp)   :: P444,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P444_int,method,epsrel,epsabs,ans,error,prob)
     P444 = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P444
   
   subroutine P444_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang444(xsq,ysq,zsq,vsq) *power_gal_ell4_novir_at(q) & !
              *sphbe4(k*r)*sphbe4(q*r) 
              

   end subroutine P444_int

     function P466(arg_k)
     implicit none
     real(kind=dp)   :: P466,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P466_int,method,epsrel,epsabs,ans,error,prob)
     P466 = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P466
   
   subroutine P466_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang466(xsq,ysq,zsq,vsq) *power_gal_ell6_novir_at(q) & !
              *sphbe6(k*r)*sphbe6(q*r) 
              
   end subroutine P466_int

     function P402(arg_k)
     implicit none
     real(kind=dp)   :: P402,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P402_int,method,epsrel,epsabs,ans,error,prob)
     P402 = - ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P402
   
   subroutine P402_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang402(xsq,ysq,zsq) *power_gal_ell2_novir_at(q) & !
              *sphbe0(k*r)*sphbe2(q*r) 
              

   end subroutine P402_int

     function P420(arg_k)
     implicit none
     real(kind=dp)   :: P420,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P420_int,method,epsrel,epsabs,ans,error,prob)
     P420 = - ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P420
   
   subroutine P420_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang420(xsq,ysq,zsq) *power_gal_ell0_novir_at(q) & !
              *sphbe2(k*r)*sphbe0(q*r) 
              

   end subroutine P420_int



     function P440(arg_k)
     implicit none
     real(kind=dp)   :: P440,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P440_int,method,epsrel,epsabs,ans,error,prob)
     P440 =  ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P440
   
   subroutine P440_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang440(xsq,ysq,zsq) *power_gal_ell0_novir_at(q) & !
              *sphbe4(k*r)*sphbe0(q*r) 
              

   end subroutine P440_int

     function P442(arg_k)
     implicit none
     real(kind=dp)   :: P442,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P442_int,method,epsrel,epsabs,ans,error,prob)
     P442 =  - ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P442
   
   subroutine P442_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang442(xsq,ysq,zsq) *power_gal_ell2_novir_at(q) & !
              *sphbe4(k*r)*sphbe2(q*r) 
              

   end subroutine P442_int


     function P404(arg_k)
     implicit none
     real(kind=dp)   :: P404,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P404_int,method,epsrel,epsabs,ans,error,prob)
     P404 =  ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P404
   
   subroutine P404_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang404(xsq,ysq,zsq) *power_gal_ell4_novir_at(q) & !
              *sphbe0(k*r)*sphbe4(q*r) 
              

   end subroutine P404_int

     function P424(arg_k)
     implicit none
     real(kind=dp)   :: P424,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P424_int,method,epsrel,epsabs,ans,error,prob)
     P424 =  - ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P424
   
   subroutine P424_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang424(xsq,ysq,zsq) *power_gal_ell4_novir_at(q) & !
              *sphbe2(k*r)*sphbe4(q*r) 
              

   end subroutine P424_int


     function P460(arg_k)
     implicit none
     real(kind=dp)   :: P460,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P460_int,method,epsrel,epsabs,ans,error,prob)
     P460 =  - ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P460
   
   subroutine P460_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang460(xsq,ysq,zsq) *power_gal_ell0_novir_at(q) & !
              *sphbe6(k*r)*sphbe0(q*r) 
              

   end subroutine P460_int

     function P462(arg_k)
     implicit none
     real(kind=dp)   :: P462,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P462_int,method,epsrel,epsabs,ans,error,prob)
     P462 =  ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P462
   
   subroutine P462_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang462(xsq,ysq,zsq) *power_gal_ell2_novir_at(q) & !
              *sphbe6(k*r)*sphbe2(q*r) 
              

   end subroutine P462_int

     function P464(arg_k)
     implicit none
     real(kind=dp)   :: P464,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P464_int,method,epsrel,epsabs,ans,error,prob)
     P464 =  - ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P464
   
   subroutine P464_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     I0r = I0(r)
     I2r = I2(r)
     psiperp = (I0r-I2r)/3._dp 
     psipar = (I0r+2._dp*I2r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     zsq = (k*f)**2 *dabs(psiperp-psipar) !take abs value to fix small-r noise
          
     h(1)   = (q*r)**2 *Fang464(xsq,ysq,zsq) *power_gal_ell4_novir_at(q) & !
              *sphbe6(k*r)*sphbe4(q*r) 
              

   end subroutine P464_int




! - iso

     function P400iso(arg_k)
     implicit none
     real(kind=dp)   :: P400iso,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P400iso_int,method,epsrel,epsabs,ans,error,prob)
     P400iso = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
     end function P400iso
   
   subroutine P400iso_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     psiperp = I0(r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang400iso(xsq,ysq,vsq) *power_gal_ell0_novir_at(q) & !
              *sphbe0(k*r)*sphbe0(q*r) 
              

   end subroutine P400iso_int

     function P422iso(arg_k)
     implicit none
     real(kind=dp)   :: P422iso,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P422iso_int,method,epsrel,epsabs,ans,error,prob)
     P422iso = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P422iso
   
   subroutine P422iso_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     psiperp = I0(r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang422iso(xsq,ysq,vsq) *power_gal_ell2_novir_at(q) & !
              *sphbe2(k*r)*sphbe2(q*r) 
              

   end subroutine P422iso_int

     function P444iso(arg_k)
     implicit none
     real(kind=dp)   :: P444iso,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/safety
     epsrel    = relacc/safety
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(P444iso_int,method,epsrel,epsabs,ans,error,prob)
     P444iso = ans(1) *(2._dp/pi) *(kmax-kmin)*(rmax-rmin)
!     write(*,*)ans,error,prob
     
   end function P444iso
   
   subroutine P444iso_int(ndim,x,ncomp,h)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),h(ncomp)
     real(kind=dp)     :: psiperp,psipar,I2r,I0r,q,r
 
     r    = rmax*x(1)+rmin*(1._dp-x(1))
     q    = kmax*x(2)+kmin*(1._dp-x(2))
     psiperp = I0(r)/3._dp 
     
     xsq = k**2 * dabs(sigv_P**2 + f**2 *sigv**2 - f**2 *psiperp) !take abs value to fix small-r noise
     ysq = (k*a_roman)**2
     vsq = k**2 * (sigv_P**2 + f**2 *sigv**2)
          
     h(1)   = (q*r)**2 *Fang444iso(xsq,ysq,vsq) *power_gal_ell4_novir_at(q) & !
              *sphbe4(k*r)*sphbe4(q*r) 
              

   end subroutine P444iso_int

!************************************************************
!  F_ell_ell1_ell2[x^2,y^2,z^2], integral over real and Fourier angles (includes 2ell+1 factor)
!************************************************************

!diagonal terms

      function Fang000(arg_xsq,arg_ysq,arg_zsq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang000
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq !in reality these are not needed since they are already defined above
      ysq = arg_ysq
      zsq = arg_zsq 
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc !/10._dp
      epsabs = absacc !/10._dp
      ans = 0._dp
 
      call qage (Fang000_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang000 = ans   
      if(ier /= 0) call warning('Fang000',ier)
    end function Fang000

      function Fang000_int(mu)        
      implicit none
      real(kind=dp)   :: Fang000_int, mu 
      
      Fang000_int = exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *Hkernel00(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
!                *HkernelNoAniso00(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
                exp(-mu**2 *vsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq)

      end function Fang000_int

      function Fang000iso(arg_xsq,arg_ysq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_vsq,Fang000iso
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq !in reality these are not needed since they are already defined above
      ysq = arg_ysq
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc !/10._dp
      epsabs = absacc !/10._dp
      ans = 0._dp
 
      call qage (Fang000iso_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang000iso = ans   
      if(ier /= 0) call warning('Fang000iso',ier)
    end function Fang000iso

      function Fang000iso_int(mu)        
      implicit none
      real(kind=dp)   :: Fang000iso_int, mu 
      
      Fang000iso_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq)) & 
                     -  exp(-mu**2 *vsq/(1._dp + mu**2 *ysq)) )/dsqrt(1._dp + mu**2 *ysq)

      end function Fang000iso_int


      function Fang022(arg_xsq,arg_ysq,arg_zsq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang022
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang022_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang022 = ans   
      if(ier /= 0) call warning('Fang022',ier)
    end function Fang022

      function Fang022_int(mu)        
      implicit none
      real(kind=dp)   :: Fang022_int, mu 
      
      Fang022_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *5._dp *Hkernel22(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
!                *5._dp *HkernelNoAniso22(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
                exp(-mu**2 *vsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq)) &
                * LegendreP2(mu)

      end function Fang022_int



      function Fang022iso(arg_xsq,arg_ysq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_vsq,Fang022iso
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang022iso_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang022iso = ans   
      if(ier /= 0) call warning('Fang022iso',ier)
    end function Fang022iso

      function Fang022iso_int(mu)        
      implicit none
      real(kind=dp)   :: Fang022iso_int, mu 
      
      Fang022iso_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq)) & 
                     -  exp(-mu**2 *vsq/(1._dp + mu**2 *ysq)) )/dsqrt(1._dp + mu**2 *ysq) &
                     * LegendreP2(mu)

      end function Fang022iso_int


      function Fang044(arg_xsq,arg_ysq,arg_zsq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang044
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang044_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang044 = ans   
      if(ier /= 0) call warning('Fang044',ier)
    end function Fang044

      function Fang044_int(mu)        
      implicit none
      real(kind=dp)   :: Fang044_int, mu 
      
      Fang044_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *9._dp *Hkernel44(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
!                *9._dp *HkernelNoAniso44(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
                exp(-mu**2 *vsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq)) &
                * LegendreP4(mu)

      end function Fang044_int


      function Fang044iso(arg_xsq,arg_ysq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_vsq,Fang044iso
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang044iso_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang044iso = ans   
      if(ier /= 0) call warning('Fang044iso',ier)
    end function Fang044iso

      function Fang044iso_int(mu)        
      implicit none
      real(kind=dp)   :: Fang044iso_int, mu 
      
      Fang044iso_int =  (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq)) & 
                     -  exp(-mu**2 *vsq/(1._dp + mu**2 *ysq)) )/dsqrt(1._dp + mu**2 *ysq) &
                     * LegendreP4(mu)

      end function Fang044iso_int

      function Fang200(arg_xsq,arg_ysq,arg_zsq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang200
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang200_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang200 = ans  *5._dp
      if(ier /= 0) call warning('Fang200',ier)
    end function Fang200

      function Fang200_int(mu)        
      implicit none
      real(kind=dp)   :: Fang200_int, mu 
      
      Fang200_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *Hkernel00(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
!                *HkernelNoAniso00(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
                exp(-mu**2 *vsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq)) &
                * LegendreP2(mu)
                
      end function Fang200_int

      function Fang200iso(arg_xsq,arg_ysq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_vsq,Fang200iso
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang200iso_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang200iso = ans  *5._dp
      if(ier /= 0) call warning('Fang200iso',ier)
    end function Fang200iso

      function Fang200iso_int(mu)        
      implicit none
      real(kind=dp)   :: Fang200iso_int, mu 
      
      Fang200iso_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq)) & 
                     -  exp(-mu**2 *vsq/(1._dp + mu**2 *ysq)) )/dsqrt(1._dp + mu**2 *ysq) &
                * LegendreP2(mu)
                
      end function Fang200iso_int


      function Fang222(arg_xsq,arg_ysq,arg_zsq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang222
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang222_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang222 = ans  *5._dp
      if(ier /= 0) call warning('Fang222',ier)
    end function Fang222

      function Fang222_int(mu)        
      implicit none
      real(kind=dp)   :: Fang222_int, mu 
      
      Fang222_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *5._dp *Hkernel22(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
!                *5._dp *HkernelNoAniso22(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
                exp(-mu**2 *vsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq)) &
                * LegendreP2(mu)**2

      end function Fang222_int

      function Fang222iso(arg_xsq,arg_ysq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang222iso
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang222iso_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang222iso = ans  *5._dp
      if(ier /= 0) call warning('Fang222iso',ier)
    end function Fang222iso

      function Fang222iso_int(mu)        
      implicit none
      real(kind=dp)   :: Fang222iso_int, mu 
      
      Fang222iso_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq)) & 
                     -  exp(-mu**2 *vsq/(1._dp + mu**2 *ysq)) )/dsqrt(1._dp + mu**2 *ysq) &
                * LegendreP2(mu)**2

      end function Fang222iso_int

      function Fang244(arg_xsq,arg_ysq,arg_zsq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang244
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang244_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang244 = ans  *5._dp
      if(ier /= 0) call warning('Fang244',ier)
    end function Fang244

      function Fang244_int(mu)        
      implicit none
      real(kind=dp)   :: Fang244_int, mu 
      
      Fang244_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *9._dp *Hkernel44(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
!                *9._dp *HkernelNoAniso44(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
                exp(-mu**2 *vsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq)) &
                * LegendreP2(mu)*LegendreP4(mu)

      end function Fang244_int

      function Fang244iso(arg_xsq,arg_ysq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_vsq,Fang244iso
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang244iso_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang244iso = ans  *5._dp
      if(ier /= 0) call warning('Fang244iso',ier)
    end function Fang244iso

      function Fang244iso_int(mu)        
      implicit none
      real(kind=dp)   :: Fang244iso_int, mu 
      
      Fang244iso_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq)) & 
                     -  exp(-mu**2 *vsq/(1._dp + mu**2 *ysq)) )/dsqrt(1._dp + mu**2 *ysq) &
                     * LegendreP2(mu)*LegendreP4(mu)

      end function Fang244iso_int

      function Fang400(arg_xsq,arg_ysq,arg_zsq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang400
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang400_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang400 = ans  *9._dp
      if(ier /= 0) call warning('Fang400',ier)
    end function Fang400

      function Fang400_int(mu)        
      implicit none
      real(kind=dp)   :: Fang400_int, mu 
      
      Fang400_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *Hkernel00(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
!                *HkernelNoAniso00(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
                exp(-mu**2 *vsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq)) &
                * LegendreP4(mu)
                
      end function Fang400_int

      function Fang400iso(arg_xsq,arg_ysq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_vsq,Fang400iso
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang400iso_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang400iso = ans  *9._dp
      if(ier /= 0) call warning('Fang400iso',ier)
    end function Fang400iso

      function Fang400iso_int(mu)        
      implicit none
      real(kind=dp)   :: Fang400iso_int, mu 
      
      Fang400iso_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq)) & 
                     -  exp(-mu**2 *vsq/(1._dp + mu**2 *ysq)) )/dsqrt(1._dp + mu**2 *ysq) &
                * LegendreP4(mu)
                
      end function Fang400iso_int

      function Fang422(arg_xsq,arg_ysq,arg_zsq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang422
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang422_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang422 = ans  *9._dp
      if(ier /= 0) call warning('Fang422',ier)
    end function Fang422

      function Fang422_int(mu)        
      implicit none
      real(kind=dp)   :: Fang422_int, mu 
      
      Fang422_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *5._dp *Hkernel22(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
!                *5._dp *HkernelNoAniso22(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
                exp(-mu**2 *vsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq)) &
                * LegendreP2(mu)*LegendreP4(mu)

      end function Fang422_int

      function Fang422iso(arg_xsq,arg_ysq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_vsq,Fang422iso
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang422iso_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang422iso = ans  *9._dp
      if(ier /= 0) call warning('Fang422iso',ier)
    end function Fang422iso

      function Fang422iso_int(mu)        
      implicit none
      real(kind=dp)   :: Fang422iso_int, mu 
      
      Fang422iso_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq)) & 
                     -  exp(-mu**2 *vsq/(1._dp + mu**2 *ysq)) )/dsqrt(1._dp + mu**2 *ysq) &
                * LegendreP2(mu)*LegendreP4(mu)

      end function Fang422iso_int

      function Fang444(arg_xsq,arg_ysq,arg_zsq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang444
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang444_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang444 = ans  *9._dp
      if(ier /= 0) call warning('Fang444',ier)
    end function Fang444

      function Fang444_int(mu)        
      implicit none
      real(kind=dp)   :: Fang444_int, mu 
      
      Fang444_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *9._dp *Hkernel44(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
!                *9._dp *HkernelNoAniso44(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
                exp(-mu**2 *vsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq)) &
                * LegendreP4(mu)**2

      end function Fang444_int

      function Fang444iso(arg_xsq,arg_ysq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_vsq,Fang444iso
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang444iso_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang444iso = ans  *9._dp
      if(ier /= 0) call warning('Fang444iso',ier)
    end function Fang444iso

      function Fang444iso_int(mu)        
      implicit none
      real(kind=dp)   :: Fang444iso_int, mu 
      
      Fang444iso_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq)) & 
                     -  exp(-mu**2 *vsq/(1._dp + mu**2 *ysq)) )/dsqrt(1._dp + mu**2 *ysq) &
                * LegendreP4(mu)**2

      end function Fang444iso_int

      function Fang466(arg_xsq,arg_ysq,arg_zsq,arg_vsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang466
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 
      vsq = arg_vsq

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang466_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang466 = ans  *9._dp
      if(ier /= 0) call warning('Fang466',ier)
    end function Fang466

      function Fang466_int(mu)        
      implicit none
      real(kind=dp)   :: Fang466_int, mu 
      
      Fang466_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *13._dp *Hkernel66(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
!                *13._dp *HkernelNoAniso66(mu**2 *zsq/(1._dp + mu**2 *ysq)) - &
                exp(-mu**2 *vsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq)) &
                * LegendreP6(mu)*LegendreP4(mu)

      end function Fang466_int


!extra-diagonal terms

!ell = 0

      function Fang002(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang002
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang002_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang002 = ans   
      if(ier /= 0) call warning('Fang002',ier)
    end function Fang002

      function Fang002_int(mu)        
      implicit none
      real(kind=dp)   :: Fang002_int, mu 
      
      Fang002_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *1._dp *Hkernel02(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) 

      end function Fang002_int

      function Fang004(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang004
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang004_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang004 = ans   
      if(ier /= 0) call warning('Fang004',ier)
    end function Fang004

      function Fang004_int(mu)        
      implicit none
      real(kind=dp)   :: Fang004_int, mu 
      
      Fang004_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *1._dp *Hkernel04(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) 

      end function Fang004_int


      function Fang020(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang020
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang020_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang020 = ans   
      if(ier /= 0) call warning('Fang020',ier)
    end function Fang020

      function Fang020_int(mu)        
      implicit none
      real(kind=dp)   :: Fang020_int, mu 
      
      Fang020_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *5._dp *Hkernel02(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP2(mu)

      end function Fang020_int


      function Fang040(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang040
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang040_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang040 = ans   
      if(ier /= 0) call warning('Fang040',ier)
    end function Fang040

      function Fang040_int(mu)        
      implicit none
      real(kind=dp)   :: Fang040_int, mu 
      
      Fang040_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *9._dp *Hkernel04(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP4(mu)

      end function Fang040_int

      function Fang024(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang024
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang024_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang024 = ans   
      if(ier /= 0) call warning('Fang024',ier)
    end function Fang024

      function Fang024_int(mu)        
      implicit none
      real(kind=dp)   :: Fang024_int, mu 
      
      Fang024_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *5._dp *Hkernel04(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP2(mu)

      end function Fang024_int

      function Fang042(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang042
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang042_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang042 = ans   
      if(ier /= 0) call warning('Fang042',ier)
    end function Fang042

      function Fang042_int(mu)        
      implicit none
      real(kind=dp)   :: Fang042_int, mu 
      
      Fang042_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *9._dp *Hkernel24(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP4(mu)

      end function Fang042_int



!ell = 2

      function Fang202(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang202
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang202_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang202 = ans  *5._dp
      if(ier /= 0) call warning('Fang202',ier)
    end function Fang202

      function Fang202_int(mu)        
      implicit none
      real(kind=dp)   :: Fang202_int, mu 
      
      Fang202_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *1._dp *Hkernel02(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP2(mu)

      end function Fang202_int


      function Fang220(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang220
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang220_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang220 = ans    *5._dp
      if(ier /= 0) call warning('Fang220',ier)
    end function Fang220

      function Fang220_int(mu)        
      implicit none
      real(kind=dp)   :: Fang220_int, mu 
      
      Fang220_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *5._dp *Hkernel02(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP2(mu)**2

      end function Fang220_int

      function Fang224(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang224
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang224_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang224 = ans    *5._dp
      if(ier /= 0) call warning('Fang224',ier)
    end function Fang224

      function Fang224_int(mu)        
      implicit none
      real(kind=dp)   :: Fang224_int, mu 
      
      Fang224_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *5._dp *Hkernel24(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP2(mu)**2

      end function Fang224_int

      function Fang242(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang242
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang242_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang242 = ans    *5._dp
      if(ier /= 0) call warning('Fang242',ier)
    end function Fang242

      function Fang242_int(mu)        
      implicit none
      real(kind=dp)   :: Fang242_int, mu 
      
      Fang242_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *9._dp *Hkernel24(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP2(mu)*LegendreP4(mu)

      end function Fang242_int


      function Fang240(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang240
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang240_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang240 = ans    *5._dp
      if(ier /= 0) call warning('Fang240',ier)
    end function Fang240

      function Fang240_int(mu)        
      implicit none
      real(kind=dp)   :: Fang240_int, mu 
      
      Fang240_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *9._dp *Hkernel04(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP2(mu)*LegendreP4(mu)

      end function Fang240_int

      function Fang204(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang204
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang204_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang204 = ans  *5._dp
      if(ier /= 0) call warning('Fang204',ier)
    end function Fang204

      function Fang204_int(mu)        
      implicit none
      real(kind=dp)   :: Fang204_int, mu 
      
      Fang204_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *1._dp *Hkernel04(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP2(mu)

      end function Fang204_int

      function Fang260(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang260
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang260_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang260 = ans    *5._dp
      if(ier /= 0) call warning('Fang260',ier)
    end function Fang260

      function Fang260_int(mu)        
      implicit none
      real(kind=dp)   :: Fang260_int, mu 
      
      Fang260_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *13._dp *Hkernel06(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP2(mu)*LegendreP6(mu)

      end function Fang260_int

      function Fang262(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang262
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang262_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang262 = ans    *5._dp
      if(ier /= 0) call warning('Fang262',ier)
    end function Fang262

      function Fang262_int(mu)        
      implicit none
      real(kind=dp)   :: Fang262_int, mu 
      
      Fang262_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *13._dp *Hkernel26(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP2(mu)*LegendreP6(mu)

      end function Fang262_int


!ell = 4

      function Fang402(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang402
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang402_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang402 = ans  *9._dp
      if(ier /= 0) call warning('Fang402',ier)
    end function Fang402

      function Fang402_int(mu)        
      implicit none
      real(kind=dp)   :: Fang402_int, mu 
      
      Fang402_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *1._dp *Hkernel02(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP4(mu)

      end function Fang402_int


      function Fang420(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang420
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang420_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang420 = ans    *9._dp
      if(ier /= 0) call warning('Fang420',ier)
    end function Fang420

      function Fang420_int(mu)        
      implicit none
      real(kind=dp)   :: Fang420_int, mu 
      
      Fang420_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *5._dp *Hkernel02(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP4(mu)*LegendreP2(mu)

      end function Fang420_int


      function Fang440(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang440
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang440_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang440 = ans    *9._dp
      if(ier /= 0) call warning('Fang440',ier)
    end function Fang440

      function Fang440_int(mu)        
      implicit none
      real(kind=dp)   :: Fang440_int, mu 
      
      Fang440_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *9._dp *Hkernel04(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP4(mu)*LegendreP4(mu)

      end function Fang440_int

      function Fang442(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang442
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang442_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang442 = ans    *9._dp
      if(ier /= 0) call warning('Fang442',ier)
    end function Fang442

      function Fang442_int(mu)        
      implicit none
      real(kind=dp)   :: Fang442_int, mu 
      
      Fang442_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *9._dp *Hkernel24(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP4(mu)*LegendreP4(mu)

      end function Fang442_int

      function Fang404(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang404
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang404_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang404 = ans    *9._dp
      if(ier /= 0) call warning('Fang404',ier)
    end function Fang404

      function Fang404_int(mu)        
      implicit none
      real(kind=dp)   :: Fang404_int, mu 
      
      Fang404_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *1._dp *Hkernel04(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP4(mu) 

      end function Fang404_int

      function Fang424(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang424
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang424_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang424 = ans    *9._dp
      if(ier /= 0) call warning('Fang424',ier)
    end function Fang424

      function Fang424_int(mu)        
      implicit none
      real(kind=dp)   :: Fang424_int, mu 
      
      Fang424_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *5._dp *Hkernel24(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP4(mu)*LegendreP2(mu)

      end function Fang424_int

      function Fang460(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang460
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang460_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang460 = ans    *9._dp
      if(ier /= 0) call warning('Fang460',ier)
    end function Fang460

      function Fang460_int(mu)        
      implicit none
      real(kind=dp)   :: Fang460_int, mu 
      
      Fang460_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *13._dp *Hkernel06(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP4(mu)*LegendreP6(mu)

      end function Fang460_int

      function Fang462(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang462
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang462_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang462 = ans    *9._dp
      if(ier /= 0) call warning('Fang462',ier)
    end function Fang462

      function Fang462_int(mu)        
      implicit none
      real(kind=dp)   :: Fang462_int, mu 
      
      Fang462_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *13._dp *Hkernel26(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP4(mu)*LegendreP6(mu)

      end function Fang462_int

      function Fang464(arg_xsq,arg_ysq,arg_zsq) 
      implicit none
      real(kind=dp)   :: arg_xsq,arg_ysq,arg_zsq,arg_vsq,Fang464
      real(kind=dp)   :: epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit) 
      real(kind=dp)   :: ans
 
      xsq = arg_xsq
      ysq = arg_ysq
      zsq = arg_zsq 

      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (Fang464_int,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Fang464 = ans    *9._dp
      if(ier /= 0) call warning('Fang464',ier)
    end function Fang464

      function Fang464_int(mu)        
      implicit none
      real(kind=dp)   :: Fang464_int, mu 
      
      Fang464_int = (exp(-mu**2 *xsq/(1._dp + mu**2 *ysq))/dsqrt(1._dp + mu**2 *ysq) & 
                *13._dp *Hkernel46(mu**2 *zsq/(1._dp + mu**2 *ysq)) ) &
                * LegendreP4(mu)*LegendreP6(mu)

      end function Fang464_int

!************************************************************
! Hkernel_ell1_ell2 is the integral of exp(-z *mu^2) with Legendre @ ell1 and @ ell2
!************************************************************

      function Hkernel00(z) 
      implicit none
      real(kind=dp)   :: z,Hkernel00
      
      if (z.gt.0.1_dp) then
         Hkernel00 = (dsqrt(pi)*erf(dsqrt(z)))/(2._dp*dsqrt(z))
      else
         Hkernel00 = 1._dp - z/3._dp + z**2/10._dp
      endif
      
      end function Hkernel00
      
      function Hkernel02(z)
      implicit none
      real(kind=dp)   :: z,Hkernel02
      
      if (z.gt.0.1_dp) then
         Hkernel02 = -3._dp/(4._dp*dexp(z)*z) + (dsqrt(pi)*(3._dp - 2._dp*z)*erf(dsqrt(z)))/(8._dp*z**1.5_dp)
      else
         Hkernel02 = (-2._dp*z)/15._dp + (2._dp*z**2)/35._dp
      endif
      
      end function Hkernel02
      
      function Hkernel04(z)
      implicit none
      real(kind=dp)   :: z,Hkernel04
      
      if (z.gt.0.1_dp) then
         Hkernel04 = (-5._dp*(21._dp+2._dp*z))/(32._dp*dexp(z)*z**2) + &
         (3._dp*dsqrt(pi)*(35._dp+ 4._dp*(-5._dp+z)*z)*erf(dsqrt(z)))/(64._dp*z**2.5_dp)
      else
         Hkernel04 = (4._dp*z**2)/315._dp
      endif
                  
      end function Hkernel04
      
      function Hkernel06(z)
      implicit none
      real(kind=dp)   :: z,Hkernel06
      
      if (z.gt.0.1_dp) then
         Hkernel06 = (-42._dp*dsqrt(z)*(165._dp + 20._dp*z + 4._dp*z**2) - &
         5*dexp(z)*dsqrt(pi)*(-693._dp+378._dp*z - 84._dp*z**2 + 8._dp*z**3)*erf(dsqrt(z)))/ &
         (256._dp*dexp(z)*z**3.5_dp)
      else
         Hkernel06 = (-8._dp*z**3)/9009._dp
      endif
                  
      end function Hkernel06
      
      function Hkernel22(z)
      implicit none
      real(kind=dp)   :: z,Hkernel22
      
      if (z.gt.0.1_dp) then
         Hkernel22 = ((-3._dp*(9._dp+2._dp*z))/(2._dp*dexp(z)*z**2) + &
                  (dsqrt(pi)*(27._dp+4._dp*(-3._dp+z)*z)*erf(dsqrt(z)))/(4._dp*z**2.5_dp))/8._dp
      else
         Hkernel22 = 0.2_dp - (11._dp*z)/105._dp + (3._dp*z**2)/70._dp
      endif
      
      end function Hkernel22
      
      function Hkernel24(z)
      implicit none
      real(kind=dp)   :: z,Hkernel24
      
      if (z.gt.0.1_dp) then
         Hkernel24 = ((-2._dp*dsqrt(z)*(1575._dp+ 4._dp*z*(75._dp+19._dp*z)))/dexp(z) + &
          3._dp*dsqrt(pi)*(525._dp-250._dp*z+52._dp*z**2-8._dp*z**3)*erf(dsqrt(z)))/(256._dp*z**3.5_dp)
      else
         Hkernel24 = (-4._dp*z)/105._dp + (68._dp*z**2)/3465._dp
      endif
                  
      end function Hkernel24
      
      function Hkernel26(z)
      implicit none
      real(kind=dp)   :: z,Hkernel26
      
      if (z.gt.0.1_dp) then
         Hkernel26 =   (-6._dp*dsqrt(z)*(24255._dp + 2._dp*z*(2205._dp + 574._dp*z + 36._dp*z**2)) + &
         5._dp*dexp(z)*dsqrt(pi)*(14553._dp + 8._dp*z*(-882._dp + z*(189._dp + 2._dp*(-12._dp + z)*z))) &
         *erf(dsqrt(z)))/(1024._dp*dexp(z)*z**4.5_dp)
      else
         Hkernel26 = (4._dp*z**2)/1001._dp - (92._dp*z**3)/45045._dp
      endif
                  
      end function Hkernel26
      
      function Hkernel44(z)
      implicit none
      real(kind=dp)   :: z,Hkernel44
      
      if (z.gt.0.1_dp) then
         Hkernel44 = dexp(-z)*(-5._dp*(21._dp+2._dp*z)*(1225._dp+4._dp*z*(25._dp+11._dp*z))) & 
         /(1024._dp*z**4) +  erf(dsqrt(z))*(3._dp*dsqrt(pi)*(42875._dp+24._dp*z*(-875._dp+z* &
         (185._dp+2._dp*(-10._dp+z)*z))))/(2048._dp*z**4.5_dp)
      else
         Hkernel44 = 1._dp/9._dp - (13._dp*z)/231._dp + (643._dp*z**2)/30030._dp
      endif
                  
      end function Hkernel44
      
      function Hkernel46(z)
      implicit none
      real(kind=dp)   :: z,Hkernel46
      
      if (z.gt.0.1_dp) then
         Hkernel46 = (-2._dp*dsqrt(z)*(7640325._dp + 1323000._dp*z + 352800._dp*z**2 + &
         26880._dp*z**3 + 2288._dp*z**4) - 15._dp*dexp(z)*dsqrt(pi)* &
         (-509355._dp + 251370._dp*z - 55272._dp*z**2 + 6832._dp*z**3 - 496._dp*z**4 + & 
         32._dp*z**5)*erf(dsqrt(z)))/(8192._dp*dexp(z)*z**5.5_dp)
      else
         Hkernel46 = (-10._dp*z)/429._dp + (106._dp*z**2)/9009._dp
      endif
                  
      end function Hkernel46
      
      function Hkernel66(z)
      implicit none
      real(kind=dp)   :: z,Hkernel66
      
      if (z.gt.0.1_dp) then
         Hkernel66 = (-42._dp*dsqrt(z)*(26413695._dp + 2._dp*z* &
         (2255715._dp + 4._dp*z*(158319._dp+2._dp*z*(6147._dp+z*(551._dp+22._dp*z))))) + & 
         5._dp*dexp(z)*dsqrt(pi)*(110937519._dp + 20._dp*z*(-2750517._dp+z*(620487._dp + &
         4._dp*z*(-20538._dp+z*(1701._dp+ 4._dp*(-21._dp + z)*z)))))* &
         erf(dsqrt(z)))/(32768._dp*dexp(z)*z**6.5_dp)
      else
         Hkernel66 = 0.07692307692307693_dp-(83._dp*z)/2145._dp+(71._dp*z**2)/4862._dp
      endif
                  
      end function Hkernel66

!************************************************************
! Hkernel_ell1_ell2 for no anisotropy
!************************************************************

      function HkernelNoAniso00(z) 
      implicit none
      real(kind=dp)   :: z,HkernelNoAniso00
      
         HkernelNoAniso00 = 1._dp

      end function HkernelNoAniso00

      function HkernelNoAniso22(z)
      implicit none
      real(kind=dp)   :: z,HkernelNoAniso22
      
         HkernelNoAniso22 = 0.2_dp
      
      end function HkernelNoAniso22
            
      function HkernelNoAniso44(z)
      implicit none
      real(kind=dp)   :: z,HkernelNoAniso44
      
         HkernelNoAniso44 = 1._dp/9._dp 
         
      end function HkernelNoAniso44
      
      function HkernelNoAniso66(z)
      implicit none
      real(kind=dp)   :: z,HkernelNoAniso66
      
         HkernelNoAniso66 = 1._dp/13._dp
         
      end function HkernelNoAniso66
!************************************************************
!   tabulated functions
!************************************************************

! Kaiser-like terms:
!************************************************************

! first separate the straight *matter* pdd,pdv,pvv

!!$    function gRPT_Pdd(arg_k)
!!$      implicit none
!!$      real(kind=dp)   gRPT_Pdd, pow_dd, arg_k
!!$
!!$      if (arg_k.ge.kmin .and. arg_k.le.kmax) then 
!!$         call splint(k_table,gRPTdd_table,gRPTdd_d_table,nk_table,arg_k,pow_dd)
!!$      else
!!$         pow_dd = 0._dp
!!$      endif
!!$      gRPT_Pdd = pow_dd
!!$    end function

    function gRPT_Pdd(arg_k)
    implicit none
    real(kind=dp)   gRPT_Pdd,  arg_k, P0, root, dP13_dd,P22inv_dd

      if (arg_k.ge.kmin .and. arg_k.le.kmax) then 

         if (arg_k > k1loop) then
         
         call splint(k_table, gRPT_root_table, gRPT_root_d_table, nk_table, arg_k, root)
         call splint(k_table, dP13_dd_table,   dP13_dd_d_table,   nk_table, arg_k, dP13_dd)
         call splint(k_table, P22inv_dd_table, P22inv_dd_d_table, nk_table, arg_k, P22inv_dd)

         P0        = A_fast*MatterPowerAt(arg_k)/fac_norm
         root      = A_fast*root
         dP13_dd   = A_fast*dP13_dd
         P22inv_dd = A_fast**2._dp*P22inv_dd
 
         gRPT_Pdd  = exp(root)*(1._dp + dP13_dd)*(P0*(1._dp - root) + P22inv_dd)
         gRPT_Pdd  = fac_norm*gRPT_Pdd 

         else
 
         gRPT_Pdd  = A_fast*MatterPowerAt(arg_k)

         endif

      else
         gRPT_Pdd  = 0._dp
      endif

    end function

!!$    function gRPT_Pdv(arg_k)
!!$      implicit none
!!$      real(kind=dp)   gRPT_Pdv, pow_dv, arg_k
!!$
!!$      if (arg_k.ge.kmin .and. arg_k.le.kmax) then 
!!$         call splint(k_table,gRPTdv_table,gRPTdv_d_table,nk_table,arg_k,pow_dv)
!!$      else
!!$         pow_dv = 0._dp
!!$      endif
!!$
!!$      gRPT_Pdv = pow_dv
!!$    end function

    function gRPT_Pdv(arg_k)
    implicit none
    real(kind=dp)   gRPT_Pdv,  arg_k, root, P0, dP13_dv,P22inv_dv

      if (arg_k.ge.kmin .and. arg_k.le.kmax) then 

         if (arg_k > k1loop) then
         
         call splint(k_table, gRPT_root_table, gRPT_root_d_table, nk_table, arg_k, root)
         call splint(k_table, dP13_dv_table,   dP13_dv_d_table,   nk_table, arg_k, dP13_dv)
         call splint(k_table, P22inv_dv_table, P22inv_dv_d_table, nk_table, arg_k, P22inv_dv)

         P0        = A_fast*MatterPowerAt(arg_k)/fac_norm
         root      = A_fast*root
         dP13_dv   = A_fast*dP13_dv
         P22inv_dv = A_fast**2._dp*P22inv_dv
 
         gRPT_Pdv  = exp(root)*(1._dp + dP13_dv)*(P0*(1._dp - root) + P22inv_dv)
         gRPT_Pdv  = fac_norm*gRPT_Pdv

         else
 
         gRPT_Pdv  = A_fast*MatterPowerAt(arg_k)

         endif

      else
         gRPT_Pdv  = 0._dp
      endif

    end function gRPT_Pdv

!!$    function gRPT_Pvv(arg_k)
!!$      implicit none
!!$      real(kind=dp)   gRPT_Pvv, pow_vv, arg_k
!!$
!!$      if (arg_k.ge.kmin .and. arg_k.le.kmax) then 
!!$         call splint(k_table,gRPTvv_table,gRPTvv_d_table,nk_table,arg_k,pow_vv)
!!$      else
!!$         pow_vv = 0._dp
!!$      endif
!!$
!!$      gRPT_Pvv = pow_vv
!!$    end function

    function gRPT_Pvv(arg_k)
    implicit none
    real(kind=dp)   gRPT_Pvv,  arg_k, root, P0, dP13_vv,P22inv_vv

      if (arg_k.ge.kmin .and. arg_k.le.kmax) then 
 
        if (arg_k > k1loop) then
         
         call splint(k_table, gRPT_root_table, gRPT_root_d_table, nk_table, arg_k, root)
         call splint(k_table, dP13_vv_table,   dP13_vv_d_table,   nk_table, arg_k, dP13_vv)
         call splint(k_table, P22inv_vv_table, P22inv_vv_d_table, nk_table, arg_k, P22inv_vv)

         P0        = A_fast*MatterPowerAt(arg_k)/fac_norm
         root      = A_fast*root
         dP13_vv   = A_fast*dP13_vv
         P22inv_vv = A_fast**2._dp*P22inv_vv
 
         gRPT_Pvv  = exp(root)*(1._dp + dP13_vv)*(P0*(1._dp - root) + P22inv_vv)
         gRPT_Pvv  = fac_norm*gRPT_Pvv

         else

         gRPT_Pvv  = A_fast*MatterPowerAt(arg_k)

         endif

      else
         gRPT_Pvv  = 0._dp
      endif

    end function gRPT_Pvv

    function gRPT_Pgg(arg_k)   
      implicit none
      real(kind=dp)   gRPT_Pgg, pow_dd, pow_dv, arg_k, b1sq
      real(kind=dp)   pow_b1b2, pow_b1g2, pow_b2b2, pow_b2g2, pow_g2g2, pow_b1g3m

      pow_dd = gRPT_Pdd(arg_k)

      call splint(k_table, b1b2_table, b1b2_d_table, nk_table, arg_k, pow_b1b2)
      call splint(k_table, b1g2_table, b1g2_d_table, nk_table, arg_k, pow_b1g2)
      call splint(k_table, b2b2_table, b2b2_d_table, nk_table, arg_k, pow_b2b2)
      call splint(k_table, b2g2_table, b2g2_d_table, nk_table, arg_k, pow_b2g2)
      call splint(k_table, g2g2_table, g2g2_d_table, nk_table, arg_k, pow_g2g2)
      call splint(k_table, b1g3m_table, b1g3m_d_table, nk_table, arg_k, pow_b1g3m)

      b1sq = b1**2 !*(1._dp + 2._dp * (arg_k*sigv_B)**2 )

      gRPT_Pgg = b1sq*pow_dd + A_fast**2*(b1*b2*pow_b1b2 + b1*gam2*pow_b1g2  &
      + b2*b2*pow_b2b2 + b2*gam2*pow_b2g2 + gam2*gam2*pow_g2g2 + b1*gam3minus*pow_b1g3m)

    end function


    function gRPT_Pgv(arg_k)  
      implicit none
      real(kind=dp)   gRPT_Pgv,pow_dv,arg_k, b1run
      real(kind=dp)   pow_b2,pow_g2,pow_g3m

      pow_dv = gRPT_Pdv(arg_k)

      call splint(k_table, b2_table, b2_d_table, nk_table, arg_k, pow_b2)
      call splint(k_table, g2_table, g2_d_table, nk_table, arg_k, pow_g2)
      call splint(k_table, g3m_table, g3m_d_table, nk_table, arg_k, pow_g3m)

      b1run    = b1 !*(1._dp + (arg_k*sigv_B)**2 )
      gRPT_Pgv = b1run*pow_dv + A_fast**2*(b2*pow_b2 + gam2*pow_g2 + gam3minus*pow_g3m)

    end function


    function gRPT_Pgd(arg_k) 
      implicit none
      real(kind=dp)   gRPT_Pgd,pow_dd,arg_k, b1run
      real(kind=dp)   pow_b1b2,pow_b1g2,pow_b1g3m

      pow_dd = gRPT_Pdd(arg_k)

      call splint(k_table, b1b2_table, b1b2_d_table, nk_table, arg_k, pow_b1b2)
      call splint(k_table, b1g2_table, b1g2_d_table, nk_table, arg_k, pow_b1g2)
      call splint(k_table, b1g3m_table, b1g3m_d_table, nk_table, arg_k, pow_b1g3m)

      b1run = b1 !*(1._dp + (arg_k*sigv_B)**2 )
      gRPT_Pgd = b1run*pow_dd + A_fast**2*0.5_dp*(b2*pow_b1b2 + gam2*pow_b1g2 + gam3minus*pow_b1g3m)
 
    end function


!!$    function gRPT_roots(arg_k)
!!$      implicit none
!!$      real(kind=dp)   gRPT_roots,root,arg_k
!!$
!!$      call splint(k_table,gRPTx_table,gRPTx_d_table,nk_table,arg_k,root)
!!$
!!$      gRPT_roots = root 
!!$
!!$    end function

    function gRPT_roots(arg_k) !this we are not using, so no scaling is needed
    implicit none
    real(kind=dp)   gRPT_roots,root,arg_k

      call splint(k_table,gRPT_root_table,gRPT_root_d_table,nk_table,arg_k,root)

      gRPT_roots = A_fast*root 

    end function
    
    function Plin_DW(arg_k) !desperate attempt
    implicit none
    real(kind=dp)   Plin_DW,root,arg_k,res,Ginv,x
    
!    call splint(q_table, PlinDW_table, PlinDW_d_table, nq_table, arg_k, res)
      
!    x = gRPT_roots(arg_k)
      
!    Ginv = dexp(x)!*(1._dp-x)
    
!    write(*,*)arg_k,x,Ginv
      
!   Plin_DW = matterpowerat(arg_k) *Ginv + res *(1._dp - Ginv)
    
    Plin_DW = gRPT_Pdv(arg_k)/A_fast

!    Plin_DW = matterpowerat(arg_k)
    
    end function Plin_DW
        
    
!    A_TNS or P12RS terms:
!************************************************************
!1.
    function RSDskew_fmu2v1sq(arg_k) !adding up 3 bias sq term that go as f mu^2 v1^2
      implicit none
      real(kind=dp)   RSDskew_fmu2v1sq, arg_k
      real(kind=dp)   pow_b1g2, pow_b1b2, pow_b1sq

      !! RSD squared bias terms !!
      call splint(k_table, fmu2v1sqb1sq_table, fmu2v1sqb1sq_d_table, nk_table, arg_k, pow_b1sq)
      call splint(k_table, fmu2v1sqb1b2_table, fmu2v1sqb1b2_d_table, nk_table, arg_k, pow_b1b2)
      call splint(k_table, fmu2v1sqb1g2_table, fmu2v1sqb1g2_d_table, nk_table, arg_k, pow_b1g2)

      if(tns_full==1)then
        RSDskew_fmu2v1sq = A_fast**2*(b1**2*pow_b1sq + b1*b2*pow_b1b2 + b1*gam2*pow_b1g2)
      else 
        RSDskew_fmu2v1sq = A_fast**2*(b1**2*pow_b1sq ) !eTNS 
      end if
      
    end function

!2.
    function RSDskew_f2mu4v1v2(arg_k) !adding up 3 biases term that go as f^2 mu^4 v1v2
      implicit none
      real(kind=dp)   RSDskew_f2mu4v1v2, arg_k
      real(kind=dp)   pow_b1, pow_b2, pow_g2

      call splint(k_table, f2mu4v1v2b1_table, f2mu4v1v2b1_d_table, nk_table, arg_k, pow_b1)
      call splint(k_table, f2mu4v1v2b2_table, f2mu4v1v2b2_d_table, nk_table, arg_k, pow_b2)
      call splint(k_table, f2mu4v1v2g2_table, f2mu4v1v2g2_d_table, nk_table, arg_k, pow_g2)

      if(tns_full==1)then
        RSDskew_f2mu4v1v2 = A_fast**2*(b1*pow_b1 + b2*pow_b2 + gam2*pow_g2)
      else
        RSDskew_f2mu4v1v2 = A_fast**2*(b1*pow_b1) !eTNS
      end if
      
    end function

!3.
    function RSDskew_f2mu2v1v2(arg_k) !adding up 3 biases term that go as f^2 mu^2 v1v2
      implicit none
      real(kind=dp)   RSDskew_f2mu2v1v2, arg_k
      real(kind=dp)   pow_b1, pow_b2, pow_g2

      call splint(k_table, f2mu2v1v2b1_table, f2mu2v1v2b1_d_table, nk_table, arg_k, pow_b1)
      call splint(k_table, f2mu2v1v2b2_table, f2mu2v1v2b2_d_table, nk_table, arg_k, pow_b2)
      call splint(k_table, f2mu2v1v2g2_table, f2mu2v1v2g2_d_table, nk_table, arg_k, pow_g2)

      if(tns_full==1)then
        RSDskew_f2mu2v1v2 = A_fast**2*(b1*pow_b1 + b2*pow_b2 + gam2*pow_g2) 
      else
        RSDskew_f2mu2v1v2 = A_fast**2*(b1*pow_b1) !eTNS 
      end if
      
    end function

!4.
    function RSDskew_f3mu4v2sq(arg_k) ! f^3 mu^4 v2^2 vel-vel term
      implicit none
      real(kind=dp)   RSDskew_f3mu4v2sq, arg_k
      real(kind=dp)   pow_vv

      call splint(k_table, f3mu4v2sq_table, f3mu4v2sq_d_table, nk_table, arg_k, pow_vv)

      RSDskew_f3mu4v2sq = A_fast**2*pow_vv
      
    end function

!5.
    function RSDskew_f3mu6v2sq(arg_k) ! f^3 mu^6 v2^2 vel-vel term
      implicit none
      real(kind=dp)   RSDskew_f3mu6v2sq, arg_k
      real(kind=dp)   pow_vv

      call splint(k_table, f3mu6v2sq_table, f3mu6v2sq_d_table, nk_table, arg_k, pow_vv)

      RSDskew_f3mu6v2sq = A_fast**2*pow_vv
      
    end function


! B_TNS or P22RS terms:
!************************************************************
!1.
    function RSDkurt_f2mu2v1sq(arg_k) ! f^2 mu^2 v1^2 b1^2 term
      implicit none
      real(kind=dp)   RSDkurt_f2mu2v1sq, arg_k
      real(kind=dp)   pow_vv

      call splint(k_table, f2mu2v1sqb1sq_table, f2mu2v1sqb1sq_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f2mu2v1sq = A_fast**2*b1**2*pow_vv
      
    end function

!2.
    function RSDkurt_f2mu4v1sq(arg_k) ! f^2 mu^4 v1^2 b1^2 term
      implicit none
      real(kind=dp)   RSDkurt_f2mu4v1sq, arg_k
      real(kind=dp)   pow_vv

      call splint(k_table, f2mu4v1sqb1sq_table, f2mu4v1sqb1sq_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f2mu4v1sq = A_fast**2*(b1**2 *pow_vv)
      
    end function

!3.
    function RSDkurt_f3mu4v1v2(arg_k) ! f^3 mu^4 v1v2 b1 term
      implicit none
      real(kind=dp)   RSDkurt_f3mu4v1v2, arg_k
      real(kind=dp)   pow_vv

      call splint(k_table, f3mu4v1v2b1_table, f3mu4v1v2b1_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f3mu4v1v2 = A_fast**2*b1*pow_vv
      
    end function

!4.
    function RSDkurt_f3mu6v1v2(arg_k) ! f^3 mu^6 v1v2 b1 term
      implicit none
      real(kind=dp)   RSDkurt_f3mu6v1v2, arg_k
      real(kind=dp)   pow_vv

      call splint(k_table, f3mu6v1v2b1_table, f3mu6v1v2b1_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f3mu6v1v2 = A_fast**2*b1*pow_vv
      
    end function

!5.
    function RSDkurt_f4mu4v2sq(arg_k) ! f^4 mu^4 v2^2 term
      implicit none
      real(kind=dp)   RSDkurt_f4mu4v2sq, arg_k
      real(kind=dp)   pow_vv

      call splint(k_table, f4mu4v2sq_table, f4mu4v2sq_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f4mu4v2sq = A_fast**2*pow_vv
      
    end function

!6.
    function RSDkurt_f4mu6v2sq(arg_k) ! f^4 mu^6 v2^2 term
      implicit none
      real(kind=dp)   RSDkurt_f4mu6v2sq, arg_k
      real(kind=dp)   pow_vv

      call splint(k_table, f4mu6v2sq_table, f4mu6v2sq_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f4mu6v2sq = A_fast**2*pow_vv
      
    end function

!7.
    function RSDkurt_f4mu8v2sq(arg_k) ! f^4 mu^8 v2^2 term
      implicit none
      real(kind=dp)   RSDkurt_f4mu8v2sq, arg_k
      real(kind=dp)   pow_vv

      call splint(k_table, f4mu8v2sq_table, f4mu8v2sq_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f4mu8v2sq = A_fast**2*pow_vv
      
    end function

!8.
    function RSDkurt_f4mu2v2sq(arg_k) ! f^4 mu^2 v2^2 term
      implicit none
      real(kind=dp)   RSDkurt_f4mu2v2sq, arg_k
      real(kind=dp)   pow_vv

      call splint(k_table, f4mu2v2sq_table, f4mu2v2sq_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f4mu2v2sq = A_fast**2*pow_vv
      
    end function

!9.
    function RSDkurt_f3mu2v1v2(arg_k) ! f^3 mu^2 v1v2 b1 term
      implicit none
      real(kind=dp)   RSDkurt_f3mu2v1v2, arg_k
      real(kind=dp)   pow_vv

      call splint(k_table, f3mu2v1v2b1_table, f3mu2v1v2b1_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f3mu2v1v2 = A_fast**2*b1*pow_vv
      
    end function


!************************************************************
! old version of multipoles with lorentzian or gaussian FOG
!************************************************************

    subroutine RSD_monopole(arg_k, w_dd, w_dv, w_vv)
      implicit   none
      real(kind=dp)    :: arg_k, w_dd, w_dv, w_vv, x, y, s, z, a !a is a new parameter !

      x = sigv_P*arg_k
!      y = atan(x/sqrt(2._dp))
!      y = asinh(x) ! for sqrt lorentz
!      s = sqrt(1._dp+x**2._dp)
      if(fog_model==1)then
        y = erf(x)*sqrt(pi) !for gaussian
        z = exp(-x**2) !for gaussian
      elseif(fog_model==2)then
        y = sqrt(2._dp)*atan(x/sqrt(2._dp))
      end if
      a = a_roman

!!  if(FoG.eq.'Lorentz') !!
  
      if(sigv_P /= 0._dp) then
       ! w_dd = sqrt(2._dp)*y/x
       ! w_dv = (4._dp*f*(x - sqrt(2._dp)*y))/x**3._dp
       ! w_vv = 2._dp*f**2._dp*(-6._dp*x + x**3 + 6._dp*sqrt(2._dp)*y)/(3._dp*x**5._dp)  
       ! w_dd = y/x
       ! w_dv = f*(s*x-y)/x**3
       ! w_vv = f**2*(-3._dp*s*x+2._dp*s*x**3+3._dp*y)/(8._dp*x**5)
       ! write(*,*)w_dd, w_dv, w_vv
       if(fog_model==1)then
         ! a-model with gaussian: (uses x,y,z,a
          w_dd = y/(2._dp*x)
          w_dv = -f*(3._dp*a*f**2*y-2._dp*sigv_P**2*y-6._dp*a*f**2*x*z+4._dp*sigv_P**2*x*z-4._dp*a*f**2*x**3*z)&
                 &/(4._dp*sigv_P**2*x**3)
          w_vv = -f**2*(15._dp*a*f**2*y-6._dp*sigv_P**2*y-30._dp*a*f**2*x*z+12._dp*sigv_P**2*x*z-&
                 &20._dp*a*f**2*x**3*z+8._dp*sigv_P**2*x**3*z-8._dp*a*f**2*x**5*z)/(16._dp*sigv_P**2*x**5)
       elseif(fog_model==2)then
         !a-model with lorentzian:
         w_dd = y/x
         w_dv = 4._dp*f*(6._dp*a**2*f**2*x + 3._dp*sigv_P**2*x - a**2*f**2*x**3 - 6._dp*a**2*f**2*y - &
                & 3._dp*sigv_P**2*y)/(3._dp*sigv_P**2*x**3)
         w_vv = 2._dp*f**2*(-60._dp*a**2*f**2*x - 30._dp*sigv_P**2*x + 10._dp*a**2*f**2*x**3 + &
                &5._dp*sigv_P**2*x**3 - 3._dp*a**2*f**2*x**5 + 60._dp*a**2*f**2*y + 30._dp*sigv_P**2*y)/(15._dp*sigv_P**2*x**5)
       end if

        !write(*,*)w_dd, w_dv, w_vv
        !read(*,*)
      else 
        w_dd = 1._dp
        w_dv = f*2._dp/3._dp
        w_vv = f**2/5._dp
      endif

    end subroutine RSD_monopole
 
    subroutine RSD_quadrupole(arg_k,w_dd,w_dv,w_vv)
      implicit   none
      real(kind=dp)     arg_k,w_dd,w_dv,w_vv,x,y, s, z, a !a is a new parameter ! 
     
      x = sigv_P*arg_k
!      y = atan(x/sqrt(2._dp))
!      y = asinh(x) ! for sqrt lorentz
!      s = sqrt(1._dp+x**2._dp)
      if(fog_model==1)then
        y = erf(x)*sqrt(pi)
        z = exp(-x**2)
      elseif(fog_model==2)then
        y = sqrt(2._dp)*atan(x/sqrt(2._dp))
      end if
      a = a_roman
     
     !  if(FoG.eq.'Lorentz') !!
      
      if(sigv_P /= 0._dp .and. x > 0.014_dp ) then
!        w_dd = (-5._dp*(-6._dp*x + sqrt(2._dp)*(6._dp + x**2)*y))/(2._dp*x**3)
!        w_dv = (10._dp*f*(-6._dp*x + sqrt(2._dp)*(6._dp + x**2)*y))/x**5
!        w_vv = (2._dp*f**2*(2._dp*x*(45._dp + x**4) - 15._dp*sqrt(2._dp)*(6._dp + x**2)*y))/(3._dp*x**7)
!        w_dd = -5._dp*(-3._dp*s*x+3._dp*y+2._dp*x**2*y)/(4._dp*x**3)
!        w_dv = 5._dp*f*(-9._dp*s*x+2._dp*s*x**3+9._dp*y+4._dp*x**2*y)/(8._dp*x**5)
!        w_vv = 5._dp*f**2*(15._dp*s*x-4._dp*s*x**3+4._dp*s*x**5-15._dp*y-6._dp*x**2*y)/(32._dp*x**7)
        if(fog_model==1)then
          w_dd = -5._dp*(-3._dp*y+2._dp*x**2*y+6._dp*x*z)/(8._dp*x**3)
          w_dv = -5._dp*f*(45._dp*a*f**2*y-18._dp*sigv_P**2*y-6._dp*a*f**2*x**2*y+4._dp*sigv_P**2*x**2*y- &  
                 90._dp*a*f**2*x*z+36._dp*sigv_P**2*x*z-48._dp*a*f**2*x**3*z+16._dp*sigv_P**2*x**3*z - &
                 16._dp*a*f**2*x**5*z)/(16._dp*sigv_P**2*x**5)
          w_vv = -5._dp*f**2*(315._dp*a*f**2*y-90._dp*sigv_P**2*y-30._dp*a*f**2*x**2*y+12._dp*sigv_P**2*x**2*y - &
                 630._dp*a*f**2*x*z+180._dp*sigv_P**2*x*z-360._dp*a*f**2*x**3*z+96._dp*sigv_P**2*x**3*z - &  
                 128._dp*a*f**2*x**5*z+32._dp*sigv_P**2*x**5*z-32._dp*a*f**2*x**7*z)/(64._dp*sigv_P**2*x**7)
        elseif(fog_model==2)then
          w_dd = -5._dp*(-6._dp*x + 6._dp*y + x**2*y)/(2._dp*x**3)
          w_dv = 2._dp*f*(-180._dp*a**2*f**2*x - 90._dp*sigv_P**2*x - 4._dp*a**2*f**2*x**5 + 180._dp*a**2*f**2*y + & 
                 90._dp*sigv_P**2*y + 30._dp*a**2*f**2*x**2*y + 15._dp*sigv_P**2*x**2*y)/(3._dp*sigv_P**2*x**5)
          w_vv = 2._dp*f**2*(1260._dp*a**2*f**2*x + 630._dp*sigv_P**2*x + 28._dp*a**2*f**2*x**5 + 14._dp*sigv_P**2*x**5 - & 
                 12._dp*a**2*f**2*x**7 - 1260._dp*a**2*f**2*y - 630._dp*sigv_P**2*y - 210._dp*a**2*f**2*x**2*y - &
                 105._dp*sigv_P**2*x**2*y)/(21._dp*sigv_P**2*x**7)
        endif
      else 
        w_dd = 0._dp
        w_dv = (4._dp*f)/3._dp
        w_vv = (4._dp*f**2)/7._dp
      endif

    end subroutine RSD_quadrupole

    subroutine RSD_hexadecapole(arg_k,w_dd,w_dv,w_vv)
      implicit   none
      real(kind=dp)     arg_k,w_dd,w_dv,w_vv,x,y,s, z, a !a is a new parameter ! 
 
      x = sigv_P*arg_k
!      y = atan(x/sqrt(2._dp))
!      y = asinh(x) ! for sqrt lorentz
!      s = sqrt(1._dp+x**2._dp)
      if(fog_model==1)then
        y = erf(x)*sqrt(pi)
        z = exp(-x**2)
      elseif(fog_model==2)then
        y = sqrt(2._dp)*atan(x/sqrt(2._dp))
      end if
      a = a_roman
 
      !!  if(FoG.eq.'Lorentz') !!
   
      if(sigv_P /= 0._dp .and. x > 0.055) then
!        w_dd = (3._dp*(-10._dp*x*(42._dp + 11._dp*x**2) + &
!               3._dp*sqrt(2._dp)*(140._dp + 60._dp*x**2 + 3._dp*x**4)*y))/(8._dp*x**5)
!        w_dv = (-3._dp*f*(-10._dp*x*(42._dp + 11._dp*x**2) + &
!               3._dp*sqrt(2._dp)*(140._dp + 60._dp*x**2 + 3._dp*x**4)*y))/(2._dp*x**7)
!        w_vv = (3._dp*f**2*(-10._dp*x*(42._dp + 11._dp*x**2) + &
!               3._dp*sqrt(2._dp)*(140._dp + 60._dp*x**2 + 3._dp*x**4)*y))/(2._dp*x**9)
!        w_dd = 9._dp*(-105._dp*s*x-50._dp*s*x**3+105._dp*y+120._dp*x**2*y+24._dp*x**4*y)/(64._dp*x**5)
!        w_dv = -3._dp*f*(-525._dp*s*x-190._dp*s*x**3+8._dp*s*x**5+525._dp*y+540._dp*x**2*y + &
!               72._dp*x**4*y)/(64._dp*x**7)
!        w_vv = 3._dp*f**2*(-3675._dp*s*x-1150._dp*s*x**3+8._dp*s*x**5+48._dp*s*x**7+3675._dp*y + & 
!               3600._dp*x**2*y+432._dp*x**4*y)/(1024._dp*x**9)
        if(fog_model==1)then
          w_dd = 9._dp*(105._dp*y-60._dp*x**2*y+12._dp*x**4*y-210._dp*x*z-20._dp*x**3*z)/(64._dp*x**5)
          w_dv = -9._dp*f*(3675._dp*a*f**2*y-1050._dp*sigv_P**2*y-900._dp*a*f**2*x**2*y+360._dp*sigv_P**2*x**2*y+ &  
                36._dp*a*f**2*x**4*y-24._dp*sigv_P**2*x**4*y-7350._dp*a*f**2*x*z+2100._dp*sigv_P**2*x*z - &  
                3100._dp*a*f**2*x**3*z+680._dp*sigv_P**2*x**3*z-832._dp*a*f**2*x**5*z + & 
                128._dp*sigv_P**2*x**5*z-128._dp*a*f**2*x**7*z)/(128._dp*sigv_P**2*x**7)
          w_vv = -9._dp*f**2*(33075._dp*a*f**2*y-7350._dp*sigv_P**2*y-6300._dp*a*f**2*x**2*y + &  
                1800._dp*sigv_P**2*x**2*y+180._dp*a*f**2*x**4*y-72._dp*sigv_P**2*x**4*y-66150._dp*a*f**2*x*z + &
                14700._dp*sigv_P**2*x*z-31500._dp*a*f**2*x**3*z+6200._dp*sigv_P**2*x**3*z - & 
                9600._dp*a*f**2*x**5*z+1664._dp*sigv_P**2*x**5*z-1920._dp*a*f**2*x**7*z + &  
                256._dp*sigv_P**2*x**7*z-256._dp*a*f**2*x**9*z)/(512._dp*sigv_P**2*x**9)
        elseif(fog_model==2)then
          w_dd = 3._dp*(-420._dp*x - 110._dp*x**3 + 420._dp*y + 180._dp*x**2*y + 9._dp*x**4*y)/(8._dp*x**5)
          w_dv = -3._dp*f*(2._dp*a**2*f**2 + sigv_P**2)*(-420._dp*x - 110._dp*x**3 + 420._dp*y + 180._dp*x**2*y + &
               9._dp*x**4*y)/(2._dp*sigv_P**2*x**7)
          w_vv = f**2*(-88200._dp*a**2*f**2*x - 44100._dp*sigv_P**2*x - 23100._dp*a**2*f**2*x**3 - &
               11550._dp*sigv_P**2*x**3 - 32._dp*a**2*f**2*x**9 + 88200._dp*a**2*f**2*y + 44100._dp*sigv_P**2*y + & 
               37800._dp*a**2*f**2*x**2*y + 18900._dp*sigv_P**2*x**2*y + 1890._dp*a**2*f**2*x**4*y + &
               945._dp*sigv_P**2*x**4*y)/(70._dp*sigv_P**2*x**9)
        endif
      else 
        w_dd = 0._dp
        w_dv = 0._dp
        w_vv = (8._dp*f**2)/35._dp
      endif
 
    end subroutine RSD_hexadecapole
    
!************************************************************
! multipoles with flexible FOG factor that gets integrated separately
!************************************************************
    subroutine RSD_Multipole_novir(arg_k, w_dd, w_dv, w_vv, wskew_fmu2v1sq, wskew_f2mu2v1v2, &
               wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq, wkurt_f2mu2v1sq, &
               wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f4mu4v2sq, &
               wkurt_f4mu6v2sq, wkurt_f4mu8v2sq, wkurt_f4mu2v2sq, wkurt_f3mu2v1v2, l)
      implicit   none
      real(kind=dp)    :: arg_k, w_dd, w_dv, w_vv
      real(kind=dp)    :: wskew_fmu2v1sq, wskew_f2mu2v1v2, wskew_f2mu4v1v2
      real(kind=dp)    :: wskew_f3mu4v2sq, wskew_f3mu6v2sq, wkurt_f2mu2v1sq
      real(kind=dp)    :: wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2
      real(kind=dp)    :: wkurt_f4mu4v2sq, wkurt_f4mu6v2sq, wkurt_f4mu8v2sq
      real(kind=dp)    :: wkurt_f4mu2v2sq, wkurt_f3mu2v1v2

      integer :: l

        if (l.eq.0) then
           w_dd = 1._dp
           w_dv = f*2._dp/3._dp
           w_vv = f**2/5._dp
           wskew_fmu2v1sq =  f/3._dp
           wskew_f2mu2v1v2 = f**2/3._dp
           wskew_f2mu4v1v2 = f**2/5._dp
           wskew_f3mu4v2sq = f**3/5._dp
           wskew_f3mu6v2sq = f**3/7._dp
           wkurt_f2mu2v1sq = f**2/3._dp
           wkurt_f2mu4v1sq = f**2/5._dp
           wkurt_f3mu4v1v2 = f**3/5._dp
           wkurt_f3mu6v1v2 = f**3/7._dp
           wkurt_f4mu4v2sq = f**4/5._dp
           wkurt_f4mu6v2sq = f**4/7._dp
           wkurt_f4mu8v2sq = f**4/9._dp
           wkurt_f4mu2v2sq = f**4/3._dp
           wkurt_f3mu2v1v2 = f**3/3._dp
        elseif (l.eq.2) then
           w_dd = 0._dp
           w_dv = (4._dp*f)/3._dp
           w_vv = (4._dp*f**2)/7._dp
           wskew_fmu2v1sq =  f    *(2._dp/3._dp)
           wskew_f2mu2v1v2 = f**2 *(2._dp/3._dp)
           wskew_f2mu4v1v2 = f**2 *(4._dp/7._dp)
           wskew_f3mu4v2sq = f**3 *(4._dp/7._dp)
           wskew_f3mu6v2sq = f**3 *(10._dp/21._dp)
           wkurt_f2mu2v1sq = f**2 *(2._dp/3._dp)
           wkurt_f2mu4v1sq = f**2 *(4._dp/7._dp)
           wkurt_f3mu4v1v2 = f**3 *(4._dp/7._dp)
           wkurt_f3mu6v1v2 = f**3 *(10._dp/21._dp)
           wkurt_f4mu4v2sq = f**4 *(4._dp/7._dp)
           wkurt_f4mu6v2sq = f**4 *(10._dp/21._dp)
           wkurt_f4mu8v2sq = f**4 *(40._dp/99._dp)
           wkurt_f4mu2v2sq = f**4 *(2._dp/3._dp)
           wkurt_f3mu2v1v2 = f**3 *(2._dp/3._dp)
        elseif (l.eq.4) then
           w_dd = 0._dp
           w_dv = 0._dp
           w_vv = (8._dp*f**2)/35._dp
           wskew_fmu2v1sq =  0._dp
           wskew_f2mu2v1v2 = 0._dp
           wskew_f2mu4v1v2 = f**2 *(8._dp/35._dp)
           wskew_f3mu4v2sq = f**3 *(8._dp/35._dp)
           wskew_f3mu6v2sq = f**3 *(24._dp/77._dp)
           wkurt_f2mu2v1sq = 0._dp
           wkurt_f2mu4v1sq = f**2 *(8._dp/35._dp)
           wkurt_f3mu4v1v2 = f**3 *(8._dp/35._dp)
           wkurt_f3mu6v1v2 = f**3 *(24._dp/77._dp)
           wkurt_f4mu4v2sq = f**4 *(8._dp/35._dp)
           wkurt_f4mu6v2sq = f**4 *(24._dp/77._dp)
           wkurt_f4mu8v2sq = f**4 *(48._dp/143._dp)
           wkurt_f4mu2v2sq = 0._dp
           wkurt_f3mu2v1v2 = 0._dp
        elseif (l.eq.6) then
           w_dd = 0._dp
           w_dv = 0._dp
           w_vv = 0._dp
           wskew_fmu2v1sq =  0._dp
           wskew_f2mu2v1v2 = 0._dp
           wskew_f2mu4v1v2 = 0._dp
           wskew_f3mu4v2sq = 0._dp
           wskew_f3mu6v2sq = f**3 *(16._dp/231._dp)
           wkurt_f2mu2v1sq = 0._dp
           wkurt_f2mu4v1sq = 0._dp
           wkurt_f3mu4v1v2 = 0._dp
           wkurt_f3mu6v1v2 = f**3 *(16._dp/231._dp)
           wkurt_f4mu4v2sq = 0._dp
           wkurt_f4mu6v2sq = f**4 *(16._dp/231._dp)
           wkurt_f4mu8v2sq = f**4 *(64._dp/495._dp)
           wkurt_f4mu2v2sq = 0._dp
           wkurt_f3mu2v1v2 = 0._dp
        elseif (l.eq.8) then
           w_dd = 0._dp
           w_dv = 0._dp
           w_vv = 0._dp
           wskew_fmu2v1sq =  0._dp
           wskew_f2mu2v1v2 = 0._dp
           wskew_f2mu4v1v2 = 0._dp
           wskew_f3mu4v2sq = 0._dp
           wskew_f3mu6v2sq = 0._dp
           wkurt_f2mu2v1sq = 0._dp
           wkurt_f2mu4v1sq = 0._dp
           wkurt_f3mu4v1v2 = 0._dp
           wkurt_f3mu6v1v2 = 0._dp
           wkurt_f4mu4v2sq = 0._dp
           wkurt_f4mu6v2sq = 0._dp
           wkurt_f4mu8v2sq = f**4 *(128._dp/6435._dp)
           wkurt_f4mu2v2sq = 0._dp
           wkurt_f3mu2v1v2 = 0._dp
        endif 


!           w_dd = 0._dp   !PROBLEMA
!           w_dv = 0._dp
!           w_vv = 0._dp
!           wskew_fmu2v1sq =  0._dp
!           wskew_f2mu2v1v2 = 0._dp
!           wskew_f2mu4v1v2 = 0._dp
!           wskew_f3mu4v2sq = 0._dp
!          ! wskew_f3mu6v2sq = 0._dp   
!           wkurt_f2mu2v1sq = 0._dp
!           wkurt_f2mu4v1sq = 0._dp
!           wkurt_f3mu4v1v2 = 0._dp
!           wkurt_f3mu6v1v2 = 0._dp
!           wkurt_f4mu4v2sq = 0._dp
!           wkurt_f4mu6v2sq = 0._dp
!           wkurt_f4mu8v2sq = 0._dp
!           wkurt_f4mu2v2sq = 0._dp
!           wkurt_f3mu2v1v2 = 0._dp


    end subroutine RSD_Multipole_novir

    subroutine RSD_Multipole(arg_k, w_dd, w_dv, w_vv, wskew_fmu2v1sq, wskew_f2mu2v1v2, &
               wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq, wkurt_f2mu2v1sq, &
               wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f4mu4v2sq, &
               wkurt_f4mu6v2sq, wkurt_f4mu8v2sq, wkurt_f4mu2v2sq, wkurt_f3mu2v1v2, l)
      implicit   none
      real(kind=dp)    :: arg_k, w_dd, w_dv, w_vv
      real(kind=dp)    :: wskew_fmu2v1sq, wskew_f2mu2v1v2, wskew_f2mu4v1v2
      real(kind=dp)    :: wskew_f3mu4v2sq, wskew_f3mu6v2sq, wkurt_f2mu2v1sq
      real(kind=dp)    :: wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2
      real(kind=dp)    :: wkurt_f4mu4v2sq, wkurt_f4mu6v2sq, wkurt_f4mu8v2sq
      real(kind=dp)    :: wkurt_f4mu2v2sq, wkurt_f3mu2v1v2
      real(kind=dp)    :: rsdv3
      integer :: l
 
         w_dd = RSD_virial3(arg_k,0,l) !*0._dp

         rsdv3 = RSD_virial3(arg_k,1,l) 
         w_dv            = rsdv3 *2._dp*f  !*0._dp 
         wskew_fmu2v1sq  = rsdv3 *f        !*0._dp
         wskew_f2mu2v1v2 = rsdv3 *f**2     !*0._dp
         wkurt_f2mu2v1sq = rsdv3 *f**2     !*0._dp
         wkurt_f3mu2v1v2 = rsdv3 *f**3     !*0._dp
         wkurt_f4mu2v2sq = rsdv3 *f**4     !*0._dp

         rsdv3 = RSD_virial3(arg_k,2,l) 
         w_vv            = rsdv3 *f**2  !*0._dp
         wskew_f2mu4v1v2 = rsdv3 *f**2  !*0._dp
         wkurt_f2mu4v1sq = rsdv3 *f**2  !*0._dp
         wskew_f3mu4v2sq = rsdv3 *f**3  !*0._dp
         wkurt_f3mu4v1v2 = rsdv3 *f**3  !*0._dp
         wkurt_f4mu4v2sq = rsdv3 *f**4  !*0._dp

         rsdv3 = RSD_virial3(arg_k,3,l) 
         wskew_f3mu6v2sq = rsdv3 *f**3  !*0._dp         
         wkurt_f3mu6v1v2 = rsdv3 *f**3  !*0._dp
         wkurt_f4mu6v2sq = rsdv3 *f**4  !*0._dp

         wkurt_f4mu8v2sq = RSD_virial3(arg_k,4,l) *f**4  !*0._dp

    end subroutine RSD_Multipole

    function RSD_virial3(arg_k,m,l) ! m for mu^2 power, ell for multipole

      implicit none
      real(kind=dp)   :: RSD_virial3,arg_k,epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit), l, m
      real(kind=dp)   :: ans
 
      ell = l ! multipole desired
      nmu = m ! power of mu^2
      k = arg_k
      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (virial3_dmu,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      RSD_virial3 = ans  *(2._dp*real(ell,kind=dp) + 1._dp) 
      if(ier /= 0) call warning('RSD_virial3',ier)
      if(ier /= 0) write(*,*)k        
    end function RSD_virial3


    function virial3_dmu(mu)        
      implicit none
      real(kind=dp)   :: virial3_dmu, mu, a, sigv_eff, b
      a = a_roman
      !sigv_eff = dsqrt(sigv_P**2 + f**2 *sigv**2)  
      !new definition
      sigv_eff = dsqrt(A_fast*f**2*sigv**2 + sigv_P)  
!      b = (k*mu*sigv_B)**3/(1._dp+a**2*k**2*mu**2)
      
      virial3_dmu = (mu**2)**nmu *LegendreP(mu) /dsqrt(1._dp+a**2*k**2*mu**2)& 
                *exp(-mu**2 *k**2 *sigv_eff**2/(1._dp+a**2*k**2*mu**2) ) !* sphbe0(b)

    end function virial3_dmu


    function RSD_virial4(arg_k,arg_z,l,lp)  

      implicit none
      real(kind=dp)   :: RSD_virial4,arg_k,epsabs,epsrel,a,b,arg_z
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit), n, l, m, lp
      real(kind=dp)   :: ans
 
      ell = l ! multipole desired (external)
      ellp = lp ! multipole desired (internal)
      vzdisp = arg_z ! scale-dependent velocity dispersion
      k = arg_k
      a = 0._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (virial4_dmu,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      RSD_virial4 = ans  !*(2._dp*real(ell,kind=dp) + 1._dp) we moved this factor to K(ell,ellp) !!!
      if(ier /= 0) call warning('RSD_virial4',ier)
      if(ier /= 0) write(*,*)k        
    end function RSD_virial4


    function virial4_dmu(mu)        
      implicit none
      real(kind=dp)   :: virial4_dmu, mu, a
      a = a_roman

      virial4_dmu = LegendrePp(mu) *LegendreP(mu) /dsqrt(1._dp+a**2*k**2*mu**2)& 
                *exp(-mu**2 *k**2 *vzdisp/(1._dp+a**2*k**2*mu**2) )

    end function virial4_dmu


    function RSD_virial2(arg_k,n,m,l) ! n=0,1,2 for v2 power, m for mu^2 power, ell for multipole

      implicit none
      real(kind=dp)   :: RSD_virial2,arg_k,epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit), n, l, m
      real(kind=dp)   :: ans
 
      nv2 = n ! power of v2
      ell = l ! multipole desired
      nmu = m ! power of mu^2
      k = arg_k
      a = -1._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (virial2_dmu,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      RSD_virial2 = ans *0.5_dp *(2._dp*real(ell,kind=dp) + 1._dp) 
      if(ier /= 0) call warning('RSD_virial2',ier)
      if(ier /= 0) write(*,*)k        
    end function RSD_virial2

    function virial2_dmu(mu)        
      implicit none
      real(kind=dp)   :: virial2_dmu, mu, a
      a = a_roman

!      virial2_dmu = (mu**2)**nmu *(1._dp/(1._dp-0.5_dp*a**2 *mu**2 *gRPT_roots(k)/20._dp))**nv2 & 
!                *exp(mu**2 *gRPT_roots(k) *sigv_P**2/20._dp) *LegendreP(mu)

      virial2_dmu = (mu**2)**nmu *(1._dp/(1._dp + 0.5_dp*a**2 *mu**2 *k**2))**nv2 & 
                *exp(-mu**2 *k**2 *sigv_P**2) *LegendreP(mu)

!      virial2_dmu = (mu**2)**nmu *(1._dp/(1._dp + 0.5_dp*a**2 *mu**2 *k**2)) & 
!                *exp(-mu**2 *k**2 *sigv_P**2)**nv2 *LegendreP(mu)

! the factor of 20 above normalizes roots so needed values of sigv_P are reasonable
! -k^2 -> + roots/20

    end function virial2_dmu

    function RSD_virial(arg_k,n,l) !n n=0,1,2 for w_dd w_dv w_vv, ell for multipole

      implicit none
      real(kind=dp)   :: RSD_virial,arg_k,epsabs,epsrel,a,b
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit), n, l
      real(kind=dp)   :: ans
 
      ndenvel = n !defines which virial factors we need to average
      ell = l !multipole desired
      k = arg_k
      a = -1._dp      ! lower limit of integration
      b = 1._dp       ! upper limit of integration
      epsrel = relacc
      epsabs = absacc
      ans = 0._dp
 
      call qage (virial_dmu,a,b,epsabs,epsrel,key,limit,ans, &
      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      RSD_virial = ans *0.5_dp *(2._dp*real(ell,kind=dp) + 1._dp) 
      if(ier /= 0) call warning('RSD_virial',ier)
        
    end function RSD_virial

    function virial_dmu(mu)        
      implicit none
      real(kind=dp)   :: virial_dmu, mu, a
      a = a_roman
!      virial_dmu = (mu**2/(1._dp+0.5_dp*a**2 *mu**2 *k**2))**ndenvel & 
      virial_dmu = (mu**2/(1._dp-0.5_dp*a**2 *mu**2 *gRPT_roots(k)/20._dp))**ndenvel & 
                *exp(mu**2 *gRPT_roots(k) *sigv_P**2/20._dp) *LegendreP(mu)
!                *exp(-mu**2 *k**2 *sigv_P**2) *LegendreP(mu)
! the factor of 20 above normalizes roots so needed values of sigv_P are reasonable
    end function virial_dmu
 
    function LegendreP(mu)        
      implicit none
      real(kind=dp)   :: mu, LegendreP
      
      if (ell == 0) then
         LegendreP = 1._dp
      elseif (ell == 2) then
         LegendreP=(3._dp*mu**2 - 1._dp)/2._dp
      elseif (ell == 4) then
         LegendreP=(35._dp*mu**4-30._dp*mu**2+3._dp)/8._dp
      elseif (ell == 6) then
         LegendreP=(231._dp*mu**6-315._dp*mu**4+105._dp*mu**2-5._dp)/16._dp
      elseif (ell == 8) then
         LegendreP=(6435._dp*mu**8-12012._dp*mu**6+6930._dp*mu**4-1260._dp*mu**2+35._dp)/128._dp
      else
         write(*,*)'this Legendre Polynomial is not coded up',ell
         stop
      endif   
      
    end function LegendreP

    function LegendrePp(mu)        
      implicit none
      real(kind=dp)   :: mu, LegendrePp
      
      if (ellp == 0) then
         LegendrePp = 1._dp
      elseif (ellp == 2) then
         LegendrePp=(3._dp*mu**2 - 1._dp)/2._dp
      elseif (ellp == 4) then
         LegendrePp=(35._dp*mu**4-30._dp*mu**2+3._dp)/8._dp
      elseif (ellp == 6) then
         LegendrePp=(231._dp*mu**6-315._dp*mu**4+105._dp*mu**2-5._dp)/16._dp
      elseif (ellp == 8) then
         LegendrePp=(6435._dp*mu**8-12012._dp*mu**6+6930._dp*mu**4-1260._dp*mu**2+35._dp)/128._dp
      else
         write(*,*)'this Legendre Polynomial is not coded up',ellp
         stop
      endif   
      
    end function LegendrePp

    function LegendreP1(mu)        
      implicit none
      real(kind=dp)   :: mu, LegendreP1
      LegendreP1=mu
    end function LegendreP1

    function LegendreP2(mu)        
      implicit none
      real(kind=dp)   :: mu, LegendreP2
      LegendreP2=(3._dp*mu**2 - 1._dp)/2._dp
    end function LegendreP2

    function LegendreP3(mu)        
      implicit none
      real(kind=dp)   :: mu, LegendreP3
      LegendreP3=(5._dp*mu**3 - 3._dp*mu)/2._dp
    end function LegendreP3

    function LegendreP4(mu)        
      implicit none
      real(kind=dp)   :: mu, LegendreP4
      LegendreP4=(35._dp*mu**4-30._dp*mu**2+3._dp)/8._dp
    end function LegendreP4

    function LegendreP6(mu)        
      implicit none
      real(kind=dp)   :: mu, LegendreP6
      LegendreP6=(231._dp*mu**6-315._dp*mu**4+105._dp*mu**2-5._dp)/16._dp
    end function LegendreP6
      
    function LegendreP8(mu)        
      implicit none
      real(kind=dp)   :: mu, LegendreP8
      LegendreP8=(6435._dp*mu**8-12012._dp*mu**6+6930._dp*mu**4-1260._dp*mu**2+35._dp)/128._dp
    end function LegendreP8
      
    function sphbessel(l,x)
      implicit none
      real(kind=dp)   :: x, sphbessel
      integer :: l  
      
      if (ell == 0) then
         sphbessel = sphbe0(x)
      elseif (ell == 1) then
         sphbessel = sphbe1(x)
      elseif (ell == 2) then
         sphbessel = sphbe2(x)
      elseif (ell == 3) then
         sphbessel = sphbe3(x)
      elseif (ell == 4) then
         sphbessel = sphbe4(x)
      elseif (ell == 5) then
         sphbessel = sphbe5(x)
      elseif (ell == 6) then
         sphbessel = sphbe6(x)
      elseif (ell == 7) then
         sphbessel = sphbe7(x)
      elseif (ell == 8) then
         sphbessel = sphbe8(x)
      else
         write(*,*)'this spherical bessel is not coded up',ell
         stop
      endif   
      
      end function sphbessel


!************************************************************
! power spectrum multipoles
!************************************************************

    function power_gal_ell0_novir_at(arg_k)
      implicit none
      real(kind=dp)   :: power_gal_ell0_novir_at,arg_k
      real(kind=dp)   :: w_dd, w_dv, w_vv, pkdd, pkdv, pkvv
      real(kind=dp)   :: wskew_fmu2v1sq, wskew_f2mu2v1v2, wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq
      real(kind=dp)   :: wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f2mu2v1sq
      real(kind=dp)   :: wkurt_f4mu4v2sq, wkurt_f4mu6v2sq, wkurt_f4mu8v2sq
      real(kind=dp)   :: wkurt_f4mu2v2sq, wkurt_f3mu2v1v2

      if(do_init_model)then
        call init_model
      end if
        
      call RSD_Multipole_novir(arg_k, w_dd, w_dv, w_vv, wskew_fmu2v1sq, wskew_f2mu2v1v2, &
               wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq, wkurt_f2mu2v1sq, &
               wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f4mu4v2sq, &
               wkurt_f4mu6v2sq, wkurt_f4mu8v2sq, wkurt_f4mu2v2sq, wkurt_f3mu2v1v2, 0)


      pkdd = gRPT_Pgg(arg_k)
      pkdv = gRPT_Pgv(arg_k)
      pkvv = gRPT_Pvv(arg_k)
  
      power_gal_ell0_novir_at = w_dd*pkdd + w_dv*pkdv + w_vv*pkvv &
      + wskew_fmu2v1sq* RSDskew_fmu2v1sq(arg_k) &
      + wskew_f2mu2v1v2* RSDskew_f2mu2v1v2(arg_k) &
      + wskew_f2mu4v1v2* RSDskew_f2mu4v1v2(arg_k) &
      + wskew_f3mu4v2sq* RSDskew_f3mu4v2sq(arg_k) &
      + wskew_f3mu6v2sq* RSDskew_f3mu6v2sq(arg_k) &
      + wkurt_f2mu2v1sq* RSDkurt_f2mu2v1sq(arg_k) &
      + wkurt_f2mu4v1sq* RSDkurt_f2mu4v1sq(arg_k) &
      + wkurt_f3mu4v1v2* RSDkurt_f3mu4v1v2(arg_k) &
      + wkurt_f3mu6v1v2* RSDkurt_f3mu6v1v2(arg_k) &
      + wkurt_f4mu4v2sq* RSDkurt_f4mu4v2sq(arg_k) &
      + wkurt_f4mu6v2sq* RSDkurt_f4mu6v2sq(arg_k) &
      + wkurt_f4mu8v2sq* RSDkurt_f4mu8v2sq(arg_k) &
      + wkurt_f4mu2v2sq* RSDkurt_f4mu2v2sq(arg_k) &
      + wkurt_f3mu2v1v2* RSDkurt_f3mu2v1v2(arg_k)   
    
    end function power_gal_ell0_novir_at
  
    function power_gal_ell2_novir_at(arg_k)
      implicit none
      real(kind=dp)   :: power_gal_ell2_novir_at, arg_k
      real(kind=dp)   :: w_dd, w_dv, w_vv, pkdd, pkdv, pkvv
      real(kind=dp)   :: wskew_fmu2v1sq, wskew_f2mu2v1v2, wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq
      real(kind=dp)    :: wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f2mu2v1sq
      real(kind=dp)    :: wkurt_f4mu4v2sq, wkurt_f4mu6v2sq, wkurt_f4mu8v2sq
      real(kind=dp)   :: wkurt_f4mu2v2sq, wkurt_f3mu2v1v2

      if(do_init_model)then
        call init_model
      end if
 
      call RSD_Multipole_novir(arg_k, w_dd, w_dv, w_vv, wskew_fmu2v1sq, wskew_f2mu2v1v2, &
               wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq, wkurt_f2mu2v1sq, &
               wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f4mu4v2sq, &
               wkurt_f4mu6v2sq, wkurt_f4mu8v2sq, wkurt_f4mu2v2sq, wkurt_f3mu2v1v2, 2)
 
      pkdd = gRPT_Pgg(arg_k)
      pkdv = gRPT_Pgv(arg_k)
      pkvv = gRPT_Pvv(arg_k)
 
      power_gal_ell2_novir_at = w_dd*pkdd + w_dv*pkdv + w_vv*pkvv &
      + wskew_fmu2v1sq* RSDskew_fmu2v1sq(arg_k) &
      + wskew_f2mu2v1v2* RSDskew_f2mu2v1v2(arg_k) &
      + wskew_f2mu4v1v2* RSDskew_f2mu4v1v2(arg_k) &
      + wskew_f3mu4v2sq* RSDskew_f3mu4v2sq(arg_k) &
      + wskew_f3mu6v2sq* RSDskew_f3mu6v2sq(arg_k) &
      + wkurt_f2mu2v1sq* RSDkurt_f2mu2v1sq(arg_k) &
      + wkurt_f2mu4v1sq* RSDkurt_f2mu4v1sq(arg_k) &
      + wkurt_f3mu4v1v2* RSDkurt_f3mu4v1v2(arg_k) &
      + wkurt_f3mu6v1v2* RSDkurt_f3mu6v1v2(arg_k) &
      + wkurt_f4mu4v2sq* RSDkurt_f4mu4v2sq(arg_k) &
      + wkurt_f4mu6v2sq* RSDkurt_f4mu6v2sq(arg_k) &
      + wkurt_f4mu8v2sq* RSDkurt_f4mu8v2sq(arg_k) &
      + wkurt_f4mu2v2sq* RSDkurt_f4mu2v2sq(arg_k) &
      + wkurt_f3mu2v1v2* RSDkurt_f3mu2v1v2(arg_k)   

    end function power_gal_ell2_novir_at
 
    function power_gal_ell4_novir_at(arg_k)
      implicit none
      real(kind=dp)  :: power_gal_ell4_novir_at, arg_k
      real(kind=dp)  :: w_dd, w_dv, w_vv, pkdd, pkdv, pkvv
      real(kind=dp)   :: wskew_fmu2v1sq, wskew_f2mu2v1v2, wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq
      real(kind=dp)    :: wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f2mu2v1sq
      real(kind=dp)    :: wkurt_f4mu4v2sq, wkurt_f4mu6v2sq, wkurt_f4mu8v2sq
      real(kind=dp)   :: wkurt_f4mu2v2sq, wkurt_f3mu2v1v2

      if(do_init_model)then
        call init_model
      end if
 
      call RSD_Multipole_novir(arg_k, w_dd, w_dv, w_vv, wskew_fmu2v1sq, wskew_f2mu2v1v2, &
               wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq, wkurt_f2mu2v1sq, &
               wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f4mu4v2sq, &
               wkurt_f4mu6v2sq, wkurt_f4mu8v2sq, wkurt_f4mu2v2sq, wkurt_f3mu2v1v2, 4)
 
      pkdd = gRPT_Pgg(arg_k)
      pkdv = gRPT_Pgv(arg_k)
      pkvv = gRPT_Pvv(arg_k)
 
      power_gal_ell4_novir_at = w_dd*pkdd + w_dv*pkdv + w_vv*pkvv &
      + wskew_fmu2v1sq* RSDskew_fmu2v1sq(arg_k) &
      + wskew_f2mu2v1v2* RSDskew_f2mu2v1v2(arg_k) &
      + wskew_f2mu4v1v2* RSDskew_f2mu4v1v2(arg_k) &
      + wskew_f3mu4v2sq* RSDskew_f3mu4v2sq(arg_k) &
      + wskew_f3mu6v2sq* RSDskew_f3mu6v2sq(arg_k) &
      + wkurt_f2mu2v1sq* RSDkurt_f2mu2v1sq(arg_k) &
      + wkurt_f2mu4v1sq* RSDkurt_f2mu4v1sq(arg_k) &
      + wkurt_f3mu4v1v2* RSDkurt_f3mu4v1v2(arg_k) &
      + wkurt_f3mu6v1v2* RSDkurt_f3mu6v1v2(arg_k) &
      + wkurt_f4mu4v2sq* RSDkurt_f4mu4v2sq(arg_k) &
      + wkurt_f4mu6v2sq* RSDkurt_f4mu6v2sq(arg_k) &
      + wkurt_f4mu8v2sq* RSDkurt_f4mu8v2sq(arg_k) &
      + wkurt_f4mu2v2sq* RSDkurt_f4mu2v2sq(arg_k) &
      + wkurt_f3mu2v1v2* RSDkurt_f3mu2v1v2(arg_k)   

    end function power_gal_ell4_novir_at
 
    function power_gal_ell6_novir_at(arg_k)
      implicit none
      real(kind=dp)  :: power_gal_ell6_novir_at, arg_k
      real(kind=dp)  :: w_dd, w_dv, w_vv, pkdd, pkdv, pkvv
      real(kind=dp)   :: wskew_fmu2v1sq, wskew_f2mu2v1v2, wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq
      real(kind=dp)    :: wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f2mu2v1sq
      real(kind=dp)    :: wkurt_f4mu4v2sq, wkurt_f4mu6v2sq, wkurt_f4mu8v2sq
      real(kind=dp)   :: wkurt_f4mu2v2sq, wkurt_f3mu2v1v2

      if(do_init_model)then
        call init_model
      end if
 
      call RSD_Multipole_novir(arg_k, w_dd, w_dv, w_vv, wskew_fmu2v1sq, wskew_f2mu2v1v2, &
               wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq, wkurt_f2mu2v1sq, &
               wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f4mu4v2sq, &
               wkurt_f4mu6v2sq, wkurt_f4mu8v2sq, wkurt_f4mu2v2sq, wkurt_f3mu2v1v2, 6)
 
      pkdd = gRPT_Pgg(arg_k)
      pkdv = gRPT_Pgv(arg_k)
      pkvv = gRPT_Pvv(arg_k)
 
      power_gal_ell6_novir_at = w_dd*pkdd + w_dv*pkdv + w_vv*pkvv &
      + wskew_fmu2v1sq* RSDskew_fmu2v1sq(arg_k) &
      + wskew_f2mu2v1v2* RSDskew_f2mu2v1v2(arg_k) &
      + wskew_f2mu4v1v2* RSDskew_f2mu4v1v2(arg_k) &
      + wskew_f3mu4v2sq* RSDskew_f3mu4v2sq(arg_k) &
      + wskew_f3mu6v2sq* RSDskew_f3mu6v2sq(arg_k) &
      + wkurt_f2mu2v1sq* RSDkurt_f2mu2v1sq(arg_k) &
      + wkurt_f2mu4v1sq* RSDkurt_f2mu4v1sq(arg_k) &
      + wkurt_f3mu4v1v2* RSDkurt_f3mu4v1v2(arg_k) &
      + wkurt_f3mu6v1v2* RSDkurt_f3mu6v1v2(arg_k) &
      + wkurt_f4mu4v2sq* RSDkurt_f4mu4v2sq(arg_k) &
      + wkurt_f4mu6v2sq* RSDkurt_f4mu6v2sq(arg_k) &
      + wkurt_f4mu8v2sq* RSDkurt_f4mu8v2sq(arg_k) &
      + wkurt_f4mu2v2sq* RSDkurt_f4mu2v2sq(arg_k) &
      + wkurt_f3mu2v1v2* RSDkurt_f3mu2v1v2(arg_k)   
 
    end function power_gal_ell6_novir_at
    
    function power_gal_ell8_novir_at(arg_k)
      implicit none
      real(kind=dp)  :: power_gal_ell8_novir_at, arg_k
      real(kind=dp)  :: w_dd, w_dv, w_vv, pkdd, pkdv, pkvv
      real(kind=dp)   :: wskew_fmu2v1sq, wskew_f2mu2v1v2, wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq
      real(kind=dp)    :: wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f2mu2v1sq
      real(kind=dp)    :: wkurt_f4mu4v2sq, wkurt_f4mu6v2sq, wkurt_f4mu8v2sq
      real(kind=dp)   :: wkurt_f4mu2v2sq, wkurt_f3mu2v1v2

      if(do_init_model)then
        call init_model
      end if
 
      call RSD_Multipole_novir(arg_k, w_dd, w_dv, w_vv, wskew_fmu2v1sq, wskew_f2mu2v1v2, &
               wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq, wkurt_f2mu2v1sq, &
               wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f4mu4v2sq, &
               wkurt_f4mu6v2sq, wkurt_f4mu8v2sq, wkurt_f4mu2v2sq, wkurt_f3mu2v1v2, 8)
 
      pkdd = gRPT_Pgg(arg_k)
      pkdv = gRPT_Pgv(arg_k)
      pkvv = gRPT_Pvv(arg_k)
 
      power_gal_ell8_novir_at = w_dd*pkdd + w_dv*pkdv + w_vv*pkvv &
      + wskew_fmu2v1sq* RSDskew_fmu2v1sq(arg_k) &
      + wskew_f2mu2v1v2* RSDskew_f2mu2v1v2(arg_k) &
      + wskew_f2mu4v1v2* RSDskew_f2mu4v1v2(arg_k) &
      + wskew_f3mu4v2sq* RSDskew_f3mu4v2sq(arg_k) &
      + wskew_f3mu6v2sq* RSDskew_f3mu6v2sq(arg_k) &
      + wkurt_f2mu2v1sq* RSDkurt_f2mu2v1sq(arg_k) &
      + wkurt_f2mu4v1sq* RSDkurt_f2mu4v1sq(arg_k) &
      + wkurt_f3mu4v1v2* RSDkurt_f3mu4v1v2(arg_k) &
      + wkurt_f3mu6v1v2* RSDkurt_f3mu6v1v2(arg_k) &
      + wkurt_f4mu4v2sq* RSDkurt_f4mu4v2sq(arg_k) &
      + wkurt_f4mu6v2sq* RSDkurt_f4mu6v2sq(arg_k) &
      + wkurt_f4mu8v2sq* RSDkurt_f4mu8v2sq(arg_k) &
      + wkurt_f4mu2v2sq* RSDkurt_f4mu2v2sq(arg_k) &
      + wkurt_f3mu2v1v2* RSDkurt_f3mu2v1v2(arg_k)   
 
    end function power_gal_ell8_novir_at
    

    function SIDpower_gal_monopole_at(arg_k)
      implicit none
      real(kind=dp)   :: SIDpower_gal_monopole_at,arg_k
      real(kind=dp)   :: w_dd, w_dv, w_vv, pkdd, pkdv, pkvv
      real(kind=dp)   :: wskew_fmu2v1sq, wskew_f2mu2v1v2, wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq
      real(kind=dp)   :: wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f2mu2v1sq
      real(kind=dp)   :: wkurt_f4mu4v2sq, wkurt_f4mu6v2sq, wkurt_f4mu8v2sq
      real(kind=dp)   :: wkurt_f4mu2v2sq, wkurt_f3mu2v1v2

      if(do_init_model)then
        call init_model
      end if
  
      call RSD_Multipole(arg_k, w_dd, w_dv, w_vv, wskew_fmu2v1sq, wskew_f2mu2v1v2, &
               wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq, wkurt_f2mu2v1sq, &
               wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f4mu4v2sq, &
               wkurt_f4mu6v2sq, wkurt_f4mu8v2sq, wkurt_f4mu2v2sq, wkurt_f3mu2v1v2, 0)


      pkdd = gRPT_Pgg(arg_k)
      pkdv = gRPT_Pgv(arg_k)
      pkvv = gRPT_Pvv(arg_k)
  
      SIDpower_gal_monopole_at = w_dd*pkdd + w_dv*pkdv + w_vv*pkvv &
      + wskew_fmu2v1sq* RSDskew_fmu2v1sq(arg_k) &
      + wskew_f2mu2v1v2* RSDskew_f2mu2v1v2(arg_k) &
      + wskew_f2mu4v1v2* RSDskew_f2mu4v1v2(arg_k) &
      + wskew_f3mu4v2sq* RSDskew_f3mu4v2sq(arg_k) &
      + wskew_f3mu6v2sq* RSDskew_f3mu6v2sq(arg_k) &
      + wkurt_f2mu2v1sq* RSDkurt_f2mu2v1sq(arg_k) &
      + wkurt_f2mu4v1sq* RSDkurt_f2mu4v1sq(arg_k) &
      + wkurt_f3mu4v1v2* RSDkurt_f3mu4v1v2(arg_k) &
      + wkurt_f3mu6v1v2* RSDkurt_f3mu6v1v2(arg_k) &
      + wkurt_f4mu4v2sq* RSDkurt_f4mu4v2sq(arg_k) &
      + wkurt_f4mu6v2sq* RSDkurt_f4mu6v2sq(arg_k) &
      + wkurt_f4mu8v2sq* RSDkurt_f4mu8v2sq(arg_k) &
      + wkurt_f4mu2v2sq* RSDkurt_f4mu2v2sq(arg_k) &
      + wkurt_f3mu2v1v2* RSDkurt_f3mu2v1v2(arg_k)   
    
    end function SIDpower_gal_monopole_at
  
    function SIDpower_gal_quadrupole_at(arg_k)
      implicit none
      real(kind=dp)   :: SIDpower_gal_quadrupole_at, arg_k
      real(kind=dp)   :: w_dd, w_dv, w_vv, pkdd, pkdv, pkvv
      real(kind=dp)   :: wskew_fmu2v1sq, wskew_f2mu2v1v2, wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq
      real(kind=dp)    :: wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f2mu2v1sq
      real(kind=dp)    :: wkurt_f4mu4v2sq, wkurt_f4mu6v2sq, wkurt_f4mu8v2sq
      real(kind=dp)   :: wkurt_f4mu2v2sq, wkurt_f3mu2v1v2

      if(do_init_model)then
        call init_model
      end if
 
      call RSD_Multipole(arg_k, w_dd, w_dv, w_vv, wskew_fmu2v1sq, wskew_f2mu2v1v2, &
               wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq, wkurt_f2mu2v1sq, &
               wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f4mu4v2sq, &
               wkurt_f4mu6v2sq, wkurt_f4mu8v2sq, wkurt_f4mu2v2sq, wkurt_f3mu2v1v2, 2)
 
      pkdd = gRPT_Pgg(arg_k)
      pkdv = gRPT_Pgv(arg_k)
      pkvv = gRPT_Pvv(arg_k)
 
      SIDpower_gal_quadrupole_at = w_dd*pkdd + w_dv*pkdv + w_vv*pkvv &
      + wskew_fmu2v1sq* RSDskew_fmu2v1sq(arg_k) &
      + wskew_f2mu2v1v2* RSDskew_f2mu2v1v2(arg_k) &
      + wskew_f2mu4v1v2* RSDskew_f2mu4v1v2(arg_k) &
      + wskew_f3mu4v2sq* RSDskew_f3mu4v2sq(arg_k) &
      + wskew_f3mu6v2sq* RSDskew_f3mu6v2sq(arg_k) &
      + wkurt_f2mu2v1sq* RSDkurt_f2mu2v1sq(arg_k) &
      + wkurt_f2mu4v1sq* RSDkurt_f2mu4v1sq(arg_k) &
      + wkurt_f3mu4v1v2* RSDkurt_f3mu4v1v2(arg_k) &
      + wkurt_f3mu6v1v2* RSDkurt_f3mu6v1v2(arg_k) &
      + wkurt_f4mu4v2sq* RSDkurt_f4mu4v2sq(arg_k) &
      + wkurt_f4mu6v2sq* RSDkurt_f4mu6v2sq(arg_k) &
      + wkurt_f4mu8v2sq* RSDkurt_f4mu8v2sq(arg_k) &
      + wkurt_f4mu2v2sq* RSDkurt_f4mu2v2sq(arg_k) &
      + wkurt_f3mu2v1v2* RSDkurt_f3mu2v1v2(arg_k)   

 
    end function SIDpower_gal_quadrupole_at
 
    function SIDpower_gal_hexadecapole_at(arg_k)
      implicit none
      real(kind=dp)  :: SIDpower_gal_hexadecapole_at, arg_k
      real(kind=dp)  :: w_dd, w_dv, w_vv, pkdd, pkdv, pkvv
      real(kind=dp)  :: wskew_fmu2v1sq, wskew_f2mu2v1v2, wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq
      real(kind=dp)  :: wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f2mu2v1sq
      real(kind=dp)  :: wkurt_f4mu4v2sq, wkurt_f4mu6v2sq, wkurt_f4mu8v2sq
      real(kind=dp)  :: wkurt_f4mu2v2sq, wkurt_f3mu2v1v2

      if(do_init_model)then
        call init_model
      end if
 
      call RSD_Multipole(arg_k, w_dd, w_dv, w_vv, wskew_fmu2v1sq, wskew_f2mu2v1v2, &
               wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq, wkurt_f2mu2v1sq, &
               wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f4mu4v2sq, &
               wkurt_f4mu6v2sq, wkurt_f4mu8v2sq, wkurt_f4mu2v2sq, wkurt_f3mu2v1v2, 4)
 
      pkdd = gRPT_Pgg(arg_k)
      pkdv = gRPT_Pgv(arg_k)
      pkvv = gRPT_Pvv(arg_k)
 
      SIDpower_gal_hexadecapole_at = w_dd*pkdd + w_dv*pkdv + w_vv*pkvv &
      + wskew_fmu2v1sq* RSDskew_fmu2v1sq(arg_k) &
      + wskew_f2mu2v1v2* RSDskew_f2mu2v1v2(arg_k) &
      + wskew_f2mu4v1v2* RSDskew_f2mu4v1v2(arg_k) &
      + wskew_f3mu4v2sq* RSDskew_f3mu4v2sq(arg_k) &
      + wskew_f3mu6v2sq* RSDskew_f3mu6v2sq(arg_k) &
      + wkurt_f2mu2v1sq* RSDkurt_f2mu2v1sq(arg_k) &
      + wkurt_f2mu4v1sq* RSDkurt_f2mu4v1sq(arg_k) &
      + wkurt_f3mu4v1v2* RSDkurt_f3mu4v1v2(arg_k) &
      + wkurt_f3mu6v1v2* RSDkurt_f3mu6v1v2(arg_k) &
      + wkurt_f4mu4v2sq* RSDkurt_f4mu4v2sq(arg_k) &
      + wkurt_f4mu6v2sq* RSDkurt_f4mu6v2sq(arg_k) &
      + wkurt_f4mu8v2sq* RSDkurt_f4mu8v2sq(arg_k) &
      + wkurt_f4mu2v2sq* RSDkurt_f4mu2v2sq(arg_k) &
      + wkurt_f3mu2v1v2* RSDkurt_f3mu2v1v2(arg_k)   

 
    end function SIDpower_gal_hexadecapole_at


    function SDcorr_power_gal_monopole_at(arg_k)  
      implicit none
      real(kind=dp)   SDcorr_power_gal_monopole_at, arg_k
      real(kind=dp)   pow_conv0

      call splint(q_table, rsd_convol0_table, rsd_convol0_d_table, nq_table, arg_k, pow_conv0)

      SDcorr_power_gal_monopole_at = pow_conv0
      
    end function


    function SDcorr_power_gal_quadrupole_at(arg_k)  
      implicit none
      real(kind=dp)   SDcorr_power_gal_quadrupole_at, arg_k
      real(kind=dp)   pow_conv2

      call splint(q_table, rsd_convol2_table, rsd_convol2_d_table, nq_table, arg_k, pow_conv2)

      SDcorr_power_gal_quadrupole_at = pow_conv2
      
    end function


    function SDcorr_power_gal_hexadecapole_at(arg_k)  
      implicit none
      real(kind=dp)   SDcorr_power_gal_hexadecapole_at, arg_k
      real(kind=dp)   pow_conv4

      call splint(q_table, rsd_convol4_table, rsd_convol4_d_table, nq_table, arg_k, pow_conv4)

      SDcorr_power_gal_hexadecapole_at = pow_conv4
      
    end function



    function gRPTpower_gal_monopole_at(arg_k)
      implicit none
      real(kind=dp)   :: gRPTpower_gal_monopole_at,arg_k
      
      if(do_init_model)then
        call init_model
      end if
  
      gRPTpower_gal_monopole_at = SIDpower_gal_monopole_at(arg_k) + &
                                0._dp* SDcorr_power_gal_monopole_at(arg_k)  
      
    
    end function gRPTpower_gal_monopole_at
  
    function gRPTpower_gal_quadrupole_at(arg_k)
      implicit none
      real(kind=dp)   :: gRPTpower_gal_quadrupole_at, arg_k
      
      if(do_init_model)then
        call init_model
      end if
 
      if (f>0._dp) then  
      gRPTpower_gal_quadrupole_at = SIDpower_gal_quadrupole_at(arg_k) + &
                                  0._dp*  SDcorr_power_gal_quadrupole_at(arg_k)  
      else !return Pgd
         gRPTpower_gal_quadrupole_at = gRPT_Pgd(arg_k) 
      endif
      
    end function gRPTpower_gal_quadrupole_at
 
    function gRPTpower_gal_hexadecapole_at(arg_k)
      implicit none
      real(kind=dp)  :: gRPTpower_gal_hexadecapole_at, arg_k
      
      if(do_init_model)then
        call init_model
      end if
 
      gRPTpower_gal_hexadecapole_at = SIDpower_gal_hexadecapole_at(arg_k) + &
                                    0._dp*  SDcorr_power_gal_hexadecapole_at(arg_k)
 
    end function gRPTpower_gal_hexadecapole_at

!!!!!!!!!!!!!!!!!!!!!!!!!
! THE FOLLOWING USED TO BE RSD.F90
!!!!!!!!!!!!!!!!!!!!!!!!!

!************************************************************
! Write down the 5 P12RS integrals that decompose into MC and PROP
!************************************************************
! in the following instead of Pdv we should use dewiggled linear power
! remember to do the same in Pg2 and Pb1g2 prop integrals

   function Pk_fmu2v1sqb1sq(arg_k)
     implicit none
     real(kind=dp)   :: Pk_fmu2v1sqb1sq,arg_k

     Pk_fmu2v1sqb1sq   =  Pk_fmu2v1sqb1sq_prop_old(arg_k) + Pk_fmu2v1sqb1sq_mc(arg_k)

   end function Pk_fmu2v1sqb1sq

   function Pk_f3mu4v2sq(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f3mu4v2sq,arg_k

     Pk_f3mu4v2sq   =  Pk_f3mu4v2sq_prop(arg_k) + Pk_f3mu4v2sq_mc(arg_k)

   end function Pk_f3mu4v2sq

   function Pk_f3mu6v2sq(arg_k) 
     implicit none
     real(kind=dp)   :: Pk_f3mu6v2sq,arg_k

     Pk_f3mu6v2sq   =   Pk_f3mu6v2sq_prop(arg_k) + Pk_f3mu6v2sq_mc(arg_k)

   end function Pk_f3mu6v2sq

   function Pk_f2mu2v1v2b1(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f2mu2v1v2b1,arg_k

     Pk_f2mu2v1v2b1   =  Pk_f2mu2v1v2b1_prop(arg_k) + Pk_f2mu2v1v2b1_mc(arg_k)

   end function Pk_f2mu2v1v2b1
   
   function Pk_f2mu4v1v2b1(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f2mu4v1v2b1,arg_k

     Pk_f2mu4v1v2b1   =   Pk_f2mu4v1v2b1_prop_old(arg_k) + Pk_f2mu4v1v2b1_mc(arg_k)

   end function Pk_f2mu4v1v2b1


!************************************************************
! 5 PROP integrals for P12RS
!************************************************************

   function Pk_fmu2v1sqb1sq_prop_old(arg_k)
     implicit none
     real(kind=dp)   :: Pk_fmu2v1sqb1sq_prop_old,arg_k,epsabs,epsrel,a,b
     integer, parameter :: limit=1000, key = 2 !! local integration rule !!
     real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
     integer  :: neval,ier,last,iord(limit)
     real(kind=dp)   :: ans

     k = arg_k
     a = kminRSD      ! lower limit of integration
     b = kmaxRSD      ! upper limit of integration
     epsrel = relaccRSD/10._dp
     epsabs = absaccRSD/10._dp
     ans = 0._dp

     call qage (fmu2v1sqb1sq_dq,a,b,epsabs,epsrel,key,limit,ans, &
     abserr,neval,ier,alist,blist,rlist,elist,iord,last)
     Pk_fmu2v1sqb1sq_prop_old = 4._dp*pi*ans/fac_norm*Plin_DW(k) !matterpowerat(k) !*gRPT_Pdv(k)
     if(ier /= 0)  call warning('Pk_fmu2v1sqb1sq_prop_old',ier)  

   end function Pk_fmu2v1sqb1sq_prop_old

   function fmu2v1sqb1sq_dq(q)        
     implicit none
     real(kind=dp)   :: fmu2v1sqb1sq_dq,q,kernel

     kernel = -(76._dp*k**5*q-96._dp*k**3*q**3+36._dp*k*q**5+9._dp*(k**2-q**2)**3*dlog((k + q)**2/(k - q)**2))/ &
       (168._dp*k**3*q**3) 
     fmu2v1sqb1sq_dq = kernel*q**2 *Plin_DW(q) !matterpowerat(q)

   end function fmu2v1sqb1sq_dq
   
   
   function Pk_fmu2v1sqb1sq_prop(arg_k)
     implicit none
     real(kind=dp)   :: Pk_fmu2v1sqb1sq_prop,arg_k,epsabs,epsrel,a,b
     integer, parameter :: limit=1000, key = 2 !! local integration rule !!
     real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
     integer  :: neval,ier,last,iord(limit)
     real(kind=dp)   :: ans,ans1,ans2

     k = arg_k
     a = kminRSD      ! lower limit of integration
     b = kmaxRSD      ! upper limit of integration
     epsrel = relaccRSD/10._dp
     epsabs = absaccRSD/10._dp
     ans = 0._dp

     call qage (fmu2v1sqb1sq_dd,a,b,epsabs,epsrel,key,limit,ans, &
     abserr,neval,ier,alist,blist,rlist,elist,iord,last)
     if(ier /= 0)  call warning('Pk_fmu2v1sqb1sq_dd',ier)  
     ans1 = 4._dp*pi*ans*gRPT_Pdd(k)/fac_norm !*matterpowerat(k)/fac_norm !

     call qage (fmu2v1sqb1sq_dv,a,b,epsabs,epsrel,key,limit,ans, &
     abserr,neval,ier,alist,blist,rlist,elist,iord,last)
     if(ier /= 0)  call warning('Pk_fmu2v1sqb1sq_dv',ier)  
     ans2 = 4._dp*pi*ans*gRPT_Pdv(k)/fac_norm !*matterpowerat(k)/fac_norm !
     
     Pk_fmu2v1sqb1sq_prop = ans1 + ans2

   end function Pk_fmu2v1sqb1sq_prop


   function fmu2v1sqb1sq_dd(q)        !ExorciseWiggles.nb
     implicit none
     real(kind=dp)   :: fmu2v1sqb1sq_dd,q,kernel

     kernel =  (-2._dp*k**2)/(3._dp*q**2)
     fmu2v1sqb1sq_dd = kernel*q**2*matterpowerat(q) !*gRPT_Pdv(q)

   end function fmu2v1sqb1sq_dd
   
   function fmu2v1sqb1sq_dv(q)         !ExorciseWiggles.nb
     implicit none
     real(kind=dp)   :: fmu2v1sqb1sq_dv,q,kernel

     kernel = (4._dp*(3._dp*k**5*q+8._dp*k**3*q**3-3._dp*k*q**5)-3._dp*(k**2-q**2)**3 &
     *dlog((k + q)**2/(k - q)**2))/(56._dp*k**3*q**3)
     fmu2v1sqb1sq_dv = kernel*q**2*matterpowerat(q) !*gRPT_Pdv(q)

   end function fmu2v1sqb1sq_dv
   
   
!******************************

   function Pk_f3mu4v2sq_prop(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f3mu4v2sq_prop,arg_k,epsabs,epsrel,a,b
     integer, parameter :: limit=1000, key = 2 !! local integration rule !!
     real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
     integer  :: neval,ier,last,iord(limit)
     real(kind=dp)   :: ans

     k = arg_k
     a = kminRSD      ! lower limit of integration
     b = kmaxRSD      ! upper limit of integration
     epsrel = relaccRSD/10._dp
     epsabs = absaccRSD/10._dp
     ans = 0._dp

     call qage (f3mu4v2sq_dq,a,b,epsabs,epsrel,key,limit,ans, &
     abserr,neval,ier,alist,blist,rlist,elist,iord,last)
     Pk_f3mu4v2sq_prop = 4._dp*pi*ans/fac_norm*Plin_DW(k) !matterpowerat(k) !gRPT_Pvv(k)
     if(ier /= 0)  call warning('Pk_f3mu4v2sq_prop',ier)  

   end function Pk_f3mu4v2sq_prop

   function f3mu4v2sq_dq(q)         !ExorciseWiggles.nb
     implicit none
     real(kind=dp)   :: f3mu4v2sq_dq,q,kernel

     kernel = -(4._dp*k*q*(-9._dp*k**6+33._dp*k**4*q**2+33._dp*k**2*q**4-9._dp*q**6) + &
          9._dp*(k**2-q**2)**4*dlog((k + q)**2/(k - q)**2))/(672._dp*k**3*q**5)
     f3mu4v2sq_dq = kernel*q**2 *Plin_DW(q) !matterpowerat(q) !*gRPT_Pvv(q)*q**2

   end function f3mu4v2sq_dq

!******************************

  function Pk_f3mu6v2sq_prop(arg_k)  
     implicit none
     real(kind=dp)   :: Pk_f3mu6v2sq_prop,arg_k,epsabs,epsrel,a,b
     integer, parameter :: limit=1000, key = 2 !! local integration rule !!
     real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
     integer  :: neval,ier,last,iord(limit)
     real(kind=dp)   :: ans

     k = arg_k
     a = kminRSD      ! lower limit of integration
     b = kmaxRSD      ! upper limit of integration
     epsrel = relaccRSD/10._dp
     epsabs = absaccRSD/10._dp
     ans = 0._dp

     call qage (f3mu6v2sq_dq,a,b,epsabs,epsrel,key,limit,ans, &
     abserr,neval,ier,alist,blist,rlist,elist,iord,last)
     Pk_f3mu6v2sq_prop = 4._dp*pi*ans/fac_norm*Plin_DW(k) !matterpowerat(k) !*gRPT_Pvv(k)
     if(ier /= 0)  call warning('Pk_f3mu6v2sq_prop',ier)  
     
!     write(*,*)k,gRPT_Pdd(k),gRPT_Pdv(k),gRPT_Pvv(k),matterpowerat(k)

   end function Pk_f3mu6v2sq_prop

   function f3mu6v2sq_dq(q)         !ExorciseWiggles.nb
     implicit none
     real(kind=dp)   :: f3mu6v2sq_dq,q,kernel

     kernel =  (4._dp*(9._dp*k**7*q - 109._dp*k**5*q**3 + 63._dp*k**3*q**5 - 27._dp*k*q**7) - &
         9._dp*(k**2 - q**2)**3*(k**2 + 3._dp*q**2)*dlog((k + q)**2/(k - q)**2))/(672._dp*k**3*q**5)
     f3mu6v2sq_dq = kernel*q**2 *Plin_DW(q) !matterpowerat(q) !*gRPT_Pvv(q)*q**2

   end function f3mu6v2sq_dq

!******************************

  function Pk_f2mu2v1v2b1_prop(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f2mu2v1v2b1_prop,arg_k,epsabs,epsrel,a,b
     integer, parameter :: limit=1000, key = 2 !! local integration rule !!
     real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
     integer  :: neval,ier,last,iord(limit)
     real(kind=dp)   :: ans

     k = arg_k
     a = kminRSD      ! lower limit of integration
     b = kmaxRSD      ! upper limit of integration
     epsrel = relaccRSD/10._dp
     epsabs = absaccRSD/10._dp
     ans = 0._dp

     call qage (f2mu2v1v2b1_dq,a,b,epsabs,epsrel,key,limit,ans, &
     abserr,neval,ier,alist,blist,rlist,elist,iord,last)
     Pk_f2mu2v1v2b1_prop = 4._dp*pi*ans/fac_norm*Plin_DW(k) !matterpowerat(k) !*gRPT_Pdv(k)
     if(ier /= 0)  call warning('Pk_f2mu2v1v2b1_prop',ier)  

   end function Pk_f2mu2v1v2b1_prop


   function f2mu2v1v2b1_dq(q)        !ExorciseWiggles.nb
     implicit none
     real(kind=dp)   :: f2mu2v1v2b1_dq,q,kernel

     kernel = -(4._dp*k*q*(-9._dp*k**6 + 33._dp*k**4*q**2 + 33._dp*k**2*q**4 - 9._dp*q**6) + &
           9._dp*(k**2 - q**2)**4*dlog((k + q)**2/(k - q)**2))/(672._dp*k**3*q**5)
     f2mu2v1v2b1_dq = kernel*q**2 *Plin_DW(q) !matterpowerat(q)

   end function f2mu2v1v2b1_dq


!******************************

  function Pk_f2mu4v1v2b1_prop_old(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f2mu4v1v2b1_prop_old,arg_k,epsabs,epsrel,a,b
     integer, parameter :: limit=1000, key = 2 !! local integration rule !!
     real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
     integer  :: neval,ier,last,iord(limit)
     real(kind=dp)   :: ans

     k = arg_k
     a = kminRSD      ! lower limit of integration
     b = kmaxRSD      ! upper limit of integration
     epsrel = relaccRSD/10._dp
     epsabs = absaccRSD/10._dp
     ans = 0._dp

     call qage (f2mu4v1v2b1_dq,a,b,epsabs,epsrel,key,limit,ans, &
     abserr,neval,ier,alist,blist,rlist,elist,iord,last)
     Pk_f2mu4v1v2b1_prop_old = 4._dp*pi*ans/fac_norm*Plin_DW(k) !matterpowerat(k) !*gRPT_Pdv(k) !
     if(ier /= 0)  call warning('Pk_f2mu4v1v2b1_prop_old',ier)  

   end function Pk_f2mu4v1v2b1_prop_old

   function f2mu4v1v2b1_dq(q)        
     implicit none
     real(kind=dp)   :: f2mu4v1v2b1_dq,q,kernel

     kernel =   (4._dp*k*q*(9._dp*k**6 - 185._dp*k**4*q**2 + 159._dp*k**2*q**4 - 63._dp*q**6) - &
         9._dp*(k**2 - q**2)**3*(k**2 + 7._dp*q**2)*dlog((k + q)**2/(k - q)**2))/(672._dp*k**3*q**5)
     f2mu4v1v2b1_dq = kernel*q**2 *Plin_DW(q) !matterpowerat(q)

     end function f2mu4v1v2b1_dq


   function Pk_f2mu4v1v2b1_prop(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f2mu4v1v2b1_prop,arg_k,epsabs,epsrel,a,b
     integer, parameter :: limit=1000, key = 2 !! local integration rule !!
     real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
     integer  :: neval,ier,last,iord(limit)
     real(kind=dp)   :: ans,ans1,ans2,ans3

     k = arg_k
     a = kminRSD      ! lower limit of integration
     b = kmaxRSD      ! upper limit of integration
     epsrel = relaccRSD/10._dp
     epsabs = absaccRSD/10._dp
     ans = 0._dp

! use 2 previous integrals
     call qage (fmu2v1sqb1sq_dd,a,b,epsabs,epsrel,key,limit,ans, &
     abserr,neval,ier,alist,blist,rlist,elist,iord,last)
     if(ier /= 0)  call warning('Pk_fmu2v1sqb1sq_dd',ier)  
     ans1 = 4._dp*pi*ans*gRPT_Pdv(k)/fac_norm !*matterpowerat(k)/fac_norm !

     call qage (fmu2v1sqb1sq_dv,a,b,epsabs,epsrel,key,limit,ans, &
     abserr,neval,ier,alist,blist,rlist,elist,iord,last)
     if(ier /= 0)  call warning('Pk_fmu2v1sqb1sq_dv',ier)  
     ans2 = 4._dp*pi*ans*gRPT_Pvv(k)/fac_norm !*matterpowerat(k)/fac_norm !

     call qage (f2mu4v1v2b1_dv,a,b,epsabs,epsrel,key,limit,ans, &
     abserr,neval,ier,alist,blist,rlist,elist,iord,last)
     if(ier /= 0)  call warning('Pk_f2mu4v1v2b1_dv',ier)  
     ans3 = 4._dp*pi*ans*gRPT_Pdv(k)/fac_norm !*matterpowerat(k)/fac_norm !
     
     Pk_f2mu4v1v2b1_prop = ans1+ans2+ans3

   end function Pk_f2mu4v1v2b1_prop


   function f2mu4v1v2b1_dv(q)         !ExorciseWiggles.nb
     implicit none
     real(kind=dp)   :: f2mu4v1v2b1_dv,q,kernel

     kernel = (4._dp*(9._dp*k**7*q - 109._dp*k**5*q**3 + 63._dp*k**3*q**5 - 27._dp*k*q**7) - & 
      9._dp*(k**2 - q**2)**3*(k**2 + 3._dp*q**2)*dlog((k + q)**2/(k - q)**2))/(672._dp*k**3*q**5)
     f2mu4v1v2b1_dv = kernel*q**2*matterpowerat(q) !*gRPT_Pvv(q)*q**2

   end function f2mu4v1v2b1_dv


!************************************************************
! now the 5 MC + 6 MC (no PROP counterpart) integrals for P12RS
!************************************************************


   function Pk_fmu2v1sqb1sq_mc(arg_k)
     implicit none
     real(kind=dp)   :: Pk_fmu2v1sqb1sq_mc,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(fmu2v1sqb1sq_int,method,epsrel,epsabs,ans,error,prob)
     Pk_fmu2v1sqb1sq_mc = ans(1)/fac_norm

   end function Pk_fmu2v1sqb1sq_mc
   
   subroutine fmu2v1sqb1sq_int(ndim,x,ncomp,f)  !MC counterpart
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
!        y    = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = Plin_DW(q) !matterpowerat(q)
        powp = Plin_DW(p) !matterpowerat(p)        
!        kernel = (-2._dp*(-k**2 + p*(p + q*y))*(-10._dp*p*(p + q*y)**2 + k**2*(10._dp*p + 7._dp*q*y)))/(7._dp*p*q**4)
        kernel = ((k**2 - p**2 + q**2)*(2*k**4 - 5*(p**2 - q**2)**2 + 3*k**2*(p**2 + q**2)))/(14.*p**2*q**4)
!        kernel = ((-(p**2 - q**2)**2 + k**2*(p**2 + q**2))* &
!         (2*k**4 - 5*(p**2 - q**2)**2 + 3*k**2*(p**2 + q**2)))/(28.*p**4*q**4)
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin) !*dexp(y)*(1._dp-y)
     else
        f(1)=0._dp
     endif

   end subroutine fmu2v1sqb1sq_int


!******************************

   function Pk_f3mu4v2sq_mc(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f3mu4v2sq_mc,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f3mu4v2sq_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f3mu4v2sq_mc = ans(1)/fac_norm

   end function Pk_f3mu4v2sq_mc

   subroutine f3mu4v2sq_int(ndim,x,ncomp,f) !MC counterpart
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
!        y    = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = Plin_DW(q) !matterpowerat(q)
        powp = Plin_DW(p) !matterpowerat(p)        
!        kernel = (-18._dp*p**2*(p + q*y)**5 + k**6*(6._dp*p + 7._dp*q*y) + 3._dp*k**2*p*(p + q*y)**3*(14._dp*p + 9._dp*q*y) - &
!         k**4*(30._dp*p**3 + 70._dp*p**2*q*y + 47._dp*p*q**2*y**2 + 7._dp*q**3*y**3))/(7._dp*k**2*p*q**4)
        kernel =((k - p - q)*(k + p - q)*(k - p + q)*(k + p + q)*(k**2 + 3*(p - q)*(p + q))* & 
         (4*k**4 - 3*(p**2 - q**2)**2 - k**2*(p**2 + q**2)))/(112.*k**2*p**4*q**4)
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f3mu4v2sq_int

!******************************

   function Pk_f3mu6v2sq_mc(arg_k) 
     implicit none
     real(kind=dp)   :: Pk_f3mu6v2sq_mc,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f3mu6v2sq_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f3mu6v2sq_mc = ans(1)/fac_norm

   end function Pk_f3mu6v2sq_mc

   subroutine f3mu6v2sq_int(ndim,x,ncomp,f) !MC counterpart
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
!        y    = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = Plin_DW(q) !matterpowerat(q)
        powp = Plin_DW(p) !matterpowerat(p)        
!        kernel = -(-30._dp*p**2*(p + q*y)**5 + k**6*(6._dp*p + 7._dp*q*y) + k**2*p*(p + q*y)**3*(66._dp*p + 53._dp*q*y) - &
!          3._dp*k**4*(14._dp*p**3 + 36._dp*p**2*q*y + 29._dp*p*q**2*y**2 + 7._dp*q**3*y**3))/(7._dp*k**2*p*q**4)  
        kernel = ((4*k**4 - 3*(p**2 - q**2)**2 - k**2*(p**2 + q**2))* &
        (k**6 - 5*(p**2 - q**2)**3 + k**4*(p**2 + 3*q**2) + 3*k**2*(p**4 + 2*p**2*q**2 - 3*q**4)))/ &
        (112.*k**2*p**4*q**4)
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f3mu6v2sq_int

!******************************

   function Pk_f2mu4v1v2b1_mc(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f2mu4v1v2b1_mc,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f2mu4v1v2b1_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f2mu4v1v2b1_mc = ans(1)/fac_norm

   end function Pk_f2mu4v1v2b1_mc

   subroutine f2mu4v1v2b1_int(ndim,x,ncomp,f) !MC counterpart
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
!        y    = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = Plin_DW(q) !matterpowerat(q)
        powp = Plin_DW(p) !matterpowerat(p)        
!        kernel = (50._dp*p**2*(p + q*y)**5 + k**6*(2._dp*p + 7._dp*q*y) - k**2*p*(p + q*y)**3*(98._dp*p + 65._dp*q*y) + &
!          k**4*(46._dp*p**3 + 102._dp*p**2*q*y + 77._dp*p*q**2*y**2 + 21._dp*q**3*y**3))/(7._dp*k**2*p*q**4)  
        kernel = (2*k**10 + 25*(p**2 - q**2)**5 + k**8*(37*p**2 + 9*q**2) - 6*k**2*(p**2 - q**2)**3*(p**2 + 10*q**2) - &
         2*k**6*(18*p**4 - 29*p**2*q**2 + 7*q**4) - 2*k**4*(11*p**6 - 50*p**4*q**2 + 23*p**2*q**4 + 16*q**6)) &
         /(112.*k**2*p**4*q**4)
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f2mu4v1v2b1_int

!******************************

   function Pk_f2mu2v1v2b1_mc(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f2mu2v1v2b1_mc,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f2mu2v1v2b1_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f2mu2v1v2b1_mc = ans(1)/fac_norm

   end function Pk_f2mu2v1v2b1_mc

   subroutine f2mu2v1v2b1_int(ndim,x,ncomp,f) !MC counterpart
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
!        y    = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = Plin_DW(q) !matterpowerat(q)
        powp = Plin_DW(p) !matterpowerat(p)        
!        kernel =  (-30._dp*p**2*(p + q*y)**5 + k**6*(10._dp*p + 7._dp*q*y) + k**2*p*(p + q*y)**3*(70._dp*p + 31._dp*q*y) - &
!          k**4*(50._dp*p**3 + 98._dp*p**2*q*y + 55._dp*p*q**2*y**2 + 7._dp*q**3*y**3))/(7._dp*k**2*p*q**4)
        kernel = ((k - p - q)*(k + p - q)*(k - p + q)*(k + p + q)*(k**2 + 3*(p - q)*(p + q))* &
          (2*k**4 - 5*(p**2 - q**2)**2 + 3*k**2*(p**2 + q**2)))/(112.*k**2*p**4*q**4) 
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f2mu2v1v2b1_int

!******************************

   function Pk_fmu2v1sqb1b2(arg_k)
     implicit none
     real(kind=dp)   :: Pk_fmu2v1sqb1b2,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(fmu2v1sqb1b2_int,method,epsrel,epsabs,ans,error,prob)
     Pk_fmu2v1sqb1b2 = 1._dp*ans(1)/fac_norm

   end function Pk_fmu2v1sqb1b2

   subroutine fmu2v1sqb1b2_int(ndim,x,ncomp,f) !use Plin here
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
!        y    = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = matterpowerat(q) !Plin_DW(q) !
        powp = matterpowerat(p) !Plin_DW(p) !      
!        kernel = (2._dp*(-k**2 + p*(p + q*y))*(-k**2 + p*(p + 2._dp*q*y)))/q**4
        kernel = 1._dp + ((k - p)*(k + p))/q**2
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine fmu2v1sqb1b2_int

!******************************

   function Pk_fmu2v1sqb1g2(arg_k)
     implicit none
     real(kind=dp)   :: Pk_fmu2v1sqb1g2,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(fmu2v1sqb1g2_int,method,epsrel,epsabs,ans,error,prob)
     Pk_fmu2v1sqb1g2 = 1._dp*ans(1)/fac_norm

   end function Pk_fmu2v1sqb1g2

   subroutine fmu2v1sqb1g2_int(ndim,x,ncomp,f) !use Plin here
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
!        y    = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = matterpowerat(q) !Plin_DW(q) !
        powp = matterpowerat(p) !Plin_DW(p) !        
!        kernel = (-4._dp*(k**4 + p*(p + q*y)**3 - k**2*(2._dp*p**2 + 3._dp*p*q*y + q**2*y**2)))/q**4
        kernel = ((k - p - q)*(k + p - q)*(k - p + q)*(k + p + q)*(k**2 - p**2 + q**2))/(2.*p**2*q**4)
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine fmu2v1sqb1g2_int

!******************************

   function Pk_f2mu4v1v2b2(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f2mu4v1v2b2,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f2mu4v1v2b2_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f2mu4v1v2b2 = 1._dp*ans(1)/fac_norm

   end function Pk_f2mu4v1v2b2

   subroutine f2mu4v1v2b2_int(ndim,x,ncomp,f) !use Plin here
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
!        y    = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = matterpowerat(q) !Plin_DW(q) !
        powp = matterpowerat(p) !Plin_DW(p) !      
!        kernel = (-k**6 + 5._dp*p**2*(p + q*y)**3*(p + 2._dp*q*y) + k**4*(7._dp*p**2 + 11._dp*p*q*y + 3._dp*q**2*y**2) - & 
!          k**2*p*(11._dp*p**3 + 36._dp*p**2*q*y + 36._dp*p*q**2*y**2 + 11._dp*q**3*y**3))/(k**2*q**4)  
        kernel = (k**6 - 5*(p**2 - q**2)**3 + k**4*(p**2 + 3*q**2) + 3*k**2*(p**4 + 2*p**2*q**2 - 3*q**4))/(8.*k**2*p**2*q**2)
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f2mu4v1v2b2_int

!******************************

   function Pk_f2mu4v1v2g2(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f2mu4v1v2g2,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f2mu4v1v2g2_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f2mu4v1v2g2 = 1._dp*ans(1)/fac_norm

   end function Pk_f2mu4v1v2g2

   subroutine f2mu4v1v2g2_int(ndim,x,ncomp,f) !use Plin here
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
!        y    = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = matterpowerat(q) !Plin_DW(q) !
        powp = matterpowerat(p) !Plin_DW(p) !       
!        kernel =  (2._dp*(k**6 - 5._dp*p*(p + q*y)**5 + k**2*(p + q*y)**3*(11._dp*p + 3._dp*q*y) - &
!            k**4*(7._dp*p**2 + 11._dp*p*q*y + 4._dp*q**2*y**2)))/(k**2*q**4)
        kernel = ((k - p - q)*(k + p - q)*(k - p + q)*(k + p + q)* &
         (k**6 - 5*(p**2 - q**2)**3 + k**4*(p**2 + 3*q**2) + 3*k**2*(p**4 + 2*p**2*q**2 - 3*q**4)))/ &
         (16.*k**2*p**4*q**4)    
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f2mu4v1v2g2_int

!******************************

   function Pk_f2mu2v1v2b2(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f2mu2v1v2b2,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f2mu2v1v2b2_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f2mu2v1v2b2 = 1._dp*ans(1)/fac_norm

   end function Pk_f2mu2v1v2b2

   subroutine f2mu2v1v2b2_int(ndim,x,ncomp,f) !use Plin here
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
!        y    = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = matterpowerat(q) !Plin_DW(q) !
        powp = matterpowerat(p) !Plin_DW(p) !       
!        kernel = (k**6 - 3._dp*p**2*(p + q*y)**3*(p + 2._dp*q*y) - k**4*(5._dp*p**2 + 7._dp*p*q*y + q**2*y**2) + &
!          k**2*p*(7._dp*p**3 + 22._dp*p**2*q*y + 20._dp*p*q**2*y**2 + 5._dp*q**3*y**3))/(k**2*q**4)
        kernel = ((k-p-q)**2*(k+p-q)**2*(k-p+q)**2*(k+p+q)**2*(k**2 + 3*(p - q)*(p + q))**2)/ &
        (64.*k**4*p**4*q**4)  
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f2mu2v1v2b2_int

!******************************

   function Pk_f2mu2v1v2g2(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f2mu2v1v2g2,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f2mu2v1v2g2_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f2mu2v1v2g2 = 1._dp*ans(1)/fac_norm

   end function Pk_f2mu2v1v2g2

   subroutine f2mu2v1v2g2_int(ndim,x,ncomp,f) !use Plin here
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
!        y    = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = matterpowerat(q) !Plin_DW(q) !
        powp = matterpowerat(p) !Plin_DW(p) !       
!        kernel = (-2._dp*(k**2 - 3._dp*p*(p + q*y))*(k**2 - (p + q*y)**2)**2)/(k**2*q**4)
        kernel = ((k**2 + 3*p**2 - 3*q**2)*(k**4 + (p**2 - q**2)**2 - 2*k**2*(p**2 + q**2))**2)/ &
        (16.*k**2*p**4*q**4)
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f2mu2v1v2g2_int


!************************************************************
! Now all 9 P22RS (MC) integrals: here I'm using gRPT Pdd,Pdv,Pvv as needed
!************************************************************


   function Pk_f2mu2v1sqb1sq(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f2mu2v1sqb1sq,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f2mu2v1sqb1sq_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f2mu2v1sqb1sq = 1._dp*ans(1)/fac_norm

   end function Pk_f2mu2v1sqb1sq

   subroutine f2mu2v1sqb1sq_int(ndim,x,ncomp,f) !RSDTensorDecomposition.nb 
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        powq = matterpowerat(q) ! gRPT_Pdv(q) !Plin_DW(q) !  
        powp = matterpowerat(p) ! gRPT_Pdv(p) !Plin_DW(p) !       
        kernel = ((-k + p - q)*(k + p - q)*(-k + p + q)*(k + p + q))/(8._dp*p**2*q**2)
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f2mu2v1sqb1sq_int


!******************************

   function Pk_f2mu4v1sqb1sq(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f2mu4v1sqb1sq,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f2mu4v1sqb1sq_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f2mu4v1sqb1sq = 1._dp*ans(1)/fac_norm

   end function Pk_f2mu4v1sqb1sq

   subroutine f2mu4v1sqb1sq_int(ndim,x,ncomp,f) !RSDTensorDecomposition.nb 
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        powq = matterpowerat(q) !  gRPT_Pdv(q) !Plin_DW(q) !
        powp = matterpowerat(p) !gRPT_Pdv(p) !Plin_DW(p) !        
        kernel = (k**4 - 3._dp*(p**2 - q**2)**2 + 2._dp*k**2*(p**2 + q**2))/(8._dp*p**2*q**2)
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f2mu4v1sqb1sq_int

!******************************

   function Pk_f3mu4v1v2b1(arg_k) !this one is quite wiggly (but not noisy), not sure why
     implicit none
     real(kind=dp)   :: Pk_f3mu4v1v2b1,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f3mu4v1v2b1_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f3mu4v1v2b1 = 1._dp*ans(1)/fac_norm

   end function Pk_f3mu4v1v2b1

   subroutine f3mu4v1v2b1_int(ndim,x,ncomp,f) !RSDTensorDecomposition.nb 
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp
     real(kind=dp)     :: pdvq,pvvp,pvvq,pdvp      

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        pdvq = matterpowerat(q) !gRPT_Pdv(q) !Plin_DW(q) !
        pvvp = matterpowerat(p) !gRPT_Pvv(p) !Plin_DW(p) !       
        pvvq = matterpowerat(q) !gRPT_Pvv(q) !Plin_DW(q) !
        pdvp = matterpowerat(p) !gRPT_Pdv(p) !Plin_DW(p) !       
        kernel = (3._dp*(k**4 + (p**2 - q**2)**2 - 2._dp*k**2*(p**2 + q**2))* & 
        (q**2*(k**4 + 2._dp*k**2*(p**2 - 3._dp*q**2) + 5._dp*(p**2 - q**2)**2)*Pdvq*Pvvp + &
        p**2*(k**4 + 5*(p**2 - q**2)**2 + k**2*(-6*p**2 + 2*q**2))*Pdvp*Pvvq))/ &
        (64.*k**2*p**4*q**4)
        f(1)   = (2._dp*pi/k)*q*p*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f3mu4v1v2b1_int

!******************************

   function Pk_f3mu6v1v2b1(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f3mu6v1v2b1,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f3mu6v1v2b1_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f3mu6v1v2b1 = 1._dp*ans(1)/fac_norm

   end function Pk_f3mu6v1v2b1

   subroutine f3mu6v1v2b1_int(ndim,x,ncomp,f) !RSDTensorDecomposition.nb 
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp
     real(kind=dp)     :: pdvq,pvvp,pvvq,pdvp      

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        pdvq = matterpowerat(q) !gRPT_Pdv(q) !Plin_DW(q) !
        pvvp = matterpowerat(p) !  gRPT_Pvv(p) !Plin_DW(p) !    
        pvvq = matterpowerat(q) !gRPT_Pvv(q) !Plin_DW(q) !
        pdvp = matterpowerat(p) !gRPT_Pdv(p) !Plin_DW(p) !     
        kernel = (q**2*(5._dp*k**8 - 35._dp*(p**2 - q**2)**4 + 4._dp*k**6*(p**2 + 5*q**2) + &
         20._dp*k**2*(p**2-q**2)**2*(p**2+5._dp*q**2)+6._dp*k**4*(p**4+6._dp*p**2*q**2-15._dp*q**4))&
        *Pdvq*Pvvp + p**2*(5._dp*k**8-35._dp*(p**2-q**2)**4+ 4._dp*k**6*(5._dp*p**2 + q**2) + &
        20._dp*k**2*(p**2-q**2)**2*(5._dp*p**2+q**2)+6._dp*k**4*(-15._dp*p**4+6._dp*p**2*q**2+q**4))&
        *Pdvp*Pvvq)/(128._dp*k**2*p**4*q**4)         
        f(1)   = (2._dp*pi/k)*q*p*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f3mu6v1v2b1_int

!******************************


   function Pk_f4mu4v2sq(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f4mu4v2sq,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f4mu4v2sq_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f4mu4v2sq = 1._dp*ans(1)/fac_norm

   end function Pk_f4mu4v2sq

   subroutine f4mu4v2sq_int(ndim,x,ncomp,f) !RSDTensorDecomposition.nb 
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        powq = matterpowerat(q) !gRPT_Pvv(q) !Plin_DW(q) !
        powp = matterpowerat(p) !    gRPT_Pvv(p) !Plin_DW(p) !  
        kernel = (3._dp*(k**4 + (p**2 - q**2)**2 - 2._dp*k**2*(p**2 + q**2))**2* &
         (k**4 - 35._dp*(p**2 - q**2)**2 + 10._dp*k**2*(p**2 + q**2)))/(1024._dp*k**4*p**4*q**4) 
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f4mu4v2sq_int

!******************************

   function Pk_f4mu6v2sq(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f4mu6v2sq,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f4mu6v2sq_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f4mu6v2sq = 1._dp*ans(1)/fac_norm

   end function Pk_f4mu6v2sq

   subroutine f4mu6v2sq_int(ndim,x,ncomp,f) !RSDTensorDecomposition.nb 
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        powq = matterpowerat(q) !gRPT_Pvv(q) !Plin_DW(q) !
        powp = matterpowerat(p) !gRPT_Pvv(p) !Plin_DW(p) !        
        kernel = (3._dp*(k**12 + 105._dp*(p**2 - q**2)**6 + 2._dp*k**10*(p**2 + q**2) - &
         350._dp*k**2*(p**2 - q**2)**4*(p**2 + q**2) + k**8*(23._dp*p**4 + 2._dp*p**2*q**2 + 23._dp*q**4) + &
         5._dp*k**4*(p**2 - q**2)**2*(83._dp*p**4 + 74._dp*p**2*q**2 + 83._dp*q**4) - &
         4._dp*k**6*(49._dp*p**6 - 9._dp*p**4*q**2 - 9._dp*p**2*q**4 + 49._dp*q**6)))/ &
         (1024._dp*k**4*p**4*q**4)    
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f4mu6v2sq_int

!******************************

   function Pk_f4mu8v2sq(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f4mu8v2sq,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f4mu8v2sq_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f4mu8v2sq = 1._dp*ans(1)/fac_norm

   end function Pk_f4mu8v2sq

   subroutine f4mu8v2sq_int(ndim,x,ncomp,f) !RSDTensorDecomposition.nb 
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        powq = matterpowerat(q) !gRPT_Pvv(q) !Plin_DW(q) !
        powp = matterpowerat(p) !gRPT_Pvv(p) !Plin_DW(p) !      
        kernel = ((5._dp*k**12 - 231._dp*(p**2 - q**2)**6 + 6._dp*k**10*(p**2 + q**2) + &
         630._dp*k**2*(p**2 - q**2)**4*(p**2 + q**2) + 3._dp*k**8*(5._dp*p**4 + 6._dp*p**2*q**2 + 5._dp*q**4) - & 
         105._dp*k**4*(p**2 - q**2)**2*(5._dp*p**4 + 6._dp*p**2*q**2 + 5._dp*q**4) + & 
         20._dp*k**6*(5._dp*p**6 + 3._dp*p**4*q**2 + 3._dp*p**2*q**4 + 5._dp*q**6)))/ &
        (1024._dp*k**4*p**4*q**4)
        f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f4mu8v2sq_int

!******************************

   function Pk_f4mu2v2sq(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f4mu2v2sq,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f4mu2v2sq_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f4mu2v2sq = 1._dp*ans(1)/fac_norm

   end function Pk_f4mu2v2sq

   subroutine f4mu2v2sq_int(ndim,x,ncomp,f) !RSDTensorDecomposition.nb 
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        powq = matterpowerat(q) !gRPT_Pvv(q) !Plin_DW(q) !
        powp = matterpowerat(p) !gRPT_Pvv(p) !Plin_DW(p) !       
        kernel = (5._dp*(k**4 +(p**2-q**2)**2-2._dp*k**2*(p**2+q**2))**3)/(1024._dp*k**4*p**4*q**4)
     f(1)   = (2._dp*pi/k)*q*p*powq*powp*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f4mu2v2sq_int

!******************************

   function Pk_f3mu2v1v2b1(arg_k)
     implicit none
     real(kind=dp)   :: Pk_f3mu2v1v2b1,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absaccRSD/10._dp
     epsrel    = relaccRSD/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(f3mu2v1v2b1_int,method,epsrel,epsabs,ans,error,prob)
     Pk_f3mu2v1v2b1 = 1._dp*ans(1)/fac_norm

   end function Pk_f3mu2v1v2b1

   subroutine f3mu2v1v2b1_int(ndim,x,ncomp,f) !RSDTensorDecomposition.nb 
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,y,kernel,powq,powp
     real(kind=dp)     :: pdvq,pvvp,pvvq,pdvp      

     q    = kmaxRSD*x(1)+kminRSD*(1._dp-x(1))
     pmin = max(kminRSD,dabs(k-q))
     pmax = min(kmaxRSD,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        pdvq = matterpowerat(q) !gRPT_Pdv(q) !Plin_DW(q) !
        pvvp = matterpowerat(p) !gRPT_Pvv(p) !Plin_DW(p) !     
        pvvq = matterpowerat(q) !gRPT_Pvv(q) !Plin_DW(q) !
        pdvp = matterpowerat(p) ! gRPT_Pdv(p) !Plin_DW(p) !      
        kernel = (-3._dp*(k**4 + (p**2 - q**2)**2 - 2._dp*k**2*(p**2 + q**2))**2* &
        (q**2*Pdvq*Pvvp + p**2*Pdvp*Pvvq))/(128._dp*k**2*p**4*q**4)
        f(1)   = (2._dp*pi/k)*q*p*kernel*(kmaxRSD-kminRSD)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine f3mu2v1v2b1_int

!************************************************************
   subroutine warning(name,ier)
   implicit   none
   integer    :: ier
   character  :: name*(*)

     if (ier .eq. 1) then
       write(*,*) trim(name)//' achieved Maximum number of subdivisions'
     elseif (ier .eq. 2) then
       write(*,*) trim(name)//' Round off error detected'
     elseif (ier .eq. 3) then
       write(*,*) trim(name)//' Extremely bad integrand behaviour.'
     elseif (ier .eq. 4) then
       write(*,*) trim(name)//' algorithm does not converge. Accur. imposs'
     elseif (ier .eq. 5) then
       write(*,*) trim(name)//' Integrand is probably divergent.'
     elseif (ier .eq. 6) then
       write(*,*) trim(name)//' input is invalid'
     elseif (ier .eq. 7) then
       write(*,*) trim(name)//' bad integrand occurs in one cycle/abnormal termination'
     elseif (ier .eq. 8) then
       write(*,*) trim(name)//' maximum number of cycles allowed has been achieved'
     elseif (ier .eq. 9) then
       write(*,*) trim(name)//' the extrapol table for convergence does not converge'
     endif

   end subroutine warning


end module model_power
