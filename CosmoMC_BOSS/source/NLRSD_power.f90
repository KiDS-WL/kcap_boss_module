module NLRSD_power
  use settings,        only : dp => mcp, feedback
  use CosmologyTypes,  only : nk_table, TNLRSD_tables
  use numrec
  use quadpack 
  use cubaint
  use NLRSD_settings,  only :  Matterpowerat, fac_norm, tns_full, & 
                              b1, b2, gam2, gam3minus, sigv_P, f, a_vir, A_fast, shot_noise, &
                              do_init_model , init_plin
  use nonlinear_power, only : gRPTpower, gRPTpower_22inv, init_za, get_root => root, &
                              I0, I2, sphbe0, sphbe1, sphbe2, sphbe3, &
                              sphbe4, sphbe5, sphbe6, sphbe7, sphbe8, sigv
  use bias,            only : Pk_b1b2, Pk_b1g2_prop, Pk_b1g2_mc, Pk_b2b2, Pk_b2g2, &
                              Pk_g2g2, Pk_b1g3m_prop, Pk_b2, Pk_g2_prop, Pk_g2_mc

  implicit none
  integer, private                             :: ell, nmu 
  real(kind=dp),  private, parameter           :: kmin = 0.0001_dp, kmax = 2.0_dp
  real(kind=dp),  private, parameter           :: relacc = 0.01_dp, absacc = 0.01_dp
  real(kind=dp),  private, parameter           :: kminRSD = 0.0001_dp, kmaxRSD = 10._dp
  real(kind=dp),  private, parameter           :: relaccRSD = 0.001_dp, absaccRSD = 0.01_dp
  real(kind=dp),  private, parameter           :: k1loop = 0.01_dp !kmin at which we start calculating 1loop
  real(kind=dp),  private                      :: k, q

  type(TNLRSD_tables), pointer                 :: NLRSD_tab
  

  contains

!************************************************************
! build tables
!************************************************************

    subroutine init_NLRSD_model()
      implicit none 
      real(kind=dp)                :: oldtime, tl, x
      call init_kvec
      call init_gRPT
      call init_bias_galgal
      if (f>0._dp) then  
        call init_bias_galvel
        call init_rsd_velvel
        call init_rsd_galvel
        call init_rsd_galgal
!         call init_VDskew
        if(feedback > 1)write(*,'(a)')'> Code Initialized'
      else
        if(feedback > 1)write(*,*)'running in real-space mode: quadrupole becomes Pgm'
      endif             
       do_init_model = .false.      
    end subroutine init_NLRSD_model

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
        NLRSD_tab%k_table(i) = (qmax - qmin)*dble(i-1)/dble(nk_table - 1) + qmin  
        NLRSD_tab%k_table(i) = sign(abs(NLRSD_tab%k_table(i))**pow,NLRSD_tab%k_table(i)) - kcenter 
      enddo
      NLRSD_tab%k_table = 10._dp**NLRSD_tab%k_table
      
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
        ak = NLRSD_tab%k_table(i)
        if (ak < k1loop) then
          !root             = get_root(ak)
          NLRSD_tab%dP13_dd_table(i) = 0._dp
          NLRSD_tab%dP13_dv_table(i) = 0._dp
          NLRSD_tab%dP13_vv_table(i) = 0._dp
          !NLRSD_tab%P22inv_dd_table(i) = 0._dp
          !NLRSD_tab%P22inv_dv_table(i) = 0._dp
          !NLRSD_tab%P22inv_vv_table(i) = 0._dp

          call gRPTpower_22inv(ak, root, P22inv_dd, P22inv_dv, P22inv_vv)

          NLRSD_tab%P22inv_dd_table(i) = P22inv_dd
          NLRSD_tab%P22inv_dv_table(i) = P22inv_dv
          NLRSD_tab%P22inv_vv_table(i) = P22inv_vv

          NLRSD_tab%gRPT_root_table(i) = root
        else
          call gRPTpower(ak, root, dP13_dd, dP13_dv, dP13_vv, P22inv_dd, P22inv_dv, P22inv_vv)

          NLRSD_tab%dP13_dd_table(i) = dP13_dd
          NLRSD_tab%dP13_dv_table(i) = dP13_dv
          NLRSD_tab%dP13_vv_table(i) = dP13_vv

          NLRSD_tab%P22inv_dd_table(i) = P22inv_dd
          NLRSD_tab%P22inv_dv_table(i) = P22inv_dv
          NLRSD_tab%P22inv_vv_table(i) = P22inv_vv

          NLRSD_tab%gRPT_root_table(i) = root
        end if
    
      enddo 
      if(feedback > 1)write(*,'(A)')'> Done with init gRPT'
      !close(93)

      call spline(NLRSD_tab%k_table, NLRSD_tab%dP13_dd_table, nk_table, 3d30,3d30, NLRSD_tab%dP13_dd_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%dP13_dv_table, nk_table, 3d30,3d30, NLRSD_tab%dP13_dv_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%dP13_vv_table, nk_table, 3d30,3d30, NLRSD_tab%dP13_vv_d_table)

      call spline(NLRSD_tab%k_table, NLRSD_tab%P22inv_dd_table, nk_table, 3d30,3d30, NLRSD_tab%P22inv_dd_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%P22inv_dv_table, nk_table, 3d30,3d30, NLRSD_tab%P22inv_dv_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%P22inv_vv_table, nk_table, 3d30,3d30, NLRSD_tab%P22inv_vv_d_table)

      call spline(NLRSD_tab%k_table, NLRSD_tab%gRPT_root_table,  nk_table, 3d30,3d30, NLRSD_tab%gRPT_root_d_table)

    end subroutine init_gRPT

    subroutine init_bias_galgal() !checked
      implicit   none
      integer         :: i
      real(kind=dp)   :: ak

      do i = 1, nk_table
        ak = NLRSD_tab%k_table(i)

        if (ak >= k1loop) then !tabulate 1loop integrals
          NLRSD_tab%b1b2_table(i)  = Pk_b1b2(ak)
          NLRSD_tab%b1g2_table(i)  = Pk_b1g2(ak)
          NLRSD_tab%b2b2_table(i)  = Pk_b2b2(ak)
          NLRSD_tab%b2g2_table(i)  = Pk_b2g2(ak)
          NLRSD_tab%g2g2_table(i)  = Pk_g2g2(ak) !was b1g2, bug!!!!
          NLRSD_tab%b1g3m_table(i) = Pk_b1g3m(ak)
        else
          NLRSD_tab%b1b2_table(i)  = 0._dp
          NLRSD_tab%b1g2_table(i)  = 0._dp
          NLRSD_tab%b2b2_table(i)  = 0._dp
          NLRSD_tab%b2g2_table(i)  = 0._dp
          NLRSD_tab%g2g2_table(i)  = 0._dp
          NLRSD_tab%b1g3m_table(i) = 0._dp
        endif        
      enddo 

      call spline(NLRSD_tab%k_table, NLRSD_tab%b1b2_table, nk_table, 3d30, 3d30, NLRSD_tab%b1b2_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%b1g2_table, nk_table, 3d30, 3d30, NLRSD_tab%b1g2_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%b2b2_table, nk_table, 3d30, 3d30, NLRSD_tab%b2b2_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%b2g2_table, nk_table, 3d30, 3d30, NLRSD_tab%b2g2_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%g2g2_table, nk_table, 3d30, 3d30, NLRSD_tab%g2g2_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%b1g3m_table, nk_table, 3d30, 3d30, NLRSD_tab%b1g3m_d_table)

    end subroutine init_bias_galgal
    
    subroutine init_bias_galvel() !checked
      implicit none
      integer       :: i
      real(kind=dp) :: ak

      do i = 1,nk_table
        ak = NLRSD_tab%k_table(i)
        if (ak >= k1loop) then !tabulate 1loop integrals
          NLRSD_tab%b2_table(i)   =  Pk_b2(ak)
          NLRSD_tab%g2_table(i)   =  Pk_g2(ak)
          NLRSD_tab%g3m_table(i)  =  Pk_g3m(ak)
        else
          NLRSD_tab%b2_table(i)  = 0._dp
          NLRSD_tab%g2_table(i)  = 0._dp
          NLRSD_tab%g3m_table(i) = 0._dp
        endif
      enddo 

      call spline(NLRSD_tab%k_table, NLRSD_tab%b2_table, nk_table, 3d30, 3d30, NLRSD_tab%b2_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%g2_table, nk_table, 3d30, 3d30, NLRSD_tab%g2_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%g3m_table, nk_table, 3d30, 3d30, NLRSD_tab%g3m_d_table)

    end subroutine init_bias_galvel

    subroutine init_rsd_velvel()
      implicit   none
      integer         :: i
      real(kind=dp)   :: ak

      do i = 1, nk_table
        ak = NLRSD_tab%k_table(i)
        if (ak >= k1loop) then !tabulate 1loop integrals
          NLRSD_tab%f3mu4v2sq_table(i) = Pk_f3mu4v2sq(ak)
          NLRSD_tab%f3mu6v2sq_table(i) = Pk_f3mu6v2sq(ak)
        
          NLRSD_tab%f4mu2v2sq_table(i) = Pk_f4mu2v2sq(ak)
          NLRSD_tab%f4mu4v2sq_table(i) = Pk_f4mu4v2sq(ak)
          NLRSD_tab%f4mu6v2sq_table(i) = Pk_f4mu6v2sq(ak)
          NLRSD_tab%f4mu8v2sq_table(i) = Pk_f4mu8v2sq(ak)
        else
          NLRSD_tab%f3mu4v2sq_table(i) = 0._dp
          NLRSD_tab%f3mu6v2sq_table(i) = 0._dp
                    
          NLRSD_tab%f4mu2v2sq_table(i) = 0._dp
          NLRSD_tab%f4mu4v2sq_table(i) = 0._dp
          NLRSD_tab%f4mu6v2sq_table(i) = 0._dp
          NLRSD_tab%f4mu8v2sq_table(i) = 0._dp
        endif
      enddo 

      call spline(NLRSD_tab%k_table, NLRSD_tab%f3mu4v2sq_table, nk_table, 3d30, 3d30, NLRSD_tab%f3mu4v2sq_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%f3mu6v2sq_table, nk_table, 3d30, 3d30, NLRSD_tab%f3mu6v2sq_d_table)
                                                                                                          
      call spline(NLRSD_tab%k_table, NLRSD_tab%f4mu2v2sq_table, nk_table, 3d30, 3d30, NLRSD_tab%f4mu2v2sq_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%f4mu4v2sq_table, nk_table, 3d30, 3d30, NLRSD_tab%f4mu4v2sq_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%f4mu6v2sq_table, nk_table, 3d30, 3d30, NLRSD_tab%f4mu6v2sq_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%f4mu8v2sq_table, nk_table, 3d30, 3d30, NLRSD_tab%f4mu8v2sq_d_table)
    end subroutine init_rsd_velvel

    subroutine init_rsd_galgal()
      implicit   none
      integer         :: i
      real(kind=dp)   :: ak

      do i = 1, nk_table
        ak = NLRSD_tab%k_table(i)

        if (ak >= k1loop) then !tabulate 1loop integrals
          NLRSD_tab%fmu2v1sqb1sq_table(i) = Pk_fmu2v1sqb1sq(ak)  
          NLRSD_tab%fmu2v1sqb1b2_table(i) = Pk_fmu2v1sqb1b2(ak)
          NLRSD_tab%fmu2v1sqb1g2_table(i) = Pk_fmu2v1sqb1g2(ak)

          NLRSD_tab%f2mu2v1sqb1sq_table(i) = Pk_f2mu2v1sqb1sq(ak)
          NLRSD_tab%f2mu4v1sqb1sq_table(i) = Pk_f2mu4v1sqb1sq(ak)
        else
          NLRSD_tab%fmu2v1sqb1sq_table(i) = 0._dp
          NLRSD_tab%fmu2v1sqb1b2_table(i) = 0._dp
          NLRSD_tab%fmu2v1sqb1g2_table(i) = 0._dp
                    
          NLRSD_tab%f2mu2v1sqb1sq_table(i) = 0._dp
          NLRSD_tab%f2mu4v1sqb1sq_table(i) = 0._dp
        endif
      enddo 

      call spline(NLRSD_tab%k_table, NLRSD_tab%fmu2v1sqb1sq_table, nk_table, 3d30, 3d30, NLRSD_tab%fmu2v1sqb1sq_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%fmu2v1sqb1b2_table, nk_table, 3d30, 3d30, NLRSD_tab%fmu2v1sqb1b2_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%fmu2v1sqb1g2_table, nk_table, 3d30, 3d30, NLRSD_tab%fmu2v1sqb1g2_d_table)
                                                                                                             
      call spline(NLRSD_tab%k_table, NLRSD_tab%f2mu2v1sqb1sq_table, nk_table, 3d30, 3d30, NLRSD_tab%f2mu2v1sqb1sq_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%f2mu4v1sqb1sq_table, nk_table, 3d30, 3d30, NLRSD_tab%f2mu4v1sqb1sq_d_table)
    end subroutine init_rsd_galgal

    subroutine init_rsd_galvel()
      implicit   none
      integer         :: i
      real(kind=dp)   :: ak

      do i = 1, nk_table
        ak = NLRSD_tab%k_table(i)

        if (ak >= k1loop) then !tabulate 1loop integrals
          NLRSD_tab%f2mu4v1v2b1_table(i)  = Pk_f2mu4v1v2b1(ak)
          NLRSD_tab%f2mu2v1v2b1_table(i)  = Pk_f2mu2v1v2b1(ak)
          NLRSD_tab%f2mu4v1v2b2_table(i)  = Pk_f2mu4v1v2b2(ak)
          NLRSD_tab%f2mu2v1v2b2_table(i)  = Pk_f2mu2v1v2b2(ak)
          NLRSD_tab%f2mu4v1v2g2_table(i)  = Pk_f2mu4v1v2g2(ak)
          NLRSD_tab%f2mu2v1v2g2_table(i)  = Pk_f2mu2v1v2g2(ak)
  
          NLRSD_tab%f3mu4v1v2b1_table(i)  = Pk_f3mu4v1v2b1(ak)
          NLRSD_tab%f3mu6v1v2b1_table(i)  = Pk_f3mu6v1v2b1(ak)

          NLRSD_tab%f3mu2v1v2b1_table(i)  = Pk_f3mu2v1v2b1(ak)
        else
          NLRSD_tab%f2mu4v1v2b1_table(i)  = 0._dp
          NLRSD_tab%f2mu2v1v2b1_table(i)  = 0._dp
          NLRSD_tab%f2mu4v1v2b2_table(i)  = 0._dp
          NLRSD_tab%f2mu2v1v2b2_table(i)  = 0._dp
          NLRSD_tab%f2mu4v1v2g2_table(i)  = 0._dp
          NLRSD_tab%f2mu2v1v2g2_table(i)  = 0._dp

          NLRSD_tab%f3mu4v1v2b1_table(i)  = 0._dp
          NLRSD_tab%f3mu6v1v2b1_table(i)  = 0._dp

          NLRSD_tab%f3mu2v1v2b1_table(i)  = 0._dp
        endif        
      enddo 

      call spline(NLRSD_tab%k_table, NLRSD_tab%f2mu4v1v2b1_table, nk_table, 3d30, 3d30, NLRSD_tab%f2mu4v1v2b1_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%f2mu2v1v2b1_table, nk_table, 3d30, 3d30, NLRSD_tab%f2mu2v1v2b1_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%f2mu4v1v2b2_table, nk_table, 3d30, 3d30, NLRSD_tab%f2mu4v1v2b2_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%f2mu2v1v2b2_table, nk_table, 3d30, 3d30, NLRSD_tab%f2mu2v1v2b2_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%f2mu4v1v2g2_table, nk_table, 3d30, 3d30, NLRSD_tab%f2mu4v1v2g2_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%f2mu2v1v2g2_table, nk_table, 3d30, 3d30, NLRSD_tab%f2mu2v1v2g2_d_table)
                                                                                                            
      call spline(NLRSD_tab%k_table, NLRSD_tab%f3mu4v1v2b1_table, nk_table, 3d30, 3d30, NLRSD_tab%f3mu4v1v2b1_d_table)
      call spline(NLRSD_tab%k_table, NLRSD_tab%f3mu6v1v2b1_table, nk_table, 3d30, 3d30, NLRSD_tab%f3mu6v1v2b1_d_table)
                                                                                                            
      call spline(NLRSD_tab%k_table, NLRSD_tab%f3mu2v1v2b1_table, nk_table, 3d30, 3d30, NLRSD_tab%f3mu2v1v2b1_d_table)
    end subroutine init_rsd_galvel
    
    
!************************************************************
!  Bias integrals that reqire wiggle exorcism (it turns out, not)
!************************************************************

   function Pk_b1g2(arg_k)
     implicit none
     real(kind=dp)   :: Pk_b1g2,arg_k

      Pk_b1g2   = Pk_b1g2_prop(arg_k)*matterpowerat(arg_k) + Pk_b1g2_mc(arg_k)

   end function Pk_b1g2
    
   function Pk_b1g3m(arg_k)
     implicit none
     real(kind=dp)   :: Pk_b1g3m,arg_k

     Pk_b1g3m   = Pk_b1g3m_prop(arg_k)*matterpowerat(arg_k)

   end function Pk_b1g3m

   function Pk_g2(arg_k)
   implicit none
   real(kind=dp)   :: Pk_g2,arg_k

    Pk_g2   = Pk_g2_prop(arg_k)*matterpowerat(arg_k) + Pk_g2_mc(arg_k)

   end function Pk_g2

   function Pk_g3m(arg_k)
     implicit none
     real(kind=dp)   :: Pk_g3m,arg_k

     Pk_g3m = Pk_b1g3m(arg_k)/2._dp

   end function Pk_g3m

!************************************************************
!   tabulated functions
!************************************************************

! Kaiser-like terms:
!************************************************************

! first separate the straight *matter* pdd,pdv,pvv

    function gRPT_Pdd(arg_k)
      implicit none
      real(kind=dp)   gRPT_Pdd,  arg_k, P0, root, dP13_dd,P22inv_dd

      if (arg_k >= kmin .and. arg_k <= kmax) then 

        if (arg_k > k1loop) then
         
          call splint(NLRSD_tab%k_table, NLRSD_tab%gRPT_root_table, NLRSD_tab%gRPT_root_d_table, nk_table, arg_k, root)
          call splint(NLRSD_tab%k_table, NLRSD_tab%dP13_dd_table,   NLRSD_tab%dP13_dd_d_table,   nk_table, arg_k, dP13_dd)
          call splint(NLRSD_tab%k_table, NLRSD_tab%P22inv_dd_table, NLRSD_tab%P22inv_dd_d_table, nk_table, arg_k, P22inv_dd)

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

    end function gRPT_Pdd

    function gRPT_Pdv(arg_k)
      implicit none
      real(kind=dp)   gRPT_Pdv,  arg_k, root, P0, dP13_dv,P22inv_dv

      if (arg_k >= kmin .and. arg_k <= kmax) then 

        if (arg_k > k1loop) then
         
          call splint(NLRSD_tab%k_table, NLRSD_tab%gRPT_root_table, NLRSD_tab%gRPT_root_d_table, nk_table, arg_k, root)
          call splint(NLRSD_tab%k_table, NLRSD_tab%dP13_dv_table,   NLRSD_tab%dP13_dv_d_table,   nk_table, arg_k, dP13_dv)
          call splint(NLRSD_tab%k_table, NLRSD_tab%P22inv_dv_table, NLRSD_tab%P22inv_dv_d_table, nk_table, arg_k, P22inv_dv)

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

    function gRPT_Pvv(arg_k)
    implicit none
    real(kind=dp)   gRPT_Pvv,  arg_k, root, P0, dP13_vv,P22inv_vv

      if (arg_k >= kmin .and. arg_k <=kmax) then 
        if (arg_k > k1loop) then
          call splint(NLRSD_tab%k_table, NLRSD_tab%gRPT_root_table, NLRSD_tab%gRPT_root_d_table, nk_table, arg_k, root)
          call splint(NLRSD_tab%k_table, NLRSD_tab%dP13_vv_table,   NLRSD_tab%dP13_vv_d_table,   nk_table, arg_k, dP13_vv)
          call splint(NLRSD_tab%k_table, NLRSD_tab%P22inv_vv_table, NLRSD_tab%P22inv_vv_d_table, nk_table, arg_k, P22inv_vv)

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

      call splint(NLRSD_tab%k_table, NLRSD_tab%b1b2_table, NLRSD_tab%b1b2_d_table, nk_table, arg_k, pow_b1b2)
      call splint(NLRSD_tab%k_table, NLRSD_tab%b1g2_table, NLRSD_tab%b1g2_d_table, nk_table, arg_k, pow_b1g2)
      call splint(NLRSD_tab%k_table, NLRSD_tab%b2b2_table, NLRSD_tab%b2b2_d_table, nk_table, arg_k, pow_b2b2)
      call splint(NLRSD_tab%k_table, NLRSD_tab%b2g2_table, NLRSD_tab%b2g2_d_table, nk_table, arg_k, pow_b2g2)
      call splint(NLRSD_tab%k_table, NLRSD_tab%g2g2_table, NLRSD_tab%g2g2_d_table, nk_table, arg_k, pow_g2g2)
      call splint(NLRSD_tab%k_table, NLRSD_tab%b1g3m_table, NLRSD_tab%b1g3m_d_table, nk_table, arg_k, pow_b1g3m)

      b1sq = b1**2 

      gRPT_Pgg = b1sq*pow_dd + A_fast**2*(b1*b2*pow_b1b2 + b1*gam2*pow_b1g2  &
      + b2*b2*pow_b2b2 + b2*gam2*pow_b2g2 + gam2*gam2*pow_g2g2 + b1*gam3minus*pow_b1g3m)

    end function gRPT_Pgg


    function gRPT_Pgv(arg_k)  
      implicit none
      real(kind=dp)   gRPT_Pgv,pow_dv,arg_k, b1run
      real(kind=dp)   pow_b2,pow_g2,pow_g3m

      pow_dv = gRPT_Pdv(arg_k)

      call splint(NLRSD_tab%k_table, NLRSD_tab%b2_table, NLRSD_tab%b2_d_table, nk_table, arg_k, pow_b2)
      call splint(NLRSD_tab%k_table, NLRSD_tab%g2_table, NLRSD_tab%g2_d_table, nk_table, arg_k, pow_g2)
      call splint(NLRSD_tab%k_table, NLRSD_tab%g3m_table, NLRSD_tab%g3m_d_table, nk_table, arg_k, pow_g3m)

      b1run    = b1 
      gRPT_Pgv = b1run*pow_dv + A_fast**2*(b2*pow_b2 + gam2*pow_g2 + gam3minus*pow_g3m)

    end function gRPT_Pgv


    function gRPT_Pgd(arg_k) 
      implicit none
      real(kind=dp)   gRPT_Pgd,pow_dd,arg_k, b1run
      real(kind=dp)   pow_b1b2,pow_b1g2,pow_b1g3m

      pow_dd = gRPT_Pdd(arg_k)

      call splint(NLRSD_tab%k_table, NLRSD_tab%b1b2_table,  NLRSD_tab%b1b2_d_table,  nk_table, arg_k, pow_b1b2)
      call splint(NLRSD_tab%k_table, NLRSD_tab%b1g2_table,  NLRSD_tab%b1g2_d_table,  nk_table, arg_k, pow_b1g2)
      call splint(NLRSD_tab%k_table, NLRSD_tab%b1g3m_table, NLRSD_tab%b1g3m_d_table, nk_table, arg_k, pow_b1g3m)

      b1run = b1 
      gRPT_Pgd = b1run*pow_dd + A_fast**2*0.5_dp*(b2*pow_b1b2 + gam2*pow_b1g2 + gam3minus*pow_b1g3m)
 
    end function


    function gRPT_roots(arg_k) !this we are not using, so no scaling is needed
    implicit none
    real(kind=dp)   gRPT_roots,root,arg_k

      call splint(NLRSD_tab%k_table,NLRSD_tab%gRPT_root_table,NLRSD_tab%gRPT_root_d_table,nk_table,arg_k,root)

      gRPT_roots = A_fast*root 

    end function
    
    function Plin_DW(arg_k) !desperate attempt
    implicit none
    real(kind=dp)   Plin_DW,root,arg_k,Ginv,x
      
    Plin_DW = gRPT_Pdv(arg_k)/A_fast
    
    end function Plin_DW
        
    
!    A_TNS or P12RS terms:
!************************************************************
!1.
    function RSDskew_fmu2v1sq(arg_k) !adding up 3 bias sq term that go as f mu^2 v1^2
      implicit none
      real(kind=dp)   RSDskew_fmu2v1sq, arg_k
      real(kind=dp)   pow_b1g2, pow_b1b2, pow_b1sq

      !! RSD squared bias terms !!
      call splint(NLRSD_tab%k_table, NLRSD_tab%fmu2v1sqb1sq_table, NLRSD_tab%fmu2v1sqb1sq_d_table, nk_table, arg_k, pow_b1sq)
      call splint(NLRSD_tab%k_table, NLRSD_tab%fmu2v1sqb1b2_table, NLRSD_tab%fmu2v1sqb1b2_d_table, nk_table, arg_k, pow_b1b2)
      call splint(NLRSD_tab%k_table, NLRSD_tab%fmu2v1sqb1g2_table, NLRSD_tab%fmu2v1sqb1g2_d_table, nk_table, arg_k, pow_b1g2)

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

      call splint(NLRSD_tab%k_table, NLRSD_tab%f2mu4v1v2b1_table, NLRSD_tab%f2mu4v1v2b1_d_table, nk_table, arg_k, pow_b1)
      call splint(NLRSD_tab%k_table, NLRSD_tab%f2mu4v1v2b2_table, NLRSD_tab%f2mu4v1v2b2_d_table, nk_table, arg_k, pow_b2)
      call splint(NLRSD_tab%k_table, NLRSD_tab%f2mu4v1v2g2_table, NLRSD_tab%f2mu4v1v2g2_d_table, nk_table, arg_k, pow_g2)

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

      call splint(NLRSD_tab%k_table, NLRSD_tab%f2mu2v1v2b1_table, NLRSD_tab%f2mu2v1v2b1_d_table, nk_table, arg_k, pow_b1)
      call splint(NLRSD_tab%k_table, NLRSD_tab%f2mu2v1v2b2_table, NLRSD_tab%f2mu2v1v2b2_d_table, nk_table, arg_k, pow_b2)
      call splint(NLRSD_tab%k_table, NLRSD_tab%f2mu2v1v2g2_table, NLRSD_tab%f2mu2v1v2g2_d_table, nk_table, arg_k, pow_g2)

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

      call splint(NLRSD_tab%k_table, NLRSD_tab%f3mu4v2sq_table, NLRSD_tab%f3mu4v2sq_d_table, nk_table, arg_k, pow_vv)

      RSDskew_f3mu4v2sq = A_fast**2*pow_vv
      
    end function

!5.
    function RSDskew_f3mu6v2sq(arg_k) ! f^3 mu^6 v2^2 vel-vel term
      implicit none
      real(kind=dp)   RSDskew_f3mu6v2sq, arg_k
      real(kind=dp)   pow_vv

      call splint(NLRSD_tab%k_table, NLRSD_tab%f3mu6v2sq_table, NLRSD_tab%f3mu6v2sq_d_table, nk_table, arg_k, pow_vv)

      RSDskew_f3mu6v2sq = A_fast**2*pow_vv
      
    end function


! B_TNS or P22RS terms:
!************************************************************
!1.
    function RSDkurt_f2mu2v1sq(arg_k) ! f^2 mu^2 v1^2 b1^2 term
      implicit none
      real(kind=dp)   RSDkurt_f2mu2v1sq, arg_k
      real(kind=dp)   pow_vv

      call splint(NLRSD_tab%k_table, NLRSD_tab%f2mu2v1sqb1sq_table, NLRSD_tab%f2mu2v1sqb1sq_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f2mu2v1sq = A_fast**2*b1**2*pow_vv
      
    end function

!2.
    function RSDkurt_f2mu4v1sq(arg_k) ! f^2 mu^4 v1^2 b1^2 term
      implicit none
      real(kind=dp)   RSDkurt_f2mu4v1sq, arg_k
      real(kind=dp)   pow_vv

      call splint(NLRSD_tab%k_table, NLRSD_tab%f2mu4v1sqb1sq_table, NLRSD_tab%f2mu4v1sqb1sq_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f2mu4v1sq = A_fast**2*(b1**2 *pow_vv)
      
    end function

!3.
    function RSDkurt_f3mu4v1v2(arg_k) ! f^3 mu^4 v1v2 b1 term
      implicit none
      real(kind=dp)   RSDkurt_f3mu4v1v2, arg_k
      real(kind=dp)   pow_vv

      call splint(NLRSD_tab%k_table, NLRSD_tab%f3mu4v1v2b1_table, NLRSD_tab%f3mu4v1v2b1_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f3mu4v1v2 = A_fast**2*b1*pow_vv
      
    end function

!4.
    function RSDkurt_f3mu6v1v2(arg_k) ! f^3 mu^6 v1v2 b1 term
      implicit none
      real(kind=dp)   RSDkurt_f3mu6v1v2, arg_k
      real(kind=dp)   pow_vv

      call splint(NLRSD_tab%k_table, NLRSD_tab%f3mu6v1v2b1_table, NLRSD_tab%f3mu6v1v2b1_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f3mu6v1v2 = A_fast**2*b1*pow_vv
      
    end function

!5.
    function RSDkurt_f4mu4v2sq(arg_k) ! f^4 mu^4 v2^2 term
      implicit none
      real(kind=dp)   RSDkurt_f4mu4v2sq, arg_k
      real(kind=dp)   pow_vv

      call splint(NLRSD_tab%k_table, NLRSD_tab%f4mu4v2sq_table, NLRSD_tab%f4mu4v2sq_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f4mu4v2sq = A_fast**2*pow_vv
      
    end function

!6.
    function RSDkurt_f4mu6v2sq(arg_k) ! f^4 mu^6 v2^2 term
      implicit none
      real(kind=dp)   RSDkurt_f4mu6v2sq, arg_k
      real(kind=dp)   pow_vv

      call splint(NLRSD_tab%k_table, NLRSD_tab%f4mu6v2sq_table, NLRSD_tab%f4mu6v2sq_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f4mu6v2sq = A_fast**2*pow_vv
      
    end function

!7.
    function RSDkurt_f4mu8v2sq(arg_k) ! f^4 mu^8 v2^2 term
      implicit none
      real(kind=dp)   RSDkurt_f4mu8v2sq, arg_k
      real(kind=dp)   pow_vv

      call splint(NLRSD_tab%k_table, NLRSD_tab%f4mu8v2sq_table, NLRSD_tab%f4mu8v2sq_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f4mu8v2sq = A_fast**2*pow_vv
      
    end function

!8.
    function RSDkurt_f4mu2v2sq(arg_k) ! f^4 mu^2 v2^2 term
      implicit none
      real(kind=dp)   RSDkurt_f4mu2v2sq, arg_k
      real(kind=dp)   pow_vv

      call splint(NLRSD_tab%k_table, NLRSD_tab%f4mu2v2sq_table, NLRSD_tab%f4mu2v2sq_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f4mu2v2sq = A_fast**2*pow_vv
      
    end function

!9.
    function RSDkurt_f3mu2v1v2(arg_k) ! f^3 mu^2 v1v2 b1 term
      implicit none
      real(kind=dp)   RSDkurt_f3mu2v1v2, arg_k
      real(kind=dp)   pow_vv

      call splint(NLRSD_tab%k_table, NLRSD_tab%f3mu2v1v2b1_table, NLRSD_tab%f3mu2v1v2b1_d_table, nk_table, arg_k, pow_vv)

      !! power contribution !!
      RSDkurt_f3mu2v1v2 = A_fast**2*b1*pow_vv
      
    end function


!************************************************************
! multipoles with flexible FOG factor that gets integrated separately
!************************************************************

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
      a = a_vir
      sigv_eff = dsqrt(A_fast*f**2*sigv**2 + sigv_P)  
      
      virial3_dmu = (mu**2)**nmu *LegendreP(mu) /dsqrt(1._dp+a**2*k**2*mu**2)& 
                *exp(-mu**2 *k**2 *sigv_eff**2/(1._dp+a**2*k**2*mu**2) ) !* sphbe0(b)

    end function virial3_dmu


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


!************************************************************
! power spectrum multipoles
!************************************************************

    function SIDpower_gal_monopole_at(arg_k)
      implicit none
      real(kind=dp)   :: SIDpower_gal_monopole_at,arg_k
      real(kind=dp)   :: w_dd, w_dv, w_vv, pkdd, pkdv, pkvv
      real(kind=dp)   :: wskew_fmu2v1sq, wskew_f2mu2v1v2, wskew_f2mu4v1v2, wskew_f3mu4v2sq, wskew_f3mu6v2sq
      real(kind=dp)   :: wkurt_f2mu4v1sq, wkurt_f3mu4v1v2, wkurt_f3mu6v1v2, wkurt_f2mu2v1sq
      real(kind=dp)   :: wkurt_f4mu4v2sq, wkurt_f4mu6v2sq, wkurt_f4mu8v2sq
      real(kind=dp)   :: wkurt_f4mu2v2sq, wkurt_f3mu2v1v2

      integer         :: i !~borrar
      real(kind=dp)   :: dlk,kkk !borrar

      if(do_init_model)then
        call init_NLRSD_model
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
        call init_NLRSD_model
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
        call init_NLRSD_model
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


    function gRPTpower_gal_monopole_at(arg_k)
      implicit none
      real(kind=dp)   :: gRPTpower_gal_monopole_at,arg_k
      if(do_init_model)then
        call init_NLRSD_model
      end if
      gRPTpower_gal_monopole_at = SIDpower_gal_monopole_at(arg_k) !+ &
                             !   0._dp* SDcorr_power_gal_monopole_at(arg_k)  
    end function gRPTpower_gal_monopole_at
  
    function gRPTpower_gal_quadrupole_at(arg_k)
      implicit none
      real(kind=dp)   :: gRPTpower_gal_quadrupole_at, arg_k
      if(do_init_model)then
        call init_NLRSD_model
      end if
      !if (f>0._dp) then  
      gRPTpower_gal_quadrupole_at = SIDpower_gal_quadrupole_at(arg_k)! + &
                              !    0._dp*  SDcorr_power_gal_quadrupole_at(arg_k)  
     ! else !return Pgd
     !   gRPTpower_gal_quadrupole_at = gRPT_Pgd(arg_k) 
     ! endif
    end function gRPTpower_gal_quadrupole_at
 
    function gRPTpower_gal_hexadecapole_at(arg_k)
      implicit none
      real(kind=dp)  :: gRPTpower_gal_hexadecapole_at, arg_k
      if(do_init_model)then
        call init_NLRSD_model
      end if
      gRPTpower_gal_hexadecapole_at = SIDpower_gal_hexadecapole_at(arg_k) !+ &
                              !      0._dp*  SDcorr_power_gal_hexadecapole_at(arg_k)
    end function gRPTpower_gal_hexadecapole_at

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


   function fmu2v1sqb1sq_dd(q)        
     implicit none
     real(kind=dp)   :: fmu2v1sqb1sq_dd,q,kernel

     kernel =  (-2._dp*k**2)/(3._dp*q**2)
     fmu2v1sqb1sq_dd = kernel*q**2*matterpowerat(q) !*gRPT_Pdv(q)

   end function fmu2v1sqb1sq_dd
   
   function fmu2v1sqb1sq_dv(q)       
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

   function f3mu4v2sq_dq(q)         
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
     

   end function Pk_f3mu6v2sq_prop

   function f3mu6v2sq_dq(q)         
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


   function f2mu2v1v2b1_dq(q)      
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


   function f2mu4v1v2b1_dv(q)     
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


end module NLRSD_power
