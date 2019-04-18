!pqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpqpq
!    settings for mcmc_engine   - scaling                     !
!bdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbd

module NLRSD_settings
  use settings, only : dp => mcp
  use CosmoTheory   !for the linear theory power spectrum
  use numrec
  implicit none
  !a few constants...
  real(kind=dp),parameter             :: pi       = 4._dp*atan(1._dp)
  real(kind=dp),parameter             :: twopi    = 2._dp*pi
  real(kind=dp),parameter             :: pisqr    = pi*pi
  real(kind=dp),parameter             :: twopisqr = 2._dp*pi*pi
  !Set different options of the model of non-linearities and RSD
  integer                             :: tns_full
  real(kind=dp)                       :: fac_norm 
  !redshift for the evaluation of the linear matter power spectrum
  !the initial value is set to -1 to check that it is initialised correctly
  real(kind=dp)                               :: NLRSD_redshift = 0._dp
  real(kind=dp)                               :: num_NLRSD_redshift = 0._dp
  real(kind=dp)                               :: tot_NLRSD_redshift = 0._dp
  !model parameters                    
  real(kind=dp)                               :: alpha_tr, alpha_lo, alpha
  real(kind=dp)                               :: b1, b2, gam2, gam3minus 
  real(kind=dp)                               :: f, sigv_P, a_vir
  real(kind=dp)                               :: A_fast, shot_noise
  !recompute tables of the NL model
  logical                                     :: do_init_model = .true.
  !pointer to linear theory power spectrum
  Type(TCosmoTheoryPK), pointer               :: PK_lin
  integer, parameter, private                 :: nlin = 800
  real(kind=dp), parameter, private           :: akmin = 0.0001_dp, akmax = 1.8_dp
  real(kind=dp), dimension(nlin), private     :: klin, plin, plind
  !integer,  private                 :: nlin = 800
  !real(kind=dp), private           :: akmin = 0.0001_dp, akmax = 10._dp
  !real(kind=dp), dimension(:), allocatable    :: klin, plin, plind
  real(kind=dp), private                      :: n_eff_1, n_eff_2


  contains

    subroutine set_NLRSD_settings()
      implicit none
      integer         :: i
      real(kind=dp)   :: dklin

      !use original TNS or eTNS
      tns_full =  1
      fac_norm = (2._dp*pi)**3
      !use true amplitude of model power spectrum
      A_fast = 1._dp

      !define k values where linear theory P(k) will be computed
      dklin = (log10(akmax)-log10(akmin))/real(nlin-1,kind=dp)
      do i = 1,nlin
        klin(i) = log10(akmin) + dklin*real(i-1,kind=dp)
      end do
      klin = 10._dp**klin
    end subroutine set_NLRSD_settings

    !wraper for the linear theory power spectrum
    function matterpowerat2(kval)
      real(kind=dp)                 :: matterpowerat2, kval
      matterpowerat2 = PK_lin%PowerAt(kval, NLRSD_redshift)
    end function matterpowerat2

    subroutine init_plin()
      implicit   none
      integer          :: i, iwr
      real(kind=dp)    :: dlp, dlk
      integer  :: check
    
      do i=1,nlin
        plin(i) = PK_lin%PowerAt(klin(i), NLRSD_redshift)
      enddo
      !compute neff at high_k
      dlp = log10(plin(3)) - log10(plin(1))
      dlk = log10(klin(3)) - log10(klin(1))
      n_eff_1 = dlp/dlk
      dlp = log10(plin(nlin)) - log10(plin(nlin-3))
      dlk = log10(klin(nlin)) - log10(klin(nlin-3))
      n_eff_2 = dlp/dlk
      call spline(klin, plin, nlin, 3d30, 3d30, plind)
    
    end subroutine init_plin

    function Matterpowerat(QQ)   !linear power
      implicit  none
      real(kind=dp)    :: Matterpowerat, QQ, lQQ
      if (QQ >= akmin .and. QQ <= akmax) then 
        call splint(klin, plin, plind, nlin, QQ, Matterpowerat)
        !rescale power spectrum amplitude
      elseif(QQ < akmin)then
        Matterpowerat = plin(1)*(QQ/klin(1))**n_eff_1
      elseif(QQ > akmax)then
        Matterpowerat = plin(nlin)*(QQ/klin(nlin))**n_eff_2
      endif
    end function Matterpowerat

    subroutine write_plin()
      implicit   none
      integer          :: i
      open(99,file = 'matterpower_lin_test.dat')
      do i=1,nlin
        write(99,*)klin(i),plin(i)
      enddo
      close(99) 
      open(99,file = 'nuisance_params_test.dat')
      write(99,*)'A_fast',A_fast
      write(99,*)'alpha_tr',alpha_tr
      write(99,*)'alpha_lo',alpha_lo
      write(99,*)'f',f
      write(99,*)'b1',b1
      write(99,*)'b2',b2
      write(99,*)'gam2',gam2
      write(99,*)'gam3minus',gam3minus
      write(99,*)'a_vir',a_vir
      close(99)
    end subroutine write_plin

end module NLRSD_settings


