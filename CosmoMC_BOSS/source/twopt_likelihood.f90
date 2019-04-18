
!================================================================================
  module twopt_model_def
    use camb,           only: nthermo_derived !: number used in 'get_DV'
    use settings,       only: mcp, feedback  
    use cosmologytypes, only: CMBParams, TCosmoTheoryParams, CosmoSettings !: cosmological parameters: H0, w0, wa, Omegas
    use cosmotheory,    only: TCosmoTheoryPredictions!, & !: derived parameters (H, DA)
    use twopt_Data_def
    use fiducial,      only : compute_scaling  !: get D_V(z), D_A(z) and E(z)
    use NLRSD_settings
    use NLRSD_power
    use compute_twopt, only : get_model_pk_multi , get_model_pk_wedges, &
                              get_model_xi_multi, get_model_xi_wedges
    implicit none

    type :: twopt_model
      !type of measurement (multipoles or wedges in conf. or Fourier space)
      integer                :: twopt_type
      integer                :: num_ell  !number of wedges or multipoles
      real(mcp)              :: redshift !effective redshift 
      real(mcp)              :: dm_fid, hub_fid!to scale distances to the fiducial cosmology
      contains
        !methods
        procedure                  :: init_model           !init the model
        procedure                  :: get_model            !computes the model
        procedure, private         :: set_fiducial_scaling !compute the fiducial Dm(zeff) and H(zeff)
        procedure, private         :: fill_model_par       !assign values of internal model parameters
    end type twopt_model

    real(kind=mcp), parameter, private      :: epsabs = 1.e-2_mcp 
    real(kind=mcp), parameter, private      :: epsrel = 1.e-2_mcp 
    real(kind=mcp), private                 :: gamma_grav
    type(CMBParams), private                :: CMB2         !cosmological parameters

  contains !: method definitions
    !intialise the model
    subroutine init_model(model, dataset)
     implicit none
      class(twopt_model)                              :: model  !object reference
      type(twopt_dataset)                             :: dataset
      model%twopt_type      = dataset%twopt_type
      model%num_ell         = dataset%nsize%num_ell
      model%redshift        = dataset%zm
      !if(model%redshift > NLRSD_redshift) NLRSD_redshift = model%redshift
      !define redshift for non-linear model as the mean of all twopt datasets
      num_NLRSD_redshift = num_NLRSD_redshift + 1._mcp
      tot_NLRSD_redshift = tot_NLRSD_redshift + model%redshift
      NLRSD_redshift = tot_NLRSD_redshift/num_NLRSD_redshift
      z_nonlin = NLRSD_redshift 
      !Compute Dm and H for the fiducial cosmology
      call model%set_fiducial_scaling(dataset%om_fid, dataset%h0_fid) 
      if(feedback > 1) write(*,*) 'Two-point model initialised'
      !fix settings of the model. this needs more flexibility.
      call set_NLRSD_settings
    end subroutine init_model
    
    !: compute the fiducial D_v(z_eff)
    subroutine set_fiducial_scaling(model, om_fid, h0_fid)
     implicit none
      class(twopt_model) :: model    !object reference
      real(mcp)        :: om_fid, h0_fid   !fiducial values of omega_m and H_0/100
      real(mcp)        :: dm, hub  !distances to scale the measurements acording to the fiducial cosmology
      call compute_scaling(om_fid, h0_fid, model%redshift, dm, hub)
      model%dm_fid = dm 
      model%hub_fid  = hub
      if(feedback > 1)write(*,*)'fiducial distances at redshift ',model%redshift
      if(feedback > 1)write(*,*)'Dm (Mpc/h) =', model%dm_fid
      if(feedback > 1)write(*,*)'100E =', model%hub_fid*100._dp
    end subroutine set_fiducial_scaling

    !returns the model evaluated in k, given the parameters and theory
    subroutine get_model( model, z_index, bands, vtheo, CMB, Theory, DataParams)
      implicit none
      class(twopt_model)                         :: model     !object reference
      integer, intent(in)                        :: z_index   !index of the redshift for linear P(k)
      real(mcp), intent(in),dimension(:)         :: bands     !where to evaluate the model
      real(mcp), intent(out),dimension(:)        :: vtheo     !model
      class(CMBParams)                           :: CMB       !cosmological parameters
      class(TCosmoTheoryPredictions), intent(in) :: Theory  !power spectrum and derived parameters (H and D_A)
      real(mcp),dimension(:), intent(in)         :: DataParams  !model parameters in the same order in paramnames
      call fill_model_par( model, z_index, CMB, Theory, DataParams)
      select case(model%twopt_type)
        case(1) !P_ell
          call get_model_pk_multi(bands, model%num_ell, vtheo)
        case(2) !P_wed
          call get_model_pk_wedges(bands, model%num_ell, vtheo)
        case(3) !xi_ell
          call get_model_xi_multi(bands, model%num_ell, vtheo)
        case(4) !xi_wed
          call get_model_xi_wedges(bands, model%num_ell, vtheo)
      end select
    end subroutine get_model

    subroutine fill_model_par(model, z_index, CMB, Theory, DataParams)
      use camb
      implicit none
      class(twopt_model)                         :: model   !object reference
      integer, intent(in)                        :: z_index
      class(CMBParams),intent(in)                :: CMB     !cosmological parameters
      class(TCosmoTheoryPredictions), target     :: Theory  !power spectrum and derived parameters (H and D_A)
      real(mcp),dimension(:), intent(in)         :: DataParams  !model parameters order the same as in paramnames
      real(kind=mcp)                             :: dm , hub, ee, dm_rat, hub_rat, z_index_out
      real(mcp)                                  :: om_z, f_growth     ! Omega_m(z), f(z)=d ln(D)/d ln(a)
      if(.not. allocated(Theory%MPK))then
         write(*,*) 'ERROR: Your Theory%MPK derived type is not initialized. Make sure you are'
         write(*,*) '       calling a SetPk routine and filling your power spectra.'
         call MPIstop()
      end if
      if(Theory%init_NLRSD)then
        do_init_model = .true.
        Theory%init_NLRSD = .false.
      end if
      z_index_out = z_index - 1 !skip index for z=0
      PK_lin    => Theory%MPK
      NLRSD_tab => Theory%NLRSD
      call init_plin

      !specify redshift index of linear theory power spectrum arrays
      !scalings for fiducial cosmology
      dm   = Theory%derived_parameters(nthermo_derived+(z_index_out-1)*npar_atz+3)*(1._mcp + model%redshift)*CMB%h
      ee   = Theory%derived_parameters(nthermo_derived+(z_index_out-1)*npar_atz+2)/100._mcp!
      hub  = ee/CMB%h
      dm_rat   = dm/model%dm_fid
      hub_rat  = hub/model%hub_fid
      alpha_tr = dm_rat
      alpha_lo = 1._mcp/hub_rat
      alpha    = (alpha_tr**2*alpha_lo)**(1._mcp/3._mcp)
      !bias parameters:
      b1 = DataParams(1)
      b2 = DataParams(2)
      if(CosmoSettings%local_lag_g2)then
        gam2 = -2._mcp/7._mcp*(b1 - 1._mcp)
      else
        gam2 = DataParams(3)
      end if
      if(CosmoSettings%local_lag_g3)then
        gam3minus = 11._dp/63._dp*(b1 - 1._dp)*1.5_dp
      else
        gam3minus = DataParams(4)
      end if
      !RSD parameters:
      sigv_P = 0._mcp
      a_vir = DataParams(5)
      f  = Theory%growth_z%value(model%redshift)/Theory%sigma8_z%Value(model%redshift)
      A_fast = (Theory%sigma8_z%Value(model%redshift)/Theory%sigma8_z%Value(NLRSD_redshift))**2
      if(CosmoSettings%use_growth)then
         CMB2 = CMB
         !redefine f as omega_m(z)**gamma
         om_z = (CMB%omdm+CMB%omb)*(1._mcp + model%redshift)**3/hub**2 
         f = om_z**DataParams(6)
         write(*,*)'antes A_fast',A_fast,f
         A_fast = A_fast*(Dgrowth_gamma(model%redshift,DataParams(6))/Dgrowth_gamma(model%redshift,gamma_gr(CMB%w)))**2 
         write(*,*)'A_fast',A_fast
      end if
      shot_noise = DataParams(9)
    end subroutine fill_model_par

    function gamma_gr(w)
      real(kind=mcp)   :: gamma_gr, w
      if(w < -1._mcp)then
        gamma_gr = 0.55_mcp + 0.02*(1._mcp + w)
      else
        gamma_gr = 0.55_mcp + 0.05*(1._mcp + w)
      end if
    end function gamma_gr

    function Dgrowth_gamma(z, gamma_val)
      real(kind=mcp), intent(in)     :: z, gamma_val
      real(kind=mcp)                 :: Dgrowth_gamma, lnD, a
      real(kind=mcp)                 :: lna0 = -4._mcp
      real(kind=mcp)                 :: abserr
      integer                       :: neval, ierr
      gamma_grav = gamma_val
      a = 1._mcp/(1._mcp+z)
      call qags(ker_Dgrowth,exp(lna0),a,epsabs,epsrel,lnD,abserr,neval,ierr)
      lnD = lnD + lna0              !Border condition: D(a0)/a0 \equiv 1 => lnDa0 = lna0
      Dgrowth_gamma = exp(lnD)
    end function Dgrowth_gamma

    function ker_Dgrowth(a)
      implicit none
      real(kind=mcp), intent(in)     :: a
      real(kind=mcp)                 :: z, ker_Dgrowth, om_z
      z = (1._mcp/a) - 1._mcp
      om_z = (CMB2%omdm+CMB2%omb)*(1._mcp + z)**3/e2(CMB2,z)
      ker_Dgrowth = om_z**gamma_grav
      ker_Dgrowth = ker_Dgrowth/a
    end function ker_Dgrowth

    function e2(CMB,z)    !E(z)=H(z)/H_0
      implicit none
      class(CMBParams),intent(in)   :: CMB     !cosmological parameters
      real(kind=mcp), intent(in)    :: z
      real(kind=mcp)                :: e2, zplus, expon, a, fac, weff, arg, Om_r
      Om_r  = 2.4969d-5/CMB%h**2
      zplus  = 1._mcp + z
      a     = 1._mcp / zplus
      fac   = 1._mcp + (1._mcp - a)/log(a)
      weff  = CMB%w + CMB%wa*fac
      expon = 3._mcp*(1._mcp+weff)
      arg = 187000._mcp*a*CMB%omnuh2
      e2 = (CMB%omb+CMB%omdm)/a**3 + CMB%omv/a**expon + CMB%omk/a**2 + Om_r/a**4*(1._mcp+0.2271_mcp*CMB%nnu*func(arg))
    end function e2

    function func(y)
      !equation (26) of Komatsu et al. (2009)
      implicit none
      real(kind=mcp)             :: func, y
      real(kind=mcp), parameter  :: a = 0.3173_mcp, p = 1.83_mcp
      real(kind=mcp), parameter  :: pinv = 1._mcp/p
      func = (1._mcp + (a*y)**p)**pinv
    end function func

end module twopt_model_def

!================================================================================
!The type LSSlikelihood is filled and added  to LikeList.
!LSSlikelihood contains a type with the data, one with the model and
!defines the appropriate likelihood
module twopt_likelihood_def
  use settings!, only: mcp, feedback, TSettingIni   !real precision and feedback level
  use CosmologyTypes!, only: CMBParams !: cosmological parameters: H0, w0, wa, Omegas
  use Likelihood_Cosmology!, only: TCosmoCalcLikelihood
  use CosmoTheory!,    only: TCosmoTheoryPredictions!, & !: derived parameters (H, DA)
  use likelihood
  use twopt_Data_def, only: twopt_dataset     !object that read and stores the data
  use twopt_model_def, only: twopt_model      !object that define the model
  use compute_twopt, only: test_found_nan     !object that define the model
  implicit none

  type, extends(TCosmoCalcLikelihood) :: TtwoptLikelihood
    type(twopt_dataset), pointer, private  :: dataset    ! read and store measurement
    type(twopt_model), pointer, private    :: model      ! implent the model
    real(mcp), allocatable, dimension(:)   :: twopt_th   ! model LSS measurement
    real(mcp), allocatable, dimension(:)   :: twopt_wth  ! convolved LSS measurement
    real(mcp), allocatable, dimension(:)   :: covdat, covth 
    character(LEN=20)                      :: twopt_version =  'November_2018'
    contains
      !likelihood
      procedure             :: LogLike => twopt_LnLike  
      procedure             :: ReadIni => init_likelihood
  end type TtwoptLikelihood

  contains

    !twopt likelihood
    function twopt_LnLike(this, CMB, Theory, DataParams)
     implicit none
      class(TtwoptLikelihood)         :: this           ! object reference 
      class(CMBParams)                :: CMB            ! cosmological parameters
      class(TCosmoTheoryPredictions), target :: Theory ! power spectrum and derived parameters
      real(mcp)                       :: DataParams(:)  ! model parameters
      real(mcp)                       :: twopt_LNLike
      real(mcp)                       :: chisq          ! chi square
      integer                         :: i,nd,iw,nw              ! delete this
      !borrame!
      real(kind=mcp)                :: oldtime, tl, x

      if(abs(this%exact_z(1)-Theory%MPK%y(this%exact_z_index(1)))>1.d-3)then
        write(*,*)'ERROR: two-point redshift does not match the value stored'
        write(*,*)'in the PK%y array.'
        write(*,*)'for dataset ', this%dataset%name
        call MpiStop()
      end if

      call cpu_time(oldtime)
      call this%model%get_model(this%exact_z_index(1), this%dataset%bands, &
                                this%twopt_th, CMB, Theory, DataParams)
      !convolve with the window function
      call this%dataset%convolve(this%twopt_th, this%twopt_wth)

      !compute chi2
      this%covdat = this%dataset%measur - this%twopt_wth       !: data minus theory
      this%covth  = matmul(this%dataset%invcov, this%covdat)   !: C^-1 * (d-t)
      chisq       = sum(this%covdat * this%covth)  
      !log(likelihood)
      twopt_LnLike = 0.5_mcp*chisq
      if(chisq /= chisq)then
        write(*,*)'chi2 is NaN!!'
        call test_found_nan
        chisq = 1.e20
      end if
         
      if(feedback>1) write(*,*) 'two-point clustering '//this%name//' chi-sq:',chisq

      call cpu_time(x)
      tl = x-oldtime
      if(feedback > 1)write(*,*)'time=',tl
      
      if (twopt_LnLike > 1e8) then
        write(*,*) 'Chisq is huge, maybe there is a problem? chisq=',chisq
      end if
    end function twopt_LnLike

    subroutine init_likelihood(this, Ini)
      use IniObjects
      implicit none
      class (Ttwoptlikelihood)                  :: this  !object reference
      class(TSettingIni)                        :: Ini
      integer                                   :: nl, nb, np

      ! some common settings
      this%LikelihoodType             = 'twopt' 
      this%version                    = this%twopt_version     
      this%needs_background_functions = .true.  
      this%needs_powerspectra         = .true.  
      this%needs_exact_z              = .true.
      this%needs_nonlinear_pk         = .false.
      this%num_z                      = 1
      this%speed                      = 2     

      ! read data and model
      allocate(this%dataset)
      call this%dataset%read_data(Ini)
      this%name = this%dataset%name

      allocate(this%exact_z(this%num_z))
      this%exact_z(1) = this%dataset%zm

      !: set up the model
      allocate(this%model)
      call this%model%init_model(this%dataset)

      !allocate only once arrays used to store the model, convolved model
      !and other quantities needed when computing the likelihood
      nl = this%dataset%nsize%num_ell
      np = this%dataset%nsize%num_points_use
      nb = this%dataset%nsize%num_bands_use
      allocate(this%twopt_th(nl*nb))   !theory
      allocate(this%twopt_wth(nl*np))  !convolved theory
      allocate(this%covdat(nl*np))   !matmul storage
      allocate(this%covth(nl*np))    !matmul storage

      if(feedback>1) write(*,*) 'Likelihood ', trim(this%name), " initialised"

    end subroutine init_likelihood

end module twopt_likelihood_def

!: This module contains the subroutine that adds the likelihood to cosmomc.
module twoptlikelihood
  use settings,             only: mcp, feedback !: real precision and feedback level
  use StringUtils,          only: numcat !: unify a string and a number to get a file name
  use likelihood,           only: TLikelihoodList
  use twopt_likelihood_def, only: Ttwoptlikelihood
  use IniObjects
  use CosmologyTypes
  use CosmoTheory
  implicit none

  logical                              :: use_twopt         !use LSS information
  integer                              :: twopt_numdatasets !number of datasets 
  contains

    !: call this function in DataLikelihoods.f90
    subroutine twopt_Likelihood_Add(LikeList, Ini)
      use IniObjects
      implicit none
      class(TLikelihoodList)                         :: LikeList
      !type(TIniFile)                                 :: Ini
      class(TSettingIni)                             :: ini
      type(Ttwoptlikelihood), pointer, dimension(:)  :: like ! hope that works also without pointer
      integer                                        :: i              !loop integer
      character(LEN=180)                             :: twopt_filename   !dataset file
      character(LEN=180)                             :: twopt_paramnames !parameter names file
      !indices to order 'like' according to redshifts. Feeded back to
      !likelihoods
      integer, dimension(:), allocatable          :: redshifts_indices 
      !store the size of 'z_outputs' before adding the redshifts in this subroutine
      integer                                     :: size_z_outputs 

      use_twopt = (Ini%Read_Logical('use_twopt',.false.))
      !if twopt is not required, do nothing and return
      if (.not. use_twopt) then
        if(feedback > 1) write(*,*) 'two-point clustering measurements not used'
        return
      end if

      !get the number of datasets
      twopt_numdatasets = Ini%Read_Int('twopt_numdatasets',0)
      if(feedback > 0) write(*,*) 'setting up ', twopt_numdatasets, " likelihoods"
      if(twopt_numdatasets < 1) then
        write(*,*) 'At least one dataset must be given'
        stop
      end if

      !allocate likelihoods
      allocate(like(twopt_numdatasets))

      !read the datasets
      do i=1, twopt_numdatasets
        !allocate and initilise the likelihood
        call like(i)%ReadDatasetFile(Ini%ReadFileName(numcat('twopt_dataset',i)) )
        !add the model parameters 
        twopt_paramnames = Ini%ReadFileName(numcat('twopt_paramnames',i))
        call like(i)%loadParamNames(twopt_paramnames) 

        call LikeList%add(like(i))
      end do

      if(feedback > 0) write(*,*) 'Two-point clustering likelihoods added to "LikeList"'
    end subroutine twopt_Likelihood_Add

end module twoptlikelihood

