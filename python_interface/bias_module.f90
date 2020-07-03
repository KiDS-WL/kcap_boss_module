module bias_module_config
    use twopt_model_def

    type wedges_config
        type(twopt_dataset)      :: dataset
        type(twopt_model)        :: model
    end type wedges_config

    public wedges_config

end module bias_module_config

module bias_module
    use, intrinsic :: ieee_arithmetic

    use iso_c_binding,  only : c_int, c_double, c_float, c_bool, c_loc, c_ptr, c_f_pointer
    use CosmologyTypes, only: CMBParams, TCosmoTheoryParams, CosmoSettings, npar_atz, max_derived_parameters
    use CosmoTheory,    only: TCosmoTheoryPredictions
    use NLRSD_settings, only: NLRSD_redshift, tot_NLRSD_redshift, num_NLRSD_redshift
    use twopt_model_def
    use FileUtils
    use CAMB,           only : nthermo_derived!derived_age, derived_zstar, derived_rstar, derived_thetastar, derived_DAstar, &
    ! derived_zdrag, derived_rdrag, derived_kD, derived_thetaD, derived_zEQ, derived_keq , &
    ! derived_thetaEQ, derived_theta_rs_EQ

    use bias_module_config

    implicit none

    !power spectrum and derived parameters (H and D_A)
    !Holds the NLRSD tables and therefore needs to be persistent between calls
    type(TCosmoTheoryPredictions)             :: Theory    
    type(CMBParams)                           :: CMB       !cosmological parameters

    contains
        function initialize_wedges(twopt_type, num_ell, num_points_use, num_bands_use, zm, om_fid, h0_fid, &
                                   window, n_window_x, n_window_y, &
                                   verbose) result(config_ptr) bind(c, name="initialize_wedges")            
            integer(kind=c_int), intent(in)   :: twopt_type, num_ell, num_points_use, num_bands_use
            real(kind=c_double), intent(in)   :: zm, om_fid, h0_fid

            integer(kind=c_int), intent(in) :: n_window_x, n_window_y
            real(kind=c_double), intent(in) :: window(n_window_x, n_window_y)
            
            integer(kind=c_int), intent(in) :: verbose

            type(c_ptr)                       :: config_ptr
            type(wedges_config), pointer      :: config
            
            allocate(config)

            if(verbose > 0) then
                feedback = 2
            else
                feedback = 0
            endif

            !Set dataset params
            config%dataset%twopt_type = twopt_type
            config%dataset%nsize%num_ell = num_ell
            config%dataset%nsize%num_points_use = num_points_use
            config%dataset%nsize%num_bands_use = num_bands_use
            config%dataset%zm = real(zm, kind=mcp)
            config%dataset%om_fid = real(om_fid, kind=mcp)
            config%dataset%h0_fid = real(h0_fid, kind=mcp)

            allocate(config%dataset%window, source=window)

            if(verbose > 0) write(*,*) "Calling init_model."
            call init_model(config%model, config%dataset)
            if(verbose > 0) write(*,*) "NLRSD_redshift", NLRSD_redshift, "tot_NLRSD_redshift", tot_NLRSD_redshift

            config_ptr = c_loc(config)
        end function initialize_wedges

        subroutine cleanup_wedges(config_ptr) bind(c, name="cleanup_wedges")
            type(c_ptr), value            :: config_ptr
            type(wedges_config), pointer  :: config

            call c_f_pointer(config_ptr, config)
            if(allocated(config%dataset%window)) deallocate(config%dataset%window)
            deallocate(config)
        end subroutine cleanup_wedges

        subroutine cleanup_cosmology() bind(c, name="cleanup_cosmology")
            if(allocated(Theory%MPK)) then
                call Theory%MPK%Clear()
                deallocate(Theory%MPK)
            end if
            if(allocated(Theory%growth_z)) then
                call Theory%growth_z%Clear()
                deallocate(Theory%growth_z)
            endif
            if(allocated(Theory%sigma8_z)) then
                call Theory%sigma8_z%Clear()
                deallocate(Theory%sigma8_z)
            end if

            NLRSD_redshift = 0.0
            tot_NLRSD_redshift = 0.0
            num_NLRSD_redshift = 0.0
        end subroutine cleanup_cosmology

        subroutine setup_cosmology(h, omdm, omb, omv, omk, omnuh2, nnu, w, wa, & !CMB params
                                   use_growth, local_lag_g2, local_lag_g3, &
                                   Pk_z, n_Pk_z, Pk_log_k, n_Pk_log_k, Pk, n_Pk_x, n_Pk_y, &
                                   growth_z, n_growth_z, sigma8_z, n_sigma8_z, &
                                   verbose) bind(c, name="setup_cosmology")
            real(kind=c_double), intent(in)  :: h, omdm, omb, omv, omk, omnuh2, nnu, w, wa
            
            logical(kind=c_bool), intent(in) :: use_growth, local_lag_g2, local_lag_g3

            integer(kind=c_int), intent(in)  :: n_Pk_z, n_Pk_log_k, n_Pk_x, n_Pk_y
            real(kind=c_double), intent(in)  :: Pk_z(n_Pk_z), Pk_log_k(n_Pk_log_k)
            real(kind=c_double), intent(in)  :: Pk(n_Pk_x, n_Pk_y)

            integer(kind=c_int), intent(in)  :: n_growth_z, n_sigma8_z
            real(kind=c_double), intent(in)  :: growth_z(n_growth_z), sigma8_z(n_sigma8_z)

            integer(kind=c_int), intent(in)  :: verbose

            !Set cosmological parameters
            CMB%h = real(h, kind=mcp)
            CMB%omdm = real(omdm, kind=mcp)
            CMB%omb = real(omb, kind=mcp)
            CMB%omv = real(omv, kind=mcp)
            CMB%omk = real(omk, kind=mcp)
            CMB%omnuh2 = real(omnuh2, kind=mcp)
            CMB%nnu = real(nnu, kind=mcp)
            CMB%w = real(w, kind=mcp)
            CMB%wa = real(wa, kind=mcp)

            Theory%init_NLRSD = .true.

            if(n_Pk_z /= n_Pk_y) then
                stop "Size of z array and second dimension of Pk do not match."
            else if(n_Pk_log_k /= n_Pk_x) then
                stop "Size of k array and first dimension of Pk do not match."
            else if(n_Pk_x /= size(Pk, dim=1)) then
                stop "n_Px_x and first dimension of Pk do not match."
            end if

            !set power spectrum (see CAMBCalc_SetDerived)
            ! call Theory%FreePK()
            if(.not. allocated(Theory%MPK)) allocate(Theory%MPK)
            call Theory%MPK%Init(Pk_log_k, Pk_z, Pk)

            if(verbose > 2) then
                write(*,*) "z: ", Pk_z
                write(*,*) "log_k: ", Pk_log_k(:4)
                write(*,*) "Pk: ", Pk(:,1)
                write(*,*) Pk_log_k(1), exp(Pk_log_k(1)), Pk_z(1)
                write(*,*) Theory%MPK%PowerAt(exp(Pk_log_k(1)), Pk_z(1)), log(Theory%MPK%PowerAt(exp(Pk_log_k(1)), Pk_z(1)))
            end if

            if(n_Pk_z /= n_growth_z) then
                stop "Size of z array and growth array do not match."
            else if(n_Pk_z /= n_sigma8_z) then
                stop "Size of z array and sigma8 array do not match."
            end if

            if(.not. allocated(Theory%growth_z)) allocate(Theory%growth_z)
            call Theory%growth_z%Clear()
            Theory%growth_z%n = n_Pk_z
            allocate(Theory%growth_z%x, source=Pk_z)
            allocate(Theory%growth_z%F, source=growth_z)

            if(.not. allocated(Theory%sigma8_z))  allocate(Theory%sigma8_z)
            call Theory%sigma8_z%Clear()
            Theory%sigma8_z%n = n_Pk_z
            allocate(Theory%sigma8_z%x, source=Pk_z)
            allocate(Theory%sigma8_z%F, source=sigma8_z)

            ! CosmoSettings (need to set z_outputs?)
            CosmoSettings%use_growth = use_growth
            CosmoSettings%local_lag_g2 = local_lag_g2
            CosmoSettings%local_lag_g3 = local_lag_g3
        end subroutine setup_cosmology

        subroutine compute_wedges(config_ptr, &
                                  b1, b2, gamma2, gamma3, a_vir, gamma, & ! Bias parameters
                                  z_index, &
                                  H_z, n_H_z, DA_z, n_DA_z, &
                                  bands, n_bands, &
                                  vtheo, n_vtheo, vtheo_convolved, n_vtheo_convolved, &
                                  Pk_mm, n_Pk_mm, Pk_gm, n_Pk_gm, Pk_gg, n_Pk_gg, &
                                  verbose) bind(c, name="compute_wedges")
            type(c_ptr), value               :: config_ptr
            
            real(kind=c_double), intent(in)  :: b1, b2, gamma2, gamma3, a_vir, gamma

            integer(kind=c_int), intent(in)  :: z_index

            integer(kind=c_int), intent(in)  :: n_H_z, n_DA_z
            real(kind=c_double), intent(in)  :: H_z(n_H_z), DA_z(n_DA_z)

            integer(kind=c_int), intent(in)  :: n_bands
            real(kind=c_double), intent(in)  :: bands(n_bands)

            integer(kind=c_int), intent(in)  :: n_vtheo, n_vtheo_convolved
            real(kind=c_double), intent(inout) :: vtheo(n_vtheo), vtheo_convolved(n_vtheo_convolved)

            integer(kind=c_int), intent(in)  :: n_Pk_mm, n_Pk_gm, n_Pk_gg
            real(kind=c_double), intent(inout) :: Pk_mm(n_Pk_mm), Pk_gm(n_Pk_gm), Pk_gg(n_Pk_gg)

            integer(kind=c_int), intent(in)  :: verbose

            !Variables
            type(wedges_config), pointer              :: config
            real(mcp), dimension(6)                   :: DataParams

            integer                                   :: i

            call c_f_pointer(config_ptr, config)

            DataParams(1) = b1
            DataParams(2) = b2
            DataParams(3) = gamma2
            DataParams(4) = gamma3
            DataParams(5) = a_vir
            DataParams(6) = gamma

            if(n_H_z /= n_DA_z) then
                stop "Sizes of H_z and DA_z do not match."
            end if
            ! Set everything to nan by default so in anything gets accessed that hasn't been set, it will hopefully break
            Theory%derived_parameters = ieee_value(Theory%derived_parameters, ieee_quiet_nan)

            if(z_index-1 > n_H_z) then
                stop "z_index larger than available H(z) redshifts"
            else if(z_index-1 > (max_derived_parameters-nthermo_derived)/npar_atz) then
                stop "z_index larger than Theroy%derived_parameters array."
            end if
            do i=1,min(n_H_z, (max_derived_parameters-nthermo_derived)/npar_atz)
                Theory%derived_parameters(nthermo_derived+(i-1)*npar_atz+2) = H_z(i)
                Theory%derived_parameters(nthermo_derived+(i-1)*npar_atz+3) = DA_z(i)
            end do
            if(verbose > 0) then
                write(*,*) "n_H_z", n_H_z
                write(*,*) "H(~z_index): ", H_z(max(1,z_index-4):min(n_H_z, z_index+4))
                write(*,*) "DA(~z_index): ", DA_z(max(1,z_index-4):min(n_DA_z, z_index+4))
                write(*,*) "z_index", z_index
                write(*,*) "H at z_index", Theory%derived_parameters(nthermo_derived+(z_index-2)*npar_atz+2)
                write(*,*) "DA at z_index", Theory%derived_parameters(nthermo_derived+(z_index-2)*npar_atz+3)
            end if

            if(n_vtheo /= config%dataset%nsize%num_bands_use * config%dataset%nsize%num_ell) then
                stop "Size of vtheo array doesn't match"
            end if
            
            if(verbose > 0) write(*,*) "Calling get_model."
            call get_model(config%model, z_index, bands, vtheo, CMB, Theory, DataParams)
            
            if(n_vtheo_convolved /= config%dataset%nsize%num_points_use * config%dataset%nsize%num_ell) then
                stop "Size of vtheo_convolved array doesn't match"
            end if
            call config%dataset%convolve(vtheo, vtheo_convolved)

            if(n_Pk_mm /= Theory%MPK%nx .or. n_Pk_gm /= n_Pk_mm .or. n_Pk_gg /= n_Pk_mm) then
                stop "Size of output power spectrum and input k array do not match."
            end if

            if(verbose > 0) write(*,*) "Getting power spectra."
            do i=1,n_Pk_mm
                Pk_mm(i) = gRPT_Pdd(exp(Theory%MPK%x(i)))
                Pk_gm(i) = gRPT_Pgd(exp(Theory%MPK%x(i)))
                Pk_gg(i) = gRPT_Pgg(exp(Theory%MPK%x(i)))
            end do

        end subroutine compute_wedges
end module bias_module