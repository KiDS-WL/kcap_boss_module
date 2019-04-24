module bias_module
    use, intrinsic :: ieee_arithmetic

    use iso_c_binding,  only : c_int, c_double, c_float, c_bool
    use CosmologyTypes, only: CMBParams, TCosmoTheoryParams, CosmoSettings, npar_atz, max_derived_parameters
    use CosmoTheory,    only: TCosmoTheoryPredictions
    use twopt_model_def
    use FileUtils
    use CAMB,           only : nthermo_derived!derived_age, derived_zstar, derived_rstar, derived_thetastar, derived_DAstar, &
    ! derived_zdrag, derived_rdrag, derived_kD, derived_thetaD, derived_zEQ, derived_keq , &
    ! derived_thetaEQ, derived_theta_rs_EQ
    implicit none

    contains
        subroutine compute_wedges(h, omdm, omb, omv, omk, omnuh2, nnu, w, wa, & !CMB params
                                  b1, b2, gamma2, gamma3, a_vir, gamma, & ! Bias parameters
                                  H_z, n_H_z, DA_z, n_DA_z, &
                                  twopt_type, num_ell, num_points_use, num_bands_use, z_index, zm, om_fid, h0_fid, & !dataset params
                                  use_growth, local_lag_g2, local_lag_g3, &
                                  window, n_window_x, n_window_y, &
                                  bands, n_bands, &
                                  Pk_z, n_Pk_z, Pk_log_k, n_Pk_log_k, Pk, n_Pk_x, n_Pk_y, &
                                  growth_z, n_growth_z, sigma8_z, n_sigma8_z, &
                                  vtheo, n_vtheo, vtheo_convolved, n_vtheo_convolved, &
                                  Pk_mm, n_Pk_mm, Pk_gm, n_Pk_gm, Pk_gg, n_Pk_gg, &
                                  verbose) bind(c, name="compute_wedges")
            real(kind=c_double), intent(in) :: h, omdm, omb, omv, omk, omnuh2, nnu, w, wa
            real(kind=c_double), intent(in) :: b1, b2, gamma2, gamma3, a_vir, gamma

            integer(kind=c_int), intent(in) :: n_H_z, n_DA_z
            real(kind=c_double), intent(in) :: H_z(n_H_z), DA_z(n_DA_z)

            integer(kind=c_int), intent(in) :: twopt_type, num_ell, num_points_use, num_bands_use, z_index
            real(kind=c_double), intent(in) :: zm, om_fid, h0_fid
            logical(kind=c_bool), intent(in) :: use_growth, local_lag_g2, local_lag_g3

            integer(kind=c_int), intent(in) :: n_window_x, n_window_y, n_bands
            real(kind=c_double), intent(in) :: window(n_window_x, n_window_y), bands(n_bands)

            integer(kind=c_int), intent(in) :: n_Pk_z, n_Pk_log_k, n_Pk_x, n_Pk_y
            real(kind=c_double), intent(in) :: Pk_z(n_Pk_z), Pk_log_k(n_Pk_log_k)
            real(kind=c_double), intent(in) :: Pk(n_Pk_x, n_Pk_y)

            integer(kind=c_int), intent(in) :: n_growth_z, n_sigma8_z
            real(kind=c_double), intent(in) :: growth_z(n_growth_z), sigma8_z(n_sigma8_z)

            integer(kind=c_int), intent(in) :: n_vtheo, n_vtheo_convolved
            real(kind=c_double), intent(inout) :: vtheo(n_vtheo), vtheo_convolved(n_vtheo_convolved)

            integer(kind=c_int), intent(in) :: n_Pk_mm, n_Pk_gm, n_Pk_gg
            real(kind=c_double), intent(inout) :: Pk_mm(n_Pk_mm), Pk_gm(n_Pk_gm), Pk_gg(n_Pk_gg)

            integer(kind=c_int), intent(in) :: verbose

            !Variables
            type(twopt_model)                         :: model     !object reference
            type(twopt_dataset)                       :: dataset

            type(CMBParams)                           :: CMB       !cosmological parameters
            type(TCosmoTheoryPredictions)             :: Theory    !power spectrum and derived parameters (H and D_A)
            real(mcp), dimension(:), allocatable      :: DataParams

            integer                                   :: i

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

            if(n_H_z /= n_DA_z) then
                stop "Sizes of H_z and DA_z do not match."
            end if
            ! Set everything to nan by default so in anything gets accessed that hasn't been set, it will hopefully break
            Theory%derived_parameters = ieee_value(Theory%derived_parameters, ieee_quiet_nan)

            if(z_index > n_H_z) then
                stop "z_index larger than available H(z) redshifts"
            else if(z_index > (max_derived_parameters-nthermo_derived)/npar_atz) then
                stop "z_index larger than Theroy%derived_parameters array."
            end if
            do i=1,min(n_H_z, (max_derived_parameters-nthermo_derived)/npar_atz)
                Theory%derived_parameters(nthermo_derived+(i-1)*npar_atz+2) = H_z(i)
                Theory%derived_parameters(nthermo_derived+(i-1)*npar_atz+3) = DA_z(i)
            end do
            if(verbose > 0) then
                write(*,*) "H(~z_index): ", H_z(z_index-4:z_index+4)
                write(*,*) "DA(~z_index): ", DA_z(z_index-4:z_index+4)
                write(*,*) "z_index", z_index
                write(*,*) "H at z_index", Theory%derived_parameters(nthermo_derived+(z_index-2)*npar_atz+2)
                write(*,*) "DA at z_index", Theory%derived_parameters(nthermo_derived+(z_index-2)*npar_atz+3)
            end if

            !Set dataset params
            dataset%twopt_type = twopt_type
            dataset%nsize%num_ell = num_ell
            dataset%nsize%num_points_use = num_points_use
            dataset%nsize%num_bands_use = num_bands_use
            dataset%zm = real(zm, kind=mcp)
            dataset%om_fid = real(om_fid, kind=mcp)
            dataset%h0_fid = real(h0_fid, kind=mcp)
    
            allocate(dataset%window, source=window)
            
            if(n_Pk_z /= n_Pk_y) then
                stop "Size of z array and second dimension of Pk do not match."
            else if(n_Pk_log_k /= n_Pk_x) then
                stop "Size of k array and first dimension of Pk do not match."
            else if(n_Pk_x /= size(Pk, dim=1)) then
                stop "n_Px_x and first dimension of Pk do not match."
            end if
            ! allocate(z(n_Pk_z), source=Pk_z)
            ! allocate(log_k(n_Pk_log_k), source=Pk_log_k)

            !set power spectrum (see CAMBCalc_SetDerived)
            ! call Theory%FreePK()
            allocate(Theory%MPK)
            call Theory%MPK%Init(Pk_log_k, Pk_z, Pk)

            if(verbose > 0) then
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

            allocate(Theory%growth_z)
            call Theory%growth_z%Clear()
            Theory%growth_z%n = n_Pk_z
            allocate(Theory%growth_z%x, source=Pk_z)
            allocate(Theory%growth_z%F, source=growth_z)

            allocate(Theory%sigma8_z)
            call Theory%sigma8_z%Clear()
            Theory%sigma8_z%n = n_Pk_z
            allocate(Theory%sigma8_z%x, source=Pk_z)
            allocate(Theory%sigma8_z%F, source=sigma8_z)

            allocate(DataParams(6))
            DataParams(1) = b1
            DataParams(2) = b2
            DataParams(3) = gamma2
            DataParams(4) = gamma3
            DataParams(5) = a_vir
            DataParams(6) = gamma

            ! CosmoSettings (need to set z_outputs?)
            CosmoSettings%use_growth = use_growth
            CosmoSettings%local_lag_g2 = local_lag_g2
            CosmoSettings%local_lag_g3 = local_lag_g3

            feedback = 2
            write(*,*) "Calling init_model."
            call init_model(model, dataset)
            write(*,*) "Calling get_model."

            if(n_vtheo /= num_bands_use*num_ell) then
                stop "Size of vtheo array doesn't match"
            end if
            call get_model(model, z_index, bands, vtheo, CMB, Theory, DataParams)
            
            if(n_vtheo_convolved /= num_points_use*num_ell) then
                stop "Size of vtheo_convolved array doesn't match"
            end if
            call dataset%convolve(vtheo, vtheo_convolved)

            if(n_Pk_mm /= n_Pk_log_k .or. n_Pk_gm /= n_Pk_log_k .or. n_Pk_gg /= n_Pk_log_k) then
                stop "Size of output power spectrum and input k array do not match."
            end if

            write(*,*) "Getting power spectra."
            do i=1,n_Pk_log_k
                Pk_mm(i) = gRPT_Pdd(exp(Pk_log_k(i)))
                Pk_gm(i) = gRPT_Pgd(exp(Pk_log_k(i)))
                Pk_gg(i) = gRPT_Pgg(exp(Pk_log_k(i)))
            end do

        end subroutine compute_wedges
end module bias_module