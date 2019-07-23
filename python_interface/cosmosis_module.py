import os
import collections.abc

import numpy as np
import scipy.interpolate

from cosmosis.datablock import option_section, names

from cosmosis.datablock.cosmosis_py import errors

import bias_module

def setup(options):
    config = {}
    
    config["module"] = bias_module.BiasModule()
    
    here = os.path.dirname(__file__)
    window_file = options[option_section, "window_file"]
    bands_file = options[option_section, "bands_file"]
    config["window"] = np.loadtxt(window_file)
    config["bands"] = np.loadtxt(bands_file)

    config["verbose"] = options.get_bool(option_section, "verbose", False)

    try:
        bands_range = options[option_section, "bands_range"]
    except errors.BlockNameNotFound:
        bands_range = [20, 160]

    try:
        points_range = options[option_section, "points_range"]
    except errors.BlockNameNotFound:
        points_range = [4, 32]
    
    n_points, n_bands = config["window"].shape
    
    cut_points_idx = list(range(points_range[0])) + list(range(points_range[1], n_points))
    cut_bands_idx = list(range(bands_range[0])) + list(range(bands_range[1], n_bands))

    if config["verbose"]:
        print("Cutting bands", cut_bands_idx)
        print("Cutting points", cut_points_idx)
    if len(cut_points_idx) > 0:
        config["window"] = np.delete(config["window"], cut_points_idx, axis=0)
    if len(cut_bands_idx) > 0:
        config["bands"] = np.delete(config["bands"], cut_bands_idx)
        config["window"] = np.delete(config["window"], cut_bands_idx, axis=1)

    config["output_section_wedges"] = options.get_string(option_section, "output_section_wedges", "xi_wedges")
    config["output_section_pk_mm"] = options.get_string(option_section, "output_section_pk_mm", "matter_matter_power_spectrum_pt")
    config["output_section_pk_gm"] = options.get_string(option_section, "output_section_pk_gm", "galaxy_matter_power_spectrum_pt")
    config["output_section_pk_gg"] = options.get_string(option_section, "output_section_pk_gg", "galaxy_galaxy_power_spectrum_pt")

    config["compute_lss_parameters"] = options.get_bool(option_section, "compute_lss_parameters", True)

    config["twopt_type"] = 4
    config["num_ell"] = 3
    config["num_points_use"] = config["window"].shape[0]
    config["num_bands_use"] = config["window"].shape[1]
    config["z_index"] = 4
    config["zm"] = options[option_section, "z_eff"]
    if not isinstance(config["zm"], collections.abc.Iterable):
        config["zm"] = [config["zm"]]
    config["om_fid"] = 0.31
    config["h0_fid"] = 0.7

    config["use_growth"] = False
    config["local_lag_g2"] = True
    config["local_lag_g3"] = False

    # config["derived_parameters"] = np.loadtxt(os.path.join(here, "../output/derived_params.txt"))
    # config["data_parameters"] = np.loadtxt(os.path.join(here, "../output/data_params.txt"))

    # config["z_pk"] = np.loadtxt(os.path.join(here, "../output/z_p_k.txt"))
    # config["log_k_pk"] = np.loadtxt(os.path.join(here, "../output/log_k_h.txt"))
    # config["log_pk"] = np.loadtxt(os.path.join(here, "../output/log_p_k.txt"))

    # if config["verbose"]:
    #     print(config)
    return config

def execute(block, config):
    h = block[names.cosmological_parameters, "h0"]
    omdm = block[names.cosmological_parameters, "omega_c"]
    omb = block[names.cosmological_parameters, "omega_b"]
    omv = block[names.cosmological_parameters, "omega_lambda"]
    omk = block[names.cosmological_parameters, "omega_k"]
    omnuh2 = block[names.cosmological_parameters, "omnuh2"]
    nnu = block.get_double(names.cosmological_parameters, "nnu", default=3.046)
    w = block[names.cosmological_parameters, "w"]
    wa = block[names.cosmological_parameters, "wa"]

    #print("h omdm omb omv omk omnuh2 nnu w wa", h, omdm, omb, omv, omk, omnuh2, nnu, w, wa)
    if config["use_growth"]:
        gamma = block["bias_parameters", f"gamma"]
    else:
        gamma = 1.0
        
    log_Pk = np.log(block[names.matter_power_lin, "p_k"])
    log_k_h = np.log(block[names.matter_power_lin, "k_h"])
    z_Pk = block[names.matter_power_lin, "z"]

    Pk_mm_pt = np.zeros((len(config["zm"]), len(log_k_h)))
    Pk_gm_pt = np.zeros((len(config["zm"]), len(log_k_h)))
    Pk_gg_pt = np.zeros((len(config["zm"]), len(log_k_h)))

    z_growth = block[names.growth_parameters, "z"]
    if not np.allclose(z_Pk, z_growth):
        raise ValueError("Redshifts of power spectrum and growth do not match.")

    sigma8 = block[names.growth_parameters, "SIGMA_8"]
    if block.has_value(names.growth_parameters, "fsigma_8"):
        growth = block[names.growth_parameters, "fsigma_8"]
    else:
        sigma2_vdelta_8 = block[names.growth_parameters, "SIGMA2_VDELTA_8"]
        growth = sigma2_vdelta_8/sigma8
    #print("z_pk", z_Pk, flush=True)
    #print("Pk shape:", log_Pk.shape, "Pk(fixed k):", log_Pk[:,log_Pk.shape[1]//2], flush=True)
    #print("z_growth:", z_growth, flush=True)
    #print("sigma_8", sigma8, flush=True)
    #print("fsigma_8:", growth, flush=True)

    for i, zm in enumerate(config["zm"]):
        b = i + 1

        z = block[names.distances, "z"]
        H_z = block[names.distances, "h"]*2.99792458e8/1e3
        DA_z = block[names.distances, "d_a"]
        z_index = np.argmin(np.abs(z-zm))

        if z_index >= (200-13)//4:
            # z_index too large to fit H_z and DA_z into derived_parameters array
            # so we take a range around z[z_index]
            s = slice(max(0, z_index-4), min(len(z), z_index+4))
            z = z[s]
            H_z = H_z[s]
            DA_z = DA_z[s]
            z_index = np.argmin(np.abs(z-zm))

        # print(f"z_index: {z_index} (z={z[z_index]})")
        # print(f"H(z_index) = {H_z[z_index]}, DA(z_index) = {DA_z[z_index]}")

        # log_Pk_intp = scipy.interpolate.RectBivariateSpline(z, log_k_h, log_Pk)
        # log_k_h = config["log_k_pk"]
        # z = config["z_pk"]
        # log_Pk = log_Pk_intp(z, log_k_h, grid=True)

        # log_Pk = config["log_pk"].T
        # log_k_h = config["log_k_pk"]
        # z = config["z_pk"]

        # sigma8 = scipy.interpolate.InterpolatedUnivariateSpline(z_growth, sigma8)(z)
        # growth = scipy.interpolate.InterpolatedUnivariateSpline(z_growth, growth)(z)
        # z_growth = z

        # Bias parameters
        b1 = block["bias_parameters", f"b1_bin_{b}"]
        b2 = block["bias_parameters", f"b2_bin_{b}"]
        if not config["local_lag_g2"]:
            gamma2 = block["bias_parameters", f"gamma2_bin_{b}"]
        else:
            gamma2 = 1.0
        if not config["local_lag_g3"]:
            gamma3 = block["bias_parameters", f"gamma3_bin_{b}"]
        else:
            gamma3 = 1.0
        a_vir = block["bias_parameters", f"a_vir_bin_{b}"]

        vtheo, vtheo_convolved, Pk_mm, Pk_gm, Pk_gg = \
            config["module"].compute_wedges(
                                h, omdm, omb, omv, omk, omnuh2, nnu, w, wa, 
                                b1, b2, gamma2, gamma3, a_vir, gamma,
                                H_z, DA_z,
                                config["twopt_type"], config["num_ell"], 
                                config["num_points_use"], config["num_bands_use"], 
                                z_index+2, zm, config["om_fid"], config["h0_fid"], 
                                config["use_growth"], config["local_lag_g2"], config["local_lag_g3"],
                                config["window"], 
                                config["bands"], 
                                z_Pk, log_k_h, log_Pk.T,
                                growth, sigma8,
                                verbose=config["verbose"])

        n = len(vtheo)//config["num_ell"]
        vtheo = np.array([vtheo[i*n:(i+1)*n] for i in range(config["num_ell"])])
        n = len(vtheo_convolved)//config["num_ell"]
        vtheo_convolved = np.array([vtheo_convolved[i*n:(i+1)*n] for i in range(config["num_ell"])])

        block[config["output_section_wedges"], f"vtheo_bin_{b}"] = vtheo
        block[config["output_section_wedges"], f"bin_{b}"] = vtheo_convolved
        Pk_mm_pt[i] = Pk_mm
        Pk_gm_pt[i] = Pk_gm
        Pk_gg_pt[i] = Pk_gg

    block[config["output_section_wedges"], "n_wedge"] = config["num_ell"]
    block[config["output_section_wedges"], "bands"] = config["bands"]
    block[config["output_section_wedges"], "z"] = config["zm"]

    block[config["output_section_pk_mm"], "z"] = config["zm"]
    block[config["output_section_pk_gm"], "z"] = config["zm"]
    block[config["output_section_pk_gg"], "z"] = config["zm"]

    block[config["output_section_pk_mm"], "k_h"] = np.exp(log_k_h)
    block[config["output_section_pk_gm"], "k_h"] = np.exp(log_k_h)
    block[config["output_section_pk_gg"], "k_h"] = np.exp(log_k_h)

    block[config["output_section_pk_mm"], "p_k"] = Pk_mm_pt
    block[config["output_section_pk_gm"], "p_k"] = Pk_gm_pt
    block[config["output_section_pk_gg"], "p_k"] = Pk_gg_pt

    if config["compute_lss_parameters"]:
        z = block["distances", "z"]
        D_m = block["distances", "d_m"]
        H = block["distances", "h"]
        r_d = block["distances", "rs_zdrag"]

        for i, zm in enumerate(config["zm"]):
            b = i + 1

            z_index = np.argmin(np.abs(z-zm))
            D_v = ((D_m**2 *z/H)**(1/3))[z_index]
            F_AP = (D_m*H)[z_index]

            z_index = np.argmin(np.abs(z_growth-zm))
            f_sigma_8 = growth[z_index]
            
            block["lss_parameters", f"d_v_over_r_d_bin_{b}"] = D_v/r_d
            block["lss_parameters", f"F_AP_bin_{b}"] = F_AP
            block["lss_parameters", f"f_sigma_8_bin_{b}"] = f_sigma_8

    return 0

def cleanup(config):
    pass
