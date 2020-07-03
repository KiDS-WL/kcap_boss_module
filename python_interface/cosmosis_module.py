import os
import sys
import collections.abc
import multiprocessing
import queue
import signal

import numpy as np
import scipy.interpolate

from cosmosis.datablock import option_section, names

from cosmosis.datablock.cosmosis_py import errors

import bias_module

def setup(options):
    config = {}

    here = os.path.dirname(__file__)
    window_file = options[option_section, "window_file"]
    bands_file = options[option_section, "bands_file"]
    config["window"] = np.loadtxt(window_file)
    config["bands"] = np.loadtxt(bands_file)

    config["verbose"] = options.get_bool(option_section, "verbose", False)
    config["timeout"] = options.get_double(option_section, "timeout", default=-1.0)
    if config["timeout"] <= 0:
        config["timeout"] = None

    try:
        bands_range = options[option_section, "bands_range"]
    except errors.BlockNameNotFound:
        bands_range = [20, 160]
    if len(bands_range) != 2 or not all(np.issubdtype(type(b), np.integer) for b in bands_range):
        raise ValueError(f"bands_range needs to be two integers, got {bands_range}, {[type(b) for b in bands_range]}")

    try:
        points_range = options[option_section, "points_range"]
    except errors.BlockNameNotFound:
        points_range = [4, 32]
    if len(points_range) != 2 or not all(np.issubdtype(type(b), np.integer) for b in points_range):
        raise ValueError(f"points_range needs to be two integers, got {points_range}")
    
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

    config["use_growth"] = options.get_bool(option_section, "use_growth", False)
    config["local_lag_g2"] = options.get_bool(option_section, "local_lag_g2", True)
    config["local_lag_g3"] = options.get_bool(option_section, "local_lag_g3", False)

    config["no_interpolation"] = options.get_bool(option_section, "no_interpolation", False)

    # config["derived_parameters"] = np.loadtxt(os.path.join(here, "../output/derived_params.txt"))
    # config["data_parameters"] = np.loadtxt(os.path.join(here, "../output/data_params.txt"))

    # config["z_pk"] = np.loadtxt(os.path.join(here, "../output/z_p_k.txt"))
    # config["log_k_pk"] = np.loadtxt(os.path.join(here, "../output/log_k_h.txt"))
    # config["log_pk"] = np.loadtxt(os.path.join(here, "../output/log_p_k.txt"))

    # if config["verbose"]:
    #     print(config)
    return config

def run_wedges(q, config, 
               h, omdm, omb, omv, omk, omnuh2, nnu, w, wa, 
               z_Pk, log_k_h, log_Pk,
               growth, sigma8,
               params):

    module = bias_module.BiasModule()
    wedges_config = []

    for z in config["zm"]:
        wedges_config.append(module.initialize_wedges(config["twopt_type"], config["num_ell"], 
                                                         config["num_points_use"], config["num_bands_use"], 
                                                         z, config["om_fid"], config["h0_fid"], 
                                                         config["window"], verbose=config["verbose"]))

    # Setup cosmology
    module.setup_cosmology(h, omdm, omb, omv, omk, omnuh2, nnu, w, wa, 
                                     config["use_growth"], config["local_lag_g2"], config["local_lag_g3"],
                                     z_Pk, log_k_h, log_Pk.T,
                                     growth, sigma8, verbose=config["verbose"])

    results = []
    #print("Running compute_wedges", flush=True)
    for i, zm in enumerate(config["zm"]):
        b = i + 1
        b1, b2, gamma2, gamma3, a_vir, gamma, z_index, H_z, DA_z = params[i]

        vtheo, vtheo_convolved, Pk_mm, Pk_gm, Pk_gg = module.compute_wedges(
                                                        wedges_config[i],
                                                        b1, b2, gamma2, gamma3, a_vir, gamma,
                                                        z_index+2,
                                                        H_z, DA_z,
                                                        config["bands"],
                                                        verbose=config["verbose"])
        results.append((vtheo, vtheo_convolved, Pk_mm, Pk_gm, Pk_gg))
    
    q.put(results)

    for c in wedges_config:
        module.cleanup_wedges(c)
    
    module.cleanup_cosmology()

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

    if config["use_growth"]:
        gamma = block["bias_parameters", f"gamma"]
    else:
        gamma = 1.0
        
    log_Pk = np.log(block[names.matter_power_lin, "p_k"])
    log_k_h = np.log(block[names.matter_power_lin, "k_h"])
    z_Pk = block[names.matter_power_lin, "z"]

    z_growth = block[names.growth_parameters, "z"]
    if not np.allclose(z_Pk, z_growth):
        raise ValueError("Redshifts of power spectrum and growth do not match.")

    sigma8 = block[names.growth_parameters, "SIGMA_8"]
    if block.has_value(names.growth_parameters, "fsigma_8"):
        growth = block[names.growth_parameters, "fsigma_8"]
    else:
        sigma2_vdelta_8 = block[names.growth_parameters, "SIGMA2_VDELTA_8"]
        growth = sigma2_vdelta_8/sigma8

    
    Pk_mm_pt = np.zeros((len(config["zm"]), len(log_k_h)))
    Pk_gm_pt = np.zeros((len(config["zm"]), len(log_k_h)))
    Pk_gg_pt = np.zeros((len(config["zm"]), len(log_k_h)))

    params = []
    for i, zm in enumerate(config["zm"]):
        b = i + 1

        z = block[names.distances, "z"]
        H_z = block[names.distances, "h"]*2.99792458e8/1e3
        DA_z = block[names.distances, "d_a"]

        if config["no_interpolation"]:
            # Original implementation taylored to the way the CosmoMC module reads
            # these parameters. Sensitive to z resolution
            z_index = np.argmin(np.abs(z-zm))

            if z_index >= (200-13)//4:
                # z_index too large to fit H_z and DA_z into derived_parameters array
                # so we take a range around z[z_index]
                s = slice(max(0, z_index-4), min(len(z), z_index+4))
                z = z[s]
                H_z = H_z[s]
                DA_z = DA_z[s]
                z_index = np.argmin(np.abs(z-zm))
        else:
            # More stable and sane.
            H_z_intp = scipy.interpolate.InterpolatedUnivariateSpline(z, H_z, ext=2)
            DA_z_intp = scipy.interpolate.InterpolatedUnivariateSpline(z, DA_z, ext=2)

            H_z = np.atleast_1d(H_z_intp(zm))
            DA_z = np.atleast_1d(DA_z_intp(zm))
            z_index = 0

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

        params.append((b1, b2, gamma2, gamma3, a_vir, gamma, z_index, H_z, DA_z))
        

    # Run wedges
    # Need fork. Using spawn launches the whole cosmosis process
    mp_context = multiprocessing.get_context("fork")

    # Queue to get results back from the process
    q = multiprocessing.Queue()

    # Create the process
    proc = mp_context.Process(target=run_wedges, 
                              args=(q, config, h, omdm, omb, omv, omk, omnuh2, nnu, w, wa, 
                                    z_Pk, log_k_h, log_Pk,
                                    growth, sigma8,
                                    params))
    proc.start()
    try:
        # Wait for results to show up in the queue
        result = q.get(block=True, timeout=config["timeout"])
    except queue.Empty:
        print(f"wedges module timed out after {config['timeout']} s. Attempting to tell process to stop.", file=sys.stderr, flush=True)
        # os.kill(proc.pid, signal.SIGTERM)
        proc.join(0.5)
        if proc.is_alive():
            # print("Process is not cooperating. Killing it.")
            proc.terminate()
        return 1

    proc.join(0.5)
    proc.close()

    # Put the results into the datablock
    for i, zm in enumerate(config["zm"]):
        b = i + 1

        vtheo, vtheo_convolved, Pk_mm, Pk_gm, Pk_gg = result[i]

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
        z = block[names.growth_parameters, "z"]
        fsigma_8 = scipy.interpolate.InterpolatedUnivariateSpline(z, block[names.growth_parameters, "fsigma_8"])
        
        z_background = block[names.distances, "z"]
        F_AP = scipy.interpolate.InterpolatedUnivariateSpline(z_background, block[names.distances, "F_AP"])
        rs_DV = scipy.interpolate.InterpolatedUnivariateSpline(z_background[1:], block[names.distances, "rs_DV"][1:])

        for i, zm in enumerate(config["zm"]):
            b = i + 1

            block["lss_parameters", f"rs_DV_bin_{b}"] = float(rs_DV(zm))
            block["lss_parameters", f"F_AP_bin_{b}"] = float(F_AP(zm))
            block["lss_parameters", f"fsigma_8_bin_{b}"] = float(fsigma_8(zm))

    return 0

def cleanup(config):
    pass


# def execute(block, config):
#     h = block[names.cosmological_parameters, "h0"]
#     omdm = block[names.cosmological_parameters, "omega_c"]
#     omb = block[names.cosmological_parameters, "omega_b"]
#     omv = block[names.cosmological_parameters, "omega_lambda"]
#     omk = block[names.cosmological_parameters, "omega_k"]
#     omnuh2 = block[names.cosmological_parameters, "omnuh2"]
#     nnu = block.get_double(names.cosmological_parameters, "nnu", default=3.046)
#     w = block[names.cosmological_parameters, "w"]
#     wa = block[names.cosmological_parameters, "wa"]

#     if config["use_growth"]:
#         gamma = block["bias_parameters", f"gamma"]
#     else:
#         gamma = 1.0
        
#     log_Pk = np.log(block[names.matter_power_lin, "p_k"])
#     log_k_h = np.log(block[names.matter_power_lin, "k_h"])
#     z_Pk = block[names.matter_power_lin, "z"]

#     z_growth = block[names.growth_parameters, "z"]
#     if not np.allclose(z_Pk, z_growth):
#         raise ValueError("Redshifts of power spectrum and growth do not match.")

#     sigma8 = block[names.growth_parameters, "SIGMA_8"]
#     if block.has_value(names.growth_parameters, "fsigma_8"):
#         growth = block[names.growth_parameters, "fsigma_8"]
#     else:
#         sigma2_vdelta_8 = block[names.growth_parameters, "SIGMA2_VDELTA_8"]
#         growth = sigma2_vdelta_8/sigma8

#     # Setup cosmology
#     config["module"].setup_cosmology(h, omdm, omb, omv, omk, omnuh2, nnu, w, wa, 
#                                      config["use_growth"], config["local_lag_g2"], config["local_lag_g3"],
#                                      z_Pk, log_k_h, log_Pk.T,
#                                      growth, sigma8, verbose=config["verbose"])

    
#     Pk_mm_pt = np.zeros((len(config["zm"]), len(log_k_h)))
#     Pk_gm_pt = np.zeros((len(config["zm"]), len(log_k_h)))
#     Pk_gg_pt = np.zeros((len(config["zm"]), len(log_k_h)))

#     for i, zm in enumerate(config["zm"]):
#         b = i + 1

#         z = block[names.distances, "z"]
#         H_z = block[names.distances, "h"]*2.99792458e8/1e3
#         DA_z = block[names.distances, "d_a"]
#         z_index = np.argmin(np.abs(z-zm))

#         if z_index >= (200-13)//4:
#             # z_index too large to fit H_z and DA_z into derived_parameters array
#             # so we take a range around z[z_index]
#             s = slice(max(0, z_index-4), min(len(z), z_index+4))
#             z = z[s]
#             H_z = H_z[s]
#             DA_z = DA_z[s]
#             z_index = np.argmin(np.abs(z-zm))

#         # Bias parameters
#         b1 = block["bias_parameters", f"b1_bin_{b}"]
#         b2 = block["bias_parameters", f"b2_bin_{b}"]
#         if not config["local_lag_g2"]:
#             gamma2 = block["bias_parameters", f"gamma2_bin_{b}"]
#         else:
#             gamma2 = 1.0
#         if not config["local_lag_g3"]:
#             gamma3 = block["bias_parameters", f"gamma3_bin_{b}"]
#         else:
#             gamma3 = 1.0
#         a_vir = block["bias_parameters", f"a_vir_bin_{b}"]

#         vtheo, vtheo_convolved, Pk_mm, Pk_gm, Pk_gg = config["module"].compute_wedges(
#                                                         config["wedges_config"][i],
#                                                         b1, b2, gamma2, gamma3, a_vir, gamma,
#                                                         z_index+2,
#                                                         H_z, DA_z,
#                                                         config["bands"],
#                                                         verbose=config["verbose"])

#         n = len(vtheo)//config["num_ell"]
#         vtheo = np.array([vtheo[i*n:(i+1)*n] for i in range(config["num_ell"])])
#         n = len(vtheo_convolved)//config["num_ell"]
#         vtheo_convolved = np.array([vtheo_convolved[i*n:(i+1)*n] for i in range(config["num_ell"])])

#         block[config["output_section_wedges"], f"vtheo_bin_{b}"] = vtheo
#         block[config["output_section_wedges"], f"bin_{b}"] = vtheo_convolved
#         Pk_mm_pt[i] = Pk_mm
#         Pk_gm_pt[i] = Pk_gm
#         Pk_gg_pt[i] = Pk_gg

#     block[config["output_section_wedges"], "n_wedge"] = config["num_ell"]
#     block[config["output_section_wedges"], "bands"] = config["bands"]
#     block[config["output_section_wedges"], "z"] = config["zm"]

#     block[config["output_section_pk_mm"], "z"] = config["zm"]
#     block[config["output_section_pk_gm"], "z"] = config["zm"]
#     block[config["output_section_pk_gg"], "z"] = config["zm"]

#     block[config["output_section_pk_mm"], "k_h"] = np.exp(log_k_h)
#     block[config["output_section_pk_gm"], "k_h"] = np.exp(log_k_h)
#     block[config["output_section_pk_gg"], "k_h"] = np.exp(log_k_h)

#     block[config["output_section_pk_mm"], "p_k"] = Pk_mm_pt
#     block[config["output_section_pk_gm"], "p_k"] = Pk_gm_pt
#     block[config["output_section_pk_gg"], "p_k"] = Pk_gg_pt

#     if config["compute_lss_parameters"]:
#         z = block[names.growth_parameters, "z"]
#         fsigma_8 = scipy.interpolate.InterpolatedUnivariateSpline(z, block[names.growth_parameters, "fsigma_8"])
        
#         z_background = block[names.distances, "z"]
#         F_AP = scipy.interpolate.InterpolatedUnivariateSpline(z_background, block[names.distances, "F_AP"])
#         rs_DV = scipy.interpolate.InterpolatedUnivariateSpline(z_background[1:], block[names.distances, "rs_DV"][1:])

#         for i, zm in enumerate(config["zm"]):
#             b = i + 1

#             block["lss_parameters", f"rs_DV_bin_{b}"] = float(rs_DV(zm))
#             block["lss_parameters", f"F_AP_bin_{b}"] = float(F_AP(zm))
#             block["lss_parameters", f"fsigma_8_bin_{b}"] = float(fsigma_8(zm))

#     return 0

# def cleanup(config):
#     for c in config["wedges_config"]:
#         config["module"].cleanup_wedges(c)
    
#     config["module"].cleanup_cosmology()
